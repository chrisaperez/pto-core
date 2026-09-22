// SPDX-License-Identifier: MIT
//
// Core value types for the fragment pipeline.
//
// Everything downstream of the readers speaks `Fragment`, a 16-byte
// coordinate-only record with no owned storage. That is the load-bearing
// decision in this module: a conventional pipeline materialises a fragment BED
// between every stage (bedtools bamtobed | awk | sort | bedtools intersect),
// so an 800M-read experiment is written and re-read as text four times. Here a
// fragment never becomes text unless the user asks for text, and the record is
// small enough that a 64 KiB batch is 4096 fragments -- one L2-resident block
// handed to each consumer.
//
// No htslib types appear in this header, so the public surface stays a plain
// C++ one and the module still compiles with BAM support switched off.
#pragma once

#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace toolkit {

inline constexpr int kCacheLineSize = 64;

// Per-contig coordinates are int32, as in every mainstream genome format's
// on-disk representation (BAM's pos is int32-derived, BED is parsed as such by
// most tools). Human chr1 is 2.5e8, three orders of magnitude below the
// ceiling. Coordinates beyond it are rejected by name at parse time rather
// than truncated -- a truncated coordinate silently places a read on the wrong
// part of the chromosome, which is worse than refusing to run.
using Coord = std::int32_t;

inline constexpr std::int64_t kMaxCoord =
    static_cast<std::int64_t>(std::numeric_limits<Coord>::max());

// Sentinel contig id, used wherever "no contig" has to be representable.
inline constexpr std::int32_t kNoTid = -1;

// Fragment-level flag bits. Distinct from SAM flags on purpose: these are
// facts this module established, not facts it read off the input.
namespace fragflags {
inline constexpr std::uint16_t kReverse = 0x1;    // 5' end is on the minus strand
inline constexpr std::uint16_t kDuplicate = 0x2;  // marked by DuplicateMarker
inline constexpr std::uint16_t kInPeak = 0x4;     // set by the FRiP pass
inline constexpr std::uint16_t kSingleEnd = 0x8;  // inferred span, not a real pair
}  // namespace fragflags

// A sequenced fragment: the half-open span [start, end) on contig `tid`.
//
// Exactly 16 bytes, and asserted so below. For a paired-end BAM this is the
// full TLEN-derived insert, not either read's aligned span -- the fragment is
// the biologically meaningful unit for size distributions, duplicate identity
// and FRiP alike, and reconstructing it once here means the three consumers do
// not each re-derive it from mate fields.
struct Fragment {
    std::int32_t tid = kNoTid;
    Coord start = 0;
    Coord end = 0;
    std::uint8_t mapq = 0;
    std::uint8_t pad = 0;  // explicit, so the layout is not compiler-chosen
    std::uint16_t flags = 0;

    [[nodiscard]] Coord length() const noexcept { return end - start; }
    [[nodiscard]] bool reverse() const noexcept {
        return (flags & fragflags::kReverse) != 0;
    }
    [[nodiscard]] bool duplicate() const noexcept {
        return (flags & fragflags::kDuplicate) != 0;
    }
};

static_assert(sizeof(Fragment) == 16, "Fragment must stay 16 bytes: the batch "
                                      "size and cache-residency argument in "
                                      "this header depend on it");

// Fragments are handed to consumers in batches rather than one at a time. The
// callback is a std::function in the public API, and calling it per fragment
// would put an indirect call in front of every record; per batch the cost is
// amortised over kBatchFragments and the consumer gets a contiguous span it
// can vectorise or hand to OpenMP.
inline constexpr std::size_t kBatchFragments = 4096;  // 64 KiB

// Adds `v` to a 64-bit running total, throwing std::overflow_error instead of
// wrapping (unsigned) or overflowing (signed, which is undefined behaviour).
//
// Reachable, not theoretical: the fragment-BED copy column lets one line emit a
// million fragments, and --max-length 0 admits a 2 Gbp span. 4,612 such lines
// push FRiP's signed base total past INT64_MAX, and 9,224 wrap the size
// histogram's unsigned length sum. Measured 2026-09-11 on an M4 Pro, Release:
// `frip` over 4,700 lines stops with this error after 142 s, `sizes` over
// 9,300 lines after 41 s. Checked once per fragment or once per batch, never
// inside a vectorised or OpenMP loop.
inline void checked_add(std::uint64_t& acc, std::uint64_t v, const char* what) {
    if (v > std::numeric_limits<std::uint64_t>::max() - acc) {
        throw std::overflow_error(std::string(what) + " exceeds the 64-bit range");
    }
    acc += v;
}

inline void checked_add(std::int64_t& acc, std::int64_t v, const char* what) {
    const bool over = v >= 0 ? acc > std::numeric_limits<std::int64_t>::max() - v
                             : acc < std::numeric_limits<std::int64_t>::min() - v;
    if (over) throw std::overflow_error(std::string(what) + " exceeds the 64-bit range");
    acc += v;
}

// Alignment-level filters applied while fragments are being reconstructed, so
// a rejected read never reaches a batch.
struct FragmentFilter {
    int min_mapq = 0;
    // Fragments shorter than this are dropped. The ATAC/CUT&Tag convention of
    // >= 100 bp (drop sub-nucleosomal noise) is expressed as min_length = 100;
    // 0 disables.
    Coord min_length = 0;
    // Fragments longer than this are dropped (0 disables). 1000 bp is the
    // usual gate for chimeric inserts.
    Coord max_length = 1000;
    // SAM bits that disqualify an alignment. Matches the conventional
    // pipeline: no unmapped, secondary, supplementary or vendor-failed reads.
    // Note kDuplicate is NOT excluded by default -- this module marks
    // duplicates itself, and silently honouring an upstream tool's marks would
    // make its own duplicate rate uninterpretable.
    std::uint16_t exclude_flags = 0x4 | 0x100 | 0x200 | 0x800;
    // Keep only concordantly mapped pairs. Set to 0 to include singletons.
    std::uint16_t require_flags = 0x2;
    // Include single-end reads, using the aligned span as the fragment.
    bool allow_single_end = false;
    // Read column 5 of a fragment BED as the number of read pairs supporting
    // the fragment, and emit it that many times -- the 10x fragments.tsv
    // layout (`--with-counts`). OFF by default: in a BED5/BED6 column 5 is a
    // SCORE, and reading it as a count multiplied every statistic by it (a
    // record scored 60 was 60 fragments, exit 0). A standard BED is one fragment
    // per interval. Ignored by the BEDPE and BAM readers.
    bool use_count_column = false;

    void validate() const {
        if (min_length < 0 || max_length < 0) {
            throw std::invalid_argument("fragment length bounds must not be negative");
        }
        if (max_length > 0 && min_length > max_length) {
            throw std::invalid_argument(
                "min fragment length (" + std::to_string(min_length) +
                ") exceeds max (" + std::to_string(max_length) + ")");
        }
    }

    [[nodiscard]] bool length_ok(Coord len) const noexcept {
        if (len <= 0) return false;
        if (len < min_length) return false;
        if (max_length > 0 && len > max_length) return false;
        return true;
    }
};

// Fragment size distribution, accumulated as a histogram rather than a vector
// of lengths.
//
// The distribution is the QC artefact everyone actually wants (the
// nucleosomal ladder), and a histogram answers every question asked of it --
// mean, median, any percentile, the mono/di-nucleosome fractions -- in
// kMaxTrackedLength counters instead of one int per fragment. At 800M
// fragments that is 8 KiB against 3.2 GB, which is the difference between a
// streaming pass and a sort to disk.
class SizeHistogram {
public:
    // Fragments longer than this are counted in `overflow_`, not binned.
    // 2000 bp is past the tri-nucleosome peak; nothing above it informs QC.
    static constexpr Coord kMaxTrackedLength = 2000;

    // Not noexcept: the length sum is checked, see checked_add. `total_` and
    // the bins cannot overflow first -- every added length is >= 1.
    void add(Coord length) {
        if (length <= 0) return;
        checked_add(sum_, static_cast<std::uint64_t>(length), "summed fragment length");
        ++total_;
        if (length > kMaxTrackedLength) {
            ++overflow_;
            return;
        }
        ++counts_[static_cast<std::size_t>(length)];
    }

    // Checked before anything is modified, so a throw leaves *this unchanged.
    void merge(const SizeHistogram& other) {
        std::uint64_t sum = sum_;
        checked_add(sum, other.sum_, "summed fragment length");
        for (std::size_t i = 0; i <= static_cast<std::size_t>(kMaxTrackedLength); ++i) {
            counts_[i] += other.counts_[i];
        }
        total_ += other.total_;
        sum_ = sum;
        overflow_ += other.overflow_;
    }

    [[nodiscard]] std::uint64_t total() const noexcept { return total_; }
    [[nodiscard]] std::uint64_t overflow() const noexcept { return overflow_; }
    [[nodiscard]] std::uint64_t count(Coord length) const noexcept {
        if (length < 0 || length > kMaxTrackedLength) return 0;
        return counts_[static_cast<std::size_t>(length)];
    }
    [[nodiscard]] double mean() const noexcept {
        return total_ == 0 ? 0.0 : static_cast<double>(sum_) / static_cast<double>(total_);
    }

    // Nearest-rank percentile over the tracked range. `q` in [0, 1].
    //
    // Returns kMaxTrackedLength when the rank falls inside the overflow bucket:
    // the exact value is not recoverable from a histogram, so the answer is
    // reported as saturated rather than invented. Callers that care can check
    // overflow() -- it is carried in the report for exactly this reason.
    [[nodiscard]] Coord percentile(double q) const noexcept {
        if (total_ == 0) return 0;
        if (q < 0.0) q = 0.0;
        if (q > 1.0) q = 1.0;
        const auto rank = static_cast<std::uint64_t>(
            q * static_cast<double>(total_ - 1));
        std::uint64_t seen = 0;
        for (std::size_t i = 0; i <= static_cast<std::size_t>(kMaxTrackedLength); ++i) {
            seen += counts_[i];
            if (seen > rank) return static_cast<Coord>(i);
        }
        return kMaxTrackedLength;
    }

    [[nodiscard]] Coord median() const noexcept { return percentile(0.5); }

    // Fraction of fragments in [lo, hi]. The mononucleosome window
    // [180, 247] and the sub-nucleosomal window [0, 100) are the two ratios
    // used to call an ATAC/CUT&Tag library good or bad.
    [[nodiscard]] double fraction_in(Coord lo, Coord hi) const noexcept {
        if (total_ == 0) return 0.0;
        std::uint64_t n = 0;
        for (Coord i = std::max<Coord>(lo, 0);
             i <= std::min(hi, kMaxTrackedLength); ++i) {
            n += counts_[static_cast<std::size_t>(i)];
        }
        return static_cast<double>(n) / static_cast<double>(total_);
    }

private:
    // +1 so kMaxTrackedLength itself is a valid index.
    alignas(kCacheLineSize) std::uint64_t
        counts_[static_cast<std::size_t>(kMaxTrackedLength) + 1] = {};
    std::uint64_t total_ = 0;
    std::uint64_t sum_ = 0;
    std::uint64_t overflow_ = 0;
};

// Why records were discarded. Reported rather than summed into one "skipped"
// number, because "90% dropped" means something different for each reason and
// a single counter has repeatedly hidden a chromosome-naming mismatch.
struct StreamStats {
    std::uint64_t records_read = 0;      // input records (reads or BEDPE lines)
    std::uint64_t fragments_emitted = 0;
    std::uint64_t dropped_flags = 0;     // failed exclude/require flag test
    std::uint64_t dropped_mapq = 0;
    std::uint64_t dropped_length = 0;    // outside [min_length, max_length]
    std::uint64_t dropped_unpaired = 0;  // single-end and not allowed
    std::uint64_t dropped_interchrom = 0;// mates on different contigs
    std::uint64_t malformed_lines = 0;   // text input only
    SizeHistogram sizes;

    void merge(const StreamStats& o) {
        records_read += o.records_read;
        fragments_emitted += o.fragments_emitted;
        dropped_flags += o.dropped_flags;
        dropped_mapq += o.dropped_mapq;
        dropped_length += o.dropped_length;
        dropped_unpaired += o.dropped_unpaired;
        dropped_interchrom += o.dropped_interchrom;
        malformed_lines += o.malformed_lines;
        sizes.merge(o.sizes);
    }
};

// A half-open genomic interval, used only at the parsing boundary. The hot
// paths use Fragment and the contig dictionary below; this type owns a string
// and is not allowed anywhere per-record.
struct Region {
    std::string chrom;
    std::int64_t start = 0;
    std::int64_t end = 0;
    std::string name;

    [[nodiscard]] std::int64_t length() const noexcept { return end - start; }
};

}  // namespace toolkit
