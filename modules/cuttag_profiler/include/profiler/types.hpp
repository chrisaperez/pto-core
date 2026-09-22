// SPDX-License-Identifier: MIT
//
// Core value types shared by the BAM reader, the signal calculator and the
// embedded HTTP layer. Deliberately free of any htslib include so that the
// public surface of the library does not leak C headers into consumers.
#pragma once

#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace profiler {

inline constexpr int kCacheLineSize = 64;

enum class Strand : char { kForward = '+', kReverse = '-', kUnknown = '.' };

// A half-open genomic interval [start, end) on `chrom`, 0-based (BED
// convention). `name` and `strand` are optional BED columns 4 and 6.
struct Region {
    std::string chrom;
    std::int64_t start = 0;
    std::int64_t end = 0;
    std::string name;
    Strand strand = Strand::kUnknown;

    [[nodiscard]] std::int64_t length() const noexcept { return end - start; }
};

// Where the fixed-width window is anchored, mirroring the vocabulary of
// `deeptools computeMatrix reference-point`.
enum class ReferencePoint { kTSS, kCenter, kTES };

// Per-bin signal normalisation. Definitions are documented in
// docs/APPLICATION_NOTE.tex and in README.md; `kRaw` leaves the counts
// untouched so downstream tools can apply their own scaling.
enum class Normalization { kRaw, kCPM, kRPKM, kBPM };

// How a bin's signal is accumulated from the alignments overlapping it.
//   kReads - number of alignments (or fragments) overlapping the bin. This is
//            what deeptools reports and is the default for parity.
//   kDepth - mean per-base sequencing depth, honouring the CIGAR string so
//            that deletions and reference skips (introns) do not contribute.
enum class CountMode { kReads, kDepth };

// SAM flag bits, redeclared here to keep htslib out of the public headers.
// Values are fixed by the SAM specification.
namespace flags {
inline constexpr std::uint16_t kPaired = 0x1;
inline constexpr std::uint16_t kProperPair = 0x2;
inline constexpr std::uint16_t kUnmapped = 0x4;
inline constexpr std::uint16_t kMateUnmapped = 0x8;
inline constexpr std::uint16_t kReverse = 0x10;
inline constexpr std::uint16_t kRead1 = 0x40;
inline constexpr std::uint16_t kRead2 = 0x80;
inline constexpr std::uint16_t kSecondary = 0x100;
inline constexpr std::uint16_t kQCFail = 0x200;
inline constexpr std::uint16_t kDuplicate = 0x400;
inline constexpr std::uint16_t kSupplementary = 0x800;
}  // namespace flags

// Narrows a parsed integer option to int after checking it against [lo, hi].
// Both entry points parse integers wider than int -- the CLI as int64 via
// from_chars, the HTTP handler as JSON numbers -- and both used to static_cast
// straight to int, so `--min-mapq 4294967296` became 0 and `--threads
// 4294967297` became 1 before any validation could see the real value.
// `hi` must itself fit in int.
[[nodiscard]] inline int narrow_option(std::int64_t value, std::int64_t lo, std::int64_t hi,
                                       const char* name) {
    if (value < lo || value > hi) {
        throw std::invalid_argument(std::string(name) + " must be between " +
                                    std::to_string(lo) + " and " + std::to_string(hi) +
                                    ", got " + std::to_string(value));
    }
    return static_cast<int>(value);
}

// Alignment-level filters applied before an alignment contributes signal.
struct FilterOptions {
    int min_mapq = 0;
    // Alignments with any of these bits set are discarded. The default matches
    // the conventional CUT&Tag pipeline: no unmapped, secondary, supplementary,
    // vendor-failed or optically duplicated reads.
    std::uint16_t exclude_flags = flags::kUnmapped | flags::kSecondary |
                                  flags::kQCFail | flags::kDuplicate |
                                  flags::kSupplementary;
    // Alignments must have all of these bits set. Set to `kProperPair` to keep
    // only concordantly mapped pairs.
    std::uint16_t require_flags = 0;

    // Treat a properly paired alignment as the full sequenced fragment
    // (TLEN-derived) rather than the aligned read span. Standard practice for
    // CUT&Tag, where the fragment marks the protected footprint.
    bool extend_to_fragment = true;
    // Fragments longer than this are dropped (0 disables). CUT&Tag protocols
    // typically gate at 1000 bp to exclude chimeric fragments.
    std::int64_t max_fragment_length = 1000;
    // For single-end data, extend each read to this many bases (0 disables).
    std::int64_t extend_reads_to = 0;

    [[nodiscard]] bool passes_flags(std::uint16_t flag) const noexcept {
        return (flag & exclude_flags) == 0 &&
               (flag & require_flags) == require_flags;
    }

    // Largest fragment length or read extension accepted, in bp.
    // alignment_span computes `pos + extend_reads_to`, which overflows for a
    // value near INT64_MAX; 2^40 bp is three orders of magnitude past any
    // assembled chromosome.
    static constexpr std::int64_t kMaxExtension = std::int64_t{1} << 40;

    // One validation point for every entry path, as ProfileOptions::validate
    // is. Before 2026-09-11 nothing checked these, and each wrong value changed
    // the answer rather than failing: --min-mapq 300 filtered out every
    // alignment (MAPQ is 0-255), --max-fragment -5 silently turned the
    // fragment gate off, and --extend-reads -100 was accepted.
    void validate() const {
        if (min_mapq < 0 || min_mapq > 255) {
            throw std::invalid_argument("minimum MAPQ must be between 0 and 255, got " +
                                        std::to_string(min_mapq));
        }
        if (max_fragment_length < 0 || max_fragment_length > kMaxExtension) {
            throw std::invalid_argument(
                "maximum fragment length must be between 0 (off) and " +
                std::to_string(kMaxExtension) + ", got " +
                std::to_string(max_fragment_length));
        }
        if (extend_reads_to < 0 || extend_reads_to > kMaxExtension) {
            throw std::invalid_argument("read extension must be between 0 (off) and " +
                                        std::to_string(kMaxExtension) + ", got " +
                                        std::to_string(extend_reads_to));
        }
    }
};

// Parameters describing the matrix to compute.
struct ProfileOptions {
    ReferencePoint reference_point = ReferencePoint::kTSS;
    std::int64_t upstream = 2000;
    std::int64_t downstream = 2000;
    std::int64_t bin_size = 50;
    Normalization normalization = Normalization::kCPM;
    CountMode count_mode = CountMode::kReads;
    // Flip the profile for regions on the minus strand so that "upstream"
    // always means 5' of the feature.
    bool respect_strand = true;
    // Value substituted where a window runs off the end of a chromosome.
    double missing_value = 0.0;
    int threads = 0;  // 0 => hardware_concurrency()

    // Largest matrix width we will allocate, per region. The regions x bins
    // matrix of double is allocated up front, so an unbounded bin count is an
    // OOM-kill with no diagnostic (AUDIT.md S6).
    static constexpr std::int64_t kMaxBins = 10'000'000;

    // Largest window accepted, in bp. Every bin boundary is computed as
    // `width * b / nbins` (compute_matrix; BamReader::query_bins and add_bases)
    // with b <= nbins <= kMaxBins, so a width of at most INT64_MAX / kMaxBins
    // keeps every such product in range. The bin-count ceiling alone did not:
    // upstream 4e18, downstream 1, bin size 1e12 is 4e6 bins, passed
    // validate(), and UBSan then reported "signed integer overflow:
    // 4000000000000000001 * 3" at signal_calc.cpp:196 -- in Release it exited 0
    // writing bin offsets 2.3e12 bp wrong (2026-09-11). The cap is ~922 Gbp,
    // still far past any chromosome.
    static constexpr std::int64_t kMaxWindow =
        std::numeric_limits<std::int64_t>::max() / kMaxBins;

    // compute_matrix clamps its worker count to 256 anyway; this bounds what a
    // caller may ask for, and refuses a negative count rather than reading it
    // as "use every core".
    static constexpr int kMaxThreads = 4096;

    // Checked: `upstream + downstream` on int64 overflows for coordinates near
    // the type limit, and the overflow happens *before* any caller's `> 0` test
    // can see it (AUDIT.md S5, same class). Saturating to INT64_MAX keeps the
    // function total, and validate() rejects the saturated value below.
    [[nodiscard]] std::int64_t window_length() const noexcept {
        std::int64_t w = 0;
        if (__builtin_add_overflow(upstream, downstream, &w)) {
            return std::numeric_limits<std::int64_t>::max();
        }
        return w;
    }

    // int64 throughout: the previous `static_cast<int>` narrowed a 4e9 window at
    // bin size 1 into a negative int, which then failed a downstream `<= 0`
    // check with a misleading "window is empty" message.
    //
    // The round-up is CHECKED, for the same reason window_length() above is.
    // `w + bin_size - 1` overflows whenever w is near the type limit, which is
    // reachable precisely because window_length() saturates to INT64_MAX on a
    // hostile `{"upstream": 9223372036854775807, "downstream": 1}`. That is UB,
    // and it did not merely produce a wrong number -- it produced a large
    // NEGATIVE one, so `bins > kMaxBins` was false and validate()'s 10M-bin
    // ceiling was bypassed completely, re-opening the unauthenticated
    // memory-exhaustion path that ceiling exists to close.
    //
    // Measured before this fix, with upstream=INT64_MAX, downstream=1,
    // bin_size=50: bin_count() = -184467440737095515 and validate() did not
    // throw, on BOTH Apple Clang -O2 and GCC 15 -O2; UBSan reported
    // "signed integer overflow: 9223372036854775807 + 50" here. Being UB, it
    // was also unstable across optimisation levels, which is why the Clang
    // Release build's test run passed while the GCC one failed -- fixing the
    // overflow in window_length() alone left this sibling live, exactly the
    // incomplete-propagation pattern REVIEW_2026-08-15 is about.
    //
    // Saturating keeps the function total and noexcept; validate() then rejects
    // the saturated value against kMaxBins.
    [[nodiscard]] std::int64_t bin_count() const noexcept {
        const std::int64_t w = window_length();
        if (bin_size <= 0 || w <= 0) return 0;
        std::int64_t numerator = 0;
        // bin_size >= 1 here, so `bin_size - 1` cannot itself overflow.
        if (__builtin_add_overflow(w, bin_size - 1, &numerator)) {
            return std::numeric_limits<std::int64_t>::max();
        }
        return numerator / bin_size;
    }

    // Single validation point for EVERY entry path -- CLI, HTTP, and any future
    // caller. These checks previously lived only in main.cpp, so the
    // POST /api/profile handler, which builds this same struct from untrusted
    // JSON, reached compute_matrix with only a subset of them applied
    // (REVIEW_2026-08-15 finding 2).
    void validate() const {
        if (bin_size <= 0) {
            throw std::invalid_argument("bin size must be greater than zero");
        }
        if (upstream < 0 || downstream < 0) {
            // window_length() alone does not catch this: upstream -100 with
            // downstream 2000 still sums to a positive 1900, and the resulting
            // profile is silently mislabelled on its x-axis.
            throw std::invalid_argument(
                "upstream and downstream must not be negative");
        }
        if (window_length() <= 0) {
            throw std::invalid_argument(
                "upstream + downstream must be greater than zero");
        }
        if (window_length() > kMaxWindow) {
            throw std::invalid_argument(
                "window of " + std::to_string(window_length()) + " bp exceeds the " +
                std::to_string(kMaxWindow) + " bp limit");
        }
        if (threads < 0 || threads > kMaxThreads) {
            throw std::invalid_argument("threads must be between 0 (all cores) and " +
                                        std::to_string(kMaxThreads) + ", got " +
                                        std::to_string(threads));
        }
        const std::int64_t bins = bin_count();
        if (bins > kMaxBins) {
            throw std::invalid_argument(
                "window of " + std::to_string(window_length()) +
                " bp at bin size " + std::to_string(bin_size) + " needs " +
                std::to_string(bins) + " bins (limit " +
                std::to_string(kMaxBins) +
                "); increase the bin size or narrow the window");
        }
    }
};

// A dense regions x bins signal matrix in row-major order.
struct SignalMatrix {
    std::size_t rows = 0;
    std::size_t cols = 0;
    std::vector<double> values;         // rows * cols
    std::vector<std::string> row_names; // size == rows
    // Column-wise mean across all regions - the "meta-profile" curve.
    std::vector<double> column_mean;
    // Bin centre offsets relative to the reference point, in base pairs.
    std::vector<std::int64_t> bin_offsets;

    // Diagnostics carried through to the API response.
    std::uint64_t total_mapped_reads = 0;
    std::uint64_t regions_skipped = 0;
    double scale_factor = 1.0;
    double elapsed_seconds = 0.0;

    [[nodiscard]] double at(std::size_t r, std::size_t c) const {
        return values[r * cols + c];
    }
};

}  // namespace profiler
