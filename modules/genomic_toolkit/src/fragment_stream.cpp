// SPDX-License-Identifier: MIT

#include "toolkit/fragment_stream.hpp"

#include <algorithm>
#include <cctype>
#include <charconv>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include "line_reader.hpp"

#ifdef GTK_HAVE_HTSLIB
#include <htslib/bgzf.h>
#include <htslib/hts.h>
#include <htslib/sam.h>
#endif

namespace toolkit {
namespace {

// Splits `line` on tabs into at most `max_fields` views. Returns how many were
// filled. std::from_chars and string_view throughout: no substr, no stoi, no
// allocation. On a 200M-line fragment file that difference is the whole
// runtime.
std::size_t split_tabs(std::string_view line, std::string_view* out,
                       std::size_t max_fields) {
    std::size_t n = 0;
    std::size_t start = 0;
    while (n < max_fields) {
        const std::size_t tab = line.find('\t', start);
        if (tab == std::string_view::npos) {
            out[n++] = line.substr(start);
            break;
        }
        out[n++] = line.substr(start, tab - start);
        start = tab + 1;
    }
    return n;
}

// Parses a whole field as a base-10 integer. A leading '-' IS accepted, so
// every caller range-checks the result. Returns false on anything else,
// including a partially numeric field ("100abc") -- from_chars alone accepts
// that prefix, and accepting it silently is how a mis-delimited file becomes a
// plausible wrong answer instead of an error.
bool parse_int(std::string_view s, std::int64_t& out) {
    if (s.empty()) return false;
    const char* first = s.data();
    const char* last = s.data() + s.size();
    const auto res = std::from_chars(first, last, out);
    return res.ec == std::errc() && res.ptr == last;
}

bool is_comment(std::string_view line) {
    return line.empty() || line.front() == '#' || line.starts_with("track") ||
           line.starts_with("browser");
}

// Accumulates fragments and flushes to the callback every kBatchFragments.
class BatchEmitter {
public:
    BatchEmitter(const BatchCallback& cb, StreamStats& stats) : cb_(cb), stats_(stats) {
        buffer_.reserve(kBatchFragments);
    }

    void emit(const Fragment& f) {
        stats_.sizes.add(f.length());
        ++stats_.fragments_emitted;
        buffer_.push_back(f);
        if (buffer_.size() == kBatchFragments) flush();
    }

    void flush() {
        if (buffer_.empty()) return;
        cb_(std::span<Fragment>(buffer_.data(), buffer_.size()));
        buffer_.clear();
    }

private:
    const BatchCallback& cb_;
    StreamStats& stats_;
    std::vector<Fragment> buffer_;
};

// Shared coordinate validation for the text paths. Rejecting an out-of-range
// coordinate rather than truncating it into int32 is the same rule PeakSet
// applies, for the same reason.
bool coords_ok(std::int64_t start, std::int64_t end) {
    return start >= 0 && end >= start && end <= kMaxCoord;
}

bool ends_with_ci(const std::string& s, std::string_view suffix) {
    if (s.size() < suffix.size()) return false;
    return std::equal(suffix.rbegin(), suffix.rend(), s.rbegin(),
                      [](char a, char b) { return std::tolower(a) == std::tolower(b); });
}

}  // namespace

StreamStats stream_bedpe(const std::string& path, const FragmentFilter& filter,
                         ContigDict& dict, const BatchCallback& on_batch) {
    filter.validate();
    StreamStats stats;
    detail::LineReader reader(path);
    BatchEmitter emitter(on_batch, stats);

    std::string_view fields[10];
    std::string_view line;
    while (reader.next(line)) {
        if (is_comment(line)) continue;
        ++stats.records_read;

        const std::size_t n = split_tabs(line, fields, 10);
        if (n < 6) {
            ++stats.malformed_lines;
            continue;
        }

        std::int64_t s1 = 0, e1 = 0, s2 = 0, e2 = 0;
        if (!parse_int(fields[1], s1) || !parse_int(fields[2], e1) ||
            !parse_int(fields[4], s2) || !parse_int(fields[5], e2)) {
            ++stats.malformed_lines;
            continue;
        }
        // BEDPE uses -1 in every coordinate column for an unmapped mate.
        if (s1 < 0 || s2 < 0) {
            ++stats.dropped_unpaired;
            continue;
        }
        // A mate whose own interval ends before it starts is corrupt. The
        // min/max below would otherwise build a plausible span out of it.
        if (e1 < s1 || e2 < s2) {
            ++stats.malformed_lines;
            continue;
        }
        if (fields[0] != fields[3]) {
            ++stats.dropped_interchrom;
            continue;
        }

        const std::int64_t start = std::min(s1, s2);
        const std::int64_t end = std::max(e1, e2);
        if (!coords_ok(start, end)) {
            ++stats.malformed_lines;
            continue;
        }

        // Column 8 is optional. An ABSENT column means "this file does not
        // record MAPQ", and 255 (the SAM "unavailable" value) lets the record
        // through -- there is nothing to filter on. A PRESENT but unparseable
        // column is a different thing entirely: the record is corrupt, and the
        // line is rejected.
        //
        // This previously discarded parse_int's result, which failed open in
        // two ways, both reproduced against a hand-written BEDPE:
        //
        //   "...\tname\tXX"     -> mapq stayed 255 and the fragment passed
        //                          --min-mapq 30 as if it were a perfect
        //                          alignment.
        //   "...\tname\t60abc"  -> std::from_chars writes the parsed prefix
        //                          before reporting the trailing junk, so the
        //                          discarded `false` left mapq = 60.
        //
        // Silently promoting a corrupt field to maximum confidence inverts the
        // filter on exactly the input it exists to catch. The `copies` column
        // in stream_fragment_bed() below already checks its parse; this is the
        // same check, in the sibling function that was missed.
        std::int64_t mapq = 255;
        if (n >= 8 && fields[7] != ".") {
            // Outside 0-255 is not a MAPQ. "-7" used to be counted as a
            // low-quality alignment (dropped_mapq), which reports a corrupt
            // record as a filtering decision.
            if (!parse_int(fields[7], mapq) || mapq < 0 || mapq > 255) {
                ++stats.malformed_lines;
                continue;
            }
        }
        if (mapq < filter.min_mapq) {
            ++stats.dropped_mapq;
            continue;
        }

        const auto len = static_cast<Coord>(end - start);
        if (!filter.length_ok(len)) {
            ++stats.dropped_length;
            continue;
        }

        Fragment f;
        f.tid = dict.intern(fields[0]);
        f.start = static_cast<Coord>(start);
        f.end = static_cast<Coord>(end);
        f.mapq = static_cast<std::uint8_t>(std::clamp<std::int64_t>(mapq, 0, 255));
        if (n >= 9 && fields[8] == "-") f.flags |= fragflags::kReverse;
        emitter.emit(f);
    }

    emitter.flush();
    return stats;
}

StreamStats stream_fragment_bed(const std::string& path, const FragmentFilter& filter,
                                ContigDict& dict, const BatchCallback& on_batch) {
    filter.validate();
    StreamStats stats;
    detail::LineReader reader(path);
    BatchEmitter emitter(on_batch, stats);

    std::string_view fields[5];
    std::string_view line;
    while (reader.next(line)) {
        if (is_comment(line)) continue;
        ++stats.records_read;

        const std::size_t n = split_tabs(line, fields, 5);
        std::int64_t start = 0, end = 0;
        if (n < 3 || !parse_int(fields[1], start) || !parse_int(fields[2], end) ||
            !coords_ok(start, end)) {
            ++stats.malformed_lines;
            continue;
        }

        const auto len = static_cast<Coord>(end - start);
        if (!filter.length_ok(len)) {
            ++stats.dropped_length;
            continue;
        }

        // Column 5 of a 10x fragments.tsv is the number of read pairs
        // supporting the fragment. Emitting one fragment per supporting pair
        // keeps counts comparable with a BAM of the same library; treating the
        // line as a single fragment would deflate every count by the mean
        // duplication rate.
        //
        // An ABSENT column means one copy. A PRESENT column must be a count in
        // [1, kMaxCopiesPerLine], or the record is malformed -- the same
        // present-but-corrupt rule as the BEDPE MAPQ column above. This used to
        // fail open: "abc", "0" and "-3" each became 1 copy, and 99999999999999
        // was clamped to a million, so one corrupt line added a million
        // fragments to every statistic with no signal. The ceiling bounds what
        // one line can emit; a real per-fragment read count is orders of
        // magnitude below it.
        static constexpr std::int64_t kMaxCopiesPerLine = 1'000'000;
        //
        // All of that applies only when the caller says column 5 IS a count
        // (FragmentFilter::use_count_column, `--with-counts`). By default it is
        // not read at all: in a BED5/BED6 it is a score, and a record scored 60
        // used to be emitted 60 times.
        std::int64_t copies = 1;
        if (filter.use_count_column && n >= 5 &&
            (!parse_int(fields[4], copies) || copies < 1 || copies > kMaxCopiesPerLine)) {
            ++stats.malformed_lines;
            continue;
        }

        Fragment f;
        f.tid = dict.intern(fields[0]);
        f.start = static_cast<Coord>(start);
        f.end = static_cast<Coord>(end);
        f.mapq = 255;
        for (std::int64_t c = 0; c < copies; ++c) emitter.emit(f);
    }

    emitter.flush();
    return stats;
}

#ifdef GTK_HAVE_HTSLIB
StreamStats stream_bam(const std::string& path, const FragmentFilter& filter,
                       ContigDict& dict, const BatchCallback& on_batch,
                       const std::vector<std::string>& regions) {
    filter.validate();
    StreamStats stats;

    htsFile* fp = hts_open(path.c_str(), "r");
    if (fp == nullptr) throw std::runtime_error("cannot open '" + path + "'");
    sam_hdr_t* hdr = sam_hdr_read(fp);
    if (hdr == nullptr) {
        hts_close(fp);
        throw std::runtime_error("cannot read a SAM/BAM header from '" + path + "'");
    }
    // A BGZF stream ends in a 28-byte empty block. Without it the file was cut
    // short, and cut at a block boundary it decodes cleanly into a prefix: on
    // ENCFF121ZQX's chr22 slice, 66,296 of 132,934 records, exit 0, with only
    // htslib's stderr warning to say so -- an interrupted copy reported as
    // half a library. fastq_stream and pto-peaks already refuse this
    // (docs/TORTURE_2026-09-10.md). Checked before any record is read, so no
    // partial report is printed. The format's compression is asked rather than
    // `is_bgzf`, which htslib documents as not implying BGZF. A pipe cannot be
    // checked (2) and is let through; CRAM carries its own EOF container.
    if (hts_get_format(fp)->compression == bgzf) {
        const int eof = bgzf_check_EOF(fp->fp.bgzf);
        if (eof == 0 || eof < 0) {
            sam_hdr_destroy(hdr);
            hts_close(fp);
            throw std::runtime_error(
                eof == 0 ? "'" + path + "' has no BGZF end-of-file marker, so it is "
                                        "truncated; statistics from it would describe only "
                                        "the part that survived"
                         : "cannot check '" + path + "' for truncation (I/O error)");
        }
    }
    // A Fragment's tid is the DICTIONARY's id, never the raw BAM tid. The two
    // coincide only when `dict` arrives empty, and `frip` deliberately
    // pre-seeds it with the peak file's contigs first -- so on a multi-contig
    // header chr22's reads (BAM tid 21) were looked up as whatever contig held
    // dictionary id 21, and FRiP was silently wrong with exit 0 whenever the
    // peak file did not list contigs in header order. Single-contig fixtures
    // cannot show it. `validate_bam_contig_order.sh` pins it.
    std::vector<std::int32_t> dict_id_of_tid(static_cast<std::size_t>(sam_hdr_nref(hdr)));
    for (int i = 0; i < sam_hdr_nref(hdr); ++i) {
        dict_id_of_tid[static_cast<std::size_t>(i)] = dict.intern(sam_hdr_tid2name(hdr, i));
    }

    // Region setup. Everything here is torn down on every throw path below,
    // which is why it is raw rather than a guard type: the two htslib handles
    // have distinct destroy functions and no shared base.
    hts_idx_t* idx = nullptr;
    hts_itr_t* itr = nullptr;
    std::int64_t region_beg = 0;  // 0-based inclusive; ownership boundary
    std::size_t next_region = 0;

    if (!regions.empty()) {
        idx = sam_index_load(fp, path.c_str());
        if (idx == nullptr) {
            sam_hdr_destroy(hdr);
            hts_close(fp);
            throw std::runtime_error(
                "--region needs an index for '" + path +
                "'; run `samtools index` on it. Without one every shard would "
                "read the whole file.");
        }
        itr = sam_itr_querys(idx, hdr, regions[0].c_str());
        if (itr == nullptr) {
            hts_idx_destroy(idx);
            sam_hdr_destroy(hdr);
            hts_close(fp);
            throw std::runtime_error("cannot parse region '" + regions[0] +
                                     "' against this header");
        }
        region_beg = itr->beg;
        next_region = 1;
    }

    BatchEmitter emitter(on_batch, stats);
    bam1_t* rec = bam_init1();
    int ret = 0;

    while (true) {
        ret = itr != nullptr ? sam_itr_next(fp, itr, rec) : sam_read1(fp, hdr, rec);
        if (ret < 0) {
            // End of one region is not end of input: advance to the next.
            // `ret` is reset to 0 so the truncation check below sees a clean
            // exhaustion rather than the last iterator's -1.
            if (itr != nullptr && ret == -1 && next_region < regions.size()) {
                hts_itr_destroy(itr);
                itr = sam_itr_querys(idx, hdr, regions[next_region].c_str());
                if (itr == nullptr) {
                    bam_destroy1(rec);
                    hts_idx_destroy(idx);
                    sam_hdr_destroy(hdr);
                    hts_close(fp);
                    throw std::runtime_error("cannot parse region '" +
                                             regions[next_region] +
                                             "' against this header");
                }
                region_beg = itr->beg;
                ++next_region;
                continue;
            }
            break;
        }
        // A record starting before the region belongs to the shard on the
        // other side of the boundary. Skipped BEFORE records_read so the
        // per-shard counters sum to the whole-file counters -- counting it
        // here would make `records_read` overlap between shards.
        if (itr != nullptr && rec->core.pos < region_beg) continue;

        ++stats.records_read;
        const std::uint16_t flag = rec->core.flag;

        if ((flag & filter.exclude_flags) != 0) {
            ++stats.dropped_flags;
            continue;
        }
        if (rec->core.qual < filter.min_mapq) {
            ++stats.dropped_mapq;
            continue;
        }

        std::int64_t start = 0;
        std::int64_t end = 0;
        std::uint16_t frag_flags = 0;

        if ((flag & 0x1) != 0) {  // paired
            if ((flag & filter.require_flags) != filter.require_flags) {
                ++stats.dropped_flags;
                continue;
            }
            if (rec->core.tid != rec->core.mtid) {
                ++stats.dropped_interchrom;
                continue;
            }
            // Exactly one read of a proper pair has a positive TLEN, so
            // emitting only on that read yields one fragment per pair with no
            // mate lookup and no read-name hash table. A pair with TLEN 0 or a
            // mate-before-self orientation is handled by its partner.
            if (rec->core.isize <= 0) continue;
            start = rec->core.pos;
            end = start + rec->core.isize;
            if ((flag & 0x10) != 0) frag_flags |= fragflags::kReverse;
        } else {
            if (!filter.allow_single_end) {
                ++stats.dropped_unpaired;
                continue;
            }
            start = rec->core.pos;
            end = bam_endpos(rec);
            frag_flags |= fragflags::kSingleEnd;
            if ((flag & 0x10) != 0) frag_flags |= fragflags::kReverse;
        }

        if (!coords_ok(start, end)) {
            ++stats.malformed_lines;
            continue;
        }
        const auto len = static_cast<Coord>(end - start);
        if (!filter.length_ok(len)) {
            ++stats.dropped_length;
            continue;
        }

        // An unplaced record (tid -1) or one naming a contig past the header
        // is malformed, not tid 0: fail closed rather than index out of range.
        if (rec->core.tid < 0 ||
            static_cast<std::size_t>(rec->core.tid) >= dict_id_of_tid.size()) {
            ++stats.malformed_lines;
            continue;
        }

        Fragment f;
        f.tid = dict_id_of_tid[static_cast<std::size_t>(rec->core.tid)];
        f.start = static_cast<Coord>(start);
        f.end = static_cast<Coord>(end);
        f.mapq = rec->core.qual;
        f.flags = frag_flags;
        emitter.emit(f);
    }

    emitter.flush();
    bam_destroy1(rec);
    if (itr != nullptr) hts_itr_destroy(itr);
    if (idx != nullptr) hts_idx_destroy(idx);
    sam_hdr_destroy(hdr);
    hts_close(fp);

    // -1 is a clean EOF; anything below it is a truncated or corrupt file, and
    // returning the partial statistics as if they were complete is how a
    // truncated BAM becomes a published number.
    if (ret < -1) throw std::runtime_error("truncated or corrupt record in '" + path + "'");
    return stats;
}
#endif

StreamStats stream_fragments(const std::string& path, const FragmentFilter& filter,
                             ContigDict& dict, const BatchCallback& on_batch,
                             const std::vector<std::string>& regions) {
    std::string stem = path;
    if (ends_with_ci(stem, ".gz")) stem.resize(stem.size() - 3);

    const bool indexable = ends_with_ci(stem, ".bam") || ends_with_ci(stem, ".cram") ||
                           ends_with_ci(stem, ".sam");
    // Fail closed. Accepting --region on a text format and silently ignoring
    // it would give every shard of a map-reduce the WHOLE file: N times the
    // work, and an answer that is N times too large rather than an error.
    if (!regions.empty() && !indexable) {
        throw std::runtime_error(
            "--region needs an indexed BAM/CRAM; '" + path +
            "' is a text fragment file, which has no index to seek with. "
            "Split text inputs by byte range instead.");
    }

    if (ends_with_ci(stem, ".bedpe")) return stream_bedpe(path, filter, dict, on_batch);
    if (indexable) {
#ifdef GTK_HAVE_HTSLIB
        return stream_bam(path, filter, dict, on_batch, regions);
#else
        throw std::runtime_error(
            "'" + path +
            "' needs BAM/CRAM support, which this build does not have. "
            "Reconfigure with -DGTK_WITH_HTSLIB=ON, or convert the input to a "
            "BEDPE/fragment BED first.");
#endif
    }
    // .bed, .tsv, .fragments and anything unrecognised: the fragment-BED
    // parser is the permissive one, and its first malformed line is reported
    // rather than guessed at.
    return stream_fragment_bed(path, filter, dict, on_batch);
}

}  // namespace toolkit
