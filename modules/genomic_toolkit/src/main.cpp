// SPDX-License-Identifier: MIT
//
// genomic_toolkit CLI.
//
// Three subcommands over one streaming pass each:
//   sizes   -- fragment size distribution and QC ratios
//   markdup -- duplicate rate and library-complexity estimate
//   frip    -- fraction of reads in peaks
//
// `frip` runs duplicate marking in the same pass when asked, which is the
// point of the module: the shell version of that combination reads the data
// three times and writes two intermediate files.

#include <cerrno>
#include <cstdio>
#include <cstring>
#include <exception>
#include <string>
#include <vector>

#include "toolkit/cli_args.hpp"
#include "toolkit/contig_dict.hpp"
#include "toolkit/dup_marker.hpp"
#include "toolkit/fragment_stream.hpp"
#include "toolkit/frip.hpp"
#include "toolkit/file_identity.hpp"
#include "toolkit/peak_set.hpp"

namespace {

using namespace toolkit;

void print_usage() {
    std::printf(
        "genomic_toolkit " GTK_VERSION " -- streaming fragment QC\n"
        "\n"
        "usage:\n"
        "  genomic_toolkit sizes   <input> [options]\n"
        "  genomic_toolkit markdup <input> [options]\n"
        "  genomic_toolkit frip    <input> --peaks <peaks.bed> [options]\n"
        "\n"
        "input: .bam/.cram%s, .bedpe, or a 3+ column fragment .bed/.tsv\n"
        "       (.gz accepted where the build has htslib)\n"
        "\n"
        "options:\n"
        "  --peaks <file>     peak BED, required by `frip`\n"
        "  --region <str>     restrict to an htslib region (\"chr1\",\n"
        "                     \"chr1:1-1000000\"); repeatable, and the regions\n"
        "                     are read as one stream with summed statistics.\n"
        "                     \"chr1:1-1000000\"); needs an indexed BAM/CRAM.\n"
        "                     A fragment belongs to the region containing its\n"
        "                     START, so tiled regions partition a file exactly\n"
        "                     and counts from the parts sum to the whole.\n"
        "  --min-length <n>   drop fragments shorter than n bp (default 0;\n"
        "                     use 100 for the usual ATAC/CUT&Tag gate)\n"
        "  --max-length <n>   drop fragments longer than n bp (default 1000, 0 = off)\n"
        "  --min-mapq <n>     drop alignments below MAPQ n (default 0)\n"
        "  --single-end       include unpaired reads, using the aligned span\n"
        "  --keep-dups        do not mark duplicates (frip only)\n"
        "  --with-counts      fragment BED: read column 5 as a read-pair count and\n"
        "                     emit the fragment that many times (10x fragments.tsv).\n"
        "                     Off by default: each interval is one fragment, so a\n"
        "                     BED5/BED6 score column is never a multiplier\n"
        "  --unsorted         duplicate mode for non-coordinate-sorted input\n"
        "                     (O(distinct fragments) memory instead of O(pile-up))\n"
        "  --genome-size <n>  effective genome size, enables FRiP enrichment\n"
        "  --json             emit JSON instead of a text report\n"
        "  --histogram <file> write the size histogram as TSV\n"
        "\n"
        "exit status is 2 when the peak file and the input disagree on contig\n"
        "naming, because that failure otherwise reports a FRiP of exactly zero.\n",
        kHaveBam ? "" : " (unavailable: built without htslib)");
}

struct Options {
    std::string command;
    std::string input;
    std::string peaks;
    std::string histogram;
    std::vector<std::string> regions;
    FragmentFilter filter;
    bool keep_dups = false;
    bool unsorted = false;
    bool json = false;
    std::int64_t genome_size = 0;
};

// Returns false and prints a message on a bad argument, rather than throwing:
// argument errors are the one class of failure where a stack unwind buys
// nothing and the user wants the usage line.
bool parse_args(int argc, char** argv, Options& opt) {
    if (argc < 3) return false;
    opt.command = argv[1];
    opt.input = argv[2];

    auto need_value = [&](int& i) -> const char* {
        if (i + 1 >= argc) {
            std::fprintf(stderr, "error: %s needs a value\n", argv[i]);
            return nullptr;
        }
        return argv[++i];
    };

    // Every numeric option goes through one checked parse; toolkit/cli_args.hpp
    // lists what std::stol did instead. The range is the option's own, so a
    // value that would have been narrowed is refused by name.
    auto number = [&](int& i, std::int64_t lo, std::int64_t hi, std::int64_t& out) -> bool {
        const char* v = need_value(i);
        if (v == nullptr) return false;
        if (!parse_bounded_int(v, lo, hi, out)) {
            std::fprintf(stderr, "error: %s expects an integer in [%lld, %lld], got '%s'\n",
                         argv[i - 1], static_cast<long long>(lo), static_cast<long long>(hi),
                         v);
            return false;
        }
        return true;
    };
    std::int64_t n = 0;

    for (int i = 3; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--peaks") {
            const char* v = need_value(i);
            if (v == nullptr) return false;
            opt.peaks = v;
        } else if (a == "--region") {
            const char* v = need_value(i);
            if (v == nullptr) return false;
            opt.regions.emplace_back(v);
        } else if (a == "--histogram") {
            const char* v = need_value(i);
            if (v == nullptr) return false;
            opt.histogram = v;
        } else if (a == "--min-length") {
            if (!number(i, 0, kMaxCoord, n)) return false;
            opt.filter.min_length = static_cast<Coord>(n);
        } else if (a == "--max-length") {
            if (!number(i, 0, kMaxCoord, n)) return false;
            opt.filter.max_length = static_cast<Coord>(n);
        } else if (a == "--min-mapq") {
            if (!number(i, 0, 255, n)) return false;
            opt.filter.min_mapq = static_cast<int>(n);
        } else if (a == "--genome-size") {
            if (!number(i, 0, std::numeric_limits<std::int64_t>::max(), n)) return false;
            opt.genome_size = n;
        } else if (a == "--single-end") {
            opt.filter.allow_single_end = true;
            opt.filter.require_flags = 0;
        } else if (a == "--keep-dups") {
            opt.keep_dups = true;
        } else if (a == "--with-counts") {
            opt.filter.use_count_column = true;
        } else if (a == "--unsorted") {
            opt.unsorted = true;
        } else if (a == "--json") {
            opt.json = true;
        } else {
            std::fprintf(stderr, "error: unknown option '%s'\n", a.c_str());
            return false;
        }
    }
    return true;
}

// Returns false, having said why, when the histogram was not written.
//
// An unopenable path used to be a "warning" and the run exited 0 without the
// file it was asked for; a write that failed after the open (ENOSPC, EFBIG)
// was not looked at at all, and fclose -- where buffered rows actually reach
// the file -- was not checked either. A requested output that does not exist
// is a failed run.
bool write_histogram(const std::string& path, const SizeHistogram& h) {
    std::FILE* fp = std::fopen(path.c_str(), "w");
    if (fp == nullptr) {
        std::fprintf(stderr, "error: cannot write histogram to %s: %s\n", path.c_str(),
                     std::strerror(errno));
        return false;
    }
    bool ok = std::fprintf(fp, "length\tcount\n") >= 0;
    for (Coord i = 0; i <= SizeHistogram::kMaxTrackedLength && ok; ++i) {
        const std::uint64_t c = h.count(i);
        if (c != 0) {
            ok = std::fprintf(fp, "%d\t%llu\n", i, static_cast<unsigned long long>(c)) >= 0;
        }
    }
    const bool closed = std::fclose(fp) == 0;
    if (!ok || !closed) {
        std::fprintf(stderr, "error: writing histogram to %s failed: %s\n", path.c_str(),
                     std::strerror(errno));
        return false;
    }
    return true;
}

void report_sizes(const StreamStats& s, bool json) {
    const SizeHistogram& h = s.sizes;
    if (json) {
        std::printf(
            "{\"records_read\":%llu,\"fragments\":%llu,\"mean\":%.2f,\"median\":%d,"
            "\"p10\":%d,\"p90\":%d,\"sub_nucleosomal\":%.4f,\"mononucleosome\":%.4f,"
            "\"over_%d_bp\":%llu,\"dropped\":{\"flags\":%llu,\"mapq\":%llu,"
            "\"length\":%llu,\"unpaired\":%llu,\"interchrom\":%llu,\"malformed\":%llu}}\n",
            static_cast<unsigned long long>(s.records_read),
            static_cast<unsigned long long>(s.fragments_emitted), h.mean(), h.median(),
            h.percentile(0.10), h.percentile(0.90), h.fraction_in(1, 99),
            h.fraction_in(180, 247), SizeHistogram::kMaxTrackedLength,
            static_cast<unsigned long long>(h.overflow()),
            static_cast<unsigned long long>(s.dropped_flags),
            static_cast<unsigned long long>(s.dropped_mapq),
            static_cast<unsigned long long>(s.dropped_length),
            static_cast<unsigned long long>(s.dropped_unpaired),
            static_cast<unsigned long long>(s.dropped_interchrom),
            static_cast<unsigned long long>(s.malformed_lines));
        return;
    }
    std::printf("records read        %llu\n", static_cast<unsigned long long>(s.records_read));
    std::printf("fragments kept      %llu\n",
                static_cast<unsigned long long>(s.fragments_emitted));
    std::printf("  dropped: flags %llu  mapq %llu  length %llu  unpaired %llu  "
                "interchrom %llu  malformed %llu\n",
                static_cast<unsigned long long>(s.dropped_flags),
                static_cast<unsigned long long>(s.dropped_mapq),
                static_cast<unsigned long long>(s.dropped_length),
                static_cast<unsigned long long>(s.dropped_unpaired),
                static_cast<unsigned long long>(s.dropped_interchrom),
                static_cast<unsigned long long>(s.malformed_lines));
    std::printf("mean length         %.1f bp\n", h.mean());
    std::printf("median length       %d bp\n", h.median());
    std::printf("p10 / p90           %d / %d bp\n", h.percentile(0.10), h.percentile(0.90));
    std::printf("sub-nucleosomal     %.1f%% (< 100 bp)\n", 100.0 * h.fraction_in(1, 99));
    std::printf("mononucleosome      %.1f%% (180-247 bp)\n", 100.0 * h.fraction_in(180, 247));
    if (h.overflow() != 0) {
        std::printf("over %d bp          %llu (not binned; percentiles saturate)\n",
                    SizeHistogram::kMaxTrackedLength,
                    static_cast<unsigned long long>(h.overflow()));
    }
}

// A fragment file tolerates a handful of malformed lines -- read_peaks_bed's
// comment gives the asymmetry with peak files. A file in which EVERY record is
// malformed is not that: it is the wrong file, or the wrong delimiter, and a
// report of zero fragments with exit 0 is indistinguishable from an empty
// library. The report has already been printed; this says why it is empty.
bool nothing_parsed(const StreamStats& s, const std::string& input) {
    if (s.records_read == 0 || s.malformed_lines != s.records_read) return false;
    std::fprintf(stderr,
                 "error: none of the %llu records in '%s' parsed as a fragment -- every one "
                 "is malformed. Check the format and the delimiter (tabs).\n",
                 static_cast<unsigned long long>(s.records_read), input.c_str());
    return true;
}

int run(const Options& opt) {
    ContigDict dict;
    StreamStats stats;

    if (opt.command == "sizes") {
        stats = stream_fragments(opt.input, opt.filter, dict, [](std::span<Fragment>) {},
                                 opt.regions);
        report_sizes(stats, opt.json);
        if (nothing_parsed(stats, opt.input)) return 2;
        if (!opt.histogram.empty() && !write_histogram(opt.histogram, stats.sizes)) return 1;
        return 0;
    }

    if (opt.command == "markdup") {
        DuplicateMarker marker(opt.unsorted ? DuplicateMarker::Mode::kHashAll
                                            : DuplicateMarker::Mode::kSortedStream);
        stats = stream_fragments(opt.input, opt.filter, dict,
                                 [&](std::span<Fragment> b) { marker.observe_batch(b); },
                                 opt.regions);
        const DuplicateStats& d = marker.stats();
        if (opt.json) {
            std::printf("{\"examined\":%llu,\"duplicates\":%llu,\"unique\":%llu,"
                        "\"duplicate_rate\":%.6f,\"estimated_library_size\":%.0f,"
                        "\"peak_live_keys\":%llu,\"out_of_order\":%llu}\n",
                        static_cast<unsigned long long>(d.examined),
                        static_cast<unsigned long long>(d.duplicates),
                        static_cast<unsigned long long>(d.unique), d.duplicate_rate(),
                        marker.estimate_library_size(),
                        static_cast<unsigned long long>(d.peak_live_keys),
                        static_cast<unsigned long long>(d.out_of_order));
        } else {
            std::printf("fragments examined  %llu\n",
                        static_cast<unsigned long long>(d.examined));
            std::printf("duplicates          %llu (%.2f%%)\n",
                        static_cast<unsigned long long>(d.duplicates),
                        100.0 * d.duplicate_rate());
            std::printf("distinct            %llu\n",
                        static_cast<unsigned long long>(d.unique));
            std::printf("est. library size   %.0f molecules\n",
                        marker.estimate_library_size());
            std::printf("peak live keys      %llu (resident duplicate state)\n",
                        static_cast<unsigned long long>(d.peak_live_keys));
        }
        if (nothing_parsed(stats, opt.input)) return 2;
        if (d.out_of_order != 0) {
            std::fprintf(stderr,
                         "error: %llu records arrived out of coordinate order, so the "
                         "duplicate rate above is an undercount. Coordinate-sort the "
                         "input, or re-run with --unsorted.\n",
                         static_cast<unsigned long long>(d.out_of_order));
            return 2;
        }
        return 0;
    }

    if (opt.command == "frip") {
        if (opt.peaks.empty()) {
            std::fprintf(stderr, "error: frip needs --peaks\n");
            return 1;
        }
        PeakSet peaks = PeakSet::build(read_peaks_bed(opt.peaks));
        // PeakSet drops zero-length intervals, which cover no bases. A file of
        // nothing else used to reach the naming diagnostic below -- "chr1 is
        // absent from the peak file" about a file that names chr1 -- because
        // the dropped intervals never interned their contig.
        if (peaks.size() == 0) {
            std::fprintf(stderr,
                         "error: every interval in '%s' is zero-length (start == end), so no "
                         "base is inside a peak. BED is 0-based, half-open: a 1 bp peak at "
                         "position p is 'p-1  p'.\n",
                         opt.peaks.c_str());
            return 1;
        }

        // The peak set is bound to the stream dictionary *after* the stream
        // has filled it, which works for BAM (header interned up front) but
        // not for text input, where contigs appear as they are read. Both
        // cases are handled by pre-seeding the stream dictionary with the peak
        // contigs: peak tids then exist before the first fragment, and any
        // extra contig the stream introduces resolves to "no peaks", which is
        // the correct answer for it.
        for (const auto& name : peaks.dict().names()) dict.intern(name);
        peaks.bind(dict);

        DuplicateMarker marker(opt.unsorted ? DuplicateMarker::Mode::kHashAll
                                            : DuplicateMarker::Mode::kSortedStream);
        FripCounter counter(peaks, !opt.keep_dups);

        stats = stream_fragments(opt.input, opt.filter, dict, [&](std::span<Fragment> b) {
            if (!opt.keep_dups) marker.observe_batch(b);
            counter.add_batch(b);
        }, opt.regions);
        const FripResult r = counter.finish(dict);

        if (opt.json) {
            // `bases_in_peaks` and `bases_total` are the two terms of
            // base_frip. Emitting the ratio alone made base FRiP impossible to
            // merge across a sharded run -- the numerator and denominator
            // cannot be recovered from it, and weighting shard ratios by
            // fragment count is wrong because fragments differ in length.
            // They cost two integers and make the statistic decomposable.
            std::printf("{\"fragments\":%llu,\"in_peaks\":%llu,\"frip\":%.6f,"
                        "\"base_frip\":%.6f,\"bases_in_peaks\":%lld,"
                        "\"bases_total\":%lld,\"duplicates_excluded\":%llu,"
                        "\"peaks_input\":%zu,\"peaks_merged\":%zu,\"peak_bases\":%lld,"
                        "\"enrichment\":%.3f,\"unmatched_contigs\":%zu}\n",
                        static_cast<unsigned long long>(r.fragments_total),
                        static_cast<unsigned long long>(r.fragments_in_peaks),
                        r.fragment_frip(), r.base_frip(),
                        static_cast<long long>(r.bases_in_peaks),
                        static_cast<long long>(r.bases_total),
                        static_cast<unsigned long long>(r.duplicates_excluded),
                        r.peaks_input, r.peaks_merged,
                        static_cast<long long>(r.peak_bases),
                        r.enrichment(opt.genome_size), r.unmatched_contigs.size());
        } else {
            std::printf("peaks               %zu input -> %zu merged, %lld bp covered\n",
                        r.peaks_input, r.peaks_merged,
                        static_cast<long long>(r.peak_bases));
            std::printf("fragments counted   %llu\n",
                        static_cast<unsigned long long>(r.fragments_total));
            if (!opt.keep_dups) {
                std::printf("duplicates excluded %llu (%.2f%%)\n",
                            static_cast<unsigned long long>(r.duplicates_excluded),
                            100.0 * marker.stats().duplicate_rate());
            }
            std::printf("in peaks            %llu\n",
                        static_cast<unsigned long long>(r.fragments_in_peaks));
            std::printf("FRiP (fragments)    %.4f\n", r.fragment_frip());
            std::printf("FRiP (bases)        %.4f\n", r.base_frip());
            if (opt.genome_size > 0) {
                std::printf("enrichment          %.2fx over uniform\n",
                            r.enrichment(opt.genome_size));
            }
        }

        if (nothing_parsed(stats, opt.input)) return 2;

        // `markdup` has refused this since it was written; `frip` runs the same
        // marker in the same sorted-stream mode and did not look. An unsorted
        // input then left duplicates unmarked and counted, exit 0 -- the FRiP
        // above includes fragments it claims to have excluded.
        if (!opt.keep_dups && marker.stats().out_of_order != 0) {
            std::fprintf(stderr,
                         "error: %llu records arrived out of coordinate order, so duplicates "
                         "were undercounted and the FRiP above counts some of them. "
                         "Coordinate-sort the input, or re-run with --unsorted (or "
                         "--keep-dups).\n",
                         static_cast<unsigned long long>(marker.stats().out_of_order));
            return 2;
        }

        // A handful of unmatched contigs is normal (decoys, scaffolds, chrM),
        // so the diagnostic fires on the combination that is never innocent:
        // nothing overlapped anything AND the naming does not line up. A zero
        // FRiP is otherwise indistinguishable from a genuinely dead library.
        if (r.fragments_in_peaks == 0 && !r.unmatched_contigs.empty()) {
            std::fprintf(stderr,
                         "error: no fragment overlapped a peak, and %zu stream contigs "
                         "(e.g. '%s') are absent from the peak file. This is almost "
                         "certainly a chromosome-naming mismatch ('chr1' vs '1'), not a "
                         "FRiP of zero.\n",
                         r.unmatched_contigs.size(), r.unmatched_contigs.front().c_str());
            return 2;
        }
        return 0;
    }

    std::fprintf(stderr, "error: unknown command '%s'\n", opt.command.c_str());
    return 1;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc >= 2 && (std::strcmp(argv[1], "-h") == 0 ||
                      std::strcmp(argv[1], "--help") == 0)) {
        print_usage();
        return 0;
    }
    if (argc >= 2 && std::strcmp(argv[1], "--version") == 0) {
        std::printf("genomic_toolkit " GTK_VERSION "\n");
        return 0;
    }

    Options opt;
    if (!parse_args(argc, argv, opt)) {
        print_usage();
        return 1;
    }

    // --histogram is written AFTER the input is streamed, so naming an input
    // replaced it with a two-column table and exited 0 -- the BAM, the peak
    // BED, or the BAM's index, which the next --region query needs.
    // toolkit/file_identity.hpp holds the semantics (links, "-").
    if (!opt.histogram.empty()) {
        const std::vector<std::string> inputs = {opt.input, opt.input + ".bai",
                                                 opt.input + ".csi", opt.input + ".crai",
                                                 opt.peaks};
        if (const auto clash = output_overwrites_input(opt.histogram, inputs)) {
            std::fprintf(stderr,
                         "error: --histogram '%s' is the same file as input '%s'; writing "
                         "it would destroy the input\n",
                         opt.histogram.c_str(), clash->c_str());
            return 1;
        }
    }

    try {
        return run(opt);
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
}
