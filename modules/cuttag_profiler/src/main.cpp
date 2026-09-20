// SPDX-License-Identifier: MIT
//
// Entry point. Two modes share one core:
//   `profile` - headless batch run, for pipelines and for benchmarking.
//   `serve`   - embedded HTTP server + local browser dashboard.
//
// `serve` exists only when PTO_CLOUD_BUILD is OFF. Under the cloud-worker
// build the whole subcommand is compiled out -- not gated at runtime, not
// hidden behind a flag -- along with the server itself and the vendored
// httplib it needs. See SECURITY_HTTP_2026-08-15 H1 and the option's comment
// in CMakeLists.txt: the dashboard's session token makes it safe for one
// operator on one workstation, and a listening port inside a multi-tenant
// worker is a way around the per-job task isolation that IS the boundary
// there. pto-cloud refuses to dispatch `serve` in two places already; this is
// the third and the only one that survives a bad deploy.
#include <algorithm>
#include <charconv>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <thread>
#include <vector>

#include "profiler/bam_reader.hpp"
#include "profiler/bed_reader.hpp"
#include "profiler/file_identity.hpp"
#include "profiler/signal_calc.hpp"
#include "profiler/types.hpp"

#ifndef PTO_CLOUD_BUILD
#include "profiler/http_server.hpp"
#endif

#ifndef PROFILER_VERSION
#define PROFILER_VERSION "0.0.0"
#endif

namespace {

using namespace profiler;

void print_usage() {
    std::puts(
        "cuttag_profiler " PROFILER_VERSION
        " - fast, air-gapped CUT&Tag / CUT&RUN signal profiling\n"
        "\n"
        "USAGE\n"
        "  cuttag_profiler profile --bam FILE --regions FILE [options]\n"
#ifndef PTO_CLOUD_BUILD
        "  cuttag_profiler serve [options]\n"
#endif
        "\n"
        "PROFILE OPTIONS\n"
        "  --bam FILE               indexed BAM/CRAM (a .bai/.csi must exist)\n"
        "  --regions FILE           BED or GTF of loci (.gz accepted)\n"
        "  --out-matrix FILE        write the regions x bins matrix as TSV\n"
        "  --out-profile FILE       write the meta-profile curve as TSV\n"
        "  --reference-point WHAT   TSS | center | TES        [TSS]\n"
        "  --upstream BP            window before the anchor  [2000]\n"
        "  --downstream BP          window after the anchor   [2000]\n"
        "  --bin-size BP            bin width                 [50]\n"
        "  --normalization WHAT     CPM | RPKM | BPM | raw    [CPM]\n"
        "  --count-mode WHAT        reads | depth             [reads]\n"
        "  --min-mapq N             discard alignments below  [0]\n"
        "  --max-fragment BP        drop longer fragments     [1000]\n"
        "  --extend-reads BP        single-end extension      [0 = off]\n"
        "  --no-extend-fragment     score aligned spans, not fragments\n"
        "  --proper-pairs-only      require the SAM proper-pair flag\n"
        "  --ignore-strand          do not orient by the BED strand column\n"
        "  -p, --threads N          worker threads            [all cores]\n"
        "  -q, --quiet              suppress the progress meter\n"
        "\n"
#ifndef PTO_CLOUD_BUILD
        "SERVE OPTIONS\n"
        "  --port N                 listen port               [8080]\n"
        "  --host ADDR              bind address              [127.0.0.1]\n"
        "  --data-root DIR          confine file access to DIR (REQUIRED)\n"
        "  --no-confinement         run with no data root; single-user hosts only\n"
        "  --no-browser             do not launch a browser\n"
        "  --no-auth                do not require the session token; local scripts only\n"
        "  --allow-remote           permit a non-loopback bind (off by default)\n"
        "\n"
#else
        "This binary was built with PTO_CLOUD_BUILD=ON: it contains no HTTP\n"
        "server, and `serve` is not a subcommand.\n"
        "\n"
#endif
        "  -h, --help               this message\n"
        "  -V, --version            print the version and exit");
}

// Reads the value that follows `--flag`, erroring out when it is missing.
std::string take_value(int argc, char** argv, int& i, const char* flag) {
    if (i + 1 >= argc) {
        throw std::runtime_error(std::string("missing value for ") + flag);
    }
    return argv[++i];
}

std::int64_t to_int(const std::string& text, const char* flag) {
    std::int64_t value = 0;
    const auto* last = text.data() + text.size();
    const auto result = std::from_chars(text.data(), last, value);
    if (result.ec != std::errc{} || result.ptr != last) {
        throw std::runtime_error(std::string("expected an integer for ") + flag +
                                 ", got '" + text + "'");
    }
    return value;
}

int run_profile(int argc, char** argv) {
    std::string bam_path;
    std::string region_path;
    std::string out_matrix;
    std::string out_profile;
    ProfileOptions opts;
    FilterOptions filters;
    bool quiet = false;

    for (int i = 2; i < argc; ++i) {
        const std::string_view arg = argv[i];
        if (arg == "--bam") {
            bam_path = take_value(argc, argv, i, "--bam");
        } else if (arg == "--regions") {
            region_path = take_value(argc, argv, i, "--regions");
        } else if (arg == "--out-matrix") {
            out_matrix = take_value(argc, argv, i, "--out-matrix");
        } else if (arg == "--out-profile") {
            out_profile = take_value(argc, argv, i, "--out-profile");
        } else if (arg == "--reference-point") {
            const std::string value = take_value(argc, argv, i, "--reference-point");
            if (value == "center") {
                opts.reference_point = ReferencePoint::kCenter;
            } else if (value == "TES" || value == "tes") {
                opts.reference_point = ReferencePoint::kTES;
            } else if (value == "TSS" || value == "tss") {
                opts.reference_point = ReferencePoint::kTSS;
            } else {
                throw std::runtime_error("unknown reference point: " + value);
            }
        } else if (arg == "--upstream") {
            opts.upstream = to_int(take_value(argc, argv, i, "--upstream"), "--upstream");
        } else if (arg == "--downstream") {
            opts.downstream = to_int(take_value(argc, argv, i, "--downstream"), "--downstream");
        } else if (arg == "--bin-size") {
            opts.bin_size = to_int(take_value(argc, argv, i, "--bin-size"), "--bin-size");
        } else if (arg == "--normalization") {
            const std::string value = take_value(argc, argv, i, "--normalization");
            if (value == "CPM") opts.normalization = Normalization::kCPM;
            else if (value == "RPKM") opts.normalization = Normalization::kRPKM;
            else if (value == "BPM") opts.normalization = Normalization::kBPM;
            else if (value == "raw" || value == "none") opts.normalization = Normalization::kRaw;
            else throw std::runtime_error("unknown normalization: " + value);
        } else if (arg == "--count-mode") {
            const std::string value = take_value(argc, argv, i, "--count-mode");
            if (value == "depth") opts.count_mode = CountMode::kDepth;
            else if (value == "reads") opts.count_mode = CountMode::kReads;
            else throw std::runtime_error("unknown count mode: " + value);
        } else if (arg == "--min-mapq") {
            filters.min_mapq = narrow_option(
                to_int(take_value(argc, argv, i, "--min-mapq"), "--min-mapq"), 0, 255,
                "--min-mapq");
        } else if (arg == "--max-fragment") {
            filters.max_fragment_length = to_int(take_value(argc, argv, i, "--max-fragment"), "--max-fragment");
        } else if (arg == "--extend-reads") {
            filters.extend_reads_to = to_int(take_value(argc, argv, i, "--extend-reads"), "--extend-reads");
        } else if (arg == "--no-extend-fragment") {
            filters.extend_to_fragment = false;
        } else if (arg == "--proper-pairs-only") {
            filters.require_flags |= flags::kProperPair;
        } else if (arg == "--ignore-strand") {
            opts.respect_strand = false;
        } else if (arg == "-p" || arg == "--threads") {
            opts.threads = narrow_option(
                to_int(take_value(argc, argv, i, "--threads"), "--threads"), 0,
                ProfileOptions::kMaxThreads, "--threads");
        } else if (arg == "-q" || arg == "--quiet") {
            quiet = true;
        } else {
            throw std::runtime_error("unknown option: " + std::string(arg));
        }
    }

    if (bam_path.empty() || region_path.empty()) {
        throw std::runtime_error("both --bam and --regions are required");
    }
    // Outputs are written after the whole matrix is computed, so an output
    // that names an input replaced it and exited 0: the region file, the BAM,
    // or its index, which every later run needs. profiler/file_identity.hpp
    // holds the semantics (links, "-"). htslib looks for the index as x.bam.bai
    // AND as x.bai, so both are inputs here.
    {
        std::vector<std::string> inputs = {bam_path, bam_path + ".bai", bam_path + ".csi",
                                           bam_path + ".crai", region_path};
        if (bam_path.size() > 4 && bam_path.compare(bam_path.size() - 4, 4, ".bam") == 0) {
            inputs.push_back(bam_path.substr(0, bam_path.size() - 4) + ".bai");
        }
        for (const std::string* out : {&out_matrix, &out_profile}) {
            if (const auto clash = output_overwrites_input(*out, inputs)) {
                throw std::runtime_error("output '" + *out + "' is the same file as input '" +
                                         *clash + "'; writing it would destroy the input");
            }
        }
        // The two outputs over each other. Neither need exist yet, so identity
        // cannot be asked of the filesystem; compare the normalised absolute
        // spellings, and inodes when they do exist.
        if (!out_matrix.empty() && !out_profile.empty()) {
            std::error_code ec_m, ec_p;
            const auto norm_m = std::filesystem::absolute(out_matrix, ec_m).lexically_normal();
            const auto norm_p = std::filesystem::absolute(out_profile, ec_p).lexically_normal();
            if ((!ec_m && !ec_p && norm_m == norm_p) || same_regular_file(out_matrix, out_profile)) {
                throw std::runtime_error("--out-matrix and --out-profile name the same file '" +
                                         out_matrix + "'; one would overwrite the other");
            }
        }
    }
    // Every window/bin guard now lives on ProfileOptions itself, so the CLI and
    // the POST /api/profile handler cannot drift apart again -- they did, and
    // the network path was the one missing the checks (REVIEW_2026-08-15
    // finding 2). Adding a third entry point picks the guards up for free.
    opts.validate();
    filters.validate();

    const auto regions = read_regions(region_path);
    if (regions.empty()) {
        throw std::runtime_error("no regions parsed from " + region_path);
    }
    if (!quiet) {
        std::fprintf(stderr, "loaded %zu regions from %s\n", regions.size(),
                     region_path.c_str());
    }

    const int worker_threads =
        opts.threads > 0 ? opts.threads
                         : static_cast<int>(std::thread::hardware_concurrency());
    // See choose_bgzf_thread_count: a shared BGZF pool does not help this
    // workload at any worker count, single-threaded included, because every
    // region query calls bgzf_seek and invalidates whatever block a pool
    // would otherwise have kept inflating ahead of.
    const int bgzf_threads = choose_bgzf_thread_count(worker_threads);

    BamReader reader(bam_path, bgzf_threads);
    if (!quiet) {
        std::fprintf(stderr, "%s: %zu contigs, %llu mapped reads\n", bam_path.c_str(),
                     reader.chromosomes().size(),
                     static_cast<unsigned long long>(reader.total_mapped_reads()));
    }

    std::function<void(double)> progress;
    if (!quiet) {
        progress = [](double fraction) {
            const int width = 34;
            const int filled = static_cast<int>(fraction * width);
            std::fprintf(stderr, "\r  [%.*s%*s] %5.1f%%", filled,
                         "##################################", width - filled, "",
                         fraction * 100.0);
            std::fflush(stderr);
        };
    }

    const SignalMatrix matrix =
        compute_matrix(reader, regions, opts, filters, progress);
    if (!quiet) std::fputc('\n', stderr);

    if (!out_matrix.empty()) write_matrix_tsv(matrix, out_matrix);
    if (!out_profile.empty()) write_profile_tsv(matrix, out_profile);

    if (out_matrix.empty() && out_profile.empty()) {
        // Nothing requested on disk: emit the curve on stdout so the tool is
        // still useful in a shell pipeline.
        std::printf("offset\tmean_signal\n");
        for (std::size_t c = 0; c < matrix.cols; ++c) {
            std::printf("%lld\t%.6f\n",
                        static_cast<long long>(matrix.bin_offsets[c]),
                        matrix.column_mean[c]);
        }
    }

    std::fprintf(stderr,
                 "profiled %zu regions x %zu bins in %.2f s (%.0f regions/s)"
                 "%s\n",
                 matrix.rows, matrix.cols, matrix.elapsed_seconds,
                 matrix.rows / std::max(matrix.elapsed_seconds, 1e-9),
                 matrix.regions_skipped > 0
                     ? ("; " + std::to_string(matrix.regions_skipped) +
                        " skipped (contig missing or window off-chromosome)")
                           .c_str()
                     : "");

    // Every region skipped is not a flat profile; it is no profile. It used to
    // exit 0 over an all-zero curve, reported only as a count at the end of
    // the summary line above -- the chromosome-naming trap genomic_toolkit's
    // frip already refuses. The outputs have been written, so the evidence is
    // on disk; the exit code is 2, the repository's "the input's shape made the
    // answer wrong" status, which pto-cloud and Nextflow never retry.
    if (matrix.rows > 0 && matrix.regions_skipped == matrix.rows) {
        std::string absent;
        for (const Region& region : regions) {
            if (reader.chromosome_length(region.chrom) < 0) {
                absent = region.chrom;
                break;
            }
        }
        if (!absent.empty()) {
            std::fprintf(stderr,
                         "error: all %zu regions were skipped, and contig '%s' (among "
                         "others, possibly) is not in the BAM. This is almost certainly "
                         "a chromosome-naming mismatch ('chr1' vs '1'), not an empty "
                         "profile.\n",
                         matrix.rows, absent.c_str());
        } else {
            std::fprintf(stderr,
                         "error: all %zu regions were skipped: every window lies off "
                         "its contig\n",
                         matrix.rows);
        }
        return 2;
    }
    return 0;
}

#ifndef PTO_CLOUD_BUILD

int run_serve(int argc, char** argv) {
    ServerOptions options;
    for (int i = 2; i < argc; ++i) {
        const std::string_view arg = argv[i];
        if (arg == "--port") {
            options.port = narrow_option(to_int(take_value(argc, argv, i, "--port"), "--port"),
                                         0, 65535, "--port");
        } else if (arg == "--host") {
            options.host = take_value(argc, argv, i, "--host");
        } else if (arg == "--data-root") {
            options.data_root = take_value(argc, argv, i, "--data-root");
        } else if (arg == "--no-confinement") {
            options.confinement_disabled = true;
        } else if (arg == "--no-browser") {
            options.open_browser = false;
        } else if (arg == "--no-auth") {
            // For a local script that would otherwise have to scrape the token
            // out of the startup banner. validate() refuses to pair it with
            // --allow-remote; there is no configuration in which publishing an
            // unauthenticated read to the network is the thing someone meant.
            options.require_token = false;
        } else if (arg == "--allow-remote") {
            options.loopback_only = false;
        } else {
            throw std::runtime_error("unknown option: " + std::string(arg));
        }
    }

    // Validate before constructing, so a misconfiguration is a clean CLI error
    // rather than an exception surfacing from inside the server. HttpServer's
    // constructor validates again -- this is the friendly front end, not the
    // control. (The control has to be on the type: see REVIEW finding 2, where
    // guards living only in main.cpp left the HTTP path unprotected.)
    options.validate();

    // The --no-confinement + --allow-remote and --no-auth + --allow-remote
    // refusals are on validate(), above, so they are enforced for every
    // construction path rather than only for this one. What is left here is
    // the part that is genuinely presentation: warnings.
    if (options.confinement_disabled) {
        std::fprintf(stderr,
                     "\n"
                     "  ############################################################\n"
                     "  #  WARNING: running with NO path confinement.              #\n"
                     "  #                                                          #\n"
                     "  #  Every file this process can read is reachable through   #\n"
                     "  #  the API. The session token is the only thing standing   #\n"
                     "  #  in front of it. Use --data-root unless this is a        #\n"
                     "  #  single-user machine you trust.                          #\n"
                     "  ############################################################\n"
                     "\n");
    }
    if (!options.require_token) {
        std::fprintf(stderr,
                     "WARNING: --no-auth disables the session token. Any process "
                     "on this machine can query the API and read anything inside "
                     "the data root.\n");
    }
    if (!options.loopback_only) {
        std::fprintf(stderr,
                     "WARNING: --allow-remote binds %s:%d beyond this machine. "
                     "Genomic data will be reachable over the network.\n",
                     options.host.c_str(), options.port);
    }

    HttpServer server(std::move(options));
    server.run();
    return 0;
}

#else  // PTO_CLOUD_BUILD

// The cloud build still RECOGNISES `serve`, and refuses it by name. Falling
// through to "unknown command" would be misleading -- it would read as a typo
// rather than as a deliberate property of this artefact -- and an operator who
// reaches for the dashboard on a worker image needs to be told which binary
// they have, not left guessing.
int run_serve(int, char**) {
    throw std::runtime_error(
        "this binary was built without the HTTP server (PTO_CLOUD_BUILD=ON), "
        "so `serve` does not exist here. The dashboard authenticates one "
        "operator's browser against one workstation; inside a multi-tenant "
        "worker the isolation boundary is the per-job task, and a listening "
        "port is a way around it. Use `cuttag_profiler profile "
        "--out-matrix/--out-profile` and render the TSV through the API, or "
        "rebuild with -DPTO_CLOUD_BUILD=OFF for local use.");
}

#endif  // PTO_CLOUD_BUILD

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        print_usage();
        return 1;
    }

    const std::string_view command = argv[1];
    try {
        if (command == "-h" || command == "--help" || command == "help") {
            print_usage();
            return 0;
        }
        if (command == "-V" || command == "--version" || command == "version") {
            std::puts(PROFILER_VERSION);
            return 0;
        }
        if (command == "profile") return run_profile(argc, argv);
        if (command == "serve") return run_serve(argc, argv);

        std::fprintf(stderr, "unknown command: %s\n\n", argv[1]);
        print_usage();
        return 1;
    } catch (const std::exception& error) {
        std::fprintf(stderr, "error: %s\n", error.what());
        return 1;
    }
}
