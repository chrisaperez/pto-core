// SPDX-License-Identifier: MIT
//
// Covers windowing, strand orientation, normalisation and the parallel matrix
// builder, plus the BED/GTF front ends.
#include "test_util.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iterator>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <thread>
#include <vector>

#include "profiler/bed_reader.hpp"
#include "profiler/signal_calc.hpp"
#include "profiler/tsv_format.hpp"
#include "profiler/worker_abort.hpp"

using namespace profiler;
using testing::cigar_op;
using testing::ReadSpec;
using testing::TempBam;
using testing::TempFile;

namespace {

constexpr std::int64_t kChromLength = 100000;

std::vector<std::pair<std::string, std::int64_t>> contigs() {
    return {{"chr1", kChromLength}};
}

// A sharp pile-up of `count` 100 bp reads all starting at `pos`.
void add_pileup(std::vector<ReadSpec>& reads, std::int64_t pos, int count,
                const std::string& tag) {
    for (int i = 0; i < count; ++i) {
        ReadSpec spec;
        spec.qname = tag + std::to_string(i);
        spec.pos = pos;
        spec.cigar = {cigar_op(BAM_CMATCH, 100)};
        reads.push_back(std::move(spec));
    }
}

// True when `body` throws std::invalid_argument. The refusals added on
// 2026-09-11 are all of that type, and a refusal that silently stopped
// throwing is the regression these tests exist for.
bool throws_invalid_argument(const std::function<void()>& body) {
    try {
        body();
    } catch (const std::invalid_argument&) {
        return true;
    } catch (...) {
    }
    return false;
}

}  // namespace

int main() {
    testing::Suite suite{"signal_calc", {}};

    suite.add("make_window anchors on the requested reference point", [] {
        ProfileOptions opts;
        opts.upstream = 1000;
        opts.downstream = 500;

        Region gene{"chr1", 5000, 8000, "g1", Strand::kForward};
        Region window;

        opts.reference_point = ReferencePoint::kTSS;
        CHECK(make_window(gene, opts, kChromLength, window));
        CHECK_EQ(window.start, std::int64_t{4000});
        CHECK_EQ(window.end, std::int64_t{5500});

        opts.reference_point = ReferencePoint::kTES;
        CHECK(make_window(gene, opts, kChromLength, window));
        CHECK_EQ(window.start, std::int64_t{7000});
        CHECK_EQ(window.end, std::int64_t{8500});

        opts.reference_point = ReferencePoint::kCenter;
        CHECK(make_window(gene, opts, kChromLength, window));
        CHECK_EQ(window.start, std::int64_t{5500});
        CHECK_EQ(window.end, std::int64_t{7000});
    });

    suite.add("minus-strand windows are mirrored", [] {
        ProfileOptions opts;
        opts.upstream = 1000;
        opts.downstream = 500;
        opts.reference_point = ReferencePoint::kTSS;

        Region gene{"chr1", 5000, 8000, "g1", Strand::kReverse};
        Region window;
        CHECK(make_window(gene, opts, kChromLength, window));
        // The TSS of a minus-strand gene is its right-hand end, and upstream
        // runs toward higher coordinates.
        CHECK_EQ(window.start, std::int64_t{7500});
        CHECK_EQ(window.end, std::int64_t{9000});

        opts.respect_strand = false;
        CHECK(make_window(gene, opts, kChromLength, window));
        CHECK_EQ(window.start, std::int64_t{4000});
        CHECK_EQ(window.end, std::int64_t{5500});
    });

    suite.add("windows entirely off the contig are rejected", [] {
        ProfileOptions opts;
        opts.upstream = 2000;
        opts.downstream = 2000;
        Region window;

        // Anchor past the right-hand end of a short contig: [3000, 7000) has
        // no overlap with [0, 400).
        CHECK(!make_window({"chr1", 5000, 5100, "", Strand::kUnknown}, opts, 400,
                           window));
        // Anchor far enough left that the window closes before the contig opens.
        CHECK(!make_window({"chr1", -5000, -4900, "", Strand::kUnknown}, opts,
                           kChromLength, window));

        // A window that merely *straddles* an edge is kept; compute_matrix
        // fills the out-of-bounds bins with the missing value rather than
        // discarding an otherwise usable region.
        CHECK(make_window({"chr1", 500, 600, "", Strand::kUnknown}, opts, 400,
                          window));
        CHECK_EQ(window.start, std::int64_t{-1500});
        CHECK(make_window({"chr1", 100, 200, "", Strand::kUnknown}, opts,
                          kChromLength, window));
        CHECK_EQ(window.start, std::int64_t{-1900});
    });

    suite.add("normalization factors match their definitions", [] {
        CHECK_NEAR(normalization_factor(Normalization::kRaw, 1000000, 50, 0), 1.0, 1e-12);
        CHECK_NEAR(normalization_factor(Normalization::kCPM, 2000000, 50, 0), 0.5, 1e-12);
        CHECK_NEAR(normalization_factor(Normalization::kRPKM, 1000000, 50, 0),
                   1.0e9 / (1.0e6 * 50.0), 1e-9);
        CHECK_NEAR(normalization_factor(Normalization::kBPM, 0, 50, 250.0),
                   1.0e6 / 250.0, 1e-9);
        // A BAM with no mapped reads must not produce a division by zero.
        CHECK_NEAR(normalization_factor(Normalization::kCPM, 0, 50, 0), 1.0, 1e-12);
    });

    suite.add("the meta-profile peaks at the reference point", [] {
        std::vector<ReadSpec> reads;
        // 20 reads centred on each of three TSSs at 10000, 20000, 30000.
        for (std::int64_t tss : {10000, 20000, 30000}) {
            add_pileup(reads, tss - 50, 20, "p" + std::to_string(tss));
        }
        TempBam bam(contigs(), reads);
        BamReader reader(bam.path(), 2);

        std::vector<Region> regions{
            {"chr1", 10000, 12000, "a", Strand::kForward},
            {"chr1", 20000, 22000, "b", Strand::kForward},
            {"chr1", 30000, 32000, "c", Strand::kForward}};

        ProfileOptions opts;
        opts.upstream = 1000;
        opts.downstream = 1000;
        opts.bin_size = 100;  // 20 bins
        opts.normalization = Normalization::kRaw;
        opts.threads = 2;

        const SignalMatrix matrix =
            compute_matrix(reader, regions, opts, FilterOptions{});

        CHECK_EQ(matrix.rows, std::size_t{3});
        CHECK_EQ(matrix.cols, std::size_t{20});
        CHECK_EQ(matrix.regions_skipped, std::uint64_t{0});

        // Bins 9 and 10 straddle the anchor and hold the pile-up.
        CHECK_NEAR(matrix.column_mean[9], 20.0, 1e-9);
        CHECK_NEAR(matrix.column_mean[10], 20.0, 1e-9);
        CHECK_NEAR(matrix.column_mean[0], 0.0, 1e-9);
        CHECK_NEAR(matrix.column_mean[19], 0.0, 1e-9);

        // Bin offsets are symmetric about the anchor.
        CHECK_EQ(matrix.bin_offsets.front(), std::int64_t{-950});
        CHECK_EQ(matrix.bin_offsets.back(), std::int64_t{950});
    });

    suite.add("an asymmetric signal is flipped for minus-strand regions", [] {
        std::vector<ReadSpec> reads;
        // Signal sits 500 bp to the right of position 10000 only.
        add_pileup(reads, 10500, 30, "right");
        TempBam bam(contigs(), reads);
        BamReader reader(bam.path(), 1);

        ProfileOptions opts;
        opts.upstream = 1000;
        opts.downstream = 1000;
        opts.bin_size = 100;
        opts.normalization = Normalization::kRaw;
        opts.threads = 1;

        // A plus-strand feature at 10000 sees the peak downstream (bin 15);
        // a minus-strand feature whose TSS is also 10000 must see it upstream.
        std::vector<Region> plus{{"chr1", 10000, 10200, "p", Strand::kForward}};
        std::vector<Region> minus{{"chr1", 9800, 10000, "m", Strand::kReverse}};

        const SignalMatrix forward =
            compute_matrix(reader, plus, opts, FilterOptions{});
        const SignalMatrix reverse =
            compute_matrix(reader, minus, opts, FilterOptions{});

        CHECK_NEAR(forward.column_mean[15], 30.0, 1e-9);
        CHECK_NEAR(reverse.column_mean[4], 30.0, 1e-9);
        CHECK_NEAR(reverse.column_mean[15], 0.0, 1e-9);
    });

    suite.add("CPM scaling divides by the library size", [] {
        std::vector<ReadSpec> reads;
        add_pileup(reads, 10000, 40, "x");
        TempBam bam(contigs(), reads);
        BamReader reader(bam.path(), 1);

        std::vector<Region> regions{{"chr1", 10000, 10100, "r", Strand::kForward}};
        ProfileOptions opts;
        opts.upstream = 0;
        opts.downstream = 100;
        opts.bin_size = 100;
        opts.threads = 1;

        opts.normalization = Normalization::kRaw;
        const SignalMatrix raw = compute_matrix(reader, regions, opts, FilterOptions{});
        CHECK_NEAR(raw.column_mean[0], 40.0, 1e-9);

        opts.normalization = Normalization::kCPM;
        const SignalMatrix cpm = compute_matrix(reader, regions, opts, FilterOptions{});
        // 40 reads out of a 40-read library => 1e6 CPM.
        CHECK_NEAR(cpm.column_mean[0], 1.0e6, 1e-3);
        CHECK_EQ(cpm.total_mapped_reads, std::uint64_t{40});
    });

    suite.add("results do not depend on the thread count", [] {
        std::vector<ReadSpec> reads;
        for (int i = 0; i < 60; ++i) {
            add_pileup(reads, 10000 + i * 500, 1 + (i % 7), "r" + std::to_string(i));
        }
        TempBam bam(contigs(), reads);
        BamReader reader(bam.path(), 4);

        std::vector<Region> regions;
        for (int i = 0; i < 60; ++i) {
            regions.push_back({"chr1", 10000 + i * 500, 10000 + i * 500 + 200,
                               "r" + std::to_string(i),
                               i % 2 == 0 ? Strand::kForward : Strand::kReverse});
        }

        ProfileOptions opts;
        opts.upstream = 500;
        opts.downstream = 500;
        opts.bin_size = 50;
        opts.normalization = Normalization::kCPM;

        opts.threads = 1;
        const SignalMatrix serial = compute_matrix(reader, regions, opts, FilterOptions{});
        opts.threads = 8;
        const SignalMatrix parallel = compute_matrix(reader, regions, opts, FilterOptions{});

        CHECK_EQ(serial.rows, parallel.rows);
        CHECK_EQ(serial.cols, parallel.cols);
        int differences = 0;
        for (std::size_t i = 0; i < serial.values.size(); ++i) {
            if (std::fabs(serial.values[i] - parallel.values[i]) > 1e-9) ++differences;
        }
        CHECK_EQ(differences, 0);
        for (std::size_t c = 0; c < serial.cols; ++c) {
            CHECK_NEAR(serial.column_mean[c], parallel.column_mean[c], 1e-9);
        }
    });

    suite.add("unknown contigs are skipped, not fatal", [] {
        std::vector<ReadSpec> reads;
        add_pileup(reads, 10000, 5, "x");
        TempBam bam(contigs(), reads);
        BamReader reader(bam.path(), 1);

        std::vector<Region> regions{
            {"chr1", 10000, 10100, "good", Strand::kForward},
            {"chrUnplaced", 10000, 10100, "bad", Strand::kForward}};

        ProfileOptions opts;
        opts.upstream = 100;
        opts.downstream = 100;
        opts.bin_size = 100;
        opts.normalization = Normalization::kRaw;

        const SignalMatrix matrix = compute_matrix(reader, regions, opts, FilterOptions{});
        CHECK_EQ(matrix.rows, std::size_t{2});
        CHECK_EQ(matrix.regions_skipped, std::uint64_t{1});
        CHECK_NEAR(matrix.at(1, 0), 0.0, 1e-9);
    });

    suite.add("BED parsing handles 3, 4 and 6 column files", [] {
        TempFile bed(".bed",
                     "# a comment\n"
                     "track name=\"demo\"\n"
                     "chr1\t100\t200\n"
                     "chr1\t300\t400\tpeak2\n"
                     "chr2\t500\t600\tpeak3\t0\t-\n"
                     "chr2\t700\t800\t.\t0\t+\n");
        const auto regions = read_bed(bed.path());

        CHECK_EQ(regions.size(), std::size_t{4});
        CHECK_EQ(regions[0].chrom, std::string("chr1"));
        CHECK_EQ(regions[0].start, std::int64_t{100});
        CHECK_EQ(regions[0].end, std::int64_t{200});
        CHECK(regions[0].name.empty());
        CHECK(regions[0].strand == Strand::kUnknown);
        CHECK_EQ(regions[1].name, std::string("peak2"));
        CHECK(regions[2].strand == Strand::kReverse);
        CHECK(regions[3].name.empty());  // "." means unnamed
        CHECK(regions[3].strand == Strand::kForward);
    });

    suite.add("malformed BED coordinates raise an error", [] {
        TempFile bed(".bed", "chr1\tnotanumber\t200\n");
        bool threw = false;
        try {
            read_bed(bed.path());
        } catch (const std::exception&) {
            threw = true;
        }
        CHECK(threw);
    });

    suite.add("GTF parsing converts to 0-based half-open coordinates", [] {
        TempFile gtf(".gtf",
                     "#!genome-build test\n"
                     "chr1\tsrc\tgene\t101\t200\t.\t+\t.\t"
                     "gene_id \"ENSG1\"; gene_name \"ACTB\";\n"
                     "chr1\tsrc\texon\t101\t150\t.\t+\t.\t"
                     "gene_id \"ENSG1\"; gene_name \"ACTB\";\n"
                     "chr2\tsrc\tgene\t501\t600\t.\t-\t.\t"
                     "gene_id \"ENSG2\"; gene_name \"GAPDH\";\n");
        const auto genes = read_gtf(gtf.path());

        CHECK_EQ(genes.size(), std::size_t{2});  // the exon row is filtered out
        CHECK_EQ(genes[0].start, std::int64_t{100});
        CHECK_EQ(genes[0].end, std::int64_t{200});
        CHECK_EQ(genes[0].name, std::string("ACTB"));
        CHECK(genes[1].strand == Strand::kReverse);
        CHECK_EQ(genes[1].name, std::string("GAPDH"));

        const auto exons = read_gtf(gtf.path(), "exon");
        CHECK_EQ(exons.size(), std::size_t{1});
    });

    suite.add("read_regions dispatches on the file extension", [] {
        TempFile bed(".bed", "chr1\t100\t200\tpeak\n");
        CHECK_EQ(read_regions(bed.path()).size(), std::size_t{1});

        TempFile gtf(".gtf",
                     "chr1\tsrc\tgene\t101\t200\t.\t+\t.\tgene_name \"X\";\n");
        const auto regions = read_regions(gtf.path());
        CHECK_EQ(regions.size(), std::size_t{1});
        CHECK_EQ(regions[0].start, std::int64_t{100});
    });

    // --------------------------------------------------- security regressions
    //
    // REVIEW_2026-08-15 finding 2. AUDIT.md S4/S5/S6 were fixed in main.cpp
    // only; POST /api/profile builds the same ProfileOptions from untrusted
    // JSON and applied just a subset of the checks. The fix moved every guard
    // onto the type, so these tests exercise the *one* thing both entry points
    // now call -- which is the point of the fix: an HTTP-level test would only
    // prove the handler still forwards to validate(), while these prove the
    // guards themselves, for every present and future caller.
    suite.add("ProfileOptions::validate rejects an unbounded bin count", [] {
        // The S6 scenario, reachable unauthenticated over HTTP before the fix:
        // bin_count() = 80,000,000 and compute_matrix then allocates 640 MB of
        // double per region, so the kernel SIGKILLs the server with no
        // diagnostic and every other in-flight request dies with it.
        ProfileOptions opts;
        opts.upstream = 2000000000;
        opts.downstream = 2000000000;
        opts.bin_size = 50;
        CHECK_EQ(opts.bin_count(), std::int64_t{80000000});
        CHECK(opts.bin_count() > ProfileOptions::kMaxBins);

        bool threw = false;
        try {
            opts.validate();
        } catch (const std::invalid_argument&) {
            threw = true;
        }
        CHECK(threw);

        // Exactly at the ceiling is allowed; one bin past it is not.
        ProfileOptions at_limit;
        at_limit.bin_size = 1;
        at_limit.upstream = ProfileOptions::kMaxBins;
        at_limit.downstream = 0;
        CHECK_EQ(at_limit.bin_count(), ProfileOptions::kMaxBins);
        bool at_limit_threw = false;
        try {
            at_limit.validate();
        } catch (const std::invalid_argument&) {
            at_limit_threw = true;
        }
        CHECK(!at_limit_threw);

        at_limit.downstream = 1;
        bool over_threw = false;
        try {
            at_limit.validate();
        } catch (const std::invalid_argument&) {
            over_threw = true;
        }
        CHECK(over_threw);
    });

    suite.add("window_length does not overflow on extreme coordinates", [] {
        // {"upstream": 9223372036854775807, "downstream": 1} overflowed *inside*
        // window_length(), so the `<= 0` guard was evaluating an already-UB
        // value. UBSan flagged it as signed integer overflow in types.hpp.
        ProfileOptions opts;
        opts.upstream = std::numeric_limits<std::int64_t>::max();
        opts.downstream = 1;
        // Saturates instead of wrapping, and validate() then rejects it.
        CHECK_EQ(opts.window_length(), std::numeric_limits<std::int64_t>::max());
        bool threw = false;
        try {
            opts.validate();
        } catch (const std::invalid_argument&) {
            threw = true;
        }
        CHECK(threw);

        // Ordinary sums are untouched by the checked add.
        ProfileOptions ok;
        ok.upstream = 2000;
        ok.downstream = 3000;
        CHECK_EQ(ok.window_length(), std::int64_t{5000});
    });

    suite.add("validate rejects negative half-windows and bad bin sizes", [] {
        // A negative half-window sums positive, so it passed the old
        // window_length() > 0 check and reached make_window with a negative
        // left edge -- producing a profile whose x-axis is mislabelled, i.e. a
        // wrong scientific result reported as success.
        ProfileOptions negative;
        negative.upstream = -100000;
        negative.downstream = 200000;
        CHECK(negative.window_length() > 0);  // why the old guard missed it
        bool neg_threw = false;
        try {
            negative.validate();
        } catch (const std::invalid_argument&) {
            neg_threw = true;
        }
        CHECK(neg_threw);

        // binSize: 0 must be refused, not silently rewritten to 1. The HTTP
        // handler used to clamp with std::max(1, ...), returning a
        // plausible-looking matrix at a resolution the client never asked for.
        ProfileOptions zero_bin;
        zero_bin.bin_size = 0;
        CHECK_EQ(zero_bin.bin_count(), std::int64_t{0});  // no division by zero
        bool zero_threw = false;
        try {
            zero_bin.validate();
        } catch (const std::invalid_argument&) {
            zero_threw = true;
        }
        CHECK(zero_threw);

        ProfileOptions neg_bin;
        neg_bin.bin_size = -50;
        bool neg_bin_threw = false;
        try {
            neg_bin.validate();
        } catch (const std::invalid_argument&) {
            neg_bin_threw = true;
        }
        CHECK(neg_bin_threw);
    });

    suite.add("bin_count does not narrow a wide window to int", [] {
        // static_cast<int> on a 4e9 bin count is implementation-defined and on
        // x86-64 clang yielded a negative value, which then failed the
        // `nbins <= 0` test with a misleading "profile window is empty" error
        // rather than the real cause.
        ProfileOptions opts;
        opts.upstream = 4000000000;
        opts.downstream = 0;
        opts.bin_size = 1;
        CHECK_EQ(opts.bin_count(), std::int64_t{4000000000});
        CHECK(opts.bin_count() > 0);
    });

    suite.add("compute_matrix re-validates rather than trusting its caller", [] {
        std::vector<ReadSpec> reads;
        add_pileup(reads, 10000, 5, "x");
        TempBam bam(contigs(), reads);
        BamReader reader(bam.path(), 1);
        std::vector<Region> regions{{"chr1", 10000, 10100, "r", Strand::kForward}};

        // Straight to the allocating function with guards never applied --
        // the position a third entry point would be in.
        ProfileOptions opts;
        opts.upstream = 2000000000;
        opts.downstream = 2000000000;
        opts.bin_size = 50;
        bool threw = false;
        try {
            (void)compute_matrix(reader, regions, opts, FilterOptions{});
        } catch (const std::invalid_argument&) {
            threw = true;
        }
        CHECK(threw);
    });

    // ------------------------------------------------- worker cancellation
    //
    // compute_matrix used to capture the first worker exception under a mutex
    // and then join every worker unconditionally. The failure was recorded
    // promptly and acted on only after every sibling had finished its entire
    // stripe -- on a large region set, the whole run -- producing a matrix
    // already guaranteed to be discarded.

    suite.add("WorkerAbort signals siblings and keeps the first error", [] {
        WorkerAbort abort;
        CHECK(!abort.aborted());
        CHECK(!abort.has_error());
        abort.rethrow_if_failed();  // no-op when nothing failed

        // Siblings spinning on aborted() must all observe the stop and exit.
        // This is the property the fix turns on: without it they run to
        // completion regardless.
        constexpr int kSiblings = 8;
        std::atomic<int> exited{0};
        std::vector<std::thread> spinners;
        for (int i = 0; i < kSiblings; ++i) {
            spinners.emplace_back([&] {
                while (!abort.aborted()) {
                    std::this_thread::sleep_for(std::chrono::milliseconds(1));
                }
                exited.fetch_add(1, std::memory_order_relaxed);
            });
        }

        abort.fail(std::make_exception_ptr(std::runtime_error("first")));
        for (auto& spinner : spinners) spinner.join();

        CHECK_EQ(exited.load(), kSiblings);
        CHECK(abort.aborted());
        CHECK(abort.has_error());

        // A later failure must not displace the one the user is shown: the
        // first is the cause, the rest are usually consequences of cancelling.
        abort.fail(std::make_exception_ptr(std::runtime_error("second")));
        std::string message;
        try {
            abort.rethrow_if_failed();
        } catch (const std::runtime_error& error) {
            message = error.what();
        }
        CHECK_EQ(message, std::string("first"));
    });

    suite.add("cancel() stops workers without inventing an error", [] {
        // A caller-requested stop and a failure must stay distinguishable --
        // otherwise a cancelled run reports a spurious exception, and a fix
        // that conflated the two would silently truncate good results.
        WorkerAbort abort;
        abort.cancel();
        CHECK(abort.aborted());
        CHECK(!abort.has_error());
        bool threw = false;
        try {
            abort.rethrow_if_failed();
        } catch (...) {
            threw = true;
        }
        CHECK(!threw);
    });

    suite.add("a failing worker stops the run instead of draining its stripe", [] {
        // A real exception out of a real worker, through the public API, with
        // the wasted work measured rather than asserted about.
        //
        // The trigger: BamReader pools exactly ONE handle (the bootstrap one
        // that read the header), so with several threads every worker but the
        // first calls acquire() -> hts_open(path), which throws once the file
        // is gone. The header and index are already resident, so the reader
        // stays usable for everything except opening a NEW handle. This is not
        // a contrived setup -- a scratch directory reaped mid-run, or an NFS
        // mount dropping, arrives at exactly this path.
        //
        // Note the two readers below are deliberately separate. Reusing one
        // would leave four handles in its pool from the baseline run, so no
        // worker would ever need to open the file again and the failure would
        // not trigger at all -- the run would quietly succeed against the
        // unlinked-but-still-open descriptors and this test would pass while
        // measuring nothing.
        constexpr int kRegions = 12000;
        std::vector<Region> regions;
        for (int i = 0; i < kRegions; ++i) {
            regions.push_back({"chr1", 10000 + (i % 40) * 500,
                               10000 + (i % 40) * 500 + 200,
                               "r" + std::to_string(i), Strand::kForward});
        }

        ProfileOptions opts;
        opts.upstream = 500;
        opts.downstream = 500;
        opts.bin_size = 10;
        opts.normalization = Normalization::kRaw;
        opts.threads = 4;

        auto build_reads = [] {
            std::vector<ReadSpec> reads;
            for (int i = 0; i < 40; ++i) {
                add_pileup(reads, 10000 + i * 500, 5, "r" + std::to_string(i));
            }
            return reads;
        };
        auto elapsed_ms = [](auto start) {
            return std::chrono::duration<double, std::milli>(
                       std::chrono::steady_clock::now() - start)
                .count();
        };

        // Best-of-N, matching the convention in PERFORMANCE_2026-08-15: this
        // measures a ratio between two runs on one machine, and a single
        // scheduling hiccup in either would otherwise decide the result.
        constexpr int kRepeats = 3;
        double failing = std::numeric_limits<double>::max();
        bool threw = false;
        for (int rep = 0; rep < kRepeats; ++rep) {
            auto bam = std::make_unique<TempBam>(contigs(), build_reads());
            BamReader reader(bam->path(), 1);
            std::error_code ec;
            std::filesystem::remove(bam->path(), ec);
            CHECK(!ec);

            const auto start = std::chrono::steady_clock::now();
            try {
                (void)compute_matrix(reader, regions, opts, FilterOptions{});
            } catch (const std::exception&) {
                threw = true;
            }
            failing = std::min(failing, elapsed_ms(start));
            bam.reset();  // its destructor tolerates the already-removed file
        }

        // First contract: a catchable std::exception on the CALLER's thread.
        // Before the per-worker tally allocation was moved inside the try, that
        // same allocation could throw where nothing was catching, and an
        // exception escaping a thread's entry point is std::terminate -- an
        // outcome no test can observe, only a crashed suite.
        CHECK(threw);

        double baseline = std::numeric_limits<double>::max();
        for (int rep = 0; rep < kRepeats; ++rep) {
            TempBam bam(contigs(), build_reads());
            BamReader reader(bam.path(), 1);
            const auto start = std::chrono::steady_clock::now();
            (void)compute_matrix(reader, regions, opts, FilterOptions{});
            baseline = std::min(baseline, elapsed_ms(start));
        }

        // Second contract, and the actual regression: a failed run must not
        // cost what a successful one costs.
        //
        // The ratio is self-calibrating, so this does not encode any machine's
        // speed. Measured on the reference host: before the fix, 0.47-0.51 --
        // worker 0 held the one pooled handle and drained its entire quarter of
        // the regions after its three siblings had already failed. After, 0.031
        // -0.033. The 0.25 threshold sits an order of magnitude clear of the
        // fixed behaviour and well below the broken one.
        //
        // Skipped rather than asserted when the baseline is too short to
        // measure: on a machine fast enough to finish 12k regions in under
        // 50 ms, thread start-up dominates both numbers and the ratio stops
        // meaning anything. Reported, not silently passed.
        if (baseline < 50.0) {
            std::printf(
                "    (skipped ratio check: baseline %.1f ms is too short to "
                "measure against)\n",
                baseline);
        } else {
            CHECK(failing < 0.25 * baseline);
        }
    });

    // -------------------------------------------- reduction thread bounding
    //
    // The column-mean reduction was hardcoded to `min(nthreads, 8)` with no
    // stated reason. The 8 was doing real work -- each reducer owns a private
    // cols-wide accumulator, so scratch memory scales with thread count TIMES
    // matrix width -- so it is now three named bounds instead of one magic
    // number. These pin the arithmetic the replacement must not break.

    suite.add("column means are exact regardless of the reducer count", [] {
        // A wide matrix and a high thread request together drive the memory
        // bound, which is the bound that did not exist before. The reduction
        // must still be exact, and identical to the single-threaded answer.
        std::vector<ReadSpec> reads;
        for (int i = 0; i < 50; ++i) {
            add_pileup(reads, 10000 + i * 300, 1 + (i % 5), "w" + std::to_string(i));
        }
        TempBam bam(contigs(), reads);
        BamReader reader(bam.path(), 2);

        std::vector<Region> regions;
        for (int i = 0; i < 200; ++i) {
            regions.push_back({"chr1", 10000 + (i % 50) * 300,
                               10000 + (i % 50) * 300 + 200,
                               "w" + std::to_string(i),
                               i % 2 == 0 ? Strand::kForward : Strand::kReverse});
        }

        ProfileOptions opts;
        opts.upstream = 2000;
        opts.downstream = 2000;
        opts.bin_size = 10;  // 400 bins: wide enough to bind the memory bound
        opts.normalization = Normalization::kRaw;

        opts.threads = 1;
        const SignalMatrix serial = compute_matrix(reader, regions, opts, FilterOptions{});
        opts.threads = 64;  // far above the old hardcoded 8
        const SignalMatrix wide = compute_matrix(reader, regions, opts, FilterOptions{});

        CHECK_EQ(serial.cols, wide.cols);
        for (std::size_t c = 0; c < serial.cols; ++c) {
            CHECK_NEAR(serial.column_mean[c], wide.column_mean[c], 1e-9);
        }

        // And the means really are the means -- a reduction that silently
        // dropped a stripe would still be self-consistent across thread counts.
        for (std::size_t c = 0; c < wide.cols; ++c) {
            double expected = 0.0;
            for (std::size_t r = 0; r < wide.rows; ++r) expected += wide.at(r, c);
            expected /= static_cast<double>(wide.rows);
            CHECK_NEAR(wide.column_mean[c], expected, 1e-9);
        }
    });

    suite.add("the reduction covers every row when rows are scarce", [] {
        // The rows-per-reducer bound must never round down to zero threads, and
        // the stripe arithmetic must still cover row 0 through row n-1 when
        // there are fewer rows than requested threads.
        std::vector<ReadSpec> reads;
        add_pileup(reads, 10000, 7, "solo");
        TempBam bam(contigs(), reads);
        BamReader reader(bam.path(), 1);

        std::vector<Region> regions{{"chr1", 10000, 10100, "only", Strand::kForward}};

        ProfileOptions opts;
        opts.upstream = 100;
        opts.downstream = 100;
        opts.bin_size = 50;
        opts.normalization = Normalization::kRaw;
        opts.threads = 32;  // many more threads than rows

        const SignalMatrix matrix = compute_matrix(reader, regions, opts, FilterOptions{});
        CHECK_EQ(matrix.rows, std::size_t{1});
        for (std::size_t c = 0; c < matrix.cols; ++c) {
            CHECK_NEAR(matrix.column_mean[c], matrix.at(0, c), 1e-9);
        }
    });

    suite.add("center anchor does not overflow on extreme coordinates", [] {
        // REGRESSION. `anchor_for` computed the centre as
        // `(region.start + region.end) / 2`, which is undefined behaviour once
        // the two coordinates sum past INT64_MAX. bed_reader accepts any
        // non-negative pair with end >= start, so this was reachable from a
        // plain BED line, and UBSan reported it directly:
        //
        //   signed integer overflow: 9223372036854775000 +
        //   9223372036854775007 cannot be represented in type 'int64_t'
        //
        // The wrapped sum is always negative, so for a narrow band of inputs
        // the anchor landed just below zero and make_window's off-contig test
        // PASSED -- the region was profiled at coordinate ~0 rather than
        // skipped, contributing signal from the wrong locus to the profile.
        // The same region under kTSS was correctly skipped, which is the
        // asymmetry that gives the bug away.
        //
        // Both halves are asserted: the extreme cases are rejected, and the
        // ordinary midpoint is unchanged.
        constexpr std::int64_t kMax = std::numeric_limits<std::int64_t>::max();
        ProfileOptions opts;
        opts.reference_point = ReferencePoint::kCenter;
        opts.upstream = 2000;
        opts.downstream = 2000;
        opts.bin_size = 50;

        Region window;
        // Sums past INT64_MAX. Every one of these must be refused, not
        // silently relocated to the start of the contig.
        CHECK(!make_window({"chr1", kMax, kMax, "", Strand::kUnknown}, opts,
                           kChromLength, window));
        CHECK(!make_window({"chr1", kMax - 807, kMax - 800, "", Strand::kUnknown},
                           opts, kChromLength, window));
        CHECK(!make_window({"chr1", kMax / 2 + 1, kMax, "", Strand::kUnknown},
                           opts, kChromLength, window));

        // ...and the ordinary case is bit-identical to what it always was.
        CHECK(make_window({"chr1", 5000, 6000, "", Strand::kUnknown}, opts,
                          kChromLength, window));
        CHECK_EQ(window.start, std::int64_t{3500});   // 5500 - 2000
        CHECK_EQ(window.end, std::int64_t{7500});     // 5500 + 2000

        // An odd span still truncates toward the lower coordinate, as integer
        // division always did.
        CHECK(make_window({"chr1", 5000, 5001, "", Strand::kUnknown}, opts,
                          kChromLength, window));
        CHECK_EQ(window.start, std::int64_t{3000});   // 5000 - 2000
        CHECK_EQ(window.end, std::int64_t{7000});     // 5000 + 2000

        // A span covering the whole coordinate range. `end - start` cannot
        // overflow here because start >= 0, which is what the readers
        // guarantee, so the midpoint is computed exactly -- and against a real
        // contig length it is then rejected as off-chromosome, which is the
        // right answer for a locus at 4.6e18.
        CHECK(!make_window({"chr1", 0, kMax, "", Strand::kUnknown}, opts,
                           kChromLength, window));

        // Same input against a contig big enough to hold it, so the midpoint
        // itself is asserted rather than inferred from a rejection. This is the
        // arithmetic the old form could not do at all.
        CHECK(make_window({"chr1", 0, kMax, "", Strand::kUnknown}, opts, kMax,
                          window));
        CHECK_EQ(window.start, kMax / 2 - 2000);
        CHECK_EQ(window.end, kMax / 2 + 2000);
    });

    // ---- 2026-09-11 audit (docs/AUDIT_2026-09-11_cuttag_profiler.md) ----

    suite.add("skipped regions do not dilute the meta-profile", [] {
        // Half the regions on a contig absent from the BAM used to halve the
        // profile: the skipped rows were filled with the missing value and
        // still divided into the column mean (demo data: 5.000363 -> 2.500181),
        // exit 0.
        std::vector<ReadSpec> reads;
        add_pileup(reads, 10000 - 50, 20, "p");
        TempBam bam(contigs(), reads);
        BamReader reader(bam.path(), 1);

        ProfileOptions opts;
        opts.upstream = 1000;
        opts.downstream = 1000;
        opts.bin_size = 100;
        opts.normalization = Normalization::kRaw;
        opts.threads = 1;

        const std::vector<Region> good{{"chr1", 10000, 12000, "a", Strand::kForward}};
        std::vector<Region> mixed = good;
        mixed.push_back({"chrAbsent", 10000, 12000, "b", Strand::kForward});
        mixed.push_back({"chrAbsent", 20000, 22000, "", Strand::kForward});

        const SignalMatrix alone = compute_matrix(reader, good, opts, FilterOptions{});
        const SignalMatrix with_skips = compute_matrix(reader, mixed, opts, FilterOptions{});
        CHECK_EQ(with_skips.rows, std::size_t{3});
        CHECK_EQ(with_skips.regions_skipped, std::uint64_t{2});
        CHECK_NEAR(alone.column_mean[9], 20.0, 1e-9);
        for (std::size_t c = 0; c < alone.cols; ++c) {
            CHECK_NEAR(with_skips.column_mean[c], alone.column_mean[c], 1e-9);
        }

        // The skipped rows are still in the matrix, and still identifiable:
        // an unnamed region is named by its coordinates, as a kept one is.
        CHECK_EQ(with_skips.row_names[1], std::string("b"));
        CHECK_EQ(with_skips.row_names[2], std::string("chrAbsent:20000-22000"));
        CHECK_NEAR(with_skips.at(1, 9), opts.missing_value, 1e-12);

        // BPM sums signal across regions; a skipped row must not enter it.
        opts.normalization = Normalization::kBPM;
        const SignalMatrix bpm_alone = compute_matrix(reader, good, opts, FilterOptions{});
        const SignalMatrix bpm_mixed = compute_matrix(reader, mixed, opts, FilterOptions{});
        CHECK_NEAR(bpm_mixed.scale_factor, bpm_alone.scale_factor, 1e-12);
        for (std::size_t c = 0; c < bpm_alone.cols; ++c) {
            CHECK_NEAR(bpm_mixed.column_mean[c], bpm_alone.column_mean[c], 1e-6);
        }
    });

    suite.add("with every region skipped the profile is the missing value, not zero", [] {
        std::vector<ReadSpec> reads;
        add_pileup(reads, 10000, 5, "x");
        TempBam bam(contigs(), reads);
        BamReader reader(bam.path(), 1);

        ProfileOptions opts;
        opts.upstream = 100;
        opts.downstream = 100;
        opts.bin_size = 100;
        opts.normalization = Normalization::kRaw;
        opts.missing_value = -1.0;  // distinguishable from a real zero

        const std::vector<Region> regions{
            {"chrAbsent", 10000, 10100, "a", Strand::kForward},
            {"chrAlsoAbsent", 500, 600, "b", Strand::kForward}};
        const SignalMatrix matrix = compute_matrix(reader, regions, opts, FilterOptions{});
        CHECK_EQ(matrix.regions_skipped, std::uint64_t{2});
        for (double v : matrix.column_mean) CHECK_NEAR(v, -1.0, 1e-12);
    });

    suite.add("validate caps the window so bin geometry cannot overflow", [] {
        // Passed the bin-count ceiling (4e6 bins) and overflowed `width * b`
        // in compute_matrix: UBSan "signed integer overflow:
        // 4000000000000000001 * 3" at signal_calc.cpp:196.
        ProfileOptions huge;
        huge.upstream = 4000000000000000000;
        huge.downstream = 1;
        huge.bin_size = 1000000000000;
        CHECK(huge.bin_count() <= ProfileOptions::kMaxBins);  // why the old guard passed it
        CHECK(throws_invalid_argument([&] { huge.validate(); }));

        // The function that did the overflowing arithmetic refuses it too.
        std::vector<ReadSpec> reads;
        add_pileup(reads, 10000, 1, "x");
        TempBam bam(contigs(), reads);
        BamReader reader(bam.path(), 1);
        const std::vector<Region> regions{{"chr1", 10000, 10100, "r", Strand::kForward}};
        CHECK(throws_invalid_argument(
            [&] { (void)compute_matrix(reader, regions, huge, FilterOptions{}); }));

        // At the cap, with a bin count inside kMaxBins: accepted. One base past
        // it: refused.
        ProfileOptions at_cap;
        at_cap.upstream = ProfileOptions::kMaxWindow;
        at_cap.downstream = 0;
        at_cap.bin_size = ProfileOptions::kMaxWindow / ProfileOptions::kMaxBins + 1;
        CHECK(at_cap.bin_count() <= ProfileOptions::kMaxBins);
        CHECK(!throws_invalid_argument([&] { at_cap.validate(); }));
        at_cap.downstream = 1;
        CHECK(throws_invalid_argument([&] { at_cap.validate(); }));

        // The property the cap exists for: the largest product fits.
        CHECK(ProfileOptions::kMaxWindow <=
              std::numeric_limits<std::int64_t>::max() / ProfileOptions::kMaxBins);
    });

    suite.add("FilterOptions::validate refuses filters that silently change the answer", [] {
        const auto refused = [](const std::function<void(FilterOptions&)>& set) {
            FilterOptions f;
            set(f);
            return throws_invalid_argument([&] { f.validate(); });
        };
        CHECK(!refused([](FilterOptions&) {}));
        CHECK(refused([](FilterOptions& f) { f.min_mapq = 256; }));
        CHECK(refused([](FilterOptions& f) { f.min_mapq = -1; }));
        CHECK(refused([](FilterOptions& f) { f.max_fragment_length = -5; }));
        CHECK(refused([](FilterOptions& f) { f.extend_reads_to = -100; }));
        CHECK(refused([](FilterOptions& f) {
            f.extend_reads_to = FilterOptions::kMaxExtension + 1;
        }));
        CHECK(!refused([](FilterOptions& f) {
            f.min_mapq = 255;
            f.max_fragment_length = 0;
            f.extend_reads_to = FilterOptions::kMaxExtension;
        }));

        // compute_matrix re-checks, as it does for ProfileOptions.
        std::vector<ReadSpec> reads;
        add_pileup(reads, 10000, 1, "x");
        TempBam bam(contigs(), reads);
        BamReader reader(bam.path(), 1);
        const std::vector<Region> regions{{"chr1", 10000, 10100, "r", Strand::kForward}};
        ProfileOptions opts;
        opts.upstream = 100;
        opts.downstream = 100;
        FilterOptions bad;
        bad.min_mapq = 300;
        CHECK(throws_invalid_argument(
            [&] { (void)compute_matrix(reader, regions, opts, bad); }));
    });

    suite.add("narrow_option refuses what an int cast would have wrapped", [] {
        CHECK_EQ(narrow_option(255, 0, 255, "--min-mapq"), 255);
        CHECK_EQ(narrow_option(0, 0, ProfileOptions::kMaxThreads, "--threads"), 0);
        CHECK(throws_invalid_argument(
            [] { (void)narrow_option(4294967296LL, 0, 255, "--min-mapq"); }));
        CHECK(throws_invalid_argument([] { (void)narrow_option(300, 0, 255, "--min-mapq"); }));
        CHECK(throws_invalid_argument([] { (void)narrow_option(-1, 0, 255, "--min-mapq"); }));
        CHECK(throws_invalid_argument([] {
            (void)narrow_option(4294967297LL, 0, ProfileOptions::kMaxThreads, "--threads");
        }));

        ProfileOptions negative_threads;
        negative_threads.threads = -3;
        CHECK(throws_invalid_argument([&] { negative_threads.validate(); }));
    });

    suite.add("a GTF name attribute is matched as a whole key, not a substring", [] {
        TempFile gtf(".gtf",
                     "chr1\tsrc\tgene\t101\t200\t.\t+\t.\t"
                     "gene_id \"G1\"; havana_gene_name \"WRONG\"; gene_name \"RIGHT\";\n"
                     "chr1\tsrc\tgene\t301\t400\t.\t+\t.\t"
                     "gene_name \"FIRST\"; other \"x\";\n"
                     "chr1\tsrc\tgene\t501\t600\t.\t+\t.\t"
                     "ID=g3;Name=n3;gene_name=GFF3\n"
                     "chr1\tsrc\tgene\t701\t800\t.\t+\t.\t"
                     "gene_id \"G4\"; havana_gene_name \"ONLY_HAVANA\";\n");
        const auto genes = read_gtf(gtf.path());
        CHECK_EQ(genes.size(), std::size_t{4});
        if (genes.size() == 4) {
            CHECK_EQ(genes[0].name, std::string("RIGHT"));
            CHECK_EQ(genes[1].name, std::string("FIRST"));
            CHECK_EQ(genes[2].name, std::string("GFF3"));
            CHECK(genes[3].name.empty());  // no gene_name key at all
        }
    });

    // write_matrix_tsv's per-value formatting (currently std::ostream with
    // `fixed`/`precision(6)`) sits inside the timed region of every `profile`
    // run and dominates it at real matrix sizes (docs/AUDIT... finding 4).
    // Nothing here pins its output today, which is exactly the risk in
    // replacing it: these tests exist BEFORE that replacement, against a
    // helper (append_fixed6) the replacement will introduce, so the two are
    // compared byte-for-byte rather than one being trusted to restate the
    // other.
    suite.add("append_fixed6 matches ostream fixed/precision(6), including its edge cases", [] {
        auto ostream_fmt = [](double v) {
            std::ostringstream oss;
            oss.setf(std::ios::fixed);
            oss.precision(6);
            oss << v;
            return oss.str();
        };
        auto chars_fmt = [](double v) {
            std::string out;
            append_fixed6(out, v);
            return out;
        };

        const double values[] = {
            0.0, -0.0, 1.0, -1.0, 123.456789, -123.456789,
            1e-10, -1e-10, 1e20, -1e20,
            0.0000005, 0.00000049999,  // rounds up / down at the 6th decimal
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::lowest(),
            std::numeric_limits<double>::min(),
            std::numeric_limits<double>::denorm_min(),
            std::numeric_limits<double>::infinity(),
            -std::numeric_limits<double>::infinity(),
        };
        for (double v : values) {
            CHECK_EQ(chars_fmt(v), ostream_fmt(v));
        }

        // NaN is checked separately, and NOT against this platform's ostream:
        // whether operator<< signs a negative NaN's mantissa bit is a libc++
        // vs. libstdc++ difference, not a standard-mandated one. libc++
        // (macOS, where this was written) never signs it, matching the
        // comment this replaced; libstdc++ (glibc, Linux) does, printing
        // "-nan" for -std::nan(""). append_fixed6 special-cases NaN and
        // always emits unsigned "nan" BEFORE reaching std::to_chars -- see
        // tsv_format.hpp -- specifically so the matrix TSV doesn't carry a
        // sign that means nothing downstream and varies by build platform.
        // The contract that matters is chars_fmt's own: both signs of NaN
        // produce the identical, unsigned string, everywhere.
        CHECK_EQ(chars_fmt(std::nan("")), std::string("nan"));
        CHECK_EQ(chars_fmt(-std::nan("")), std::string("nan"));
    });

    suite.add("write_matrix_tsv is byte-identical to the ostream version it replaces", [] {
        SignalMatrix m;
        m.rows = 2;
        m.cols = 3;
        m.row_names = {"geneA", "geneB"};
        m.bin_offsets = {-100, 0, 100};
        m.values = {
            0.0, 123.456789, -1e-10,
            std::numeric_limits<double>::max(), std::nan(""), -0.0,
        };

        const testing::TempFile out_new("", "");
        write_matrix_tsv(m, out_new.path());

        // The exact formatter write_matrix_tsv used before this change,
        // reproduced verbatim so the two write functions can be diffed
        // without needing the pre-change binary on hand.
        std::ostringstream expected;
        expected << "region";
        for (std::int64_t offset : m.bin_offsets) expected << '\t' << offset;
        expected << '\n';
        expected.setf(std::ios::fixed);
        expected.precision(6);
        for (std::size_t r = 0; r < m.rows; ++r) {
            expected << m.row_names[r];
            for (std::size_t c = 0; c < m.cols; ++c) expected << '\t' << m.at(r, c);
            expected << '\n';
        }

        std::ifstream got_stream(out_new.path());
        const std::string got((std::istreambuf_iterator<char>(got_stream)),
                              std::istreambuf_iterator<char>());
        CHECK_EQ(got, expected.str());
    });

    // /api/matrix (http_server.cpp) builds this same TSV shape as an
    // in-memory HTTP response body rather than a file, and nothing tests its
    // body today. format_matrix_tsv is the extracted, HTTP-infrastructure-free
    // unit that change is built on: same byte layout as write_matrix_tsv's
    // file output (verified directly against it, not just against ostream
    // again), so a body downloaded from the dashboard and a file from
    // --out-matrix for the same run are identical either way.
    suite.add("format_matrix_tsv matches write_matrix_tsv's file output exactly", [] {
        SignalMatrix m;
        m.rows = 2;
        m.cols = 3;
        m.row_names = {"geneA", "geneB"};
        m.bin_offsets = {-100, 0, 100};
        m.values = {
            0.0, 123.456789, -1e-10,
            std::numeric_limits<double>::max(), std::nan(""), -0.0,
        };

        const testing::TempFile out_file("", "");
        write_matrix_tsv(m, out_file.path());
        std::ifstream file_stream(out_file.path());
        const std::string from_file((std::istreambuf_iterator<char>(file_stream)),
                                    std::istreambuf_iterator<char>());

        CHECK_EQ(format_matrix_tsv(m), from_file);
    });

    return suite.run();
}
