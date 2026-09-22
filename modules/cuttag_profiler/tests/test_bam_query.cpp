// SPDX-License-Identifier: MIT
//
// Exercises the indexed-query layer against BAMs synthesised at run time, so
// the suite needs no checked-in fixtures and no reference genome.
#include "test_util.hpp"

#include <algorithm>
#include <atomic>
#include <string>
#include <thread>
#include <vector>

#include "profiler/bam_reader.hpp"

using namespace profiler;
using testing::cigar_op;
using testing::ReadSpec;
using testing::TempBam;

namespace {

constexpr std::int64_t kChromLength = 100000;

// Ten 100 bp reads tiling chr1:1000-2000, one per 100 bp bin.
std::vector<ReadSpec> tiling_reads() {
    std::vector<ReadSpec> reads;
    for (int i = 0; i < 10; ++i) {
        ReadSpec spec;
        spec.qname = "tile" + std::to_string(i);
        spec.pos = 1000 + i * 100;
        spec.cigar = {cigar_op(BAM_CMATCH, 100)};
        reads.push_back(std::move(spec));
    }
    return reads;
}

std::vector<std::pair<std::string, std::int64_t>> contigs() {
    return {{"chr1", kChromLength}, {"chr2", kChromLength}};
}

}  // namespace

int main() {
    testing::Suite suite{"bam_query", {}};

    // A shared BGZF pool sized to hardware_concurrency() for the single
    // region-query-worker case measured 3.4x slower wall time and ~14
    // core-seconds of system time per wall second on 100k indexed region
    // queries (docs/AUDIT_2026-09-11_cuttag_profiler.md finding 1): every
    // such query calls bgzf_seek, which invalidates the currently-decoded
    // block, so there is no run of sequential blocks under one query for a
    // pool to overlap ahead of, only a per-block queue hand-off cost paid a
    // hundred thousand times. This never scales with worker count.
    suite.add("no BGZF pool is built for the region-query workload, at any worker count", [] {
        for (int worker_threads : {1, 2, 4, 8, 12, 32}) {
            CHECK_EQ(choose_bgzf_thread_count(worker_threads), 0);
        }
    });

    suite.add("header exposes contigs and lengths", [] {
        TempBam bam(contigs(), tiling_reads());
        BamReader reader(bam.path(), 2);

        CHECK_EQ(reader.chromosomes().size(), std::size_t{2});
        CHECK_EQ(reader.chromosomes()[0], std::string("chr1"));
        CHECK_EQ(reader.chromosome_length("chr1"), kChromLength);
        CHECK_EQ(reader.chromosome_length("chr2"), kChromLength);
        CHECK_EQ(reader.chromosome_length("chrX"), std::int64_t{-1});
        CHECK_EQ(reader.tid_for("chr2"), 1);
        CHECK_EQ(reader.tid_for("nonexistent"), -1);
    });

    suite.add("index metadata yields the mapped-read total", [] {
        TempBam bam(contigs(), tiling_reads());
        BamReader reader(bam.path(), 1);
        // Read from the .bai statistics, not by streaming the records.
        CHECK_EQ(reader.total_mapped_reads(), std::uint64_t{10});
    });

    suite.add("count_reads honours interval boundaries", [] {
        TempBam bam(contigs(), tiling_reads());
        BamReader reader(bam.path(), 1);
        FilterOptions filters;

        CHECK_EQ(reader.count_reads({"chr1", 1000, 2000}, filters), std::uint64_t{10});
        CHECK_EQ(reader.count_reads({"chr1", 1000, 1100}, filters), std::uint64_t{1});
        // A read ending exactly at the interval start must not be counted.
        CHECK_EQ(reader.count_reads({"chr1", 1100, 1200}, filters), std::uint64_t{1});
        CHECK_EQ(reader.count_reads({"chr1", 5000, 6000}, filters), std::uint64_t{0});
        CHECK_EQ(reader.count_reads({"chr2", 1000, 2000}, filters), std::uint64_t{0});
        CHECK_EQ(reader.count_reads({"chrX", 1000, 2000}, filters), std::uint64_t{0});
    });

    suite.add("query_bins places each read in its own bin", [] {
        TempBam bam(contigs(), tiling_reads());
        BamReader reader(bam.path(), 1);
        FilterOptions filters;

        std::vector<BinTally> bins(10);
        reader.query_bins({"chr1", 1000, 2000}, bins, filters);

        for (std::size_t i = 0; i < bins.size(); ++i) {
            CHECK_NEAR(bins[i].read_count, 1.0, 1e-9);
            CHECK_NEAR(bins[i].base_count, 100.0, 1e-9);
        }
    });

    suite.add("a read spanning two bins is counted in both", [] {
        std::vector<ReadSpec> reads;
        ReadSpec spec;
        spec.qname = "straddle";
        spec.pos = 1050;  // spans 1050-1150, i.e. bins 0 and 1
        spec.cigar = {cigar_op(BAM_CMATCH, 100)};
        reads.push_back(spec);

        TempBam bam(contigs(), reads);
        BamReader reader(bam.path(), 1);
        FilterOptions filters;

        std::vector<BinTally> bins(10);
        reader.query_bins({"chr1", 1000, 2000}, bins, filters);

        CHECK_NEAR(bins[0].read_count, 1.0, 1e-9);
        CHECK_NEAR(bins[1].read_count, 1.0, 1e-9);
        CHECK_NEAR(bins[2].read_count, 0.0, 1e-9);
        // Bases split evenly across the boundary.
        CHECK_NEAR(bins[0].base_count, 50.0, 1e-9);
        CHECK_NEAR(bins[1].base_count, 50.0, 1e-9);
    });

    suite.add("MAPQ filter discards low-confidence alignments", [] {
        std::vector<ReadSpec> reads = tiling_reads();
        for (std::size_t i = 0; i < reads.size(); i += 2) reads[i].mapq = 5;

        TempBam bam(contigs(), reads);
        BamReader reader(bam.path(), 1);

        FilterOptions permissive;
        CHECK_EQ(reader.count_reads({"chr1", 1000, 2000}, permissive), std::uint64_t{10});

        FilterOptions strict;
        strict.min_mapq = 30;
        CHECK_EQ(reader.count_reads({"chr1", 1000, 2000}, strict), std::uint64_t{5});
    });

    suite.add("flag filters drop duplicates and secondaries", [] {
        std::vector<ReadSpec> reads = tiling_reads();
        reads[0].flag = flags::kDuplicate;
        reads[1].flag = flags::kSecondary;
        reads[2].flag = flags::kQCFail;
        reads[3].flag = flags::kSupplementary;

        TempBam bam(contigs(), reads);
        BamReader reader(bam.path(), 1);

        FilterOptions defaults;  // excludes all four by default
        CHECK_EQ(reader.count_reads({"chr1", 1000, 2000}, defaults), std::uint64_t{6});

        FilterOptions keep_everything;
        keep_everything.exclude_flags = flags::kUnmapped;
        CHECK_EQ(reader.count_reads({"chr1", 1000, 2000}, keep_everything),
                 std::uint64_t{10});
    });

    suite.add("depth accounting is CIGAR-aware", [] {
        std::vector<ReadSpec> reads;
        ReadSpec spec;
        spec.qname = "spliced";
        spec.pos = 1000;
        // 40M 20D 40M covers 100 bp of reference but only 80 aligned bases.
        spec.cigar = {cigar_op(BAM_CMATCH, 40), cigar_op(BAM_CDEL, 20),
                      cigar_op(BAM_CMATCH, 40)};
        reads.push_back(spec);

        ReadSpec clipped;
        clipped.qname = "clipped";
        clipped.pos = 2000;
        // Soft-clipped bases consume query but not reference: 50 bp covered.
        clipped.cigar = {cigar_op(BAM_CSOFT_CLIP, 20), cigar_op(BAM_CMATCH, 50)};
        reads.push_back(clipped);

        TempBam bam(contigs(), reads);
        BamReader reader(bam.path(), 1);
        FilterOptions filters;

        std::vector<BinTally> bins(10);
        reader.query_bins({"chr1", 1000, 2000}, bins, filters);
        double total_bases = 0.0;
        for (const auto& bin : bins) total_bases += bin.base_count;
        CHECK_NEAR(total_bases, 80.0, 1e-9);

        std::vector<BinTally> clipped_bins(1);
        reader.query_bins({"chr1", 2000, 2100}, clipped_bins, filters);
        CHECK_NEAR(clipped_bins[0].base_count, 50.0, 1e-9);
    });

    suite.add("paired fragments are counted once, end to end", [] {
        // A 300 bp fragment from two 100 bp mates at 3000 and 3200.
        std::vector<ReadSpec> reads;
        ReadSpec left;
        left.qname = "pair1";
        left.pos = 3000;
        left.flag = flags::kPaired | flags::kProperPair | flags::kRead1;
        left.cigar = {cigar_op(BAM_CMATCH, 100)};
        left.mtid = 0;
        left.mpos = 3200;
        left.isize = 300;
        reads.push_back(left);

        ReadSpec right;
        right.qname = "pair1";
        right.pos = 3200;
        right.flag = flags::kPaired | flags::kProperPair | flags::kRead2 | flags::kReverse;
        right.cigar = {cigar_op(BAM_CMATCH, 100)};
        right.mtid = 0;
        right.mpos = 3000;
        right.isize = -300;
        reads.push_back(right);

        TempBam bam(contigs(), reads);
        BamReader reader(bam.path(), 1);

        FilterOptions fragment_mode;  // extend_to_fragment is on by default
        CHECK_EQ(reader.count_reads({"chr1", 3000, 3300}, fragment_mode),
                 std::uint64_t{1});

        // The gap between the mates is covered when the fragment is extended.
        std::vector<BinTally> bins(3);
        reader.query_bins({"chr1", 3000, 3300}, bins, fragment_mode);
        CHECK_NEAR(bins[0].base_count, 100.0, 1e-9);
        CHECK_NEAR(bins[1].base_count, 100.0, 1e-9);  // the unsequenced middle
        CHECK_NEAR(bins[2].base_count, 100.0, 1e-9);

        // Without extension both mates score, and the middle stays empty.
        FilterOptions read_mode;
        read_mode.extend_to_fragment = false;
        CHECK_EQ(reader.count_reads({"chr1", 3000, 3300}, read_mode), std::uint64_t{2});
        std::vector<BinTally> read_bins(3);
        reader.query_bins({"chr1", 3000, 3300}, read_bins, read_mode);
        CHECK_NEAR(read_bins[1].base_count, 0.0, 1e-9);
    });

    suite.add("fragments starting before the window still contribute", [] {
        // The mate pair starts at 2900 but reaches into a window opening at
        // 3100; a naive iterator query would miss it entirely.
        std::vector<ReadSpec> reads;
        ReadSpec left;
        left.qname = "reach";
        left.pos = 2900;
        left.flag = flags::kPaired | flags::kProperPair | flags::kRead1;
        left.cigar = {cigar_op(BAM_CMATCH, 50)};
        left.mtid = 0;
        left.mpos = 3150;
        left.isize = 300;  // fragment 2900-3200
        reads.push_back(left);

        ReadSpec right;
        right.qname = "reach";
        right.pos = 3150;
        right.flag = flags::kPaired | flags::kProperPair | flags::kRead2 | flags::kReverse;
        right.cigar = {cigar_op(BAM_CMATCH, 50)};
        right.mtid = 0;
        right.mpos = 2900;
        right.isize = -300;
        reads.push_back(right);

        TempBam bam(contigs(), reads);
        BamReader reader(bam.path(), 1);
        FilterOptions filters;

        std::vector<BinTally> bins(1);
        reader.query_bins({"chr1", 3100, 3200}, bins, filters);
        CHECK_NEAR(bins[0].read_count, 1.0, 1e-9);
        CHECK_NEAR(bins[0].base_count, 100.0, 1e-9);
    });

    suite.add("over-long fragments are rejected", [] {
        std::vector<ReadSpec> reads;
        ReadSpec left;
        left.qname = "long";
        left.pos = 4000;
        left.flag = flags::kPaired | flags::kProperPair | flags::kRead1;
        left.cigar = {cigar_op(BAM_CMATCH, 100)};
        left.mtid = 0;
        left.mpos = 6000;
        left.isize = 2100;
        reads.push_back(left);

        TempBam bam(contigs(), reads);
        BamReader reader(bam.path(), 1);

        FilterOptions gated;  // max_fragment_length defaults to 1000
        CHECK_EQ(reader.count_reads({"chr1", 4000, 4200}, gated), std::uint64_t{0});

        FilterOptions ungated;
        ungated.max_fragment_length = 0;
        CHECK_EQ(reader.count_reads({"chr1", 4000, 4200}, ungated), std::uint64_t{1});
    });

    suite.add("concurrent queries agree with serial ones", [] {
        TempBam bam(contigs(), tiling_reads());
        BamReader reader(bam.path(), 4);
        FilterOptions filters;

        std::vector<BinTally> expected(10);
        reader.query_bins({"chr1", 1000, 2000}, expected, filters);

        constexpr int kThreads = 8;
        constexpr int kIterations = 40;
        std::vector<std::thread> workers;
        std::atomic<int> mismatches{0};

        for (int t = 0; t < kThreads; ++t) {
            workers.emplace_back([&] {
                for (int n = 0; n < kIterations; ++n) {
                    std::vector<BinTally> bins(10);
                    reader.query_bins({"chr1", 1000, 2000}, bins, filters);
                    for (std::size_t i = 0; i < bins.size(); ++i) {
                        if (bins[i].read_count != expected[i].read_count ||
                            bins[i].base_count != expected[i].base_count) {
                            mismatches.fetch_add(1);
                        }
                    }
                }
            });
        }
        for (auto& worker : workers) worker.join();
        CHECK_EQ(mismatches.load(), 0);
    });

    // Regression: BamReader::tid_for() used to call sam_hdr_name2tid() directly
    // on the shared header. That function lazily builds the header's hrecs index
    // on first use, without synchronisation, so concurrent callers each built one
    // and orphaned all but the last -- an unfreeable ~160 KB per losing thread
    // that ASan reported as a leak out of sam_hdr_fill_hrecs. Names that hit the
    // constructor-built table never reach htslib now; names that miss (unknown
    // contigs, AN: aliases) still do, but only under a lock. Hammer both kinds
    // at once. This asserts nothing a leak checker would not; it exists so the
    // ASan and TSan jobs have something to catch.
    suite.add("concurrent name lookups do not race on the lazy header index", [] {
        TempBam bam({testing::Contig{"chr1", kChromLength, {}},
                     testing::Contig{"chr2", kChromLength, {}}},
                    tiling_reads());
        BamReader reader(bam.path(), 1);

        constexpr int kThreads = 8;
        constexpr int kIterations = 50;
        std::vector<std::thread> workers;
        std::atomic<int> wrong{0};

        // Rendezvous before the first lookup. Without it this test is worthless:
        // spawning eight threads takes long enough that the first one is done
        // filling hrecs before the last one starts, so the defect goes
        // unobserved. Verified by reverting the fix -- unbarriered, this suite
        // passed; barriered, it aborts on the leak.
        std::atomic<int> arrived{0};
        std::atomic<bool> go{false};

        for (int t = 0; t < kThreads; ++t) {
            workers.emplace_back([&] {
                arrived.fetch_add(1, std::memory_order_release);
                while (!go.load(std::memory_order_acquire)) {
                    std::this_thread::yield();
                }
                for (int n = 0; n < kIterations; ++n) {
                    if (reader.tid_for("chr1") != 0) wrong.fetch_add(1);
                    if (reader.tid_for("chr2") != 1) wrong.fetch_add(1);
                    // Misses share one name so every thread contends on the
                    // fallback path, which is where htslib is still reached.
                    if (reader.tid_for("chrUnplaced") >= 0) wrong.fetch_add(1);
                    if (reader.chromosome_length("chr1") != kChromLength) {
                        wrong.fetch_add(1);
                    }
                }
            });
        }
        while (arrived.load(std::memory_order_acquire) < kThreads) {
            std::this_thread::yield();
        }
        go.store(true, std::memory_order_release);

        for (auto& worker : workers) worker.join();
        CHECK_EQ(wrong.load(), 0);
    });

    // The reason tid_for() still falls through to htslib on a miss: an @SQ line
    // may carry AN: alternative names, and sam_hdr_name2tid() resolves those.
    // A BED written against Ensembl or RefSeq contig names is the real case.
    suite.add("@SQ AN: alternative names resolve to their contig", [] {
        TempBam bam({testing::Contig{"chr1", kChromLength, "1,NC_000001.11"},
                     testing::Contig{"chr2", kChromLength, {}}},
                    tiling_reads());
        BamReader reader(bam.path(), 1);

        CHECK_EQ(reader.tid_for("chr1"), 0);
        CHECK_EQ(reader.tid_for("1"), 0);
        CHECK_EQ(reader.tid_for("NC_000001.11"), 0);
        CHECK_EQ(reader.tid_for("chr2"), 1);
        CHECK(reader.tid_for("chrUnplaced") < 0);

        // Memoised aliases must keep answering, and keep answering the same.
        CHECK_EQ(reader.tid_for("1"), 0);
        CHECK_EQ(reader.chromosome_length("1"), kChromLength);
    });

    suite.add("a missing index is reported, not ignored", [] {
        bool threw = false;
        try {
            BamReader reader("/nonexistent/path/to/sample.bam", 1);
        } catch (const std::exception&) {
            threw = true;
        }
        CHECK(threw);
    });

    return suite.run();
}
