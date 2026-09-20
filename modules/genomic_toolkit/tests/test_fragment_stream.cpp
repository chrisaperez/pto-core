// SPDX-License-Identifier: MIT
//
// Text fragment sources, the size histogram, and the filters.
//
// The BAM path is exercised only when the build has htslib; the text paths are
// always exercised, which is why they carry the parser edge cases (CRLF, a
// final line with no newline, a line longer than the read buffer, partially
// numeric fields).

#include "toolkit/fragment_stream.hpp"

#include <string>
#include <vector>

#include "test_util.hpp"

using namespace toolkit;

namespace {

// Collects every emitted fragment. Real consumers never do this -- that is the
// point of the batch interface -- but a test needs the fragments to assert on.
struct Collector {
    std::vector<Fragment> fragments;
    std::vector<std::size_t> batch_sizes;

    BatchCallback callback() {
        return [this](std::span<Fragment> b) {
            batch_sizes.push_back(b.size());
            fragments.insert(fragments.end(), b.begin(), b.end());
        };
    }
};

}  // namespace

int main() {
    testing::Suite suite{"fragment_stream", {}};

    suite.add("bedpe pairs become one fragment spanning both mates", [] {
        const testing::TempFile f(".bedpe",
                                  "chr1\t100\t150\tchr1\t300\t350\tpair1\t60\t+\t-\n"
                                  "chr1\t500\t550\tchr1\t560\t610\tpair2\t60\t-\t+\n");
        ContigDict dict;
        Collector c;
        FragmentFilter filter;
        const StreamStats s = stream_bedpe(f.path(), filter, dict, c.callback());

        CHECK_EQ(s.records_read, std::uint64_t{2});
        CHECK_EQ(s.fragments_emitted, std::uint64_t{2});
        CHECK_EQ(c.fragments.size(), std::size_t{2});
        CHECK_EQ(c.fragments[0].start, Coord{100});
        CHECK_EQ(c.fragments[0].end, Coord{350});   // min(start), max(end)
        CHECK_EQ(c.fragments[0].length(), Coord{250});
        // Column 9 is strand1: '+' on the first pair, '-' on the second.
        CHECK(!c.fragments[0].reverse());
        CHECK(c.fragments[1].reverse());
        CHECK_EQ(c.fragments[1].start, Coord{500});
        CHECK_EQ(c.fragments[1].end, Coord{610});
        CHECK_EQ(dict.size(), std::size_t{1});
        CHECK_EQ(dict.name(c.fragments[0].tid), std::string("chr1"));
    });

    suite.add("inter-chromosomal and unmapped mates are counted, not emitted", [] {
        const testing::TempFile f(".bedpe",
                                  "chr1\t100\t150\tchr2\t300\t350\tx\t60\n"
                                  "chr1\t-1\t-1\tchr1\t300\t350\ty\t60\n"
                                  "chr1\t100\t150\tchr1\t300\t350\tz\t60\n");
        ContigDict dict;
        Collector c;
        const StreamStats s = stream_bedpe(f.path(), FragmentFilter{}, dict, c.callback());
        CHECK_EQ(s.dropped_interchrom, std::uint64_t{1});
        CHECK_EQ(s.dropped_unpaired, std::uint64_t{1});
        CHECK_EQ(s.fragments_emitted, std::uint64_t{1});
    });

    suite.add("length and mapq filters", [] {
        const testing::TempFile f(".bedpe",
                                  "chr1\t0\t10\tchr1\t40\t50\ta\t60\n"    // 50 bp
                                  "chr1\t0\t10\tchr1\t190\t200\tb\t60\n"  // 200 bp
                                  "chr1\t0\t10\tchr1\t1990\t2000\tc\t60\n" // 2000 bp
                                  "chr1\t0\t10\tchr1\t190\t200\td\t5\n");  // low mapq
        ContigDict dict;
        Collector c;
        FragmentFilter filter;
        filter.min_length = 100;   // the >= 100 bp gate from the ATAC recipe
        filter.max_length = 1000;
        filter.min_mapq = 30;
        const StreamStats s = stream_bedpe(f.path(), filter, dict, c.callback());

        CHECK_EQ(s.fragments_emitted, std::uint64_t{1});
        CHECK_EQ(s.dropped_length, std::uint64_t{2});  // 50 bp and 2000 bp
        CHECK_EQ(s.dropped_mapq, std::uint64_t{1});
        CHECK_EQ(c.fragments.size(), std::size_t{1});
        CHECK_EQ(c.fragments[0].length(), Coord{200});
    });

    suite.add("malformed lines are counted, not fatal", [] {
        const testing::TempFile f(".bedpe",
                                  "# a comment\n"
                                  "track name=whatever\n"
                                  "chr1\t100\n"                              // too few columns
                                  "chr1\tabc\t150\tchr1\t300\t350\tx\t60\n"  // non-numeric
                                  "chr1\t100x\t150\tchr1\t300\t350\ty\t60\n" // numeric prefix
                                  "chr1\t100\t150\tchr1\t300\t350\tz\t60\n");
        ContigDict dict;
        Collector c;
        const StreamStats s = stream_bedpe(f.path(), FragmentFilter{}, dict, c.callback());
        CHECK_EQ(s.malformed_lines, std::uint64_t{3});
        CHECK_EQ(s.fragments_emitted, std::uint64_t{1});
    });

    suite.add("fragment BED honours the count column only when asked", [] {
        // 10x fragments.tsv column 5 is the number of supporting read pairs.
        const testing::TempFile f(".bed",
                                  "chr1\t100\t300\tAAACGA-1\t3\n"
                                  "chr1\t400\t600\tAAACGA-1\t1\n"
                                  "chr1\t700\t900\tAAACGA-1\n");  // absent => 1
        ContigDict dict;
        Collector c;
        FragmentFilter with_counts;
        with_counts.use_count_column = true;
        const StreamStats s = stream_fragment_bed(f.path(), with_counts, dict, c.callback());
        CHECK_EQ(s.records_read, std::uint64_t{3});
        CHECK_EQ(s.fragments_emitted, std::uint64_t{5});
        CHECK_EQ(c.fragments.size(), std::size_t{5});

        // The same 10x file without the flag counts unique fragments.
        ContigDict d2;
        Collector c2;
        const StreamStats plain = stream_fragment_bed(f.path(), FragmentFilter{}, d2,
                                                      c2.callback());
        CHECK_EQ(plain.fragments_emitted, std::uint64_t{3});
        CHECK_EQ(plain.malformed_lines, std::uint64_t{0});
    });

    suite.add("a BED5/BED6 score column is one fragment per interval by default", [] {
        // Before 2026-09-11 column 5 was always a count: the first record below
        // was 60 fragments and every size, duplicate and FRiP statistic was
        // multiplied by the score, exit 0.
        const testing::TempFile f(".bed",
                                  "chr1\t100\t300\tread1\t60\t+\n"   // BED6, MAPQ-like score
                                  "chr1\t400\t600\tread2\t0\t-\n"    // score 0
                                  "chr1\t700\t900\tread3\t.\t+\n"    // score absent
                                  "chr1\t1000\t1200\tpeak\t1000\n"); // BED5
        ContigDict dict;
        Collector c;
        const StreamStats s = stream_fragment_bed(f.path(), FragmentFilter{}, dict,
                                                  c.callback());
        CHECK_EQ(s.records_read, std::uint64_t{4});
        CHECK_EQ(s.fragments_emitted, std::uint64_t{4});
        CHECK_EQ(s.malformed_lines, std::uint64_t{0});
        CHECK_EQ(s.sizes.total(), std::uint64_t{4});

        // The dispatching entry point every subcommand uses behaves the same.
        ContigDict d3;
        Collector c3;
        const StreamStats via_dispatch = stream_fragments(f.path(), FragmentFilter{}, d3,
                                                          c3.callback());
        CHECK_EQ(via_dispatch.fragments_emitted, std::uint64_t{4});

        // Asking for counts on this file is a claim that column 5 IS a count,
        // and the count rule then applies: "0" and "." are corrupt counts.
        FragmentFilter with_counts;
        with_counts.use_count_column = true;
        ContigDict d2;
        Collector c2;
        const StreamStats counted = stream_fragment_bed(f.path(), with_counts, d2,
                                                        c2.callback());
        CHECK_EQ(counted.fragments_emitted, std::uint64_t{60 + 1000});
        CHECK_EQ(counted.malformed_lines, std::uint64_t{2});
    });

    suite.add("size histogram statistics", [] {
        std::string body;
        // 100 fragments of 200 bp, 100 of 50 bp, 1 of 5000 bp.
        for (int i = 0; i < 100; ++i) {
            body += "chr1\t" + std::to_string(i * 10) + "\t" +
                    std::to_string(i * 10 + 200) + "\n";
        }
        for (int i = 0; i < 100; ++i) {
            body += "chr2\t" + std::to_string(i * 10) + "\t" +
                    std::to_string(i * 10 + 50) + "\n";
        }
        body += "chr3\t0\t5000\n";

        const testing::TempFile f(".bed", body);
        ContigDict dict;
        Collector c;
        FragmentFilter filter;
        filter.max_length = 0;  // keep the 5 kb outlier
        const StreamStats s = stream_fragment_bed(f.path(), filter, dict, c.callback());

        CHECK_EQ(s.sizes.total(), std::uint64_t{201});
        CHECK_EQ(s.sizes.count(200), std::uint64_t{100});
        CHECK_EQ(s.sizes.count(50), std::uint64_t{100});
        CHECK_EQ(s.sizes.overflow(), std::uint64_t{1});  // beyond kMaxTrackedLength
        CHECK_EQ(s.sizes.median(), Coord{200});
        CHECK_NEAR(s.sizes.fraction_in(180, 247), 100.0 / 201.0, 1e-9);
        CHECK_NEAR(s.sizes.fraction_in(1, 99), 100.0 / 201.0, 1e-9);
        // The mean includes the outlier's true length even though it is not
        // binned -- otherwise the overflow bucket would quietly bias it down.
        CHECK_NEAR(s.sizes.mean(), (100.0 * 200 + 100.0 * 50 + 5000.0) / 201.0, 1e-6);
    });

    suite.add("batching: full batches then a short tail", [] {
        std::string body;
        const int n = static_cast<int>(kBatchFragments) + 7;
        for (int i = 0; i < n; ++i) {
            body += "chr1\t" + std::to_string(i) + "\t" + std::to_string(i + 100) + "\n";
        }
        const testing::TempFile f(".bed", body);
        ContigDict dict;
        Collector c;
        stream_fragment_bed(f.path(), FragmentFilter{}, dict, c.callback());
        CHECK_EQ(c.fragments.size(), static_cast<std::size_t>(n));
        CHECK_EQ(c.batch_sizes.size(), std::size_t{2});
        CHECK_EQ(c.batch_sizes[0], kBatchFragments);
        CHECK_EQ(c.batch_sizes[1], std::size_t{7});
    });

    suite.add("CRLF, a missing final newline, and a very long line", [] {
        // The long line is padded past the reader's 1 MiB chunk, so the buffer
        // has to grow mid-line. A fixed buffer would split it into two
        // malformed records and still report a plausible-looking count.
        const std::string padding(1'200'000, 'N');
        const testing::TempFile f(".bed",
                                  "chr1\t100\t300\r\n"
                                  "chr1\t400\t600\tname\t1\t" + padding + "\n"
                                  "chr1\t700\t900");  // no trailing newline
        ContigDict dict;
        Collector c;
        const StreamStats s = stream_fragment_bed(f.path(), FragmentFilter{}, dict,
                                                  c.callback());
        CHECK_EQ(s.malformed_lines, std::uint64_t{0});
        CHECK_EQ(s.fragments_emitted, std::uint64_t{3});
        CHECK_EQ(c.fragments[0].end, Coord{300});  // \r did not leak into the field
        CHECK_EQ(c.fragments[2].start, Coord{700});
    });

    suite.add("coordinates past the int32 ceiling are rejected, not truncated", [] {
        const testing::TempFile f(".bed", "chr1\t0\t3000000000\nchr1\t10\t210\n");
        ContigDict dict;
        Collector c;
        FragmentFilter filter;
        filter.max_length = 0;
        const StreamStats s = stream_fragment_bed(f.path(), filter, dict, c.callback());
        CHECK_EQ(s.malformed_lines, std::uint64_t{1});
        CHECK_EQ(s.fragments_emitted, std::uint64_t{1});
    });

    suite.add("dispatch by extension, and a clear error for unsupported input", [] {
        const testing::TempFile bed(".bed", "chr1\t100\t300\n");
        ContigDict dict;
        Collector c;
        const StreamStats s = stream_fragments(bed.path(), FragmentFilter{}, dict,
                                               c.callback());
        CHECK_EQ(s.fragments_emitted, std::uint64_t{1});

        if constexpr (!kHaveBam) {
            // Without htslib a .bam input must fail loudly rather than being
            // handed to the text parser, which would report every line as
            // malformed and a fragment count of zero.
            ContigDict d2;
            Collector c2;
            CHECK_THROWS(std::runtime_error,
                         stream_fragments("nonexistent.bam", FragmentFilter{}, d2,
                                          c2.callback()));
        }
    });

    suite.add("a missing file is an error", [] {
        ContigDict dict;
        Collector c;
        CHECK_THROWS(std::runtime_error,
                     stream_fragment_bed("/nonexistent/path/x.bed", FragmentFilter{},
                                         dict, c.callback()));
    });

    suite.add("contradictory filter bounds are rejected up front", [] {
        const testing::TempFile f(".bed", "chr1\t100\t300\n");
        ContigDict dict;
        Collector c;
        FragmentFilter filter;
        filter.min_length = 500;
        filter.max_length = 100;
        CHECK_THROWS(std::invalid_argument,
                     stream_fragment_bed(f.path(), filter, dict, c.callback()));
    });

    // ---- 2026-09-11 audit (docs/AUDIT_2026-09-11_genomic_toolkit.md) ----

    suite.add("a BEDPE mate with an inverted interval is malformed, not a fragment", [] {
        // min(start)/max(end) used to build a plausible span out of a mate
        // whose own end precedes its start: "150 100" became part of a 150 bp
        // fragment with 0 malformed.
        const testing::TempFile f(".bedpe",
                                  "chr1\t150\t100\tchr1\t250\t300\ta\t60\n"   // mate 1 inverted
                                  "chr1\t100\t150\tchr1\t300\t250\tb\t60\n"   // mate 2 inverted
                                  "chr1\t100\t150\tchr1\t250\t300\tc\t60\n"); // valid
        ContigDict dict;
        Collector c;
        const StreamStats s = stream_bedpe(f.path(), FragmentFilter{}, dict, c.callback());
        CHECK_EQ(s.malformed_lines, std::uint64_t{2});
        CHECK_EQ(s.fragments_emitted, std::uint64_t{1});
    });

    suite.add("a BEDPE MAPQ outside 0-255 is malformed, not a filtering decision", [] {
        // "-7" used to land in dropped_mapq, reporting a corrupt record as a
        // low-quality alignment; "256" passed any threshold.
        const testing::TempFile f(".bedpe",
                                  "chr1\t100\t150\tchr1\t250\t300\ta\t-7\n"
                                  "chr1\t100\t150\tchr1\t250\t300\tb\t256\n"
                                  "chr1\t100\t150\tchr1\t250\t300\tc\t255\n");
        ContigDict dict;
        Collector c;
        FragmentFilter filter;
        filter.min_mapq = 30;
        const StreamStats s = stream_bedpe(f.path(), filter, dict, c.callback());
        CHECK_EQ(s.malformed_lines, std::uint64_t{2});
        CHECK_EQ(s.dropped_mapq, std::uint64_t{0});
        CHECK_EQ(s.fragments_emitted, std::uint64_t{1});
    });

    suite.add("BED is 0-based half-open: 1 bp, empty and inverted intervals", [] {
        const testing::TempFile f(".bed",
                                  "chr1\t0\t1\n"       // one base, at the very start
                                  "chr1\t100\t100\n"   // empty: no base, not corrupt
                                  "chr1\t100\t99\n");  // inverted: corrupt
        ContigDict dict;
        Collector c;
        const StreamStats s = stream_fragment_bed(f.path(), FragmentFilter{}, dict,
                                                  c.callback());
        CHECK_EQ(s.records_read, std::uint64_t{3});
        CHECK_EQ(s.fragments_emitted, std::uint64_t{1});
        CHECK_EQ(s.dropped_length, std::uint64_t{1});
        CHECK_EQ(s.malformed_lines, std::uint64_t{1});
        if (c.fragments.size() == 1) {
            CHECK_EQ(c.fragments[0].start, Coord{0});
            CHECK_EQ(c.fragments[0].length(), Coord{1});
        }
    });

    suite.add("unsorted and overlapping fragments are all emitted, in input order", [] {
        // The stream imposes no order; order is the duplicate marker's
        // business, and the CLI refuses on its behalf (scripts/test_cli.sh).
        const testing::TempFile f(".bed",
                                  "chr2\t500\t700\n"
                                  "chr1\t100\t400\n"
                                  "chr1\t50\t450\n"    // contains the previous one
                                  "chr1\t100\t400\n"); // identical to the second
        ContigDict dict;
        Collector c;
        const StreamStats s = stream_fragment_bed(f.path(), FragmentFilter{}, dict,
                                                  c.callback());
        CHECK_EQ(s.fragments_emitted, std::uint64_t{4});
        CHECK_EQ(s.malformed_lines, std::uint64_t{0});
        if (c.fragments.size() == 4) {
            CHECK_EQ(dict.name(c.fragments[0].tid), std::string("chr2"));
            CHECK_EQ(c.fragments[2].start, Coord{50});
        }
    });

    suite.add("the size histogram's length sum throws instead of wrapping", [] {
        // Reachable through the copy column: `sizes` over 9,300 crafted lines
        // reaches it in 41 s (see checked_add). Driven here by doubling one
        // max-length fragment with
        // merge(): after k merges the sum is (2^31 - 1) * 2^k, which fits for
        // k = 33 and not for k = 34. Before the fix the 34th merge wrapped to a
        // small number and mean() reported a plausible length.
        SizeHistogram h;
        h.add(static_cast<Coord>(kMaxCoord));
        int merges = 0;
        bool threw = false;
        for (; merges < 64; ++merges) {
            const SizeHistogram copy = h;
            const std::uint64_t total_before = h.total();
            const double mean_before = h.mean();
            try {
                h.merge(copy);
            } catch (const std::overflow_error&) {
                threw = true;
                CHECK_EQ(h.total(), total_before);  // a refused merge changes nothing
                CHECK_NEAR(h.mean(), mean_before, 1e-3);
                break;
            }
        }
        CHECK(threw);
        CHECK_EQ(merges, 33);
    });

    return suite.run();
}
