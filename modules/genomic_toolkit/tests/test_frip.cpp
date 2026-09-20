// SPDX-License-Identifier: MIT
//
// FRiP end to end: peak BED -> merged cover -> streaming count, plus the two
// failure modes that produce a wrong number silently rather than an error
// (double counting against overlapping peaks, and a chromosome-naming
// mismatch).

#include "toolkit/frip.hpp"

#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "toolkit/dup_marker.hpp"
#include "toolkit/fragment_stream.hpp"
#include "test_util.hpp"

using namespace toolkit;

namespace {

FripResult run_frip(const std::string& frag_path, const std::string& peak_path,
                    bool skip_duplicates, StreamStats* out_stats = nullptr) {
    PeakSet peaks = PeakSet::build(read_peaks_bed(peak_path));
    ContigDict dict;
    for (const auto& name : peaks.dict().names()) dict.intern(name);
    peaks.bind(dict);

    DuplicateMarker marker;
    FripCounter counter(peaks, skip_duplicates);
    const StreamStats s =
        stream_fragments(frag_path, FragmentFilter{}, dict, [&](std::span<Fragment> b) {
            if (skip_duplicates) marker.observe_batch(b);
            counter.add_batch(b);
        });
    if (out_stats != nullptr) *out_stats = s;
    return counter.finish(dict);
}

}  // namespace

int main() {
    testing::Suite suite{"frip", {}};

    suite.add("counts fragments touching a peak", [] {
        const testing::TempFile peaks(".bed", "chr1\t1000\t2000\nchr1\t5000\t6000\n");
        const testing::TempFile frags(".bed",
                                      "chr1\t1100\t1300\n"   // inside peak 1
                                      "chr1\t900\t1050\n"    // straddles peak 1's start
                                      "chr1\t3000\t3200\n"   // between peaks
                                      "chr1\t5900\t6100\n"   // straddles peak 2's end
                                      "chr1\t900\t1000\n");  // ends exactly at the start
        const FripResult r = run_frip(frags.path(), peaks.path(), false);

        CHECK_EQ(r.fragments_total, std::uint64_t{5});
        CHECK_EQ(r.fragments_in_peaks, std::uint64_t{3});
        CHECK_NEAR(r.fragment_frip(), 3.0 / 5.0, 1e-12);

        // Base FRiP counts only the bases actually inside a peak: 200 from the
        // contained fragment, 50 from the left straddle, 100 from the right.
        CHECK_EQ(r.bases_in_peaks, std::int64_t{200 + 50 + 100});
        CHECK_EQ(r.bases_total, std::int64_t{200 + 150 + 200 + 200 + 100});
        CHECK_NEAR(r.base_frip(), 350.0 / 850.0, 1e-12);
    });

    suite.add("a fragment spanning overlapping peaks counts once", [] {
        // The bedtools-and-wc recipe over an unmerged peak file counts this
        // fragment three times and reports a FRiP above 1. Merging at build
        // time is what makes the statistic well defined.
        const testing::TempFile peaks(".bed",
                                      "chr1\t1000\t1500\n"
                                      "chr1\t1200\t1700\n"
                                      "chr1\t1600\t1800\n");
        const testing::TempFile frags(".bed", "chr1\t1000\t1800\n");
        const FripResult r = run_frip(frags.path(), peaks.path(), false);

        CHECK_EQ(r.peaks_input, std::size_t{3});
        CHECK_EQ(r.peaks_merged, std::size_t{1});
        CHECK_EQ(r.peak_bases, std::int64_t{800});
        CHECK_EQ(r.fragments_in_peaks, std::uint64_t{1});
        CHECK_NEAR(r.fragment_frip(), 1.0, 1e-12);
        CHECK_EQ(r.bases_in_peaks, std::int64_t{800});
    });

    suite.add("duplicates are excluded from both numerator and denominator", [] {
        const testing::TempFile peaks(".bed", "chr1\t1000\t2000\n");
        const testing::TempFile frags(".bed",
                                      "chr1\t1100\t1300\n"
                                      "chr1\t1100\t1300\n"   // duplicate, in peak
                                      "chr1\t1100\t1300\n"   // duplicate, in peak
                                      "chr1\t3000\t3200\n"
                                      "chr1\t3000\t3200\n"); // duplicate, out of peak
        const FripResult kept = run_frip(frags.path(), peaks.path(), false);
        CHECK_EQ(kept.fragments_total, std::uint64_t{5});
        CHECK_EQ(kept.fragments_in_peaks, std::uint64_t{3});

        const FripResult deduped = run_frip(frags.path(), peaks.path(), true);
        CHECK_EQ(deduped.duplicates_excluded, std::uint64_t{3});
        CHECK_EQ(deduped.fragments_total, std::uint64_t{2});
        CHECK_EQ(deduped.fragments_in_peaks, std::uint64_t{1});
        CHECK_NEAR(deduped.fragment_frip(), 0.5, 1e-12);
    });

    suite.add("a naming mismatch is visible, not a silent zero", [] {
        const testing::TempFile peaks(".bed", "1\t1000\t2000\n");      // Ensembl style
        const testing::TempFile frags(".bed", "chr1\t1100\t1300\n");   // UCSC style
        const FripResult r = run_frip(frags.path(), peaks.path(), false);

        CHECK_EQ(r.fragments_in_peaks, std::uint64_t{0});
        CHECK_EQ(r.fragment_frip(), 0.0);
        // The diagnostic that separates this from a genuinely dead library.
        CHECK_EQ(r.unmatched_contigs.size(), std::size_t{1});
        CHECK_EQ(r.unmatched_contigs.front(), std::string("chr1"));
    });

    suite.add("a contig with no peaks is not a mismatch", [] {
        const testing::TempFile peaks(".bed", "chr1\t1000\t2000\n");
        const testing::TempFile frags(".bed", "chr1\t1100\t1300\nchrM\t50\t250\n");
        const FripResult r = run_frip(frags.path(), peaks.path(), false);
        CHECK_EQ(r.fragments_total, std::uint64_t{2});
        CHECK_EQ(r.fragments_in_peaks, std::uint64_t{1});
        // chrM is genuinely absent from the peak file, which is normal; the
        // count is still right, so this is reported without being an error.
        CHECK_EQ(r.unmatched_contigs.size(), std::size_t{1});
    });

    suite.add("enrichment over uniform", [] {
        const testing::TempFile peaks(".bed", "chr1\t0\t100\n");
        const testing::TempFile frags(".bed", "chr1\t10\t60\nchr1\t500\t550\n");
        const FripResult r = run_frip(frags.path(), peaks.path(), false);
        CHECK_NEAR(r.fragment_frip(), 0.5, 1e-12);
        // Peaks cover 100 of 1000 bp, so uniform scatter would give 0.1.
        CHECK_NEAR(r.enrichment(1000), 5.0, 1e-9);
        CHECK_EQ(r.enrichment(0), 0.0);  // no genome size => no claim
    });

    suite.add("multi-contig, many fragments, against a brute-force count", [] {
        std::string peak_body;
        std::string frag_body;
        // Peaks every 1000 bp, 200 bp wide, on two contigs.
        for (int i = 0; i < 50; ++i) {
            peak_body += "chr1\t" + std::to_string(i * 1000) + "\t" +
                         std::to_string(i * 1000 + 200) + "\n";
            peak_body += "chr2\t" + std::to_string(i * 1000) + "\t" +
                         std::to_string(i * 1000 + 200) + "\n";
        }
        // Fragments every 250 bp: exactly one in four starts inside a peak.
        std::uint64_t expected = 0;
        for (int i = 0; i < 400; ++i) {
            const int start = i * 250;
            const int end = start + 100;
            frag_body += "chr1\t" + std::to_string(start) + "\t" + std::to_string(end) + "\n";
            // A peak covers [k*1000, k*1000+200) for k < 50, so a 100 bp
            // fragment lands in one only when its start is a multiple of 1000
            // below 50000. Fragments past the last peak are the tail that
            // makes this a real check rather than a restatement of the loop.
            if (start < 50 * 1000 && start % 1000 < 200) ++expected;
        }
        const testing::TempFile peaks(".bed", peak_body);
        const testing::TempFile frags(".bed", frag_body);
        const FripResult r = run_frip(frags.path(), peaks.path(), false);

        CHECK_EQ(r.fragments_total, std::uint64_t{400});
        CHECK_EQ(r.fragments_in_peaks, expected);
        CHECK_EQ(r.peaks_merged, std::size_t{100});
        CHECK_EQ(r.peak_bases, std::int64_t{100 * 200});
    });

    suite.add("a per-thread overlap cursor survives a batch boundary mid-contig", [] {
        // FripCounter's OpenMP loop hands each thread a persistent
        // OverlapCursor (toolkit/peak_set.hpp) that must reset itself
        // whenever the resolved ContigPeaks changes -- and a batch is a fixed
        // number of LINES (kBatchFragments == 4096), with no relationship to
        // contig boundaries, so a real stream routinely hands one batch two
        // different contigs' worth of fragments, and the next batch a third.
        // Three contigs, 8000 fragments total (two full 4096-line batches),
        // sized so the chr1->chr2 transition falls inside batch 1 (line 3000
        // of 4096) and the chr2->chr3 transition falls inside batch 2 (line
        // 6000, i.e. line 1904 of the second batch).
        std::string peak_body;
        std::string frag_body;
        for (const char* chrom : {"chr1", "chr2", "chr3"}) {
            for (int i = 0; i < 50; ++i) {
                peak_body += std::string(chrom) + "\t" + std::to_string(i * 1000) + "\t" +
                             std::to_string(i * 1000 + 200) + "\n";
            }
        }
        std::uint64_t expected = 0;
        for (const auto& [chrom, n] :
             std::vector<std::pair<std::string, int>>{{"chr1", 3000}, {"chr2", 3000}, {"chr3", 2000}}) {
            for (int i = 0; i < n; ++i) {
                const int start = i * 250;
                const int end = start + 100;
                frag_body += chrom + "\t" + std::to_string(start) + "\t" + std::to_string(end) + "\n";
                if (start < 50 * 1000 && start % 1000 < 200) ++expected;
            }
        }
        const testing::TempFile peaks(".bed", peak_body);
        const testing::TempFile frags(".bed", frag_body);
        const FripResult r = run_frip(frags.path(), peaks.path(), false);

        CHECK_EQ(r.fragments_total, std::uint64_t{8000});
        CHECK_EQ(r.fragments_in_peaks, expected);
        CHECK_EQ(expected, std::uint64_t{150});  // 50 per contig, spot-checked by hand
        CHECK_EQ(r.peaks_merged, std::size_t{150});
        CHECK_EQ(r.peak_bases, std::int64_t{150 * 200});
    });

    suite.add("a malformed or empty peak file is fatal", [] {
        const testing::TempFile bad(".bed", "chr1\t1000\n");
        CHECK_THROWS(std::invalid_argument, read_peaks_bed(bad.path()));

        const testing::TempFile empty(".bed", "# only a comment\n");
        CHECK_THROWS(std::invalid_argument, read_peaks_bed(empty.path()));

        CHECK_THROWS(std::runtime_error, read_peaks_bed("/nonexistent/peaks.bed"));
    });

    // ---- 2026-09-11 audit (docs/AUDIT_2026-09-11_genomic_toolkit.md) ----

    suite.add("checked_add refuses to wrap or overflow, at the exact boundary", [] {
        // FRiP's base totals are signed; past INT64_MAX that is undefined
        // behaviour, reachable from 4,612 crafted lines (`frip` over 4,700 now
        // stops with the error after 142 s). The accumulator wiring in
        // FripCounter::add_batch cannot be driven that far in a unit test, so
        // the boundary is pinned on the helper itself and the wiring by that
        // measured run (docs/AUDIT_2026-09-11_genomic_toolkit.md).
        constexpr std::uint64_t kU = std::numeric_limits<std::uint64_t>::max();
        constexpr std::int64_t kI = std::numeric_limits<std::int64_t>::max();
        constexpr std::int64_t kIMin = std::numeric_limits<std::int64_t>::min();

        std::uint64_t u = kU - 5;
        checked_add(u, 5, "u");
        CHECK_EQ(u, kU);  // reaching the maximum is fine
        CHECK_THROWS(std::overflow_error, checked_add(u, 1, "u"));
        CHECK_EQ(u, kU);  // and a refused add leaves the total alone

        std::int64_t s = kI - 5;
        checked_add(s, 5, "s");
        CHECK_EQ(s, kI);
        CHECK_THROWS(std::overflow_error, checked_add(s, 1, "s"));
        CHECK_EQ(s, kI);

        std::int64_t neg = kIMin + 5;
        checked_add(neg, -5, "neg");
        CHECK_EQ(neg, kIMin);
        CHECK_THROWS(std::overflow_error, checked_add(neg, -1, "neg"));
        CHECK_EQ(neg, kIMin);

        std::int64_t mixed = kI;
        checked_add(mixed, -1, "mixed");  // moving away from the edge is not overflow
        CHECK_EQ(mixed, kI - 1);
    });

    suite.add("unsorted and overlapping peak records give the sorted, merged answer", [] {
        const testing::TempFile sorted(".bed", "chr1\t1000\t1800\nchr1\t5000\t6000\n");
        const testing::TempFile shuffled(".bed",
                                         "chr1\t5000\t6000\n"
                                         "chr1\t1600\t1800\n"
                                         "chr1\t1000\t1500\n"
                                         "chr1\t1200\t1700\n"
                                         "chr1\t5500\t5600\n");  // nested
        const testing::TempFile frags(".bed",
                                      "chr1\t900\t1100\nchr1\t1700\t1900\n"
                                      "chr1\t3000\t3100\nchr1\t5550\t5560\n");
        const FripResult a = run_frip(frags.path(), sorted.path(), false);
        const FripResult b = run_frip(frags.path(), shuffled.path(), false);
        CHECK_EQ(b.peaks_merged, std::size_t{2});
        CHECK_EQ(b.peak_bases, a.peak_bases);
        CHECK_EQ(b.fragments_in_peaks, a.fragments_in_peaks);
        CHECK_EQ(b.bases_in_peaks, a.bases_in_peaks);
        CHECK_EQ(a.fragments_in_peaks, std::uint64_t{3});
    });

    suite.add("a CRLF peak file gives the LF answer", [] {
        const testing::TempFile lf(".bed", "chr1\t1000\t2000\nchr1\t5000\t6000\n");
        const testing::TempFile crlf(".bed", "chr1\t1000\t2000\r\nchr1\t5000\t6000\r\n");
        const testing::TempFile frags(".bed", "chr1\t1900\t2100\nchr1\t5999\t6001\n");
        const FripResult a = run_frip(frags.path(), lf.path(), false);
        const FripResult b = run_frip(frags.path(), crlf.path(), false);
        CHECK_EQ(b.peak_bases, a.peak_bases);
        CHECK_EQ(b.fragments_in_peaks, a.fragments_in_peaks);
        CHECK_EQ(b.bases_in_peaks, a.bases_in_peaks);
        CHECK_EQ(b.bases_in_peaks, std::int64_t{100 + 1});
        CHECK_EQ(b.unmatched_contigs.size(), std::size_t{0});  // "\r" did not join a name
    });

    return suite.run();
}
