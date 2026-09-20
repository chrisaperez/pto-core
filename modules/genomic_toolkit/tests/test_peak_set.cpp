// SPDX-License-Identifier: MIT
//
// PeakSet: merging, overlap, base coverage.
//
// The core of this suite is a differential test against a brute-force
// reference that states the overlap rule the other way round -- as a bitmap
// over positions -- so the comparison is between two independent
// formulations rather than one expression against a restatement of itself.
// Every interval defect this module could have (the half-open off-by-one, a
// missed book-ended merge, a nested interval, a query starting inside a peak)
// is a disagreement between those two.

#include "toolkit/peak_set.hpp"

#include <random>
#include <vector>

#include "test_util.hpp"

using namespace toolkit;

namespace {

Region reg(std::string chrom, std::int64_t s, std::int64_t e) {
    Region r;
    r.chrom = std::move(chrom);
    r.start = s;
    r.end = e;
    return r;
}

// Reference: a per-base bitmap. Slow and obviously correct.
struct Bitmap {
    std::vector<char> covered;

    explicit Bitmap(std::size_t size) : covered(size, 0) {}

    void add(std::int64_t s, std::int64_t e) {
        for (std::int64_t i = s; i < e && i < static_cast<std::int64_t>(covered.size()); ++i) {
            covered[static_cast<std::size_t>(i)] = 1;
        }
    }
    [[nodiscard]] std::int64_t bases(std::int64_t s, std::int64_t e) const {
        std::int64_t n = 0;
        for (std::int64_t i = s; i < e && i < static_cast<std::int64_t>(covered.size()); ++i) {
            n += covered[static_cast<std::size_t>(i)];
        }
        return n;
    }
    [[nodiscard]] bool any(std::int64_t s, std::int64_t e) const { return bases(s, e) > 0; }
    [[nodiscard]] std::int64_t total() const {
        std::int64_t n = 0;
        for (char c : covered) n += c;
        return n;
    }
};

}  // namespace

int main() {
    testing::Suite suite{"peak_set", {}};

    suite.add("merges overlapping, nested and book-ended intervals", [] {
        const std::vector<Region> peaks = {
            reg("chr1", 100, 200),
            reg("chr1", 150, 250),  // overlaps the first
            reg("chr1", 160, 170),  // nested inside the union
            reg("chr1", 250, 300),  // book-ended: shares the boundary
            reg("chr1", 400, 500),  // disjoint
            reg("chr2", 0, 10),
        };
        const PeakSet ps = PeakSet::build(peaks);
        CHECK_EQ(ps.input_count(), std::size_t{6});
        CHECK_EQ(ps.size(), std::size_t{3});  // [100,300) [400,500) chr2:[0,10)
        CHECK_EQ(ps.covered_bases(), std::int64_t{200 + 100 + 10});

        // The merge is what stops a fragment spanning 190-260 counting twice.
        CHECK(ps.overlaps_named("chr1", 190, 260));
        CHECK(!ps.overlaps_named("chr1", 300, 400));  // the gap
        CHECK(!ps.overlaps_named("chr3", 0, 1000));   // unknown contig
    });

    suite.add("half-open boundaries", [] {
        const PeakSet ps = PeakSet::build({reg("chr1", 100, 200)});
        CHECK(!ps.overlaps_named("chr1", 0, 100));    // ends exactly at the start
        CHECK(ps.overlaps_named("chr1", 99, 101));
        CHECK(ps.overlaps_named("chr1", 199, 200));
        CHECK(!ps.overlaps_named("chr1", 200, 300));  // starts exactly at the end
        CHECK(!ps.overlaps_named("chr1", 150, 150));  // empty query covers nothing
    });

    suite.add("zero-length peaks cover no bases", [] {
        // Legal in BED as insertion points. They must not merge a gap shut and
        // must not make an empty query report a hit.
        const PeakSet ps = PeakSet::build({reg("chr1", 100, 100), reg("chr1", 200, 210)});
        CHECK_EQ(ps.size(), std::size_t{1});
        CHECK_EQ(ps.covered_bases(), std::int64_t{10});
        CHECK(!ps.overlaps_named("chr1", 90, 110));
    });

    suite.add("a set of only zero-length intervals is empty but counts its input", [] {
        // The library reports this; the CLI refuses it by name (scripts/
        // test_cli.sh). Before 2026-09-11 it reached frip's naming diagnostic,
        // because the dropped intervals never interned their contig.
        const PeakSet ps = PeakSet::build({reg("chr1", 100, 100), reg("chr1", 200, 200)});
        CHECK_EQ(ps.size(), std::size_t{0});
        CHECK_EQ(ps.input_count(), std::size_t{2});
        CHECK_EQ(ps.covered_bases(), std::int64_t{0});
        CHECK(!ps.overlaps_named("chr1", 0, 1000));
    });

    suite.add("rejects invalid intervals", [] {
        CHECK_THROWS(std::invalid_argument, (void)PeakSet::build({reg("chr1", -1, 10)}));
        CHECK_THROWS(std::invalid_argument, (void)PeakSet::build({reg("chr1", 50, 10)}));
        CHECK_THROWS(std::invalid_argument,
                     (void)PeakSet::build({reg("chr1", 0, kMaxCoord + 1)}));
    });

    suite.add("binds contigs and reports mismatches", [] {
        PeakSet ps = PeakSet::build({reg("chr1", 100, 200), reg("chr2", 0, 50)});
        ContigDict stream;
        const std::int32_t t1 = stream.intern("chr1");
        const std::int32_t tX = stream.intern("1");  // wrong convention
        ps.bind(stream);
        CHECK(ps.overlaps(t1, 150, 160));
        CHECK(!ps.overlaps(tX, 150, 160));
        const auto unmatched = stream.unmatched_against(ps.dict());
        CHECK_EQ(unmatched.size(), std::size_t{1});
        CHECK_EQ(unmatched.front(), std::string("1"));
    });

    suite.add("differential vs per-base bitmap", [] {
        constexpr std::int64_t kSpan = 4000;
        std::mt19937 rng(20260817);
        std::uniform_int_distribution<std::int64_t> pos(0, kSpan - 1);
        std::uniform_int_distribution<std::int64_t> len(0, 120);

        for (int trial = 0; trial < 40; ++trial) {
            std::vector<Region> peaks;
            Bitmap ref(static_cast<std::size_t>(kSpan));
            for (int i = 0; i < 60; ++i) {
                const std::int64_t s = pos(rng);
                const std::int64_t e = std::min(kSpan, s + len(rng));
                peaks.push_back(reg("chr1", s, e));
                ref.add(s, e);
            }
            const PeakSet ps = PeakSet::build(peaks);
            CHECK_EQ(ps.covered_bases(), ref.total());

            for (int q = 0; q < 400; ++q) {
                const std::int64_t s = pos(rng);
                const std::int64_t e = std::min(kSpan, s + len(rng));
                const bool got = ps.overlaps_named("chr1", static_cast<Coord>(s),
                                                   static_cast<Coord>(e));
                CHECK_EQ(got, ref.any(s, e));
            }
        }
    });

    // FRiP's fragment stream is coordinate-sorted, so overlap_bases() is
    // called with a monotone `qs` almost always -- the property an
    // OverlapCursor exploits to answer in O(1) amortised instead of a cold
    // O(log n) binary search per fragment. The cursor must be exactly as
    // correct as the plain overload in every query order, not merely the
    // sorted one: a driver can still go backwards (an out-of-order input,
    // which genomic_toolkit detects and reports elsewhere, not silently
    // mis-answers here), and correctness must not depend on which order
    // arrives. Checked against the SAME independent bitmap oracle the
    // existing tests above use, in three orders.
    suite.add("overlap_bases with a cursor matches the bitmap exactly, in every query order", [] {
        constexpr std::int64_t kSpan = 2000;
        std::mt19937 rng(99);
        std::uniform_int_distribution<std::int64_t> pos(0, kSpan - 1);
        std::uniform_int_distribution<std::int64_t> len(0, 300);

        std::vector<Region> peaks;
        Bitmap ref(static_cast<std::size_t>(kSpan));
        for (int i = 0; i < 30; ++i) {
            const std::int64_t s = pos(rng);
            const std::int64_t e = std::min(kSpan, s + len(rng));
            peaks.push_back(reg("chr1", s, e));
            ref.add(s, e);
        }
        PeakSet ps = PeakSet::build(peaks);
        ContigDict stream;
        const std::int32_t tid = stream.intern("chr1");
        ps.bind(stream);

        std::vector<std::pair<std::int64_t, std::int64_t>> queries;
        for (int q = 0; q < 2000; ++q) {
            const std::int64_t s = pos(rng);
            const std::int64_t e = std::min(kSpan, s + len(rng));
            queries.emplace_back(s, e);
        }

        for (int order = 0; order < 3; ++order) {
            auto qs_order = queries;
            if (order == 1) {
                std::sort(qs_order.begin(), qs_order.end());
            } else if (order == 2) {
                std::sort(qs_order.rbegin(), qs_order.rend());
            }
            OverlapCursor cursor;
            for (const auto& [s, e] : qs_order) {
                const std::int64_t got = ps.overlap_bases(
                    tid, static_cast<Coord>(s), static_cast<Coord>(e), cursor);
                CHECK_EQ(got, ref.bases(s, e));
            }
        }
    });

    suite.add("the cursor resets across a contig change, forwards or backwards", [] {
        // Two contigs, one shared cursor. A cursor that carried the wrong
        // contig's `i`/`last_q` forward would index into the wrong ContigPeaks
        // entirely -- a silently wrong answer, not a crash, since both arrays
        // are in bounds for their own contig.
        std::vector<Region> peaks = {
            reg("chr1", 100, 200),
            reg("chr1", 500, 600),
            reg("chr2", 50, 90),
            reg("chr2", 300, 999),
        };
        PeakSet ps = PeakSet::build(peaks);
        ContigDict stream;
        const std::int32_t t1 = stream.intern("chr1");
        const std::int32_t t2 = stream.intern("chr2");
        ps.bind(stream);

        OverlapCursor cursor;
        // Advance the cursor deep into chr1's interval list...
        CHECK_EQ(ps.overlap_bases(t1, 550, 600, cursor), 50);
        // ...then jump to chr2, at coordinates that would be nonsense if the
        // cursor incorrectly carried chr1's position/index forward.
        CHECK_EQ(ps.overlap_bases(t2, 0, 100, cursor), 40);  // [50,90) clipped to [0,100)
        CHECK_EQ(ps.overlap_bases(t2, 300, 400, cursor), 100);
        // Back to chr1, at a coordinate EARLIER than the cursor's last chr1
        // use -- this must be treated as a fresh contig, not as "chr1 going
        // backwards" (which would also reset, but via the wrong branch).
        CHECK_EQ(ps.overlap_bases(t1, 0, 150, cursor), 50);
    });

    suite.add("overlap_bases matches the bitmap exactly", [] {
        constexpr std::int64_t kSpan = 2000;
        std::mt19937 rng(7);
        std::uniform_int_distribution<std::int64_t> pos(0, kSpan - 1);
        std::uniform_int_distribution<std::int64_t> len(0, 300);

        std::vector<Region> peaks;
        Bitmap ref(static_cast<std::size_t>(kSpan));
        for (int i = 0; i < 30; ++i) {
            const std::int64_t s = pos(rng);
            const std::int64_t e = std::min(kSpan, s + len(rng));
            peaks.push_back(reg("chr1", s, e));
            ref.add(s, e);
        }
        PeakSet ps = PeakSet::build(peaks);
        ContigDict stream;
        const std::int32_t tid = stream.intern("chr1");
        ps.bind(stream);

        for (int q = 0; q < 3000; ++q) {
            const std::int64_t s = pos(rng);
            const std::int64_t e = std::min(kSpan, s + len(rng));
            CHECK_EQ(ps.overlap_bases(tid, static_cast<Coord>(s), static_cast<Coord>(e)),
                     ref.bases(s, e));
        }
    });

    return suite.run();
}
