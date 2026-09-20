// SPDX-License-Identifier: MIT
//
// DuplicateMarker: both modes, and the agreement between them.
//
// The central claim of the sorted-stream mode is that keeping state for one
// start position at a time loses nothing on coordinate-sorted input. That is
// checked here the only way it can be: by running the O(1)-memory mode and the
// O(n)-memory mode over the same shuffled-then-sorted fragment set and
// requiring identical duplicate counts.

#include "toolkit/dup_marker.hpp"

#include <algorithm>
#include <random>
#include <vector>

#include "test_util.hpp"

using namespace toolkit;

namespace {

Fragment frag(std::int32_t tid, Coord start, Coord end, bool reverse = false) {
    Fragment f;
    f.tid = tid;
    f.start = start;
    f.end = end;
    if (reverse) f.flags |= fragflags::kReverse;
    return f;
}

void sort_by_coord(std::vector<Fragment>& v) {
    std::sort(v.begin(), v.end(), [](const Fragment& a, const Fragment& b) {
        if (a.tid != b.tid) return a.tid < b.tid;
        if (a.start != b.start) return a.start < b.start;
        return a.end < b.end;
    });
}

}  // namespace

int main() {
    testing::Suite suite{"dup_marker", {}};

    suite.add("identical fragments after the first are duplicates", [] {
        DuplicateMarker m;
        Fragment a = frag(0, 100, 300);
        Fragment b = frag(0, 100, 300);
        Fragment c = frag(0, 100, 300);
        CHECK(!m.observe(a));
        CHECK(m.observe(b));
        CHECK(m.observe(c));
        CHECK(!a.duplicate());
        CHECK(b.duplicate());
        CHECK_EQ(m.stats().duplicates, std::uint64_t{2});
        CHECK_EQ(m.stats().unique, std::uint64_t{1});
        CHECK_NEAR(m.stats().duplicate_rate(), 2.0 / 3.0, 1e-12);
    });

    suite.add("same start, different end or strand is not a duplicate", [] {
        DuplicateMarker m;
        Fragment a = frag(0, 100, 300);
        Fragment b = frag(0, 100, 301);         // different 3' end
        Fragment c = frag(0, 100, 300, true);   // different orientation
        CHECK(!m.observe(a));
        CHECK(!m.observe(b));
        CHECK(!m.observe(c));
        CHECK_EQ(m.stats().duplicates, std::uint64_t{0});
        // Three distinct keys alive at one position; this is the number that
        // bounds the mode's memory.
        CHECK_EQ(m.stats().peak_live_keys, std::uint64_t{3});
    });

    suite.add("state is released when the start position advances", [] {
        DuplicateMarker m;
        Fragment a = frag(0, 100, 300);
        Fragment b = frag(0, 101, 301);
        Fragment c = frag(0, 102, 302);
        m.observe(a);
        m.observe(b);
        m.observe(c);
        // One key alive at a time: this is the whole memory argument, so it is
        // asserted rather than described.
        CHECK_EQ(m.stats().peak_live_keys, std::uint64_t{1});
        CHECK_EQ(m.stats().duplicates, std::uint64_t{0});
    });

    suite.add("detects unsorted input instead of silently undercounting", [] {
        DuplicateMarker m;
        Fragment a = frag(0, 500, 600);
        Fragment b = frag(0, 100, 200);  // goes backwards
        m.observe(a);
        m.observe(b);
        CHECK_EQ(m.stats().out_of_order, std::uint64_t{1});
    });

    suite.add("detects a contig revisited after being left", [] {
        // Ascending starts throughout, so a naive monotonicity check sees
        // nothing wrong -- but chr0's earlier state is gone by the time its
        // second block arrives, so duplicates across the interruption are lost.
        DuplicateMarker m;
        Fragment a = frag(0, 100, 200);
        Fragment b = frag(1, 100, 200);
        Fragment c = frag(0, 300, 400);
        m.observe(a);
        m.observe(b);
        m.observe(c);
        CHECK_EQ(m.stats().out_of_order, std::uint64_t{1});
    });

    suite.add("hash mode handles arbitrary input order", [] {
        DuplicateMarker m(DuplicateMarker::Mode::kHashAll);
        Fragment a = frag(0, 500, 600);
        Fragment b = frag(0, 100, 200);
        Fragment c = frag(0, 500, 600);  // duplicate of `a`, far away in the stream
        CHECK(!m.observe(a));
        CHECK(!m.observe(b));
        CHECK(m.observe(c));
        CHECK_EQ(m.stats().out_of_order, std::uint64_t{0});
    });

    suite.add("both modes agree on coordinate-sorted input", [] {
        std::mt19937 rng(4242);
        std::uniform_int_distribution<int> tid(0, 2);
        std::uniform_int_distribution<Coord> start(0, 2000);
        std::uniform_int_distribution<Coord> len(50, 400);
        std::bernoulli_distribution rev(0.5);

        std::vector<Fragment> frags;
        for (int i = 0; i < 20000; ++i) {
            const Coord s = start(rng);
            frags.push_back(frag(tid(rng), s, s + len(rng), rev(rng)));
        }

        // Hash mode over the shuffled order.
        DuplicateMarker hashed(DuplicateMarker::Mode::kHashAll);
        for (Fragment f : frags) hashed.observe(f);

        // Sorted mode over the same set, coordinate-sorted.
        sort_by_coord(frags);
        DuplicateMarker sorted(DuplicateMarker::Mode::kSortedStream);
        for (Fragment& f : frags) sorted.observe(f);

        CHECK_EQ(sorted.stats().examined, hashed.stats().examined);
        CHECK_EQ(sorted.stats().duplicates, hashed.stats().duplicates);
        CHECK_EQ(sorted.stats().unique, hashed.stats().unique);
        CHECK_EQ(sorted.stats().out_of_order, std::uint64_t{0});
        // And the memory claim: pile-up depth, not file size.
        CHECK(sorted.stats().peak_live_keys < hashed.stats().peak_live_keys / 10);
    });

    suite.add("library size estimate is sane", [] {
        DuplicateMarker m;
        // 1000 distinct positions, each observed twice: heavy duplication, so
        // the library must be estimated as only modestly larger than what was
        // seen.
        for (Coord i = 0; i < 1000; ++i) {
            Fragment a = frag(0, i * 10, i * 10 + 100);
            Fragment b = frag(0, i * 10, i * 10 + 100);
            m.observe(a);
            m.observe(b);
        }
        CHECK_EQ(m.stats().duplicates, std::uint64_t{1000});
        const double L = m.estimate_library_size();
        CHECK(L >= static_cast<double>(m.stats().unique));
        CHECK(L < 10.0 * static_cast<double>(m.stats().unique));

        // No duplicates at all leaves the estimate unbounded above; the
        // function reports 0 rather than an invented number.
        DuplicateMarker clean;
        for (Coord i = 0; i < 100; ++i) {
            Fragment f = frag(0, i * 10, i * 10 + 100);
            clean.observe(f);
        }
        CHECK_EQ(clean.estimate_library_size(), 0.0);
    });

    suite.add("observe_batch marks in place", [] {
        DuplicateMarker m;
        std::vector<Fragment> batch = {frag(0, 10, 110), frag(0, 10, 110), frag(0, 20, 120)};
        CHECK_EQ(m.observe_batch(batch), std::uint64_t{1});
        CHECK(!batch[0].duplicate());
        CHECK(batch[1].duplicate());
        CHECK(!batch[2].duplicate());
    });

    return suite.run();
}
