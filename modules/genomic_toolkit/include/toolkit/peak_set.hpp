// SPDX-License-Identifier: MIT
//
// Merged, non-overlapping peak intervals with O(log n) overlap and
// base-coverage queries.
//
// WHY NOT AN INTERVAL TREE. cuttag_profiler has one (an implicit tree over
// sorted arrays with a max-end column) and it is the right structure there,
// because it answers "give me every annotation overlapping this window" and
// the answer is a set. FRiP asks a strictly weaker question -- "does this
// fragment touch any peak", and "how many of its bases do" -- for which the
// set of hits is not needed and, worse, is actively wrong to enumerate: peak
// callers routinely emit overlapping or bookended intervals, and counting a
// fragment once per overlapping peak inflates FRiP without bound. So the peaks
// are *merged* at build time into a disjoint, ascending cover, and the query
// collapses to one binary search:
//
//     the only candidate is the last interval whose start <= query start,
//     plus the run of intervals starting inside the query.
//
// Merging is not an optimisation, it is the definition of the statistic. This
// is also what `bedtools intersect -u` does internally, and the reason a naive
// `-wa` pipeline over unmerged peaks reports a different FRiP than the paper
// it is being compared to.
//
// LAYOUT. Per contig, two parallel int32 arrays plus a prefix-sum column of
// covered bases. Structure-of-arrays because the binary search only ever
// touches `starts`: at 1M peaks the descent walks ~20 cache lines of starts
// instead of ~20 lines of interleaved (start, end, sum) triples.
//
// THREADING. Immutable once built, so it is shared by const reference across
// FRiP workers with no synchronisation.
#pragma once

#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include "toolkit/contig_dict.hpp"
#include "toolkit/types.hpp"

namespace toolkit {

// Per-call-site state for a coordinate-sorted sweep of overlap_bases()
// queries. FRiP's fragment stream is sorted, so `qs` is almost always
// non-decreasing across calls -- the property this exploits to answer in O(1)
// amortised instead of a cold O(log n) binary search every fragment: the
// candidate interval only ever moves forward, so it gallops from where the
// last query left it instead of re-deriving it from scratch.
//
// Correctness never depends on the caller actually being sorted: a query that
// moves backwards, or against a different contig than this cursor's last use
// (detected by PeakSet::overlap_bases via `last_contig`), resets to the exact
// same binary search overlap_bases() always did. The cursor can only make an
// answer cheaper, never different -- see floor_index(Coord, OverlapCursor&).
//
// One cursor per logical sweep, not shared across concurrent sweeps: it is
// mutated on every call and carries no synchronisation of its own, matching
// FripCounter's per-thread ownership of one cursor each.
struct OverlapCursor {
    std::size_t i = 0;
    Coord last_q = -1;
    const void* last_contig = nullptr;
};

// One contig's disjoint intervals, ascending and non-adjacent.
struct ContigPeaks {
    alignas(kCacheLineSize) std::vector<Coord> starts;
    alignas(kCacheLineSize) std::vector<Coord> ends;
    // covered_before[i] = total peak bases in intervals [0, i). Lets a
    // multi-interval base-overlap query skip the interior run entirely.
    alignas(kCacheLineSize) std::vector<std::int64_t> covered_before;
    std::int64_t covered_total = 0;

    [[nodiscard]] std::size_t size() const noexcept { return starts.size(); }

    // Index of the last interval with start <= q, or npos if none.
    [[nodiscard]] std::size_t floor_index(Coord q) const noexcept {
        // upper_bound gives the first start > q; step back one.
        const auto it = std::upper_bound(starts.begin(), starts.end(), q);
        if (it == starts.begin()) return kNpos;
        return static_cast<std::size_t>(it - starts.begin()) - 1;
    }

    static constexpr std::size_t kNpos = static_cast<std::size_t>(-1);

    // Cursor-aware equivalent of floor_index(Coord). `cursor` must already
    // belong to THIS ContigPeaks (PeakSet::overlap_bases enforces that via
    // `last_contig`; calling this directly on a cursor from a different
    // ContigPeaks is a wrong answer, not a memory-safety issue -- `i` would
    // simply index this array using a position derived from another one).
    [[nodiscard]] std::size_t floor_index(Coord q, OverlapCursor& cursor) const noexcept {
        if (q < cursor.last_q) {
            cursor.i = 0;
            cursor.last_q = -1;  // not monotone: fall through to the cold path below
        }
        if (cursor.last_q < 0) {
            const std::size_t result = floor_index(q);
            cursor.last_q = q;
            cursor.i = (result == kNpos) ? 0 : result;
            return result;
        }
        cursor.last_q = q;
        while (cursor.i + 1 < starts.size() && starts[cursor.i + 1] <= q) ++cursor.i;
        return starts[cursor.i] > q ? kNpos : cursor.i;
    }

    // True if [qs, qe) touches any interval. Empty queries touch nothing.
    [[nodiscard]] bool overlaps(Coord qs, Coord qe) const noexcept {
        if (qs >= qe || starts.empty()) return false;
        // Case 1: the interval covering (or preceding) qs reaches past qs.
        const std::size_t i = floor_index(qs);
        if (i != kNpos && ends[i] > qs) return true;
        // Case 2: some interval starts strictly inside [qs, qe). Because the
        // cover is disjoint and ascending, that is exactly the next index.
        const std::size_t next = (i == kNpos) ? 0 : i + 1;
        return next < starts.size() && starts[next] < qe;
    }

    // Number of bases of [qs, qe) that fall inside a peak.
    [[nodiscard]] std::int64_t overlap_bases(Coord qs, Coord qe) const noexcept {
        if (qs >= qe || starts.empty()) return 0;
        std::int64_t total = 0;
        std::size_t i = floor_index(qs);
        if (i == kNpos) {
            i = 0;
        } else {
            total += std::max<std::int64_t>(
                0, std::min(static_cast<std::int64_t>(ends[i]), static_cast<std::int64_t>(qe)) -
                       static_cast<std::int64_t>(qs));
            ++i;
        }
        // Whole intervals contained in the query. `ends` is ascending because
        // the cover is disjoint, so the far edge is a binary search and the
        // interior run is then a prefix-sum subtraction -- a fragment spanning
        // a peak-dense region costs two searches, never a scan.
        const std::size_t j = static_cast<std::size_t>(
            std::upper_bound(ends.begin() + static_cast<std::ptrdiff_t>(i), ends.end(), qe) -
            ends.begin());
        if (j > i) total += covered_before[j] - covered_before[i];
        // Partial interval at the far edge.
        if (j < starts.size() && starts[j] < qe) {
            total += static_cast<std::int64_t>(qe) - static_cast<std::int64_t>(starts[j]);
        }
        return total;
    }

    // Cursor-aware equivalent of overlap_bases(Coord, Coord). Identical body
    // to the plain overload above but for the floor_index call -- kept as a
    // full copy rather than routed through a shared template so the plain,
    // hot, cursor-free path (used by overlaps_named() and anything without a
    // sweep) carries no cursor-branch overhead at all.
    [[nodiscard]] std::int64_t overlap_bases(Coord qs, Coord qe,
                                             OverlapCursor& cursor) const noexcept {
        if (qs >= qe || starts.empty()) return 0;
        std::int64_t total = 0;
        std::size_t i = floor_index(qs, cursor);
        if (i == kNpos) {
            i = 0;
        } else {
            total += std::max<std::int64_t>(
                0, std::min(static_cast<std::int64_t>(ends[i]), static_cast<std::int64_t>(qe)) -
                       static_cast<std::int64_t>(qs));
            ++i;
        }
        // A forward scan rather than upper_bound: since `ends` is ascending,
        // "advance while ends[j] <= qe" reaches exactly the same first index
        // past qe that upper_bound would, and for a sorted sweep the cursor
        // has already put `i` next to the answer -- a real fragment overlaps
        // a small, near-constant number of peaks, so this is typically 0-2
        // comparisons against upper_bound's guaranteed log2(n).
        std::size_t j = i;
        while (j < ends.size() && ends[j] <= qe) ++j;
        if (j > i) total += covered_before[j] - covered_before[i];
        if (j < starts.size() && starts[j] < qe) {
            total += static_cast<std::int64_t>(qe) - static_cast<std::int64_t>(starts[j]);
        }
        return total;
    }
};

class PeakSet {
public:
    PeakSet() = default;

    // Builds from arbitrary regions: sorts, merges overlapping *and*
    // book-ended intervals, and interns contig names.
    //
    // Throws std::invalid_argument on a negative or inverted interval, or a
    // coordinate past the int32 ceiling.
    [[nodiscard]] static PeakSet build(const std::vector<Region>& regions) {
        PeakSet ps;
        std::vector<std::vector<std::pair<Coord, Coord>>> groups;

        for (std::size_t i = 0; i < regions.size(); ++i) {
            const Region& r = regions[i];
            validate(r, i);
            if (r.end == r.start) continue;  // zero-length: covers no bases
            const std::int32_t tid = ps.dict_.intern(r.chrom);
            if (static_cast<std::size_t>(tid) >= groups.size()) groups.resize(tid + 1);
            groups[static_cast<std::size_t>(tid)].emplace_back(
                static_cast<Coord>(r.start), static_cast<Coord>(r.end));
        }

        ps.by_tid_.resize(groups.size());
        for (std::size_t t = 0; t < groups.size(); ++t) {
            auto& g = groups[t];
            std::sort(g.begin(), g.end());
            ContigPeaks& cp = ps.by_tid_[t];
            for (const auto& [s, e] : g) {
                // `s <= cp.ends.back()` merges book-ended intervals too
                // ([100,200) and [200,300) become [100,300)). They cover a
                // contiguous stretch of bases, so leaving them split would let
                // a fragment spanning the join be counted against two peaks.
                if (!cp.starts.empty() && s <= cp.ends.back()) {
                    cp.ends.back() = std::max(cp.ends.back(), e);
                } else {
                    cp.starts.push_back(s);
                    cp.ends.push_back(e);
                }
            }
            cp.covered_before.resize(cp.starts.size() + 1, 0);
            for (std::size_t i = 0; i < cp.starts.size(); ++i) {
                cp.covered_before[i + 1] =
                    cp.covered_before[i] +
                    (static_cast<std::int64_t>(cp.ends[i]) - cp.starts[i]);
            }
            cp.covered_total = cp.covered_before.back();
            ps.merged_count_ += cp.starts.size();
            ps.covered_bases_ += cp.covered_total;
        }
        ps.input_count_ = regions.size();
        return ps;
    }

    // Resolves this peak set's contigs against the fragment stream's
    // dictionary, so queries take the stream's tid directly. Call once, before
    // the streaming pass.
    void bind(const ContigDict& stream_dict) {
        local_of_stream_tid_ = stream_dict.map_onto(dict_);
    }

    // Overlap test by *stream* tid. An unbound or unknown tid answers false --
    // a BAM legitimately contains contigs a peak file does not cover -- which
    // is why bind()'s naming mismatches must be reported separately by the
    // caller rather than inferred from a zero result.
    [[nodiscard]] bool overlaps(std::int32_t stream_tid, Coord qs, Coord qe) const noexcept {
        const ContigPeaks* cp = resolve(stream_tid);
        return cp != nullptr && cp->overlaps(qs, qe);
    }

    [[nodiscard]] std::int64_t overlap_bases(std::int32_t stream_tid, Coord qs,
                                             Coord qe) const noexcept {
        const ContigPeaks* cp = resolve(stream_tid);
        return cp == nullptr ? 0 : cp->overlap_bases(qs, qe);
    }

    // Cursor-aware equivalent, for a coordinate-sorted sweep of queries (a
    // fragment stream). `cursor` is reset automatically the first time it is
    // used and whenever the resolved contig changes from its last use, so one
    // cursor safely follows a stream across contig boundaries -- it must
    // still be one cursor per concurrent sweep, never shared between threads.
    [[nodiscard]] std::int64_t overlap_bases(std::int32_t stream_tid, Coord qs, Coord qe,
                                             OverlapCursor& cursor) const noexcept {
        const ContigPeaks* cp = resolve(stream_tid);
        if (cp == nullptr) return 0;
        if (cp != cursor.last_contig) cursor = OverlapCursor{.last_contig = cp};
        return cp->overlap_bases(qs, qe, cursor);
    }

    // Query by name, for callers with no stream dictionary in hand.
    [[nodiscard]] bool overlaps_named(std::string_view chrom, Coord qs, Coord qe) const {
        const std::int32_t tid = dict_.lookup(chrom);
        if (tid == kNoTid) return false;
        return by_tid_[static_cast<std::size_t>(tid)].overlaps(qs, qe);
    }

    [[nodiscard]] const ContigDict& dict() const noexcept { return dict_; }
    // Intervals after merging; compare against input_count() to see how much
    // the peak file overlapped itself.
    [[nodiscard]] std::size_t size() const noexcept { return merged_count_; }
    [[nodiscard]] std::size_t input_count() const noexcept { return input_count_; }
    [[nodiscard]] std::int64_t covered_bases() const noexcept { return covered_bases_; }

private:
    [[nodiscard]] const ContigPeaks* resolve(std::int32_t stream_tid) const noexcept {
        if (stream_tid < 0 ||
            static_cast<std::size_t>(stream_tid) >= local_of_stream_tid_.size()) {
            return nullptr;
        }
        const std::int32_t local = local_of_stream_tid_[static_cast<std::size_t>(stream_tid)];
        if (local == kNoTid) return nullptr;
        return &by_tid_[static_cast<std::size_t>(local)];
    }

    static void validate(const Region& r, std::size_t i) {
        const std::string where = "peak " + std::to_string(i) + " (" + r.chrom + ":" +
                                  std::to_string(r.start) + "-" + std::to_string(r.end) + ")";
        if (r.start < 0 || r.end < 0) {
            throw std::invalid_argument(where + " has a negative coordinate");
        }
        if (r.end < r.start) throw std::invalid_argument(where + " ends before it starts");
        if (r.end > kMaxCoord) {
            throw std::invalid_argument(
                where + " exceeds the " + std::to_string(kMaxCoord) +
                " bp coordinate ceiling of the int32 layout");
        }
    }

    std::vector<ContigPeaks> by_tid_;
    ContigDict dict_;
    std::vector<std::int32_t> local_of_stream_tid_;
    std::size_t merged_count_ = 0;
    std::size_t input_count_ = 0;
    std::int64_t covered_bases_ = 0;
};

}  // namespace toolkit
