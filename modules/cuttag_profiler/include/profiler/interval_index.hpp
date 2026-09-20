// SPDX-License-Identifier: MIT
//
// Genomic interval indexing: "which annotations overlap this range?"
//
// The profiler's existing matrix path anchors a fixed-width window on a
// reference point (TSS/centre/TES). That cannot express a variable-length
// feature -- a gene body, an enhancer-promoter link, a multi-region
// intersection -- because the answer is a set whose size depends on the query.
// This header is the structure for those questions.
//
// LAYOUT. Per contig, three parallel arrays sorted by start, plus a payload
// column. A pointer-based red-black interval tree is the textbook answer and
// the wrong one here: one heap node per interval is one cache miss per level of
// the descent. Instead the tree is *implicit* in the sorted array -- the node
// for the index range [lo, hi) is at (lo + hi) / 2, its children are the ranges
// either side -- so navigating it is index arithmetic over contiguous memory
// and the whole contig is four allocations regardless of interval count.
//
//   starts_   [ s0 s1 s2 ... ]  ascending
//   ends_     [ e0 e1 e2 ... ]
//   max_end_  [ ... ]           max end over the implicit subtree rooted here
//   payload_  [ ... ]           index into the caller's Region vector
//
// max_end_ is what makes the descent O(log n + m) rather than a prefix scan: a
// subtree whose maximum end lies at or before the query start cannot contain an
// overlap and is pruned whole. It is also the only thing that handles fully
// nested intervals correctly -- [100, 200) containing [120, 130) -- which is
// precisely the case a "binary search the start, then scan forward" gets wrong.
//
// ZERO COPY. A query appends *payload ids* to caller-owned scratch and returns
// a span over it. Regions own std::strings; copying them per hit would put an
// allocation on the hot path. Callers resolve regions[id] only for the hits
// they actually use.
//
// THREADING. Immutable once built, therefore shareable as a const reference
// across the profiler's workers with no synchronisation -- the same contract
// BamReader's index and header already follow. The scratch vector is the only
// mutable state and it belongs to the calling thread.
#pragma once

#include <algorithm>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

#include "profiler/types.hpp"

namespace profiler {

// Per-contig coordinates are int32, not the int64 used on the Region API.
//
// Halving the hot arrays is the entire point of this structure: at 1M
// annotations, int32 keeps starts_ + ends_ + max_end_ at 12 MB instead of 24,
// which is the difference between fitting a chromosome's index in L3 and not.
// Human chr1 is 2.5e8, three orders of magnitude below the limit.
//
// Contigs at or beyond INT32_MAX -- some conifer and amphibian assemblies get
// there -- are rejected by name at build time. Truncating silently would place
// annotations at wrong coordinates, which is worse than not running.
using Coord = std::int32_t;

inline constexpr std::int64_t kMaxIndexableCoord =
    static_cast<std::int64_t>(std::numeric_limits<Coord>::max());

// One contig's intervals. Public members because this is a data structure, not
// an abstraction: signal_calc.cpp reads these arrays directly.
class ContigIntervals {
public:
    // Appends the payload id of every interval overlapping [qs, qe) to `out`,
    // in ascending start order. Returns the number appended. Allocates only if
    // `out` has to grow.
    std::size_t overlap(Coord qs, Coord qe, std::vector<std::uint32_t>& out) const {
        if (qs >= qe || starts_.empty()) return 0;
        const std::size_t before = out.size();
        descend(0, starts_.size(), qs, qe, out);
        return out.size() - before;
    }

    [[nodiscard]] std::size_t size() const noexcept { return starts_.size(); }

    [[nodiscard]] std::span<const Coord> starts() const noexcept { return starts_; }
    [[nodiscard]] std::span<const Coord> ends() const noexcept { return ends_; }
    [[nodiscard]] std::span<const std::uint32_t> payloads() const noexcept {
        return payload_;
    }

private:
    friend class IntervalIndex;

    // Half-open overlap. Spelled out once, here, and used nowhere else in this
    // file: BED is half-open and the off-by-one in this predicate is the single
    // most common defect in interval code.
    //
    // The `as < ae` term is not redundant. The familiar two-term form
    // `as < be && bs < ae` is only equivalent to "the intersection is
    // non-empty" when BOTH intervals are non-empty; for a zero-length interval
    // [x, x) strictly inside the query it evaluates true, reporting an overlap
    // with a feature that covers no bases. Zero-length records are legal in BED
    // (insertion points), so this is reachable from real input, and the
    // corresponding degenerate *query* is short-circuited in overlap().
    //
    // The rule, applied to both sides: an interval with start >= end is empty
    // and overlaps nothing. test_interval_index.cpp pins both halves; its
    // brute-force reference states the same rule the other way round, as
    // max(as, bs) < min(ae, be), so the differential test compares two
    // independent formulations rather than one expression against itself.
    static bool overlaps(Coord as, Coord ae, Coord bs, Coord be) noexcept {
        return as < ae && as < be && bs < ae;
    }

    // Recursive descent over the implicit tree on [lo, hi).
    //
    // Depth is ceil(log2(n)) -- 25 for 30M intervals, 31 at the int32 ceiling
    // -- so the C++ stack is in no danger and an explicit stack would only add
    // bookkeeping. Two prunes carry the complexity bound:
    //
    //   * max_end_[mid] <= qs: nothing in this subtree reaches the query.
    //   * starts_[mid] >= qe: the node and everything right of it starts at or
    //     after the query ends (the array is sorted by start), so only the left
    //     child can still contribute.
    void descend(std::size_t lo, std::size_t hi, Coord qs, Coord qe,
                 std::vector<std::uint32_t>& out) const {
        if (lo >= hi) return;
        const std::size_t mid = lo + (hi - lo) / 2;
        if (max_end_[mid] <= qs) return;

        descend(lo, mid, qs, qe, out);

        if (starts_[mid] >= qe) return;
        if (overlaps(starts_[mid], ends_[mid], qs, qe)) out.push_back(payload_[mid]);
        descend(mid + 1, hi, qs, qe, out);
    }

    // Post-order fill of max_end_ over the same implicit tree. Returns the
    // subtree maximum so the parent can fold it in without a second pass.
    Coord build_max(std::size_t lo, std::size_t hi) {
        if (lo >= hi) return std::numeric_limits<Coord>::min();
        const std::size_t mid = lo + (hi - lo) / 2;
        Coord m = ends_[mid];
        m = std::max(m, build_max(lo, mid));
        m = std::max(m, build_max(mid + 1, hi));
        max_end_[mid] = m;
        return m;
    }

    alignas(kCacheLineSize) std::vector<Coord> starts_;
    alignas(kCacheLineSize) std::vector<Coord> ends_;
    alignas(kCacheLineSize) std::vector<Coord> max_end_;
    alignas(kCacheLineSize) std::vector<std::uint32_t> payload_;
};

class IntervalIndex {
public:
    IntervalIndex() = default;

    // Builds an index over `regions`. Payload ids are indices into that vector,
    // which the caller must keep alive for as long as it resolves hits -- the
    // index stores no copy of it.
    //
    // Throws std::invalid_argument on a negative coordinate, an inverted
    // interval, or a coordinate past the int32 ceiling. Zero-length intervals
    // are kept: they are legal in BED as insertion points and simply overlap
    // nothing.
    [[nodiscard]] static IntervalIndex build(const std::vector<Region>& regions) {
        IntervalIndex index;
        if (regions.size() > std::numeric_limits<std::uint32_t>::max()) {
            throw std::invalid_argument(
                "interval index supports at most 2^32 regions, got " +
                std::to_string(regions.size()));
        }

        // Group by contig first so each contig is sorted independently. One
        // global sort by (chrom, start) would work too, but grouping keeps the
        // comparison on integers instead of strings.
        std::vector<std::vector<std::uint32_t>> groups;
        for (std::size_t i = 0; i < regions.size(); ++i) {
            const Region& r = regions[i];
            validate(r, i);
            auto it = index.local_of_name_.find(r.chrom);
            if (it == index.local_of_name_.end()) {
                it = index.local_of_name_.emplace(r.chrom, index.contig_names_.size()).first;
                index.contig_names_.push_back(r.chrom);
                groups.emplace_back();
            }
            groups[it->second].push_back(static_cast<std::uint32_t>(i));
        }

        index.by_contig_.resize(groups.size());
        for (std::size_t c = 0; c < groups.size(); ++c) {
            auto& ids = groups[c];
            // Ties on start broken by end, then by original index, so the
            // layout -- and therefore the order hits come back in -- is
            // reproducible for a given input rather than dependent on the
            // sort's internal choices.
            std::sort(ids.begin(), ids.end(), [&](std::uint32_t a, std::uint32_t b) {
                if (regions[a].start != regions[b].start) {
                    return regions[a].start < regions[b].start;
                }
                if (regions[a].end != regions[b].end) return regions[a].end < regions[b].end;
                return a < b;
            });

            ContigIntervals& ci = index.by_contig_[c];
            const std::size_t n = ids.size();
            ci.starts_.resize(n);
            ci.ends_.resize(n);
            ci.max_end_.resize(n);
            ci.payload_.resize(n);
            for (std::size_t j = 0; j < n; ++j) {
                const Region& r = regions[ids[j]];
                ci.starts_[j] = static_cast<Coord>(r.start);
                ci.ends_[j] = static_cast<Coord>(r.end);
                ci.payload_[j] = ids[j];
            }
            ci.build_max(0, n);
        }
        index.total_ = regions.size();
        return index;
    }

    // Resolves this index's contigs against a BAM header's ordering, so the hot
    // path takes a tid and never hashes a chromosome name. Contigs absent from
    // the header stay unreachable by tid (and are reported by
    // unbound_contigs()); a header contig with no annotations maps to nothing
    // and answers every query with an empty span.
    void bind_contigs(const std::vector<std::string>& header_order) {
        local_of_tid_.assign(header_order.size(), kNoContig);
        bound_ = 0;
        for (std::size_t tid = 0; tid < header_order.size(); ++tid) {
            auto it = local_of_name_.find(header_order[tid]);
            if (it != local_of_name_.end()) {
                local_of_tid_[tid] = it->second;
                ++bound_;
            }
        }
    }

    // Annotation contigs that bind_contigs() could not match to the header.
    // Usually a naming mismatch ("chr1" vs "1"), which otherwise shows up as a
    // silently empty result -- so callers should report this, not ignore it.
    [[nodiscard]] std::vector<std::string> unbound_contigs() const {
        std::vector<std::string> out;
        if (local_of_tid_.empty()) return contig_names_;
        std::vector<bool> seen(contig_names_.size(), false);
        for (const std::size_t local : local_of_tid_) {
            if (local != kNoContig) seen[local] = true;
        }
        for (std::size_t c = 0; c < contig_names_.size(); ++c) {
            if (!seen[c]) out.push_back(contig_names_[c]);
        }
        return out;
    }

    // Overlap query by BAM tid. Requires bind_contigs() first; an unknown or
    // unbound tid yields an empty span rather than an error, because a BAM
    // legitimately contains contigs the annotation does not cover.
    [[nodiscard]] std::span<const std::uint32_t> overlap(
        int tid, Coord qs, Coord qe, std::vector<std::uint32_t>& scratch) const {
        scratch.clear();
        if (tid < 0 || static_cast<std::size_t>(tid) >= local_of_tid_.size()) return {};
        const std::size_t local = local_of_tid_[static_cast<std::size_t>(tid)];
        if (local == kNoContig) return {};
        by_contig_[local].overlap(qs, qe, scratch);
        return scratch;
    }

    // Overlap query by contig name, for callers with no BAM header in hand
    // (the CLI's --annotation path over two BED files, say).
    [[nodiscard]] std::span<const std::uint32_t> overlap(
        std::string_view chrom, Coord qs, Coord qe,
        std::vector<std::uint32_t>& scratch) const {
        scratch.clear();
        const auto it = local_of_name_.find(std::string(chrom));
        if (it == local_of_name_.end()) return {};
        by_contig_[it->second].overlap(qs, qe, scratch);
        return scratch;
    }

    [[nodiscard]] std::size_t size() const noexcept { return total_; }
    [[nodiscard]] std::size_t contig_count() const noexcept { return by_contig_.size(); }
    [[nodiscard]] const std::vector<std::string>& contigs() const noexcept {
        return contig_names_;
    }
    // Contigs successfully matched by the last bind_contigs() call.
    [[nodiscard]] std::size_t bound_contig_count() const noexcept { return bound_; }

private:
    static constexpr std::size_t kNoContig = static_cast<std::size_t>(-1);

    static void validate(const Region& r, std::size_t i) {
        const std::string where =
            "region " + std::to_string(i) + " (" + r.chrom + ":" +
            std::to_string(r.start) + "-" + std::to_string(r.end) + ")";
        if (r.start < 0 || r.end < 0) {
            throw std::invalid_argument(where + " has a negative coordinate");
        }
        if (r.end < r.start) {
            throw std::invalid_argument(where + " ends before it starts");
        }
        if (r.end > kMaxIndexableCoord) {
            throw std::invalid_argument(
                where + " exceeds the " + std::to_string(kMaxIndexableCoord) +
                " bp coordinate ceiling of the interval index. Contig '" + r.chrom +
                "' is too long for the int32 layout this structure uses.");
        }
    }

    std::vector<ContigIntervals> by_contig_;
    std::vector<std::string> contig_names_;
    std::unordered_map<std::string, std::size_t> local_of_name_;
    std::vector<std::size_t> local_of_tid_;
    std::size_t bound_ = 0;
    std::size_t total_ = 0;
};

}  // namespace profiler
