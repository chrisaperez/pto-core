// SPDX-License-Identifier: MIT
//
// Streaming duplicate identification.
//
// THE PROBLEM WITH THE INCUMBENT. Picard MarkDuplicates is a sort in disguise:
// it buffers read ends in a JVM heap, spills to a temp directory when the heap
// runs out, and its runtime on a large BAM is dominated by that spill, not by
// the duplicate logic. The usual pipeline pays for it twice, because
// coordinate-sorting for the aligner and coordinate-sorting for MarkDuplicates
// are the same order and it still round-trips through disk between them.
//
// THE OBSERVATION. Two fragments are PCR/optical duplicates when their 5' ends
// agree -- same contig, same start, same end, same orientation. On a
// coordinate-sorted input, every fragment sharing a start arrives
// *contiguously*. So duplicate detection needs no global structure at all: it
// needs a set keyed on (end, orientation) that lives exactly as long as one
// start position. Memory is O(deepest pile-up), typically a few hundred
// entries, never a function of file size, and there is no spill because there
// is nothing to spill.
//
// That is kSortedStream below, and it is the mode to use.
//
// kHashAll is the fallback for input that is not coordinate-sorted (a
// name-sorted BAM, a concatenation of fragment files). It keeps a 64-bit
// fingerprint of every distinct fragment in an open-addressed table, so memory
// is O(distinct fragments) -- still far below Picard's per-read Java objects,
// but no longer constant. It exists so that "your input is not sorted" is a
// slower answer rather than no answer.
//
// TRADE-OFFS, STATED PLAINLY. This is not a drop-in MarkDuplicates:
//
//   * Duplicate identity uses the *aligned* fragment span. Picard uses the
//     unclipped 5' position, so a duplicate pair whose copies were soft-clipped
//     differently is one duplicate to Picard and two distinct fragments here.
//     For assays where clipping is rare (ATAC, CUT&Tag, WGS) the rates agree
//     closely; for amplicon panels with primer clipping they will not.
//   * No optical/sequencer duplicate distinction: read names are not parsed for
//     tile coordinates, so `estimate_library_size()` treats all duplicates as
//     PCR duplicates and is therefore a lower bound on library complexity.
//   * No UMI awareness. A UMI-aware run should key on the UMI instead, which
//     this structure supports in principle but the CLI does not expose yet.
//   * The representative kept is the first fragment seen at a position, not
//     the highest-quality one. For counting and for FRiP this is immaterial;
//     for variant calling on the marked output it is not, which is why this
//     tool reports duplicates and does not rewrite a BAM.
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <span>
#include <vector>

#include "toolkit/types.hpp"

namespace toolkit {

struct DuplicateStats {
    std::uint64_t examined = 0;
    std::uint64_t duplicates = 0;
    std::uint64_t unique = 0;
    // Largest number of distinct fragments alive at one start position
    // (kSortedStream) or table occupancy (kHashAll). The headline number for
    // "how much memory did this actually need".
    std::uint64_t peak_live_keys = 0;
    // Fragments that arrived before their predecessor in kSortedStream mode.
    // Non-zero means the input was not coordinate-sorted and the duplicate
    // rate is an undercount; the CLI turns this into an error.
    std::uint64_t out_of_order = 0;

    [[nodiscard]] double duplicate_rate() const noexcept {
        return examined == 0 ? 0.0
                             : static_cast<double>(duplicates) / static_cast<double>(examined);
    }
};

class DuplicateMarker {
public:
    enum class Mode {
        kSortedStream,  // O(pile-up depth) memory; requires coordinate order
        kHashAll,       // O(distinct fragments) memory; any input order
    };

    explicit DuplicateMarker(Mode mode = Mode::kSortedStream) : mode_(mode) {
        if (mode_ == Mode::kHashAll) rehash(kInitialSlots);
    }

    // Returns true if `f` duplicates a fragment already seen. Sets
    // fragflags::kDuplicate on `f` as a side effect, so a caller streaming a
    // batch can filter in place.
    bool observe(Fragment& f) {
        const bool dup = (mode_ == Mode::kSortedStream) ? observe_sorted(f) : observe_hashed(f);
        ++stats_.examined;
        if (dup) {
            ++stats_.duplicates;
            f.flags |= fragflags::kDuplicate;
        } else {
            ++stats_.unique;
        }
        return dup;
    }

    // Convenience for the batch callback: marks every fragment in `batch` and
    // returns how many were duplicates.
    std::uint64_t observe_batch(std::span<Fragment> batch) {
        std::uint64_t n = 0;
        for (Fragment& f : batch) n += observe(f) ? 1 : 0;
        return n;
    }

    [[nodiscard]] const DuplicateStats& stats() const noexcept { return stats_; }
    [[nodiscard]] Mode mode() const noexcept { return mode_; }

    // Lander-Waterman library size estimate: the number of distinct molecules
    // in the library that would produce `unique` distinct observations out of
    // `examined` draws. Solved by bisection on
    //
    //     unique / examined == (1 - exp(-examined / L)) / (examined / L)
    //
    // Returns 0 when no duplicates were seen (the estimate is unbounded above)
    // or when every read was a duplicate (degenerate).
    [[nodiscard]] double estimate_library_size() const noexcept {
        const double n = static_cast<double>(stats_.examined);
        const double u = static_cast<double>(stats_.unique);
        if (n <= 0.0 || u <= 0.0 || u >= n) return 0.0;
        // f(L) = u/n - (1 - e^{-n/L}) / (n/L) is strictly DECREASING in L: a
        // larger library means a larger expected distinct fraction, so the
        // residual falls. Bisection must bracket accordingly -- with the
        // inequalities the other way round the search runs away from the root
        // and returns an absurd library size.
        auto f = [&](double L) {
            const double x = n / L;
            return u / n - (1.0 - std::exp(-x)) / x;
        };
        double lo = u;  // at least as many molecules as distinct observations
        double hi = u * 2.0;
        for (int i = 0; i < 64 && f(hi) > 0.0; ++i) hi *= 2.0;
        for (int i = 0; i < 200; ++i) {
            const double mid = 0.5 * (lo + hi);
            if (f(mid) > 0.0) {
                lo = mid;  // still short of the root; move right
            } else {
                hi = mid;
            }
        }
        return 0.5 * (lo + hi);
    }

private:
    // --- kSortedStream --------------------------------------------------
    //
    // `live_` holds the (end, orientation) keys seen at the current
    // (tid, start). It is a small sorted-insert vector, not a hash set: pile-up
    // depth is typically < 100, where a linear scan over a contiguous array
    // beats hashing and allocates nothing after the first position.
    bool observe_sorted(const Fragment& f) {
        if (f.tid != cur_tid_ || f.start != cur_start_) {
            if (f.tid == cur_tid_) {
                if (f.start < cur_start_) ++stats_.out_of_order;
            } else {
                // A contig we already finished has come back. Its earlier
                // fragments are gone from `live_`, so duplicates spanning the
                // interruption are missed -- exactly the undercount
                // out_of_order exists to flag. Checked as well as the
                // within-contig regression, because a BAM interleaved by
                // contig has ascending starts throughout and would otherwise
                // look perfectly sorted.
                if (std::find(seen_tids_.begin(), seen_tids_.end(), f.tid) !=
                    seen_tids_.end()) {
                    ++stats_.out_of_order;
                } else {
                    seen_tids_.push_back(f.tid);
                }
            }
            cur_tid_ = f.tid;
            cur_start_ = f.start;
            live_.clear();
        }
        const std::uint64_t key = end_key(f);
        for (const std::uint64_t k : live_) {
            if (k == key) return true;
        }
        live_.push_back(key);
        stats_.peak_live_keys = std::max<std::uint64_t>(stats_.peak_live_keys, live_.size());
        return false;
    }

    // --- kHashAll -------------------------------------------------------
    //
    // Open addressing with linear probing over a power-of-two table of 64-bit
    // fingerprints. Linear probing rather than chaining because the whole point
    // is cache behaviour: a probe sequence walks one cache line, a chain walks
    // one pointer dereference per step.
    //
    // A fingerprint collision reports a false duplicate. At 64 bits and 1e9
    // fragments the expected number of collisions is ~2.7e-2, i.e. this
    // essentially never fires; it is a real (if negligible) source of
    // over-counting and is recorded here rather than left implicit.
    bool observe_hashed(const Fragment& f) {
        if (used_ * 10 >= slots_.size() * 7) rehash(slots_.size() * 2);
        const std::uint64_t h = mix(full_key(f));
        const std::uint64_t mask = slots_.size() - 1;
        std::uint64_t i = h & mask;
        while (slots_[i] != kEmpty) {
            if (slots_[i] == h) return true;
            i = (i + 1) & mask;
        }
        slots_[i] = h;
        ++used_;
        stats_.peak_live_keys = used_;
        return false;
    }

    void rehash(std::size_t new_slots) {
        if (new_slots < kInitialSlots) new_slots = kInitialSlots;
        std::vector<std::uint64_t> next(new_slots, kEmpty);
        const std::uint64_t mask = new_slots - 1;
        for (const std::uint64_t h : slots_) {
            if (h == kEmpty) continue;
            std::uint64_t i = h & mask;
            while (next[i] != kEmpty) i = (i + 1) & mask;
            next[i] = h;
        }
        slots_.swap(next);
    }

    // (end, orientation) -- identity within one already-fixed start position.
    static std::uint64_t end_key(const Fragment& f) noexcept {
        return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(f.end)) << 1) |
               (f.reverse() ? 1ULL : 0ULL);
    }

    // (tid, start, end, orientation) -- full identity, for unsorted input.
    static std::uint64_t full_key(const Fragment& f) noexcept {
        std::uint64_t k = static_cast<std::uint32_t>(f.tid);
        k = k * 0x100000001b3ULL + static_cast<std::uint32_t>(f.start);
        k = k * 0x100000001b3ULL + static_cast<std::uint32_t>(f.end);
        k = k * 0x100000001b3ULL + (f.reverse() ? 1ULL : 0ULL);
        return k;
    }

    // splitmix64 finaliser. The FNV-ish accumulation above leaves structure in
    // the low bits, and linear probing is unusually sensitive to that.
    static std::uint64_t mix(std::uint64_t x) noexcept {
        x += 0x9e3779b97f4a7c15ULL;
        x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
        x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
        x ^= x >> 31;
        // kEmpty is a legal hash value; nudge it so the sentinel stays unique.
        return x == kEmpty ? 1ULL : x;
    }

    static constexpr std::uint64_t kEmpty = 0;
    static constexpr std::size_t kInitialSlots = 1024;

    Mode mode_;
    DuplicateStats stats_;

    std::int32_t cur_tid_ = kNoTid;
    Coord cur_start_ = -1;
    std::vector<std::uint64_t> live_;
    std::vector<std::int32_t> seen_tids_;

    std::vector<std::uint64_t> slots_;
    std::size_t used_ = 0;
};

}  // namespace toolkit
