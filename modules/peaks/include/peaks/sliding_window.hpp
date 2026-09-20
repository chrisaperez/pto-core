// SPDX-License-Identifier: MIT
//
// pto-peaks: the multi-scale sliding background window.
//
// What it is
// ----------
// A per-base pileup and three centred background windows (1kb, 5kb, 10kb) over
// a stream of coordinate-sorted fragments, in O(1) per base and in memory that
// does not depend on the length of the chromosome.
//
// This is the structure MACS2 spends most of its wall time working around. It
// writes a bedGraph of the pileup to disk, then re-scans that file once per
// background scale, which is four passes over a genome-sized intermediate. Here
// the pileup exists only as a running integer and the three background sums are
// running differences, so the whole chromosome costs one pass and 256 KiB.
//
// Two rings, and why there are two
// --------------------------------
// `delta_` holds the +1/-1 events of fragment starts and ends. Walking it
// forward and accumulating gives the pileup, and each slot is cleared as the
// frontier passes it, so the ring is self-cleaning.
//
// `pileup_` holds the resulting per-base values, because a CENTRED window needs
// to subtract a base that is now 5000 bp behind the evaluation position, and
// that value is no longer recoverable from a running total. Both rings are
// kRingCapacity slots of int64 and are allocated exactly once, in the
// constructor; reset() only zeroes them.
//
// Centred, not trailing
// ---------------------
// The background window for base p is [p - W/2, p + W/2). Trailing windows are
// cheaper -- no lookahead, one ring -- and they are wrong in a way that is hard
// to see in the output: the background a peak is judged against would include
// the peak's own rising edge, so every call is biased shorter at its 3' end.
//
// Centring costs a lag. Fragments arrive sorted by start, so the pileup at base
// q is final once a fragment with start > q has been seen; the evaluation
// position therefore trails the frontier by the widest half-window, 5000 bp.
// That is what `evaluable_limit()` reports and what `step()` refuses to go past.
//
// Coordinate arithmetic
// ---------------------
// Every position computation here is int64, including the ones whose operands
// are int32 coordinates. `pos + 5000` on an int32 near INT32_MAX is undefined
// behaviour, and it happens BEFORE any range check can look at the result --
// which is exactly the defect docs/REVIEW_2026-08-15.md finding (b) records in
// cuttag_profiler's `window_length()`. The accumulators additionally go through
// `__builtin_add_overflow`; they cannot realistically overflow an int64, but a
// saturating accumulator that says so beats one that wraps quietly, and
// `saturated()` is wired to a lambda of infinity, which is p = 1, which is the
// fail-closed direction.
#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

namespace pto::peaks {

// Per-contig coordinates are int32, as in every mainstream genome format's
// on-disk representation. Human chr1 is 2.5e8, an order of magnitude below the
// ceiling; coordinates past it are refused rather than truncated, because a
// truncated coordinate silently places signal on the wrong part of a
// chromosome.
using Coord = std::int32_t;

inline constexpr std::int32_t kNoTid = -1;

// MACS2's three local background scales, in this order.
inline constexpr int kScaleCount = 3;
inline constexpr Coord kScaleWidth[kScaleCount] = {1000, 5000, 10000};
inline constexpr Coord kMaxScaleWidth = 10000;
inline constexpr Coord kMaxHalfWidth = kMaxScaleWidth / 2;

// The longest fragment the delta ring can represent. Real paired-end inserts
// are under 1 kb and genomic_toolkit's default filter caps them at 1 kb; four
// times that leaves room for long-fragment protocols without letting a
// malformed record wrap the ring onto live data.
inline constexpr Coord kMaxFragmentSpan = 4096;

// Power of two, so the ring index is a mask rather than a division, and large
// enough that the widest window and the longest fragment never overlap in it.
inline constexpr std::int64_t kRingCapacity = 16384;

static_assert(kRingCapacity > kMaxScaleWidth + kMaxFragmentSpan,
              "the ring must hold the widest window and the longest fragment");
static_assert((kRingCapacity & (kRingCapacity - 1)) == 0,
              "kRingCapacity must be a power of two");

// Why a fragment was refused. Anything other than kAccepted leaves the window
// exactly as it was: a record that cannot be placed correctly is not placed at
// all, and the count is what lets a CLI exit 2 rather than print a track that
// looks fine.
//
// kNeedsDrain is not a rejection and is not counted as one. It is back
// pressure: the ring holds a bounded lookahead, and a caller that added a
// whole chromosome before evaluating any of it would lap the reader. The
// protocol is to step() until the window has caught up and then offer the same
// fragment again -- which is what the natural streaming loop does anyway:
//
//     while (next_fragment(f)) {
//       while (w.add_fragment(f.start, f.end) == AddResult::kNeedsDrain) {
//         if (!w.step()) break;
//         visit(w);
//       }
//       while (w.step()) visit(w);
//     }
//
// Progress is guaranteed: a drain always frees at least
// kRingCapacity - kMaxFragmentSpan - kMaxScaleWidth bases of lookahead, so the
// frontier can always advance, whatever the gap between two fragments.
enum class AddResult {
  kAccepted = 0,
  kOutOfOrder = 1,
  kOutOfBounds = 2,
  kBadSpan = 3,
  kNeedsDrain = 4,
};

class MultiScaleWindow {
 public:
  // The only allocation this class ever performs.
  MultiScaleWindow()
      : delta_(static_cast<std::size_t>(kRingCapacity), 0),
        pileup_(static_cast<std::size_t>(kRingCapacity), 0) {}

  // Start a contig. Allocates nothing; the rings are zeroed in place, which is
  // what keeps the previous contig's tail from wrapping onto this one's head.
  void reset(std::int32_t tid, Coord contig_length) noexcept {
    std::fill(delta_.begin(), delta_.end(), 0);
    std::fill(pileup_.begin(), pileup_.end(), 0);
    tid_ = tid;
    contig_len_ = contig_length < 0 ? 0 : static_cast<std::int64_t>(contig_length);
    frontier_ = -1;
    running_ = 0;
    // Start far enough behind that every window is empty, so the rolling sums
    // can be primed by the same update that maintains them. step() walks
    // through the negative positions without reporting them.
    pos_ = -static_cast<std::int64_t>(kMaxHalfWidth) - 1;
    last_start_ = -1;
    max_end_ = 0;
    finished_ = false;
    saturated_ = false;
    for (int i = 0; i < kScaleCount; ++i) sum_[i] = 0;
    for (auto& r : rejected_) r = 0;
  }

  [[nodiscard]] std::int32_t tid() const noexcept { return tid_; }

  // Fragments must arrive with non-decreasing `start`. Equal starts are in
  // order -- a pile of identical fragments is what a real peak looks like.
  AddResult add_fragment(Coord start, Coord end) noexcept {
    const std::int64_t s = start;
    const std::int64_t e = end;
    if (e <= s || e - s > kMaxFragmentSpan) return reject(AddResult::kBadSpan);
    if (s < 0 || e > contig_len_) return reject(AddResult::kOutOfBounds);
    if (finished_) return reject(AddResult::kOutOfOrder);
    if (s < last_start_) return reject(AddResult::kOutOfOrder);

    // The ring cannot hold an unbounded distance between the writer and the
    // reader. Advance as far as is safe, then ask the caller to drain and
    // offer this fragment again; the partial advance is what guarantees the
    // retry makes progress rather than spinning.
    if (s > safe_start_limit()) {
      advance_frontier_to(safe_start_limit() - 1);
      return AddResult::kNeedsDrain;
    }

    // Every delta at a position strictly below `s` is now known, because no
    // later fragment can start there. Deltas AT `s` are not: the next fragment
    // may share this start.
    advance_frontier_to(s - 1);
    delta_[ring(s)] += 1;
    delta_[ring(e)] -= 1;
    last_start_ = s;
    if (e > max_end_) max_end_ = e;
    return AddResult::kAccepted;
  }

  // Fast-forward across a region with no coverage.
  //
  // A chromosome is not uniformly sequenced: centromeres, satellite arrays and
  // unmappable stretches are megabases with no fragments at all, and stepping
  // through them one base at a time is work that cannot produce a peak. This
  // skips them.
  //
  // It is only legal when the jump provably lands in empty space -- every
  // fragment seen so far must end at least a widest-half-window behind the
  // target, so all three rolling sums are zero on the far side. Otherwise it
  // refuses and changes nothing, because a seek that silently discarded live
  // coverage would erase a peak rather than mis-score it.
  bool seek(Coord target) noexcept {
    const std::int64_t t = target;
    if (finished_ || t >= contig_len_ || t <= pos_) return false;
    if (t - max_end_ < kMaxHalfWidth) return false;

    advance_frontier_to(std::min<std::int64_t>(max_end_, t - 1));
    if (running_ != 0) return false;  // something is still open; refuse

    std::fill(delta_.begin(), delta_.end(), 0);
    std::fill(pileup_.begin(), pileup_.end(), 0);
    pos_ = t;
    frontier_ = t - 1;
    running_ = 0;
    for (int i = 0; i < kScaleCount; ++i) sum_[i] = 0;
    return true;
  }

  // No more fragments. Every remaining base becomes evaluable, but the
  // frontier is NOT run to the end of the contig here: it advances lazily in
  // step(), staying one half-window ahead of the evaluation position.
  //
  // Racing it to the end instead is the mistake this comment exists to record.
  // The frontier writes into the same pileup ring the rolling sums read from,
  // so on a contig longer than kRingCapacity it laps the reader and every
  // window sum comes back zero -- with no crash, no diagnostic, and a
  // perfectly plausible empty peak list.
  void finish() noexcept { finished_ = true; }

  // The highest base whose widest window is fully known.
  [[nodiscard]] Coord evaluable_limit() const noexcept {
    const std::int64_t by_contig = contig_len_ - 1;
    if (finished_) return static_cast<Coord>(by_contig);
    return static_cast<Coord>(
        std::min(frontier_ - kMaxHalfWidth + 1, by_contig));
  }

  // Advance one base. False once the evaluable limit is reached.
  bool step() noexcept {
    const std::int64_t limit = evaluable_limit();
    if (limit < 0 || pos_ >= limit) return false;
    do {
      advance_one();
    } while (pos_ < 0);
    return true;
  }

  [[nodiscard]] Coord position() const noexcept {
    return static_cast<Coord>(pos_ < 0 ? -1 : pos_);
  }

  [[nodiscard]] std::int64_t pileup() const noexcept { return pileup_at(pos_); }

  [[nodiscard]] std::int64_t window_sum(int scale) const noexcept {
    return sum_[scale];
  }

  // Truncated at both contig ends, so a base 100 bp from a telomere is divided
  // by the sequence that exists rather than by 10 kb of nothing.
  [[nodiscard]] std::int64_t window_width(int scale) const noexcept {
    const std::int64_t h = kScaleWidth[scale] / 2;
    const std::int64_t lo = std::max<std::int64_t>(0, pos_ - h);
    const std::int64_t hi = std::min<std::int64_t>(contig_len_, pos_ + h);
    return std::max<std::int64_t>(0, hi - lo);
  }

  [[nodiscard]] double background_rate(int scale) const noexcept {
    const std::int64_t w = window_width(scale);
    if (w <= 0) return 0.0;
    return static_cast<double>(sum_[scale]) / static_cast<double>(w);
  }

  // MACS2's lambda_local: the largest of the genome background and the three
  // local rates. Taking the max is what stops a copy-number amplification or a
  // mapping artefact from being called as a peak -- the background rises with
  // it.
  //
  // A saturated accumulator returns infinity, which log_poisson_sf treats as a
  // degenerate lambda and answers p = 1. Not a peak: the fail-closed direction.
  [[nodiscard]] double local_lambda(double lambda_bg) const noexcept {
    if (saturated_) return std::numeric_limits<double>::infinity();
    double lam = lambda_bg;
    for (int i = 0; i < kScaleCount; ++i) {
      const double r = background_rate(i);
      if (r > lam) lam = r;
    }
    return lam;
  }

  [[nodiscard]] std::int64_t rejected(AddResult why) const noexcept {
    return rejected_[static_cast<std::size_t>(why)];
  }

  [[nodiscard]] std::int64_t rejected_total() const noexcept {
    return rejected_[1] + rejected_[2] + rejected_[3];
  }

  // The highest fragment start the ring can accept right now. Bases below the
  // contig do not occupy ring slots -- pileup_at answers zero for them without
  // a load -- so the live low end is clamped at zero, which is what gives the
  // very first fragments enough room to push the frontier past the first
  // half-window and make stepping possible at all.
  [[nodiscard]] std::int64_t safe_start_limit() const noexcept {
    const std::int64_t low = std::max<std::int64_t>(0, pos_ - kMaxHalfWidth - 1);
    return low + kRingCapacity - kMaxFragmentSpan - 1;
  }

  [[nodiscard]] bool saturated() const noexcept { return saturated_; }

  [[nodiscard]] std::int64_t contig_length() const noexcept { return contig_len_; }

  // The last base whose window could still see a nonzero contribution from a
  // fragment already added -- past it, every window is provably empty until a
  // closer fragment arrives, which is exactly seek()'s own precondition.
  // -1 means no fragment has been added to this contig at all, so nothing is
  // protected and any target at or past base 0 is eligible.
  [[nodiscard]] std::int64_t signal_limit() const noexcept {
    if (max_end_ == 0) return -1;  // max_end_ == 0 is otherwise unreachable: every
                                   // accepted fragment has e >= 1 (add_fragment
                                   // requires e > s >= 0)
    return std::min<std::int64_t>(contig_len_ - 1, max_end_ + kMaxHalfWidth - 1);
  }

 private:
  [[nodiscard]] static std::size_t ring(std::int64_t p) noexcept {
    // p is never negative where this is called, but the mask is applied to the
    // unsigned form so a stray negative wraps into the ring instead of
    // indexing before it.
    return static_cast<std::size_t>(static_cast<std::uint64_t>(p) &
                                    static_cast<std::uint64_t>(kRingCapacity - 1));
  }

  // Positions outside the contig hold no signal. Returning zero for them is
  // what makes the rolling difference correct at both boundaries: a base
  // enters and leaves the window with the same value either way.
  [[nodiscard]] std::int64_t pileup_at(std::int64_t p) const noexcept {
    if (p < 0 || p >= contig_len_) return 0;
    return pileup_[ring(p)];
  }

  void accumulate(std::int64_t& acc, std::int64_t addend) noexcept {
    std::int64_t out = 0;
    if (__builtin_add_overflow(acc, addend, &out)) {
      saturated_ = true;
      acc = addend >= 0 ? std::numeric_limits<std::int64_t>::max()
                        : std::numeric_limits<std::int64_t>::min();
      return;
    }
    acc = out;
  }

  void advance_frontier_to(std::int64_t target) noexcept {
    while (frontier_ < target) {
      ++frontier_;
      const std::size_t idx = ring(frontier_);
      accumulate(running_, delta_[idx]);
      delta_[idx] = 0;
      pileup_[idx] = running_;
    }
  }

  // One base of the centred rolling update. Moving from p-1 to p, each window
  // gains the base at p + h - 1 and loses the one at p - h - 1.
  void advance_one() noexcept {
    ++pos_;
    // The widest window needs the base at pos + 5000 - 1. Before finish() the
    // frontier is already past it and this is a no-op; after finish() it is
    // what walks the pileup through the tail of the contig, one base at a
    // time, so the ring is never lapped.
    advance_frontier_to(pos_ + kMaxHalfWidth - 1);
    for (int i = 0; i < kScaleCount; ++i) {
      const std::int64_t h = kScaleWidth[i] / 2;
      accumulate(sum_[i], pileup_at(pos_ + h - 1));
      accumulate(sum_[i], -pileup_at(pos_ - h - 1));
    }
  }

  AddResult reject(AddResult why) noexcept {
    ++rejected_[static_cast<std::size_t>(why)];
    return why;
  }

  std::vector<std::int64_t> delta_;
  std::vector<std::int64_t> pileup_;
  std::int32_t tid_ = kNoTid;
  std::int64_t contig_len_ = 0;
  std::int64_t frontier_ = -1;
  std::int64_t running_ = 0;
  std::int64_t pos_ = -static_cast<std::int64_t>(kMaxHalfWidth) - 1;
  std::int64_t last_start_ = -1;
  std::int64_t max_end_ = 0;
  std::int64_t sum_[kScaleCount] = {};
  std::int64_t rejected_[5] = {};
  bool finished_ = false;
  bool saturated_ = false;
};

// The position a driver should seek() to in order to evaluate `candidate` --
// an incoming fragment's start minus the half-window, or a contig's own last
// base -- given `live`, or a negative value when no seek is available and the
// driver should drain normally instead. One decision, used by both the CLI
// driver (src/main.cpp) and its differential test (tests/test_caller.cpp),
// rather than two copies of the same checks.
//
// `live` is signal_limit() as it was captured by the CALLER at the right
// moment, not necessarily what MultiScaleWindow::signal_limit() reports right
// now -- on a fragment arrival, the caller must capture it BEFORE adding the
// fragment that triggered this decision, since that add already extends
// max_end_ (and so signal_limit()) past any real target, which would make
// this refuse unconditionally. Taking it as a parameter rather than deriving
// it from a MultiScaleWindow reference keeps that ordering a single explicit
// fact at the call site instead of a hidden requirement on when this runs.
//
// Callers own a second obligation this function cannot see or enforce: it
// answers "is `candidate` far enough to be worth a jump", not "is it SAFE to
// jump there right now". seek()'s own ring-wipe makes that a real question
// whenever `live >= 0` (something has already been added to this contig) --
// the caller must first drain up to `live` using whatever frontier is
// ALREADY available, and only proceed with seek() if that fully succeeds
// (MultiScaleWindow::position() reaches `live`). If it does not (the usual
// case for a fragment arriving after a real, narrow peak: nothing but that
// peak's own start has pushed the frontier, nowhere near `live`), seeking now
// would silently discard whatever of that peak was never stepped through --
// the caller must fall back to draining normally instead. This is why
// interior gaps between two dense, narrow peaks are not skipped by this
// design: only a LEADING gap (live < 0, nothing added yet, so there is
// nothing to lose) or a TRAILING one (this same reasoning, applied once at
// end of contig, without needing seek() at all -- see src/main.cpp) can
// satisfy that obligation for certain.
//
// `score_cutoff` must be strictly positive: at cutoff == 0, a fully empty
// (k=0) base's own score is exactly 0 (see detail::provably_below_cutoff),
// which meets "at or above the cutoff" -- observe()'s "a skipped base is
// below the cutoff by construction" contract does not hold there, so no base
// may ever be skipped.
[[nodiscard]] inline std::int64_t choose_seek_target(std::int64_t live, std::int64_t candidate,
                                                     std::int64_t contig_length,
                                                     std::int64_t current_position,
                                                     double score_cutoff) noexcept {
  if (!(score_cutoff > 0.0)) return -1;
  // seek()'s own precondition is candidate - max_end_ >= kMaxHalfWidth. When
  // live >= 0 that is exactly `candidate > live` (live == max_end_ +
  // kMaxHalfWidth - 1 by construction). When live < 0 (signal_limit()'s
  // sentinel for "nothing added to this contig yet", max_end_ == 0), the same
  // bound is candidate >= kMaxHalfWidth, i.e. candidate > kMaxHalfWidth - 1.
  const std::int64_t floor = (live < 0) ? (kMaxHalfWidth - 1) : live;
  if (candidate <= floor) return -1;             // seek() would refuse
  if (candidate >= contig_length) return -1;     // seek()'s own bound
  if (candidate <= current_position) return -1;  // already there
  return candidate;
}

}  // namespace pto::peaks
