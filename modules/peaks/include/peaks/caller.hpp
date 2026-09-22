// SPDX-License-Identifier: MIT
//
// pto-peaks: the peak state machine, the summit detector and Benjamini-Hochberg.
//
// The caller turns a per-base stream of (pileup, local lambda) into candidate
// peaks. Almost every subtlety in peak calling lives in the edges between its
// four states rather than in the p-value:
//
//   IDLE           nothing open. A base at or above the cutoff opens a peak.
//   IN_PEAK        above the cutoff and the summit is still improving.
//   SUMMIT_SEARCH  above the cutoff but past the current summit. The peak is
//                  still growing; the machine is now looking for either its
//                  end or a second, higher summit. A bimodal enhancer has two,
//                  and the summit must move to the taller one.
//   CLOSING        under the cutoff, but not yet for max_gap bases. Real
//                  signal is noisy and a single sub-threshold base in the
//                  middle of a peak is not the end of it. Once the drop
//                  persists, the peak is finalised and the machine returns to
//                  IDLE.
//
// The gap is measured in COORDINATES, not in calls to observe(). A driver may
// skip a stretch with no coverage at all -- MultiScaleWindow::seek() exists to
// do exactly that across a centromere -- and a peak must not be stitched
// across the skipped region merely because no base was reported inside it.
//
// Memory
// ------
// Nothing here allocates per base. The candidate vector is reserved once, and
// it grows only when a peak is FINALISED -- tens of thousands of times per
// genome, not billions. That is the bound the plan means by "peak memory
// remains bounded by O(pileup depth) rather than chromosome length": the ring
// in sliding_window.hpp is the part that must not grow with the chromosome,
// and the peak list is bounded by the number of peaks.
//
// Why the correction is not streaming
// -----------------------------------
// Benjamini-Hochberg is a rank statistic: the q-value of the best peak depends
// on how many peaks there are in total, so it cannot be computed until the
// last one has been seen. finalize() is therefore a second pass, over the
// candidate list only. Any "streaming FDR" would have to guess the final
// count, and guessing it low is the direction that invents peaks.
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

#include "peaks/poisson_model.hpp"
#include "peaks/sliding_window.hpp"

namespace pto::peaks {

// Exactly 32 bytes, and asserted below. A run of candidates stays in L1 while
// the coordinate sweeps past them, which is what the plan's "peak merging and
// summit tracking in L1 cache" is about.
//
// The q-value is deliberately NOT a member: it does not exist until every peak
// has been called, and reserving a slot for it here would put a field in the
// hot struct that is meaningless for the entire streaming pass. It lives in a
// parallel array, produced by finalize().
struct CandidatePeak {
  std::int32_t tid = kNoTid;
  Coord start = 0;
  Coord end = 0;          // half-open
  Coord summit_pos = 0;
  float summit_score = 0.0f;   // -log10 p at the summit
  float summit_lambda = 0.0f;  // the local background it was scored against
  std::int64_t total_signal = 0;
};

static_assert(sizeof(CandidatePeak) == 32,
              "CandidatePeak must stay 32 bytes: the cache-residency argument "
              "in this header and in implementation_plan.md depends on it");

enum class CallerState { kIdle, kInPeak, kSummitSearch, kClosing };

struct PeakCallOptions {
  // -log10 p required to open or extend a peak. 2.0 is p = 0.01, MACS2's
  // default for the initial scan; the FDR gate below is what actually decides
  // what survives.
  double score_cutoff = 2.0;
  // Calls shorter than this are dropped. Below about a fragment length a
  // "peak" is a single read.
  Coord min_length = 50;
  // How long a sub-threshold stretch may be before it ends a peak.
  Coord max_gap = 30;
  // Symmetric padding added to each finalised peak's REPORTED interval, in bp.
  // 0 (the default) is the sharp, high-resolution call: the enrichment core and
  // nothing more, which is what transcription-factor footprinting wants. A
  // positive value is a compatibility mode for pipelines that expect the wider
  // spans MACS2/MACS3 produce (`pto-core/modules/peaks/CONCORDANCE.md` records
  // the ~3x median-width gap this closes). It moves the bounds and NOTHING
  // else: the summit, its score, its background and the summed signal are the
  // measurement and are left exactly as called, and the padding is applied
  // AFTER the min_length gate so a short core is never padded into existence.
  // The 5' side is clamped at 0 here; the 3' side is clamped against the contig
  // length by the narrowPeak writer, which is the layer that knows it.
  Coord extend_peaks = 0;
  // Benjamini-Hochberg level applied by finalize().
  double qvalue_cutoff = 0.05;
  // Candidates reserved up front, so the streaming pass never allocates.
  std::size_t reserve_peaks = 1u << 16;
  // How many bases are scored per vectorised call. This is a throughput knob
  // and nothing else: every depth must produce the same peaks, which
  // `test_every_queue_depth_gives_the_same_peaks` checks. 256 keeps the queue
  // (three arrays of 8 bytes) inside L1 while being long enough that the
  // vector kernels reach steady state.
  std::size_t score_batch = 256;

  void validate() const {
    if (!(score_cutoff >= 0.0)) {  // also catches NaN
      throw std::invalid_argument("score cutoff must be a non-negative number");
    }
    if (min_length < 0) {
      throw std::invalid_argument("min peak length must not be negative");
    }
    if (max_gap < 0) {
      throw std::invalid_argument("max gap must not be negative");
    }
    if (extend_peaks < 0) {
      throw std::invalid_argument("peak extension must not be negative");
    }
    if (!(qvalue_cutoff > 0.0) || !(qvalue_cutoff <= 1.0)) {
      throw std::invalid_argument("q-value cutoff must be in (0, 1]");
    }
    if (score_batch == 0 || score_batch > (1u << 20)) {
      throw std::invalid_argument("score batch must be in [1, 1048576]");
    }
    // Passed straight to vector::reserve, so an absurd value is a length_error
    // out of a constructor rather than the argument error it actually is.
    if (reserve_peaks > (1u << 28)) {
      throw std::invalid_argument("reserved peak count must be at most 268435456");
    }
  }
};

namespace detail {

// True when this base cannot possibly clear the cutoff, so it need not be
// scored at all.
//
// The bound is a theorem, not a heuristic, which is the only kind of prefilter
// worth having: Choi (1994) shows the median m of a Poisson(lambda) satisfies
// m >= lambda - ln 2, so for k <= lambda - ln 2 at least half the distribution
// lies at or above k, P(X >= k) >= 1/2, and -log10 p <= 0.30103. A cutoff
// above that can never be met. Below it the filter refuses to fire.
//
// NaN falls through to the full evaluation on purpose: every comparison
// against it is false, so a degenerate lambda is not silently classified here
// -- log_poisson_sf's own fail-closed guard answers it.
[[nodiscard]] inline bool provably_below_cutoff(std::uint32_t k, double lambda,
                                                double cutoff) noexcept {
  constexpr double kLn2 = 0.6931471805599453;
  constexpr double kScoreAtMedian = 0.3010299956639812;  // -log10(1/2)
  if (k == 0) return cutoff > 0.0;  // P(X >= 0) = 1 exactly
  if (!(cutoff > kScoreAtMedian)) return false;
  return static_cast<double>(k) <= lambda - kLn2;
}

// Choi's bound above fires below the background (k <= lambda - ln2); a
// background track's modal state is k ~= lambda, exactly where it never
// fires, so every such base pays a full Poisson evaluation just to be told
// it scores close to 0.9. This is the complementary prefilter: an open
// interval of lambda, per k, where poisson_pmf_upper_bound (poisson_model.hpp)
// already proves the exact score cannot reach the cutoff -- and that bound is
// smallest exactly at lambda = k, precisely where Choi's is loosest.
//
// lo > hi (the default) means empty: no lambda is ever admitted. Exactly 16
// bytes so a table of these (kSkipBandMaxK + 1 entries) stays a few KiB.
struct SkipBand {
  double lo = 1.0, hi = 0.0;
};

// Bases past this k always fall through to Choi's bound and the exact
// evaluation; U's minimum (at lambda = k) grows with k roughly like
// k*ln(k)/2, so realistic background pileup depths are covered well below it,
// and the per-run table this bounds (PeakCaller's skip_bands_) stays a few
// KiB and a sub-10ms one-time build regardless.
inline constexpr std::uint32_t kSkipBandMaxK = 1024;

// Solves poisson_pmf_upper_bound(k, lambda) == cutoff for lambda, bisecting
// between `below` (where U(k, below) <= cutoff) and `above` (where
// U(k, above) > cutoff) -- named for which side of cutoff they are on, not
// for numeric order, since U(k, .) is convex with one minimum (poisson_model.hpp)
// and this is used on BOTH sides of it: `above` is smaller than `below` for
// the left root (U decreasing there) and larger for the right one (U
// increasing). Classifying `mid` and keeping the invariant on whichever
// endpoint it replaces works identically either way; assuming a fixed
// direction here would not -- the first version of this function did, and
// silently returned a right root a few hundred ULPs on the WRONG side of
// lambda = k, caught by the differential test below.
//
// 200 iterations halves the bracket past 2^-200 of its start, converging to
// the same double from either side long before that budget is spent -- the
// same fixed-iteration approach genomic_toolkit's estimate_library_size()
// uses for the same reason: provably enough beats a tolerance check that
// could spin on a pathological bracket.
[[nodiscard]] inline double solve_u_equals_cutoff(std::uint32_t k, double cutoff, double below,
                                                  double above) noexcept {
  for (int i = 0; i < 200; ++i) {
    const double mid = 0.5 * (below + above);
    if (poisson_pmf_upper_bound(k, mid) <= cutoff) {
      below = mid;
    } else {
      above = mid;
    }
  }
  return 0.5 * (below + above);
}

// U(k, .) is strictly convex with its minimum at lambda = k (poisson_model.hpp).
// Finds the two roots of U(k, lambda) == cutoff on either side of that
// minimum, or an empty band if even the minimum does not clear the cutoff.
// Called once per k, at PeakCaller construction, for the run's own cutoff --
// never on the per-base path.
[[nodiscard]] inline SkipBand compute_skip_band(std::uint32_t k, double cutoff) noexcept {
  if (k == 0) return SkipBand{};  // provably_below_cutoff's own case, exactly, not a bound
  const double lambda_min = static_cast<double>(k);
  const double u_min = poisson_pmf_upper_bound(k, lambda_min);
  if (!(u_min < cutoff)) return SkipBand{};  // even the best case does not clear it

  // Left root: U -> +infinity as lambda -> 0+ (the -k*ln(lambda) term), so
  // any lambda approaching 0 brackets it together with lambda_min (where
  // U == u_min <= cutoff). 1e-300 is comfortably inside log_approx's domain
  // and gives U a value in the billions for any k this table reaches.
  const double lo = solve_u_equals_cutoff(k, cutoff, lambda_min, 1e-300);

  // Right root: U -> +infinity as lambda -> infinity too (the +lambda term),
  // so doubling the bracket must eventually clear the cutoff again.
  double hi_bracket = lambda_min * 2.0 + 1.0;
  while (poisson_pmf_upper_bound(k, hi_bracket) < cutoff) hi_bracket *= 2.0;
  const double hi = solve_u_equals_cutoff(k, cutoff, lambda_min, hi_bracket);

  // solve_u_equals_cutoff's 200 iterations converge lo/hi to within a few
  // ULPs of the true root, where floating-point error alone could put the
  // COMPUTED U on either side of cutoff -- shrinking the interval very
  // slightly toward its own interior guarantees every lambda strictly inside
  // (lo, hi) scores under the cutoff bound, never merely at it. The margin
  // this trades away is nowhere near load-bearing: the tightest
  // cutoff-minus-score gap measured across the whole (k, lambda) grid this
  // table covers is twelve orders of magnitude above double rounding.
  SkipBand band;
  band.lo = lo * (1.0 + 1e-9);
  band.hi = hi * (1.0 - 1e-9);
  return band;
}

}  // namespace detail

// Benjamini-Hochberg, in -log10 space throughout.
//
// `scores` are -log10 p; `out` receives -log10 q for the same positions. The
// whole procedure is done on logarithms because a real summit reaches
// -log10 p of several thousand, and p itself is then zero in double -- every
// such peak would get the same adjusted value, and their ranking, which is the
// only thing a q-value is for, would be destroyed.
//
// The step-up is the standard one: with m tests and the j-th smallest p-value,
// the adjustment is p*m/j, made monotone by taking the running minimum from
// the largest rank down (a running MAXIMUM here, since the sign is flipped),
// and capped at q = 1, which is a floor of zero on this scale.
inline void benjamini_hochberg(std::span<const float> scores, std::vector<float>& out) {
  const std::size_t m = scores.size();
  out.assign(m, 0.0f);
  if (m == 0) return;

  std::vector<std::uint32_t> order(m);
  std::iota(order.begin(), order.end(), 0u);
  // Descending score is ascending p-value, and equal scores keep their genomic
  // order so a run of tied peaks is ranked the way it was called.
  //
  // The tie-break is in the COMPARATOR, not in the algorithm, and that is the
  // whole reason this is `std::sort`. Ordering by index when the scores are
  // equal makes the order total, which is what `stable_sort` was being asked
  // for -- and libstdc++'s stable_sort buys that stability with a
  // `_Temporary_buffer`, an allocation it takes from
  // `operator new(size_t, nothrow_t)` and returns through
  // `operator delete(void*)`. A translation unit that replaces one of those
  // and not the other -- which is every test binary in this module, and which
  // ASan reports as `alloc-dealloc-mismatch` -- then aborts here rather than
  // in the code that made the mistake. Sorting indices with a total order
  // needs no buffer and no allocator at all.
  std::sort(order.begin(), order.end(), [&](std::uint32_t a, std::uint32_t b) {
    if (scores[a] != scores[b]) return scores[a] > scores[b];
    return a < b;
  });

  const double log10_m = std::log10(static_cast<double>(m));
  double running = 0.0;  // -log10 of the running minimum adjusted p
  for (std::size_t rank = m; rank-- > 0;) {
    const std::uint32_t i = order[rank];
    // rank is 0-based here; the procedure's j is rank + 1.
    const double adjusted = static_cast<double>(scores[i]) - log10_m +
                            std::log10(static_cast<double>(rank + 1));
    const double capped = adjusted > 0.0 ? adjusted : 0.0;  // q <= 1
    if (capped > running) running = capped;
    out[i] = static_cast<float>(running);
  }
}

class PeakCaller {
 public:
  explicit PeakCaller(const PeakCallOptions& opts) : opts_(opts) {
    opts_.validate();
    peaks_.reserve(opts_.reserve_peaks);
    // Every queue is sized once, here. Nothing below this line allocates.
    q_pos_.resize(opts_.score_batch);
    q_pile_.resize(opts_.score_batch);
    q_lambda_.resize(opts_.score_batch);
    q_score_.resize(opts_.score_batch);
    cand_k_.resize(opts_.score_batch);
    // One-time, sub-10ms table build (see kSkipBandMaxK) for this run's own
    // cutoff -- never touched again on the per-base path.
    skip_bands_.resize(detail::kSkipBandMaxK + 1);
    for (std::uint32_t k = 0; k <= detail::kSkipBandMaxK; ++k) {
      skip_bands_[k] = detail::compute_skip_band(k, opts_.score_cutoff);
    }
  }

  // Starts a contig, flushing anything still open on the previous one. A peak
  // may never span a boundary however the signal looks: the coordinates on
  // either side are unrelated.
  void begin_contig(std::int32_t tid) {
    end_contig();
    tid_ = tid;
    last_pos_ = -1;
  }

  // Enqueues a base. The queue is scored in one vectorised call when it fills,
  // and the state machine then walks the results in coordinate order.
  //
  // The machine is NOT vectorised and cannot be: the state at base i+1 depends
  // on the state at base i, and no masking changes that. What vectorises is
  // the Poisson evaluation, which is where the time goes, and running it ahead
  // of the machine over a fixed queue is how the vector units stay busy
  // without changing a single transition.
  //
  // `pos` must be non-decreasing within a contig. A driver is free to SKIP
  // bases -- that is what the gap accounting in step_machine is for -- but not
  // to go backwards. MultiScaleWindow::step() and seek() both satisfy this;
  // extend() clamps rather than trusting it, because a peak whose end fell
  // below its start would be silently swallowed by flush()'s length test and
  // counted as too short.
  //
  // A base that either prefilter can PROVE cannot clear the cutoff is not
  // queued at all, exactly as if the driver itself had skipped that
  // coordinate -- step_machine's gap accounting (see the comment there) treats
  // an unobserved coordinate range as a below-cutoff run already, which is
  // provably what this base's own exact score would have been, so the two are
  // indistinguishable to the machine. This is the same argument seek() relies
  // on for coverage deserts, applied one base at a time instead of one region
  // at a time: it removes both the queue slot AND the eventual step_machine
  // call, where the earlier scheme (score it as 0.0, still walk the machine
  // over it) removed only the Poisson evaluation.
  void observe(Coord pos, std::int64_t pileup, double lambda) {
    const std::uint32_t k = clamp_count(pileup);
    if (is_provably_below_cutoff(k, lambda)) return;
    q_pos_[queued_] = pos;
    q_pile_[queued_] = pileup;
    q_lambda_[queued_] = lambda;
    cand_k_[queued_] = k;
    if (++queued_ == opts_.score_batch) drain_queue();
  }

  void drain_queue() {
    if (queued_ == 0) return;
    // Every queued base already survived both prefilters in observe(), so
    // the whole queue is scored in one vectorised call -- no compaction step
    // is needed here any more.
    neg_log10_poisson_sf_batch(cand_k_.data(), q_lambda_.data(), q_score_.data(), queued_);

    // The sequential machine, over precomputed scores, in order.
    const std::size_t n = queued_;
    queued_ = 0;  // cleared first: step_machine can flush a peak, never re-enter
    for (std::size_t i = 0; i < n; ++i) {
      step_machine(q_pos_[i], q_pile_[i], q_lambda_[i], q_score_[i]);
    }
  }

 private:
  void step_machine(Coord pos, std::int64_t pileup, double lambda, double score) {
    // A driver may skip bases with no coverage rather than reporting each of
    // them. Those bases are under the cutoff by construction, so the skip is a
    // gap and is measured the same way -- but the every-base stream this must
    // reproduce checks `above` BEFORE it checks gap length (the switch below,
    // CLOSING branch), so a real base landing exactly on the coordinate where
    // the gap would have closed still gets to resume the peak if it clears the
    // cutoff. Force-flushing here on coordinate arithmetic alone, with no
    // chance for that base's own score to be seen first, stitches the peak
    // short of where an every-base stream would have left it.
    //
    // The fix is to flush now only when the closing coordinate falls STRICTLY
    // before `pos` -- i.e. among bases that were actually skipped, which are
    // below the cutoff by the contract above, so the every-base stream
    // provably flushed there and stayed IDLE the rest of the way to `pos`.
    // When `pos` IS the closing coordinate (or earlier), `pos`'s real score
    // has not been consulted yet: fall through to the switch below with
    // gap_start_ set correctly, and let its own `above`-first check decide,
    // exactly as it would have for that base in an every-base stream.
    if (last_pos_ >= 0 && pos > last_pos_ + 1 && state_ != CallerState::kIdle) {
      const bool closing = (state_ == CallerState::kClosing);
      const std::int64_t gap_from = closing ? gap_start_ : last_pos_ + 1;
      // The base that OPENS a gap (sets gap_start_) is never itself checked
      // against max_gap -- only the CLOSING branch does that, and it is first
      // reached on the base after. So the earliest coordinate the every-base
      // stream could flush at is one past a fresh transition's gap_from, not
      // gap_from + max_gap directly: at max_gap == 0 those differ, since
      // gap_from + max_gap (== gap_from) is never reachable at all -- pos is
      // already >= last_pos_ + 2 here, one past gap_from = last_pos_ + 1. An
      // already-CLOSING run has no such gap: gap_start_ was set on a real,
      // already-checked base, so gap_start_ + max_gap always dominates.
      const std::int64_t earliest = last_pos_ + (closing ? 1 : 2);
      const std::int64_t closes_at = std::max(earliest, gap_from + opts_.max_gap);
      if (pos > closes_at) {
        flush();
      } else if (!closing) {
        gap_start_ = gap_from;
        state_ = CallerState::kClosing;
      }
    }
    last_pos_ = pos;

    const bool above = score >= opts_.score_cutoff;

    switch (state_) {
      case CallerState::kIdle:
        if (above) open(pos, score, lambda, pileup);
        break;
      case CallerState::kInPeak:
      case CallerState::kSummitSearch:
        if (above) {
          extend(pos, score, lambda, pileup);
        } else {
          gap_start_ = pos;
          state_ = CallerState::kClosing;
        }
        break;
      case CallerState::kClosing:
        if (above) {
          // The sub-threshold stretch was short enough to be interior: the
          // peak resumes and swallows it.
          extend(pos, score, lambda, pileup);
        } else if (static_cast<std::int64_t>(pos) - gap_start_ >= opts_.max_gap) {
          flush();
        }
        break;
    }
  }

 public:
  // Flushes an open peak. A peak that runs to the last base of a chromosome is
  // still a peak, and it is the one a "close when the coordinate passes
  // end + gap" rule drops, because that coordinate never arrives.
  void end_contig() {
    drain_queue();
    flush();
    state_ = CallerState::kIdle;
    last_pos_ = -1;
  }

  // The machine's state after the bases it has PROCESSED, which is not
  // necessarily every base observed: enqueued bases have not reached it yet.
  // Call drain_queue() first to make an observation per-base meaningful, or
  // set score_batch to 1. Nothing in production reads this per base; it is
  // here because the four-state walk is worth pinning.
  [[nodiscard]] CallerState state() const noexcept { return state_; }
  [[nodiscard]] const std::vector<CandidatePeak>& peaks() const noexcept { return peaks_; }
  [[nodiscard]] const std::vector<float>& qscores() const noexcept { return qscores_; }
  [[nodiscard]] std::int64_t discarded_short() const noexcept { return discarded_short_; }
  [[nodiscard]] bool saturated() const noexcept { return saturated_; }

  // Multiple-testing correction over every candidate, then the FDR gate.
  // Both vectors come back the same length and in genomic order.
  void finalize() {
    drain_queue();
    flush();
    scratch_.clear();
    scratch_.reserve(peaks_.size());
    for (const CandidatePeak& p : peaks_) scratch_.push_back(p.summit_score);
    benjamini_hochberg(scratch_, qscores_);

    const float gate = static_cast<float>(-std::log10(opts_.qvalue_cutoff));
    std::size_t keep = 0;
    for (std::size_t i = 0; i < peaks_.size(); ++i) {
      if (qscores_[i] >= gate) {
        peaks_[keep] = peaks_[i];
        qscores_[keep] = qscores_[i];
        ++keep;
      }
    }
    peaks_.resize(keep);
    qscores_.resize(keep);
  }

 private:
  // The Poisson model is indexed by an unsigned 32-bit count. A pileup past
  // that is four billion fragments on one base -- not data, but clamping is
  // still the right answer: the score saturates at "certain", which is what
  // such a position is, and `saturated()` records that a number was pinned.
  std::uint32_t clamp_count(std::int64_t k) noexcept {
    if (k <= 0) return 0u;
    constexpr std::int64_t kMax = std::numeric_limits<std::uint32_t>::max();
    if (k >= kMax) {
      saturated_ = true;
      return static_cast<std::uint32_t>(kMax);
    }
    return static_cast<std::uint32_t>(k);
  }

  // True when either prefilter can PROVE this base cannot clear the cutoff:
  // Choi's bound (k depressed relative to lambda) or the skip band (lambda
  // close enough to k that even the tightest single-term bound on the exact
  // score stays under the cutoff). Checked in observe(), before a base ever
  // reaches the queue.
  [[nodiscard]] bool is_provably_below_cutoff(std::uint32_t k, double lambda) const noexcept {
    if (detail::provably_below_cutoff(k, lambda, opts_.score_cutoff)) return true;
    if (k <= detail::kSkipBandMaxK) {
      const detail::SkipBand& band = skip_bands_[k];
      if (lambda > band.lo && lambda < band.hi) return true;
    }
    return false;
  }

  void open(Coord pos, double score, double lambda, std::int64_t pileup) noexcept {
    cur_ = CandidatePeak{};
    cur_.tid = tid_;
    cur_.start = pos;
    cur_.end = next_coord(pos);
    cur_.summit_pos = pos;
    cur_.summit_score = static_cast<float>(score);
    cur_.summit_lambda = static_cast<float>(lambda);
    cur_.total_signal = 0;
    summit_lo_ = pos;
    summit_hi_ = pos;
    add_signal(pileup);
    open_ = true;
    state_ = CallerState::kInPeak;
  }

  // The comparison is done on the FLOAT, deliberately.
  //
  // Comparing the incoming double against `static_cast<double>(summit_score)`
  // looks equivalent and is not: the stored value has been rounded to float,
  // so a genuinely equal score compares GREATER than its own rounded copy and
  // every tie reads as an improvement. On a flat peak -- which is what a
  // uniform pileup over a fragment-length plateau is, and what most real peaks
  // look like at their top -- the summit then walks to the last base of the
  // plateau, every time. It was found by running the CLI on planted data and
  // reading a summit offset of 199 on a 200 bp peak.
  //
  // Ties instead extend a plateau, and the summit reported is its midpoint.
  // That is what a symmetric pileau deserves, and it is stable: adding one
  // base to either end of the plateau moves the summit by half a base, not by
  // the whole width.
  void extend(Coord pos, double score, double lambda, std::int64_t pileup) noexcept {
    cur_.end = std::max(cur_.end, next_coord(pos));
    add_signal(pileup);
    const float f = static_cast<float>(score);
    if (f > cur_.summit_score) {
      cur_.summit_score = f;
      cur_.summit_lambda = static_cast<float>(lambda);
      summit_lo_ = pos;
      summit_hi_ = pos;
      state_ = CallerState::kInPeak;
    } else if (f == cur_.summit_score &&
               summit_hi_ == static_cast<std::int64_t>(pos) - 1) {
      summit_hi_ = pos;  // the plateau continues
      state_ = CallerState::kInPeak;
    } else {
      state_ = CallerState::kSummitSearch;
    }
  }

  // `pos + 1` where pos is a Coord is an int32 addition, and at the top of the
  // range it is undefined behaviour that happens before any test of the result
  // can look at it -- sliding_window.hpp's header names this as the defect
  // docs/REVIEW_2026-08-15.md finding (b) records. The window clamps position()
  // to contig_len - 1 so the CLI cannot reach it, but observe() is public and
  // this header is installed, so the arithmetic is done in int64 and saturated
  // rather than left to a precondition nothing enforces.
  [[nodiscard]] static Coord next_coord(Coord pos) noexcept {
    const std::int64_t next = static_cast<std::int64_t>(pos) + 1;
    constexpr std::int64_t kCeiling = std::numeric_limits<Coord>::max();
    return static_cast<Coord>(next > kCeiling ? kCeiling : next);
  }

  // Checked, for the reason docs/REVIEW_2026-08-15.md finding (b) gives: a
  // signed overflow is undefined behaviour that happens before any test of the
  // result can run, so it is caught at the addition rather than after it.
  void add_signal(std::int64_t pileup) noexcept {
    std::int64_t out = 0;
    if (__builtin_add_overflow(cur_.total_signal, pileup, &out)) {
      saturated_ = true;
      cur_.total_signal = std::numeric_limits<std::int64_t>::max();
      return;
    }
    cur_.total_signal = out;
  }

  // Widens a finalised peak's reported interval symmetrically by
  // opts_.extend_peaks on each side, clamping the 5' end at 0. Called only from
  // flush(), only after the min_length gate, and only on a peak that is about
  // to be kept. Touches the bounds and nothing else -- see PeakCallOptions.
  void apply_extension(CandidatePeak& p) const noexcept {
    const std::int64_t ext = static_cast<std::int64_t>(opts_.extend_peaks);
    const std::int64_t lo = static_cast<std::int64_t>(p.start) - ext;
    const std::int64_t hi = static_cast<std::int64_t>(p.end) + ext;
    constexpr std::int64_t kCeiling = std::numeric_limits<Coord>::max();
    p.start = static_cast<Coord>(lo < 0 ? 0 : lo);
    p.end = static_cast<Coord>(hi > kCeiling ? kCeiling : hi);
  }

  void flush() {
    if (open_) {
      // The midpoint of the highest plateau. For a single maximum the two
      // ends coincide and this is that base.
      cur_.summit_pos = static_cast<Coord>(
          summit_lo_ + (summit_hi_ - summit_lo_) / 2);
      if (cur_.end - cur_.start >= opts_.min_length) {
        // The length gate has passed on the enrichment core; only now does the
        // optional compatibility padding move the reported bounds.
        if (opts_.extend_peaks > 0) apply_extension(cur_);
        peaks_.push_back(cur_);
      } else {
        ++discarded_short_;
      }
      open_ = false;
    }
    state_ = CallerState::kIdle;
  }

  PeakCallOptions opts_;
  // The score queue. Three inputs, one output, and a compacted candidate view
  // of them; all sized once in the constructor.
  std::vector<Coord> q_pos_;
  std::vector<std::int64_t> q_pile_;
  std::vector<double> q_lambda_;
  std::vector<double> q_score_;
  // k, precomputed in observe() (clamp_count is already paid there to run the
  // prefilters), aligned 1:1 with q_pos_/q_pile_/q_lambda_.
  std::vector<std::uint32_t> cand_k_;
  // One entry per k in [0, kSkipBandMaxK], built once at construction for
  // this run's score_cutoff. See detail::SkipBand.
  std::vector<detail::SkipBand> skip_bands_;
  std::size_t queued_ = 0;
  std::vector<CandidatePeak> peaks_;
  std::vector<float> qscores_;
  std::vector<float> scratch_;
  CandidatePeak cur_{};
  CallerState state_ = CallerState::kIdle;
  std::int32_t tid_ = kNoTid;
  std::int64_t last_pos_ = -1;
  std::int64_t gap_start_ = 0;
  std::int64_t summit_lo_ = 0;
  std::int64_t summit_hi_ = 0;
  std::int64_t discarded_short_ = 0;
  bool open_ = false;
  bool saturated_ = false;
};

}  // namespace pto::peaks
