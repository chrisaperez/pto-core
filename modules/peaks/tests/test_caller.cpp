// SPDX-License-Identifier: MIT
//
// pto-peaks: the peak state machine, the summit detector and Benjamini-Hochberg.
//
// The caller turns a per-base stream of (pileup, local lambda) into a list of
// candidate peaks. It is a four-state machine -- IDLE, IN_PEAK, SUMMIT_SEARCH,
// CLOSING -- and almost every subtlety in peak calling lives in the edges
// between those states rather than in the p-value:
//
//   * a peak is not closed by the first base that drops under the cutoff, or
//     every noisy summit becomes a string of fragments;
//   * it IS closed once the drop persists for max_gap, or two genes' promoters
//     merge into one call;
//   * the summit is the highest-scoring base, and when a second, higher one
//     appears later inside the same peak the summit must move to it;
//   * a peak may never cross a contig boundary, however the signal looks.
//
// The q-value half is tested against a direct implementation of the textbook
// step-up procedure rather than against itself, and separately for the
// property that survives extreme input: monotonicity.
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <random>
#include <stdexcept>
#include <vector>

#include "peaks_test_util.hpp"

#include "peaks/caller.hpp"
#include "peaks/poisson_model.hpp"
#include "peaks/sliding_window.hpp"

namespace {

using pto::peaks::CallerState;
using pto::peaks::CandidatePeak;
using pto::peaks::Coord;
using pto::peaks::PeakCallOptions;
using pto::peaks::PeakCaller;

// Options with the gates opened, so a case can exercise one rule at a time
// without the others discarding its peaks first.
PeakCallOptions permissive() {
    PeakCallOptions o;
    o.score_cutoff = 2.0;
    o.min_length = 1;
    o.max_gap = 10;
    o.qvalue_cutoff = 1.0;
    o.reserve_peaks = 256;
    return o;
}

// A pileup that scores just above / just below the cutoff against lambda = 1.
// Chosen by search rather than by hand so the intent survives a change to the
// scoring function.
std::int64_t pileup_scoring(double target, double lambda) {
    for (std::int64_t k = 1; k < 100000; ++k) {
        if (pto::peaks::neg_log10_poisson_sf(static_cast<std::uint32_t>(k), lambda) >= target) {
            return k;
        }
    }
    return 100000;
}

// ---------------------------------------------------------------------------

void test_candidate_peak_is_thirty_two_bytes() {
    // The plan sizes this struct so a run of candidates stays in L1 while the
    // coordinate sweeps past them. Growing it is allowed, but not by accident.
    CHECK_EQ(sizeof(CandidatePeak), 32u);
    CHECK(std::is_trivially_copyable_v<CandidatePeak>);
}

void test_idle_stays_idle_below_the_cutoff() {
    PeakCaller c(permissive());
    c.begin_contig(0);
    for (Coord p = 0; p < 500; ++p) c.observe(p, 1, 10.0);
    c.end_contig();
    CHECK(c.state() == CallerState::kIdle);
    CHECK_EQ(c.peaks().size(), 0u);
}

void test_the_state_machine_walks_all_four_states() {
    // IDLE -> IN_PEAK (score rising) -> SUMMIT_SEARCH (score past its maximum)
    // -> CLOSING (under the cutoff, gap not yet spent) -> IDLE (peak emitted).
    auto opts = permissive();
    opts.max_gap = 5;
    // Depth 1: state() reports the machine after the bases it has processed,
    // and this case is about the transitions between them, so each observe()
    // has to reach the machine before the next assertion.
    opts.score_batch = 1;
    // Below kScoreAtMedian (0.30103, caller.hpp) Choi's bound can never fire
    // (its own guard refuses any cutoff at or under it), and the skip band's
    // own minimum for k=1 (~0.434, poisson_pmf_upper_bound(1,1)) sits above
    // this cutoff too, so it is empty. The below-cutoff sentinel below (k=1,
    // lambda=1.0) is real background under EITHER cutoff; this one just keeps
    // it out of both prefilters, so it still reaches step_machine directly
    // and the per-call state() transitions below stay observable one base at
    // a time, exactly as they were before caller.hpp started dropping
    // provably-below-cutoff bases before they are ever queued. The "hi"/"mid"
    // targets are picked by pileup_scoring against a fixed score target
    // (8.0/4.0), not against opts.score_cutoff, so lowering the cutoff here
    // does not change what they mean.
    opts.score_cutoff = 0.25;
    PeakCaller c(opts);
    c.begin_contig(0);

    const std::int64_t hi = pileup_scoring(8.0, 1.0);
    const std::int64_t mid = pileup_scoring(4.0, 1.0);

    Coord p = 0;
    c.observe(p++, 1, 1.0);
    CHECK(c.state() == CallerState::kIdle);

    c.observe(p++, mid, 1.0);
    CHECK(c.state() == CallerState::kInPeak);

    c.observe(p++, hi, 1.0);  // still improving
    CHECK(c.state() == CallerState::kInPeak);

    c.observe(p++, mid, 1.0);  // above cutoff, past the summit
    CHECK(c.state() == CallerState::kSummitSearch);

    c.observe(p++, 1, 1.0);  // under the cutoff, gap opens
    CHECK(c.state() == CallerState::kClosing);

    for (int i = 0; i < 5; ++i) c.observe(p++, 1, 1.0);
    CHECK(c.state() == CallerState::kIdle);
    CHECK_EQ(c.peaks().size(), 1u);
}

void test_summit_is_the_highest_scoring_base() {
    PeakCaller c(permissive());
    c.begin_contig(7);
    const std::int64_t lo = pileup_scoring(3.0, 1.0);
    const std::int64_t hi = pileup_scoring(9.0, 1.0);
    Coord p = 1000;
    c.observe(p++, lo, 1.0);
    c.observe(p++, lo, 1.0);
    const Coord summit = p;
    c.observe(p++, hi, 1.0);
    c.observe(p++, lo, 1.0);
    c.end_contig();

    CHECK_EQ(c.peaks().size(), 1u);
    if (!c.peaks().empty()) {
        const CandidatePeak& pk = c.peaks()[0];
        CHECK_EQ(pk.tid, 7);
        CHECK_EQ(pk.start, 1000);
        CHECK_EQ(pk.end, p);
        CHECK_EQ(pk.summit_pos, summit);
        CHECK(pk.summit_score > 8.0f);
    }
}

void test_a_later_higher_summit_moves_the_summit() {
    // The reason SUMMIT_SEARCH exists as a state rather than as a comparison:
    // once the score has come off its maximum the machine is looking for
    // either the end of the peak or a better summit, and both must be possible.
    auto opts = permissive();
    opts.max_gap = 20;
    opts.score_batch = 1;  // see the four-state walk: state() is per processed base
    PeakCaller c(opts);
    c.begin_contig(0);
    const std::int64_t a = pileup_scoring(5.0, 1.0);
    const std::int64_t b = pileup_scoring(12.0, 1.0);

    Coord p = 0;
    c.observe(p++, a, 1.0);
    const Coord first_summit = 0;
    c.observe(p++, a - 1, 1.0);
    CHECK(c.state() == CallerState::kSummitSearch);
    const Coord second_summit = p;
    c.observe(p++, b, 1.0);
    CHECK(c.state() == CallerState::kInPeak);
    c.end_contig();

    CHECK_EQ(c.peaks().size(), 1u);
    if (!c.peaks().empty()) {
        CHECK_EQ(c.peaks()[0].summit_pos, second_summit);
        CHECK(c.peaks()[0].summit_pos != first_summit);
    }
}

void test_a_flat_plateau_puts_the_summit_at_its_middle() {
    // Most real peaks are flat at the top: a uniform pileup over a
    // fragment-length plateau. Comparing the incoming double score against the
    // float it was stored as makes every tie look like an improvement, and the
    // summit then walks to the LAST base of the plateau. This is the case that
    // catches it -- found by running the CLI on planted data and reading a
    // summit offset of 199 on a 200 bp peak.
    PeakCaller c(permissive());
    c.begin_contig(0);
    const std::int64_t flat = pileup_scoring(9.0, 1.0);
    const Coord first = 1000;
    const Coord last = 1099;
    for (Coord p = first; p <= last; ++p) c.observe(p, flat, 1.0);
    c.end_contig();

    CHECK_EQ(c.peaks().size(), 1u);
    if (!c.peaks().empty()) {
        const Coord summit = c.peaks()[0].summit_pos;
        CHECK_MSG(summit == first + (last - first) / 2,
                  "summit at " + std::to_string(summit) + ", plateau [" +
                      std::to_string(first) + "," + std::to_string(last) + "]");
        CHECK(summit != last);
        CHECK(summit != first);
    }
}

void test_a_short_gap_does_not_split_a_peak() {
    auto opts = permissive();
    opts.max_gap = 30;
    opts.score_batch = 1;  // see the four-state walk: state() is per processed base
    // See the four-state walk's comment: below kScoreAtMedian neither
    // prefilter can fire, so the (k=1, lambda=1.0) sentinel below still
    // reaches step_machine directly and CLOSING stays observable per base.
    opts.score_cutoff = 0.25;
    PeakCaller c(opts);
    c.begin_contig(0);
    const std::int64_t hi = pileup_scoring(6.0, 1.0);

    Coord p = 100;
    for (int i = 0; i < 10; ++i) c.observe(p++, hi, 1.0);
    for (int i = 0; i < 20; ++i) c.observe(p++, 1, 1.0);  // 20 < max_gap
    CHECK(c.state() == CallerState::kClosing);
    for (int i = 0; i < 10; ++i) c.observe(p++, hi, 1.0);
    CHECK(c.state() != CallerState::kIdle);
    c.end_contig();

    CHECK_EQ(c.peaks().size(), 1u);
    if (!c.peaks().empty()) {
        CHECK_EQ(c.peaks()[0].start, 100);
        CHECK_EQ(c.peaks()[0].end, p);  // the gap is interior to the peak
    }
}

void test_a_long_gap_splits_into_two_peaks() {
    auto opts = permissive();
    opts.max_gap = 30;
    PeakCaller c(opts);
    c.begin_contig(0);
    const std::int64_t hi = pileup_scoring(6.0, 1.0);

    Coord p = 100;
    for (int i = 0; i < 10; ++i) c.observe(p++, hi, 1.0);
    const Coord first_end = p;
    for (int i = 0; i < 40; ++i) c.observe(p++, 1, 1.0);  // 40 >= max_gap
    const Coord second_start = p;
    for (int i = 0; i < 10; ++i) c.observe(p++, hi, 1.0);
    c.end_contig();

    CHECK_EQ(c.peaks().size(), 2u);
    if (c.peaks().size() == 2) {
        CHECK_EQ(c.peaks()[0].start, 100);
        CHECK_EQ(c.peaks()[0].end, first_end);
        CHECK_EQ(c.peaks()[1].start, second_start);
    }
}

void test_a_skipped_coordinate_range_counts_as_a_gap() {
    // The driver may jump over a stretch with no coverage at all rather than
    // calling observe() a million times. The gap that closes a peak is
    // therefore measured in COORDINATES, not in calls.
    auto opts = permissive();
    opts.max_gap = 30;
    PeakCaller c(opts);
    c.begin_contig(0);
    const std::int64_t hi = pileup_scoring(6.0, 1.0);
    for (Coord p = 100; p < 110; ++p) c.observe(p, hi, 1.0);
    c.observe(100000, hi, 1.0);  // an enormous jump
    c.end_contig();
    CHECK_EQ(c.peaks().size(), 2u);
}

// A driver that skips a sub-threshold run must produce EXACTLY the peaks an
// every-base stream would have, not merely "a gap of the same length". The
// every-base machine checks `above` before it checks gap length (step_machine's
// CLOSING branch), so a base that lands exactly on gap_from + max_gap and is
// itself above the cutoff RESUMES the peak instead of closing it -- and a
// driver that skipped straight to that base must resume it too, not flush
// early because the raw coordinate arithmetic alone reached max_gap.
//
// Both sub-cases of step_machine's gap preamble are exercised: `gap_from` is
// last_pos_ + 1 when the gap is opened by the jump itself (the prior real base
// was IN_PEAK/SUMMIT_SEARCH), and gap_start_ when a real sub-threshold base was
// already observed before the jump (the prior real base was already CLOSING).
void run_every_base_vs_skipped(PeakCallOptions opts, Coord run_len,
                                bool explicit_low_base_first, std::int64_t hi,
                                std::int64_t lo) {
    // Reference: the ground truth, one call per base, nothing skipped.
    PeakCaller ref(opts);
    ref.begin_contig(0);
    for (Coord p = 0; p < run_len; ++p) ref.observe(p, hi, 1.0);
    const Coord gap_begin = run_len;  // first base below the cutoff
    // last_pos_ + 1 (fresh transition) or gap_start_ (already closing) is
    // gap_begin either way -- the every-base walk cannot tell them apart, only
    // the skipping driver's bookkeeping can, which is exactly what is under
    // test.
    // The base that sets gap_start_ (the fresh-transition base at gap_begin
    // itself) is never flush-checked -- only the CLOSING branch does that, and
    // it is reached starting from the base AFTER. So the earliest coordinate
    // the every-base stream could flush at is gap_begin + max(1, max_gap), not
    // gap_begin + max_gap: at max_gap == 0 those differ, and gap_begin +
    // max_gap (== gap_begin + 0 == gap_begin) is not even reachable through a
    // jump, since the caller below always lands at gap_begin or later.
    const Coord gap_len = static_cast<Coord>(opts.max_gap);
    const Coord boundary = gap_begin + std::max<Coord>(1, gap_len);
    for (Coord p = gap_begin; p < boundary; ++p) ref.observe(p, lo, 1.0);
    ref.observe(boundary, hi, 1.0);  // above cutoff, exactly at the boundary
    ref.end_contig();

    // Skipping driver: the same real observations, but the sub-threshold run
    // is never actually walked -- only jumped over, exactly as a driver using
    // MultiScaleWindow::seek() across a coverage desert would.
    PeakCaller skip(opts);
    skip.begin_contig(0);
    for (Coord p = 0; p < run_len; ++p) skip.observe(p, hi, 1.0);
    if (explicit_low_base_first) {
        // One real sub-threshold base lands the machine in CLOSING with
        // gap_start_ == gap_begin before the jump, exercising that branch of
        // the gap_from ternary rather than the "fresh transition" one.
        skip.observe(gap_begin, lo, 1.0);
        skip.observe(boundary, hi, 1.0);
    } else {
        skip.observe(boundary, hi, 1.0);  // straight to the boundary base
    }
    skip.end_contig();

    CHECK_MSG(skip.peaks().size() == ref.peaks().size(),
              "explicit_low_base_first=" + std::to_string(explicit_low_base_first) +
                  ": " + std::to_string(skip.peaks().size()) + " peaks, expected " +
                  std::to_string(ref.peaks().size()));
    if (skip.peaks().size() == ref.peaks().size() && !ref.peaks().empty()) {
        const auto& r = ref.peaks()[0];
        const auto& s = skip.peaks()[0];
        CHECK_MSG(s.start == r.start && s.end == r.end &&
                      s.summit_pos == r.summit_pos && s.summit_score == r.summit_score,
                  "explicit_low_base_first=" + std::to_string(explicit_low_base_first) +
                      ": skip=[" + std::to_string(s.start) + "," + std::to_string(s.end) +
                      ") ref=[" + std::to_string(r.start) + "," + std::to_string(r.end) + ")");
    }
}

void test_a_gap_boundary_base_above_cutoff_resumes_the_skipped_peak() {
    auto opts = permissive();
    opts.max_gap = 30;
    opts.min_length = 1;
    const std::int64_t hi = pileup_scoring(9.0, 1.0);
    const std::int64_t lo = pileup_scoring(0.1, 1.0);
    // Fresh transition: the gap is opened by the jump itself (no real
    // sub-threshold base observed before it).
    run_every_base_vs_skipped(opts, /*run_len=*/10, /*explicit_low_base_first=*/false, hi, lo);
    // Already closing: one real sub-threshold base is observed first, so
    // gap_start_ (not last_pos_ + 1) is the ternary's other branch.
    run_every_base_vs_skipped(opts, /*run_len=*/10, /*explicit_low_base_first=*/true, hi, lo);
}

void test_a_gap_boundary_base_above_cutoff_resumes_the_skipped_peak_max_gap_zero() {
    // The degenerate case: zero tolerance. The every-base machine still
    // reserves one base to notice the drop before it can flush (step_machine's
    // CLOSING branch is reached for the first time only on the base AFTER the
    // one that set gap_start_), so the earliest reachable boundary is
    // gap_from + 1, not gap_from + 0 -- a plain `pos - gap_from >= max_gap`
    // gets this wrong even at max_gap == 0, where naively gap_from + max_gap
    // IS gap_from, a coordinate no jump can land on since the smallest
    // reachable pos is already gap_from + 1.
    //
    // Only the fresh-transition sub-case is meaningful here: "already
    // closing" requires one real low base to be observed at gap_from first,
    // which already consumes the only base max_gap == 0 leaves before the
    // boundary, so there is no room left for an actual jump in that sub-case.
    auto opts = permissive();
    opts.max_gap = 0;
    opts.min_length = 1;
    const std::int64_t hi = pileup_scoring(9.0, 1.0);
    const std::int64_t lo = pileup_scoring(0.1, 1.0);
    run_every_base_vs_skipped(opts, /*run_len=*/10, /*explicit_low_base_first=*/false, hi, lo);
}

void test_a_gap_boundary_base_below_cutoff_still_flushes() {
    // Companion regression pin: when the boundary base is NOT above the
    // cutoff, both sub-cases must still flush exactly where the every-base
    // stream would -- this rules out a fix that swings the other way and
    // never flushes at the boundary at all.
    auto opts = permissive();
    opts.max_gap = 30;
    opts.min_length = 1;
    const std::int64_t hi = pileup_scoring(9.0, 1.0);
    const std::int64_t lo = pileup_scoring(0.1, 1.0);

    for (bool explicit_low_base_first : {false, true}) {
        PeakCaller ref(opts);
        ref.begin_contig(0);
        for (Coord p = 0; p < 10; ++p) ref.observe(p, hi, 1.0);
        for (Coord p = 10; p <= 40; ++p) ref.observe(p, lo, 1.0);  // boundary at 40, below cutoff
        ref.observe(41, hi, 1.0);  // a second, separate peak
        ref.end_contig();

        PeakCaller skip(opts);
        skip.begin_contig(0);
        for (Coord p = 0; p < 10; ++p) skip.observe(p, hi, 1.0);
        if (explicit_low_base_first) {
            skip.observe(10, lo, 1.0);
            skip.observe(40, lo, 1.0);
        } else {
            skip.observe(40, lo, 1.0);
        }
        skip.observe(41, hi, 1.0);
        skip.end_contig();

        CHECK_MSG(skip.peaks().size() == ref.peaks().size() && ref.peaks().size() == 2u,
                  "explicit_low_base_first=" + std::to_string(explicit_low_base_first) +
                      ": skip=" + std::to_string(skip.peaks().size()) +
                      " ref=" + std::to_string(ref.peaks().size()));
    }
}

void test_peaks_shorter_than_min_length_are_discarded() {
    auto opts = permissive();
    opts.min_length = 50;
    PeakCaller c(opts);
    c.begin_contig(0);
    const std::int64_t hi = pileup_scoring(6.0, 1.0);
    Coord p = 0;
    for (int i = 0; i < 10; ++i) c.observe(p++, hi, 1.0);  // 10 bp, too short
    for (int i = 0; i < 40; ++i) c.observe(p++, 1, 1.0);
    const Coord keep_start = p;
    for (int i = 0; i < 80; ++i) c.observe(p++, hi, 1.0);  // 80 bp, kept
    c.end_contig();

    CHECK_EQ(c.peaks().size(), 1u);
    CHECK_EQ(c.discarded_short(), 1);
    if (!c.peaks().empty()) CHECK_EQ(c.peaks()[0].start, keep_start);
}

void test_total_signal_sums_the_pileup_over_the_peak() {
    PeakCaller c(permissive());
    c.begin_contig(0);
    const std::int64_t hi = pileup_scoring(6.0, 1.0);
    Coord p = 0;
    std::int64_t want = 0;
    for (int i = 0; i < 25; ++i) {
        c.observe(p++, hi + i, 1.0);
        want += hi + i;
    }
    c.end_contig();
    CHECK_EQ(c.peaks().size(), 1u);
    if (!c.peaks().empty()) CHECK_EQ(c.peaks()[0].total_signal, want);
}

void test_end_contig_flushes_an_open_peak() {
    // A peak that runs to the last base of a chromosome is still a peak. It is
    // also the one an "emit when the coordinate passes end + gap" rule drops
    // on the floor, because that coordinate never arrives.
    PeakCaller c(permissive());
    c.begin_contig(3);
    const std::int64_t hi = pileup_scoring(6.0, 1.0);
    for (Coord p = 900; p < 1000; ++p) c.observe(p, hi, 1.0);
    CHECK_EQ(c.peaks().size(), 0u);  // still open
    c.end_contig();
    CHECK_EQ(c.peaks().size(), 1u);
    if (!c.peaks().empty()) {
        CHECK_EQ(c.peaks()[0].end, 1000);
        CHECK_EQ(c.peaks()[0].tid, 3);
    }
}

void test_peaks_never_span_contigs() {
    auto opts = permissive();
    opts.max_gap = 1000;
    PeakCaller c(opts);
    const std::int64_t hi = pileup_scoring(6.0, 1.0);
    c.begin_contig(0);
    for (Coord p = 0; p < 20; ++p) c.observe(p, hi, 1.0);
    c.begin_contig(1);  // a new contig, without an explicit end_contig
    for (Coord p = 0; p < 20; ++p) c.observe(p, hi, 1.0);
    c.end_contig();

    CHECK_EQ(c.peaks().size(), 2u);
    if (c.peaks().size() == 2) {
        CHECK_EQ(c.peaks()[0].tid, 0);
        CHECK_EQ(c.peaks()[1].tid, 1);
        CHECK(c.peaks()[1].start == 0);
    }
}

void test_options_reject_nonsense() {
    PeakCallOptions o;
    o.max_gap = -1;
    CHECK_THROWS(std::invalid_argument, o.validate());
    o = PeakCallOptions{};
    o.min_length = -5;
    CHECK_THROWS(std::invalid_argument, o.validate());
    o = PeakCallOptions{};
    o.qvalue_cutoff = 1.5;
    CHECK_THROWS(std::invalid_argument, o.validate());
    o = PeakCallOptions{};
    o.qvalue_cutoff = 0.0;
    CHECK_THROWS(std::invalid_argument, o.validate());
    o = PeakCallOptions{};
    o.score_cutoff = -1.0;
    CHECK_THROWS(std::invalid_argument, o.validate());
    CHECK_THROWS(std::invalid_argument, [] {
        PeakCallOptions bad;
        bad.score_cutoff = std::numeric_limits<double>::quiet_NaN();
        bad.validate();
    }());
}

// --- Benjamini-Hochberg ----------------------------------------------------

// The textbook step-up, in linear probability space, so the header is checked
// against the procedure rather than against a rearrangement of itself.
std::vector<double> reference_bh(std::vector<double> p) {
    const std::size_t m = p.size();
    std::vector<std::size_t> order(m);
    for (std::size_t i = 0; i < m; ++i) order[i] = i;
    std::sort(order.begin(), order.end(),
              [&](std::size_t a, std::size_t b) { return p[a] < p[b]; });
    std::vector<double> q(m, 1.0);
    double running = 1.0;
    for (std::size_t rank = m; rank-- > 0;) {
        const std::size_t i = order[rank];
        const double adj = p[i] * static_cast<double>(m) / static_cast<double>(rank + 1);
        running = std::min(running, std::min(1.0, adj));
        q[i] = running;
    }
    return q;
}

void test_bh_matches_the_textbook_procedure() {
    std::mt19937_64 rng(20260922);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    for (int trial = 0; trial < 40; ++trial) {
        const std::size_t m = 1 + static_cast<std::size_t>(rng() % 200);
        std::vector<double> p(m);
        std::vector<float> scores(m);
        for (std::size_t i = 0; i < m; ++i) {
            p[i] = std::pow(10.0, -u(rng) * 12.0);  // p in [1e-12, 1]
            scores[i] = static_cast<float>(-std::log10(p[i]));
        }
        const std::vector<double> want = reference_bh(p);
        std::vector<float> got;
        pto::peaks::benjamini_hochberg(scores, got);
        CHECK_EQ(got.size(), m);
        for (std::size_t i = 0; i < m && i < got.size(); ++i) {
            const double got_q = std::pow(10.0, -static_cast<double>(got[i]));
            // q comes back as -log10 q in a float, and that storage is what
            // limits agreement with a double reference: one float ulp of the
            // logarithm is a relative error of ln(10) x ulp in q itself. At
            // q = 1e-11 that is about 2e-6, so a flat 1e-6 would be demanding
            // more than the type holds. The bound below is that error model
            // with four ulp of headroom, not a round number.
            const double logq = -std::log10(std::max(want[i], 1e-300));
            const double tol = want[i] * 2.302585092994046 *
                               std::max(1.0, logq) * 4.0 *
                               static_cast<double>(std::numeric_limits<float>::epsilon());
            CHECK_MSG(std::fabs(got_q - want[i]) <= tol,
                      "i=" + std::to_string(i) + " got q " + testing::describe(got_q) +
                          " want " + testing::describe(want[i]));
        }
    }
}

void test_bh_is_monotone_even_for_extreme_scores() {
    // Real summits reach -log10 p of several thousand, which is a p-value of
    // zero in double. The adjustment must stay ordered anyway: a q-value that
    // inverts the p-value ranking would reorder the output track.
    std::vector<float> scores{5000.0f, 4000.0f, 1e-9f, 0.0f, 300.0f, 12.5f, 3000.0f};
    std::vector<float> q;
    pto::peaks::benjamini_hochberg(scores, q);
    CHECK_EQ(q.size(), scores.size());
    for (std::size_t i = 0; i < q.size(); ++i) {
        CHECK(std::isfinite(q[i]));
        CHECK(q[i] >= 0.0f);
        CHECK_MSG(q[i] <= scores[i] + 1e-3f,
                  "q must never be more significant than p at i=" + std::to_string(i));
    }
    // Ordering is preserved: a higher score never gets a lower q-score.
    for (std::size_t i = 0; i < q.size(); ++i) {
        for (std::size_t j = 0; j < q.size(); ++j) {
            if (scores[i] > scores[j]) {
                CHECK_MSG(q[i] >= q[j] - 1e-3f,
                          "BH inverted the ranking at " + std::to_string(i) + "," +
                              std::to_string(j));
            }
        }
    }
}

void test_bh_degenerate_sizes() {
    std::vector<float> q;
    pto::peaks::benjamini_hochberg({}, q);
    CHECK_EQ(q.size(), 0u);
    // With one test, BH is the identity.
    std::vector<float> one{7.0f};
    pto::peaks::benjamini_hochberg(one, q);
    CHECK_EQ(q.size(), 1u);
    if (!q.empty()) CHECK_NEAR(static_cast<double>(q[0]), 7.0, 1e-5);
}

void test_finalize_filters_at_the_q_cutoff() {
    auto opts = permissive();
    opts.min_length = 1;
    opts.max_gap = 5;
    // Stricter than the score cutoff implies, and that is not incidental.
    // Every candidate already has p <= 10^-score_cutoff, and the q-value of
    // the WEAKEST candidate is exactly its own p-value -- the step-up's
    // running minimum reaches it at rank m, where the m/j factor is one. So an
    // FDR level at or above 10^-score_cutoff can never remove anything, and a
    // test that used one would be asserting that filtering happened while
    // measuring a procedure that cannot filter.
    opts.qvalue_cutoff = 1e-4;
    PeakCaller c(opts);
    c.begin_contig(0);

    // One unambiguous peak and a run of marginal ones. After correcting for
    // the number of tests, only the strong one should survive.
    const std::int64_t strong = pileup_scoring(40.0, 1.0);
    const std::int64_t weak = pileup_scoring(2.1, 1.0);
    Coord p = 0;
    for (int i = 0; i < 5; ++i) c.observe(p++, strong, 1.0);
    for (int i = 0; i < 10; ++i) c.observe(p++, 1, 1.0);
    for (int rep = 0; rep < 60; ++rep) {
        for (int i = 0; i < 3; ++i) c.observe(p++, weak, 1.0);
        for (int i = 0; i < 10; ++i) c.observe(p++, 1, 1.0);
    }
    c.end_contig();
    const std::size_t before = c.peaks().size();
    CHECK(before > 50u);
    c.finalize();
    CHECK(c.peaks().size() < before);
    CHECK(c.peaks().size() >= 1u);
    CHECK_EQ(c.qscores().size(), c.peaks().size());
    if (!c.peaks().empty()) {
        CHECK(c.peaks()[0].summit_score > 30.0f);
        for (float q : c.qscores()) CHECK(q >= 4.0f);  // -log10(1e-4)
    }
}

// ===========================================================================
// SIMD batching
//
// The state machine itself cannot be vectorised: the state at base i+1 depends
// on the state at base i, and no amount of masking changes that. What CAN be
// vectorised is the expensive part -- the Poisson survival evaluation -- and
// the way to do it is to run it AHEAD of the machine, over a fixed queue of
// pending bases, then walk the machine over the precomputed scores in order.
//
// So the invariant these cases pin is that batching changes nothing. The queue
// depth is a performance knob; every one of them must produce the same peaks,
// the same summits and the same q-values.
// ===========================================================================

void test_a_coordinate_at_the_top_of_the_range_does_not_overflow() {
    // `pos + 1` in a Coord is an int32 addition, and at the top of the range it
    // is undefined behaviour that happens BEFORE any test of the result can
    // look at it -- the defect sliding_window.hpp's header names, out of
    // docs/REVIEW_2026-08-15 finding (b). MultiScaleWindow clamps position()
    // to contig_len - 1 so the CLI cannot reach it, but observe() is public and
    // this header is installed. Run this under UBSan: without a sanitizer a
    // pass proves only that the machine wrapped quietly.
    PeakCallOptions opts;
    opts.score_cutoff = 0.0;   // every base clears it
    opts.min_length = 0;
    opts.score_batch = 1;
    PeakCaller caller(opts);
    caller.begin_contig(0);
    caller.observe(std::numeric_limits<Coord>::max(), 5, 0.01);
    caller.end_contig();
    for (const CandidatePeak& p : caller.peaks()) {
        CHECK(p.end >= p.start);
        CHECK(p.summit_pos >= p.start);
    }
}

void test_a_position_that_goes_backwards_cannot_invert_a_peak() {
    // observe() requires non-decreasing positions and every driver in this
    // repository supplies them. The clamp is here because the failure mode of
    // trusting it is invisible: extend() would write an `end` BELOW `start`,
    // flush()'s length test would read the negative width as too short, and a
    // real peak would vanish into discarded_short() with nothing reported.
    PeakCallOptions opts;
    opts.score_cutoff = 0.0;
    opts.min_length = 0;
    opts.score_batch = 1;
    PeakCaller caller(opts);
    caller.begin_contig(0);
    for (Coord p : {100, 101, 102, 103, 50}) caller.observe(p, 5, 0.01);
    caller.end_contig();

    CHECK_EQ(caller.peaks().size(), 1u);
    if (!caller.peaks().empty()) {
        const CandidatePeak& p = caller.peaks()[0];
        CHECK_EQ(p.start, 100);
        CHECK_MSG(p.end >= 104, "a backwards position pulled the peak's end below its start");
    }
    CHECK_EQ(caller.discarded_short(), 0);
}

void test_validate_rejects_an_absurd_reservation() {
    // reserve_peaks goes straight to vector::reserve in the constructor, so an
    // absurd value arrives as a length_error out of a constructor rather than
    // as the argument error it is. Every other option in this struct is
    // bounded; this one was not.
    PeakCallOptions opts;
    opts.reserve_peaks = std::numeric_limits<std::size_t>::max();
    bool threw = false;
    try {
        opts.validate();
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    CHECK(threw);
}

void test_the_prefilter_never_skips_a_real_candidate() {
    // Bases that cannot clear the cutoff bypass the queue entirely -- that is
    // what "batch the candidates" means. The skip has to be provable, not
    // heuristic: the median of a Poisson exceeds lambda - ln 2, so a pileup at
    // or below that has P(X >= k) >= 1/2 and a score of at most 0.302.
    // Anything else must go to the queue.
    //
    // A sweep, because a prefilter that is wrong for one (k, lambda) in ten
    // thousand drops one peak in ten thousand and nothing reports it.
    std::mt19937_64 rng(20260923);
    std::uniform_real_distribution<double> lu(0.01, 3000.0);
    long long skipped = 0;
    long long checked = 0;
    for (int i = 0; i < 200000; ++i) {
        const double lam = lu(rng);
        const auto k = static_cast<std::uint32_t>(
            std::min(4e9, lam * std::uniform_real_distribution<double>(0.0, 2.5)(rng)));
        for (double cutoff : {0.5, 1.0, 2.0, 5.0}) {
            ++checked;
            if (!pto::peaks::detail::provably_below_cutoff(k, lam, cutoff)) continue;
            ++skipped;
            const double truth = pto::peaks::neg_log10_poisson_sf(k, lam);
            CHECK_MSG(truth < cutoff,
                      "prefilter skipped k=" + std::to_string(k) + " lam=" +
                          testing::describe(lam) + " cutoff=" + testing::describe(cutoff) +
                          " but its score is " + testing::describe(truth));
        }
    }
    // A prefilter that never fires is safe and useless; say so rather than
    // reporting a pass over nothing.
    CHECK_MSG(skipped > checked / 20,
              "the prefilter fired on only " + std::to_string(skipped) + " of " +
                  std::to_string(checked) + " cases");
}

// Choi's bound (provably_below_cutoff, above) fires below the background
// (k <= lambda - ln2); a background track's modal state is k ~= lambda,
// exactly where it never fires. The skip band (poisson_pmf_upper_bound's
// U(k, lambda) < cutoff, bisected once per k in the PeakCaller constructor)
// is smallest exactly there, so it is meant to cover what Choi's bound
// misses -- checked directly against the ACTUAL band compute_skip_band
// returns, not hand-picked points, so a bisection bug in one direction or a
// shrink applied to the wrong side shows up here as a false skip.
void test_the_skip_band_never_skips_a_real_candidate() {
    std::mt19937_64 rng(20260916);
    std::uniform_int_distribution<std::uint32_t> kd(1, 5000);
    long long skipped = 0, checked = 0;
    for (int i = 0; i < 5000; ++i) {
        const std::uint32_t k = kd(rng);
        for (double cutoff : {0.30103, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0}) {
            ++checked;
            const auto band = pto::peaks::detail::compute_skip_band(k, cutoff);
            if (band.lo >= band.hi) continue;  // empty band: nothing to check
            // Several points across the band's width, including close to its
            // edges, where a rounding mistake in the bisection or the shrink
            // would show first.
            for (double frac : {0.001, 0.05, 0.25, 0.5, 0.75, 0.95, 0.999}) {
                const double lam = band.lo + frac * (band.hi - band.lo);
                ++skipped;
                const double truth = pto::peaks::neg_log10_poisson_sf(k, lam);
                CHECK_MSG(truth < cutoff,
                          "skip band admitted k=" + std::to_string(k) + " lambda=" +
                              testing::describe(lam) + " cutoff=" + testing::describe(cutoff) +
                              " but its score is " + testing::describe(truth));
            }
            // The band's own claimed centre (lambda == k) must be inside it
            // whenever the band is non-empty -- that is the whole point of
            // U's convexity argument (minimised at lambda = k).
            CHECK_MSG(static_cast<double>(k) > band.lo && static_cast<double>(k) < band.hi,
                      "k=" + std::to_string(k) + " cutoff=" + testing::describe(cutoff) +
                          " band=[" + testing::describe(band.lo) + "," +
                          testing::describe(band.hi) + ") excludes its own lambda=k centre");
        }
    }
    // A band that never forms is safe and useless; say so rather than
    // reporting a pass over nothing.
    CHECK_MSG(skipped > checked,
              "the skip band fired on too few of " + std::to_string(checked) + " cases (" +
                  std::to_string(skipped) + " points sampled)");
}

// The two prefilters are complementary by design: Choi covers k well below
// lambda, the skip band covers k near lambda. A background track's bases
// (k ~= lambda) must still be correctly scored -- not silently classified as
// "above cutoff" by a band that reached too far -- and must still fail to
// open a peak, exactly as they did before either prefilter existed.
void test_background_track_bases_are_unaffected_by_either_prefilter() {
    // cutoff=8.0 means a truly Poisson(20)-distributed base clears it with
    // probability ~1e-8; over 200000 iid draws the expected false-positive
    // count is ~0.002, so an empty peak list is the correct null result, not
    // a statistical coin flip. (cutoff=3.0 expects ~200 false positives from
    // pure sampling noise alone -- that was this test's original, wrong,
    // invariant, caught by re-deriving the math rather than by a code bug.)
    auto opts = permissive();
    opts.score_cutoff = 8.0;
    opts.min_length = 1;
    PeakCaller c(opts);
    c.begin_contig(0);
    std::mt19937_64 rng(20260916);
    std::poisson_distribution<int> pois(20.0);
    for (Coord p = 0; p < 200000; ++p) {
        const auto k = static_cast<std::int64_t>(std::max(0, pois(rng)));
        c.observe(p, k, 20.0);  // k centred on lambda: Choi never fires here
    }
    c.end_contig();
    CHECK_MSG(c.peaks().empty(),
              "a background track (k ~ Poisson(20), lambda=20) produced " +
                  std::to_string(c.peaks().size()) + " peaks at cutoff 8.0");
}

// Finding 9 (2026-09-16 audit): observe() now drops a base the instant either
// prefilter proves it, before it is ever queued, rather than queueing it and
// letting drain_queue() score it as 0.0. The argument that this cannot move a
// peak boundary is two already-tested facts composed: the prefilters are
// SOUND (test_the_prefilter_never_skips_a_real_candidate,
// test_the_skip_band_never_skips_a_real_candidate -- a dropped base's real
// score, if computed, would provably be under the cutoff), and step_machine
// already treats an unobserved coordinate range exactly like a run of
// below-cutoff bases (test_a_skipped_coordinate_range_counts_as_a_gap, and
// seek() itself). This test exercises the composition directly rather than
// re-deriving it: two peaks planted in heavy Poisson background at a cutoff
// where most background bases ARE dropped (unlike the test above, which picks
// a cutoff specifically so nothing is), and the caller must still localise
// both at their planted coordinates.
void test_planted_peaks_survive_heavy_prefilter_dropping() {
    auto opts = permissive();
    // As in the background-track test above: cutoff=8.0 against lambda=15
    // makes a background false positive a ~1e-8 event per base, so ~59500
    // background bases here contribute an expected ~0.0006 spurious peaks --
    // the two planted peaks below are the only ones that should ever appear.
    opts.score_cutoff = 8.0;
    opts.min_length = 20;
    opts.max_gap = 10;
    PeakCaller c(opts);
    c.begin_contig(0);

    constexpr Coord kContigLen = 60000;
    constexpr Coord kPeakA_Start = 10000, kPeakA_End = 10200;    // [start, end)
    constexpr Coord kPeakB_Start = 40000, kPeakB_End = 40300;
    const std::int64_t peak_pileup = pileup_scoring(40.0, 15.0);  // clears cutoff=8.0 by a wide margin

    std::mt19937_64 rng(20260916);
    std::poisson_distribution<int> pois(15.0);
    for (Coord p = 0; p < kContigLen; ++p) {
        const bool in_a = p >= kPeakA_Start && p < kPeakA_End;
        const bool in_b = p >= kPeakB_Start && p < kPeakB_End;
        const std::int64_t k =
            (in_a || in_b) ? peak_pileup : static_cast<std::int64_t>(std::max(0, pois(rng)));
        c.observe(p, k, 15.0);
    }
    c.end_contig();

    CHECK_EQ(c.peaks().size(), 2u);
    if (c.peaks().size() == 2) {
        CHECK_EQ(c.peaks()[0].start, kPeakA_Start);
        CHECK_EQ(c.peaks()[0].end, kPeakA_End);
        CHECK_EQ(c.peaks()[1].start, kPeakB_Start);
        CHECK_EQ(c.peaks()[1].end, kPeakB_End);
    }
}

void test_every_queue_depth_gives_the_same_peaks() {
    // Well-separated scores, so the comparison is about the batching and not
    // about a base sitting one ulp from the cutoff.
    std::mt19937_64 rng(20260924);
    std::vector<std::int64_t> pile;
    std::vector<double> lam;
    const std::int64_t hi = pileup_scoring(20.0, 1.0);
    for (int i = 0; i < 4000; ++i) {
        const bool enriched = (i / 37) % 3 == 0;
        pile.push_back(enriched ? hi + static_cast<std::int64_t>(rng() % 5) : 0);
        lam.push_back(1.0);
    }

    std::vector<CandidatePeak> reference;
    std::vector<float> reference_q;
    for (std::size_t depth : {std::size_t{1}, std::size_t{2}, std::size_t{3},
                              std::size_t{8}, std::size_t{64}, std::size_t{256},
                              std::size_t{4096}}) {
        auto opts = permissive();
        opts.min_length = 5;
        opts.max_gap = 10;
        opts.score_batch = depth;
        PeakCaller c(opts);
        c.begin_contig(0);
        for (std::size_t i = 0; i < pile.size(); ++i) {
            c.observe(static_cast<Coord>(i), pile[i], lam[i]);
        }
        c.end_contig();
        c.finalize();

        if (reference.empty()) {
            reference = c.peaks();
            reference_q = c.qscores();
            CHECK(reference.size() > 10u);
            continue;
        }
        CHECK_MSG(c.peaks().size() == reference.size(),
                  "depth " + std::to_string(depth) + ": " +
                      std::to_string(c.peaks().size()) + " peaks vs " +
                      std::to_string(reference.size()));
        for (std::size_t i = 0; i < c.peaks().size() && i < reference.size(); ++i) {
            const CandidatePeak& a = c.peaks()[i];
            const CandidatePeak& b = reference[i];
            const std::string at = "depth " + std::to_string(depth) + " peak " +
                                   std::to_string(i);
            CHECK_MSG(a.start == b.start && a.end == b.end, at + ": bounds");
            CHECK_MSG(a.summit_pos == b.summit_pos, at + ": summit");
            CHECK_MSG(a.total_signal == b.total_signal, at + ": signal");
            CHECK_MSG(a.summit_score == b.summit_score, at + ": score");
        }
        for (std::size_t i = 0; i < c.qscores().size() && i < reference_q.size(); ++i) {
            CHECK_MSG(c.qscores()[i] == reference_q[i],
                      "depth " + std::to_string(depth) + " q " + std::to_string(i));
        }
    }
}

void test_a_partial_queue_is_flushed_at_every_boundary() {
    // A peak that ends mid-queue must still be emitted. The queue is drained
    // at end_contig and at begin_contig, and both are tested because the
    // second is the one a caller forgets.
    for (std::size_t depth : {std::size_t{1}, std::size_t{8}, std::size_t{64}}) {
        auto opts = permissive();
        opts.score_batch = depth;
        opts.min_length = 1;
        PeakCaller c(opts);
        const std::int64_t hi = pileup_scoring(9.0, 1.0);

        c.begin_contig(0);
        for (Coord p = 0; p < 3; ++p) c.observe(p, hi, 1.0);  // fewer bases than a full queue
        c.begin_contig(1);                                     // no end_contig first
        for (Coord p = 0; p < 3; ++p) c.observe(p, hi, 1.0);
        c.end_contig();

        CHECK_MSG(c.peaks().size() == 2u,
                  "depth " + std::to_string(depth) + ": " +
                      std::to_string(c.peaks().size()) + " peaks");
        if (c.peaks().size() == 2) {
            CHECK(c.peaks()[0].tid == 0);
            CHECK(c.peaks()[1].tid == 1);
        }
    }
}

void test_the_queue_preserves_coordinate_order_across_a_jump() {
    // A coordinate jump is a gap, and the gap logic reads positions off the
    // queue. If the queue reordered or coalesced anything, this stitches two
    // peaks into one.
    for (std::size_t depth : {std::size_t{1}, std::size_t{16}, std::size_t{512}}) {
        auto opts = permissive();
        opts.score_batch = depth;
        opts.max_gap = 30;
        opts.min_length = 1;
        PeakCaller c(opts);
        c.begin_contig(0);
        const std::int64_t hi = pileup_scoring(9.0, 1.0);
        for (Coord p = 100; p < 110; ++p) c.observe(p, hi, 1.0);
        c.observe(500000, hi, 1.0);
        c.end_contig();
        CHECK_MSG(c.peaks().size() == 2u,
                  "depth " + std::to_string(depth) + ": " +
                      std::to_string(c.peaks().size()) + " peaks");
    }
}

void test_observe_allocates_nothing_once_reserved() {
    auto opts = permissive();
    opts.reserve_peaks = 4096;
    PeakCaller c(opts);
    c.begin_contig(0);
    const std::int64_t hi = pileup_scoring(6.0, 1.0);
    {
        testing::AllocationGuard guard;
        Coord p = 0;
        for (int rep = 0; rep < 500; ++rep) {
            for (int i = 0; i < 5; ++i) c.observe(p++, hi, 1.0);
            for (int i = 0; i < 30; ++i) c.observe(p++, 1, 1.0);
        }
        c.end_contig();
        // Read the counter BEFORE building any message: std::to_string
        // allocates, and the order in which a function's arguments are
        // evaluated is unspecified, so on libstdc++ the message's own
        // allocation was being counted against the loop.
        const long long allocs = guard.count();
        if (testing::AllocationGuard::counting()) {
            CHECK_MSG(allocs == 0,
                      std::to_string(allocs) + " allocation(s) in the observe loop");
        } else {
            testing::skip("allocation counting is off under sanitizers");
        }
    }
    CHECK(c.peaks().size() >= 400u);
}

void test_total_signal_saturates_rather_than_wrapping() {
    PeakCaller c(permissive());
    c.begin_contig(0);
    const std::int64_t huge = std::numeric_limits<std::int64_t>::max() / 4;
    for (Coord p = 0; p < 10; ++p) c.observe(p, huge, 1.0);
    c.end_contig();
    CHECK_EQ(c.peaks().size(), 1u);
    if (!c.peaks().empty()) {
        CHECK(c.peaks()[0].total_signal > 0);
        CHECK(c.saturated());
    }
}

void test_end_to_end_window_and_caller_localise_a_planted_peak() {
    // The only case that runs the two halves against each other: a flat
    // background with one dense pile, driven through the real window.
    constexpr Coord kLen = 60000;
    constexpr Coord kPeakAt = 30000;
    pto::peaks::MultiScaleWindow w;
    w.reset(0, kLen);

    std::vector<std::pair<Coord, Coord>> frags;
    for (Coord s = 0; s + 200 < kLen; s += 200) frags.emplace_back(s, s + 200);  // background
    for (int i = 0; i < 400; ++i) frags.emplace_back(kPeakAt, kPeakAt + 200);    // the peak
    std::sort(frags.begin(), frags.end());

    auto opts = permissive();
    opts.min_length = 20;
    opts.max_gap = 50;
    opts.score_cutoff = 5.0;
    opts.qvalue_cutoff = 0.05;
    PeakCaller c(opts);
    c.begin_contig(0);

    // The real streaming loop: add, drain, repeat. The ring holds a bounded
    // lookahead, so the two halves have to advance together.
    auto drain = [&] {
        while (w.step()) c.observe(w.position(), w.pileup(), w.local_lambda(0.5));
    };
    for (auto [fs, fe] : frags) {
        auto r = w.add_fragment(fs, fe);
        while (r == pto::peaks::AddResult::kNeedsDrain) {
            if (!w.step()) break;
            c.observe(w.position(), w.pileup(), w.local_lambda(0.5));
            r = w.add_fragment(fs, fe);
        }
        CHECK(r == pto::peaks::AddResult::kAccepted);
        drain();
    }
    w.finish();
    drain();
    c.end_contig();
    c.finalize();

    CHECK_MSG(c.peaks().size() >= 1u,
              "found " + std::to_string(c.peaks().size()) + " peaks");
    if (!c.peaks().empty()) {
        const CandidatePeak& pk = c.peaks()[0];
        CHECK_MSG(pk.summit_pos >= kPeakAt && pk.summit_pos < kPeakAt + 200,
                  "summit at " + std::to_string(pk.summit_pos) + ", planted at " +
                      std::to_string(kPeakAt));
        CHECK(pk.start >= kPeakAt - 500);
        CHECK(pk.end <= kPeakAt + 700);
    }
    // The flat background must not produce a second call.
    CHECK_MSG(c.peaks().size() <= 2u,
              "flat background produced " + std::to_string(c.peaks().size()) + " peaks");
}

// ===========================================================================
// Wiring seek() into the driver: a fast-forward across a coverage desert must
// never change which peaks are called, only how many bases it takes to get
// there. main.cpp's Driver uses the SAME choose_seek_target() predicate this
// test does (peaks/sliding_window.hpp) -- one decision, not a restatement of
// it here and another one there.
// ===========================================================================

// Mimics main.cpp's Driver::add()/end_contig(), with or without the seek
// fast-forward. `use_seek=false` is exactly
// test_end_to_end_window_and_caller_localise_a_planted_peak's own streaming
// loop.
//
// `use_seek=true` reproduces the driver's actual protocol. seek()'s
// unconditional ring wipe makes safety, not just eligibility, the deciding
// question -- see choose_seek_target's own comment for why that means only a
// LEADING gap (nothing added to the contig yet) or a TRAILING one (handled
// once, at end of contig) are ever skipped; an interior gap between two real,
// narrow peaks is not, because bounded-draining the first peak's own region
// needs frontier that only a fragment far enough past it can provide, and by
// the time such a fragment exists it has already been added -- seeking to a
// point BEFORE it at that stage would erase ITS delta events too.
//
// 1. On a fragment arrival, BEFORE adding it: if nothing has been added to
//    this contig yet (signal_limit() < 0), there is nothing to lose, so a
//    leading seek is attempted directly. Otherwise, bounded-drain up to
//    signal_limit() using whatever frontier earlier fragments already
//    provided; only if that FULLY succeeds (position reaches it) is seeking
//    attempted -- if it does not, no not-yet-observed signal has been
//    discarded, and the fragment is simply added and drained normally, same
//    as use_seek=false.
// 2. At end of contig: finish() unlocks evaluable_limit() to the contig's end
//    regardless of the frontier, which is what lets a bounded drain reach the
//    last fragment's own signal_limit() at all (the frontier alone is stuck
//    at that fragment's start, not its end). seek() cannot be called once
//    finished_ is set, but nothing needs it to be: everything past
//    signal_limit() is provably k=0 background, and PeakCaller::flush()
//    finalises an open peak from its own already-stored fields regardless of
//    how much of that background is ever stepped through -- so the bounded
//    drain simply stops, and end_contig() is called directly. Still gated by
//    choose_seek_target (score_cutoff > 0, a real gap, past the current
//    position): at cutoff == 0 every base is "above" regardless of pileup, so
//    the open peak keeps growing through the background too, and stopping
//    early would wrongly cut it short.
std::vector<CandidatePeak> run_driver(const PeakCallOptions& opts, Coord contig_len,
                                      const std::vector<std::pair<Coord, Coord>>& frags,
                                      double lambda_bg, bool use_seek) {
    pto::peaks::MultiScaleWindow w;
    w.reset(0, contig_len);
    PeakCaller c(opts);
    c.begin_contig(0);

    auto visit = [&] { c.observe(w.position(), w.pileup(), w.local_lambda(lambda_bg)); };
    auto drain = [&] {
        while (w.step()) visit();
    };
    auto bounded_drain = [&](std::int64_t up_to) {
        while (w.position() < up_to && w.step()) visit();
    };

    for (auto [fs, fe] : frags) {
        if (use_seek) {
            const std::int64_t live = w.signal_limit();  // before THIS fragment
            bool safe_to_seek;
            if (live < 0) {
                safe_to_seek = true;  // nothing added yet: nothing to lose
            } else {
                bounded_drain(live);
                safe_to_seek = (w.position() >= live);
            }
            if (safe_to_seek) {
                const std::int64_t candidate =
                    static_cast<std::int64_t>(fs) - pto::peaks::kMaxHalfWidth;
                const std::int64_t t = pto::peaks::choose_seek_target(
                    live, candidate, w.contig_length(), w.position(), opts.score_cutoff);
                if (t >= 0) CHECK(w.seek(static_cast<Coord>(t)));
            }
        }

        auto r = w.add_fragment(fs, fe);
        while (r == pto::peaks::AddResult::kNeedsDrain) {
            if (!w.step()) break;
            visit();
            r = w.add_fragment(fs, fe);
        }
        CHECK(r == pto::peaks::AddResult::kAccepted);
        drain();
    }

    const std::int64_t final_live = w.signal_limit();
    w.finish();
    if (use_seek) {
        const std::int64_t stop_at = pto::peaks::choose_seek_target(
            final_live, w.contig_length() - 1, w.contig_length(), w.position(), opts.score_cutoff);
        if (stop_at >= 0) {
            bounded_drain(final_live);  // the last peak's own region only
        } else {
            drain();  // no exploitable gap (cutoff <= 0 included): walk normally
        }
    } else {
        drain();
    }
    c.end_contig();
    c.finalize();
    return c.peaks();
}

void check_seek_matches_no_seek(const char* label, const PeakCallOptions& opts, Coord contig_len,
                                const std::vector<std::pair<Coord, Coord>>& frags,
                                double lambda_bg = 0.5) {
    const auto ref = run_driver(opts, contig_len, frags, lambda_bg, /*use_seek=*/false);
    const auto got = run_driver(opts, contig_len, frags, lambda_bg, /*use_seek=*/true);
    CHECK_MSG(got.size() == ref.size(),
              std::string(label) + ": " + std::to_string(got.size()) + " peaks, expected " +
                  std::to_string(ref.size()));
    for (std::size_t i = 0; i < got.size() && i < ref.size(); ++i) {
        CHECK_MSG(got[i].start == ref[i].start && got[i].end == ref[i].end &&
                      got[i].summit_pos == ref[i].summit_pos &&
                      got[i].summit_score == ref[i].summit_score &&
                      got[i].total_signal == ref[i].total_signal,
                  std::string(label) + " peak " + std::to_string(i) + ": got=[" +
                      std::to_string(got[i].start) + "," + std::to_string(got[i].end) +
                      ") ref=[" + std::to_string(ref[i].start) + "," +
                      std::to_string(ref[i].end) + ")");
    }
}

// One dense pile (a real peak), optionally a second, with a huge gap of no
// coverage at all in between and/or after -- the shape every case below
// reuses.
std::vector<std::pair<Coord, Coord>> peak_then_desert(Coord peak_at, Coord second_peak_at = -1) {
    std::vector<std::pair<Coord, Coord>> frags;
    for (int i = 0; i < 100; ++i) frags.emplace_back(peak_at, peak_at + 200);
    if (second_peak_at >= 0) {
        for (int i = 0; i < 100; ++i) frags.emplace_back(second_peak_at, second_peak_at + 200);
    }
    return frags;
}

void test_seek_never_changes_the_peak_list() {
    auto opts = permissive();
    opts.min_length = 20;
    opts.score_cutoff = 5.0;

    // A real desert: one peak, then nothing, then the contig ends far away.
    check_seek_matches_no_seek("desert then end of contig", opts, 300000,
                               peak_then_desert(1000), 0.5);

    // A desert between two real peaks.
    check_seek_matches_no_seek("desert between two peaks", opts, 300000,
                               peak_then_desert(1000, 250000), 0.5);

    // Large max_gap: far bigger than any true desert here would need, but the
    // caller's own gap accounting must still agree exactly either way.
    auto big_gap = opts;
    big_gap.max_gap = 100000;
    check_seek_matches_no_seek("large max_gap across a desert", big_gap, 300000,
                               peak_then_desert(1000, 250000), 0.5);

    // Empty contig tail: nothing at all after the one peak, out to a contig
    // long enough that an every-base walk of the tail alone would be the
    // dominant cost -- exactly the case seek() exists for.
    check_seek_matches_no_seek("empty contig tail, very long contig", opts, 5000000,
                               peak_then_desert(1000), 0.5);

    // cutoff == 0: a k=0 base's own score is exactly 0, which meets "at or
    // above the cutoff" -- see detail::provably_below_cutoff -- so
    // choose_seek_target must refuse unconditionally and seek must never be
    // taken, checked both directly and through the same end-to-end
    // comparison as every other case above.
    auto zero_cutoff = opts;
    zero_cutoff.score_cutoff = 0.0;
    zero_cutoff.qvalue_cutoff = 1.0;
    check_seek_matches_no_seek("cutoff == 0", zero_cutoff, 300000,
                               peak_then_desert(1000, 250000), 0.5);
    {
        pto::peaks::MultiScaleWindow w;
        w.reset(0, 300000);
        w.add_fragment(1000, 1200);
        // A target that WOULD be eligible at any positive cutoff (a huge gap
        // past the fragment's own end, well inside the contig, well past the
        // window's own current position) must still be refused at cutoff 0.
        CHECK_EQ(pto::peaks::choose_seek_target(w.signal_limit(), 200000, w.contig_length(),
                                                w.position(), 0.0),
                 std::int64_t(-1));
    }
}

// ===========================================================================
// --extend-peaks: the dual-mode boundary architecture.
//
// The default is a sharp, high-resolution call -- the enrichment core and
// nothing more, which is what precise transcription-factor footprinting wants.
// `extend_peaks > 0` pads the REPORTED interval outward, symmetrically, for
// pipelines that expect MACS-style spans. It touches the bounds and nothing
// else: the summit, its score, its background and the summed signal are the
// measurement and are left exactly as called. It is applied AFTER the
// min_length gate, so a short core is never padded into existence.
// ===========================================================================

// Drives the machine over one flat plateau [first, last] and returns its peak.
CandidatePeak call_one_plateau(const PeakCallOptions& opts, Coord first, Coord last) {
    PeakCaller c(opts);
    c.begin_contig(0);
    const std::int64_t flat = pileup_scoring(9.0, 1.0);
    for (Coord p = first; p <= last; ++p) c.observe(p, flat, 1.0);
    c.end_contig();
    CHECK_EQ(c.peaks().size(), 1u);
    return c.peaks().empty() ? CandidatePeak{} : c.peaks()[0];
}

void test_extend_peaks_defaults_to_sharp_bounds() {
    auto opts = permissive();
    CHECK_EQ(opts.extend_peaks, 0);  // Option A is the default
    const CandidatePeak pk = call_one_plateau(opts, 1000, 1099);
    // The enrichment core is [1000, 1100): 100 bp, and that is what comes back.
    CHECK_EQ(pk.start, 1000);
    CHECK_EQ(pk.end, 1100);
}

void test_extend_peaks_widens_the_interval_symmetrically() {
    auto opts = permissive();
    opts.extend_peaks = 250;
    const CandidatePeak pk = call_one_plateau(opts, 1000, 1099);
    CHECK_EQ(pk.start, 750);   // 1000 - 250
    CHECK_EQ(pk.end, 1350);    // 1100 + 250
}

void test_extend_peaks_preserves_the_summit_and_every_score() {
    const Coord first = 4000;
    const Coord last = 4199;  // a 200 bp plateau
    auto sharp = permissive();
    auto padded = permissive();
    padded.extend_peaks = 500;
    const CandidatePeak a = call_one_plateau(sharp, first, last);
    const CandidatePeak b = call_one_plateau(padded, first, last);

    CHECK_EQ(a.summit_pos, b.summit_pos);      // the absolute summit does not move
    CHECK_EQ(a.summit_score, b.summit_score);
    CHECK_EQ(a.summit_lambda, b.summit_lambda);
    CHECK_EQ(a.total_signal, b.total_signal);  // padding adds no signal
    CHECK_EQ(a.tid, b.tid);
    // Only the bounds moved, and by exactly the extension on each side.
    CHECK_EQ(b.start, a.start - 500);
    CHECK_EQ(b.end, a.end + 500);
    // The summit is still an in-bounds position of the widened peak.
    CHECK(b.summit_pos >= b.start);
    CHECK(b.summit_pos < b.end);
}

void test_extend_peaks_clamps_the_start_at_zero() {
    auto opts = permissive();
    opts.extend_peaks = 5000;
    // A plateau near the origin: 5000 bp of padding would take the start
    // negative, and Coord is signed, so this is the difference between a valid
    // BED and a file every downstream tool rejects. The 3' side is left for the
    // writer to clamp against the contig length, which is the layer that has it.
    const CandidatePeak pk = call_one_plateau(opts, 40, 139);
    CHECK_EQ(pk.start, 0);
    CHECK_EQ(pk.end, 140 + 5000);
    CHECK(pk.summit_pos >= 0);
    CHECK(pk.summit_pos < pk.end);
}

void test_extend_peaks_does_not_rescue_a_call_shorter_than_min_length() {
    auto opts = permissive();
    opts.min_length = 200;
    opts.extend_peaks = 1000;
    PeakCaller c(opts);
    c.begin_contig(0);
    const std::int64_t flat = pileup_scoring(9.0, 1.0);
    for (Coord p = 1000; p < 1030; ++p) c.observe(p, flat, 1.0);  // 30 bp core
    c.end_contig();
    // The core is 30 bp and the gate is 200; padding is applied after the gate,
    // so this stays discarded rather than becoming a 2030 bp "peak".
    CHECK_EQ(c.peaks().size(), 0u);
    CHECK_EQ(c.discarded_short(), 1);
}

void test_extend_peaks_is_deterministic_across_queue_depths() {
    // Padding is a per-peak transform, so it must not interact with the score
    // batching -- the same invariant test_every_queue_depth_gives_the_same_peaks
    // pins for the sharp path.
    std::vector<CandidatePeak> reference;
    for (std::size_t depth : {std::size_t{1}, std::size_t{7}, std::size_t{256}}) {
        auto opts = permissive();
        opts.min_length = 5;
        opts.extend_peaks = 120;
        opts.score_batch = depth;
        PeakCaller c(opts);
        c.begin_contig(0);
        const std::int64_t hi = pileup_scoring(9.0, 1.0);
        Coord p = 500;
        for (int rep = 0; rep < 6; ++rep) {
            for (int i = 0; i < 30; ++i) c.observe(p++, hi, 1.0);
            for (int i = 0; i < 60; ++i) c.observe(p++, 1, 1.0);
        }
        c.end_contig();
        if (reference.empty()) {
            reference = c.peaks();
            CHECK(reference.size() >= 3u);
            continue;
        }
        CHECK_EQ(c.peaks().size(), reference.size());
        for (std::size_t i = 0; i < c.peaks().size() && i < reference.size(); ++i) {
            CHECK_MSG(c.peaks()[i].start == reference[i].start &&
                          c.peaks()[i].end == reference[i].end,
                      "depth " + std::to_string(depth) + " peak " + std::to_string(i));
        }
    }
}

void test_options_reject_a_negative_extension() {
    PeakCallOptions o;
    o.extend_peaks = -1;
    CHECK_THROWS(std::invalid_argument, o.validate());
}

}  // namespace

int main() {
    testing::Suite suite{"peaks/caller", {}};
    suite.add("CandidatePeak is 32 bytes", test_candidate_peak_is_thirty_two_bytes);
    suite.add("idle stays idle below the cutoff", test_idle_stays_idle_below_the_cutoff);
    suite.add("the state machine walks all four states", test_the_state_machine_walks_all_four_states);
    suite.add("summit is the highest-scoring base", test_summit_is_the_highest_scoring_base);
    suite.add("a later higher summit moves the summit", test_a_later_higher_summit_moves_the_summit);
    suite.add("a flat plateau puts the summit at its middle", test_a_flat_plateau_puts_the_summit_at_its_middle);
    suite.add("a short gap does not split a peak", test_a_short_gap_does_not_split_a_peak);
    suite.add("a long gap splits into two peaks", test_a_long_gap_splits_into_two_peaks);
    suite.add("a skipped coordinate range counts as a gap", test_a_skipped_coordinate_range_counts_as_a_gap);
    suite.add("a gap boundary base above cutoff resumes the skipped peak", test_a_gap_boundary_base_above_cutoff_resumes_the_skipped_peak);
    suite.add("a gap boundary base above cutoff resumes the skipped peak, max_gap=0", test_a_gap_boundary_base_above_cutoff_resumes_the_skipped_peak_max_gap_zero);
    suite.add("a gap boundary base below cutoff still flushes", test_a_gap_boundary_base_below_cutoff_still_flushes);
    suite.add("peaks shorter than min_length are discarded", test_peaks_shorter_than_min_length_are_discarded);
    suite.add("total_signal sums the pileup", test_total_signal_sums_the_pileup_over_the_peak);
    suite.add("end_contig flushes an open peak", test_end_contig_flushes_an_open_peak);
    suite.add("peaks never span contigs", test_peaks_never_span_contigs);
    suite.add("options reject nonsense", test_options_reject_nonsense);
    suite.add("BH matches the textbook procedure", test_bh_matches_the_textbook_procedure);
    suite.add("BH is monotone for extreme scores", test_bh_is_monotone_even_for_extreme_scores);
    suite.add("BH degenerate sizes", test_bh_degenerate_sizes);
    suite.add("finalize filters at the q cutoff", test_finalize_filters_at_the_q_cutoff);
    suite.add("a coordinate at the top of the range does not overflow", test_a_coordinate_at_the_top_of_the_range_does_not_overflow);
    suite.add("a backwards position cannot invert a peak", test_a_position_that_goes_backwards_cannot_invert_a_peak);
    suite.add("validate rejects an absurd reservation", test_validate_rejects_an_absurd_reservation);
    suite.add("the prefilter never skips a real candidate", test_the_prefilter_never_skips_a_real_candidate);
    suite.add("the skip band never skips a real candidate", test_the_skip_band_never_skips_a_real_candidate);
    suite.add("background track bases are unaffected by either prefilter", test_background_track_bases_are_unaffected_by_either_prefilter);
    suite.add("planted peaks survive heavy prefilter dropping", test_planted_peaks_survive_heavy_prefilter_dropping);
    suite.add("every queue depth gives the same peaks", test_every_queue_depth_gives_the_same_peaks);
    suite.add("a partial queue is flushed at every boundary", test_a_partial_queue_is_flushed_at_every_boundary);
    suite.add("the queue preserves coordinate order", test_the_queue_preserves_coordinate_order_across_a_jump);
    suite.add("observe allocates nothing once reserved", test_observe_allocates_nothing_once_reserved);
    suite.add("total_signal saturates rather than wrapping", test_total_signal_saturates_rather_than_wrapping);
    suite.add("end to end: window and caller localise a peak", test_end_to_end_window_and_caller_localise_a_planted_peak);
    suite.add("seek never changes the peak list", test_seek_never_changes_the_peak_list);
    suite.add("extend_peaks defaults to sharp bounds", test_extend_peaks_defaults_to_sharp_bounds);
    suite.add("extend_peaks widens the interval symmetrically", test_extend_peaks_widens_the_interval_symmetrically);
    suite.add("extend_peaks preserves the summit and every score", test_extend_peaks_preserves_the_summit_and_every_score);
    suite.add("extend_peaks clamps the start at zero", test_extend_peaks_clamps_the_start_at_zero);
    suite.add("extend_peaks does not rescue a short call", test_extend_peaks_does_not_rescue_a_call_shorter_than_min_length);
    suite.add("extend_peaks is deterministic across queue depths", test_extend_peaks_is_deterministic_across_queue_depths);
    suite.add("options reject a negative extension", test_options_reject_a_negative_extension);
    return suite.run();
}
