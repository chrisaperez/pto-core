// SPDX-License-Identifier: MIT
//
// pto-peaks: the multi-scale sliding background window.
//
// This is the structure MACS2 spends most of its runtime avoiding. It keeps a
// per-base pileup and three centred background windows (1kb, 5kb, 10kb) over a
// stream of coordinate-sorted fragments, in O(1) per base and in memory that
// does not grow with the chromosome.
//
// The properties worth testing here are the ones an eyeball cannot check on a
// real chromosome:
//
//   * the rolling sums equal a direct recomputation, at every base, for a
//     random pileup -- an O(1) update that drifts is the whole failure mode;
//   * the windows are CENTRED, because a trailing window would include a
//     peak's own rising edge in the background it is judged against, and
//     systematically truncate its 3' end;
//   * they truncate at contig boundaries, so a peak near a telomere is not
//     divided by 10kb of chromosome that does not exist;
//   * an out-of-order fragment is refused rather than silently mis-binned;
//   * coordinate arithmetic near INT32_MAX does not overflow. The window adds
//     half-widths to positions, and `pos + 5000` on an int32 near the ceiling
//     is undefined behaviour BEFORE any range check can look at it. That is
//     the exact defect REVIEW_2026-08-15 finding (b) records in
//     cuttag_profiler's window_length(), reappearing in a new module;
//   * and that the streaming loop allocates nothing at all.
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <cstdint>
#include <limits>
#include <random>
#include <vector>

#include "peaks_test_util.hpp"

#include "peaks/sliding_window.hpp"

namespace {

using pto::peaks::AddResult;
using pto::peaks::Coord;
using pto::peaks::MultiScaleWindow;

constexpr Coord kMaxCoord = std::numeric_limits<Coord>::max();

// Drives the window to `target`, returning false if it could not get there.
bool run_to(MultiScaleWindow& w, Coord target) {
    while (w.position() < target) {
        if (!w.step()) return false;
    }
    return w.position() == target;
}

// The real streaming loop, written once. The ring holds a bounded lookahead,
// so a caller that added a whole chromosome before evaluating any of it is
// told to drain; every test that spans more than a few kb goes through here,
// because that is how the CLI drives it too.
template <class Visit>
void drive(MultiScaleWindow& w, const std::vector<std::pair<Coord, Coord>>& frags,
           Visit&& visit) {
    for (const auto& f : frags) {
        AddResult r = w.add_fragment(f.first, f.second);
        while (r == AddResult::kNeedsDrain) {
            if (!w.step()) break;
            visit(w);
            r = w.add_fragment(f.first, f.second);
        }
        CHECK_MSG(r == AddResult::kAccepted,
                  "fragment [" + std::to_string(f.first) + "," +
                      std::to_string(f.second) + ") refused");
        while (w.step()) visit(w);
    }
    w.finish();
    while (w.step()) visit(w);
}

// The definition the window must match, computed the slow way from a dense
// pileup array. Centred, half-open, truncated at both contig ends.
struct Reference {
    std::vector<std::int64_t> pileup;
    Coord length = 0;

    explicit Reference(Coord len) : pileup(static_cast<std::size_t>(len), 0), length(len) {}

    void add(Coord s, Coord e) {
        for (Coord q = std::max<Coord>(0, s); q < std::min(e, length); ++q) {
            ++pileup[static_cast<std::size_t>(q)];
        }
    }
    [[nodiscard]] std::int64_t at(Coord p) const {
        if (p < 0 || p >= length) return 0;
        return pileup[static_cast<std::size_t>(p)];
    }
    [[nodiscard]] std::int64_t sum(Coord p, Coord width) const {
        const Coord half = width / 2;
        const Coord lo = std::max<Coord>(0, p - half);
        const Coord hi = std::min<Coord>(length, p + half);
        std::int64_t total = 0;
        for (Coord q = lo; q < hi; ++q) total += at(q);
        return total;
    }
    [[nodiscard]] std::int64_t width(Coord p, Coord w) const {
        const Coord half = w / 2;
        return std::max<Coord>(0, std::min<Coord>(length, p + half) -
                                      std::max<Coord>(0, p - half));
    }
};

// ---------------------------------------------------------------------------

void test_scales_are_the_macs_widths() {
    // 1kb / 5kb / 10kb, in that order, and the ring must be big enough to hold
    // the largest of them plus the lookahead the centring needs.
    CHECK_EQ(pto::peaks::kScaleCount, 3);
    CHECK_EQ(pto::peaks::kScaleWidth[0], 1000);
    CHECK_EQ(pto::peaks::kScaleWidth[1], 5000);
    CHECK_EQ(pto::peaks::kScaleWidth[2], 10000);
    CHECK_EQ(pto::peaks::kMaxScaleWidth, 10000);
    CHECK(pto::peaks::kRingCapacity > pto::peaks::kMaxScaleWidth +
                                          pto::peaks::kMaxFragmentSpan);
    // A power of two, so the ring index is a mask rather than a division.
    CHECK((pto::peaks::kRingCapacity & (pto::peaks::kRingCapacity - 1)) == 0);
}

void test_the_ring_index_uses_masking_not_division() {
    // The ring index is on the hottest path in the module -- once per base per
    // scale -- and an integer division there costs 20-40 cycles against one
    // for an AND. The capacity is a power of two so that the mask is exact,
    // and this scans the header to keep it that way: a later edit that reaches
    // for `%` would be correct and several times slower, and no output value
    // would show it.
    CHECK((pto::peaks::kRingCapacity & (pto::peaks::kRingCapacity - 1)) == 0);

    const std::filesystem::path header =
        std::filesystem::path(PTO_PEAKS_SOURCE_DIR) / "include" / "peaks" /
        "sliding_window.hpp";
    std::ifstream in(header);
    if (!in) {
        std::fprintf(stderr, "  FAIL cannot open %s for the division scan\n",
                     header.string().c_str());
        ++testing::g_failures;
        ++testing::g_checks;
        return;
    }
    std::string line;
    std::size_t lineno = 0;
    std::size_t scanned = 0;
    while (std::getline(in, line)) {
        ++lineno;
        const std::size_t comment = line.find("//");
        if (comment != std::string::npos) line.erase(comment);
        ++scanned;
        if (line.find('%') != std::string::npos) {
            testing::report(false, "modulo in the sliding window", __FILE__,
                            static_cast<int>(lineno), line);
        } else {
            ++testing::g_checks;
        }
    }
    CHECK(scanned > 200);
}

void test_empty_contig_is_all_zero() {
    MultiScaleWindow w;
    w.reset(0, 5000);
    w.finish();
    int steps = 0;
    while (w.step()) {
        ++steps;
        CHECK_EQ(w.pileup(), 0);
        for (int s = 0; s < pto::peaks::kScaleCount; ++s) CHECK_EQ(w.window_sum(s), 0);
    }
    CHECK_EQ(steps, 5000);
}

void test_single_fragment_pileup_profile() {
    MultiScaleWindow w;
    w.reset(0, 4000);
    CHECK(w.add_fragment(100, 200) == AddResult::kAccepted);
    w.finish();
    while (w.step()) {
        const Coord p = w.position();
        const std::int64_t want = (p >= 100 && p < 200) ? 1 : 0;
        CHECK_MSG(w.pileup() == want,
                  "pos " + std::to_string(p) + ": " + std::to_string(w.pileup()) +
                      " vs " + std::to_string(want));
    }
}

void test_overlapping_fragments_stack() {
    MultiScaleWindow w;
    w.reset(0, 4000);
    CHECK(w.add_fragment(100, 200) == AddResult::kAccepted);
    CHECK(w.add_fragment(150, 250) == AddResult::kAccepted);
    CHECK(w.add_fragment(150, 160) == AddResult::kAccepted);
    w.finish();
    while (w.step()) {
        const Coord p = w.position();
        std::int64_t want = 0;
        if (p >= 100 && p < 200) ++want;
        if (p >= 150 && p < 250) ++want;
        if (p >= 150 && p < 160) ++want;
        CHECK_MSG(w.pileup() == want, "pos " + std::to_string(p));
    }
}

void test_abutting_fragments_do_not_double_count() {
    // [100,200) and [200,300) share the coordinate 200 in the sense that one
    // ends where the other begins. Half-open means the pileup is 1 throughout
    // and never 2 -- an off-by-one in the delta encoding shows up here and
    // nowhere else, and it would inflate every background estimate slightly.
    MultiScaleWindow w;
    w.reset(0, 1000);
    CHECK(w.add_fragment(100, 200) == AddResult::kAccepted);
    CHECK(w.add_fragment(200, 300) == AddResult::kAccepted);
    w.finish();
    while (w.step()) {
        const Coord p = w.position();
        const std::int64_t want = (p >= 100 && p < 300) ? 1 : 0;
        CHECK_MSG(w.pileup() == want, "pos " + std::to_string(p));
    }
}

void test_rolling_sums_match_a_direct_recomputation() {
    // The load-bearing test in this file. An O(1) rolling update is a running
    // difference, and a running difference that is wrong by one at any base
    // stays wrong for the rest of the chromosome.
    constexpr Coord kLen = 30000;
    std::mt19937_64 rng(20260921);
    std::uniform_int_distribution<int> gap(1, 60);
    std::uniform_int_distribution<int> span(40, 400);

    Reference ref(kLen);
    MultiScaleWindow w;
    w.reset(0, kLen);
    std::vector<std::pair<Coord, Coord>> frags;
    Coord at = 0;
    while (at < kLen) {
        const Coord s = at;
        const Coord e = std::min<Coord>(kLen, s + span(rng));
        if (s < e) {
            frags.emplace_back(s, e);
            ref.add(s, e);
        }
        at += gap(rng);
    }

    long long compared = 0;
    drive(w, frags, [&](MultiScaleWindow& win) {
        const Coord p = win.position();
        CHECK_MSG(win.pileup() == ref.at(p), "pileup at " + std::to_string(p));
        for (int sc = 0; sc < pto::peaks::kScaleCount; ++sc) {
            const std::int64_t want = ref.sum(p, pto::peaks::kScaleWidth[sc]);
            CHECK_MSG(win.window_sum(sc) == want,
                      "scale " + std::to_string(pto::peaks::kScaleWidth[sc]) + " at " +
                          std::to_string(p) + ": " + std::to_string(win.window_sum(sc)) +
                          " vs " + std::to_string(want));
            CHECK_MSG(win.window_width(sc) == ref.width(p, pto::peaks::kScaleWidth[sc]),
                      "width scale " + std::to_string(sc) + " at " + std::to_string(p));
        }
        ++compared;
    });
    CHECK_EQ(compared, static_cast<long long>(kLen));
}

void test_windows_are_centred_not_trailing() {
    // One fragment pile at 20000..20100 on a long contig. A CENTRED 1kb window
    // sees it from 19600 through 20500; a trailing one would see it from 20100
    // through 21100. The distinction decides whether a peak is judged against
    // the background beside it or against itself.
    constexpr Coord kLen = 40000;
    MultiScaleWindow w;
    w.reset(0, kLen);
    std::vector<std::pair<Coord, Coord>> frags(10, {20000, 20100});

    std::int64_t at_19600 = -1, at_19400 = -1, at_20600 = -1, at_20400 = -1;
    drive(w, frags, [&](MultiScaleWindow& win) {
        switch (win.position()) {
            case 19400: at_19400 = win.window_sum(0); break;
            case 19600: at_19600 = win.window_sum(0); break;
            case 20400: at_20400 = win.window_sum(0); break;
            case 20600: at_20600 = win.window_sum(0); break;
            default: break;
        }
    });
    // Inside half a window of the pile on BOTH sides, the signal is visible.
    CHECK(at_19600 > 0);
    CHECK(at_20400 > 0);
    // Beyond half a window on both sides, it is not.
    CHECK_EQ(at_19400, 0);
    CHECK_EQ(at_20600, 0);
}

void test_windows_truncate_at_the_contig_start() {
    // Position 100 on a fresh contig has no 10kb of upstream sequence to
    // average over. Dividing by 10000 anyway understates the background by two
    // orders of magnitude and makes every telomeric position a peak.
    constexpr Coord kLen = 30000;
    MultiScaleWindow w;
    w.reset(0, kLen);
    for (int i = 0; i < 50; ++i) CHECK(w.add_fragment(0, 200) == AddResult::kAccepted);
    w.finish();
    CHECK(run_to(w, 100));
    CHECK_EQ(w.window_width(0), 600);    // [0, 600)
    CHECK_EQ(w.window_width(1), 2600);   // [0, 2600)
    CHECK_EQ(w.window_width(2), 5100);   // [0, 5100)
    // The rate uses the truncated width, so it stays a real per-base rate.
    CHECK_NEAR(w.background_rate(0),
               static_cast<double>(w.window_sum(0)) / 600.0, 1e-12);
}

void test_windows_truncate_at_the_contig_end() {
    constexpr Coord kLen = 20000;
    MultiScaleWindow w;
    w.reset(0, kLen);
    std::int64_t w0 = -1, w2 = -1;
    Coord last = -1;
    drive(w, {{19000, 19100}}, [&](MultiScaleWindow& win) {
        last = win.position();
        if (win.position() != kLen - 1) return;
        w0 = win.window_width(0);
        w2 = win.window_width(2);
    });
    CHECK_EQ(w0, 500 + 1);    // [19499, 20000)
    CHECK_EQ(w2, 5000 + 1);   // [14999, 20000)
    CHECK_EQ(last, kLen - 1);
    CHECK(!w.step());         // nothing past the contig end
}

void test_seek_skips_a_coverage_desert() {
    // Centromeres and satellite arrays are megabases with no fragments. Walking
    // them one base at a time cannot produce a peak, and on a real genome it is
    // a measurable fraction of the run.
    constexpr Coord kLen = 200000;
    MultiScaleWindow w;
    w.reset(0, kLen);
    CHECK(w.add_fragment(100, 200) == AddResult::kAccepted);
    // Too close to the coverage just added: the windows would not be empty on
    // the far side, so the jump is refused rather than silently erasing it.
    CHECK(!w.seek(1000));
    CHECK(w.seek(100000));
    CHECK_EQ(w.position(), 100000);
    for (int s = 0; s < pto::peaks::kScaleCount; ++s) CHECK_EQ(w.window_sum(s), 0);
    // Streaming resumes normally on the far side.
    CHECK(w.add_fragment(100500, 100600) == AddResult::kAccepted);
    w.finish();
    CHECK(run_to(w, 100550));
    CHECK_EQ(w.pileup(), 1);
    CHECK(w.seek(150000) == false);  // refused after finish()
}

void test_a_contig_shorter_than_the_widest_window() {
    // A 500 bp contig -- scaffolds like this are all over a draft assembly.
    // Every window is truncated to the contig, and all three rates coincide.
    constexpr Coord kLen = 500;
    MultiScaleWindow w;
    w.reset(0, kLen);
    for (int i = 0; i < 5; ++i) CHECK(w.add_fragment(200, 300) == AddResult::kAccepted);
    w.finish();
    CHECK(run_to(w, 250));
    for (int s = 0; s < pto::peaks::kScaleCount; ++s) {
        CHECK_MSG(w.window_width(s) <= kLen, "scale " + std::to_string(s));
        CHECK(w.window_width(s) > 0);
    }
    CHECK_EQ(w.window_sum(2), 500);  // all 5 fragments x 100 bp, whole contig
}

void test_multi_scale_rates_separate_local_from_broad() {
    // 1000 fragments packed into 500 bp of a 30 kb contig. The 1kb window sees
    // the pile at nearly full density; the 10kb window averages it over twenty
    // times the sequence, so its rate must be an order of magnitude lower.
    // This is the entire reason MACS keeps three scales rather than one.
    constexpr Coord kLen = 30000;
    MultiScaleWindow w;
    w.reset(0, kLen);
    std::vector<std::pair<Coord, Coord>> frags(1000, {15000, 15500});
    double r1 = 0, r5 = 0, r10 = 0, lam0 = 0, lam_big = 0;
    drive(w, frags, [&](MultiScaleWindow& win) {
        if (win.position() != 15250) return;
        r1 = win.background_rate(0);
        r5 = win.background_rate(1);
        r10 = win.background_rate(2);
        lam0 = win.local_lambda(0.0);
        lam_big = win.local_lambda(1e9);
    });
    CHECK_MSG(r1 > r5, "1kb " + testing::describe(r1) + " vs 5kb " + testing::describe(r5));
    CHECK_MSG(r5 > r10, "5kb " + testing::describe(r5) + " vs 10kb " + testing::describe(r10));
    CHECK(r1 / r10 > 5.0);
    // local_lambda takes the largest of the three and the genome background.
    CHECK_NEAR(lam0, r1, 1e-12);
    CHECK_NEAR(lam_big, 1e9, 1e-6);
}

void test_out_of_order_fragments_are_refused_not_mis_binned() {
    // A coordinate-sorted stream is the module's precondition. Accepting a
    // record that arrives late would place its signal at whatever base the
    // ring happens to be pointing at -- a wrong answer reported as success.
    // Refusing it is the fail-closed direction, and the count is what lets the
    // CLI exit 2 instead of printing a plausible track.
    MultiScaleWindow w;
    w.reset(0, 4000);
    CHECK(w.add_fragment(1000, 1100) == AddResult::kAccepted);
    CHECK(w.add_fragment(1200, 1300) == AddResult::kAccepted);
    CHECK(w.add_fragment(1100, 1150) == AddResult::kOutOfOrder);
    CHECK_EQ(w.rejected(AddResult::kOutOfOrder), 1);

    // Equal starts are in order: a pile-up of identical fragments is normal.
    CHECK(w.add_fragment(1200, 1400) == AddResult::kAccepted);

    w.finish();
    while (w.step()) {
        const Coord p = w.position();
        std::int64_t want = 0;
        if (p >= 1000 && p < 1100) ++want;
        if (p >= 1200 && p < 1300) ++want;
        if (p >= 1200 && p < 1400) ++want;
        CHECK_MSG(w.pileup() == want,
                  "rejected fragment leaked into pos " + std::to_string(p));
    }
}

void test_malformed_fragments_are_refused() {
    MultiScaleWindow w;
    w.reset(0, 10000);
    CHECK(w.add_fragment(500, 500) == AddResult::kBadSpan);     // empty
    CHECK(w.add_fragment(600, 500) == AddResult::kBadSpan);     // inverted
    CHECK(w.add_fragment(-5, 100) == AddResult::kOutOfBounds);  // before the contig
    CHECK(w.add_fragment(9990, 10100) == AddResult::kOutOfBounds);  // past the end
    CHECK(w.add_fragment(100, 100 + pto::peaks::kMaxFragmentSpan + 1) ==
          AddResult::kBadSpan);  // longer than the ring can represent
    CHECK_EQ(w.rejected(AddResult::kBadSpan), 3);
    CHECK_EQ(w.rejected(AddResult::kOutOfBounds), 2);
    CHECK_EQ(w.rejected(AddResult::kAccepted), 0);

    w.finish();
    while (w.step()) CHECK_EQ(w.pileup(), 0);
}

void test_coordinates_near_int32_max_do_not_overflow() {
    // REVIEW_2026-08-15 finding (b), in a new module. The window adds a 5000 bp
    // half-width to the evaluation position; on an int32 near the ceiling that
    // is undefined behaviour, and it happens before any bounds test can look
    // at the value. Run this under UBSan -- a pass here with no sanitizer
    // proves only that the machine wrapped quietly.
    const Coord len = kMaxCoord;
    MultiScaleWindow w;
    w.reset(0, len);
    // Two billion empty bases are not worth stepping through, and seek() is
    // what a centromere needs anyway.
    CHECK(w.seek(kMaxCoord - 5000));
    const Coord s = kMaxCoord - 300;
    CHECK(w.add_fragment(s, s + 100) == AddResult::kAccepted);
    w.finish();
    CHECK(run_to(w, kMaxCoord - 250));
    CHECK_EQ(w.pileup(), 1);
    for (int sc = 0; sc < pto::peaks::kScaleCount; ++sc) {
        CHECK(w.window_width(sc) > 0);
        CHECK(w.window_sum(sc) >= 0);
        CHECK(w.background_rate(sc) >= 0.0);
    }
    // The right edge is clamped to the contig, not wrapped to a negative.
    CHECK_EQ(w.window_width(2), 5000 + 250);
}

void test_deep_pileup_does_not_saturate_the_sums() {
    // 200k fragments stacked on one base. The per-base pileup and the 10kb sum
    // are both int64 accumulations; a 32-bit counter would wrap here and a
    // saturating one must say so rather than report a small number.
    MultiScaleWindow w;
    w.reset(0, 20000);
    for (int i = 0; i < 200000; ++i) CHECK(w.add_fragment(10000, 10001) == AddResult::kAccepted);
    w.finish();
    CHECK(run_to(w, 10000));
    CHECK_EQ(w.pileup(), 200000);
    CHECK_EQ(w.window_sum(0), 200000);
    CHECK(!w.saturated());
}

void test_reset_clears_state_between_contigs() {
    MultiScaleWindow w;
    w.reset(0, 4000);
    for (int i = 0; i < 20; ++i) CHECK(w.add_fragment(1000, 1500) == AddResult::kAccepted);
    w.finish();
    while (w.step()) {
    }

    // A second contig must start from nothing: leftover deltas in the ring
    // would appear as signal at whatever base they wrap onto.
    w.reset(1, 4000);
    CHECK_EQ(w.tid(), 1);
    CHECK_EQ(w.rejected(AddResult::kOutOfOrder), 0);
    // Coordinates that were "behind" on the previous contig are fine now.
    CHECK(w.add_fragment(10, 20) == AddResult::kAccepted);
    w.finish();
    while (w.step()) {
        const Coord p = w.position();
        const std::int64_t want = (p >= 10 && p < 20) ? 1 : 0;
        CHECK_MSG(w.pileup() == want, "leak at " + std::to_string(p));
    }
}

void test_streaming_allocates_nothing() {
    // The ring buffers are sized once. Everything after that -- every
    // add_fragment, every step, every query -- runs out of storage that
    // already exists.
    MultiScaleWindow w;
    w.reset(0, 60000);
    std::vector<std::pair<Coord, Coord>> frags;  // built before the guard
    for (Coord s = 0; s + 200 < 60000; s += 37) frags.emplace_back(s, s + 200);
    {
        testing::AllocationGuard guard;
        std::int64_t acc = 0;
        for (const auto& f : frags) {
            AddResult r = w.add_fragment(f.first, f.second);
            while (r == AddResult::kNeedsDrain) {
                if (!w.step()) break;
                acc += w.pileup();
                r = w.add_fragment(f.first, f.second);
            }
            while (w.step()) {
                acc += w.pileup() + w.window_sum(0) + w.window_sum(1) + w.window_sum(2);
                acc += static_cast<std::int64_t>(w.local_lambda(0.5));
            }
        }
        w.finish();
        while (w.step()) acc += w.pileup();
        const long long allocs = guard.count();  // read before any message allocates
        CHECK(acc > 0);
        if (testing::AllocationGuard::counting()) {
            CHECK_MSG(allocs == 0,
                      std::to_string(allocs) + " allocation(s) in the streaming loop");
        } else {
            testing::skip("allocation counting is off under sanitizers");
        }
    }
    // reset() on an already-constructed window must not reallocate either.
    {
        testing::AllocationGuard guard;
        w.reset(2, 1000);
        const long long allocs = guard.count();
        if (testing::AllocationGuard::counting()) {
            CHECK_MSG(allocs == 0, std::to_string(allocs) + " allocation(s) in reset()");
        } else {
            testing::skip("allocation counting is off under sanitizers");
        }
    }
}

void test_step_stops_at_the_evaluable_limit() {
    // Positions whose 10kb window is not yet fully known must not be handed
    // out: their background would be computed from a half-filled window and
    // come back too low.
    MultiScaleWindow w;
    w.reset(0, 100000);
    CHECK(w.add_fragment(0, 100) == AddResult::kAccepted);
    CHECK(w.add_fragment(10000, 10100) == AddResult::kAccepted);
    // The frontier is at 9999, so only bases up to 9999 - 4999 are complete.
    CHECK(w.evaluable_limit() <= 5000);
    while (w.step()) {
    }
    CHECK(w.position() <= 5000);
    CHECK(w.position() > 4000);
    w.finish();
    CHECK(run_to(w, 99999));

    // A fragment far past the ring's lookahead is back pressure, not an error.
    MultiScaleWindow w2;
    w2.reset(0, 100000);
    CHECK(w2.add_fragment(0, 100) == AddResult::kAccepted);
    CHECK(w2.add_fragment(90000, 90100) == AddResult::kNeedsDrain);
    CHECK_EQ(w2.rejected_total(), 0);  // not counted as a rejection
}

// Torture 2026-09-10: one window object reused across many contigs, the way the
// CLI drives it, against a brute-force reference at EVERY base.
//
// test_rolling_sums_match_a_direct_recomputation covers one 30 kb contig with
// dense, short-gapped fragments. What it does not reach is where a ring buffer
// goes wrong: contig lengths either side of the ring capacity, the widest
// half-window and the longest fragment; a contig of 1 bp; fragments hugging
// both contig ends, so the previous contig's last deltas sit in the slots the
// next contig's first bases map onto; and gaps longer than the ring, which
// force the back-pressure drain in the middle of live signal. The reference
// uses prefix sums so every base of ~1.4 Mb of contigs can be checked.
void test_one_window_reused_across_contigs_matches_the_reference() {
    using pto::peaks::kMaxFragmentSpan;
    using pto::peaks::kRingCapacity;
    using pto::peaks::kScaleCount;
    using pto::peaks::kScaleWidth;

    struct Prefix {
        std::vector<std::int64_t> pile;
        std::vector<std::int64_t> cum;  // cum[i] = sum of pile[0..i)
        Coord len;
        explicit Prefix(Coord l) : pile(static_cast<std::size_t>(l) + 1, 0), cum(), len(l) {}
        void add(Coord s, Coord e) {  // difference array, finished by build()
            pile[static_cast<std::size_t>(s)] += 1;
            pile[static_cast<std::size_t>(e)] -= 1;
        }
        void build() {
            std::int64_t run = 0;
            cum.assign(static_cast<std::size_t>(len) + 1, 0);
            for (Coord i = 0; i < len; ++i) {
                run += pile[static_cast<std::size_t>(i)];
                pile[static_cast<std::size_t>(i)] = run;
                cum[static_cast<std::size_t>(i) + 1] = cum[static_cast<std::size_t>(i)] + run;
            }
        }
        [[nodiscard]] std::int64_t sum(Coord p, Coord width) const {
            const std::int64_t half = width / 2;
            const std::int64_t lo = std::max<std::int64_t>(0, p - half);
            const std::int64_t hi = std::min<std::int64_t>(len, p + half);
            return hi <= lo ? 0 : cum[static_cast<std::size_t>(hi)] - cum[static_cast<std::size_t>(lo)];
        }
    };

    const Coord ring = static_cast<Coord>(kRingCapacity);
    const std::vector<Coord> lengths = {
        1, 2, 3, 999, 1000, 1001, 4095, kMaxFragmentSpan, 4097, 4999, 5000, 5001, 9999, 10000,
        10001, ring - kMaxFragmentSpan - 1, ring - 1, ring, ring + 1, 2 * ring - 1, 2 * ring,
        3 * ring + 7, 131072};

    std::mt19937_64 rng(20260910);
    MultiScaleWindow w;
    long long compared = 0;
    long long mismatches = 0;
    std::string first_mismatch;
    std::int32_t tid = 0;

    for (int round = 0; round < 3; ++round) {
        for (const Coord len : lengths) {
            Prefix ref(len);
            std::vector<std::pair<Coord, Coord>> frags;
            std::uniform_int_distribution<int> coin(0, 99);
            Coord at = 0;
            while (at < len) {
                const Coord max_span = std::min<Coord>(kMaxFragmentSpan, len - at);
                std::uniform_int_distribution<Coord> span(1, max_span);
                // A pile of identical fragments now and then: a real peak.
                const int copies = coin(rng) < 5 ? 50 : 1;
                const Coord e = at + span(rng);
                for (int c = 0; c < copies; ++c) frags.emplace_back(at, e);
                const int roll = coin(rng);
                if (roll < 3) {
                    at += static_cast<Coord>(ring + 1 + coin(rng) * 300);  // longer than the ring
                } else if (roll < 60) {
                    at += 0;  // equal starts are in order
                    at += static_cast<Coord>(coin(rng) % 3);
                } else {
                    at += static_cast<Coord>(1 + coin(rng) * 4);
                }
            }
            // Signal right up to both ends of the contig.
            const Coord tail_span = std::min<Coord>(kMaxFragmentSpan, len);
            frags.emplace_back(len - tail_span, len);
            frags.emplace_back(len - 1, len);
            frags.emplace_back(0, std::min<Coord>(len, 7));
            std::stable_sort(frags.begin(), frags.end(),
                             [](const auto& a, const auto& b) { return a.first < b.first; });
            for (const auto& fr : frags) ref.add(fr.first, fr.second);
            ref.build();

            w.reset(tid++, len);
            long long seen = 0;
            drive(w, frags, [&](MultiScaleWindow& win) {
                const Coord p = win.position();
                ++seen;
                bool ok = win.pileup() == ref.pile[static_cast<std::size_t>(p)];
                for (int sc = 0; sc < kScaleCount; ++sc) {
                    ok = ok && win.window_sum(sc) == ref.sum(p, kScaleWidth[sc]);
                }
                if (!ok) {
                    if (mismatches == 0) {
                        first_mismatch = "contig length " + std::to_string(len) + " round " +
                                         std::to_string(round) + " base " + std::to_string(p) +
                                         ": pileup " + std::to_string(win.pileup()) + " vs " +
                                         std::to_string(ref.pile[static_cast<std::size_t>(p)]);
                    }
                    ++mismatches;
                }
            });
            CHECK_MSG(seen == len, "contig length " + std::to_string(len) + ": visited " +
                                       std::to_string(seen) + " bases");
            compared += seen;
        }
    }
    CHECK_MSG(mismatches == 0, first_mismatch);
    CHECK(compared > 1000000);
}

}  // namespace

int main() {
    testing::Suite suite{"peaks/sliding_window", {}};
    suite.add("scales are the MACS widths", test_scales_are_the_macs_widths);
    suite.add("the ring index uses masking not division", test_the_ring_index_uses_masking_not_division);
    suite.add("empty contig is all zero", test_empty_contig_is_all_zero);
    suite.add("single fragment pileup profile", test_single_fragment_pileup_profile);
    suite.add("overlapping fragments stack", test_overlapping_fragments_stack);
    suite.add("abutting fragments do not double count", test_abutting_fragments_do_not_double_count);
    suite.add("rolling sums match a direct recomputation", test_rolling_sums_match_a_direct_recomputation);
    suite.add("windows are centred not trailing", test_windows_are_centred_not_trailing);
    suite.add("windows truncate at the contig start", test_windows_truncate_at_the_contig_start);
    suite.add("windows truncate at the contig end", test_windows_truncate_at_the_contig_end);
    suite.add("seek skips a coverage desert", test_seek_skips_a_coverage_desert);
    suite.add("a contig shorter than the widest window", test_a_contig_shorter_than_the_widest_window);
    suite.add("multi-scale rates separate local from broad", test_multi_scale_rates_separate_local_from_broad);
    suite.add("out-of-order fragments are refused", test_out_of_order_fragments_are_refused_not_mis_binned);
    suite.add("malformed fragments are refused", test_malformed_fragments_are_refused);
    suite.add("coordinates near INT32_MAX do not overflow", test_coordinates_near_int32_max_do_not_overflow);
    suite.add("deep pileup does not saturate", test_deep_pileup_does_not_saturate_the_sums);
    suite.add("reset clears state between contigs", test_reset_clears_state_between_contigs);
    suite.add("streaming allocates nothing", test_streaming_allocates_nothing);
    suite.add("step stops at the evaluable limit", test_step_stops_at_the_evaluable_limit);
    suite.add("one window reused across contigs matches the reference",
              test_one_window_reused_across_contigs_matches_the_reference);
    return suite.run();
}
