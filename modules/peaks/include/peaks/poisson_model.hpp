// SPDX-License-Identifier: MIT
//
// pto-peaks: Poisson tail probabilities, vectorised.
//
// What this header computes
// -------------------------
// ln P(X >= k) for X ~ Poisson(lambda), which is the p-value a peak caller
// assigns to a pileup of k reads over a local background of lambda. It is the
// same quantity MACS2 computes, and the same regularised lower incomplete
// gamma:  P(X >= k) = gamma(k, lambda) / Gamma(k).
//
// Everything is returned as a natural logarithm. Not as a convenience -- a
// significant peak has a p-value of 1e-300 or smaller, which is not a double,
// and the caller wants -log10(p) in the output column anyway. There is no
// entry point here that returns a linear probability, deliberately.
//
// Why there is no call to libm in it
// ----------------------------------
// A peak call evaluates this once per candidate position -- hundreds of
// millions of times per genome -- so the evaluation has to vectorise, and
// `std::log` does not: there is no libm entry point that takes a __m256d, and
// a scalar call inside the loop serialises every lane and clobbers the vector
// registers around it. So the three transcendentals this needs (log, exp,
// log1p) are open-coded from their series, once per ISA, and
// tests/test_poisson_model.cpp scans this file to keep them that way.
//
// Classification (isfinite, isnan) and fabs are not calls and are used freely.
//
// Every polynomial is evaluated as two interleaved Horner chains rather than
// one, because a profile of the finished kernel found both the scalar and the
// vector paths waiting on latency rather than throughput: a 20-term Horner is
// 20 dependent fused multiply-adds, and at four cycles each that dominated
// everything around it. Splitting into even- and odd-indexed coefficients
// halves the depth for one extra multiply. The logarithm additionally uses a
// shorter coefficient count than log1p (kAtanhShortTerms), because its own
// range reduction guarantees a smaller argument.
//
// The shape of the computation
// ----------------------------
// Three paths, and which one answers is decided by the numbers rather than by
// a rule about the inputs.
//
// FAST PATH -- the finite lower sum. P(X >= k) = 1 - sum_{i<k} pmf(i), and
// that sum has exactly k terms, each the previous times lambda/i. For
// k <= kFiniteMaxK it is the cheapest thing available and it has the property
// the vector kernels care about most: a trip count known before the loop
// starts, no division (reciprocals are tabulated), no convergence test, and
// no cross-lane operation of any kind. The lanes stay in the vector registers
// from the first term to the last. It also does not care which side of
// lambda = k the input falls on, so it replaces BOTH branches below for the
// large majority of positions on a real track -- measured at 92% of a
// background track, and it is the reason this file is roughly three times
// faster than the version that only had the two exact branches.
//
// What it cannot do is reach the deep tail: it produces Q and the answer is
// 1 - Q, so once Q is within a few eps of 1 there is nothing left to subtract.
// So it is checked AFTER the fact -- the sum hands back the tail it could
// resolve, and whether that tail survived the subtraction is a property of the
// number, not something predicted from k and lambda. kFiniteMinTail is where
// it stops being allowed to answer, and that constant is measured by
// test_finite_path_gate_is_conservative rather than argued for here.
//
// The other two paths are exact everywhere and are what the enriched positions
// -- the peaks, the reason any of this exists -- actually run. Both are built
// on one leading log-term,
//
//     lead = -lambda + k*ln(lambda) - lnGamma(k)
//
// and differ in the correction, because the two series available converge on
// opposite sides of lambda = k+1:
//
//   lambda < k+1   The tail is a convergent power series and the answer is
//                  lead + ln(series), all of it in log space, so a p-value of
//                  1e-3000 costs nothing and underflows nothing.
//
//   lambda >= k+1  P is near 1. The tail is computed as its complement
//                  Q = exp(lead)*CF through Lentz's continued fraction, and
//                  the answer is log1p(-Q) -- log1p and not log, because a
//                  p-value of 1 - 3.8e-42 is not the same number as 1, and
//                  log(1 - 3.8e-42) is exactly the number that loses it.
//
// lnGamma(k) comes from a 513-entry table of log-factorials for k <= 512 and
// from the Stirling asymptotic series above it. The table is not a
// micro-optimisation: it is a load where the series is a division, a
// logarithm and five fused multiply-adds, and pileup depths in real ChIP and
// ATAC data are almost entirely under 512.
//
// What a Gaussian approximation would buy, and why there isn't one
// ---------------------------------------------------------------
// The obvious way to make a Poisson tail cheap for large k and lambda is to
// stop computing it and use a Normal approximation with a continuity
// correction instead. It is fixed-cost, branchless and vectorises perfectly,
// and it is wrong by amounts that have no place in a p-value: against the
// golden table in the test suite it misses ln P by 25% at k = 100, by a
// factor of 13 at k = 1e5 with lambda = 1e3, and by a factor of 84 at
// k = 1000 with lambda = 1. It converges only where P is already near 1 --
// which is exactly where the answer does not matter. The finite lower sum
// above is the fixed-cost branchless path this file has instead, and unlike
// the Gaussian it is EXACT wherever it is used.
//
// The accurate version of that idea is Temme's uniform asymptotic expansion,
// which is a genuine Poisson-to-Normal transform with a correction series.
// It would help the one regime nothing else makes cheap -- lambda near k with
// k large, where both exact series need O(sqrt(k)) terms. It is not here
// because that regime is a p-value of about 0.5, i.e. not a peak, and because
// its coefficient scheme is several hundred lines that would have to be right
// in four kernels.
//
// Accuracy, stated rather than implied
// ------------------------------------
// `lead` is a difference of three terms that each reach O(k ln k) and then
// cancel, so its ABSOLUTE error is eps times the largest term, not eps times
// the answer. At k = 1e5, lambda = 1e3 the terms are ~1e6 and the answer is
// ~-3.6e5: an absolute error near 1e-10 is the best a double holds, and that
// is what this returns.
//
// Where that lands depends on the branch, and the difference is worth knowing
// before trusting a digit:
//
//   * deep tail (the significant peaks): the error is ADDITIVE in a large
//     answer, so the relative accuracy is ~1e-15. This is the regime the
//     whole file exists for.
//   * P near 1 with a very large k: the error is RELATIVE, because an
//     absolute error in `lead` is a relative error in Q. At k = 1e6 it is
//     ~1e-9. These are p-values of 0.8; nothing downstream can tell.
//
// Closing that second gap means computing ln(pmf) from Loader's saddle-point
// form -- `stirlerr` plus a cancellation-free `bd0` -- instead of from three
// logarithms. It is the known upgrade and it is not free: bd0 needs its own
// series and its own branch in each of the four kernels below. It has not
// been done because no consumer of this header can observe the difference.
//
// Degenerate input fails CLOSED
// -----------------------------
// A lambda that is zero, negative, NaN or infinite returns ln p = 0, i.e.
// p = 1, i.e. not a peak. This is the repository's rule about filters, and
// the direction matters more here than almost anywhere else: the naive
// evaluation of -lambda + k*ln(lambda) as lambda -> 0 tends to -infinity,
// which is a p-value of ZERO, which would report every unmeasured window --
// every centromere, every assembly gap, every region with no background
// signal at all -- as the most significant peak in the genome. A window with
// no background is not evidence of enrichment; it is an absence of evidence.
//
// The same rule covers a series or continued fraction that fails to converge
// inside its iteration cap. Returning the partial sum would understate the
// tail, and understating the tail overstates significance.
#pragma once

#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <limits>

// ---------------------------------------------------------------------------
// ISA gating
//
// Same arrangement as modules/scrna_matrix/include/matrix/simd_math.hpp, and
// for the same reason: the x86 kernels carry per-function `target` attributes
// rather than relying on -mavx2 being applied to the translation unit, so the
// rest of the binary stays baseline x86-64 and a CPU without AVX2 cannot
// execute an AVX2 instruction because some unrelated loop got vectorised.
// PEAKS_HAVE_AVX2 / PEAKS_HAVE_AVX512 are set by CMake after probing that the
// compiler can actually build such a function.
//
// NEON is gated differently and that is not an oversight: Advanced SIMD is
// mandatory in ARMv8-A, so there is no capability to probe, no target
// attribute to apply, and no way to build a binary whose NEON path the CPU
// cannot run. It keeps its place in the dispatch enum only so that a single
// machine can be made to run every compiled path.
// ---------------------------------------------------------------------------
#if defined(__x86_64__) || defined(__i386__) || defined(_M_X64) || defined(_M_IX86)
#define PEAKS_ARCH_X86 1
#endif

#if defined(PEAKS_ARCH_X86) && (defined(PEAKS_HAVE_AVX2) || defined(PEAKS_HAVE_AVX512))
#include <immintrin.h>
#define PEAKS_BUILD_X86_SIMD 1
#endif

#if defined(__aarch64__) || defined(_M_ARM64)
#define PEAKS_ARCH_ARM64 1
#endif

#if defined(PEAKS_ARCH_ARM64) && defined(__ARM_NEON)
#include <arm_neon.h>
#define PEAKS_BUILD_ARM_NEON 1
#endif

namespace pto::peaks {

// ===========================================================================
// Constants and series coefficients
//
// Every polynomial below is the EXACT Taylor series of the function it
// approximates, truncated where the next term falls under 1e-18 relative over
// the argument range the reduction guarantees. None of them is a minimax fit.
// That costs a few terms against a fitted polynomial of the same accuracy and
// buys an error bound anyone can re-derive from the range and the term index,
// which is what makes the tolerances in the test suite arguable rather than
// observed.
// ===========================================================================
namespace detail {

inline constexpr double kInfinity = std::numeric_limits<double>::infinity();
inline constexpr double kQNaN = std::numeric_limits<double>::quiet_NaN();
inline constexpr double kEpsilon = std::numeric_limits<double>::epsilon();

// ln 2, split so that n*kLn2Hi is exact for every n a range reduction can
// produce (kLn2Hi has 21 significant bits).
inline constexpr double kLn2Hi = 6.93147180369123816490e-01;
inline constexpr double kLn2Lo = 1.90821492927058770002e-10;
inline constexpr double kInvLn2 = 1.44269504088896338700e+00;
inline constexpr double kLn10 = 2.30258509299404568402e+00;
inline constexpr double kInvLn10 = 4.34294481903251827651e-01;
inline constexpr double kSqrtHalf = 7.07106781186547524401e-01;
inline constexpr double kHalfLogTwoPi = 9.18938533204672741781e-01;

// Round-to-nearest-integer by addition: x + kShift discards everything below
// the units place, and the low bits of the result hold the integer. Used
// instead of std::nearbyint, which is a call on baseline x86-64 and has no
// vector form.
inline constexpr double kShift = 6755399441055744.0;  // 0x1.8p52
inline constexpr std::uint64_t kShiftBits = 0x4338000000000000ull;

// Overflow and underflow gates for exp(). Outside them the range reduction
// would produce an exponent that does not fit an int, so these are load
// bearing rather than an optimisation.
inline constexpr double kExpOverflow = 709.782712893384;
inline constexpr double kExpUnderflow = -746.0;

// 2*atanh(f) = ln((1+f)/(1-f)) = 2f * (1 + f^2/3 + f^4/5 + ...).
//
// 20 terms, because log1p() drives this at |f| <= 1/3 (f = y/(2+y) with
// |y| <= 1/2), where the first omitted term is 0.111^20/43 ~ 7e-21. The
// logarithm's own range reduction only reaches |f| <= 0.172, where a dozen
// terms would do; it shares the longer polynomial rather than carrying a
// second one, because two nearly-identical coefficient tables in four kernels
// apiece is eight chances for one of them to drift.
inline constexpr int kAtanhTerms = 20;

// The logarithm's own range reduction only ever reaches |f| <= 0.172, where
// the first omitted term at 12 is 0.0294^12/27 ~ 1.6e-20 -- so log pays for
// eight terms it does not need if it shares log1p's polynomial. It did, until
// the profile showed both kernels latency-bound on exactly these chains.
inline constexpr int kAtanhShortTerms = 12;

struct AtanhCoeffs {
  double c[kAtanhTerms];
};

inline constexpr AtanhCoeffs make_atanh_coeffs() {
  AtanhCoeffs a{};
  for (int j = 0; j < kAtanhTerms; ++j) {
    a.c[j] = 1.0 / static_cast<double>(2 * j + 3);
  }
  return a;
}

inline constexpr AtanhCoeffs kAtanh = make_atanh_coeffs();

// exp(r) = 1 + r + r^2 * (1/2! + r/3! + ...) for |r| <= ln(2)/2 = 0.3466.
// 12 terms: the first omitted is 0.3466^14/14! ~ 4e-18, and pulling r^2 out
// front damps the polynomial's own error by a further factor of 8.
inline constexpr int kExpTerms = 12;

struct ExpCoeffs {
  double c[kExpTerms];
};

inline constexpr ExpCoeffs make_exp_coeffs() {
  ExpCoeffs e{};
  double fact = 2.0;  // 2!
  for (int j = 0; j < kExpTerms; ++j) {
    e.c[j] = 1.0 / fact;
    fact *= static_cast<double>(j + 3);
  }
  return e;
}

inline constexpr ExpCoeffs kExp = make_exp_coeffs();

// Stirling: lnGamma(z) = (z-1/2)ln z - z + ln(2pi)/2 + sum B_2n/(2n(2n-1)z^(2n-1)).
// Five correction terms, valid from z = 16 up: the first omitted term is
// 691/(360360 * 16^11) = 1.1e-16, which is 0.03 of an ulp of lnGamma(16).
inline constexpr double kStirlingMin = 16.0;
inline constexpr double kS0 = 1.0 / 12.0;
inline constexpr double kS1 = -1.0 / 360.0;
inline constexpr double kS2 = 1.0 / 1260.0;
inline constexpr double kS3 = -1.0 / 1680.0;
inline constexpr double kS4 = 1.0 / 1188.0;

// The finite lower sum: how far it reaches, and when it is allowed to answer.
//
// P(X >= k) = 1 - sum_{i<k} pmf(i), and that sum is FINITE -- exactly k terms,
// each the previous one times lambda/i. Below kFiniteMaxK it is cheaper than
// either infinite series, it needs no convergence test (so a vector loop runs
// a fixed number of iterations in lockstep), and with reciprocals tabulated it
// needs no division either.
//
// What it cannot do is reach the deep tail. The sum produces Q, and the answer
// is 1 - Q; once Q is within a few eps of 1 there is no tail left to subtract.
// The error in ln P is about eps*Q/P, so admissibility is a condition on P
// itself: kFiniteMinTail is the smallest upper tail for which the subtraction
// still carries the answer. It is a power of two so the comparison is exact,
// and it is validated by a dense sweep rather than argued -- see
// test_finite_path_gate_is_conservative, which walks the whole (k, lambda)
// grid and asserts the margin against the error model the suite already uses.
//
// Above kFiniteMaxK, or below kFiniteMinTail, the exact series and continued
// fraction below answer instead. On a real background track that is under 1%
// of positions -- and they are the enriched ones, where the series converges
// fastest anyway, because lambda/k is then well under 1.
inline constexpr std::uint32_t kFiniteMaxK = 64;
inline constexpr double kFiniteMinTail = 0.03125;  // 2^-5

// The finite sum walks UP from pmf(0) = exp(-lambda), so every term inherits
// that first one's precision. Past lambda = -ln(DBL_MIN) it is subnormal and
// past 745 it is zero outright, and the sum then reports a tail of zero where
// the exact path -- which works in logs and never materialises exp(-lambda) --
// still resolves one. The values involved are around 1e-280, so no consumer
// could tell; the bound is here anyway, because a fast path that agrees with
// the reference everywhere it runs is one fewer thing to reason about, and
// this costs a comparison in a regime no background track reaches (708 reads
// of local background at a single base).
inline constexpr double kFiniteMaxLambda = 708.0;

// Iteration caps. Both series converge in O(sqrt(k)) terms when lambda is
// near k and in a handful of terms anywhere else. The cap exists to bound the
// loop, not to be reached, and exhausting it fails closed (p = 1).
//
// It has to cover the largest k the caller can pass, which is the clamp in
// PeakCaller::clamp_count: 2^32 - 1. Near lambda = k the series' terms fall as
// exp(-n^2 / 2k), so reaching kSeriesEps takes about sqrt(2 k ln(1e18)) terms
// -- roughly 600,000 at the ceiling. This constant used to be 100,000 with a
// comment saying that covered "every k a pileup can hold"; it covered k up to
// about 1e8. Measured against mpmath at 40 digits: a base 4 sigma above its
// background at k = 1e9 (exact -log10 p = 4.4996) and 4 and 9 sigma at
// k = 2^32 - 1 (4.4995, 9.0063) all came back p = 1 -- significant bases
// silently scored as nothing, in every kernel. The cost is paid only by the
// positions that need it, because every loop exits at convergence: measured on
// Apple M4, the worst case (k = lambda = 2^32 - 1) is 0.96 ms scalar and
// 0.33 ms per lane batched on NEON, k = 1e6 is ~10 us, and a typical peak base
// is unchanged.
inline constexpr int kMaxIterations = 1 << 20;

// How often the vector kernels test their convergence masks. Between checks
// the lanes run in lockstep with no cross-lane operation at all; the check is
// a horizontal reduction, and doing one per iteration is what serialises an
// otherwise vector loop. Overshooting a converged lane by up to a stride is
// harmless in both loops: the series' terms are still shrinking, and the
// continued fraction's delta is already within eps of 1.
inline constexpr int kConvergenceStride = 8;
inline constexpr double kSeriesEps = 1e-18;
inline constexpr double kCfEps = 4.0 * kEpsilon;
// Lentz's method renormalises any partial denominator that approaches zero.
inline constexpr double kCfTiny = 1e-300;

}  // namespace detail

// The largest k with a tabulated log-factorial. Public because the test suite
// checks the seam on both sides of it, and because a caller sizing a batch may
// reasonably care where the cheap path ends.
inline constexpr std::uint32_t kLogFactorialTableMax = 512;

// ===========================================================================
// Scalar kernels
//
// These are the correctness contract. Every vector kernel below computes the
// same expressions in the same order and is required to agree with these to
// within the tolerance test_poisson_model.cpp derives -- so when a parity test
// fails, this is the side that is right.
// ===========================================================================

// 2*atanh(f), for |f| <= 0.35.
//
// Written as 2f + 2f*r*P(r) rather than as 2f*(1 + r*P(r)): the leading 2f is
// exact (a multiply by a power of two), and everything the polynomial can get
// wrong is confined to a correction that is at most 1.3% of the result.
//
// P is evaluated as two interleaved Horner chains -- the even-indexed
// coefficients and the odd ones, both in r^2 -- recombined once at the end.
// Plain Horner is a strictly serial chain of N fused multiply-adds, and at
// four cycles each that is what both the scalar and the vector kernels were
// waiting on: splitting it halves the depth for one extra multiply, and the
// two chains issue in parallel on every machine this targets.
template <int Terms>
[[nodiscard]] constexpr double atanh2_n(double f) noexcept {
  static_assert(Terms >= 4 && Terms <= detail::kAtanhTerms);
  const double r = f * f;
  const double r2 = r * r;
  constexpr int last_even = ((Terms - 1) / 2) * 2;
  constexpr int last_odd = ((Terms - 2) / 2) * 2 + 1;
  double even = detail::kAtanh.c[last_even];
  for (int j = last_even - 2; j >= 0; j -= 2) even = even * r2 + detail::kAtanh.c[j];
  double odd = detail::kAtanh.c[last_odd];
  for (int j = last_odd - 2; j >= 1; j -= 2) odd = odd * r2 + detail::kAtanh.c[j];
  const double acc = even + odd * r;
  const double two_f = f + f;
  return two_f + two_f * r * acc;
}

// The full-range form, for log1p, whose argument reaches |f| = 1/3.
[[nodiscard]] constexpr double atanh2(double f) noexcept {
  return atanh2_n<detail::kAtanhTerms>(f);
}

// Natural logarithm. Domain and special values match libm exactly: log(0) is
// -inf, log of a negative is NaN, log(+inf) is +inf, log(1) is +0.
[[nodiscard]] constexpr double log_approx(double x) noexcept {
  if (x != x) return x;
  if (x < 0.0) return detail::kQNaN;
  if (x == 0.0) return -detail::kInfinity;
  if (x == detail::kInfinity) return x;

  int e = 0;
  std::uint64_t bits = std::bit_cast<std::uint64_t>(x);
  if (bits < 0x0010000000000000ull) {
    // Subnormal: the exponent field is zero and the mantissa is not
    // normalised, so the split below would read a garbage exponent. Scale into
    // the normal range first and pay it back in `e`.
    x = x * 0x1p54;
    bits = std::bit_cast<std::uint64_t>(x);
    e = -54;
  }
  e += static_cast<int>(bits >> 52) - 1022;
  double m = std::bit_cast<double>((bits & 0x000FFFFFFFFFFFFFull) |
                                   0x3FE0000000000000ull);  // m in [0.5, 1)
  if (m < detail::kSqrtHalf) {
    // Centre the mantissa on 1, so that |f| <= 0.172 and the leading e*ln2
    // term can never cancel against ln(m): |e*ln2| >= 0.693 whenever e != 0.
    m = m + m;
    e -= 1;
  }
  const double f = (m - 1.0) / (m + 1.0);
  const double de = static_cast<double>(e);
  return de * detail::kLn2Hi + (de * detail::kLn2Lo + atanh2_n<detail::kAtanhShortTerms>(f));
}

// 2^n as a bit pattern, for |n| <= 1023.
[[nodiscard]] constexpr double pow2i(int n) noexcept {
  return std::bit_cast<double>(static_cast<std::uint64_t>(n + 1023) << 52);
}

// e^x. Gradual underflow is preserved: exp(-745) is a denormal, not zero.
[[nodiscard]] constexpr double exp_approx(double x) noexcept {
  if (x != x) return x;
  if (x > detail::kExpOverflow) return detail::kInfinity;
  if (x < detail::kExpUnderflow) return 0.0;

  const double t = x * detail::kInvLn2 + detail::kShift;
  const double dn = t - detail::kShift;  // = round(x / ln 2)
  const int n = static_cast<int>(dn);
  // Two-step subtraction: dn*kLn2Hi is exact, so the residual keeps every bit
  // it has, and kLn2Lo carries the rest of ln 2.
  const double r = (x - dn * detail::kLn2Hi) - dn * detail::kLn2Lo;

  // Even/odd split, for the reason atanh2_n gives.
  const double r2 = r * r;
  double even = detail::kExp.c[detail::kExpTerms - 2];
  for (int j = detail::kExpTerms - 4; j >= 0; j -= 2) even = even * r2 + detail::kExp.c[j];
  double odd = detail::kExp.c[detail::kExpTerms - 1];
  for (int j = detail::kExpTerms - 3; j >= 1; j -= 2) odd = odd * r2 + detail::kExp.c[j];
  const double p = even + odd * r;
  const double er = 1.0 + (r + r2 * p);

  // Scaled in two halves. 2^n is not representable for |n| > 1023, and n
  // reaches 1024 for arguments just under the overflow threshold whose true
  // value is still finite -- splitting the scale keeps both ends exact and
  // needs no branch, which is what the vector kernels want.
  const int n1 = n >> 1;
  return (er * pow2i(n1)) * pow2i(n - n1);
}

// ln(1 + y), accurate as y -> 0 where ln(1 + y) is not.
[[nodiscard]] constexpr double log1p_approx(double y) noexcept {
  if (y != y) return y;
  if (y < -1.0) return detail::kQNaN;
  if (y == -1.0) return -detail::kInfinity;
  if (y >= -0.5 && y <= 0.5) {
    // f = y/(2+y) maps [-1/2, 1/2] into [-1/3, 1/5], inside the polynomial's
    // range, and 2+y is exact for every y here.
    return atanh2(y / (2.0 + y));
  }
  return log_approx(1.0 + y);
}

namespace detail {

// The branch-free half of log1p, for the continued-fraction path where the
// argument is known to be in [-1/2, 0]. Q > 1/2 is unreachable there: the
// branch is taken only when lambda >= k+1, and the median of a Poisson exceeds
// lambda - ln 2, so P(X <= k-1) < 1/2 by construction.
[[nodiscard]] constexpr double log1p_bounded(double y) noexcept {
  return atanh2(y / (2.0 + y));
}

// lnGamma by the Stirling asymptotic series. Valid for z >= kStirlingMin;
// callers are responsible for shifting smaller arguments up.
[[nodiscard]] constexpr double lgamma_stirling(double z) noexcept {
  const double w = 1.0 / z;
  const double w2 = w * w;
  const double s = w * (kS0 + w2 * (kS1 + w2 * (kS2 + w2 * (kS3 + w2 * kS4))));
  return (z - 0.5) * log_approx(z) - z + kHalfLogTwoPi + s;
}

// ln(k!) for k = 0 .. kLogFactorialTableMax, built at compile time.
//
// Up to 20 the factorial itself is an exact double, so one logarithm of an
// exact integer is the shortest path to a correctly-rounded answer. Above it
// the Stirling series is already good to a fraction of an ulp at z = 21, and
// summing 500 logarithms instead would accumulate more error than either.
struct LogFactorialTable {
  double v[kLogFactorialTableMax + 1];
};

inline constexpr LogFactorialTable make_log_factorial_table() {
  LogFactorialTable t{};
  t.v[0] = 0.0;
  double fact = 1.0;
  for (std::uint32_t k = 1; k <= 20; ++k) {
    fact *= static_cast<double>(k);
    t.v[k] = log_approx(fact);
  }
  for (std::uint32_t k = 21; k <= kLogFactorialTableMax; ++k) {
    t.v[k] = lgamma_stirling(static_cast<double>(k) + 1.0);
  }
  return t;
}

inline constexpr LogFactorialTable kLogFactorial = make_log_factorial_table();

// 1/i for the finite lower sum, so its inner loop is two multiplies rather
// than a divide. Each entry is a correctly-rounded quotient, and the extra
// rounding per term is worth it: FDIV is ~10 cycles and does not pipeline on
// either NEON or AVX2, where FMUL is 3-4 and does.
struct ReciprocalTable {
  double v[kFiniteMaxK + 1];
};

inline constexpr ReciprocalTable make_reciprocal_table() {
  ReciprocalTable t{};
  t.v[0] = 0.0;
  for (std::uint32_t i = 1; i <= kFiniteMaxK; ++i) {
    t.v[i] = 1.0 / static_cast<double>(i);
  }
  return t;
}

inline constexpr ReciprocalTable kReciprocal = make_reciprocal_table();

// Q = P(X <= k-1) = sum_{i=0}^{k-1} e^-lambda lambda^i / i!, for k <= kFiniteMaxK.
//
// Exactly k terms, no division, no convergence test, no early exit -- which is
// what lets the vector form run a fixed trip count with the lanes in lockstep.
// The recurrence cannot overflow: t starts at exp(-lambda) <= 1 and every term
// is a Poisson probability, so the running term is bounded by 1 whatever
// lambda is. For lambda > 745 the first term underflows to zero and the whole
// sum is zero, which is the right answer -- P(X >= k) is 1 there.
[[nodiscard]] inline double finite_lower_sum(std::uint32_t k, double lambda) noexcept {
  double t = ::pto::peaks::exp_approx(-lambda);
  double q = t;
  for (std::uint32_t i = 1; i < k; ++i) {
    t *= lambda * kReciprocal.v[i];
    q += t;
  }
  return q;
}

}  // namespace detail

// ln(k!). A load for every pileup depth a sequencing run realistically
// produces, and the Stirling series above that.
[[nodiscard]] inline double log_factorial(std::uint32_t k) noexcept {
  if (k <= kLogFactorialTableMax) return detail::kLogFactorial.v[k];
  return detail::lgamma_stirling(static_cast<double>(k) + 1.0);
}

// lnGamma(z) for z > 0.
//
// Non-positive arguments return +infinity rather than a reflected value. Every
// non-positive integer is a pole, and the alternative -- computing the
// reflection formula for the non-integers in between -- would answer a
// question nothing in a peak caller asks, in a way that turns a bad argument
// into a plausible number.
[[nodiscard]] inline double lgamma_approx(double z) noexcept {
  if (z != z) return z;
  if (!(z > 0.0)) return detail::kInfinity;
  if (z == detail::kInfinity) return z;
  // The two exact zeros. Reached by the shift path they would come out as a
  // difference of two numbers near 28, i.e. as some 1e-14; libm returns the
  // zero, and so does this, because a downstream sum that accumulates 1e-14
  // forever is a worse answer than one that accumulates nothing.
  if (z == 1.0 || z == 2.0) return 0.0;
  if (z >= detail::kStirlingMin) return detail::lgamma_stirling(z);

  // Shift up to the Stirling range: lnGamma(z) = lnGamma(z+n) - ln(z(z+1)...).
  // At most 16 steps, and the product cannot overflow (it is bounded by 16!)
  // nor underflow (the smallest denormal times 15! is still representable).
  double p = 1.0;
  double y = z;
  while (y < detail::kStirlingMin) {
    p = p * y;
    y = y + 1.0;
  }
  return detail::lgamma_stirling(y) - log_approx(p);
}

namespace detail {

// sum_{n>=0} x^n * Gamma(a)/Gamma(a+n+1) = (1/a)(1 + x/(a+1) + ...).
//
// Converges for x < a+1; the ratio between successive terms is x/(a+n) < 1 and
// falling. `converged` reports whether the cap was reached, because a truncated
// sum understates the tail, and understating the tail overstates significance.
[[nodiscard]] inline double gamma_series(double a, double x, bool& converged) noexcept {
  double ap = a;
  double del = 1.0 / a;  // a >= 1 here; no division by zero is reachable
  double sum = del;
  converged = false;
  for (int i = 0; i < kMaxIterations; ++i) {
    ap += 1.0;
    del *= x / ap;
    sum += del;
    if (del <= sum * kSeriesEps) {
      converged = true;
      break;
    }
  }
  return sum;
}

// Lentz's modified continued fraction for the upper incomplete gamma, giving h
// with Q(a,x) = exp(lead) * h. Converges for x >= a+1, where the first partial
// denominator b = x+1-a is at least 2.
[[nodiscard]] inline double gamma_cf(double a, double x, bool& converged) noexcept {
  double b = x + 1.0 - a;
  double c = 1.0 / kCfTiny;
  double d = 1.0 / b;
  double h = d;
  converged = false;
  for (int i = 1; i <= kMaxIterations; ++i) {
    const double an = -static_cast<double>(i) * (static_cast<double>(i) - a);
    b += 2.0;
    d = an * d + b;
    if (std::fabs(d) < kCfTiny) d = kCfTiny;  // renormalise; never divide by 0
    c = b + an / c;
    if (std::fabs(c) < kCfTiny) c = kCfTiny;
    d = 1.0 / d;
    const double del = d * c;
    h *= del;
    if (std::fabs(del - 1.0) <= kCfEps) {
      converged = true;
      break;
    }
  }
  return h;
}

// The one place the final answer is sanitised, shared by every kernel so that
// scalar and vector agree bit for bit on the degenerate cases. A p-value lies
// in [0, 1], so its logarithm lies in [-inf, 0]; anything else -- a NaN from a
// pathological input, a +1e-17 from rounding at p = 1 -- becomes 0, which is
// p = 1, which is not a peak.
[[nodiscard]] inline double finalise(double out) noexcept {
  return (out < 0.0 && out > -kInfinity) ? out : 0.0;
}

}  // namespace detail

namespace detail {

// The exact path: the infinite upper-tail series below lambda = k+1 and
// Lentz's continued fraction above it. Assumes k >= 1 and a finite positive
// lambda; log_poisson_sf applies those guards.
//
// Exposed under detail:: rather than kept inside log_poisson_sf so that the
// suite can drive it directly. That is what makes the finite path's
// admissibility gate testable: the two paths are compared against each other
// over the whole (k, lambda) grid, rather than the constant being argued for
// in a comment.
[[nodiscard]] inline double sf_exact(std::uint32_t k, double lambda) noexcept {
  const double a = static_cast<double>(k);
  const double lead =
      -lambda + a * ::pto::peaks::log_approx(lambda) -
      ::pto::peaks::log_factorial(k - 1);  // lnGamma(k)

  bool converged = false;
  double out = 0.0;
  if (lambda < a + 1.0) {
    const double sum = gamma_series(a, lambda, converged);
    out = lead + ::pto::peaks::log_approx(sum);
  } else {
    const double h = gamma_cf(a, lambda, converged);
    double q = ::pto::peaks::exp_approx(lead) * h;
    // Clamped, not corrected: Q < 1/2 is a property of this branch, so the
    // clamp only ever absorbs rounding noise at the extremes. It also keeps
    // the argument inside log1p_bounded's polynomial range unconditionally.
    if (!(q > 0.0)) q = 0.0;
    if (q > 0.5) q = 0.5;
    out = log1p_bounded(-q);
  }
  if (!converged) return 0.0;
  return finalise(out);
}

}  // namespace detail

// ln P(X >= k) for X ~ Poisson(lambda).
//
// Returns 0 (p = 1) for k = 0, which is exact, and for any lambda that is not
// finite and positive, which is the fail-closed rule.
[[nodiscard]] inline double log_poisson_sf(std::uint32_t k, double lambda) noexcept {
  // Rejects negatives, zero and NaN in one comparison, and infinity in the
  // second. Written as `!(lambda > 0)` rather than `lambda <= 0` because NaN
  // compares false against both and would fall through the second form.
  if (!(lambda > 0.0) || !(lambda < detail::kInfinity)) return 0.0;
  if (k == 0) return 0.0;

  // Fast path: the finite lower sum. Handles both sides of lambda = k with one
  // loop and no branch, and answers for the large majority of positions on a
  // real track. It is checked AFTER the fact rather than predicted: the sum
  // hands back the tail it could resolve, and whether that tail survived the
  // subtraction is a property of the number, not of a rule about the inputs.
  if (k <= detail::kFiniteMaxK && lambda <= detail::kFiniteMaxLambda) {
    const double q = detail::finite_lower_sum(k, lambda);
    if (1.0 - q >= detail::kFiniteMinTail) {
      return detail::finalise(log1p_approx(-q));
    }
  }
  return detail::sf_exact(k, lambda);
}

// -log10 p, the column MACS2 writes and the one a browser track displays.
[[nodiscard]] inline double neg_log10_poisson_sf(std::uint32_t k, double lambda) noexcept {
  return -log_poisson_sf(k, lambda) * detail::kInvLn10;
}

// An upper bound on neg_log10_poisson_sf, cheap enough to solve for lambda
// (peaks/caller.hpp's skip bands do exactly that) where the exact score has
// no closed form to solve at all.
//
// P(X >= k) >= P(X = k), the single term the sum's own first element already
// is, so -log10 P(X >= k) <= -log10 P(X = k) = (lambda - k*ln(lambda) +
// ln(k!)) / ln(10). Undefined at k = 0 (ln(0) is -infinity, and correctly so:
// P(X=0) = P(X>=0) = 1 exactly there, so the "bound" would be the exact
// answer already) -- callers needing k = 0 have detail::provably_below_cutoff
// for that, which is exact rather than a bound.
//
// Strictly convex in lambda: d/dlambda = (1 - k/lambda) / ln(10), zero
// exactly at lambda = k (where the bound is tightest -- P(X=k) is the mode's
// own probability) and positive on either side, so for any fixed k the set
// {lambda : this is below some cutoff} is a single open interval, not a
// half-line or a union of two.
[[nodiscard]] inline double poisson_pmf_upper_bound(std::uint32_t k, double lambda) noexcept {
  return (lambda - static_cast<double>(k) * log_approx(lambda) + log_factorial(k)) *
         detail::kInvLn10;
}

// ===========================================================================
// Vector kernels
//
// One per ISA, written out rather than shared behind a traits template. That
// is not a preference: a `target`-attributed helper cannot be inlined into a
// function without the same attribute, so a single algorithm template over a
// per-ISA operation traits class either fails to compile or degrades every
// vector operation into a call. modules/scrna_matrix/include/matrix/simd_math.hpp
// reached the same arrangement for the same reason.
//
// What keeps the copies honest is the parity suite: every kernel compiled into
// the binary is driven explicitly through *_batch_isa and compared against the
// scalar reference on a workload built to reach both branches, both sides of
// the log-factorial table seam, and every degenerate lambda.
//
// Three properties are shared by all of them and worth stating once:
//
//   * Lanes that would take the other branch are given BENIGN inputs, not
//     masked out. The series is run with x = 0 (converges on the first term)
//     and the continued fraction with a = 1, x = 64 (converges in about five).
//     Masking the arithmetic instead would leave the inactive lanes carrying
//     whatever the other branch's inputs produce -- for the series that is a
//     divergent term that reaches infinity in a few hundred iterations and
//     then contaminates the horizontal convergence test.
//   * A degenerate lambda is replaced by 1.0 before any arithmetic and the
//     lane is zeroed at the end. A NaN allowed through the kernel would come
//     out as a NaN p-value, and a NaN compares false against every threshold:
//     the position is silently skipped rather than visibly wrong.
//   * The table index is clamped in DOUBLE before it is converted, so the
//     gather cannot address past the table however large k is. A gather with
//     an unclamped count reads whatever follows the table and returns a
//     plausible double.
// ===========================================================================

#if defined(PEAKS_BUILD_ARM_NEON)
namespace detail {
namespace neon {

using V = float64x2_t;
using M = uint64x2_t;

// Even/odd split Horner; see the scalar atanh2_n for why the chain is halved.
template <int Terms>
inline V vatanh2_n(V f) noexcept {
  const V r = vmulq_f64(f, f);
  const V r2 = vmulq_f64(r, r);
  constexpr int last_even = ((Terms - 1) / 2) * 2;
  constexpr int last_odd = ((Terms - 2) / 2) * 2 + 1;
  V even = vdupq_n_f64(kAtanh.c[last_even]);
  for (int j = last_even - 2; j >= 0; j -= 2) {
    even = vfmaq_f64(vdupq_n_f64(kAtanh.c[j]), even, r2);
  }
  V odd = vdupq_n_f64(kAtanh.c[last_odd]);
  for (int j = last_odd - 2; j >= 1; j -= 2) {
    odd = vfmaq_f64(vdupq_n_f64(kAtanh.c[j]), odd, r2);
  }
  const V acc = vfmaq_f64(even, odd, r);
  const V two_f = vaddq_f64(f, f);
  return vfmaq_f64(two_f, vmulq_f64(two_f, r), acc);
}

inline V vatanh2(V f) noexcept { return vatanh2_n<kAtanhTerms>(f); }

// log for arguments the caller has already established are positive and
// finite. Subnormals are handled; +-inf, zero and NaN are not, because no
// caller in this file can present one.
inline V vlog(V x) noexcept {
  const M sub = vcltq_f64(x, vdupq_n_f64(2.2250738585072014e-308));
  const V xs = vbslq_f64(sub, vmulq_f64(x, vdupq_n_f64(0x1p54)), x);
  const int64x2_t adj = vbslq_s64(sub, vdupq_n_s64(-54), vdupq_n_s64(0));

  const uint64x2_t bits = vreinterpretq_u64_f64(xs);
  int64x2_t e = vaddq_s64(vreinterpretq_s64_u64(vshrq_n_u64(bits, 52)),
                          vdupq_n_s64(-1022));
  e = vaddq_s64(e, adj);
  V m = vreinterpretq_f64_u64(
      vorrq_u64(vandq_u64(bits, vdupq_n_u64(0x000FFFFFFFFFFFFFull)),
                vdupq_n_u64(0x3FE0000000000000ull)));
  const M lo = vcltq_f64(m, vdupq_n_f64(kSqrtHalf));
  m = vbslq_f64(lo, vaddq_f64(m, m), m);
  e = vaddq_s64(e, vbslq_s64(lo, vdupq_n_s64(-1), vdupq_n_s64(0)));

  const V f = vdivq_f64(vsubq_f64(m, vdupq_n_f64(1.0)), vaddq_f64(m, vdupq_n_f64(1.0)));
  const V de = vcvtq_f64_s64(e);
  return vfmaq_f64(vfmaq_f64(vatanh2_n<kAtanhShortTerms>(f), de, vdupq_n_f64(kLn2Lo)),
                   de, vdupq_n_f64(kLn2Hi));
}

inline V vpow2(int64x2_t n) noexcept {
  return vreinterpretq_f64_s64(vshlq_n_s64(vaddq_s64(n, vdupq_n_s64(1023)), 52));
}

// exp, clamped to the range over which the reduction is exact rather than
// branching to +-inf. Every argument this file presents is a `lead` term
// bounded by ln(sqrt(k/2pi)) < 11, so the clamp is unreachable here; it exists
// so the exponent arithmetic below cannot overflow an int64 on any input.
inline V vexp(V x0) noexcept {
  const V x = vminq_f64(vmaxq_f64(x0, vdupq_n_f64(kExpUnderflow)),
                        vdupq_n_f64(kExpOverflow));
  const V t = vfmaq_f64(vdupq_n_f64(kShift), x, vdupq_n_f64(kInvLn2));
  const V dn = vsubq_f64(t, vdupq_n_f64(kShift));
  const int64x2_t n = vsubq_s64(vreinterpretq_s64_f64(t),
                                vdupq_n_s64(static_cast<std::int64_t>(kShiftBits)));

  V r = vfmsq_f64(x, dn, vdupq_n_f64(kLn2Hi));
  r = vfmsq_f64(r, dn, vdupq_n_f64(kLn2Lo));

  const V r2 = vmulq_f64(r, r);
  V ev = vdupq_n_f64(kExp.c[kExpTerms - 2]);
  for (int j = kExpTerms - 4; j >= 0; j -= 2) ev = vfmaq_f64(vdupq_n_f64(kExp.c[j]), ev, r2);
  V od = vdupq_n_f64(kExp.c[kExpTerms - 1]);
  for (int j = kExpTerms - 3; j >= 1; j -= 2) od = vfmaq_f64(vdupq_n_f64(kExp.c[j]), od, r2);
  const V p = vfmaq_f64(ev, od, r);
  const V er = vaddq_f64(vdupq_n_f64(1.0), vfmaq_f64(r, r2, p));

  const int64x2_t n1 = vshrq_n_s64(n, 1);
  return vmulq_f64(vmulq_f64(er, vpow2(n1)), vpow2(vsubq_s64(n, n1)));
}

inline V vlog1p_bounded(V y) noexcept {
  return vatanh2(vdivq_f64(y, vaddq_f64(vdupq_n_f64(2.0), y)));
}

// ln(1 + y) over the whole range the finite sum can produce, in ONE logarithm.
//
// The bounded form above needs |y| <= 0.35 and the finite sum reaches -0.97, so
// a select between two evaluations would be the obvious construction -- and it
// would put a second 20-term polynomial in the hottest loop in the file. This
// is the standard alternative: take log of the rounded u = 1 + y and add back
// the error that rounding introduced, c = (y - (u-1))/u. When y is small
// enough that u rounds to 1, log(u) is 0 and c is exactly y, which is the case
// the whole function exists for -- a p-value of 1 - 3.8e-42 must not come back
// as 1. u is clamped away from zero for the lanes that are about to be
// discarded anyway; an admitted lane always has u >= kFiniteMinTail.
inline V vlog1p_full(V y) noexcept {
  const V one = vdupq_n_f64(1.0);
  const V u = vmaxq_f64(vaddq_f64(one, y), vdupq_n_f64(kCfTiny));
  const V c = vdivq_f64(vsubq_f64(y, vsubq_f64(u, one)), u);
  return vaddq_f64(vlog(u), c);
}

// Horizontal max, by lane extraction rather than a reduction intrinsic: two
// lanes, and this runs once per vector to choose a trip count, never inside a
// loop.
inline double hmax(V v) noexcept {
  const double a0 = vgetq_lane_f64(v, 0);
  const double a1 = vgetq_lane_f64(v, 1);
  return a0 > a1 ? a0 : a1;
}

inline V vlgamma_stirling(V z) noexcept {
  const V w = vdivq_f64(vdupq_n_f64(1.0), z);
  const V w2 = vmulq_f64(w, w);
  V s = vdupq_n_f64(kS4);
  s = vfmaq_f64(vdupq_n_f64(kS3), s, w2);
  s = vfmaq_f64(vdupq_n_f64(kS2), s, w2);
  s = vfmaq_f64(vdupq_n_f64(kS1), s, w2);
  s = vfmaq_f64(vdupq_n_f64(kS0), s, w2);
  s = vmulq_f64(s, w);
  const V t = vfmaq_f64(vsubq_f64(vdupq_n_f64(kHalfLogTwoPi), z),
                        vsubq_f64(z, vdupq_n_f64(0.5)), vlog(z));
  return vaddq_f64(t, s);
}

inline bool all_true(M mask) noexcept {
  return vminvq_u32(vreinterpretq_u32_u64(mask)) == 0xFFFFFFFFu;
}

inline bool any_true(M mask) noexcept {
  return vmaxvq_u32(vreinterpretq_u32_u64(mask)) != 0u;
}

inline V vlgamma(V z) noexcept {
  const V one = vdupq_n_f64(1.0);
  const M ok = vandq_u64(vcgtq_f64(z, vdupq_n_f64(0.0)),
                         vcltq_f64(z, vdupq_n_f64(kInfinity)));
  const V ze = vbslq_f64(ok, z, one);

  V p = one;
  V y = ze;
  for (int step = 0; step < 16; ++step) {
    const M below = vcltq_f64(y, vdupq_n_f64(kStirlingMin));
    if (!any_true(below)) break;
    p = vmulq_f64(p, vbslq_f64(below, y, one));
    y = vaddq_f64(y, vbslq_f64(below, one, vdupq_n_f64(0.0)));
  }
  V out = vsubq_f64(vlgamma_stirling(y), vlog(p));
  // The two exact zeros, matched to the scalar so a parity failure elsewhere
  // is never confused with this.
  const M zero_pt = vorrq_u64(vceqq_f64(ze, one), vceqq_f64(ze, vdupq_n_f64(2.0)));
  out = vbslq_f64(zero_pt, vdupq_n_f64(0.0), out);
  // z <= 0 and z = +inf are poles or overflow; NaN propagates.
  out = vbslq_f64(ok, out, vdupq_n_f64(kInfinity));
  return vbslq_f64(vceqq_f64(z, z), out, z);
}

inline void lgamma_batch(const double* z, double* out, std::size_t n) noexcept {
  std::size_t i = 0;
  for (; i + 2 <= n; i += 2) vst1q_f64(out + i, vlgamma(vld1q_f64(z + i)));
  for (; i < n; ++i) out[i] = ::pto::peaks::lgamma_approx(z[i]);
}

inline void poisson_batch(const std::uint32_t* ks, const double* lams, double* out,
                          std::size_t n) noexcept {
  const V one = vdupq_n_f64(1.0);
  const V zero = vdupq_n_f64(0.0);
  const V half = vdupq_n_f64(0.5);
  const V tiny = vdupq_n_f64(kCfTiny);
  const V fmaxk = vdupq_n_f64(static_cast<double>(kFiniteMaxK));
  const V fmaxlam = vdupq_n_f64(kFiniteMaxLambda);
  const V fmintail = vdupq_n_f64(kFiniteMinTail);
  std::size_t i = 0;
  for (; i + 2 <= n; i += 2) {
    const V raw_lam = vld1q_f64(lams + i);
    const V raw_k = vcvtq_f64_u64(vmovl_u32(vld1_u32(ks + i)));

    const M valid = vandq_u64(vandq_u64(vcgtq_f64(raw_lam, zero),
                                        vcltq_f64(raw_lam, vdupq_n_f64(kInfinity))),
                              vcgtq_f64(raw_k, zero));
    const V lam = vbslq_f64(valid, raw_lam, one);
    const V a = vbslq_f64(valid, raw_k, one);

    // ---- fast path: the finite lower sum -----------------------------------
    //
    // This loop is the point of the whole kernel. It has a trip count fixed
    // before it starts, no division, no convergence test and no cross-lane
    // operation of any kind -- the lanes stay in the vector registers from the
    // first term to the last. The trip count is one horizontal max, taken once
    // per vector, of the counts that are actually on this path.
    M admit = vandq_u64(valid, vandq_u64(vcleq_f64(a, fmaxk), vcleq_f64(lam, fmaxlam)));
    const int trip = static_cast<int>(hmax(vbslq_f64(admit, a, one)));

    V t = vexp(vnegq_f64(lam));  // pmf(0)
    V qf = t;
    for (int j = 1; j < trip; ++j) {
      t = vmulq_f64(t, vmulq_f64(lam, vdupq_n_f64(kReciprocal.v[j])));
      // Lanes past their own k contribute nothing. Masking the ADDEND rather
      // than the term keeps the recurrence itself branch-free, and the term
      // cannot run away: every value it takes is a Poisson probability.
      qf = vaddq_f64(qf, vreinterpretq_f64_u64(vandq_u64(
                             vreinterpretq_u64_f64(t),
                             vcgtq_f64(a, vdupq_n_f64(static_cast<double>(j))))));
    }
    admit = vandq_u64(admit, vcgeq_f64(vsubq_f64(one, qf), fmintail));
    const V fast_res = vlog1p_full(vnegq_f64(qf));

    V res;
    if (all_true(admit)) {
      res = fast_res;
    } else {
      // ---- exact path: series below lambda = k+1, continued fraction above.
      // Reached only when a lane needs the deep tail or a count past the
      // table; on a background track that is well under one vector in a
      // hundred, and on an enriched one it is the whole vector, which is why
      // it stays vectorised rather than falling back to scalar.
      const V idxd = vsubq_f64(
          vminq_f64(a, vdupq_n_f64(static_cast<double>(kLogFactorialTableMax) + 1.0)),
          one);
      const int64x2_t idx = vcvtq_s64_f64(idxd);
      const V tabv = vsetq_lane_f64(kLogFactorial.v[vgetq_lane_s64(idx, 1)],
                                    vdupq_n_f64(kLogFactorial.v[vgetq_lane_s64(idx, 0)]), 1);
      const M tabled =
          vcleq_f64(a, vdupq_n_f64(static_cast<double>(kLogFactorialTableMax) + 1.0));
      const V lngam = vbslq_f64(tabled, tabv, vlgamma_stirling(a));

      const V lead = vsubq_f64(vfmaq_f64(vnegq_f64(lam), a, vlog(lam)), lngam);
      const M is_series = vcltq_f64(lam, vaddq_f64(a, one));

      V ln_series = lead;
      M sdone = vdupq_n_u64(0xFFFFFFFFFFFFFFFFull);
      if (any_true(is_series)) {
        const V xs = vbslq_f64(is_series, lam, zero);
        V ap = a;
        V del = vdivq_f64(one, a);
        V sum = del;
        sdone = vdupq_n_u64(0);
        // Strided: the horizontal convergence test is what serialises a vector
        // loop, so it runs once per kConvergenceStride terms instead of once
        // per term. Overshooting is harmless -- the terms are still shrinking.
        for (int it = 0; it < kMaxIterations; it += kConvergenceStride) {
          for (int s2 = 0; s2 < kConvergenceStride; ++s2) {
            ap = vaddq_f64(ap, one);
            del = vmulq_f64(del, vdivq_f64(xs, ap));
            sum = vaddq_f64(sum, del);
          }
          sdone = vcleq_f64(del, vmulq_f64(sum, vdupq_n_f64(kSeriesEps)));
          if (all_true(sdone)) break;
        }
        ln_series = vaddq_f64(lead, vlog(sum));
      }

      V ln_cf = zero;
      M cdone = vdupq_n_u64(0xFFFFFFFFFFFFFFFFull);
      if (!all_true(is_series)) {
        const V ac = vbslq_f64(is_series, one, a);
        const V xc = vbslq_f64(is_series, vdupq_n_f64(64.0), lam);
        V b = vaddq_f64(vsubq_f64(xc, ac), one);
        V c = vdupq_n_f64(1.0 / kCfTiny);
        V d = vdivq_f64(one, b);
        V h = d;
        cdone = vdupq_n_u64(0);
        int it = 1;
        while (it <= kMaxIterations) {
          for (int s2 = 0; s2 < kConvergenceStride; ++s2, ++it) {
            const V di = vdupq_n_f64(static_cast<double>(it));
            const V an = vnegq_f64(vmulq_f64(di, vsubq_f64(di, ac)));
            b = vaddq_f64(b, vdupq_n_f64(2.0));
            d = vfmaq_f64(b, an, d);
            d = vbslq_f64(vcltq_f64(vabsq_f64(d), tiny), tiny, d);
            c = vaddq_f64(b, vdivq_f64(an, c));
            c = vbslq_f64(vcltq_f64(vabsq_f64(c), tiny), tiny, c);
            d = vdivq_f64(one, d);
            const V dl = vmulq_f64(d, c);
            // The freeze stays per iteration: the continued fraction's delta
            // oscillates around 1, so a lane that has converged must stop
            // updating h then, not at the next stride boundary.
            h = vbslq_f64(cdone, h, vmulq_f64(h, dl));
            cdone = vorrq_u64(cdone, vcleq_f64(vabsq_f64(vsubq_f64(dl, one)),
                                               vdupq_n_f64(kCfEps)));
          }
          if (all_true(cdone)) break;
        }
        V q = vmulq_f64(vexp(lead), h);
        q = vbslq_f64(vcgtq_f64(q, zero), q, zero);  // NaN and negatives -> 0
        q = vminq_f64(q, half);
        ln_cf = vlog1p_bounded(vnegq_f64(q));
      }

      V exact_res = vbslq_f64(is_series, ln_series, ln_cf);
      exact_res = vbslq_f64(vbslq_u64(is_series, sdone, cdone), exact_res, zero);
      res = vbslq_f64(admit, fast_res, exact_res);
    }

    const M good = vandq_u64(vcltq_f64(res, zero),
                             vcgtq_f64(res, vdupq_n_f64(-kInfinity)));
    res = vbslq_f64(vandq_u64(good, valid), res, zero);
    vst1q_f64(out + i, res);
  }
  for (; i < n; ++i) out[i] = ::pto::peaks::log_poisson_sf(ks[i], lams[i]);
}

}  // namespace neon
}  // namespace detail
#endif  // PEAKS_BUILD_ARM_NEON

#if defined(PEAKS_BUILD_X86_SIMD) && defined(PEAKS_HAVE_AVX2)
namespace detail {
namespace avx2 {

// Per-function, not per-translation-unit. A -mavx2 on the target would let the
// compiler emit AVX2 in code no runtime check guards, and the binary would
// fault on a CPU without it.
#define PEAKS_AVX2_FN __attribute__((target("avx2,fma"))) inline

using V = __m256d;

PEAKS_AVX2_FN V vsel(V m, V t, V f) noexcept { return _mm256_blendv_pd(f, t, m); }
PEAKS_AVX2_FN bool all_true(V m) noexcept { return _mm256_movemask_pd(m) == 0xF; }
PEAKS_AVX2_FN bool any_true(V m) noexcept { return _mm256_movemask_pd(m) != 0; }
PEAKS_AVX2_FN V vabs(V a) noexcept {
  return _mm256_andnot_pd(_mm256_set1_pd(-0.0), a);
}

// Even/odd split Horner; see the scalar atanh2_n for why the chain is halved.
template <int Terms>
__attribute__((target("avx2,fma"))) inline V vatanh2_n(V f) noexcept {
  const V r = _mm256_mul_pd(f, f);
  const V r2 = _mm256_mul_pd(r, r);
  constexpr int last_even = ((Terms - 1) / 2) * 2;
  constexpr int last_odd = ((Terms - 2) / 2) * 2 + 1;
  V even = _mm256_set1_pd(kAtanh.c[last_even]);
  for (int j = last_even - 2; j >= 0; j -= 2) {
    even = _mm256_fmadd_pd(even, r2, _mm256_set1_pd(kAtanh.c[j]));
  }
  V odd = _mm256_set1_pd(kAtanh.c[last_odd]);
  for (int j = last_odd - 2; j >= 1; j -= 2) {
    odd = _mm256_fmadd_pd(odd, r2, _mm256_set1_pd(kAtanh.c[j]));
  }
  const V acc = _mm256_fmadd_pd(odd, r, even);
  const V two_f = _mm256_add_pd(f, f);
  return _mm256_fmadd_pd(_mm256_mul_pd(two_f, r), acc, two_f);
}

PEAKS_AVX2_FN V vatanh2(V f) noexcept { return vatanh2_n<kAtanhTerms>(f); }

PEAKS_AVX2_FN V vlog(V x) noexcept {
  const V sub = _mm256_cmp_pd(x, _mm256_set1_pd(2.2250738585072014e-308), _CMP_LT_OQ);
  const V xs = vsel(sub, _mm256_mul_pd(x, _mm256_set1_pd(0x1p54)), x);
  const __m256i adj = _mm256_castpd_si256(
      _mm256_and_pd(sub, _mm256_castsi256_pd(_mm256_set1_epi64x(-54))));

  const __m256i bits = _mm256_castpd_si256(xs);
  __m256i e = _mm256_add_epi64(_mm256_srli_epi64(bits, 52), _mm256_set1_epi64x(-1022));
  e = _mm256_add_epi64(e, adj);
  V m = _mm256_castsi256_pd(
      _mm256_or_si256(_mm256_and_si256(bits, _mm256_set1_epi64x(0x000FFFFFFFFFFFFFll)),
                      _mm256_set1_epi64x(0x3FE0000000000000ll)));
  const V lo = _mm256_cmp_pd(m, _mm256_set1_pd(kSqrtHalf), _CMP_LT_OQ);
  m = vsel(lo, _mm256_add_pd(m, m), m);
  // A true lane of `lo` is all ones, which as an int64 is -1.
  e = _mm256_add_epi64(e, _mm256_castpd_si256(
                              _mm256_and_pd(lo, _mm256_castsi256_pd(
                                                    _mm256_set1_epi64x(-1)))));

  const V f = _mm256_div_pd(_mm256_sub_pd(m, _mm256_set1_pd(1.0)),
                            _mm256_add_pd(m, _mm256_set1_pd(1.0)));
  // int64 -> double without AVX-512: adding the round-to-integer constant's
  // bit pattern puts the value in the mantissa, and one subtraction extracts it.
  const V de = _mm256_sub_pd(
      _mm256_castsi256_pd(_mm256_add_epi64(
          e, _mm256_set1_epi64x(static_cast<long long>(kShiftBits)))),
      _mm256_set1_pd(kShift));
  return _mm256_fmadd_pd(
      de, _mm256_set1_pd(kLn2Hi),
      _mm256_fmadd_pd(de, _mm256_set1_pd(kLn2Lo), vatanh2_n<kAtanhShortTerms>(f)));
}

PEAKS_AVX2_FN V vpow2(__m256i n) noexcept {
  return _mm256_castsi256_pd(
      _mm256_slli_epi64(_mm256_add_epi64(n, _mm256_set1_epi64x(1023)), 52));
}

PEAKS_AVX2_FN V vexp(V x0) noexcept {
  const V x = _mm256_min_pd(_mm256_max_pd(x0, _mm256_set1_pd(kExpUnderflow)),
                            _mm256_set1_pd(kExpOverflow));
  const V t = _mm256_fmadd_pd(x, _mm256_set1_pd(kInvLn2), _mm256_set1_pd(kShift));
  const V dn = _mm256_sub_pd(t, _mm256_set1_pd(kShift));
  const __m256i n = _mm256_sub_epi64(
      _mm256_castpd_si256(t),
      _mm256_set1_epi64x(static_cast<long long>(kShiftBits)));

  V r = _mm256_fnmadd_pd(dn, _mm256_set1_pd(kLn2Hi), x);
  r = _mm256_fnmadd_pd(dn, _mm256_set1_pd(kLn2Lo), r);

  const V r2 = _mm256_mul_pd(r, r);
  V ev = _mm256_set1_pd(kExp.c[kExpTerms - 2]);
  for (int j = kExpTerms - 4; j >= 0; j -= 2) ev = _mm256_fmadd_pd(ev, r2, _mm256_set1_pd(kExp.c[j]));
  V od = _mm256_set1_pd(kExp.c[kExpTerms - 1]);
  for (int j = kExpTerms - 3; j >= 1; j -= 2) od = _mm256_fmadd_pd(od, r2, _mm256_set1_pd(kExp.c[j]));
  const V p = _mm256_fmadd_pd(od, r, ev);
  const V er = _mm256_add_pd(_mm256_set1_pd(1.0), _mm256_fmadd_pd(r2, p, r));

  // n >> 1, arithmetic. AVX2 has no 64-bit arithmetic right shift, and the
  // bias makes the logical one exact: floor((n + 2048)/2) - 1024 = floor(n/2)
  // for every n the reduction can produce.
  const __m256i u = _mm256_add_epi64(n, _mm256_set1_epi64x(2048));
  const __m256i n1 =
      _mm256_sub_epi64(_mm256_srli_epi64(u, 1), _mm256_set1_epi64x(1024));
  return _mm256_mul_pd(_mm256_mul_pd(er, vpow2(n1)),
                       vpow2(_mm256_sub_epi64(n, n1)));
}

PEAKS_AVX2_FN V vlog1p_bounded(V y) noexcept {
  return vatanh2(_mm256_div_pd(y, _mm256_add_pd(_mm256_set1_pd(2.0), y)));
}

// ln(1+y) over the finite sum's full range in one logarithm; see the NEON
// vlog1p_full for the derivation of the correction term.
PEAKS_AVX2_FN V vlog1p_full(V y) noexcept {
  const V one = _mm256_set1_pd(1.0);
  const V u = _mm256_max_pd(_mm256_add_pd(one, y), _mm256_set1_pd(kCfTiny));
  const V c = _mm256_div_pd(_mm256_sub_pd(y, _mm256_sub_pd(u, one)), u);
  return _mm256_add_pd(vlog(u), c);
}

// Horizontal max of four lanes; runs once per vector to choose a trip count.
PEAKS_AVX2_FN double hmax(V v) noexcept {
  __m128d m = _mm_max_pd(_mm256_castpd256_pd128(v), _mm256_extractf128_pd(v, 1));
  m = _mm_max_pd(m, _mm_unpackhi_pd(m, m));
  return _mm_cvtsd_f64(m);
}

PEAKS_AVX2_FN V vlgamma_stirling(V z) noexcept {
  const V w = _mm256_div_pd(_mm256_set1_pd(1.0), z);
  const V w2 = _mm256_mul_pd(w, w);
  V s = _mm256_set1_pd(kS4);
  s = _mm256_fmadd_pd(s, w2, _mm256_set1_pd(kS3));
  s = _mm256_fmadd_pd(s, w2, _mm256_set1_pd(kS2));
  s = _mm256_fmadd_pd(s, w2, _mm256_set1_pd(kS1));
  s = _mm256_fmadd_pd(s, w2, _mm256_set1_pd(kS0));
  s = _mm256_mul_pd(s, w);
  const V t = _mm256_fmadd_pd(_mm256_sub_pd(z, _mm256_set1_pd(0.5)), vlog(z),
                              _mm256_sub_pd(_mm256_set1_pd(kHalfLogTwoPi), z));
  return _mm256_add_pd(t, s);
}

PEAKS_AVX2_FN V vlgamma(V z) noexcept {
  const V one = _mm256_set1_pd(1.0);
  const V ok = _mm256_and_pd(_mm256_cmp_pd(z, _mm256_setzero_pd(), _CMP_GT_OQ),
                             _mm256_cmp_pd(z, _mm256_set1_pd(kInfinity), _CMP_LT_OQ));
  const V ze = vsel(ok, z, one);

  V p = one;
  V y = ze;
  for (int step = 0; step < 16; ++step) {
    const V below = _mm256_cmp_pd(y, _mm256_set1_pd(kStirlingMin), _CMP_LT_OQ);
    if (!any_true(below)) break;
    p = _mm256_mul_pd(p, vsel(below, y, one));
    y = _mm256_add_pd(y, vsel(below, one, _mm256_setzero_pd()));
  }
  V out = _mm256_sub_pd(vlgamma_stirling(y), vlog(p));
  const V zero_pt = _mm256_or_pd(_mm256_cmp_pd(ze, one, _CMP_EQ_OQ),
                                 _mm256_cmp_pd(ze, _mm256_set1_pd(2.0), _CMP_EQ_OQ));
  out = vsel(zero_pt, _mm256_setzero_pd(), out);
  out = vsel(ok, out, _mm256_set1_pd(kInfinity));
  return vsel(_mm256_cmp_pd(z, z, _CMP_EQ_OQ), out, z);
}

PEAKS_AVX2_FN void lgamma_batch(const double* z, double* out, std::size_t n) noexcept {
  std::size_t i = 0;
  for (; i + 4 <= n; i += 4) _mm256_storeu_pd(out + i, vlgamma(_mm256_loadu_pd(z + i)));
  for (; i < n; ++i) out[i] = ::pto::peaks::lgamma_approx(z[i]);
}

PEAKS_AVX2_FN void poisson_batch(const std::uint32_t* ks, const double* lams,
                                 double* out, std::size_t n) noexcept {
  const V one = _mm256_set1_pd(1.0);
  const V zero = _mm256_setzero_pd();
  const V tiny = _mm256_set1_pd(kCfTiny);
  const V tabmax = _mm256_set1_pd(static_cast<double>(kLogFactorialTableMax) + 1.0);
  const V fmaxk = _mm256_set1_pd(static_cast<double>(kFiniteMaxK));
  const V fmaxlam = _mm256_set1_pd(kFiniteMaxLambda);
  const V fmintail = _mm256_set1_pd(kFiniteMinTail);
  std::size_t i = 0;
  for (; i + 4 <= n; i += 4) {
    const V raw_lam = _mm256_loadu_pd(lams + i);
    // uint32 -> double: OR into the mantissa of 2^52 and subtract it. A signed
    // convert would turn every count past 2^31 into a negative.
    const __m256i k64 = _mm256_cvtepu32_epi64(
        _mm_loadu_si128(reinterpret_cast<const __m128i*>(ks + i)));
    const V raw_k = _mm256_sub_pd(
        _mm256_castsi256_pd(
            _mm256_or_si256(k64, _mm256_set1_epi64x(0x4330000000000000ll))),
        _mm256_set1_pd(4503599627370496.0));

    const V valid = _mm256_and_pd(
        _mm256_and_pd(_mm256_cmp_pd(raw_lam, zero, _CMP_GT_OQ),
                      _mm256_cmp_pd(raw_lam, _mm256_set1_pd(kInfinity), _CMP_LT_OQ)),
        _mm256_cmp_pd(raw_k, zero, _CMP_GT_OQ));
    const V lam = vsel(valid, raw_lam, one);
    const V a = vsel(valid, raw_k, one);

    // ---- fast path: the finite lower sum. Fixed trip count, no division, no
    // convergence test, no cross-lane operation inside the loop.
    V admit = _mm256_and_pd(
        valid, _mm256_and_pd(_mm256_cmp_pd(a, fmaxk, _CMP_LE_OQ),
                             _mm256_cmp_pd(lam, fmaxlam, _CMP_LE_OQ)));
    const int trip = static_cast<int>(hmax(vsel(admit, a, one)));

    V t = vexp(_mm256_sub_pd(zero, lam));  // pmf(0)
    V qf = t;
    for (int j = 1; j < trip; ++j) {
      t = _mm256_mul_pd(t, _mm256_mul_pd(lam, _mm256_set1_pd(kReciprocal.v[j])));
      qf = _mm256_add_pd(
          qf, _mm256_and_pd(t, _mm256_cmp_pd(a, _mm256_set1_pd(static_cast<double>(j)),
                                             _CMP_GT_OQ)));
    }
    admit = _mm256_and_pd(
        admit, _mm256_cmp_pd(_mm256_sub_pd(one, qf), fmintail, _CMP_GE_OQ));
    const V fast_res = vlog1p_full(_mm256_sub_pd(zero, qf));

    V res;
    if (all_true(admit)) {
      res = fast_res;
    } else {
      // ---- exact path: series below lambda = k+1, continued fraction above.
      const __m128i idx =
          _mm256_cvtpd_epi32(_mm256_sub_pd(_mm256_min_pd(a, tabmax), one));
      const V tabv = _mm256_i32gather_pd(kLogFactorial.v, idx, 8);
      const V lngam =
          vsel(_mm256_cmp_pd(a, tabmax, _CMP_LE_OQ), tabv, vlgamma_stirling(a));

      const V lead = _mm256_sub_pd(
          _mm256_fmadd_pd(a, vlog(lam), _mm256_sub_pd(zero, lam)), lngam);
      const V is_series = _mm256_cmp_pd(lam, _mm256_add_pd(a, one), _CMP_LT_OQ);

      const V all_ones = _mm256_castsi256_pd(_mm256_set1_epi64x(-1));
      V ln_series = lead;
      V sdone = all_ones;
      if (any_true(is_series)) {
        const V xs = _mm256_and_pd(is_series, lam);  // inactive lanes get x = 0
        V ap = a;
        V del = _mm256_div_pd(one, a);
        V sum = del;
        sdone = zero;
        // Strided: the horizontal convergence test runs once per
        // kConvergenceStride terms rather than once per term.
        for (int it = 0; it < kMaxIterations; it += kConvergenceStride) {
          for (int s2 = 0; s2 < kConvergenceStride; ++s2) {
            ap = _mm256_add_pd(ap, one);
            del = _mm256_mul_pd(del, _mm256_div_pd(xs, ap));
            sum = _mm256_add_pd(sum, del);
          }
          sdone = _mm256_cmp_pd(del, _mm256_mul_pd(sum, _mm256_set1_pd(kSeriesEps)),
                                _CMP_LE_OQ);
          if (all_true(sdone)) break;
        }
        ln_series = _mm256_add_pd(lead, vlog(sum));
      }

      V ln_cf = zero;
      V cdone = all_ones;
      if (!all_true(is_series)) {
        const V ac = vsel(is_series, one, a);
        const V xc = vsel(is_series, _mm256_set1_pd(64.0), lam);
        V b = _mm256_add_pd(_mm256_sub_pd(xc, ac), one);
        V c = _mm256_set1_pd(1.0 / kCfTiny);
        V d = _mm256_div_pd(one, b);
        V h = d;
        cdone = zero;
        int it = 1;
        while (it <= kMaxIterations) {
          for (int s2 = 0; s2 < kConvergenceStride; ++s2, ++it) {
            const V di = _mm256_set1_pd(static_cast<double>(it));
            const V an = _mm256_sub_pd(zero, _mm256_mul_pd(di, _mm256_sub_pd(di, ac)));
            b = _mm256_add_pd(b, _mm256_set1_pd(2.0));
            d = _mm256_fmadd_pd(an, d, b);
            d = vsel(_mm256_cmp_pd(vabs(d), tiny, _CMP_LT_OQ), tiny, d);
            c = _mm256_add_pd(b, _mm256_div_pd(an, c));
            c = vsel(_mm256_cmp_pd(vabs(c), tiny, _CMP_LT_OQ), tiny, c);
            d = _mm256_div_pd(one, d);
            const V dl = _mm256_mul_pd(d, c);
            // The freeze stays per iteration: the delta oscillates around 1.
            h = vsel(cdone, h, _mm256_mul_pd(h, dl));
            cdone = _mm256_or_pd(cdone, _mm256_cmp_pd(vabs(_mm256_sub_pd(dl, one)),
                                                      _mm256_set1_pd(kCfEps),
                                                      _CMP_LE_OQ));
          }
          if (all_true(cdone)) break;
        }
        V q = _mm256_mul_pd(vexp(lead), h);
        q = vsel(_mm256_cmp_pd(q, zero, _CMP_GT_OQ), q, zero);
        q = _mm256_min_pd(q, _mm256_set1_pd(0.5));
        ln_cf = vlog1p_bounded(_mm256_sub_pd(zero, q));
      }

      V exact_res = vsel(is_series, ln_series, ln_cf);
      exact_res = vsel(vsel(is_series, sdone, cdone), exact_res, zero);
      res = vsel(admit, fast_res, exact_res);
    }

    const V good = _mm256_and_pd(
        _mm256_cmp_pd(res, zero, _CMP_LT_OQ),
        _mm256_cmp_pd(res, _mm256_set1_pd(-kInfinity), _CMP_GT_OQ));
    res = vsel(_mm256_and_pd(good, valid), res, zero);
    _mm256_storeu_pd(out + i, res);
  }
  for (; i < n; ++i) out[i] = ::pto::peaks::log_poisson_sf(ks[i], lams[i]);
}

#undef PEAKS_AVX2_FN

}  // namespace avx2
}  // namespace detail
#endif  // PEAKS_HAVE_AVX2

#if defined(PEAKS_BUILD_X86_SIMD) && defined(PEAKS_HAVE_AVX512)
namespace detail {
namespace avx512 {

#define PEAKS_AVX512_FN __attribute__((target("avx512f,avx512dq"))) inline

using V = __m512d;
using M = __mmask8;

PEAKS_AVX512_FN V vsel(M m, V t, V f) noexcept { return _mm512_mask_blend_pd(m, f, t); }
PEAKS_AVX512_FN bool all_true(M m) noexcept { return m == 0xFF; }
PEAKS_AVX512_FN bool any_true(M m) noexcept { return m != 0; }
PEAKS_AVX512_FN V vabs(V a) noexcept {
  return _mm512_castsi512_pd(_mm512_and_si512(
      _mm512_castpd_si512(a), _mm512_set1_epi64(0x7FFFFFFFFFFFFFFFll)));
}

// Even/odd split Horner; see the scalar atanh2_n for why the chain is halved.
template <int Terms>
__attribute__((target("avx512f,avx512dq"))) inline V vatanh2_n(V f) noexcept {
  const V r = _mm512_mul_pd(f, f);
  const V r2 = _mm512_mul_pd(r, r);
  constexpr int last_even = ((Terms - 1) / 2) * 2;
  constexpr int last_odd = ((Terms - 2) / 2) * 2 + 1;
  V even = _mm512_set1_pd(kAtanh.c[last_even]);
  for (int j = last_even - 2; j >= 0; j -= 2) {
    even = _mm512_fmadd_pd(even, r2, _mm512_set1_pd(kAtanh.c[j]));
  }
  V odd = _mm512_set1_pd(kAtanh.c[last_odd]);
  for (int j = last_odd - 2; j >= 1; j -= 2) {
    odd = _mm512_fmadd_pd(odd, r2, _mm512_set1_pd(kAtanh.c[j]));
  }
  const V acc = _mm512_fmadd_pd(odd, r, even);
  const V two_f = _mm512_add_pd(f, f);
  return _mm512_fmadd_pd(_mm512_mul_pd(two_f, r), acc, two_f);
}

PEAKS_AVX512_FN V vatanh2(V f) noexcept { return vatanh2_n<kAtanhTerms>(f); }

PEAKS_AVX512_FN V vlog(V x) noexcept {
  const M sub = _mm512_cmp_pd_mask(x, _mm512_set1_pd(2.2250738585072014e-308),
                                   _CMP_LT_OQ);
  const V xs = vsel(sub, _mm512_mul_pd(x, _mm512_set1_pd(0x1p54)), x);
  const __m512i bits = _mm512_castpd_si512(xs);
  __m512i e = _mm512_add_epi64(_mm512_srli_epi64(bits, 52), _mm512_set1_epi64(-1022));
  e = _mm512_mask_add_epi64(e, sub, e, _mm512_set1_epi64(-54));

  V m = _mm512_castsi512_pd(
      _mm512_or_si512(_mm512_and_si512(bits, _mm512_set1_epi64(0x000FFFFFFFFFFFFFll)),
                      _mm512_set1_epi64(0x3FE0000000000000ll)));
  const M lo = _mm512_cmp_pd_mask(m, _mm512_set1_pd(kSqrtHalf), _CMP_LT_OQ);
  m = vsel(lo, _mm512_add_pd(m, m), m);
  e = _mm512_mask_add_epi64(e, lo, e, _mm512_set1_epi64(-1));

  const V f = _mm512_div_pd(_mm512_sub_pd(m, _mm512_set1_pd(1.0)),
                            _mm512_add_pd(m, _mm512_set1_pd(1.0)));
  const V de = _mm512_sub_pd(
      _mm512_castsi512_pd(_mm512_add_epi64(
          e, _mm512_set1_epi64(static_cast<long long>(kShiftBits)))),
      _mm512_set1_pd(kShift));
  return _mm512_fmadd_pd(
      de, _mm512_set1_pd(kLn2Hi),
      _mm512_fmadd_pd(de, _mm512_set1_pd(kLn2Lo), vatanh2_n<kAtanhShortTerms>(f)));
}

PEAKS_AVX512_FN V vpow2(__m512i n) noexcept {
  return _mm512_castsi512_pd(
      _mm512_slli_epi64(_mm512_add_epi64(n, _mm512_set1_epi64(1023)), 52));
}

PEAKS_AVX512_FN V vexp(V x0) noexcept {
  const V x = _mm512_min_pd(_mm512_max_pd(x0, _mm512_set1_pd(kExpUnderflow)),
                            _mm512_set1_pd(kExpOverflow));
  const V t = _mm512_fmadd_pd(x, _mm512_set1_pd(kInvLn2), _mm512_set1_pd(kShift));
  const V dn = _mm512_sub_pd(t, _mm512_set1_pd(kShift));
  const __m512i n = _mm512_sub_epi64(
      _mm512_castpd_si512(t), _mm512_set1_epi64(static_cast<long long>(kShiftBits)));

  V r = _mm512_fnmadd_pd(dn, _mm512_set1_pd(kLn2Hi), x);
  r = _mm512_fnmadd_pd(dn, _mm512_set1_pd(kLn2Lo), r);

  const V r2 = _mm512_mul_pd(r, r);
  V ev = _mm512_set1_pd(kExp.c[kExpTerms - 2]);
  for (int j = kExpTerms - 4; j >= 0; j -= 2) ev = _mm512_fmadd_pd(ev, r2, _mm512_set1_pd(kExp.c[j]));
  V od = _mm512_set1_pd(kExp.c[kExpTerms - 1]);
  for (int j = kExpTerms - 3; j >= 1; j -= 2) od = _mm512_fmadd_pd(od, r2, _mm512_set1_pd(kExp.c[j]));
  const V p = _mm512_fmadd_pd(od, r, ev);
  const V er = _mm512_add_pd(_mm512_set1_pd(1.0), _mm512_fmadd_pd(r2, p, r));

  const __m512i n1 = _mm512_srai_epi64(n, 1);
  return _mm512_mul_pd(_mm512_mul_pd(er, vpow2(n1)),
                       vpow2(_mm512_sub_epi64(n, n1)));
}

PEAKS_AVX512_FN V vlog1p_bounded(V y) noexcept {
  return vatanh2(_mm512_div_pd(y, _mm512_add_pd(_mm512_set1_pd(2.0), y)));
}

// ln(1+y) over the finite sum's full range in one logarithm; see the NEON
// vlog1p_full for the derivation of the correction term.
PEAKS_AVX512_FN V vlog1p_full(V y) noexcept {
  const V one = _mm512_set1_pd(1.0);
  const V u = _mm512_max_pd(_mm512_add_pd(one, y), _mm512_set1_pd(kCfTiny));
  const V c = _mm512_div_pd(_mm512_sub_pd(y, _mm512_sub_pd(u, one)), u);
  return _mm512_add_pd(vlog(u), c);
}

PEAKS_AVX512_FN double hmax(V v) noexcept { return _mm512_reduce_max_pd(v); }

PEAKS_AVX512_FN V vlgamma_stirling(V z) noexcept {
  const V w = _mm512_div_pd(_mm512_set1_pd(1.0), z);
  const V w2 = _mm512_mul_pd(w, w);
  V s = _mm512_set1_pd(kS4);
  s = _mm512_fmadd_pd(s, w2, _mm512_set1_pd(kS3));
  s = _mm512_fmadd_pd(s, w2, _mm512_set1_pd(kS2));
  s = _mm512_fmadd_pd(s, w2, _mm512_set1_pd(kS1));
  s = _mm512_fmadd_pd(s, w2, _mm512_set1_pd(kS0));
  s = _mm512_mul_pd(s, w);
  const V t = _mm512_fmadd_pd(_mm512_sub_pd(z, _mm512_set1_pd(0.5)), vlog(z),
                              _mm512_sub_pd(_mm512_set1_pd(kHalfLogTwoPi), z));
  return _mm512_add_pd(t, s);
}

PEAKS_AVX512_FN V vlgamma(V z) noexcept {
  const V one = _mm512_set1_pd(1.0);
  const M ok = static_cast<M>(
      _mm512_cmp_pd_mask(z, _mm512_setzero_pd(), _CMP_GT_OQ) &
      _mm512_cmp_pd_mask(z, _mm512_set1_pd(kInfinity), _CMP_LT_OQ));
  const V ze = vsel(ok, z, one);

  V p = one;
  V y = ze;
  for (int step = 0; step < 16; ++step) {
    const M below = _mm512_cmp_pd_mask(y, _mm512_set1_pd(kStirlingMin), _CMP_LT_OQ);
    if (!any_true(below)) break;
    p = _mm512_mul_pd(p, vsel(below, y, one));
    y = _mm512_add_pd(y, vsel(below, one, _mm512_setzero_pd()));
  }
  V out = _mm512_sub_pd(vlgamma_stirling(y), vlog(p));
  const M zero_pt = static_cast<M>(_mm512_cmp_pd_mask(ze, one, _CMP_EQ_OQ) |
                                   _mm512_cmp_pd_mask(ze, _mm512_set1_pd(2.0),
                                                      _CMP_EQ_OQ));
  out = vsel(zero_pt, _mm512_setzero_pd(), out);
  out = vsel(ok, out, _mm512_set1_pd(kInfinity));
  return vsel(_mm512_cmp_pd_mask(z, z, _CMP_EQ_OQ), out, z);
}

PEAKS_AVX512_FN void lgamma_batch(const double* z, double* out, std::size_t n) noexcept {
  std::size_t i = 0;
  for (; i + 8 <= n; i += 8) _mm512_storeu_pd(out + i, vlgamma(_mm512_loadu_pd(z + i)));
  for (; i < n; ++i) out[i] = ::pto::peaks::lgamma_approx(z[i]);
}

PEAKS_AVX512_FN void poisson_batch(const std::uint32_t* ks, const double* lams,
                                   double* out, std::size_t n) noexcept {
  const V one = _mm512_set1_pd(1.0);
  const V zero = _mm512_setzero_pd();
  const V tiny = _mm512_set1_pd(kCfTiny);
  const V tabmax = _mm512_set1_pd(static_cast<double>(kLogFactorialTableMax) + 1.0);
  const V fmaxk = _mm512_set1_pd(static_cast<double>(kFiniteMaxK));
  const V fmaxlam = _mm512_set1_pd(kFiniteMaxLambda);
  const V fmintail = _mm512_set1_pd(kFiniteMinTail);
  std::size_t i = 0;
  for (; i + 8 <= n; i += 8) {
    const V raw_lam = _mm512_loadu_pd(lams + i);
    const __m512i k64 = _mm512_cvtepu32_epi64(
        _mm256_loadu_si256(reinterpret_cast<const __m256i*>(ks + i)));
    const V raw_k = _mm512_sub_pd(
        _mm512_castsi512_pd(
            _mm512_or_si512(k64, _mm512_set1_epi64(0x4330000000000000ll))),
        _mm512_set1_pd(4503599627370496.0));

    const M valid = static_cast<M>(
        _mm512_cmp_pd_mask(raw_lam, zero, _CMP_GT_OQ) &
        _mm512_cmp_pd_mask(raw_lam, _mm512_set1_pd(kInfinity), _CMP_LT_OQ) &
        _mm512_cmp_pd_mask(raw_k, zero, _CMP_GT_OQ));
    const V lam = vsel(valid, raw_lam, one);
    const V a = vsel(valid, raw_k, one);

    // ---- fast path: the finite lower sum. Fixed trip count, no division, no
    // convergence test, no cross-lane operation inside the loop.
    M admit = static_cast<M>(valid & _mm512_cmp_pd_mask(a, fmaxk, _CMP_LE_OQ) &
                             _mm512_cmp_pd_mask(lam, fmaxlam, _CMP_LE_OQ));
    const int trip = static_cast<int>(hmax(vsel(admit, a, one)));

    V t = vexp(_mm512_sub_pd(zero, lam));  // pmf(0)
    V qf = t;
    for (int j = 1; j < trip; ++j) {
      t = _mm512_mul_pd(t, _mm512_mul_pd(lam, _mm512_set1_pd(kReciprocal.v[j])));
      qf = _mm512_mask_add_pd(
          qf,
          _mm512_cmp_pd_mask(a, _mm512_set1_pd(static_cast<double>(j)), _CMP_GT_OQ),
          qf, t);
    }
    admit = static_cast<M>(
        admit & _mm512_cmp_pd_mask(_mm512_sub_pd(one, qf), fmintail, _CMP_GE_OQ));
    const V fast_res = vlog1p_full(_mm512_sub_pd(zero, qf));

    V res;
    if (all_true(admit)) {
      res = fast_res;
    } else {
      // ---- exact path: series below lambda = k+1, continued fraction above.
      const __m256i idx =
          _mm512_cvtpd_epi32(_mm512_sub_pd(_mm512_min_pd(a, tabmax), one));
      const V tabv = _mm512_i32gather_pd(idx, kLogFactorial.v, 8);
      const V lngam = vsel(_mm512_cmp_pd_mask(a, tabmax, _CMP_LE_OQ), tabv,
                           vlgamma_stirling(a));

      const V lead = _mm512_sub_pd(
          _mm512_fmadd_pd(a, vlog(lam), _mm512_sub_pd(zero, lam)), lngam);
      const M is_series = _mm512_cmp_pd_mask(lam, _mm512_add_pd(a, one), _CMP_LT_OQ);

      V ln_series = lead;
      M sdone = 0xFF;
      if (any_true(is_series)) {
        const V xs = _mm512_maskz_mov_pd(is_series, lam);
        V ap = a;
        V del = _mm512_div_pd(one, a);
        V sum = del;
        sdone = 0;
        for (int it = 0; it < kMaxIterations; it += kConvergenceStride) {
          for (int s2 = 0; s2 < kConvergenceStride; ++s2) {
            ap = _mm512_add_pd(ap, one);
            del = _mm512_mul_pd(del, _mm512_div_pd(xs, ap));
            sum = _mm512_add_pd(sum, del);
          }
          sdone = _mm512_cmp_pd_mask(del, _mm512_mul_pd(sum, _mm512_set1_pd(kSeriesEps)),
                                     _CMP_LE_OQ);
          if (all_true(sdone)) break;
        }
        ln_series = _mm512_add_pd(lead, vlog(sum));
      }

      V ln_cf = zero;
      M cdone = 0xFF;
      if (!all_true(is_series)) {
        const V ac = vsel(is_series, one, a);
        const V xc = vsel(is_series, _mm512_set1_pd(64.0), lam);
        V b = _mm512_add_pd(_mm512_sub_pd(xc, ac), one);
        V c = _mm512_set1_pd(1.0 / kCfTiny);
        V d = _mm512_div_pd(one, b);
        V h = d;
        cdone = 0;
        int it = 1;
        while (it <= kMaxIterations) {
          for (int s2 = 0; s2 < kConvergenceStride; ++s2, ++it) {
            const V di = _mm512_set1_pd(static_cast<double>(it));
            const V an = _mm512_sub_pd(zero, _mm512_mul_pd(di, _mm512_sub_pd(di, ac)));
            b = _mm512_add_pd(b, _mm512_set1_pd(2.0));
            d = _mm512_fmadd_pd(an, d, b);
            d = vsel(_mm512_cmp_pd_mask(vabs(d), tiny, _CMP_LT_OQ), tiny, d);
            c = _mm512_add_pd(b, _mm512_div_pd(an, c));
            c = vsel(_mm512_cmp_pd_mask(vabs(c), tiny, _CMP_LT_OQ), tiny, c);
            d = _mm512_div_pd(one, d);
            const V dl = _mm512_mul_pd(d, c);
            h = vsel(cdone, h, _mm512_mul_pd(h, dl));
            cdone = static_cast<M>(
                cdone | _mm512_cmp_pd_mask(vabs(_mm512_sub_pd(dl, one)),
                                           _mm512_set1_pd(kCfEps), _CMP_LE_OQ));
          }
          if (all_true(cdone)) break;
        }
        V q = _mm512_mul_pd(vexp(lead), h);
        q = vsel(_mm512_cmp_pd_mask(q, zero, _CMP_GT_OQ), q, zero);
        q = _mm512_min_pd(q, _mm512_set1_pd(0.5));
        ln_cf = vlog1p_bounded(_mm512_sub_pd(zero, q));
      }

      V exact_res = vsel(is_series, ln_series, ln_cf);
      const M converged = static_cast<M>((is_series & sdone) |
                                         (static_cast<M>(~is_series) & cdone));
      exact_res = vsel(converged, exact_res, zero);
      res = vsel(admit, fast_res, exact_res);
    }

    const M good = static_cast<M>(
        _mm512_cmp_pd_mask(res, zero, _CMP_LT_OQ) &
        _mm512_cmp_pd_mask(res, _mm512_set1_pd(-kInfinity), _CMP_GT_OQ) & valid);
    res = vsel(good, res, zero);
    _mm512_storeu_pd(out + i, res);
  }
  for (; i < n; ++i) out[i] = ::pto::peaks::log_poisson_sf(ks[i], lams[i]);
}

#undef PEAKS_AVX512_FN

}  // namespace avx512
}  // namespace detail
#endif  // PEAKS_HAVE_AVX512

// ===========================================================================
// Runtime dispatch
//
// Resolved once, lazily, on first use. Compile-time gating is not enough on
// its own: the build machine's capabilities are not the run machine's.
// ===========================================================================

enum class IsaLevel { kScalar, kNeon, kAvx2, kAvx512 };

// Whether a path was compiled into THIS binary. The parity suite iterates over
// every level and skips the ones this returns false for, so a test run says
// what it actually exercised rather than what the enum lists.
[[nodiscard]] inline bool isa_compiled(IsaLevel level) noexcept {
  switch (level) {
    case IsaLevel::kScalar:
      return true;
    case IsaLevel::kNeon:
#if defined(PEAKS_BUILD_ARM_NEON)
      return true;
#else
      return false;
#endif
    case IsaLevel::kAvx2:
#if defined(PEAKS_BUILD_X86_SIMD) && defined(PEAKS_HAVE_AVX2)
      return true;
#else
      return false;
#endif
    case IsaLevel::kAvx512:
#if defined(PEAKS_BUILD_X86_SIMD) && defined(PEAKS_HAVE_AVX512)
      return true;
#else
      return false;
#endif
  }
  return false;
}

[[nodiscard]] inline std::size_t isa_lanes(IsaLevel level) noexcept {
  switch (level) {
    case IsaLevel::kAvx512: return 8;
    case IsaLevel::kAvx2:   return 4;
    case IsaLevel::kNeon:   return 2;
    case IsaLevel::kScalar: break;
  }
  return 1;
}

[[nodiscard]] inline const char* isa_name(IsaLevel level) noexcept {
  switch (level) {
    case IsaLevel::kAvx512: return "avx512";
    case IsaLevel::kAvx2:   return "avx2";
    case IsaLevel::kNeon:   return "neon";
    case IsaLevel::kScalar: break;
  }
  return "scalar";
}

namespace detail {

[[nodiscard]] inline bool env_is(const char* name, const char* value) noexcept {
  const char* v = std::getenv(name);
  if (v == nullptr) return false;
  const char* p = value;
  while (*v != '\0' && *p != '\0' && *v == *p) {
    ++v;
    ++p;
  }
  return *v == '\0' && *p == '\0';
}

// PTO_PEAKS_FORCE_ISA=scalar|neon|avx2|avx512 pins the choice so one machine
// can run every compiled path, which is what a CI sweep needs. Requesting a
// path that was not compiled in is ignored; requesting one the hardware cannot
// execute will fault, and that is the caller's explicit choice -- which is why
// it is opt-in and never a default.
[[nodiscard]] inline IsaLevel detect_isa() noexcept {
  if (env_is("PTO_PEAKS_FORCE_ISA", "scalar")) return IsaLevel::kScalar;
#if defined(PEAKS_BUILD_ARM_NEON)
  if (env_is("PTO_PEAKS_FORCE_ISA", "neon")) return IsaLevel::kNeon;
#endif
#if defined(PEAKS_BUILD_X86_SIMD) && defined(PEAKS_HAVE_AVX512)
  if (env_is("PTO_PEAKS_FORCE_ISA", "avx512")) return IsaLevel::kAvx512;
#endif
#if defined(PEAKS_BUILD_X86_SIMD) && defined(PEAKS_HAVE_AVX2)
  if (env_is("PTO_PEAKS_FORCE_ISA", "avx2")) return IsaLevel::kAvx2;
#endif

#if defined(PEAKS_BUILD_X86_SIMD)
#if defined(PEAKS_HAVE_AVX512)
  if (__builtin_cpu_supports("avx512f") && __builtin_cpu_supports("avx512dq")) {
    return IsaLevel::kAvx512;
  }
#endif
#if defined(PEAKS_HAVE_AVX2)
  if (__builtin_cpu_supports("avx2") && __builtin_cpu_supports("fma")) {
    return IsaLevel::kAvx2;
  }
#endif
#endif

  // NEON is architecturally guaranteed on AArch64: if it compiled in, the CPU
  // has it. There is no probe and nothing to fall back from.
#if defined(PEAKS_BUILD_ARM_NEON)
  return IsaLevel::kNeon;
#else
  return IsaLevel::kScalar;
#endif
}

}  // namespace detail

// Thread-safe through the magic-static rule, and evaluated after static
// initialisation, so there is no __builtin_cpu_init ordering question.
[[nodiscard]] inline IsaLevel isa_level() noexcept {
  static const IsaLevel level = detail::detect_isa();
  return level;
}

// ===========================================================================
// Batch entry points
//
// `out` must not alias `k` or `lambda`. The *_isa forms exist so a test can
// drive a specific kernel rather than whatever the host CPU dispatches to;
// asking for a level this binary does not carry runs the scalar path, which is
// the same answer, only slower.
// ===========================================================================

inline void lgamma_batch_isa(IsaLevel isa, const double* z, double* out,
                             std::size_t n) noexcept {
  switch (isa) {
#if defined(PEAKS_BUILD_ARM_NEON)
    case IsaLevel::kNeon: detail::neon::lgamma_batch(z, out, n); return;
#endif
#if defined(PEAKS_BUILD_X86_SIMD) && defined(PEAKS_HAVE_AVX2)
    case IsaLevel::kAvx2: detail::avx2::lgamma_batch(z, out, n); return;
#endif
#if defined(PEAKS_BUILD_X86_SIMD) && defined(PEAKS_HAVE_AVX512)
    case IsaLevel::kAvx512: detail::avx512::lgamma_batch(z, out, n); return;
#endif
    default: break;
  }
  for (std::size_t i = 0; i < n; ++i) out[i] = lgamma_approx(z[i]);
}

inline void lgamma_batch(const double* z, double* out, std::size_t n) noexcept {
  lgamma_batch_isa(isa_level(), z, out, n);
}

inline void log_poisson_sf_batch_isa(IsaLevel isa, const std::uint32_t* k,
                                     const double* lambda, double* out,
                                     std::size_t n) noexcept {
  switch (isa) {
#if defined(PEAKS_BUILD_ARM_NEON)
    case IsaLevel::kNeon: detail::neon::poisson_batch(k, lambda, out, n); return;
#endif
#if defined(PEAKS_BUILD_X86_SIMD) && defined(PEAKS_HAVE_AVX2)
    case IsaLevel::kAvx2: detail::avx2::poisson_batch(k, lambda, out, n); return;
#endif
#if defined(PEAKS_BUILD_X86_SIMD) && defined(PEAKS_HAVE_AVX512)
    case IsaLevel::kAvx512: detail::avx512::poisson_batch(k, lambda, out, n); return;
#endif
    default: break;
  }
  // The scalar batch is the elementwise scalar function, and the test suite
  // pins that identity: it is what makes this a usable reference for the
  // vector kernels rather than a fourth implementation to keep in step.
  for (std::size_t i = 0; i < n; ++i) out[i] = log_poisson_sf(k[i], lambda[i]);
}

inline void log_poisson_sf_batch(const std::uint32_t* k, const double* lambda,
                                 double* out, std::size_t n) noexcept {
  log_poisson_sf_batch_isa(isa_level(), k, lambda, out, n);
}

// -log10 p for a whole batch, in place over the natural logarithms so there is
// one pass and no temporary.
inline void neg_log10_poisson_sf_batch(const std::uint32_t* k, const double* lambda,
                                       double* out, std::size_t n) noexcept {
  log_poisson_sf_batch(k, lambda, out, n);
  for (std::size_t i = 0; i < n; ++i) out[i] = -out[i] * detail::kInvLn10;
}

}  // namespace pto::peaks
