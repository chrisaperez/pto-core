#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <span>
#include <type_traits>
#include <vector>

// ---------------------------------------------------------------------------
// ISA gating
//
// x86 SIMD is compiled in only when (a) the target really is x86 and (b) the
// build system verified the compiler can emit that ISA through a function-level
// `target` attribute (SCRNA_HAVE_AVX2 / SCRNA_HAVE_AVX512, set by CMake). On a
// non-x86 target -- notably arm64 / Apple Silicon -- every AVX path below
// compiles out and the scalar path is the only one built.
//
// The ISA-specific kernels carry their own `target` attributes instead of
// relying on -mavx2 / -mavx512f being applied to the whole translation unit.
// That is what makes the runtime dispatch below actually safe: the rest of the
// binary stays baseline x86-64, so a CPU without AVX2 cannot execute an AVX2
// instruction merely because some unrelated loop got auto-vectorized. Applying
// the ISA flag target-wide (the previous arrangement) meant a binary built on
// an AVX2 host would SIGILL on a host without it.
// ---------------------------------------------------------------------------
#if defined(__x86_64__) || defined(__i386__) || defined(_M_X64) || defined(_M_IX86)
#define SCRNA_ARCH_X86 1
#endif

#if defined(SCRNA_ARCH_X86) && (defined(SCRNA_HAVE_AVX2) || defined(SCRNA_HAVE_AVX512))
#include <immintrin.h>
#define SCRNA_BUILD_X86_SIMD 1
#endif

// ---------------------------------------------------------------------------
// ARM64 NEON
//
// Gated differently from the x86 paths, and the difference is not an oversight:
//
//   * NEON (Advanced SIMD) is MANDATORY in ARMv8-A, so unlike AVX2 there is no
//     such thing as an AArch64 CPU without it. There is therefore no CPUID
//     check, no `target` attribute, and none of the SIGILL hazard the x86
//     comment above describes -- the ISA-flag-target-wide problem cannot arise
//     because no flag is needed.
//   * `vaddvq_f32` (horizontal add) is AArch64-only, not ARMv7 NEON, so this
//     keys off `__aarch64__` rather than `__ARM_NEON` alone.
//
// It still participates in the runtime dispatch enum. Not for capability
// detection, but so `SCRNA_FORCE_ISA=scalar` can reach the scalar reference at
// run time -- which is what lets one machine test both paths, and what the
// parity tests rely on.
#if defined(__aarch64__) || defined(_M_ARM64)
#define SCRNA_ARCH_ARM64 1
#endif

#if defined(SCRNA_ARCH_ARM64) && defined(__ARM_NEON)
#include <arm_neon.h>
#define SCRNA_BUILD_ARM_NEON 1
#endif

namespace scrna::simd {

// ---------------------------------------------------------------------------
// Wide accumulation: the module's ONE way to write a floating-point reduction
//
// Every reduction in scrna_matrix accumulates in `WideAcc<T>` and narrows
// exactly once, through `checked_narrow`. This is not a style preference. Hand
// -rolled `float` accumulation silently saturating to +/-inf is the single most
// repeated defect class in this module's audit trail -- found four separate
// times, by four different techniques (differential testing against a double
// reference, boundary review, a defence-in-depth read, and a regression test
// written for an unrelated finding), in `sparse_norm_sq`, in three HNSW norm
// loops, in `cosine_from_dot`'s narrowing, and in the dot kernels themselves.
//
// Four independent discoveries of one root cause is a fact about the
// abstraction, not about any one call site. So the primitives come first, the
// reductions are written in terms of them, and `tests/test_wide_accumulation.cpp`
// greps the module for a bare `float` accumulator so a fifth instance fails the
// build rather than waiting for a fifth audit.
//
// The failure mode is worth stating precisely, because "loses precision" is the
// wrong intuition. `float` saturation does not degrade an answer, it inverts
// one: two IDENTICAL rows scaled to 1e30 produce a dot product of +inf, the
// downstream `!isfinite -> 0` guard turns that into a similarity of ZERO, and
// the row that should be every other row's nearest neighbour matches nothing at
// all. Finite input, no diagnostic, silently wrong index.
// ---------------------------------------------------------------------------

// Widest sensible accumulator for a value type: float squares in double.
template <typename T>
using WideAcc = std::conditional_t<std::is_same_v<T, float>, double, T>;

// The one narrowing step, and the only place a wide accumulator becomes a T.
//
// Checked AFTER the cast, never only before it: `WideAcc<float>` is double, so
// a value can be perfectly finite as a double and still overflow float on the
// way down. That exact mistake is finding T5 -- `cosine_from_dot` guarded the
// wide quotient, cast it, and returned +inf.
//
// Non-finite becomes T{0}. That is a deliberate, documented degradation, and it
// is only correct because it is the LAST step: a quantity that cannot be
// represented is reported as absent, never as a number. It is emphatically not
// a licence to narrow early -- see `cosine_from_dot_wide`, which divides in the
// wide type precisely so that a 1e60 dot over a 1e60 denominator narrows to
// 1.0 instead of narrowing to zero twice and reporting no similarity.
template <typename T>
[[nodiscard]] inline T checked_narrow(WideAcc<T> wide) noexcept {
  const T narrowed = static_cast<T>(wide);
  return std::isfinite(narrowed) ? narrowed : T{0};
}

// ---------------------------------------------------------------------------
// Scalar reference kernels
//
// These are always compiled, on every target, and define the correctness
// contract that the vector kernels must match to within floating-point
// tolerance. Nothing here is ISA-gated.
// ---------------------------------------------------------------------------

// Sparse dot product, accumulated in the wide type. This is the primitive; the
// narrowing `sparse_dot` below is a convenience over it.
//
// BOTH the multiply and the accumulation are widened, and the multiply is the
// half that is easy to miss. Widening only the accumulator fixes nothing here:
// two finite floats of 1e30 multiply to 1e60, which is +inf BEFORE the addition
// ever runs, so `double acc += a[i] * b[i]` accumulates infinities in double
// just as faithfully as float did. The cast has to be on the operands.
//
// Callers that go on to divide by a norm should stay wide -- see
// `cosine_from_dot_wide`. Narrowing here and dividing there throws away the
// headroom this function exists to provide.
template <typename T, typename Index>
[[nodiscard]] inline WideAcc<T> sparse_dot_wide(std::span<const Index> idx_a,
                                                std::span<const T> val_a,
                                                std::span<const Index> idx_b,
                                                std::span<const T> val_b) {
  assert(idx_a.size() == val_a.size() && "sparse_dot_wide: idx_a/val_a length mismatch");
  assert(idx_b.size() == val_b.size() && "sparse_dot_wide: idx_b/val_b length mismatch");

  WideAcc<T> acc{0};
  std::size_t i = 0, j = 0;
  while (i < idx_a.size() && j < idx_b.size()) {
    const Index a = idx_a[i];
    const Index b = idx_b[j];
    if (a == b) {
      acc += static_cast<WideAcc<T>>(val_a[i]) * static_cast<WideAcc<T>>(val_b[j]);
      ++i;
      ++j;
    } else if (a < b) {
      ++i;
    } else {
      ++j;
    }
  }
  return acc;
}

// Sparse dot product between two sorted-index sparse rows via a merge-join on
// column indices (the standard approach for unaligned sparsity patterns where a
// dense gather isn't profitable). Runs at O(nnz_a + nnz_b).
//
// Narrows through `checked_narrow`, so a dot product too large to be a float is
// reported as 0 rather than as +inf. If you are about to divide the result by
// something equally large, call `sparse_dot_wide` instead and keep it wide --
// that is what the k-NN sweep does.
//
// Precondition: idx/val spans are pairwise equal-length. This is asserted, not
// thrown: the function sits inside an O(n^2) k-NN sweep where an always-on
// branch is measurable, and every in-tree caller obtains both spans from the
// same Block_CSR::row_range(), which makes a mismatch a programming error
// rather than a data error. Data-level invariants are validated once, and do
// throw, in Block_CSR's constructor.
template <typename T, typename Index>
[[nodiscard]] inline T sparse_dot(std::span<const Index> idx_a, std::span<const T> val_a,
                                   std::span<const Index> idx_b, std::span<const T> val_b) {
  assert(idx_a.size() == val_a.size() && "sparse_dot: idx_a/val_a length mismatch");
  assert(idx_b.size() == val_b.size() && "sparse_dot: idx_b/val_b length mismatch");

  return checked_narrow<T>(sparse_dot_wide<T, Index>(idx_a, val_a, idx_b, val_b));
}

// Squared L2 norm of a sparse row's value array. Contiguous and branch-free, so
// compilers auto-vectorize this well at -O2/-O3; it is deliberately left as
// plain scalar code rather than hand-written intrinsics.
//
// Accumulates in the wide type. Squaring is where finite input turns infinite:
// a value of 1e30 is a perfectly ordinary float, but 1e30^2 = 1e60 is not, and
// the float sum saturates to inf. That inf then divides an inf dot product and
// yields NaN -- a NaN manufactured from input that passed every finiteness
// check at the boundary. Found by fuzzing at 300k iterations; 20k did not reach it.
template <typename T>
[[nodiscard]] inline WideAcc<T> sparse_norm_sq(std::span<const T> val) {
  WideAcc<T> acc{0};
  for (const T v : val) acc += static_cast<WideAcc<T>>(v) * static_cast<WideAcc<T>>(v);
  return acc;
}

// L2 norm, with the square root taken before narrowing back to T. sqrt(1e60)
// is 1e30, which is representable -- so the norm survives values whose square
// does not.
template <typename T>
[[nodiscard]] inline T sparse_norm(std::span<const T> val) {
  return static_cast<T>(std::sqrt(sparse_norm_sq<T>(val)));
}

// Combine a precomputed dot product with precomputed row norms.
//
// The denominator is formed in the wide type: two finite float norms of 1e30
// multiply to 1e60, which overflows float but not double.
//
// The final finiteness guard is the invariant that matters. A dot product can
// still saturate to inf for extreme magnitudes no matter how the norms are
// computed, and inf/inf is NaN -- which, reaching the top-k comparator, is
// undefined behaviour in std::sort_heap rather than merely a bad score. Yielding
// 0 for such a pair is a deliberate, documented degradation: a similarity that
// cannot be represented is reported as no similarity, never as NaN.
template <typename T>
[[nodiscard]] inline T cosine_from_dot_wide(WideAcc<T> dot, T norm_a, T norm_b) noexcept {
  if (dot == WideAcc<T>{0}) return T{0};
  const WideAcc<T> denom =
      static_cast<WideAcc<T>>(norm_a) * static_cast<WideAcc<T>>(norm_b);
  if (denom == WideAcc<T>{0}) return T{0};
  const WideAcc<T> r = dot / denom;
  if (!std::isfinite(r)) return T{0};

  // Finiteness must be re-checked AFTER the narrowing cast, not only on the
  // wide quotient.
  //
  // WideAcc<float> is double, so `r` can be finite in double and still
  // overflow float: cosine_from_dot(1.0f, FLT_DENORM_MIN, 1.0f) computes
  // r = 7.1e44, passes isfinite(r), and returns +inf. Other reachable
  // combinations do the same -- (1, 1e-30, 1e-20) and (1e30, 1e-20, 1e-20)
  // among them.
  //
  // That defeats the entire purpose of the guard directly above, whose
  // documented contract is that a similarity which cannot be represented is
  // reported as no similarity. An infinite score is worse than a NaN for the
  // caller: it is a valid, totally-ordered value that compares greater than
  // every real similarity, so it does not trip the sort_heap precondition --
  // it silently wins the top-k and makes a degenerate row the nearest
  // neighbour of everything it touches.
  //
  // Not reachable through a Block_CSR built by this library (from_raw()
  // rejects non-finite values and caps ||row||^2 at 1e34, and for consistent
  // dot/norm triples Cauchy-Schwarz bounds the ratio at 1). This is
  // defence-in-depth for direct callers of a public header, and it costs one
  // predictable compare on a path that already branches twice.
  return checked_narrow<T>(r);
}

// Convenience overload for a dot product that is already narrow.
//
// Kept because callers outside the k-NN sweep have a T in hand, but it is the
// lossy entry point BY CONSTRUCTION: whatever headroom the wide accumulator
// bought was spent before this was called. `build_knn_graph_bruteforce` uses
// the wide form above, and so should anything else computing dot and norms
// together -- a dot of 1e60 over a denominator of 1e60 is 1.0 in the wide form
// and is a float `+inf` narrowed to 0, i.e. "matches nothing", in this one.
template <typename T>
[[nodiscard]] inline T cosine_from_dot(T dot, T norm_a, T norm_b) noexcept {
  return cosine_from_dot_wide<T>(static_cast<WideAcc<T>>(dot), norm_a, norm_b);
}

// Cosine similarity between two sparse rows: dot(a, b) / (||a|| * ||b||).
// Norms are computed directly from the sparse values (no zero contributions),
// so this is exact, not an approximation.
//
// Note for callers doing many-vs-many work: this recomputes both norms on every
// call, which is O(nnz) of redundant work per pair. build_knn_graph_bruteforce
// precomputes row norms once and uses cosine_from_dot() below instead.
template <typename T, typename Index>
[[nodiscard]] inline T sparse_cosine_similarity(std::span<const Index> idx_a,
                                                 std::span<const T> val_a,
                                                 std::span<const Index> idx_b,
                                                 std::span<const T> val_b) {
  assert(idx_a.size() == val_a.size() && "sparse_cosine_similarity: a length mismatch");
  assert(idx_b.size() == val_b.size() && "sparse_cosine_similarity: b length mismatch");

  // Wide from the dot straight into the division. Narrowing in between is what
  // made two identical 1e30-scaled rows score 0 instead of 1.
  const WideAcc<T> dot = sparse_dot_wide<T, Index>(idx_a, val_a, idx_b, val_b);
  if (dot == WideAcc<T>{0}) return T{0};

  return cosine_from_dot_wide<T>(dot, sparse_norm<T>(val_a), sparse_norm<T>(val_b));
}

// ---------------------------------------------------------------------------
// Dense scratch buffer for the scatter/gather kernels
//
// The gather kernels compute dot(dense_row, sparse_row) where the dense side is
// one expanded sparse row. In a many-vs-many sweep this is a large win: the
// "query" row is scattered into the scratch once, then every candidate row is a
// single O(nnz_candidate) gather pass with no merge-join branching.
//
// Cost is O(n_cols) floats of residency per thread (n_cols is the gene count for
// scRNA-seq, so typically 2k-30k -> 8-120 KB, comfortably L2-resident). Clearing
// touches only the scattered positions, never the whole buffer, keeping the
// per-row cost O(nnz) rather than O(n_cols).
// ---------------------------------------------------------------------------
class DenseScratch {
 public:
  // Sizes the buffer to hold n_cols floats, ALL ZERO.
  //
  // The zeroing is unconditional. Previously this only allocated when growing,
  // so re-sizing an already-large buffer left whatever the last user wrote in
  // place -- and the gather kernels read every position a row's column indices
  // point at, so stale values silently contaminate similarity scores rather than
  // failing. Reusing one scratch across two matrices was enough to trigger it.
  void resize(std::size_t n_cols) {
    if (buf_.size() < n_cols) {
      buf_.assign(n_cols, 0.0f);
    } else {
      std::fill(buf_.begin(), buf_.begin() + static_cast<std::ptrdiff_t>(n_cols), 0.0f);
    }
    n_cols_ = n_cols;
  }

  // Expands a sparse row into the dense buffer.
  void scatter(std::span<const std::int32_t> idx, std::span<const float> val) {
    assert(idx.size() == val.size() && "DenseScratch::scatter: length mismatch");
    for (std::size_t i = 0; i < idx.size(); ++i) {
      assert(idx[i] >= 0 && static_cast<std::size_t>(idx[i]) < n_cols_ &&
             "DenseScratch::scatter: column index out of range");
      buf_[static_cast<std::size_t>(idx[i])] = val[i];
    }
  }

  // Restores the buffer to all-zero by clearing only the positions scatter()
  // touched. Must be called with the same index span that was scattered.
  void clear(std::span<const std::int32_t> idx) {
    for (const std::int32_t c : idx) {
      assert(c >= 0 && static_cast<std::size_t>(c) < n_cols_ &&
             "DenseScratch::clear: column index out of range");
      buf_[static_cast<std::size_t>(c)] = 0.0f;
    }
  }

  [[nodiscard]] const float* data() const noexcept { return buf_.data(); }
  [[nodiscard]] std::size_t size() const noexcept { return n_cols_; }

 private:
  std::vector<float> buf_;
  std::size_t n_cols_ = 0;
};

// ---------------------------------------------------------------------------
// Gather kernels: dot(dense, sparse)
//
// Every variant computes sum over nz of dense[idx[nz]] * val[nz]. `dense` must
// have at least max(idx)+1 elements; Block_CSR validates that column indices are
// in [0, n_cols) at construction time, which is what makes the unchecked
// hardware gather below memory-safe. The same guarantee is what makes the
// NEON path's lane-wise `dense[idx[i]]` loads safe -- NEON has no gather
// instruction, but an out-of-range index is just as out of bounds there.
// ---------------------------------------------------------------------------

// Reference gather, accumulated wide. This is also the EXACT fallback the
// dispatcher below drops to when a vector kernel saturates, which is why it is
// worth having a widened version even though it is not the fast path -- and,
// under SCRNA_FORCE_ISA=scalar, it IS the fast path.
//
// Four independent `double` accumulators rather than one, for the same
// latency-vs-throughput reason as the vector kernels' own unroll
// (dense_dot_avx2's comment): a single running `acc +=` chain cannot start
// its next multiply-add until the previous one has fully landed, so it runs
// at FMA latency, not at the throughput several independent chains let the
// core reach. Combined pairwise rather than in declaration order --
// (acc0+acc1)+(acc2+acc3) is a shorter dependency chain for the final
// reduction than a strictly left-to-right fold. All four stay `double`; this
// is the widened path, not the narrow one `checked_narrow` guards below.
[[nodiscard]] inline double gather_dot_scalar_wide(const float* dense,
                                                   const std::int32_t* idx,
                                                   const float* val, std::size_t nnz) {
  double acc0 = 0.0, acc1 = 0.0, acc2 = 0.0, acc3 = 0.0;
  std::size_t i = 0;
  for (; i + 4 <= nnz; i += 4) {
    acc0 += static_cast<double>(dense[idx[i]]) * static_cast<double>(val[i]);
    acc1 += static_cast<double>(dense[idx[i + 1]]) * static_cast<double>(val[i + 1]);
    acc2 += static_cast<double>(dense[idx[i + 2]]) * static_cast<double>(val[i + 2]);
    acc3 += static_cast<double>(dense[idx[i + 3]]) * static_cast<double>(val[i + 3]);
  }
  double acc = (acc0 + acc1) + (acc2 + acc3);
  for (; i < nnz; ++i) {
    acc += static_cast<double>(dense[idx[i]]) * static_cast<double>(val[i]);
  }
  return acc;
}

[[nodiscard]] inline float gather_dot_scalar(const float* dense, const std::int32_t* idx,
                                              const float* val, std::size_t nnz) {
  return checked_narrow<float>(gather_dot_scalar_wide(dense, idx, val, nnz));
}

#if defined(SCRNA_BUILD_X86_SIMD) && defined(SCRNA_HAVE_AVX2)

// AVX2 horizontal sum of 8 packed floats.
[[nodiscard]] __attribute__((target("avx2,fma")))
inline float hsum256_ps(__m256 v) {
  __m128 lo = _mm256_castps256_ps128(v);
  __m128 hi = _mm256_extractf128_ps(v, 1);
  lo = _mm_add_ps(lo, hi);
  __m128 shuf = _mm_movehdup_ps(lo);
  __m128 sums = _mm_add_ps(lo, shuf);
  shuf = _mm_movehl_ps(shuf, sums);
  sums = _mm_add_ss(sums, shuf);
  return _mm_cvtss_f32(sums);
}

// AVX2 gather-based dot(dense, sparse). Processes 8 lanes per iteration.
// AVX2 has no lane-masked load, so the < 8 element remainder is a scalar tail.
[[nodiscard]] __attribute__((target("avx2,fma")))
inline float gather_dot_avx2(const float* dense, const std::int32_t* idx,
                              const float* val, std::size_t nnz) {
  constexpr std::size_t kLanes = 8;  // AVX2: 8 floats per 256-bit register
  __m256 acc = _mm256_setzero_ps();
  std::size_t i = 0;
  for (; i + kLanes <= nnz; i += kLanes) {
    const __m256i vidx = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(idx + i));
    // NB: AVX2 gather argument order is (base, index, scale) -- the AVX-512
    // form below is (index, base, scale). Easy to transpose by accident.
    const __m256 vgat = _mm256_i32gather_ps(dense, vidx, 4);
    const __m256 vval = _mm256_loadu_ps(val + i);
    acc = _mm256_fmadd_ps(vgat, vval, acc);
  }
  float sum = hsum256_ps(acc);
  // The tail accumulates into the register the horizontal reduction just
  // produced, so it is narrow by necessity; the dispatcher recomputes the
  // whole kernel in double if this comes back non-finite.
  // SCRNA_FLOAT_ACC_OK: vector-kernel scalar tail.
  for (; i < nnz; ++i) sum += dense[idx[i]] * val[i];
  return sum;
}

// Dense dot product over 8-float vectors, 4-way unrolled: four independent
// accumulators instead of one. A single accumulator makes every `fmadd` wait
// on the PREVIOUS iteration's result -- the FMA's own latency (multiple
// cycles on every relevant microarchitecture), not its throughput, is then
// the loop's speed limit, however many FMA ports the core has. Four
// independent chains give the core four iterations' worth of independent
// work in flight, which is what lets it issue at throughput instead of
// stalling on latency; they are combined into one register only once, after
// the unrolled loop, so the loop itself never re-serialises them.
//
// AVX2 has no lane-masked load, so the remainder (< 4*8 = 32 elements) falls
// through a single-vector loop and then a scalar tail -- note the 8 here
// against the AVX-512 version's 16.
[[nodiscard]] __attribute__((target("avx2,fma")))
inline float dense_dot_avx2(const float* a, const float* b, std::size_t n) {
  constexpr std::size_t kLanes = 8;
  __m256 acc0 = _mm256_setzero_ps();
  __m256 acc1 = _mm256_setzero_ps();
  __m256 acc2 = _mm256_setzero_ps();
  __m256 acc3 = _mm256_setzero_ps();
  std::size_t i = 0;
  for (; i + 4 * kLanes <= n; i += 4 * kLanes) {
    acc0 = _mm256_fmadd_ps(_mm256_loadu_ps(a + i), _mm256_loadu_ps(b + i), acc0);
    acc1 = _mm256_fmadd_ps(_mm256_loadu_ps(a + i + kLanes), _mm256_loadu_ps(b + i + kLanes), acc1);
    acc2 = _mm256_fmadd_ps(_mm256_loadu_ps(a + i + 2 * kLanes), _mm256_loadu_ps(b + i + 2 * kLanes),
                            acc2);
    acc3 = _mm256_fmadd_ps(_mm256_loadu_ps(a + i + 3 * kLanes), _mm256_loadu_ps(b + i + 3 * kLanes),
                            acc3);
  }
  for (; i + kLanes <= n; i += kLanes) {
    acc0 = _mm256_fmadd_ps(_mm256_loadu_ps(a + i), _mm256_loadu_ps(b + i), acc0);
  }
  float sum =
      hsum256_ps(_mm256_add_ps(_mm256_add_ps(acc0, acc1), _mm256_add_ps(acc2, acc3)));
  // The tail accumulates into the register the horizontal reduction just
  // produced, so it is narrow by necessity; the dispatcher recomputes the
  // whole kernel in double if this comes back non-finite.
  // SCRNA_FLOAT_ACC_OK: vector-kernel scalar tail.
  for (; i < n; ++i) sum += a[i] * b[i];
  return sum;
}

#endif  // SCRNA_HAVE_AVX2

#if defined(SCRNA_BUILD_ARM_NEON)

// ---------------------------------------------------------------------------
// ARM64 NEON kernels
//
// THE CENTRAL CONSTRAINT: NEON has no gather instruction.
//
// AVX2 has `_mm256_i32gather_ps`; AVX-512 has `_mm512_i32gather_ps`; SVE has
// `svld1_gather_*`. Baseline ARMv8-A Advanced SIMD has nothing equivalent, and
// Apple Silicon does not expose SVE. So the "gather" below is four ordinary
// scalar loads assembled into a vector lane by lane.
//
// That bounds the achievable speed-up, and it is worth being explicit about why
// rather than discovering it in a profile. Per 4 elements:
//
//   scalar : 4 loads + 4 multiplies + 4 adds, each add serially dependent on
//            the last (one long dependency chain through `acc`)
//   NEON   : 4 loads (unchanged -- no gather to save them) + 1 vector FMA
//
// The loads do not get cheaper. What improves is the arithmetic and, above all,
// the dependency chain: the scalar accumulator serialises every iteration at
// FP-add latency, whereas two independent vector accumulators keep several FMAs
// in flight.
//
// MEASURED on an M4 Pro (best of 7, 2000 nnz, GCC 16 -O3 -mcpu=native):
//
//     gather_dot  scalar  1414 ns/call   1.41 G elem/s
//     gather_dot  NEON     314 ns/call   6.36 G elem/s    4.50x
//
// 4.5x is more than the load-bound reasoning above would predict, which says
// the scalar version was dependency-chain bound rather than load bound -- the
// loads were already overlapping, and what was costing time was the serialised
// accumulate. Do not read 4.5x as "the gather got vectorised"; it did not.
//
// End-to-end k-NN gains less (~2x) because scatter, clear, the top-k heap and
// the norm pass do not vectorise. See docs/PERFORMANCE_2026-08-15.md.
// ---------------------------------------------------------------------------

// Assembles four gathered floats into one vector.
//
// `vsetq_lane_f32` with a literal lane index compiles to a single `ld1 {v.s}[n]`
// on AArch64 -- the load targets the lane directly, with no round trip through
// a GPR and no stack spill. Writing this as a `float tmp[4]` staging array
// instead costs a store-forwarding stall per group of four, which is precisely
// the wrong trade in a kernel that is already load-bound.
//
// Bounds are the caller's contract: Block_CSR::validate() proves every column
// index lies in [0, n_cols) at construction, and DenseScratch is sized to
// n_cols. That is what makes an unchecked load safe here, exactly as it is for
// the x86 hardware gathers above.
[[nodiscard]] inline float32x4_t gather4_neon(const float* dense,
                                              const std::int32_t* idx) noexcept {
  float32x4_t v = vdupq_n_f32(0.0f);
  v = vsetq_lane_f32(dense[idx[0]], v, 0);
  v = vsetq_lane_f32(dense[idx[1]], v, 1);
  v = vsetq_lane_f32(dense[idx[2]], v, 2);
  v = vsetq_lane_f32(dense[idx[3]], v, 3);
  return v;
}

// NEON gather-based dot(dense, sparse).
//
// Two independent accumulators over 8 elements per iteration. Two rather than
// one because a single accumulator serialises on FMA latency (~4 cycles on
// M-series) and would leave most of the pipeline idle; two rather than four
// because at four the kernel is load-port bound and the extra registers buy
// nothing. Eight elements per iteration also matches the AVX2 kernel's
// granularity, so the tail-boundary tests exercise the same residue classes on
// both architectures.
[[nodiscard]] inline float gather_dot_neon(const float* dense, const std::int32_t* idx,
                                           const float* val, std::size_t nnz) noexcept {
  constexpr std::size_t kLanes = 4;
  float32x4_t acc0 = vdupq_n_f32(0.0f);
  float32x4_t acc1 = vdupq_n_f32(0.0f);

  std::size_t i = 0;
  for (; i + 2 * kLanes <= nnz; i += 2 * kLanes) {
    const float32x4_t g0 = gather4_neon(dense, idx + i);
    const float32x4_t g1 = gather4_neon(dense, idx + i + kLanes);
    acc0 = vfmaq_f32(acc0, g0, vld1q_f32(val + i));
    acc1 = vfmaq_f32(acc1, g1, vld1q_f32(val + i + kLanes));
  }
  // One more 4-wide step when 4 <= remainder < 8.
  for (; i + kLanes <= nnz; i += kLanes) {
    acc0 = vfmaq_f32(acc0, gather4_neon(dense, idx + i), vld1q_f32(val + i));
  }

  // vaddvq_f32 is a single AArch64 instruction (`faddp`-tree), unlike the
  // shuffle sequence the x86 paths need.
  float sum = vaddvq_f32(vaddq_f32(acc0, acc1));

  // NEON has no lane-masked load, so the < 4 remainder is a scalar tail --
  // same structure as the AVX2 kernel, different width.
  // The tail accumulates into the register the horizontal reduction just
  // produced, so it is narrow by necessity; the dispatcher recomputes the
  // whole kernel in double if this comes back non-finite.
  // SCRNA_FLOAT_ACC_OK: vector-kernel scalar tail.
  for (; i < nnz; ++i) sum += dense[idx[i]] * val[i];
  return sum;
}

// Dense dot product. No gather, so this is the kernel that actually saturates:
// two contiguous vector loads and an FMA per 4 elements, nothing scalar in the
// steady state. Used by the hnswlib space adapter, which is the HNSW hot path.
//
// Four independent accumulators rather than two (see dense_dot_avx2's comment
// for why this is a latency-vs-throughput fix, not a cosmetic one): a single
// `vfmaq_f32` chain is bound by FMA latency regardless of how many FMA
// pipelines the core has, and doubling to two chains only halves that
// bottleneck. Combined into one vector only after the unrolled loop ends.
[[nodiscard]] inline float dense_dot_neon(const float* a, const float* b,
                                          std::size_t n) noexcept {
  constexpr std::size_t kLanes = 4;
  float32x4_t acc0 = vdupq_n_f32(0.0f);
  float32x4_t acc1 = vdupq_n_f32(0.0f);
  float32x4_t acc2 = vdupq_n_f32(0.0f);
  float32x4_t acc3 = vdupq_n_f32(0.0f);

  std::size_t i = 0;
  for (; i + 4 * kLanes <= n; i += 4 * kLanes) {
    acc0 = vfmaq_f32(acc0, vld1q_f32(a + i), vld1q_f32(b + i));
    acc1 = vfmaq_f32(acc1, vld1q_f32(a + i + kLanes), vld1q_f32(b + i + kLanes));
    acc2 = vfmaq_f32(acc2, vld1q_f32(a + i + 2 * kLanes), vld1q_f32(b + i + 2 * kLanes));
    acc3 = vfmaq_f32(acc3, vld1q_f32(a + i + 3 * kLanes), vld1q_f32(b + i + 3 * kLanes));
  }
  for (; i + kLanes <= n; i += kLanes) {
    acc0 = vfmaq_f32(acc0, vld1q_f32(a + i), vld1q_f32(b + i));
  }

  float sum = vaddvq_f32(vaddq_f32(vaddq_f32(acc0, acc1), vaddq_f32(acc2, acc3)));
  // The tail accumulates into the register the horizontal reduction just
  // produced, so it is narrow by necessity; the dispatcher recomputes the
  // whole kernel in double if this comes back non-finite.
  // SCRNA_FLOAT_ACC_OK: vector-kernel scalar tail.
  for (; i < n; ++i) sum += a[i] * b[i];
  return sum;
}

#endif  // SCRNA_BUILD_ARM_NEON

#if defined(SCRNA_BUILD_X86_SIMD) && defined(SCRNA_HAVE_AVX512)

// AVX-512 gather-based dot(dense, sparse). Processes 16 lanes per iteration --
// twice the AVX2 width. The tail is handled with a lane mask rather than a
// scalar loop, which sidesteps the classic 8-vs-16 remainder-bound bug that
// copy-pasting the AVX2 version would introduce.
[[nodiscard]] __attribute__((target("avx512f,avx512bw")))
inline float gather_dot_avx512(const float* dense, const std::int32_t* idx,
                                const float* val, std::size_t nnz) {
  constexpr std::size_t kLanes = 16;  // AVX-512: 16 floats per 512-bit register
  __m512 acc = _mm512_setzero_ps();
  std::size_t i = 0;

  for (; i + kLanes <= nnz; i += kLanes) {
    const __m512i vidx = _mm512_loadu_si512(reinterpret_cast<const void*>(idx + i));
    const __m512 vgat = _mm512_i32gather_ps(vidx, dense, 4);
    const __m512 vval = _mm512_loadu_ps(val + i);
    acc = _mm512_fmadd_ps(vgat, vval, acc);
  }

  if (i < nnz) {
    const unsigned remaining = static_cast<unsigned>(nnz - i);
    const __mmask16 m = static_cast<__mmask16>((1u << remaining) - 1u);
    const __m512i vidx = _mm512_maskz_loadu_epi32(m, idx + i);
    const __m512 vval = _mm512_maskz_loadu_ps(m, val + i);
    // Masked-off lanes read nothing (no fault) and contribute zero.
    const __m512 vgat = _mm512_mask_i32gather_ps(_mm512_setzero_ps(), m, vidx, dense, 4);
    acc = _mm512_fmadd_ps(vgat, vval, acc);
  }

  return _mm512_reduce_add_ps(acc);
}

// Dense dot product over 16-float vectors, with a masked tail.
[[nodiscard]] __attribute__((target("avx512f,avx512bw")))
inline float dense_dot_avx512(const float* a, const float* b, std::size_t n) {
  constexpr std::size_t kLanes = 16;
  __m512 acc = _mm512_setzero_ps();
  std::size_t i = 0;
  for (; i + kLanes <= n; i += kLanes) {
    acc = _mm512_fmadd_ps(_mm512_loadu_ps(a + i), _mm512_loadu_ps(b + i), acc);
  }
  if (i < n) {
    const __mmask16 m = static_cast<__mmask16>((1u << static_cast<unsigned>(n - i)) - 1u);
    acc = _mm512_fmadd_ps(_mm512_maskz_loadu_ps(m, a + i), _mm512_maskz_loadu_ps(m, b + i), acc);
  }
  return _mm512_reduce_add_ps(acc);
}

#endif  // SCRNA_HAVE_AVX512

// Four independent accumulators, same reasoning as gather_dot_scalar_wide
// above -- and the same note that this is dense_dot's actual fast path under
// SCRNA_FORCE_ISA=scalar, not merely the saturation fallback.
[[nodiscard]] inline double dense_dot_scalar_wide(const float* a, const float* b,
                                                  std::size_t n) {
  double acc0 = 0.0, acc1 = 0.0, acc2 = 0.0, acc3 = 0.0;
  std::size_t i = 0;
  for (; i + 4 <= n; i += 4) {
    acc0 += static_cast<double>(a[i]) * static_cast<double>(b[i]);
    acc1 += static_cast<double>(a[i + 1]) * static_cast<double>(b[i + 1]);
    acc2 += static_cast<double>(a[i + 2]) * static_cast<double>(b[i + 2]);
    acc3 += static_cast<double>(a[i + 3]) * static_cast<double>(b[i + 3]);
  }
  double acc = (acc0 + acc1) + (acc2 + acc3);
  for (; i < n; ++i) {
    acc += static_cast<double>(a[i]) * static_cast<double>(b[i]);
  }
  return acc;
}

[[nodiscard]] inline float dense_dot_scalar(const float* a, const float* b, std::size_t n) {
  return checked_narrow<float>(dense_dot_scalar_wide(a, b, n));
}

// ---------------------------------------------------------------------------
// Runtime ISA dispatch
//
// Resolved once, lazily, on first use. Compile-time gating alone is not enough:
// the build machine's capabilities are not the run machine's.
// ---------------------------------------------------------------------------

enum class IsaLevel { kScalar, kNeon, kAvx2, kAvx512 };

// Test/CI override. SCRNA_FORCE_ISA=scalar|avx2|avx512 pins the kernel choice,
// so a single machine can exercise every compiled path rather than only the one
// its CPUID happens to select. Requesting a path that was not compiled in is
// ignored. Requesting one the hardware cannot execute will fault -- that is the
// caller's explicit choice, which is why this is opt-in and never a default.
// SCRNA_FORCE_SCALAR_DISPATCH is kept as a shorthand for the common CI case.
[[nodiscard]] inline bool env_is(const char* name, const char* value) noexcept {
  const char* v = std::getenv(name);
  if (v == nullptr) return false;
  const char* p = value;
  while (*v != '\0' && *p != '\0' && *v == *p) { ++v; ++p; }
  return *v == '\0' && *p == '\0';
}

[[nodiscard]] inline IsaLevel detect_isa() noexcept {
  if (const char* forced = std::getenv("SCRNA_FORCE_SCALAR_DISPATCH");
      forced != nullptr && forced[0] != '\0' && forced[0] != '0') {
    return IsaLevel::kScalar;
  }
  if (env_is("SCRNA_FORCE_ISA", "scalar")) return IsaLevel::kScalar;
#if defined(SCRNA_BUILD_ARM_NEON)
  if (env_is("SCRNA_FORCE_ISA", "neon")) return IsaLevel::kNeon;
#endif
#if defined(SCRNA_BUILD_X86_SIMD) && defined(SCRNA_HAVE_AVX512)
  if (env_is("SCRNA_FORCE_ISA", "avx512")) return IsaLevel::kAvx512;
#endif
#if defined(SCRNA_BUILD_X86_SIMD) && defined(SCRNA_HAVE_AVX2)
  if (env_is("SCRNA_FORCE_ISA", "avx2")) return IsaLevel::kAvx2;
#endif

#if defined(SCRNA_BUILD_X86_SIMD)
#if defined(SCRNA_HAVE_AVX512)
  if (__builtin_cpu_supports("avx512f") && __builtin_cpu_supports("avx512bw")) {
    return IsaLevel::kAvx512;
  }
#endif
#if defined(SCRNA_HAVE_AVX2)
  if (__builtin_cpu_supports("avx2") && __builtin_cpu_supports("fma")) {
    return IsaLevel::kAvx2;
  }
#endif
#endif

  // NEON is architecturally guaranteed on AArch64 -- if it compiled in, the
  // CPU has it. No probe, and nothing to fall back from.
#if defined(SCRNA_BUILD_ARM_NEON)
  return IsaLevel::kNeon;
#else
  return IsaLevel::kScalar;
#endif
}

// Thread-safe (C++11 magic statics) and evaluated after static init, so no
// __builtin_cpu_init() ordering concerns.
[[nodiscard]] inline IsaLevel isa_level() noexcept {
  static const IsaLevel level = detect_isa();
  return level;
}

[[nodiscard]] inline const char* isa_name() noexcept {
  switch (isa_level()) {
    case IsaLevel::kAvx512: return "avx512";
    case IsaLevel::kAvx2:   return "avx2";
    case IsaLevel::kNeon:   return "neon";
    case IsaLevel::kScalar: break;
  }
  return "scalar";
}

// ---------------------------------------------------------------------------
// Dispatch, and where the vector kernels' float accumulators are made safe
//
// The ISA kernels above accumulate in float, and that is deliberate rather than
// an unfixed instance of the defect this header opens by describing. Widening
// them means halving the lane count -- the NEON gather kernel is 4.5x the
// scalar one and is the module's headline benchmark -- to buy correctness on
// inputs that `Block_CSR::validate()` already rejects (it caps ||row||^2 at
// 1e34 precisely so this accumulator provably cannot saturate).
//
// The exposure is direct callers of this public header, who have no Block_CSR
// and no such guarantee. They are covered here instead, for the cost of one
// perfectly-predicted branch per row rather than 2x on every FMA:
//
//   run the vector kernel; if the answer is not finite, it saturated, so
//   recompute it exactly with the widened scalar kernel.
//
// That guard is complete, not a heuristic. Saturation cannot hide: once any
// partial sum or any single product reaches +/-inf, every subsequent operation
// propagates it (inf + x = inf, inf - inf = NaN, inf * 0 = NaN), a horizontal
// reduction over a vector containing inf yields inf or NaN, and both fail
// isfinite. There is no input for which a float accumulator overflows and the
// final result comes back finite and wrong.
//
// What it does not catch is ordinary rounding drift, which is a precision
// question rather than a saturation one and is what the parity tests bound.
// ---------------------------------------------------------------------------

[[nodiscard]] inline double gather_dot_wide(const float* dense, const std::int32_t* idx,
                                            const float* val, std::size_t nnz) {
  float fast;
  switch (isa_level()) {
#if defined(SCRNA_BUILD_X86_SIMD) && defined(SCRNA_HAVE_AVX512)
    case IsaLevel::kAvx512: fast = gather_dot_avx512(dense, idx, val, nnz); break;
#endif
#if defined(SCRNA_BUILD_X86_SIMD) && defined(SCRNA_HAVE_AVX2)
    case IsaLevel::kAvx2:   fast = gather_dot_avx2(dense, idx, val, nnz); break;
#endif
#if defined(SCRNA_BUILD_ARM_NEON)
    case IsaLevel::kNeon:   fast = gather_dot_neon(dense, idx, val, nnz); break;
#endif
    default:
      // The scalar path is already wide; there is nothing to fall back from.
      return gather_dot_scalar_wide(dense, idx, val, nnz);
  }
  if (!std::isfinite(fast)) return gather_dot_scalar_wide(dense, idx, val, nnz);
  return static_cast<double>(fast);
}

// Dispatching gather dot product. This is the hot inner kernel of the k-NN
// sweep; the branch resolves to a cached static and predicts perfectly.
//
// Prefer `gather_dot_wide` when the result is about to be divided by a norm:
// this one narrows, and a dot product of 1e60 narrows to 0 however carefully it
// was computed.
[[nodiscard]] inline float gather_dot(const float* dense, const std::int32_t* idx,
                                       const float* val, std::size_t nnz) {
  return checked_narrow<float>(gather_dot_wide(dense, idx, val, nnz));
}

[[nodiscard]] inline double gather_dot_wide(const DenseScratch& scratch,
                                            std::span<const std::int32_t> idx,
                                            std::span<const float> val) {
  assert(idx.size() == val.size() && "gather_dot_wide: idx/val length mismatch");
  return gather_dot_wide(scratch.data(), idx.data(), val.data(), idx.size());
}

[[nodiscard]] inline float gather_dot(const DenseScratch& scratch,
                                       std::span<const std::int32_t> idx,
                                       std::span<const float> val) {
  assert(idx.size() == val.size() && "gather_dot: idx/val length mismatch");
  return gather_dot(scratch.data(), idx.data(), val.data(), idx.size());
}

// Dispatching dense dot product, for consumers that hand us two contiguous
// vectors rather than a sparse/dense pair -- notably the hnswlib space adapter,
// whose distance-function signature is dense-only.
[[nodiscard]] inline double dense_dot_wide(const float* a, const float* b, std::size_t n) {
  float fast;
  switch (isa_level()) {
#if defined(SCRNA_BUILD_X86_SIMD) && defined(SCRNA_HAVE_AVX512)
    case IsaLevel::kAvx512: fast = dense_dot_avx512(a, b, n); break;
#endif
#if defined(SCRNA_BUILD_X86_SIMD) && defined(SCRNA_HAVE_AVX2)
    case IsaLevel::kAvx2:   fast = dense_dot_avx2(a, b, n); break;
#endif
#if defined(SCRNA_BUILD_ARM_NEON)
    case IsaLevel::kNeon:   fast = dense_dot_neon(a, b, n); break;
#endif
    default: return dense_dot_scalar_wide(a, b, n);
  }
  if (!std::isfinite(fast)) return dense_dot_scalar_wide(a, b, n);
  return static_cast<double>(fast);
}

[[nodiscard]] inline float dense_dot(const float* a, const float* b, std::size_t n) {
  return checked_narrow<float>(dense_dot_wide(a, b, n));
}

}  // namespace scrna::simd
