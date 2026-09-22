// Minimal, dependency-free sanity tests (no gtest/Catch2 required so the
// project stays "pure C++20, static dependencies"). Exits non-zero on
// failure so it composes with `ctest`.

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <random>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

#include "matrix/block_csr.hpp"
#include "matrix/hnsw_index.hpp"
#include "matrix/knn_graph.hpp"
#include "matrix/simd_math.hpp"

using scrna::AlignedVector;
using scrna::Block_CSR;
using scrna::kCacheLineBytes;

namespace {

int g_failures = 0;

#define SCRNA_CHECK(cond)                                                   \
  do {                                                                      \
    if (!(cond)) {                                                          \
      std::fprintf(stderr, "FAILED: %s (%s:%d)\n", #cond, __FILE__, __LINE__); \
      ++g_failures;                                                         \
    }                                                                       \
  } while (0)

// Like SCRNA_CHECK, but with a caller-supplied std::string detail (e.g. the
// specific values that made this instance of a property-style test fail) --
// unlike a bare condition, a loop over seeds/sizes gives no other way to tell
// which iteration failed.
#define SCRNA_CHECK_MSG(cond, msg)                                          \
  do {                                                                      \
    if (!(cond)) {                                                          \
      std::fprintf(stderr, "FAILED: %s -- %s (%s:%d)\n", #cond,             \
                    std::string(msg).c_str(), __FILE__, __LINE__);          \
      ++g_failures;                                                         \
    }                                                                       \
  } while (0)

// Asserts that `expr` throws std::invalid_argument.
#define SCRNA_CHECK_THROWS(expr)                                            \
  do {                                                                      \
    bool threw = false;                                                     \
    try {                                                                   \
      (void)(expr);                                                         \
    } catch (const std::invalid_argument&) {                                \
      threw = true;                                                         \
    } catch (...) {                                                         \
    }                                                                       \
    if (!threw) {                                                           \
      std::fprintf(stderr, "FAILED (no throw): %s (%s:%d)\n", #expr, __FILE__, __LINE__); \
      ++g_failures;                                                         \
    }                                                                       \
  } while (0)

constexpr float kEps = 1e-5f;

bool near(float a, float b, float eps = kEps) { return std::fabs(a - b) <= eps; }

// Build a tiny 3x4 matrix:
// row0: [1, 0, 2, 0]
// row1: [0, 3, 0, 0]
// row2: [0, 0, 0, 4]
Block_CSR<float> make_toy_matrix() {
  AlignedVector<float> values = {1.0f, 2.0f, 3.0f, 4.0f};
  AlignedVector<std::int32_t> col_idx = {0, 2, 1, 3};
  AlignedVector<std::int32_t> row_ptr = {0, 2, 3, 4};
  return Block_CSR<float>(3, 4, std::move(values), std::move(col_idx), std::move(row_ptr));
}

void test_block_csr_alignment() {
  auto mat = make_toy_matrix();
  SCRNA_CHECK(mat.rows() == 3);
  SCRNA_CHECK(mat.cols() == 4);
  SCRNA_CHECK(mat.nnz() == 4);
  SCRNA_CHECK(scrna::is_cache_aligned(mat.values_data()));
  SCRNA_CHECK(scrna::is_cache_aligned(mat.col_idx_data()));
  // row_ptr was previously unchecked anywhere.
  SCRNA_CHECK(scrna::is_cache_aligned(mat.row_ptr_data()));
  SCRNA_CHECK(!mat.is_zero_copy());

  auto row0_vals = mat.row_values(0);
  auto row0_idx = mat.row_indices(0);
  SCRNA_CHECK(row0_vals.size() == 2);
  SCRNA_CHECK(row0_idx[0] == 0 && row0_idx[1] == 2);
  SCRNA_CHECK(row0_vals[0] == 1.0f && row0_vals[1] == 2.0f);
}

// Structural validation must reject malformed CSR. The column-index bound in
// particular is what keeps the hardware gather kernels in-bounds.
void test_block_csr_validation() {
  // Column index >= n_cols.
  SCRNA_CHECK_THROWS(Block_CSR<float>(1, 2,
                                       AlignedVector<float>{1.0f},
                                       AlignedVector<std::int32_t>{7},
                                       AlignedVector<std::int32_t>{0, 1}));
  // Negative column index.
  SCRNA_CHECK_THROWS(Block_CSR<float>(1, 4,
                                       AlignedVector<float>{1.0f},
                                       AlignedVector<std::int32_t>{-1},
                                       AlignedVector<std::int32_t>{0, 1}));
  // row_ptr not starting at 0.
  SCRNA_CHECK_THROWS(Block_CSR<float>(1, 4,
                                       AlignedVector<float>{1.0f},
                                       AlignedVector<std::int32_t>{0},
                                       AlignedVector<std::int32_t>{1, 1}));
  // row_ptr decreasing.
  SCRNA_CHECK_THROWS(Block_CSR<float>(2, 4,
                                       AlignedVector<float>{1.0f, 2.0f},
                                       AlignedVector<std::int32_t>{0, 1},
                                       AlignedVector<std::int32_t>{0, 2, 1}));
  // row_ptr[n_rows] != nnz.
  SCRNA_CHECK_THROWS(Block_CSR<float>(1, 4,
                                       AlignedVector<float>{1.0f, 2.0f},
                                       AlignedVector<std::int32_t>{0, 1},
                                       AlignedVector<std::int32_t>{0, 1}));
  // values/col_idx length mismatch.
  SCRNA_CHECK_THROWS(Block_CSR<float>(1, 4,
                                       AlignedVector<float>{1.0f, 2.0f},
                                       AlignedVector<std::int32_t>{0},
                                       AlignedVector<std::int32_t>{0, 2}));
  // Wrong row_ptr length.
  SCRNA_CHECK_THROWS(Block_CSR<float>(2, 4,
                                       AlignedVector<float>{1.0f},
                                       AlignedVector<std::int32_t>{0},
                                       AlignedVector<std::int32_t>{0, 1}));
}

// adopt() must borrow in place (no copy) and must refuse misaligned input.
void test_block_csr_adopt() {
  AlignedVector<float> values = {1.0f, 2.0f, 3.0f, 4.0f};
  AlignedVector<std::int32_t> col_idx = {0, 2, 1, 3};
  AlignedVector<std::int32_t> row_ptr = {0, 2, 3, 4};

  const float* vptr = values.data();
  const std::int32_t* cptr = col_idx.data();
  const std::int32_t* rptr = row_ptr.data();

  struct DummyOwner {
    int tag = 0;
  };

  auto mat = Block_CSR<float>::adopt(3, 4, std::span<const float>(values),
                                      std::span<const std::int32_t>(col_idx),
                                      std::span<const std::int32_t>(row_ptr), DummyOwner{});

  // The defining property: the matrix points at the caller's memory.
  SCRNA_CHECK(mat.values_data() == vptr);
  SCRNA_CHECK(mat.col_idx_data() == cptr);
  SCRNA_CHECK(mat.row_ptr_data() == rptr);
  SCRNA_CHECK(mat.is_zero_copy());
  SCRNA_CHECK(mat.nnz() == 4);
  SCRNA_CHECK(near(mat.row_values(0)[1], 2.0f));

  // A copied matrix must NOT alias the source.
  auto copied = Block_CSR<float>::from_raw(3, 4, std::span<const float>(values),
                                            std::span<const std::int32_t>(col_idx),
                                            std::span<const std::int32_t>(row_ptr));
  SCRNA_CHECK(copied.values_data() != vptr);
  SCRNA_CHECK(!copied.is_zero_copy());

  // Misaligned adoption must be rejected, not silently accepted.
  std::vector<float> unaligned_pool(16, 1.0f);
  const float* misaligned = unaligned_pool.data() + 1;  // deliberately off-boundary
  if (!scrna::is_cache_aligned(misaligned)) {
    AlignedVector<std::int32_t> one_idx = {0};
    AlignedVector<std::int32_t> one_ptr = {0, 1};
    SCRNA_CHECK_THROWS(
        Block_CSR<float>::adopt(1, 4, std::span<const float>(misaligned, 1),
                                 std::span<const std::int32_t>(one_idx),
                                 std::span<const std::int32_t>(one_ptr), DummyOwner{}));
  }

  // The owner token must keep foreign memory alive past the caller's scope.
  auto heap_vals = std::make_shared<AlignedVector<float>>(AlignedVector<float>{5.0f, 6.0f});
  auto heap_idx = std::make_shared<AlignedVector<std::int32_t>>(AlignedVector<std::int32_t>{0, 1});
  auto heap_ptr = std::make_shared<AlignedVector<std::int32_t>>(AlignedVector<std::int32_t>{0, 2});
  struct HeapOwner {
    std::shared_ptr<AlignedVector<float>> v;
    std::shared_ptr<AlignedVector<std::int32_t>> c;
    std::shared_ptr<AlignedVector<std::int32_t>> r;
  };
  Block_CSR<float> adopted_heap = Block_CSR<float>::adopt(
      1, 4, std::span<const float>(*heap_vals),
      std::span<const std::int32_t>(*heap_idx), std::span<const std::int32_t>(*heap_ptr),
      HeapOwner{heap_vals, heap_idx, heap_ptr});
  const long use_after = heap_vals.use_count();
  SCRNA_CHECK(use_after > 1);  // matrix holds a reference
  SCRNA_CHECK(near(adopted_heap.row_values(0)[0], 5.0f));
}

// Column indices must be strictly increasing within a row. Unsorted or
// duplicated indices previously produced two *different* wrong answers
// depending on which kernel ran: merge-join scored two identical vectors at
// 0.80, and gather scored a duplicated index at 1.20 -- a cosine above 1.0.
void test_rejects_unsorted_and_duplicate_columns() {
  // Unsorted within row 0.
  SCRNA_CHECK_THROWS(Block_CSR<float>(1, 4,
                                       AlignedVector<float>{2.0f, 1.0f},
                                       AlignedVector<std::int32_t>{2, 0},
                                       AlignedVector<std::int32_t>{0, 2}));
  // Duplicate column within row 0.
  SCRNA_CHECK_THROWS(Block_CSR<float>(1, 4,
                                       AlignedVector<float>{1.0f, 3.0f},
                                       AlignedVector<std::int32_t>{1, 1},
                                       AlignedVector<std::int32_t>{0, 2}));
  // Unsorted in the *second* row only -- the check must be per-row, and must
  // not be fooled by the row boundary resetting to a lower index.
  SCRNA_CHECK_THROWS(Block_CSR<float>(2, 4,
                                       AlignedVector<float>{1.0f, 2.0f, 1.0f, 2.0f},
                                       AlignedVector<std::int32_t>{0, 2, 3, 1},
                                       AlignedVector<std::int32_t>{0, 2, 4}));
  // A descending index across a row boundary is legal: row 1 restarts at 0.
  auto ok = Block_CSR<float>(2, 4,
                              AlignedVector<float>{1.0f, 2.0f, 1.0f, 2.0f},
                              AlignedVector<std::int32_t>{2, 3, 0, 1},
                              AlignedVector<std::int32_t>{0, 2, 4});
  SCRNA_CHECK(ok.nnz() == 4);

  // With the invariant enforced, both kernels must now agree on every matrix
  // they accept. Cross-check them on identical rows.
  auto m = Block_CSR<float>(2, 4,
                             AlignedVector<float>{1.0f, 2.0f, 1.0f, 2.0f},
                             AlignedVector<std::int32_t>{0, 2, 0, 2},
                             AlignedVector<std::int32_t>{0, 2, 4});
  const float merge = scrna::simd::sparse_cosine_similarity<float, std::int32_t>(
      m.row_indices(0), m.row_values(0), m.row_indices(1), m.row_values(1));
  const auto g = scrna::build_knn_graph_bruteforce(m, 1);
  SCRNA_CHECK(near(merge, 1.0f, 1e-4f));
  SCRNA_CHECK(near(g.neighbor_sim[0], 1.0f, 1e-4f));
  SCRNA_CHECK(near(merge, g.neighbor_sim[0], 1e-4f));
  SCRNA_CHECK(g.neighbor_sim[0] <= 1.0f + 1e-4f);  // cosine cannot exceed 1
}

// The span factories must reject mismatched extents outright. This is the case
// the old (pointer, nnz) signature could not express: one length was applied to
// all three buffers, so a short `indices` array was read to the length of
// `values` (heap-buffer-overflow, confirmed under ASan). The mismatch is now a
// property of the arguments themselves.
void test_span_factories_reject_mismatched_extents() {
  AlignedVector<float> values(64, 1.0f);
  AlignedVector<std::int32_t> short_idx = {0, 1};     // must be 64 long
  AlignedVector<std::int32_t> good_idx;               // strictly increasing 0..63
  for (std::int32_t c = 0; c < 64; ++c) good_idx.push_back(c);
  AlignedVector<std::int32_t> good_ptr = {0, 64};     // n_rows = 1
  AlignedVector<std::int32_t> short_ptr = {0};        // must be n_rows + 1

  struct DummyOwner {};

  // indices shorter than values
  SCRNA_CHECK_THROWS(Block_CSR<float>::from_raw(1, 64, std::span<const float>(values),
                                                 std::span<const std::int32_t>(short_idx),
                                                 std::span<const std::int32_t>(good_ptr)));
  SCRNA_CHECK_THROWS(Block_CSR<float>::adopt(1, 64, std::span<const float>(values),
                                              std::span<const std::int32_t>(short_idx),
                                              std::span<const std::int32_t>(good_ptr),
                                              DummyOwner{}));

  // indptr not n_rows + 1
  SCRNA_CHECK_THROWS(Block_CSR<float>::from_raw(1, 64, std::span<const float>(values),
                                                 std::span<const std::int32_t>(good_idx),
                                                 std::span<const std::int32_t>(short_ptr)));
  SCRNA_CHECK_THROWS(Block_CSR<float>::adopt(5, 64, std::span<const float>(values),
                                              std::span<const std::int32_t>(good_idx),
                                              std::span<const std::int32_t>(good_ptr),
                                              DummyOwner{}));

  // Non-empty span with a null pointer.
  SCRNA_CHECK_THROWS(Block_CSR<float>::from_raw(
      1, 64, std::span<const float>(static_cast<const float*>(nullptr), 4),
      std::span<const std::int32_t>(good_idx), std::span<const std::int32_t>(good_ptr)));

  // The well-formed case still succeeds.
  auto ok = Block_CSR<float>::from_raw(1, 64, std::span<const float>(values),
                                        std::span<const std::int32_t>(good_idx),
                                        std::span<const std::int32_t>(good_ptr));
  SCRNA_CHECK(ok.nnz() == 64);
}

// Round-tripping a matrix through its own span accessors must reproduce it
// exactly -- the accessors and the factories agree on extents by construction.
void test_span_accessors_round_trip() {
  auto mat = make_toy_matrix();
  auto copy = Block_CSR<float>::from_raw(mat.rows(), mat.cols(), mat.values_span(),
                                          mat.col_idx_span(), mat.row_ptr_span());
  SCRNA_CHECK(copy.rows() == mat.rows());
  SCRNA_CHECK(copy.nnz() == mat.nnz());
  SCRNA_CHECK(copy.values_span().size() == mat.nnz());
  SCRNA_CHECK(copy.col_idx_span().size() == mat.nnz());
  SCRNA_CHECK(copy.row_ptr_span().size() == mat.rows() + 1);
  for (std::size_t i = 0; i < mat.nnz(); ++i) {
    SCRNA_CHECK(copy.values_span()[i] == mat.values_span()[i]);
    SCRNA_CHECK(copy.col_idx_span()[i] == mat.col_idx_span()[i]);
  }
}

void test_sparse_dot_and_cosine() {
  // a = [1, 0, 2, 0], b = [1, 0, 2, 0] -> identical, cosine == 1
  std::vector<std::int32_t> idx_a = {0, 2};
  std::vector<float> val_a = {1.0f, 2.0f};
  std::vector<std::int32_t> idx_b = {0, 2};
  std::vector<float> val_b = {1.0f, 2.0f};

  const float dot = scrna::simd::sparse_dot<float, std::int32_t>(idx_a, val_a, idx_b, val_b);
  SCRNA_CHECK(dot == 5.0f);  // 1*1 + 2*2

  const float cos = scrna::simd::sparse_cosine_similarity<float, std::int32_t>(idx_a, val_a, idx_b, val_b);
  SCRNA_CHECK(cos > 0.999f && cos < 1.001f);

  // Disjoint support -> orthogonal -> cosine == 0
  std::vector<std::int32_t> idx_c = {1, 3};
  std::vector<float> val_c = {5.0f, 6.0f};
  const float cos_orth = scrna::simd::sparse_cosine_similarity<float, std::int32_t>(idx_a, val_a, idx_c, val_c);
  SCRNA_CHECK(cos_orth == 0.0f);
}

// The dense scratch must round-trip: scatter then clear must leave it all-zero,
// otherwise state leaks between rows of the k-NN sweep and corrupts results.
void test_dense_scratch_roundtrip() {
  scrna::simd::DenseScratch scratch;
  const std::size_t n_cols = 128;
  scratch.resize(n_cols);

  std::vector<std::int32_t> idx = {3, 17, 64, 127};
  std::vector<float> val = {1.5f, -2.0f, 0.25f, 9.0f};

  scratch.scatter(idx, val);
  for (std::size_t i = 0; i < idx.size(); ++i) {
    SCRNA_CHECK(near(scratch.data()[idx[i]], val[i]));
  }
  scratch.clear(idx);
  for (std::size_t c = 0; c < n_cols; ++c) {
    SCRNA_CHECK(scratch.data()[c] == 0.0f);
  }
}

// The dispatched gather kernel must agree with (a) the scalar gather and (b)
// the independent merge-join implementation, across lengths that straddle both
// the 8-lane (AVX2) and 16-lane (AVX-512) boundaries and their remainders.
void test_gather_dot_parity() {
  std::mt19937 rng(12345);
  std::uniform_real_distribution<float> vdist(-2.0f, 2.0f);

  const std::size_t n_cols = 512;
  const std::size_t lengths[] = {0, 1, 2, 7, 8, 9, 15, 16, 17, 31, 32, 33, 63, 64, 65, 200};

  for (const std::size_t nnz : lengths) {
    // Strictly increasing column indices, as CSR requires for the merge-join.
    std::vector<std::int32_t> idx_a;
    std::vector<float> val_a;
    std::vector<std::int32_t> idx_b;
    std::vector<float> val_b;

    std::int32_t c = 0;
    for (std::size_t i = 0; i < nnz && c < static_cast<std::int32_t>(n_cols); ++i) {
      idx_a.push_back(c);
      val_a.push_back(vdist(rng));
      c += 1 + static_cast<std::int32_t>(rng() % 3);
    }
    c = 0;
    for (std::size_t i = 0; i < nnz && c < static_cast<std::int32_t>(n_cols); ++i) {
      idx_b.push_back(c);
      val_b.push_back(vdist(rng));
      c += 1 + static_cast<std::int32_t>(rng() % 2);
    }

    scrna::simd::DenseScratch scratch;
    scratch.resize(n_cols);
    scratch.scatter(idx_a, val_a);

    const float via_dispatch = scrna::simd::gather_dot(scratch, idx_b, val_b);
    const float via_scalar = scrna::simd::gather_dot_scalar(
        scratch.data(), idx_b.data(), val_b.data(), idx_b.size());
    const float via_mergejoin =
        scrna::simd::sparse_dot<float, std::int32_t>(idx_a, val_a, idx_b, val_b);

    SCRNA_CHECK(near(via_dispatch, via_scalar, 1e-4f));
    // Cross-check against a structurally different algorithm.
    SCRNA_CHECK(near(via_dispatch, via_mergejoin, 1e-4f));

    scratch.clear(idx_a);
  }
}

// ---------------------------------------------------------------------------
// SIMD parity, done against an exact reference rather than against another
// approximation.
//
// WHY NOT BITWISE PARITY. It is not achievable and never was, for two
// independent reasons, both of which apply to every vector path here (AVX2,
// AVX-512, NEON):
//
//   1. FUSION. `vfmaq_f32` / `_mm256_fmadd_ps` round once; the scalar
//      `acc += a * b` rounds twice. Different results by construction.
//   2. REASSOCIATION. Vector kernels keep several partial accumulators and sum
//      them at the end. Floating-point addition is not associative, so a
//      different summation order is a different answer.
//
// Demanding bitwise equality would therefore mean forbidding FMA and forbidding
// multiple accumulators -- i.e. forbidding the optimisation. The meaningful
// question is not "are they identical" but "are they both close to the truth,
// and is the vector path no worse".
//
// So these helpers compute the exact-ish value in double and bound each path's
// error against it, using the standard forward error bound for a dot product:
//
//      |computed - exact|  <=  C * n * u * SUM|a_i * b_i|
//
// with u = 2^-24 (float unit roundoff). The bound scales with SUM|a_i*b_i|, not
// with |exact|, which is what makes it correct under cancellation: a dot
// product that sums to ~0 out of large terms is ill-conditioned, and a relative
// tolerance would spuriously fail on it.
// ---------------------------------------------------------------------------

constexpr double kFloatUnitRoundoff = 5.9604644775390625e-08;  // 2^-24

struct DotReference {
  double exact = 0.0;    // sum in double
  double sum_abs = 0.0;  // sum of |a_i * b_i|, the conditioning term
};

// Reference for dot(dense_gathered, val).
DotReference gather_reference(const float* dense, const std::int32_t* idx,
                              const float* val, std::size_t nnz) {
  DotReference r;
  for (std::size_t i = 0; i < nnz; ++i) {
    const double term = static_cast<double>(dense[idx[i]]) * static_cast<double>(val[i]);
    r.exact += term;
    r.sum_abs += std::fabs(term);
  }
  return r;
}

DotReference dense_reference(const float* a, const float* b, std::size_t n) {
  DotReference r;
  for (std::size_t i = 0; i < n; ++i) {
    const double term = static_cast<double>(a[i]) * static_cast<double>(b[i]);
    r.exact += term;
    r.sum_abs += std::fabs(term);
  }
  return r;
}

// C = 8 leaves headroom over the strict n*u*sum_abs bound for the extra
// rounding FMA-free tails and the final horizontal reduction contribute. It is
// still tight enough to catch a real defect: a transposed gather index or a
// dropped tail element moves the result by orders of magnitude more than this.
double dot_error_budget(const DotReference& ref, std::size_t n) {
  constexpr double kC = 8.0;
  // The +1 keeps the budget non-zero for n == 0 and for exactly-representable
  // inputs, where the bound would otherwise be 0 and demand bitwise equality.
  return kC * static_cast<double>(n + 1) * kFloatUnitRoundoff * ref.sum_abs +
         std::numeric_limits<float>::min();
}

// Exercises every tail-boundary residue class for a kernel whose main loop is
// 8 wide with a 4-wide step and a <4 scalar tail. Sweeping 0..72 contiguously
// covers every (main-loop iterations, 4-step taken, tail length) combination,
// which a hand-picked list of lengths does not guarantee.
std::vector<std::size_t> tail_boundary_lengths() {
  std::vector<std::size_t> lengths;
  for (std::size_t n = 0; n <= 72; ++n) lengths.push_back(n);
  for (std::size_t n : {100u, 127u, 128u, 129u, 200u, 255u, 256u, 257u, 1000u}) {
    lengths.push_back(n);
  }
  return lengths;
}

// Every compiled kernel, so a test failure names the ISA rather than "the SIMD
// path". The dispatcher is checked separately, in test_isa_dispatch_routes.
void check_gather_kernels(const float* dense, const std::int32_t* idx,
                          const float* val, std::size_t nnz, const char* what) {
  const DotReference ref = gather_reference(dense, idx, val, nnz);
  const double budget = dot_error_budget(ref, nnz);

  const double scalar_err =
      std::fabs(static_cast<double>(
                    scrna::simd::gather_dot_scalar(dense, idx, val, nnz)) - ref.exact);
  if (!(scalar_err <= budget)) {
    std::fprintf(stderr, "  scalar gather_dot exceeded budget (%s, nnz=%zu): "
                         "err=%.3e budget=%.3e\n", what, nnz, scalar_err, budget);
    ++g_failures;
  }

#if defined(SCRNA_BUILD_ARM_NEON)
  const double neon_err =
      std::fabs(static_cast<double>(
                    scrna::simd::gather_dot_neon(dense, idx, val, nnz)) - ref.exact);
  if (!(neon_err <= budget)) {
    std::fprintf(stderr, "  NEON gather_dot exceeded budget (%s, nnz=%zu): "
                         "err=%.3e budget=%.3e\n", what, nnz, neon_err, budget);
    ++g_failures;
  }
#endif
}

// NEON gather kernel: exhaustive tail coverage against the double reference.
void test_neon_gather_tail_boundaries() {
  std::mt19937 rng(918273);
  std::uniform_real_distribution<float> vdist(-2.0f, 2.0f);
  const std::size_t n_cols = 4096;

  std::vector<float> dense(n_cols);
  for (float& d : dense) d = vdist(rng);

  for (const std::size_t nnz : tail_boundary_lengths()) {
    if (nnz > n_cols) continue;
    // Strictly increasing indices, as Block_CSR guarantees.
    std::vector<std::int32_t> idx;
    std::vector<float> val;
    std::int32_t c = 0;
    for (std::size_t i = 0; i < nnz; ++i) {
      idx.push_back(c);
      val.push_back(vdist(rng));
      c += 1 + static_cast<std::int32_t>(rng() % 3);
      if (c >= static_cast<std::int32_t>(n_cols)) break;
    }
    check_gather_kernels(dense.data(), idx.data(), val.data(), idx.size(), "uniform");
  }
}

// Numerically hostile inputs. Well-conditioned random data hides sign and
// ordering bugs that these surface.
void test_neon_gather_adversarial_values() {
  const std::size_t n_cols = 512;
  std::vector<float> dense(n_cols, 0.0f);

  // (a) Massive cancellation: alternating +/- terms summing to ~0. A dropped
  //     or duplicated tail element shows up enormously here.
  for (std::size_t nnz : {7u, 8u, 9u, 15u, 16u, 17u, 63u, 64u, 65u}) {
    std::vector<std::int32_t> idx(nnz);
    std::vector<float> val(nnz);
    for (std::size_t i = 0; i < nnz; ++i) {
      idx[i] = static_cast<std::int32_t>(i);
      dense[i] = 1.0e6f;
      val[i] = (i % 2 == 0) ? 1.0f : -1.0f;
    }
    check_gather_kernels(dense.data(), idx.data(), val.data(), nnz, "cancellation");
  }

  // (b) Wide dynamic range: one huge term among many tiny ones, so the small
  //     contributions are at the edge of being absorbed.
  for (std::size_t nnz : {8u, 12u, 16u, 33u}) {
    std::vector<std::int32_t> idx(nnz);
    std::vector<float> val(nnz);
    for (std::size_t i = 0; i < nnz; ++i) {
      idx[i] = static_cast<std::int32_t>(i);
      dense[i] = (i == 0) ? 1.0e8f : 1.0e-4f;
      val[i] = (i == 0) ? 1.0e8f : 1.0e-4f;
    }
    check_gather_kernels(dense.data(), idx.data(), val.data(), nnz, "dynamic range");
  }

  // (c) Subnormals and zeros. These take a different path in some FPUs, and
  //     Apple Silicon flushes-to-zero in some modes -- worth pinning that both
  //     kernels agree about it rather than assuming.
  for (std::size_t nnz : {4u, 8u, 13u}) {
    std::vector<std::int32_t> idx(nnz);
    std::vector<float> val(nnz);
    for (std::size_t i = 0; i < nnz; ++i) {
      idx[i] = static_cast<std::int32_t>(i);
      dense[i] = (i % 3 == 0) ? 0.0f : std::numeric_limits<float>::denorm_min();
      val[i] = (i % 2 == 0) ? 1.0f : std::numeric_limits<float>::denorm_min();
    }
    check_gather_kernels(dense.data(), idx.data(), val.data(), nnz, "subnormal");
  }

  // (d) Gather ORDER. Every value equal makes a permutation invisible; distinct
  //     values with shuffled (still strictly increasing) indices do not. This is
  //     what catches a transposed gather, which is the single easiest mistake to
  //     make in a hand-written gather kernel.
  {
    const std::size_t nnz = 32;
    std::vector<std::int32_t> idx(nnz);
    std::vector<float> val(nnz);
    for (std::size_t i = 0; i < nnz; ++i) {
      idx[i] = static_cast<std::int32_t>(i * 7);          // strictly increasing
      dense[idx[i]] = static_cast<float>(i + 1);          // distinct per lane
      val[i] = static_cast<float>((i % 5) + 1);
    }
    check_gather_kernels(dense.data(), idx.data(), val.data(), nnz, "gather order");
  }
}

// NEON dense kernel (the hnswlib space adapter's hot path).
void test_neon_dense_tail_boundaries() {
  std::mt19937 rng(5150);
  std::uniform_real_distribution<float> vdist(-1.5f, 1.5f);

  for (const std::size_t n : tail_boundary_lengths()) {
    std::vector<float> a(n), b(n);
    for (std::size_t i = 0; i < n; ++i) { a[i] = vdist(rng); b[i] = vdist(rng); }

    const DotReference ref = dense_reference(a.data(), b.data(), n);
    const double budget = dot_error_budget(ref, n);

    const double scalar_err = std::fabs(
        static_cast<double>(scrna::simd::dense_dot_scalar(a.data(), b.data(), n)) - ref.exact);
    SCRNA_CHECK(scalar_err <= budget);

#if defined(SCRNA_BUILD_ARM_NEON)
    const double neon_err = std::fabs(
        static_cast<double>(scrna::simd::dense_dot_neon(a.data(), b.data(), n)) - ref.exact);
    SCRNA_CHECK(neon_err <= budget);
#endif
  }
}

// A single sequential FLOAT accumulator -- what gather_dot_scalar was before
// the wide-accumulation refactor, kept here as the peer the NEON kernel is
// compared against.
//
// It has to live in the test rather than in the header now. The library's
// scalar kernel accumulates in double, so measuring NEON against it no longer
// answers "does vectorising cost accuracy?" -- it answers "is float worse than
// double?", to which the answer is yes, always, by construction, and the test
// would be pinning nothing. The property below is about SUMMATION TREE SHAPE at
// a fixed precision, so both sides of the comparison must be at that precision.
float gather_dot_float_acc(const float* dense, const std::int32_t* idx,
                           const float* val, std::size_t nnz) {
  float acc = 0.0f;  // NOLINT(scrna-wide-accumulation) -- see above
  for (std::size_t i = 0; i < nnz; ++i) acc += dense[idx[i]] * val[i];
  return acc;
}

// Vectorising should not COST accuracy. Two accumulators plus a horizontal
// reduction is a shallower summation tree than one sequential accumulator, so
// the vector path is expected to be at least as accurate on average.
//
// Asserted in aggregate over many trials rather than per-trial: any individual
// case can go either way on rounding luck, but a systematic regression (say,
// accumulating the tail into the wrong place) moves the mean.
void test_neon_accuracy_not_worse_than_scalar() {
#if defined(SCRNA_BUILD_ARM_NEON)
  std::mt19937 rng(31337);
  std::uniform_real_distribution<float> vdist(-1.0f, 1.0f);
  const std::size_t n_cols = 2048;
  std::vector<float> dense(n_cols);
  for (float& d : dense) d = vdist(rng);

  double scalar_total = 0.0, neon_total = 0.0, wide_total = 0.0;
  const int kTrials = 200;
  for (int t = 0; t < kTrials; ++t) {
    const std::size_t nnz = 64 + (rng() % 512);
    std::vector<std::int32_t> idx;
    std::vector<float> val;
    std::int32_t c = 0;
    for (std::size_t i = 0; i < nnz && c < static_cast<std::int32_t>(n_cols); ++i) {
      idx.push_back(c);
      val.push_back(vdist(rng));
      c += 1 + static_cast<std::int32_t>(rng() % 3);
    }
    const DotReference ref = gather_reference(dense.data(), idx.data(), val.data(), idx.size());
    scalar_total += std::fabs(
        static_cast<double>(gather_dot_float_acc(
            dense.data(), idx.data(), val.data(), idx.size())) - ref.exact);
    neon_total += std::fabs(
        static_cast<double>(scrna::simd::gather_dot_neon(
            dense.data(), idx.data(), val.data(), idx.size())) - ref.exact);
    wide_total += std::fabs(
        static_cast<double>(scrna::simd::gather_dot_scalar(
            dense.data(), idx.data(), val.data(), idx.size())) - ref.exact);
  }
  // Generous factor: the claim is "not systematically worse", not "always
  // better". A real regression blows past 1.5x.
  SCRNA_CHECK(neon_total <= scalar_total * 1.5);

  // And the other direction, which is what the wide-accumulation refactor
  // bought: the library's scalar kernel now accumulates in double, so it must
  // be strictly better than either float path. Without this assertion a
  // revert to `float acc` would sail through -- the test above would simply
  // start comparing two float kernels again and still pass.
  SCRNA_CHECK(wide_total < scalar_total);
  SCRNA_CHECK(wide_total < neon_total);

  std::printf("  accumulated |error| over %d trials: float-acc=%.3e neon=%.3e "
              "wide-acc=%.3e\n", kTrials, scalar_total, neon_total, wide_total);
#endif
}

// The dispatcher must actually route to the kernel it names. Without this, a
// mis-wired switch would leave every kernel test passing while production ran
// the scalar path -- which is exactly the failure mode that made the NEON gap
// invisible in the first place.
void test_isa_dispatch_routes() {
  std::mt19937 rng(24680);
  std::uniform_real_distribution<float> vdist(-2.0f, 2.0f);
  const std::size_t n_cols = 256, nnz = 37;  // 37 = 4*8 + 5: main loop, 4-step, tail

  std::vector<float> dense(n_cols);
  for (float& d : dense) d = vdist(rng);
  std::vector<std::int32_t> idx(nnz);
  std::vector<float> val(nnz);
  for (std::size_t i = 0; i < nnz; ++i) {
    idx[i] = static_cast<std::int32_t>(i * 3);
    val[i] = vdist(rng);
  }

  const DotReference ref = gather_reference(dense.data(), idx.data(), val.data(), nnz);
  const double budget = dot_error_budget(ref, nnz);
  const double dispatched = static_cast<double>(
      scrna::simd::gather_dot(dense.data(), idx.data(), val.data(), nnz));
  SCRNA_CHECK(std::fabs(dispatched - ref.exact) <= budget);

  // isa_name() must describe the path that was actually taken. Run the suite
  // under SCRNA_FORCE_ISA=scalar and =neon to cover both directions; CI does.
  const std::string isa = scrna::simd::isa_name();
#if defined(SCRNA_BUILD_ARM_NEON)
  SCRNA_CHECK(isa == "neon" || isa == "scalar");
#endif
  std::printf("  dispatch reports isa=%s\n", isa.c_str());
}

// The dense dot dispatcher (used by the hnswlib space adapter) must agree with
// its scalar reference across the 8-lane and 16-lane boundaries.
void test_dense_dot_parity() {
  std::mt19937 rng(2024);
  std::uniform_real_distribution<float> vdist(-1.5f, 1.5f);
  for (const std::size_t n : {std::size_t{0}, std::size_t{1}, std::size_t{7}, std::size_t{8},
                               std::size_t{9}, std::size_t{15}, std::size_t{16}, std::size_t{17},
                               std::size_t{31}, std::size_t{32}, std::size_t{33},
                               std::size_t{50}, std::size_t{257}}) {
    std::vector<float> a(n), b(n);
    for (std::size_t i = 0; i < n; ++i) { a[i] = vdist(rng); b[i] = vdist(rng); }
    const float got = scrna::simd::dense_dot(a.data(), b.data(), n);
    const float want = scrna::simd::dense_dot_scalar(a.data(), b.data(), n);
    SCRNA_CHECK(near(got, want, 1e-4f));
  }
}

// Reference implementation: full scoring + full sort, no heap. The production
// path must match it exactly, including tie-breaking.
std::vector<std::pair<float, std::int32_t>> reference_topk(const Block_CSR<float>& mat,
                                                            std::size_t i, std::size_t k) {
  std::vector<std::pair<float, std::int32_t>> scored;
  for (std::size_t j = 0; j < mat.rows(); ++j) {
    if (j == i) continue;
    const float sim = scrna::simd::sparse_cosine_similarity<float, std::int32_t>(
        mat.row_indices(i), mat.row_values(i), mat.row_indices(j), mat.row_values(j));
    scored.emplace_back(sim, static_cast<std::int32_t>(j));
  }
  std::sort(scored.begin(), scored.end(), scrna::detail::BetterNeighbor{});
  if (scored.size() > k) scored.resize(k);
  return scored;
}

Block_CSR<float> make_random_matrix(std::size_t n_rows, std::size_t n_cols,
                                     std::size_t nnz_per_row, std::mt19937& rng) {
  AlignedVector<float> values;
  AlignedVector<std::int32_t> col_idx;
  AlignedVector<std::int32_t> row_ptr;
  row_ptr.push_back(0);
  std::uniform_real_distribution<float> vdist(0.1f, 3.0f);

  for (std::size_t r = 0; r < n_rows; ++r) {
    // Sample strictly increasing column indices.
    std::vector<std::int32_t> cols;
    for (std::size_t c = 0; c < n_cols && cols.size() < nnz_per_row; ++c) {
      if (rng() % n_cols < nnz_per_row) cols.push_back(static_cast<std::int32_t>(c));
    }
    if (cols.empty()) cols.push_back(static_cast<std::int32_t>(rng() % n_cols));
    for (const std::int32_t c : cols) {
      col_idx.push_back(c);
      values.push_back(vdist(rng));
    }
    row_ptr.push_back(static_cast<std::int32_t>(values.size()));
  }
  return Block_CSR<float>(n_rows, n_cols, std::move(values), std::move(col_idx),
                           std::move(row_ptr));
}

// Finding 10 (2026-09-16 audit): build_knn_graph_bruteforce now sweeps only
// the upper triangle of the n x n similarity matrix, computing each unordered
// pair once and applying it to BOTH rows' top-k heaps, scheduled by the
// classic round-robin ("circle method") tournament algorithm over row tiles
// so that concurrent tile-pairs never touch the same row. This is the
// scheduler's own correctness property, independent of the similarity math:
// every unordered pair of tiles appears in exactly one round, and no tile
// repeats within a round.
void test_round_robin_schedule_covers_every_pair_exactly_once() {
  for (const std::size_t num_tiles : {std::size_t{0}, std::size_t{1}, std::size_t{2},
                                        std::size_t{3}, std::size_t{4}, std::size_t{5},
                                        std::size_t{6}, std::size_t{7}, std::size_t{8},
                                        std::size_t{11}, std::size_t{16}}) {
    const auto rounds = scrna::detail::round_robin_schedule(num_tiles);

    std::vector<std::vector<bool>> seen(num_tiles, std::vector<bool>(num_tiles, false));
    std::size_t pair_count = 0;
    for (const auto& pairs : rounds) {
      std::vector<bool> touched(num_tiles, false);
      for (const auto& [a, b] : pairs) {
        SCRNA_CHECK(a < b);
        SCRNA_CHECK(b < num_tiles);
        // No tile repeats within a round: the property that makes a round
        // safe to run fully in parallel with no lock on per-tile state.
        SCRNA_CHECK(!touched[a]);
        SCRNA_CHECK(!touched[b]);
        touched[a] = touched[b] = true;
        // No pair is scheduled twice, in this round or any other.
        SCRNA_CHECK(!seen[a][b]);
        seen[a][b] = true;
        ++pair_count;
      }
    }
    // Every one of the C(num_tiles, 2) unordered pairs was scheduled exactly
    // once: pair_count already rules out duplicates (checked above), so
    // matching the total also rules out omissions.
    const std::size_t expected = num_tiles < 2 ? 0 : num_tiles * (num_tiles - 1) / 2;
    SCRNA_CHECK(pair_count == expected);
  }
}

// The scheduler's parallel-safety property translates, through the sweep, into
// an exact-equality guarantee the old two-computations-per-pair arrangement
// never had: sim(i,j) and sim(j,i) are not merely close, they are the SAME
// float, because there is only one computation. Checked directly on the graph
// output rather than assumed from the scheduling test above, since it is the
// property callers (and the ULP/tie-break guardrails around this change)
// actually depend on.
void test_symmetric_sweep_is_strictly_symmetric() {
  std::mt19937 rng(20260916);
  for (const std::size_t n : {std::size_t{2}, std::size_t{3}, std::size_t{9}, std::size_t{50}}) {
    auto mat = make_random_matrix(n, 40, 10, rng);
    const std::size_t k = std::min<std::size_t>(6, n - 1);
    auto graph = scrna::build_knn_graph_bruteforce(mat, k);

    std::size_t mutual_pairs = 0;
    for (std::size_t i = 0; i < n; ++i) {
      for (std::size_t r = 0; r < k; ++r) {
        const std::int32_t j = graph.neighbor_idx[i * k + r];
        if (j < 0) continue;
        const float sim_ij = graph.neighbor_sim[i * k + r];
        // Is i ALSO in j's own list? Cosine similarity is symmetric, but
        // being each other's top-k neighbour is not -- j may have k closer
        // rows than i. Only check the pairs where the relation IS mutual.
        for (std::size_t s = 0; s < k; ++s) {
          if (graph.neighbor_idx[static_cast<std::size_t>(j) * k + s] ==
              static_cast<std::int32_t>(i)) {
            const float sim_ji = graph.neighbor_sim[static_cast<std::size_t>(j) * k + s];
            SCRNA_CHECK_MSG(sim_ij == sim_ji,
                            "n=" + std::to_string(n) + " i=" + std::to_string(i) +
                                " j=" + std::to_string(j) + " sim_ij=" + std::to_string(sim_ij) +
                                " sim_ji=" + std::to_string(sim_ji));
            ++mutual_pairs;
            break;
          }
        }
      }
    }
    // The check above is vacuous if no pair ever turns out mutual -- with
    // random dense-ish rows and k this close to n-1, that should not happen.
    SCRNA_CHECK(mutual_pairs > 0);
  }
}

// Extends test_knn_matches_reference (below) across several seeds, sizes and
// k, against the SAME independent naive reference (full score + full sort,
// no heap, no tiling). The round-robin sweep changes summation order
// relative to the old per-row arrangement (a tile-pair's gather order is
// fixed by whichever side is scattered, not by "the row currently being
// resolved"), so this is where a genuine near-tie rank swap would show up;
// none does across this sweep of seeds; see this change's commit message for
// the larger-scale before/after comparison against the pre-change binary.
void test_knn_matches_reference_multi_seed() {
  const std::size_t seeds[] = {1, 2, 3, 42, 987, 20260916};
  for (const std::size_t seed : seeds) {
    std::mt19937 rng(static_cast<unsigned>(seed));
    const std::size_t n = 30 + (seed % 40);
    const std::size_t n_cols = 50;
    const std::size_t k = 5;
    auto mat = make_random_matrix(n, n_cols, 10, rng);
    auto graph = scrna::build_knn_graph_bruteforce(mat, k);
    for (std::size_t i = 0; i < n; ++i) {
      const auto expected = reference_topk(mat, i, k);
      for (std::size_t r = 0; r < k; ++r) {
        const std::int32_t got_idx = graph.neighbor_idx[i * k + r];
        const float got_sim = graph.neighbor_sim[i * k + r];
        if (r < expected.size()) {
          SCRNA_CHECK(got_idx == expected[r].second);
          SCRNA_CHECK(near(got_sim, expected[r].first, 1e-4f));
        } else {
          SCRNA_CHECK(got_idx == -1);
        }
      }
    }
  }
}

void test_knn_graph_bruteforce() {
  auto mat = make_toy_matrix();
  auto graph = scrna::build_knn_graph_bruteforce(mat, /*k=*/2);
  SCRNA_CHECK(graph.k == 2);
  SCRNA_CHECK(graph.neighbor_idx.size() == 3 * 2);
  // Every row should get at least one non-negative neighbor id.
  for (std::size_t i = 0; i < 3; ++i) {
    SCRNA_CHECK(graph.neighbor_idx[i * 2] != -1);
  }
}

// The bounded-heap rewrite must produce exactly the same graph the naive
// full-sort reference produces.
void test_knn_matches_reference() {
  std::mt19937 rng(987);
  const std::size_t n = 60, n_cols = 80, k = 7;
  auto mat = make_random_matrix(n, n_cols, 12, rng);
  auto graph = scrna::build_knn_graph_bruteforce(mat, k);

  SCRNA_CHECK(graph.neighbor_idx.size() == n * k);

  for (std::size_t i = 0; i < n; ++i) {
    const auto expected = reference_topk(mat, i, k);
    for (std::size_t r = 0; r < k; ++r) {
      const std::int32_t got_idx = graph.neighbor_idx[i * k + r];
      const float got_sim = graph.neighbor_sim[i * k + r];
      if (r < expected.size()) {
        SCRNA_CHECK(got_idx == expected[r].second);
        SCRNA_CHECK(near(got_sim, expected[r].first, 1e-4f));
      } else {
        SCRNA_CHECK(got_idx == -1);
      }
    }
  }
}

// k larger than the number of available neighbors must leave sentinels, not
// garbage or out-of-range writes.
void test_knn_k_exceeds_rows() {
  auto mat = make_toy_matrix();
  const std::size_t k = 10;  // only 2 other rows exist
  auto graph = scrna::build_knn_graph_bruteforce(mat, k);
  SCRNA_CHECK(graph.neighbor_idx.size() == 3 * k);
  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t r = 2; r < k; ++r) {
      SCRNA_CHECK(graph.neighbor_idx[i * k + r] == -1);
      SCRNA_CHECK(graph.neighbor_sim[i * k + r] == 0.0f);
    }
  }
}

// Repeated runs must agree exactly, including under OpenMP scheduling.
void test_knn_deterministic() {
  std::mt19937 rng(555);
  auto mat = make_random_matrix(40, 50, 8, rng);
  auto a = scrna::build_knn_graph_bruteforce(mat, 5);
  auto b = scrna::build_knn_graph_bruteforce(mat, 5);
  SCRNA_CHECK(a.neighbor_idx == b.neighbor_idx);
  for (std::size_t i = 0; i < a.neighbor_sim.size(); ++i) {
    SCRNA_CHECK(a.neighbor_sim[i] == b.neighbor_sim[i]);
  }
}

void test_knn_edge_cases() {
  auto mat = make_toy_matrix();
  auto g0 = scrna::build_knn_graph_bruteforce(mat, 0);
  SCRNA_CHECK(g0.neighbor_idx.empty());

  Block_CSR<float> empty;
  auto ge = scrna::build_knn_graph_bruteforce(empty, 5);
  SCRNA_CHECK(ge.neighbor_idx.empty());
}

// A k-NN graph built through the adopted (zero-copy) path must match one built
// through the copying path.
void test_knn_on_adopted_matrix() {
  std::mt19937 rng(77);
  auto owned = make_random_matrix(30, 40, 6, rng);

  struct NoopOwner {};
  auto adopted = Block_CSR<float>::adopt(owned.rows(), owned.cols(), owned.values_span(),
                                          owned.col_idx_span(), owned.row_ptr_span(),
                                          NoopOwner{});
  SCRNA_CHECK(adopted.is_zero_copy());

  auto ga = scrna::build_knn_graph_bruteforce(owned, 4);
  auto gb = scrna::build_knn_graph_bruteforce(adopted, 4);
  SCRNA_CHECK(ga.neighbor_idx == gb.neighbor_idx);
}

#if defined(SCRNA_ENABLE_HNSW)

// HNSW is approximate, so it is checked against exact brute force by recall
// rather than by equality.
void test_hnsw_recall_vs_bruteforce() {
  std::mt19937 rng(4242);
  const std::size_t n = 400, n_cols = 64, k = 10;
  auto mat = make_random_matrix(n, n_cols, 20, rng);

  const auto exact = scrna::build_knn_graph_bruteforce(mat, k);

  scrna::HnswParams p;
  p.ef_construction = 400;
  p.ef_search = 200;
  const auto approx = scrna::hnsw_build_and_search(mat, k, p);

  SCRNA_CHECK(approx.k == k);
  SCRNA_CHECK(approx.neighbor_idx.size() == n * k);

  std::size_t hits = 0, total = 0;
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t r = 0; r < k; ++r) {
      const std::int32_t got = approx.neighbor_idx[i * k + r];
      // Never return the query point itself.
      SCRNA_CHECK(got != static_cast<std::int32_t>(i));
      if (got < 0) continue;
      ++total;
      for (std::size_t s = 0; s < k; ++s) {
        if (exact.neighbor_idx[i * k + s] == got) { ++hits; break; }
      }
    }
  }
  const double recall = total ? static_cast<double>(hits) / static_cast<double>(total) : 0.0;
  std::printf("  hnsw recall@%zu vs exact: %.3f\n", k, recall);
  SCRNA_CHECK(recall > 0.90);

  // Similarities must be sane and in descending order per row.
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t r = 1; r < k; ++r) {
      if (approx.neighbor_idx[i * k + r] < 0) continue;
      SCRNA_CHECK(approx.neighbor_sim[i * k + r] <= approx.neighbor_sim[i * k + r - 1] + 1e-4f);
    }
    if (approx.neighbor_idx[i * k] >= 0) {
      const float s = approx.neighbor_sim[i * k];
      SCRNA_CHECK(s >= -1.01f && s <= 1.01f);
    }
  }
}

// The dense entry point (the one that scales) must agree with the sparse one.
void test_hnsw_dense_matches_sparse() {
  std::mt19937 rng(99);
  const std::size_t n = 200, dim = 32, k = 8;
  auto mat = make_random_matrix(n, dim, 12, rng);

  // Materialise the same data densely.
  std::vector<float> dense(n * dim, 0.0f);
  for (std::size_t i = 0; i < n; ++i) {
    const auto idx = mat.row_indices(i);
    const auto val = mat.row_values(i);
    for (std::size_t t = 0; t < idx.size(); ++t) dense[i * dim + idx[t]] = val[t];
  }

  scrna::HnswParams p;
  p.ef_construction = 300;
  p.ef_search = 150;
  p.num_threads = 1;  // determinism for a like-for-like comparison

  const auto from_sparse = scrna::hnsw_build_and_search(mat, k, p);
  const auto from_dense = scrna::hnsw_build_and_search_dense(dense.data(), n, dim, k, p);

  SCRNA_CHECK(from_sparse.neighbor_idx == from_dense.neighbor_idx);
}

// Finding 19 (2026-09-16 audit): hnsw_build_and_search{,_dense} used to build
// and sentinel-fill a full n*k KnnGraph unconditionally, then discard it on
// every call that reached index.query_{dense,csr}() -- which builds and
// returns its OWN complete graph. The early-exit graph now exists only inside
// the n==0/k==0/dim==0 branch; this pins that branch still returns exactly
// the same sentinel-filled shape for hnsw_build_and_search_dense specifically
// (test_hnsw_edge_cases above already covers the CSR entry point's k==0 case).
void test_hnsw_dense_early_exit_shapes() {
  std::vector<float> dense(30, 1.0f);  // 3 rows x 10 dims, content irrelevant
  scrna::HnswParams p;

  const auto g_n0 = scrna::hnsw_build_and_search_dense(dense.data(), 0, 10, 5, p);
  SCRNA_CHECK(g_n0.k == 5);
  SCRNA_CHECK(g_n0.neighbor_idx.empty());
  SCRNA_CHECK(g_n0.neighbor_sim.empty());

  const auto g_k0 = scrna::hnsw_build_and_search_dense(dense.data(), 3, 10, 0, p);
  SCRNA_CHECK(g_k0.k == 0);
  SCRNA_CHECK(g_k0.neighbor_idx.empty());
  SCRNA_CHECK(g_k0.neighbor_sim.empty());

  const auto g_dim0 = scrna::hnsw_build_and_search_dense(dense.data(), 3, 0, 5, p);
  SCRNA_CHECK(g_dim0.k == 5);
  SCRNA_CHECK(g_dim0.neighbor_idx.size() == 3 * 5);
  for (const std::int32_t v : g_dim0.neighbor_idx) SCRNA_CHECK(v == -1);
  for (const float v : g_dim0.neighbor_sim) SCRNA_CHECK(v == 0.0f);
}

// The memory guard must refuse rather than attempt a huge allocation.
void test_hnsw_memory_guard() {
  auto mat = make_toy_matrix();
  scrna::HnswParams p;
  p.max_index_bytes = 16;  // absurdly small
  SCRNA_CHECK_THROWS(scrna::hnsw_build_and_search(mat, 2, p));

  // And the estimate should scale with dimensionality as documented.
  const std::size_t small = scrna::hnsw_index_bytes(1000, 50);
  const std::size_t big = scrna::hnsw_index_bytes(1000, 2000);
  SCRNA_CHECK(big > small * 10);
}

void test_hnsw_edge_cases() {
  auto mat = make_toy_matrix();
  scrna::HnswParams p;
  auto g0 = scrna::hnsw_build_and_search(mat, 0, p);
  SCRNA_CHECK(g0.neighbor_idx.empty());

  // k larger than the population: extra slots stay sentinel.
  auto gk = scrna::hnsw_build_and_search(mat, 10, p);
  SCRNA_CHECK(gk.neighbor_idx.size() == 3 * 10);
  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t r = 2; r < 10; ++r) {
      SCRNA_CHECK(gk.neighbor_idx[i * 10 + r] == -1);
    }
  }
}

// Both spaces must produce equivalent results; they differ only in which SIMD
// implementation computes the inner product.
void test_hnsw_space_equivalence() {
  std::mt19937 rng(31337);
  const std::size_t n = 150, dim = 48, k = 6;
  auto mat = make_random_matrix(n, dim, 15, rng);

  scrna::HnswParams a;
  a.num_threads = 1; a.ef_construction = 300; a.ef_search = 150; a.use_scrna_space = true;
  scrna::HnswParams b = a;
  b.use_scrna_space = false;

  const auto ga = scrna::hnsw_build_and_search(mat, k, a);
  const auto gb = scrna::hnsw_build_and_search(mat, k, b);
  SCRNA_CHECK(ga.neighbor_idx == gb.neighbor_idx);
}

#endif  // SCRNA_ENABLE_HNSW

}  // namespace

// ------------------------------------------------------- security regressions
// One test per closed finding in docs/REVIEW_2026-08-15.md.

// REVIEW finding 4 (Medium): `n * k` sizes the result and indexes every write.
// Unchecked it wraps -- k = 2**62 over 3 rows gives n*k == 0, both assigns
// succeed against an empty vector, and the writes at i*k+r then land outside it.
void test_knn_k_overflow_is_rejected() {
  auto mat = make_toy_matrix();  // 3 rows

  // 3 * 2**63 wraps to 2**63; must throw rather than allocate against it.
  // (2**62 is deliberately not used here: 3 * 2**62 does NOT overflow, so it is
  // an honest allocation failure, not this defect.)
  SCRNA_CHECK_THROWS(scrna::build_knn_graph_bruteforce(mat, std::size_t{1} << 63));
  SCRNA_CHECK_THROWS(
      scrna::build_knn_graph_bruteforce(mat, std::numeric_limits<std::size_t>::max()));
  SCRNA_CHECK_THROWS(
      scrna::build_knn_graph_bruteforce(mat, std::numeric_limits<std::size_t>::max() / 2));

  // k = 0 is a legal no-op and must NOT be caught by the overflow guard (the
  // division by k in the check has to be predicated on k != 0).
  auto empty = scrna::build_knn_graph_bruteforce(mat, 0);
  SCRNA_CHECK(empty.k == 0);
  SCRNA_CHECK(empty.neighbor_idx.empty());

  // The largest k that does not overflow must still be accepted, so the guard
  // is not merely rejecting everything large.
  auto ok = scrna::build_knn_graph_bruteforce(mat, 2);
  SCRNA_CHECK(ok.neighbor_idx.size() == 6);
}

// REVIEW finding 5 (Medium): three norm loops accumulated in float, so any
// element around 1e19 or above saturated norm_sq to +inf. inf > 0 passes the
// zero-guard, 1/sqrt(inf) is 0, and the row was written out as an all-zero
// vector that matches nothing -- from finite input, with no diagnostic.
void test_norms_do_not_saturate_on_large_magnitudes() {
  // The residual this test used to pin -- two identical 1e30-scaled rows
  // silently scoring 0 because simd::gather_dot accumulates in float and
  // saturates to +inf -- is now unreachable: Block_CSR::validate() rejects any
  // row whose squared norm exceeds 1e34, which is the point below which the
  // float accumulator provably cannot saturate (Cauchy-Schwarz, four orders of
  // headroom under FLT_MAX). Measured cliff: exact top-15 agreement with a
  // double reference at value scale 1e17, every row wrong at 1e18.
  {
    const float big = 1e30f;
    AlignedVector<float> values = {big, big, big, big, 1.0f};
    AlignedVector<std::int32_t> col_idx = {0, 2, 0, 2, 1};
    AlignedVector<std::int32_t> row_ptr = {0, 2, 4, 5};
    bool threw = false;
    try {
      Block_CSR<float> mat(3, 4, std::move(values), std::move(col_idx), std::move(row_ptr));
    } catch (const std::invalid_argument&) {
      threw = true;
    }
    SCRNA_CHECK(threw);
  }

  // Large-but-representable magnitudes stay exact and stay ranked correctly.
  // 1e15 gives a squared norm of 2e30, comfortably under the ceiling.
  {
    const float big = 1e15f;
    AlignedVector<float> values = {big, big, big, big, 1.0f};
    AlignedVector<std::int32_t> col_idx = {0, 2, 0, 2, 1};
    AlignedVector<std::int32_t> row_ptr = {0, 2, 4, 5};
    Block_CSR<float> mat(3, 4, std::move(values), std::move(col_idx), std::move(row_ptr));

    auto graph = scrna::build_knn_graph_bruteforce(mat, 2);
    for (float s : graph.neighbor_sim) SCRNA_CHECK(std::isfinite(s));
    for (std::int32_t v : graph.neighbor_idx) SCRNA_CHECK(v >= -1 && v < 3);
    // Identical directions now score 1.0 instead of silently scoring 0.
    SCRNA_CHECK(graph.neighbor_idx[0] == 1);
    SCRNA_CHECK(near(graph.neighbor_sim[0], 1.0f, 1e-4f));
    SCRNA_CHECK(graph.neighbor_idx[1] == 2);
    SCRNA_CHECK(near(graph.neighbor_sim[1], 0.0f, 1e-4f));
  }

  // The ceiling must be a ceiling, not a blanket ban on large values: a single
  // row at the boundary is accepted.
  {
    const float v = 1e16f;  // norm^2 = 1e32 < 1e34
    AlignedVector<float> values = {v};
    AlignedVector<std::int32_t> col_idx = {0};
    AlignedVector<std::int32_t> row_ptr = {0, 1};
    bool threw = false;
    try {
      Block_CSR<float> mat(1, 4, std::move(values), std::move(col_idx), std::move(row_ptr));
    } catch (const std::invalid_argument&) {
      threw = true;
    }
    SCRNA_CHECK(!threw);
  }

  const float big = 1e30f;  // still exercised through the dense helpers below

#if defined(SCRNA_ENABLE_HNSW)
  // The shared helper directly: a 1e30-scaled row must normalise to a unit
  // vector, not to all zeros. Squaring in float gives 1e60 -> +inf -> inv 0.
  {
    const std::int32_t idx[2] = {0, 2};
    const float val[2] = {big, big};
    float out[4] = {0.0f, 0.0f, 0.0f, 0.0f};
    scrna::detail::densify_normalized(std::span<const std::int32_t>(idx, 2),
                                      std::span<const float>(val, 2), out);
    SCRNA_CHECK(out[0] != 0.0f);
    SCRNA_CHECK(near(out[0], 1.0f / std::sqrt(2.0f)));
    SCRNA_CHECK(near(out[0] * out[0] + out[2] * out[2], 1.0f));
  }

  // Dense path, both check_finite settings.
  {
    const float row[3] = {big, 0.0f, big};
    float out[3] = {0.0f, 0.0f, 0.0f};
    scrna::detail::normalize_dense_row(row, 3, out, /*check_finite=*/true, 0);
    SCRNA_CHECK(near(out[0] * out[0] + out[2] * out[2], 1.0f));
    scrna::detail::normalize_dense_row(row, 3, out, /*check_finite=*/false, 0);
    SCRNA_CHECK(near(out[0] * out[0] + out[2] * out[2], 1.0f));
  }

  // A genuinely all-zero row still normalises to zero rather than NaN, and a
  // non-finite row is still rejected on the build pass.
  {
    const float zero[3] = {0.0f, 0.0f, 0.0f};
    float out[3] = {1.0f, 1.0f, 1.0f};
    scrna::detail::normalize_dense_row(zero, 3, out, /*check_finite=*/true, 0);
    SCRNA_CHECK(out[0] == 0.0f && out[1] == 0.0f && out[2] == 0.0f);

    const float nan_row[2] = {1.0f, std::numeric_limits<float>::quiet_NaN()};
    float out2[2] = {0.0f, 0.0f};
    SCRNA_CHECK_THROWS(
        (scrna::detail::normalize_dense_row(nan_row, 2, out2, /*check_finite=*/true, 7), 0));
  }
#endif  // SCRNA_ENABLE_HNSW
}


// Regression: cosine_from_dot could return +inf.
//
// The guard checked isfinite() on the WIDE quotient and then returned it cast
// to T. WideAcc<float> is double, so a quotient that is perfectly finite in
// double can overflow on the narrowing cast:
//
//     cosine_from_dot(1.0f, FLT_DENORM_MIN, 1.0f)
//         r = 1.0 / 1.4e-45 = 7.1e44   -- finite as a double
//         (float)r          = +inf     -- returned to the caller
//
// That defeats the guard's documented contract, and an infinite score is worse
// for the caller than a NaN: it is a valid, totally-ordered value that compares
// greater than every real similarity, so instead of tripping a sort_heap
// precondition it silently wins the top-k and makes a degenerate row the
// nearest neighbour of everything it touches.
//
// Not reachable through a Block_CSR built by this library -- from_raw() rejects
// non-finite values and caps ||row||^2, and Cauchy-Schwarz bounds a consistent
// dot/norm triple at 1 -- so this is defence-in-depth for direct callers of a
// public header.
void test_cosine_from_dot_never_returns_non_finite() {
  struct Case { float dot, norm_a, norm_b; const char* label; };
  const Case cases[] = {
    {1.0f, std::numeric_limits<float>::denorm_min(), 1.0f, "subnormal norm"},
    {1.0f, 1e-30f, 1e-20f, "two tiny norms"},
    {1e30f, 1e-20f, 1e-20f, "large dot, tiny norms"},
    {std::numeric_limits<float>::max(), 1e-20f, 1.0f, "FLT_MAX dot"},
    {std::numeric_limits<float>::infinity(), 1.0f, 1.0f, "infinite dot"},
    {std::numeric_limits<float>::quiet_NaN(), 1.0f, 1.0f, "NaN dot"},
    {1.0f, std::numeric_limits<float>::infinity(), 1.0f, "infinite norm"},
    {1.0f, 0.0f, 1.0f, "zero norm a"},
    {1.0f, 1.0f, 0.0f, "zero norm b"},
    {0.0f, 0.0f, 0.0f, "all zero"},
    {1.0f, std::numeric_limits<float>::denorm_min(),
     std::numeric_limits<float>::denorm_min(), "two subnormal norms"},
  };

  for (const Case& c : cases) {
    const float r = scrna::simd::cosine_from_dot(c.dot, c.norm_a, c.norm_b);
    SCRNA_CHECK(std::isfinite(r));
    // A representable similarity is a cosine and must lie in [-1, 1]; an
    // unrepresentable one is reported as no similarity.
    SCRNA_CHECK(r >= -1.0f && r <= 1.0f);
  }

  // The ordinary path is untouched: identical unit vectors still score 1.
  SCRNA_CHECK(std::fabs(scrna::simd::cosine_from_dot(1.0f, 1.0f, 1.0f) - 1.0f) < 1e-6f);
  SCRNA_CHECK(std::fabs(scrna::simd::cosine_from_dot(-1.0f, 1.0f, 1.0f) + 1.0f) < 1e-6f);
  SCRNA_CHECK(std::fabs(scrna::simd::cosine_from_dot(6.0f, 2.0f, 6.0f) - 0.5f) < 1e-6f);
}

// SM6 (audit 2026-09-11). `indptr.size() != n_rows + 1` wraps: n_rows is
// caller-supplied, SIZE_MAX + 1 is 0, and an EMPTY but non-null indptr then
// passed the extent check for SIZE_MAX rows. validate() went on to read
// row_ptr[0], row_ptr[1], ... past the end of the buffer -- ASan
// heap-buffer-overflow at block_csr.hpp:302.
//
// The allocation below is posix_memalign'd to EXACTLY four bytes on purpose:
// AlignedAllocator rounds every allocation up to 64, which is what masked
// AUDIT.md finding 1 from an earlier ASan run.
void test_row_count_overflow_is_refused() {
  void* raw = nullptr;
  if (::posix_memalign(&raw, kCacheLineBytes, sizeof(std::int32_t)) != 0) {
    std::fprintf(stderr, "posix_memalign failed; skipping\n");
    return;
  }
  std::memset(raw, 0, sizeof(std::int32_t));
  const AlignedVector<float> no_values;
  const AlignedVector<std::int32_t> no_indices;
  const std::span<const std::int32_t> empty_indptr(static_cast<const std::int32_t*>(raw), 0);

  SCRNA_CHECK_THROWS(Block_CSR<float>::adopt(
      SIZE_MAX, 4, std::span<const float>(no_values),
      std::span<const std::int32_t>(no_indices), empty_indptr, 0));
  SCRNA_CHECK_THROWS(Block_CSR<float>::from_raw(
      SIZE_MAX, 4, std::span<const float>(no_values),
      std::span<const std::int32_t>(no_indices), empty_indptr));
  std::free(raw);

  // An ordinary short indptr is still refused, and an empty matrix still builds.
  const AlignedVector<std::int32_t> one_entry{0};
  SCRNA_CHECK_THROWS(Block_CSR<float>::from_raw(
      3, 4, std::span<const float>(no_values),
      std::span<const std::int32_t>(no_indices), one_entry));
  const auto empty = Block_CSR<float>::from_raw(
      0, 4, std::span<const float>(no_values),
      std::span<const std::int32_t>(no_indices), one_entry);
  SCRNA_CHECK(empty.rows() == 0 && empty.nnz() == 0);
}

// The int32 index ceiling, checked from the extents alone. A matrix with more
// than INT32_MAX nonzeros cannot be addressed by `index_type`, and the refusal
// must happen before any element is read -- so this hands in spans that CLAIM
// 2^31 elements over a one-element buffer. If the guard moved after the first
// read, ASan would report it here rather than the test merely failing.
void test_nnz_beyond_int32_is_refused_before_any_read() {
  const AlignedVector<float> one_value{1.0f};
  const AlignedVector<std::int32_t> one_index{0};
  const AlignedVector<std::int32_t> indptr{0, 0};
  const std::size_t past_int32 = std::size_t{1} << 31;
  SCRNA_CHECK_THROWS(Block_CSR<float>::adopt(
      1, 4, std::span<const float>(one_value.data(), past_int32),
      std::span<const std::int32_t>(one_index.data(), past_int32), indptr, 0));
}

// A row of stored zeros -- an empty cell, which scRNA-seq produces and nothing
// upstream guarantees was filtered. Its norm is zero, so every cosine against
// it divides by zero. The contract is that this is a DEFINED zero, never a
// NaN: cosine_from_dot_wide returns 0 when either norm is 0, so the row scores
// 0 against everything rather than poisoning the top-k comparator.
void test_a_zero_norm_row_scores_zero_not_nan() {
  const AlignedVector<float> values{1.0f, 2.0f, 0.0f, 3.0f, 1.0f};
  const AlignedVector<std::int32_t> cols{0, 1, 2, 0, 2};
  const AlignedVector<std::int32_t> ptr{0, 2, 3, 5};  // row 1 is a stored zero
  const auto mat = Block_CSR<float>::from_raw(3, 4, values, cols, ptr);

  const auto graph = scrna::build_knn_graph_bruteforce(mat, 2);
  for (const float s : graph.neighbor_sim) SCRNA_CHECK(std::isfinite(s));
  SCRNA_CHECK(graph.neighbor_sim[2] == 0.0f);  // row 1's neighbours
  SCRNA_CHECK(graph.neighbor_sim[3] == 0.0f);
  // The other rows still score normally against each other.
  SCRNA_CHECK(graph.neighbor_sim[0] > 0.0f);

  // And the primitive says so directly, for both orders of the zero operand.
  SCRNA_CHECK(scrna::simd::cosine_from_dot_wide<float>(0.0, 0.0f, 1.0f) == 0.0f);
  SCRNA_CHECK(scrna::simd::cosine_from_dot_wide<float>(5.0, 1.0f, 0.0f) == 0.0f);
}

// Empty and zero-dimension matrices are legal shapes, not errors: a filtered
// AnnData object can legitimately have no cells left, or no genes. They must
// construct and search without indexing a zero-byte allocation -- note
// AlignedAllocator returns nullptr for n == 0, so a span over one is null with
// extent 0 and every loop below must simply not run.
void test_empty_and_zero_dimension_matrices() {
  const AlignedVector<float> no_vals;
  const AlignedVector<std::int32_t> no_idx;
  const AlignedVector<std::int32_t> just_zero{0};

  const auto zero_by_zero = Block_CSR<float>::from_raw(0, 0, no_vals, no_idx, just_zero);
  SCRNA_CHECK(zero_by_zero.rows() == 0 && zero_by_zero.cols() == 0);
  SCRNA_CHECK(zero_by_zero.nnz() == 0);
  SCRNA_CHECK(scrna::build_knn_graph_bruteforce(zero_by_zero, 5).neighbor_idx.empty());

  const auto no_rows = Block_CSR<float>::from_raw(0, 8, no_vals, no_idx, just_zero);
  SCRNA_CHECK(no_rows.cols() == 8 && no_rows.nnz() == 0);

  // Three cells with no counts at all: every row empty, so every similarity is
  // a defined zero and the neighbour slots are filled by ordering alone.
  const AlignedVector<std::int32_t> all_zero_ptr{0, 0, 0, 0};
  const auto no_counts = Block_CSR<float>::from_raw(3, 4, no_vals, no_idx, all_zero_ptr);
  const auto graph = scrna::build_knn_graph_bruteforce(no_counts, 2);
  SCRNA_CHECK(graph.neighbor_idx.size() == 6);
  for (const float s : graph.neighbor_sim) SCRNA_CHECK(s == 0.0f);

  // Zero columns is likewise a shape, not a failure.
  const AlignedVector<std::int32_t> two_rows_ptr{0, 0, 0};
  const auto no_cols = Block_CSR<float>::from_raw(2, 0, no_vals, no_idx, two_rows_ptr);
  SCRNA_CHECK(no_cols.cols() == 0);
  SCRNA_CHECK(scrna::build_knn_graph_bruteforce(no_cols, 1).neighbor_idx.size() == 2);

  // And the adopt path over empty spans, which must not trip the null-with-
  // extent guard (the extent is genuinely zero here).
  const auto adopted = Block_CSR<float>::adopt(
      0, 0, std::span<const float>(no_vals), std::span<const std::int32_t>(no_idx),
      std::span<const std::int32_t>(just_zero), 0);
  SCRNA_CHECK(adopted.rows() == 0);
}

int main() {
  std::printf("SIMD dispatch: %s\n", scrna::simd::isa_name());

  test_block_csr_alignment();
  test_block_csr_validation();
  test_block_csr_adopt();
  test_rejects_unsorted_and_duplicate_columns();
  test_span_factories_reject_mismatched_extents();
  test_span_accessors_round_trip();
  test_sparse_dot_and_cosine();
  test_dense_scratch_roundtrip();
  test_gather_dot_parity();
  test_dense_dot_parity();
  test_round_robin_schedule_covers_every_pair_exactly_once();
  test_knn_graph_bruteforce();
  test_knn_matches_reference();
  test_knn_matches_reference_multi_seed();
  test_symmetric_sweep_is_strictly_symmetric();
  test_knn_k_exceeds_rows();
  test_knn_deterministic();
  test_knn_edge_cases();
  test_knn_on_adopted_matrix();
  test_neon_gather_tail_boundaries();
  test_neon_gather_adversarial_values();
  test_neon_dense_tail_boundaries();
  test_neon_accuracy_not_worse_than_scalar();
  test_isa_dispatch_routes();
  test_knn_k_overflow_is_rejected();
  test_norms_do_not_saturate_on_large_magnitudes();
  test_cosine_from_dot_never_returns_non_finite();
  test_row_count_overflow_is_refused();
  test_nnz_beyond_int32_is_refused_before_any_read();
  test_a_zero_norm_row_scores_zero_not_nan();
  test_empty_and_zero_dimension_matrices();

#if defined(SCRNA_ENABLE_HNSW)
  std::printf("HNSW enabled; running approximate k-NN tests\n");
  test_hnsw_recall_vs_bruteforce();
  test_hnsw_dense_matches_sparse();
  test_hnsw_dense_early_exit_shapes();
  test_hnsw_memory_guard();
  test_hnsw_edge_cases();
  test_hnsw_space_equivalence();
#endif

  if (g_failures == 0) {
    std::printf("All tests passed.\n");
    return 0;
  }
  std::fprintf(stderr, "%d test(s) failed.\n", g_failures);
  return 1;
}
