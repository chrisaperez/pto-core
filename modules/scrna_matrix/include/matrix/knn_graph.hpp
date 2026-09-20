#pragma once

// `<atomic>`, `<limits>`, `<stdexcept>` and `<string>` are unconditional: the
// brute-force path also validates its arguments and runs through
// detail::parallel_for_checked, neither of which depends on hnswlib.
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "matrix/block_csr.hpp"
#include "matrix/simd_math.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

#if defined(SCRNA_ENABLE_HNSW)
#include <memory>
#include <mutex>
#include <queue>
#include <span>
#include "hnswlib/hnswlib.h"
#endif

namespace scrna {

// Dense k-NN graph in CSR-like flat form: for row i, its k neighbors are
// `neighbor_idx[i * k .. i * k + k)` sorted by descending cosine similarity,
// with matching scores in `neighbor_sim`. Rows with fewer than k available
// neighbors keep the -1 / 0.0f sentinel in the unused slots.
struct KnnGraph {
  std::size_t k = 0;
  std::vector<std::int32_t> neighbor_idx;
  std::vector<float> neighbor_sim;
};

namespace detail {

// (similarity, row). Ordered by descending similarity, ties broken by ascending
// row index so results are deterministic regardless of thread scheduling.
using Scored = std::pair<float, std::int32_t>;

struct BetterNeighbor {
  [[nodiscard]] bool operator()(const Scored& a, const Scored& b) const noexcept {
    // NaN must be handled explicitly. `a.first != b.first` is true when either
    // operand is NaN and every subsequent comparison is false, which makes NaN
    // compare "equivalent" to every other value while those values remain
    // ordered among themselves. That breaks transitivity of equivalence, and a
    // comparator that is not a strict weak ordering is undefined behaviour for
    // std::push_heap / std::sort_heap -- not merely a wrong answer.
    //
    // Block_CSR now rejects non-finite input, so a NaN score should be
    // unreachable; this keeps the ordering total regardless, since the cost of
    // being wrong here is UB rather than a bad neighbour list.
    const bool a_nan = std::isnan(a.first);
    const bool b_nan = std::isnan(b.first);
    if (a_nan || b_nan) {
      if (a_nan != b_nan) return b_nan;  // any real score beats NaN
      return a.second < b.second;        // both NaN: order by index, deterministically
    }
    if (a.first != b.first) return a.first > b.first;
    return a.second < b.second;
  }
};

// Bounded top-k accumulator.
//
// Uses BetterNeighbor as the heap comparator, which puts the *worst* current
// candidate at the root, so a new candidate is accepted in O(log k) by
// displacing that root. Memory is O(k) and never O(n).
//
// std::sort_heap with this comparator already emits descending-by-similarity
// order, so no reverse pass is needed afterwards.
class BoundedTopK {
 public:
  explicit BoundedTopK(std::size_t k) : k_(k) { heap_.reserve(k); }

  void clear() noexcept { heap_.clear(); }

  void offer(float sim, std::int32_t row) {
    const Scored cand{sim, row};
    if (heap_.size() < k_) {
      heap_.push_back(cand);
      std::push_heap(heap_.begin(), heap_.end(), BetterNeighbor{});
    } else if (k_ > 0 && BetterNeighbor{}(cand, heap_.front())) {
      std::pop_heap(heap_.begin(), heap_.end(), BetterNeighbor{});
      heap_.back() = cand;
      std::push_heap(heap_.begin(), heap_.end(), BetterNeighbor{});
    }
  }

  // Sorts in place (descending) and returns the result.
  [[nodiscard]] const std::vector<Scored>& finalize() {
    std::sort_heap(heap_.begin(), heap_.end(), BetterNeighbor{});
    return heap_;
  }

  [[nodiscard]] std::size_t size() const noexcept { return heap_.size(); }

 private:
  std::size_t k_;
  std::vector<Scored> heap_;
};

// Runs a body over [0, n) in parallel, converting any exception into a
// rethrow AFTER the parallel region.
//
// An exception that escapes an OpenMP structured block is undefined behaviour
// (OpenMP 5.2 §2.1: a throw must be caught within the same structured block),
// and in practice libomp calls std::terminate -- the process dies with no
// catchable Python error. hnswlib throws std::runtime_error from inside
// addPoint and searchKnn, including "Not enough memory" on allocation failure,
// which is precisely the failure mode expected when indexing a million cells.
//
// This lives outside the SCRNA_ENABLE_HNSW guard because it has no hnswlib
// dependency and the brute-force region needs the same protection: its
// per-thread state and its BoundedTopK both allocate, so it, too, can throw
// from inside a structured block (REVIEW_2026-08-15 finding 4).
//
// `Chunk` is the `schedule(dynamic)` chunk size. The HNSW regions were tuned at
// 256 and keep it by default; the brute-force loop wants 64, since its per-row
// cost is O(n * nnz) rather than a single index operation.
template <std::size_t Chunk = 256, typename Body>
inline void parallel_for_checked(std::size_t n, int num_threads, const char* what, Body&& body) {
  std::atomic<bool> failed{false};
  std::string message;

  #pragma omp parallel num_threads(num_threads)
  {
    #pragma omp for schedule(dynamic, Chunk)
    for (std::size_t i = 0; i < n; ++i) {
      if (failed.load(std::memory_order_relaxed)) continue;  // drain quickly, do not break out
      try {
        body(i);
      } catch (const std::exception& e) {
        bool expected = false;
        if (failed.compare_exchange_strong(expected, true)) {
          #pragma omp critical(scrna_parallel_error)
          message = e.what();
        }
      } catch (...) {
        failed.store(true, std::memory_order_relaxed);
      }
    }
  }

  if (failed.load(std::memory_order_relaxed)) {
    throw std::runtime_error(std::string(what) + " failed: " +
                              (message.empty() ? "unknown error" : message));
  }
}

// Slots in a KnnGraph of `n` rows at stride `k`, i.e. n * k -- checked.
//
// That product sizes both result vectors and then offsets every write into
// them, so a wrapped value is not a small allocation, it is a write outside
// one: n=4, k=2^62 gives 0, both assigns succeed against empty vectors, and
// row i's neighbours land at data() + i*k. Brute force had this check since
// REVIEW_2026-08-15 finding 4; the four HNSW producers did not, and a
// UBSan-confirmed store through a null pointer from query_dense() was the
// result (docs/AUDIT_2026-09-11_scrna_matrix.md, SM1). One definition now,
// called by every producer, so a fifth cannot be written without it.
[[nodiscard]] inline std::size_t checked_graph_slots(std::size_t n, std::size_t k) {
  if (k != 0 && n > std::numeric_limits<std::size_t>::max() / k) {
    throw std::invalid_argument(
        "k=" + std::to_string(k) + " with " + std::to_string(n) +
        " rows overflows the result size (n_rows * k exceeds SIZE_MAX)");
  }
  return n * k;
}

// Row-tile size for the symmetric brute-force sweep below. Large enough that
// scattering one row (the T=float path) is amortised over many gather
// partners rather than paid once per comparison -- DenseScratch's own
// comment is why that scatter exists at all -- small enough that a round has
// at least `nt` independent tile-pairs to keep every thread busy. Aiming for
// 2*nt tiles gives exactly nt pairs in a typical (non-bye) round; for n
// comparable to or smaller than that, the ceiling below collapses tile_size
// to 1, degenerating gracefully to one tile per row.
[[nodiscard]] inline std::size_t choose_tile_size(std::size_t n, int nt) noexcept {
  const std::size_t target_tiles = static_cast<std::size_t>(nt < 1 ? 2 : 2 * nt);
  const std::size_t size = (n + target_tiles - 1) / target_tiles;
  return size == 0 ? 1 : size;
}

// The classic round-robin ("circle method") tournament schedule over
// `num_tiles` slots: one round per entry of the returned vector, each round
// a list of (a, b) slot pairs with a < b. Every unordered pair of DISTINCT
// slots appears in EXACTLY one round, over num_tiles - 1 rounds (num_tiles
// even) or num_tiles rounds (odd, see below) -- and within a single round
// every slot appears in AT MOST one pair. That second property is what makes
// a round parallel-safe: two threads processing two different pairs from the
// SAME round touch disjoint slots, so the per-slot state a pair writes to
// (here, two tiles' worth of row_topk entries) needs no lock.
//
// Odd `num_tiles` gets a virtual "bye" slot (index num_tiles) so the
// even-slot algorithm below stays applicable; any pair touching it is
// dropped from its round, so that slot's real tile simply sits out one round
// in rotation -- exactly a tournament bye.
[[nodiscard]] inline std::vector<std::vector<std::pair<std::size_t, std::size_t>>>
round_robin_schedule(std::size_t num_tiles) {
  std::vector<std::vector<std::pair<std::size_t, std::size_t>>> rounds;
  if (num_tiles < 2) return rounds;

  const bool has_bye = (num_tiles % 2) != 0;
  const std::size_t m = has_bye ? num_tiles + 1 : num_tiles;  // always even
  const std::size_t bye = num_tiles;

  std::vector<std::size_t> ring(m);
  for (std::size_t i = 0; i < m; ++i) ring[i] = i;

  rounds.reserve(m - 1);
  for (std::size_t r = 0; r + 1 < m; ++r) {
    std::vector<std::pair<std::size_t, std::size_t>> pairs;
    pairs.reserve(m / 2);
    for (std::size_t i = 0; i < m / 2; ++i) {
      std::size_t a = ring[i];
      std::size_t b = ring[m - 1 - i];
      if (has_bye && (a == bye || b == bye)) continue;
      if (a > b) std::swap(a, b);
      pairs.emplace_back(a, b);
    }
    rounds.push_back(std::move(pairs));

    // Rotate every slot but ring[0] by one position: new[1] = old[m-1],
    // new[j] = old[j-1] for j in [2, m). ring[0] stays fixed throughout,
    // which is what guarantees every OTHER pair of slots meets exactly once
    // as the rotation cycles through m-1 positions.
    const std::size_t last = ring[m - 1];
    for (std::size_t j = m - 1; j >= 2; --j) ring[j] = ring[j - 1];
    ring[1] = last;
  }
  return rounds;
}

}  // namespace detail

// Brute-force exact k-NN over all rows of a Block_CSR using cosine similarity.
//
// O(n_rows^2 * avg_nnz) time. That is a correctness baseline and a practical
// tool up to the low tens of thousands of rows; it is NOT a path to 1M+ cells,
// which needs an approximate index (HNSW/Annoy-style) that this library does not
// yet provide. SIMD and memory work here improves the constant factor and the
// memory ceiling, not the asymptotic complexity.
//
// Memory: O(n_rows * k) for the result, plus O(k) heap and O(n_cols) dense
// scratch per thread. Nothing scales with n_rows per thread.
//
// For T = float the inner loop takes the scatter/gather path: row i is expanded
// once into a per-thread dense scratch buffer, then each candidate row is a
// single gather pass dispatched to AVX-512 / AVX2 / scalar at runtime. Other
// value types use the portable merge-join kernel.
template <typename T>
[[nodiscard]] inline KnnGraph build_knn_graph_bruteforce(const Block_CSR<T>& mat,
                                                           std::size_t k) {
  using Index = typename Block_CSR<T>::index_type;

  const std::size_t n = mat.rows();

  // `n * k` sizes the result and indexes every write below. Unchecked it wraps:
  // n=4, k=2**62 gives 0, both assigns succeed against an empty vector, and the
  // writes at i*k+r then land outside it. Reject before allocating anything.
  // (REVIEW_2026-08-15 finding 4.)
  const std::size_t slots = detail::checked_graph_slots(n, k);
  // Deliberately NOT clamped to n. REVIEW_2026-08-15 suggests `if (k > n) k = n`
  // as a companion to the check above, but that would change a documented,
  // tested contract: k is the row stride of the result, and a row with fewer
  // than k available neighbours keeps -1/0.0f sentinels in the unused slots
  // (see KnnGraph, and test_knn_k_exceeds_rows). The overflow check alone
  // closes the defect; clamping would silently reshape the caller's output.
  KnnGraph graph;
  graph.k = k;
  graph.neighbor_idx.assign(slots, -1);
  graph.neighbor_sim.assign(slots, 0.0f);
  if (n == 0 || k == 0) return graph;

  // Row norms, computed once. The previous arrangement recomputed both operands'
  // norms inside every one of the O(n^2) similarity calls, which is O(n^2 * nnz)
  // of redundant work; hoisting it out makes it O(n * nnz).
  std::vector<T> norms(n, T{0});
  #pragma omp parallel for schedule(static)
  for (std::size_t i = 0; i < n; ++i) {
    norms[i] = simd::sparse_norm<T>(mat.row_values(i));
  }

  int nt = 1;
#ifdef _OPENMP
  nt = omp_get_max_threads();
  if (nt < 1) nt = 1;
#endif

  // Persistent, ONE HEAP PER ROW rather than per thread. cosine similarity is
  // symmetric (sim(i,j) == sim(j,i)), so the sweep below computes each
  // unordered pair exactly once and applies the SAME score to both rows'
  // heaps -- which means a row's heap is written across many rounds,
  // potentially by a DIFFERENT thread in each one (see round_robin_schedule),
  // and can no longer be thread-local the way the old per-row-owned-by-one-
  // thread arrangement allowed. Built up front, like the old per-thread
  // state, so every throwing allocation happens on the calling thread rather
  // than inside a parallel region (REVIEW_2026-08-15 finding 4).
  std::vector<detail::BoundedTopK> row_topk;
  row_topk.reserve(n);
  for (std::size_t i = 0; i < n; ++i) row_topk.emplace_back(k);

  // Per-thread scratch for the T=float scatter/gather path only; the other
  // value types use the portable merge-join kernel directly (see
  // score_row_against_range) and need none.
  std::vector<simd::DenseScratch> thread_scratch;
  if constexpr (std::is_same_v<T, float>) {
    thread_scratch.resize(static_cast<std::size_t>(nt));
    for (auto& s : thread_scratch) s.resize(mat.cols());
  }

  const std::size_t tile_size = detail::choose_tile_size(n, nt);
  const std::size_t num_tiles = (n + tile_size - 1) / tile_size;
  auto tile_range = [tile_size, n](std::size_t t) -> std::pair<std::size_t, std::size_t> {
    const std::size_t lo = t * tile_size;
    return {lo, std::min(n, lo + tile_size)};
  };

  // Scores row `a` against every row in [b_lo, b_hi), b != a, applying each
  // result to BOTH sides. `a` is scattered (T=float path) ONCE for the whole
  // partner range rather than once per partner -- DenseScratch's own comment
  // is why that scatter exists at all, and a version of this sweep that
  // scattered per-pair measured SLOWER than the code it replaced despite
  // halving the dot-product count, because it paid an O(nnz_a) scatter
  // n_partners times over instead of once. Each pair is computed exactly
  // once and applied to both row_topk[a] and row_topk[b], so the two
  // neighbour-list entries this produces -- b in a's list, a in b's -- carry
  // the exact same float, not merely two that agree to a few ULPs: there is
  // only one computation, so there is nothing for the two to disagree about.
  auto score_row_against_range = [&](std::size_t a, std::size_t b_lo, std::size_t b_hi, int tid) {
    if constexpr (std::is_same_v<T, float>) {
      simd::DenseScratch& scratch = thread_scratch[static_cast<std::size_t>(tid)];
      const auto idx_a = mat.row_indices(a);
      scratch.scatter(idx_a, mat.row_values(a));
      for (std::size_t b = b_lo; b < b_hi; ++b) {
        if (b == a) continue;
        // Wide dot straight into the wide division -- no narrowing between
        // them. Both operands of the ratio can legitimately be ~1e60 for
        // large-magnitude rows, and each is representable only as a double;
        // narrowing the dot first turns a similarity of 1.0 into a float
        // +inf, which the guard then reports as 0. That is the brute-force
        // k-NN baseline silently deciding a row matches nothing.
        const double dot = simd::gather_dot_wide(scratch, mat.row_indices(b), mat.row_values(b));
        const float sim = simd::cosine_from_dot_wide<float>(dot, norms[a], norms[b]);
        row_topk[a].offer(sim, static_cast<std::int32_t>(b));
        row_topk[b].offer(sim, static_cast<std::int32_t>(a));
      }
      // Restore the scratch to all-zero by touching only what was scattered.
      scratch.clear(idx_a);
    } else {
      (void)tid;
      for (std::size_t b = b_lo; b < b_hi; ++b) {
        if (b == a) continue;
        const simd::WideAcc<T> dot = simd::sparse_dot_wide<T, Index>(
            mat.row_indices(a), mat.row_values(a), mat.row_indices(b), mat.row_values(b));
        const float sim = static_cast<float>(simd::cosine_from_dot_wide<T>(dot, norms[a], norms[b]));
        row_topk[a].offer(sim, static_cast<std::int32_t>(b));
        row_topk[b].offer(sim, static_cast<std::int32_t>(a));
      }
    }
  };

  // Phase 1: diagonal blocks. Tile t's own internal (i, j) pairs, i < j both
  // in t's row range -- fully independent across tiles (disjoint row
  // ranges), so this runs as one flat parallel loop with no round needed.
  // Row i is scored against (i, hi) -- everything AFTER it in the tile --
  // which is exactly the upper triangle of the tile's own n_tile x n_tile
  // block; rows before i were already covered when THEY were the "a" side.
  detail::parallel_for_checked<1>(num_tiles, nt, "brute-force k-NN (diagonal)",
                                  [&](std::size_t t) {
    int tid = 0;
#ifdef _OPENMP
    tid = omp_get_thread_num();
#endif
    const auto [lo, hi] = tile_range(t);
    for (std::size_t i = lo; i < hi; ++i) {
      score_row_against_range(i, i + 1, hi, tid);
    }
  });

  // Phase 2: off-diagonal tile pairs, one round-robin round at a time. Within
  // a round every tile appears in at most one pair (round_robin_schedule's
  // own invariant), so the pairs in one round touch disjoint row ranges and
  // run in parallel with no lock; rounds run one after another --
  // parallel_for_checked's implicit barrier at the end of each -- which is
  // what lets a row's heap move safely from one round's writer thread to a
  // different round's without a data race.
  const auto rounds = detail::round_robin_schedule(num_tiles);
  for (const auto& pairs : rounds) {
    if (pairs.empty()) continue;
    detail::parallel_for_checked<1>(pairs.size(), nt, "brute-force k-NN (round)",
                                    [&](std::size_t p) {
      int tid = 0;
#ifdef _OPENMP
      tid = omp_get_thread_num();
#endif
      const auto [ta, tb] = pairs[p];
      const auto [a_lo, a_hi] = tile_range(ta);
      const auto [b_lo, b_hi] = tile_range(tb);
      for (std::size_t a = a_lo; a < a_hi; ++a) {
        score_row_against_range(a, b_lo, b_hi, tid);
      }
    });
  }

  // Finalise: sort each row's heap into the flat result. Parallel because the
  // writes are disjoint across rows, same as the norm loop above.
  #pragma omp parallel for schedule(static)
  for (std::size_t i = 0; i < n; ++i) {
    const auto& best = row_topk[i].finalize();
    for (std::size_t r = 0; r < best.size(); ++r) {
      graph.neighbor_idx[i * k + r] = best[r].second;
      graph.neighbor_sim[i * k + r] = best[r].first;
    }
  }

  return graph;
}

// ===========================================================================
// Approximate k-NN via HNSW (hnswlib), for the scale brute force cannot reach.
//
// Enabled with -DSCRNA_ENABLE_HNSW=ON. Header-only, vendored in
// third_party/hnswlib (Apache-2.0).
//
// Two things about hnswlib shape this integration and are worth stating plainly:
//
// 1. It indexes DENSE, fixed-dimension vectors. There is no sparse
//    representation, so every row must be expanded to `dim` floats inside the
//    index. Residency is therefore ~n_rows * dim * 4 bytes plus link overhead,
//    independent of how sparse the input was. At 1.3M cells this is ~530 MB for
//    a 50-dim PCA embedding but ~10.7 GB in 2000-gene space. Use
//    hnsw_index_bytes() to check before building, and prefer a reduced
//    representation for large datasets. Construction itself densifies one row at
//    a time into per-thread scratch, so it adds only O(dim) per thread on top of
//    the index -- the full dense matrix is never materialised.
//
// 2. Its own AVX/AVX-512 kernels are gated on compile-time __AVX__/__AVX512F__
//    with no per-function target attributes, so under this project's
//    deliberately-baseline ISA flags hnswlib compiles down to SSE on x86. The
//    ScrnaCosineSpace adapter below routes distances through simd::dense_dot
//    instead, which IS runtime-dispatched to AVX2/AVX-512 via target attributes.
//    That is why it is the default: it is the only way to get AVX-512 distances
//    without putting ISA flags on the whole target and reintroducing the SIGILL
//    hazard those attributes exist to prevent.
//
// Unlike build_knn_graph_bruteforce, results are APPROXIMATE and not
// deterministic across runs when built multi-threaded: HNSW graph structure
// depends on insertion order. Raise ef_construction/ef_search for better recall.
// ===========================================================================
#if defined(SCRNA_ENABLE_HNSW)

struct HnswParams {
  std::size_t M = 16;                  // graph degree; memory scales with this
  std::size_t ef_construction = 200;   // build-time candidate breadth
  std::size_t ef_search = 0;           // query breadth; 0 -> max(2k, 64)
  std::size_t random_seed = 100;
  int num_threads = 0;                 // 0 -> OpenMP default (respects OMP_NUM_THREADS)
  // Refuse to build an index larger than this. Guards against silently
  // allocating tens of GB when handed high-dimensional input.
  std::size_t max_index_bytes = 8ull << 30;
  // Use our target-attributed SIMD kernels rather than hnswlib's own space.
  bool use_scrna_space = true;
};

// Estimated resident bytes for an index of `n` points in `dim` dimensions.
// Covers the per-element payload, level-0 links, label, per-element mutex and
// level counter -- the terms that actually dominate.
[[nodiscard]] inline std::size_t hnsw_index_bytes(std::size_t n, std::size_t dim,
                                                   std::size_t M = 16) {
  // Saturating arithmetic, not wrapping. This figure gates an allocation, so an
  // overflow here would silently wrap to a small number, wave the request past
  // the max_index_bytes guard, and hand the impossible allocation straight to
  // hnswlib. Returning SIZE_MAX makes an unrepresentable index fail the guard,
  // which is the correct direction to err in.
  constexpr std::size_t kMax = std::numeric_limits<std::size_t>::max();

  if (dim != 0 && dim > kMax / sizeof(float)) return kMax;
  const std::size_t data = dim * sizeof(float);

  if (M != 0 && M > (kMax - sizeof(unsigned int)) / (2 * sizeof(unsigned int))) return kMax;
  const std::size_t links0 = 2 * M * sizeof(unsigned int) + sizeof(unsigned int);

  const std::size_t label = sizeof(hnswlib::labeltype);
  const std::size_t per_elem_locks = sizeof(std::mutex) + sizeof(int);

  if (data > kMax - links0 - label - per_elem_locks) return kMax;
  const std::size_t per_element = data + links0 + label + per_elem_locks;

  if (per_element != 0 && n > kMax / per_element) return kMax;
  return n * per_element;
}

namespace detail {

// hnswlib space backed by this project's runtime-dispatched dense kernels.
// Vectors are unit-normalised on insertion, so inner product is cosine
// similarity and the returned distance is 1 - cosine, matching hnswlib's own
// InnerProductSpace convention.
//
// Finding 16 (2026-09-16 audit): the straightforward way to write this,
// `1.0f - simd::dense_dot(a, b, n)`, re-resolves WHICH kernel to call on
// every single invocation -- simd::isa_level()'s result is cached (a magic
// static), so this is not re-running CPUID, but it is still a function call
// into dense_dot, a static-guard check and a switch, on a callback hnswlib
// invokes through its own function-pointer indirection for every edge
// comparison during construction and every candidate during search: tens of
// millions of calls on a dataset this module's audit already sizes at 1.3M
// cells, none of which can see the ISA level change mid-run. get_dist_func()
// is called ONCE per HnswIndex (hnswlib caches the pointer it returns), so
// resolving which kernel to bind happens once here too, in the constructor,
// and every actual distance call after that goes straight to the bound
// kernel with no intervening dispatch at all.
class ScrnaCosineSpace : public hnswlib::SpaceInterface<float> {
 public:
  explicit ScrnaCosineSpace(std::size_t dim)
      : dim_(dim), data_size_(dim * sizeof(float)), dist_func_(resolve_dist_func()) {}

  size_t get_data_size() override { return data_size_; }
  hnswlib::DISTFUNC<float> get_dist_func() override { return dist_func_; }
  void* get_dist_func_param() override { return &dim_; }

 private:
  // Each wrapper below reproduces simd::dense_dot_wide's own saturation
  // fallback (run the vector kernel; a non-finite result means it saturated,
  // so recompute exactly in double) -- see that function's comment for why
  // the check is complete rather than a heuristic -- but calls its OWN kernel
  // directly, with no ISA switch. One of these is picked once, in
  // resolve_dist_func(), and that choice cannot change for the lifetime of
  // the space.
  static float distance_scalar(const void* a, const void* b, const void* dim_ptr) {
    const std::size_t n = *static_cast<const std::size_t*>(dim_ptr);
    const double wide = simd::dense_dot_scalar_wide(static_cast<const float*>(a),
                                                     static_cast<const float*>(b), n);
    return 1.0f - static_cast<float>(wide);
  }

#if defined(SCRNA_BUILD_ARM_NEON)
  static float distance_neon(const void* a, const void* b, const void* dim_ptr) {
    const std::size_t n = *static_cast<const std::size_t*>(dim_ptr);
    const auto* pa = static_cast<const float*>(a);
    const auto* pb = static_cast<const float*>(b);
    const float fast = simd::dense_dot_neon(pa, pb, n);
    if (!std::isfinite(fast)) {
      return 1.0f - static_cast<float>(simd::dense_dot_scalar_wide(pa, pb, n));
    }
    return 1.0f - fast;
  }
#endif

#if defined(SCRNA_BUILD_X86_SIMD) && defined(SCRNA_HAVE_AVX2)
  static float distance_avx2(const void* a, const void* b, const void* dim_ptr) {
    const std::size_t n = *static_cast<const std::size_t*>(dim_ptr);
    const auto* pa = static_cast<const float*>(a);
    const auto* pb = static_cast<const float*>(b);
    const float fast = simd::dense_dot_avx2(pa, pb, n);
    if (!std::isfinite(fast)) {
      return 1.0f - static_cast<float>(simd::dense_dot_scalar_wide(pa, pb, n));
    }
    return 1.0f - fast;
  }
#endif

#if defined(SCRNA_BUILD_X86_SIMD) && defined(SCRNA_HAVE_AVX512)
  static float distance_avx512(const void* a, const void* b, const void* dim_ptr) {
    const std::size_t n = *static_cast<const std::size_t*>(dim_ptr);
    const auto* pa = static_cast<const float*>(a);
    const auto* pb = static_cast<const float*>(b);
    const float fast = simd::dense_dot_avx512(pa, pb, n);
    if (!std::isfinite(fast)) {
      return 1.0f - static_cast<float>(simd::dense_dot_scalar_wide(pa, pb, n));
    }
    return 1.0f - fast;
  }
#endif

  [[nodiscard]] static hnswlib::DISTFUNC<float> resolve_dist_func() noexcept {
    switch (simd::isa_level()) {
#if defined(SCRNA_BUILD_X86_SIMD) && defined(SCRNA_HAVE_AVX512)
      case simd::IsaLevel::kAvx512: return &distance_avx512;
#endif
#if defined(SCRNA_BUILD_X86_SIMD) && defined(SCRNA_HAVE_AVX2)
      case simd::IsaLevel::kAvx2: return &distance_avx2;
#endif
#if defined(SCRNA_BUILD_ARM_NEON)
      case simd::IsaLevel::kNeon: return &distance_neon;
#endif
      default: return &distance_scalar;
    }
  }

  std::size_t dim_;
  std::size_t data_size_;
  hnswlib::DISTFUNC<float> dist_func_;
};

// Turn a squared norm into the reciprocal scale factor for unit-normalisation.
//
// Shared by all three normalisation sites so the saturation reasoning lives in
// one place. `norm_sq` must already have been accumulated in double: a float
// accumulator saturates to +inf for any element around 1e19 or above, and since
// inf > 0 the zero-guard passes, 1/sqrt(inf) is 0, and the row is written out
// as an all-zero vector that matches nothing. Finite input, silently wrong
// index, no diagnostic. Same defect as AUDIT.md S2, which fixed
// simd::sparse_norm_sq but not these loops (REVIEW_2026-08-15 finding 5).
//
// sqrt before the reciprocal: sqrt(1e60) = 1e30 is representable as a float
// even though 1e60 is not, so taking the root first keeps rows that a
// float-squared accumulator would have destroyed.
[[nodiscard]] inline float inv_norm_from_sq(double norm_sq) noexcept {
  const double norm = std::sqrt(norm_sq);
  return (norm > 0.0 && std::isfinite(norm)) ? static_cast<float>(1.0 / norm) : 0.0f;
}

// Unit-normalise `row` (length dim) into `out`, accumulating in double.
// `check_finite` rejects NaN/Inf input on the build pass; the query pass skips
// it because the build pass already vetted every row.
inline void normalize_dense_row(const float* row, std::size_t dim, float* out,
                                bool check_finite, std::size_t row_index) {
  double norm_sq = 0.0;
  for (std::size_t d = 0; d < dim; ++d) {
    if (check_finite && !std::isfinite(row[d])) {
      throw std::invalid_argument("row " + std::to_string(row_index) +
                                  " contains a non-finite value");
    }
    norm_sq += static_cast<double>(row[d]) * static_cast<double>(row[d]);
  }
  const float inv = inv_norm_from_sq(norm_sq);
  for (std::size_t d = 0; d < dim; ++d) out[d] = row[d] * inv;
}

// Expand one sparse row into `out` (length dim) and unit-normalise it.
// `out` must be zeroed on entry and is left zeroed for the caller to reuse only
// after clear_dense_row() -- we zero just the touched positions, never all dim.
inline void densify_normalized(std::span<const std::int32_t> idx,
                                std::span<const float> val, float* out) {
  const float inv = inv_norm_from_sq(simd::sparse_norm_sq<float>(val));
  for (std::size_t t = 0; t < idx.size(); ++t) {
    out[idx[t]] = val[t] * inv;
  }
}

inline void clear_dense_row(std::span<const std::int32_t> idx, float* out) {
  for (const std::int32_t c : idx) out[c] = 0.0f;
}

// Drain hnswlib's max-heap result (farthest first) into nearest-first order,
// dropping the query point itself.
inline void drain_into(std::priority_queue<std::pair<float, hnswlib::labeltype>>& pq,
                        std::size_t self_label, std::size_t k,
                        std::int32_t* out_idx, float* out_sim) {
  std::vector<std::pair<float, hnswlib::labeltype>> tmp;
  tmp.reserve(pq.size());
  while (!pq.empty()) {
    tmp.push_back(pq.top());
    pq.pop();
  }
  std::size_t w = 0;
  for (auto it = tmp.rbegin(); it != tmp.rend() && w < k; ++it) {
    if (it->second == self_label) continue;  // exclude self
    out_idx[w] = static_cast<std::int32_t>(it->second);
    // Space distance is 1 - cosine, so recover similarity directly.
    out_sim[w] = 1.0f - it->first;
    ++w;
  }
}

}  // namespace detail

// The HNSW index object itself -- construction, query, and save/load -- lives
// in matrix/hnsw_index.hpp, which includes this header for the primitives
// above. The one-shot hnsw_build_and_search{,_dense} entry points moved there
// with it, so include that header rather than this one if you call them.

#endif  // SCRNA_ENABLE_HNSW

}  // namespace scrna
