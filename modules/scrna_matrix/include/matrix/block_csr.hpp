#pragma once

#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <memory>
#include <new>
#include <span>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace scrna {

inline constexpr std::size_t kCacheLineBytes = 64;

// Minimal STL-compatible allocator that guarantees kCacheLineBytes alignment
// for every allocation, so `values_` / `col_idx_` arrays start on (and are
// padded to) an L1/L2 cache-line boundary. This is what lets a row's nonzero
// run be touched with a single cache-line-aligned prefetch/gather instead of
// straddling two lines.
template <typename T, std::size_t Align = kCacheLineBytes>
struct AlignedAllocator {
  using value_type = T;

  template <typename U>
  struct rebind {
    using other = AlignedAllocator<U, Align>;
  };

  AlignedAllocator() noexcept = default;
  template <typename U>
  AlignedAllocator(const AlignedAllocator<U, Align>&) noexcept {}

  [[nodiscard]] T* allocate(std::size_t n) {
    if (n == 0) return nullptr;
    const std::size_t bytes = ((n * sizeof(T) + Align - 1) / Align) * Align;
    void* p = ::operator new(bytes, std::align_val_t{Align});
    return static_cast<T*>(p);
  }

  void deallocate(T* p, std::size_t) noexcept {
    ::operator delete(p, std::align_val_t{Align});
  }

  template <typename U>
  bool operator==(const AlignedAllocator<U, Align>&) const noexcept { return true; }
  template <typename U>
  bool operator!=(const AlignedAllocator<U, Align>&) const noexcept { return false; }
};

template <typename T>
using AlignedVector = std::vector<T, AlignedAllocator<T, kCacheLineBytes>>;

[[nodiscard]] inline bool is_cache_aligned(const void* p) noexcept {
  return (reinterpret_cast<std::uintptr_t>(p) % kCacheLineBytes) == 0;
}

// ---------------------------------------------------------------------------
// BufferHandle: storage that is either owned (allocated through
// AlignedAllocator) or adopted (borrowed from a caller who guarantees both the
// alignment and the lifetime).
//
// The adopted case exists so a NumPy/SciPy buffer that already satisfies the
// alignment invariant can be used in place instead of being copied. The
// lifetime of the foreign memory is pinned by a type-erased owner token held
// here, so the Block_CSR cannot outlive the buffer it borrows.
//
// data()/size() are computed rather than cached so that copying or moving a
// handle can never leave a stale pointer behind.
// ---------------------------------------------------------------------------
template <typename T>
class BufferHandle {
 public:
  BufferHandle() = default;

  // Owning: takes possession of cache-aligned storage.
  explicit BufferHandle(AlignedVector<T> owned) : owned_(std::move(owned)), owns_(true) {}

  // Adopting: borrows `ptr` for `size` elements. `owner` is any type whose
  // continued existence keeps that memory valid; it is kept alive here.
  template <typename Owner>
  static BufferHandle adopt(const T* ptr, std::size_t size, Owner owner) {
    BufferHandle h;
    h.owns_ = false;
    h.ext_ptr_ = ptr;
    h.ext_size_ = size;
    h.ext_owner_ = std::make_shared<Owner>(std::move(owner));
    return h;
  }

  [[nodiscard]] const T* data() const noexcept { return owns_ ? owned_.data() : ext_ptr_; }
  [[nodiscard]] std::size_t size() const noexcept { return owns_ ? owned_.size() : ext_size_; }
  [[nodiscard]] bool owns() const noexcept { return owns_; }

 private:
  AlignedVector<T> owned_;
  const T* ext_ptr_ = nullptr;
  std::size_t ext_size_ = 0;
  // Type-erased keep-alive for the adopted case. NOTE for the pybind11 binding:
  // when Owner holds Python objects, this handle must be destroyed with the GIL
  // held, which is the case for any Block_CSR owned by a py::class_ instance.
  std::shared_ptr<void> ext_owner_;
  bool owns_ = true;
};

// Row-major, cache-aligned Compressed Sparse Row matrix.
//
// Layout mirrors scipy.sparse.csr_matrix (data / indices / indptr) so it can
// be built directly from, or exported zero-copy to, a SciPy CSR buffer. The
// difference is that `values_` and `col_idx_` sit on 64-byte boundaries, either
// because AlignedAllocator put them there or because an adopted foreign buffer
// was checked to already satisfy that.
template <typename T>
class alignas(kCacheLineBytes) Block_CSR {
 public:
  using value_type = T;
  using index_type = std::int32_t;

  Block_CSR() = default;

  Block_CSR(std::size_t n_rows, std::size_t n_cols,
            AlignedVector<T> values, AlignedVector<index_type> col_idx,
            AlignedVector<index_type> row_ptr)
      : n_rows_(n_rows), n_cols_(n_cols),
        values_(std::move(values)), col_idx_(std::move(col_idx)),
        row_ptr_(std::move(row_ptr)) {
    validate();
  }

  // Copies the CSR triple into cache-aligned storage. Use when the source
  // buffers are not already aligned, are the wrong dtype, or must not be pinned
  // for the lifetime of this matrix.
  //
  // Each buffer carries its own extent, so there is no separate `nnz` argument
  // that could disagree with the data it describes.
  static Block_CSR from_raw(std::size_t n_rows, std::size_t n_cols,
                             std::span<const T> values,
                             std::span<const index_type> indices,
                             std::span<const index_type> indptr) {
    check_extents(n_rows, values, indices, indptr, "from_raw");
    return Block_CSR(n_rows, n_cols,
                      AlignedVector<T>(values.begin(), values.end()),
                      AlignedVector<index_type>(indices.begin(), indices.end()),
                      AlignedVector<index_type>(indptr.begin(), indptr.end()));
  }

  // Adopts the CSR triple in place, with no copy.
  //
  // Requirements, all enforced here rather than trusted:
  //   * the three extents agree with each other and with `n_rows`;
  //   * every pointer is already kCacheLineBytes-aligned (throws otherwise --
  //     silently accepting a misaligned buffer would void the invariant the
  //     rest of the engine is built on);
  //   * `owner` keeps that memory alive and unmodified for as long as the
  //     returned matrix exists.
  // The usual structural validation runs too, so an adopted buffer is no less
  // checked than a copied one.
  template <typename Owner>
  static Block_CSR adopt(std::size_t n_rows, std::size_t n_cols,
                          std::span<const T> values,
                          std::span<const index_type> indices,
                          std::span<const index_type> indptr,
                          Owner owner) {
    check_extents(n_rows, values, indices, indptr, "adopt");
    if (!is_cache_aligned(values.data()) || !is_cache_aligned(indices.data()) ||
        !is_cache_aligned(indptr.data())) {
      throw std::invalid_argument(
          "Block_CSR::adopt requires all three buffers to be 64-byte aligned");
    }
    Block_CSR m;
    m.n_rows_ = n_rows;
    m.n_cols_ = n_cols;
    // One shared owner token per buffer; they typically alias the same source.
    m.values_ = BufferHandle<T>::adopt(values.data(), values.size(), owner);
    m.col_idx_ = BufferHandle<index_type>::adopt(indices.data(), indices.size(), owner);
    m.row_ptr_ = BufferHandle<index_type>::adopt(indptr.data(), indptr.size(), std::move(owner));
    m.validate();
    return m;
  }

  [[nodiscard]] std::size_t rows() const noexcept { return n_rows_; }
  [[nodiscard]] std::size_t cols() const noexcept { return n_cols_; }
  [[nodiscard]] std::size_t nnz() const noexcept { return values_.size(); }

  // True when the value/index storage was adopted in place rather than copied.
  [[nodiscard]] bool is_zero_copy() const noexcept { return !values_.owns(); }

  [[nodiscard]] std::span<const T> row_values(std::size_t row) const {
    const auto [begin, end] = row_range(row);
    return {values_.data() + begin, end - begin};
  }

  [[nodiscard]] std::span<const index_type> row_indices(std::size_t row) const {
    const auto [begin, end] = row_range(row);
    return {col_idx_.data() + begin, end - begin};
  }

  // Whole-buffer views. Prefer these to the raw pointers below: they carry
  // their extent, so they can be handed straight back to from_raw()/adopt()
  // or to any bounds-aware algorithm without a separate length travelling
  // alongside and drifting out of sync.
  [[nodiscard]] std::span<const T> values_span() const noexcept {
    return {values_.data(), values_.size()};
  }
  [[nodiscard]] std::span<const index_type> col_idx_span() const noexcept {
    return {col_idx_.data(), col_idx_.size()};
  }
  [[nodiscard]] std::span<const index_type> row_ptr_span() const noexcept {
    return {row_ptr_.data(), row_ptr_.size()};
  }

  // Raw pointers, retained for interop with C APIs and for identity/alignment
  // checks. They carry no extent; use the span accessors above for anything
  // that reads through them.
  [[nodiscard]] const T* values_data() const noexcept { return values_.data(); }
  [[nodiscard]] const index_type* col_idx_data() const noexcept { return col_idx_.data(); }
  [[nodiscard]] const index_type* row_ptr_data() const noexcept { return row_ptr_.data(); }

 private:
  // Cross-checks the three extents against each other and against n_rows,
  // before any of them is used to read memory.
  //
  // This is the check that a (pointer, length) signature could not express.
  // Previously the caller passed one `nnz` that was applied to all three
  // buffers, so a short `indices` array was read to the length of `values` --
  // an out-of-bounds read that validate() could not catch, because by then the
  // only sizes it could see were the ones the caller had asserted. A span
  // carries its own extent, so the mismatch is visible here and rejected before
  // a single element is touched.
  static void check_extents(std::size_t n_rows,
                             std::span<const T> values,
                             std::span<const index_type> indices,
                             std::span<const index_type> indptr,
                             const char* who) {
    if (indices.size() != values.size()) {
      throw std::invalid_argument(
          std::string("Block_CSR::") + who + ": indices.size()=" +
          std::to_string(indices.size()) + " must equal values.size()=" +
          std::to_string(values.size()));
    }
    // `indptr.size() - 1 != n_rows`, never `indptr.size() != n_rows + 1`.
    // n_rows is caller-supplied, and SIZE_MAX + 1 wraps to 0: an EMPTY but
    // non-null indptr then passed for SIZE_MAX rows, and validate() read
    // row_ptr[0], row_ptr[1], ... straight past a 4-byte buffer (ASan
    // heap-buffer-overflow; docs/AUDIT_2026-09-11_scrna_matrix.md, SM6).
    if (indptr.empty() || indptr.size() - 1 != n_rows) {
      throw std::invalid_argument(
          std::string("Block_CSR::") + who + ": indptr.size()=" +
          std::to_string(indptr.size()) + " must equal n_rows+1 (n_rows=" +
          std::to_string(n_rows) + ")");
    }
    // A non-empty span must actually point somewhere. std::span does not
    // enforce this, and a null data pointer with a nonzero extent would
    // otherwise reach the element loops in validate().
    if (!values.empty() && values.data() == nullptr) {
      throw std::invalid_argument(std::string("Block_CSR::") + who +
                                   ": values span is non-empty but null");
    }
    if (!indices.empty() && indices.data() == nullptr) {
      throw std::invalid_argument(std::string("Block_CSR::") + who +
                                   ": indices span is non-empty but null");
    }
    if (!indptr.empty() && indptr.data() == nullptr) {
      throw std::invalid_argument(std::string("Block_CSR::") + who +
                                   ": indptr span is non-empty but null");
    }
  }

  // Structural validation, run for both the copying and the adopting path.
  //
  // The column-index range check is deliberately a throw and not an assert: the
  // AVX2/AVX-512 gather kernels index a dense scratch buffer with these values
  // directly, with no per-lane bounds check available in hardware. An
  // out-of-range column index would become an out-of-bounds read at vector
  // speed, so it has to be rejected here, once, in every build configuration.
  void validate() const {
    // Same wrap as check_extents(): the public constructor reaches here directly.
    if (row_ptr_.size() == 0 || row_ptr_.size() - 1 != n_rows_) {
      throw std::invalid_argument("row_ptr size must equal n_rows + 1");
    }
    if (values_.size() != col_idx_.size()) {
      throw std::invalid_argument("values and col_idx must have equal length");
    }

    const std::size_t nnz = values_.size();
    if (nnz > static_cast<std::size_t>(std::numeric_limits<index_type>::max())) {
      // int32 indptr cannot address the nonzeros. Reachable at the scale this
      // library targets: ~1.3M cells x ~2k nonzeros/cell overflows int32.
      throw std::invalid_argument(
          "nnz exceeds int32 index range; a 64-bit index_type is required for this matrix");
    }

    const index_type* rp = row_ptr_.data();
    if (rp == nullptr) {
      if (n_rows_ != 0) throw std::invalid_argument("row_ptr must not be null");
      return;
    }
    if (rp[0] != 0) {
      throw std::invalid_argument("row_ptr[0] must be 0");
    }
    for (std::size_t r = 0; r < n_rows_; ++r) {
      if (rp[r + 1] < rp[r]) {
        throw std::invalid_argument("row_ptr must be non-decreasing");
      }
    }
    if (static_cast<std::size_t>(rp[n_rows_]) != nnz) {
      throw std::invalid_argument("row_ptr[n_rows] must equal nnz");
    }

    const index_type* ci = col_idx_.data();
    const T* vals = values_.data();
    const auto n_cols_i = static_cast<std::int64_t>(n_cols_);

    // One pass per row, checking range, ordering and finiteness together.
    //
    // Ordering is not cosmetic. The merge-join kernel (simd::sparse_dot) walks
    // two rows in lockstep and *requires* strictly increasing column indices;
    // the scatter/gather kernel does not. Feeding unsorted or duplicated indices
    // therefore produced two different wrong answers depending on which kernel
    // ran -- measured: two identical vectors scored 0.80 through merge-join, and
    // a duplicated index scored 1.20 through gather, a cosine similarity above
    // 1.0 that is not a possible value. Neither path complained.
    //
    // SciPy does not guarantee sorted indices (hence has_sorted_indices and
    // sort_indices()), so this is reachable from ordinary input, and rejecting
    // is the only safe response: the buffers may be adopted read-only, so
    // sorting them in place is not ours to do.
    for (std::size_t r = 0; r < n_rows_; ++r) {
      const auto begin = static_cast<std::size_t>(rp[r]);
      const auto end = static_cast<std::size_t>(rp[r + 1]);
      index_type prev = -1;
      double row_norm_sq = 0.0;
      for (std::size_t i = begin; i < end; ++i) {
        if (ci[i] < 0 || static_cast<std::int64_t>(ci[i]) >= n_cols_i) {
          throw std::invalid_argument("column index out of range [0, n_cols): index " +
                                       std::to_string(i));
        }
        if (ci[i] <= prev) {
          throw std::invalid_argument(
              std::string("column indices must be strictly increasing within each row; row ") +
              std::to_string(r) + " has " + std::to_string(ci[i]) + " after " +
              std::to_string(prev) +
              (ci[i] == prev ? " (duplicate column)" : " (unsorted)") +
              ". For a SciPy matrix call .sum_duplicates() then .sort_indices().");
        }
        prev = ci[i];

        // Reject NaN/Inf at the boundary. A NaN value propagates through the
        // norm and the dot product into the similarity score, where it reaches
        // the top-k comparator; a comparator fed NaN cannot be a strict weak
        // ordering, which makes std::sort_heap undefined behaviour rather than
        // merely inaccurate. Infinities are equally corrosive.
        if constexpr (std::is_floating_point_v<T>) {
          if (!std::isfinite(vals[i])) {
            throw std::invalid_argument(
                "values must be finite; found NaN or Inf at nonzero index " + std::to_string(i));
          }
          row_norm_sq += static_cast<double>(vals[i]) * static_cast<double>(vals[i]);
        }
      }

      // Magnitude ceiling.
      //
      // Finiteness is necessary but not sufficient. The SIMD dot kernels
      // accumulate in `float`, and by Cauchy-Schwarz every partial sum is
      // bounded by ||a||*||b||. Once that product can exceed FLT_MAX (3.4e38)
      // the accumulator saturates to +inf, and the downstream
      // `!isfinite -> 0` guard then reports a similarity of 0 for what may be
      // two IDENTICAL rows. No NaN, no exception, no warning -- just a
      // silently wrong neighbour graph.
      //
      // The kernels are no longer the last line of defence: the dispatchers in
      // simd_math.hpp recompute in double whenever a vector kernel returns a
      // non-finite result, so a caller who reaches those functions directly
      // gets the right answer regardless of this ceiling. Keeping the ceiling
      // anyway is deliberate -- it keeps the fast path fast by making the
      // fallback unreachable for anything this constructor admits, and a
      // squared norm above 1e34 is not expression data whatever the kernels do
      // with it. What changed is that it is now an optimisation and a sanity
      // check rather than the only thing between a user and a wrong graph.
      //
      // Measured on a 300x20000 matrix at 2000 nnz/row against an exact double
      // reference: at value scale 1e17 and below, max |similarity error| is
      // 7e-8 and the ordered top-15 matches exactly for every row. At 1e18 the
      // error jumps to 1.1e-1 and EVERY row gets a different neighbour set.
      // The failure has a cliff, not a slope, and the cliff is silent.
      //
      // Requiring ||row||^2 <= 1e34 gives ||a||*||b|| <= 1e34, four orders of
      // magnitude below FLT_MAX. That is unreachable for any real expression
      // matrix (raw counts, CPM, log1p and z-scores are all far below it) and
      // it is checked here, where every value is already being touched, rather
      // than by widening the O(n^2) kernel that is this library's hot path.
      if constexpr (std::is_floating_point_v<T>) {
        constexpr double kMaxRowNormSq = 1e34;
        if (row_norm_sq > kMaxRowNormSq) {
          throw std::invalid_argument(
              "row " + std::to_string(r) + " has squared norm " +
              std::to_string(row_norm_sq) +
              ", above the 1e34 ceiling, beyond which the SIMD dot kernels' "
              "float accumulator saturates and the similarity of two identical "
              "rows falls back to an exact but slow double recomputation. "
              "Rescale the matrix (these magnitudes do not occur in expression "
              "data).");
        }
      }
    }

    // Regression guard on the alignment invariant itself. Never expected to
    // fire: owned storage comes from AlignedAllocator, adopted storage is
    // checked in adopt() before we get here.
    assert((values_.size() == 0 || is_cache_aligned(values_.data())) &&
           "Block_CSR: values buffer is not cache-line aligned");
    assert((col_idx_.size() == 0 || is_cache_aligned(col_idx_.data())) &&
           "Block_CSR: col_idx buffer is not cache-line aligned");
    assert((row_ptr_.size() == 0 || is_cache_aligned(row_ptr_.data())) &&
           "Block_CSR: row_ptr buffer is not cache-line aligned");
  }

  [[nodiscard]] std::pair<std::size_t, std::size_t> row_range(std::size_t row) const {
    if (row >= n_rows_) throw std::out_of_range("row index out of range");
    return {static_cast<std::size_t>(row_ptr_.data()[row]),
            static_cast<std::size_t>(row_ptr_.data()[row + 1])};
  }

  std::size_t n_rows_ = 0;
  std::size_t n_cols_ = 0;
  BufferHandle<T> values_;
  BufferHandle<index_type> col_idx_;
  BufferHandle<index_type> row_ptr_;
};

}  // namespace scrna
