#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <cmath>
#include <cstring>
#include <limits>
#include <memory>
#include <string>
#include <type_traits>

#include "matrix/block_csr.hpp"
#include "matrix/hnsw_index.hpp"
#include "matrix/knn_graph.hpp"
#include "matrix/simd_math.hpp"

namespace py = pybind11;
using scrna::Block_CSR;
using scrna::KnnGraph;

namespace {

// Keeps the source NumPy arrays alive for as long as a Block_CSR borrows their
// memory. py::object is refcounted, so holding one here is all that is needed --
// no manual capsule bookkeeping, unlike the export direction below where the
// memory is C++-owned and NumPy is the borrower.
//
// The destructor acquires the GIL rather than merely documenting that callers
// must hold it. Releasing a py::object without the GIL is a decref race against
// the interpreter, i.e. silent heap corruption inside CPython. This object is
// reachable from a Block_CSR, which is an ordinary C++ value that any future
// binding could copy or destroy inside a py::call_guard<py::gil_scoped_release>
// region -- none does today, and that is exactly the kind of invariant that
// holds until someone adds one entry point. gil_scoped_acquire is reentrant, so
// the common already-held case costs a TLS check.
struct NumpyOwner {
  py::object data;
  py::object indices;
  py::object indptr;

  NumpyOwner(py::object d, py::object i, py::object p)
      : data(std::move(d)), indices(std::move(i)), indptr(std::move(p)) {}

  // Copying incref's three objects, so it needs the GIL for the same reason the
  // destructor does. Block_CSR::adopt takes its owner by value per buffer, so
  // this runs on the ordinary construction path, not just in exotic cases.
  NumpyOwner(const NumpyOwner& other) {
    py::gil_scoped_acquire gil;
    data = other.data;
    indices = other.indices;
    indptr = other.indptr;
  }
  NumpyOwner& operator=(const NumpyOwner&) = delete;
  // Moving only transfers pointers -- no refcount traffic, so no GIL needed.
  // The moved-from object is left null and its destructor is a no-op, which is
  // safe because pybind11's dec_ref() is Py_XDECREF.
  NumpyOwner(NumpyOwner&&) noexcept = default;
  NumpyOwner& operator=(NumpyOwner&&) = delete;

  ~NumpyOwner() {
    // Py_IsInitialized() first. After Py_Finalize() every Python object this
    // holds has already been torn down, and PyGILState_Ensure() -- which
    // gil_scoped_acquire calls -- is undefined against a dead interpreter; it
    // typically crashes or hangs. A Block_CSR outliving the interpreter is
    // reachable through a static, a leaked module reference, or an embedding
    // host that finalises while C++ objects are still alive.
    //
    // In that state the correct action is to leak the handles rather than touch
    // them: the memory is already reclaimed by interpreter shutdown, so there is
    // nothing to free and nothing that can safely be decremented.
    if (Py_IsInitialized() == 0) {
      (void)data.release();
      (void)indices.release();
      (void)indptr.release();
      return;
    }
    // Reentrant: costs a TLS check when the GIL is already held, and correctly
    // acquires it when a caller released it (e.g. a future binding that takes a
    // Block_CSR by value under py::call_guard<py::gil_scoped_release>).
    py::gil_scoped_acquire gil;
    // Released explicitly inside the guarded scope rather than left to the
    // members' own destructors: those run after this body returns, by which
    // point `gil` has already been destroyed and the GIL dropped again.
    data.release().dec_ref();
    indices.release().dec_ref();
    indptr.release().dec_ref();
  }
};

// Name stamped on the capsules to_aligned() creates. Adoption is gated on this
// (see is_adoptable), so it is a security boundary, not a debugging aid.
constexpr const char* kAlignedCapsuleName = "scrna_matrix.aligned_buffer.v1";

// Clears NPY_ARRAY_WRITEABLE on an array we own.
void freeze(py::array& arr) {
  py::detail::array_proxy(arr.ptr())->flags &=
      ~static_cast<int>(py::detail::npy_api::NPY_ARRAY_WRITEABLE_);
}

// Whether a buffer may be borrowed in place rather than copied.
//
// SECURITY. Block_CSR::adopt validates structure once, and the vector gather
// kernels then index a dense scratch buffer with the adopted column indices
// using *unchecked loads* -- hardware gathers on AVX2/AVX-512, and lane-wise
// scalar loads in the NEON path (which has no gather instruction, but is no
// less unchecked for it). That is memory-safe only while the
// validated bytes cannot change. Borrowing anything else is an arbitrary-write
// primitive; four escapes were confirmed under ASan, each from pure Python:
//
//   * `indices[0] = n_cols + 4` after adoption  -> heap-buffer-overflow WRITE
//     at an attacker-chosen offset with an attacker-chosen value;
//   * `indptr[n_rows] = 1 << 20` after adoption -> heap-buffer-overflow READ;
//   * `data.resize(k, refcheck=False)`          -> heap-use-after-free (the
//     keep-alive pins the object; NumPy reallocates the buffer beneath it, and
//     a live buffer export does not prevent it either);
//   * the first, from a second thread, while build_knn_graph holds
//     gil_scoped_release -- so revalidating on entry is not sufficient, the
//     TOCTOU window stays open for the whole kernel.
//
// This test is deliberately about PROVENANCE, not about flags. An earlier
// version required `!arr.writeable()`, and that was defeated in one line:
//
//     v = writeable_array[:]; v.flags.writeable = False
//
// gives a read-only *view* whose base is still writeable, which passed the flag
// check and restored the arbitrary write in full (re-confirmed under ASan).
// Flags describe a handle; they say nothing about the other handles onto the
// same bytes. Chasing the base chain instead would mean re-deriving NumPy's
// aliasing rules and getting them right for every future NumPy release.
//
// So: adopt only buffers this module allocated, identified by our own named
// capsule, and additionally require that the array still spans exactly that
// allocation. Those buffers are provably immutable and non-reallocatable --
// NumPy refuses `resize` ("does not own its data") and refuses
// `flags.writeable = True` ("cannot set WRITEABLE flag to True of this array"),
// both measured -- and no aliasing writeable handle to them can exist, because
// the only reference to the underlying AlignedVector is the capsule itself.
//
// The consequence is deliberate: a stock, writeable `csr.data.astype("float32")`
// takes the copy path. Call to_aligned() to get the in-place path back.
template <typename T>
[[nodiscard]] bool is_adoptable(const py::array& arr) {
  if (!arr.dtype().is(py::dtype::of<T>())) return false;
  if ((arr.flags() & py::array::c_style) == 0) return false;
  if (arr.ndim() != 1) return false;
  if (arr.writeable()) return false;  // cheap reject; not relied upon

  // The base must be OUR capsule -- not an ndarray, not a foreign capsule.
  PyObject* base = py::detail::array_proxy(arr.ptr())->base;
  if (base == nullptr) return false;  // owns its data => resizable => refuse
  if (!PyCapsule_CheckExact(base)) return false;
  if (PyCapsule_IsValid(base, kAlignedCapsuleName) != 1) return false;

  void* raw = PyCapsule_GetPointer(base, kAlignedCapsuleName);
  if (raw == nullptr) {
    PyErr_Clear();
    return false;
  }
  const auto* held = static_cast<const std::shared_ptr<scrna::AlignedVector<T>>*>(raw);
  if (held == nullptr || !*held) return false;

  // Reject any reshaped/offset/truncated view that merely inherited our capsule
  // as its base: the array must cover exactly the allocation we made.
  if (arr.data() != (*held)->data()) return false;
  if (static_cast<std::size_t>(arr.size()) != (*held)->size()) return false;

  // Follows from the allocator, re-checked because adopt() depends on it.
  return scrna::is_cache_aligned(arr.data());
}

// Views a 1-D NumPy array as a span, taking the extent from the array itself.
//
// This is the whole point of the span-based factories: the pointer and the
// length come from one object and cannot drift apart. Callers can no longer
// hand three buffers and a single length that only describes one of them.
template <typename T>
[[nodiscard]] std::span<const T> as_span(const py::array& arr) {
  return {static_cast<const T*>(arr.data()), static_cast<std::size_t>(arr.size())};
}

// Builds a Block_CSR<float> from a scipy.sparse.csr_matrix's raw
// (data, indices, indptr) buffers.
//
// Takes the in-place path when all three arrays are float32/int32,
// C-contiguous, 64-byte aligned AND non-writeable; otherwise copies into
// aligned storage. In practice that means the in-place path is reached by
// passing to_aligned() output, which is frozen and capsule-backed.
// Deliberately does NOT declare its parameters with `forcecast`: that would let
// pybind11 silently materialise a converted copy before this function can
// inspect the caller's actual buffers, which both hides a copy and destroys the
// information needed to decide whether adoption is possible.
//
// NOTE (behaviour change): a stock `csr.data.astype("float32")` is NOT adopted
// any more, even when NumPy happens to have aligned it. Such an array stays
// writeable and resizable from Python, and borrowing it was memory-unsafe --
// see the escapes enumerated on is_adoptable(). Route through to_aligned() for
// the in-place path, and read .is_zero_copy to see which path was taken.
// Defined below, beside to_aligned()'s narrowing check, which they share.
void check_value_dtype(const py::array& arr, const char* what);
void check_index_dtype(const py::array& arr, const char* what);

Block_CSR<float> from_scipy_csr(py::array data, py::array indices, py::array indptr,
                                 std::size_t n_rows, std::size_t n_cols) {
  // Dimensionality is a NumPy-level concern a span cannot express, so it is
  // still checked here. Extent agreement is not: each span below carries the
  // array's own length, and Block_CSR::check_extents cross-checks the three of
  // them before any element is read.
  if (data.ndim() != 1 || indices.ndim() != 1 || indptr.ndim() != 1) {
    throw std::invalid_argument("from_scipy_csr: data, indices and indptr must each be 1-D");
  }

  if (is_adoptable<float>(data) && is_adoptable<std::int32_t>(indices) &&
      is_adoptable<std::int32_t>(indptr)) {
    return Block_CSR<float>::adopt(
        n_rows, n_cols,
        as_span<float>(data), as_span<std::int32_t>(indices), as_span<std::int32_t>(indptr),
        NumpyOwner{data, indices, indptr});  // pins all three past the caller
  }

  // Fallback: convert (if needed) and copy into 64-byte-aligned storage.
  //
  // Checked BEFORE the forcecast, which performs the conversion these checks
  // are about and cannot report on it afterwards.
  check_value_dtype(data, "from_scipy_csr: data");
  check_index_dtype(indices, "from_scipy_csr: indices");
  check_index_dtype(indptr, "from_scipy_csr: indptr");

  auto data_c = py::array_t<float, py::array::c_style | py::array::forcecast>::ensure(data);
  auto indices_c =
      py::array_t<std::int32_t, py::array::c_style | py::array::forcecast>::ensure(indices);
  auto indptr_c =
      py::array_t<std::int32_t, py::array::c_style | py::array::forcecast>::ensure(indptr);
  if (!data_c || !indices_c || !indptr_c) {
    throw std::invalid_argument("from_scipy_csr: could not interpret inputs as float32/int32 arrays");
  }

  return Block_CSR<float>::from_raw(
      n_rows, n_cols,
      as_span<float>(data_c), as_span<std::int32_t>(indices_c),
      as_span<std::int32_t>(indptr_c));
}

// Returns a 64-byte-aligned copy of `src`, so that a subsequent
// from_scipy_csr() can adopt it in place instead of copying again. One copy
// here, then none afterwards -- worthwhile when the same matrix is used for
// repeated queries.
template <typename T>
py::array_t<T> to_aligned_impl(py::array src) {
  auto typed = py::array_t<T, py::array::c_style | py::array::forcecast>::ensure(src);
  if (!typed) throw std::invalid_argument("to_aligned: unsupported input array");

  auto buf = typed.request();
  const auto n = static_cast<std::size_t>(buf.size);

  auto owned = std::make_shared<scrna::AlignedVector<T>>(n);
  if (n > 0) std::memcpy(owned->data(), buf.ptr, n * sizeof(T));

  // A NAMED capsule. The name is what is_adoptable() authenticates, so this is
  // the token that distinguishes "memory this module allocated and nobody else
  // holds a writeable handle to" from "some array that currently looks aligned".
  //
  // The PyObject* destructor form is used (rather than the void* one) because
  // only it can carry a name. It must not throw or touch the interpreter beyond
  // the capsule itself: it runs during deallocation, possibly while an
  // exception is in flight, so the pending error state is saved and restored.
  auto* held = new std::shared_ptr<scrna::AlignedVector<T>>(owned);
  py::capsule owner(static_cast<void*>(held), kAlignedCapsuleName,
                    [](PyObject* o) {
                      PyObject *etype = nullptr, *evalue = nullptr, *etb = nullptr;
                      PyErr_Fetch(&etype, &evalue, &etb);
                      void* p = PyCapsule_GetPointer(o, kAlignedCapsuleName);
                      if (p != nullptr) {
                        delete static_cast<std::shared_ptr<scrna::AlignedVector<T>>*>(p);
                      } else {
                        PyErr_Clear();
                      }
                      PyErr_Restore(etype, evalue, etb);
                    });
  if (!owner) {
    delete held;
    throw std::runtime_error("to_aligned: could not create the owning capsule");
  }

  // `owner` is passed as the array's base, so NumPy holds a reference to the
  // capsule, the capsule holds a shared_ptr to the AlignedVector, and the
  // AlignedVector owns the bytes. The array cannot outlive its memory even if
  // every Python-side name for it is dropped.
  py::array_t<T> out({n}, {sizeof(T)}, owned->data(), owner);

  // Frozen on the way out. Combined with the capsule base this makes the buffer
  // immutable AND non-reallocatable: NumPy refuses `resize(refcheck=False)`
  // ("does not own its data") and refuses `flags.writeable = True` ("cannot set
  // WRITEABLE flag to True of this array"). Both measured.
  freeze(out);
  return out;
}

// Narrowing conversions are range-checked before they happen.
//
// `forcecast` will happily wrap int64 -> int32 (3e9 becomes -1294967296) and turn
// float64 1e300 into inf, with no exception. Either silently corrupts a matrix:
// a wrapped indptr yields a structurally invalid CSR, and an inf value poisons
// every similarity it touches. Refusing is the only safe behaviour.
// `what` is the whole message prefix (e.g. "to_aligned: int64 input"), because
// from_scipy_csr now shares this check and its messages name a different call.
template <typename Src, typename Dst>
void check_narrowing(const py::array& src, const char* what) {
  auto typed = py::array_t<Src, py::array::c_style | py::array::forcecast>::ensure(src);
  if (!typed) throw std::invalid_argument(std::string(what) + ": could not read input array");
  const auto* p = static_cast<const Src*>(typed.request().ptr);
  const auto n = static_cast<std::size_t>(typed.size());
  for (std::size_t i = 0; i < n; ++i) {
    if constexpr (std::is_floating_point_v<Src>) {
      if (!std::isfinite(p[i])) {
        throw std::invalid_argument(std::string(what) +
                                     " contains a non-finite value at index " + std::to_string(i));
      }
      if (static_cast<double>(p[i]) > static_cast<double>(std::numeric_limits<Dst>::max()) ||
          static_cast<double>(p[i]) < static_cast<double>(std::numeric_limits<Dst>::lowest())) {
        throw std::invalid_argument(std::string(what) +
                                     " value at index " + std::to_string(i) +
                                     " overflows float32");
      }
    } else {
      if (p[i] > static_cast<Src>(std::numeric_limits<Dst>::max()) ||
          p[i] < static_cast<Src>(std::numeric_limits<Dst>::lowest())) {
        throw std::invalid_argument(std::string(what) +
                                     " value at index " + std::to_string(i) +
                                     " does not fit in int32; this matrix needs 64-bit indices");
      }
    }
  }
}

// The dtype rules for from_scipy_csr's COPY path.
//
// to_aligned() has range-checked its narrowing since REVIEW_2026-08-15
// finding 3. from_scipy_csr's fallback did not: it handed the caller's arrays
// straight to py::array_t<...>::ensure with forcecast, which is that same
// conversion with no check. Measured on the unpatched module: an int64 column
// index of 2**32 + 5 -- a legal index in a matrix 2**33 columns wide -- wrapped
// to 5, which is IN range, so validate() accepted the matrix and every
// similarity computed from it silently used the wrong gene. A float64 index
// array truncated 1.9 to column 1 the same way.
// (docs/AUDIT_2026-09-11_scrna_matrix.md, SM7.)
//
// Accepted dtypes are the ones SciPy actually produces, and nothing else is
// converted on the caller's behalf: a silent conversion is what this closes.
void check_index_dtype(const py::array& arr, const char* what) {
  if (arr.dtype().is(py::dtype::of<std::int32_t>())) return;
  if (arr.dtype().is(py::dtype::of<std::int64_t>())) {
    check_narrowing<std::int64_t, std::int32_t>(arr, what);
    return;
  }
  throw std::invalid_argument(
      std::string(what) + " must be an int32 or int64 array (SciPy's own index dtypes); "
      "convert it yourself if a different dtype is what you mean");
}

// Values: float32 as-is, float64 range-checked. Integer counts are accepted
// because they cannot overflow float (int64's maximum is ~9.2e18, far below
// FLT_MAX) -- only values past 2**24 lose precision, exactly as float64 does.
void check_value_dtype(const py::array& arr, const char* what) {
  const auto& dt = arr.dtype();
  if (dt.is(py::dtype::of<float>())) {
    // Finiteness, not narrowing: float32 needs no conversion, but a NaN or Inf
    // here would be refused by Block_CSR::validate() a moment later with an
    // index into the CSR rather than into the caller's own array. Same
    // reasoning as SM11 for to_aligned; the two entry points now agree.
    check_narrowing<float, float>(arr, what);
    return;
  }
  if (dt.is(py::dtype::of<double>())) {
    check_narrowing<double, float>(arr, what);
    return;
  }
  if (dt.is(py::dtype::of<std::int32_t>()) || dt.is(py::dtype::of<std::int64_t>())) return;
  throw std::invalid_argument(
      std::string(what) + " must be a float32, float64, int32 or int64 array; "
      "convert it yourself if a different dtype is what you mean");
}

py::object to_aligned(py::array src) {
  if (src.dtype().is(py::dtype::of<std::int32_t>())) return to_aligned_impl<std::int32_t>(src);
  if (src.dtype().is(py::dtype::of<float>())) {
    // float32 was the one dtype that reached the aligned buffer unchecked: the
    // narrowing check below covers float64, so NaN and Inf were refused from a
    // float64 array and accepted from a float32 one. Block_CSR::validate()
    // still caught them at construction, so nothing wrong was ever computed --
    // but this function is documented as the checked path, and an array that
    // will be rejected later should be rejected here, where the index that
    // names the offending value is still the caller's own index.
    // docs/AUDIT_2026-09-11_scrna_matrix.md, SM11.
    check_narrowing<float, float>(src, "to_aligned: float32 input");
    return to_aligned_impl<float>(src);
  }
  if (src.dtype().is(py::dtype::of<double>())) {
    check_narrowing<double, float>(src, "to_aligned: float64 input");
    return to_aligned_impl<float>(src);
  }
  if (src.dtype().is(py::dtype::of<std::int64_t>())) {
    check_narrowing<std::int64_t, std::int32_t>(src, "to_aligned: int64 input");
    return to_aligned_impl<std::int32_t>(src);
  }
  throw std::invalid_argument("to_aligned: expected a float32/float64 or int32/int64 array");
}

// Zero-copy export: wrap KnnGraph's owned std::vector buffers directly as
// NumPy arrays via a py::capsule that keeps the KnnGraph alive for the
// array's lifetime, instead of copying into a fresh buffer. (The mirror image of
// from_scipy_csr's adoption path, which borrows in the other direction.)
py::tuple knn_graph_to_numpy(std::shared_ptr<KnnGraph> graph) {
  const std::size_t n = graph->neighbor_idx.size();
  KnnGraph* raw = graph.get();

  py::capsule idx_owner(new std::shared_ptr<KnnGraph>(graph),
                         [](void* p) { delete static_cast<std::shared_ptr<KnnGraph>*>(p); });
  py::capsule sim_owner(new std::shared_ptr<KnnGraph>(graph),
                         [](void* p) { delete static_cast<std::shared_ptr<KnnGraph>*>(p); });

  py::array_t<std::int32_t> idx_arr({n}, {sizeof(std::int32_t)}, raw->neighbor_idx.data(), idx_owner);
  py::array_t<float> sim_arr({n}, {sizeof(float)}, raw->neighbor_sim.data(), sim_owner);

  return py::make_tuple(idx_arr, sim_arr, graph->k);
}

}  // namespace

PYBIND11_MODULE(scrna_matrix_py, m) {
  m.doc() = "Cache-aligned Block_CSR sparse matrix engine for scRNA-seq (zero-copy pybind11 bindings)";

  py::class_<Block_CSR<float>>(m, "BlockCSR")
      .def(py::init<>())
      .def_static("from_scipy_csr", &from_scipy_csr,
                  py::arg("data"), py::arg("indices"), py::arg("indptr"),
                  py::arg("n_rows"), py::arg("n_cols"),
                  "Build from scipy CSR buffers. Adopts them in place when they are "
                  "float32/int32, C-contiguous, 64-byte aligned and non-writeable; "
                  "copies otherwise. In practice that means pass to_aligned() output "
                  "to get the in-place path -- a writeable array is always copied, "
                  "because borrowing memory Python can still mutate or reallocate is "
                  "not memory-safe. Check .is_zero_copy to see which happened.\n\n"
                  "The copy path narrows like to_aligned(): float64 values and int64 "
                  "indices are range-checked, and any other dtype is refused rather "
                  "than converted silently.")
      .def_property_readonly("shape", [](const Block_CSR<float>& m) {
        return py::make_tuple(m.rows(), m.cols());
      })
      .def_property_readonly("nnz", &Block_CSR<float>::nnz)
      .def_property_readonly("is_zero_copy", &Block_CSR<float>::is_zero_copy,
                             "True if the source buffers were adopted rather than copied.");

  py::class_<KnnGraph, std::shared_ptr<KnnGraph>>(m, "KnnGraph")
      .def_readonly("k", &KnnGraph::k)
      .def("to_numpy", [](std::shared_ptr<KnnGraph> self) { return knn_graph_to_numpy(self); });

  m.def(
      "build_knn_graph",
      [](const Block_CSR<float>& mat, std::size_t k, const std::string& method,
         std::size_t M, std::size_t ef_construction, std::size_t ef_search,
         int num_threads, std::size_t max_index_bytes, bool use_scrna_space) {
        if (method == "brute") {
          return std::make_shared<KnnGraph>(scrna::build_knn_graph_bruteforce(mat, k));
        }
        if (method == "hnsw") {
#if defined(SCRNA_ENABLE_HNSW)
          scrna::HnswParams p;
          p.M = M;
          p.ef_construction = ef_construction;
          p.ef_search = ef_search;
          p.num_threads = num_threads;
          p.max_index_bytes = max_index_bytes;
          p.use_scrna_space = use_scrna_space;
          return std::make_shared<KnnGraph>(scrna::hnsw_build_and_search(mat, k, p));
#else
          throw std::invalid_argument(
              "method='hnsw' requires building with -DSCRNA_ENABLE_HNSW=ON");
#endif
        }
        throw std::invalid_argument("method must be 'brute' or 'hnsw', got '" + method + "'");
      },
      py::arg("matrix"), py::arg("k"), py::arg("method") = "brute",
      py::arg("M") = 16, py::arg("ef_construction") = 200, py::arg("ef_search") = 0,
      py::arg("num_threads") = 0, py::arg("max_index_bytes") = (8ull << 30),
      py::arg("use_scrna_space") = true,
      // The C++ side is OpenMP-parallel numeric work that touches no Python
      // objects, so holding the GIL across it only blocks other Python threads.
      py::call_guard<py::gil_scoped_release>(),
      "Cosine k-NN graph over a BlockCSR.\n\n"
      "method='brute': exact, O(n^2), deterministic.\n"
      "method='hnsw':  approximate, scales to 1M+ cells, NOT deterministic when\n"
      "                built multi-threaded. Note hnswlib stores DENSE vectors, so\n"
      "                memory is n_cells * n_genes * 4 bytes regardless of input\n"
      "                sparsity -- for wide gene-space matrices prefer\n"
      "                build_knn_graph_dense() on a PCA embedding.");

  m.def(
      "build_knn_graph_dense",
      [](py::array_t<float, py::array::c_style | py::array::forcecast> data, std::size_t k,
         const std::string& method, std::size_t M, std::size_t ef_construction,
         std::size_t ef_search, int num_threads, std::size_t max_index_bytes,
         bool use_scrna_space) {
        auto buf = data.request();
        if (buf.ndim != 2) throw std::invalid_argument("expected a 2-D (n_obs, n_dims) array");
        const auto n = static_cast<std::size_t>(buf.shape[0]);
        const auto dim = static_cast<std::size_t>(buf.shape[1]);
        const float* ptr = static_cast<const float*>(buf.ptr);

        if (method != "hnsw") {
          throw std::invalid_argument("build_knn_graph_dense currently supports method='hnsw'");
        }
#if defined(SCRNA_ENABLE_HNSW)
        scrna::HnswParams p;
        p.M = M;
        p.ef_construction = ef_construction;
        p.ef_search = ef_search;
        p.num_threads = num_threads;
        p.max_index_bytes = max_index_bytes;
        p.use_scrna_space = use_scrna_space;

        std::shared_ptr<KnnGraph> out;
        {
          // `ptr` is plain memory; no Python API is touched inside.
          py::gil_scoped_release release;
          out = std::make_shared<KnnGraph>(
              scrna::hnsw_build_and_search_dense(ptr, n, dim, k, p));
        }
        return out;
#else
        throw std::invalid_argument(
            "build_knn_graph_dense requires building with -DSCRNA_ENABLE_HNSW=ON");
#endif
      },
      py::arg("data"), py::arg("k"), py::arg("method") = "hnsw",
      py::arg("M") = 16, py::arg("ef_construction") = 200, py::arg("ef_search") = 0,
      py::arg("num_threads") = 0, py::arg("max_index_bytes") = (8ull << 30),
      py::arg("use_scrna_space") = true,
      "Approximate cosine k-NN over a dense (n_obs, n_dims) float32 array, e.g. a\n"
      "PCA embedding. This is the path that scales to 1M+ cells.");

  m.def(
      "hnsw_index_bytes",
      [](std::size_t n, std::size_t dim, std::size_t M) {
#if defined(SCRNA_ENABLE_HNSW)
        return scrna::hnsw_index_bytes(n, dim, M);
#else
        (void)n; (void)dim; (void)M;
        throw std::invalid_argument("built without -DSCRNA_ENABLE_HNSW=ON");
#endif
      },
      py::arg("n_obs"), py::arg("n_dims"), py::arg("M") = 16,
      "Estimated resident bytes for an HNSW index. Check this before building on\n"
      "wide matrices: hnswlib is dense, so sparsity does not reduce it.");

#if defined(SCRNA_ENABLE_HNSW)
  // ------------------------------------------------------------------------
  // HnswIndex — a reference atlas you build once and query many times.
  //
  // build_knn_graph_dense() above is still the right call for "give me the
  // graph of this one dataset". This class exists for the workflow it cannot
  // express: build an index over a reference (minutes), persist it, and later
  // map new cells against it (milliseconds) in a different process.
  // ------------------------------------------------------------------------
  py::class_<scrna::HnswIndex>(m, "HnswIndex", R"doc(
A persistable approximate nearest-neighbour index over dense float32 vectors.

Vectors are unit-normalised on insertion, so the reported similarity is cosine.

    idx = scrna_matrix.HnswIndex(dim=50, capacity=1_300_000)
    idx.add_dense(atlas_pca)          # float32 (n, 50)
    idx.set_labels(barcodes)          # optional
    idx.save_index("pbmc1M.hnsw")     # writes pbmc1M.hnsw + pbmc1M.hnsw.meta

    idx = scrna_matrix.HnswIndex.load_index("pbmc1M.hnsw")
    graph = idx.query_dense(new_cells_pca, k=15)

SECURITY: load_index() trusts the file it is given. It validates the header,
geometry and checksum, but the graph body is parsed by vendored hnswlib code.
Load index files only from sources you would trust to run code, and never wire
this into a network-facing service. See docs/AUDIT.md, 2026-08-16.
)doc")
      .def(py::init([](std::size_t dim, std::size_t capacity, std::size_t M,
                       std::size_t ef_construction, std::size_t ef_search,
                       int num_threads, std::size_t max_index_bytes,
                       bool use_scrna_space) {
             scrna::HnswParams p;
             p.M = M;
             p.ef_construction = ef_construction;
             p.ef_search = ef_search;
             p.num_threads = num_threads;
             p.max_index_bytes = max_index_bytes;
             p.use_scrna_space = use_scrna_space;
             return std::make_unique<scrna::HnswIndex>(dim, capacity, p);
           }),
           py::arg("dim"), py::arg("capacity"), py::arg("M") = 16,
           py::arg("ef_construction") = 200, py::arg("ef_search") = 0,
           py::arg("num_threads") = 0, py::arg("max_index_bytes") = (8ull << 30),
           py::arg("use_scrna_space") = true)

      .def("add_dense",
           [](scrna::HnswIndex& self,
              py::array_t<float, py::array::c_style | py::array::forcecast> data,
              int num_threads) {
             auto buf = data.request();
             if (buf.ndim != 2) {
               throw std::invalid_argument("expected a 2-D (n_obs, n_dims) array");
             }
             if (static_cast<std::size_t>(buf.shape[1]) != self.dim()) {
               throw std::invalid_argument(
                   "array has " + std::to_string(buf.shape[1]) +
                   " columns but the index has dim " + std::to_string(self.dim()));
             }
             const auto n = static_cast<std::size_t>(buf.shape[0]);
             const float* ptr = static_cast<const float*>(buf.ptr);
             // `ptr` is plain memory and the body touches no Python objects.
             // Construction is the multi-threaded, minutes-long part, so the
             // GIL must not be held across it.
             py::gil_scoped_release release;
             self.add_dense(ptr, n, num_threads);
           },
           py::arg("data"), py::arg("num_threads") = 0,
           "Append rows of a dense (n_obs, n_dims) float32 array to the index.")

      .def("add_csr",
           [](scrna::HnswIndex& self, const Block_CSR<float>& mat, int num_threads) {
             self.add_csr(mat, num_threads);
           },
           py::arg("matrix"), py::arg("num_threads") = 0,
           py::call_guard<py::gil_scoped_release>(),
           "Append every row of a BlockCSR, densifying on the fly. The index "
           "stores dense vectors, so sparsity does not reduce its memory.")

      .def("query_dense",
           [](const scrna::HnswIndex& self,
              py::array_t<float, py::array::c_style | py::array::forcecast> data,
              std::size_t k, int num_threads) {
             auto buf = data.request();
             if (buf.ndim != 2) {
               throw std::invalid_argument("expected a 2-D (n_obs, n_dims) array");
             }
             if (static_cast<std::size_t>(buf.shape[1]) != self.dim()) {
               throw std::invalid_argument(
                   "array has " + std::to_string(buf.shape[1]) +
                   " columns but the index has dim " + std::to_string(self.dim()));
             }
             const auto n = static_cast<std::size_t>(buf.shape[0]);
             const float* ptr = static_cast<const float*>(buf.ptr);
             std::shared_ptr<KnnGraph> out;
             {
               py::gil_scoped_release release;
               // No self-exclusion: query rows are foreign to the index, which
               // is the entire point of a saved reference. Excluding label i
               // from row i's results here would silently drop a legitimate
               // neighbour whose label happens to match the query's position.
               out = std::make_shared<KnnGraph>(self.query_dense(ptr, n, k, num_threads));
             }
             return out;
           },
           py::arg("data"), py::arg("k"), py::arg("num_threads") = 0,
           "k nearest index members for each row of a dense query array.\n\n"
           "Neighbour ids are index labels (insertion order). Query rows are "
           "treated as foreign to the index, so a query identical to an indexed "
           "vector will return that vector as its own nearest neighbour.")

      .def("query_csr",
           [](const scrna::HnswIndex& self, const Block_CSR<float>& mat, std::size_t k,
              int num_threads) {
             return std::make_shared<KnnGraph>(self.query_csr(mat, k, num_threads));
           },
           py::arg("matrix"), py::arg("k"), py::arg("num_threads") = 0,
           py::call_guard<py::gil_scoped_release>(),
           "k nearest index members for each row of a BlockCSR query matrix.")

      .def("save_index", &scrna::HnswIndex::save_index, py::arg("path"),
           py::call_guard<py::gil_scoped_release>(),
           "Write the index as a pair of files: <path> and <path>.meta. Both are "
           "required to load it again.")

      .def_static("load_index",
                  [](const std::string& path, std::size_t ef_search, int num_threads,
                     std::size_t max_index_bytes, bool use_scrna_space,
                     bool verify_checksum) {
                    scrna::HnswParams p;
                    p.ef_search = ef_search;
                    p.num_threads = num_threads;
                    p.max_index_bytes = max_index_bytes;
                    p.use_scrna_space = use_scrna_space;
                    py::gil_scoped_release release;
                    return std::make_unique<scrna::HnswIndex>(
                        scrna::HnswIndex::load_index(path, p, verify_checksum));
                  },
                  py::arg("path"), py::arg("ef_search") = 0, py::arg("num_threads") = 0,
                  py::arg("max_index_bytes") = (8ull << 30),
                  py::arg("use_scrna_space") = true, py::arg("verify_checksum") = true,
                  "Load an index written by save_index().\n\n"
                  "verify_checksum re-reads the payload to check its CRC-32. Leave "
                  "it on unless you just wrote the file yourself: it is what catches "
                  "a payload truncated or corrupted AFTER it was written. A save that "
                  "runs out of disk is refused by save_index() itself, which compares "
                  "the bytes on disk against the size hnswlib computes for the graph.")

      .def("set_labels", &scrna::HnswIndex::set_labels, py::arg("names"),
           "Attach per-element names (e.g. cell barcodes), in insertion order. "
           "Serialized with the index.")
      .def_property_readonly("labels", &scrna::HnswIndex::labels)
      .def_property_readonly("dim", &scrna::HnswIndex::dim)
      .def_property_readonly("capacity", &scrna::HnswIndex::capacity)
      .def("set_ef", &scrna::HnswIndex::set_ef, py::arg("ef"),
           "Query-time candidate breadth. Higher is more accurate and slower; "
           "0 restores the default of max(2k, 64).")
      .def("__len__", &scrna::HnswIndex::size);
#endif  // SCRNA_ENABLE_HNSW

  m.def("has_hnsw", []() {
#if defined(SCRNA_ENABLE_HNSW)
    return true;
#else
    return false;
#endif
  }, "Whether approximate k-NN (hnswlib) was compiled in.");

  m.def("to_aligned", &to_aligned, py::arg("array"),
        "Return a 64-byte-aligned, READ-ONLY copy of a 1-D array so from_scipy_csr "
        "can adopt it without copying. float64 is narrowed to float32 and int64 to "
        "int32.\n\n"
        "The result is deliberately immutable: it is capsule-backed, so NumPy "
        "refuses both in-place writes and resize(), which is what makes it safe for "
        "the engine to borrow. Mutate your source array before calling this, not the "
        "array it returns.");

  m.def("simd_isa", []() { return std::string(scrna::simd::isa_name()); },
        "Name of the SIMD kernel selected at runtime: 'avx512', 'avx2' (x86-64), "
        "'neon' (arm64) or 'scalar'.\n\n"
        "Override with SCRNA_FORCE_ISA=scalar|neon|avx2|avx512 to pin the kernel, "
        "e.g. to compare a vector path against the scalar reference on one host. "
        "Requesting a path that was not compiled in is ignored.");
}
