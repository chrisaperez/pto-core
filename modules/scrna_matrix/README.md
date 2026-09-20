# scrna_matrix

A cache-aligned sparse matrix engine for high-dimensional scRNA-seq.

Header-only C++20, built around a 64-byte cache-line-aligned block-compressed
sparse row format (`Block_CSR<T>`), with runtime-dispatched SIMD sparse
dot-product and cosine-similarity kernels and zero-copy `pybind11`
interoperability with NumPy and SciPy `csr_matrix` buffers.

## Motivation

scRNA-seq datasets regularly exceed 1 million cells. Standard analysis
frameworks store expression data as CSC or CSR matrices, but k-NN graph
construction, cell-cell cosine distance and library-size normalization access
that memory non-contiguously. The result is cache thrashing, and at scale an
out-of-memory failure on hardware that should have been sufficient.

`Block_CSR` aligns its value and index arrays to cache-line boundaries and
exposes SIMD kernels tuned for sparse-vector operations, so these workloads stay
in cache. Equally important, it avoids duplicating data when crossing the
C++/Python boundary: a buffer that already satisfies the alignment invariant is
adopted in place rather than copied.

## Directory structure

```
scrna_matrix/
├── CMakeLists.txt
├── docs/
│   ├── AUDIT.md            # static-analysis and security review log
│   ├── BENCHMARKS.md
│   └── MANUSCRIPT.tex
├── include/matrix/
│   ├── block_csr.hpp       # alignas(64) Block_CSR<T>, AlignedAllocator, adopt path
│   ├── knn_graph.hpp       # exact brute-force k-NN and the HNSW pathway
│   └── simd_math.hpp       # scalar/NEON/AVX2/AVX-512 kernels, runtime dispatch,
│                           #   WideAcc and checked_narrow
├── third_party/
│   └── hnswlib/            # vendored hnswlib v0.9.0 (Apache-2.0), locally patched
├── src/
│   └── python_bindings.cpp # pybind11 zero-copy bindings
├── tests/
│   ├── test_matrix_ops.cpp
│   ├── test_wide_accumulation.cpp  # source scan; fails the build on a bare float acc
│   ├── test_hnsw_persist.cpp
│   ├── test_python_bindings.py
│   └── benchmark_1m_cells.py
└── examples/
    ├── scrna_anndata.py           # scanpy-compatible neighbors()
    └── scanpy_integration.ipynb
```

The vendored hnswlib is **patched, not pristine**. The patches are recorded in
`third_party/hnswlib/VERSION.txt`; the most important is patch 4, which threads
the `ef` search parameter through `searchKnnWithEf` so a query never writes
shared index state.

## Build

Requires a C++20 compiler and CMake 3.20 or newer, and optionally OpenMP, Eigen3
and pybind11.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```

This produces `build/test_matrix_ops` and, when pybind11 is found, a
`scrna_matrix_py*.so` extension module importable from Python by adding `build/`
to `PYTHONPATH`. `ctest` runs both the C++ tests and, when pytest is available,
the binding tests.

Build a `Debug` configuration too when changing the engine. The alignment and
span-length invariants are `assert`-based and compile out under `NDEBUG`.

AVX-512 kernels are opt-in at configure time with `-DSCRNA_ENABLE_AVX512=ON`,
and AVX2 is compiled in automatically wherever the toolchain supports it. Both
are selected at runtime through CPUID, with per-function `target` attributes
rather than a whole-binary ISA flag, so enabling them does not make the binary
unrunnable on older CPUs.

**OpenMP is mandatory for this module and configuration fails without it.** This
is the one module in `pto-core` where a missing dependency is a hard error, and
it is deliberate. It used to be optional with a warning fallback, and the result
was that `_OPENMP` went undefined, every `#pragma omp` compiled to serial code,
and the suite reported 100% passed while never executing a parallel iteration. A
silently single-threaded engine that CI signs off on is worse than a build that
stops.

### OpenMP on macOS

Apple's clang does not ship an OpenMP runtime, so `find_package(OpenMP)` fails
by default. With `libomp` installed from Homebrew or conda, point CMake at it:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DOpenMP_CXX_FLAGS="-Xpreprocessor -fopenmp -I$PREFIX/include" \
  -DOpenMP_CXX_LIB_NAMES=omp \
  -DOpenMP_omp_LIBRARY=$PREFIX/lib/libomp.dylib \
  -DCMAKE_BUILD_RPATH=$PREFIX/lib
```

Note `-Xpreprocessor` rather than `-Xclang`. The latter passes CMake's own
detection test and then fails to link, which is exactly the kind of
false-positive detection the hard-fail probe exists to catch. The configure
output reports which case you are in.

## Python usage

```python
import scrna_matrix_py as scrna
import scipy.sparse as sp

csr = sp.load_npz("cells.npz")

# to_aligned guarantees the 64-byte alignment that lets BlockCSR adopt these
# buffers in place. Without it, adoption still happens whenever NumPy's own
# allocation is already aligned (usual at realistic sizes) but is not guaranteed.
mat = scrna.BlockCSR.from_scipy_csr(
    scrna.to_aligned(csr.data.astype("float32")),
    scrna.to_aligned(csr.indices.astype("int32")),
    scrna.to_aligned(csr.indptr.astype("int32")),
    csr.shape[0], csr.shape[1],
)
print(mat.is_zero_copy, scrna.simd_isa())   # True scalar|avx2|avx512

graph = scrna.build_knn_graph(mat, k=15)          # exact, O(n^2)
neighbor_idx, neighbor_sim, k = graph.to_numpy()
```

**Do not `.astype()` your indices before handing them over.** The narrowing
conversion at the NumPy boundary is range-checked in both entry points, and an
int64 column index that wraps during your own cast lands *in range*, so nothing
downstream can tell it was ever wrong.

A buffer is adopted only when it is provably immutable, and that is not
fastidiousness. The AVX2 and AVX-512 gather kernels index a dense scratch buffer
with the adopted column indices using unchecked hardware gathers. This means the
bounds validation performed once at construction is the only thing standing
between the inner loop and an arbitrary memory write, which is exactly why a
buffer the caller can still mutate afterward may not be adopted at all.

### Scaling past roughly 50k cells

`method="hnsw"` swaps the exact search for an HNSW index. The caveat that
matters is that hnswlib indexes **dense** vectors and has no sparse
representation, so it belongs on a reduced representation rather than raw gene
space:

```python
# Recommended: approximate search on a PCA embedding.
graph = scrna.build_knn_graph_dense(X_pca.astype("float32"), k=15, method="hnsw")

# Check the cost before building; sparsity does not reduce it.
scrna.hnsw_index_bytes(1_300_000, 50)    # 0.49 GB
scrna.hnsw_index_bytes(1_300_000, 2000)  # 9.94 GB
```

Measured: **1.3M cells by 50 dimensions in 37.1 s at 1.01 GB peak RSS**,
including index construction and all 1.3M queries, against roughly 1.6 hours
extrapolated for the exact path at about 3.5 ns per n-squared.

In 1,000-dimensional gene space the ordering inverts and HNSW is **slower** than
exact, by 0.05x. Each HNSW distance touches 1,000 dense dimensions where the
exact path touches roughly 60 non-zeros, so a smaller number of comparisons does
not pay for their individual cost. On a 50-dimensional PCA embedding with real
cluster structure, HNSW wins by 11.2x at 50,000 cells at a recall of 0.997. In
other words, the representation and not the algorithm decides which backend
wins, and this module states that rather than leaving you to rediscover it by
benchmarking both.

HNSW results are approximate and not reproducible across multi-threaded runs.
Use `method="brute"` when you need ground truth or determinism.

### Parameter bounds that are not arbitrary

| Parameter | Bound | Why |
|:--|:--|:--|
| `HnswParams::M` | 2 to 10,000 | The range hnswlib honours. Below it the library is undefined behaviour, since `1/log(1)` is infinity cast to `int`. Above it hnswlib clamps silently and the index saves but cannot load. |
| `capacity` | at most 2,147,483,647 | `KnnGraph` neighbour ids are int32. Refused in both the constructor and the loader. |
| `k` | any | It is a row stride and is never clamped. Short rows keep the `-1` sentinel, even though the candidate count asked of hnswlib is clamped. |
| `n_cols` | unbounded | A stated residual. A huge declared width fails as `std::bad_alloc` rather than a named error. |

### Scanpy and AnnData

`examples/scrna_anndata.py` is a drop-in replacement for `sc.pp.neighbors`:

```python
import scanpy as sc
from scrna_anndata import neighbors

adata = sc.datasets.pbmc3k()
# ... standard preprocessing ...
neighbors(adata, n_neighbors=15)   # fills obsp/uns the way scanpy expects
sc.tl.leiden(adata)                # works unchanged
```

See [`examples/scanpy_integration.ipynb`](examples/scanpy_integration.ipynb) for
the full worked example with live outputs.

### Selecting the SIMD kernel

Dispatch is automatic. `SCRNA_FORCE_ISA=scalar|avx2|avx512` overrides it for
testing, which is how every compiled path gets exercised on one machine.
Requesting a kernel the CPU cannot execute will fault, which is why it is
opt-in.

## What this module does not do

Read this before planning work here, because the boundary is not obvious from
the name.

There is **no Matrix Market reader, no barcode or feature TSV reader, and no
CLI**. That parsing belongs to scanpy and anndata upstream of this module. There
is also no per-cell or per-gene accumulator to overflow, which rules out a whole
class of question people arrive with.

## Status

Implemented and tested:

* `Block_CSR<T>` cache-aligned layout, with structural validation (column
  indices bounds-checked at construction, which is what makes the gather kernels
  memory-safe) and an adopt-in-place path for buffers that already satisfy the
  alignment invariant.
* Scalar, NEON, AVX2 and AVX-512 kernels for sparse cosine similarity, selected
  by runtime CPU dispatch.
* Exact brute-force k-NN with O(k) per-thread memory through a bounded top-k
  heap, and row norms hoisted out of the O(n^2) loop.
* Approximate k-NN via vendored hnswlib, with a memory guard that refuses
  oversized indices rather than thrashing, densification done one row at a time
  so the full dense matrix is never materialised, and distances routed through
  the runtime-dispatched kernels.
* Zero-copy Python interop in both directions, with `.is_zero_copy` reporting
  which path was taken, and the GIL released during graph construction.
* End-to-end AnnData integration writing `obsp['distances']`,
  `obsp['connectivities']` and `uns['neighbors']` in scanpy's own layout, so
  `sc.tl.leiden`, `umap` and `paga` consume it directly, matching
  `sc.pp.neighbors` neighbour for neighbour on pbmc3k.

Not done, and not to be quoted until it is:

* **A like-for-like 1.3M-cell comparison against `scanpy.pp.neighbors`.** The
  scalability goal is met. The 20x-speedup and 60%-RSS claims are comparisons
  that have not been run. Note that `sc.pp.neighbors` on PCA is itself
  approximate and near-linear, so it is a peer of `method="hnsw"` rather than a
  baseline the module should categorically dominate.
* **HNSW recall above 50k cells.** It is 0.997 at 50k; beyond that there is no
  feasible exact baseline to measure against. There has been no HNSW parameter
  sweep at all, and every number uses the defaults.
* **Throughput numbers for the AVX2 and AVX-512 kernels.** They are compile- and
  correctness-verified but have not run on hardware that natively supports them,
  so every benchmark on record is from the scalar path.
* **A structure-aware fuzzer over `from_scipy_csr`.** That is the entry point at
  which one of the audit's eleven findings, an unchecked narrowing conversion,
  was found by manual probing rather than by fuzzing.

## Correctness rules specific to this module

Two of these are enforced by tests that read the source rather than run it,
because a rule a human has to remember is a rule that gets broken in the next
sibling function.

**Accumulate floating-point reductions in `WideAcc<T>` and narrow once through
`checked_narrow()`.** The intuition about why is wrong in an instructive way. A
bare `float` accumulator does not gradually lose precision on large input. It
saturates to infinity, the downstream non-finite guard maps that to a similarity
of zero, and two identical rows become each other's worst possible match. The
input was finite, no exception was raised, no diagnostic was printed, and the
neighbour graph is wrong. This is the module's most repeated defect: four
instances found by four different techniques across three review passes, none of
which was looking for it.

Two corollaries. Stay wide until the last operation, because
`cosine_from_dot_wide(dot, na, nb)` divides 1e60 by 1e60 and gets 1.0 where
narrowing the dot first gets infinity and reports 0. And do not add a
`float acc = 0` at all: `test_wide_accumulation.cpp` scans the module's sources
and fails the build. A deliberate narrow accumulator marks itself
`// SCRNA_FLOAT_ACC_OK: <why>`, which is a grep target rather than a rubber
stamp. The four vector-kernel scalar tails are the only ones today, and they are
safe because the dispatchers recompute in double whenever a kernel returns a
non-finite result.

**Size every `KnnGraph` through `detail::checked_graph_slots(n, k)`.** That
product also offsets every write, so an unchecked `n * k` is a wild store rather
than a small allocation. Brute force had the check. The four HNSW producers did
not.

**Serialize index mutations.** `add_*`, `set_labels` and `save_index` run under
`write_mutex_`, because `add_*` reads `size()` as its first label and hnswlib's
`addPoint` *updates* an existing label instead of failing. Two concurrent adds
silently kept half the cells. Queries deliberately do not take that lock:
hnswlib's search is safe against a concurrent insertion, and serializing queries
would forfeit exactly the parallel throughput that releasing the GIL exists to
provide.

This module's record is **18 confirmed defects, 4 of them High or Critical**,
across the reviews and the 2026-09-11 audit and its second round. The two High
ones were silent: the unchecked `n * k` storing through a null pointer, and the
concurrent-mutation race leaving a perfectly valid index holding half its rows.
Round 2 also established by measurement rather than assertion that a zero-norm
vector produces no NaN anywhere (it scores a defined 0 and never displaces a
real neighbour), and that every CSR structural invariant already failed closed.

Two sanitizer notes for reproducing any of that: a UBSan suppression is inert
under `-fno-sanitize-recover=all`, so use `-fno-sanitize=alignment` for
hnswlib's known misaligned store; and `AlignedAllocator` rounds allocations up
to 64 bytes, so an ASan overread reproduction needs `posix_memalign` at the
exact size.

The Python binding layer's sanitizer coverage comes from Linux CI rather than
from a macOS host, and the reason is checkable rather than assumed: macOS System
Integrity Protection strips `DYLD_*` variables from a protected process before
it starts, and those variables are exactly what preload-based instrumentation of
a dynamically loaded Python extension requires.

## Benchmarking

See [`docs/BENCHMARKS.md`](docs/BENCHMARKS.md) for datasets, metrics, the
baseline comparison protocol, and the two methodology pitfalls that make naive
comparisons misleading: numba JIT warm-up dominating a single first call, and a
graph built from `X` being a different computation from scanpy's PCA default.
[`tests/benchmark_1m_cells.py`](tests/benchmark_1m_cells.py) runs it.
[`docs/AUDIT.md`](docs/AUDIT.md) is the module-level review log, with findings,
the evidence for each, and the items still open.

## License

MIT, under the root [`../../LICENSE`](../../LICENSE). Vendored hnswlib is
Apache-2.0 and its grant is reproduced in
[`../../THIRD_PARTY_LICENSES.md`](../../THIRD_PARTY_LICENSES.md).
