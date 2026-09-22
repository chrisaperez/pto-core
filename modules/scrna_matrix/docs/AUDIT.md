# Audit Trail

Static analysis and security review log for `scrna_matrix`. Newest entry first.
Each entry records what was examined, what was found, the evidence that
established it, and what was left open.

Convention: a finding is only recorded here if it was **reproduced** — with a
sanitizer trace, a failing assertion, or a measured wrong result — before a fix
was written. Suspicions that did not reproduce are listed under *Examined and
cleared* so a later reviewer does not re-derive them.

---

## 2026-09-11 — Adversarial audit (SM1–SM8)

Full record: [`../../../docs/AUDIT_2026-09-11_scrna_matrix.md`](../../../docs/AUDIT_2026-09-11_scrna_matrix.md).
Eight findings, all reproduced before patching and all closed with regression
tests; two High.

| # | Severity | Finding |
|---|---|---|
| SM1 | High | `n * k` unchecked in all four HNSW producers → store through a null/wild pointer (UBSan) |
| SM2 | High | Two concurrent `add_*` calls share a label base; hnswlib updates in place, so half the rows vanish silently |
| SM3 | Medium | `M` of 0 or 1 is UB inside hnswlib; `M > 10000` is clamped silently, giving an index that saves and will not load |
| SM4 | Medium | `save_index` checksummed whatever reached the disk, sealing a short write as intact |
| SM5 | Medium | Labels set before a later `add_*` were saved as-is, and the pair is refused at load |
| SM6 | Low | `n_rows + 1` wraps, so an empty `indptr` passed for `SIZE_MAX` rows → OOB read (ASan) |
| SM7 | Medium | `from_scipy_csr`'s copy path forcecast unchecked: an out-of-range int64 column index wrapped to an **in-range** one |
| SM8 | Low | `scrna_anndata.py` pre-cast with `.astype`, so `to_aligned`'s range check never saw the real value |
| SM9 | Medium | Query breadth was applied by writing hnswlib's shared `ef_` before each search — a data race between concurrent queries (TSan) |
| SM10 | Low | `capacity > INT32_MAX` was accepted, though `drain_into` narrows labels into an `int32` neighbour id |
| SM11 | Low | `to_aligned` checked float64 for NaN/Inf and not float32 |

Three of the eight are an already-closed defect recurring in a sibling path —
SM1 is REVIEW_2026-08-15 finding 4 in the four functions that lacked it; SM7
and SM8 are finding 3 of the entry below, in the copy path and in its caller.

SM9–SM11 are round 2 of the same pass, and two of them were **residuals of
round 1 that turned out to be closable**: the `ef_` race did not need the
OpenMP TSan container (building the probe without `-fopenmp` leaves
`std::thread`s TSan can see, and the race is between two queries, not inside a
parallel region), and "unreachable under the default memory budget" was the
wrong reason to leave `capacity` unguarded, because that budget is a parameter
the caller sets. Closing SM9 required the first *functional* patch to vendored
hnswlib — `searchKnnWithEf`, patch 4 in `third_party/hnswlib/VERSION.txt`.

Round 2 also answered three questions without finding a defect, and those
answers are in the main document rather than here: zero-norm vectors produce no
NaN anywhere (they score a defined 0 and can never displace a real neighbour),
every CSR structural invariant already failed closed through the Python
boundary, and `k` behaves as the documented row stride at both edges.

Two things about running the sanitizers here, which cost time to rediscover:
a UBSan **suppression does nothing** under `-fno-sanitize-recover=all`, so
`docker/ubsan.supp` cannot silence hnswlib's known misaligned store in that
mode — use `-fno-sanitize=alignment`. And `AlignedAllocator` rounds every
allocation up to 64 bytes, so an ASan reproduction of an overread must use
`posix_memalign` at the exact size; that is what masked finding 1 below from an
earlier run, and it masked SM6 on the first attempt too.

---

## 2026-08-15 — Monorepo stress test: fuzzing, sanitizers, thread churn

**Scope:** `scrna_matrix`, `fastq_stream`, `cuttag_profiler`.
**Method:** ASan+UBSan and TSan builds; 1,000,000-iteration structural fuzzer;
36 adversarial FASTQ/BGZF inputs; 48 adversarial BED/BAM/parameter inputs;
thread counts 1, 2, 3, 8, 16, 32; 1.3M-cell scale run.

Six defects found, all fixed. Two were latent in code this audit trail had
already reviewed — the earlier passes were not sufficient.

### S1 (scrna_matrix, high) — unsorted/duplicate column indices silently mis-scored

`sparse_dot` merge-joins two rows in lockstep and requires strictly increasing
column indices. Nothing enforced that. The two kernels then disagreed on the
same matrix:

| input | merge-join | gather | correct |
|---|---|---|---|
| two identical rows, unsorted indices | **0.80** | 1.00 | 1.00 |
| duplicated column index | 1.00 | **1.20** | 1.00 |

A cosine similarity of 1.20 is not a possible value. SciPy does not guarantee
sorted indices (hence `has_sorted_indices` / `sort_indices()`), so this was
reachable from ordinary input.

Fixed by extending `validate()` to walk per row and require strictly increasing
columns, sharing the pass that already checks range and finiteness. The error
names the row and points at `.sum_duplicates()` / `.sort_indices()`.
`_as_csr()` in `scrna_anndata.py` now normalises on the way in. Rejecting is the
only safe option: adopted buffers are read-only, so sorting them is not ours to do.

### S2 (scrna_matrix, high) — NaN manufactured from finite input

Found only at 300k fuzz iterations; 20k did not reach it. Input of 1e30 is a
perfectly ordinary float and passed the finiteness check added earlier, but
`1e30² = 1e60` saturates float to `inf`, and `inf/inf` is NaN — which then
reached the top-k comparator, the exact UB the earlier NaN work existed to
prevent. Boundary finiteness checking was necessary but not sufficient.

Fixed by accumulating norms in a wide type (`WideAcc<float> = double`), taking
the square root before narrowing (`sqrt(1e60) = 1e30` is representable), forming
the denominator wide, and adding a final finiteness guard in `cosine_from_dot`
that returns 0 rather than NaN. A similarity that cannot be represented is now
reported as no similarity — a documented degradation, never undefined behaviour.

### S3 (fastq_stream, medium) — unsigned underflow in mean quality

`phred_calculator.hpp:127` computed
`qsum_in += simd::sum_u8(qual, len) - len * kPhredOffset` in `uint64_t`.
Quality bytes below the +33 offset make the subtraction wrap: four bytes of
value 1..4 give `sum_u8 = 10` against `len*33 = 132`, and the reported mean
quality was **Q4611686018427387904**. The clamp three lines below already
guarded the histogram against exactly this input; only the aggregate was
exposed. Folded the sum into the loop that already walks every base, using the
same clamped value. Now reports Q0.00.

### S4 (cuttag_profiler, medium) — division by zero on `--bin-size 0`

`types.hpp:111` divides by `bin_size` in `bin_count()`; nothing validated it.
UBSan: "division by zero". Added the missing option check alongside the
existing `window_length() > 0` one.

### S5 (cuttag_profiler, medium) — signed overflow from BED coordinates

`signal_calc.cpp:79`, `anchor + right` with a coordinate near `INT64_MAX`.
UBSan: "9223372036854775806 + 2000 cannot be represented in type 'int64_t'".
The off-contig rejection immediately below could not help — the overflow
happened first and left it comparing a wrapped negative. Replaced with checked
arithmetic (`__builtin_add_overflow` / `__builtin_sub_overflow`) that saturates.

Note: the first fix here was itself wrong, and the fuzzer caught it on re-run.
A hand-rolled bound (`anchor > kMax - right`) overflows when `right` is
negative, which `--upstream -100` makes possible. The builtins are sign-agnostic;
negative half-windows are now also rejected outright.

### S6 (cuttag_profiler, medium) — OOM-kill from an unvalidated window

`--upstream 2000000000` asks for ~40M bins per region. The regions × bins
matrix of double is allocated up front, and the kernel SIGKILLed the process
before any diagnostic could be printed. Added a 10M-bin ceiling with an error
naming the window, the bin size, and the resulting bin count.

### Not ours: hnswlib misaligned store

`third_party/hnswlib/hnswalg.h:1202` stores an 8-byte `labeltype` to an address
that is not 8-byte aligned; UBSan flags it on every HNSW build. It is upstream's
packing scheme: `size_links_level0_` is 132 bytes, so the label offset lands on
a 4-byte boundary whenever the vector dimension is even (which includes both the
4-dim test fixtures and 50-dim PCA). Benign on x86-64 and arm64, where unaligned
scalar access is supported in hardware, but it is UB by the standard and it is
the only UBSan finding in the tree. Not patched — vendored code is kept
unmodified. It is worth an upstream issue.

### Sign-off

| Check | Result |
|---|---|
| Structural fuzz, 1,000,000 iterations, ASan+UBSan | clean; 1,030,870 invariant assertions |
| FASTQ/BGZF adversarial inputs (36) | no crash, hang, or sanitizer hit |
| BED/BAM/parameter adversarial inputs (48) | no crash, hang, or sanitizer hit |
| ASan+UBSan, threads 1/2/3/8/16/32 | 0 errors, 0 UBSan findings in our code |
| TSan, `std::thread` hammer, threads 2/3/8/16/32 | 0 warnings |
| 1.3M-cell HNSW scale run | 35.5 s, 1.02 GB, all similarities finite and in range |
| Notebook end-to-end | executes clean |
| Test suites, all three modules | pass |

**On TSan and OpenMP:** the OpenMP build reports races across
`#pragma omp parallel for` boundaries. These are artifacts — the mambaforge
libomp has zero `__tsan` undefined symbols, so TSan cannot see its barriers or
fork/join. Rather than assert that, the `std::thread` hammer removes OpenMP
entirely and exercises the same shared state (shared const matrix, magic-static
ISA dispatch, per-thread scratch isolation, `adopt()` refcount churn, concurrent
rejection); it is clean to 32 threads, where every report would be real.

**ASan and TSan cannot be combined** — clang rejects
`-fsanitize=address,thread`, as their shadow-memory models are incompatible.
They were run as separate passes; ASan+UBSan do combine.

---

## 2026-08-15 — Remediation: span-based factories (closes open item 1)

**Scope:** `include/matrix/block_csr.hpp`, `src/python_bindings.cpp`,
`tests/test_matrix_ops.cpp`.

Finding 1 of the audit below was fixed at the pybind11 layer, which left the
C++ API able to reproduce it. `from_raw` and `adopt` took a raw pointer per
buffer plus **one** `nnz` that was applied to all three — a signature that
cannot express "these three buffers have different lengths", so a short
`indices` array was read to the length of `values`.

Both factories now take `std::span`. The extent travels with the pointer, so
the three lengths are checkable against each other, and `check_extents()`
rejects any disagreement before an element is touched.

```diff
- static Block_CSR from_raw(std::size_t n_rows, std::size_t n_cols,
-                           const T* data, const index_type* indices,
-                           const index_type* indptr, std::size_t nnz);
+ static Block_CSR from_raw(std::size_t n_rows, std::size_t n_cols,
+                           std::span<const T> values,
+                           std::span<const index_type> indices,
+                           std::span<const index_type> indptr);

- template <typename Owner>
- static Block_CSR adopt(std::size_t n_rows, std::size_t n_cols,
-                        const T* data, const index_type* indices,
-                        const index_type* indptr, std::size_t nnz, Owner owner);
+ template <typename Owner>
+ static Block_CSR adopt(std::size_t n_rows, std::size_t n_cols,
+                        std::span<const T> values,
+                        std::span<const index_type> indices,
+                        std::span<const index_type> indptr, Owner owner);
```

New private helper, run first by both factories:

```cpp
static void check_extents(std::size_t n_rows,
                          std::span<const T> values,
                          std::span<const index_type> indices,
                          std::span<const index_type> indptr,
                          const char* who);
```

It enforces `indices.size() == values.size()`, `indptr.size() == n_rows + 1`,
and rejects a non-empty span with a null `data()` — which `std::span` itself
permits and which would otherwise reach the element loops in `validate()`.

Also added `values_span()`, `col_idx_span()` and `row_ptr_span()` accessors.
The raw `*_data()` accessors remain for C interop and pointer-identity checks,
but callers that read through them now have a bounds-carrying alternative, and
a matrix can be round-tripped back into the factories without a loose length.

**Call sites updated:** `from_scipy_csr` builds spans via a new `as_span<T>()`
helper that takes the extent from the NumPy array itself; its two hand-written
length checks were removed as now-redundant, leaving one enforcement point. The
`ndim() == 1` check stays — dimensionality is a NumPy concern a span cannot
express. Tests updated to pass spans throughout.

**Verification** — the original ASan reproduction, unchanged apart from the
call signature:

```
adopt: values=64, indices=2 (spans carry true extents)
REJECTED before any read: Block_CSR::adopt: indices.size()=2 must equal values.size()=64
```

No sanitizer report; the mismatch is now rejected at the API boundary. Two C++
regressions added (`test_span_factories_reject_mismatched_extents`,
`test_span_accessors_round_trip`) covering short `indices`, wrong-length
`indptr`, null-with-extent spans, and accessor round-trip. Full sweep passes:
arm64 Release/Debug, OpenMP at 1/4/12 threads, x86-64 cross with AVX2 forced,
`SCRNA_ENABLE_HNSW=OFF`, and 24 Python tests.

**Residual, stated precisely:** `std::span` constructed from `{ptr, len}` is
still only as truthful as `len`. What changed is that the three extents are now
mutually checked and no separate `nnz` can contradict `values.size()` — the
specific defect is gone, and a caller must now actively construct a lying span
rather than merely pass a plausible signature wrongly. Callers building spans
from a container or from `as_span()` cannot get it wrong at all.

---

## 2026-08-15 — Static analysis of engine, bindings, and AnnData layer

**Scope:** `include/matrix/*.hpp`, `src/python_bindings.cpp`,
`examples/scrna_anndata.py`, `tests/`, `CMakeLists.txt`.
Excludes `third_party/hnswlib` (vendored upstream, unmodified).

**Focus:** silent failures and error swallowing; memory, concurrency and
resource handling; numerical and boundary edge cases; type-safety and
lifetime bugs.

**Method:** manual review, AddressSanitizer, targeted reproduction harnesses,
and differential testing against reference implementations. The full report is
reproduced inline below rather than linked externally, so this document stays
self-contained.

**Outcome:** 9 findings, all fixed and verified in the same pass.

| # | Severity | Finding | Location |
|---|---|---|---|
| 1 | Critical | Out-of-bounds read from unvalidated buffer lengths | `src/python_bindings.cpp:60-89` |
| 2 | High | NaN breaks strict weak ordering → UB in `std::sort_heap` | `include/matrix/knn_graph.hpp:46-48` |
| 3 | High | Silent narrowing corruption in `to_aligned` | `src/python_bindings.cpp:113-118` |
| 4 | High | Exceptions escaping OpenMP regions → `std::terminate` | `include/matrix/knn_graph.hpp:343-368,425-448` |
| 5 | Medium | `eliminate_zeros()` deletes real nearest neighbours | `examples/scrna_anndata.py:119-120` |
| 6 | Medium | Integer overflow defeats the memory guard | `include/matrix/knn_graph.hpp:231` |
| 7 | Medium | Broad `except Exception` swallowed computation errors | `examples/scrna_anndata.py:149` |
| 8 | Low | Stale scratch contamination on non-growing resize | `include/matrix/simd_math.hpp:140-143` |
| 9 | Low | Degenerate neighbour padding hidden from the caller | `examples/scrna_anndata.py:138-142` |

Line numbers refer to the code **as it stood before the fix**.

---

### 1. Critical — Out-of-bounds read from unvalidated buffer lengths

`src/python_bindings.cpp:60-89`, reaching `include/matrix/block_csr.hpp:173,242`

`from_scipy_csr` derived `nnz` from `data.size()`, then read `indices` and
`indptr` using that length without checking those arrays were that long.

The failure mode worth remembering: this *looked* defended.
`Block_CSR::validate()` compares `values_.size()` against `col_idx_.size()`, but
on the adopt path those are the numbers just passed in — the comparison is
tautological and inspects no memory. Validation sat at a layer that could not
see the invariant, which reads as safety while providing none. The pybind layer
is the only place the true extents are known.

**Evidence** — AddressSanitizer, exactly-sized 64-byte-aligned allocations
(matching what NumPy hands over):

```
ERROR: AddressSanitizer: heap-buffer-overflow
READ of size 4 at 0x609000000048 thread T0
    #0 scrna::Block_CSR<float>::validate() const   block_csr.hpp:242
    #1 scrna::Block_CSR<float>::adopt<Owner>(...)  block_csr.hpp:173
0x609000000048 is located 0 bytes after 8-byte region
```

Reachable from Python. Whether it then throws or silently yields a corrupt
matrix depends on adjacent heap contents; both were observed. Note an earlier
ASan run on `AlignedVector` buffers stayed *quiet* — `AlignedAllocator` rounds
allocations up to 64 bytes and masked the overread. NumPy provides no such
padding, which is why the reproduction had to use `posix_memalign`.

**Fix**

```diff
+  if (data.ndim() != 1 || indices.ndim() != 1 || indptr.ndim() != 1) {
+    throw std::invalid_argument("from_scipy_csr: data, indices and indptr must each be 1-D");
+  }
+  const auto nnz = static_cast<std::size_t>(data.size());
+  if (static_cast<std::size_t>(indices.size()) != nnz) {
+    throw std::invalid_argument(
+        "from_scipy_csr: len(indices)=" + std::to_string(indices.size()) +
+        " must equal len(data)=" + std::to_string(nnz));
+  }
+  if (static_cast<std::size_t>(indptr.size()) != n_rows + 1) {
+    throw std::invalid_argument(
+        "from_scipy_csr: len(indptr)=" + std::to_string(indptr.size()) +
+        " must equal n_rows+1=" + std::to_string(n_rows + 1));
+  }
```

Regression: `tests/test_python_bindings.py::test_rejects_mismatched_buffer_lengths`

---

### 2. High — NaN breaks strict weak ordering → UB in `std::sort_heap`

`include/matrix/knn_graph.hpp:46-48`

The top-k comparator was `if (a.first != b.first) return a.first > b.first;`.
With NaN, `a != b` is true and every following comparison false, so `comp(a,b)`
and `comp(b,a)` are both false. NaN becomes "equivalent" to every value while
those values stay ordered among themselves — equivalence is non-transitive, so
this is not a strict weak ordering, and `std::push_heap` / `std::sort_heap`
become **undefined behaviour** rather than merely inaccurate.

**Evidence** — NaN/Inf input was accepted and propagated to the comparator:

```
accepted NaN/Inf input: nnz= 960
similarities NaN count: 64   inf count: 0
sample sims: [nan nan nan nan nan 0.4231826]
```

**Fix** — reject at the boundary, and make the ordering total anyway. The
finiteness check folds into the existing O(nnz) column-index pass, so it costs
no extra traversal.

```diff
   // block_csr.hpp, inside validate()
+  if constexpr (std::is_floating_point_v<T>) {
+    if (!std::isfinite(vals[i])) {
+      throw std::invalid_argument(
+          "values must be finite; found NaN or Inf at nonzero index " + std::to_string(i));
+    }
+  }
```

```diff
   // knn_graph.hpp, BetterNeighbor::operator()
+  const bool a_nan = std::isnan(a.first), b_nan = std::isnan(b.first);
+  if (a_nan || b_nan) {
+    if (a_nan != b_nan) return b_nan;   // any real score beats NaN
+    return a.second < b.second;         // both NaN: deterministic by index
+  }
   if (a.first != b.first) return a.first > b.first;
   return a.second < b.second;
```

Regression: `test_rejects_non_finite_values`, `test_hnsw_dense_rejects_non_finite`

---

### 3. High — Silent narrowing corruption in `to_aligned`

`src/python_bindings.cpp:113-118`

`forcecast` wraps without raising. Measured:

| Input | Produced |
|---|---|
| `int64  3_000_000_000` | `int32 -1294967296` |
| `float64 1e300` | `float32 inf` |

SciPy switches to int64 `indptr` whenever nnz exceeds 2³¹ — exactly this
library's stated target scale. A wrapped `indptr` is a structurally invalid CSR;
an `inf` value poisons every similarity it participates in.

**Fix** — range-check before converting, and refuse rather than wrap.

```diff
-  if (src.dtype().is(py::dtype::of<std::int64_t>())) return to_aligned_impl<std::int32_t>(src);
+  if (src.dtype().is(py::dtype::of<std::int64_t>())) {
+    check_narrowing<std::int64_t, std::int32_t>(src, "int64 input");
+    return to_aligned_impl<std::int32_t>(src);
+  }
```

`check_narrowing<Src, Dst>` rejects non-finite values and anything outside the
destination range, with an actionable message
(`"...does not fit in int32; this matrix needs 64-bit indices"`). In-range
narrowing still succeeds.

Regression: `test_to_aligned_refuses_lossy_narrowing`

---

### 4. High — Exceptions escaping OpenMP regions → `std::terminate`

`include/matrix/knn_graph.hpp:343-368, 425-448`

`addPoint` and `searchKnn` sat inside `#pragma omp parallel` with no `try`.
hnswlib throws `std::runtime_error` from roughly a dozen sites reachable from
those calls, including `"Not enough memory"` at `hnswalg.h:128,140` — precisely
the failure expected when indexing 1.3M cells, which is the workload this path
exists to serve.

An exception leaving an OpenMP structured block is undefined behaviour
(OpenMP 5.2 §2.1). libomp calls `std::terminate`: the process dies with no
catchable Python error.

**Fix** — capture the first exception, drain the loop, rethrow outside the
region.

```diff
+  template <typename Body>
+  inline void parallel_for_checked(std::size_t n, int num_threads,
+                                   const char* what, Body&& body) {
+    std::atomic<bool> failed{false};
+    std::string message;
+    #pragma omp parallel num_threads(num_threads)
+    {
+      #pragma omp for schedule(dynamic, 256)
+      for (std::size_t i = 0; i < n; ++i) {
+        if (failed.load(std::memory_order_relaxed)) continue;  // drain, don't break out
+        try { body(i); }
+        catch (const std::exception& e) {
+          bool expected = false;
+          if (failed.compare_exchange_strong(expected, true)) {
+            #pragma omp critical(scrna_hnsw_error)
+            message = e.what();
+          }
+        } catch (...) { failed.store(true, std::memory_order_relaxed); }
+      }
+    }
+    if (failed.load(std::memory_order_relaxed))
+      throw std::runtime_error(std::string(what) + " failed: " +
+                               (message.empty() ? "unknown error" : message));
+  }
```

Scratch buffers are additionally cleared on the throw path, so a partial failure
cannot leak one row's values into later iterations on that thread.

---

### 5. Medium — `eliminate_zeros()` deletes real nearest neighbours

`examples/scrna_anndata.py:119-120`

A cosine distance of exactly `0.0` is legitimate: it means two **identical
cells**. `setdiag(0.0)` followed by `eliminate_zeros()` cannot distinguish that
from the zeroed diagonal, so it silently deleted the strongest edges in the
graph. Duplicate expression profiles are routine in scRNA-seq — doublets,
low-complexity cells, technical replicates.

**Evidence** — 4 identical cells, `k=10`:

| | Before | After |
|---|---|---|
| Neighbours stored for cells 0–3 | 7 of 10 | 10 of 10 |
| Rows with fewer than k | 4 of 60 | 0 of 60 |
| Zero-distance edges preserved | 0 | 12 |
| Self-loops | none | none |

**Fix** — filter self-edges structurally, never by value.

```diff
-  keep = cols >= 0
-  matrix = sp.csr_matrix((dist[keep], (rows[keep], cols[keep])), shape=(n_obs, n_obs))
-  matrix.setdiag(0.0)
-  matrix.eliminate_zeros()
+  keep = (cols >= 0) & (cols != rows)
+  matrix = sp.csr_matrix((dist[keep], (rows[keep], cols[keep])), shape=(n_obs, n_obs))
+  # coo->csr sums duplicate (row, col) pairs; top-k emits distinct neighbours, so
+  # this asserts a future approximate backend cannot silently double an edge.
+  assert matrix.nnz == int(keep.sum()), "duplicate neighbour entries were summed"
```

---

### 6. Medium — Integer overflow defeats the memory guard

`include/matrix/knn_graph.hpp:231`

`hnsw_index_bytes` ended in `return n * (data + links0 + ...)`, which wraps. That
figure gates an allocation, so overflow waves the request straight past the
`max_index_bytes` check.

**Evidence**

```
hnsw_index_bytes(2**60, 1024)  ->  0          # guard passes, allocation attempted
hnsw_index_bytes(2**58,   64)  ->  wrapped
```

An overflowing bounds check is worse than no bounds check: it reads as
protection while removing it.

**Fix** — saturating arithmetic, so an unrepresentable index always *fails* the
guard.

```diff
+  constexpr std::size_t kMax = std::numeric_limits<std::size_t>::max();
+  if (dim != 0 && dim > kMax / sizeof(float)) return kMax;
   ...
+  if (per_element != 0 && n > kMax / per_element) return kMax;
   return n * per_element;
```

Regression: `test_hnsw_index_bytes_saturates_instead_of_wrapping`

---

### 7. Medium — Broad `except Exception` swallowed computation errors

`examples/scrna_anndata.py:149`

One `except Exception` wrapped both the optional import *and* the call to
scanpy's fuzzy-simplicial-set routine. A genuine failure inside that call — bad
shapes, non-finite distances — silently substituted a **different algorithm**,
surfacing as subtly wrong connectivities rather than an error.

**Fix** — only the import is optional; computation errors propagate.

```diff
-  try:
-      from scanpy.neighbors._connectivity import umap as _umap_conn
-      conn = _umap_conn(knn_idx, knn_dist, n_obs=n_obs, n_neighbors=k + 1)
-      return sp.csr_matrix(conn), "scanpy.umap"
-  except Exception:
-      ...fallback...
+  try:
+      from scanpy.neighbors._connectivity import umap as _umap_conn
+  except (ImportError, AttributeError):
+      _umap_conn = None
+
+  if _umap_conn is not None:
+      conn = _umap_conn(knn_idx, knn_dist, n_obs=n_obs, n_neighbors=k + 1)
+      return sp.csr_matrix(conn), "scanpy.umap"
+
+  warnings.warn("...NOT the UMAP connectivity...", RuntimeWarning, stacklevel=3)
```

---

### 8. Low — Stale scratch contamination on non-growing resize

`include/matrix/simd_math.hpp:140-143`

`DenseScratch::resize()` documented "all zero" but only zeroed when *growing*.
Reusing one scratch across two matrices left stale values in place, and the
gather kernels read every position a row's column indices point at — so residue
is consumed as data rather than caught.

```diff
-  if (buf_.size() < n_cols) buf_.assign(n_cols, 0.0f);
+  if (buf_.size() < n_cols) {
+    buf_.assign(n_cols, 0.0f);
+  } else {
+    std::fill(buf_.begin(), buf_.begin() + static_cast<std::ptrdiff_t>(n_cols), 0.0f);
+  }
```

Same complexity; the active region is now always zero on entry.

---

### 9. Low — Degenerate neighbour padding hidden from the caller

`examples/scrna_anndata.py:138-142`

When `k >= n_obs - 1`, unused slots were padded with repeated self-references.
Those duplicated zero-distance entries distort UMAP's local-connectivity
estimate, and nothing signalled the result was degenerate. Now emits a
`RuntimeWarning` naming the affected cell count.

---

### Examined and cleared

Checked and found sound; recorded so a later reviewer does not repeat the work.

- **`TODO` / `FIXME` / stub / mock values** — none in any production path.
- **`catch(...)` swallowing, bare `except:`** — none.
- **`BufferHandle` lifetime** — `data()` and `size()` are computed, not cached,
  so copying or moving a handle cannot leave a stale pointer behind.
- **`ScrnaCosineSpace` lifetime** — declared before the `HierarchicalNSW` that
  stores its `dist_func_param_`, so destruction order is correct; the index dies
  first.
- **`to_numpy` capsule keep-alive** — correctly pins the `KnnGraph` behind the
  exported NumPy arrays; verified by dropping the graph and re-reading.
- **GIL handling** — `py::call_guard<gil_scoped_release>` is instantiated after
  argument casting, so `std::string` and array arguments are safe; return-value
  conversion happens after the guard is destroyed. Verified by a thread-progress
  test.
- **`row_values()` throwing inside an OpenMP loop** — bounded by `i < mat.rows()`
  at every call site, so unreachable in practice.

### Open items

1. ~~**The fix for finding 1 lives in the pybind11 layer.**~~
   **Closed 2026-08-15** — see the entry below. The factories now take
   `std::span`, so the C++ API no longer admits the mismatch either.
2. **The ASan reproduction and the x86-64 cross-build are not reproducible from
   the repo.** Both were run from throwaway scratch files. Folding the ASan probe
   into the test suite and adding a CI-shaped cross-build script would keep those
   checks alive.
3. **No fuzzing.** `from_scipy_csr` is the untrusted-input boundary and would
   suit a structure-aware fuzzer over (data, indices, indptr, n_rows, n_cols).

### Verification

Every fix verified across all build configurations:

| Configuration | Result |
|---|---|
| arm64 Release | pass |
| arm64 Debug (assertions live) | pass |
| OpenMP, 1 / 4 / 12 threads | pass |
| x86-64 cross-build, `SCRNA_FORCE_ISA=avx2` | pass |
| `SCRNA_ENABLE_HNSW=OFF` | pass |

Python suite grew 19 → 24 tests; the five additions are regressions pinned to
findings 1, 2, 3 and 6. The integration notebook re-executes clean end to end.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```

Note that every benchmark and test on record ran the **scalar** SIMD path; this
host is arm64 and has no AVX2/AVX-512. The AVX2 kernel is correctness-verified
under Rosetta emulation, the AVX-512 kernel is compile-verified only.

---

## 2026-08-16 — New trust boundary: `HnswIndex::load_index`

Index serialization (`include/matrix/hnsw_index.hpp`) introduces the first API
in this module that takes a **file** and lets its contents decide how much
memory to allocate and how to stride through a buffer. Recorded here as a
deliberate boundary rather than discovered later as a finding.

### The mechanism

`hnswlib::HierarchicalNSW::loadIndex` (third_party/hnswlib/hnswalg.h:755) takes
a `SpaceInterface*` from the caller and the element geometry from the file, and
never checks that they agree. `data_size_` is overwritten from our space while
`size_data_per_element_`, `offsetData_` and `label_offset_` come from the
payload. A payload written at dim 8, loaded with a space constructed for dim
4096, therefore reads at strides the allocation does not cover. No error, no
assertion — out-of-bounds reads inside the loaded blob.

### What is checked, and where

`read_and_validate_meta()` runs before hnswlib is handed the file at all:

| Check | Catches |
|---|---|
| magic + `format_version` | wrong file entirely; version skew |
| `metric`, `normalized` | an index whose distances mean something else |
| `scalar_layout` | non-IEEE754 or big-endian producer |
| `dim` in (0, 2^20] | garbage or hostile dimension |
| `element_count <= capacity` | internally inconsistent header |
| `hnsw_index_bytes(capacity, dim, M) <= max_index_bytes` | a header claiming 2^32 elements, i.e. the same memory-exhaustion shape finding 6 closed for the *built* path |
| `validate_payload_header()` | **the dim/payload mismatch above** |
| CRC-32 + byte count of the payload | truncated or corrupted save |

`validate_payload_header()` is the one worth singling out. Every other check
only establishes that the sidecar is *self-consistent and plausible* — a
sidecar declaring dim 4096 for an 8-dim payload passes all of them. It reads
the leading 96 bytes of the hnswlib payload and requires
`label_offset_ - offsetData_ == dim * sizeof(float)`, which is the single
relation that pins the payload's true stride to our declared dimension, plus
agreement on `max_elements_`, `cur_element_count` and `M_`.

This was found by `test_dim_mismatch_is_refused` in
`tests/test_hnsw_persist.cpp` failing against the first implementation, which
had the plausibility checks but not the cross-check. Exactly the
REVIEW_2026-08-15 pattern: the obvious validation was written, and the sibling
path it did not cover was live.

### Residual (open)

**Loading an index file is equivalent to trusting whoever produced it.** The
checks above cover accident — truncation, corruption, version skew, a
mismatched pair of files — and the specific geometry confusion described here.
They do **not** make the load safe against a *deliberately malicious* payload:
the graph body is consumed by vendored hnswlib code that we have not audited,
and a hostile link array can send its traversal out of bounds regardless of
what our header says. The offsets `validate_payload_header()` reads are
themselves attacker-controlled; they are cross-checked against each other, not
independently established.

Consequences to hold to:

1. Do not make `load_index` reachable from any unauthenticated entry point.
   In particular it must not be wired into `cuttag_profiler`'s embedded HTTP
   server, whose finding H1c (no authentication) is still open.
2. Treat `.hnsw` files like any other executable artefact from a third party:
   provenance is the control, not the format.
3. A structure-aware fuzzer over the payload would be the way to move this from
   "trusted producer" to "hardened", and belongs with open item 3 above
   (`from_scipy_csr` fuzzing) as the same piece of work.

### Verification

`tests/test_hnsw_persist.cpp`, 12 cases, arm64 Release and Debug, scalar SIMD
path. Round-trip equality is asserted on **exact** float comparison, not a
tolerance: the same graph over the same vectors must yield bit-identical
distances, and a tolerance there would hide a partially-written payload.
