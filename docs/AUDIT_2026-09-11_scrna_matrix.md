# Audit — `scrna_matrix`, 2026-09-11

Adversarial pass over `modules/scrna_matrix`, following the playbook of
`AUDIT_2026-09-11_genomic_toolkit.md` (G2–G10),
`AUDIT_2026-09-11_cuttag_profiler.md` (C1–C6) and
`AUDIT_2026-09-11_fastq_stream.md` (A1–A8): reproduce before patching, add a
regression test per finding, and state the residuals rather than rounding them
off.

Findings are numbered **SM1–SM8** to avoid colliding with
`modules/scrna_matrix/docs/AUDIT.md`'s S1–S6 and with `SECURITY_HTTP`'s M-series.

**Convention, unchanged from the earlier passes:** nothing is recorded as a
finding unless it was reproduced first — a sanitizer trace, a wrong number, or
a file that would not load. Suspicions that did not reproduce are under
*Examined and cleared*, so a later reviewer does not re-derive them.

---

## 0. Scope correction, before anything else

The brief for this pass listed Matrix Market parsing, barcode/feature TSV
parsing, and a standalone CLI as the primary attack surface. **None of the
three exists in this module**, and no patch here invents them:

```
$ grep -rniI "matrixmarket|\.mtx|barcodes\.tsv|features\.tsv|mmread" \
      modules/scrna_matrix --exclude-dir=third_party
(no match outside a comment about HNSW label strings)
```

`scrna_matrix` is a header-only engine plus a pybind11 extension. It has no
`main()`, installs no binary, and is the one module the `pto-core` wheel builds
**OFF** (`pyproject.toml`, `PTO_BUILD_SCRNA_MATRIX = "OFF"`). Reading a `.mtx`
triple, a `barcodes.tsv` or a `features.tsv` is done by the caller — `scanpy` /
`anndata` in every in-tree example — and the module is handed the CSR buffers
that come out. So the parsing checks in the brief have no target here, and
asserting them "clean" would be asserting something about scanpy.

What the module does have, and what this pass attacked instead:

| Trust boundary | What crosses it |
|---|---|
| `Block_CSR::from_raw` / `::adopt` | a CSR triple plus caller-declared `n_rows` / `n_cols` |
| `from_scipy_csr` (pybind11) | three NumPy arrays of arbitrary dtype, and two dimensions |
| `to_aligned` (pybind11) | one NumPy array of arbitrary dtype |
| `HnswIndex` construction/query | `HnswParams` (M, ef, capacity, budget) and `k` |
| `HnswIndex::load_index` | a `.hnsw` payload and its `.meta` sidecar |
| `HnswIndex::save_index` | the filesystem's answer to a write |
| concurrency | Python threads, which the bindings release the GIL for |

The brief's *other* three themes do map onto this module, and are answered
below: 64-bit index bounds (SM6 and *Examined and cleared*), accumulator
wraparound (cleared — and why), and the C++/Python binding boundary (SM7, SM8).

---

## 1. Baseline

| Configuration | Result |
|---|---|
| Release, `ctest` (3 suites) | 3/3 pass |
| Release, Python bindings (`test_python_bindings.py`) | 44 passed |
| ASan+UBSan, `-fno-sanitize-recover=all`, `ctest` | **1/3 pass** — `test_matrix_ops` and `test_hnsw_persist` abort |
| ASan+UBSan, `-fno-sanitize=alignment`, `ctest` | 3/3 pass |

The two strict-mode aborts are **not findings**. Both are the already-recorded
misaligned `labeltype` store in vendored hnswlib
(`third_party/hnswlib/hnswalg.h:1308`), the NOT-PATCHED item of that
directory's `VERSION.txt` (numbered 5 after this pass added patch 4): fixing it would change
`size_data_per_element_` and invalidate every `.hnsw` file ever written.

Worth knowing for the next pass: the project suppresses it via
`docker/ubsan.supp`, and **a UBSan suppression does nothing under
`-fno-sanitize-recover=all`** — suppressions apply only to recoverable checks.
So the brief's exact flag combination cannot use the suppression file, and the
way to run this module strictly is `-fno-sanitize=alignment`, which disables
that one check and leaves the other ~20 UBSan checks armed. Every sanitizer
result below is from that configuration.

---

## 2. Findings

| # | Severity | Finding | Location |
|---|---|---|---|
| SM1 | High | `n * k` unchecked in all four HNSW producers → store through a null/wild pointer | `hnsw_index.hpp` (`query_dense`, `query_csr`, both one-shots) |
| SM2 | High | Two concurrent `add_*` calls take the same label base; half the rows are silently overwritten | `hnsw_index.hpp:294,317` |
| SM3 | Medium | `M` outside hnswlib's honoured range: UB at 0 and 1, silent clamp above 10000 → an index that saves and cannot load | `hnsw_index.hpp` ctor, `read_and_validate_meta` |
| SM4 | Medium | `save_index` checksums whatever reached the disk, sealing a short write as intact | `hnsw_index.hpp:437` |
| SM5 | Medium | Labels set before a later `add_*` are saved as-is; the pair is refused at load | `hnsw_index.hpp:437` |
| SM6 | Low | `n_rows + 1` wraps, so an empty `indptr` passes for `SIZE_MAX` rows → OOB read | `block_csr.hpp:247,278` |
| SM7 | Medium | `from_scipy_csr`'s copy path forcecasts unchecked → an out-of-range column index wraps to an **in-range** one | `src/python_bindings.cpp:216` |
| SM8 | Low | The AnnData helper pre-casts with `.astype`, bypassing `to_aligned`'s range check | `examples/scrna_anndata.py:75` |
| SM9 | Medium | Query breadth applied by writing hnswlib's shared `ef_` before each search → data race (TSan) | `hnsw_index.hpp:385` (`apply_query_ef`) |
| SM10 | Low | `capacity > INT32_MAX` accepted, though neighbour ids are `int32` | `hnsw_index.hpp` ctor / `read_and_validate_meta` |
| SM11 | Low | `to_aligned` checked float64 for NaN/Inf and not float32 | `src/python_bindings.cpp:324` |

SM9–SM11 are round 2 of this pass; the first two were residuals 1 and 2 of the
first round, and closing them is what §6 records.

Line numbers are as the code stood **before** the fix.

Three of the eight (SM1, SM7, SM8) are the same defect as an earlier,
already-closed finding, recurring in a sibling path — the pattern
`REVIEW_2026-08-15`'s retrospective names and this repository keeps rediscovering.
SM1 is REVIEW finding 4 (`n * k` overflow) in the four functions that did not
have it; SM7 and SM8 are AUDIT.md finding 3 (unchecked narrowing) in the copy
path and in the caller that feeds it.

---

### SM1 (High) — `n * k` unchecked in every HNSW producer

`KnnGraph` is a flat table of `n * k` slots, and row `i`'s neighbours are
written at `data() + i * k`. `build_knn_graph_bruteforce` has checked that
product since `REVIEW_2026-08-15` finding 4. `HnswIndex::query_dense`,
`HnswIndex::query_csr`, `hnsw_build_and_search` and
`hnsw_build_and_search_dense` each wrote `graph.neighbor_idx.assign(n * k, -1)`
with no check at all.

With 4 rows and `k = 2^62`, `n * k` is exactly 0: both `assign` calls succeed
against empty vectors, the search then runs normally, and `drain_into` writes
its results to `data() + i * k`.

**Evidence** — ASan+UBSan, `query_dense(data, n=4, k=2**62)`:

```
knn_graph.hpp:462:5: runtime error: store to null pointer of type 'std::int32_t'
SUMMARY: UndefinedBehaviorSanitizer: undefined-behavior knn_graph.hpp:462:5
```

Row 0 hits a null pointer; rows 1 and above are a write at an arbitrary
multiple of `k` from it. Reachable from Python in one line —
`scrna.build_knn_graph_dense(x, k=2**62)` — and the realistic form is a `k`
read from a config file or a notebook cell, not a hostile one.

**Fix.** One checked helper, in `knn_graph.hpp`, called by **every** producer
including brute force, so the definition exists once:

```cpp
[[nodiscard]] inline std::size_t checked_graph_slots(std::size_t n, std::size_t k) {
  if (k != 0 && n > std::numeric_limits<std::size_t>::max() / k) {
    throw std::invalid_argument(
        "k=" + std::to_string(k) + " with " + std::to_string(n) +
        " rows overflows the result size (n_rows * k exceeds SIZE_MAX)");
  }
  return n * k;
}
```

A `k` larger than the index is still legal and still fills the unused slots
with the `-1` sentinel; only the overflow is refused.

Regressions: `test_a_k_that_overflows_the_result_is_refused`
(`test_hnsw_persist.cpp`, all four entry points) and the Python test of the
same name (`build_knn_graph_dense`, `build_knn_graph` in both methods,
`HnswIndex.query_dense`).

---

### SM2 (High) — concurrent `add_*` silently discards half the index

`add_dense` reads `size()` as the label of its first new element, then inserts
`base + i` for each row. hnswlib's `addPoint` on a label that **already exists**
updates that element in place rather than failing
(`hnswalg.h:1232-1249`). Nothing serialised the read of `size()` against
another thread's insertions, and both binding entry points release the GIL.

**Evidence** — two `std::thread`s, 2,000 rows each, into a capacity of 4,000:

```
r7: added 2000 + 2000 rows; index size=2000 (expected 4000)
```

No exception, no warning, and the surviving index is a perfectly valid HNSW
graph over half the cells. The same shape reproduces from Python threads.

**Fix.** A `write_mutex_` held across base + capacity check + insertion in
`add_dense`, `add_csr`, `set_labels` and `save_index`. Queries deliberately do
**not** take it — hnswlib's `searchKnn` is safe against a concurrent
`addPoint`, and serialising queries would forfeit parallel query throughput,
which is the reason the bindings release the GIL in the first place. The mutex
is held behind a `unique_ptr` so the class stays movable (`load_index` returns
by value).

Regressions: `test_concurrent_adds_do_not_share_labels` (C++),
`test_concurrent_add_dense_keeps_every_row` (Python).

---

### SM3 (Medium) — `M` outside the range hnswlib honours

`HnswParams::M` went to hnswlib unvalidated. hnswlib validates it by *clamping*
at the top and not at all at the bottom:

| `M` | What happened |
|---|---|
| 0 | `maxM0_ = 0`; `mutuallyConnectNewElement` indexes a vector at offset −4 |
| 1 | `mult_ = 1 / log(1) = +inf`; `getRandomLevel`'s `(int) r` is UB |
| 20000 | capped to 10000 with a stderr warning, while our meta and budget keep 20000 |

**Evidence** — ASan+UBSan:

```
hnswalg.h:284:22: runtime error: inf is outside the range of representable values of type 'int'    # M=1
vector.h:434:27: runtime error: applying non-zero offset 18446744073709551612 to null pointer      # M=0
```

and for the clamp, a completed build whose own file will not load:

```
r4_mbig: save succeeded, meta M=20000
r4_mbig: load FAILED after a successful save:
         index metadata declares M=20000 but the payload was built with M=10000
```

That last one is the expensive shape: the atlas build finishes, the files are
written, and the failure appears when someone tries to use them.

**Fix.** `validate_params()` requires `M` in `[2, 10000]` — hnswlib's own
honoured range — at construction, and `read_and_validate_meta` applies the same
bound to a value read from a file, because `loadIndex` takes `mult_` from the
payload rather than recomputing it, so a poked `M` reaches the same undefined
level draw on the next `add_*`. `ef_construction` and `ef_search` are also
checked against the `uint32` header fields that were silently truncating them.

Regressions: `test_M_outside_the_range_hnswlib_honours_is_refused` (C++,
including a poked sidecar), `test_hnsw_index_refuses_M_outside_the_supported_range`
(Python).

---

### SM4 (Medium) — a short write was sealed as a valid index

`save_index` wrote the payload through hnswlib's `saveIndex`, which uses an
**unchecked** `std::ofstream`, and then computed the CRC and byte count *from
the file it had just written*. A write that failed partway therefore produced a
sidecar that faithfully described the truncated payload — and the checksum
whose documented job is to catch truncation matched it.

The binding docstring said this check was "the only check that catches a save
truncated by a full disk or a killed job". It was not: for a full disk it
sealed the damage instead.

**Evidence** — `RLIMIT_FSIZE` of 256 KiB in a forked child, with `SIGXFSZ`
ignored first (which is what CPython does at startup, so a Python caller gets
`EFBIG` rather than a killed process — the same shape as `ENOSPC` or a quota):

```
r2 child : save_index returned SUCCESS
r2 parent: payload on disk 262144 bytes; meta exists=1
r2 parent: load refused: Index seems to be corrupted or unsupported
```

hnswlib's own `loadIndex` catches it at load time by comparing structure to
file length, so this is a fail-late, not a wrong answer — but the failure is
reported hours later on a machine that no longer holds the source data, and
the repository's rule is that every write is checked.

**Fix.** After `saveIndex`, compare the bytes on disk against
`index_->indexFileSize()` — hnswlib's own count for this graph, computed
independently of the file — and refuse. The stale sidecar is removed *before*
the payload is written, so every failure path now leaves a payload with no
meta, which `load_index` already refuses. `out.flush()` became `out.close()`,
so a failure to flush the sidecar is caught too.

Regression: `test_a_short_write_is_not_sealed_as_a_valid_index`.

---

### SM5 (Medium) — labels that no longer match, saved anyway

`set_labels` requires one name per element *at the time it is called*. Adding
more elements afterwards left the labels describing a prefix. `save_index`
wrote `label_count = labels_.size()` beside `element_count = size()`, and
`load_index` refuses that pair.

**Evidence**

```
r3: save_index succeeded (size=20 labels=10)
r3: load FAILED after a successful save: index carries 10 labels for 20 elements
```

**Fix.** `save_index` refuses when labels are present and do not match the
element count, naming both numbers and saying how to fix it — at the point
where `set_labels` can still be called, rather than after the build is gone.
Nothing is written when it refuses.

Regressions: `test_labels_stale_after_a_later_add_are_refused_at_save` (C++),
`test_save_refuses_labels_that_no_longer_match` (Python).

---

### SM6 (Low) — `n_rows + 1` wraps, defeating the extent check

`check_extents` required `indptr.size() == n_rows + 1`. `n_rows` is
caller-supplied, and `SIZE_MAX + 1` is 0 — so an **empty but non-null** `indptr`
satisfied the check for `SIZE_MAX` rows, and `validate()` then read
`row_ptr[0]`, `row_ptr[1]`, … off the end of the buffer.

**Evidence** — ASan, with the buffer `posix_memalign`'d to exactly 4 bytes
(`AlignedAllocator` rounds every allocation up to 64, which is what masked
AUDIT.md finding 1 from an earlier ASan run):

```
ERROR: AddressSanitizer: heap-buffer-overflow
READ of size 4 at 0x609000000044 thread T0
    #0 scrna::Block_CSR<float>::validate() const   block_csr.hpp:302
    #1 scrna::Block_CSR<float>::adopt<int>(...)    block_csr.hpp:183
0x609000000044 is located 0 bytes after 4-byte region
```

Low, because the reachable callers cannot supply it: `to_aligned` returns a
null `data()` for an empty array, so the Python path throws "row_ptr must not
be null" first. It is the **C++** API — a public header — that admits it.

**Fix.** Both the factory check and `validate()` now test
`indptr.size() - 1 != n_rows` with an empty-span guard, an expression with no
addition to overflow.

Regression: `test_row_count_overflow_is_refused`.

---

### SM7 (Medium) — the copy path narrowed silently

`to_aligned` has range-checked its narrowing since AUDIT.md finding 3.
`from_scipy_csr`'s **copy** path did not: it handed the caller's arrays
straight to `py::array_t<...>::ensure` with `forcecast`, which is the same
conversion with no check.

**Evidence** — unpatched module, one row, one nonzero:

```
1 int64 column wrap  : ACCEPTED, shape=(1, 8589934592) nnz=1 (column 2**32+5 became 5)
2 same via to_aligned: REFUSED (to_aligned: int64 input value at index 0 does not fit in int32)
4 float index array  : ACCEPTED nnz=2 (0.9 and 1.9 truncated to columns 0 and 1)
```

Case 1 is the one that matters. The wrapped index is **in range**, so
`validate()` accepts the matrix and nothing downstream can tell: every
similarity involving that row is computed against a different gene. A matrix
2³³ columns wide is contrived, but `csr.indices` is int64 for any matrix past
2³¹ nonzeros, which is this library's stated target scale, and the two
entry points disagreeing about the same value is the defect regardless.

(Case 3, float64 values of 1e300, was already refused — but by `validate()`'s
finiteness check *after* `forcecast` turned them into `inf`, not by the
conversion.)

**Fix.** `check_index_dtype` / `check_value_dtype` run **before** the
forcecast: int64 indices are range-checked with the same `check_narrowing` that
`to_aligned` uses, float64 values likewise, integer values are accepted (int64's
maximum is ~9.2e18, far below `FLT_MAX`, so they cannot overflow — only lose
precision past 2²⁴, exactly as float64 does), and any other dtype is refused
with a message rather than converted on the caller's behalf.

Regression: `test_from_scipy_csr_refuses_lossy_narrowing`.

---

### SM8 (Low) — the AnnData helper cast before the check could run

`blockcsr_from_anndata` called `csr.indices.astype("int32")` and
`csr.data.astype("float32")` and passed the results to `to_aligned`. NumPy's
`astype` wraps, so the check existed and the caller stepped around it:

```
5 anndata helper pre-cast: astype gave 5 from 4294967301;
                           to_aligned on the original raised
```

**Fix.** The index arrays go to `to_aligned` at their own width; the values go
as float64 (lossless for every dtype AnnData stores there) so the one narrowing
happens on the checked path.

---

### SM9 (Medium) — query breadth was published through a shared member

`ef` — the search breadth — depends on `k`: the default is `max(2k, 64)`.
hnswlib v0.9.0's `searchKnn` takes it from the shared `ef_` member
(`std::max(ef_, k)`) and offers no per-query argument, so `apply_query_ef(k)`
called `index_->setEf()` immediately before each search. Two queries with
different `k` on one index therefore wrote and read that member concurrently.

**Evidence** — ThreadSanitizer, two `std::thread`s querying one 3,000-element
index with `k=2` and `k=400` (built without OpenMP so TSan can see the
threads; `docker/Dockerfile.tsan` is the OpenMP route and was not needed):

```
WARNING: ThreadSanitizer: data race (pid=16415)
  Write of size 8 at 0x000109d039d0 by thread T2:
    #0 scrna::HnswIndex::query_dense(...) const  hnsw_index.hpp:385
  Previous read of size 8 by thread T1:
    #0 hnswlib::HierarchicalNSW<float>::searchKnn(...) const  hnswalg.h:1414
```

The effect is mild — the loser searches at the other thread's breadth, so
recall changes rather than memory being corrupted — but it is a data race and
therefore UB, and it makes a query's result depend on what another thread
happened to be asking for.

**Fix.** A per-query `ef`, which needs the vendored library to accept one:
`third_party/hnswlib/VERSION.txt` **patch 4** moves `searchKnn`'s body into
`searchKnnWithEf(query, k, ef_override, isIdAllowed)` and leaves `searchKnn` as
a one-line delegation with `ef_override = 0`. `ef_` is read only when the
override is zero, so both the signature and the behaviour of the upstream API
are unchanged. (Adding the parameter to `searchKnn` itself does not compile:
it is a pure virtual of `AlgorithmInterface`, so the class becomes abstract —
recorded in VERSION.txt so the next person does not retry it.)

`HnswIndex` now computes `query_ef(k)` and passes it; nothing writes `ef_`
during a query. The configured value moved into an
`std::atomic<std::uint32_t>` beside the write mutex, so `set_ef()` concurrent
with a query is a defined read rather than a race.

**Verification.** The same TSan binary, same workload: **0 warnings**, exit 0.
Regression: `test_concurrent_queries_do_not_share_search_breadth` pins the
observable half — each thread's results now equal the answers it gets with
nothing else running.

---

### SM10 (Low) — a capacity larger than an `int32` neighbour id

`drain_into` narrows hnswlib's `labeltype` into the `std::int32_t` a `KnnGraph`
stores. An index holding more than `INT32_MAX` elements could therefore answer
with a wrapped, negative id — which is exactly the `-1` value that means *no
neighbour*, so the failure would read as a missing result rather than a wrong
one.

This was residual 2 of round 1, left unguarded because it is unreachable under
the default 8 GB budget (~300 GB of index) and could not be reproduced without
thrashing the host. That reasoning was wrong in one respect: `max_index_bytes`
is a caller-settable parameter, so "unreachable" is a property of a default
rather than of the code. The guard is now explicit:

```cpp
static constexpr std::size_t kMaxCapacity =
    static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max());
```

applied in the constructor **before any allocation** and to the `capacity` read
from a `.meta` sidecar. Verified with the budget raised to `SIZE_MAX`, so the
capacity check is demonstrably what refuses it:

```
capacity: threw HnswIndex: capacity 2147483649 exceeds 2147483647; neighbour ids
are int32 in a KnnGraph, so a larger index could return a wrapped id. Shard the
atlas instead.
```

Honest note: this is a guard on a path that still has no end-to-end
reproduction — the wrap itself was never observed, because observing it needs
~300 GB. What is tested is that the guard fires, in both the constructor and
the loader, before anything is allocated.

Regression: `test_capacity_beyond_int32_is_refused` (C++, constructor and a
poked sidecar), `test_hnsw_index_refuses_a_capacity_past_int32` (Python).

---

### SM11 (Low) — `to_aligned` checked float64 for NaN, not float32

`check_narrowing<double, float>` rejects non-finite values on the float64 path.
A float32 array went straight to `to_aligned_impl<float>` with no check, so
`to_aligned(np.array([1.0, np.nan], dtype="float32"))` returned an aligned
buffer containing NaN:

```
to_aligned NaN (float32): ACCEPTED, returned [1.0, nan, 3.0]
to_aligned Inf (float32): ACCEPTED, returned [inf]
```

Nothing wrong was ever *computed* — `Block_CSR::validate()` refuses non-finite
values on both the adopt and copy paths, which is the check that keeps NaN away
from the top-k comparator, and it did. But `to_aligned` is documented as the
checked path, and an array that will certainly be rejected later should be
rejected here, where the index in the message is still the caller's own.

**Fix.** `check_narrowing<float, float>` on the float32 path: same finiteness
test, no narrowing to perform. One extra pass over a buffer this function is
already copying.

**Test churn, stated plainly:** this broke the existing
`test_rejects_non_finite_values`, which built its NaN matrix *through*
`to_aligned` and expected the refusal to come from `BlockCSR.from_scipy_csr`.
The rejection now happens a line earlier. That test was rewritten to assert the
refusal at **both** boundaries — the aligned path and the copy path — rather
than at whichever one happened to catch it first. It is the only pre-existing
test this audit changed.

---

## 3. Examined and cleared

Recorded so a later reviewer does not repeat the work.

- **Integer bounds on `indptr` / `nnz`.** `index_type` is `int32`, and
  `validate()` refuses `nnz > INT32_MAX` with a message naming the 64-bit
  requirement — *before* any element is read, which
  `test_nnz_beyond_int32_is_refused_before_any_read` pins by handing in spans
  that claim 2³¹ elements over a one-element buffer (ASan would report a read
  if the guard ever moved after the first dereference). There is no silent
  32-bit truncation anywhere in the CSR path; the ceiling is a stated
  capability limit, not a wraparound.
- **UMI / expression accumulators.** The brief asked about cell-wise and
  gene-wise summation wrapping. **This module computes no such totals** — it has
  no per-cell or per-gene sum, only norms and dot products, which accumulate in
  `WideAcc<T>` (double) and narrow once through `checked_narrow`, enforced at
  build time by `test_wide_accumulation.cpp`'s source scan.
- **Duplicate barcodes in `set_labels`.** Accepted, deliberately. Labels are
  positional metadata, never a lookup key, and duplicate barcodes across
  libraries are ordinary in a concatenated atlas; refusing them would reject
  valid input to protect nothing. The count is enforced, which is the property
  the index depends on.
- **`load_index` header validation.** Re-read against the brief's "oversized
  declared dimensions" case: `dim` is bounded at 2²⁰, `capacity` is passed
  through the same `hnsw_index_bytes` budget as a built index, the
  saturating estimator refuses rather than wraps, and `validate_payload_header`
  pins the declared `dim` to the payload's true stride. A label length is
  capped at 1 MiB. All still hold; SM3 adds `M` to that list.
- **`BufferHandle` / capsule provenance.** The adoption gate (`is_adoptable`)
  was re-read against the four escapes it was built for; it still requires our
  own named capsule and exact extent, and a read-only view over a writeable
  base is still refused.
- **Non-finite input.** NaN/Inf reach neither the comparator nor the index:
  refused in `validate()`, in `normalize_dense_row`, and by `to_aligned`
  (float32 included, as of SM11).
- **Zero-norm vectors (an empty cell).** Measured, not assumed: a row of all
  zeros produces **no NaN anywhere** — not in the brute-force path, not in the
  HNSW index, not in a query. `inv_norm_from_sq` returns 0 for a non-positive
  norm, so the row is stored as a defined zero vector; `cosine_from_dot_wide`
  returns 0 when either norm is 0, so it scores 0 rather than `0/0`. Its
  distance to everything is `1 - 0 = 1`, the maximum, so it can never displace
  a real neighbour in anyone else's list, and its own list comes back filled at
  similarity 0. The degradation is "matches nothing", which is the right answer
  for a cell with no counts and is what `test_a_zero_norm_row_scores_zero_not_nan`
  (C++), `test_a_zero_norm_vector_produces_no_nan` (C++) and
  `test_an_empty_cell_scores_zero_not_nan` (Python) pin.
- **CSR structural invariants at the Python boundary.** `indptr[0] == 0`,
  `indptr` non-decreasing, `indptr[n_rows] == nnz`, `0 <= indices[j] < n_cols`,
  strictly increasing within a row, and `len(indices) == len(data)` were all
  already enforced in `validate()` and covered in C++ by
  `test_block_csr_validation` — and all nine malformed shapes were confirmed to
  fail closed through `from_scipy_csr` with messages that name the invariant.
  What was missing was a test saying so at that boundary, since the gather
  kernels index a dense scratch buffer with those values using unchecked
  hardware gathers; `test_malformed_csr_structures_fail_closed` is now that
  test.
- **`k` as a row stride.** `k = 0` gives an empty graph and `k` beyond the
  index gives `-1` sentinels in the unused slots — the documented contract,
  previously pinned only for brute force. It is now pinned for the HNSW paths
  too, and the *request* made to hnswlib is clamped to the number of elements
  held, so an oversized `k` no longer sizes a candidate heap for results that
  cannot exist. `k` itself is deliberately **not** clamped: it is the caller's
  output stride, and silently reshaping the result would be worse than the
  sentinels it replaces.

## 4. Residuals — what this pass did **not** close

1. ~~**`apply_query_ef` mutates shared state from a `const` query.**~~
   **Closed in round 2 — SM9.** It was reproduced after all: TSan does not need
   the OpenMP container, because building the probe without `-fopenmp` leaves
   `std::thread`s TSan can see, and the race is between two *queries*, not
   inside a parallel region. The fix was not the expensive one feared here
   (locking the query path); it was a per-query `ef`, which hnswlib v0.9.0
   genuinely does not offer and which patch 4 adds.
2. ~~**`capacity > INT32_MAX` is accepted.**~~ **Closed in round 2 — SM10.**
   The reasoning for leaving it ("unreachable under the default budget") was
   wrong in one respect: the budget is a caller-settable parameter. The guard
   is explicit now, though the wrap itself still has no end-to-end
   reproduction — see SM10.
3. **`n_cols` is not bounded.** A matrix may declare 2⁴⁰ columns with valid
   int32 indices; the brute-force path then allocates an `n_cols`-float dense
   scratch **per thread**. It fails as an allocation failure rather than
   silently, which is why it is not a finding, but the diagnostic is
   `std::bad_alloc` rather than a message naming the dimension. Not measured:
   the run would thrash the host.
4. **The vendored graph reader is still trusted.** Unchanged from AUDIT.md
   2026-08-16: the header checks catch accident (truncation, corruption,
   version skew, geometry mismatch) and now a bad `M`, but a deliberately
   hostile payload is consumed by unaudited hnswlib code. Provenance remains
   the control.
5. **The misaligned `labeltype` store is still UB.** See §1; fixing it breaks
   the on-disk format.
6. **No fuzzer.** AUDIT.md's open item 3 (a structure-aware fuzzer over
   `from_scipy_csr`) is still open. SM7 would have been found by one.

---

## 5. Verification

Every configuration was run after the patches; the first two were also run
before them, which is where the baseline in §1 comes from.

| Configuration | Before | After |
|---|---|---|
| Release, `ctest` | 3/3 | **3/3** |
| ASan+UBSan, `-fno-sanitize-recover=all -fno-sanitize=alignment` | 3/3 | **3/3** |
| ASan+UBSan, strict (alignment on) | 1/3 (vendored hnswlib) | 1/3 (unchanged, not ours) |
| `SCRNA_ENABLE_HNSW=OFF`, `ctest` | — | **3/3** |
| `SCRNA_BUILD_PYTHON=ON`, `ctest` | — | **4/4**, incl. `test_python_bindings` |
| CI's exact ASan+UBSan flags + `docker/ubsan.supp` + CI's final grep | — | **3/3**, no unsuppressed reports (§7.1) |
| `test_python_bindings.py` under pytest | 44 passed | **65 passed** (21 new, 1 rewritten) |
| `ruff check --select F` (the rule set `ci-cloud.yml` gates on) | — | **clean** |
| `mypy --ignore-missing-imports` | — | **clean** |
| `SCRNA_FORCE_ISA=scalar test_matrix_ops` | pass | **pass** (HNSW recall@10 vs exact 1.000) |
| TSan, two concurrent queries at different `k` | **data race** (SM9) | **0 warnings** |
| First-party warnings (`modules/*/(src\|include\|tests)/`) | 0 | **0** |

The Python-enabled tree is listed separately on purpose: `ctest` reports
"100% tests passed" over three suites when pybind11 is missing, so the run that
matters is the one whose list has **four** entries. This one does.

Each reproduction was re-run against the patched headers and now fails closed:

```
r1_query_dense   : k=4611686018427387904 with 4 rows overflows the result size
r3_stale_labels  : save_index: the index holds 20 elements but 10 labels ...
r4_m1 / m0 / mbig: HnswIndex: M=1 / 0 / 20000 is outside [2, 10000]
r6_rows_wrap     : Block_CSR::adopt: indptr.size()=0 must equal n_rows+1 (n_rows=18446744073709551615)
r7_concurrent    : added 2000 + 2000 rows; index size=4000 (expected 4000)
r2_short_write   : save_index: wrote 262144 of 809412 bytes to '...' (disk full, quota or
                   file-size limit?). No .meta was written ...  [no sidecar left behind]
```

Round 2, same form:

```
ef_race (TSan)   : 0 warnings, exit 0        [was: race, hnsw_index.hpp:385 vs hnswalg.h:1414]
capacity         : HnswIndex: capacity 2147483649 exceeds 2147483647; neighbour ids are int32
                   in a KnnGraph ... Shard the atlas instead.   [budget raised to SIZE_MAX, so
                   the capacity check is demonstrably what refuses it]
to_aligned NaN   : to_aligned: float32 input contains a non-finite value at index 1
to_aligned Inf   : to_aligned: float32 input contains a non-finite value at index 0
zero_norm_dense  : non-finite sims=0; empty cell's neighbours all at similarity 0.000
zero_norm_sparse : non-finite=0 on both the brute-force and HNSW paths
k_edges          : k=0 -> 0 slots; k=50 over 8 elements -> 50 slots, 7 filled, tail=-1
malformed CSR    : 9 of 9 shapes refused through from_scipy_csr, each naming its invariant
```

### Throughput

Interleaved A/B, 14 runs of each binary in both orders (`base` first for seven,
`fixed` first for seven, to cancel drift), Release `-O3`, nothing else running.
Two measurement sessions: round 1 after SM1–SM8, round 2 after SM9–SM11.

| Stage | Round 1 Δ | Round 2 Δ |
|---|---|---|
| `from_raw` validate, 200k rows × 1000 nnz (200M nnz) | +0.0% | −1.6% |
| brute-force k-NN, 4000 × 20000, 1500 nnz/row, k=15 | +0.4% | +0.4% |
| HNSW build, 100k × 50 dense | +0.9% | +1.2% |
| HNSW query, 100k × 50, k=15 | +0.9% | +0.3% |
| `save_index`, 100k × 50 | **+5.1%** | **+6.3%** |

Medians of 14 per binary per session. **Only the deltas are comparable between
the two sessions, not the absolute times**, and the reason is worth recording
rather than hiding: `perf-base` is a byte-identical binary in both, and its own
median went from 0.625 s to 2.009 s on the brute-force stage between sessions.
A 3.2× change in the unchanged control is the machine, not the code — measured
during the round-2 session as a load average of 8.9 from unrelated desktop
processes (a browser and the window server), with `pmset -g therm` recording no
thermal warning, so it is contention rather than throttling. What survives that
is the paired comparison, which is exactly why the harness interleaves the two
binaries instead of running all of one and then all of the other; a
"before, then after" script on this host would have reported the patches as a
3× slowdown.

Everything except `save_index` is inside run-to-run spread, in both sessions
and in both directions. The hot paths were not touched: SM1 adds one division
per *call* rather than per row, SM2's mutex is taken once per `add_*` rather
than per insertion, and SM9's per-query `ef` replaced a `setEf()` write with an
argument — if anything marginally less work, which is consistent with HNSW
query moving from +0.9% to +0.3%.

`save_index` is a real regression of 5–6%, reproduced in both sessions, and it
is the intended cost of SM4: `indexFileSize()` walks `cur_element_count`
link-list lengths to compute the expected size independently of the file. It is
6–9 ms on a 100k-element index, paid once per save, and it is the difference
between a save that reports success on a truncated file and one that does not.
Deliberately not optimised away by deriving the size from our own bookkeeping:
the whole point is that the figure comes from hnswlib's view of the graph
rather than from the same assumptions that wrote it.

### Files changed

```
modules/scrna_matrix/include/matrix/block_csr.hpp      SM6
modules/scrna_matrix/include/matrix/knn_graph.hpp      SM1 (shared helper)
modules/scrna_matrix/include/matrix/hnsw_index.hpp     SM1-SM5, SM9, SM10
modules/scrna_matrix/src/python_bindings.cpp           SM7, SM11 (+ two docstrings)
modules/scrna_matrix/examples/scrna_anndata.py         SM8
modules/scrna_matrix/third_party/hnswlib/hnswalg.h     SM9 (local patch 4)
modules/scrna_matrix/third_party/hnswlib/VERSION.txt   SM9 (patch 4 recorded)
modules/scrna_matrix/tests/test_matrix_ops.cpp         SM6, zero-norm + the int32 nnz ceiling
modules/scrna_matrix/tests/test_hnsw_persist.cpp       SM1-SM5, SM9, SM10, k edges, zero-norm
modules/scrna_matrix/tests/test_python_bindings.py     SM1-SM3, SM5, SM7, SM10, SM11,
                                                       CSR invariants, k edges, empty cell
modules/scrna_matrix/docs/AUDIT.md                     pointer entry
docs/AUDIT_2026-09-11_scrna_matrix.md                  this file
CLAUDE.md                                              input rules + trail entry
```

`hnswalg.h` is the one change outside first-party code, and the first
*functional* patch this tree has applied to hnswlib — the three before it fixed
races and an overread without touching an API. It adds `searchKnnWithEf` and
reduces `searchKnn` to a delegation; upstream's signature and behaviour are
unchanged, and `VERSION.txt` records it for the next version bump.

---

## 6. Round 2 — the five checks requested before sign-off

Recorded as asked-and-answered, because two of the five produced no defect and
that is a result worth being able to find later.

**1. Parameter bounds, no untested residuals.**
`capacity > INT32_MAX` is now refused in the constructor *and* in
`read_and_validate_meta`, with the memory budget raised to `SIZE_MAX` in the
test so the capacity check is demonstrably what refuses it — **SM10**. The
honest remainder is stated there: the guard is tested, the wrap it prevents is
not, because observing it needs ~300 GB.

For `k`: **rejecting or clamping it would have broken a documented, tested
contract**, so neither was done, and this is the one place I did not implement
the request as written. `k` is the row STRIDE of the result — `k = 0` yields an
empty graph and `k` beyond the index yields `-1` sentinels in the unused slots
(`KnnGraph`'s definition, `test_knn_k_exceeds_rows`, `test_knn_edge_cases`).
Clamping `k` would silently reshape the caller's output array; rejecting it
would break callers who legitimately ask for more neighbours than a small
cluster contains. What *was* wrong is that the number asked of hnswlib was the
unclamped `k`, sizing a candidate heap for results that cannot exist — that
request is now clamped to the element count (`fetch_count`), which is the
"cleanly clamps" half of the ask, applied where it changes no answer. Both
edges are now pinned at the HNSW entry points as well as for brute force
(`test_k_at_the_edges_of_the_index`, `test_hnsw_k_at_the_edges_of_the_index`),
and `n * k` overflow is refused outright (SM1).

**2. Floating-point poisoning.** No NaN is produced by a zero-norm vector, in
any path — measured, and now pinned by three tests. See *Examined and cleared*
for the mechanism and why the resulting "matches nothing" is the correct
degradation rather than a defect. The one real gap was `to_aligned` accepting
non-finite **float32** while checking float64 — **SM11**. `from_scipy_csr` was
already refusing NaN/Inf through `validate()` on both its paths, which is what
kept such values away from the top-k comparator; it now also refuses them with
a message naming the entry point.

**3. Sparse CSR structural invariants.** All four named invariants — plus
`indptr[0] == 0`, strict ordering within a row, and length agreement — were
already enforced in `Block_CSR::validate()` and covered in C++. Nothing
verified they survived the pybind11 layer, which is the boundary real input
crosses; all nine malformed shapes were confirmed to fail closed there with
messages naming the invariant, and `test_malformed_csr_structures_fail_closed`
now says so. **No code change was needed, and none was made** — the tests are
the deliverable.

**4. Thread safety in `apply_query_ef`.** Eliminated, not mitigated: the
function is gone, and query breadth is passed per call rather than published
through hnswlib's shared `ef_` — **SM9**, TSan-confirmed before (data race) and
after (0 warnings). This required local patch 4 to vendored hnswlib, which is
the trade-off worth reviewing explicitly: the alternative that needs no patch
is serialising queries behind the write mutex, which would forfeit the parallel
query throughput the bindings release the GIL to provide.

**5. Final verification.** §5, re-run in full after these changes.

---

## 7. Round 3 — CI parity and the dimension checks

### 7.1 The sanitizer job needs no CMake change, and should not get one

The concern raised was that `CI / ASan + UBSan` runs without this pass's local
`-fno-sanitize=alignment`, so the vendored hnswlib misaligned store would abort
the job unless CMake suppressed it on the affected targets. It does not, and a
CMake-level suppression would be the wrong repair. `.github/workflows/ci.yml`
already configures:

```
-fsanitize=address,undefined -fno-sanitize-recover=all \
-fsanitize-recover=alignment,pointer-overflow
```

with `UBSAN_OPTIONS=...suppressions=pto-core/docker/ubsan.supp` and a final
"No unsuppressed UBSan reports" step. That is a *better* arrangement than the
local one this audit used: the alignment check stays compiled in and armed, the
two upstream classes are made recoverable so the suppression file can apply to
them at all (a suppression cannot silence a compiled-in trap — the same
property §1 records), and anything the file does not cover still fails the job.
Disabling the check in CMake would remove it for first-party code too, which
is exactly what `docker/ubsan.supp`'s scope discipline forbids.

Two things were verified rather than assumed:

* **`ubsan.supp` matches by FILE, not line** — `alignment:hnswalg.h` and
  `pointer-overflow:hnswalg.h`. Local patch 4 (SM9) shifted line numbers inside
  that file but not its name, so the suppressions still apply. A line-anchored
  suppression would have silently stopped matching, and the job would have
  failed on the next push rather than here.
* **The module builds and passes under CI's exact flag string**, not just the
  local approximation: `RelWithDebInfo`, the flags above, `ASAN_OPTIONS` and
  `UBSAN_OPTIONS` as CI sets them, and CI's final grep. 3/3 suites pass and the
  gate reports no unsuppressed reports.

One local-only trap worth recording, because it cost a run and produced a
failure that looks like a code fault: **the sanitizer option parser splits on
spaces**, so `suppressions=/path/to/pto-core repo/...` silently breaks the
whole option string and every test aborts with
`AddressSanitizer: ERROR: expected '=' in UBSAN_OPTIONS` — no test having run.
This repository's checkout path contains a space; CI's does not. Copy the
suppression file somewhere space-free before pointing `UBSAN_OPTIONS` at it.

### 7.2 The three input and dimension checks

All three were probed against the current code. **None produced a defect**, and
the deliverable is the tests that now pin them.

| Check | Result |
|---|---|
| Dense query/add dimension vs index `dim` | already enforced at every entry point: `shape[1] != dim()` in the `add_dense`/`query_dense` bindings, `mat.cols() != dim()` in `add_csr`/`query_csr`. 1-D input is refused separately. |
| `dim == 0`, `capacity == 0` | refused in the `HnswIndex` constructor, before any allocation |
| `n_rows == 0`, `n_cols == 0`, `nnz == 0` | legal shapes: they construct and search, returning an empty or all-sentinel graph. No zero-byte allocation is indexed — `AlignedAllocator` returns `nullptr` for `n == 0` and every loop over it simply does not run. |
| Null `data` with a positive count | refused (`add_dense: null data`, `query_dense: null data`) |
| `NaN`/`Inf` in the sparse `data` buffer | refused on both paths — `validate()` catches it for adopted and copied buffers alike, and `to_aligned` catches it before that (SM11) |

One tightening followed from the last row rather than from a defect:
`check_value_dtype` now runs the finiteness check for **float32** values too, so
`from_scipy_csr` reports `from_scipy_csr: data contains a non-finite value at
index 1` — an index into the caller's own array — instead of deferring to
`validate()`'s `values must be finite; found NaN or Inf at nonzero index 1`,
which counts in CSR terms. Same refusal, earlier and in the caller's
coordinates.

New tests: `test_degenerate_dimensions_are_refused_or_empty` and
`test_empty_and_zero_dimension_matrices` (C++),
`test_query_and_add_require_the_index_dimension`,
`test_empty_and_zero_dimension_matrices_are_defined` and
`test_hnsw_refuses_zero_dim_and_zero_capacity` (Python).

### 7.3 The one CI step that cannot be run on this host

CI runs the Python binding suite under ASan with the runtime `LD_PRELOAD`ed
(`ci.yml`, "Test (Python bindings, runtime preloaded)"). **That step has no
local equivalent on this machine, and the reason is the host, not the code:**

* macOS System Integrity Protection strips `DYLD_*` from the environment before
  a protected binary starts. Verified from inside the process — both the venv's
  `python3` and its real base executable report
  `os.environ.get("DYLD_INSERT_LIBRARIES")` as absent, with
  `csrutil status: enabled`.
* Without the preload, ASan is loaded late via `dlopen` and aborts at import:
  `ERROR: Interceptors are not working ... Please launch the executable with
  DYLD_INSERT_LIBRARIES=...` — the very variable SIP has removed.

So the binding layer's ASan coverage comes from CI (ubuntu-24.04, clang-18)
and not from this host. What *was* verified locally is the whole C++ engine
under CI's exact sanitizer flags, and the binding suite un-instrumented (65
passed). If local coverage of the binding layer is wanted later, the way in is
an embedded-interpreter harness built **as** an ASan executable — the runtime is
then linked into the main binary and no preload is needed — rather than any
amount of work on the preload.

No public signature changed. The behaviour changes a caller can observe are all
refusals of input that was previously accepted and wrong: a `k` that overflows
`n * k`, an `M` outside `[2, 10000]`, a save whose labels no longer match, a
save that ran short, `SIZE_MAX` rows, and a NumPy dtype that would have been
converted silently.
