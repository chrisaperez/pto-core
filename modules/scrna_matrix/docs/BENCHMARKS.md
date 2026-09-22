# Benchmarking Protocol

## Choosing a method

Two k-NN backends, with genuinely different operating ranges:

| | `method="brute"` | `method="hnsw"` |
|---|---|---|
| Result | Exact | Approximate |
| Deterministic | Yes | No (insertion-order dependent; pin `num_threads=1`) |
| Complexity | O(n_obs²) | ~O(n_obs log n_obs) |
| Input | Sparse CSR, native | Dense only — densified internally |
| Memory | O(n_obs·k) + O(n_genes)/thread | n_obs × n_dims × 4 B, sparsity-independent |
| Use when | ≲ 20-50k cells, or ground truth needed | ≳ 50k cells **on a reduced representation** |

The decisive property of hnswlib is that it indexes **dense fixed-dimension
vectors**. It has no sparse representation, so gene-space input must be
densified — which costs memory proportional to `n_genes` and makes every
distance an `n_genes`-wide dense operation instead of a ~70-element sparse
gather. Measured below, HNSW is **4-20x slower than exact brute force in
1,000-dimensional gene space**, and 11x faster on a 50-dimensional PCA
embedding. The representation, not the algorithm, decides which wins.

Rule of thumb: exact up to a few tens of thousands of cells; HNSW on a PCA
embedding beyond that; never HNSW on wide gene-space input.

## Datasets

| Name | Source | Cells | Notes |
|---|---|---|---|
| `pbmc3k` | 10x Genomics (via `sc.datasets.pbmc3k()`) | ~2,700 | Downloaded on demand by scanpy; used by `examples/scanpy_integration.ipynb` |
| `synthetic` | `benchmark_1m_cells.py --make-synthetic N` | any | Generated locally, not committed. scRNA-seq-shaped sparsity |
| `1k_hgmm_v3` | 10x Genomics | ~1,000 | Fast smoke-test dataset for CI-scale benchmarking |
| `1.3M_cells_neuron` | 10x Genomics "1.3 Million Brain Cells" | ~1,300,000 | Aspirational target; see the scope note above |

The two 10x datasets must be staged locally (not committed to this repo) before
running `tests/benchmark_1m_cells.py`.

## Metrics

- **Peak RSS memory (GB)** — `resource.getrusage(RUSAGE_SELF).ru_maxrss` per run, or `/usr/bin/time -v` wrapping the process.
- **k-NN graph construction wall-clock time (s)** — measured strictly around the `build_knn_graph` / `sc.pp.neighbors` call, excluding I/O.
- **Cache miss rate** — `perf stat -e cache-misses,cache-references` (Linux only; not available on macOS/M-series, use `xctrace`/Instruments as a substitute there).

### Measurement pitfalls found in practice

1. **numba JIT warm-up dominates a single first call.** scanpy's
   fuzzy-simplicial-set routine and `pynndescent` both JIT-compile on first use —
   several seconds, unrelated to the search itself. Warm both backends before
   timing, or the numbers are meaningless. An unwarmed run of the notebook
   attributed 5.1 s to a search that actually takes 0.037 s.
2. **`sc.pp.neighbors` defaults to a 50-component PCA**, not `X`. Comparing a
   graph built from `X` against scanpy's default compares two different
   representations: measured neighbour overlap was 0.26 that way versus 1.00 with
   `use_rep="X"`. Always pass `use_rep="X"` for a like-for-like comparison.
3. **Exact vs approximate.** This engine is exact; `sc.pp.neighbors` is not.
   Disagreement is a statement about the approximation's recall, not necessarily
   an error on either side.

## Baselines

- `scanpy.pp.neighbors()` (default UMAP-style approximate k-NN via `pynndescent`), invoked with `use_rep="X"`, `metric="cosine"`
- `scipy.sparse` default linear algebra routines (dense dot product loop over `csr_matrix` rows)

## Procedure

1. Load the dataset into both representations (`AnnData`/`scipy.sparse.csr_matrix` and `Block_CSR`).
2. Warm each backend with one discarded call.
3. Run each backend 5 times; report median wall-clock and peak RSS.
4. Run once under `perf stat` for cache-miss counts (separate invocation — `perf` overhead skews timing).
5. Record hardware: CPU model, core count, RAM, whether OpenMP thread count was pinned (`OMP_NUM_THREADS`), and which SIMD kernel was selected (`scrna_matrix_py.simd_isa()`).

## Measured results

Hardware: Apple M4 Pro (12 cores), macOS 15. Built Release with OpenMP
(12 threads). **SIMD kernel: `scalar`** — arm64 has no AVX2/AVX-512, so these
numbers reflect the portable path only. x86 numbers with the gather kernels
active are not yet collected; see "Unverified" below.

### Real data — `pbmc3k`, 2,700 cells x 1,000 HVGs, k=15

Reproduced by `examples/scanpy_integration.ipynb` (committed with live outputs).

| Backend | Wall clock | Neighbour overlap |
|---|---|---|
| `scrna_matrix` (exact, search only) | 0.037 s | — |
| `scrna_matrix` (search + AnnData conversion) | 0.040 s | — |
| `sc.pp.neighbors` (approximate, `use_rep="X"`) | 0.244 s | 1.0000 vs ours |

Brute force wins here because 2,700 cells is small enough that the quadratic term
is cheap, the working set is L2-resident, and the search parallelises perfectly.
This does not extrapolate — see below.

### Scaling — synthetic, 1,000 genes, 6% density (~60 nnz/cell), k=15

| Cells | nnz | Wall clock | vs previous | ns per n² |
|---|---|---|---|---|
| 2,500 | 150,000 | 0.027 s | — | 4.36 |
| 5,000 | 300,000 | 0.096 s | 3.51x | 3.82 |
| 10,000 | 600,000 | 0.351 s | 3.67x | 3.51 |
| 20,000 | 1,200,000 | 1.465 s | 4.17x | 3.66 |
| 40,000 | 2,400,000 | 5.454 s | 3.72x | 3.41 |

The near-constant final column is the quadratic law holding across a 16x range.
Extrapolating at ~3.5 ns/n²: **1.3M cells ≈ 1.6 hours** for the exact path at
this density. Not an out-of-memory failure — an unacceptable runtime, and the
reason the HNSW backend exists.

### HNSW — 1,000-dim gene space, k=15 (the case where it loses)

Synthetic, 6% density. Same data as the exact table above.

| Cells | Exact | HNSW | Ratio | HNSW recall |
|---|---|---|---|---|
| 2,500 | 0.027 s | 0.518 s | **0.05x** | 0.988 |
| 5,000 | 0.092 s | 1.597 s | 0.06x | 0.960 |
| 10,000 | 0.366 s | 4.441 s | 0.08x | 0.885 |
| 20,000 | 1.399 s | 9.744 s | 0.14x | 0.767 |
| 40,000 | 5.618 s | 22.312 s | 0.25x | 0.617 |

HNSW is *slower everywhere here*, and recall degrades as n grows. Each HNSW
distance touches 1,000 dense dimensions where the exact path touches ~60
non-zeros, so HNSW's smaller number of comparisons does not pay for their cost.
Recall degradation is partly the data: random sparse vectors have no cluster
structure for an ANN graph to exploit.

### HNSW — 50-dim embedding with cluster structure, k=15 (the case it is for)

| Cells | Exact | HNSW | Speed-up | HNSW recall |
|---|---|---|---|---|
| 10,000 | 0.25 s | 0.11 s | 2.4x | 1.000 |
| 25,000 | 1.57 s | 0.25 s | 6.4x | 0.999 |
| 50,000 | 6.32 s | 0.57 s | **11.2x** | 0.997 |
| 100,000 | ~25 s (extrapolated) | 1.40 s | ~18x | not measured |
| 200,000 | ~100 s (extrapolated) | 3.65 s | ~28x | not measured |
| **1,300,000** | **~1.6 h (extrapolated)** | **37.1 s** | **~155x** | not measured |

The 1.3M-cell row is a real run, not an extrapolation: **37.1 s wall clock at
1.01 GB peak RSS**, including index construction and all 1.3M queries. Recall is
marked "not measured" above 50k because the exact baseline needed to measure it
is precisely what is infeasible at that size.

## Target Results (per manuscript claims, §1.4)

- 20x wall-clock speedup over `scanpy.pp.neighbors()` / `scipy.sparse` baselines on the 1.3M-cell dataset.
- 60% peak RSS reduction over the same baselines.

**Status: the 1M+ scalability goal is met; the two comparative claims are still
unmeasured.**

What is now demonstrated: 1.3M cells in 37.1 s at 1.01 GB, which was previously
projected at ~1.6 hours. That closes the feasibility gap.

What is *not* demonstrated: neither headline number is a comparison against
`scanpy.pp.neighbors` at 1.3M cells, and that comparison has not been run. Note
also that `sc.pp.neighbors` on PCA is itself approximate and near-linear, so it
is a peer of `method="hnsw"` rather than a baseline it should be expected to beat
by 20x. Until a like-for-like 1.3M-cell run against scanpy exists — same
representation, same k, both JIT-warmed, RSS captured — neither the 20x nor the
60% figure should be quoted.

## Unverified

- **AVX2 and AVX-512 gather kernels have no benchmark numbers.** They compile
  (verified by x86-64 cross-compilation) and the AVX2 kernel is numerically
  correct (verified by forcing it under Rosetta, which emulates AVX2), but
  neither has been run on hardware that natively supports it, so no throughput
  claim is made. `SCRNA_FORCE_ISA=scalar|avx2|avx512` pins the kernel for
  like-for-like measurement once such hardware is available. Every number on
  this page is therefore from the **scalar** path; the SIMD kernels can only
  improve them.
- **HNSW recall above 50k cells**, for the reason given above. Validate on a
  subsample of your own data before trusting the graph at scale.
- **HNSW parameter tuning.** All runs use the defaults (M=16,
  ef_construction=200, ef_search=max(2k,64)). No sweep has been done, so the
  speed/recall trade-off here is a single point, not a frontier.
- **hnswlib's own AVX kernels.** By default they compile at baseline ISA (SSE on
  x86) because hnswlib gates them on compile-time `__AVX__`/`__AVX512F__` without
  target attributes; distances instead run through `ScrnaCosineSpace`, which is
  runtime-dispatched. `-DSCRNA_HNSW_NATIVE_ISA=ON` enables hnswlib's own kernels
  at the cost of a binary that requires those extensions. The two have not been
  compared.
- **Peak RSS is not yet tabulated.** `benchmark_1m_cells.py` reports it, but no
  systematic baseline comparison has been run.
- **Cache-miss counts** require Linux `perf`; not collected on the macOS box used
  for the numbers above.
