# pto-core — performance report

**Date:** 2026-08-15
**Host:** Apple M4 Pro, 12 cores, 24 GB, macOS 26.5.2
**Toolchain:** Homebrew GCC 16.1.0, CMake 4.3.4, Ninja
**Flags:** `-O3 -DNDEBUG -std=gnu++20 -arch arm64 -fopenmp` (`-mcpu=native` via `FQ_NATIVE=ON` / `PROFILER_NATIVE_ARCH=ON`)
**Baselines:** deeptools 3.5.6, scanpy 1.11.5 / anndata 0.11.4 / numpy 2.2.6 / scipy 1.15.3, samtools 1.22

Every number below was measured on this host in this session. Correctness was
verified first — `ctest` is 6/6 on the same binaries — because a fast wrong
answer is not a benchmark result.

---

## Reproducing

```bash
# 1. Build (Release, GCC 16, native tuning)
cmake -S . -B build-bench -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER=/opt/homebrew/bin/gcc-16 \
  -DCMAKE_CXX_COMPILER=/opt/homebrew/bin/g++-16 \
  -DFQ_NATIVE=ON -DPROFILER_NATIVE_ARCH=ON -DSCRNA_HNSW_NATIVE_ISA=ON \
  -Dpybind11_DIR="$(python3 -c 'import pybind11; print(pybind11.get_cmake_dir())')"
cmake --build build-bench -j12
ctest --test-dir build-bench --output-on-failure     # verify before timing

# 2. fastq_stream — component + end-to-end microbenchmarks
./build-bench/modules/fastq_stream/benchmark

# 3. scrna_matrix — scaling
cd modules/scrna_matrix
export PYTHONPATH=../../build-bench/modules/scrna_matrix
python tests/benchmark_1m_cells.py --make-synthetic 200000 --output synth.h5ad
python tests/benchmark_1m_cells.py --input synth.h5ad --backend scrna --method hnsw --n-pcs 50 --k 15
python tests/benchmark_1m_cells.py --input synth.h5ad --backend scanpy       --n-pcs 50 --k 15

# 4. cuttag_profiler — head-to-head against deeptools
cd modules/cuttag_profiler
pip install deeptools
scripts/benchmark.sh \
  --bam data/bench.bam --regions data/bench_promoters.bed \
  --profiler ../../build-bench/modules/cuttag_profiler/cuttag_profiler \
  --threads 12 --replicates 3
```

---

## 1. fastq_stream

### Component kernels (best of 5)

| Kernel | Throughput | Run-to-run spread |
|---|---:|---:|
| `sum_u8` (mean Phred) | 82.30 GB/s | 21% |
| `trim_sliding_window` | 80.47 GB/s | 10% |
| `count_ge` (Q30 fraction) | 59.15 GB/s | 19% |
| `add_widen` (per-cycle quality) | 13.87 GB/s | 15% |
| adapter seed+extend, 2 adapters | **1.63 GB/s** | 1% |

SIMD backend: **NEON**. Adapter matching is ~50× slower than every other
kernel and is the pipeline's compute bottleneck — `--no-adapter` is
substantially faster when you do not need it.

Lock-free SPSC ring: **54.36 M items/s**, mean hand-off latency **18.39 ns**.

### End-to-end, 4M reads × 150 bp

Measured through the CLI (`--qc-only`), best of 3 per configuration.

| Input | Threads | Reads/s | MB/s | Peak RSS |
|---|---:|---:|---:|---:|
| uncompressed (1197 MiB) | 4 | 2,514,580 | 752 | 17 MB |
| uncompressed | 8 | 11,697,638 | 3,500 | 29 MB |
| uncompressed | 12 | **16,556,711** | **4,954** | 42 MB |
| gzip (583 MiB) | 4 | 2,246,549 | 327 | 8 MB |
| gzip | 8 | 2,204,618 | 321 | 9 MB |
| gzip | 12 | 2,233,370 | 325 | 8 MB |
| BGZF (581 MiB) | 4 | 2,388,179 | 347 | 8 MB |
| BGZF | 8 | 2,463,038 | 358 | 9 MB |
| BGZF | 12 | **7,139,031** | **1,037** | 9 MB |

MB/s is measured against the **on-disk** size, so compressed rows are not
comparable to the uncompressed row — compare reads/s across rows instead.

> ### Finding: plain gzip does not scale with threads
>
> gzip input is flat at ~2.2 M reads/s from 4 to 12 threads — adding 3× the
> threads buys **1%**. BGZF on the same data reaches 7.1 M reads/s at 12
> threads, **3.2× faster**.
>
> This is structural, not a bug: a single-member gzip stream cannot be
> block-decoded, so `reader.hpp` inflates it serially through zlib and the
> decompressor becomes the bottleneck before the worker pool is saturated.
> BGZF is a series of independent blocks, so it parallelises.
>
> **Recommendation:** store FASTQ as BGZF (`bgzip`), not plain gzip. If inputs
> arrive as gzip and are read more than once, a one-off `zcat | bgzip` pays for
> itself immediately. At 12 threads the uncompressed path is faster still
> (16.6 M reads/s) when the I/O budget allows it.

Memory is flat and small throughout — 42 MB peak at the highest throughput
configuration, and under 10 MB on every compressed run. The pool ceiling is
structural and does not grow with input size.

---

## 2. scrna_matrix

Synthetic datasets, 2000 genes, ~100 nnz/cell, k=15.

### Exact brute force — O(n²), as designed

| Cells | Wall clock | Peak RSS | Scaling vs previous |
|---:|---:|---:|---|
| 5,000 | 0.146 s | 0.175 GB | — |
| 10,000 | 0.560 s | 0.166 GB | 3.8× for 2× cells |
| 20,000 | 2.321 s | 0.191 GB | 4.1× for 2× cells |

Quadratic, as documented. This is a correctness baseline, not a path to 1M
cells: extrapolating, 1M cells would be ~1.6 hours of similarity evaluation.

### Approximate HNSW — the path that scales

| Cells | PCs | Wall clock | Peak RSS | Est. index |
|---:|---:|---:|---:|---:|
| 50,000 | 50 | 1.195 s | 0.240 GB | 0.02 GB |
| 200,000 | 50 | 6.383 s | 0.615 GB | 0.08 GB |

4× the cells for 5.3× the time — near-linear, and memory grows with the index
rather than with n².

### Head-to-head vs scanpy, 50,000 cells, identical input

| Backend | Wall clock | Peak RSS |
|---|---:|---:|
| `scanpy.pp.neighbors` | 21.092 s | 0.881 GB |
| `scrna_matrix` (HNSW, 50 PCs) | **1.195 s** | **0.240 GB** |
| **Advantage** | **17.6× faster** | **3.7× less memory** |

> ### RESOLVED: NEON kernels added
>
> The numbers above were measured before `simd_math.hpp` had any ARM64 path —
> every run reported `simd=scalar`, because the AVX2/AVX-512 kernels are gated
> on `SCRNA_BUILD_X86_SIMD`. NEON `gather_dot` and `dense_dot` kernels have
> since been added; see [NEON results](#4-neon-arm64-kernels) below.
>
> Treat the `scrna_matrix` figures in this section as the **scalar baseline**.

---

## 3. cuttag_profiler

100,000 promoter windows, ±2 kb, 50 bp bins, 87 MB BAM, 12 threads, 3 replicates.

| Stage | Mean wall clock | Peak RSS |
|---|---:|---:|
| **`cuttag_profiler`** (BAM → matrix) | **1.59 s** | **88.9 MB** |
| `deeptools bamCoverage` (BAM → bigWig) | 7.49 s | 334.6 MB |
| `deeptools computeMatrix` (bigWig → matrix) | 44.21 s | 421.5 MB |
| deeptools pipeline total | 51.70 s | 421.5 MB |

| Metric | Value |
|---|---:|
| Speed-up, end to end | **32.6×** |
| Speed-up vs `computeMatrix` alone | 27.9× |
| Peak memory reduction | **4.7×** |

Per-replicate spread was 1.45–1.82 s for the profiler and 43.73–44.81 s for
`computeMatrix`, so the ratio is not an artefact of one lucky run.

### Accuracy against the baseline

Speed means nothing if the matrices disagree, so the harness compares them
bin-for-bin:

| Metric | Value |
|---|---:|
| Meta-profile Pearson r | **0.999991** |
| Max residual (% of peak) | 0.384% |
| Per-region r, median (n=2000) | 0.988098 |
| Per-region r, 10th percentile | 0.977268 |
| Fraction of regions with r > 0.95 | **99.80%** |

`PASS: meta-profile correlation >= 0.99`.

A fitted global scale factor of **1.997** is applied before comparison. That
near-exact 2× is a normalisation-convention difference between the two tools,
not disagreement about shape — the correlation after scaling is 0.999991. Worth
knowing if you compare absolute values across tools; irrelevant if you compare
profiles.

---

## 4. NEON (ARM64) kernels

Added after the measurements above. `gather_dot` and `dense_dot` now have
ARM64 NEON implementations; the dispatch enum gained `kNeon`, and
`SCRNA_FORCE_ISA=neon|scalar` selects between them at run time on one host.

`sparse_dot` was deliberately **not** vectorised — see [Why not
`sparse_dot`](#why-not-sparse_dot).

### Isolated kernels (M4 Pro, GCC 16 `-O3 -mcpu=native`, best of 7)

| Kernel | Scalar | NEON | Speed-up |
|---|---:|---:|---:|
| `gather_dot` (2000 nnz) | 1414 ns/call · 1.41 G elem/s | 314 ns/call · **6.36 G elem/s** | **4.50×** |
| `dense_dot` (512 dims) | 174 ns/call · 2.94 G elem/s | 29 ns/call · **17.49 G elem/s** | **5.94×** |

### End-to-end

| Workload | Scalar | NEON | Speed-up |
|---|---:|---:|---:|
| Brute-force k-NN, 20k cells (`gather_dot` path) | 2.249 s | **1.116 s** | **2.02×** |
| HNSW, 50k cells (`dense_dot` path) | 1.155 s | **0.918 s** | 1.26× |

The gap between kernel and end-to-end speed-up is Amdahl, not a defect: the
scatter, the `clear`, the bounded top-K heap and the norm pass do not
vectorise, so they become the larger share of the remaining time.

### Why `gather_dot` gained 4.5× when NEON has no gather

NEON has no gather instruction — no `_mm256_i32gather_ps` equivalent, and
Apple Silicon does not expose SVE. The four "gathered" loads are still four
scalar loads assembled lane-wise with `vsetq_lane_f32`.

So 4.5× is **not** "the gather got vectorised". It says the scalar version was
**dependency-chain bound, not load bound**: the loads were already overlapping,
and what cost time was `acc += ...` serialising every iteration at FP-add
latency. Two independent `float32x4_t` accumulators break that chain and keep
several FMAs in flight. `dense_dot`, which has no gather at all, gains more
(5.94×) — consistent with that reading.

### Why not `sparse_dot`

1. **It is not in the float hot path.** `knn_graph.hpp:265` uses `gather_dot`
   for `T = float`; `sparse_dot` is reached only at `:274`, the non-float
   branch. Vectorising it would not have touched anything the benchmarks above
   execute.
2. **A merge-join is the wrong shape for NEON.** It is a data-dependent
   branch per element (`a == b` / `a < b`), inherently sequential. SIMD set
   intersection needs either conflict-detection instructions NEON lacks, or a
   different algorithm entirely.

### Accuracy

NEON is **more** accurate than scalar, not less — two accumulators plus a
horizontal reduction is a shallower summation tree than one sequential
accumulator. Accumulated absolute error against a double reference over 200
random trials:

| Path | Accumulated \|error\| |
|---|---:|
| scalar | 2.373e-04 |
| **NEON** | **1.048e-04** |

## Summary

| Module | Headline | vs baseline |
|---|---|---|
| `fastq_stream` | 16.6 M reads/s, 4,954 MB/s (12 threads, uncompressed) | not measured — `fastp`/`Trimmomatic` not installed |
| `scrna_matrix` | 200k cells in 6.4 s (HNSW, scalar); NEON adds **2.02×** on brute-force k-NN | **17.6× faster, 3.7× less memory** than scanpy |
| `cuttag_profiler` | 100k windows in 1.59 s, 88.9 MB | **32.6× faster, 4.7× less memory** than deeptools, r=0.999991 |

### Two actionable findings

1. **Plain gzip caps `fastq_stream` at ~2.2 M reads/s** regardless of thread
   count. BGZF gives 3.2× at 12 threads. Change the storage format, not the
   flags.
2. ~~`scrna_matrix` has no NEON kernel~~ — **done.** NEON `gather_dot` and
   `dense_dot` landed: 4.50× and 5.94× on the isolated kernels, 2.02× on
   end-to-end brute-force k-NN. Section 4 has the detail.

### What was not measured, and why

- **`fastq_stream` vs `fastp` / `Trimmomatic`.** Neither is installed on this
  host. `docs/PAPER_DRAFT.tex` Table 2 deliberately reads `TODO` rather than
  carrying estimates; that stays true. To fill it in:
  `brew install fastp`, then `scripts/run_benchmarks.py` on GIAB HG002.
- **The 1M-cell claim.** The largest run here is 200k cells. `docs/BENCHMARKS.md`
  is explicit that the 1M+ figure is not yet substantiated, and nothing here
  changes that. HNSW's near-linear scaling makes it plausible; plausible is not
  measured.
- **Real sequencing data for `fastq_stream`.** Reads are synthetic (adapter
  read-through, 3′ quality decay, N-rich reads). Kernel throughput is
  content-independent, but filter-pass rates on real data will differ.
- **Thermals.** Single-host laptop-class numbers with a hybrid CPU. Consecutive
  runs of an unchanged kernel varied up to 21% (see the spread column), which is
  why every component figure is a best-of-N. Treat these as ranking evidence,
  not absolutes.
