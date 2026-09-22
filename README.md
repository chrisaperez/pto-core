# pto-core

[![CI](https://github.com/chrisaperez/pto-core/actions/workflows/ci.yml/badge.svg)](https://github.com/chrisaperez/pto-core/actions/workflows/ci.yml)
[![Wheels](https://github.com/chrisaperez/pto-core/actions/workflows/wheels.yml/badge.svg)](https://github.com/chrisaperez/pto-core/actions/workflows/wheels.yml)
[![PyPI](https://img.shields.io/pypi/v/pto-core.svg)](https://pypi.org/project/pto-core/)
[![Python 3.9+](https://img.shields.io/badge/python-3.9%2B-blue.svg)](pyproject.toml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Five C++20 engines for single-cell and epigenomic data, built so the analysis
runs on the machine that already holds the sequencing data.

## Quickstart

```bash
pip install pto-core
```

```bash
pto-core doctor
```

`doctor` never touches your data. It reports which of the vendored engines
this wheel actually resolves, where each one came from (bundled inside the
wheel versus a `PATH` or `PTO_BIN_DIR` override), and the version each one
reports back — a wheel that can't clear this is broken, not slow, and this is
the one command to run before spending evaluator time on the rest.

```bash
pto-core qc your_reads.fastq.gz
```

That runs the streaming QC path (`fastq_stream --qc-only`) against a real
file and prints the summary — read counts, Q30 rate, mean quality — with no
build step, no CMake, and no cluster job. `pto-core man` lists every other
command the same way; `pto-core man frip` (or any other) prints one command's
parameters.

## `scrna_matrix` is a separate install

`pto-core` ships four of the five engines below — `fastq_stream`,
`genomic_toolkit`, `cuttag_profiler`, and `pto-peaks` — as one wheel with zero
Python runtime dependencies, and nothing left for the host to provide at the
system level either: htslib and libdeflate are already linked into the
binaries. The install this package exists to make frictionless is landing on
a login node or a shared workstation without touching whatever numpy or
pydantic version the site has already pinned, and without a compiler in
sight.

`scrna_matrix` is not in that wheel, and the reason is the same one that
makes its own CMake configure step a hard failure without OpenMP: a silently
single-threaded fallback that still reports success is worse than a build
that stops (see [OpenMP is mandatory](#openmp-is-mandatory-for-scrna_matrix),
below). Bundling that requirement into `pto-core` would mean a macOS host
without `libomp` could not install the wheel at all, just to get a FRiP
number out of `genomic_toolkit`. So it is published on its own, as a compiled
extension module rather than a CLI:

```bash
pip install pto-core scrna-matrix
```

Reach for both when the workflow needs a k-NN graph over an expression
matrix; `pto-core` by itself is enough for QC, fragment metrics, and CUT&Tag
profiling.

## Modules

The modules were developed independently and are vendored here so they can be
built, tested and reviewed together. Each remains a self-contained CMake project
with its own dependencies, options, tests and install rules, and each must stay
independently buildable with the other four deleted.

| Module | What it does | Language surface |
|:--|:--|:--|
| [`modules/scrna_matrix`](modules/scrna_matrix) | Cache-line-aligned sparse matrix engine (`Block_CSR<T>`), SIMD sparse-dot and cosine kernels, exact and HNSW k-NN | Header-only C++20 plus zero-copy `pybind11` bindings |
| [`modules/fastq_stream`](modules/fastq_stream) | Streaming FASTQ QC, adapter and quality trimming, piped straight into an aligner | C++20 library and CLI |
| [`modules/cuttag_profiler`](modules/cuttag_profiler) | CUT&Tag and CUT&RUN reference-point signal matrices, with an embedded loopback dashboard | Single C++20 binary |
| [`modules/genomic_toolkit`](modules/genomic_toolkit) | Fragment size distributions, duplicate rate and FRiP in one streaming pass, with no sort and no intermediate files | C++20 library and CLI |
| [`modules/peaks`](modules/peaks) | `pto-peaks`: single-pass pileup, multi-scale rolling background, SIMD Poisson tails, Benjamini-Hochberg, narrowPeak out | C++20 headers and CLI |

What unifies these is a deployment constraint rather than a shared abstraction:
each tool has to run on the cluster or workstation where the data already
lives, with no service dependency and no data leaving the host. There is no
shared runtime between the modules, and none is intended. See
[Repository layout](#repository-layout) for why that absence is a decision
rather than an omission.

There is no repository-wide document explaining how any one engine works
internally, in keeping with the "no shared layer" decision above — that lives
in each module's own README and, where one exists, its `docs/ARCHITECTURE.md`
(`fastq_stream` has one). [Architecture & System
Design](#architecture--system-design), further down, covers the other kind of
architecture: how the five are built, gated and packaged together.

## Demo

One command runs three of the five modules end to end and serves the dashboard:

```bash
./scripts/run_demo.sh
```

It is offline by default, with fixtures checked in or synthesised locally. See
[`docs/DEMO.md`](docs/DEMO.md) for options, deployment and troubleshooting.

The deployment kit reaches the dashboard over an SSH tunnel rather than a public
port. `cuttag_profiler` requires a per-run session token, but that token is one
shared secret for the run and not a login system; see
[`deploy/README.md`](deploy/README.md) for the actual security model. An
internet-facing port 8080 is unsupported.

## Tests

`ctest` at the top level runs thirty-two suites, or thirty-three when
`pybind11` is available and the Python binding tests register.

| Test | Module | Covers |
|:--|:--|:--|
| `test_matrix_ops` | scrna_matrix | `Block_CSR` invariants, SIMD parity against a double reference across every vector-tail boundary, exact k-NN, HNSW recall |
| `test_wide_accumulation` | scrna_matrix | scans the module's own sources and fails the build on a bare `float` reduction accumulator |
| `test_hnsw_persist` | scrna_matrix | HNSW serialization is lossless: a reloaded index returns byte-identical neighbour lists, not merely approximately equal ones |
| `test_python_bindings`\* | scrna_matrix | buffer adoption, lifetime and keepalive, GIL behaviour, NumPy and SciPy round-trips |
| `trimmer` | fastq_stream | SIMD kernels, lock-free ring, quality and adapter trimming, record framing, end-to-end pipeline over raw, gzip and BGZF |
| `fastq_stream_cli` | fastq_stream | the binary as a black box: argument validation, report files, non-regular-file inputs |
| `bam_query` | cuttag_profiler | BAM iteration, filtering, index handling |
| `signal_calc` | cuttag_profiler | window construction, normalisation, strand orientation, thread-count invariance, BED and GTF parsing |
| `matrix_cache` | cuttag_profiler | the cached-matrix memory bound: a cache capped at a fixed entry count is not a memory bound when nothing caps a region count |
| `interval_index` | cuttag_profiler | the pruning descent, differentially checked against an O(n) scan over 20k random queries |
| `secrets` | cuttag_profiler | session-token minting from the CSPRNG and constant-time comparison |
| `region_names` | cuttag_profiler | region-name handling from the HTTP security audit |
| `bed_reader` | cuttag_profiler | bounded line growth and other parsing behaviour |
| `file_identity` | cuttag_profiler | an output path may not name an input |
| `profiler_cli` | cuttag_profiler | the binary as a black box |
| `server_linkage` | cuttag_profiler | asserts against the linked artefact, in both directions, whether the embedded HTTP server is present |
| `http_security` | cuttag_profiler | `Host` and `Origin` validation: loopback names, IPv6 literals, port mismatch, and the substring tricks a `find`-based check would accept |
| `gtk_peak_set` | genomic_toolkit | interval merging, half-open boundaries, coordinate ceiling, differential test against a per-base bitmap |
| `gtk_dup_marker` | genomic_toolkit | sorted-stream and hash duplicate modes agreeing on the same set, the O(pile-up) memory claim, out-of-order and contig-revisit detection |
| `gtk_fragment_stream` | genomic_toolkit | BEDPE pairing, size and MAPQ filters, malformed-line accounting, CRLF, over-long lines, batching boundaries |
| `gtk_frip` | genomic_toolkit | end-to-end FRiP, the overlapping-peak double-count trap, duplicate exclusion, chromosome-naming mismatch |
| `gtk_fuzz_parsers` | genomic_toolkit | seeded mutation fuzzing of the BEDPE and fragment-BED parsers, plus emitted-fragment invariants |
| `gtk_file_identity` | genomic_toolkit | `--histogram` may not name an input, its index, or the peak file |
| `gtk_cli_args` | genomic_toolkit | the checked integer parse behind every numeric option: partial numbers, int32 narrowing, int64 overflow, per-option ranges |
| `gtk_cli` | genomic_toolkit | the binary as a black box: exit codes, refusals, truncated-BAM cases |
| `gtk_region_partition` | genomic_toolkit | htslib only: `--region` tiles partition a file exactly |
| `gtk_bam_contig_order` | genomic_toolkit | htslib only: FRiP on a multi-contig BAM does not depend on the peak file's contig order |
| `peaks_poisson_model` | peaks | Poisson tails against mpmath golden values, the degenerate-background fail-closed rule, the finite-lower-sum fast-path gate, convergence at the iteration ceiling |
| `peaks_sliding_window` | peaks | the centred multi-scale background rings, contig-end truncation, the coverage-desert skip, and back pressure |
| `peaks_caller` | peaks | the state machine, the float-comparison summit-plateau fix, Benjamini-Hochberg's cutoff-relative floor, batch-size invariance |
| `peaks_bam_streamer` | peaks | native BGZF and BAM decoding by byte shift, per-block CRC32, an empty BGZF block mid-stream staying transparent |
| `peaks_file_identity` | peaks | an output path may not name an input, including chrom-sizes and index files |
| `peaks_cli` | peaks | the binary as a black box: argument parsing, the BED reader's fail-closed checks, the narrowPeak writer, and the exit codes the job contract depends on |

\* Registered only when CMake finds `pybind11`, and skipped with a
configure-time message otherwise.

The C++ suites are dependency-free by design, with no GoogleTest and no Catch2,
for the same reason the tools are: they must run on a cluster that has a
compiler and nothing else.

### Sanitizers

Sanitizer builds are per-module rather than wired into the default build. For
`fastq_stream`:

```bash
cmake -S modules/fastq_stream -B build/fq-asan -DFQ_SANITIZE=ON -DFQ_NATIVE=OFF
cmake --build build/fq-asan -j && ./build/fq-asan/test_trimmer
```

No macOS toolchain combines ThreadSanitizer with OpenMP: Apple Clang has no
OpenMP, and Homebrew GCC has no sanitizer runtimes on arm64. Use the container,
which pins clang, `libomp-dev` and libarcher:

```bash
docker build -f docker/Dockerfile.tsan -t pto-tsan . && docker run --rm pto-tsan
```

`scripts/run_tsan.sh` sweeps thread counts and repeats, and refuses to run
rather than "verifying" a serial build when the toolchain cannot supply
`-fopenmp -fsanitize=thread`. `fastq_stream` and `cuttag_profiler` are
`std::thread`-based and need no container.

Two sanitizer notes worth knowing before you reproduce an audit. A UBSan
suppression is inert under `-fno-sanitize-recover=all`, so use
`-fno-sanitize=alignment` for hnswlib's known misaligned store. And
`AlignedAllocator` rounds allocations up to 64 bytes, so an ASan overread
reproduction needs `posix_memalign` at the exact size.

### Public-data validation

The suites above prove each module does what its author intended on fixtures
its author wrote. [`tests/validation/`](tests/validation) answers the other
question: on real public sequencing runs, does `pto-core` produce the same
answers as the tools it replaces?

Those are not the same question. A synthetic FASTQ has no adapter dimers, no
quality cliff at cycle 140 and no duplicate structure. A synthetic BAM has no
soft-clipping and no unplaced contigs. Every interesting disagreement between
two tools lives in exactly that material.

```bash
tests/validation/run_validation.py list                    # what can run here
tests/validation/run_validation.py fetch --allow-download  # stage inputs
tests/validation/run_validation.py run
```

Eleven cases across four modules: TSS meta-profiles against deeptools
(Pearson r at least 0.99), FRiP against samtools and bedtools (exact count),
fragment sizes against samtools TLEN (KS at most 0.02), QC counters against
seqkit (exact), trimming against fastp, alignment rate through bowtie2,
duplicate rate against Picard, and HNSW recall against `scrna_matrix`'s own
exact path. Datasets are real pinned accessions (ENCODE `ENCSR000AKA`, ENA
`PRJNA512492`) verified by checksum.

Two properties are load-bearing rather than stylistic. Network access is
opt-in: without `--allow-download` nothing touches the internet, and a missing
input is a clean skip carrying the command to stage it. These tools run where
the data is, and that host is the last place a test harness should be reaching
out from. And skips are not failures, but a case with no metrics is an error.
Silently passing because every comparison was skipped is the most dangerous
thing a validation suite can do, and it is the default if you only check for
failures.

Where exact equality is not expected, the reason is documented next to the
threshold rather than absorbed into a comfortable tolerance, and both observed
values are printed on every run so that the gap changing is visible.
`markdup` against Picard compares aligned span to unclipped 5' end;
`cuttag_profiler` against deeptools compares span-overlap to midpoint binning.

## Performance

Measured on an Apple M4 Pro with 12 cores and Homebrew GCC 16. See
[`docs/PERFORMANCE_2026-08-15.md`](docs/PERFORMANCE_2026-08-15.md) for the full
report, the exact commands, and what was not measured. The cross-tool Pareto
table, which is the number to read if you are deciding whether to adopt any of
this, is the one directly below.

| Module | Headline | Against the baseline |
|:--|:--|:--|
| `fastq_stream` | 16.6 M reads/s, 4,954 MB/s at 12 threads; constant 10 MB peak RSS | Peak RSS does not grow with input size at all |
| `scrna_matrix` | 1.3M cells in 37.1 s at 1.01 GB peak RSS (HNSW); NEON gives 2.02x on brute-force k-NN | See the caveat below |
| `cuttag_profiler` | 100k windows in 1.59 s at 88.9 MB | 32.6x faster and 4.7x less memory than deeptools, at r = 0.999991 |
| `genomic_toolkit`\* | 400k fragments scored in 0.05 s, no intermediate files | 28x faster than the samtools-sort, bamtobed, sort, intersect chain, with an identical in-peak count of 148,680 |

\* `genomic_toolkit` is not part of the 2026-08-15 performance run. It postdates
it, and its row was measured separately with Apple Clang on the repository's own
`demo.bam`. Reproduce both sides, counts included, with
`modules/genomic_toolkit/scripts/validate_against_bedtools.sh`.

**The `scrna_matrix` row is feasibility, not superiority.** The 1.3M-cell figure
is a real run including index construction and all 1.3M queries, and the exact
path extrapolates to roughly 1.6 hours on the same data at about 3.5 ns per
n-squared, which is a projected 155x. But there is no like-for-like
`scanpy.pp.neighbors` run at 1.3M cells, so no comparative claim is made at that
scale. A 50,000-cell head-to-head does exist in the performance report, and even
that carries two caveats stated there rather than buried: numba JIT warm-up
dominates a single first call, and a graph built from `X` against scanpy's
default compares two different computations. `sc.pp.neighbors` on a PCA
representation is itself approximate and near-linear, which makes it a peer of
the HNSW backend rather than a baseline it should categorically beat.

**`peaks` has no row here at all.** It has not been benchmarked for speed
against MACS3, only for concordance, and on that axis this project's own 0.98
parity gates are not met. Summit location agrees to a 1 bp median offset with
97.6% of matched summits inside 50 bp, and peak strength correlates at Pearson
r 0.990, but `pto-peaks` calls intervals roughly 3x narrower, so mean IoU is
0.458. Full numbers and the open options for closing the gap:
[`modules/peaks/CONCORDANCE.md`](modules/peaks/CONCORDANCE.md).

### Reproducing the numbers yourself

The table above is a record of one run on one machine.
[`scripts/benchmark_audit.py`](scripts/benchmark_audit.py) sweeps a
configuration matrix of input-size tiers by thread counts, and with
`--baselines` runs the standard pipeline each tool replaces on the same input at
the same core count in the same session:

```bash
scripts/benchmark_audit.py --list
scripts/benchmark_audit.py --inputs small --threads 1 2 4 8 --baselines --repeats 3
```

Four choices in it are deliberately unflattering, because a benchmark that only
reports favourable numbers is marketing:

* A speedup is printed only when both sides ran in the same session. Comparing a
  fresh run against a number quoted from a README is how benchmarks start lying.
* Peak RSS comes from the kernel's `ru_maxrss` via `wait4` and never from a
  sampling thread. A peak is a spike, and a sampler is asleep during it.
* For a multi-stage baseline, wall time sums across stages but peak RSS takes
  the maximum, since the stages do not overlap. Summing would overstate the
  baseline and flatter `pto-core`.
* Parallel efficiency above 1.0 is flagged rather than celebrated. It almost
  always means the 1-thread baseline was contended, and here specifically it
  means htslib decompressed BGZF on its own threads: a 1-thread baseline on the
  `cuttag/profile` case shows 4.32 cores busy, so a figure above 100% is a
  property of the baseline and not superlinear scaling.

Skipped configurations are listed explicitly, with the note that any speedup
claim covers only the rows that actually ran.

One finding still worth acting on: plain gzip caps `fastq_stream` at roughly
2.2 M reads/s regardless of thread count, and BGZF is 3.2x faster at 12 threads.
That is the framing result from
[`fastq_stream`'s architecture doc](modules/fastq_stream/docs/ARCHITECTURE.md)
observed from the other end.

## Security and review trail

Findings are logged before they are fixed, and residuals are recorded rather
than quietly dropped. Across seven adversarial passes between 2026-08-15 and
2026-09-11 these documents record **74 confirmed defects, 24 of them rated High
or Critical**, every one carrying a regression test added in the commit that
fixed it.

These are living documents. Extend them rather than replacing them when you
find or fix something in their scope.

| Document | Scope | Outcome |
|:--|:--|:--|
| [`REVIEW_2026-08-15.md`](docs/REVIEW_2026-08-15.md) | cross-module architecture and security | 5 findings, all closed with regression tests |
| [`SECURITY_HTTP_2026-08-15.md`](docs/SECURITY_HTTP_2026-08-15.md) | `cuttag_profiler`'s embedded HTTP server and path confinement | all 9 closed, one stated residual |
| [`TORTURE_2026-08-17.md`](docs/TORTURE_2026-08-17.md) | adversarial fuzzing, sanitizer and concurrency sweep | 5 defects fixed, 1 false positive left unsuppressed |
| [`TORTURE_2026-09-10.md`](docs/TORTURE_2026-09-10.md) | `fastq_stream`, then `pto-peaks` | 14 and 8 findings respectively |
| [`AUDIT_2026-09-11_fastq_stream.md`](docs/AUDIT_2026-09-11_fastq_stream.md) | Phred offsets, record framing, per-cycle truncation | 8 findings, all closed |
| [`AUDIT_2026-09-11_genomic_toolkit.md`](docs/AUDIT_2026-09-11_genomic_toolkit.md) | `sizes`, `markdup`, `frip` | 9 findings, all closed |
| [`AUDIT_2026-09-11_cuttag_profiler.md`](docs/AUDIT_2026-09-11_cuttag_profiler.md) | `profile`, bin geometry, skipped regions | 6 findings, all closed |
| [`AUDIT_2026-09-11_scrna_matrix.md`](docs/AUDIT_2026-09-11_scrna_matrix.md) | HNSW producers, index concurrency, NumPy boundary | 11 findings, all closed |

The headline findings are worth naming, because in every case the program ran to
completion and printed a plausible number:

* A BAM cut at a BGZF block boundary reported half a library at exit 0.
* Regions on a contig absent from the BAM were averaged into the meta-profile as
  zeros, diluting the curve in direct proportion to how many were skipped.
* An unchecked `n * k` in every HNSW producer stored through a null pointer.
* Two concurrent index mutations left an index holding half its rows, with no
  exception raised and a perfectly valid file on disk.
* An iteration cap silently scored bases 4 to 19 orders of magnitude past the
  significance cutoff as p = 1 for counts above roughly 1e9.
* `-i X -o X` truncated the input.
* Phred+64 quality strings were read as Phred+33, moving mean Q from 19.85 to
  50.85 and the Q30 rate from 0.265 to 1.000.

**The recurrence is the methodological finding.** Four of the review's five
findings were a previously-closed defect surviving in a sibling function or a
second entry point the earlier pass did not walk. In `scrna_matrix`, the same
unchecked narrowing and the same unchecked `n * k` were each found twice. When
you fix one of these, grep for its shape across the tree before calling it done.
Several of the guards here are source-level scans precisely because a rule a
human has to remember is a rule that gets broken in the next sibling function.

`pto-peaks` has no separate audit document of its own. Its 8 hardening findings
are the second half of the 2026-09-10 pass, and its MACS3 concordance status is
a different question tracked in
[`modules/peaks/CONCORDANCE.md`](modules/peaks/CONCORDANCE.md) rather than as a
security finding. `modules/scrna_matrix/docs/AUDIT.md` is the earlier
module-level equivalent and the model the others follow: reproduction before
fix, honest residuals, and a section for what turned out not to be ours.

## Architecture & System Design

### Building from source

The top-level `CMakeLists.txt` is an aggregate. It adds each module as a
subdirectory and registers every module's tests with one CTest project, so the
whole tree builds and tests in a single pass.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```

Requires CMake 3.20 or newer, a C++20 compiler, and a working OpenMP toolchain
for `scrna_matrix`.

#### OpenMP is mandatory for `scrna_matrix`

Configuration fails without it, and this is deliberate. It used to be optional
with a `message(WARNING)` fallback, and the result was that `_OPENMP` went
undefined, every `#pragma omp` compiled to serial code, and the test suite
reported 100% passed while never executing a parallel iteration. A silently
single-threaded engine that CI signs off on is worse than a build that stops.
It is also the reason `scrna_matrix` ships as [its own
wheel](#scrna_matrix-is-a-separate-install) rather than inside `pto-core`.

| Platform | What to do |
|:--|:--|
| macOS, Apple Clang | `brew install libomp`, then plain `cmake -S . -B build` works. The keg is auto-detected. |
| macOS, Homebrew GCC | `brew install gcc`, then `-DCMAKE_C_COMPILER=gcc-15 -DCMAKE_CXX_COMPILER=g++-15`. |
| Debian/Ubuntu, Clang | `sudo apt-get install libomp-dev` |
| Debian/Ubuntu, GCC | Already present as libgomp. |

Apple Clang needs `-Xclang -fopenmp` and Homebrew installs `libomp` keg-only,
so CMake's `FindOpenMP` reports NOT FOUND on a stock Mac even after
`brew install libomp`. [`cmake/PtoOpenMP.cmake`](cmake/PtoOpenMP.cmake) closes
that gap: it locates the keg and hands `FindOpenMP` the flags it cannot derive.
For a non-standard Homebrew prefix:

```bash
cmake -S . -B build -DLIBOMP_ROOT=$(brew --prefix libomp)
```

That makes the requirement satisfiable. It does not relax it. Detection is
followed by a `check_cxx_source_compiles` probe that fails the build unless
`_OPENMP` is actually defined and a `#pragma omp parallel` region compiles,
because `find_package` succeeding proves only that the compiler accepts the
flag.

`fastq_stream` and `cuttag_profiler` contain no OpenMP code and build without
any of this, under `-DPTO_BUILD_SCRNA_MATRIX=OFF`.

#### Module switches

A host missing one module's dependencies can still build and test the rest.

| Option | Default | Effect |
|:--|:--|:--|
| `PTO_BUILD_SCRNA_MATRIX` | ON | build `modules/scrna_matrix` |
| `PTO_BUILD_FASTQ_STREAM` | ON | build `modules/fastq_stream` |
| `PTO_BUILD_CUTTAG_PROFILER` | ON | build `modules/cuttag_profiler` |
| `PTO_BUILD_GENOMIC_TOOLKIT` | ON | build `modules/genomic_toolkit` |
| `PTO_BUILD_PEAKS` | ON | build `modules/peaks` |

```bash
cmake -S . -B build -DPTO_BUILD_CUTTAG_PROFILER=OFF
```

Disabling every module is a configure-time error rather than an empty build,
and CI asserts that it stays one.

#### Dependencies

| Module | Required | Optional |
|:--|:--|:--|
| `scrna_matrix` | none (header-only) | `pybind11` for Python bindings, OpenMP for parallel k-NN; `hnswlib` is vendored |
| `fastq_stream` | `zlib`; `libdeflate` (found if installed, otherwise fetched and built statically) | none |
| `cuttag_profiler` | `htslib` via `pkg-config` | none |
| `genomic_toolkit` | none | `htslib` for BAM/CRAM/`.gz` input (text paths and the whole test suite work without it), OpenMP for the FRiP reduction |
| `peaks` | none | `libdeflate`, else `zlib`, for BAM/BGZF input. Discovered, never fetched, via `PEAKS_WITH_BGZF=AUTO`/`ON`/`OFF`. Takes no htslib dependency at all. |

On macOS with Homebrew, add `-DCMAKE_PREFIX_PATH=/opt/homebrew`. `httplib` and
`nlohmann/json` are vendored under `modules/cuttag_profiler/third_party/`.

Without `pybind11` the `scrna_matrix` Python module and its binding tests are
skipped with a configure-time message, and the C++ suite still runs. To include
them, point CMake at an interpreter that has it:

```bash
cmake -S . -B build -Dpybind11_DIR=$(python -c 'import pybind11; print(pybind11.get_cmake_dir())')
```

#### Building a single module

Nothing at the top level is required to use a module, and a module's build may
not depend on being invoked from there:

```bash
cmake -S modules/fastq_stream -B build/fastq_stream -DCMAKE_BUILD_TYPE=Release
cmake --build build/fastq_stream -j
```

Module-specific options (`FQ_*`, `SCRNA_*`, `CUTTAG_*`, `GTK_*`, `PEAKS_*`) are
documented in each module's README and work identically in both modes.

### Packaging and distribution

The wheel-facing decisions live in [`pyproject.toml`](pyproject.toml) and
[`modules/scrna_matrix/pyproject.toml`](modules/scrna_matrix/pyproject.toml),
each with the reasoning kept inline rather than assumed. The parts worth
knowing without reading both files end to end:

* **Two wheels, not one.** `pto-core` configures with
  `PTO_BUILD_SCRNA_MATRIX=OFF` — the split explained
  [above](#scrna_matrix-is-a-separate-install) — and its own distribution is
  ABI-neutral (`py3-none-<platform>`), because it contains three ordinary
  executables and a pure-Python wrapper, nothing bound to an interpreter's
  ABI. `scrna-matrix` is the opposite: a compiled extension module, tagged
  normally per interpreter.
* **htslib is statically linked on Linux, dynamically on macOS, and the
  difference is not a style choice.** manylinux has no htslib package, so
  `packaging/scripts/build_htslib_manylinux.sh` builds it from source and
  links it in statically; `ci.yml`'s `static-link` job asserts the resulting
  binary carries no `libhts`-family shared-object dependency. The same static
  link succeeds against Homebrew htslib on macOS too, because its remaining
  need — libcurl — is the system dylib at `/usr/lib/libcurl.4.dylib`; that is
  why `cuttag_profiler` links htslib statically there while `fastq_stream` and
  `genomic_toolkit` still link it dynamically, with `delocate-wheel` vendoring
  `libhts`, `libomp`, `libdeflate` and `liblzma` into `pto/.dylibs/` and
  rewriting the two binaries' install names.
* **No `universal2` macOS wheel.** Homebrew's htslib and libomp are
  per-architecture; an arm64 runner cannot link an x86_64 `libhts.dylib`, and
  the cross build either fails at link time or produces a wheel `dyld`
  rejects on the user's machine. Both architectures ship as native wheels
  instead, one runner per architecture (macos-13 for x86_64, macos-14 for
  arm64).
* **Every native-ISA flag defaults OFF for a wheel build**
  (`PROFILER_NATIVE_ARCH`, `FQ_NATIVE`, `GTK_NATIVE_ARCH`,
  `PEAKS_ENABLE_AVX512`). These flags put the *build host's* instruction set
  onto the binary; a wheel installs on machines it never saw, and an AVX-512
  instruction reaching a CPU without one is a `SIGILL` with no diagnostic the
  user can act on — and not even at import, but partway through someone's
  batch.
* **Third-party code is vendored, not fetched at install time.**
  `cpp-httplib` and `nlohmann/json` live under
  `modules/cuttag_profiler/third_party/`; `hnswlib`, with the local patches
  recorded in `third_party/hnswlib/VERSION.txt`, lives under
  `modules/scrna_matrix/third_party/`. Full license text for all three is in
  [`THIRD_PARTY_LICENSES.md`](THIRD_PARTY_LICENSES.md).

## Repository layout

```
pto-core/
├── CMakeLists.txt          # aggregate build; owns no sources of its own
├── pyproject.toml          # the `pto-core` wheel; see packaging/
├── cmake/                  # shared CMake modules (OpenMP toolchain resolution)
├── deploy/                 # container entrypoint and demo provisioning
├── docker/                 # TSan verification image and suppressions
├── docs/                   # repository-wide review and audit documents
├── packaging/              # the wheel's Python source (`pto`) and its own suite
├── scripts/                # verification harnesses
├── tests/validation/       # public-data differential suite; drives the CLIs, links nothing
└── modules/
    ├── scrna_matrix/       # header-only engine and pybind11 bindings
    ├── fastq_stream/       # streaming FASTQ QC and trimming
    ├── cuttag_profiler/    # CUT&Tag profiling with embedded dashboard
    ├── genomic_toolkit/    # streaming fragment QC
    └── peaks/              # pto-peaks: streaming peak caller
```

CI is not in this directory. GitHub reads workflows only from the repository
root, one level up, and every `run` step there executes with
`working-directory: pto-core`. That is load-bearing: the monorepo restructure
moved `CMakeLists.txt` down a level while the workflow still ran `cmake -S .`
from the repository root, and all seven jobs failed with an error that reads
like a build failure and was a path failure.

Every directory here holds files. Earlier revisions carried empty `include/`,
`python/`, `tests/` and `third_party/` directories at the root, anticipating a
shared top-level layer: a common headers namespace, a unified Python package,
cross-module integration tests. That layer does not exist and is not planned.
The modules share no code and have disjoint dependency sets, so introducing one
would be a design decision to take deliberately rather than a gap to fill in.

`tests/validation/` is not a reinstatement of that layer, and the distinction is
worth being precise about. It contains no C++, links against no module, and no
module's build or test suite references it; nothing under `modules/` can break
because of it. It is a Python harness that runs the shipped binaries as black
boxes against public data. The rule it must keep obeying: a module must remain
independently buildable and testable with this directory deleted.

## License

MIT. [`LICENSE`](LICENSE) governs the whole of `pto-core`: the aggregate build
and all five modules, including `fastq_stream`, which also carries its own
identical copy from before that file existed. First-party sources carry
`SPDX-License-Identifier: MIT`; that header is a per-file marker rather than a
substitute for the root grant, and no module needs a license file of its own for
the grant to apply.

Vendored third-party code keeps its own license, reproduced in full in
[`THIRD_PARTY_LICENSES.md`](THIRD_PARTY_LICENSES.md): `hnswlib` (Apache-2.0),
`cpp-httplib` (MIT), `nlohmann/json` (MIT).
