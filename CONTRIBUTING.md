# Contributing to pto-core

This repository is five independently-developed C++20 tools, vendored
together because they share a deployment constraint (must run on the machine
that already holds the sequencing data) rather than a codebase. See the
[README](README.md) — especially "Modules" and "Repository layout" — before
your first change: the "no shared layer" rule there is a deliberate,
previously-reverted decision, not an oversight, and it shapes what a PR here
is allowed to touch.

## Setting up a development environment

You need:

* **CMake 3.20+** and a **C++20 compiler** (Clang or GCC; see the README's
  [Building from source](README.md#building-from-source) for the exact
  per-platform OpenMP setup, which is the one dependency worth reading about
  before you configure anything).
* **A Python virtual environment** for anything under `packaging/`,
  `scripts/`, or `tests/validation/` — never install into system Python.

```bash
python3 -m venv .venv
source .venv/bin/activate          # .venv\Scripts\activate on Windows
pip install --upgrade pip
```

### Building the C++ engines

Either the aggregate build or a single module — both are required to keep
working, and a PR that only exercises one should test that module standalone:

```bash
# everything
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j

# one module, in isolation (must always work with the other four deleted)
cmake -S modules/fastq_stream -B build/fastq_stream -DCMAKE_BUILD_TYPE=Release
cmake --build build/fastq_stream -j
```

### Building an editable `pto-core` wheel

For working on `packaging/pto/` (the CLI and Python API) against a real
build of the engines, without reinstalling on every change:

```bash
pip install scikit-build-core[pyproject] ninja pybind11
pip install -e . --no-build-isolation
```

`--no-build-isolation` is required, not optional flourish: `pto-core`'s
[build backend is `scikit-build-core`](pyproject.toml), which runs CMake
under the hood, and PEP 517's default isolated build environment would
re-download and reconfigure the whole C++ tree into a throwaway venv on every
`pip install`. Skipping isolation means pip reuses your active environment
and, more importantly, your existing `build/` directory — an edit under
`packaging/pto/` reinstalls in seconds; an edit under `modules/` still costs
a real (if incremental) CMake rebuild, because scikit-build-core does not
know your source tree well enough to skip that safely. Building the
scikit-build-core codepath's tests, if you're changing `pyproject.toml`
itself: `python -m pytest packaging/tests -q` after the editable install.

### The zero-runtime-dependency rule

`pto-core`'s `pyproject.toml` declares `dependencies = []`, and that line is
enforced by design, not left empty by omission. The install this package
exists to make frictionless is landing on a login node or a shared
workstation without resolving against whatever `numpy` or `pydantic` version
the site has already pinned — every dependency added here is a chance for
that to fail before a user gets a FRiP number.

What follows from that for a PR:

* **A new Python runtime dependency is very unlikely to be accepted.** The
  Python layer (`packaging/pto/`) returns each engine's own JSON as plain
  dicts specifically so nothing needs a validation or schema library.
* **A new C++ dependency is vendored, not resolved from PyPI or pip.**
  `cuttag_profiler` vendors `cpp-httplib` and `nlohmann/json` under
  `modules/cuttag_profiler/third_party/`; `scrna_matrix` vendors `hnswlib`
  under `modules/scrna_matrix/third_party/` (with local patches recorded in
  its `VERSION.txt`). A C++ library a module genuinely needs goes there, with
  its license added to [`THIRD_PARTY_LICENSES.md`](THIRD_PARTY_LICENSES.md),
  not into a Python wheel dependency that would need to resolve at install
  time on a machine you don't control.
* **System libraries a module discovers at configure time (htslib, zlib,
  libdeflate, OpenMP) are a different category** — those are fine, and
  already how `cuttag_profiler`, `fastq_stream`, `genomic_toolkit` and
  `scrna_matrix` work. What's off the table is a dependency the *Python
  package* resolves at `pip install` time.

If a module you're changing genuinely needs something new, open an issue
first (the feature-request template asks for exactly this) rather than
landing it in a PR — this is one of the places where "no" is the default
answer and the reasoning above is the bar it has to clear.

## Running the test suite

```bash
ctest --test-dir build --output-on-failure   # all 32+ suites, one pass
```

Or run one suite directly while iterating:

```bash
./build/modules/peaks/test_peaks_poisson_model
ctest --test-dir build -R test_http_security --output-on-failure
```

All C++ suites are dependency-free by design — no GoogleTest, no Catch2 —
for the same reason the tools themselves are: they have to run on a cluster
that has a compiler and nothing else. Don't introduce a test framework
dependency to make a new suite more convenient to write.

### Sanitizers

Sanitizer builds are per-module, not part of the default build, and a PR
touching concurrency, pointer arithmetic, or byte-level parsing (BGZF/BAM
decoding, the lock-free queue, HNSW index mutation) should be run through
the relevant one before review:

```bash
# ASan + UBSan (fastq_stream)
cmake -S modules/fastq_stream -B build/fq-asan -DFQ_SANITIZE=ON -DFQ_NATIVE=OFF
cmake --build build/fq-asan -j && ./build/fq-asan/test_trimmer

# TSan (fastq_stream / cuttag_profiler are std::thread-based; no container needed)
cmake -S . -B build-tsan -DPTO_BUILD_SCRNA_MATRIX=OFF -DCMAKE_CXX_FLAGS=-fsanitize=thread

# TSan for scrna_matrix specifically needs OpenMP + TSan together, which no
# macOS toolchain currently provides both halves of:
docker build -f docker/Dockerfile.tsan -t pto-tsan . && docker run --rm pto-tsan
```

See the README's [Sanitizers](README.md#sanitizers) section for the two
non-obvious traps (a UBSan suppression file path containing a space, and
`AlignedAllocator`'s 64-byte rounding) before you spend time debugging a
sanitizer run that looks like a false positive.

### Python tests

```bash
python -m pytest packaging/tests -q
```

covers `packaging/pto/` — the CLI, the API wrapper, and the binary-resolution
logic. It does not exercise the C++ engines' own correctness; that's what
the `ctest` suites above are for.

[`tests/validation/`](tests/validation) is a third, separate thing: a
public-data differential harness that runs the shipped binaries as black
boxes against real ENCODE/ENA accessions and diffs them against
samtools/bedtools/deeptools/Picard/fastp/seqkit. It needs `--allow-download`
to touch the network at all and is not part of the PR-time test loop for most
changes — see the README's [Public-data
validation](README.md#public-data-validation) section for when it matters.

## A note on `scrna_matrix`

`scrna_matrix`'s CMake configuration **hard-fails without OpenMP**, on
purpose: it used to fall back to a `message(WARNING)` and a silently serial
build that the test suite reported as 100% passing. If your build is failing
to configure specifically on `scrna_matrix`, that's almost certainly missing
OpenMP, not a bug — see the README's [OpenMP is
mandatory](README.md#openmp-is-mandatory-for-scrna_matrix) section for the
per-platform fix. If you don't need this module for the change you're
making, `-DPTO_BUILD_SCRNA_MATRIX=OFF` builds and tests the other four
without it.

This is also why `scrna_matrix` ships as the separate `scrna-matrix` PyPI
package rather than inside `pto-core`: bundling a hard OpenMP requirement
into the CLI wheel would mean a macOS host without `libomp` couldn't install
`pto-core` at all, just to get a FRiP number out of `genomic_toolkit`. See
the README's ["`scrna_matrix` is a separate
install"](README.md#scrna_matrix-is-a-separate-install) section for the full
reasoning.

## Code style and pre-commit

```bash
pip install pre-commit
pre-commit install
```

installs the hooks in [`.pre-commit-config.yaml`](.pre-commit-config.yaml) to
run on every commit: trailing-whitespace/EOF/YAML hygiene, `ruff` for Python,
and `clang-format` for C++. Each module has its own `.clang-format`
inheriting shared conservative settings from the repo root's — the five
modules were developed independently and don't share an indent convention,
and the hook respects that rather than imposing one. Vendored third-party
code under `modules/*/third_party/` is excluded entirely; it isn't ours to
reformat.

Run it against everything at once with:

```bash
pre-commit run --all-files
```

## Security-sensitive changes

If you're touching `cuttag_profiler`'s embedded HTTP server
(`http_server.cpp`, `http_security.hpp`, `secrets.hpp`), BAM/BGZF decoding in
any module, or anything doing byte-level parsing of untrusted input, read
`docs/SECURITY_HTTP_2026-08-15.md` and the relevant `docs/AUDIT_*.md` first.
These are living documents — extend them rather than replacing them when you
find or fix something in their scope, and if you're fixing a defect that
looks like a pattern (an unchecked narrowing cast, a filter that fails open),
grep for the same shape elsewhere before calling the fix done. Four of five
findings in `docs/REVIEW_2026-08-15.md` were exactly that: the same defect
surviving in a sibling function the first pass didn't walk.
