# pto-core architecture diagram suite

Five diagrams of how `pto-core` is built and how data moves through it, in two forms: terminal-grade
Unicode flowcharts, and an editable graphic specification for a conference deck. Every claim on them was
read from the source at commit `0ba20d1` and, where it is behavioural, run against binaries built fresh from
that commit. Where a document and the code disagreed, the code won; the disagreements are listed below, and
as of this update they are also fixed at the source, so the docs and the diagrams now agree with the code in
the same 18 places. The suite is scoped to `pto-core` itself: Nextflow and `pto-cloud`, a separate proprietary
project, are not named on any diagram. Diagrams 1 and 4 used to draw them as dashed "outside this repo" boxes;
they are omitted entirely now, and the two boxes still drawn dashed (a shell, a browser) name no product.

| File | What it is |
|:--|:--|
| [`ascii-diagrams.md`](ascii-diagrams.md) | the five flowcharts in box-drawing characters, each with its source-of-truth list |
| [`graphic-spec.md`](graphic-spec.md) | the design system (colours, type, shapes, edges) and every node, edge and container of the five slides, with coordinates |
| [`pto-core-diagrams.drawio`](pto-core-diagrams.drawio) | the five slides as a multi-page draw.io file; containers are real parents |
| [`svg/`](svg) | the same slides as SVG, for Canva, PowerPoint, Keynote or Figma |
| [`diagram-spec.json`](diagram-spec.json) | the same model as data, for Excalidraw or scripting |
| [`src/`](src) | the generators (Python standard library only) |

| # | Diagram | Answers |
|:-:|:--|:--|
| 1 | System overview | which engines exist, who calls them, and where the process boundary is |
| 2 | `fastq_stream` streaming | how bytes become chunks, chunks become records, and what is copied |
| 3 | The four CLIs | the control path and the core loop of each |
| 4 | Integration | how Python, shells, orchestrators and a browser reach the binaries |
| 5 | Concurrency, ISA, memory | threads, SIMD dispatch, bounds and alignment |

## What was verified, and how

1. **Read, not recalled.** Every module's headers, sources, CMake, the Python wrapper, `pyproject.toml`, both
   workflows and the deploy kit were read in full or in the parts each claim comes from. Line references are on
   the diagrams. The detached worktree under `.claude/worktrees/` is a different commit and was excluded.
2. **Built fresh.** Out-of-tree, from `main`, with the wheel's own flags (Apple clang 21, htslib 1.24, arm64):

   ```bash
   cmake -S . -B /tmp/pto-verify -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=/opt/homebrew \
     -DPTO_BUILD_SCRNA_MATRIX=OFF -DFQ_ENABLE_HTS=ON -DFQ_NATIVE=OFF -DPROFILER_NATIVE_ARCH=OFF \
     -DGTK_NATIVE_ARCH=OFF -DPEAKS_WITH_BGZF=ON -DPROFILER_STATIC_HTSLIB=OFF \
     -DFQ_BUILD_TESTS=OFF -DPROFILER_BUILD_TESTS=OFF -DGTK_BUILD_TESTS=OFF -DPEAKS_BUILD_TESTS=OFF
   cmake --build /tmp/pto-verify -j
   ```

3. **Run.** The results that back the diagrams:

   | Area | Result |
   |:--|:--|
   | thread split | 14 configurations (`-t 1..12`, plain gzip and BGZF): observed inflater and worker counts match `pipeline.cpp:279-318` in every one |
   | determinism | `fastq_stream` counters identical across raw, gzip, BGZF and across NEON and scalar; `pto-peaks` output byte-identical for BAM, BED, stdin BED, stdin BAM and forced-scalar; profiler matrix and profile byte-identical at 1 and 8 threads; FRiP identical at `OMP_NUM_THREADS` 1 and 8 |
   | numbers | `frip --keep-dups` gives 148,680 in-peak fragments, the README's figure; region tiles sum to the whole file |
   | exit codes | every code on Diagram 4 reproduced, including `fastq_stream` exiting 2 for usage errors |
   | I/O guards | output equal to input refused for the same file, a hard link and a symlink; a FIFO output blocks in `open` until a reader attaches and then delivers all 239,984 lines |
   | HTTP | 14 requests against a live loopback server: 401, 421, 403, 415, 413, 400, 404, 200; Host is checked before the token; `/api/matrix` needs both tokens |
   | wrapper | `doctor`, `versions`, `bin-dir`, argv construction, exit-code mapping, `exec` pass-through |
   | artifacts | layout, tags and entry point of the built `dist/pto_core-0.1.6-*.whl` files |

4. **Not verified here**, and marked `(static)` where it appears: the `scrna_matrix` build and tests (needs
   pybind11 and an OpenMP toolchain), the x86-64 AVX2 and AVX-512 kernels (the host is arm64), the sanitizer and
   TSan jobs, CRAM input, the Docker image, and Nextflow and `pto-cloud`, which are not in this repository. The
   project's own `ctest` suites were not run; the count of 32 suites (33 with pybind11) comes from enumerating
   the `add_test` registrations, and matches the top-level README.

## Where documents and code disagree

Diagrams follow the right-hand column. None of these changes a result; each would have put a wrong number,
name or claim on a slide had the diagram been copied from the document.

| # | A document says | The code does | Evidence | Slide |
|:-:|:--|:--|:--|:-:|
| 1 | `fastq_stream/docs/ARCHITECTURE.md:180`: "There is no runtime dispatch; build once per target ISA" | Kernels carry per-function `target` attributes and are chosen once by `__builtin_cpu_supports`; `FQ_FORCE_ISA` overrides (`simd.hpp:4-40,514-568`). The README is right | ran | 5 |
| 2 | `ARCHITECTURE.md:55`: `D = spare/3` inflaters | `D = max(1, spare/3)` (`pipeline.cpp:317`): BGZF at `-t 1..5` runs one inflater, not zero, so `-t 1` spawns 5 threads | ran | 2 |
| 3 | `ARCHITECTURE.md:96`: ring padding uses `hardware_destructive_interference_size`, falling back to 128 B on arm64 and 64 B elsewhere | Deliberately not used (ODR hazard); a fixed 128 B on arm64 **and** x86-64, 64 B elsewhere (`lockfree_queue.hpp:19-42`) | static | 2 |
| 4 | `buffer.hpp:9`: free list "behind a spinlock" | `std::mutex` (`buffer.hpp:118,133`); `ARCHITECTURE.md` says mutex | static | 2 |
| 5 | `ARCHITECTURE.md:288` and the README: "No BAM input" | BAM, CRAM and SAM are read through `HtsInput` when built with `FQ_ENABLE_HTS` (off in the module, on in the wheel, `pyproject.toml:229`); a BAM is recognised only in a regular file, never on stdin | ran | 1, 2 |
| 6 | `ARCHITECTURE.md` §5 lists `sum_u8` as the mean-Phred kernel and bounds 32-bit sums at ~650 reads per chunk | `sum_u8` has no production caller (`simd.hpp:666-675`); the real bound is ~43,000 records per chunk (`phred_calculator.hpp:6-12`). The conclusion, no overflow, holds | static | 5 |
| 7 | `Buffer::seq` (`buffer.hpp:80`) reads as the ordering key | It is written (`reader.cpp:237,290,301,378`) and never read; order comes from lane round-robin | grep | 2 |
| 8 | `cuttag_profiler/README.md:321-323`: a shared BGZF pool is created for a single-threaded caller | `choose_bgzf_thread_count()` returns 0 unconditionally (`bam_reader.hpp:65-67`); no pool is ever built | static | 3, 5 |
| 9 | `deploy/Dockerfile:12-16`, `deploy/provision_aws.sh:15-17`, `scripts/run_demo.sh:335`: no authentication, session token "unimplemented" | A 256-bit session token is minted per run and enforced on every `/api/` request (`http_server.cpp:430-494`); the README and `entrypoint.sh` already say so | ran | 4 |
| 10 | `scrna_matrix/README.md:229` and the function header: O(k) per-thread memory | One `BoundedTopK` per row (O(n·k)) and a symmetric tile sweep in circle-method rounds (`knn_graph.hpp:187-241,299-413`) | static | 5 |
| 11 | `peaks/README.md:12`: `pto-peaks atac.bam peaks.narrowPeak` | Rejected: "at most one input file", exit 1 (`main.cpp:839-842`); the output needs `-o` | ran | 3 |
| 12 | `CMakeLists.txt:49-51` and the README: under `PTO_CLOUD_BUILD`, `serve` "is not a subcommand" | `serve` is recognised and refused by name, exit 1 (`main.cpp:399-415`); `--help` omits it; an `nm`-based test checks the linked binary | static | 4 |
| 13 | `packaging/README.md:250-255`, `api.py:272-277`: the dashboard is "not reachable from this package" | True for `pto.profile()` and `pto-core profile`. But nothing sets `PTO_CLOUD_BUILD` for the wheel, the shipped binary contains the server, and `pto-core exec cuttag_profiler serve --data-root D` reaches it through `os.execv` with no guard (`cli.py:83-92`) | ran | 4 |
| 14 | `pyproject.toml:1` "three command-line engines"; `:329` cites `wheels-cli.yml` and macOS 13/14 runners; `CMakeLists.txt:3` "three modules" | Four CLIs are wrapped (`_binaries.py:44`); there is no `wheels-cli.yml`, `wheels.yml` holds both jobs on `macos-15-intel` and `macos-15` with a 15.0 target; there are five modules | static | 1, 4 |
| 15 | `_cpu.py` docstring and the README: every binary defaults to `hardware_concurrency` and the wrapper passes an explicit count | Only `fastq_stream` (`-t`) and `cuttag_profiler` (`--threads`) take one from the wrapper; `genomic_toolkit` has no flag (OpenMP, `OMP_NUM_THREADS` left unset); `pto-peaks` is single-threaded | ran | 4, 5 |
| 16 | "Exit 2 means the input's shape made the answer wrong" | True for `genomic_toolkit`, `cuttag_profiler` and `pto-peaks`. `fastq_stream` also exits 2 for CLI usage errors (`main.cpp:82-85`), so `pto.run_qc(min_len=-1)` raises `InvalidInput` | ran | 1, 4 |
| 17 | Comments and docs cite `nextflow/conf/pto.config`, `nextflow/tests`, `pto-cloud` (`job_specs.py`, `TOOL_EXIT_CODES`, `RETRYABLE_EXIT_CODES`, `FORBIDDEN_SUBCOMMANDS`), `pto-cli`, `pto-studio`, `benchmarks/_measure.py`, `examples/`, `CLAUDE.md` | None of these paths exist in this repository. Diagrams 1 and 4 used to draw them as dashed "outside this repo" boxes; per an explicit decision to keep the public suite scoped to `pto-core`, they now name no external product at all | grep | 1, 4 |
| 18 | `deploy/Dockerfile:52-58`: "Only cuttag_profiler is needed" | The build disables `scrna_matrix` and `fastq_stream` only, so `genomic_toolkit` and `pto-peaks` are compiled too; only `cuttag_profiler` is copied to the runtime stage | static | n/a |

Claims that held, and are drawn as documented: the thread formula for plain gzip, the framing rules, the record
limit, the exit codes of the other three CLIs, the 148,680 FRiP figure, the suite counts, the ISA override
variables, the absence of any shared layer (four private `file_identity.hpp` copies), and that `pto-peaks` holds
its peaks until Benjamini-Hochberg runs, so only its input streams.

## Regenerating

```bash
python3 docs/diagrams/src/make_ascii.py docs/diagrams/ascii-diagrams.md
python3 docs/diagrams/src/build_graphics.py docs/diagrams
```

Edit the wording in `src/d1.py` to `src/d5.py` (ASCII) and `src/pages.py` (graphics). The ASCII canvas draws
every border and junction itself and refuses a line that would overrun its box; the graphics build refuses a
label that does not fit its node, an overlap, or a dangling edge. Both use only the Python standard library.
