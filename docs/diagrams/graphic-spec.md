# pto-core: editable graphic specification

Five 16:9 slides on a 1920 x 1080 canvas. Everything here is generated from one model,
`src/pages.py`, so the machine files and this document cannot disagree:

| File | Use it to |
|:--|:--|
| `pto-core-diagrams.drawio` | open all five pages in draw.io / diagrams.net (File > Open from > Device). Containers are real parent shapes, so moving one moves its children; edges stay connected to their nodes. |
| `svg/*.svg` | upload to Canva (Uploads, drag the five files), or place in PowerPoint, Keynote or Figma. Text stays live text. |
| `diagram-spec.json` | rebuild in Excalidraw or any scripting tool: `pages[].nodes` are rectangles (x, y, w, h, kind, text) and `pages[].edges` are arrows (from, to, exit side, entry side, via points). |

Regenerate after editing `src/pages.py`:

```bash
python3 docs/diagrams/src/build_graphics.py docs/diagrams
```

The build fails, and writes nothing, if a node leaves the page, overlaps a sibling or sits outside its
parent, if an edge names a missing node, or if a label needs more height than its box has.

## Canvas, grid and type

| Item | Value |
|:--|:--|
| Slide | 1920 x 1080 px, 16:9; outer margin 64 px; snap to an 8 px grid |
| Title band | y 28-104: title 32 px bold `#0F2A43`, subtitle 17 px `#546E7A` (one sentence) |
| Footer | y 1050, 11 px `#78909C`: source files and the commit the claims were read at |
| Text face | Inter (Helvetica or Arial as fallback); code identifiers may use JetBrains Mono |
| Sizes | container and node titles 16 px bold; body 13 px; table cells 12 px; edge labels 12 px |
| Node anatomy | title top-left, 10 px inset; body lines below, no bullets; 1.5 px stroke; radius 8 px |
| Container anatomy | title top-left in the container's own colour; children inset at least 24 px |

## Colour language

| kind | fill | stroke | text | outline | contrast | meaning |
|:--|:--|:--|:--|:--|--:|:--|
| `cpp` | `#E3F0FD` | `#1565C0` | `#0D3B72` | solid, rounded | 9.6:1 | C++20 engine or core component |
| `py` | `#E6F4EA` | `#2E7D32` | `#1B4D20` | solid, rounded | 8.7:1 | Python orchestration (stdlib only) |
| `build` | `#EFE6FA` | `#6A3FA0` | `#3F2266` | solid, rounded | 10.6:1 | build, packaging, CI, deploy |
| `data` | `#FFF4DB` | `#B7791F` | `#6B4A0E` | solid, document | 7.4:1 | artifact, format, stream |
| `sync` | `#DDF5F1` | `#00796B` | `#004D40` | solid, process | 8.6:1 | concurrency primitive: ring, pool, lock |
| `gate` | `#FDE8E8` | `#C62828` | `#7F1D1D` | solid, rounded | 8.5:1 | security gate |
| `ext` | `#F4F4F4` | `#757575` | `#424242` | dashed, rounded | 9.1:1 | outside this repository (contract only) |
| `panel` | `#FFFFFF` | `#B0BEC5` | `#37474F` | solid, container | 9.6:1 | neutral grouping container |
| `zone_cpp` | `#F5F9FF` | `#1565C0` | `#0D3B72` | solid, container | 10.5:1 | container: C++ engines |
| `zone_py` | `#F4FAF5` | `#2E7D32` | `#1B4D20` | solid, container | 9.3:1 | container: Python layer |
| `zone_gate` | `#FFF6F6` | `#C62828` | `#7F1D1D` | dashed, container | 9.4:1 | container: trust boundary |
| `boundary` | `#EFE6FA` | `#6A3FA0` | `#3F2266` | solid, band | 10.6:1 | the process boundary |
| `callout` | `#FFFBEA` | `#E0B84C` | `#5B4300` | solid, rounded | 9.0:1 | callout: a finding worth reading |

## Edge language

| kind | colour | width px | dash | meaning |
|:--|:--|--:|:--|:--|
| `payload` | `#37474F` | 2.5 | solid | bytes, buffer handles, artifacts |
| `control` | `#1565C0` | 1.6 | 6,4 | call, argv, configuration |
| `error` | `#C62828` | 1.6 | 2,4 | abort, refusal, failure exit |
| `contract` | `#6A3FA0` | 3.0 | solid | crosses the process boundary |

Arrowheads are filled triangles, 7 px. Edge labels are 12 px in the edge colour on a white halo, used only where the gap is wider than the label.

## Shape language

| Shape | Meaning |
|:--|:--|
| rounded rectangle | a stage, component or gate |
| folded-corner document | an artifact, file format or stream |
| rectangle with two vertical bars | a concurrency primitive: ring, pool, lock |
| wide horizontal band | the process boundary |
| dashed outline | outside this repository: its contract is defined here, its code is not |
| large container | grouping only; carries no behaviour |

## Accessibility

Meaning never rides on colour alone: each kind also differs by outline (solid, dashed) or shape (document,
double bar, band), and every engine node states "C++20" or "CLI" in its own text. Text-on-fill contrast is
checked at build time; the lowest pair in this palette is 7.4:1 (WCAG AAA needs 7:1). The smallest
text is 11 px and is used only for footers.

## Page 1: 1 Overview

**Slide title:** Five independent C++20 engines behind one process boundary  
**Subtitle:** No shared runtime and no shared code: callers reach an engine only through argv, files, stdio and an exit code.  
**Footer:** Sources: packaging/pto/*, pyproject.toml, CMakeLists.txt, modules/*/src/main.cpp (commit 0ba20d1)

**Talk track**

- Five engines: four command-line tools and one Python extension. Each is its own CMake project and builds and tests with the other four deleted.
- The only coupling is the process boundary: argv, files, stdio and an exit code. There is no shared runtime, and the output-aliasing guard is deliberately copied four times.
- What's left dashed here is generic -- a shell and a browser -- not any named external product; this slide stays inside what pto-core itself ships.

**Containers (grouping hierarchy)**

| id | kind | x | y | w | h | title | children |
|:--|:--|--:|--:|--:|--:|:--|--:|
| `z_call` | panel | 64 | 124 | 1792 | 140 | CALLERS | 5 |
| `z_py` | zone_py | 64 | 292 | 912 | 270 | PYTHON ORCHESTRATION   packaging/pto   stdlib only, zero runtime dependencies | 6 |
| `band` | boundary | 64 | 590 | 1792 | 48 | PROCESS BOUNDARY   argv + files + stdio + exit code   no shared memory   no IPC | 0 |
| `z_eng` | zone_cpp | 64 | 666 | 1792 | 250 | C++20 ENGINES | 5 |

**Nodes**

| id | kind | parent | x | y | w | h | text (title / body lines) |
|:--|:--|:--|--:|--:|--:|--:|:--|
| `c_cli` | py | z_call | 88 | 164 | 272 | 88 | **pto-core CLI** / console script pto.cli:main |
| `c_api` | py | z_call | 384 | 164 | 272 | 88 | **import pto** / Python API, 7 functions |
| `c_path` | ext | z_call | 1000 | 164 | 196 | 88 | **shell / SLURM** / bare tool names on PATH |
| `c_brw` | ext | z_call | 1318 | 164 | 196 | 88 | **browser** / dashboard, loopback only |
| `c_scp` | ext | z_call | 1636 | 164 | 196 | 88 | **scanpy / AnnData** / import scrna_matrix |
| `py_cli` | py | z_py | 88 | 344 | 272 | 96 | **cli.py** / argparse, 12 subcommands / exec uses os.execv |
| `py_api` | py | z_py | 384 | 344 | 272 | 96 | **api.py** / 7 functions / None kwargs are omitted |
| `py_cpu` | py | z_py | 680 | 344 | 272 | 96 | **_cpu.py** / min(affinity, cgroup quota, / cpu_count), never below 1 |
| `py_bin` | py | z_py | 88 | 452 | 272 | 96 | **_binaries.py** / $PTO_BIN_DIR > pto/_vendor/bin / > $PATH |
| `py_run` | py | z_py | 384 | 452 | 272 | 96 | **_run.py** / Popen, shell=False / os.wait4: wall time + peak RSS |
| `py_err` | py | z_py | 680 | 452 | 272 | 96 | **errors.py** / rc 2 = InvalidInput / other rc = ToolError |
| `e_fq` | cpp | z_eng | 88 | 712 | 328 | 190 | **fastq_stream 0.1.0** / CLI: streaming QC and trimming / in: FASTQ raw/gz/BGZF; BAM/CRAM/SAM / out: QC JSON + text; FASTQ to file/FIFO/stdout / threads 3+D+W; SIMD runtime / exit 0 \| 1 \| 2 (usage errors too) |
| `e_gtk` | cpp | z_eng | 440 | 712 | 328 | 190 | **genomic_toolkit 0.1.0** / CLI: sizes \| markdup \| frip / in: BAM/CRAM, BEDPE, fragment BED / out: JSON on stdout (+ histogram TSV) / threads 1 + OpenMP (FRiP); no hand SIMD / exit 0 \| 1 \| 2 |
| `e_pk` | cpp | z_eng | 792 | 712 | 328 | 190 | **pto-peaks 0.1.0** / CLI: streaming peak caller / in: fragment BED or sorted BAM (stdin ok) / out: narrowPeak; BAM decoded natively / threads 1; SIMD runtime (Poisson) / exit 0 \| 1 \| 2 |
| `e_cp` | cpp | z_eng | 1144 | 712 | 328 | 190 | **cuttag_profiler 1.0.0** / CLI: profile \| serve / in: indexed BAM + BED/GTF / out: matrix.tsv + profile.tsv, or HTTP / threads N workers; no hand SIMD / exit 0 \| 1 \| 2 |
| `e_sc` | cpp | z_eng | 1496 | 712 | 328 | 190 | **scrna_matrix** / header-only C++20 + pybind11 module / in: NumPy/SciPy CSR buffers / out: KnnGraph as NumPy views / threads OpenMP; SIMD runtime / no CLI, no exit code; own wheel |
| `l_fq` | ext |  | 88 | 928 | 328 | 56 | zlib + libdeflate / htslib optional (ON in the wheel) |
| `l_gtk` | ext |  | 440 | 928 | 328 | 56 | htslib (AUTO-detected) / OpenMP optional |
| `l_pk` | ext |  | 792 | 928 | 328 | 56 | libdeflate or zlib (BAM/BGZF) / no htslib |
| `l_cp` | ext |  | 1144 | 928 | 328 | 56 | htslib (required) / httplib + json (server only) |
| `l_sc` | ext |  | 1496 | 928 | 328 | 56 | OpenMP (required, hard fail) / hnswlib (vendored, patched) |
| `build` | build |  | 64 | 996 | 1792 | 44 | cmake -S . : PTO_BUILD_* x5 (PTO_CLOUD_BUILD off)  ->  wheel pto-core py3-none-<platform>, 4 CLIs in pto/_vendor/bin  \|  wheel scrna-matrix cp39-cp313  \|  four private copies of file_identity.hpp |

**Edges (directed)**

| id | from | exit | to | entry | kind | label | via |
|:--|:--|:-:|:--|:-:|:--|:--|:--|
| e1 | `c_cli` | b | `py_cli` | t | control |  |  |
| e2 | `c_api` | b | `py_api` | t | control |  |  |
| e3 | `py_cli` | r | `py_api` | l | control |  |  |
| e4 | `py_cpu` | l | `py_api` | r | control |  |  |
| e5 | `py_api` | b | `py_run` | t | control |  |  |
| e6 | `py_run` | l | `py_bin` | r | control |  |  |
| e7 | `py_run` | r | `py_err` | l | error |  |  |
| e8 | `py_run` | b | `band` | t | contract |  |  |
| e9 | `c_path` | b | `band` | t | contract | bare tool names |  |
| e10 | `band` | b | `e_fq` | t | contract |  |  |
| e11 | `band` | b | `e_gtk` | t | contract |  |  |
| e12 | `band` | b | `e_pk` | t | contract |  |  |
| e13 | `band` | b | `e_cp` | t | contract |  |  |
| e14 | `c_brw` | b | `e_cp` | t | payload | HTTP 127.0.0.1 | (1416,690) > (1450,690) |
| e15 | `c_scp` | b | `e_sc` | t | payload | in-process import |  |

## Page 2: 2 Streaming

**Slide title:** fastq_stream moves buffer handles, not bytes  
**Subtitle:** 3 + D + W threads joined by lock-free SPSC rings; a 256 KiB payload is copied at four points only.  
**Footer:** Sources: fastq_stream/src/pipeline.cpp:263-613, reader.cpp:32-380, buffer.hpp, lockfree_queue.hpp (commit 0ba20d1); thread math checked at 14 configurations

**Talk track**

- Buffers, not bytes, move between threads: 16-slot single-producer single-consumer rings hand Buffer handles across, and the payload is copied at four points only.
- Order comes from topology, not sequence numbers: chunk k rides lane k mod D, and the assembler and writer read lanes in the same round-robin order.
- Memory is bounded by the number of lanes, not the input: 15.8 MB of pool for a 12.8 MB file at -t 6, measured.

**Containers (grouping hierarchy)**

| id | kind | x | y | w | h | title | children |
|:--|:--|--:|--:|--:|--:|:--|--:|
| `p_buf` | panel | 64 | 632 | 880 | 184 | BUFFER ANATOMY   buffer.hpp:28-35 | 3 |
| `p_ring` | panel | 976 | 632 | 880 | 184 | SPSC RING   lockfree_queue.hpp:47 | 6 |
| `p_copy` | panel | 64 | 832 | 880 | 184 | THE FOUR COPY POINTS | 4 |
| `p_ctl` | panel | 976 | 832 | 880 | 184 | CONTROL AND BOUNDS | 3 |

**Nodes**

| id | kind | parent | x | y | w | h | text (title / body lines) |
|:--|:--|:--|--:|--:|--:|--:|:--|
| `src` | data |  | 64 | 136 | 232 | 200 | **SOURCE** / path \| - \| FIFO / BAM/CRAM/SAM by content / (regular files only) / framing: 1f 8b + BC = BGZF / 1f 8b alone = gzip; else raw |
| `read` | cpp |  | 336 | 136 | 360 | 200 | **READER  1 thread** / FdSource 1 MiB stage + fadvise / BGZF: <= 4 whole members per Buffer (c1) / gzip: zlib inflate into the Buffer / raw: memcpy stage -> Buffer (c1) / BAM: htslib decode pool, FASTQ text / EOF: BGZF must end in its 28 B block |
| `ringA` | sync |  | 736 | 136 | 200 | 200 | **lane x D** / SpscRing<Buffer,16> / 15 usable slots / closed flag / BGZF only |
| `infl` | cpp |  | 976 | 136 | 320 | 200 | **INFLATERS  D threads** / libdeflate, 1 decompressor each / 4 members -> one 256 KiB Buffer / compressed Buffer -> in_pool / D = max(1, spare/3) |
| `ringB` | sync |  | 1336 | 136 | 200 | 200 | **lane x D** / pop order = push order / chunk k on lane k mod D |
| `asm` | cpp |  | 1526 | 400 | 330 | 200 | **ASSEMBLER  1 thread** / pop lane[k mod D] or the direct lane / carry spliced into 64 KiB prefix (c2) / cut at the 4th newline (memchr) / whole records -> lane j mod W / partial tail -> carry (max 64 KiB) |
| `ringC` | sync |  | 1290 | 400 | 200 | 200 | **lane x W** / input of worker i / chunk j on lane j mod W |
| `work` | cpp |  | 890 | 400 | 360 | 200 | **WORKERS  W threads** / out Buffer from out_pool (320 KiB) / checks -> Phred+64 rebase -> QC on RAW / trim -> adapters -> filters / survivors memcpy'd into out (c3) / SIMD: runtime ISA kernels / always ONE out Buffer per chunk |
| `ringD` | sync |  | 650 | 400 | 200 | 200 | **lane x W** / output of worker i / one Buffer per input chunk |
| `writ` | cpp |  | 364 | 400 | 250 | 200 | **WRITER  1 thread** / pop lane k mod W in order / write(2), EINTR retried (c4) / SIGPIPE ignored; close checked |
| `sink` | data |  | 64 | 400 | 260 | 200 | **SINK** / stdout \| file \| FIFO / --mkfifo: open blocks / until a reader attaches / refused: -i X -o X |
| `b_pre` | data | p_buf | 88 | 676 | 280 | 96 | **kPrefix 64 KiB** / headroom in front of data() / extend_front() splices the carry |
| `b_pay` | data | p_buf | 368 | 676 | 552 | 96 | **payload 256 KiB (in_pool)   320 KiB (out_pool)** / data() ... size() ... capacity() / 64 B aligned by aligned_alloc |
| `b_note` | note | p_buf | 88 | 780 | 832 | 28 | handle = {pool*, base*, cap, size, front_} + seq (written, never read); free list under std::mutex |
| `r0` | sync | p_ring | 1000 | 676 | 160 | 96 | **head_** / consumer |
| `r1` | sync | p_ring | 1172 | 676 | 160 | 96 | **cached_tail_** / consumer-private |
| `r2` | sync | p_ring | 1344 | 676 | 160 | 96 | **tail_** / producer |
| `r3` | sync | p_ring | 1516 | 676 | 160 | 96 | **cached_head_** / producer-private |
| `r4` | sync | p_ring | 1688 | 676 | 160 | 96 | **slots_[16]** / Buffer handles / 15 usable |
| `r_note` | note | p_ring | 1000 | 780 | 832 | 28 | each member alignas(128), so producer and consumer share no line; kCacheLine = 128 B on arm64 and x86-64 |
| `c1` | callout | p_copy | 88 | 876 | 196 | 124 | **(c1) stage -> Buffer** / raw and BGZF compressed batch; gzip inflates in place |
| `c2` | callout | p_copy | 300 | 876 | 196 | 124 | **(c2) carry splice** / only the straddling partial record, at most 64 KiB |
| `c3` | callout | p_copy | 512 | 876 | 196 | 124 | **(c3) serialise** / worker writes surviving records into the out Buffer |
| `c4` | callout | p_copy | 724 | 876 | 196 | 124 | **(c4) write(2)** / the kernel's copy to file, FIFO or stdout |
| `k_gate` | sync | p_ctl | 1000 | 876 | 272 | 124 | **ErrorGate** / first error wins / abort flag polled in every ring wait / then join, then throw |
| `k_term` | sync | p_ctl | 1284 | 876 | 272 | 124 | **Termination** / producer release-stores closed / consumer exits on empty AND closed |
| `k_mem` | data | p_ctl | 1568 | 876 | 272 | 124 | **Memory bound** / (lanes x 15 + in-flight) x <= 384 KiB / measured 47 blocks = 15.8 MB, 12.8 MB FASTQ, -t 6 |

**Edges (directed)**

| id | from | exit | to | entry | kind | label | via |
|:--|:--|:-:|:--|:-:|:--|:--|:--|
| e1 | `src` | r | `read` | l | payload |  |  |
| e2 | `read` | r | `ringA` | l | payload |  |  |
| e3 | `ringA` | r | `infl` | l | payload |  |  |
| e4 | `infl` | r | `ringB` | l | payload |  |  |
| e5 | `ringB` | r | `asm` | t | payload |  |  |
| e6 | `read` | b | `asm` | t | payload | direct lane (gzip, raw, BAM) | (516,376) > (1608,376) |
| e7 | `asm` | l | `ringC` | r | payload |  |  |
| e8 | `ringC` | l | `work` | r | payload |  |  |
| e9 | `work` | l | `ringD` | r | payload |  |  |
| e10 | `ringD` | l | `writ` | r | payload |  |  |
| e11 | `writ` | l | `sink` | r | payload |  |  |

## Page 3: 3 Engines

**Slide title:** Each engine is one pass: a control path and a core loop  
**Subtitle:** Top row of each panel is main(); the row below is the loop that touches every record or region.  
**Footer:** Sources: modules/*/src/main.cpp and the files named in each panel (commit 0ba20d1)

**Talk track**

- Every engine is one streaming pass with no sort and no temporary file. The one place results are held is pto-peaks' candidate list, because Benjamini-Hochberg needs every p-value.
- Inputs that would otherwise print a plausible wrong number exit 2 instead: naming mismatch, unsorted input, every region skipped, Phred+64.

**Containers (grouping hierarchy)**

| id | kind | x | y | w | h | title | children |
|:--|:--|--:|--:|--:|--:|:--|--:|
| `pa` | panel | 64 | 124 | 880 | 436 | 1  fastq_stream   qc | trim   (main.cpp, pipeline.cpp) | 11 |
| `pb` | panel | 976 | 124 | 880 | 436 | 2  genomic_toolkit   sizes | markdup | frip | 11 |
| `pc` | panel | 64 | 580 | 880 | 436 | 3  cuttag_profiler   profile   (main.cpp, signal_calc.cpp, bam_reader.cpp) | 11 |
| `pd` | panel | 976 | 580 | 880 | 436 | 4  pto-peaks   (main.cpp, sliding_window / caller / poisson_model / bam_streamer .hpp) | 9 |

**Nodes**

| id | kind | parent | x | y | w | h | text (title / body lines) |
|:--|:--|:--|--:|--:|--:|--:|:--|
| `a0` | cpp | pa | 88 | 188 | 150 | 100 | **argv** / -i -o -j -t / --qc-only ... |
| `a1` | cpp | pa | 258 | 188 | 150 | 100 | **validate** / ranges, adapters; / usage error = rc 2 |
| `a2` | cpp | pa | 428 | 188 | 150 | 100 | **probe** / BAM/CRAM by content / (regular files) |
| `a3` | cpp | pa | 598 | 188 | 150 | 100 | **run_pipeline** / see page 2: / 3 + D + W threads |
| `a4` | cpp | pa | 768 | 188 | 150 | 100 | **report + rc** / text, JSON (-j) / rc 0 \| 1 \| 2 |
| `b0` | cpp | pa | 88 | 344 | 150 | 110 | **checks** / header, SEQ, QUAL / printable; / Phred+64 rebase |
| `b1` | cpp | pa | 258 | 344 | 150 | 110 | **QC on RAW** / counts, Q20/Q30 / GC, N, per cycle |
| `b2` | cpp | pa | 428 | 344 | 150 | 110 | **trim** / front, tail, / sliding window |
| `b3` | cpp | pa | 598 | 344 | 150 | 110 | **adapters** / seed + extend / TruSeq, Nextera |
| `b4` | cpp | pa | 768 | 344 | 150 | 110 | **filters** / min_len, max_n, / min_mean_q |
| `pa_n` | note | pa | 88 | 478 | 832 | 40 | rc 2 also when the qualities fit only Phred+64 and no --phred-offset was given (the report is printed first) |
| `c0` | cpp | pb | 1000 | 188 | 190 | 100 | **argv** / sizes\|markdup\|frip / checked int parse |
| `c1` | cpp | pb | 1210 | 188 | 190 | 100 | **peaks [frip]** / read_peaks_bed / malformed = fatal |
| `c2` | cpp | pb | 1420 | 188 | 190 | 100 | **PeakSet::build** / sort, merge overlap / and book-ended; SoA |
| `c3` | cpp | pb | 1630 | 188 | 190 | 100 | **bind** / pre-seed contigs / peak tid -> stream tid |
| `d0` | cpp | pb | 1000 | 304 | 190 | 110 | **stream_fragments** / .bam .cram .sam: htslib; / else .bedpe; else BED |
| `d1` | cpp | pb | 1210 | 304 | 190 | 110 | **reader** / BAM: TLEN>0 = 1 pair / text: 1 MiB, cap 64 MiB/line |
| `d2` | cpp | pb | 1420 | 304 | 190 | 110 | **BatchEmitter** / Fragment = 16 B / 4096 per batch = 64 KiB / size histogram here |
| `d3` | cpp | pb | 1630 | 304 | 190 | 110 | **consumers** / sizes: histogram / markdup: marker / frip: marker + count |
| `dm` | cpp | pb | 1000 | 434 | 408 | 96 | **DuplicateMarker** / sorted: keys (end, strand) of ONE start / hash (--unsorted): 64-bit fingerprints / sets kDuplicate in place |
| `fc` | cpp | pb | 1424 | 434 | 408 | 96 | **FripCounter.add_batch** / omp for if batch >= 512; per-thread cursor / 5 integer sums; checked_add outside |
| `pb_n` | note | pb | 1000 | 538 | 832 | 20 | exit 2: nothing parsed \| out of order (markdup; frip without --keep-dups) \| no overlap and bad naming |
| `f0` | cpp | pc | 88 | 644 | 150 | 100 | **guards** / output != input / limits checked |
| `f1` | cpp | pc | 258 | 644 | 150 | 100 | **regions** / BED/GTF (.gz) / names sanitised |
| `f2` | cpp | pc | 428 | 644 | 150 | 100 | **BamReader** / header, index; / mapped total from / index metadata |
| `f3` | cpp | pc | 598 | 644 | 150 | 100 | **compute** / compute_matrix: / N<=256, N<=rows; / row blocks |
| `f4` | cpp | pc | 768 | 644 | 150 | 100 | **finish** / scale, means, / write TSV; / all skipped: rc 2 |
| `g0` | cpp | pc | 88 | 800 | 150 | 110 | **window** / TSS/center/TES / anchor; minus / strand mirrored |
| `g1` | cpp | pc | 258 | 800 | 150 | 110 | **skip?** / no such contig, / or off the contig: / row = missing |
| `g2` | cpp | pc | 428 | 800 | 150 | 110 | **query_bins** / widen by max fragment; / sam_itr_queryi |
| `g3` | cpp | pc | 598 | 800 | 150 | 110 | **alignment** / flags + MAPQ; / leftmost mate: / TLEN span; / reads + bases |
| `g4` | cpp | pc | 768 | 800 | 150 | 110 | **row** / fill bins; / minus strand reversed: / 5' to 3' |
| `pc_n` | note | pc | 88 | 934 | 832 | 60 | matrix = regions x bins doubles, allocated up front (bins <= 10,000,000) / means: min(N, rows/64, 64 MiB/(cols*8)) reducers, cache-line aligned / skipped rows stay in the matrix but leave means, BPM total and scaling |
| `h0` | cpp | pd | 1000 | 644 | 190 | 100 | **detect** / peek first byte: / 1f = BAM, else BED / (stdin ok) |
| `h1` | cpp | pd | 1210 | 644 | 190 | 100 | **records** / BAM: BGZF + CRC32, / leftmost mate TLEN>0 / BED: 4 MiB blocks |
| `h2` | cpp | pd | 1420 | 644 | 190 | 100 | **Driver.add** / order and bounds / checked; seek across / coverage deserts |
| `h3` | cpp | pd | 1630 | 644 | 190 | 100 | **finalize** / BH over ALL / candidates, then / narrowPeak |
| `i0` | cpp | pd | 1000 | 800 | 190 | 110 | **MultiScaleWindow** / delta + pileup rings, / 2 x 16384 x int64 / = 256 KiB |
| `i1` | cpp | pd | 1210 | 800 | 190 | 110 | **3 backgrounds** / centred 1k, 5k, 10k / running sums; / lag = 5000 bases |
| `i2` | cpp | pd | 1420 | 800 | 190 | 110 | **PeakCaller** / prefilter, queue of 256, / SIMD Poisson -log10 p |
| `i3` | cpp | pd | 1630 | 800 | 190 | 110 | **state machine** / IDLE, IN_PEAK, / SUMMIT_SEARCH, / CLOSING (max_gap) |
| `pd_n` | note | pd | 1000 | 934 | 832 | 60 | lambda = max(genome floor, 3 local rates); a degenerate lambda scores p = 1 (fails closed) / peaks are held (32 B each) until finalize(): BH needs every p-value, so only the INPUT streams / exit 2: unsorted, contig not in sizes, bad field, BED fragment > 4096 bp (BAM: counted skipped) |

**Edges (directed)**

| id | from | exit | to | entry | kind | label | via |
|:--|:--|:-:|:--|:-:|:--|:--|:--|
| e1 | `a0` | r | `a1` | l | payload |  |  |
| e2 | `a1` | r | `a2` | l | payload |  |  |
| e3 | `a2` | r | `a3` | l | payload |  |  |
| e4 | `a3` | r | `a4` | l | payload |  |  |
| e5 | `b0` | r | `b1` | l | payload |  |  |
| e6 | `b1` | r | `b2` | l | payload |  |  |
| e7 | `b2` | r | `b3` | l | payload |  |  |
| e8 | `b3` | r | `b4` | l | payload |  |  |
| e9 | `a3` | b | `b0` | t | control | per record, in each of W workers |  |
| e10 | `c0` | r | `c1` | l | payload |  |  |
| e11 | `c1` | r | `c2` | l | payload |  |  |
| e12 | `c2` | r | `c3` | l | payload |  |  |
| e13 | `d0` | r | `d1` | l | payload |  |  |
| e14 | `d1` | r | `d2` | l | payload |  |  |
| e15 | `d2` | r | `d3` | l | payload |  |  |
| e16 | `c1` | b | `d0` | t | control |  |  |
| e17 | `d3` | b | `dm` | t | payload |  |  |
| e18 | `d3` | b | `fc` | t | payload |  |  |
| e19 | `f0` | r | `f1` | l | payload |  |  |
| e20 | `f1` | r | `f2` | l | payload |  |  |
| e21 | `f2` | r | `f3` | l | payload |  |  |
| e22 | `f3` | r | `f4` | l | payload |  |  |
| e23 | `g0` | r | `g1` | l | payload |  |  |
| e24 | `g1` | r | `g2` | l | payload |  |  |
| e25 | `g2` | r | `g3` | l | payload |  |  |
| e26 | `g3` | r | `g4` | l | payload |  |  |
| e27 | `f3` | b | `g0` | t | control | per region, in each of N workers |  |
| e28 | `h0` | r | `h1` | l | payload |  |  |
| e29 | `h1` | r | `h2` | l | payload |  |  |
| e30 | `h2` | r | `h3` | l | payload |  |  |
| e31 | `i0` | r | `i1` | l | payload |  |  |
| e32 | `i1` | r | `i2` | l | payload |  |  |
| e33 | `i2` | r | `i3` | l | payload |  |  |
| e34 | `h2` | b | `i0` | t | control | per base, one thread |  |

## Page 4: 4 Integration

**Slide title:** Everything outside the engines meets them at one contract  
**Subtitle:** argv in; JSON, files and an exit code out. The wrapper adds resolution, thread budgeting and error typing.  
**Footer:** Sources: packaging/pto/*, pyproject.toml, ci.yml, wheels.yml, http_server.cpp:401-803, deploy/* (commit 0ba20d1); dashboard probed live

**Talk track**

- The Python wrapper is thin: it resolves a binary, budgets threads, runs argv without a shell, and types the exit code.
- The dashboard is the only network path: Host, then Origin, then a 256-bit token, all in one pre-routing handler, then a body cap and path confinement.
- A finding: the wheel's cuttag_profiler still contains serve, reachable through pto-core exec, although the Python API never builds it.

**Containers (grouping hierarchy)**

| id | kind | x | y | w | h | title | children |
|:--|:--|--:|--:|--:|--:|:--|--:|
| `pA` | panel | 64 | 124 | 1100 | 480 | A   FOUR WAYS IN | 15 |
| `pB` | panel | 1196 | 124 | 660 | 480 | B   RESOLUTION AND THREADS | 6 |
| `tbl` | panel | 64 | 620 | 1100 | 274 |  | 0 |
| `pD` | panel | 1196 | 620 | 660 | 274 | D   BUILD, PACKAGE, PUBLISH | 6 |
| `pE` | zone_gate | 64 | 908 | 1792 | 140 | TRUST BOUNDARY   127.0.0.1 only   Host, then Origin, then token, in one pre-routing handler for /api/* | 8 |

**Nodes**

| id | kind | parent | x | y | w | h | text (title / body lines) |
|:--|:--|:--|--:|--:|--:|--:|:--|
| `lab0` | note | pA | 80 | 202 | 112 | 24 | PYTHON API |
| `l00` | py | pA | 200 | 168 | 210 | 96 | **pto.frip(...)** / 7 functions / None kwargs omitted |
| `l01` | py | pA | 430 | 168 | 210 | 96 | **run_tool** / resolve binary / argv, shell=False |
| `l02` | py | pA | 660 | 168 | 210 | 96 | **Popen + wait4** / os.wait4: rc + max RSS / stdio to temp files |
| `l03` | py | pA | 890 | 168 | 210 | 96 | **result** / rc 0: Result / rc 2: InvalidInput / else ToolError |
| `lab1` | note | pA | 80 | 362 | 112 | 24 | pto-core |
| `l10` | py | pA | 200 | 328 | 210 | 96 | **pto-core frip** / argparse subcommand / --json optional |
| `l11` | py | pA | 430 | 328 | 210 | 96 | **cli._cmd_frip** / calls pto.api.frip |
| `l12` | py | pA | 660 | 328 | 210 | 96 | **output** / metrics: stdout / cost: stderr |
| `l13` | py | pA | 890 | 328 | 210 | 96 | **exit code** / 0 ok \| 1 error / 2 InvalidInput / 141 broken pipe |
| `lab2` | note | pA | 80 | 522 | 112 | 24 | exec / PATH |
| `l20` | py | pA | 200 | 488 | 210 | 96 | **pto-core exec** / os.execv: no Python / frame remains |
| `l21` | py | pA | 430 | 488 | 210 | 96 | **bin-dir** / prints the engine / directory |
| `l22` | py | pA | 660 | 488 | 210 | 96 | **PATH** / export PATH= / $(pto-core bin-dir):$PATH |
| `l23` | ext | pA | 890 | 488 | 210 | 96 | **shell, SLURM, Nextflow** / bare tool names |
| `r1` | py | pB | 1220 | 176 | 190 | 88 | **1  $PTO_BIN_DIR** / override; no reinstall |
| `r2` | py | pB | 1430 | 176 | 190 | 88 | **2  pto/_vendor/bin** / wheel-vendored |
| `r3` | py | pB | 1640 | 176 | 190 | 88 | **3  $PATH** / shutil.which |
| `r_err` | gate | pB | 1220 | 296 | 610 | 84 | **none found: BinaryNotFound** / lists the searched directories and prints the cmake build recipe |
| `r_tools` | callout | pB | 1220 | 400 | 610 | 76 | **TOOLS = the four engines** / a test fails if the wheel vendors a fifth binary that TOOLS lacks |
| `r_thr` | callout | pB | 1220 | 492 | 610 | 96 | **thread budget** / min(sched_getaffinity, cgroup quota, cpu_count), at least 1 / passed as -t (fastq_stream) and --threads (cuttag_profiler) only / genomic_toolkit and pto-peaks get none; OMP_NUM_THREADS is not set |
| `d1` | build | pD | 1220 | 668 | 196 | 96 | **cmake -S .** / PTO_BUILD_* x5 / PTO_CLOUD_BUILD off |
| `d2` | build | pD | 1428 | 668 | 196 | 96 | **scikit-build-core** / install-dir pto/_vendor / headers dropped |
| `d3` | build | pD | 1636 | 668 | 196 | 96 | **cibuildwheel** / cp312 pin, py3-none-* / doctor + pytest |
| `d4` | build | pD | 1220 | 784 | 196 | 96 | **repair** / delocate (macOS) / auditwheel (Linux) |
| `d5` | build | pD | 1428 | 784 | 196 | 96 | **publish** / v* tag -> PyPI / OIDC, no token |
| `d6` | callout | pD | 1636 | 784 | 196 | 96 | **vendored libs** / macOS: libomp, libhts / Linux: libgomp, static htslib |
| `w0` | ext | pE | 88 | 952 | 200 | 88 | **client** / browser or curl / URL carries ?t=TOKEN |
| `w1` | gate | pE | 308 | 952 | 200 | 88 | **Host  421** / loopback name + port / defeats DNS rebinding |
| `w2` | gate | pE | 528 | 952 | 200 | 88 | **Origin  403** / absent ok; present must / be this server |
| `w3` | gate | pE | 748 | 952 | 200 | 88 | **token  401** / ?t= or Bearer / constant-time compare |
| `w4` | gate | pE | 968 | 952 | 200 | 88 | **body gates** / <= 1 MiB: 413 / JSON type: 415 |
| `w5` | gate | pE | 1188 | 952 | 200 | 88 | **resolve_input** / under data root: 400 / open once, on the fd |
| `w6` | cpp | pE | 1408 | 952 | 200 | 88 | **compute** / compute_matrix / HTTP pool = 16 threads |
| `w7` | cpp | pE | 1628 | 952 | 200 | 88 | **cache + fetch** / MatrixCache 512 MiB, 4 / GET /api/matrix?token=M |

**Table `tbl`** (columns [190, 220, 200, 290, 200] px)

| engine | machine-readable | artifacts | exit codes | thread control |
|:--|:--|:--|:--|:--|
| fastq_stream<br>qc \| trim | JSON file: -j PATH<br>text: --report PATH<br>wrapper reads the JSON | trimmed FASTQ:<br>-o FILE \| FIFO \| - | 0 ok, 1 runtime<br>2 usage error, or<br>Phred+64 suspected | -t N: wrapper passes it |
| genomic_toolkit<br>sizes markdup frip | JSON on stdout (--json) | histogram TSV:<br>--histogram FILE | 0 ok, 1 usage or system<br>2 input shape: naming,<br>order, nothing parsed | no flag; OpenMP reads<br>OMP_NUM_THREADS<br>(wrapper leaves it unset) |
| cuttag_profiler<br>profile \| serve | none: summary on stderr | matrix.tsv, profile.tsv<br>(or HTTP: serve) | 0 ok, 1 error<br>2 every region skipped | --threads N: wrapper<br>passes it |
| pto-peaks | none: counts on stderr | narrowPeak: --out FILE<br>(default stdout) | 0 ok, 1 usage or I/O<br>2 unsorted, contig missing,<br>bad field | single-threaded;<br>no thread option |

**Edges (directed)**

| id | from | exit | to | entry | kind | label | via |
|:--|:--|:-:|:--|:-:|:--|:--|:--|
| e1 | `l00` | r | `l01` | l | control |  |  |
| e2 | `l01` | r | `l02` | l | control |  |  |
| e3 | `l02` | r | `l03` | l | control |  |  |
| e4 | `l10` | r | `l11` | l | control |  |  |
| e5 | `l11` | r | `l12` | l | control |  |  |
| e6 | `l12` | r | `l13` | l | control |  |  |
| e7 | `l20` | r | `l21` | l | control |  |  |
| e8 | `l21` | r | `l22` | l | control |  |  |
| e9 | `l22` | r | `l23` | l | control |  |  |
| e10 | `r1` | r | `r2` | l | control |  |  |
| e11 | `r2` | r | `r3` | l | control |  |  |
| e12 | `r3` | b | `r_err` | t | error |  |  |
| e13 | `d1` | r | `d2` | l | control |  |  |
| e14 | `d2` | r | `d3` | l | control |  |  |
| e15 | `d3` | b | `d4` | t | control |  |  |
| e16 | `d4` | r | `d5` | l | control |  |  |
| e17 | `w0` | r | `w1` | l | control |  |  |
| e18 | `w1` | r | `w2` | l | control |  |  |
| e19 | `w2` | r | `w3` | l | control |  |  |
| e20 | `w3` | r | `w4` | l | control |  |  |
| e21 | `w4` | r | `w5` | l | control |  |  |
| e22 | `w5` | r | `w6` | l | control |  |  |
| e23 | `w6` | r | `w7` | l | control |  |  |

## Page 5: 5 Concurrency and ISA

**Slide title:** Determinism is tested, not assumed: threads, ISA and memory  
**Subtitle:** (ran) = executed in this review on an arm64 host; (static) = read from source or CI configuration only.  
**Footer:** Sources: simd.hpp:4-568, poisson_model.hpp:1715-1889, simd_math.hpp:704-770, knn_graph.hpp:131-427, signal_calc.cpp, frip.hpp, ci.yml (commit 0ba20d1)

**Talk track**

- One dispatch pattern, three engines: fastq_stream, pto-peaks and scrna_matrix. genomic_toolkit and cuttag_profiler contain no hand-written SIMD.
- Determinism was checked by running: identical results across framings, ISAs and thread counts.
- There is no mmap anywhere in first-party code; memory is bounded by structure, not by a mapping.

**Containers (grouping hierarchy)**

| id | kind | x | y | w | h | title | children |
|:--|:--|--:|--:|--:|--:|:--|--:|
| `tA` | panel | 64 | 124 | 1792 | 224 |  | 0 |
| `tB` | panel | 64 | 372 | 1240 | 198 |  | 0 |
| `mech` | panel | 1336 | 372 | 520 | 268 | ONE DISPATCH MECHANISM | 4 |
| `tC` | panel | 64 | 664 | 1792 | 186 |  | 0 |

**Nodes**

| id | kind | parent | x | y | w | h | text (title / body lines) |
|:--|:--|:--|--:|--:|--:|--:|:--|
| `m0` | sync | mech | 1360 | 420 | 472 | 44 | **1  FORCE_ISA env, if that ISA is compiled in** |
| `m1` | sync | mech | 1360 | 476 | 472 | 44 | **2  __builtin_cpu_supports (also OS save)** |
| `m2` | sync | mech | 1360 | 532 | 472 | 44 | **3  NEON on AArch64: no probe** |
| `m3` | sync | mech | 1360 | 588 | 472 | 44 | **4  scalar reference** |

**Table `tA`** (columns [170, 520, 560, 542] px)

| engine | threads and primitive | shared state and synchronisation | ordering and determinism evidence |
|:--|:--|:--|:--|
| fastq_stream | 3+D+W std::thread: reader, D inflaters, assembler, W workers, writer; BAM adds an htslib decode pool | SPSC rings with acquire/release atomics, no CAS; buffer pools under std::mutex; ErrorGate = mutex + abort atomic | lane order = stream order; raw, gzip, BGZF and NEON vs scalar give 1 identical result (ran) |
| genomic_toolkit | 1 streaming thread; OpenMP parallel for in FRiP when the batch has >= 512 fragments, schedule(static) | PeakSet is immutable; per-thread OverlapCursor; integer reductions; checked_add outside the parallel region | identical JSON at OMP_NUM_THREADS=1 and 8; region tiles sum to the whole file (ran) |
| cuttag_profiler | profile: N std::thread workers, block-partitioned rows, then min(N, rows/64, memory) reducers. serve: httplib pool of 16 | one leased htsFile per worker (mutex-guarded list); index and header shared read-only; alias lookup under a mutex; reader map and MatrixCache under mutexes | profile: matrix and profile byte-identical at 1 and 8 threads (sha256, ran); serve gates probed live (ran) |
| pto-peaks | 1 thread; SIMD lanes only inside the Poisson evaluation | nothing shared | scalar and NEON outputs byte-identical (ran) |
| scrna_matrix | OpenMP: norms static; tile pairs in circle-method rounds (dynamic); HNSW build and query in parallel | one heap per row; a round's tile pairs are disjoint, so no locks; write_mutex_ guards index mutation | brute force deterministic (tie-break by row); HNSW approximate and not reproducible when threaded (static, not built here) |

**Table `tB`** (columns [170, 400, 300, 370] px)

| engine | kernels | ISAs | override and evidence |
|:--|:--|:--|:--|
| fastq_stream | all_printable, count_ge, add_widen_u8_to_u32, mismatches, find_seed4 (sum_u8: tests only) | scalar, NEON, AVX2, AVX-512BW (opt-in FQ_ENABLE_AVX512) | FQ_FORCE_ISA. ran: NEON by default, scalar on request, avx2 ignored on arm64 |
| pto-peaks | Poisson -log10 p in batches of 256 bases; lanes 8/4/2/1 doubles | AVX-512F+DQ (opt-in), AVX2+FMA, NEON, scalar | PTO_PEAKS_FORCE_ISA. ran: scalar output == NEON output, byte for byte |
| scrna_matrix | gather_dot, dense_dot; HNSW distance via ScrnaCosineSpace; sparse merge-join is scalar | AVX-512F+BW (opt-in), AVX2+FMA, NEON, scalar | SCRNA_FORCE_ISA; Python simd_isa(). static + CI sweep only |
| genomic_toolkit<br>cuttag_profiler | none: no intrinsics or CPU probes in first-party sources (grep) | compiler auto-vectorisation only | -march=native is OFF by default in both |

**Table `tC`** (columns [170, 800, 822] px)

| engine | resident bound | layout and I/O primitive |
|:--|:--|:--|
| fastq_stream | (lanes x 15 slots + in-flight) x <= 384 KiB; not a function of input size (47 blocks = 15.8 MB, ran) | Buffer 64 B aligned; ring cursors 128 B apart; read(2) into a 1 MiB stage with posix_fadvise(SEQUENTIAL) |
| genomic_toolkit | histogram 2,001 counters; dup marker O(deepest pile-up) or O(distinct) with --unsorted; batch 64 KiB | Fragment 16 B; alignas(64) arrays; LineReader 1 MiB growing to a 64 MiB cap; htslib for BAM/CRAM |
| cuttag_profiler | rows x bins x 8 B allocated up front; one tally per worker; matrix cache 512 MiB / 4 entries (server) | per-reducer aligned accumulators; indexed queries through htslib; server parses regions from the fd |
| pto-peaks | 2 rings x 128 KiB; candidates 32 B each + ~12 B BH scratch; not a function of contig length | 4 MiB block reads, 1 MiB line cap; BGZF decoded natively in place, no per-record allocation |
| scrna_matrix | CSR + KnnGraph n*k*8 B; per-row heaps O(n*k); HNSW ~ n*dim*4 B + links (hnsw_index_bytes guards it) | AlignedAllocator 64 B; NumPy buffers adopted in place if aligned and immutable; no file I/O |
| all engines | memory-mapping: NONE. No mmap, madvise or MAP_* in first-party code (grep) | the only page-cache hint is one posix_fadvise(SEQUENTIAL), fastq_stream reader.cpp:76 |

**Edges (directed)**

| id | from | exit | to | entry | kind | label | via |
|:--|:--|:-:|:--|:-:|:--|:--|:--|
| e1 | `m0` | b | `m1` | t | control |  |  |
| e2 | `m1` | b | `m2` | t | control |  |  |
| e3 | `m2` | b | `m3` | t | control |  |  |
