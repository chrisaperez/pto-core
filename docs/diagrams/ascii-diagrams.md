# pto-core: terminal-grade flowcharts

Five diagrams, drawn in box-drawing characters only (`┌ ─ ┐ │ └ ┘ ├ ┤ ┬ ┴ ┼`, the double-line
set `═ ║ ╔ ╗ ╚ ╝ ╠ ╣ ╦ ╩ ╬`, the dashed `┄ ┆`, and the arrowheads `► ◄ ▲ ▼`). Every line is at
most 100 columns: it fits a wide terminal, a monospace code block on a 16:9 slide, or a
monospace text frame in Canva, draw.io or Excalidraw.

**Conventions**

| Mark | Meaning |
|:--|:--|
| solid box | a component that exists in this repository |
| dashed box (`┄ ┆`) | outside this repository: its contract is defined here, its code is not |
| `──►` light arrow | control: a call, argv, an exit code |
| `══►` double arrow | payload: a buffer handle, a stream of bytes, a boundary crossing |
| `(c1)`..`(c4)` | the four places `fastq_stream` copies payload bytes (Diagram 2) |
| `(ran)` / `(static)` | a claim executed in this review, or one read from source / CI config only |

**How these were checked.** Every number, constant and step was read from the source at commit
`0ba20d1` (not from the READMEs), and the behavioural claims were then run against binaries
built fresh from that commit. The ledger of places where a document disagreed with the code is
in [`README.md`](README.md); nothing on a diagram is taken from a document that the code
contradicts. File and line references are given so any claim can be re-checked.

Regenerate with `python3 docs/diagrams/src/make_ascii.py docs/diagrams/ascii-diagrams.md`. The
diagrams are drawn by a small canvas (`src/canvas.py`), so borders and junctions cannot drift when
wording changes; edit the text in `src/d1.py` .. `src/d5.py`, not the glyphs in this file.


## Diagram 1: system overview and tool boundaries

```text
DIAGRAM 1 / 5   pto-core: system overview and tool boundaries
solid = component in this repo    dashed = outside it: a caller, or code whose contract is set here
══╦══ the process boundary: argv + files + stdio + exit code, no shared memory, no IPC
CALLERS
┌─ pto-core ──┐     ┌─ import pto ┐     ┌┄ PATH ┄┄┄┄┄┄┐      ┌┄ browser ┄┄┄┐      ┌┄ Python ┄┄┄┄┐
│console cmd  │     │Python API   │     ┆SLURM, Make, ┆      ┆dashboard UI ┆      ┆import       ┆
│pto.cli:main │     │7 functions  │     ┆bare names   ┆      ┆loopback only┆      ┆scrna_matrix ┆
│             │     │             │     ┆             ┆      ┆             ┆      ┆             ┆
└──────┬──────┘     └──────┬──────┘     └┄┄┄┄┄┄┬┄┄┄┄┄┄┘      └┄┄┄┄┄┄┬┄┄┄┄┄┄┘      └┄┄┄┄┄┄┬┄┄┄┄┄┄┘
       │                   │                   ┆                    ┆                    ┆
       ▼                   ▼                   ┆                    ┆                    ┆
┌─ packaging/pto (stdlib only) ┐               ┆                    ┆                    ┆
│ cli.py > api.py > _run.py    │               ┆                    ┆                    ┆
│ resolve: $PTO_BIN_DIR >      │               ┆                    ┆                    ┆
│   pto/_vendor/bin > $PATH    │               ┆                    ┆                    ┆
│ argv list, shell=False       │       bare tool names       HTTP 127.0.0.1           in-proc
│ threads = min(affinity,      │               ┆                    ┆                    ┆
│   cgroup quota, cpu_count)   │               ┆                    ┆                    ┆
│ wait4: wall + peak RSS       │               ┆                    ┆                    ┆
│ rc 2 -> InvalidInput         │               ┆                    ┆                    ┆
└───────────────┬──────────────┘               ┆                    ┆                    ┆
                │                              ┆                    ┆                    ┆
                │                              ┆                    ┆                    ┆
PROCESS         ▼                              ▼                    ┆                    ┆
BOUNDARY ╔═══════════════════╦═══════════════════╦═════════════╗    ┆                    ┆
         ║                   ║                   ║             ║    ┆                    ┆
┌────────▼────────┐ ┌────────▼────────┐ ┌────────▼────────┐ ┌──▼────▼─────────┐ ┌────────▼────────┐
│fastq_stream     │ │genomic_toolkit  │ │pto-peaks        │ │cuttag_profiler  │ │scrna_matrix     │
│C++20 CLI  0.1.0 │ │C++20 CLI  0.1.0 │ │C++20 CLI  0.1.0 │ │C++20 CLI  1.0.0 │ │C++20 header-only│
│in : FASTQ raw/gz│ │sizes, markdup,  │ │in : frag BED or │ │profile | serve  │ │+ pybind11 module│
│     BGZF, BAM(c)│ │frip             │ │     sorted BAM, │ │in : indexed BAM+│ │in : NumPy/SciPy │
│out: QC JSON+text│ │in : BAM/CRAM,   │ │     (stdin ok)  │ │     BED or GTF  │ │     CSR buffers │
│     FASTQ: file,│ │     BEDPE, BED  │ │out: narrowPeak  │ │out: matrix.tsv, │ │out: KnnGraph as │
│     FIFO, stdout│ │out: JSON stdout │ │BAM: own decoder │ │     profile.tsv │ │     NumPy views │
│threads 3+D+W    │ │     (+ hist TSV)│ │threads 1        │ │     or HTTP (d) │ │threads OpenMP   │
│SIMD runtime ISA │ │threads 1+OpenMP │ │SIMD runtime ISA │ │threads N workers│ │SIMD runtime ISA │
│rc 0 | 1 | 2 (a) │ │no hand SIMD     │ │  (Poisson only) │ │no hand SIMD     │ │no CLI, no rc (b)│
│                 │ │rc 0 | 1 | 2     │ │rc 0 | 1 | 2     │ │rc 0 | 1 | 2     │ │                 │
└─────────────────┘ └─────────────────┘ └─────────────────┘ └─────────────────┘ └─────────────────┘

LINKS
┌┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┐ ┌┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┐ ┌┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┐ ┌┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┐ ┌┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┐
┆zlib + libdeflate┆ ┆htslib (AUTO)    ┆ ┆libdeflate | zlib┆ ┆htslib (required)┆ ┆OpenMP (required)┆
┆htslib (opt) (c) ┆ ┆OpenMP (opt)     ┆ ┆no htslib        ┆ ┆httplib+json (d) ┆ ┆hnswlib (vendor) ┆
└┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┘ └┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┘ └┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┘ └┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┘ └┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┘


┌─ build and distribution ────────────────────────────────────────────────────────────────────────┐
│ cmake -S .  ->  PTO_BUILD_{SCRNA_MATRIX,FASTQ_STREAM,CUTTAG_PROFILER,GENOMIC_TOOLKIT,PEAKS}=ON  │
│ wheel pto-core     py3-none-<platform>: 4 CLIs in pto/_vendor/bin (scrna_matrix excluded)       │
│ wheel scrna-matrix cp39..cp313:        scrna_matrix_py extension (OpenMP hard-fails configure)  │
│ no shared layer: 4 CLIs each carry a private file_identity.hpp; no common headers or runtime    │
└─────────────────────────────────────────────────────────────────────────────────────────────────┘

(a) fastq_stream also exits 2 on CLI usage errors; the other CLIs use 1 there and reserve 2
    for 'the input's shape made the answer wrong' (naming mismatch, unsorted, all regions skipped)
(b) no wrapper in pto-core: installed with `pip install scrna-matrix`, imported in-process
(c) FQ_ENABLE_HTS defaults OFF in the module, ON in the wheel   (d) compiled out by PTO_CLOUD_BUILD
```

**Verified against**

- `packaging/pto/_binaries.py:44-84` (TOOLS, search order), `_run.py:119-200`, `_cpu.py:75-83`, `api.py:53-55`, `cli.py:229-385`
- `pyproject.toml:87-88,119,131,186-249` (console script, ABI tag, install dir, `PTO_BUILD_*`, `*_NATIVE` OFF)
- `CMakeLists.txt:28-58` (module switches, `PTO_CLOUD_BUILD`)
- exit codes: `fastq_stream/src/main.cpp:82-85,338-351`, `genomic_toolkit/src/main.cpp:253-423`, `cuttag_profiler/src/main.cpp:322,421-447`, `peaks/src/main.cpp:63-65`
- no shared layer: four `modules/*/include/*/file_identity.hpp` copies; SIMD only in `fastq_stream/.../simd.hpp`, `peaks/.../poisson_model.hpp`, `scrna_matrix/.../simd_math.hpp` (grep)
- ran: `--version` of the four binaries (0.1.0, 0.1.0, 1.0.0, 0.1.0); layout of the built `dist/pto_core-0.1.6-*.whl` files

## Diagram 2: streaming data flow and buffer management

```text
DIAGRAM 2 / 5   fastq_stream: streaming data flow and buffer management
══► a Buffer HANDLE moves by pointer through a lock-free SPSC ring; payload bytes are not copied
(c1)-(c4): the only four places payload bytes ARE copied.   T=--threads  D=inflaters  W=workers

┌─ SOURCE ────────────────────────────────────────────────────────────────────────────────────────┐
│ path | "-" (stdin) | FIFO       BAM/CRAM/SAM recognised by content (htslib), regular files only │
│ framing, first 18 bytes: 1f 8b + BGZF 'BC' subfield = BGZF; 1f 8b alone = gzip; else raw        │
└─────────────────────────────┬───────────────────────────────────────────────────────────────────┘
                              ▼
┌─ READER  1 thread  pipeline.cpp:364  reader.cpp:169-303 ────────────────────────────────────────┐
│ FdSource: 1 MiB stage (doubles on demand), posix_fadvise(SEQUENTIAL), read(2), EINTR retried    │
│ BGZF: memcpy up to 4 whole compressed members into one pooled Buffer; lane = k mod D      (c1)  │
│ gzip: zlib inflate (15+16 bits, inflateReset per member) into a pooled Buffer; serial           │
│ raw : memcpy stage -> Buffer                                                              (c1)  │
│ BAM : htslib sam_read1, decode pool of max(1, spare/3) threads, records rewritten as FASTQ      │
│ EOF: BGZF must end with its 28-byte empty block; a missing one, or trailing bytes, is an error  │
└─────────────────────────────┬──────────────────────────────────────────╥────────────────────────┘
                              │                                          ║
                              ▼                                          ║
┌─ lane x D = SpscRing<Buffer,16> ────────────────────────────┐          ║ direct lane (1 ring)
│ 15 usable slots + closed flag; BGZF only                    │          ║ used by gzip, raw, BAM:
└─────────────────────────────╥───────────────────────────────┘          ║ no inflate rank, so the
                              ▼                                          ║ reader feeds assembler
┌─ INFLATERS  D threads  pipeline.cpp:385 ────────────────────┐          ║
│ libdeflate_gzip_decompress, one decompressor per thread     │          ║
│ 4 members -> one 256 KiB Buffer from in_pool                │          ║
│ the compressed Buffer returns to in_pool                    │          ║
│ D = max(1, spare/3), spare = max(1, T-3); BGZF only         │          ║
└─────────────────────────────╥───────────────────────────────┘          ║
                              ▼                                          ║
┌─ lane x D = SpscRing<Buffer,16> ────────────────────────────┐          ║
│ pop order = push order: chunk k on lane k mod D             │          ║
└─────────────────────────────╥───────────────────────────────┘          ║
                              ▼                                          ▼
┌─ ASSEMBLER  1 thread  pipeline.cpp:409  record.hpp:36 ──────────────────────────────────────────┐
│ pop lane[k mod D] (or the direct lane) in the reader's order: stream order, no reorder buffer   │
│ carry: extend_front(carry) into the 64 KiB prefix ahead of data(), an in-place splice     (c2)  │
│ cut = last_record_boundary(): memchr '\n', every 4th newline ends a record ('@' can be QUAL)    │
│ whole records -> work_in[j mod W]; trailing partial record -> carry (over 64 KiB = error)       │
│ EOF: a last record without its final '\n' is accepted; any other short tail is an error         │
│ one thread on purpose: memchr runs at memory bandwidth, well above the inflate rate             │
└─────────────────────────────╥───────────────────────────────────────────────────────────────────┘
                              ▼
┌─ lane x W = SpscRing<Buffer,16> ────────────────────────────┐
│ chunk j on lane j mod W; input of worker i                  │
└─────────────────────────────╥───────────────────────────────┘
                              ▼
┌─ WORKERS  W threads  pipeline.cpp:488  process() :152 ──────────────────────────────────────────┐
│ W = max(1, spare-D). out = out_pool.acquire(): cap 320 KiB = 256 KiB chunk + 64 KiB carry       │
│ per record, in order: header/'+' printable -> Phred+64 rebase (only --phred-offset 64)          │
│   -> accumulate_read on the RAW read -> trim_front -> trim_back -> sliding window               │
│   -> adapter seed-and-extend (find_seed4 + mismatches) -> min_len, max_n_rate, min_mean_q       │
│   -> surviving records memcpy'd into out                                                  (c3)  │
│ SIMD (runtime ISA): all_printable, count_ge x2, add_widen_u8_to_u32, mismatches, find_seed4     │
│ per-cycle sums: 32-bit CycleScratch folded into 64-bit QcStats once per chunk; merged at join   │
│ ALWAYS pushes exactly ONE out Buffer per input chunk, even if empty: writer order needs it      │
│ --qc-only: identical work, serialisation skipped, an empty out Buffer is still pushed           │
└─────────────────────────────╥───────────────────────────────────────────────────────────────────┘
                              ▼
┌─ lane x W = SpscRing<Buffer,16> ────────────────────────────┐
│ worker i's output: one Buffer per input chunk               │
└─────────────────────────────╥───────────────────────────────┘
                              ▼
┌─ WRITER  1 thread  pipeline.cpp:566 ────────────────────────────────────────────────────────────┐
│ pop work_out[k mod W] round-robin = original order; write(2) loop, EINTR retried          (c4)  │
│ SIGPIPE ignored: a closed FIFO becomes EPIPE -> error gate; close() checked (NFS/FUSE ENOSPC)   │
└─────────────────────────────┬───────────────────────────────────────────────────────────────────┘
                              ▼
┌─ SINK ──────────────────────────────────────────────────────────────────────────────────────────┐
│ stdout | file (O_TRUNC) | FIFO: --mkfifo, open(2) BLOCKS until a reader attaches (verified)     │
│ refused before anything is opened: -i X -o X, also through hard links and symlinks              │
└─────────────────────────────────────────────────────────────────────────────────────────────────┘

BUFFER ANATOMY   buffer.hpp:28-35   BufferPool uses aligned_alloc(64 B)
┌──────────────────────────────────────┬──────────────────────────────────────────────────────────┐
│ kPrefix 64 KiB: headroom             │ payload 256 KiB (in_pool) / 320 KiB (out_pool)           │
│ extend_front() splices carry here    │ data() .... size() .............. capacity()             │
└──────────────────────────────────────┴──────────────────────────────────────────────────────────┘
handle = {pool*, base*, cap, size, front_} + seq; seq is written by reader/inflater, never read
pool free list: std::mutex (buffer.hpp:118), one lock pair per 256 KiB of stream, off the data path

SPSC RING  lockfree_queue.hpp:47  kCacheLine = 128 B on arm64 and x86-64, 64 B elsewhere
┌─ one producer, one consumer, no CAS ────────────────────────────────────────────────────────────┐
│ head_, cached_tail_, tail_, cached_head_, slots_[16]: each alignas(128); no line is shared      │
│ try_push: relaxed load of tail; reload head (acquire) only when full; release-store tail        │
│ try_pop : relaxed load of head; reload tail (acquire) only when empty; release-store head       │
│ Backoff pipeline.hpp:72: 128 spins, 16 yields, then 100 us sleeps (yield-only burned ~11 cores) │
└─────────────────────────────────────────────────────────────────────────────────────────────────┘

CONTROL AND BOUNDS
┌─ error gate, termination, memory ───────────────────────────────────────────────────────────────┐
│ ErrorGate (pipeline.cpp:26): first error wins; abort flag polled in every ring wait; then join  │
│ close: the producer release-stores closed; the consumer exits when the ring is empty AND closed │
│ a thread that cannot start trips the gate and joins the rest: no terminate (pipeline.cpp:599)   │
│ resident bound: (lanes x 15 + in-flight) x <= 384 KiB per block; not a function of input size   │
│ measured: 47 blocks = 15.8 MB pool for a 12.8 MB FASTQ at -t 6 (fastq_stream 0.1.0, this build) │
└─────────────────────────────────────────────────────────────────────────────────────────────────┘

(c1) stage -> Buffer: raw, BGZF compressed batch; plain gzip inflates into the Buffer directly
(c2) carry splice: only the straddling partial record, at most 64 KiB
(c3) worker serialises survivors into the out Buffer      (c4) writer's write(2): the kernel's copy
```

**Verified against**

- `fastq_stream/src/pipeline.cpp:26` (ErrorGate), `:152` (process), `:263-345` (sizing, source choice), `:364-613` (the five stage lambdas and thread start)
- `include/fastq_stream/pipeline.hpp:46` (ring slots), `:72-98` (Backoff, Lane), `:129-172` (Config)
- `include/fastq_stream/buffer.hpp:28-35,103-153` (chunk, prefix, alignment, pool, mutex), `lockfree_queue.hpp:19-104`
- `src/reader.cpp:32-303,350-380` (stage, sniffing, BGZF/gzip/raw), `include/fastq_stream/record.hpp:36-54`, `phred_calculator.hpp:186-246`
- `src/hts_input.cpp:14-137`, `include/fastq_stream/file_identity.hpp:63-66` (device+inode identity)
- ran: 14 thread configurations (`-t 1..12`, plain and BGZF), 0 mismatches against the formula; framing detection; FIFO blocks until a reader attaches and delivers 239,984 lines; hard link and symlink outputs refused, input untouched

## Diagram 3: the four command-line engines

```text
DIAGRAM 3 / 5   the four command-line engines: control path and core loop
each panel: top = main() control flow; bottom = the loop that touches every record or region

┌─ 1  fastq_stream   qc | trim   src/main.cpp  src/pipeline.cpp ──────────────────────────────────┐
│ ┌─ argv ────────┐  ┌─ validate ────┐  ┌─ probe ───────┐  ┌─ run_pipeline ┐  ┌─ report + rc ─┐   │
│ │-i -o -j -t    │  │range checks;  │  │BAM/CRAM by    │  │see Diagram 2: │  │text, JSON (-j)│   │
│ │--qc-only ...  ├─►│usage err: rc 2├─►│content (files)├─►│3+D+W threads  ├─►│rc 0 | 1 | 2   │   │
│ │               │  │               │  │               │  │               │  │               │   │
│ └───────────────┘  └───────────────┘  └───────────────┘  └───────────────┘  └───────────────┘   │
│                                                                                                 │
│ per record, in each of W workers (RecordProcessor::process, pipeline.cpp:152):                  │
│ ┌─ checks ──────┐  ┌─ QC on RAW ───┐  ┌─ trim ────────┐  ┌─ adapters ────┐  ┌─ filters ─────┐   │
│ │header/SEQ/QUAL│  │counts, Q20/Q30│  │front, tail,   │  │seed + extend  │  │min_len, max_n,│   │
│ │printable;     ├─►│GC, N, cycles  ├─►│sliding window ├─►│TruSeq, Nextera├─►│min_mean_q     │   │
│ │Phred+64 rebase│  │               │  │               │  │               │  │               │   │
│ └───────────────┘  └───────────────┘  └───────────────┘  └───────────────┘  └───────────────┘   │
│ rc 2 also when qualities fit only Phred+64 and no --phred-offset was given (report first)       │
└─────────────────────────────────────────────────────────────────────────────────────────────────┘

┌─ 2  genomic_toolkit   sizes | markdup | frip   src/main.cpp  src/fragment_stream.cpp ───────────┐
│ ┌─ argv ─────────────┐  ┌─ peaks [frip] ─────┐  ┌─ PeakSet::build ───┐  ┌─ bind ─────────────┐  │
│ │sizes|markdup|frip  │  │read_peaks_bed      │  │sort; merge overlap │  │pre-seed contigs    │  │
│ │checked int parse   ├─►│malformed = fatal   ├─►│+ book-ended; SoA   ├─►│peak tid -> stream  │  │
│ │                    │  │                    │  │                    │  │                    │  │
│ └────────────────────┘  └────────────────────┘  └────────────────────┘  └────────────────────┘  │
│                                                                                                 │
│ stream: one pass, no sort, no temp file                                                         │
│ ┌─ stream_fragments ─┐  ┌─ reader ───────────┐  ┌─ BatchEmitter ─────┐  ┌─ consumers ────────┐  │
│ │.bam .cram .sam:    │  │BAM: TLEN>0 = 1 pair│  │Fragment = 16 B     │  │sizes: histogram    │  │
│ │htslib; else .bedpe;│  │text: 1 MiB reader  │  │4096/batch = 64 KiB │  │markdup: marker     │  │
│ │else fragment BED   ├─►│cap 64 MiB / line   ├─►│size histogram here ├─►│frip: marker+count  │  │
│ │                    │  │                    │  │                    │  │                    │  │
│ └────────────────────┘  └────────────────────┘  └────────────────────┘  └────────────────────┘  │
│                                                                                                 │
│ ┌─ DuplicateMarker ──────────────────────────┐  ┌─ FripCounter.add_batch ─────────────────────┐ │
│ │sorted: keys (end,strand) of ONE start      │  │omp for, if batch >= 512; per-thread         │ │
│ │hash (--unsorted): 64-bit fingerprints      │  │OverlapCursor; 5 integer sums; checked_add   │ │
│ │sets kDuplicate in place on the batch       │  │                                             │ │
│ └────────────────────────────────────────────┘  └─────────────────────────────────────────────┘ │
│ exit 2: nothing parsed | out of order (markdup; frip w/o --keep-dups) | no overlap + bad naming │
└─────────────────────────────────────────────────────────────────────────────────────────────────┘

┌─ 3  cuttag_profiler   profile   src/main.cpp  src/signal_calc.cpp  src/bam_reader.cpp ──────────┐
│ ┌─ guards ──────┐  ┌─ regions ─────┐  ┌─ BamReader ───┐  ┌─ compute ─────┐  ┌─ finish ──────┐   │
│ │output != input│  │BED/GTF (.gz); │  │header, index; │  │compute_matrix:│  │scale, means,  │   │
│ │limits checked ├─►│names sanitised├─►│mapped total:  ├─►│N<=256, <=rows;├─►│write TSV;     │   │
│ │               │  │               │  │index metadata │  │blocks of rows │  │all skip: rc 2 │   │
│ └───────────────┘  └───────────────┘  └───────────────┘  └───────────────┘  └───────────────┘   │
│                                                                                                 │
│ per region, in each of N workers (signal_calc.cpp:242); one leased htsFile each, no BGZF pool:  │
│ ┌─ window ──────┐  ┌─ skip? ───────┐  ┌─ query_bins ──┐  ┌─ alignment ───┐  ┌─ row ─────────┐   │
│ │TSS/center/TES │  │no such contig,│  │widen by max   │  │flags + MAPQ;  │  │fill bins;     │   │
│ │anchor; minus  │  │or off-contig: │  │frag length;   │  │leftmost mate  │  │minus reversed │   │
│ │strand mirrored├─►│row = missing  ├─►│sam_itr_queryi ├─►│= TLEN span;   ├─►│= 5' to 3'     │   │
│ │               │  │               │  │               │  │reads, bases   │  │               │   │
│ └───────────────┘  └───────────────┘  └───────────────┘  └───────────────┘  └───────────────┘   │
│                                                                                                 │
│ matrix = regions x bins doubles, allocated up front (bins <= 10,000,000, window <= ~922 Gbp)    │
│ column means: min(N, rows/64, 64 MiB/(cols*8)) reducers, cache-line-aligned accumulators        │
│ skipped rows stay in the matrix, named, but are excluded from means, BPM total and scaling      │
│                                                                                                 │
└─────────────────────────────────────────────────────────────────────────────────────────────────┘

┌─ 4  pto-peaks   main.cpp + sliding_window, caller, poisson_model, bam_streamer (.hpp) ──────────┐
│ ┌─ detect ───────────┐  ┌─ records ──────────┐  ┌─ Driver.add ───────┐  ┌─ finalize ─────────┐  │
│ │peek first byte:    │  │BAM: BGZF, CRC32,   │  │order, bounds check;│  │BH over ALL         │  │
│ │1f => BAM, else BED ├─►│leftmost mate TLEN>0├─►│seek across empty   ├─►│candidates, then    │  │
│ │(stdin ok)          │  │BED: 4 MiB blocks   │  │coverage deserts    │  │narrowPeak          │  │
│ └────────────────────┘  └────────────────────┘  └────────────────────┘  └────────────────────┘  │
│                                                                                                 │
│ per base, in one thread, memory independent of contig length (sliding_window.hpp, caller.hpp):  │
│ ┌─ MultiScaleWindow ─┐  ┌─ 3 backgrounds ────┐  ┌─ PeakCaller ───────┐  ┌─ state machine ────┐  │
│ │delta+pileup rings, │  │centred 1k, 5k, 10k │  │prefilter, queue of │  │IDLE IN_PEAK        │  │
│ │2 x 16384 x int64   │  │running sums; lag = │  │256, SIMD Poisson   │  │SUMMIT_SEARCH       │  │
│ │= 256 KiB           ├─►│5000 bases          ├─►│-log10 p            ├─►│CLOSING (max_gap)   │  │
│ │                    │  │                    │  │                    │  │                    │  │
│ └────────────────────┘  └────────────────────┘  └────────────────────┘  └────────────────────┘  │
│                                                                                                 │
│ lambda = max(genome floor, 3 local rates); a degenerate lambda scores p = 1 (fails closed)      │
│ peaks are held (32 B each) until finalize(): BH needs every p-value, so only the INPUT streams  │
│ exit 2: unsorted, contig not in sizes, bad field, BED fragment > 4096 bp (BAM: skipped)         │
│                                                                                                 │
│                                                                                                 │
└─────────────────────────────────────────────────────────────────────────────────────────────────┘
```

**Verified against**

- `fastq_stream`: `src/main.cpp`, `src/pipeline.cpp:152-253`
- `genomic_toolkit`: `src/main.cpp:253-423`, `src/fragment_stream.cpp:62-85,271-501`, `include/toolkit/types.hpp:58-84`, `dup_marker.hpp:156-212`, `frip.hpp:100-166`, `peak_set.hpp:154-235`, `src/line_reader.hpp:153-157`
- `cuttag_profiler`: `src/main.cpp:116-325`, `src/signal_calc.cpp:98-451`, `src/bam_reader.cpp:123-372`, `include/profiler/types.hpp:139-262`
- `pto-peaks`: `src/main.cpp:331-430,594-764`, `sliding_window.hpp:71-388`, `caller.hpp:57-500`, `poisson_model.hpp:1784-1889`, `bam_streamer.hpp:100-110`
- ran: FRiP with `--keep-dups` = 148,680 in-peak fragments (README figure reproduced); region tiling sums to the whole; matrix and profile byte-identical at 1 and 8 threads; `pto-peaks` output identical for BAM, BED, stdin BED and stdin BAM

## Diagram 4: integration and orchestration

```text
DIAGRAM 4 / 5   integration and orchestration: how callers reach the binaries
everything below ships in this repository: no external product or proprietary service is named

A   FOUR WAYS IN
           ┌─ pto.frip(...) ──┐  ┌─ run_tool ───────┐  ┌─ Popen+wait4 ────┐  ┌─ result ─────────┐
           │7 functions;      │  │resolve binary;   │  │spawn child;      │  │rc 0: Result      │
PYTHON API │None = omitted    ├─►│argv, shell=False ├─►│os.wait4 reaps:   ├─►│rc 2: InvalidInput│
           │existence checked │  │stdio -> tmpfiles │  │rc + ru_maxrss    │  │rc n: ToolError   │
           └──────────────────┘  └──────────────────┘  └──────────────────┘  └──────────────────┘

           ┌─ pto-core frip ──┐  ┌─ cli._cmd_frip ──┐  ┌─ output ─────────┐  ┌─ exit code ──────┐
           │argparse command  │  │calls pto.api.frip│  │metrics to stdout;│  │0 ok | 1 error    │
pto-core   │--json optional   ├─►│same path as above├─►│cost to stderr    ├─►│2 InvalidInput    │
           │also: man, doctor │  │                  │  │                  │  │141 broken pipe   │
           └──────────────────┘  └──────────────────┘  └──────────────────┘  └──────────────────┘

           ┌─ exec ───────────┐  ┌─ bin-dir ────────┐  ┌─ PATH ───────────┐  ┌─ consumer ───────┐
           │pto-core exec T   │  │prints dir of the │  │PATH=$(pto-core   │  │shell, SLURM,     │
exec / PATH│os.execv replaces ├─►│first resolvable  ├─►│ bin-dir):$PATH   ├─►│Nextflow: no      │
           │the process       │  │engine            │  │bare tool names   │  │Python in the path│
           └──────────────────┘  └──────────────────┘  └──────────────────┘  └──────────────────┘

┌─ resolution and threads   packaging/pto/_binaries.py:64-84  _cpu.py:75 ─────────────────────────┐
│ $PTO_BIN_DIR  >  pto/_vendor/bin  >  $PATH (shutil.which);  none = BinaryNotFound + build hint  │
│ TOOLS = the 4 engines (a test fails if the wheel vendors a fifth binary that TOOLS lacks)       │
│ threads = min(sched_getaffinity, cgroup v2/v1 quota, cpu_count), never below 1                  │
│ passed as -t to fastq_stream and --threads to cuttag_profiler ONLY; the other two receive none  │
│                                                                                                 │
└─────────────────────────────────────────────────────────────────────────────────────────────────┘

B   THE RESULT CONTRACT   (verified by running each binary)
┌─────────────────┬────────────────────┬──────────────────┬──────────────────────┬─────────────────┐
│ engine          │ machine-readable   │ artifacts        │ exit codes           │ thread control  │
├─────────────────┼────────────────────┼──────────────────┼──────────────────────┼─────────────────┤
│ fastq_stream    │ JSON file: -j PATH │ trimmed FASTQ:   │ 0 ok  1 runtime      │ -t N: wrapper   │
│ qc | trim       │ text: --report     │ -o FILE|FIFO|-   │ 2 usage error, or    │ passes it       │
│                 │ wrapper reads JSON │                  │ Phred+64 suspected   │                 │
├─────────────────┼────────────────────┼──────────────────┼──────────────────────┼─────────────────┤
│ genomic_toolkit │ JSON on stdout     │ histogram TSV:   │ 0 ok  1 usage/system │ no flag: OpenMP │
│ sizes markdup   │ (--json flag)      │ --histogram FILE │ 2 input shape:       │ OMP_NUM_THREADS │
│ frip            │                    │                  │ naming, order, empty │ wrapper: unset  │
├─────────────────┼────────────────────┼──────────────────┼──────────────────────┼─────────────────┤
│ cuttag_profiler │ none;              │ matrix.tsv,      │ 0 ok  1 error        │ --threads N     │
│ profile | serve │ summary on stderr  │ profile.tsv      │ 2 every region       │ wrapper passes  │
│                 │                    │ (or HTTP: serve) │ skipped              │                 │
├─────────────────┼────────────────────┼──────────────────┼──────────────────────┼─────────────────┤
│ pto-peaks       │ none;              │ narrowPeak       │ 0 ok  1 usage/IO     │ one thread;     │
│                 │ counts on stderr   │ --out FILE (or   │ 2 unsorted, contig   │ no thread flag  │
│                 │                    │ stdout)          │ missing, bad field   │                 │
└─────────────────┴────────────────────┴──────────────────┴──────────────────────┴─────────────────┘
fastq_stream exits 2 for usage errors too: pto.InvalidInput fires on a bad parameter there

C   BUILD, PACKAGE, PUBLISH
┌─ cmake -S . ────┐ ┌─ scikit-build ──┐ ┌─ cibuildwheel ──┐ ┌─ repair ────────┐ ┌─ publish ───────┐
│PTO_BUILD_* x5   │ │install-dir =    │ │cp312 pin;       │ │delocate (mac)   │ │v* tag -> PyPI   │
│PTO_CLOUD_BUILD  ├►│pto/_vendor;     ├►│py3-none-<plat>  ├►│auditwheel(lnx)  ├►│OIDC, no token   │
│default OFF      │ │no headers       │ │doctor + pytest  │ │vendors libs     │ │setuptools_scm   │
└─────────────────┘ └─────────────────┘ └─────────────────┘ └─────────────────┘ └─────────────────┘
vendored libs: macOS pto/.dylibs (omp hts deflate z lzma); Linux pto_core.libs (gomp bz2 lzma)
ci.yml jobs: build-and-test (ISA + OMP sweeps), sanitizers, config-matrix, static-link, TSan

D   THE ONE NETWORK PATH: cuttag_profiler serve (http_server.cpp; absent in cloud builds)
╔═ TRUST BOUNDARY   127.0.0.1 only; --allow-remote must be explicit; reach it by SSH tunnel ══════╗
║ one handler, /api/* only, in this order (http_server.cpp:430-494): Host before token            ║
║ ┌─ client ──────────┐   ┌─ Host  421 ───────┐   ┌─ Origin  403 ─────┐   ┌─ token  401 ──────┐   ║
║ │browser or curl;   │   │loopback name and  │   │absent is fine;    │   │?t= or Bearer;     │   ║
║ │URL has ?t=TOKEN   │   │bound port only:   │   │present must be    │   │constant_time_     │   ║
║ │(printed once)     ├──►│defeats rebinding  ├──►│this server        ├──►│equals, 256 bit    │   ║
║ │                   │   │                   │   │                   │   │                   │   ║
║ └───────────────────┘   └───────────────────┘   └───────────────────┘   └─────────┬─────────┘   ║
║           ┌───────────────────────────────────────────────────────────────────────┘             ║
║           ▼                                                                                     ║
║ ┌─ body gates ──────┐   ┌─ resolve_input ───┐   ┌─ compute ─────────┐   ┌─ cache + fetch ───┐   ║
║ │<= 1 MiB: 413      │   │under data root;   │   │compute_matrix as  │   │MatrixCache 512MiB,│   ║
║ │JSON type: 415     │   │open once, on the  │   │in Diagram 3;      │   │4 entries, FIFO;   │   ║
║ │30 s timeouts      ├──►│fd: no re-lookup   ├──►│HTTP pool = 16     ├──►│GET /api/matrix    │   ║
║ │                   │   │                   │   │                   │   │                   │   ║
║ └───────────────────┘   └───────────────────┘   └───────────────────┘   └───────────────────┘   ║
║ GET /api/matrix?token=M needs BOTH the session token and the matrix token M (256 bit)           ║
║ confinement: canonical path under --data-root (else 400); FIFO, symlink, device refused         ║
║ assets from rodata need no token (CSP default-src 'self'); API replies: no-store, frame DENY    ║
║ container: 127.0.0.1:8080 only; read_only, cap_drop ALL, no-new-privileges, 4 GiB, 2 CPUs       ║
╚═════════════════════════════════════════════════════════════════════════════════════════════════╝

E   THE serve SUBCOMMAND: REACHABILITY AND THE CLOUD-BUILD GATE
pto.profile and `pto-core profile` never build it, yet the wheel's binary still contains it,
so `pto-core exec cuttag_profiler serve --data-root D` reaches it: no guard in _cmd_exec (ran it)
PTO_CLOUD_BUILD=ON drops the server entirely; `serve` is refused by name (exit 1); nm-based test
```

**Verified against**

- `packaging/pto/_binaries.py:64-84`, `_run.py:119-200`, `_cpu.py:75-83`, `cli.py:83-92,229-385`, `errors.py`
- `pyproject.toml`, `.github/workflows/ci.yml:49-546`, `.github/workflows/wheels.yml`, built wheels in `dist/`
- `cuttag_profiler/src/http_server.cpp:208-235,277-303,401-803`, `include/profiler/matrix_cache.hpp:52-55`, `CMakeLists.txt:269,383-405`
- `deploy/Dockerfile`, `deploy/docker-compose.yml`, `deploy/entrypoint.sh`
- ran: every HTTP status on the diagram (401, 421, 403, 415, 413, 400, 404, 200) against a live loopback server; the Python wrapper's argv, exit-code mapping and `exec` passthrough; `strings` on the shipped wheel binary shows the `serve` code is present

## Diagram 5: concurrency, ISA dispatch and memory

```text
DIAGRAM 5 / 5   concurrency, ISA dispatch and memory
evidence tags: (ran) = executed in this review;  (static) = read from source or CI config only

A   CONCURRENCY MODEL PER ENGINE
┌─────────────────┬───────────────────────────┬───────────────────────────┬────────────────────────┐
│ engine          │ threads and primitive     │ shared state and sync     │ ordering / determinism │
├─────────────────┼───────────────────────────┼───────────────────────────┼────────────────────────┤
│ fastq_stream    │ 3+D+W std::thread:        │ SPSC rings: acq/rel       │ lane order = stream    │
│                 │ reader, D inflaters,      │ atomics, no CAS; pools:   │ order. raw, gzip, BGZF │
│                 │ assembler, W workers,     │ std::mutex; ErrorGate:    │ and NEON vs scalar:    │
│                 │ writer; BAM adds an       │ mutex + abort atomic      │ 1 identical result     │
│                 │ htslib decode pool        │                           │ (ran)                  │
├─────────────────┼───────────────────────────┼───────────────────────────┼────────────────────────┤
│ genomic_toolkit │ 1 streaming thread;       │ PeakSet immutable;        │ identical JSON at      │
│                 │ OpenMP parallel for in    │ per-thread cursor;        │ OMP_NUM_THREADS=1      │
│                 │ FRiP if batch >= 512,     │ integer reductions;       │ and 8; region tiles    │
│                 │ schedule(static)          │ checked_add outside       │ sum to the whole       │
│                 │                           │ the parallel region       │ (ran)                  │
├─────────────────┼───────────────────────────┼───────────────────────────┼────────────────────────┤
│ cuttag_profiler │ N std::thread workers,    │ one leased htsFile per    │ matrix and profile     │
│ profile         │ block-partitioned rows;   │ worker (mutex list);      │ byte-identical at 1    │
│                 │ then reducers =           │ index + header shared     │ and 8 threads          │
│                 │ min(N, rows/64, mem)      │ read-only; alias lookup   │ (sha256, ran)          │
│                 │                           │ under a mutex             │                        │
├─────────────────┼───────────────────────────┼───────────────────────────┼────────────────────────┤
│ cuttag_profiler │ httplib pool: 16          │ reader map: mutex;        │ not benchmarked;       │
│ serve           │ threads (CMake sets       │ MatrixCache: mutex,       │ gates unit-tested      │
│                 │ THREAD_POOL_COUNT)        │ FIFO by bytes             │ and probed live (ran)  │
│                 │                           │                           │                        │
│                 │                           │                           │                        │
├─────────────────┼───────────────────────────┼───────────────────────────┼────────────────────────┤
│ pto-peaks       │ 1 thread; SIMD lanes      │ nothing shared            │ scalar and NEON        │
│                 │ only inside the           │                           │ outputs byte-          │
│                 │ Poisson evaluation        │                           │ identical (ran)        │
│                 │                           │                           │                        │
│                 │                           │                           │                        │
├─────────────────┼───────────────────────────┼───────────────────────────┼────────────────────────┤
│ scrna_matrix    │ OpenMP: norms static;     │ one heap per row; a       │ brute force: tie-break │
│                 │ tile pairs in circle-     │ round's tile pairs are    │ by row, deterministic; │
│                 │ method rounds, dynamic;   │ disjoint (no locks);      │ HNSW approximate, not  │
│                 │ HNSW build and query      │ write_mutex_ guards       │ reproducible threaded  │
│                 │                           │ index mutation            │ (static only)          │
└─────────────────┴───────────────────────────┴───────────────────────────┴────────────────────────┘

B   ISA DISPATCH   one mechanism, three engines
per-function target attributes, never a target-wide -mavx2; level resolved once (function static)
  order: FORCE_ISA env (if compiled in) > __builtin_cpu_supports (checks OS register save too)
         > NEON on AArch64 (no probe) > scalar.  wheels build with every *_NATIVE flag OFF

┌─────────────────┬────────────────────────────┬────────────────────────┬──────────────────────────┐
│ engine          │ kernels                    │ ISAs                   │ override and evidence    │
├─────────────────┼────────────────────────────┼────────────────────────┼──────────────────────────┤
│ fastq_stream    │ all_printable, count_ge,   │ scalar, NEON, AVX2,    │ FQ_FORCE_ISA. ran:       │
│                 │ add_widen_u8_to_u32,       │ AVX-512BW (opt-in:     │ NEON default, scalar     │
│                 │ mismatches, find_seed4     │ FQ_ENABLE_AVX512)      │ on request; avx2 is      │
│                 │ (sum_u8: tests only)       │                        │ ignored on arm64         │
├─────────────────┼────────────────────────────┼────────────────────────┼──────────────────────────┤
│ pto-peaks       │ Poisson -log10 p in        │ AVX-512F+DQ (opt-in),  │ PTO_PEAKS_FORCE_ISA      │
│                 │ batches of 256 bases;      │ AVX2+FMA, NEON, scalar │ ran: scalar output ==    │
│                 │ lanes 8/4/2/1 doubles      │                        │ NEON output, byte for    │
│                 │                            │                        │ byte                     │
├─────────────────┼────────────────────────────┼────────────────────────┼──────────────────────────┤
│ scrna_matrix    │ gather_dot, dense_dot;     │ AVX-512F+BW (opt-in),  │ SCRNA_FORCE_ISA; Python  │
│                 │ HNSW distance through      │ AVX2+FMA, NEON, scalar │ simd_isa() reports it.   │
│                 │ ScrnaCosineSpace; sparse   │                        │ static + CI sweep only   │
│                 │ merge-join is scalar       │                        │ (not built here)         │
├─────────────────┼────────────────────────────┼────────────────────────┼──────────────────────────┤
│ genomic_toolkit │ none: no intrinsics or     │ compiler auto-vec      │ n/a                      │
│ cuttag_profiler │ CPU probes in first-       │ only; -march=native    │                          │
│                 │ party sources (grep)       │ OFF by default         │                          │
│                 │                            │                        │                          │
└─────────────────┴────────────────────────────┴────────────────────────┴──────────────────────────┘

C   MEMORY: BOUNDS, ALIGNMENT, I/O
┌─────────────────┬───────────────────────────────────────────┬────────────────────────────────────┐
│ engine          │ resident bound                            │ layout and I/O primitive           │
├─────────────────┼───────────────────────────────────────────┼────────────────────────────────────┤
│ fastq_stream    │ (lanes x 15 slots + in-flight) x          │ Buffer 64 B aligned; ring cursors  │
│                 │ <= 384 KiB; not a function of             │ 128 B apart; read(2) into a 1 MiB  │
│                 │ input size (47 blocks = 15.8 MB, ran)     │ stage with fadvise(SEQUENTIAL)     │
├─────────────────┼───────────────────────────────────────────┼────────────────────────────────────┤
│ genomic_toolkit │ size histogram 2,001 counters; dup        │ Fragment 16 B; alignas(64) arrays; │
│                 │ marker O(deepest pile-up), or             │ LineReader 1 MiB, grows to a       │
│                 │ O(distinct) with --unsorted; batch 64 KiB │ 64 MiB cap; htslib for BAM/CRAM    │
├─────────────────┼───────────────────────────────────────────┼────────────────────────────────────┤
│ cuttag_profiler │ rows x bins x 8 B allocated up front;     │ per-reducer aligned accumulators;  │
│                 │ one tally per worker; matrix cache        │ indexed queries through htslib;    │
│                 │ 512 MiB and 4 entries (server only)       │ server: regions via the fd         │
├─────────────────┼───────────────────────────────────────────┼────────────────────────────────────┤
│ pto-peaks       │ 2 rings x 128 KiB; candidates 32 B        │ 4 MiB block reads, 1 MiB line cap; │
│                 │ each + ~12 B BH scratch; not a            │ BGZF decoded natively in place     │
│                 │ function of contig length                 │ (no per-record allocation)         │
├─────────────────┼───────────────────────────────────────────┼────────────────────────────────────┤
│ scrna_matrix    │ CSR + KnnGraph n*k*8 B; per-row           │ AlignedAllocator 64 B; NumPy bufs  │
│                 │ heaps O(n*k) extra; HNSW ~ n*dim*4 B      │ adopted in place if aligned and    │
│                 │ + links (hnsw_index_bytes guards it)      │ immutable; no file I/O             │
└─────────────────┴───────────────────────────────────────────┴────────────────────────────────────┘

memory-mapping: NONE. No mmap, madvise or MAP_* in first-party code (grep); the one page-cache
hint is posix_fadvise(SEQUENTIAL) in fastq_stream (reader.cpp:76). Libraries do their own I/O.

CI covers it: *_FORCE_ISA sweeps, OMP_NUM_THREADS 1/8/32, ASan+UBSan, static-link, TSan container
(2/8/32 threads x2 repeats; ci.yml:184-218, 534-546)  (static)
```

**Verified against**

- `fastq_stream/include/fastq_stream/simd.hpp:4-40,514-568`, `peaks/include/peaks/poisson_model.hpp:1715-1889`, `scrna_matrix/include/matrix/simd_math.hpp:704-770`
- `scrna_matrix/include/matrix/knn_graph.hpp:131-427`, `cuttag_profiler/src/signal_calc.cpp:227-451`, `genomic_toolkit/include/toolkit/frip.hpp:100-166`, `cuttag_profiler/CMakeLists.txt:269`
- `.github/workflows/ci.yml:184-218,534-546` (ISA, thread and TSan sweeps)
- grep sweeps over first-party sources: no `mmap`/`madvise`/`MAP_*`; the only `posix_fadvise` is `fastq_stream/src/reader.cpp:76`
- ran: `FQ_FORCE_ISA` and `PTO_PEAKS_FORCE_ISA` change the selected kernel and results stay identical; `OMP_NUM_THREADS` 1 vs 8; profiler 1 vs 8 threads
- static only: `scrna_matrix` (not built here, needs pybind11), AVX2/AVX-512 paths (this host is arm64), TSan container
