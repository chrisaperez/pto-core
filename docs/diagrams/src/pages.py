"""The five slides. Coordinates are absolute pixels on a 1920x1080 canvas."""
from graphics import PageBuilder

REV = "commit 0ba20d1"


def page1():
    b = PageBuilder("p1", "1 Overview", "Five independent C++20 engines behind one process boundary",
                    "No shared runtime and no shared code: callers reach an engine only through argv, files, stdio and an exit code.",
                    f"Sources: packaging/pto/*, pyproject.toml, CMakeLists.txt, modules/*/src/main.cpp ({REV})")
    b.box("z_call", "panel", 64, 124, 1792, 140, "CALLERS")
    b.box("c_cli", "py", 88, 164, 272, 88, "pto-core CLI", ["console script pto.cli:main"], "z_call")
    b.box("c_api", "py", 384, 164, 272, 88, "import pto", ["Python API, 7 functions"], "z_call")
    b.box("c_path", "ext", 1000, 164, 196, 88, "shell / SLURM", ["bare tool names on PATH"], "z_call")
    b.box("c_brw", "ext", 1318, 164, 196, 88, "browser", ["dashboard, loopback only"], "z_call")
    b.box("c_scp", "ext", 1636, 164, 196, 88, "scanpy / AnnData", ["import scrna_matrix"], "z_call")

    b.box("z_py", "zone_py", 64, 292, 912, 270, "PYTHON ORCHESTRATION   packaging/pto   stdlib only, zero runtime dependencies")
    cols = [88, 384, 680]
    b.box("py_cli", "py", cols[0], 344, 272, 96, "cli.py", ["argparse, 12 subcommands", "exec uses os.execv"], "z_py")
    b.box("py_api", "py", cols[1], 344, 272, 96, "api.py", ["7 functions", "None kwargs are omitted"], "z_py")
    b.box("py_cpu", "py", cols[2], 344, 272, 96, "_cpu.py", ["min(affinity, cgroup quota,", "cpu_count), never below 1"], "z_py")
    b.box("py_bin", "py", cols[0], 452, 272, 96, "_binaries.py", ["$PTO_BIN_DIR > pto/_vendor/bin", "> $PATH"], "z_py")
    b.box("py_run", "py", cols[1], 452, 272, 96, "_run.py", ["Popen, shell=False", "os.wait4: wall time + peak RSS"], "z_py")
    b.box("py_err", "py", cols[2], 452, 272, 96, "errors.py", ["rc 2 = InvalidInput", "other rc = ToolError"], "z_py")

    b.box("band", "boundary", 64, 590, 1792, 48,
          "PROCESS BOUNDARY   argv + files + stdio + exit code   no shared memory   no IPC")

    b.box("z_eng", "zone_cpp", 64, 666, 1792, 250, "C++20 ENGINES")
    ex = [88, 440, 792, 1144, 1496]
    b.box("e_fq", "cpp", ex[0], 712, 328, 190, "fastq_stream 0.1.0", [
        "CLI: streaming QC and trimming", "in: FASTQ raw/gz/BGZF; BAM/CRAM/SAM", "out: QC JSON + text; FASTQ to file/FIFO/stdout",
        "threads 3+D+W; SIMD runtime", "exit 0 | 1 | 2 (usage errors too)"], "z_eng")
    b.box("e_gtk", "cpp", ex[1], 712, 328, 190, "genomic_toolkit 0.1.0", [
        "CLI: sizes | markdup | frip", "in: BAM/CRAM, BEDPE, fragment BED", "out: JSON on stdout (+ histogram TSV)",
        "threads 1 + OpenMP (FRiP); no hand SIMD", "exit 0 | 1 | 2"], "z_eng")
    b.box("e_pk", "cpp", ex[2], 712, 328, 190, "pto-peaks 0.1.0", [
        "CLI: streaming peak caller", "in: fragment BED or sorted BAM (stdin ok)", "out: narrowPeak; BAM decoded natively",
        "threads 1; SIMD runtime (Poisson)", "exit 0 | 1 | 2"], "z_eng")
    b.box("e_cp", "cpp", ex[3], 712, 328, 190, "cuttag_profiler 1.0.0", [
        "CLI: profile | serve", "in: indexed BAM + BED/GTF", "out: matrix.tsv + profile.tsv, or HTTP",
        "threads N workers; no hand SIMD", "exit 0 | 1 | 2"], "z_eng")
    b.box("e_sc", "cpp", ex[4], 712, 328, 190, "scrna_matrix", [
        "header-only C++20 + pybind11 module", "in: NumPy/SciPy CSR buffers", "out: KnnGraph as NumPy views",
        "threads OpenMP; SIMD runtime", "no CLI, no exit code; own wheel"], "z_eng")

    libs = [("l_fq", ["zlib + libdeflate", "htslib optional (ON in the wheel)"]), ("l_gtk", ["htslib (AUTO-detected)", "OpenMP optional"]),
            ("l_pk", ["libdeflate or zlib (BAM/BGZF)", "no htslib"]), ("l_cp", ["htslib (required)", "httplib + json (server only)"]),
            ("l_sc", ["OpenMP (required, hard fail)", "hnswlib (vendored, patched)"])]
    for i, (lid, body) in enumerate(libs):
        b.box(lid, "ext", ex[i], 928, 328, 56, "", body)
    b.box("build", "build", 64, 996, 1792, 44, "", [
        "cmake -S . : PTO_BUILD_* x5 (PTO_CLOUD_BUILD off)  ->  wheel pto-core py3-none-<platform>, 4 CLIs in pto/_vendor/bin  |  "
        "wheel scrna-matrix cp39-cp313  |  four private copies of file_identity.hpp"])

    # edges
    b.edge("c_cli", "py_cli", "control", s="b", d="t")
    b.edge("c_api", "py_api", "control", s="b", d="t")
    b.edge("py_cli", "py_api", "control")
    b.edge("py_cpu", "py_api", "control", s="l", d="r")
    b.edge("py_api", "py_run", "control", s="b", d="t")
    b.edge("py_run", "py_bin", "control", s="l", d="r")
    b.edge("py_run", "py_err", "error", s="r", d="l")
    b.edge("py_run", "band", "contract", s="b", d="t", d_at=(520 - 64) / 1792)
    b.edge("c_path", "band", "contract", "bare tool names", s="b", d="t", d_at=(1098 - 64) / 1792)
    for eid, cx in (("e_fq", 252), ("e_gtk", 604), ("e_pk", 956), ("e_cp", 1308)):
        b.edge("band", eid, "contract", s="b", d="t", s_at=(cx - 64) / 1792)
    b.edge("c_brw", "e_cp", "payload", "HTTP 127.0.0.1", s="b", d="t", via=[(1416, 690), (1450, 690)], d_at=(1450 - 1144) / 328)
    b.edge("c_scp", "e_sc", "payload", "in-process import", s="b", d="t", d_at=(1734 - 1496) / 328)
    b.page.talk = [
        "Five engines: four command-line tools and one Python extension. Each is its own CMake project and builds and tests with the other four deleted.",
        "The only coupling is the process boundary: argv, files, stdio and an exit code. There is no shared runtime, and the output-aliasing guard is deliberately copied four times.",
        "What's left dashed here is generic -- a shell and a browser -- not any named external product; this slide stays inside what pto-core itself ships.",
    ]
    return b.page


def page2():
    b = PageBuilder("p2", "2 Streaming", "fastq_stream moves buffer handles, not bytes",
                    "3 + D + W threads joined by lock-free SPSC rings; a 256 KiB payload is copied at four points only.",
                    f"Sources: fastq_stream/src/pipeline.cpp:263-613, reader.cpp:32-380, buffer.hpp, lockfree_queue.hpp ({REV}); thread math checked at 14 configurations")
    y1, h = 136, 200
    b.box("src", "data", 64, y1, 232, h, "SOURCE", ["path | - | FIFO", "BAM/CRAM/SAM by content", "(regular files only)", "framing: 1f 8b + BC = BGZF", "1f 8b alone = gzip; else raw"])
    b.box("read", "cpp", 336, y1, 360, h, "READER  1 thread", [
        "FdSource 1 MiB stage + fadvise", "BGZF: <= 4 whole members per Buffer (c1)", "gzip: zlib inflate into the Buffer",
        "raw: memcpy stage -> Buffer (c1)", "BAM: htslib decode pool, FASTQ text", "EOF: BGZF must end in its 28 B block"])
    b.box("ringA", "sync", 736, y1, 200, h, "lane x D", ["SpscRing<Buffer,16>", "15 usable slots", "closed flag", "BGZF only"])
    b.box("infl", "cpp", 976, y1, 320, h, "INFLATERS  D threads", ["libdeflate, 1 decompressor each", "4 members -> one 256 KiB Buffer",
                                                                    "compressed Buffer -> in_pool", "D = max(1, spare/3)"])
    b.box("ringB", "sync", 1336, y1, 200, h, "lane x D", ["pop order = push order", "chunk k on lane k mod D"])
    y2 = 400
    b.box("asm", "cpp", 1526, y2, 330, h, "ASSEMBLER  1 thread", ["pop lane[k mod D] or the direct lane", "carry spliced into 64 KiB prefix (c2)",
                                                                    "cut at the 4th newline (memchr)", "whole records -> lane j mod W", "partial tail -> carry (max 64 KiB)"])
    b.box("ringC", "sync", 1290, y2, 200, h, "lane x W", ["input of worker i", "chunk j on lane j mod W"])
    b.box("work", "cpp", 890, y2, 360, h, "WORKERS  W threads", [
        "out Buffer from out_pool (320 KiB)", "checks -> Phred+64 rebase -> QC on RAW", "trim -> adapters -> filters",
        "survivors memcpy'd into out (c3)", "SIMD: runtime ISA kernels", "always ONE out Buffer per chunk"])
    b.box("ringD", "sync", 650, y2, 200, h, "lane x W", ["output of worker i", "one Buffer per input chunk"])
    b.box("writ", "cpp", 364, y2, 250, h, "WRITER  1 thread", ["pop lane k mod W in order", "write(2), EINTR retried (c4)", "SIGPIPE ignored; close checked"])
    b.box("sink", "data", 64, y2, 260, h, "SINK", ["stdout | file | FIFO", "--mkfifo: open blocks", "until a reader attaches", "refused: -i X -o X"])
    b.chain(["src", "read", "ringA", "infl", "ringB"], "payload")
    b.edge("ringB", "asm", "payload", s="r", d="t")
    b.edge("read", "asm", "payload", "direct lane (gzip, raw, BAM)", s="b", d="t", via=[(516, 376), (1608, 376)], d_at=0.25)
    b.edge("asm", "ringC", "payload", s="l", d="r")
    b.edge("ringC", "work", "payload", s="l", d="r")
    b.edge("work", "ringD", "payload", s="l", d="r")
    b.edge("ringD", "writ", "payload", s="l", d="r")
    b.edge("writ", "sink", "payload", s="l", d="r")

    b.box("p_buf", "panel", 64, 632, 880, 184, "BUFFER ANATOMY   buffer.hpp:28-35")
    b.box("b_pre", "data", 88, 676, 280, 96, "kPrefix 64 KiB", ["headroom in front of data()", "extend_front() splices the carry"], "p_buf")
    b.box("b_pay", "data", 368, 676, 552, 96, "payload 256 KiB (in_pool)   320 KiB (out_pool)", ["data() ... size() ... capacity()", "64 B aligned by aligned_alloc"], "p_buf")
    b.box("b_note", "note", 88, 780, 832, 28, "", ["handle = {pool*, base*, cap, size, front_} + seq (written, never read); free list under std::mutex"], "p_buf")
    b.box("p_ring", "panel", 976, 632, 880, 184, "SPSC RING   lockfree_queue.hpp:47")
    for i, (lab, body) in enumerate([("head_", ["consumer"]), ("cached_tail_", ["consumer-private"]), ("tail_", ["producer"]),
                                     ("cached_head_", ["producer-private"]), ("slots_[16]", ["Buffer handles", "15 usable"])]):
        b.box(f"r{i}", "sync", 1000 + i * 172, 676, 160, 96, lab, body, "p_ring")
    b.box("r_note", "note", 1000, 780, 832, 28, "", ["each member alignas(128), so producer and consumer share no line; kCacheLine = 128 B on arm64 and x86-64"], "p_ring")
    b.box("p_copy", "panel", 64, 832, 880, 184, "THE FOUR COPY POINTS")
    for i, (lab, body) in enumerate([("(c1) stage -> Buffer", ["raw and BGZF compressed batch; gzip inflates in place"]),
                                     ("(c2) carry splice", ["only the straddling partial record, at most 64 KiB"]),
                                     ("(c3) serialise", ["worker writes surviving records into the out Buffer"]),
                                     ("(c4) write(2)", ["the kernel's copy to file, FIFO or stdout"])]):
        b.box(f"c{i + 1}", "callout", 88 + i * 212, 876, 196, 124, lab, body, "p_copy")
    b.box("p_ctl", "panel", 976, 832, 880, 184, "CONTROL AND BOUNDS")
    b.box("k_gate", "sync", 1000, 876, 272, 124, "ErrorGate", ["first error wins", "abort flag polled in every ring wait", "then join, then throw"], "p_ctl")
    b.box("k_term", "sync", 1284, 876, 272, 124, "Termination", ["producer release-stores closed", "consumer exits on empty AND closed"], "p_ctl")
    b.box("k_mem", "data", 1568, 876, 272, 124, "Memory bound", ["(lanes x 15 + in-flight) x <= 384 KiB", "measured 47 blocks = 15.8 MB, 12.8 MB FASTQ, -t 6"], "p_ctl")
    b.page.talk = [
        "Buffers, not bytes, move between threads: 16-slot single-producer single-consumer rings hand Buffer handles across, and the payload is copied at four points only.",
        "Order comes from topology, not sequence numbers: chunk k rides lane k mod D, and the assembler and writer read lanes in the same round-robin order.",
        "Memory is bounded by the number of lanes, not the input: 15.8 MB of pool for a 12.8 MB file at -t 6, measured.",
    ]
    return b.page


def page3():
    b = PageBuilder("p3", "3 Engines", "Each engine is one pass: a control path and a core loop",
                    "Top row of each panel is main(); the row below is the loop that touches every record or region.",
                    f"Sources: modules/*/src/main.cpp and the files named in each panel ({REV})")

    def panel(pid, x, y, title):
        b.box(pid, "panel", x, y, 880, 436, title)
        return x + 24, y + 64

    def row(prefix, x, y, w, h, gap, specs, kind, parent, first=None):
        ids = []
        for i, (lab, body) in enumerate(specs):
            nid = f"{prefix}{i}"
            b.box(nid, kind, x + i * (w + gap), y, w, h, lab, body, parent)
            ids.append(nid)
        b.chain(ids, "control" if kind == "py" else "payload")
        return ids

    # 1 fastq_stream
    x, y = panel("pa", 64, 124, "1  fastq_stream   qc | trim   (main.cpp, pipeline.cpp)")
    a = row("a", x, y, 150, 100, 20, [("argv", ["-i -o -j -t", "--qc-only ..."]), ("validate", ["ranges, adapters;", "usage error = rc 2"]),
                                      ("probe", ["BAM/CRAM by content", "(regular files)"]), ("run_pipeline", ["see page 2:", "3 + D + W threads"]),
                                      ("report + rc", ["text, JSON (-j)", "rc 0 | 1 | 2"])], "cpp", "pa")
    bb = row("b", x, y + 156, 150, 110, 20, [("checks", ["header, SEQ, QUAL", "printable;", "Phred+64 rebase"]), ("QC on RAW", ["counts, Q20/Q30", "GC, N, per cycle"]),
                                             ("trim", ["front, tail,", "sliding window"]), ("adapters", ["seed + extend", "TruSeq, Nextera"]),
                                             ("filters", ["min_len, max_n,", "min_mean_q"])], "cpp", "pa")
    b.edge("a3", "b0", "control", "per record, in each of W workers", s="b", d="t")
    b.box("pa_n", "note", x, y + 290, 832, 40, "", ["rc 2 also when the qualities fit only Phred+64 and no --phred-offset was given (the report is printed first)"], "pa")

    # 2 genomic_toolkit
    x, y = panel("pb", 976, 124, "2  genomic_toolkit   sizes | markdup | frip")
    row("c", x, y, 190, 100, 20, [("argv", ["sizes|markdup|frip", "checked int parse"]), ("peaks [frip]", ["read_peaks_bed", "malformed = fatal"]),
                                  ("PeakSet::build", ["sort, merge overlap", "and book-ended; SoA"]), ("bind", ["pre-seed contigs", "peak tid -> stream tid"])], "cpp", "pb")
    row("d", x, y + 116, 190, 110, 20, [("stream_fragments", [".bam .cram .sam: htslib;", "else .bedpe; else BED"]), ("reader", ["BAM: TLEN>0 = 1 pair", "text: 1 MiB, cap 64 MiB/line"]),
                                        ("BatchEmitter", ["Fragment = 16 B", "4096 per batch = 64 KiB", "size histogram here"]),
                                        ("consumers", ["sizes: histogram", "markdup: marker", "frip: marker + count"])], "cpp", "pb")
    b.edge("c1", "d0", "control", s="b", d="t")
    b.box("dm", "cpp", x, y + 246, 408, 96, "DuplicateMarker", ["sorted: keys (end, strand) of ONE start", "hash (--unsorted): 64-bit fingerprints", "sets kDuplicate in place"], "pb")
    b.box("fc", "cpp", x + 424, y + 246, 408, 96, "FripCounter.add_batch", ["omp for if batch >= 512; per-thread cursor", "5 integer sums; checked_add outside"], "pb")
    b.edge("d3", "dm", "payload", s="b", d="t", d_at=0.5)
    b.edge("d3", "fc", "payload", s="b", d="t", d_at=0.5)
    b.box("pb_n", "note", x, y + 350, 832, 20, "", ["exit 2: nothing parsed | out of order (markdup; frip without --keep-dups) | no overlap and bad naming"], "pb")

    # 3 cuttag_profiler
    x, y = panel("pc", 64, 580, "3  cuttag_profiler   profile   (main.cpp, signal_calc.cpp, bam_reader.cpp)")
    row("f", x, y, 150, 100, 20, [("guards", ["output != input", "limits checked"]), ("regions", ["BED/GTF (.gz)", "names sanitised"]),
                                  ("BamReader", ["header, index;", "mapped total from", "index metadata"]), ("compute", ["compute_matrix:", "N<=256, N<=rows;", "row blocks"]),
                                  ("finish", ["scale, means,", "write TSV;", "all skipped: rc 2"])], "cpp", "pc")
    row("g", x, y + 156, 150, 110, 20, [("window", ["TSS/center/TES", "anchor; minus", "strand mirrored"]), ("skip?", ["no such contig,", "or off the contig:", "row = missing"]),
                                        ("query_bins", ["widen by max fragment;", "sam_itr_queryi"]), ("alignment", ["flags + MAPQ;", "leftmost mate:", "TLEN span;", "reads + bases"]),
                                        ("row", ["fill bins;", "minus strand reversed:", "5' to 3'"])], "cpp", "pc")
    b.edge("f3", "g0", "control", "per region, in each of N workers", s="b", d="t")
    b.box("pc_n", "note", x, y + 290, 832, 60, "", ["matrix = regions x bins doubles, allocated up front (bins <= 10,000,000)",
                                                     "means: min(N, rows/64, 64 MiB/(cols*8)) reducers, cache-line aligned",
                                                     "skipped rows stay in the matrix but leave means, BPM total and scaling"], "pc")

    # 4 pto-peaks
    x, y = panel("pd", 976, 580, "4  pto-peaks   (main.cpp, sliding_window / caller / poisson_model / bam_streamer .hpp)")
    row("h", x, y, 190, 100, 20, [("detect", ["peek first byte:", "1f = BAM, else BED", "(stdin ok)"]), ("records", ["BAM: BGZF + CRC32,", "leftmost mate TLEN>0", "BED: 4 MiB blocks"]),
                                  ("Driver.add", ["order and bounds", "checked; seek across", "coverage deserts"]), ("finalize", ["BH over ALL", "candidates, then", "narrowPeak"])], "cpp", "pd")
    row("i", x, y + 156, 190, 110, 20, [("MultiScaleWindow", ["delta + pileup rings,", "2 x 16384 x int64", "= 256 KiB"]), ("3 backgrounds", ["centred 1k, 5k, 10k", "running sums;", "lag = 5000 bases"]),
                                        ("PeakCaller", ["prefilter, queue of 256,", "SIMD Poisson -log10 p"]), ("state machine", ["IDLE, IN_PEAK,", "SUMMIT_SEARCH,", "CLOSING (max_gap)"])], "cpp", "pd")
    b.edge("h2", "i0", "control", "per base, one thread", s="b", d="t")
    b.box("pd_n", "note", x, y + 290, 832, 60, "", ["lambda = max(genome floor, 3 local rates); a degenerate lambda scores p = 1 (fails closed)",
                                                     "peaks are held (32 B each) until finalize(): BH needs every p-value, so only the INPUT streams",
                                                     "exit 2: unsorted, contig not in sizes, bad field, BED fragment > 4096 bp (BAM: counted skipped)"], "pd")
    b.page.talk = [
        "Every engine is one streaming pass with no sort and no temporary file. The one place results are held is pto-peaks' candidate list, because Benjamini-Hochberg needs every p-value.",
        "Inputs that would otherwise print a plausible wrong number exit 2 instead: naming mismatch, unsorted input, every region skipped, Phred+64.",
    ]
    return b.page


def page4():
    b = PageBuilder("p4", "4 Integration", "Everything outside the engines meets them at one contract",
                    "argv in; JSON, files and an exit code out. The wrapper adds resolution, thread budgeting and error typing.",
                    f"Sources: packaging/pto/*, pyproject.toml, ci.yml, wheels.yml, http_server.cpp:401-803, deploy/* ({REV}); dashboard probed live")
    b.box("pA", "panel", 64, 124, 1100, 480, "A   FOUR WAYS IN")
    lanes = [
        ("PYTHON API", "py", [("pto.frip(...)", ["7 functions", "None kwargs omitted"]), ("run_tool", ["resolve binary", "argv, shell=False"]),
                              ("Popen + wait4", ["os.wait4: rc + max RSS", "stdio to temp files"]), ("result", ["rc 0: Result", "rc 2: InvalidInput", "else ToolError"])]),
        ("pto-core", "py", [("pto-core frip", ["argparse subcommand", "--json optional"]), ("cli._cmd_frip", ["calls pto.api.frip"]),
                            ("output", ["metrics: stdout", "cost: stderr"]), ("exit code", ["0 ok | 1 error", "2 InvalidInput", "141 broken pipe"])]),
        ("exec / PATH", "py", [("pto-core exec", ["os.execv: no Python", "frame remains"]), ("bin-dir", ["prints the engine", "directory"]),
                               ("PATH", ["export PATH=", "$(pto-core bin-dir):$PATH"]), ("shell, SLURM, Nextflow", ["bare tool names"])]),
    ]
    for li, (lab, kind, specs) in enumerate(lanes):
        y = 168 + li * 160
        b.box(f"lab{li}", "note", 80, y + 34, 112, 24, "", [lab], "pA", size=12)
        ids = []
        for j, (t, body) in enumerate(specs):
            k = "ext" if (kind == "ext" or t.startswith("shell")) else "py"
            nid = f"l{li}{j}"
            b.box(nid, k, 200 + j * 230, y, 210, 96, t, body, "pA")
            ids.append(nid)
        b.chain(ids, "control" if kind == "py" else "control")

    b.box("pB", "panel", 1196, 124, 660, 480, "B   RESOLUTION AND THREADS")
    b.box("r1", "py", 1220, 176, 190, 88, "1  $PTO_BIN_DIR", ["override; no reinstall"], "pB")
    b.box("r2", "py", 1430, 176, 190, 88, "2  pto/_vendor/bin", ["wheel-vendored"], "pB")
    b.box("r3", "py", 1640, 176, 190, 88, "3  $PATH", ["shutil.which"], "pB")
    b.chain(["r1", "r2", "r3"], "control")
    b.box("r_err", "gate", 1220, 296, 610, 84, "none found: BinaryNotFound", ["lists the searched directories and prints the cmake build recipe"], "pB")
    b.edge("r3", "r_err", "error", s="b", d="t", d_at=0.85)
    b.box("r_tools", "callout", 1220, 400, 610, 76, "TOOLS = the four engines", ["a test fails if the wheel vendors a fifth binary that TOOLS lacks"], "pB")
    b.box("r_thr", "callout", 1220, 492, 610, 96, "thread budget", ["min(sched_getaffinity, cgroup quota, cpu_count), at least 1", "passed as -t (fastq_stream) and --threads (cuttag_profiler) only",
                                                                  "genomic_toolkit and pto-peaks get none; OMP_NUM_THREADS is not set"], "pB")

    b.table("tbl", 64, 620, 1100, 274, [190, 220, 200, 290, 200], ["engine", "machine-readable", "artifacts", "exit codes", "thread control"], [
        (["fastq_stream", "qc | trim"], ["JSON file: -j PATH", "text: --report PATH", "wrapper reads the JSON"], ["trimmed FASTQ:", "-o FILE | FIFO | -"],
         ["0 ok, 1 runtime", "2 usage error, or", "Phred+64 suspected"], ["-t N: wrapper passes it"]),
        (["genomic_toolkit", "sizes markdup frip"], ["JSON on stdout (--json)"], ["histogram TSV:", "--histogram FILE"],
         ["0 ok, 1 usage or system", "2 input shape: naming,", "order, nothing parsed"], ["no flag; OpenMP reads", "OMP_NUM_THREADS", "(wrapper leaves it unset)"]),
        (["cuttag_profiler", "profile | serve"], ["none: summary on stderr"], ["matrix.tsv, profile.tsv", "(or HTTP: serve)"],
         ["0 ok, 1 error", "2 every region skipped"], ["--threads N: wrapper", "passes it"]),
        (["pto-peaks"], ["none: counts on stderr"], ["narrowPeak: --out FILE", "(default stdout)"],
         ["0 ok, 1 usage or I/O", "2 unsorted, contig missing,", "bad field"], ["single-threaded;", "no thread option"]),
    ], size=12)

    b.box("pD", "panel", 1196, 620, 660, 274, "D   BUILD, PACKAGE, PUBLISH")
    b.box("d1", "build", 1220, 668, 196, 96, "cmake -S .", ["PTO_BUILD_* x5", "PTO_CLOUD_BUILD off"], "pD")
    b.box("d2", "build", 1428, 668, 196, 96, "scikit-build-core", ["install-dir pto/_vendor", "headers dropped"], "pD")
    b.box("d3", "build", 1636, 668, 196, 96, "cibuildwheel", ["cp312 pin, py3-none-*", "doctor + pytest"], "pD")
    b.box("d4", "build", 1220, 784, 196, 96, "repair", ["delocate (macOS)", "auditwheel (Linux)"], "pD")
    b.box("d5", "build", 1428, 784, 196, 96, "publish", ["v* tag -> PyPI", "OIDC, no token"], "pD")
    b.box("d6", "callout", 1636, 784, 196, 96, "vendored libs", ["macOS: libomp, libhts", "Linux: libgomp, static htslib"], "pD")
    b.chain(["d1", "d2", "d3"], "control")
    b.edge("d3", "d4", "control", s="b", d="t")
    b.edge("d4", "d5", "control")

    b.box("pE", "zone_gate", 64, 908, 1792, 140, "TRUST BOUNDARY   127.0.0.1 only   Host, then Origin, then token, in one pre-routing handler for /api/*")
    specs = [("client", "ext", ["browser or curl", "URL carries ?t=TOKEN"]), ("Host  421", "gate", ["loopback name + port", "defeats DNS rebinding"]),
             ("Origin  403", "gate", ["absent ok; present must", "be this server"]), ("token  401", "gate", ["?t= or Bearer", "constant-time compare"]),
             ("body gates", "gate", ["<= 1 MiB: 413", "JSON type: 415"]), ("resolve_input", "gate", ["under data root: 400", "open once, on the fd"]),
             ("compute", "cpp", ["compute_matrix", "HTTP pool = 16 threads"]), ("cache + fetch", "cpp", ["MatrixCache 512 MiB, 4", "GET /api/matrix?token=M"])]
    ids = []
    for i, (t, k, body) in enumerate(specs):
        b.box(f"w{i}", k, 88 + i * 220, 952, 200, 88, t, body, "pE")
        ids.append(f"w{i}")
    b.chain(ids, "control")
    b.page.talk = [
        "The Python wrapper is thin: it resolves a binary, budgets threads, runs argv without a shell, and types the exit code.",
        "The dashboard is the only network path: Host, then Origin, then a 256-bit token, all in one pre-routing handler, then a body cap and path confinement.",
        "A finding: the wheel's cuttag_profiler still contains serve, reachable through pto-core exec, although the Python API never builds it.",
    ]
    return b.page


def page5():
    b = PageBuilder("p5", "5 Concurrency and ISA", "Determinism is tested, not assumed: threads, ISA and memory",
                    "(ran) = executed in this review on an arm64 host; (static) = read from source or CI configuration only.",
                    f"Sources: simd.hpp:4-568, poisson_model.hpp:1715-1889, simd_math.hpp:704-770, knn_graph.hpp:131-427, signal_calc.cpp, frip.hpp, ci.yml ({REV})")
    b.table("tA", 64, 124, 1792, 0, [170, 520, 560, 542], ["engine", "threads and primitive", "shared state and synchronisation", "ordering and determinism evidence"], [
        (["fastq_stream"], ["3+D+W std::thread: reader, D inflaters, assembler, W workers, writer; BAM adds an htslib decode pool"],
         ["SPSC rings with acquire/release atomics, no CAS; buffer pools under std::mutex; ErrorGate = mutex + abort atomic"],
         ["lane order = stream order; raw, gzip, BGZF and NEON vs scalar give 1 identical result (ran)"]),
        (["genomic_toolkit"], ["1 streaming thread; OpenMP parallel for in FRiP when the batch has >= 512 fragments, schedule(static)"],
         ["PeakSet is immutable; per-thread OverlapCursor; integer reductions; checked_add outside the parallel region"],
         ["identical JSON at OMP_NUM_THREADS=1 and 8; region tiles sum to the whole file (ran)"]),
        (["cuttag_profiler"], ["profile: N std::thread workers, block-partitioned rows, then min(N, rows/64, memory) reducers. serve: httplib pool of 16"],
         ["one leased htsFile per worker (mutex-guarded list); index and header shared read-only; alias lookup under a mutex; reader map and MatrixCache under mutexes"],
         ["profile: matrix and profile byte-identical at 1 and 8 threads (sha256, ran); serve gates probed live (ran)"]),
        (["pto-peaks"], ["1 thread; SIMD lanes only inside the Poisson evaluation"], ["nothing shared"],
         ["scalar and NEON outputs byte-identical (ran)"]),
        (["scrna_matrix"], ["OpenMP: norms static; tile pairs in circle-method rounds (dynamic); HNSW build and query in parallel"],
         ["one heap per row; a round's tile pairs are disjoint, so no locks; write_mutex_ guards index mutation"],
         ["brute force deterministic (tie-break by row); HNSW approximate and not reproducible when threaded (static, not built here)"]),
    ], size=12)
    yB = 124 + b.page.nodes[-1].h + 24
    b.table("tB", 64, yB, 1240, 0, [170, 400, 300, 370], ["engine", "kernels", "ISAs", "override and evidence"], [
        (["fastq_stream"], ["all_printable, count_ge, add_widen_u8_to_u32, mismatches, find_seed4 (sum_u8: tests only)"],
         ["scalar, NEON, AVX2, AVX-512BW (opt-in FQ_ENABLE_AVX512)"], ["FQ_FORCE_ISA. ran: NEON by default, scalar on request, avx2 ignored on arm64"]),
        (["pto-peaks"], ["Poisson -log10 p in batches of 256 bases; lanes 8/4/2/1 doubles"], ["AVX-512F+DQ (opt-in), AVX2+FMA, NEON, scalar"],
         ["PTO_PEAKS_FORCE_ISA. ran: scalar output == NEON output, byte for byte"]),
        (["scrna_matrix"], ["gather_dot, dense_dot; HNSW distance via ScrnaCosineSpace; sparse merge-join is scalar"], ["AVX-512F+BW (opt-in), AVX2+FMA, NEON, scalar"],
         ["SCRNA_FORCE_ISA; Python simd_isa(). static + CI sweep only"]),
        (["genomic_toolkit", "cuttag_profiler"], ["none: no intrinsics or CPU probes in first-party sources (grep)"], ["compiler auto-vectorisation only"],
         ["-march=native is OFF by default in both"]),
    ], size=12)
    hB = b.page.nodes[-1].h
    b.box("mech", "panel", 1336, yB, 520, max(hB, 268), "ONE DISPATCH MECHANISM")
    for i, t in enumerate(["1  FORCE_ISA env, if that ISA is compiled in", "2  __builtin_cpu_supports (also OS save)", "3  NEON on AArch64: no probe", "4  scalar reference"]):
        b.box(f"m{i}", "sync", 1360, yB + 48 + i * 56, 472, 44, t, [], "mech")
    b.chain(["m0", "m1", "m2", "m3"], "control", s="b", d="t")
    yC = yB + max(hB, 268) + 24
    b.table("tC", 64, yC, 1792, 0, [170, 800, 822], ["engine", "resident bound", "layout and I/O primitive"], [
        (["fastq_stream"], ["(lanes x 15 slots + in-flight) x <= 384 KiB; not a function of input size (47 blocks = 15.8 MB, ran)"],
         ["Buffer 64 B aligned; ring cursors 128 B apart; read(2) into a 1 MiB stage with posix_fadvise(SEQUENTIAL)"]),
        (["genomic_toolkit"], ["histogram 2,001 counters; dup marker O(deepest pile-up) or O(distinct) with --unsorted; batch 64 KiB"],
         ["Fragment 16 B; alignas(64) arrays; LineReader 1 MiB growing to a 64 MiB cap; htslib for BAM/CRAM"]),
        (["cuttag_profiler"], ["rows x bins x 8 B allocated up front; one tally per worker; matrix cache 512 MiB / 4 entries (server)"],
         ["per-reducer aligned accumulators; indexed queries through htslib; server parses regions from the fd"]),
        (["pto-peaks"], ["2 rings x 128 KiB; candidates 32 B each + ~12 B BH scratch; not a function of contig length"],
         ["4 MiB block reads, 1 MiB line cap; BGZF decoded natively in place, no per-record allocation"]),
        (["scrna_matrix"], ["CSR + KnnGraph n*k*8 B; per-row heaps O(n*k); HNSW ~ n*dim*4 B + links (hnsw_index_bytes guards it)"],
         ["AlignedAllocator 64 B; NumPy buffers adopted in place if aligned and immutable; no file I/O"]),
        (["all engines"], ["memory-mapping: NONE. No mmap, madvise or MAP_* in first-party code (grep)"],
         ["the only page-cache hint is one posix_fadvise(SEQUENTIAL), fastq_stream reader.cpp:76"]),
    ], size=12)
    b.page.talk = [
        "One dispatch pattern, three engines: fastq_stream, pto-peaks and scrna_matrix. genomic_toolkit and cuttag_profiler contain no hand-written SIMD.",
        "Determinism was checked by running: identical results across framings, ISAs and thread counts.",
        "There is no mmap anywhere in first-party code; memory is bounded by structure, not by a mapping.",
    ]
    return b.page


PAGES = [page1, page2, page3, page4, page5]
