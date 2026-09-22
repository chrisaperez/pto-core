from canvas import Canvas

W_, H_ = 100, 110
FW = 99            # full-width box (inner width 95)
NW = 63            # narrow box (leaves room for the bypass lane on the right)
BYP = 73           # x of the direct-lane bypass
MID = 30           # x of the main vertical spine


def tag(s: str, t: str, w: int = 90) -> str:
    return s.ljust(w) + t


def build() -> list[str]:
    c = Canvas(W_, H_)
    c.text(0, 0, "DIAGRAM 2 / 5   fastq_stream: streaming data flow and buffer management")
    c.text(0, 1, "══► a Buffer HANDLE moves by pointer through a lock-free SPSC ring; payload bytes are not copied")
    c.text(0, 2, "(c1)-(c4): the only four places payload bytes ARE copied.   T=--threads  D=inflaters  W=workers")

    y = 4
    # ---- SOURCE -------------------------------------------------------------
    c.box(0, y, FW, 4, "SOURCE", [
        'path | "-" (stdin) | FIFO       BAM/CRAM/SAM recognised by content (htslib), regular files only',
        "framing, first 18 bytes: 1f 8b + BGZF 'BC' subfield = BGZF; 1f 8b alone = gzip; else raw",
    ])
    c.path([(MID, y + 3), (MID, y + 4)])
    y += 5

    # ---- READER -------------------------------------------------------------
    c.box(0, y, FW, 8, "READER  1 thread  pipeline.cpp:364  reader.cpp:169-303", [
        "FdSource: 1 MiB stage (doubles on demand), posix_fadvise(SEQUENTIAL), read(2), EINTR retried",
        tag("BGZF: memcpy up to 4 whole compressed members into one pooled Buffer; lane = k mod D", "(c1)"),
        "gzip: zlib inflate (15+16 bits, inflateReset per member) into a pooled Buffer; serial",
        tag("raw : memcpy stage -> Buffer", "(c1)"),
        "BAM : htslib sam_read1, decode pool of max(1, spare/3) threads, records rewritten as FASTQ",
        "EOF: BGZF must end with its 28-byte empty block; a missing one, or trailing bytes, is an error",
    ])
    ry = y + 7
    c.path([(MID, ry), (MID, ry + 2)])
    c.path([(BYP, ry), (BYP, ry + 1)], style="d", head=None)
    y += 9

    # ---- ring A -------------------------------------------------------------
    ringA = y + 1
    c.box(0, ringA, NW, 3, "lane x D = SpscRing<Buffer,16>", ["15 usable slots + closed flag; BGZF only"])
    c.path([(MID, ringA + 2), (MID, ringA + 3)], style="d")
    y = ringA + 4

    # ---- INFLATERS ----------------------------------------------------------
    c.box(0, y, NW, 6, "INFLATERS  D threads  pipeline.cpp:385", [
        "libdeflate_gzip_decompress, one decompressor per thread",
        "4 members -> one 256 KiB Buffer from in_pool",
        "the compressed Buffer returns to in_pool",
        "D = max(1, spare/3), spare = max(1, T-3); BGZF only",
    ])
    c.path([(MID, y + 5), (MID, y + 6)], style="d")
    y += 7

    # ---- ring B -------------------------------------------------------------
    ringB = y
    c.box(0, ringB, NW, 3, "lane x D = SpscRing<Buffer,16>", ["pop order = push order: chunk k on lane k mod D"])
    c.path([(MID, ringB + 2), (MID, ringB + 3)], style="d")

    # bypass line + label
    for i, ln in enumerate(["direct lane (1 ring)", "used by gzip, raw, BAM:", "no inflate rank, so the",
                            "reader feeds assembler"]):
        c.text(BYP + 2, ry + 3 + i, ln)
    asm_top = ringB + 4
    c.path([(BYP, ry + 1), (BYP, asm_top - 1)], style="d")
    y = asm_top

    # ---- ASSEMBLER ----------------------------------------------------------
    c.box(0, y, FW, 8, "ASSEMBLER  1 thread  pipeline.cpp:409  record.hpp:36", [
        "pop lane[k mod D] (or the direct lane) in the reader's order: stream order, no reorder buffer",
        tag("carry: extend_front(carry) into the 64 KiB prefix ahead of data(), an in-place splice", "(c2)"),
        "cut = last_record_boundary(): memchr '\\n', every 4th newline ends a record ('@' can be QUAL)",
        "whole records -> work_in[j mod W]; trailing partial record -> carry (over 64 KiB = error)",
        "EOF: a last record without its final '\\n' is accepted; any other short tail is an error",
        "one thread on purpose: memchr runs at memory bandwidth, well above the inflate rate",
    ])
    c.path([(MID, y + 7), (MID, y + 8)], style="d")
    y += 9

    # ---- ring C -------------------------------------------------------------
    c.box(0, y, NW, 3, "lane x W = SpscRing<Buffer,16>", ["chunk j on lane j mod W; input of worker i"])
    c.path([(MID, y + 2), (MID, y + 3)], style="d")
    y += 4

    # ---- WORKERS ------------------------------------------------------------
    c.box(0, y, FW, 11, "WORKERS  W threads  pipeline.cpp:488  process() :152", [
        "W = max(1, spare-D). out = out_pool.acquire(): cap 320 KiB = 256 KiB chunk + 64 KiB carry",
        "per record, in order: header/'+' printable -> Phred+64 rebase (only --phred-offset 64)",
        "  -> accumulate_read on the RAW read -> trim_front -> trim_back -> sliding window",
        "  -> adapter seed-and-extend (find_seed4 + mismatches) -> min_len, max_n_rate, min_mean_q",
        tag("  -> surviving records memcpy'd into out", "(c3)"),
        "SIMD (runtime ISA): all_printable, count_ge x2, add_widen_u8_to_u32, mismatches, find_seed4",
        "per-cycle sums: 32-bit CycleScratch folded into 64-bit QcStats once per chunk; merged at join",
        "ALWAYS pushes exactly ONE out Buffer per input chunk, even if empty: writer order needs it",
        "--qc-only: identical work, serialisation skipped, an empty out Buffer is still pushed",
    ])
    c.path([(MID, y + 10), (MID, y + 11)], style="d")
    y += 12

    # ---- ring D -------------------------------------------------------------
    c.box(0, y, NW, 3, "lane x W = SpscRing<Buffer,16>", ["worker i's output: one Buffer per input chunk"])
    c.path([(MID, y + 2), (MID, y + 3)], style="d")
    y += 4

    # ---- WRITER -------------------------------------------------------------
    c.box(0, y, FW, 4, "WRITER  1 thread  pipeline.cpp:566", [
        tag("pop work_out[k mod W] round-robin = original order; write(2) loop, EINTR retried", "(c4)"),
        "SIGPIPE ignored: a closed FIFO becomes EPIPE -> error gate; close() checked (NFS/FUSE ENOSPC)",
    ])
    c.path([(MID, y + 3), (MID, y + 4)])
    y += 5

    # ---- SINK ---------------------------------------------------------------
    c.box(0, y, FW, 4, "SINK", [
        "stdout | file (O_TRUNC) | FIFO: --mkfifo, open(2) BLOCKS until a reader attaches (verified)",
        "refused before anything is opened: -i X -o X, also through hard links and symlinks",
    ])
    y += 5

    # ==================== panel B: buffer anatomy ============================
    c.text(0, y, "BUFFER ANATOMY   buffer.hpp:28-35   BufferPool uses aligned_alloc(64 B)")
    y += 1
    c.box(0, y, 40, 4, "", ["kPrefix 64 KiB: headroom", "extend_front() splices carry here"], pad=2)
    c.box(39, y, 60, 4, "", ["payload 256 KiB (in_pool) / 320 KiB (out_pool)", "data() .... size() .............. capacity()"], pad=2)
    y += 4
    c.text(0, y, "handle = {pool*, base*, cap, size, front_} + seq; seq is written by reader/inflater, never read")
    c.text(0, y + 1, "pool free list: std::mutex (buffer.hpp:118), one lock pair per 256 KiB of stream, off the data path")
    y += 3

    # ==================== panel C: ring + backoff ============================
    c.text(0, y, "SPSC RING  lockfree_queue.hpp:47  kCacheLine = 128 B on arm64 and x86-64, 64 B elsewhere")
    y += 1
    c.box(0, y, FW, 6, "one producer, one consumer, no CAS", [
        "head_, cached_tail_, tail_, cached_head_, slots_[16]: each alignas(128); no line is shared",
        "try_push: relaxed load of tail; reload head (acquire) only when full; release-store tail",
        "try_pop : relaxed load of head; reload tail (acquire) only when empty; release-store head",
        "Backoff pipeline.hpp:72: 128 spins, 16 yields, then 100 us sleeps (yield-only burned ~11 cores)",
    ])
    y += 7

    # ==================== panel D: control ====================================
    c.text(0, y, "CONTROL AND BOUNDS")
    y += 1
    c.box(0, y, FW, 7, "error gate, termination, memory", [
        "ErrorGate (pipeline.cpp:26): first error wins; abort flag polled in every ring wait; then join",
        "close: the producer release-stores closed; the consumer exits when the ring is empty AND closed",
        "a thread that cannot start trips the gate and joins the rest: no terminate (pipeline.cpp:599)",
        "resident bound: (lanes x 15 + in-flight) x <= 384 KiB per block; not a function of input size",
        "measured: 47 blocks = 15.8 MB pool for a 12.8 MB FASTQ at -t 6 (fastq_stream 0.1.0, this build)",
    ])
    y += 8

    # ---- copy legend --------------------------------------------------------
    c.text(0, y, "(c1) stage -> Buffer: raw, BGZF compressed batch; plain gzip inflates into the Buffer directly")
    c.text(0, y + 1, "(c2) carry splice: only the straddling partial record, at most 64 KiB")
    c.text(0, y + 2, "(c3) worker serialises survivors into the out Buffer      (c4) writer's write(2): the kernel's copy")
    return c.render()
