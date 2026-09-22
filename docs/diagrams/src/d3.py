from canvas import Canvas

W_, H_ = 100, 92
FW = 99
NOTE_MAX = 95      # inner text must stop before the panel's right border (x=98)


def row(c: Canvas, x0: int, y: int, specs, w: int, h: int, gap: int = 2):
    """Boxes left->right with an arrow between neighbours."""
    out, x = [], x0
    for i, (title, lines) in enumerate(specs):
        assert x + w <= 97, f"row overruns the panel: box {title!r} ends at {x + w - 1}"
        c.box(x, y, w, h, title, lines, pad=1)
        out.append((x, y, w, h))
        if i > 0:
            px, _, pw, _ = out[i - 1]
            ym = y + h // 2
            c.path([(px + pw - 1, ym), (x - 1, ym)])
        x += w + gap
    return out


def note(c: Canvas, x: int, y: int, s: str):
    if len(s) > NOTE_MAX - (x - 2):
        c._fail(f"note ({len(s)}) would cross the panel border: {s[:50]!r}")
    c.text(x, y, s)


def build() -> list[str]:
    c = Canvas(W_, H_)
    c.text(0, 0, "DIAGRAM 3 / 5   the four command-line engines: control path and core loop")
    c.text(0, 1, "each panel: top = main() control flow; bottom = the loop that touches every record or region")

    # ======================= 1. fastq_stream ================================
    y = 3
    c.box(0, y, FW, 15, "1  fastq_stream   qc | trim   src/main.cpp  src/pipeline.cpp", [])
    row(c, 2, y + 1, [
        ("argv", ["-i -o -j -t", "--qc-only ..."]),
        ("validate", ["range checks;", "usage err: rc 2"]),
        ("probe", ["BAM/CRAM by", "content (files)"]),
        ("run_pipeline", ["see Diagram 2:", "3+D+W threads"]),
        ("report + rc", ["text, JSON (-j)", "rc 0 | 1 | 2"]),
    ], w=17, h=5)
    note(c, 2, y + 7, "per record, in each of W workers (RecordProcessor::process, pipeline.cpp:152):")
    row(c, 2, y + 8, [
        ("checks", ["header/SEQ/QUAL", "printable;", "Phred+64 rebase"]),
        ("QC on RAW", ["counts, Q20/Q30", "GC, N, cycles"]),
        ("trim", ["front, tail,", "sliding window"]),
        ("adapters", ["seed + extend", "TruSeq, Nextera"]),
        ("filters", ["min_len, max_n,", "min_mean_q"]),
    ], w=17, h=5)
    note(c, 2, y + 13, "rc 2 also when qualities fit only Phred+64 and no --phred-offset was given (report first)")

    # ======================= 2. genomic_toolkit =============================
    y = 19
    c.box(0, y, FW, 22, "2  genomic_toolkit   sizes | markdup | frip   src/main.cpp  src/fragment_stream.cpp", [])
    row(c, 2, y + 1, [
        ("argv", ["sizes|markdup|frip", "checked int parse"]),
        ("peaks [frip]", ["read_peaks_bed", "malformed = fatal"]),
        ("PeakSet::build", ["sort; merge overlap", "+ book-ended; SoA"]),
        ("bind", ["pre-seed contigs", "peak tid -> stream"]),
    ], w=22, h=5)
    note(c, 2, y + 7, "stream: one pass, no sort, no temp file")
    row(c, 2, y + 8, [
        ("stream_fragments", [".bam .cram .sam:", "htslib; else .bedpe;", "else fragment BED"]),
        ("reader", ["BAM: TLEN>0 = 1 pair", "text: 1 MiB reader", "cap 64 MiB / line"]),
        ("BatchEmitter", ["Fragment = 16 B", "4096/batch = 64 KiB", "size histogram here"]),
        ("consumers", ["sizes: histogram", "markdup: marker", "frip: marker+count"]),
    ], w=22, h=6)
    c.box(2, y + 15, 46, 5, "DuplicateMarker", [
        "sorted: keys (end,strand) of ONE start",
        "hash (--unsorted): 64-bit fingerprints",
        "sets kDuplicate in place on the batch",
    ], pad=1)
    c.box(50, y + 15, 47, 5, "FripCounter.add_batch", [
        "omp for, if batch >= 512; per-thread",
        "OverlapCursor; 5 integer sums; checked_add",
    ], pad=1)
    note(c, 2, y + 20, "exit 2: nothing parsed | out of order (markdup; frip w/o --keep-dups) | no overlap + bad naming")

    # ======================= 3. cuttag_profiler =============================
    y = 42
    c.box(0, y, FW, 20, "3  cuttag_profiler   profile   src/main.cpp  src/signal_calc.cpp  src/bam_reader.cpp", [])
    row(c, 2, y + 1, [
        ("guards", ["output != input", "limits checked"]),
        ("regions", ["BED/GTF (.gz);", "names sanitised"]),
        ("BamReader", ["header, index;", "mapped total:", "index metadata"]),
        ("compute", ["compute_matrix:", "N<=256, <=rows;", "blocks of rows"]),
        ("finish", ["scale, means,", "write TSV;", "all skip: rc 2"]),
    ], w=17, h=5)
    note(c, 2, y + 7, "per region, in each of N workers (signal_calc.cpp:242); one leased htsFile each, no BGZF pool:")
    row(c, 2, y + 8, [
        ("window", ["TSS/center/TES", "anchor; minus", "strand mirrored"]),
        ("skip?", ["no such contig,", "or off-contig:", "row = missing"]),
        ("query_bins", ["widen by max", "frag length;", "sam_itr_queryi"]),
        ("alignment", ["flags + MAPQ;", "leftmost mate", "= TLEN span;", "reads, bases"]),
        ("row", ["fill bins;", "minus reversed", "= 5' to 3'"]),
    ], w=17, h=6)
    note(c, 2, y + 15, "matrix = regions x bins doubles, allocated up front (bins <= 10,000,000, window <= ~922 Gbp)")
    note(c, 2, y + 16, "column means: min(N, rows/64, 64 MiB/(cols*8)) reducers, cache-line-aligned accumulators")
    note(c, 2, y + 17, "skipped rows stay in the matrix, named, but are excluded from means, BPM total and scaling")

    # ======================= 4. pto-peaks ===================================
    y = 63
    c.box(0, y, FW, 21, "4  pto-peaks   main.cpp + sliding_window, caller, poisson_model, bam_streamer (.hpp)", [])
    row(c, 2, y + 1, [
        ("detect", ["peek first byte:", "1f => BAM, else BED", "(stdin ok)"]),
        ("records", ["BAM: BGZF, CRC32,", "leftmost mate TLEN>0", "BED: 4 MiB blocks"]),
        ("Driver.add", ["order, bounds check;", "seek across empty", "coverage deserts"]),
        ("finalize", ["BH over ALL", "candidates, then", "narrowPeak"]),
    ], w=22, h=5)
    note(c, 2, y + 7, "per base, in one thread, memory independent of contig length (sliding_window.hpp, caller.hpp):")
    row(c, 2, y + 8, [
        ("MultiScaleWindow", ["delta+pileup rings,", "2 x 16384 x int64", "= 256 KiB"]),
        ("3 backgrounds", ["centred 1k, 5k, 10k", "running sums; lag =", "5000 bases"]),
        ("PeakCaller", ["prefilter, queue of", "256, SIMD Poisson", "-log10 p"]),
        ("state machine", ["IDLE IN_PEAK", "SUMMIT_SEARCH", "CLOSING (max_gap)"]),
    ], w=22, h=6)
    note(c, 2, y + 15, "lambda = max(genome floor, 3 local rates); a degenerate lambda scores p = 1 (fails closed)")
    note(c, 2, y + 16, "peaks are held (32 B each) until finalize(): BH needs every p-value, so only the INPUT streams")
    note(c, 2, y + 17, "exit 2: unsorted, contig not in sizes, bad field, BED fragment > 4096 bp (BAM: skipped)")
    return c.render()
