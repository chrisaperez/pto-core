from canvas import Canvas
from d4 import table, note

W_, H_ = 100, 120


def build() -> list[str]:
    c = Canvas(W_, H_)
    c.text(0, 0, "DIAGRAM 5 / 5   concurrency, ISA dispatch and memory")
    c.text(0, 1, "evidence tags: (ran) = executed in this review;  (static) = read from source or CI config only")

    # ================= A. concurrency ========================================
    y = 3
    c.text(0, y, "A   CONCURRENCY MODEL PER ENGINE")
    y = table(
        c, 0, y + 1, [0, 18, 46, 74, 99],
        ["engine", "threads and primitive", "shared state and sync", "ordering / determinism"],
        [
            (["fastq_stream", "", "", "", ""],
             ["3+D+W std::thread:", "reader, D inflaters,", "assembler, W workers,", "writer; BAM adds an", "htslib decode pool"],
             ["SPSC rings: acq/rel", "atomics, no CAS; pools:", "std::mutex; ErrorGate:", "mutex + abort atomic", ""],
             ["lane order = stream", "order. raw, gzip, BGZF", "and NEON vs scalar:", "1 identical result", "(ran)"]),
            (["genomic_toolkit", "", "", "", ""],
             ["1 streaming thread;", "OpenMP parallel for in", "FRiP if batch >= 512,", "schedule(static)", ""],
             ["PeakSet immutable;", "per-thread cursor;", "integer reductions;", "checked_add outside", "the parallel region"],
             ["identical JSON at", "OMP_NUM_THREADS=1", "and 8; region tiles", "sum to the whole", "(ran)"]),
            (["cuttag_profiler", "profile", "", "", ""],
             ["N std::thread workers,", "block-partitioned rows;", "then reducers =", "min(N, rows/64, mem)", ""],
             ["one leased htsFile per", "worker (mutex list);", "index + header shared", "read-only; alias lookup", "under a mutex"],
             ["matrix and profile", "byte-identical at 1", "and 8 threads", "(sha256, ran)", ""]),
            (["cuttag_profiler", "serve", "", "", ""],
             ["httplib pool: 16", "threads (CMake sets", "THREAD_POOL_COUNT)", "", ""],
             ["reader map: mutex;", "MatrixCache: mutex,", "FIFO by bytes", "", ""],
             ["not benchmarked;", "gates unit-tested", "and probed live (ran)", "", ""]),
            (["pto-peaks", "", "", "", ""],
             ["1 thread; SIMD lanes", "only inside the", "Poisson evaluation", "", ""],
             ["nothing shared", "", "", "", ""],
             ["scalar and NEON", "outputs byte-", "identical (ran)", "", ""]),
            (["scrna_matrix", "", "", "", ""],
             ["OpenMP: norms static;", "tile pairs in circle-", "method rounds, dynamic;", "HNSW build and query", ""],
             ["one heap per row; a", "round's tile pairs are", "disjoint (no locks);", "write_mutex_ guards", "index mutation"],
             ["brute force: tie-break", "by row, deterministic;", "HNSW approximate, not", "reproducible threaded", "(static only)"]),
        ], row_h=5,
    )
    y += 2

    # ================= B. ISA dispatch =======================================
    c.text(0, y, "B   ISA DISPATCH   one mechanism, three engines")
    note(c, 0, y + 1, "per-function target attributes, never a target-wide -mavx2; level resolved once (function static)")
    note(c, 0, y + 2, "  order: FORCE_ISA env (if compiled in) > __builtin_cpu_supports (checks OS register save too)")
    note(c, 0, y + 3, "         > NEON on AArch64 (no probe) > scalar.  wheels build with every *_NATIVE flag OFF")
    y = table(
        c, 0, y + 5, [0, 18, 47, 72, 99],
        ["engine", "kernels", "ISAs", "override and evidence"],
        [
            (["fastq_stream", "", "", ""],
             ["all_printable, count_ge,", "add_widen_u8_to_u32,", "mismatches, find_seed4", "(sum_u8: tests only)"],
             ["scalar, NEON, AVX2,", "AVX-512BW (opt-in:", "FQ_ENABLE_AVX512)", ""],
             ["FQ_FORCE_ISA. ran:", "NEON default, scalar", "on request; avx2 is", "ignored on arm64"]),
            (["pto-peaks", "", "", ""],
             ["Poisson -log10 p in", "batches of 256 bases;", "lanes 8/4/2/1 doubles", ""],
             ["AVX-512F+DQ (opt-in),", "AVX2+FMA, NEON, scalar", "", ""],
             ["PTO_PEAKS_FORCE_ISA", "ran: scalar output ==", "NEON output, byte for", "byte"]),
            (["scrna_matrix", "", "", ""],
             ["gather_dot, dense_dot;", "HNSW distance through", "ScrnaCosineSpace; sparse", "merge-join is scalar"],
             ["AVX-512F+BW (opt-in),", "AVX2+FMA, NEON, scalar", "", ""],
             ["SCRNA_FORCE_ISA; Python", "simd_isa() reports it.", "static + CI sweep only", "(not built here)"]),
            (["genomic_toolkit", "cuttag_profiler", "", ""],
             ["none: no intrinsics or", "CPU probes in first-", "party sources (grep)", ""],
             ["compiler auto-vec", "only; -march=native", "OFF by default", ""],
             ["n/a", "", "", ""]),
        ], row_h=4,
    )
    y += 2

    # ================= C. memory =============================================
    c.text(0, y, "C   MEMORY: BOUNDS, ALIGNMENT, I/O")
    y = table(
        c, 0, y + 1, [0, 18, 62, 99],
        ["engine", "resident bound", "layout and I/O primitive"],
        [
            (["fastq_stream", "", ""],
             ["(lanes x 15 slots + in-flight) x", "<= 384 KiB; not a function of", "input size (47 blocks = 15.8 MB, ran)"],
             ["Buffer 64 B aligned; ring cursors", "128 B apart; read(2) into a 1 MiB", "stage with fadvise(SEQUENTIAL)"]),
            (["genomic_toolkit", "", ""],
             ["size histogram 2,001 counters; dup", "marker O(deepest pile-up), or", "O(distinct) with --unsorted; batch 64 KiB"],
             ["Fragment 16 B; alignas(64) arrays;", "LineReader 1 MiB, grows to a", "64 MiB cap; htslib for BAM/CRAM"]),
            (["cuttag_profiler", "", ""],
             ["rows x bins x 8 B allocated up front;", "one tally per worker; matrix cache", "512 MiB and 4 entries (server only)"],
             ["per-reducer aligned accumulators;", "indexed queries through htslib;", "server: regions via the fd"]),
            (["pto-peaks", "", ""],
             ["2 rings x 128 KiB; candidates 32 B", "each + ~12 B BH scratch; not a", "function of contig length"],
             ["4 MiB block reads, 1 MiB line cap;", "BGZF decoded natively in place", "(no per-record allocation)"]),
            (["scrna_matrix", "", ""],
             ["CSR + KnnGraph n*k*8 B; per-row", "heaps O(n*k) extra; HNSW ~ n*dim*4 B", "+ links (hnsw_index_bytes guards it)"],
             ["AlignedAllocator 64 B; NumPy bufs", "adopted in place if aligned and", "immutable; no file I/O"]),
        ], row_h=3,
    )
    y += 2
    note(c, 0, y, "memory-mapping: NONE. No mmap, madvise or MAP_* in first-party code (grep); the one page-cache")
    note(c, 0, y + 1, "hint is posix_fadvise(SEQUENTIAL) in fastq_stream (reader.cpp:76). Libraries do their own I/O.")
    note(c, 0, y + 3, "CI covers it: *_FORCE_ISA sweeps, OMP_NUM_THREADS 1/8/32, ASan+UBSan, static-link, TSan container")
    note(c, 0, y + 4, "(2/8/32 threads x2 repeats; ci.yml:184-218, 534-546)  (static)")
    return c.render()
