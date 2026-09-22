from canvas import Canvas

W_, H_ = 100, 64


def build() -> list[str]:
    c = Canvas(W_, H_)
    c.text(0, 0, "DIAGRAM 1 / 5   pto-core: system overview and tool boundaries")
    c.text(0, 1, "solid = component in this repo    dashed = outside it: a caller, or code whose contract is set here")
    c.text(0, 2, "══╦══ the process boundary: argv + files + stdio + exit code, no shared memory, no IPC")

    # ---- consumers -------------------------------------------------------
    xs = [0, 20, 40, 61, 82]
    cons = [
        ("pto-core", ["console cmd", "pto.cli:main"], False),
        ("import pto", ["Python API", "7 functions"], False),
        ("PATH", ["SLURM, Make,", "bare names"], True),
        ("browser", ["dashboard UI", "loopback only"], True),
        ("Python", ["import", "scrna_matrix"], True),
    ]
    for x, (t, ls, dashed) in zip(xs, cons):
        c.box(x, 4, 15, 5, t, ls, dashed=dashed, pad=1)
    c.text(0, 3, "CALLERS")

    # ---- Python orchestration band --------------------------------------
    c.box(0, 11, 32, 10, "packaging/pto (stdlib only)", [
        "cli.py > api.py > _run.py",
        "resolve: $PTO_BIN_DIR >",
        "  pto/_vendor/bin > $PATH",
        "argv list, shell=False",
        "threads = min(affinity,",
        "  cgroup quota, cpu_count)",
        "wait4: wall + peak RSS",
        "rc 2 -> InvalidInput",
    ], pad=2)
    c.path([(7, 8), (7, 10)])       # pto-core CLI -> band
    c.path([(27, 8), (27, 10)])     # import pto  -> band

    # ---- bus (process boundary) -----------------------------------------
    BUS_Y = 24
    c.hline(9, 63, BUS_Y, "d")
    for x in (9, 29, 49, 63):
        c.path([(x, BUS_Y), (x, 26)], style="d")
    c.path([(16, 20), (16, BUS_Y - 1)])                   # Python band -> bus
    c.path([(47, 8), (47, BUS_Y - 1)], dashed=True)       # shell / PATH -> bus
    c.text(39, 15, "bare tool names")
    c.text(0, BUS_Y - 1, "PROCESS")
    c.text(0, BUS_Y, "BOUNDARY")

    # browser and scanpy get their own doors
    c.path([(68, 8), (68, 26)], dashed=True)
    c.text(61, 15, "HTTP 127.0.0.1")
    c.path([(89, 8), (89, 26)], dashed=True)
    c.text(86, 15, "in-proc")

    # ---- engines ----------------------------------------------------------
    EY, EH = 27, 13
    eng = [
        (0, "fastq_stream", [
            "C++20 CLI  0.1.0",
            "in : FASTQ raw/gz",
            "     BGZF, BAM(c)",
            "out: QC JSON+text",
            "     FASTQ: file,",
            "     FIFO, stdout",
            "threads 3+D+W",
            "SIMD runtime ISA",
            "rc 0 | 1 | 2 (a)",
        ]),
        (20, "genomic_toolkit", [
            "C++20 CLI  0.1.0",
            "sizes, markdup,",
            "frip",
            "in : BAM/CRAM,",
            "     BEDPE, BED",
            "out: JSON stdout",
            "     (+ hist TSV)",
            "threads 1+OpenMP",
            "no hand SIMD",
            "rc 0 | 1 | 2",
        ]),
        (40, "pto-peaks", [
            "C++20 CLI  0.1.0",
            "in : frag BED or",
            "     sorted BAM,",
            "     (stdin ok)",
            "out: narrowPeak",
            "BAM: own decoder",
            "threads 1",
            "SIMD runtime ISA",
            "  (Poisson only)",
            "rc 0 | 1 | 2",
        ]),
        (60, "cuttag_profiler", [
            "C++20 CLI  1.0.0",
            "profile | serve",
            "in : indexed BAM+",
            "     BED or GTF",
            "out: matrix.tsv,",
            "     profile.tsv",
            "     or HTTP (d)",
            "threads N workers",
            "no hand SIMD",
            "rc 0 | 1 | 2",
        ]),
        (80, "scrna_matrix", [
            "C++20 header-only",
            "+ pybind11 module",
            "in : NumPy/SciPy",
            "     CSR buffers",
            "out: KnnGraph as",
            "     NumPy views",
            "threads OpenMP",
            "SIMD runtime ISA",
            "no CLI, no rc (b)",
        ]),
    ]
    for x, t, ls in eng:
        c.box(x, EY - 1, 19, EH, "", [t] + ls, pad=1)

    # ---- link dependencies ------------------------------------------------
    LY = EY + EH + 1
    c.text(0, LY - 1, "LINKS")
    libs = [
        (0, ["zlib + libdeflate", "htslib (opt) (c)"]),
        (20, ["htslib (AUTO)", "OpenMP (opt)"]),
        (40, ["libdeflate | zlib", "no htslib"]),
        (60, ["htslib (required)", "httplib+json (d)"]),
        (80, ["OpenMP (required)", "hnswlib (vendor)"]),
    ]
    for x, ls in libs:
        c.box(x, LY, 19, 4, "", ls, dashed=True, pad=1)

    # ---- build + distribution ---------------------------------------------
    BY = LY + 6
    c.box(0, BY, 99, 6, "build and distribution", [
        "cmake -S .  ->  PTO_BUILD_{SCRNA_MATRIX,FASTQ_STREAM,CUTTAG_PROFILER,GENOMIC_TOOLKIT,PEAKS}=ON",
        "wheel pto-core     py3-none-<platform>: 4 CLIs in pto/_vendor/bin (scrna_matrix excluded)",
        "wheel scrna-matrix cp39..cp313:        scrna_matrix_py extension (OpenMP hard-fails configure)",
        "no shared layer: 4 CLIs each carry a private file_identity.hpp; no common headers or runtime",
    ], pad=2)

    c.text(0, BY + 7, "(a) fastq_stream also exits 2 on CLI usage errors; the other CLIs use 1 there and reserve 2")
    c.text(0, BY + 8, "    for 'the input's shape made the answer wrong' (naming mismatch, unsorted, all regions skipped)")
    c.text(0, BY + 9, "(b) no wrapper in pto-core: installed with `pip install scrna-matrix`, imported in-process")
    c.text(0, BY + 10, "(c) FQ_ENABLE_HTS defaults OFF in the module, ON in the wheel   (d) compiled out by PTO_CLOUD_BUILD")
    return c.render()


if __name__ == "__main__":
    from canvas import check_render
    rows = build()
    print("\n".join(rows))
    print(check_render(rows) or "lint ok", "| rows:", len(rows), "| max width:", max(len(r) for r in rows))
