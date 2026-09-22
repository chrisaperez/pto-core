from canvas import Canvas

W_, H_ = 100, 112
FW = 99


def row(c, x0, y, specs, w, h, gap=2, dashed=False, limit=98):
    out, x = [], x0
    for i, (title, lines) in enumerate(specs):
        assert x + w <= limit, f"row overruns: {title!r} ends at {x + w - 1}"
        c.box(x, y, w, h, title, lines, pad=1, dashed=dashed)
        out.append((x, y, w, h))
        if i > 0:
            px, _, pw, _ = out[i - 1]
            c.path([(px + pw - 1, y + h // 2), (x - 1, y + h // 2)], dashed=dashed)
        x += w + gap
    return out


def note(c, x, y, s, limit=99):
    if x + len(s) > limit:
        c._fail(f"note ({len(s)}) at x={x} exceeds {limit}: {s[:50]!r}")
    c.text(x, y, s)


def table(c, x, y, xs, header, rows, row_h):
    """xs = column border x positions (first is left edge, last is right edge)."""
    total_h = 3 + len(rows) * row_h + (len(rows) - 1) + 1     # top, header, sep, rows + seps, bottom
    x2, y2 = xs[-1], y + total_h - 1
    c.hline(x, x2, y, "l"); c.hline(x, x2, y2, "l"); c.vline(x, y, y2, "l"); c.vline(x2, y, y2, "l")
    for xv in xs[1:-1]:
        c.vline(xv, y, y2, "l")
    c.hline(x, x2, y + 2, "l")
    for i, h in enumerate(header):
        c.text(xs[i] + 2, y + 1, h) if len(h) <= xs[i + 1] - xs[i] - 3 else c._fail(f"header {h!r} too wide")
    ry = y + 3
    for r, cells in enumerate(rows):
        for ci, cell_lines in enumerate(cells):
            wmax = xs[ci + 1] - xs[ci] - 3
            for li, ln in enumerate(cell_lines):
                if len(ln) > wmax:
                    c._fail(f"cell {ln!r} ({len(ln)}>{wmax}) col {ci}")
                c.text(xs[ci] + 2, ry + li, ln)
        ry += row_h
        if r < len(rows) - 1:
            c.hline(x, x2, ry, "l")
            ry += 1
    return y2


def build() -> list[str]:
    c = Canvas(W_, H_)
    c.text(0, 0, "DIAGRAM 4 / 5   integration and orchestration: how callers reach the binaries")
    c.text(0, 1, "everything below ships in this repository: no external product or proprietary service is named")

    # ================= A. invocation lanes ===================================
    y = 3
    c.text(0, y, "A   FOUR WAYS IN")
    lanes = [
        ("PYTHON API", False, [
            ("pto.frip(...)", ["7 functions;", "None = omitted", "existence checked"]),
            ("run_tool", ["resolve binary;", "argv, shell=False", "stdio -> tmpfiles"]),
            ("Popen+wait4", ["spawn child;", "os.wait4 reaps:", "rc + ru_maxrss"]),
            ("result", ["rc 0: Result", "rc 2: InvalidInput", "rc n: ToolError"]),
        ]),
        ("pto-core", False, [
            ("pto-core frip", ["argparse command", "--json optional", "also: man, doctor"]),
            ("cli._cmd_frip", ["calls pto.api.frip", "same path as above", ""]),
            ("output", ["metrics to stdout;", "cost to stderr", ""]),
            ("exit code", ["0 ok | 1 error", "2 InvalidInput", "141 broken pipe"]),
        ]),
        ("exec / PATH", False, [
            ("exec", ["pto-core exec T", "os.execv replaces", "the process"]),
            ("bin-dir", ["prints dir of the", "first resolvable", "engine"]),
            ("PATH", ["PATH=$(pto-core", " bin-dir):$PATH", "bare tool names"]),
            ("consumer", ["shell, SLURM,", "Nextflow: no", "Python in the path"]),
        ]),
    ]
    ly = y + 1
    for label, dashed, specs in lanes:
        c.text(0, ly + 2, label)
        row(c, 11, ly, specs, w=20, h=5, gap=2, dashed=dashed)
        ly += 6
    y = ly

    c.box(0, y, FW, 7, "resolution and threads   packaging/pto/_binaries.py:64-84  _cpu.py:75", [
        "$PTO_BIN_DIR  >  pto/_vendor/bin  >  $PATH (shutil.which);  none = BinaryNotFound + build hint",
        "TOOLS = the 4 engines (a test fails if the wheel vendors a fifth binary that TOOLS lacks)",
        "threads = min(sched_getaffinity, cgroup v2/v1 quota, cpu_count), never below 1",
        "passed as -t to fastq_stream and --threads to cuttag_profiler ONLY; the other two receive none",
    ])
    y += 8

    # ================= B. result contract =====================================
    c.text(0, y, "B   THE RESULT CONTRACT   (verified by running each binary)")
    y += 1
    yb = table(
        c, 0, y, [0, 18, 39, 58, 81, 99],
        ["engine", "machine-readable", "artifacts", "exit codes", "thread control"],
        [
            (["fastq_stream", "qc | trim", ""],
             ["JSON file: -j PATH", "text: --report", "wrapper reads JSON"],
             ["trimmed FASTQ:", "-o FILE|FIFO|-", ""],
             ["0 ok  1 runtime", "2 usage error, or", "Phred+64 suspected"],
             ["-t N: wrapper", "passes it", ""]),
            (["genomic_toolkit", "sizes markdup", "frip"],
             ["JSON on stdout", "(--json flag)", ""],
             ["histogram TSV:", "--histogram FILE", ""],
             ["0 ok  1 usage/system", "2 input shape:", "naming, order, empty"],
             ["no flag: OpenMP", "OMP_NUM_THREADS", "wrapper: unset"]),
            (["cuttag_profiler", "profile | serve", ""],
             ["none;", "summary on stderr", ""],
             ["matrix.tsv,", "profile.tsv", "(or HTTP: serve)"],
             ["0 ok  1 error", "2 every region", "skipped"],
             ["--threads N", "wrapper passes", ""]),
            (["pto-peaks", "", ""],
             ["none;", "counts on stderr", ""],
             ["narrowPeak", "--out FILE (or", "stdout)"],
             ["0 ok  1 usage/IO", "2 unsorted, contig", "missing, bad field"],
             ["one thread;", "no thread flag", ""]),
        ], row_h=3,
    )
    y = yb + 1
    note(c, 0, y, "fastq_stream exits 2 for usage errors too: pto.InvalidInput fires on a bad parameter there")
    y += 2

    # ================= C. build -> package -> publish =========================
    c.text(0, y, "C   BUILD, PACKAGE, PUBLISH")
    y += 1
    row(c, 0, y, [
        ("cmake -S .", ["PTO_BUILD_* x5", "PTO_CLOUD_BUILD", "default OFF"]),
        ("scikit-build", ["install-dir =", "pto/_vendor;", "no headers"]),
        ("cibuildwheel", ["cp312 pin;", "py3-none-<plat>", "doctor + pytest"]),
        ("repair", ["delocate (mac)", "auditwheel(lnx)", "vendors libs"]),
        ("publish", ["v* tag -> PyPI", "OIDC, no token", "setuptools_scm"]),
    ], w=19, h=5, gap=1, limit=99)
    y += 5
    note(c, 0, y, "vendored libs: macOS pto/.dylibs (omp hts deflate z lzma); Linux pto_core.libs (gomp bz2 lzma)")
    note(c, 0, y + 1, "ci.yml jobs: build-and-test (ISA + OMP sweeps), sanitizers, config-matrix, static-link, TSan")
    y += 3

    # ================= D. dashboard trust boundary ============================
    c.text(0, y, "D   THE ONE NETWORK PATH: cuttag_profiler serve (http_server.cpp; absent in cloud builds)")
    y += 1
    c.box(0, y, FW, 21, "TRUST BOUNDARY   127.0.0.1 only; --allow-remote must be explicit; reach it by SSH tunnel", [], style="d")
    note(c, 2, y + 1, "one handler, /api/* only, in this order (http_server.cpp:430-494): Host before token", limit=98)
    r1 = row(c, 2, y + 2, [
        ("client", ["browser or curl;", "URL has ?t=TOKEN", "(printed once)"]),
        ("Host  421", ["loopback name and", "bound port only:", "defeats rebinding"]),
        ("Origin  403", ["absent is fine;", "present must be", "this server"]),
        ("token  401", ["?t= or Bearer;", "constant_time_", "equals, 256 bit"]),
    ], w=21, h=6, gap=3)
    #, limit=98)
    # second row
    y2 = y + 10
    r2 = row(c, 2, y2, [
        ("body gates", ["<= 1 MiB: 413", "JSON type: 415", "30 s timeouts"]),
        ("resolve_input", ["under data root;", "open once, on the", "fd: no re-lookup"]),
        ("compute", ["compute_matrix as", "in Diagram 3;", "HTTP pool = 16"]),
        ("cache + fetch", ["MatrixCache 512MiB,", "4 entries, FIFO;", "GET /api/matrix"]),
    ], w=21, h=6, gap=3)
    # elbow from the token box down/left into the body-gates box
    xb4 = r1[3][0] + r1[3][2] // 2
    xb1 = r2[0][0] + r2[0][2] // 2
    c.path([(xb4, y + 7), (xb4, y + 8)], head=None)
    c.path([(xb4, y + 8), (xb1, y + 8), (xb1, y2 - 1)])
    note(c, 2, y + 16, "GET /api/matrix?token=M needs BOTH the session token and the matrix token M (256 bit)", limit=98)
    note(c, 2, y + 17, "confinement: canonical path under --data-root (else 400); FIFO, symlink, device refused", limit=98)
    note(c, 2, y + 18, "assets from rodata need no token (CSP default-src 'self'); API replies: no-store, frame DENY", limit=98)
    note(c, 2, y + 19, "container: 127.0.0.1:8080 only; read_only, cap_drop ALL, no-new-privileges, 4 GiB, 2 CPUs", limit=98)
    y += 22

    # ================= E. the serve subcommand: reachability and its gate =====
    c.text(0, y, "E   THE serve SUBCOMMAND: REACHABILITY AND THE CLOUD-BUILD GATE")
    note(c, 0, y + 1, "pto.profile and `pto-core profile` never build it, yet the wheel's binary still contains it,")
    note(c, 0, y + 2, "so `pto-core exec cuttag_profiler serve --data-root D` reaches it: no guard in _cmd_exec (ran it)")
    note(c, 0, y + 3, "PTO_CLOUD_BUILD=ON drops the server entirely; `serve` is refused by name (exit 1); nm-based test")
    return c.render()
