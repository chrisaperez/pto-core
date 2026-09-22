"""Assemble docs/diagrams/ascii-diagrams.md from the five renderers."""
import importlib
import sys
from pathlib import Path

from canvas import Canvas, check_render

OUT = Path(sys.argv[1])

INTRO = """# pto-core: terminal-grade flowcharts

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
"""

SOURCES = {
    "d1": ("Diagram 1: system overview and tool boundaries", [
        "`packaging/pto/_binaries.py:44-84` (TOOLS, search order), `_run.py:119-200`, `_cpu.py:75-83`, `api.py:53-55`, `cli.py:229-385`",
        "`pyproject.toml:87-88,119,131,186-249` (console script, ABI tag, install dir, `PTO_BUILD_*`, `*_NATIVE` OFF)",
        "`CMakeLists.txt:28-58` (module switches, `PTO_CLOUD_BUILD`)",
        "exit codes: `fastq_stream/src/main.cpp:82-85,338-351`, `genomic_toolkit/src/main.cpp:253-423`, `cuttag_profiler/src/main.cpp:322,421-447`, `peaks/src/main.cpp:63-65`",
        "no shared layer: four `modules/*/include/*/file_identity.hpp` copies; SIMD only in `fastq_stream/.../simd.hpp`, `peaks/.../poisson_model.hpp`, `scrna_matrix/.../simd_math.hpp` (grep)",
        "ran: `--version` of the four binaries (0.1.0, 0.1.0, 1.0.0, 0.1.0); layout of the built `dist/pto_core-0.1.6-*.whl` files",
    ]),
    "d2": ("Diagram 2: streaming data flow and buffer management", [
        "`fastq_stream/src/pipeline.cpp:26` (ErrorGate), `:152` (process), `:263-345` (sizing, source choice), `:364-613` (the five stage lambdas and thread start)",
        "`include/fastq_stream/pipeline.hpp:46` (ring slots), `:72-98` (Backoff, Lane), `:129-172` (Config)",
        "`include/fastq_stream/buffer.hpp:28-35,103-153` (chunk, prefix, alignment, pool, mutex), `lockfree_queue.hpp:19-104`",
        "`src/reader.cpp:32-303,350-380` (stage, sniffing, BGZF/gzip/raw), `include/fastq_stream/record.hpp:36-54`, `phred_calculator.hpp:186-246`",
        "`src/hts_input.cpp:14-137`, `include/fastq_stream/file_identity.hpp:63-66` (device+inode identity)",
        "ran: 14 thread configurations (`-t 1..12`, plain and BGZF), 0 mismatches against the formula; framing detection; FIFO blocks until a reader attaches and delivers 239,984 lines; hard link and symlink outputs refused, input untouched",
    ]),
    "d3": ("Diagram 3: the four command-line engines", [
        "`fastq_stream`: `src/main.cpp`, `src/pipeline.cpp:152-253`",
        "`genomic_toolkit`: `src/main.cpp:253-423`, `src/fragment_stream.cpp:62-85,271-501`, `include/toolkit/types.hpp:58-84`, `dup_marker.hpp:156-212`, `frip.hpp:100-166`, `peak_set.hpp:154-235`, `src/line_reader.hpp:153-157`",
        "`cuttag_profiler`: `src/main.cpp:116-325`, `src/signal_calc.cpp:98-451`, `src/bam_reader.cpp:123-372`, `include/profiler/types.hpp:139-262`",
        "`pto-peaks`: `src/main.cpp:331-430,594-764`, `sliding_window.hpp:71-388`, `caller.hpp:57-500`, `poisson_model.hpp:1784-1889`, `bam_streamer.hpp:100-110`",
        "ran: FRiP with `--keep-dups` = 148,680 in-peak fragments (README figure reproduced); region tiling sums to the whole; matrix and profile byte-identical at 1 and 8 threads; `pto-peaks` output identical for BAM, BED, stdin BED and stdin BAM",
    ]),
    "d4": ("Diagram 4: integration and orchestration", [
        "`packaging/pto/_binaries.py:64-84`, `_run.py:119-200`, `_cpu.py:75-83`, `cli.py:83-92,229-385`, `errors.py`",
        "`pyproject.toml`, `.github/workflows/ci.yml:49-546`, `.github/workflows/wheels.yml`, built wheels in `dist/`",
        "`cuttag_profiler/src/http_server.cpp:208-235,277-303,401-803`, `include/profiler/matrix_cache.hpp:52-55`, `CMakeLists.txt:269,383-405`",
        "`deploy/Dockerfile`, `deploy/docker-compose.yml`, `deploy/entrypoint.sh`",
        "ran: every HTTP status on the diagram (401, 421, 403, 415, 413, 400, 404, 200) against a live loopback server; the Python wrapper's argv, exit-code mapping and `exec` passthrough; `strings` on the shipped wheel binary shows the `serve` code is present",
    ]),
    "d5": ("Diagram 5: concurrency, ISA dispatch and memory", [
        "`fastq_stream/include/fastq_stream/simd.hpp:4-40,514-568`, `peaks/include/peaks/poisson_model.hpp:1715-1889`, `scrna_matrix/include/matrix/simd_math.hpp:704-770`",
        "`scrna_matrix/include/matrix/knn_graph.hpp:131-427`, `cuttag_profiler/src/signal_calc.cpp:227-451`, `genomic_toolkit/include/toolkit/frip.hpp:100-166`, `cuttag_profiler/CMakeLists.txt:269`",
        "`.github/workflows/ci.yml:184-218,534-546` (ISA, thread and TSan sweeps)",
        "grep sweeps over first-party sources: no `mmap`/`madvise`/`MAP_*`; the only `posix_fadvise` is `fastq_stream/src/reader.cpp:76`",
        "ran: `FQ_FORCE_ISA` and `PTO_PEAKS_FORCE_ISA` change the selected kernel and results stay identical; `OMP_NUM_THREADS` 1 vs 8; profiler 1 vs 8 threads",
        "static only: `scrna_matrix` (not built here, needs pybind11), AVX2/AVX-512 paths (this host is arm64), TSan container",
    ]),
}


def main() -> int:
    Canvas.problems = []
    parts = [INTRO]
    total_rows = 0
    for key in ("d1", "d2", "d3", "d4", "d5"):
        mod = importlib.import_module(key)
        rows = mod.build()
        lint = check_render(rows)
        if Canvas.problems or lint:
            print(key, "PROBLEMS", Canvas.problems, lint)
            return 1
        title, srcs = SOURCES[key]
        parts.append(f"\n## {title}\n")
        parts.append("```text")
        parts.extend(rows)
        parts.append("```\n")
        parts.append("**Verified against**\n")
        parts.extend(f"- {s}" for s in srcs)
        total_rows += len(rows)
        print(f"{key}: {len(rows)} rows, max width {max(len(r) for r in rows)}")
    OUT.write_text("\n".join(parts) + "\n", encoding="utf-8")
    print("wrote", OUT, "| diagram rows:", total_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
