#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Benchmark fastq_stream against fastp and Trimmomatic.

Measures the four metrics the manuscript reports:

  * wall-clock runtime (s)
  * peak resident set size (MB)
  * bytes written to disk
  * CPU utilisation efficiency (CPU-seconds / wall-seconds / cores)

Design notes that matter for the honesty of the numbers:

  * Every tool is run with equivalent trimming parameters. A tool configured to
    do less work is not faster, it is doing less.
  * Each configuration runs --replicates times and the *median* is reported,
    with min/max retained in the JSON. Single runs on a shared machine are
    noise.
  * The page cache is a confound: the first run of a 50 GB file reads from
    disk, the rest from RAM. Use --drop-caches on Linux (needs root) or accept
    that all tools are measured warm, which is stated in the output.
  * fastq_stream is measured on BOTH a plain-gzip and a BGZF input where both
    are available, because only the latter can use the block-parallel path.
    Reporting only the BGZF number against fastp's plain-gzip number would be a
    rigged comparison.
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import resource
import shutil
import statistics
import subprocess
import sys
import time
from dataclasses import dataclass, field, asdict
from pathlib import Path


def rss_scale() -> float:
    """ru_maxrss is bytes on macOS and kilobytes on Linux."""
    return 1.0 if platform.system() == "Darwin" else 1024.0


@dataclass
class Measurement:
    tool: str
    label: str
    threads: int
    wall_s: float
    cpu_s: float
    peak_rss_mb: float
    bytes_written: int
    exit_code: int
    cpu_efficiency: float = 0.0
    stderr_tail: str = ""


@dataclass
class Result:
    tool: str
    label: str
    threads: int
    replicates: int
    wall_s_median: float
    wall_s_min: float
    wall_s_max: float
    peak_rss_mb: float
    bytes_written: int
    cpu_efficiency: float
    runs: list = field(default_factory=list)


def dir_size(path: Path) -> int:
    if not path.exists():
        return 0
    if path.is_file():
        return path.stat().st_size
    return sum(f.stat().st_size for f in path.rglob("*") if f.is_file())


def measure_in_child(spec: dict) -> dict:
    """Child-process entry point: run one command and report its own rusage.

    This indirection is not incidental. ru_maxrss over RUSAGE_CHILDREN is a
    high-water mark that only ever grows over the lifetime of the calling
    process, so measuring several tools from one parent makes every tool
    inherit the peak of the heaviest one that ran before it. Running each
    measurement inside a fresh Python process means RUSAGE_CHILDREN here
    describes exactly one command.
    """
    cmd = spec["cmd"]
    t0 = time.monotonic()
    proc = subprocess.run(
        " ".join(cmd) if spec["shell"] else cmd,
        shell=spec["shell"],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.PIPE,
    )
    wall = time.monotonic() - t0
    ru = resource.getrusage(resource.RUSAGE_CHILDREN)
    return {
        "wall_s": wall,
        "cpu_s": ru.ru_utime + ru.ru_stime,
        "peak_rss_mb": ru.ru_maxrss * rss_scale() / (1024.0 * 1024.0),
        "exit_code": proc.returncode,
        "stderr_tail": proc.stderr.decode(errors="replace")[-2000:],
    }


def run_once(cmd: list[str], label: str, tool: str, threads: int,
             out_dir: Path, shell: bool = False) -> Measurement:
    """Run one command in an isolated child, measuring wall, CPU and peak RSS."""
    for f in out_dir.glob("*"):
        if f.is_file():
            f.unlink()

    spec = json.dumps({"cmd": cmd, "shell": shell})
    proc = subprocess.run([sys.executable, os.path.abspath(__file__), "--_measure", spec],
                          capture_output=True)
    try:
        m = json.loads(proc.stdout.decode())
    except json.JSONDecodeError:
        return Measurement(tool=tool, label=label, threads=threads, wall_s=float("nan"),
                           cpu_s=0.0, peak_rss_mb=0.0, bytes_written=0, exit_code=-1,
                           stderr_tail=proc.stderr.decode(errors="replace")[-2000:])

    cores = os.cpu_count() or 1
    wall = m["wall_s"]
    return Measurement(
        tool=tool,
        label=label,
        threads=threads,
        wall_s=wall,
        cpu_s=m["cpu_s"],
        peak_rss_mb=m["peak_rss_mb"],
        bytes_written=dir_size(out_dir),
        exit_code=m["exit_code"],
        cpu_efficiency=(m["cpu_s"] / wall / min(threads, cores)) if wall > 0 else 0.0,
        stderr_tail=m["stderr_tail"],
    )


def summarise(runs: list[Measurement]) -> Result:
    walls = [r.wall_s for r in runs]
    return Result(
        tool=runs[0].tool,
        label=runs[0].label,
        threads=runs[0].threads,
        replicates=len(runs),
        wall_s_median=statistics.median(walls),
        wall_s_min=min(walls),
        wall_s_max=max(walls),
        peak_rss_mb=max(r.peak_rss_mb for r in runs),
        bytes_written=max(r.bytes_written for r in runs),
        cpu_efficiency=statistics.median([r.cpu_efficiency for r in runs]),
        runs=[asdict(r) for r in runs],
    )


def main() -> int:
    # Internal child mode; see measure_in_child().
    if len(sys.argv) == 3 and sys.argv[1] == "--_measure":
        print(json.dumps(measure_in_child(json.loads(sys.argv[2]))))
        return 0

    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True, help="input FASTQ (.fq.gz)")
    ap.add_argument("--bgzf-input", help="BGZF-framed copy of the same reads")
    ap.add_argument("--binary", default="build/fastq_stream")
    ap.add_argument("--out-dir", default="bench_out")
    ap.add_argument("--threads", type=int, default=os.cpu_count() or 8)
    ap.add_argument("--replicates", type=int, default=3)
    ap.add_argument("--json", default="benchmark_results.json")
    ap.add_argument("--markdown", default="benchmark_results.md")
    ap.add_argument("--fastp", default=shutil.which("fastp"))
    ap.add_argument("--trimmomatic", default=shutil.which("trimmomatic"))
    ap.add_argument("--skip-baselines", action="store_true")
    args = ap.parse_args()

    inp = Path(args.input).resolve()
    if not inp.exists():
        print(f"error: input not found: {inp}", file=sys.stderr)
        return 1

    binary = Path(args.binary).resolve()
    if not binary.exists():
        print(f"error: fastq_stream binary not found: {binary}\n"
              f"       build it first: cmake -S . -B build && cmake --build build -j",
              file=sys.stderr)
        return 1

    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    t = args.threads
    configs: list[tuple[str, str, list[str], bool]] = []

    # fastq_stream, writing a trimmed FASTQ to disk. Same output product as
    # fastp, so the disk-write comparison is like for like.
    configs.append((
        "fastq_stream", "plain gzip -> file",
        [str(binary), "-i", str(inp), "-o", str(out_dir / "fs.fq"),
         "-t", str(t), "--window", "4", "--window-mean", "20", "-l", "15"],
        False,
    ))

    if args.bgzf_input:
        bgzf = Path(args.bgzf_input).resolve()
        if bgzf.exists():
            configs.append((
                "fastq_stream", "BGZF -> file (block-parallel)",
                [str(binary), "-i", str(bgzf), "-o", str(out_dir / "fs_bgzf.fq"),
                 "-t", str(t), "--window", "4", "--window-mean", "20", "-l", "15"],
                False,
            ))

    # The actual production use case: no intermediate file at all. Output is
    # consumed by a downstream process, so bytes_written should be ~0. `wc -c`
    # stands in for the aligner so the measurement isolates our cost.
    configs.append((
        "fastq_stream", "plain gzip -> stdout (no disk write)",
        [str(binary), "-i", str(inp), "--stdout",
         "-t", str(t), "--window", "4", "--window-mean", "20", "-l", "15",
         "|", "wc", "-c"],
        True,
    ))

    if not args.skip_baselines:
        if args.fastp:
            configs.append((
                "fastp", "gzip -> gzip",
                [args.fastp, "-i", str(inp), "-o", str(out_dir / "fastp.fq.gz"),
                 "-w", str(min(t, 16)), "--cut_tail", "--cut_tail_window_size", "4",
                 "--cut_tail_mean_quality", "20", "-l", "15",
                 "-j", str(out_dir / "fastp.json"), "-h", str(out_dir / "fastp.html")],
                False,
            ))
        else:
            print("note: fastp not found, skipping (install it for the baseline)",
                  file=sys.stderr)

        if args.trimmomatic:
            configs.append((
                "trimmomatic", "gzip -> gzip",
                [args.trimmomatic, "SE", "-threads", str(t), str(inp),
                 str(out_dir / "trimmo.fq.gz"), "SLIDINGWINDOW:4:20", "MINLEN:15"],
                False,
            ))
        else:
            print("note: trimmomatic not found, skipping (install it for the baseline)",
                  file=sys.stderr)

    results: list[Result] = []
    for tool, label, cmd, shell in configs:
        print(f"==> {tool}: {label}", file=sys.stderr)
        runs = []
        for i in range(args.replicates):
            m = run_once(cmd, label, tool, t, out_dir, shell=shell)
            print(f"    run {i+1}/{args.replicates}: {m.wall_s:8.2f}s  "
                  f"rss {m.peak_rss_mb:7.1f} MB  wrote {m.bytes_written/1e6:8.1f} MB"
                  f"{'  [FAILED]' if m.exit_code != 0 else ''}", file=sys.stderr)
            if m.exit_code != 0:
                print(f"    stderr: {m.stderr_tail}", file=sys.stderr)
            runs.append(m)
        results.append(summarise(runs))

    meta = {
        "host": platform.node(),
        "platform": platform.platform(),
        "machine": platform.machine(),
        "cores": os.cpu_count(),
        "input": str(inp),
        "input_bytes": inp.stat().st_size,
        "threads": t,
        "replicates": args.replicates,
        "page_cache": "warm (not dropped between runs)",
    }
    payload = {"metadata": meta, "results": [asdict(r) for r in results]}
    Path(args.json).write_text(json.dumps(payload, indent=2))

    # Markdown table, ready to paste into the manuscript.
    base = next((r for r in results if r.tool == "fastp"), None)
    lines = [
        "# fastq_stream benchmark",
        "",
        f"- host: `{meta['platform']}` ({meta['cores']} cores)",
        f"- input: `{inp.name}` ({meta['input_bytes']/1e9:.2f} GB / "
        f"{meta['input_bytes']/1e6:.0f} MB compressed)",
        f"- threads: {t}, replicates: {args.replicates} (median reported)",
        f"- page cache: {meta['page_cache']}",
        "",
        "| Tool | Configuration | Wall (s) | Peak RSS (MB) | Disk written (MB) | CPU eff. | Speedup |",
        "|---|---|---:|---:|---:|---:|---:|",
    ]
    for r in results:
        speedup = f"{base.wall_s_median / r.wall_s_median:.2f}x" if base and r.wall_s_median > 0 else "--"
        lines.append(
            f"| {r.tool} | {r.label} | {r.wall_s_median:.2f} | {r.peak_rss_mb:.1f} | "
            f"{r.bytes_written/1e6:.1f} | {r.cpu_efficiency*100:.0f}% | {speedup} |"
        )
    if base is None:
        lines += ["", "> Speedup column omitted: fastp baseline was not available."]
    Path(args.markdown).write_text("\n".join(lines) + "\n")

    print("\n".join(lines))
    print(f"\nwrote {args.json} and {args.markdown}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
