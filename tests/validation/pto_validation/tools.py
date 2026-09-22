# SPDX-License-Identifier: MIT
"""Subprocess execution and the file-format readers the cases share.

Deliberately self-contained: it does not import from scripts/benchmark_audit.py
and that script does not import from here, even though both time subprocesses.
The repo has no cross-module layer by design (see the README's "Repository
layout"), and the two have genuinely different jobs -- this one cares whether a
command produced the right answer, that one cares what it cost.
"""

from __future__ import annotations

import json
import os
import shlex
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, Sequence


class ToolError(RuntimeError):
    """A reference or pto command exited non-zero.

    Carries the full command line and the tail of stderr, because the useful
    half of a bedtools or deeptools failure is always in the last few lines and
    a bare CalledProcessError throws it away.
    """

    def __init__(self, command: Sequence[str], returncode: int, stderr: str):
        self.command = list(command)
        self.returncode = returncode
        self.stderr = stderr
        tail = "\n".join(stderr.strip().splitlines()[-12:])
        super().__init__(
            f"command failed with exit {returncode}:\n"
            f"  {shlex.join(self.command)}\n"
            f"{tail}"
        )


@dataclass
class CommandResult:
    command: list[str]
    returncode: int
    stdout: str
    stderr: str
    wall_seconds: float

    def json(self) -> dict:
        return json.loads(self.stdout)


def run(
    command: Sequence[str],
    *,
    cwd: Path | None = None,
    stdout_path: Path | None = None,
    check: bool = True,
    timeout: float | None = None,
    env: dict[str, str] | None = None,
) -> CommandResult:
    """Run a command, capturing output and wall time.

    `stdout_path` streams stdout straight to a file instead of through a pipe.
    That matters for the tools whose output *is* the artefact -- a genome-wide
    bedgraph is hundreds of megabytes, and buffering it into a Python string to
    immediately write it back out is the one way this harness could plausibly
    run the machine out of memory.
    """
    command = [str(part) for part in command]
    merged_env = {**os.environ, **(env or {})}

    started = time.perf_counter()
    if stdout_path is not None:
        stdout_path.parent.mkdir(parents=True, exist_ok=True)
        with stdout_path.open("wb") as sink:
            process = subprocess.run(
                command,
                cwd=cwd,
                stdout=sink,
                stderr=subprocess.PIPE,
                timeout=timeout,
                env=merged_env,
            )
        stdout_text = ""
    else:
        process = subprocess.run(
            command,
            cwd=cwd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=timeout,
            env=merged_env,
        )
        stdout_text = process.stdout.decode("utf-8", "replace")

    elapsed = time.perf_counter() - started
    stderr_text = (process.stderr or b"").decode("utf-8", "replace")

    if check and process.returncode != 0:
        raise ToolError(command, process.returncode, stderr_text)

    return CommandResult(
        command=command,
        returncode=process.returncode,
        stdout=stdout_text,
        stderr=stderr_text,
        wall_seconds=elapsed,
    )


# ---------------------------------------------------------------------------
# Format readers
# ---------------------------------------------------------------------------


def read_bed(path: Path, *, max_records: int | None = None) -> list[tuple[str, int, int]]:
    """Read (contig, start, end) from a BED, skipping track/browser/comment lines.

    Intentionally tolerant of the header junk that real UCSC and MACS2 output
    carries, and intentionally strict about the coordinate columns -- a
    non-integer start is a malformed file, not a line to skip, because silently
    dropping it would understate every downstream overlap metric.
    """
    out: list[tuple[str, int, int]] = []
    with path.open() as handle:
        for lineno, line in enumerate(handle, 1):
            if not line.strip() or line[0] == "#":
                continue
            if line.startswith(("track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            try:
                out.append((fields[0], int(fields[1]), int(fields[2])))
            except ValueError as exc:
                raise ValueError(f"{path}:{lineno}: bad BED coordinates: {line!r}") from exc
            if max_records is not None and len(out) >= max_records:
                break
    return out


def read_tsv_matrix(path: Path, *, skip_columns: int = 0) -> Iterator[list[float]]:
    """Stream a numeric TSV row by row.

    A generator rather than a list: cuttag_profiler's region x bin matrix for a
    genome-wide promoter set is comfortably larger than the machine's RAM at
    Python float sizes, and every consumer here reduces it streaming anyway.
    """
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line[0] == "#":
                continue
            fields = line.rstrip("\n").split("\t")[skip_columns:]
            row: list[float] = []
            for field in fields:
                try:
                    row.append(float(field))
                except ValueError:
                    row.append(float("nan"))
            yield row


def read_profile_tsv(path: Path) -> tuple[list[float], list[float]]:
    """Read cuttag_profiler's `--out-profile` output: offset, mean_signal.

    Returns (offsets, values). The header line is detected by trying to parse
    it rather than by assuming its presence, because the same reader is pointed
    at deeptools output further down and the two disagree about headers.
    """
    offsets: list[float] = []
    values: list[float] = []
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line[0] == "#":
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue
            try:
                offset, value = float(fields[0]), float(fields[1])
            except ValueError:
                continue  # header
            offsets.append(offset)
            values.append(value)
    return offsets, values


def read_histogram_tsv(path: Path) -> dict[int, float]:
    """Read a two-column `size<TAB>count` histogram into {size: count}."""
    hist: dict[int, float] = {}
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line[0] == "#":
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue
            try:
                hist[int(fields[0])] = hist.get(int(fields[0]), 0.0) + float(fields[1])
            except ValueError:
                continue  # header
    return hist


def deeptools_matrix_profile(path: Path) -> list[float]:
    """Collapse a deeptools computeMatrix gzip output into a mean meta-profile.

    computeMatrix writes a JSON header line prefixed with '@', then one row per
    region: six BED-ish annotation columns followed by the per-bin values, with
    'nan' for bins falling off the end of a contig. Averaging down the columns
    while skipping those nans is exactly what plotProfile shows, and it is the
    only representation the two tools can be compared in without reimplementing
    deeptools' region sorting.
    """
    import gzip

    sums: list[float] = []
    counts: list[int] = []
    opener = gzip.open if path.suffix == ".gz" else open

    with opener(path, "rt") as handle:  # type: ignore[operator]
        for line in handle:
            if not line or line[0] == "@":
                continue
            fields = line.rstrip("\n").split("\t")[6:]
            if not fields:
                continue
            if not sums:
                sums = [0.0] * len(fields)
                counts = [0] * len(fields)
            for i, field in enumerate(fields):
                if i >= len(sums):
                    break
                if field in ("nan", "NA", ""):
                    continue
                try:
                    sums[i] += float(field)
                except ValueError:
                    continue
                counts[i] += 1

    return [s / c if c else float("nan") for s, c in zip(sums, counts)]
