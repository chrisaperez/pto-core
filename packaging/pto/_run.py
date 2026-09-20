"""Run one pto-core binary and describe what it did.

TWO RULES, BOTH INHERITED RATHER THAN INVENTED
----------------------------------------------
1. ARGV, NEVER A SHELL. Every invocation is a list handed to
   `subprocess.Popen(..., shell=False)`. No caller value is ever interpolated
   into a command line. This is `api/app/job_specs.py` rule 1, and it is not
   weaker on a workstation: `cuttag_profiler`'s SECURITY_HTTP L8 finding was a
   `--host` value reaching `std::system()`, on a local tool.

2. EXIT 2 IS NOT A FAILURE TO RETRY. `genomic_toolkit` returns 2 when the
   input's own shape makes the answer wrong -- contig naming, record order --
   and it prints a valid report first. `InvalidInput` carries that report.

MEASUREMENT
-----------
Wall time and peak RSS come back with every call because a lab evaluating this
engine is evaluating exactly those two numbers, and asking them to re-run
under `/usr/bin/time` to get them is friction with no purpose.

Peak RSS comes from `os.wait4`, NOT from `resource.getrusage(RUSAGE_CHILDREN)`
as `benchmarks/_measure.py` uses. That file's header explains why it must
interpose a whole process per measurement: `RUSAGE_CHILDREN.ru_maxrss` is the
high-water mark over every child the calling process has ever reaped, so a
second call in the same interpreter reports the first call's peak. A library
cannot fork a helper per call and must not silently return that contaminated
number, and `wait4` is scoped to one child by construction -- it is the same
kernel counter, attributed correctly. The KILOBYTES-on-Linux,
BYTES-on-macOS difference is normalised here exactly as it is there.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
import threading
import time
from dataclasses import dataclass, field
from pathlib import Path

from . import _binaries
from .errors import InvalidInput, ToolError

__all__ = ["Result", "run_tool"]

# ru_maxrss units. One branch, here, mirroring benchmarks/_measure.py.
_RSS_SCALE = 1 if sys.platform == "darwin" else 1024


@dataclass(frozen=True)
class Result:
    """What one pto-core invocation produced and what it cost.

    `metrics` is the tool's own JSON report, parsed -- the same object the
    cloud API stores on a job and the same one the CLI prints under `--json`.
    It is NOT reshaped into a schema of this package's own devising: the tools
    own their report fields, a lab's downstream script reads them by name, and
    a translation layer here would be a second schema to keep in step with
    three binaries.
    """

    tool: str
    subcommand: str | None
    argv: list[str]
    returncode: int
    metrics: dict = field(default_factory=dict)
    artifacts: dict = field(default_factory=dict)   # role -> Path
    wall_seconds: float = 0.0
    peak_rss_bytes: int | None = None
    stderr: str = ""

    @property
    def peak_rss_mb(self) -> float | None:
        if self.peak_rss_bytes is None:
            return None
        return round(self.peak_rss_bytes / 1e6, 1)

    def __getitem__(self, key: str):
        """`result["frip"]` -- the metrics are the point of the object."""
        return self.metrics[key]

    def get(self, key: str, default=None):
        return self.metrics.get(key, default)

    def as_dict(self) -> dict:
        """A JSON-serialisable record of the whole run, for a lab notebook."""
        return {
            "tool": self.tool,
            "subcommand": self.subcommand,
            "argv": list(self.argv),
            "returncode": self.returncode,
            "metrics": self.metrics,
            "artifacts": {k: str(v) for k, v in self.artifacts.items()},
            "wall_seconds": round(self.wall_seconds, 4),
            "peak_rss_bytes": self.peak_rss_bytes,
        }


def _wait_measured(proc: subprocess.Popen) -> tuple[int, int | None]:
    """Reap `proc`, returning (exit code, peak RSS bytes or None).

    `os.wait4` is unavailable on some platforms; there the exit code is still
    exact and peak RSS is reported as unknown rather than as a fabricated 0.
    """
    if not hasattr(os, "wait4"):
        return proc.wait(), None

    _, status, usage = os.wait4(proc.pid, 0)
    # Tell Popen the child is already reaped, or it waits on a pid that no
    # longer exists (and, on interpreter shutdown, warns about it).
    code = os.waitstatus_to_exitcode(status)
    proc.returncode = code
    return code, max(usage.ru_maxrss, 0) * _RSS_SCALE


def run_tool(
    tool: str,
    argv_tail: list[str],
    *,
    subcommand: str | None = None,
    artifacts: dict | None = None,
    stdout_json: bool = False,
    env: dict | None = None,
    timeout: float | None = None,
    json_path: Path | None = None,
    check: bool = True,
) -> Result:
    """Execute one binary.

    `stdout_json` is `genomic_toolkit`'s shape (a JSON report on stdout);
    `json_path` is `fastq_stream`'s (`-j PATH`). Both exist because the two
    tools genuinely differ, which `job_specs.py` records as well -- the runner
    there captures stdout for one and collects a file for the other.
    """
    binary = _binaries.resolve(tool)
    argv = [str(binary)] + ([subcommand] if subcommand else []) + [str(a) for a in argv_tail]

    full_env = {**os.environ, **(env or {})}

    # Files rather than pipes: `wait4` has to be the thing that reaps the
    # child, so there is no `communicate()` draining pipes concurrently, and a
    # tool that filled a pipe buffer would deadlock against our own wait.
    with tempfile.TemporaryDirectory(prefix="pto-run-") as scratch:
        out_path = Path(scratch) / "stdout"
        err_path = Path(scratch) / "stderr"
        started = time.perf_counter()
        with out_path.open("wb") as out_fh, err_path.open("wb") as err_fh:
            proc = subprocess.Popen(
                argv, stdout=out_fh, stderr=err_fh, env=full_env, shell=False,
            )
            timer: threading.Timer | None = None
            if timeout is not None:
                timer = threading.Timer(timeout, proc.kill)
                timer.daemon = True
                timer.start()
            try:
                returncode, peak_rss = _wait_measured(proc)
            finally:
                if timer is not None:
                    timer.cancel()
        wall = time.perf_counter() - started

        stdout_text = out_path.read_text(errors="replace")
        stderr_text = err_path.read_text(errors="replace")

    metrics: dict = {}
    if stdout_json and stdout_text.strip():
        try:
            metrics = json.loads(stdout_text)
        except json.JSONDecodeError:
            # A truncated report is worth surfacing as-is rather than as a
            # JSONDecodeError from inside a QC helper; the exit code below is
            # what says whether the run was sound.
            metrics = {"_unparsed_stdout": stdout_text}
    elif json_path is not None and json_path.is_file():
        try:
            metrics = json.loads(json_path.read_text())
        except json.JSONDecodeError:
            metrics = {}

    result = Result(
        tool=tool,
        subcommand=subcommand,
        argv=argv,
        returncode=returncode,
        metrics=metrics,
        artifacts={k: Path(v) for k, v in (artifacts or {}).items()},
        wall_seconds=wall,
        peak_rss_bytes=peak_rss,
        stderr=stderr_text,
    )

    if check and returncode != 0:
        if returncode == 2:
            raise InvalidInput(tool, argv, stderr_text, result=result)
        raise ToolError(tool, argv, returncode, stderr_text)
    return result
