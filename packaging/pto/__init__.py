"""pto -- Python bindings for pto-core's genomics engines.

Five C++20 tools that run on the machine holding the sequencing data: no
service dependency, no data leaving the host. This package wraps four of them
(the fifth, `scrna_matrix`, is a compiled extension module distributed
separately as `scrna-matrix`).

    import pto

    qc = pto.run_qc("sample_R1.fastq.gz")
    print(qc["summary"]["q30_rate"], qc.wall_seconds, qc.peak_rss_mb)

    sizes = pto.fragment_sizes("atac.bam")
    called = pto.call_peaks("atac.bam", "peaks.narrowPeak")
    enrichment = pto.frip("atac.bam", called.artifacts["peaks"])
    print(enrichment["frip"], sizes["median"])

Every call returns a `Result`: the tool's own JSON report under `.metrics`
(indexable directly), the artifacts it wrote under `.artifacts`, and what the
run cost under `.wall_seconds` / `.peak_rss_mb`.

If a binary cannot be found, `pto doctor` says where it looked.
"""

from __future__ import annotations

from ._binaries import TOOLS, bin_dir, resolve, which
from ._cpu import available_cpus
from ._run import Result, run_tool
from .api import (
    call_peaks,
    fragment_sizes,
    frip,
    mark_duplicates,
    profile,
    run_qc,
    trim,
)
from .errors import BinaryNotFound, InvalidInput, PtoError, ToolError

__all__ = [
    "BinaryNotFound",
    "InvalidInput",
    "PtoError",
    "Result",
    "TOOLS",
    "ToolError",
    "available_cpus",
    "bin_dir",
    "call_peaks",
    "fragment_sizes",
    "frip",
    "mark_duplicates",
    "profile",
    "resolve",
    "run_qc",
    "run_tool",
    "tool_versions",
    "trim",
    "which",
    "__version__",
]


def _package_version() -> str:
    try:
        from importlib.metadata import PackageNotFoundError, version
    except ImportError:  # pragma: no cover - Python < 3.8
        return "0.0.0+unknown"
    try:
        return version("pto-core")
    except PackageNotFoundError:
        # Running from a source checkout with no install. Not an error; the
        # binaries are what matter and `tool_versions()` reports those.
        return "0.0.0+source"


__version__ = _package_version()


def tool_versions() -> dict:
    """`{tool: version string or None}` for each wrapped binary.

    The distribution's own `__version__` says which wrapper is installed; this
    says which *engines* it found, which is the number that belongs in a
    methods section. A tool that is absent maps to None rather than raising, so
    this stays usable as a diagnostic when only some binaries are present.
    """
    import subprocess

    versions: dict[str, str | None] = {}
    for tool in TOOLS:
        path = which(tool)
        if path is None:
            versions[tool] = None
            continue
        try:
            proc = subprocess.run(
                [str(path), "--version"], capture_output=True, text=True,
                timeout=30, shell=False, check=False,
            )
            out = (proc.stdout + proc.stderr).strip().splitlines()
            versions[tool] = out[0].strip() if out else None
        except (OSError, subprocess.SubprocessError):
            versions[tool] = None
    return versions
