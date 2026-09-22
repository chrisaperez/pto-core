# SPDX-License-Identifier: MIT
"""The case registry and the context object handed to every case.

A "case" is one comparison of one pto module against one gold-standard
pipeline on one dataset. Cases register themselves with `@case(...)`, declare
which binaries and reference tools they need, and are handed a `Context` that
has already resolved paths, staged inputs, and created a scratch directory.

Cases raise `Unavailable` to skip. They do not call sys.exit, do not assert,
and do not print -- everything they learn goes back as Metric objects, so the
same code path serves the terminal summary, the JSON artefact and CI.
"""

from __future__ import annotations

import traceback
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Iterable

from . import fetch
from .config import Dataset, find_pto_binary
from .report import CaseResult, Metric, Report, Status


class Unavailable(Exception):
    """Raised by a case that cannot run here. Becomes a SKIP, never a failure."""


@dataclass
class Context:
    dataset: Dataset
    data_dir: Path
    build_dir: Path
    work_dir: Path
    threads: int
    allow_download: bool
    reference_tools: dict[str, str | None]
    #: Populated by the case as it goes; surfaces in the report even on failure.
    artifacts: dict[str, str] = field(default_factory=dict)

    # -- resource resolution ----------------------------------------------

    def binary(self, module: str) -> Path:
        path = find_pto_binary(module, self.build_dir)
        if path is None:
            raise Unavailable(
                f"{module} binary not found under {self.build_dir}. Build it:\n"
                f"  cmake -S . -B {self.build_dir} -DCMAKE_BUILD_TYPE=Release && "
                f"cmake --build {self.build_dir} -j"
            )
        return path

    def reference(self, name: str) -> str:
        path = self.reference_tools.get(name)
        if not path:
            raise Unavailable(f"reference tool {name!r} is not on PATH")
        return path

    def input(self, role: str) -> Path:
        """Stage and return one of the dataset's files.

        Wraps the download-disabled case into `Unavailable` so an air-gapped
        run reports a tidy skip with staging instructions instead of a
        traceback -- that host is the primary deployment target, not an edge
        case, and its output should read like a normal outcome.
        """
        if role not in self.dataset.files:
            raise Unavailable(
                f"dataset {self.dataset.key!r} declares no file for role {role!r}"
            )
        try:
            return fetch.ensure_file(
                self.dataset, role, self.data_dir, allow_download=self.allow_download
            )
        except fetch.DownloadsDisabled as exc:
            raise Unavailable(str(exc)) from exc
        except fetch.FetchError as exc:
            raise Unavailable(f"could not stage {role}: {exc}") from exc

    def param(self, name: str, default=None):
        return self.dataset.params.get(name, default)

    def scratch(self, name: str) -> Path:
        path = self.work_dir / name
        path.mkdir(parents=True, exist_ok=True)
        return path

    def artifact(self, label: str, path: Path) -> Path:
        self.artifacts[label] = str(path)
        return path


@dataclass
class Case:
    name: str
    module: str
    reference: tuple[str, ...]
    summary: str
    function: Callable[[Context], list[Metric]]


REGISTRY: dict[str, Case] = {}


def case(name: str, *, module: str, reference: Iterable[str] = (), summary: str = ""):
    """Register a validation case.

    `reference` is checked before the case body runs, so a missing deeptools
    skips with a precise message instead of failing three subprocess calls
    deep with whatever error that tool's absence happens to produce.
    """

    def decorate(function: Callable[[Context], list[Metric]]) -> Callable:
        if name in REGISTRY:
            raise RuntimeError(f"duplicate validation case {name!r}")
        REGISTRY[name] = Case(
            name=name,
            module=module,
            reference=tuple(reference),
            summary=summary or (function.__doc__ or "").strip().splitlines()[0],
            function=function,
        )
        return function

    return decorate


def run_case(spec: Case, context: Context, report: Report) -> CaseResult:
    """Execute one case, converting every outcome into a CaseResult.

    Nothing escapes: an Unavailable becomes a skip, any other exception becomes
    an ERROR carrying the traceback. A case that crashes must not take the
    remaining cases down with it -- a five-hour suite that aborts on the first
    missing index file is a suite nobody runs twice.
    """
    import time

    started = time.perf_counter()
    missing = [name for name in spec.reference if not context.reference_tools.get(name)]
    if missing:
        return report.add(
            CaseResult.skipped(
                spec.name,
                context.dataset.key,
                "missing reference tool(s): " + ", ".join(missing),
            )
        )

    try:
        metrics = spec.function(context)
        result = CaseResult(
            case=spec.name,
            dataset=context.dataset.key,
            status=Status.PASS,
            metrics=metrics,
            artifacts=dict(context.artifacts),
            wall_seconds=time.perf_counter() - started,
            reference_tools={
                name: context.reference_tools[name] or ""
                for name in spec.reference
            },
        )
    except Unavailable as exc:
        result = CaseResult.skipped(spec.name, context.dataset.key, str(exc))
        result.artifacts = dict(context.artifacts)
    except Exception:  # noqa: BLE001 - deliberate catch-all, see docstring
        result = CaseResult.errored(
            spec.name, context.dataset.key, traceback.format_exc(limit=8)
        )
        result.artifacts = dict(context.artifacts)

    result.wall_seconds = time.perf_counter() - started
    return report.add(result)


def load_cases() -> dict[str, Case]:
    """Import the case modules for their registration side effects."""
    from .cases import (  # noqa: F401
        cuttag_profiler,
        fastq_stream,
        genomic_toolkit,
        scrna_matrix,
    )

    return REGISTRY
