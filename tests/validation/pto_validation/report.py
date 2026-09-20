# SPDX-License-Identifier: MIT
"""Result model and reporting.

The governing rule, from the brief and from the repo's own audit documents: a
discrepancy produces a *diff report* you can act on, never a bare assertion
failure. A metric that misses its threshold is written out with its inputs
alongside it -- the two curves, the disagreeing intervals, the counter that
drifted -- so the next person starts from evidence rather than from a
reproduction.

Metrics that pass are reported just as loudly. A validation suite that only
prints on failure teaches you nothing about the margin you are operating on,
and the margin is what tells you a regression is coming before it lands.
"""

from __future__ import annotations

import json
import platform
import sys
import time
from dataclasses import asdict, dataclass, field
from enum import Enum
from pathlib import Path


class Status(str, Enum):
    PASS = "pass"
    FAIL = "fail"
    SKIP = "skip"
    ERROR = "error"


class Direction(str, Enum):
    """Which side of the threshold is acceptable."""

    MIN = "min"  # value must be >= threshold (correlations, Jaccard, recall)
    MAX = "max"  # value must be <= threshold (KS distance, relative error)


@dataclass
class Metric:
    name: str
    value: float | None
    threshold: float | None
    direction: Direction
    unit: str = ""
    #: Set False for numbers that are reported for context but never gate the
    #: run -- observed duplicate-rate gaps against Picard, timings, counts.
    required: bool = True
    detail: str = ""

    @property
    def status(self) -> Status:
        if self.value is None:
            return Status.FAIL if self.required else Status.SKIP
        if self.threshold is None:
            return Status.PASS
        if self.direction is Direction.MIN:
            ok = self.value >= self.threshold
        else:
            ok = self.value <= self.threshold
        if ok:
            return Status.PASS
        return Status.FAIL if self.required else Status.SKIP

    def render(self) -> str:
        if self.value is None:
            value = "n/a"
        else:
            value = f"{self.value:.6g}{self.unit}"
        if self.threshold is None:
            bound = ""
        else:
            symbol = ">=" if self.direction is Direction.MIN else "<="
            bound = f"  (need {symbol} {self.threshold:g})"
        flag = "" if self.required else "  [informational]"
        return f"{value}{bound}{flag}"


@dataclass
class CaseResult:
    case: str
    dataset: str
    status: Status
    metrics: list[Metric] = field(default_factory=list)
    #: Human-readable reason for SKIP/ERROR. Always populated for those.
    message: str = ""
    #: Paths to artefacts written during the run (matrices, diffs, logs).
    artifacts: dict[str, str] = field(default_factory=dict)
    wall_seconds: float = 0.0
    reference_tools: dict[str, str] = field(default_factory=dict)

    @classmethod
    def skipped(cls, case: str, dataset: str, message: str) -> "CaseResult":
        return cls(case=case, dataset=dataset, status=Status.SKIP, message=message)

    @classmethod
    def errored(cls, case: str, dataset: str, message: str) -> "CaseResult":
        return cls(case=case, dataset=dataset, status=Status.ERROR, message=message)

    def finalize(self) -> "CaseResult":
        """Derive the case status from its metrics.

        A case with no metrics at all is an ERROR, not a PASS. Silently
        succeeding because every comparison was skipped is the single most
        dangerous outcome a validation harness can produce, and it is the
        outcome you get by default if you only check for failures.
        """
        if self.status in (Status.SKIP, Status.ERROR):
            return self
        if not self.metrics:
            self.status = Status.ERROR
            self.message = "case produced no metrics"
            return self
        if any(m.status is Status.FAIL for m in self.metrics):
            self.status = Status.FAIL
        else:
            self.status = Status.PASS
        return self


@dataclass
class Report:
    results: list[CaseResult] = field(default_factory=list)
    started_at: float = field(default_factory=time.time)
    output_dir: Path | None = None

    def add(self, result: CaseResult) -> CaseResult:
        self.results.append(result.finalize())
        return self.results[-1]

    def counts(self) -> dict[str, int]:
        counts = {status.value: 0 for status in Status}
        for result in self.results:
            counts[result.status.value] += 1
        return counts

    @property
    def ok(self) -> bool:
        """True when nothing failed or errored. Skips do not fail the run.

        Skips are the expected state on an air-gapped box with no deeptools,
        and treating them as failures would train everyone to ignore the exit
        code. The summary prints the skip count prominently instead, so a run
        that quietly validated nothing is visible rather than green.
        """
        return not any(
            r.status in (Status.FAIL, Status.ERROR) for r in self.results
        )

    # -- rendering ---------------------------------------------------------

    def write_json(self, path: Path) -> None:
        payload = {
            "schema": 1,
            "started_at": self.started_at,
            "duration_seconds": time.time() - self.started_at,
            "host": {
                "platform": platform.platform(),
                "python": platform.python_version(),
                "machine": platform.machine(),
            },
            "counts": self.counts(),
            "ok": self.ok,
            "results": [
                {**asdict(r), "status": r.status.value,
                 "metrics": [
                     {**asdict(m), "direction": m.direction.value,
                      "status": m.status.value}
                     for m in r.metrics
                 ]}
                for r in self.results
            ],
        }
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(payload, indent=2, default=str) + "\n")

    def write_markdown(self, path: Path) -> None:
        lines = [
            "# pto-core public-data validation",
            "",
            f"- host: `{platform.platform()}`",
            f"- generated: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(self.started_at))}",
            f"- duration: {time.time() - self.started_at:.1f}s",
            "",
            "| case | dataset | status | metrics |",
            "|---|---|---|---|",
        ]
        for result in self.results:
            summary = ", ".join(
                f"{m.name}={'n/a' if m.value is None else format(m.value, '.4g')}"
                for m in result.metrics
            ) or result.message
            lines.append(
                f"| `{result.case}` | `{result.dataset}` | "
                f"{_badge(result.status)} | {summary} |"
            )

        for result in self.results:
            if result.status is Status.PASS and not result.metrics:
                continue
            lines += ["", f"## {result.case} — {result.dataset}", ""]
            if result.message:
                lines += [f"> {result.message}", ""]
            if result.reference_tools:
                tools = ", ".join(f"`{k}`" for k in sorted(result.reference_tools))
                lines += [f"Compared against: {tools}", ""]
            if result.metrics:
                lines += ["| metric | value | threshold | status |", "|---|---|---|---|"]
                for m in result.metrics:
                    symbol = ">=" if m.direction is Direction.MIN else "<="
                    bound = "—" if m.threshold is None else f"{symbol} {m.threshold:g}"
                    value = "n/a" if m.value is None else f"{m.value:.6g}{m.unit}"
                    lines.append(f"| {m.name} | {value} | {bound} | {_badge(m.status)} |")
                lines.append("")
            for name, artifact in sorted(result.artifacts.items()):
                lines.append(f"- {name}: `{artifact}`")

        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("\n".join(lines) + "\n")

    def print_summary(self, stream=sys.stdout) -> None:
        width = max((len(r.case) for r in self.results), default=10) + 2
        print("\n" + "=" * 74, file=stream)
        print("VALIDATION SUMMARY", file=stream)
        print("=" * 74, file=stream)
        for result in self.results:
            print(
                f"  {_symbol(result.status)} {result.case:<{width}} "
                f"{result.dataset:<24} {result.status.value}",
                file=stream,
            )
            for metric in result.metrics:
                print(
                    f"      {_symbol(metric.status)} {metric.name:<34} {metric.render()}",
                    file=stream,
                )
            if result.message:
                for line in result.message.splitlines():
                    print(f"        {line}", file=stream)

        counts = self.counts()
        print("-" * 74, file=stream)
        print(
            f"  {counts['pass']} passed, {counts['fail']} failed, "
            f"{counts['error']} errored, {counts['skip']} skipped",
            file=stream,
        )
        if counts["skip"] and not counts["fail"]:
            print(
                "  note: skipped cases validated nothing. Install the missing "
                "reference\n        tools or stage the datasets to widen coverage.",
                file=stream,
            )
        print("=" * 74, file=stream)


def _badge(status: Status) -> str:
    return {
        Status.PASS: "✅ pass",
        Status.FAIL: "❌ **fail**",
        Status.SKIP: "⏭️ skip",
        Status.ERROR: "💥 error",
    }[status]


def _symbol(status: Status) -> str:
    return {Status.PASS: "PASS", Status.FAIL: "FAIL", Status.SKIP: "SKIP", Status.ERROR: "ERR "}[
        status
    ]


# ---------------------------------------------------------------------------
# Diff reports
# ---------------------------------------------------------------------------


def write_curve_diff(
    path: Path,
    label_a: str,
    label_b: str,
    xs: list[float],
    a: list[float],
    b: list[float],
) -> Path:
    """Write a per-position comparison of two signal curves.

    Ordered so the largest absolute divergence is trivially greppable, and
    carrying both the absolute and relative gap: a 0.01 difference means
    something entirely different at the flanks of a TSS window than it does at
    the summit, and only the relative column separates the two.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    for x, va, vb in zip(xs, a, b):
        delta = va - vb
        rel = abs(delta) / abs(vb) if vb else float("nan")
        rows.append((x, va, vb, delta, rel))

    with path.open("w") as handle:
        handle.write(f"# largest divergence first\noffset\t{label_a}\t{label_b}\tdelta\trel_delta\n")
        for x, va, vb, delta, rel in sorted(rows, key=lambda r: -abs(r[3])):
            handle.write(f"{x:g}\t{va:.6g}\t{vb:.6g}\t{delta:+.6g}\t{rel:.6g}\n")
    return path


def write_interval_diff(
    path: Path,
    label_a: str,
    label_b: str,
    a: list[tuple[str, int, int]],
    b: list[tuple[str, int, int]],
    *,
    limit: int = 5000,
) -> Path:
    """Write the intervals present in exactly one of two sets.

    Capped, because a naming-convention mismatch ('chr1' vs '1') puts *every*
    interval on both sides of the diff and a 3-million-line report is not a
    diagnostic. The cap is stated in the file so a truncated report is never
    mistaken for a complete one.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    sa, sb = set(a), set(b)
    only_a = sorted(sa - sb)
    only_b = sorted(sb - sa)

    with path.open("w") as handle:
        handle.write(
            f"# intervals unique to one side\n"
            f"# {label_a}-only: {len(only_a)}   {label_b}-only: {len(only_b)}"
            f"   shared: {len(sa & sb)}\n"
        )
        if len(only_a) > limit or len(only_b) > limit:
            handle.write(f"# TRUNCATED to {limit} rows per side\n")
        handle.write("side\tcontig\tstart\tend\n")
        for contig, start, end in only_a[:limit]:
            handle.write(f"{label_a}\t{contig}\t{start}\t{end}\n")
        for contig, start, end in only_b[:limit]:
            handle.write(f"{label_b}\t{contig}\t{start}\t{end}\n")
    return path


def write_histogram_diff(
    path: Path,
    label_a: str,
    label_b: str,
    a: dict[int, float],
    b: dict[int, float],
) -> Path:
    """Write two size distributions side by side as normalised fractions.

    Normalised, because the two tools filter differently upstream (MAPQ gates,
    proper-pair requirements) and comparing raw counts would report a filtering
    difference as a shape difference. The shape is what the KS statistic tests,
    so the diff has to show the same thing.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    total_a = sum(a.values()) or 1.0
    total_b = sum(b.values()) or 1.0

    with path.open("w") as handle:
        handle.write(f"size\t{label_a}_frac\t{label_b}_frac\tdelta\n")
        for key in sorted(set(a) | set(b)):
            fa = a.get(key, 0.0) / total_a
            fb = b.get(key, 0.0) / total_b
            handle.write(f"{key}\t{fa:.8g}\t{fb:.8g}\t{fa - fb:+.8g}\n")
    return path
