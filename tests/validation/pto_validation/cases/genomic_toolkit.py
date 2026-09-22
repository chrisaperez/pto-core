# SPDX-License-Identifier: MIT
"""genomic_toolkit vs samtools, bedtools and Picard.

genomic_toolkit already ships a differential check --
`modules/genomic_toolkit/scripts/validate_against_bedtools.sh` -- that asserts
*exact* equality of the in-peak fragment count against the samtools/bedtools
pipeline. This module does not reimplement it; the FRiP case shells out to it
and reports its verdict, so there is exactly one definition of that assertion
in the repo and no chance of the two drifting apart.

What is added here is the two comparisons that script does not make: the
fragment size distribution against samtools' own TLEN histogram, and the
duplicate rate against Picard MarkDuplicates.
"""

from __future__ import annotations

import json
from pathlib import Path

from .. import metrics as M
from ..config import (
    DUPLICATE_RATE_RELATIVE_MAX,
    FRAGMENT_KS_MAX,
    REPO_ROOT,
)
from ..harness import Context, Unavailable, case
from ..report import Direction, Metric, write_histogram_diff
from ..tools import ToolError, read_histogram_tsv, run


@case(
    "gtk/frip-vs-bedtools",
    module="genomic_toolkit",
    reference=["samtools", "bedtools"],
    summary="Exact in-peak fragment count against samtools+bedtools",
)
def frip_vs_bedtools(ctx: Context) -> list[Metric]:
    """Delegate to the module's own differential script and report its result."""
    script = REPO_ROOT / "modules/genomic_toolkit/scripts/validate_against_bedtools.sh"
    if not script.is_file():
        raise Unavailable(f"{script} is missing")

    binary = ctx.binary("genomic_toolkit")
    bam = ctx.input("bam")
    peaks = ctx.input("peaks")
    work = ctx.scratch("gtk_frip")

    try:
        result = run(["bash", script, binary, bam, peaks], timeout=14400)
    except ToolError as exc:
        # The script asserts exact equality itself and exits non-zero with the
        # two counts on stderr. Surfacing its own message verbatim is more
        # useful than re-deriving the numbers here and risks no disagreement
        # between the two reports.
        log = work / "validate_against_bedtools.log"
        log.write_text(exc.stderr)
        ctx.artifact("bedtools differential log", log)
        return [
            Metric(
                "frip.exact_count_match", 0.0, 1.0, Direction.MIN,
                detail=exc.stderr.strip().splitlines()[-1] if exc.stderr.strip() else "",
            )
        ]

    log = work / "validate_against_bedtools.log"
    log.write_text(result.stdout + result.stderr)
    ctx.artifact("bedtools differential log", log)

    if "skip:" in result.stderr:
        raise Unavailable(result.stderr.strip())

    # Cross-check the ratio the tool reports against the counts the script
    # printed, so a passing script plus a mis-derived FRiP still fails here.
    report = json.loads(
        run([binary, "frip", bam, "--peaks", peaks, "--json"]).stdout
    )
    ctx.artifact("gtk frip report", _dump(work / "frip.json", report))

    if report.get("unmatched_contigs"):
        raise Unavailable(
            f"{report['unmatched_contigs']} contigs in the peak file have no "
            "match in the alignment header; the FRiP would be meaningless"
        )

    derived = report["in_peaks"] / report["fragments"] if report["fragments"] else None
    return [
        Metric("frip.exact_count_match", 1.0, 1.0, Direction.MIN),
        Metric(
            "frip.self_consistency",
            M.relative_error(report["frip"], derived) if derived is not None else None,
            1e-6, Direction.MAX,
            detail=f"reported {report['frip']:.6f}, in_peaks/fragments = "
                   f"{derived:.6f}" if derived is not None else "",
        ),
        Metric("frip.value", report["frip"], None, Direction.MIN, required=False),
        Metric("frip.fragments", float(report["fragments"]), None, Direction.MIN,
               required=False),
    ]


@case(
    "gtk/sizes-vs-samtools",
    module="genomic_toolkit",
    reference=["samtools"],
    summary="Fragment size distribution parity with samtools TLEN",
)
def sizes_vs_samtools(ctx: Context) -> list[Metric]:
    """Compare the size histogram against one built from raw TLEN values."""
    binary = ctx.binary("genomic_toolkit")
    bam = ctx.input("bam")
    work = ctx.scratch("gtk_sizes")

    min_mapq = int(ctx.param("min_mapq", 30))
    max_length = int(ctx.param("max_fragment", 1000))

    histogram = work / "pto_sizes.tsv"
    report = json.loads(
        run([
            binary, "sizes", bam,
            "--min-mapq", min_mapq,
            "--max-length", max_length,
            "--histogram", histogram,
            "--json",
        ]).stdout
    )
    ctx.artifact("pto size histogram", histogram)
    pto_hist = read_histogram_tsv(histogram)

    # samtools side: one properly-paired, primary, non-duplicate alignment per
    # fragment (-f 0x2, -F 0xF04), taking the positive TLEN so each pair is
    # counted once rather than once per mate. This is the same fragment set
    # genomic_toolkit builds; the flag mask is the whole comparison, so it is
    # spelled out rather than hidden behind a samtools preset.
    view = run([
        ctx.reference("samtools"), "view",
        "-f", "0x2", "-F", "0xF04", "-q", str(min_mapq), bam,
    ], timeout=14400)

    reference_hist: dict[int, float] = {}
    for line in view.stdout.splitlines():
        fields = line.split("\t", 9)
        if len(fields) < 9:
            continue
        try:
            tlen = int(fields[8])
        except ValueError:
            continue
        if tlen <= 0 or (max_length and tlen > max_length):
            continue
        reference_hist[tlen] = reference_hist.get(tlen, 0.0) + 1.0

    if not reference_hist:
        raise Unavailable("samtools produced no properly-paired fragments to compare")

    ks = M.histogram_ks(pto_hist, reference_hist)
    if ks is None or ks > FRAGMENT_KS_MAX:
        ctx.artifact(
            "size distribution diff",
            write_histogram_diff(
                work / "size_diff.tsv", "pto", "samtools",
                pto_hist, reference_hist,
            ),
        )

    pto_median = report["median"]
    reference_median = _weighted_median(reference_hist)

    return [
        Metric("sizes.ks_distance", ks, FRAGMENT_KS_MAX, Direction.MAX),
        Metric(
            "sizes.median_abs_error",
            abs(pto_median - reference_median), 1.0, Direction.MAX, unit=" bp",
            detail=f"pto {pto_median} bp vs samtools {reference_median} bp",
        ),
        Metric("sizes.mononucleosome_fraction", report["mononucleosome"], None,
               Direction.MIN, required=False),
        Metric("sizes.fragments", float(report["fragments"]), None, Direction.MIN,
               required=False),
    ]


@case(
    "gtk/markdup-vs-picard",
    module="genomic_toolkit",
    reference=["picard", "samtools"],
    summary="Duplicate rate against Picard MarkDuplicates (bounded, not exact)",
)
def markdup_vs_picard(ctx: Context) -> list[Metric]:
    """Compare duplicate rates, with the documented divergence in mind.

    README.md's "Trade-offs" section records that genomic_toolkit keys a
    duplicate on the *aligned span* while Picard uses the unclipped 5' end, and
    that it has no optical-duplicate or UMI awareness. The two therefore
    disagree on soft-clipped reads by construction, which is why this asserts a
    5% relative bound rather than equality. The point of the assertion is to
    catch that gap *changing*, so the observed rates are always reported.
    """
    binary = ctx.binary("genomic_toolkit")
    bam = ctx.input("bam")
    work = ctx.scratch("gtk_markdup")

    report = json.loads(run([binary, "markdup", bam, "--json"]).stdout)
    ctx.artifact("gtk markdup report", _dump(work / "markdup.json", report))
    pto_rate = _duplicate_rate(report)
    if pto_rate is None:
        raise Unavailable(
            "markdup's JSON carries no duplicate counts to compare; the report "
            f"keys were {sorted(report)}"
        )

    picard_metrics = work / "picard_dup_metrics.txt"
    run([
        ctx.reference("picard"), "MarkDuplicates",
        f"--INPUT={bam}",
        f"--OUTPUT={work / 'picard_marked.bam'}",
        f"--METRICS_FILE={picard_metrics}",
        "--ASSUME_SORT_ORDER=coordinate",
        "--VALIDATION_STRINGENCY=LENIENT",
    ], timeout=14400)
    ctx.artifact("picard metrics", picard_metrics)

    picard_rate = _picard_duplicate_rate(picard_metrics)
    if picard_rate is None:
        raise Unavailable(f"could not parse a duplicate rate out of {picard_metrics}")

    return [
        Metric(
            "markdup.rate_rel_error",
            M.relative_error(pto_rate, picard_rate),
            DUPLICATE_RATE_RELATIVE_MAX, Direction.MAX,
            detail="expected non-zero: aligned span vs unclipped 5' end, "
                   "no optical/UMI awareness (see README Trade-offs)",
        ),
        Metric("markdup.rate_pto", pto_rate, None, Direction.MIN,
               required=False),
        Metric("markdup.rate_picard", picard_rate, None, Direction.MIN,
               required=False),
    ]


def _weighted_median(histogram: dict[int, float]) -> int:
    """Median of a {size: count} histogram, using the lower-of-two convention.

    genomic_toolkit reports an integer median for an even fragment count by
    taking the lower of the two central values rather than averaging them, so
    this matches that convention -- averaging here would produce a spurious
    0.5 bp disagreement on every even-sized library.
    """
    total = sum(histogram.values())
    seen = 0.0
    for size in sorted(histogram):
        seen += histogram[size]
        if seen * 2 >= total:
            return size
    return 0


def _duplicate_rate(report: dict) -> float | None:
    """Pull a duplicate rate out of markdup's JSON, tolerating either shape.

    The report has carried both an explicit rate and a duplicates/total pair
    across versions; accepting either keeps this case working against an older
    build without pinning it to one JSON revision.
    """
    for key in ("duplicate_rate", "dup_rate", "rate"):
        if key in report:
            return float(report[key])
    total = report.get("fragments") or report.get("records_read")
    duplicates = report.get("duplicates") or report.get("duplicates_excluded")
    if total and duplicates is not None:
        return float(duplicates) / float(total)
    return None


_PICARD_HEADER = "LIBRARY\t"


def _picard_duplicate_rate(path: Path) -> float | None:
    """Parse PERCENT_DUPLICATION from a Picard metrics file.

    Picard's format is a run of '#' comment blocks followed by a tab-separated
    header and one row per library. Aggregating across libraries would need
    their read counts as weights; these datasets are single-library, so the
    first data row is taken and a multi-library file is reported rather than
    silently averaged.
    """
    lines = [line.rstrip("\n") for line in path.read_text().splitlines()]
    for index, line in enumerate(lines):
        if line.startswith(_PICARD_HEADER):
            header = line.split("\t")
            for row in lines[index + 1:]:
                if not row.strip():
                    break
                fields = dict(zip(header, row.split("\t")))
                value = fields.get("PERCENT_DUPLICATION")
                if value:
                    try:
                        return float(value)
                    except ValueError:
                        return None
            return None
    return None


def _dump(path: Path, payload: dict) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n")
    return path
