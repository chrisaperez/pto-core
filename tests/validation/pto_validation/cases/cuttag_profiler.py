# SPDX-License-Identifier: MIT
"""cuttag_profiler vs deeptools (computeMatrix) and bedtools.

The comparison that matters for this module is the meta-profile: given the same
indexed BAM, the same reference-point window and the same bin size, does the
curve agree with what deeptools produces? Everything downstream of that -- the
heatmap, the top-regions table, the dashboard -- is a rendering of the same
matrix, so agreement on the curve is agreement on the tool.

Known, expected divergence (do not "fix" by loosening the threshold without
reading this): deeptools assigns a fragment to the bin containing its midpoint,
while cuttag_profiler's `reads` count mode adds it to every bin its span
overlaps. On a 50 bp bin against ~200 bp CUT&Tag fragments that is a mild
smoothing difference, worth ~0.005 of Pearson r and concentrated at the window
flanks where the signal is lowest. It is not a counting error in either tool.
"""

from __future__ import annotations

from pathlib import Path

from .. import metrics as M
from ..config import (
    PEAK_JACCARD_MIN,
    SIGNAL_PEARSON_MIN,
    SIGNAL_SPEARMAN_MIN,
)
from ..harness import Context, Unavailable, case
from ..report import Direction, Metric
from ..tools import (
    deeptools_matrix_profile,
    read_bed,
    read_profile_tsv,
    run,
)


def _require_index(ctx: Context, bam: Path) -> None:
    """Ensure a BAM index exists, building one with samtools if it does not.

    ENCODE distributes alignments and their indexes as separate accessions,
    and this suite pins only the alignment, so the index is normally absent
    after a fresh stage. Building it here rather than failing is worth it: it
    is a one-off cost per dataset, it is cached alongside the BAM, and the
    alternative is a htslib read error that says nothing about which of the
    two files is missing.
    """
    if any(
        bam.with_suffix(bam.suffix + suffix).exists() or bam.with_suffix(suffix).exists()
        for suffix in (".bai", ".csi", ".crai")
    ):
        return

    samtools = ctx.reference_tools.get("samtools")
    if not samtools:
        raise Unavailable(
            f"no .bai/.csi index alongside {bam.name} and samtools is not on "
            f"PATH to build one; index it with `samtools index {bam}`"
        )
    run([samtools, "index", "-@", ctx.threads, bam], timeout=14400)


@case(
    "cuttag/profile-vs-deeptools",
    module="cuttag_profiler",
    # bamCoverage is listed alongside computeMatrix because the comparison
    # needs both halves of the deeptools pipeline; a conda env with only one
    # of them installed should skip up front, not fail midway.
    reference=["computeMatrix", "bamCoverage", "samtools"],
    summary="TSS meta-profile agreement with deeptools computeMatrix",
)
def profile_vs_deeptools(ctx: Context) -> list[Metric]:
    """Compare the reference-point meta-profile against deeptools."""
    profiler = ctx.binary("cuttag_profiler")
    bam = ctx.input("bam")
    _require_index(ctx, bam)
    regions = ctx.input("regions")
    work = ctx.scratch("cuttag_deeptools")

    upstream = int(ctx.param("upstream", 2000))
    downstream = int(ctx.param("downstream", 2000))
    bin_size = int(ctx.param("bin_size", 50))
    min_mapq = int(ctx.param("min_mapq", 30))

    # pto side. `raw` normalization on both sides: CPM and deeptools' own
    # scaling are computed over different denominators (deeptools uses the
    # library size after its own filters), so comparing normalized curves would
    # test the two scale factors rather than the binning, and a scale offset is
    # exactly what Pearson r is blind to anyway.
    pto_profile = work / "pto_profile.tsv"
    run(
        [
            profiler, "profile",
            "--bam", bam,
            "--regions", regions,
            "--out-profile", pto_profile,
            "--reference-point", ctx.param("reference_point", "TSS"),
            "--upstream", upstream,
            "--downstream", downstream,
            "--bin-size", bin_size,
            "--normalization", "raw",
            "--count-mode", "reads",
            "--min-mapq", min_mapq,
            "--threads", ctx.threads,
            "--quiet",
        ]
    )
    ctx.artifact("pto profile", pto_profile)

    # deeptools side. --missingDataAsZero matches how cuttag_profiler treats an
    # uncovered bin; without it deeptools emits nan there and the two curves
    # would be compared over different supports.
    dt_matrix = work / "deeptools_matrix.gz"
    run(
        [
            ctx.reference("computeMatrix"), "reference-point",
            "--referencePoint", "TSS",
            "-S", _require_bigwig(ctx, work),
            "-R", regions,
            "--beforeRegionStartLength", upstream,
            "--afterRegionStartLength", downstream,
            "--binSize", bin_size,
            "--missingDataAsZero",
            "--numberOfProcessors", ctx.threads,
            "-o", dt_matrix,
        ],
        timeout=7200,
    )
    ctx.artifact("deeptools matrix", dt_matrix)

    offsets, pto_values = read_profile_tsv(pto_profile)
    deeptools_values = deeptools_matrix_profile(dt_matrix)

    if not pto_values or not deeptools_values:
        raise Unavailable("one side produced an empty profile")

    # Bin counts can differ by one when the window is not an exact multiple of
    # the bin size and the two tools round the last partial bin differently.
    # Truncating to the shared length compares the bins that both defined,
    # rather than silently pairing bin i against bin i+1 for the whole curve.
    n = min(len(pto_values), len(deeptools_values))
    pto_values, deeptools_values = pto_values[:n], deeptools_values[:n]
    offsets = offsets[:n]

    pearson = M.pearson(pto_values, deeptools_values)
    spearman = M.spearman(pto_values, deeptools_values)

    if pearson is None or pearson < SIGNAL_PEARSON_MIN:
        from ..report import write_curve_diff

        ctx.artifact(
            "profile diff",
            write_curve_diff(
                work / "profile_diff.tsv",
                "pto", "deeptools",
                offsets, pto_values, deeptools_values,
            ),
        )

    return [
        Metric("profile.pearson_r", pearson, SIGNAL_PEARSON_MIN, Direction.MIN),
        Metric("profile.spearman_rho", spearman, SIGNAL_SPEARMAN_MIN, Direction.MIN),
        Metric("profile.bins_compared", float(n), None, Direction.MIN, required=False),
        Metric(
            "profile.peak_offset_bp",
            float(offsets[pto_values.index(max(pto_values))]),
            None,
            Direction.MIN,
            unit=" bp",
            required=False,
            detail="offset of maximum mean signal; should sit at or near 0 for a TSS mark",
        ),
    ]


def _require_bigwig(ctx: Context, work: Path) -> Path:
    """Produce the bigWig computeMatrix needs, via deeptools bamCoverage.

    computeMatrix consumes coverage tracks, not alignments, so a BAM-to-BAM
    comparison is not available -- the conversion is part of the pipeline being
    replaced, and its cost is part of what the benchmark suite reports.
    Cached across cases in the same run: on a deep whole-genome BAM this single
    step outlasts every other operation in the suite combined.
    """
    bigwig = work.parent / "reference_coverage.bw"
    if bigwig.exists() and bigwig.stat().st_size > 0:
        return bigwig

    bam = ctx.input("bam")
    run(
        [
            ctx.reference("bamCoverage"),
            "-b", bam,
            "-o", bigwig,
            "--binSize", int(ctx.param("bin_size", 50)),
            "--minMappingQuality", int(ctx.param("min_mapq", 30)),
            "--extendReads",
            "--numberOfProcessors", ctx.threads,
        ],
        timeout=14400,
    )
    return bigwig


@case(
    "cuttag/enriched-regions-vs-bedtools",
    module="cuttag_profiler",
    reference=["bedtools", "samtools"],
    summary="Top-signal regions overlap the reference peak call",
)
def enriched_regions_vs_bedtools(ctx: Context) -> list[Metric]:
    """Check that the highest-signal windows land inside the published peaks.

    This is a different question from the profile comparison, and a weaker one
    on purpose: it does not ask whether cuttag_profiler reproduces a peak
    caller (it is not one), only whether the regions it ranks highest are the
    regions an established caller also called. A module that computed signal
    correctly but wired up its region ordering backwards would pass the
    profile check and fail this one.
    """
    profiler = ctx.binary("cuttag_profiler")
    bam = ctx.input("bam")
    _require_index(ctx, bam)
    regions = ctx.input("regions")
    peaks = ctx.input("peaks")
    work = ctx.scratch("cuttag_peaks")

    matrix = work / "matrix.tsv"
    run(
        [
            profiler, "profile",
            "--bam", bam,
            "--regions", regions,
            "--out-matrix", matrix,
            "--upstream", int(ctx.param("upstream", 2000)),
            "--downstream", int(ctx.param("downstream", 2000)),
            "--bin-size", int(ctx.param("bin_size", 50)),
            "--normalization", "CPM",
            "--min-mapq", int(ctx.param("min_mapq", 30)),
            "--threads", ctx.threads,
            "--quiet",
        ]
    )
    ctx.artifact("pto matrix", matrix)

    ranked = _top_regions_bed(matrix, work / "top_regions.bed", fraction=0.10)
    ctx.artifact("top regions", ranked)

    hits = run(
        [ctx.reference("bedtools"), "intersect", "-u", "-a", ranked, "-b", peaks]
    )
    overlapping = sum(1 for line in hits.stdout.splitlines() if line.strip())
    total = sum(1 for _ in read_bed(ranked))
    if total == 0:
        raise Unavailable("matrix produced no ranked regions")

    top_intervals = read_bed(ranked)
    peak_intervals = read_bed(peaks)
    bp_jaccard = M.interval_jaccard(top_intervals, peak_intervals)
    precision = overlapping / total

    if precision < PEAK_JACCARD_MIN:
        from ..report import write_interval_diff

        ctx.artifact(
            "region/peak diff",
            write_interval_diff(
                work / "region_peak_diff.tsv",
                "top_regions", "reference_peaks",
                top_intervals, peak_intervals,
            ),
        )

    return [
        Metric("peaks.top_region_precision", precision, PEAK_JACCARD_MIN, Direction.MIN),
        Metric(
            "peaks.bp_jaccard", bp_jaccard, None, Direction.MIN, required=False,
            detail="base-pair overlap of the top decile with the published peak set; "
                   "informational because window width and peak width differ",
        ),
        Metric("peaks.regions_ranked", float(total), None, Direction.MIN, required=False),
    ]


def _top_regions_bed(matrix: Path, destination: Path, *, fraction: float) -> Path:
    """Write the top `fraction` of matrix rows, by mean signal, as a BED.

    The matrix carries its region coordinates in the leading annotation
    columns, so this needs no second pass over the original BED and cannot
    drift out of sync with the row ordering the tool actually used.
    """
    scored: list[tuple[float, tuple[str, int, int]]] = []
    with matrix.open() as handle:
        for line in handle:
            if not line.strip() or line[0] == "#":
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                continue
            try:
                contig, start, end = fields[0], int(fields[1]), int(fields[2])
            except ValueError:
                continue  # header row
            values = []
            for field in fields[3:]:
                try:
                    values.append(float(field))
                except ValueError:
                    continue
            if not values:
                continue
            scored.append((sum(values) / len(values), (contig, start, end)))

    scored.sort(key=lambda item: -item[0])
    keep = max(1, int(len(scored) * fraction))
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w") as handle:
        for contig, start, end in sorted(item[1] for item in scored[:keep]):
            handle.write(f"{contig}\t{start}\t{end}\n")
    return destination
