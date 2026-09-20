# SPDX-License-Identifier: MIT
"""fastq_stream vs seqkit, fastp, and a bowtie2 round trip.

Three questions, in increasing order of how much they can tell you:

1. Does it count the same reads and bases as seqkit on the same file? These are
   integers off the same input, so the tolerance is exactly zero.
2. Does its QC agree with fastp's, and does its trimming remove the same
   adapter content? Trimmers legitimately differ in where they cut, so this
   compares the *effect* -- surviving reads and bases, per-cycle quality --
   rather than demanding byte-identical FASTQ.
3. Does an aligner accept its output and map at least as well as it maps the
   raw reads? That is the question the tool actually exists to answer, and it
   is the one a counter comparison cannot reach.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

from .. import metrics as M
from ..config import QC_COUNT_RELATIVE_MAX, QC_QUALITY_ABS_MAX
from ..harness import Context, Unavailable, case
from ..report import Direction, Metric
from ..tools import run


@case(
    "fastq/counts-vs-seqkit",
    module="fastq_stream",
    reference=["seqkit"],
    summary="Read/base/GC counters against seqkit stats",
)
def counts_vs_seqkit(ctx: Context) -> list[Metric]:
    """Assert exact agreement on read count, base count and GC rate."""
    binary = ctx.binary("fastq_stream")
    fastq = ctx.input("fastq_r1")
    work = ctx.scratch("fastq_seqkit")

    report = work / "qc.json"
    run([binary, "-i", fastq, "--qc-only", "-j", report, "--threads", ctx.threads])
    ctx.artifact("pto QC", report)
    qc = json.loads(report.read_text())["summary"]

    stats = run([ctx.reference("seqkit"), "stats", "-a", "-T", fastq])
    rows = [line.split("\t") for line in stats.stdout.strip().splitlines()]
    if len(rows) < 2:
        raise Unavailable("seqkit stats produced no data row")
    header, values = rows[0], rows[1]
    seqkit = dict(zip(header, values))

    reference_reads = float(seqkit["num_seqs"])
    reference_bases = float(seqkit["sum_len"])
    # seqkit reports GC as a percentage; the JSON reports a rate.
    reference_gc = float(seqkit.get("GC(%)", "nan")) / 100.0

    metrics = [
        Metric(
            "qc.reads_in.rel_error",
            M.relative_error(qc["reads_in"], reference_reads),
            QC_COUNT_RELATIVE_MAX, Direction.MAX,
            detail=f"pto {qc['reads_in']} vs seqkit {reference_reads:.0f}",
        ),
        Metric(
            "qc.bases_in.rel_error",
            M.relative_error(qc["bases_in"], reference_bases),
            QC_COUNT_RELATIVE_MAX, Direction.MAX,
            detail=f"pto {qc['bases_in']} vs seqkit {reference_bases:.0f}",
        ),
    ]
    if reference_gc == reference_gc:  # not nan
        metrics.append(
            Metric(
                "qc.gc_rate.abs_error",
                abs(qc["gc_rate"] - reference_gc),
                # seqkit rounds GC% to two decimals, so 5e-5 on a rate is the
                # tightest bound the reference itself can support.
                5e-5, Direction.MAX,
                detail=f"pto {qc['gc_rate']:.6f} vs seqkit {reference_gc:.6f}",
            )
        )
    return metrics


@case(
    "fastq/trimming-vs-fastp",
    module="fastq_stream",
    reference=["fastp"],
    summary="Adapter/quality trimming parity with fastp",
)
def trimming_vs_fastp(ctx: Context) -> list[Metric]:
    """Compare surviving reads, surviving bases and per-cycle quality."""
    binary = ctx.binary("fastq_stream")
    fastq = ctx.input("fastq_r1")
    work = ctx.scratch("fastq_fastp")

    # Both tools are configured to the same policy: TruSeq adapter, a 4 bp
    # sliding window cut at mean Q20, and a 15 bp minimum length. Any parameter
    # left at a tool's own default would show up in the metrics as a trimming
    # difference and be indistinguishable from a real divergence.
    pto_out = work / "pto_trimmed.fq"
    pto_json = work / "pto.json"
    run([
        binary,
        "-i", fastq, "-o", pto_out, "-j", pto_json,
        "--window", 4, "--window-mean", 20,
        "--min-len", 15,
        "--adapter", "AGATCGGAAGAGC",
        "--threads", ctx.threads,
    ])
    ctx.artifact("pto trimmed", pto_out)

    fastp_out = work / "fastp_trimmed.fq"
    fastp_json = work / "fastp.json"
    run([
        ctx.reference("fastp"),
        "-i", fastq, "-o", fastp_out,
        "--json", fastp_json, "--html", work / "fastp.html",
        "--cut_right", "--cut_right_window_size", 4, "--cut_right_mean_quality", 20,
        "--length_required", 15,
        "--adapter_sequence", "AGATCGGAAGAGC",
        "--thread", min(ctx.threads, 16),  # fastp caps at 16 and warns above it
    ])
    ctx.artifact("fastp trimmed", fastp_out)

    pto = json.loads(pto_json.read_text())
    fastp = json.loads(fastp_json.read_text())
    after = fastp["summary"]["after_filtering"]

    reads_rel = M.relative_error(pto["summary"]["reads_out"], after["total_reads"])
    bases_rel = M.relative_error(pto["summary"]["bases_out"], after["total_bases"])
    q30_abs = abs(pto["summary"]["q30_rate"] - after["q30_rate"])

    if reads_rel is not None and reads_rel > 0.02:
        _write_length_diff(ctx, work, pto_out, fastp_out)

    return [
        # 2%, not 0: the two use different window semantics at the read end
        # (fastp's --cut_right stops at the last full window, fastq_stream
        # evaluates the trailing partial one), which moves a small number of
        # reads across the 15 bp survival threshold. A larger gap is a real
        # difference in the trimming decision, not a boundary effect.
        Metric("trim.reads_out.rel_error", reads_rel, 0.02, Direction.MAX,
               detail=f"pto {pto['summary']['reads_out']} vs "
                      f"fastp {after['total_reads']}"),
        Metric("trim.bases_out.rel_error", bases_rel, 0.02, Direction.MAX,
               detail=f"pto {pto['summary']['bases_out']} vs "
                      f"fastp {after['total_bases']}"),
        Metric("trim.q30_rate.abs_error", q30_abs, QC_QUALITY_ABS_MAX, Direction.MAX),
        Metric("trim.adapters_trimmed", float(pto["filtering"]["adapter_trimmed"]),
               None, Direction.MIN, required=False,
               detail=f"fastp reports "
                      f"{fastp.get('adapter_cutting', {}).get('adapter_trimmed_reads', 'n/a')}"),
    ]


def _write_length_diff(ctx: Context, work: Path, a: Path, b: Path) -> None:
    """Emit the read-length distributions of two trimmed FASTQs, side by side.

    When the surviving-read counts diverge, the length histograms localise it
    immediately: a spike at the minimum length is a threshold disagreement, a
    uniform shift is an adapter-position disagreement, and the two call for
    completely different investigations.
    """
    from ..report import write_histogram_diff

    ctx.artifact(
        "trim length diff",
        write_histogram_diff(
            work / "length_diff.tsv", "pto", "fastp",
            _length_histogram(a), _length_histogram(b),
        ),
    )


def _length_histogram(path: Path) -> dict[int, float]:
    hist: dict[int, float] = {}
    with path.open() as handle:
        for index, line in enumerate(handle):
            if index % 4 == 1:  # sequence line
                length = len(line.rstrip("\n"))
                hist[length] = hist.get(length, 0.0) + 1.0
    return hist


@case(
    "fastq/alignment-rate-vs-bowtie2",
    module="fastq_stream",
    reference=["bowtie2", "samtools"],
    summary="Trimmed reads align at least as well as untrimmed reads",
)
def alignment_rate(ctx: Context) -> list[Metric]:
    """Align raw and trimmed reads with bowtie2 and compare mapping rates.

    The assertion is directional rather than exact: trimming is supposed to
    *raise* the overall alignment rate (adapter-contaminated reads start
    mapping) while lowering the total read count. A trimmer that mangles reads
    shows up here as a rate that drops, which no counter comparison detects --
    a corrupted base still counts as a base.
    """
    index_prefix = ctx.param("bowtie2_index")
    if not index_prefix:
        raise Unavailable(
            "dataset declares no `bowtie2_index` param; set it to a prebuilt "
            "bowtie2 index prefix (building one from a genome FASTA is a "
            "multi-hour job this harness deliberately does not do for you)"
        )
    index_prefix = Path(index_prefix).expanduser()
    if not Path(str(index_prefix) + ".1.bt2").exists():
        raise Unavailable(f"no bowtie2 index at {index_prefix}.*.bt2")

    binary = ctx.binary("fastq_stream")
    fastq = ctx.input("fastq_r1")
    work = ctx.scratch("fastq_bowtie2")

    trimmed = work / "trimmed.fq"
    run([binary, "-i", fastq, "-o", trimmed, "--threads", ctx.threads])

    raw_rate = _align_rate(ctx, index_prefix, fastq, work / "raw.log")
    trimmed_rate = _align_rate(ctx, index_prefix, trimmed, work / "trimmed.log")
    ctx.artifact("bowtie2 raw log", work / "raw.log")
    ctx.artifact("bowtie2 trimmed log", work / "trimmed.log")

    return [
        Metric(
            "align.rate_delta", trimmed_rate - raw_rate,
            # Allowed to be very slightly negative: on an already-clean library
            # there is no adapter to remove, and quality trimming shortens a
            # few reads past bowtie2's seed length. A drop beyond 0.5 points
            # means the trimmer is damaging reads.
            -0.005, Direction.MIN, unit=" (fraction)",
            detail=f"raw {raw_rate:.4f} -> trimmed {trimmed_rate:.4f}",
        ),
        Metric("align.rate_raw", raw_rate, None, Direction.MIN, required=False),
        Metric("align.rate_trimmed", trimmed_rate, None, Direction.MIN, required=False),
    ]


_RATE = re.compile(r"([\d.]+)% overall alignment rate")


def _align_rate(ctx: Context, index_prefix: Path, reads: Path, log: Path) -> float:
    """Run bowtie2 and parse the overall alignment rate off its stderr summary.

    Output goes to /dev/null: this case cares only about the summary line, and
    a genome-scale SAM stream is tens of gigabytes that nothing here reads.
    """
    result = run(
        [
            ctx.reference("bowtie2"),
            "-x", index_prefix, "-U", reads,
            "-p", ctx.threads, "-S", "/dev/null",
        ],
        timeout=14400,
    )
    log.write_text(result.stderr)
    match = _RATE.search(result.stderr)
    if not match:
        raise Unavailable("could not parse bowtie2's alignment rate summary")
    return float(match.group(1)) / 100.0
