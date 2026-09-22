"""The Python surface: one function per analysis a lab actually asks for.

TWO DESIGN DECISIONS WORTH THE PARAGRAPH EACH
---------------------------------------------
**Unspecified parameters are omitted from argv, not defaulted here.** Every
keyword below defaults to `None`, and a `None` contributes no flag. The binary
then applies its own documented default.

The alternative -- mirroring the defaults in these signatures here as well --
would put a second copy of every bound and default in the tree, after the
binary's own CLI parser already states them once, in the one place a change
cannot be missed.

A hosted, multi-tenant front end over these same tools is right to restate
its own bounds: it validates before provisioning compute and moving a large
input, so a bad value must fail fast and cheaply rather than after the work
has started. That economic argument does not exist here. The binary is on
this host, it validates in microseconds, and its own error message names the
flag and the range. Being a thin argv wrapper is the honest shape for this
layer.

**Reports come back parsed, files only when asked for.** `run_qc(reads)`
returns the QC report as a dict with nothing left on disk, because that is the
call an evaluator makes first. Pass `out_dir=` and the artifacts are kept and
their paths returned.
"""

from __future__ import annotations

import contextlib
import tempfile
from pathlib import Path
from typing import Iterable, Iterator, Sequence

from ._cpu import available_cpus
from ._run import Result, run_tool

__all__ = [
    "call_peaks",
    "frip",
    "fragment_sizes",
    "mark_duplicates",
    "profile",
    "run_qc",
    "trim",
]


def _threads(threads: int | None) -> str:
    """Always explicit. See `_cpu.available_cpus` for why."""
    return str(max(1, threads if threads else available_cpus()))


def _existing(path: str | Path, what: str) -> Path:
    """Fail before spawning a process, with the path the caller passed.

    Not validation of the input's *contents* -- the binary owns that. This is
    the one check worth doing here, because a typo'd path otherwise surfaces
    as a tool's own "cannot open" on stderr wrapped in a ToolError, several
    frames from the call that made it.
    """
    resolved = Path(path).expanduser()
    if not resolved.exists():
        raise FileNotFoundError(f"{what} not found: {resolved}")
    return resolved


def _flag(argv: list[str], name: str, value: object) -> None:
    """Append `--name value` unless the caller said nothing."""
    if value is None:
        return
    if isinstance(value, bool):
        raise TypeError(f"{name} is a valued flag, not a switch")
    if isinstance(value, float):
        argv += [name, f"{value:g}"]
    else:
        argv += [name, str(value)]


def _switch(argv: list[str], name: str, value: bool | None) -> None:
    """Append a bare `--name` only when the caller explicitly asked for it."""
    if value:
        argv.append(name)


@contextlib.contextmanager
def _report_dir(out_dir: str | Path | None) -> Iterator[tuple[Path, bool]]:
    """Somewhere for a tool to write its reports; (dir, keep)."""
    if out_dir is not None:
        target = Path(out_dir).expanduser()
        target.mkdir(parents=True, exist_ok=True)
        yield target, True
        return
    with tempfile.TemporaryDirectory(prefix="pto-report-") as scratch:
        yield Path(scratch), False


# ---------------------------------------------------------------------------
# fastq_stream
# ---------------------------------------------------------------------------

def _fastq_argv(
    reads: Path, *, threads, trim_front, trim_tail, window, window_mean, min_len,
    max_n_rate, min_mean_q, trim_adapters, adapters, adapter_seed,
    adapter_overlap, adapter_mismatch, reference,
) -> list[str]:
    argv = ["-i", str(reads), "-t", _threads(threads)]
    _flag(argv, "--trim-front", trim_front)
    _flag(argv, "--trim-tail", trim_tail)
    _flag(argv, "--window", window)
    _flag(argv, "--window-mean", window_mean)
    _flag(argv, "--min-len", min_len)
    _flag(argv, "--max-n-rate", max_n_rate)
    _flag(argv, "--min-mean-q", min_mean_q)

    # `trim_adapters=False` disables it; leaving it None keeps the tool's
    # built-in TruSeq/Nextera pair, which is what a first run should use.
    if trim_adapters is False:
        argv.append("--no-adapter")
    for seq in adapters or ():
        argv += ["-a", str(seq)]
    _flag(argv, "--adapter-seed", adapter_seed)
    _flag(argv, "--adapter-overlap", adapter_overlap)
    _flag(argv, "--adapter-mismatch", adapter_mismatch)

    # CRAM without a local reference makes htslib resolve it over the network.
    # On an air-gapped cluster node that stalls rather than fails.
    if reference is not None:
        argv += ["--reference", str(_existing(reference, "reference FASTA"))]
    return argv


def run_qc(
    reads: str | Path,
    *,
    out_dir: str | Path | None = None,
    threads: int | None = None,
    trim_front: int | None = None,
    trim_tail: int | None = None,
    window: int | None = None,
    window_mean: float | None = None,
    min_len: int | None = None,
    max_n_rate: float | None = None,
    min_mean_q: float | None = None,
    trim_adapters: bool | None = None,
    adapters: Sequence[str] | None = None,
    adapter_seed: int | None = None,
    adapter_overlap: int | None = None,
    adapter_mismatch: float | None = None,
    reference: str | Path | None = None,
    timeout: float | None = None,
) -> Result:
    """Streaming FASTQ QC. Writes no reads.

    >>> import pto
    >>> qc = pto.run_qc("sample_R1.fastq.gz")
    >>> qc["summary"]["reads_in"], qc["summary"]["q30_rate"]

    Accepts FASTQ (plain/gzip/BGZF) and, where `fastq_stream` was built with
    htslib, BAM/CRAM/SAM detected by content rather than extension.
    """
    reads_path = _existing(reads, "reads")
    with _report_dir(out_dir) as (target, keep):
        json_out = target / "qc_report.json"
        text_out = target / "qc_report.txt"
        argv = _fastq_argv(
            reads_path, threads=threads, trim_front=trim_front, trim_tail=trim_tail,
            window=window, window_mean=window_mean, min_len=min_len,
            max_n_rate=max_n_rate, min_mean_q=min_mean_q,
            trim_adapters=trim_adapters, adapters=adapters,
            adapter_seed=adapter_seed, adapter_overlap=adapter_overlap,
            adapter_mismatch=adapter_mismatch, reference=reference,
        )
        argv += ["--qc-only", "-j", str(json_out), "--report", str(text_out)]
        return run_tool(
            "fastq_stream", argv, json_path=json_out, timeout=timeout,
            artifacts={"qc_json": json_out, "qc_report": text_out} if keep else {},
        )


def trim(
    reads: str | Path,
    out: str | Path,
    *,
    out_dir: str | Path | None = None,
    threads: int | None = None,
    trim_front: int | None = None,
    trim_tail: int | None = None,
    window: int | None = None,
    window_mean: float | None = None,
    min_len: int | None = None,
    max_n_rate: float | None = None,
    min_mean_q: float | None = None,
    trim_adapters: bool | None = None,
    adapters: Sequence[str] | None = None,
    adapter_seed: int | None = None,
    adapter_overlap: int | None = None,
    adapter_mismatch: float | None = None,
    reference: str | Path | None = None,
    timeout: float | None = None,
) -> Result:
    """Trim and filter reads, writing cleaned FASTQ to `out`.

    `out="-"` streams to stdout, which is how the tool is meant to be used in
    an alignment pipeline -- though from Python a named FIFO is usually easier
    to reason about than plumbing a pipe through `subprocess`.

    This is the one pto-core job whose output is roughly the size of its input;
    everything else here is a reducer.
    """
    reads_path = _existing(reads, "reads")
    out_path = Path(out).expanduser() if str(out) != "-" else Path("-")
    if str(out_path) != "-":
        out_path.parent.mkdir(parents=True, exist_ok=True)

    with _report_dir(out_dir) as (target, keep):
        json_out = target / "qc_report.json"
        text_out = target / "qc_report.txt"
        argv = _fastq_argv(
            reads_path, threads=threads, trim_front=trim_front, trim_tail=trim_tail,
            window=window, window_mean=window_mean, min_len=min_len,
            max_n_rate=max_n_rate, min_mean_q=min_mean_q,
            trim_adapters=trim_adapters, adapters=adapters,
            adapter_seed=adapter_seed, adapter_overlap=adapter_overlap,
            adapter_mismatch=adapter_mismatch, reference=reference,
        )
        argv += ["-o", str(out_path), "-j", str(json_out), "--report", str(text_out)]
        artifacts = {"trimmed": out_path}
        if keep:
            artifacts |= {"qc_json": json_out, "qc_report": text_out}
        return run_tool(
            "fastq_stream", argv, json_path=json_out, timeout=timeout,
            artifacts=artifacts,
        )


# ---------------------------------------------------------------------------
# cuttag_profiler
# ---------------------------------------------------------------------------

def profile(
    alignment: str | Path,
    regions: str | Path,
    out_dir: str | Path,
    *,
    threads: int | None = None,
    reference_point: str | None = None,
    upstream: int | None = None,
    downstream: int | None = None,
    bin_size: int | None = None,
    normalization: str | None = None,
    count_mode: str | None = None,
    min_mapq: int | None = None,
    max_fragment: int | None = None,
    extend_reads: int | None = None,
    extend_fragment: bool | None = None,
    proper_pairs_only: bool | None = None,
    respect_strand: bool | None = None,
    timeout: float | None = None,
) -> Result:
    """CUT&Tag / CUT&RUN reference-point signal matrix from an indexed BAM.

    Writes `matrix.tsv` (regions x bins) and `profile.tsv` (the meta-profile
    curve) into `out_dir`. The BAM must be indexed -- `cuttag_profiler` opens
    it through the htslib index, and a missing `.bai`/`.csi` is an error rather
    than a full scan.

    `serve`, the tool's interactive dashboard, is not something this function
    builds: a library call that started a listening socket on a shared
    cluster node would be the wrong default even where the binary still has
    one compiled in. It is a single-tenant loopback tool whose session token
    is bound to one operator's browser. Render `profile.tsv` instead.

    (`pto-core exec cuttag_profiler serve ...` reaches it anyway -- `exec` has
    no subcommand allowlist -- so "unreachable" describes this function, not
    the installed binary.)
    """
    bam = _existing(alignment, "alignment")
    regions_path = _existing(regions, "regions BED/GTF")
    target = Path(out_dir).expanduser()
    target.mkdir(parents=True, exist_ok=True)
    matrix, curve = target / "matrix.tsv", target / "profile.tsv"

    argv = [
        "--bam", str(bam),
        "--regions", str(regions_path),
        "--out-matrix", str(matrix),
        "--out-profile", str(curve),
        "--threads", _threads(threads),
        "--quiet",
    ]
    _flag(argv, "--reference-point", reference_point)
    _flag(argv, "--upstream", upstream)
    _flag(argv, "--downstream", downstream)
    _flag(argv, "--bin-size", bin_size)
    _flag(argv, "--normalization", normalization)
    _flag(argv, "--count-mode", count_mode)
    _flag(argv, "--min-mapq", min_mapq)
    _flag(argv, "--max-fragment", max_fragment)
    _flag(argv, "--extend-reads", extend_reads)
    if extend_fragment is False:
        argv.append("--no-extend-fragment")
    _switch(argv, "--proper-pairs-only", proper_pairs_only)
    if respect_strand is False:
        argv.append("--ignore-strand")

    return run_tool(
        "cuttag_profiler", argv, subcommand="profile", timeout=timeout,
        artifacts={"matrix": matrix, "profile": curve},
    )


# ---------------------------------------------------------------------------
# genomic_toolkit
# ---------------------------------------------------------------------------

def _gtk_argv(
    alignment: Path, *, regions, min_length, max_length, min_mapq, single_end,
    unsorted,
) -> list[str]:
    argv = [str(alignment), "--json"]
    for region in regions or ():
        # Repeatable, and read as one stream with summed statistics. A fragment
        # belongs to the region containing its START, so a tiling partitions
        # the file exactly rather than double-counting boundary-straddling
        # fragments the way an overlap query would.
        argv += ["--region", str(region)]
    _flag(argv, "--min-length", min_length)
    _flag(argv, "--max-length", max_length)
    _flag(argv, "--min-mapq", min_mapq)
    _switch(argv, "--single-end", single_end)
    _switch(argv, "--unsorted", unsorted)
    return argv


def fragment_sizes(
    alignment: str | Path,
    *,
    histogram: str | Path | None = None,
    regions: Iterable[str] | None = None,
    min_length: int | None = None,
    max_length: int | None = None,
    min_mapq: int | None = None,
    single_end: bool | None = None,
    unsorted: bool | None = None,
    timeout: float | None = None,
) -> Result:
    """Fragment-size distribution in one streaming pass -- no sort, no temp files.

    >>> pto.fragment_sizes("atac.bam")["median"]

    Accepts `.bam`/`.cram`, `.bedpe`, or a 3+ column fragment `.bed`/`.tsv`;
    `genomic_toolkit` dispatches on the extension. Pass `histogram=` to keep
    the per-length TSV.
    """
    bam = _existing(alignment, "alignment")
    argv = _gtk_argv(
        bam, regions=regions, min_length=min_length, max_length=max_length,
        min_mapq=min_mapq, single_end=single_end, unsorted=unsorted,
    )
    artifacts = {}
    if histogram is not None:
        hist = Path(histogram).expanduser()
        hist.parent.mkdir(parents=True, exist_ok=True)
        argv += ["--histogram", str(hist)]
        artifacts["histogram"] = hist
    return run_tool(
        "genomic_toolkit", argv, subcommand="sizes", stdout_json=True,
        timeout=timeout, artifacts=artifacts,
    )


def mark_duplicates(
    alignment: str | Path,
    *,
    regions: Iterable[str] | None = None,
    min_length: int | None = None,
    max_length: int | None = None,
    min_mapq: int | None = None,
    single_end: bool | None = None,
    unsorted: bool | None = None,
    timeout: float | None = None,
) -> Result:
    """Duplicate rate and estimated library size, streaming.

    The default mode keeps duplicate state for one start position at a time,
    which is exact on coordinate-sorted input and costs O(pile-up depth) in
    memory. On unsorted input it reports `out_of_order` and exits 2 -- raising
    `InvalidInput` with the report attached, because the rate is then an
    undercount rather than nonsense. Pass `unsorted=True` for a full
    fingerprint table instead.
    """
    bam = _existing(alignment, "alignment")
    argv = _gtk_argv(
        bam, regions=regions, min_length=min_length, max_length=max_length,
        min_mapq=min_mapq, single_end=single_end, unsorted=unsorted,
    )
    return run_tool(
        "genomic_toolkit", argv, subcommand="markdup", stdout_json=True,
        timeout=timeout,
    )


def frip(
    alignment: str | Path,
    peaks: str | Path,
    *,
    genome_size: int | None = None,
    keep_dups: bool | None = None,
    regions: Iterable[str] | None = None,
    min_length: int | None = None,
    max_length: int | None = None,
    min_mapq: int | None = None,
    single_end: bool | None = None,
    unsorted: bool | None = None,
    timeout: float | None = None,
) -> Result:
    """Fraction of reads in peaks, in the same pass as size and duplicate QC.

    >>> pto.frip("atac.bam", "peaks.narrowPeak", genome_size=2_913_022_398)["frip"]

    Peaks are merged into a disjoint cover at load time: that is the definition
    of FRiP rather than an optimisation, so overlapping input peaks do not
    inflate the denominator.

    Exit 2 -- `InvalidInput` -- means no fragment overlapped any peak *and* the
    stream's contigs are absent from the peak file. That is almost always
    `chr1` vs `1` naming, and it is raised rather than returned because the
    alternative is a FRiP of exactly zero that looks like a failed experiment.

    `keep_dups=True` scores with duplicates retained. Scoring both ways is the
    dual-path mode low-input CUT&RUN work needs; call this twice and compare,
    which is exactly what the cloud API's `duplicate_handling="both"` does.
    """
    bam = _existing(alignment, "alignment")
    peaks_path = _existing(peaks, "peaks BED")
    argv = _gtk_argv(
        bam, regions=regions, min_length=min_length, max_length=max_length,
        min_mapq=min_mapq, single_end=single_end, unsorted=unsorted,
    )
    # Inserted after the input but before the shared tail, matching the CLI's
    # documented order; the parser is order-independent, the diff is not.
    argv = argv[:2] + ["--peaks", str(peaks_path)] + argv[2:]
    _flag(argv, "--genome-size", genome_size)
    _switch(argv, "--keep-dups", keep_dups)
    return run_tool(
        "genomic_toolkit", argv, subcommand="frip", stdout_json=True,
        timeout=timeout,
    )


# ---------------------------------------------------------------------------
# pto-peaks
# ---------------------------------------------------------------------------

def call_peaks(
    fragments: str | Path,
    output_path: str | Path,
    *,
    chrom_sizes: str | Path | None = None,
    cutoff: float | None = None,
    qvalue: float | None = None,
    min_length: int | None = None,
    max_gap: int | None = None,
    extend_peaks: int | None = None,
    lambda_bg: float | None = None,
    min_mapq: int | None = None,
    timeout: float | None = None,
) -> Result:
    """Call peaks in one streaming pass; writes narrowPeak to `output_path`.

    >>> pto.call_peaks("frags.bed", "peaks.narrowPeak", chrom_sizes="hg38.chrom.sizes")
    >>> pto.call_peaks("atac.bam", "peaks.narrowPeak")     # BAM carries its own

    Input is a coordinate-sorted fragment BED (chrom, start, end) or a
    coordinate-sorted BAM. `pto-peaks` tells them apart from the first byte, so
    this wrapper does not: a second copy of that detection here would be a
    second thing to keep in step, and it would have to guess on a path it never
    opens. That is also why `chrom_sizes` is optional rather than conditionally
    required -- a BAM carries its contig lengths in its header, a BED does not,
    and the binary raises the one clear error either way.

    `output_path` is required, unlike the tool's own `--out`, which defaults to
    stdout. Returning peaks through stdout would mean this package buffering an
    unbounded TSV in memory to hand back a string, and `run_tool` captures
    stdout to a scratch file it then deletes -- so an optional `output_path`
    would have a default that silently discards the answer.

    `extend_peaks` pads every reported interval symmetrically by that many bp.
    The default (`None`, equivalent to `0`) is `pto-peaks`' sharp,
    high-resolution call -- the enrichment core and nothing more. A positive
    value is a compatibility mode for pipelines that expect the wider spans
    MACS2/MACS3 produce (`pto-core/modules/peaks/CONCORDANCE.md` records the
    ~3x median-width gap it closes). It moves only the bounds: the summit, its
    -log10 p and -log10 q are unchanged, the padding is applied after the
    `min_length` gate, and the 3' edge is clamped to the contig length.

    **Paired-end only.** A BAM yields one fragment per proper pair, spanning
    TLEN from the leftmost mate; every record with `TLEN <= 0` is skipped, and
    a single-end BAM has TLEN 0 on all of them. Such a run exits 0 having
    called nothing, and reports the count under `skipped_unpaired` on
    `.stderr`. Convert single-end data to a fragment BED with the fragment
    length you intend, and pass that.

    `.metrics` is empty: `pto-peaks` prints narrowPeak, not a JSON report, and
    this package does not invent one. What the run cost is on `.wall_seconds`
    and `.peak_rss_mb`; what it skipped and why is on `.stderr`; the peaks are
    the file named by `.artifacts["peaks"]`.

    Exit 2 -- `InvalidInput` -- means the input's shape made the answer wrong:
    records out of coordinate order, a contig regrouped after being closed, or
    a contig absent from `chrom_sizes` (`chr1` vs `1`, almost always). It is
    deterministic, so it is not worth retrying.
    """
    source = _existing(fragments, "fragments BED/BAM")
    target = Path(output_path).expanduser()
    target.parent.mkdir(parents=True, exist_ok=True)

    argv: list[str] = []
    if chrom_sizes is not None:
        argv += ["--chrom-sizes", str(_existing(chrom_sizes, "chrom sizes"))]
    argv += ["--out", str(target)]
    _flag(argv, "--cutoff", cutoff)
    _flag(argv, "--qvalue", qvalue)
    _flag(argv, "--min-length", min_length)
    _flag(argv, "--max-gap", max_gap)
    _flag(argv, "--extend-peaks", extend_peaks)
    _flag(argv, "--lambda-bg", lambda_bg)
    _flag(argv, "--min-mapq", min_mapq)
    # The input is positional and last: `--chrom-sizes` before it matches the
    # tool's own usage line, and the parser is order-independent anyway.
    argv.append(str(source))

    return run_tool(
        "pto-peaks", argv, artifacts={"peaks": target}, timeout=timeout,
    )
