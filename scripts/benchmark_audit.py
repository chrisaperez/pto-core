#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Resource and performance audit for pto-core.

    scripts/benchmark_audit.py --list
    scripts/benchmark_audit.py --inputs small --threads 1 2 4 8
    scripts/benchmark_audit.py --inputs small medium --baselines --repeats 3
    scripts/benchmark_audit.py --report build-bench/audit

Measures wall-clock time, peak RSS and disk I/O for each tool across a
configuration matrix of input sizes and thread counts, and -- when
`--baselines` is given -- for the standard pipeline each tool replaces, on the
same input and the same core count.

**The point is the unvarnished number, not the flattering one.** Three rules
follow from that and are enforced in the code rather than left to discipline:

* A speedup is only reported when both sides ran on the same input in the same
  session. Comparing today's run against a number from a README is how
  benchmarks start lying.
* Peak RSS is read from the OS's own accounting for the child process tree
  (`wait4`/`getrusage`), not sampled from Python. A sampler misses the peak
  almost by definition -- it is a spike, and the sampler is asleep during it.
* Thread scaling is reported as measured parallel efficiency against this
  machine's own 1-thread run, so a tool that scales badly says so. Efficiency
  above 1.0 is flagged rather than celebrated: it usually means the 1-thread
  baseline was contended, not that the tool beat Amdahl.

Standard library only, same reasoning as the rest of the repo's tooling: this
has to run on the cluster that holds the data.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
import platform
import shutil
import statistics
import subprocess
import sys
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Callable, Sequence

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_BUILD_DIR = Path(os.environ.get("PTO_BUILD_DIR", REPO_ROOT / "build"))
DEFAULT_DATA_DIR = Path(
    os.environ.get("PTO_VALIDATION_DATA", Path.home() / ".cache" / "pto-validation")
)

#: ru_maxrss is bytes on Darwin and kibibytes on Linux. Getting this wrong
#: silently reports a 1000x memory difference between platforms, which is
#: exactly the kind of error a benchmark table makes authoritative.
_RSS_SCALE = 1 if sys.platform == "darwin" else 1024


# ---------------------------------------------------------------------------
# Measurement
# ---------------------------------------------------------------------------


@dataclass
class Measurement:
    """One execution of one command."""

    wall_seconds: float
    user_seconds: float
    system_seconds: float
    peak_rss_bytes: int
    #: Blocks the OS actually moved to/from disk. On Linux these are 512-byte
    #: block counts; on Darwin they are usually 0 because the kernel does not
    #: charge page-cache-satisfied reads to the process. Reported as-is with
    #: that caveat attached rather than converted into a fabricated MB/s.
    block_reads: int
    block_writes: int
    voluntary_switches: int
    involuntary_switches: int
    returncode: int
    output_bytes: int = 0

    @property
    def cpu_seconds(self) -> float:
        return self.user_seconds + self.system_seconds

    @property
    def cpu_utilisation(self) -> float:
        """Mean cores busy over the run. ~1.0 means it never went parallel."""
        return self.cpu_seconds / self.wall_seconds if self.wall_seconds > 0 else 0.0


def measure(
    command: Sequence[str],
    *,
    stdout_path: Path | None = None,
    env: dict[str, str] | None = None,
    cwd: Path | None = None,
) -> Measurement:
    """Run a command and return its resource usage, from the OS's accounting.

    Uses `os.wait4`, which hands back the child's `rusage` directly. That is
    the whole reason this is not `subprocess.run` plus a sampling thread: peak
    RSS is a spike, and a sampler is statistically guaranteed to miss it. The
    kernel is not.

    Note this attributes usage to the immediate child only. Every command
    driven here is a single process (or one that reaps its own workers into
    its own rusage), so that is what we want -- but a shell pipeline measured
    this way would report only the shell, which is why `pipeline()` below
    measures each stage separately instead of wrapping one in `sh -c`.
    """
    sink = None
    try:
        if stdout_path is not None:
            stdout_path.parent.mkdir(parents=True, exist_ok=True)
            sink = stdout_path.open("wb")
            stdout_fd = sink.fileno()
        else:
            stdout_fd = subprocess.DEVNULL

        started = time.perf_counter()
        process = subprocess.Popen(
            [str(part) for part in command],
            stdout=stdout_fd,
            stderr=subprocess.DEVNULL,
            env={**os.environ, **(env or {})},
            cwd=cwd,
        )
        _, status, usage = os.wait4(process.pid, 0)
        elapsed = time.perf_counter() - started
        process.returncode = os.waitstatus_to_exitcode(status)
    finally:
        if sink is not None:
            sink.close()

    return Measurement(
        wall_seconds=elapsed,
        user_seconds=usage.ru_utime,
        system_seconds=usage.ru_stime,
        peak_rss_bytes=usage.ru_maxrss * _RSS_SCALE,
        block_reads=usage.ru_inblock,
        block_writes=usage.ru_oublock,
        voluntary_switches=usage.ru_nvcsw,
        involuntary_switches=usage.ru_nivcsw,
        returncode=process.returncode,
        output_bytes=(
            stdout_path.stat().st_size
            if stdout_path is not None and stdout_path.exists()
            else 0
        ),
    )


def combine(measurements: Sequence[Measurement]) -> Measurement:
    """Fold a multi-stage pipeline into one Measurement.

    Wall time sums, because the stages here are run sequentially (materialising
    the intermediate file), which is exactly how the legacy pipeline being
    compared against is normally run. Peak RSS takes the **max**, not the sum:
    the stages do not overlap, so the machine only ever had to hold the largest
    one. Summing would overstate the baseline's cost and flatter pto, which
    is the specific direction of error this whole script exists to avoid.
    """
    return Measurement(
        wall_seconds=sum(m.wall_seconds for m in measurements),
        user_seconds=sum(m.user_seconds for m in measurements),
        system_seconds=sum(m.system_seconds for m in measurements),
        peak_rss_bytes=max((m.peak_rss_bytes for m in measurements), default=0),
        block_reads=sum(m.block_reads for m in measurements),
        block_writes=sum(m.block_writes for m in measurements),
        voluntary_switches=sum(m.voluntary_switches for m in measurements),
        involuntary_switches=sum(m.involuntary_switches for m in measurements),
        returncode=next((m.returncode for m in measurements if m.returncode != 0), 0),
        output_bytes=sum(m.output_bytes for m in measurements),
    )


# ---------------------------------------------------------------------------
# Benchmark definitions
# ---------------------------------------------------------------------------


@dataclass
class Result:
    benchmark: str
    variant: str          # "pto" or "baseline"
    input_key: str
    input_bytes: int
    threads: int
    repeats: list[Measurement] = field(default_factory=list)
    skipped: str = ""

    @property
    def best(self) -> Measurement | None:
        """The fastest repeat, not the mean.

        Deliberate: on a shared machine the slow repeats are measuring other
        people's jobs, not this one. The spread across repeats is reported
        separately so a noisy machine is visible rather than averaged into the
        headline number.
        """
        completed = [m for m in self.repeats if m.returncode == 0]
        if not completed:
            return None
        return min(completed, key=lambda m: m.wall_seconds)

    @property
    def spread(self) -> float | None:
        """Relative standard deviation of wall time across repeats."""
        times = [m.wall_seconds for m in self.repeats if m.returncode == 0]
        if len(times) < 2:
            return None
        mean = statistics.fmean(times)
        return statistics.stdev(times) / mean if mean > 0 else None


@dataclass
class Benchmark:
    name: str
    module: str
    #: (context, threads) -> list of commands to measure, or raises Skip.
    pto: Callable
    baseline: Callable | None = None
    baseline_tools: tuple[str, ...] = ()
    #: Thread-scaling is meaningless for a single-threaded tool; those
    #: benchmarks declare it so the matrix does not run five identical
    #: configurations and present them as a scaling curve.
    scales: bool = True
    summary: str = ""


class Skip(Exception):
    """This benchmark cannot run in this configuration."""


@dataclass
class Env:
    build_dir: Path
    data_dir: Path
    work_dir: Path
    inputs: dict[str, Path]
    #: The input tier this Env was built for ("small"/"medium"/"deep"). Only
    #: consulted by benchmarks with no staged file to size themselves off of
    #: (scrna_matrix's synthetic fixture) -- every file-based benchmark reads
    #: its size from `inputs` instead, which is already tier-resolved.
    tier: str = ""

    def binary(self, module: str, name: str | None = None) -> Path:
        """Path to a built CLI binary.

        `name` overrides the executable name when it differs from the module
        directory name (`pto-peaks` under `modules/peaks/`, matching the CMake
        target in modules/peaks/CMakeLists.txt) -- every other caller omits it
        and gets today's behaviour unchanged.
        """
        name = name or module
        candidates = [
            self.build_dir / f"modules/{module}/{name}",
            self.build_dir / module / name,
        ]
        for candidate in candidates:
            if candidate.is_file() and os.access(candidate, os.X_OK):
                return candidate
        raise Skip(f"{name} is not built under {self.build_dir}")

    def python_module_dir(self, module: str, name: str) -> Path:
        """Directory holding a built pybind11 extension, for sys.path.

        Mirrors binary()'s two-candidate search (standalone vs. aggregate
        build): scrna_matrix's pybind11_add_module target lands under a
        directory named for the module, the same as the compiled CLI binaries
        do, just with a platform-tagged `.so`/`.pyd` instead of a fixed name.
        """
        candidates = [
            self.build_dir / f"modules/{module}",
            self.build_dir / module,
        ]
        for directory in candidates:
            if any(directory.glob(f"{name}*.so")) or any(directory.glob(f"{name}*.pyd")):
                return directory
        raise Skip(
            f"{name} is not built under {self.build_dir} "
            f"(configure {module} with -DSCRNA_BUILD_PYTHON=ON and pybind11 available)"
        )

    def tool(self, name: str) -> str:
        path = shutil.which(name)
        if not path:
            raise Skip(f"{name} is not on PATH")
        return path

    def input(self, key: str) -> Path:
        path = self.inputs.get(key)
        if path is None or not path.exists():
            raise Skip(
                f"input {key!r} is not staged. Stage it with:\n"
                f"    tests/validation/run_validation.py fetch --allow-download"
            )
        return path

    def out(self, *parts: str) -> Path:
        path = self.work_dir.joinpath(*parts)
        path.parent.mkdir(parents=True, exist_ok=True)
        return path


BENCHMARKS: dict[str, Benchmark] = {}


def benchmark(name: str, *, module: str, baseline_tools: Sequence[str] = (),
              scales: bool = True, summary: str = ""):
    def decorate(pair):
        BENCHMARKS[name] = Benchmark(
            name=name, module=module,
            pto=pair[0], baseline=pair[1],
            baseline_tools=tuple(baseline_tools),
            scales=scales, summary=summary,
        )
        return pair

    return decorate


# -- fastq_stream -----------------------------------------------------------


def _fq_pto(env: Env, threads: int, key: str):
    binary = env.binary("fastq_stream")
    fastq = env.input(key)
    return [(
        [binary, "-i", fastq, "-o", env.out("bench", f"fq_{threads}.fq"),
         "--threads", threads],
        None,
    )]


def _fq_baseline(env: Env, threads: int, key: str):
    # fastp is the closest like-for-like: same job, same single invocation,
    # same output. It caps its own thread count at 16 and warns beyond that,
    # so asking for more would silently measure a 16-thread run and label it
    # 32 -- which would make pto look better than it is at high core
    # counts, for no reason other than a mislabelled axis.
    fastp = env.tool("fastp")
    fastq = env.input(key)
    return [(
        [fastp, "-i", fastq, "-o", env.out("bench", "fastp.fq"),
         "--json", env.out("bench", "fastp.json"),
         "--html", env.out("bench", "fastp.html"),
         "--thread", min(threads, 16)],
        None,
    )]


benchmark(
    "fastq_stream/trim",
    module="fastq_stream",
    baseline_tools=["fastp"],
    summary="FASTQ QC + adapter/quality trimming, vs fastp",
)((_fq_pto, _fq_baseline))


# -- genomic_toolkit --------------------------------------------------------


def _gtk_frip_pto(env: Env, threads: int, key: str):
    binary = env.binary("genomic_toolkit")
    return [(
        [binary, "frip", env.input(key), "--peaks", env.input("peaks"), "--json"],
        None,
    )]


def _gtk_frip_baseline(env: Env, threads: int, key: str):
    """The samtools + bedtools pipeline genomic_toolkit replaces.

    Measured as three separate processes rather than one `sh -c` pipeline for
    two reasons: `wait4` on a shell reports the shell's rusage and not its
    children's, and materialising the intermediates is what makes the
    baseline's real disk cost visible. That cost -- a name-sorted BAM and a
    fragment BED on disk -- is the specific thing `genomic_toolkit` exists to
    avoid, so hiding it inside a pipe would erase the comparison's whole point.
    """
    samtools = env.tool("samtools")
    bedtools = env.tool("bedtools")
    bam = env.input(key)
    namesorted = env.out("bench", "namesorted.bam")
    fragments = env.out("bench", "fragments.bed")

    return [
        ([samtools, "sort", "-n", "-@", threads, "-o", namesorted, bam], None),
        ([bedtools, "bamtobed", "-bedpe", "-i", namesorted], fragments),
        ([bedtools, "intersect", "-u", "-a", fragments, "-b", env.input("peaks")], None),
    ]


benchmark(
    "genomic_toolkit/frip",
    module="genomic_toolkit",
    baseline_tools=["samtools", "bedtools"],
    summary="FRiP in one streaming pass, vs sort + bamtobed + intersect",
)((_gtk_frip_pto, _gtk_frip_baseline))


def _gtk_sizes_pto(env: Env, threads: int, key: str):
    binary = env.binary("genomic_toolkit")
    return [(
        [binary, "sizes", env.input(key),
         "--histogram", env.out("bench", "sizes.tsv"), "--json"],
        None,
    )]


def _gtk_sizes_baseline(env: Env, threads: int, key: str):
    samtools = env.tool("samtools")
    return [(
        [samtools, "stats", "-@", threads, env.input(key)],
        env.out("bench", "samtools_stats.txt"),
    )]


benchmark(
    "genomic_toolkit/sizes",
    module="genomic_toolkit",
    baseline_tools=["samtools"],
    scales=False,  # the size pass is single-threaded by design; it is I/O bound
    summary="Fragment size distribution, vs samtools stats",
)((_gtk_sizes_pto, _gtk_sizes_baseline))


# -- cuttag_profiler --------------------------------------------------------


def _ctp_pto(env: Env, threads: int, key: str):
    binary = env.binary("cuttag_profiler")
    return [(
        [binary, "profile",
         "--bam", env.input(key),
         "--regions", env.input("regions"),
         "--out-matrix", env.out("bench", "matrix.tsv"),
         "--upstream", 2000, "--downstream", 2000, "--bin-size", 50,
         "--min-mapq", 30, "--threads", threads, "--quiet"],
        None,
    )]


def _ctp_baseline(env: Env, threads: int, key: str):
    """deeptools bamCoverage + computeMatrix, the pipeline being replaced.

    Both stages are measured. bamCoverage usually dominates and is easy to
    leave out of a benchmark by accident -- computeMatrix alone looks
    competitive, but you cannot run it without first producing the bigWig, so
    excluding that step would be measuring half a pipeline against a whole one.
    """
    bam_coverage = env.tool("bamCoverage")
    compute_matrix = env.tool("computeMatrix")
    bigwig = env.out("bench", "coverage.bw")

    return [
        ([bam_coverage, "-b", env.input(key), "-o", bigwig,
          "--binSize", 50, "--minMappingQuality", 30, "--extendReads",
          "--numberOfProcessors", threads], None),
        ([compute_matrix, "reference-point", "--referencePoint", "TSS",
          "-S", bigwig, "-R", env.input("regions"),
          "--beforeRegionStartLength", 2000, "--afterRegionStartLength", 2000,
          "--binSize", 50, "--missingDataAsZero",
          "--numberOfProcessors", threads,
          "-o", env.out("bench", "dt_matrix.gz")], None),
    ]


benchmark(
    "cuttag_profiler/profile",
    module="cuttag_profiler",
    baseline_tools=["bamCoverage", "computeMatrix"],
    summary="Reference-point signal matrix, vs bamCoverage + computeMatrix",
)((_ctp_pto, _ctp_baseline))


# A single, wide, real-coverage window rather than the TSS-anchored regions
# file above: samtools depth has no notion of "reference point" or binning, so
# the sound comparison is the same primitive both tools actually compute --
# per-base read depth over one interval -- not cuttag_profiler's TSS
# meta-profile machinery. chr1:1,000,000-1,100,000 in the staged H3K4me3 BAM
# has real coverage throughout (16,753 nonzero positions per `samtools depth`),
# so this measures real work, not an empty window.
_CTP_DEPTH_REGION = ("chr1", 1_000_000, 1_100_000)


def _ctp_depth_region_bed(env: Env) -> Path:
    chrom, start, end = _CTP_DEPTH_REGION
    bed = env.out("bench", "depth_region.bed")
    bed.write_text(f"{chrom}\t{start}\t{end}\tregion1\t0\t+\n")
    return bed


def _ctp_depth_pto(env: Env, threads: int, key: str):
    """cuttag_profiler's own BAM reader, in raw per-base depth mode.

    --reference-point center with upstream/downstream set to half the region
    width reconstructs exactly _CTP_DEPTH_REGION; --bin-size 1 and
    --normalization raw make the output one raw depth value per base, the
    same quantity `samtools depth` reports.
    """
    chrom, start, end = _CTP_DEPTH_REGION
    half = (end - start) // 2
    binary = env.binary("cuttag_profiler")
    return [(
        [binary, "profile",
         "--bam", env.input(key),
         "--regions", str(_ctp_depth_region_bed(env)),
         "--out-matrix", env.out("bench", "depth_matrix.tsv"),
         "--reference-point", "center", "--upstream", half, "--downstream", half,
         "--bin-size", 1, "--count-mode", "depth", "--normalization", "raw",
         "--min-mapq", 0, "--threads", threads, "--quiet"],
        None,
    )]


def _ctp_depth_baseline(env: Env, threads: int, key: str):
    """`samtools depth` over the identical interval -- the reference tool
    cuttag_profiler's README already cites its BAM reader against
    (docs/AUDIT_2026-09-11_cuttag_profiler.md's chromosome_length()/tid_for()
    work exists to serve exactly this kind of per-region query). Single
    process, single measurement: samtools depth has no separate index-build
    or intermediate-file stage the way the bedtools/deepTools baselines do.
    """
    chrom, start, end = _CTP_DEPTH_REGION
    samtools = env.tool("samtools")
    # samtools depth's -r region is 1-based inclusive; the BED region above is
    # 0-based half-open. +1 on start, unchanged end, is the standard conversion.
    region = f"{chrom}:{start + 1}-{end}"
    return [(
        [samtools, "depth", "-a", "-r", region, env.input(key)],
        env.out("bench", "samtools_depth.txt"),
    )]


benchmark(
    "cuttag_profiler/depth",
    module="cuttag_profiler",
    baseline_tools=["samtools"],
    scales=False,  # single fixed 100 kb window; not a throughput sweep
    summary="Per-base read depth over one region, vs samtools depth",
)((_ctp_depth_pto, _ctp_depth_baseline))


# -- peaks --------------------------------------------------------------


def _peaks_pto(env: Env, threads: int, key: str):
    binary = env.binary("peaks", name="pto-peaks")
    return [(
        [binary, "-o", env.out("bench", "peaks.narrowPeak"), env.input(key)],
        None,
    )]


def _peaks_baseline(env: Env, threads: int, key: str):
    """MACS2 callpeak -- the peak caller pto-peaks' own --help text already
    measures compatibility against (`--extend-peaks` documents matching the
    wider spans MACS2/MACS3 produce), and the reference tool
    tests/validation/datasets.json's `encode-atac-k562` entry stages ENCODE's
    own MACS2 narrowPeak calls against for the separate concordance harness
    (implementation_plan.md Phase 1.4). This benchmark now uses that same
    `encode-atac-k562` BAM (INPUT_TIERS' "bam_paired" role, medium tier only)
    rather than the "medium" tier's H3K4me3 BAM: that file is single-end
    alignments, so both pto-peaks and MACS2 -f BAMPE either read zero proper
    pairs or misinterpret unpaired reads as fragments -- a benchmark that
    measures ~nothing, not a fast one. ATAC-seq is also the domain pto-peaks
    is built for (narrow, punctate accessibility peaks), so this is a truer
    comparison than the ChIP-seq mark it replaced, not just a paired-end one.

    Treatment-only, paired-end, human genome size -- matches the ENCODE input
    this repo already stages. pto-peaks has no --threads flag (the state
    machine is single-base-sequential by design, see CLAUDE.md's "peaks"
    section), so this benchmark does not sweep thread counts either.
    """
    try:
        macs = env.tool("macs2")
    except Skip:
        macs = env.tool("macs3")  # same callpeak interface; lets this Skip if neither is present
    outdir = env.out("bench", "macs2")
    return [(
        [macs, "callpeak", "-t", env.input(key), "-n", "bench",
         "--outdir", outdir, "-f", "BAMPE", "-g", "hs", "-q", "0.05"],
        None,
    )]


benchmark(
    "peaks/call",
    module="peaks",
    baseline_tools=["macs2", "macs3"],
    scales=False,  # pto-peaks has no thread-count knob; see _peaks_baseline
    summary="Streaming peak calling from a coordinate-sorted BAM, vs MACS2 callpeak",
)((_peaks_pto, _peaks_baseline))


# -- scrna_matrix -------------------------------------------------------

_SCRNA_BENCH_SCRIPT = REPO_ROOT / "modules/scrna_matrix/tests/benchmark_1m_cells.py"
#: Cell counts per tier. Brute-force k-NN is O(n^2) (by design -- it is the
#: module's correctness baseline, see knn_graph.hpp), so these stay well under
#: the point where a routine `--inputs small` run stops being routine; "deep"
#: is the opt-in-only tier already, same as fastq_stream's 766 MB deep fixture.
_SCRNA_SYNTH_CELLS = {"small": 3000, "medium": 10000, "deep": 20000}
_SCRNA_K = 15


def _scrna_synth_path(env: Env) -> Path:
    """A synthetic sparse matrix, generated once per (work_dir, tier).

    scrna_matrix ships no reader for a real single-cell format -- that parsing
    is scanpy/anndata's job upstream of it (docs/AUDIT_2026-09-11_scrna_matrix.md
    Sec 0) -- and tests/validation/'s staged datasets are all FASTQ/BAM/BED, so
    there is no "real" input this benchmark could point at the way the other
    four do. tests/benchmark_1m_cells.py's own --make-synthetic already exists
    to fill exactly this gap (seeded with a fixed RNG, so the fixture -- not
    the timing -- is reproducible across runs); this just reuses it instead of
    inventing a second generator.
    """
    n_cells = _SCRNA_SYNTH_CELLS.get(env.tier, 3000)
    path = env.out("bench", f"scrna_synth_{n_cells}.h5ad")
    if not path.exists():
        result = subprocess.run(
            [sys.executable, str(_SCRNA_BENCH_SCRIPT),
             "--make-synthetic", str(n_cells), "--output", str(path)],
            capture_output=True, text=True,
        )
        if result.returncode != 0:
            raise Skip(
                "cannot generate the synthetic scRNA fixture (needs anndata, "
                f"numpy, scipy): {result.stderr.strip()[-200:]}"
            )
    return path


def _scrna_pto(env: Env, threads: int, key: str):
    py_dir = env.python_module_dir("scrna_matrix", "scrna_matrix_py")
    data = _scrna_synth_path(env)
    # Runs the script in-process via runpy rather than as `python3 script.py`
    # so the built extension's directory can be put on sys.path without an
    # env-var plumbing change to measure()/_execute() (which pass no per-stage
    # environment today) -- run_name='__main__' keeps the script's own
    # `if __name__ == "__main__"` entry point working unmodified.
    snippet = (
        "import sys, runpy; "
        f"sys.path.insert(0, {str(py_dir)!r}); "
        f"sys.argv = ['benchmark_1m_cells.py', '--input', {str(data)!r}, "
        f"'--backend', 'scrna', '--k', '{_SCRNA_K}']; "
        f"runpy.run_path({str(_SCRNA_BENCH_SCRIPT)!r}, run_name='__main__')"
    )
    return [([sys.executable, "-c", snippet], env.out("bench", "scrna_pto.txt"))]


def _scrna_baseline(env: Env, threads: int, key: str):
    """scanpy.pp.neighbors -- the call this module's brute-force and (when
    built) HNSW paths both exist to outrun for an anndata-based workflow.
    Optional like every other baseline tool here: scanpy is a heavy dependency
    this repo does not otherwise require, so its absence is a Skip, not a
    failure, exactly like a missing `fastp`/`samtools`/`bamCoverage`.
    """
    if importlib.util.find_spec("scanpy") is None or importlib.util.find_spec("anndata") is None:
        raise Skip("scanpy/anndata not installed (pip install scanpy)")
    data = _scrna_synth_path(env)
    return [(
        [sys.executable, str(_SCRNA_BENCH_SCRIPT),
         "--input", str(data), "--backend", "scanpy", "--k", str(_SCRNA_K)],
        env.out("bench", "scrna_baseline.txt"),
    )]


benchmark(
    "scrna_matrix/knn",
    module="scrna_matrix",
    baseline_tools=["scanpy"],
    # method="brute" (the default _scrna_pto passes implicitly) is
    # deterministic and needs no extra build flag, unlike method="hnsw" which
    # needs -DSCRNA_ENABLE_HNSW=ON. OpenMP threading is real here (unlike
    # peaks/genomic_toolkit's sizes pass) but is controlled by OMP_NUM_THREADS,
    # which nothing in this harness threads through to a subprocess yet --
    # scales=False states that honestly rather than re-running one
    # OMP_NUM_THREADS-less configuration under several --threads values and
    # presenting flat numbers as a scaling curve.
    scales=False,
    summary="Exact brute-force k-NN on a synthetic sparse matrix, vs scanpy.pp.neighbors",
)((_scrna_pto, _scrna_baseline))


# ---------------------------------------------------------------------------
# Input tiers
# ---------------------------------------------------------------------------

#: The size tiers the matrix runs over. Each maps to the dataset files the
#: validation harness already stages, so the two suites share one cache
#: instead of each downloading its own copy of the same 3 GB.
INPUT_TIERS = {
    "small": {
        "fastq": "cuttag-k562-igg/SRR8754612_1.fastq.gz",
        "description": "CUT&Tag IgG control, 1.2M read pairs (~23 MB)",
    },
    "medium": {
        "fastq": "encode-h3k4me3-gm12878/ENCFF000AUB.fastq.gz",
        "bam": "encode-h3k4me3-gm12878/ENCFF019VEK.bam",
        # peaks/call gets its own role rather than sharing "bam" with
        # cuttag_profiler/genomic_toolkit: ENCFF019VEK.bam above is single-end,
        # so a fragment-inferring, paired-end-only caller (see CLAUDE.md's
        # "peaks" section: "one fragment per proper pair, from the leftmost
        # mate's TLEN") correctly reads it as zero fragments and exits having
        # measured nothing. ENCFF121ZQX.bam is the paired-end ATAC-seq BAM
        # tests/validation/datasets.json's `encode-atac-k562` entry already
        # stages for exactly this module (implementation_plan.md Phase 1.4).
        "bam_paired": "encode-atac-k562/ENCFF121ZQX.bam",
        "peaks": "encode-h3k4me3-gm12878/ENCFF115YTM.bed.gz",
        "regions": "encode-h3k4me3-gm12878/ENCFF298YZW.gtf.gz",
        "description": "ENCODE H3K4me3 GM12878, whole-genome (~2.8 GB); peaks/call uses the paired-end ENCODE ATAC-seq K562 BAM (~230 MB) instead",
    },
    "deep": {
        "fastq": "cuttag-k562-h3k4me3-deep/SRR8754610_1.fastq.gz",
        "description": "CUT&Tag H3K4me3, 37M read pairs (~766 MB compressed)",
    },
}


def resolve_inputs(tier: str, data_dir: Path) -> dict[str, Path]:
    spec = INPUT_TIERS[tier]
    return {
        role: data_dir / relative
        for role, relative in spec.items()
        if role != "description"
    }


# ---------------------------------------------------------------------------
# Matrix execution
# ---------------------------------------------------------------------------


def run_matrix(args) -> list[Result]:
    results: list[Result] = []

    for tier in args.inputs:
        inputs = resolve_inputs(tier, args.data_dir)
        env = Env(
            build_dir=args.build_dir,
            data_dir=args.data_dir,
            work_dir=args.work_dir / tier,
            inputs=inputs,
            tier=tier,
        )

        for name, spec in BENCHMARKS.items():
            if args.benchmarks and name not in args.benchmarks:
                continue

            if spec.module == "fastq_stream":
                role = "fastq"
            elif spec.module == "scrna_matrix":
                # No staged file: the module reads no on-disk format of its
                # own (docs/AUDIT_2026-09-11_scrna_matrix.md Sec 0 -- that
                # parsing is scanpy/anndata's job upstream), so there is
                # nothing under INPUT_TIERS to point this role at.
                role = "synthetic"
            elif spec.module == "peaks":
                # Distinct from the generic "bam" role below: peaks/call needs
                # paired-end data specifically (see the "bam_paired" comment
                # in INPUT_TIERS), whereas cuttag_profiler/genomic_toolkit are
                # correctly measured against the H3K4me3 single-end BAM.
                role = "bam_paired"
            else:
                role = "bam"
            thread_counts = args.threads if spec.scales else [max(args.threads)]

            for threads in thread_counts:
                results.append(
                    _execute(env, spec, "pto", spec.pto, role, tier,
                             threads, args.repeats)
                )
                if args.baselines and spec.baseline is not None:
                    results.append(
                        _execute(env, spec, "baseline", spec.baseline, role, tier,
                                 threads, args.repeats)
                    )

    return results


def _execute(env: Env, spec: Benchmark, variant: str, builder, role: str,
             tier: str, threads: int, repeats: int) -> Result:
    label = f"{spec.name} [{variant}] {tier} t={threads}"
    result = Result(
        benchmark=spec.name, variant=variant, input_key=tier,
        input_bytes=0, threads=threads,
    )

    try:
        stages = builder(env, threads, role)
    except Skip as exc:
        result.skipped = str(exc)
        print(f"  SKIP {label}: {exc}")
        return result

    try:
        # A file-based role always resolves here too, since builder() above
        # already resolved the same role via env.input() to build its command
        # -- this is purely for the report's input_bytes column. A role with
        # no staged file (scrna_matrix's "synthetic") has nothing to size and
        # is left at 0 rather than turning a successful build into a skip.
        result.input_bytes = env.input(role).stat().st_size
    except Skip:
        pass

    print(f"  {label} ... ", end="", flush=True)
    for index in range(repeats):
        try:
            measurement = combine([
                measure(command, stdout_path=sink) for command, sink in stages
            ])
        except OSError as exc:
            result.skipped = f"execution failed: {exc}"
            print(f"failed ({exc})")
            return result

        if measurement.returncode != 0:
            result.skipped = f"exited {measurement.returncode}"
            print(f"exit {measurement.returncode}")
            return result
        result.repeats.append(measurement)
        if index == 0:
            print(f"{measurement.wall_seconds:.2f}s ", end="", flush=True)

    best = result.best
    if best:
        print(
            f"| best {best.wall_seconds:.2f}s, "
            f"{best.peak_rss_bytes / 1e6:.0f} MB peak, "
            f"{best.cpu_utilisation:.1f}x CPU"
        )
    return result


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------


def _key(result: Result) -> tuple:
    return (result.benchmark, result.input_key, result.threads)


def render_report(results: list[Result], args) -> str:
    lines: list[str] = [
        "# pto-core resource audit",
        "",
        f"- host: `{platform.platform()}`",
        f"- cores: {os.cpu_count()}",
        f"- generated: {time.strftime('%Y-%m-%d %H:%M:%S')}",
        f"- repeats per configuration: {args.repeats} (best reported)",
        "",
        "Peak RSS is the OS's own `ru_maxrss` for the process, not a sample. "
        "Where a baseline is a multi-stage pipeline, its wall time is the sum "
        "of the stages and its peak RSS is the **max**, since the stages do "
        "not overlap -- summing would overstate the baseline and flatter "
        "pto.",
        "",
    ]

    pto = {_key(r): r for r in results if r.variant == "pto"}
    baseline = {_key(r): r for r in results if r.variant == "baseline"}

    # -- head-to-head ------------------------------------------------------
    pairs = sorted(set(pto) & set(baseline))
    if pairs:
        lines += [
            "## Versus the standard pipeline",
            "",
            "| benchmark | input | threads | pto | baseline | speedup | "
            "pto peak RSS | baseline peak RSS |",
            "|---|---|---:|---:|---:|---:|---:|---:|",
        ]
        for key in pairs:
            a, b = pto[key].best, baseline[key].best
            if not a or not b:
                continue
            speedup = b.wall_seconds / a.wall_seconds if a.wall_seconds > 0 else None
            lines.append(
                f"| `{key[0]}` | {key[1]} | {key[2]} | "
                f"{a.wall_seconds:.2f}s | {b.wall_seconds:.2f}s | "
                f"{speedup:.2f}x | {a.peak_rss_bytes / 1e6:.0f} MB | "
                f"{b.peak_rss_bytes / 1e6:.0f} MB |"
            )
        lines.append("")

    # -- thread scaling ----------------------------------------------------
    lines += [
        "## Thread scaling",
        "",
        "Efficiency is measured speedup divided by thread count, against this "
        "machine's own 1-thread run. Values above 1.0 are flagged: they "
        "normally mean the 1-thread baseline was contended, not that the tool "
        "beat Amdahl.",
        "",
        "| benchmark | variant | input | threads | wall | speedup vs t=1 | "
        "efficiency | CPU util | peak RSS | spread |",
        "|---|---|---|---:|---:|---:|---:|---:|---:|---:|",
    ]

    for result in sorted(results, key=lambda r: (r.benchmark, r.variant,
                                                 r.input_key, r.threads)):
        best = result.best
        if not best:
            continue
        single = next(
            (r.best for r in results
             if r.benchmark == result.benchmark and r.variant == result.variant
             and r.input_key == result.input_key and r.threads == 1 and r.best),
            None,
        )
        if single and best.wall_seconds > 0:
            speedup = single.wall_seconds / best.wall_seconds
            efficiency = speedup / result.threads
            efficiency_text = f"{efficiency:.2f}" + (" ⚠" if efficiency > 1.05 else "")
            speedup_text = f"{speedup:.2f}x"
        else:
            speedup_text = efficiency_text = "—"

        spread = result.spread
        lines.append(
            f"| `{result.benchmark}` | {result.variant} | {result.input_key} | "
            f"{result.threads} | {best.wall_seconds:.2f}s | {speedup_text} | "
            f"{efficiency_text} | {best.cpu_utilisation:.1f}x | "
            f"{best.peak_rss_bytes / 1e6:.0f} MB | "
            f"{'—' if spread is None else format(spread * 100, '.1f') + '%'} |"
        )
    lines.append("")

    # -- I/O ---------------------------------------------------------------
    lines += [
        "## Disk I/O",
        "",
        "Block counts come from `ru_inblock`/`ru_oublock`, which count what "
        "the kernel actually moved. On Linux these are 512-byte blocks. On "
        "macOS they are usually **zero**, because the kernel does not charge "
        "page-cache-satisfied reads to the process -- that is reported as-is "
        "rather than converted into a MB/s figure that would be fiction.",
        "",
        "| benchmark | variant | input | threads | read blocks | write blocks | "
        "bytes written |",
        "|---|---|---|---:|---:|---:|---:|",
    ]
    for result in sorted(results, key=lambda r: (r.benchmark, r.variant, r.threads)):
        best = result.best
        if not best:
            continue
        lines.append(
            f"| `{result.benchmark}` | {result.variant} | {result.input_key} | "
            f"{result.threads} | {best.block_reads} | {best.block_writes} | "
            f"{best.output_bytes / 1e6:.1f} MB |"
        )
    lines.append("")

    skipped = [r for r in results if r.skipped]
    if skipped:
        lines += ["## Not measured", ""]
        for result in skipped:
            lines.append(
                f"- `{result.benchmark}` [{result.variant}] {result.input_key} "
                f"t={result.threads}: {result.skipped}"
            )
        lines.append("")
        lines.append(
            "A skipped configuration measured nothing. Any speedup claim above "
            "covers only the rows that actually ran."
        )
        lines.append("")

    return "\n".join(lines) + "\n"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="benchmark_audit.py",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--list", action="store_true",
                        help="show benchmarks and input tiers, then exit")
    parser.add_argument("--inputs", nargs="+", default=["small"],
                        choices=sorted(INPUT_TIERS),
                        help="input size tiers to run (default: small)")
    parser.add_argument("--threads", nargs="+", type=int,
                        default=[1, 2, 4, 8, 16, 32],
                        help="thread counts to sweep (default: 1..32)")
    parser.add_argument("--benchmarks", nargs="+", metavar="NAME",
                        help="restrict to these benchmark names")
    parser.add_argument("--baselines", action="store_true",
                        help="also measure the standard pipeline each tool "
                             "replaces, on the same input and core count")
    parser.add_argument("--repeats", type=int, default=1,
                        help="executions per configuration; the best is "
                             "reported and the spread across them is printed")
    parser.add_argument("--build-dir", type=Path, default=DEFAULT_BUILD_DIR)
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--work-dir", type=Path,
                        default=REPO_ROOT / "build-bench" / "work")
    parser.add_argument("--report", type=Path,
                        default=REPO_ROOT / "build-bench" / "audit",
                        help="output prefix; writes <prefix>.md and <prefix>.json")

    args = parser.parse_args(argv)

    if args.list:
        print("BENCHMARKS")
        for name, spec in BENCHMARKS.items():
            tools = ", ".join(spec.baseline_tools) or "none"
            print(f"  {name}")
            print(f"      {spec.summary}")
            print(f"      baseline tools: {tools}"
                  f"{'' if spec.scales else '   (single-threaded, not swept)'}")
        print("\nINPUT TIERS")
        for tier, spec in INPUT_TIERS.items():
            staged = sum(
                1 for role, rel in spec.items()
                if role != "description" and (args.data_dir / rel).exists()
            )
            total = len(spec) - 1
            print(f"  {tier:<8} {spec['description']}")
            print(f"           staged: {staged}/{total}")
        return 0

    # Ceilings above the machine's core count measure oversubscription, which
    # is a real thing to know but not what a scaling curve means. Dropping them
    # with a note beats silently plotting a flat tail as if it were scaling.
    cores = os.cpu_count() or 1
    requested = sorted(set(args.threads))
    args.threads = [t for t in requested if t <= cores] or [1]
    dropped = [t for t in requested if t > cores]
    if dropped:
        print(
            f"note: this host has {cores} cores; dropping thread counts "
            f"{dropped} rather than reporting oversubscription as scaling.",
            file=sys.stderr,
        )

    print(f"host: {platform.platform()} · {cores} cores")
    print(f"threads: {args.threads} · inputs: {args.inputs} · "
          f"repeats: {args.repeats}\n")

    results = run_matrix(args)
    if not results:
        print("no benchmarks ran", file=sys.stderr)
        return 1

    args.report.parent.mkdir(parents=True, exist_ok=True)
    markdown = render_report(results, args)
    args.report.with_suffix(".md").write_text(markdown)
    args.report.with_suffix(".json").write_text(
        json.dumps(
            {
                "host": {"platform": platform.platform(), "cores": cores},
                "generated": time.time(),
                "repeats": args.repeats,
                "results": [asdict(r) for r in results],
            },
            indent=2,
        )
        + "\n"
    )

    print(f"\nReport: {args.report.with_suffix('.md')}")
    measured = sum(1 for r in results if r.best)
    print(f"{measured}/{len(results)} configurations measured")
    return 0


if __name__ == "__main__":
    sys.exit(main())
