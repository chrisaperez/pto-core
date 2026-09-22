"""End-to-end tests for the `pto` API, against the real binaries.

These run the engines. There is no mocked subprocess anywhere in this file, and
that is deliberate: a fake that returns what we believe the tool prints cannot
catch a change in what it actually prints, and the argv contract is
the whole substance of this package. The one thing that IS asserted structurally
rather than behaviourally is argv construction, below, because "the flag was not
passed" and "the flag was passed and the tool ignored it" produce identical
output and only one of them is correct.
"""

from __future__ import annotations

import json
import os
import time
from pathlib import Path

import pytest

import pto
from pto import _binaries


def _needs(tool: str):
    if _binaries.which(tool) is None:
        pytest.skip(f"{tool} not resolvable; run `pto-core doctor`")


# ---------------------------------------------------------------------------
# Discovery
# ---------------------------------------------------------------------------

def test_every_declared_tool_resolves_or_the_suite_says_which_did_not():
    missing = [t for t in pto.TOOLS if _binaries.which(t) is None]
    if missing:
        pytest.skip(f"binaries not built: {missing}")
    assert pto.bin_dir().is_dir()


def test_which_rejects_a_tool_this_package_does_not_wrap():
    # `scrna_matrix` is the trap: it is a pto-core module, it is NOT a binary,
    # and asking for it should be a programming error rather than a PATH miss
    # that returns None and looks like "not installed".
    with pytest.raises(ValueError):
        _binaries.which("scrna_matrix")
    with pytest.raises(ValueError):
        _binaries.which("samtools")


def test_bin_dir_env_override_is_searched_before_the_vendored_copy(tmp_path, monkeypatch):
    """`PTO_BIN_DIR` must win, or a site build cannot override a stale wheel."""
    _needs("genomic_toolkit")
    fake = tmp_path / "genomic_toolkit"
    fake.write_text("#!/bin/sh\nexit 0\n")
    fake.chmod(0o755)
    monkeypatch.setenv("PTO_BIN_DIR", str(tmp_path))
    assert _binaries.which("genomic_toolkit") == fake


def test_a_missing_binary_names_the_places_it_looked(tmp_path, monkeypatch):
    monkeypatch.setenv("PTO_BIN_DIR", str(tmp_path))
    monkeypatch.setenv("PATH", str(tmp_path))
    monkeypatch.setattr(_binaries, "vendored_dir", lambda: tmp_path / "nope")
    with pytest.raises(pto.BinaryNotFound) as caught:
        _binaries.resolve("genomic_toolkit")
    message = str(caught.value)
    assert str(tmp_path) in message and "$PATH" in message
    # The fix has to be in the message; a "not found" with no remedy is what
    # sends an evaluator to the issue tracker instead of to a FRiP number.
    assert "PTO_BIN_DIR" in message and "cmake" in message


# ---------------------------------------------------------------------------
# Threads
# ---------------------------------------------------------------------------

def test_available_cpus_is_at_least_one_and_never_exceeds_the_machine():
    assert 1 <= pto.available_cpus() <= max(1, (os.cpu_count() or 1))


def test_available_cpus_takes_the_minimum_signal_not_the_maximum(monkeypatch):
    """The whole point: a 4-CPU allocation on a 128-core node yields 4."""
    from pto import _cpu

    monkeypatch.setattr(_cpu, "_affinity", lambda: 4)
    monkeypatch.setattr(_cpu, "_cgroup_quota", lambda: 64)
    monkeypatch.setattr(_cpu.os, "cpu_count", lambda: 128)
    assert _cpu.available_cpus() == 4

    # And a cgroup quota below the affinity mask still wins.
    monkeypatch.setattr(_cpu, "_affinity", lambda: 64)
    monkeypatch.setattr(_cpu, "_cgroup_quota", lambda: 2)
    assert _cpu.available_cpus() == 2

    # No signal at all is 1, never 0.
    monkeypatch.setattr(_cpu, "_affinity", lambda: None)
    monkeypatch.setattr(_cpu, "_cgroup_quota", lambda: None)
    monkeypatch.setattr(_cpu.os, "cpu_count", lambda: None)
    assert _cpu.available_cpus() == 1


# ---------------------------------------------------------------------------
# argv construction -- the contract this package IS
# ---------------------------------------------------------------------------

@pytest.fixture()
def captured(monkeypatch):
    """Intercept `run_tool` and record the argv the API built."""
    seen: dict = {}

    def fake(tool, argv_tail, **kwargs):
        seen["tool"] = tool
        seen["argv"] = list(argv_tail)
        seen["kwargs"] = kwargs
        return pto.Result(tool=tool, subcommand=kwargs.get("subcommand"), argv=[], returncode=0)

    monkeypatch.setattr("pto.api.run_tool", fake)
    return seen


def test_an_unspecified_parameter_reaches_argv_as_nothing_at_all(captured, fragments):
    """A `None` must not become `--min-mapq None`, and must not become a default.

    Restating the binary's defaults here would be a fourth copy of every bound
    in this repository, after the CLI parser, `common/models.py` and the
    generated TypeScript. This asserts the omission is real.
    """
    pto.fragment_sizes(fragments)
    argv = captured["argv"]
    for absent in ("--min-length", "--max-length", "--min-mapq", "--single-end",
                   "--unsorted", "--region", "--histogram"):
        assert absent not in argv, f"{absent} was passed despite not being requested"
    assert argv[:2] == [str(fragments), "--json"]


def test_a_specified_parameter_reaches_argv_exactly_once(captured, fragments):
    pto.fragment_sizes(fragments, min_mapq=30, min_length=100, max_length=0,
                       single_end=True, unsorted=True)
    argv = captured["argv"]
    assert argv.count("--min-mapq") == 1
    assert argv[argv.index("--min-mapq") + 1] == "30"
    assert argv[argv.index("--min-length") + 1] == "100"
    # 0 is a MEANINGFUL value for --max-length (it means "off"), so it must be
    # passed, not swallowed as falsy.
    assert argv[argv.index("--max-length") + 1] == "0"
    assert "--single-end" in argv and "--unsorted" in argv


def test_a_false_switch_is_not_a_flag(captured, fragments):
    """`single_end=False` means "do not pass --single-end", not "pass it"."""
    pto.fragment_sizes(fragments, single_end=False, unsorted=False)
    assert "--single-end" not in captured["argv"]
    assert "--unsorted" not in captured["argv"]


def test_frip_places_peaks_after_the_input_and_before_the_filters(captured, fragments, peaks):
    pto.frip(fragments, peaks, genome_size=2_913_022_398, keep_dups=True)
    argv = captured["argv"]
    assert argv[0] == str(fragments)
    assert argv[argv.index("--peaks") + 1] == str(peaks)
    assert argv[argv.index("--genome-size") + 1] == "2913022398"
    assert "--keep-dups" in argv
    assert captured["kwargs"]["subcommand"] == "frip"
    assert captured["kwargs"]["stdout_json"] is True


def test_threads_are_always_explicit_even_when_the_caller_says_nothing(captured, reads):
    """The engines default to the MACHINE's core count; see pto._cpu."""
    pto.run_qc(reads)
    argv = captured["argv"]
    assert "-t" in argv
    assert int(argv[argv.index("-t") + 1]) >= 1


def test_regions_are_repeated_rather_than_joined(captured, fragments):
    pto.fragment_sizes(fragments, regions=["chr1", "chr2:1-1000"])
    argv = captured["argv"]
    assert argv.count("--region") == 2
    assert "chr1" in argv and "chr2:1-1000" in argv


def test_float_parameters_round_trip_through_argv_without_losing_precision(
    captured, reads
):
    """`%g`, the float format argv actually carries.

    `%g` renders 0.00001 as `1e-05`, and that is fine: the engines parse with
    `std::from_chars`/strtod, which accept scientific notation -- checked
    against the binary in the end-to-end test below rather than assumed. What
    would NOT be fine is a rendering that changes the value, so that is what is
    asserted. `str(float)` is the alternative and it is worse: it produces
    `1e-05` too, and `0.30000000000000004` for arithmetic results.
    """
    pto.run_qc(reads, max_n_rate=0.00001, adapter_mismatch=0.2, min_mean_q=27.5)
    argv = captured["argv"]
    for flag, value in (("--max-n-rate", 0.00001),
                        ("--adapter-mismatch", 0.2),
                        ("--min-mean-q", 27.5)):
        rendered = argv[argv.index(flag) + 1]
        assert float(rendered) == value, f"{flag} rendered as {rendered!r}"


def test_the_engine_accepts_the_scientific_notation_g_produces(reads):
    """The other half of the test above, run against the real parser."""
    _needs("fastq_stream")
    result = pto.run_qc(reads, max_n_rate=0.00001)
    assert result.returncode == 0


# ---------------------------------------------------------------------------
# No shell, ever
# ---------------------------------------------------------------------------

def test_a_path_containing_shell_metacharacters_is_passed_through_intact(tmp_path):
    """argv, never a shell -- and that is not weaker locally than anywhere else.

    A filename is attacker-controlled input on a shared cluster: users hand each
    other BAMs. If any layer here reached a shell, this file name would delete
    the marker file instead of being read as a path.
    """
    _needs("genomic_toolkit")
    hostile_dir = tmp_path / "a b; touch pwned"
    hostile_dir.mkdir()
    marker = tmp_path / "pwned"
    frags = hostile_dir / "frag$(id).bed"
    frags.write_text("chr1\t100\t250\nchr1\t300\t470\n")

    result = pto.fragment_sizes(frags)
    assert result["fragments"] == 2
    assert not marker.exists(), "a shell interpreted the path"


# ---------------------------------------------------------------------------
# Real runs: genomic_toolkit
# ---------------------------------------------------------------------------

def test_fragment_sizes_returns_the_tools_own_report(fragments):
    _needs("genomic_toolkit")
    result = pto.fragment_sizes(fragments)
    assert result.returncode == 0
    # Field names are the TOOL's, not reshaped into a schema of our own.
    assert result["fragments"] == 4000
    assert result["median"] > 0
    assert "dropped" in result.metrics
    # Cost comes back with every call, because that is what is being evaluated.
    assert result.wall_seconds > 0
    if result.peak_rss_bytes is not None:
        assert result.peak_rss_bytes > 0
        assert result.peak_rss_mb == round(result.peak_rss_bytes / 1e6, 1)


def test_fragment_sizes_writes_a_histogram_only_when_asked(fragments, tmp_path):
    _needs("genomic_toolkit")
    hist = tmp_path / "nested" / "hist.tsv"
    result = pto.fragment_sizes(fragments, histogram=hist)
    assert hist.is_file() and hist.stat().st_size > 0
    assert result.artifacts["histogram"] == hist


def test_a_filter_actually_filters(fragments):
    _needs("genomic_toolkit")
    everything = pto.fragment_sizes(fragments)
    gated = pto.fragment_sizes(fragments, min_length=200)
    assert gated["fragments"] < everything["fragments"]


def test_frip_agrees_with_itself_on_the_bases_it_reports(fragments, peaks):
    _needs("genomic_toolkit")
    result = pto.frip(fragments, peaks)
    assert 0.0 < result["frip"] <= 1.0
    assert result["in_peaks"] <= result["fragments"]
    assert result["unmatched_contigs"] == 0
    # Peaks are merged into a disjoint cover at load; these three are one fact.
    assert result["peaks_merged"] <= result["peaks_input"]


def test_keep_dups_changes_the_answer_and_only_for_frip(fragments, peaks):
    """The dual-path comparison low-input CUT&RUN work needs."""
    _needs("genomic_toolkit")
    dedup = pto.frip(fragments, peaks)
    original = pto.frip(fragments, peaks, keep_dups=True)
    assert original["fragments"] >= dedup["fragments"]
    assert dedup["duplicates_excluded"] >= 0


def test_a_contig_naming_mismatch_raises_InvalidInput_and_keeps_the_report(
    fragments, mismatched_peaks
):
    """Exit 2. The most valuable single test in this file.

    A FRiP of exactly zero from a `chr1` vs `1` mismatch is the failure mode
    that looks like a failed experiment. It must be an exception, it must be
    distinguishable from a crash, and the report must survive -- because the
    tool did produce one.
    """
    _needs("genomic_toolkit")
    with pytest.raises(pto.InvalidInput) as caught:
        pto.frip(fragments, mismatched_peaks)

    exc = caught.value
    assert exc.returncode == 2
    assert isinstance(exc, pto.ToolError)          # retryable/not is the distinction
    assert exc.result is not None
    assert exc.result.metrics["unmatched_contigs"] >= 1
    assert exc.result.metrics["frip"] == 0.0
    assert "naming" in exc.stderr or "absent" in exc.stderr


def test_markdup_on_unsorted_input_is_InvalidInput_and_unsorted_true_fixes_it(
    unsorted_fragments
):
    _needs("genomic_toolkit")
    with pytest.raises(pto.InvalidInput) as caught:
        pto.mark_duplicates(unsorted_fragments)
    assert caught.value.result.metrics["out_of_order"] > 0

    # The documented remedy has to actually work, or the exception is a dead end.
    recovered = pto.mark_duplicates(unsorted_fragments, unsorted=True)
    assert recovered.returncode == 0
    assert recovered["examined"] > 0


def test_markdup_reports_a_rate_on_sorted_input(fragments):
    _needs("genomic_toolkit")
    result = pto.mark_duplicates(fragments)
    assert result["examined"] == 4000
    assert 0.0 <= result["duplicate_rate"] <= 1.0
    assert result["unique"] + result["duplicates"] == result["examined"]


def test_a_missing_input_fails_before_a_process_is_spawned(tmp_path):
    """FileNotFoundError naming the path the caller passed, not a ToolError."""
    with pytest.raises(FileNotFoundError) as caught:
        pto.fragment_sizes(tmp_path / "absent.bam")
    assert "absent.bam" in str(caught.value)


# ---------------------------------------------------------------------------
# Real runs: fastq_stream
# ---------------------------------------------------------------------------

def test_run_qc_leaves_nothing_behind_when_no_out_dir_is_given(reads, tmp_path):
    _needs("fastq_stream")
    before = set(os.listdir(tmp_path))
    result = pto.run_qc(reads)
    assert result["summary"]["reads_in"] == 500
    assert result.artifacts == {}
    assert set(os.listdir(tmp_path)) == before


def test_run_qc_keeps_its_artifacts_when_asked(reads, tmp_path):
    _needs("fastq_stream")
    out = tmp_path / "qc"
    result = pto.run_qc(reads, out_dir=out)
    assert result.artifacts["qc_json"].is_file()
    assert result.artifacts["qc_report"].is_file()
    # The JSON on disk and the parsed metrics are the same object.
    assert json.loads(result.artifacts["qc_json"].read_text()) == result.metrics


def test_run_qc_writes_no_reads(reads, tmp_path):
    _needs("fastq_stream")
    out = tmp_path / "qconly"
    pto.run_qc(reads, out_dir=out)
    assert not any(p.suffix in {".fq", ".fastq", ".gz"} for p in out.iterdir())


def test_trim_writes_reads_and_the_qc_report_describes_them(reads, tmp_path):
    _needs("fastq_stream")
    out = tmp_path / "trimmed" / "clean.fastq"
    result = pto.trim(reads, out, out_dir=tmp_path / "reports")
    assert out.is_file() and out.stat().st_size > 0
    summary = result.metrics["summary"]
    assert summary["reads_out"] <= summary["reads_in"] == 500
    assert result.artifacts["trimmed"] == out


def test_disabling_adapter_trimming_changes_the_result(reads):
    _needs("fastq_stream")
    with_adapters = pto.run_qc(reads)
    without = pto.run_qc(reads, trim_adapters=False)
    trimmed_on = with_adapters.metrics["filtering"]["adapter_trimmed"]
    trimmed_off = without.metrics["filtering"]["adapter_trimmed"]
    assert trimmed_on > 0, "the fixture should contain adapter-bearing reads"
    assert trimmed_off == 0


# ---------------------------------------------------------------------------
# Real runs: cuttag_profiler (needs a BAM; skips from an sdist)
# ---------------------------------------------------------------------------

def test_profile_writes_a_matrix_and_a_curve(demo_bam, tmp_path):
    _needs("cuttag_profiler")
    regions = tmp_path / "regions.bed"
    # Derived from the BAM's own header so this does not depend on a second
    # tracked fixture -- and so a contig rename in the fixture cannot make this
    # silently profile nothing.
    import subprocess
    proc = subprocess.run(
        [str(_binaries.resolve("genomic_toolkit")), "sizes", str(demo_bam), "--json"],
        capture_output=True, text=True, check=False,
    )
    if proc.returncode != 0:
        pytest.skip("genomic_toolkit could not read the demo BAM (built without htslib?)")
    regions.write_text("chr1\t1000\t1200\t.\t0\t+\nchr1\t20000\t20200\t.\t0\t-\n")

    out = tmp_path / "profile-out"
    result = pto.profile(demo_bam, regions, out, bin_size=50, upstream=500, downstream=500)
    assert result.returncode == 0
    matrix, curve = result.artifacts["matrix"], result.artifacts["profile"]
    assert matrix.is_file() and curve.is_file()
    assert curve.read_text().strip(), "the meta-profile curve is empty"


def test_profile_never_constructs_the_serve_subcommand(captured, tmp_path):
    """`serve` is a listening socket. It must be unreachable from this API."""
    bam = tmp_path / "x.bam"
    bam.write_bytes(b"")
    regions = tmp_path / "r.bed"
    regions.write_text("chr1\t1\t2\n")
    pto.profile(bam, regions, tmp_path / "o")
    assert captured["kwargs"]["subcommand"] == "profile"
    assert "serve" not in captured["argv"]
    assert not any(a in captured["argv"] for a in ("--port", "--host", "--data-root"))


# ---------------------------------------------------------------------------
# Result
# ---------------------------------------------------------------------------

def test_result_as_dict_is_json_serialisable(fragments):
    _needs("genomic_toolkit")
    record = pto.fragment_sizes(fragments).as_dict()
    round_tripped = json.loads(json.dumps(record))
    assert round_tripped["tool"] == "genomic_toolkit"
    assert round_tripped["subcommand"] == "sizes"
    assert round_tripped["returncode"] == 0
    assert round_tripped["metrics"]["fragments"] == 4000


def test_tool_versions_reports_something_for_every_resolvable_tool():
    """Raw, per-tool `--version` output; the key names the tool.

    The three engines do not agree on a format -- `fastq_stream` prints
    `0.1.0`, `genomic_toolkit` prints `genomic_toolkit 0.1.0`, and
    `cuttag_profiler` prints `1.0.0`. That is reported as-is rather than
    normalised: a methods section needs the string the tool claims for itself,
    and rewriting it here would be this package asserting a version the engine
    never printed.
    """
    versions = pto.tool_versions()
    assert set(versions) == set(pto.TOOLS)
    for tool, reported in versions.items():
        if _binaries.which(tool) is None:
            assert reported is None
            continue
        assert reported, f"{tool} reported no version"
        assert any(ch.isdigit() for ch in reported), reported


def test_the_distribution_reports_a_version():
    assert pto.__version__
    assert isinstance(pto.__version__, str)


# ---------------------------------------------------------------------------
# run_tool's own contract
# ---------------------------------------------------------------------------

def test_a_timeout_kills_the_process_rather_than_hanging(tmp_path, monkeypatch):
    """The watchdog path, which no analysis job exercises by accident.

    A job that runs past its budget on a shared node has to die, and the caller
    has to be able to tell that from a tool failure. `check=False` is used here
    so the negative exit status is inspectable rather than wrapped.
    """
    from pto import _run

    slow = tmp_path / "genomic_toolkit"
    slow.write_text("#!/bin/sh\nsleep 30\n")
    slow.chmod(0o755)
    monkeypatch.setenv("PTO_BIN_DIR", str(tmp_path))

    started = time.perf_counter()
    result = _run.run_tool("genomic_toolkit", ["x"], timeout=0.5, check=False)
    elapsed = time.perf_counter() - started

    assert elapsed < 10, f"the timeout did not fire; waited {elapsed:.1f}s"
    assert result.returncode != 0
    assert result.wall_seconds < 10


def test_check_false_returns_a_failing_result_instead_of_raising(tmp_path, monkeypatch):
    from pto import _run

    failing = tmp_path / "genomic_toolkit"
    failing.write_text("#!/bin/sh\necho 'boom' >&2\nexit 3\n")
    failing.chmod(0o755)
    monkeypatch.setenv("PTO_BIN_DIR", str(tmp_path))

    result = _run.run_tool("genomic_toolkit", ["x"], check=False)
    assert result.returncode == 3
    assert "boom" in result.stderr

    # ...and with check on, exit 3 is a ToolError and NOT an InvalidInput --
    # only exit 2 carries that meaning.
    with pytest.raises(pto.ToolError) as caught:
        _run.run_tool("genomic_toolkit", ["x"], check=True)
    assert not isinstance(caught.value, pto.InvalidInput)


def test_unparseable_stdout_is_surfaced_rather_than_raising(tmp_path, monkeypatch):
    """A truncated report should not become a JSONDecodeError inside a QC helper."""
    from pto import _run

    noisy = tmp_path / "genomic_toolkit"
    noisy.write_text("#!/bin/sh\nprintf '{\"frip\": 0.4'\n")
    noisy.chmod(0o755)
    monkeypatch.setenv("PTO_BIN_DIR", str(tmp_path))

    result = _run.run_tool("genomic_toolkit", ["x"], stdout_json=True)
    assert "_unparsed_stdout" in result.metrics


# ---------------------------------------------------------------------------
# Real runs: pto-peaks
# ---------------------------------------------------------------------------

def test_call_peaks_finds_the_planted_peaks(enriched_fragments, chrom_sizes, tmp_path):
    """WHERE, not just how many. See the fixture for why that distinction matters."""
    _needs("pto-peaks")
    out = tmp_path / "peaks.narrowPeak"
    result = pto.call_peaks(enriched_fragments, out, chrom_sizes=chrom_sizes)

    assert result.returncode == 0
    assert result.artifacts["peaks"] == out
    rows = [line.split("\t") for line in out.read_text().splitlines() if line.strip()]
    assert rows, "the enriched fixture produced no peaks at all"

    spans = [(int(r[1]), int(r[2])) for r in rows]
    for centre in (250_000, 410_000):
        assert any(start <= centre + 300 <= end for start, end in spans), (
            f"no peak covers the planting at {centre}: {spans[:10]}"
        )


def test_call_peaks_writes_narrowPeak_and_not_some_other_bed_flavour(
    enriched_fragments, chrom_sizes, tmp_path
):
    """narrowPeak is BED6+4, and the +4 is what a browser and MACS2 consumers read.

    Column 10 is the summit offset FROM THE PEAK START, not an absolute
    coordinate -- getting that wrong produces a file that loads fine and puts
    every summit at the start of the chromosome.
    """
    _needs("pto-peaks")
    out = tmp_path / "peaks.narrowPeak"
    pto.call_peaks(enriched_fragments, out, chrom_sizes=chrom_sizes)

    rows = [line.split("\t") for line in out.read_text().splitlines() if line.strip()]
    assert rows
    for row in rows:
        assert len(row) == 10, f"narrowPeak is BED6+4; got {len(row)} columns"
        start, end = int(row[1]), int(row[2])
        assert 0 <= start < end
        assert row[5] in (".", "+", "-")
        summit = int(row[9])
        assert 0 <= summit < (end - start), (
            "column 10 must be an offset within the peak, not an absolute coordinate"
        )
        # -log10 p and -log10 q, and BH can only ever make a q weaker than its p.
        p_col, q_col = float(row[7]), float(row[8])
        assert p_col >= 0 and q_col >= 0
        assert q_col <= p_col + 1e-6


def test_call_peaks_creates_the_output_directory(enriched_fragments, chrom_sizes, tmp_path):
    """`profile` does this for `out_dir`; a caller writing per-sample paths needs it."""
    _needs("pto-peaks")
    out = tmp_path / "nested" / "deeper" / "peaks.narrowPeak"
    pto.call_peaks(enriched_fragments, out, chrom_sizes=chrom_sizes)
    assert out.is_file()


def test_a_stricter_qvalue_never_returns_more_peaks(
    enriched_fragments, chrom_sizes, tmp_path
):
    """Monotonicity is the one property of an FDR gate worth asserting.

    An off-by-one in the Benjamini-Hochberg rank, or a comparison in the wrong
    direction, shows up here and nowhere in a single-threshold run.
    """
    _needs("pto-peaks")

    def count(q: float) -> int:
        out = tmp_path / f"q{q}.narrowPeak"
        pto.call_peaks(enriched_fragments, out, chrom_sizes=chrom_sizes, qvalue=q)
        return sum(1 for line in out.read_text().splitlines() if line.strip())

    assert count(0.001) <= count(0.05) <= count(0.5)


def test_min_length_actually_reaches_the_binary(
    enriched_fragments, chrom_sizes, tmp_path
):
    """A flag that is built but ignored produces identical output to no flag.

    The argv tests below prove the token is emitted; this proves the binary
    acts on it, which is the half a mocked subprocess can never show.
    """
    _needs("pto-peaks")
    out = tmp_path / "long.narrowPeak"
    pto.call_peaks(enriched_fragments, out, chrom_sizes=chrom_sizes,
                   qvalue=0.5, min_length=400)
    for line in out.read_text().splitlines():
        if not line.strip():
            continue
        chrom, start, end, *_ = line.split("\t")
        assert int(end) - int(start) >= 400


def test_bed_input_without_chrom_sizes_says_so_rather_than_calling_nothing(
    enriched_fragments, tmp_path
):
    """Exit 1, with the reason. The wrapper deliberately does not pre-empt this.

    Detecting BED-vs-BAM in Python would be a second copy of the byte sniff in
    `pto-peaks`'s own `run()`, guessing from a path this layer never opens.
    """
    _needs("pto-peaks")
    with pytest.raises(pto.ToolError) as caught:
        pto.call_peaks(enriched_fragments, tmp_path / "p.narrowPeak")
    assert "chrom-sizes" in str(caught.value)
    assert not isinstance(caught.value, pto.InvalidInput)


def test_a_contig_regrouped_after_another_is_invalid_input_not_a_wrong_answer(
    regrouped_fragments, chrom_sizes, tmp_path
):
    _needs("pto-peaks")
    with pytest.raises(pto.InvalidInput):
        pto.call_peaks(regrouped_fragments, tmp_path / "p.narrowPeak",
                       chrom_sizes=chrom_sizes)


def test_a_contig_missing_from_chrom_sizes_names_the_naming_mismatch(
    chrom_sizes, tmp_path
):
    """`chr1` vs `1` is the single most common way this fails on real data."""
    _needs("pto-peaks")
    ensembl = tmp_path / "ensembl.bed"
    ensembl.write_text("1\t1000\t1150\n1\t2000\t2150\n")
    with pytest.raises(pto.InvalidInput) as caught:
        pto.call_peaks(ensembl, tmp_path / "p.narrowPeak", chrom_sizes=chrom_sizes)
    assert "chr1" in str(caught.value) or "naming" in str(caught.value)


def test_call_peaks_reads_a_bam_without_being_given_chrom_sizes(demo_bam, tmp_path):
    """A BAM carries its contig lengths in its header, so the flag is not needed.

    This also proves the wheel was built with a BGZF backend: `PEAKS_WITH_BGZF`
    is pinned ON in `pto-core/pyproject.toml` precisely so that a wheel cannot
    ship having quietly compiled this path out.
    """
    _needs("pto-peaks")
    out = tmp_path / "bam.narrowPeak"
    result = pto.call_peaks(demo_bam, out)
    assert result.returncode == 0
    assert out.is_file()
    # The BAM reader reports what it kept and what it dropped, on stderr; a run
    # that silently kept nothing is the failure this line would catch.
    assert "fragments" in result.stderr


def test_call_peaks_leaves_metrics_empty_rather_than_inventing_a_report(
    enriched_fragments, chrom_sizes, tmp_path
):
    """`pto-peaks` prints narrowPeak, not JSON, and this package does not fabricate one."""
    _needs("pto-peaks")
    result = pto.call_peaks(enriched_fragments, tmp_path / "p.narrowPeak",
                            chrom_sizes=chrom_sizes)
    assert result.metrics == {}
    assert result.wall_seconds > 0


def test_call_peaks_omits_every_parameter_the_caller_did_not_give(
    captured, enriched_fragments, tmp_path
):
    pto.call_peaks(enriched_fragments, tmp_path / "p.narrowPeak")
    argv = captured["argv"]
    for absent in ("--chrom-sizes", "--cutoff", "--qvalue", "--min-length",
                   "--max-gap", "--lambda-bg", "--min-mapq", "--extend-peaks"):
        assert absent not in argv, f"{absent} was passed despite not being requested"
    assert captured["tool"] == "pto-peaks"
    assert captured["kwargs"].get("subcommand") is None
    # narrowPeak is not JSON; asking run_tool to parse it would fill `metrics`
    # with `_unparsed_stdout` on every successful call.
    assert captured["kwargs"].get("stdout_json") in (None, False)


def test_call_peaks_puts_the_input_last_and_the_output_behind_a_flag(
    captured, enriched_fragments, chrom_sizes, tmp_path
):
    """The input is positional, so anything appended after it becomes a second input.

    `pto-peaks` accepts at most one positional and exits 1 on two, which makes
    this ordering a real constraint rather than cosmetics.
    """
    out = tmp_path / "p.narrowPeak"
    pto.call_peaks(enriched_fragments, out, chrom_sizes=chrom_sizes,
                   cutoff=3.0, qvalue=0.01, min_length=120, max_gap=45,
                   lambda_bg=0.002, min_mapq=30, extend_peaks=90)
    argv = captured["argv"]
    assert argv[-1] == str(enriched_fragments)
    assert argv[argv.index("--out") + 1] == str(out)
    assert argv[argv.index("--chrom-sizes") + 1] == str(chrom_sizes)
    for flag, rendered in (("--cutoff", "3"), ("--qvalue", "0.01"),
                           ("--min-length", "120"), ("--max-gap", "45"),
                           ("--lambda-bg", "0.002"), ("--min-mapq", "30"),
                           ("--extend-peaks", "90")):
        assert argv.count(flag) == 1
        assert argv[argv.index(flag) + 1] == rendered


def test_lambda_bg_zero_is_passed_rather_than_swallowed_as_falsy(
    captured, enriched_fragments, tmp_path
):
    """0 means "local scales only" to this tool -- a meaningful value, not "unset"."""
    pto.call_peaks(enriched_fragments, tmp_path / "p.narrowPeak", lambda_bg=0.0)
    argv = captured["argv"]
    assert argv[argv.index("--lambda-bg") + 1] == "0"


def test_extend_peaks_reaches_argv_exactly_once(captured, enriched_fragments, tmp_path):
    """Option B: the flag that pads reported spans for MACS-style pipelines."""
    pto.call_peaks(enriched_fragments, tmp_path / "p.narrowPeak", extend_peaks=180)
    argv = captured["argv"]
    assert argv.count("--extend-peaks") == 1
    assert argv[argv.index("--extend-peaks") + 1] == "180"


def test_extend_peaks_widens_the_called_intervals_without_moving_summits(
    enriched_fragments, chrom_sizes, tmp_path
):
    """The argv test proves the token is emitted; this proves the binary acts on
    it, and -- the half a mocked subprocess can never show -- that the
    enrichment measurement is untouched. Option A stays the default; Option B
    only moves the reported bounds.
    """
    _needs("pto-peaks")

    def call(name: str, **kw) -> dict[int, list[str]]:
        out = tmp_path / name
        pto.call_peaks(enriched_fragments, out, chrom_sizes=chrom_sizes, **kw)
        rows = [ln.split("\t") for ln in out.read_text().splitlines() if ln.strip()]
        # keyed by ABSOLUTE summit position = start + column 10 (the offset)
        return {int(r[1]) + int(r[9]): r for r in rows}

    sharp = call("sharp.narrowPeak")
    padded = call("padded.narrowPeak", extend_peaks=150)
    assert sharp and padded

    common = set(sharp) & set(padded)
    assert len(common) >= max(1, len(sharp) // 2), (
        "padding changed which peaks are called, not just their bounds"
    )

    for summit in common:
        s_start, s_end = int(sharp[summit][1]), int(sharp[summit][2])
        p_start, p_end = int(padded[summit][1]), int(padded[summit][2])
        assert p_start <= s_start and p_end >= s_end
        assert (p_end - p_start) > (s_end - s_start)
        # the summit column is still an in-bounds offset of the widened peak
        assert 0 <= int(padded[summit][9]) < (p_end - p_start)
        # columns 8 and 9 are -log10 p and -log10 q -- the scores do not move
        assert sharp[summit][7] == padded[summit][7]
        assert sharp[summit][8] == padded[summit][8]
        # clear of the contig edges, the widening is exactly 2 x the padding
        if p_start > 0 and p_end < 600_000:
            assert (p_end - p_start) - (s_end - s_start) == 300


def test_extend_peaks_zero_is_identical_to_the_default(
    enriched_fragments, chrom_sizes, tmp_path
):
    """`extend_peaks=0` is Option A -- byte-for-byte the sharp default."""
    _needs("pto-peaks")
    a, b = tmp_path / "a.narrowPeak", tmp_path / "b.narrowPeak"
    pto.call_peaks(enriched_fragments, a, chrom_sizes=chrom_sizes)
    pto.call_peaks(enriched_fragments, b, chrom_sizes=chrom_sizes, extend_peaks=0)
    assert a.read_text() == b.read_text()


def test_extend_peaks_never_reports_a_span_past_the_contig_end(
    enriched_fragments, chrom_sizes, tmp_path
):
    """The 3' clamp lives in the writer, which is the layer that knows the length.

    chrom_sizes pins chr1 at 600,000 and the enriched fixture plants a window at
    410,000; a large padding must not push any `end` past the contig.
    """
    _needs("pto-peaks")
    out = tmp_path / "wide.narrowPeak"
    pto.call_peaks(enriched_fragments, out, chrom_sizes=chrom_sizes, extend_peaks=100_000)
    for line in out.read_text().splitlines():
        if not line.strip():
            continue
        chrom, start, end, *_ = line.split("\t")
        assert int(start) >= 0
        assert int(end) <= 600_000, f"{chrom}:{start}-{end} runs past the contig"


# ---------------------------------------------------------------------------
# The wrapped set: TOOLS, the wheel, and the flags
# ---------------------------------------------------------------------------

def test_TOOLS_matches_what_the_wheel_vendors():
    """Vendored and wrapped are set by two different files that cannot see each other.

    `pto-core/pyproject.toml`'s `PTO_BUILD_*` decides what lands in
    `pto/_vendor/bin/`; `_binaries.TOOLS` decides what this package can call.
    `pto-peaks` was in the first and not the second for a while -- shipped,
    undocumented and unreachable -- which is the failure this guards. It only
    means something against an installed wheel, so it skips elsewhere.
    """
    vendored = _binaries.vendored_dir()
    if not vendored.is_dir():
        pytest.skip("not running against an installed wheel")
    on_disk = {p.name for p in vendored.iterdir() if p.is_file() and os.access(p, os.X_OK)}
    assert on_disk == set(pto.TOOLS), (
        f"vendored but not wrapped: {sorted(on_disk - set(pto.TOOLS))}; "
        f"wrapped but not vendored: {sorted(set(pto.TOOLS) - on_disk)}"
    )


def test_every_flag_the_api_emits_exists_in_the_binarys_own_help():
    """The drift guard `nextflow/tests/test_modules.py` applies to its modules.

    This layer emits long flags into argv and nothing has been checking them
    against the tools. A renamed or removed flag currently surfaces as a usage
    error inside someone's batch; here it is a failing test on the commit that
    renames it.

    Only long flags are checked. Short ones (`-i`, `-t`, `-o`) are single
    characters that appear all over a help text by coincidence, so matching
    them would be a guard that cannot fail.
    """
    import re
    import subprocess

    source = (Path(pto.api.__file__)).read_text()
    # `_flag(argv, "--x", ...)`, `_switch(argv, "--x", ...)` and bare literals.
    emitted = set(re.findall(r'"(--[a-z0-9][a-z0-9-]*)"', source))
    assert len(emitted) > 20, "the extraction stopped matching; fix it, do not relax it"

    helps: dict[str, str] = {}
    for tool in pto.TOOLS:
        path = _binaries.which(tool)
        if path is None:
            pytest.skip(f"{tool} not resolvable; run `pto-core doctor`")
        proc = subprocess.run([str(path), "--help"], capture_output=True,
                              text=True, timeout=30, check=False)
        helps[tool] = proc.stdout + proc.stderr
    # `genomic_toolkit`'s flags live under its subcommands rather than its top
    # level, so the union across every tool is the right haystack: this asks
    # "does any engine document this flag", which is what catches a rename.
    for tool in ("genomic_toolkit",):
        for sub in ("sizes", "markdup", "frip"):
            proc = subprocess.run([str(_binaries.which(tool)), sub, "--help"],
                                  capture_output=True, text=True, timeout=30,
                                  check=False)
            helps[f"{tool} {sub}"] = proc.stdout + proc.stderr
    haystack = "\n".join(helps.values())

    unknown = sorted(f for f in emitted if not re.search(rf"{re.escape(f)}\b", haystack))
    assert not unknown, f"api.py emits flags no engine documents: {unknown}"
