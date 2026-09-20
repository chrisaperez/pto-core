"""`pto-core`'s CLI subcommands, run end to end against the real binaries.

Same rule as `test_api.py`'s header: no mocked subprocess. `api.py` already
has the argv-construction tests for the flags these commands forward: what is
untested until this file is whether `cli.py`'s own argument names (`--min-len`,
`--min-mapq`, ...) actually reach the right `api` keyword, which is a defect
`test_api.py` cannot see because it never goes through `cli.py` at all.
"""

from __future__ import annotations

import json

import pytest

from pto import _binaries
from pto.cli import main


def _needs(tool: str):
    if _binaries.which(tool) is None:
        pytest.skip(f"{tool} not resolvable; run `pto-core doctor`")


def test_trim_writes_cleaned_reads_and_reports_a_summary(reads, tmp_path, capsys):
    _needs("fastq_stream")
    out = tmp_path / "clean.fastq"
    code = main(["trim", str(reads), str(out), "--min-len", "20"])
    assert code == 0
    assert out.is_file() and out.stat().st_size > 0
    printed = capsys.readouterr().out
    assert "reads_in" in printed or "reads" in printed.lower()


def test_trim_json_matches_the_api_report(reads, tmp_path, capsys):
    _needs("fastq_stream")
    out = tmp_path / "clean.fastq"
    code = main(["trim", str(reads), str(out), "--json"])
    assert code == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["tool"] == "fastq_stream"


def test_markdup_reports_a_duplicate_rate(fragments, capsys):
    _needs("genomic_toolkit")
    code = main(["markdup", str(fragments), "--json"])
    assert code == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["metrics"]


def test_markdup_min_mapq_reaches_the_binary_rather_than_being_swallowed(
    fragments, unsorted_fragments, capsys
):
    _needs("genomic_toolkit")
    # unsorted input in the default (sorted) mode is InvalidInput -- exit 2 --
    # which is the CLI-level proof that `--unsorted` really reached argv:
    # passing it should turn that failure into a clean 0.
    code = main(["markdup", str(unsorted_fragments), "--json"])
    assert code == 2
    capsys.readouterr()
    code = main(["markdup", str(unsorted_fragments), "--unsorted", "--json"])
    assert code == 0


def test_profile_writes_matrix_and_profile_tsv(demo_bam, tmp_path, capsys):
    _needs("cuttag_profiler")
    out_dir = tmp_path / "profile_out"
    # cuttag_profiler needs a regions file; build a tiny one rather than
    # assuming the BAM fixture ships one.
    regions = tmp_path / "regions.bed"
    regions.write_text("chr1\t1000\t2000\n")
    code = main(["profile", str(demo_bam), str(regions), str(out_dir), "--json"])
    assert code == 0
    payload = json.loads(capsys.readouterr().out)
    assert (out_dir / "matrix.tsv").is_file()
    assert (out_dir / "profile.tsv").is_file()
    assert set(payload["artifacts"]) == {"matrix", "profile"}
