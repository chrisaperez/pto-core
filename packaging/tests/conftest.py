"""Fixtures for the `pto` wheel suite.

Every fixture here is GENERATED AS TEXT rather than shipped as a binary, and
that is what lets this suite run against an installed wheel and from an sdist:
`modules/cuttag_profiler/data` is 99 MB of benchmark BAM and is excluded from
the sdist, so a suite that depended on it would silently skip exactly where a
`pip install --no-binary` most needs checking.

It is also what pto-core's own suites do -- `genomic_toolkit`'s whole test suite
drives the BEDPE / fragment-BED parsers on purpose so it passes with or without
htslib. The BAM-only path (`cuttag_profiler profile`) skips when no BAM fixture
is present, which is the honest outcome rather than a fabricated one.
"""

from __future__ import annotations

import random
from pathlib import Path

import pytest

# The tracked demo BAM, when this runs from a checkout rather than an sdist.
_CHECKOUT_BAM = (
    Path(__file__).resolve().parents[2] / "modules" / "cuttag_profiler" / "data" / "demo.bam"
)


@pytest.fixture(scope="session")
def fragments(tmp_path_factory) -> Path:
    """A coordinate-sorted fragment BED with a realistic size distribution.

    Sorted, because `genomic_toolkit markdup`'s default duplicate mode is exact
    on coordinate-sorted input and reports `out_of_order` (exit 2) otherwise --
    both behaviours are tested, so the fixture has to be the well-formed one.
    """
    path = tmp_path_factory.mktemp("pto-fixtures") / "fragments.bed"
    rng = random.Random(20260908)
    rows = []
    for _ in range(4000):
        start = rng.randrange(1_000, 400_000)
        # Nucleosome-ish multimodality, so median/p10/p90 are not degenerate.
        length = rng.choice([70, 90, 120, 150, 170, 185, 210, 260, 340, 480])
        rows.append((start, start + length))
    rows.sort()
    path.write_text("".join(f"chr1\t{s}\t{e}\n" for s, e in rows))
    return path


@pytest.fixture(scope="session")
def unsorted_fragments(tmp_path_factory) -> Path:
    """Deliberately out of coordinate order, to reach the exit-2 markdup path."""
    path = tmp_path_factory.mktemp("pto-fixtures") / "unsorted.bed"
    rng = random.Random(11)
    rows = [(rng.randrange(1_000, 400_000), 150) for _ in range(2000)]
    path.write_text("".join(f"chr1\t{s}\t{s + n}\n" for s, n in rows))
    return path


@pytest.fixture(scope="session")
def peaks(tmp_path_factory) -> Path:
    """Peaks on `chr1`, overlapping the fragments."""
    path = tmp_path_factory.mktemp("pto-fixtures") / "peaks.bed"
    path.write_text(
        "chr1\t1000\t60000\n"
        "chr1\t90000\t150000\n"
        "chr1\t200000\t260000\n"
    )
    return path


@pytest.fixture(scope="session")
def mismatched_peaks(tmp_path_factory) -> Path:
    """Peaks whose contig naming disagrees with the fragments -- the exit-2 case."""
    path = tmp_path_factory.mktemp("pto-fixtures") / "ensembl_named.bed"
    path.write_text("1\t1000\t60000\n1\t90000\t150000\n")
    return path


@pytest.fixture(scope="session")
def reads(tmp_path_factory) -> Path:
    """A small FASTQ with adapter-bearing and low-quality reads in it."""
    path = tmp_path_factory.mktemp("pto-fixtures") / "reads.fastq"
    rng = random.Random(4242)
    adapter = "AGATCGGAAGAGC"
    records = []
    for i in range(500):
        body = "".join(rng.choice("ACGT") for _ in range(60))
        if i % 5 == 0:
            body = body[:40] + adapter + body[: 60 - 40 - len(adapter)]
        # A low-quality 3' tail on some reads so trimming has work to do.
        qual = "I" * len(body) if i % 3 else "I" * (len(body) - 12) + "#" * 12
        records.append(f"@read{i}\n{body}\n+\n{qual[:len(body)]}\n")
    path.write_text("".join(records))
    return path


@pytest.fixture(scope="session")
def demo_bam() -> Path:
    """The tracked demo BAM, or a skip. Present in a checkout, absent in an sdist."""
    if not _CHECKOUT_BAM.is_file():
        pytest.skip("no BAM fixture (expected: running from an sdist, not a checkout)")
    return _CHECKOUT_BAM


# ---------------------------------------------------------------------------
# pto-peaks
# ---------------------------------------------------------------------------
#
# `fragments` above is deliberately UNIFORM -- it exists to give
# `genomic_toolkit` a realistic size distribution, and a peak caller finds
# nothing in it, correctly. Asserting "no peaks" would pass just as well
# against a caller that never emits anything, so the peak fixtures below plant
# enrichment whose coordinates the tests then check.


@pytest.fixture(scope="session")
def chrom_sizes(tmp_path_factory) -> Path:
    """The two-column contig table BED input requires; a BAM carries its own."""
    path = tmp_path_factory.mktemp("pto-fixtures") / "test.chrom.sizes"
    path.write_text("chr1\t600000\nchr2\t120000\n")
    return path


@pytest.fixture(scope="session")
def enriched_fragments(tmp_path_factory) -> Path:
    """Coordinate-sorted fragments with two plantings of known position.

    Background is uniform over chr1 at roughly 2x; the two enriched windows sit
    at 250,000-250,600 and 410,000-410,400 and are deep enough that neither
    depends on the exact background model to clear an FDR of 0.05. The point of
    fixing the coordinates is that `test_call_peaks_finds_the_planted_peaks`
    can assert WHERE, which "some peaks came back" cannot.
    """
    path = tmp_path_factory.mktemp("pto-fixtures") / "enriched.bed"
    rng = random.Random(20260910)
    rows: list[tuple[int, int]] = []
    for _ in range(6000):
        start = rng.randrange(1_000, 500_000)
        rows.append((start, start + rng.choice([120, 150, 180, 220])))
    for centre, n in ((250_000, 1200), (410_000, 900)):
        for _ in range(n):
            start = rng.randrange(centre, centre + 600)
            rows.append((start, start + rng.choice([120, 150, 180])))
    rows.sort()
    path.write_text("".join(f"chr1\t{s}\t{e}\n" for s, e in rows))
    return path


@pytest.fixture(scope="session")
def regrouped_fragments(tmp_path_factory) -> Path:
    """A contig that reappears after another one -- `pto-peaks`'s exit-2 case.

    Sorted within each block, so this is not "unsorted input" in the ordinary
    sense; it is the specific shape a naive `cat` of two per-contig files
    produces, and the caller refuses it rather than emitting peaks for the
    first chr1 block and silently dropping the second.
    """
    path = tmp_path_factory.mktemp("pto-fixtures") / "regrouped.bed"
    path.write_text(
        "chr1\t1000\t1150\n"
        "chr1\t2000\t2150\n"
        "chr2\t1000\t1150\n"
        "chr1\t3000\t3150\n"
    )
    return path
