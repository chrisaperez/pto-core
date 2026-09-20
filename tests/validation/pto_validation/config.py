# SPDX-License-Identifier: MIT
"""Paths, dataset registry, and the assertion thresholds.

Thresholds live here rather than next to each check so the whole tolerance
budget of the suite is readable in one screen. Every number is a claim about
what these tools promise, so each one carries the reason it is set where it is.
Loosening one is a real decision -- record why in the comment when you do.
"""

from __future__ import annotations

import json
import os
import shutil
from dataclasses import dataclass, field
from pathlib import Path

# tests/validation/pto_validation/config.py -> repo root
REPO_ROOT = Path(__file__).resolve().parents[3]
VALIDATION_ROOT = REPO_ROOT / "tests" / "validation"
DATASETS_FILE = VALIDATION_ROOT / "datasets.json"

#: Downloaded inputs are cached outside the repo and shared across runs and
#: branches. These files are multi-gigabyte; putting them under the worktree
#: makes every `git status` slow and invites an accidental commit.
DEFAULT_DATA_DIR = Path(
    os.environ.get("PTO_VALIDATION_DATA", Path.home() / ".cache" / "pto-validation")
)

DEFAULT_BUILD_DIR = Path(os.environ.get("PTO_BUILD_DIR", REPO_ROOT / "build"))


# ---------------------------------------------------------------------------
# Tolerances
# ---------------------------------------------------------------------------

#: Signal-track agreement against deeptools computeMatrix.
#:
#: Not 1.0, and it never will be: deeptools bins by fragment midpoint overlap
#: while cuttag_profiler's `reads` mode counts a fragment into every bin its
#: span touches, so the two differ at bin boundaries by construction. 0.99 is
#: what that boundary effect costs on a 50 bp bin; a drop below it means
#: something other than binning changed.
SIGNAL_PEARSON_MIN = 0.99
SIGNAL_SPEARMAN_MIN = 0.99

#: Base-pair overlap agreement against the bedtools peak intersection.
PEAK_JACCARD_MIN = 0.95

#: Fragment size distribution parity against the samtools TLEN histogram.
#: KS D of 0.02 means the two ECDFs never diverge by more than 2 percentage
#: points anywhere in the support.
FRAGMENT_KS_MAX = 0.02

#: FRiP and duplicate rate are counts, not estimates. genomic_toolkit's
#: validate_against_bedtools.sh already asserts *exact* equality of the
#: in-peak fragment count, and this harness holds the derived ratio to the
#: rounding of the JSON report rather than to a statistical tolerance.
FRIP_RELATIVE_MAX = 1e-6

#: markdup vs Picard MarkDuplicates. This one is deliberately loose and is NOT
#: a defect budget: README.md "Trade-offs" records that genomic_toolkit keys
#: duplicates on the aligned span while Picard uses the unclipped 5' end, and
#: has no optical/UMI awareness. The two therefore disagree on soft-clipped
#: reads by design. The assertion exists to catch a *change* in that gap, so
#: the report prints the observed value on every run, pass or fail.
DUPLICATE_RATE_RELATIVE_MAX = 0.05

#: QC counters against seqkit/fastp. Read and base counts are exact integers
#: coming off the same input file; any difference at all is a bug.
QC_COUNT_RELATIVE_MAX = 0.0

#: Per-cycle mean quality against fastp. Both compute the same mean over the
#: same bases; the tolerance covers only fastp reporting to one decimal.
QC_QUALITY_ABS_MAX = 0.05

#: k-NN graph agreement between the exact and HNSW paths in scrna_matrix.
#: HNSW is approximate on purpose, and 0.90 recall at k=15 is the operating
#: point the module's own docs claim. This is a self-consistency check, not a
#: comparison against an external tool.
KNN_RECALL_MIN = 0.90


@dataclass(frozen=True)
class RemoteFile:
    """One downloadable input.

    `sha256` is pinned in datasets.json and verified after every download. It
    is the only thing standing between this harness and a silently-changed
    upstream file, which would show up as an unexplained metric regression
    weeks later. ENCODE entries may instead carry `encode_accession`, in which
    case the URL and checksum are resolved from the portal's own metadata API
    at fetch time -- see fetch.py for why that indirection exists.
    """

    name: str
    url: str | None = None
    sha256: str | None = None
    md5: str | None = None
    encode_accession: str | None = None
    #: Path, relative to the repo root, of a fixture already checked into the
    #: tree. Takes precedence over every remote field and never touches the
    #: network -- this is how the offline-by-default datasets are expressed.
    repo_path: str | None = None
    approx_bytes: int | None = None
    note: str = ""

    @property
    def filename(self) -> str:
        return self.name


@dataclass(frozen=True)
class Dataset:
    """A named public run plus the files it needs and the cases it feeds."""

    key: str
    title: str
    assay: str
    source: str
    description: str
    files: dict[str, RemoteFile]
    cases: list[str] = field(default_factory=list)
    params: dict = field(default_factory=dict)

    def local_dir(self, data_dir: Path) -> Path:
        return data_dir / self.key

    def path(self, data_dir: Path, role: str) -> Path:
        return self.local_dir(data_dir) / self.files[role].filename


def load_datasets(path: Path | None = None) -> dict[str, Dataset]:
    """Parse datasets.json into Dataset objects, failing loudly on typos.

    A misspelled role key is the failure mode worth catching early: it would
    otherwise surface as a KeyError deep inside a case, after a multi-gigabyte
    download has already completed.
    """
    path = path or DATASETS_FILE
    raw = json.loads(path.read_text())

    datasets: dict[str, Dataset] = {}
    for key, entry in raw["datasets"].items():
        files = {
            role: RemoteFile(**spec) for role, spec in entry.get("files", {}).items()
        }
        for role, spec in files.items():
            if not (spec.url or spec.encode_accession or spec.repo_path):
                raise ValueError(
                    f"dataset {key!r} file {role!r} has none of `url`, "
                    "`encode_accession` or `repo_path`; one is required"
                )
        datasets[key] = Dataset(
            key=key,
            title=entry["title"],
            assay=entry["assay"],
            source=entry["source"],
            description=entry.get("description", ""),
            files=files,
            cases=entry.get("cases", []),
            params=entry.get("params", {}),
        )
    return datasets


# ---------------------------------------------------------------------------
# Tool discovery
# ---------------------------------------------------------------------------

#: Reference tools, by the capability they provide. A missing tool skips the
#: cases that need it -- it never fails them. The whole point of this repo is
#: that the tools run where the data is, and the machine holding the data is
#: exactly the machine least likely to have a conda environment full of
#: deeptools on it.
REFERENCE_TOOLS = {
    "samtools": "samtools",
    "bedtools": "bedtools",
    "bowtie2": "bowtie2",
    "picard": "picard",
    "computeMatrix": "deeptools computeMatrix",
    "bamCoverage": "deeptools bamCoverage",
    "seqkit": "seqkit",
    "fastp": "fastp",
}


def find_reference_tools() -> dict[str, str | None]:
    """Resolve each reference tool to an absolute path, or None if absent."""
    return {name: shutil.which(name) for name in REFERENCE_TOOLS}


#: pto binaries, by module. Resolved against the build directory; the
#: harness reports a clear "build it first" error rather than a FileNotFound.
PTO_BINARIES = {
    "fastq_stream": "modules/fastq_stream/fastq_stream",
    "cuttag_profiler": "modules/cuttag_profiler/cuttag_profiler",
    "genomic_toolkit": "modules/genomic_toolkit/genomic_toolkit",
}


def find_pto_binary(module: str, build_dir: Path) -> Path | None:
    """Locate a module's CLI in an aggregate *or* standalone build tree.

    Both layouts are supported because both are documented: the top-level build
    puts binaries under build/modules/<module>/, while a standalone
    `cmake -S modules/<module> -B build/<module>` puts them at the root of that
    tree. A harness that only understood the first would quietly skip every
    module for anyone following the standalone instructions.
    """
    # scrna_matrix has no CLI at all -- it is a header-only engine plus
    # pybind11 bindings -- so "no binary" is its normal state rather than an
    # unbuilt one. Returning None lets callers treat both uniformly.
    if module not in PTO_BINARIES:
        return None

    candidates = [
        build_dir / PTO_BINARIES[module],
        build_dir / module / module,
        build_dir / module,
    ]
    for candidate in candidates:
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return candidate
    return None
