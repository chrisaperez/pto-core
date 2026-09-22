"""Where the pto-core binaries are, and what to say when they are not there.

Three places are searched, in this order:

  1. `$PTO_BIN_DIR`, if set. The escape hatch: a site that builds pto-core
     itself with its own htslib, or a developer pointing the wheel's Python
     API at a freshly-built tree, needs one variable and no reinstall.
  2. The wheel's vendored directory (`pto/_vendor/bin/`), populated by
     scikit-build-core from the modules' own `install(TARGETS ...)` rules --
     hence the `bin/` tail: those rules say `RUNTIME DESTINATION bin` and this
     package does not rewrite them. This is the path a `pip install` takes.
  3. `$PATH`. An sdist install, a conda environment, a module-file on an HPC
     cluster, or a source checkout whose `bin` is already exported.

Deliberately NOT searched: a hardcoded `../build/modules/...` relative to
this file. It works exactly once -- in the author's checkout -- and then
silently runs a stale binary from a build tree that has since moved on.
`PTO_BIN_DIR` covers the same need explicitly.
"""

from __future__ import annotations

import os
import shutil
from pathlib import Path

from .errors import BinaryNotFound

__all__ = ["TOOLS", "bin_dir", "resolve", "vendored_dir", "which"]

# The four CLIs this distribution wraps. `scrna_matrix` is absent on purpose:
# it is a Python extension module with its own wheel (`pip install
# scrna-matrix`), not a binary, and it is the one module whose build hard-fails
# without OpenMP -- bundling it here would make an unrelated wheel unbuildable
# on any macOS host without libomp.
#
# This tuple is the definition of "wrapped", and `pto-core/pyproject.toml`'s
# `PTO_BUILD_*` lines are the definition of "vendored". They have to agree in
# both directions and neither one enforces the other, so
# `packaging/tests/test_api.py` compares them against the binaries actually
# present. `pto-peaks` is here because it went in the wheel first and was
# reachable by nothing for a while; that is the failure the test now watches
# for.
TOOLS: tuple[str, ...] = (
    "fastq_stream", "cuttag_profiler", "genomic_toolkit", "pto-peaks",
)

_ENV_OVERRIDE = "PTO_BIN_DIR"

# Where each tool lands in a CMake build tree, for the "build it yourself"
# hint in `resolve()`. Every module but one puts its binary in a directory
# named after the binary, so the mapping is identity and the dict holds the
# exception: the target is `pto-peaks` and the module directory is `peaks`.
# Printing `build/modules/pto-peaks` would send someone to a path that does
# not exist, in the message whose entire job is to be actionable.
_BUILD_SUBDIR: dict[str, str] = {"pto-peaks": "peaks"}


def vendored_dir() -> Path:
    """The in-wheel binary directory. May not exist (sdist/source installs)."""
    return Path(__file__).resolve().parent / "_vendor" / "bin"


def _search_path() -> list[Path]:
    dirs: list[Path] = []
    override = os.environ.get(_ENV_OVERRIDE)
    if override:
        dirs.append(Path(override).expanduser())
    dirs.append(vendored_dir())
    return dirs


def which(tool: str) -> Path | None:
    """Resolve one tool, or None. Never raises."""
    if tool not in TOOLS:
        raise ValueError(f"unknown pto-core tool {tool!r}; expected one of {TOOLS}")

    for directory in _search_path():
        candidate = directory / tool
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return candidate

    found = shutil.which(tool)
    return Path(found) if found else None


def resolve(tool: str) -> Path:
    """Resolve one tool or raise `BinaryNotFound` saying how to fix it."""
    found = which(tool)
    if found is not None:
        return found

    searched = [str(d) for d in _search_path()] + ["$PATH"]
    raise BinaryNotFound(
        f"pto-core binary {tool!r} not found.\n"
        f"  searched: {', '.join(searched)}\n"
        f"\n"
        f"  A wheel installed from PyPI vendors it; a source or sdist install\n"
        f"  does not build it. Either install a wheel for this platform, or\n"
        f"  build pto-core and point pto at it:\n"
        f"\n"
        f"    cmake -S pto-core -B build -DCMAKE_BUILD_TYPE=Release\n"
        f"    cmake --build build -j\n"
        f"    export {_ENV_OVERRIDE}=$PWD/build/modules/{_BUILD_SUBDIR.get(tool, tool)}\n"
    )


def bin_dir() -> Path:
    """A directory to prepend to `$PATH` so the bare tool names resolve.

    This is what a Nextflow process or a shell pipeline wants: the binaries
    under their own names, with no Python interpreter in the call path. It
    returns the directory holding an actually-resolvable tool rather than
    asserting the vendored one exists, so it keeps working under
    `PTO_BIN_DIR` and under a PATH install.
    """
    return resolve(TOOLS[0]).parent
