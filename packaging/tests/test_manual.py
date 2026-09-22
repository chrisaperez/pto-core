"""`pto-core man` / `pto-core glossary`.

Two kinds of check, deliberately kept apart:

1. Structural, against the LIVE parser (`build_parser()`) -- so a subcommand
   added to `cli.py` without a `_manual.SUMMARIES` entry fails here instead
   of shipping a glossary with a silent gap, and a name removed from `cli.py`
   fails here instead of leaving a stale entry that `man` can never reach.
2. Cross-reference, against `examples/` and `benchmarks/compare_gold_standard.py`
   -- so a `--case gtk/frip` token in `_manual.EXAMPLES` that no longer names a
   real case is caught here, the same way `test_every_flag_the_api_emits_exists_
   in_the_binarys_own_help` catches a stale flag in `test_api.py`. These skip
   rather than fail when the sibling directories are not on disk at all (a
   detached wheel-only checkout has no `examples/`), because that is a missing
   fixture, not a defect in this package.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import pytest

from pto import _manual
from pto.cli import build_parser

# Utility commands genuinely have nothing to demonstrate in a benchmark
# comparison -- `doctor`/`bin-dir`/`versions` report on the install, `exec` is
# an escape hatch, and `man` is this command. Every command that runs an
# analysis must have a real cross-reference; this is the list of the ones
# that are exempt, not a place to add a command that was merely forgotten.
_NO_EXAMPLE_NEEDED = {"doctor", "bin-dir", "versions", "exec", "man"}

_REPO_ROOT = Path(__file__).resolve().parents[3]


def _canonical_names() -> set[str]:
    parser = build_parser()
    subparsers_action = next(
        action for action in parser._actions  # noqa: SLF001
        if hasattr(action, "choices") and isinstance(action.choices, dict)
    )
    return set(_manual.visible_parsers(subparsers_action))


def test_every_live_command_has_a_summary():
    live = _canonical_names()
    documented = set(_manual.SUMMARIES)
    missing = live - documented
    stale = documented - live
    assert not missing, f"cli.py commands with no _manual.SUMMARIES entry: {missing}"
    assert not stale, f"_manual.SUMMARIES entries for commands cli.py no longer has: {stale}"


def test_every_analysis_command_cross_references_the_examples_playground():
    live = _canonical_names()
    needs_example = live - _NO_EXAMPLE_NEEDED
    missing = needs_example - set(_manual.EXAMPLES)
    assert not missing, (
        f"commands with no examples/ cross-reference in _manual.EXAMPLES: {missing}"
    )


def test_example_case_tokens_name_a_real_benchmark_case():
    benchmarks_dir = _REPO_ROOT / "benchmarks"
    if not (benchmarks_dir / "compare_gold_standard.py").is_file():
        pytest.skip("benchmarks/ is not on disk (not a full checkout)")

    sys.path.insert(0, str(benchmarks_dir))
    sys.path.insert(0, str(_REPO_ROOT / "scripts"))
    import compare_gold_standard  # noqa: PLC0415

    cases = set(compare_gold_standard.REGISTRY)
    for name, line in _manual.EXAMPLES.items():
        match = re.search(r"--case ([a-z_]+/[a-z_]+)", line)
        if match:
            assert match.group(1) in cases, (
                f"_manual.EXAMPLES[{name!r}] points at --case {match.group(1)!r}, "
                f"which is not in compare_gold_standard.REGISTRY ({sorted(cases)})"
            )


def test_example_paths_exist_on_disk():
    if not (_REPO_ROOT / "examples").is_dir():
        pytest.skip("examples/ is not on disk (not a full checkout)")
    for name, line in _manual.EXAMPLES.items():
        match = re.search(r"(examples/[\w./-]+)", line)
        assert match, f"_manual.EXAMPLES[{name!r}] does not name a file under examples/: {line!r}"
        target = _REPO_ROOT / match.group(1)
        assert target.is_file(), f"_manual.EXAMPLES[{name!r}] points at {target}, which does not exist"


# ---------------------------------------------------------------------------
# Behaviour, through the real CLI
# ---------------------------------------------------------------------------


def test_man_with_no_argument_lists_every_command(capsys):
    from pto.cli import main

    assert main(["man"]) == 0
    out = capsys.readouterr().out
    for name in _canonical_names():
        assert name in out


def test_man_on_one_command_shows_required_and_optional_sections(capsys):
    from pto.cli import main

    assert main(["man", "frip"]) == 0
    out = capsys.readouterr().out
    assert "Required parameters:" in out
    assert "Optional / beneficial parameters:" in out
    assert "alignment" in out and "peaks" in out
    assert "See it run:" in out


def test_glossary_is_an_alias_for_man(capsys):
    from pto.cli import main

    assert main(["glossary", "qc"]) == 0
    out = capsys.readouterr().out
    assert "qc --" in out


def test_an_unknown_topic_fails_and_lists_what_exists(capsys):
    from pto.cli import main

    assert main(["man", "not-a-real-command"]) == 1
    err = capsys.readouterr().err
    assert "not-a-real-command" in err
    assert "available:" in err
