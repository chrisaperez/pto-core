"""`pto-core man` / `pto-core glossary` -- the searchable command dictionary.

The parameter lists in this file's output are read off the LIVE
`argparse` parser (`_visible_parsers`, `_describe`), not retyped here --
this package has already been bitten twice by a second hand-copied set of
flags/bounds drifting from the parser that owns them (`api.py`'s own header
names both instances). A `--flag` `build_parser()` knows about and this file
does not would be that defect a third time.

What genuinely cannot come from argparse -- a one-line plain-English summary
and which `examples/` workflow demonstrates the command -- is declared once
below, per command. `packaging/tests/test_manual.py` checks both dicts
against the live parser (so a new subcommand without an entry here fails the
suite instead of shipping a silent gap) and checks every `--case NAME` token
against `benchmarks/compare_gold_standard.py`'s actual case registry (so a
renamed case is caught here too).
"""

from __future__ import annotations

import argparse
import sys
from typing import Iterable

__all__ = ["run"]

# Command name -> one-line summary. `man` with no argument prints these.
SUMMARIES: dict[str, str] = {
    "doctor": "Report which engines were found, their versions, and where they came from.",
    "bin-dir": "Print the directory holding the resolved binaries, for putting on $PATH.",
    "versions": "Print pto-core's and every engine's version, as JSON.",
    "exec": "Hand the process over to one wrapped binary directly, with its own native flags.",
    "man": "This command: list every pto-core subcommand, or describe one in detail.",
    "qc": "Streaming FASTQ QC -- read/base counts, quality, adapter content. Writes no reads.",
    "trim": "Adapter and quality trimming, writing cleaned FASTQ (single-end only).",
    "sizes": "Fragment-size distribution from a BAM/CRAM/BEDPE, in one streaming pass.",
    "markdup": "Duplicate rate and estimated library size, streaming -- no sort required.",
    "frip": "Fraction of reads in peaks, in the same streaming pass as size/duplicate QC.",
    "peaks": "Call peaks from a coordinate-sorted fragment BED or BAM (pto-peaks).",
    "profile": "CUT&Tag / CUT&RUN reference-point signal matrix from an indexed BAM.",
}

# Command name -> where to see it run, printed as the glossary's "See it run"
# line. Kept short; `examples/README.md` carries the full walkthrough.
EXAMPLES: dict[str, str] = {
    "qc": "examples/run_benchmarks.sh --case fastq/qc       (fastq_stream vs fastp)",
    "trim": "examples/README.md, \"FASTQ QC & trimming\" -- qc and trim share fastq_stream",
    "sizes": "examples/run_benchmarks.sh --case gtk/sizes     (genomic_toolkit vs samtools)",
    "markdup": "examples/README.md, \"Duplicate marking\"",
    "frip": "examples/run_benchmarks.sh --case gtk/frip       (genomic_toolkit vs samtools+bedtools)",
    "peaks": "examples/README.md, \"Peak calling\" -- see also pto-core/modules/peaks/CONCORDANCE.md",
    "profile": "examples/run_benchmarks.sh --case cuttag/profile (cuttag_profiler vs deepTools)",
}

# Commands with no example line get this instead of a KeyError.
_NO_EXAMPLE = "(utility command; not part of the benchmark comparisons)"


def visible_parsers(
    subparsers_action: "argparse._SubParsersAction",
) -> dict[str, argparse.ArgumentParser]:
    """Name -> parser, in registration order, aliases collapsed onto the first.

    `_SubParsersAction.choices` maps every name AND alias to the same parser
    object (`man`/`glossary` both point at one parser), so this keeps only
    the first name argparse saw for each -- which is the canonical name
    because `add_parser(name, aliases=[...])` always registers `name` first.
    """
    seen: set[int] = set()
    ordered: dict[str, argparse.ArgumentParser] = {}
    for name, parser in subparsers_action.choices.items():
        if id(parser) in seen:
            continue
        seen.add(id(parser))
        ordered[name] = parser
    return ordered


class Param:
    __slots__ = ("flag", "help", "required")

    def __init__(self, flag: str, help: str, required: bool):
        self.flag = flag
        self.help = help
        self.required = required


def describe_params(parser: argparse.ArgumentParser) -> list[Param]:
    """Every user-facing argument of `parser`, positionals first, in order.

    `parser._actions` is argparse's only introspection surface for "what
    arguments does this parser take" -- there is no public equivalent in the
    standard library, which is why this is the one place in the package that
    reaches past a leading underscore.
    """
    params: list[Param] = []
    for action in parser._actions:  # noqa: SLF001 -- see docstring
        if isinstance(action, argparse._HelpAction):
            continue
        help_text = action.help or "(no description)"
        if not action.option_strings:
            # Positional. REMAINDER/'*'/'?' arguments are the one case where a
            # positional isn't mandatory (`exec`'s trailing argv).
            optional_positional = action.nargs in ("?", "*", argparse.REMAINDER)
            params.append(Param(action.dest, help_text, required=not optional_positional))
            continue
        flag = max(action.option_strings, key=len)
        default = action.default
        if default not in (None, False, argparse.SUPPRESS) and not action.required:
            help_text = f"{help_text} [default: {default}]"
        params.append(Param(flag, help_text, required=bool(action.required)))
    return params


def _print_entry(name: str, parser: argparse.ArgumentParser, *, out) -> None:
    summary = SUMMARIES.get(name, parser.description or parser.format_usage().strip())
    print(f"{name} -- {summary}", file=out)
    print(file=out)

    params = describe_params(parser)
    required = [p for p in params if p.required]
    optional = [p for p in params if not p.required]

    print("Required parameters:", file=out)
    if required:
        for p in required:
            print(f"  {p.flag:<18} {p.help}", file=out)
    else:
        print("  (none)", file=out)
    print(file=out)

    print("Optional / beneficial parameters:", file=out)
    if optional:
        for p in optional:
            print(f"  {p.flag:<18} {p.help}", file=out)
    else:
        print("  (none)", file=out)
    print(file=out)

    print(f"See it run:  {EXAMPLES.get(name, _NO_EXAMPLE)}", file=out)


def _print_index(parsers: Iterable[tuple[str, argparse.ArgumentParser]], *, out) -> None:
    print("pto-core commands  (`pto-core man <command>` for parameters and an example)", file=out)
    print(file=out)
    for name, _parser in parsers:
        print(f"  {name:<10} {SUMMARIES.get(name, '')}", file=out)


def run(subparsers_action: "argparse._SubParsersAction", topic: str | None, *, out=None) -> int:
    out = out if out is not None else sys.stdout
    parsers = visible_parsers(subparsers_action)

    if topic is None:
        _print_index(parsers.items(), out=out)
        return 0

    if topic not in parsers:
        print(f"pto-core man: no such command {topic!r}", file=sys.stderr)
        print(f"available: {', '.join(parsers)}", file=sys.stderr)
        return 1

    _print_entry(topic, parsers[topic], out=out)
    return 0
