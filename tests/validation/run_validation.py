#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Public-data end-to-end validation for pto-core.

    run_validation.py list                      # what exists, and what can run here
    run_validation.py fetch --allow-download    # stage the inputs
    run_validation.py run                       # compare against the gold standards
    run_validation.py run --case cuttag/profile-vs-deeptools

Exits 0 when nothing failed. Skips do not fail the run -- a host without
deeptools installed, or an air-gapped host with nothing staged, is an expected
state, and treating it as a failure trains everyone to ignore the exit code.
The summary reports the skip count prominently instead, so a run that quietly
validated nothing is visible rather than green.

See README.md in this directory for what each case actually compares and why
the thresholds sit where they do.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

# Allow `python tests/validation/run_validation.py` from anywhere without an
# install step. The suite has to run on a cluster login node where `pip
# install -e .` is not a thing anyone is going to do.
sys.path.insert(0, str(Path(__file__).resolve().parent))

from pto_validation import fetch  # noqa: E402
from pto_validation.config import (  # noqa: E402
    DATASETS_FILE,
    DEFAULT_BUILD_DIR,
    DEFAULT_DATA_DIR,
    find_pto_binary,
    find_reference_tools,
    load_datasets,
)
from pto_validation.harness import (  # noqa: E402
    Context,
    load_cases,
    run_case,
)
from pto_validation.report import Report  # noqa: E402


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="run_validation.py",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "command",
        choices=["list", "fetch", "run", "verify"],
        help="list: show datasets/cases and what can run here. "
             "fetch: stage inputs. run: execute cases. "
             "verify: re-hash cached inputs against their pins.",
    )
    parser.add_argument(
        "--datasets", nargs="+", metavar="KEY",
        help="restrict to these dataset keys (default: all except the "
             "deep-coverage tier, which is opt-in)",
    )
    parser.add_argument(
        "--case", dest="cases", nargs="+", metavar="NAME",
        help="restrict to these case names",
    )
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--build-dir", type=Path, default=DEFAULT_BUILD_DIR)
    parser.add_argument(
        "--work-dir", type=Path, default=None,
        help="scratch and artefact directory (default: <data-dir>/work)",
    )
    parser.add_argument(
        "--threads", type=int, default=os.cpu_count() or 4,
        help="threads handed to both pto and the reference tools, so the "
             "comparison is not confounded by one side being given more cores",
    )
    parser.add_argument(
        "--allow-download", action="store_true",
        help="permit network access. Off by default: these tools are built for "
             "hosts that hold the sequencing data, and a harness that silently "
             "reaches for the internet is the wrong default there.",
    )
    parser.add_argument(
        "--record-checksums", action="store_true",
        help="fetch only: pin missing sha256 digests in datasets.json from "
             "already-cached files, for review and commit",
    )
    parser.add_argument("--report-dir", type=Path, default=None)
    parser.add_argument(
        "--all", action="store_true",
        help="include the deep-coverage datasets that are excluded by default",
    )

    args = parser.parse_args(argv)

    datasets = load_datasets()
    if args.datasets:
        unknown = set(args.datasets) - set(datasets)
        if unknown:
            parser.error(f"unknown dataset(s): {', '.join(sorted(unknown))}")
        selected = {k: v for k, v in datasets.items() if k in args.datasets}
    elif args.all:
        selected = datasets
    else:
        # The deep tier is ~1.5 GB of download and tens of minutes of compute
        # per case. Opting into it should be a decision, not something that
        # happens because someone typed `run`.
        selected = {k: v for k, v in datasets.items() if not k.endswith("-deep")}

    registry = load_cases()
    if args.cases:
        unknown = set(args.cases) - set(registry)
        if unknown:
            parser.error(f"unknown case(s): {', '.join(sorted(unknown))}")

    if args.command == "list":
        return _list(datasets, registry, args)
    if args.command == "fetch":
        return _fetch(selected, args)
    if args.command == "verify":
        return _verify(selected, args)
    return _run(selected, registry, args)


def _list(datasets, registry, args) -> int:
    reference = find_reference_tools()

    print("DATASETS")
    for key, dataset in datasets.items():
        staged = 0
        for role in dataset.files:
            spec = dataset.files[role]
            if spec.repo_path:
                staged += 1
            elif (dataset.local_dir(args.data_dir) / spec.filename).exists():
                staged += 1
        total_bytes = sum(f.approx_bytes or 0 for f in dataset.files.values())
        print(f"  {key}")
        print(f"      {dataset.title}")
        print(f"      {dataset.assay} · {dataset.source}")
        print(
            f"      files staged: {staged}/{len(dataset.files)}"
            + (f" · ~{total_bytes / 1e9:.2f} GB total" if total_bytes else "")
        )

    print("\nCASES")
    for name, spec in sorted(registry.items()):
        missing = [t for t in spec.reference if not reference.get(t)]
        binary = find_pto_binary(spec.module, args.build_dir)
        blockers = []
        if missing:
            blockers.append("needs " + ", ".join(missing))
        if spec.module != "scrna_matrix" and binary is None:
            blockers.append(f"{spec.module} not built")
        state = "ready" if not blockers else "SKIP: " + "; ".join(blockers)
        print(f"  {name:<38} {state}")
        print(f"      {spec.summary}")

    print("\nREFERENCE TOOLS")
    for name, path in sorted(reference.items()):
        print(f"  {name:<16} {path or '-- not found --'}")
    return 0


def _fetch(datasets, args) -> int:
    if args.record_checksums:
        fetch.record_checksums(datasets, args.data_dir, DATASETS_FILE)
        return 0

    failures = 0
    for key, dataset in datasets.items():
        print(f"\n{key}: {dataset.title}")
        for role in dataset.files:
            try:
                path = fetch.ensure_file(
                    dataset, role, args.data_dir, allow_download=args.allow_download
                )
                print(f"  {role:<12} {path}")
            except fetch.DownloadsDisabled as exc:
                print(f"  {role:<12} not staged ({exc.args[0].splitlines()[0]})")
                failures += 1
            except fetch.FetchError as exc:
                print(f"  {role:<12} FAILED: {exc}", file=sys.stderr)
                failures += 1

    if failures and not args.allow_download:
        print(
            f"\n{failures} input(s) not staged. Re-run with --allow-download to "
            "fetch them.",
            file=sys.stderr,
        )
        # Not an error: refusing to download without authorisation is this
        # command working correctly, and CI on an air-gapped runner should not
        # go red for it.
        return 0
    return 1 if failures else 0


def _verify(datasets, args) -> int:
    bad = 0
    for key, dataset in datasets.items():
        for role, spec in dataset.files.items():
            path = dataset.local_dir(args.data_dir) / spec.filename
            if spec.repo_path or not path.exists():
                continue
            algorithm = "sha256" if spec.sha256 else "md5" if spec.md5 else None
            if algorithm is None:
                print(f"  {key}/{role}: no digest pinned, cannot verify")
                continue
            expected = spec.sha256 or spec.md5 or ""
            actual = fetch.digest(path, algorithm)
            ok = actual.lower() == expected.lower()
            print(f"  {'OK  ' if ok else 'BAD '} {key}/{role}  ({algorithm})")
            if not ok:
                print(f"        expected {expected}\n        actual   {actual}")
                bad += 1
    return 1 if bad else 0


def _run(datasets, registry, args) -> int:
    work_dir = args.work_dir or (args.data_dir / "work")
    report_dir = args.report_dir or (work_dir / "reports")
    reference = find_reference_tools()
    report = Report(output_dir=report_dir)

    for key, dataset in datasets.items():
        names = [n for n in dataset.cases if n in registry]
        if args.cases:
            names = [n for n in names if n in args.cases]
        if not names:
            continue

        print(f"\n=== {key}: {dataset.title} ===")
        for name in names:
            print(f"  -> {name}")
            context = Context(
                dataset=dataset,
                data_dir=args.data_dir,
                build_dir=args.build_dir,
                work_dir=work_dir / key,
                threads=args.threads,
                allow_download=args.allow_download,
                reference_tools=reference,
            )
            result = run_case(registry[name], context, report)
            print(f"     {result.status.value} ({result.wall_seconds:.1f}s)")

    if not report.results:
        print(
            "\nNo cases ran. Check `run_validation.py list` -- most likely the "
            "selected datasets declare no cases that are registered, or the "
            "--case filter matched nothing.",
            file=sys.stderr,
        )
        return 1

    report.print_summary()
    report_dir.mkdir(parents=True, exist_ok=True)
    report.write_json(report_dir / "validation.json")
    report.write_markdown(report_dir / "validation.md")
    print(f"\nReports: {report_dir / 'validation.md'}")
    return 0 if report.ok else 1


if __name__ == "__main__":
    sys.exit(main())
