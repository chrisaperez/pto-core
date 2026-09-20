"""`pto-core` -- the console script that ships with this wheel.

NAMED `pto-core`, NOT `pto`, and that is not a style choice: `pto-cli`'s
`[project.scripts]` already binds `pto` to `pto_cli.main:main`, the client for
the managed service. Two distributions installing the same script name means
whichever was installed second wins, silently, and a lab that pip-installed
both would find `pto run` either submitting a cloud job or refusing to. The
import name `pto` is free and is what the Python API uses; the script name is
not.

What this exists for is the first ten minutes of an evaluation:

    pip install pto-core
    pto-core doctor                 # did it find the engines, and where
    pto-core man                    # every command, one line each
    pto-core qc reads.fastq.gz      # one real answer, as JSON
    export PATH="$(pto-core bin-dir):$PATH"   # bare tool names, for pipelines

`exec` and `bin-dir` are the two that matter for a pipeline: `bin-dir` prints
the directory holding the vendored binaries so Nextflow, Snakemake or a shell
script can put them on `PATH` with no Python in the call path at all. `man`
(alias `glossary`) reads its parameter lists off this file's own parser --
see `_manual.py` -- so it cannot drift the way a hand-written cheat sheet does.
"""

from __future__ import annotations

import argparse
import json
import os
import sys

from . import __version__, _manual, api, tool_versions
from ._binaries import TOOLS, bin_dir, resolve, which
from ._cpu import available_cpus
from .errors import InvalidInput, PtoError, ToolError


def _cmd_doctor(args: argparse.Namespace) -> int:
    """Say what is installed and where it came from. Never raises."""
    print(f"pto-core {__version__}  (python {sys.version.split()[0]})")
    print(f"threads available to this process: {available_cpus()}")

    override = os.environ.get("PTO_BIN_DIR")
    print(f"PTO_BIN_DIR: {override or '(unset)'}")

    versions = tool_versions()
    missing = []
    print("\nengines:")
    for tool in TOOLS:
        path = which(tool)
        if path is None:
            missing.append(tool)
            print(f"  {tool:<18} NOT FOUND")
            continue
        source = "vendored" if "_vendor" in path.parts else "PATH/override"
        print(f"  {tool:<18} {versions.get(tool) or '(version unknown)'}")
        print(f"  {'':<18} {path}  [{source}]")

    if missing:
        print(
            f"\n{len(missing)} engine(s) missing. A wheel built for this platform\n"
            f"vendors them; an sdist install does not build them. Either install\n"
            f"a wheel, or build pto-core and export PTO_BIN_DIR:\n"
            f"\n"
            f"  cmake -S pto-core -B build -DCMAKE_BUILD_TYPE=Release\n"
            f"  cmake --build build -j\n"
            f"  export PTO_BIN_DIR=$PWD/build/modules/<module>\n"
            f"\n"
            f"  ...one directory per engine; each binary lands in its own.\n"
        )
        return 1

    print("\nall engines resolved.")
    return 0


def _cmd_bin_dir(args: argparse.Namespace) -> int:
    print(bin_dir())
    return 0


def _cmd_exec(args: argparse.Namespace) -> int:
    """Hand the process over to a pto-core binary.

    `os.execv`, not a subprocess: the caller gets the binary's exit code, its
    signals and its stdio with no Python frame left in between, which is what a
    pipeline stage needs.
    """
    binary = resolve(args.tool)
    os.execv(str(binary), [str(binary)] + list(args.argv))
    raise AssertionError("unreachable")  # pragma: no cover


def _emit(result, as_json: bool) -> None:
    if as_json:
        print(json.dumps(result.as_dict(), indent=2, sort_keys=True))
        return
    for key, value in result.metrics.items():
        if isinstance(value, (dict, list)):
            continue
        print(f"{key:<28} {value}")
    cost = f"{result.wall_seconds:.2f}s"
    if result.peak_rss_mb is not None:
        cost += f", peak RSS {result.peak_rss_mb} MB"
    print(f"\n[{cost}]", file=sys.stderr)


def _cmd_qc(args: argparse.Namespace) -> int:
    result = api.run_qc(args.reads, out_dir=args.out_dir, threads=args.threads)
    if args.json:
        _emit(result, True)
    else:
        summary = result.metrics.get("summary", {})
        for key, value in summary.items():
            print(f"{key:<28} {value}")
        print(f"\n[{result.wall_seconds:.2f}s, peak RSS {result.peak_rss_mb} MB]",
              file=sys.stderr)
    return 0


def _cmd_sizes(args: argparse.Namespace) -> int:
    result = api.fragment_sizes(args.alignment, histogram=args.histogram)
    _emit(result, args.json)
    return 0


def _cmd_frip(args: argparse.Namespace) -> int:
    result = api.frip(args.alignment, args.peaks, genome_size=args.genome_size)
    _emit(result, args.json)
    return 0


def _cmd_peaks(args: argparse.Namespace) -> int:
    result = api.call_peaks(
        args.fragments, args.out,
        chrom_sizes=args.chrom_sizes,
        qvalue=args.qvalue,
        cutoff=args.cutoff,
        min_length=args.min_length,
        max_gap=args.max_gap,
        extend_peaks=args.extend_peaks,
        lambda_bg=args.lambda_bg,
        min_mapq=args.min_mapq,
    )
    called = result.artifacts["peaks"]
    # Counted from the file this call just wrote, not parsed out of the tool's
    # diagnostics and not a field this package made up: `pto-peaks` emits
    # narrowPeak and no JSON report, so `result.metrics` is legitimately empty.
    # A data line in narrowPeak is a peak, and the writer emits no header.
    n = sum(1 for line in called.read_text().splitlines() if line.strip())
    if args.json:
        payload = result.as_dict()
        payload["peaks_called"] = n
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        print(f"{'peaks_called':<28} {n}")
        print(f"{'narrowPeak':<28} {called}")
        cost = f"{result.wall_seconds:.2f}s"
        if result.peak_rss_mb is not None:
            cost += f", peak RSS {result.peak_rss_mb} MB"
        # stdout is block-buffered into a pipe while stderr is not, so without
        # this the cost line and the tool's diagnostics land above the numbers
        # they describe whenever this is redirected.
        sys.stdout.flush()
        print(f"\n[{cost}]", file=sys.stderr)
    # What the BAM reader skipped, and why, is only on stderr. Passing it
    # through matters most in the case that looks like success: a single-end
    # BAM calls nothing and exits 0, and the reason is the `not the leftmost
    # mate` count in this line.
    if result.stderr.strip():
        sys.stdout.flush()
        print(result.stderr.rstrip(), file=sys.stderr)
    return 0


def _cmd_versions(args: argparse.Namespace) -> int:
    print(json.dumps({"pto-core": __version__, "engines": tool_versions()}, indent=2))
    return 0


def _cmd_man(args: argparse.Namespace) -> int:
    return _manual.run(args.subparsers, args.topic)


def _cmd_trim(args: argparse.Namespace) -> int:
    result = api.trim(
        args.reads, args.out, out_dir=args.out_dir, threads=args.threads,
        min_len=args.min_len, min_mean_q=args.min_mean_q,
    )
    if args.json:
        _emit(result, True)
        return 0
    summary = result.metrics.get("summary", {})
    for key, value in summary.items():
        print(f"{key:<28} {value}")
    print(f"\ntrimmed reads written to    {args.out}")
    print(f"[{result.wall_seconds:.2f}s, peak RSS {result.peak_rss_mb} MB]", file=sys.stderr)
    return 0


def _cmd_markdup(args: argparse.Namespace) -> int:
    result = api.mark_duplicates(
        args.alignment, min_length=args.min_length, max_length=args.max_length,
        min_mapq=args.min_mapq, single_end=args.single_end, unsorted=args.unsorted,
    )
    _emit(result, args.json)
    return 0


def _cmd_profile(args: argparse.Namespace) -> int:
    result = api.profile(
        args.alignment, args.regions, args.out_dir,
        threads=args.threads, reference_point=args.reference_point,
        upstream=args.upstream, downstream=args.downstream, bin_size=args.bin_size,
    )
    if args.json:
        print(json.dumps(result.as_dict(), indent=2, sort_keys=True))
        return 0
    for role, path in result.artifacts.items():
        print(f"{role:<28} {path}")
    cost = f"{result.wall_seconds:.2f}s"
    if result.peak_rss_mb is not None:
        cost += f", peak RSS {result.peak_rss_mb} MB"
    print(f"\n[{cost}]", file=sys.stderr)
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="pto-core",
        description="Run pto-core's genomics engines. See `import pto` for the API.",
    )
    parser.add_argument("--version", action="version", version=f"pto-core {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    sub.add_parser("doctor", help="report which engines were found and where").set_defaults(
        func=_cmd_doctor)
    sub.add_parser(
        "bin-dir",
        help="print the directory to prepend to $PATH for the bare tool names",
    ).set_defaults(func=_cmd_bin_dir)
    sub.add_parser("versions", help="engine versions as JSON").set_defaults(
        func=_cmd_versions)

    p_man = sub.add_parser(
        "man", aliases=["glossary"],
        help="describe one command (parameters + example), or list them all")
    p_man.add_argument("topic", nargs="?", default=None,
                       help="a command name; omit to list every command")
    # `_cmd_man` needs the live subparsers action to read every OTHER
    # command's own parser, which is the point: the glossary can't drift
    # from `--help` because it is read off the same object `--help` is.
    p_man.set_defaults(func=_cmd_man, subparsers=sub)

    p_exec = sub.add_parser(
        "exec", help="exec one engine directly, passing the rest of the command line")
    p_exec.add_argument("tool", choices=TOOLS, help="which vendored binary to run")
    p_exec.add_argument("argv", nargs=argparse.REMAINDER,
                        help="the binary's own flags, passed through unchanged")
    p_exec.set_defaults(func=_cmd_exec)

    p_qc = sub.add_parser("qc", help="FASTQ QC (fastq_stream --qc-only)")
    p_qc.add_argument("reads", help="FASTQ (plain/gzip/BGZF), or BAM/CRAM/SAM")
    p_qc.add_argument("--out-dir", default=None, help="keep the QC artifacts here")
    p_qc.add_argument("-t", "--threads", type=int, default=None,
                      help="[default: hardware concurrency]")
    p_qc.add_argument("--json", action="store_true", help="full report as JSON")
    p_qc.set_defaults(func=_cmd_qc)

    p_trim = sub.add_parser(
        "trim", help="adapter/quality trimming, writing cleaned FASTQ (fastq_stream)")
    p_trim.add_argument("reads", help="FASTQ (plain/gzip/BGZF), or BAM/CRAM/SAM")
    p_trim.add_argument("out", help="trimmed FASTQ output path (\"-\" for stdout)")
    p_trim.add_argument("--out-dir", default=None, help="keep the QC report here too")
    p_trim.add_argument("-t", "--threads", type=int, default=None,
                        help="[default: hardware concurrency]")
    p_trim.add_argument("--min-len", type=int, default=None,
                        help="discard reads shorter than N after trimming [default: 15]")
    p_trim.add_argument("--min-mean-q", type=float, default=None,
                        help="discard reads with mean Phred below Q [default: off]")
    p_trim.add_argument("--json", action="store_true", help="full report as JSON")
    p_trim.set_defaults(func=_cmd_trim)

    p_sizes = sub.add_parser("sizes", help="fragment-size distribution (genomic_toolkit)")
    p_sizes.add_argument("alignment", help=".bam/.cram, .bedpe, or a fragment .bed/.tsv")
    p_sizes.add_argument("--histogram", default=None, help="write the per-length TSV here")
    p_sizes.add_argument("--json", action="store_true", help="full report as JSON")
    p_sizes.set_defaults(func=_cmd_sizes)

    p_markdup = sub.add_parser(
        "markdup", help="duplicate rate and library size, streaming (genomic_toolkit)")
    p_markdup.add_argument("alignment", help=".bam/.cram, .bedpe, or a fragment .bed/.tsv")
    p_markdup.add_argument("--min-length", type=int, default=None,
                           help="drop fragments shorter than N bp [default: 0]")
    p_markdup.add_argument("--max-length", type=int, default=None,
                           help="drop fragments longer than N bp [default: 1000, 0 = off]")
    p_markdup.add_argument("--min-mapq", type=int, default=None,
                           help="drop alignments below MAPQ N [default: 0]")
    p_markdup.add_argument("--single-end", action="store_true",
                           help="include unpaired reads, using the aligned span")
    p_markdup.add_argument("--unsorted", action="store_true",
                           help="duplicate mode for non-coordinate-sorted input")
    p_markdup.add_argument("--json", action="store_true", help="full report as JSON")
    p_markdup.set_defaults(func=_cmd_markdup)

    p_frip = sub.add_parser("frip", help="fraction of reads in peaks (genomic_toolkit)")
    p_frip.add_argument("alignment", help=".bam/.cram, .bedpe, or a fragment .bed/.tsv")
    p_frip.add_argument("peaks", help="peaks BED (merged into a disjoint cover at load time)")
    p_frip.add_argument("--genome-size", type=int, default=None,
                        help="effective genome size; enables FRiP enrichment")
    p_frip.add_argument("--json", action="store_true", help="full report as JSON")
    p_frip.set_defaults(func=_cmd_frip)

    p_profile = sub.add_parser(
        "profile", help="CUT&Tag/CUT&RUN reference-point signal matrix (cuttag_profiler)")
    p_profile.add_argument("alignment", help="indexed BAM (a .bai/.csi must exist)")
    p_profile.add_argument("regions", help="BED or GTF of loci (.gz accepted)")
    p_profile.add_argument("out_dir", help="matrix.tsv and profile.tsv are written here")
    p_profile.add_argument("--reference-point", default=None, choices=("TSS", "center", "TES"),
                           help="[default: TSS]")
    p_profile.add_argument("--upstream", type=int, default=None,
                           help="window before the anchor, bp [default: 2000]")
    p_profile.add_argument("--downstream", type=int, default=None,
                           help="window after the anchor, bp [default: 2000]")
    p_profile.add_argument("--bin-size", type=int, default=None, help="bin width, bp [default: 50]")
    p_profile.add_argument("-t", "--threads", type=int, default=None,
                           help="[default: all cores]")
    p_profile.add_argument("--json", action="store_true", help="full report as JSON")
    p_profile.set_defaults(func=_cmd_profile)

    # Every default below is None so the binary applies its own, which is the
    # rule api.py's header states. argparse would otherwise print `(default:
    # None)`-free help that reads as if this layer owned the numbers.
    p_peaks = sub.add_parser(
        "peaks", help="call peaks from a fragment BED or BAM (pto-peaks)")
    p_peaks.add_argument("fragments", help="coordinate-sorted fragment BED, or BAM")
    p_peaks.add_argument("out", help="narrowPeak output path")
    p_peaks.add_argument("--chrom-sizes", default=None,
                         help="two columns, name and length; required for BED input")
    p_peaks.add_argument("--qvalue", type=float, default=None,
                         help="Benjamini-Hochberg FDR level")
    p_peaks.add_argument("--cutoff", type=float, default=None,
                         help="-log10 p to open a peak")
    p_peaks.add_argument("--min-length", type=int, default=None,
                         help="drop candidate peaks shorter than N bp")
    p_peaks.add_argument("--max-gap", type=int, default=None,
                         help="merge candidate peaks separated by at most N bp")
    p_peaks.add_argument("--extend-peaks", type=int, default=None,
                         help="pad reported spans by N bp per side for "
                              "MACS-style pipelines; 0 keeps the sharp default")
    p_peaks.add_argument("--lambda-bg", type=float, default=None,
                         help="override the local-background lambda floor")
    p_peaks.add_argument("--min-mapq", type=int, default=None,
                         help="BAM only: skip alignments below this")
    p_peaks.add_argument("--json", action="store_true", help="full report as JSON")
    p_peaks.set_defaults(func=_cmd_peaks)

    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return args.func(args)
    except InvalidInput as exc:
        # Exit 2 means the input's shape made the answer wrong, and the tool
        # printed a valid report first. Show both, and preserve the exit code:
        # a `for` loop over 300 samples needs to tell this from a crash.
        print(f"invalid input: {exc}", file=sys.stderr)
        metrics = getattr(exc.result, "metrics", None) if exc.result is not None else None
        if metrics:
            print(json.dumps(metrics, indent=2), file=sys.stderr)
        return 2
    except (ToolError, PtoError, FileNotFoundError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    except BrokenPipeError:  # pragma: no cover - `| head`
        devnull = os.open(os.devnull, os.O_WRONLY)
        os.dup2(devnull, sys.stdout.fileno())
        return 141


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
