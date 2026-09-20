#!/usr/bin/env python3
"""Generate a synthetic CUT&Tag dataset for smoke tests and benchmarking.

Emits a coordinate-sorted BAM plus a matching promoter BED. Fragments are
drawn so that a known fraction sit under a TSS-centred peak and the rest are
uniform background, which gives the profiler a ground truth to hit.

Requires only `samtools` on PATH; nothing here needs pysam or a reference
genome.

    python3 scripts/make_demo_data.py --out-dir data --regions 2000 --reads 500000
"""

from __future__ import annotations

import argparse
import os
import random
import subprocess
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out-dir", type=Path, default=Path("data"))
    parser.add_argument("--prefix", default="demo")
    parser.add_argument("--regions", type=int, default=2000,
                        help="number of promoter loci")
    parser.add_argument("--reads", type=int, default=500_000,
                        help="number of fragments (each becomes two mates)")
    parser.add_argument("--chrom-length", type=int, default=10_000_000)
    parser.add_argument("--read-length", type=int, default=50)
    parser.add_argument("--peak-fraction", type=float, default=0.55,
                        help="fraction of fragments drawn from the TSS peak")
    parser.add_argument("--peak-sigma", type=float, default=180.0,
                        help="standard deviation of the peak, in bp")
    parser.add_argument("--seed", type=int, default=42)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    rng = random.Random(args.seed)

    if not 0.0 <= args.peak_fraction <= 1.0:
        sys.exit("--peak-fraction must lie in [0, 1]")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    bam_path = args.out_dir / f"{args.prefix}.bam"
    bed_path = args.out_dir / f"{args.prefix}_promoters.bed"
    sam_path = args.out_dir / f"{args.prefix}.sam"

    chrom = "chr1"
    margin = 20_000
    usable = args.chrom_length - 2 * margin
    if usable <= args.regions:
        sys.exit("--chrom-length is too small for the requested region count")

    # Evenly spaced TSSs, jittered so the spacing is not perfectly periodic.
    spacing = usable // args.regions
    tss_sites: list[tuple[int, str]] = []
    for i in range(args.regions):
        pos = margin + i * spacing + rng.randint(0, max(1, spacing // 4))
        strand = "+" if rng.random() < 0.5 else "-"
        tss_sites.append((pos, strand))

    with bed_path.open("w") as bed:
        for i, (pos, strand) in enumerate(tss_sites):
            # A 200 bp promoter whose anchor (BED start for +, end for -) is
            # the TSS, matching how the profiler picks its reference point.
            start, end = (pos, pos + 200) if strand == "+" else (pos - 200, pos)
            bed.write(f"{chrom}\t{start}\t{end}\tgene{i:05d}\t0\t{strand}\n")

    peak_count = int(args.reads * args.peak_fraction)
    read_len = args.read_length

    fragments: list[tuple[int, int]] = []  # (start, fragment_length)
    for _ in range(peak_count):
        pos, _strand = tss_sites[rng.randrange(len(tss_sites))]
        centre = int(rng.gauss(pos, args.peak_sigma))
        length = max(2 * read_len, int(rng.gauss(160, 40)))  # nucleosome-ish
        start = centre - length // 2
        if margin <= start <= args.chrom_length - margin - length:
            fragments.append((start, length))

    for _ in range(args.reads - peak_count):
        length = max(2 * read_len, int(rng.gauss(200, 60)))
        start = rng.randint(margin, args.chrom_length - margin - length)
        fragments.append((start, length))

    fragments.sort()

    print(f"writing {len(fragments) * 2} alignments to {sam_path}", file=sys.stderr)
    with sam_path.open("w") as sam:
        sam.write("@HD\tVN:1.6\tSO:coordinate\n")
        sam.write(f"@SQ\tSN:{chrom}\tLN:{args.chrom_length}\n")
        sam.write("@RG\tID:demo\tSM:demo\tLB:cuttag\tPL:ILLUMINA\n")

        seq = "A" * read_len
        qual = "I" * read_len
        cigar = f"{read_len}M"

        # Mates are emitted together; samtools sort fixes the ordering.
        for i, (start, length) in enumerate(fragments):
            name = f"frag{i:08d}"
            mate_start = start + length - read_len
            # 99/147: proper pair, read1 forward / read2 reverse.
            sam.write(f"{name}\t99\t{chrom}\t{start + 1}\t60\t{cigar}\t=\t"
                      f"{mate_start + 1}\t{length}\t{seq}\t{qual}\tRG:Z:demo\n")
            sam.write(f"{name}\t147\t{chrom}\t{mate_start + 1}\t60\t{cigar}\t=\t"
                      f"{start + 1}\t{-length}\t{seq}\t{qual}\tRG:Z:demo\n")

    print("sorting and indexing", file=sys.stderr)
    subprocess.run(["samtools", "sort", "-o", str(bam_path), str(sam_path)],
                   check=True)
    subprocess.run(["samtools", "index", str(bam_path)], check=True)
    os.remove(sam_path)

    print(f"\nBAM     {bam_path}", file=sys.stderr)
    print(f"regions {bed_path} ({args.regions} loci)", file=sys.stderr)
    print(f"peak    {args.peak_fraction:.0%} of fragments, sigma {args.peak_sigma:.0f} bp",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
