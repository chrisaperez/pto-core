#!/usr/bin/env python3
"""Quantify agreement between a cuttag_profiler matrix and a deeptools one.

The two tools normalise on different denominators (deeptools scales a bigWig
of per-bin coverage; cuttag_profiler scales fragment counts straight off the
BAM), so an exact bit-for-bit match is not the right target. What must hold for
the tools to be interchangeable in practice is that the *shape* of the signal
is the same: the curves agree up to a single global scale factor.

Reported statistics:
  * Pearson r between the two meta-profiles;
  * Pearson r between matched region rows, summarised across regions;
  * the largest residual after fitting one global scale factor.

Depends only on the standard library, so it runs wherever the profiler builds.

    python3 scripts/compare_matrices.py --profiler p.tsv --deeptools d.tsv
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path


def read_profiler_matrix(path: Path) -> tuple[list[float], list[list[float]]]:
    """Returns (bin offsets, rows). Row 0 of the file is the offset header."""
    rows: list[list[float]] = []
    offsets: list[float] = []
    with path.open() as handle:
        header = handle.readline().rstrip("\n").split("\t")
        offsets = [float(v) for v in header[1:]]
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue
            rows.append([float(v) for v in fields[1:]])
    return offsets, rows


def read_deeptools_matrix(path: Path) -> list[list[float]]:
    """deeptools --outFileNameMatrix writes three header lines, then raw rows."""
    rows: list[list[float]] = []
    with path.open() as handle:
        for index, line in enumerate(handle):
            if index < 3:
                continue
            fields = line.rstrip("\n").split("\t")
            if not fields or fields == [""]:
                continue
            rows.append([0.0 if f in ("", "nan", "NA") else float(f) for f in fields])
    return rows


def pearson(a: list[float], b: list[float]) -> float:
    n = len(a)
    if n == 0 or n != len(b):
        return float("nan")
    mean_a = sum(a) / n
    mean_b = sum(b) / n
    num = sum((x - mean_a) * (y - mean_b) for x, y in zip(a, b))
    den_a = math.sqrt(sum((x - mean_a) ** 2 for x in a))
    den_b = math.sqrt(sum((y - mean_b) ** 2 for y in b))
    if den_a == 0.0 or den_b == 0.0:
        return float("nan")
    return num / (den_a * den_b)


def column_means(rows: list[list[float]]) -> list[float]:
    if not rows:
        return []
    cols = len(rows[0])
    sums = [0.0] * cols
    for row in rows:
        for c in range(min(cols, len(row))):
            sums[c] += row[c]
    return [s / len(rows) for s in sums]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--profiler", type=Path, required=True)
    parser.add_argument("--deeptools", type=Path, required=True)
    parser.add_argument("--row-sample", type=int, default=2000,
                        help="how many matched rows to correlate individually")
    parser.add_argument("--min-r", type=float, default=0.99,
                        help="fail the run if the profile correlation is below this")
    args = parser.parse_args()

    _offsets, p_rows = read_profiler_matrix(args.profiler)
    d_rows = read_deeptools_matrix(args.deeptools)

    if not p_rows or not d_rows:
        sys.exit("one of the matrices is empty")

    n_rows = min(len(p_rows), len(d_rows))
    n_cols = min(len(p_rows[0]), len(d_rows[0]))
    print(f"profiler matrix : {len(p_rows)} x {len(p_rows[0])}")
    print(f"deeptools matrix: {len(d_rows)} x {len(d_rows[0])}")
    print(f"compared        : {n_rows} x {n_cols}")
    if len(p_rows) != len(d_rows):
        print("  note: row counts differ; comparing the leading rows only")

    p_profile = column_means([row[:n_cols] for row in p_rows[:n_rows]])
    d_profile = column_means([row[:n_cols] for row in d_rows[:n_rows]])

    profile_r = pearson(p_profile, d_profile)

    # One global scale factor, fitted by least squares through the origin.
    num = sum(p * d for p, d in zip(p_profile, d_profile))
    den = sum(p * p for p in p_profile)
    scale = num / den if den else float("nan")

    residuals = [abs(scale * p - d) for p, d in zip(p_profile, d_profile)]
    peak = max(d_profile) if d_profile else 0.0
    max_residual_pct = 100.0 * max(residuals) / peak if peak else float("nan")

    # Per-region agreement over a sample of matched rows.
    step = max(1, n_rows // max(1, args.row_sample))
    row_rs = []
    for i in range(0, n_rows, step):
        p_row = p_rows[i][:n_cols]
        d_row = d_rows[i][:n_cols]
        r = pearson(p_row, d_row)
        if not math.isnan(r):
            row_rs.append(r)

    print()
    print(f"meta-profile Pearson r      : {profile_r:.6f}")
    print(f"fitted global scale factor  : {scale:.6g}")
    print(f"max residual (% of peak)    : {max_residual_pct:.3f}%")
    if row_rs:
        row_rs.sort()
        median = row_rs[len(row_rs) // 2]
        print(f"per-region r (n={len(row_rs)})")
        print(f"  median                    : {median:.6f}")
        print(f"  10th percentile           : {row_rs[len(row_rs) // 10]:.6f}")
        print(f"  fraction above 0.95       : "
              f"{sum(1 for r in row_rs if r > 0.95) / len(row_rs):.4f}")

    if math.isnan(profile_r) or profile_r < args.min_r:
        print(f"\nFAIL: meta-profile correlation {profile_r:.6f} < {args.min_r}")
        return 1
    print(f"\nPASS: meta-profile correlation >= {args.min_r}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
