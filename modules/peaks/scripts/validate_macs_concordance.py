#!/usr/bin/env python3
"""
MACS2/MACS3 concordance harness for pto-peaks.

Runs both callers on one BAM (or reads their narrowPeak files), pairs peaks by
best IoU, and reports recall, precision, summit agreement, mean IoU and the
Spearman correlation of -log10 p. The recorded results, and what they mean,
are in pto-core/modules/peaks/CONCORDANCE.md.

WHAT THE EXIT CODE MEANS -- READ THIS BEFORE QUOTING IT

Exit 0 means NO REGRESSION against the first real-data baseline (2026-09-11).
It does NOT mean pto-peaks is at parity with MACS3. The pre-registered parity
gate (PARITY_GATE, 0.98) is still computed, printed and written to --out-json,
and on ENCODE K562 ATAC-seq it is still not met. See SUMMIT_FLOOR below for why
the floors sit where they do.
"""

import sys
import json
import argparse
import subprocess
import os
import shutil
import tempfile

# Regression floors, re-baselined on 2026-09-11 after the first real-data run:
# ENCODE K562 ATAC-seq ENCFF121ZQX, MACS3 3.0.4 `callpeak -f BAMPE -q 0.05`,
# pto-peaks 0.1.0 at defaults. That run measured
#
#     summit within 50 bp, over matched pairs      0.9616
#     Spearman r of -log10 p, over matched pairs   0.8568
#
# The floors sit just under those numbers and were chosen AFTER seeing them.
# They are a tripwire for pto-peaks drifting from its own measured behaviour,
# not an independent biological criterion: nothing in ATAC-seq biology puts a
# line at 0.85. Two properties of the measurement limit what passing says:
#
#   * Both rates are computed over the IoU-matched pairs only -- 18,434 of
#     MACS3's 31,607 peaks (58%). Agreement on the unmatched 42% can change
#     without moving either number.
#   * Mean IoU (0.46) is reported and not gated. pto-peaks calls intervals
#     about 3x narrower (median 96 bp vs 275 bp), and a 96 bp call inside a
#     275 bp one caps IoU near 0.35, so a low IoU is what that width difference
#     predicts even where summits agree (median offset 1 bp on the 7,442 pairs
#     with >= 50% reciprocal overlap). That explains the IoU. It does not make
#     the two callers' intervals interchangeable.
#
# PARITY_GATE is the bar written before any run. It is still reported, and
# until it is met nothing should describe pto-peaks as validated against MACS3.
SUMMIT_FLOOR = 0.95
SPEARMAN_FLOOR = 0.85
PARITY_GATE = 0.98

class Peak:
    def __init__(self, chrom, start, end, pval, qval, summit_offset):
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.pval = float(pval)
        self.qval = float(qval)
        self.summit = self.start + int(summit_offset)
        self.matched = False

def load_peaks(filepath):
    peaks = {}
    with open(filepath) as f:
        for line in f:
            if not line.strip() or line.startswith('#'): continue
            parts = line.strip().split('\t')
            if len(parts) >= 10:
                p = Peak(parts[0], parts[1], parts[2], parts[7], parts[8], parts[9])
                peaks.setdefault(p.chrom, []).append(p)
    # Sort by start coordinate
    for chrom in peaks:
        peaks[chrom].sort(key=lambda x: x.start)
    return peaks

def match_peaks(macs_peaks, pto_peaks):
    matched_macs = []
    matched_pto = []
    
    for chrom in macs_peaks:
        if chrom not in pto_peaks: continue
        macs_chrom = macs_peaks[chrom]
        pto_chrom = pto_peaks[chrom]
        
        pto_idx = 0
        for mp in macs_chrom:
            # Advance until pto end is past macs start
            while pto_idx < len(pto_chrom) and pto_chrom[pto_idx].end <= mp.start:
                pto_idx += 1
            
            best_pto = None
            best_iou = -1
            
            idx = pto_idx
            while idx < len(pto_chrom) and pto_chrom[idx].start < mp.end:
                pp = pto_chrom[idx]
                overlap = min(mp.end, pp.end) - max(mp.start, pp.start)
                union = max(mp.end, pp.end) - min(mp.start, pp.start)
                iou = overlap / union if union > 0 else 0
                if iou > best_iou:
                    best_iou = iou
                    best_pto = pp
                idx += 1
                
            if best_pto is not None and not best_pto.matched:
                best_pto.matched = True
                mp.matched = True
                matched_macs.append(mp)
                matched_pto.append(best_pto)
                
    return matched_macs, matched_pto

def get_ranks(data):
    temp = [(val, i) for i, val in enumerate(data)]
    temp.sort()
    ranks = [0] * len(data)
    i = 0
    while i < len(temp):
        j = i
        while j < len(temp) and temp[j][0] == temp[i][0]:
            j += 1
        avg_rank = sum(range(i + 1, j + 1)) / (j - i)
        for k in range(i, j):
            ranks[temp[k][1]] = avg_rank
        i = j
    return ranks

def pearson(x, y):
    n = len(x)
    if n == 0: return 0.0
    mean_x = sum(x) / n
    mean_y = sum(y) / n
    num = sum((xi - mean_x) * (yi - mean_y) for xi, yi in zip(x, y))
    den = (sum((xi - mean_x)**2 for xi in x) * sum((yi - mean_y)**2 for yi in y)) ** 0.5
    return num / den if den != 0 else 0.0

def spearman(x, y):
    return pearson(get_ranks(x), get_ranks(y))

def run_macs3(bam_file, outdir):
    prefix = "macs3_test"
    cmd = ["macs3", "callpeak", "-t", bam_file, "-f", "BAMPE", "-n", prefix, "--outdir", outdir, "-q", "0.05"]
    print(f"Running MACS3: {' '.join(cmd)}")
    subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return os.path.join(outdir, f"{prefix}_peaks.narrowPeak")

def run_pto_peaks(bam_file, out_file):
    cmd = ["pto-peaks", bam_file, "-o", out_file]
    print(f"Running pto-peaks: {' '.join(cmd)}")
    subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return out_file

def main():
    parser = argparse.ArgumentParser(description="Validate pto-peaks concordance against MACS3")
    parser.add_argument("--bam", help="Test BAM file")
    parser.add_argument("--macs-peaks", help="Reference MACS3 narrowPeak file")
    parser.add_argument("--pto-peaks", help="Output pto-peaks narrowPeak file")
    parser.add_argument("--out-json", help="Path to write JSON metrics report")
    parser.add_argument("--keep-dir",
                        help="Copy both callers' narrowPeak files here (macs3_peaks.narrowPeak, "
                             "pto_peaks.narrowPeak). Without it they live in a tempdir and the "
                             "per-peak widths behind the metrics are gone when the run ends.")
    args = parser.parse_args()

    with tempfile.TemporaryDirectory() as tmpdir:
        # Resolve MACS peaks
        macs_peaks_path = args.macs_peaks
        if not macs_peaks_path:
            if not args.bam:
                sys.exit("Error: Must provide either --macs-peaks or --bam to run MACS3")
            try:
                macs_peaks_path = run_macs3(args.bam, tmpdir)
            except FileNotFoundError:
                sys.exit("Error: macs3 is not installed or not in PATH, and no --macs-peaks provided.")
            except subprocess.CalledProcessError as e:
                sys.exit(f"Error running MACS3: {e.stderr.decode()}")
                
        # Resolve PTO peaks
        pto_peaks_path = args.pto_peaks
        if not pto_peaks_path:
            if not args.bam:
                sys.exit("Error: Must provide either --pto-peaks or --bam to run pto-peaks")
            pto_peaks_path = os.path.join(tmpdir, "pto_test_peaks.narrowPeak")
            try:
                run_pto_peaks(args.bam, pto_peaks_path)
            except FileNotFoundError:
                sys.exit("Error: pto-peaks is not in PATH. Ensure the binary is built and accessible.")
            except subprocess.CalledProcessError as e:
                sys.exit(f"Error running pto-peaks: {e.stderr.decode()}")

        print("Loading peaks...")
        macs_data = load_peaks(macs_peaks_path)
        pto_data = load_peaks(pto_peaks_path)
        
        total_macs = sum(len(v) for v in macs_data.values())
        total_pto = sum(len(v) for v in pto_data.values())
        print(f"Total MACS3 peaks: {total_macs}")
        print(f"Total pto-peaks peaks: {total_pto}")

        if total_macs == 0 or total_pto == 0:
            print("Warning: One of the callers produced zero peaks. Exiting.")
            sys.exit(1)

        matched_macs, matched_pto = match_peaks(macs_data, pto_data)
        
        if not matched_macs:
            print("Error: No overlapping peaks found.")
            sys.exit(1)

        ious = []
        macs_p = []
        pto_p = []
        macs_q = []
        pto_q = []
        
        summit_within_50 = 0
        for m, p in zip(matched_macs, matched_pto):
            diff = abs(m.summit - p.summit)
            if diff <= 50:
                summit_within_50 += 1
                
            overlap = min(m.end, p.end) - max(m.start, p.start)
            union = max(m.end, p.end) - min(m.start, p.start)
            ious.append(overlap / union if union > 0 else 0)
            
            macs_p.append(m.pval)
            pto_p.append(p.pval)
            macs_q.append(m.qval)
            pto_q.append(p.qval)

        recall = len(matched_macs) / total_macs
        precision = len(matched_pto) / total_pto
        summit_rate = summit_within_50 / len(matched_macs)
        mean_iou = sum(ious) / len(ious)
        p_corr = spearman(macs_p, pto_p)
        q_corr = spearman(macs_q, pto_q)

        metrics = {
            "total_macs3_peaks": total_macs,
            "total_pto_peaks": total_pto,
            "matched_peaks": len(matched_macs),
            "recall": recall,
            "precision": precision,
            "summit_overlap_rate": summit_rate,
            "mean_iou": mean_iou,
            "spearman_p_value": p_corr,
            "spearman_q_value": q_corr
        }
        # Recorded in the JSON so a reader of the artefact alone sees both
        # verdicts, not just the one the exit code carries.
        parity_met = summit_rate >= PARITY_GATE and p_corr >= PARITY_GATE
        metrics["regression_floors"] = {"summit_overlap_rate": SUMMIT_FLOOR,
                                        "spearman_p_value": SPEARMAN_FLOOR,
                                        "set_on": "2026-09-11"}
        metrics["parity_gate"] = PARITY_GATE
        metrics["parity_met"] = parity_met

        print("\n=== Concordance Metrics ===")
        for k, v in metrics.items():
            print(f"{k}: {v:.4f}" if isinstance(v, float) else f"{k}: {v}")

        if args.out_json:
            with open(args.out_json, "w") as f:
                json.dump(metrics, f, indent=4)
            print(f"Wrote metrics to {args.out_json}")

        # Before the gates, so a failing run -- the current state -- keeps the
        # calls that explain why it failed.
        if args.keep_dir:
            os.makedirs(args.keep_dir, exist_ok=True)
            shutil.copyfile(macs_peaks_path, os.path.join(args.keep_dir, "macs3_peaks.narrowPeak"))
            shutil.copyfile(pto_peaks_path, os.path.join(args.keep_dir, "pto_peaks.narrowPeak"))
            print(f"Kept both callers' narrowPeak files in {args.keep_dir}")

        print(f"\nParity with MACS3 (pre-registered {PARITY_GATE} gates): "
              f"{'met' if parity_met else 'NOT met'}")

        errors = []
        if summit_rate < SUMMIT_FLOOR:
            errors.append(f"Summit overlap rate {summit_rate:.4f} is below the "
                          f"{SUMMIT_FLOOR} regression floor")
        if p_corr < SPEARMAN_FLOOR:
            errors.append(f"Spearman p-value correlation {p_corr:.4f} is below the "
                          f"{SPEARMAN_FLOOR} regression floor")

        if errors:
            print("\nREGRESSION against the 2026-09-11 baseline:")
            for e in errors:
                print(f" - {e}")
            sys.exit(1)
        else:
            print("\nNO REGRESSION: summit agreement and strength correlation are at or above "
                  "the 2026-09-11 floors. This is not a parity result -- see the line above.")
            sys.exit(0)

if __name__ == "__main__":
    main()
