#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Self-test for the statistics the validation suite's verdicts rest on.

    python tests/validation/test_metrics.py

Framework-free and dependency-free, in the same style as the C++ suites: a
`main` that runs checks and returns non-zero on the first failure. Pulling in
pytest here would contradict the rule that keeps the rest of the repo's tests
runnable on a bare cluster.

This exists because `metrics.py` is the layer that decides whether a run passes
or fails, and a defect in it is silent by nature -- a subtly wrong correlation
still returns a plausible number in [-1, 1], and the suite would keep printing
green. These checks run in milliseconds; the alternative is discovering the
same bug after a five-hour run against a 3 GB BAM.

The two-sample KS walk in particular had exactly that defect: advancing one
cursor at a time reports D = 1/n for two *identical* samples, because it
measures a gap in the middle of a tie that neither distribution has. Fragment
lengths are integers with enormous tie multiplicity, so that is the common
case here, not a corner. `test_ks_ties` pins it.
"""

from __future__ import annotations

import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from pto_validation import metrics as M  # noqa: E402

FAILURES: list[str] = []


def check(condition: bool, description: str) -> None:
    if not condition:
        FAILURES.append(description)


def close(a, b, tolerance: float = 1e-12) -> bool:
    return a is not None and abs(a - b) <= tolerance


def test_pearson() -> None:
    check(M.pearson([1, 2, 3], [2, 4, 6]) == 1.0, "pearson: perfect positive")
    check(close(M.pearson([1, 2, 3], [6, 4, 2]), -1.0), "pearson: perfect negative")
    # Undefined, not 1.0: two flat tracks are consistent but uncorrelated, and
    # reporting agreement there would pass a run that validated nothing.
    check(M.pearson([1, 1, 1], [1, 2, 3]) is None, "pearson: zero variance -> None")
    check(M.pearson([1.0], [1.0]) is None, "pearson: n<2 -> None")
    check(M.pearson([], []) is None, "pearson: empty -> None")
    check(
        M.pearson([1, float("nan"), 3], [2, 5, 6]) == 1.0,
        "pearson: non-finite pairs dropped, not zero-filled",
    )


def test_spearman() -> None:
    check(close(M.spearman([1, 2, 3], [10, 20, 30]), 1.0), "spearman: monotone")
    check(M.spearman([1, 2, 3], [3, 2, 1]) == -1.0, "spearman: reversed")
    # Rank correlation is scale-free where Pearson is not; this is the pair
    # that distinguishes a normalization difference from a counting bug.
    check(
        close(M.spearman([1, 2, 3, 4], [1, 4, 9, 16]), 1.0),
        "spearman: nonlinear but monotone still 1.0",
    )
    check(M._rank([5, 5, 5]) == [2.0, 2.0, 2.0], "rank: ties averaged")
    check(M._rank([10, 20, 20, 30]) == [1.0, 2.5, 2.5, 4.0], "rank: partial ties")


def test_jaccard() -> None:
    check(M.jaccard([1, 2, 3], [2, 3, 4]) == 0.5, "jaccard: half overlap")
    check(M.jaccard([1], [2]) == 0.0, "jaccard: disjoint")
    check(M.jaccard([], []) is None, "jaccard: both empty -> None, not 1.0")


def test_interval_jaccard() -> None:
    check(
        close(M.interval_jaccard([("c", 0, 100)], [("c", 50, 150)]), 50 / 150),
        "interval_jaccard: partial overlap",
    )
    check(
        M.interval_jaccard([("c", 0, 10)], [("d", 0, 10)]) == 0.0,
        "interval_jaccard: same coordinates on different contigs never overlap",
    )
    check(close(M.interval_jaccard([("c", 0, 10)], [("c", 0, 10)]), 1.0),
          "interval_jaccard: identical")
    # Overlapping inputs on one side must be merged first, or their shared
    # bases get counted twice into the union and the score sags.
    check(
        close(M.interval_jaccard([("c", 0, 10), ("c", 5, 15)], [("c", 0, 15)]), 1.0),
        "interval_jaccard: overlapping inputs merged before comparison",
    )
    check(M.interval_jaccard([], []) is None, "interval_jaccard: both empty -> None")


def test_ks_ties() -> None:
    """The regression that motivated this file."""
    check(M.ks_statistic([1, 2, 3, 4], [1, 2, 3, 4]) == 0.0,
          "ks: identical samples -> 0")
    check(M.ks_statistic([1, 1, 1], [1, 1, 1]) == 0.0,
          "ks: identical constant samples -> 0")
    check(close(M.ks_statistic([1, 2, 3, 4], [5, 6, 7, 8]), 1.0),
          "ks: disjoint supports -> 1")
    check(close(M.ks_statistic([0, 1], [1, 2]), 0.5), "ks: half shifted")
    check(M.ks_statistic([], [1]) is None, "ks: empty -> None")


def test_histogram_ks_matches_ks() -> None:
    """histogram_ks is an optimisation; it must not be a different statistic.

    It exists only so a 200-million-fragment histogram does not have to be
    expanded back into a sample list. If the two ever disagree, the binned
    path is silently reporting a different number on exactly the large inputs
    where nobody can check it by hand.
    """
    rng = random.Random(20260817)
    for trial in range(50):
        a = [rng.randint(1, 20) for _ in range(rng.randint(10, 400))]
        b = [rng.randint(1, 25) for _ in range(rng.randint(10, 300))]
        ha: dict[int, float] = {}
        hb: dict[int, float] = {}
        for v in a:
            ha[v] = ha.get(v, 0.0) + 1.0
        for v in b:
            hb[v] = hb.get(v, 0.0) + 1.0
        direct = M.ks_statistic(a, b)
        binned = M.histogram_ks(ha, hb)
        check(
            direct is not None and binned is not None and abs(direct - binned) < 1e-12,
            f"histogram_ks disagrees with ks_statistic on trial {trial}: "
            f"{direct} vs {binned}",
        )

    check(close(M.histogram_ks({1: 10}, {2: 10}), 1.0), "histogram_ks: disjoint")
    check(M.histogram_ks({}, {1: 1}) is None, "histogram_ks: empty -> None")


def test_relative_error() -> None:
    check(M.relative_error(0, 0) == 0.0, "relative_error: 0 vs 0 is exact")
    check(M.relative_error(1, 0) is None, "relative_error: nonzero vs 0 undefined")
    check(close(M.relative_error(11, 10), 0.1), "relative_error: 10%")
    check(close(M.relative_error(9, 10), 0.1), "relative_error: symmetric")
    check(M.relative_error(float("inf"), 10) is None, "relative_error: non-finite")


def test_mean() -> None:
    check(M.mean([1, 2, 3]) == 2.0, "mean")
    check(M.mean([1, float("nan"), 3]) == 2.0, "mean: skips nan")
    check(M.mean([]) is None, "mean: empty -> None")


def main() -> int:
    tests = [value for name, value in sorted(globals().items())
             if name.startswith("test_") and callable(value)]
    for test in tests:
        test()

    if FAILURES:
        print(f"FAILED ({len(FAILURES)}):", file=sys.stderr)
        for failure in FAILURES:
            print(f"  - {failure}", file=sys.stderr)
        return 1

    print(f"ok: {len(tests)} metric test groups passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
