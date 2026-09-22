# SPDX-License-Identifier: MIT
"""Statistical equivalence metrics, implemented on the standard library.

No numpy/scipy. The validation harness has to run on the same bare cluster the
tools themselves target -- a compiler and a Python interpreter, nothing else --
so every statistic here is written out longhand. The sample sizes involved
(bins in a signal matrix, peaks in a peak set, reads in a QC report) are small
enough that the pure-Python cost is irrelevant next to the tool runs being
measured.

Every function is total: it returns a float or None, and never raises on
degenerate input (empty vectors, zero variance, all-zero tracks). A None means
"undefined for this input", which the report renders as `n/a` and treats as a
failure only when the metric was declared required.
"""

from __future__ import annotations

import math
from typing import Iterable, Sequence


def _clean_pair(a: Sequence[float], b: Sequence[float]) -> tuple[list[float], list[float]]:
    """Drop index-aligned pairs where either side is NaN/inf.

    Both tracks come from real pipelines, and a bin with no coverage can land as
    a NaN in one tool and a 0.0 in the other. Dropping only the genuinely
    non-finite entries keeps the comparison honest without silently rewriting
    one tool's zeros into the other's missing values.
    """
    xs: list[float] = []
    ys: list[float] = []
    for x, y in zip(a, b):
        if math.isfinite(x) and math.isfinite(y):
            xs.append(float(x))
            ys.append(float(y))
    return xs, ys


def pearson(a: Sequence[float], b: Sequence[float]) -> float | None:
    """Product-moment correlation, or None if either side has zero variance.

    Zero variance is the common degenerate case here, not an exotic one: two
    all-zero signal tracks are perfectly consistent with each other but have no
    defined correlation. Returning None rather than 1.0 keeps the harness from
    reporting a flat-line agreement as a passing result.
    """
    xs, ys = _clean_pair(a, b)
    n = len(xs)
    if n < 2:
        return None

    mx = sum(xs) / n
    my = sum(ys) / n
    sxy = sxx = syy = 0.0
    for x, y in zip(xs, ys):
        dx = x - mx
        dy = y - my
        sxy += dx * dy
        sxx += dx * dx
        syy += dy * dy

    if sxx <= 0.0 or syy <= 0.0:
        return None
    return sxy / math.sqrt(sxx * syy)


def _rank(values: Sequence[float]) -> list[float]:
    """Fractional ranks, averaging ties.

    Ties matter for genomic tracks specifically: long runs of exactly-zero bins
    are the norm, and ranking them by position instead of averaging would
    manufacture a correlation out of nothing but input ordering.
    """
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
            j += 1
        shared = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks[order[k]] = shared
        i = j + 1
    return ranks


def spearman(a: Sequence[float], b: Sequence[float]) -> float | None:
    """Rank correlation. Reported alongside Pearson, never instead of it.

    The two disagree in a way that is diagnostic: a high Spearman with a low
    Pearson means the tools rank loci identically but disagree on scale, which
    is normally a normalization difference (CPM vs RPKM vs raw) rather than a
    defect in the counting.
    """
    xs, ys = _clean_pair(a, b)
    if len(xs) < 2:
        return None
    return pearson(_rank(xs), _rank(ys))


def jaccard(a: Iterable, b: Iterable) -> float | None:
    """|A n B| / |A u B| over two sets of hashable items.

    Used for peak-overlap agreement, where the items are merged interval keys.
    Two empty sets return None, not 1.0: "neither tool called anything" is a
    failed run to investigate, not agreement to celebrate.
    """
    sa, sb = set(a), set(b)
    if not sa and not sb:
        return None
    return len(sa & sb) / len(sa | sb)


def interval_jaccard(
    a: Sequence[tuple[str, int, int]], b: Sequence[tuple[str, int, int]]
) -> float | None:
    """Base-pair Jaccard over two interval sets, per contig.

    This is the metric that actually corresponds to what `bedtools jaccard`
    reports, and it differs from the set-of-peaks Jaccard above in the way that
    matters: two callers that agree on where the signal is but disagree on
    boundary placement score near 1.0 here and near 0.0 on peak identity. Both
    are reported, because the gap between them localises the disagreement.

    Inputs are merged first -- overlapping intervals within one set would
    otherwise double-count their shared bases into the union.
    """
    ma, mb = _merge(a), _merge(b)
    total_a = sum(e - s for _, s, e in ma)
    total_b = sum(e - s for _, s, e in mb)
    if total_a == 0 and total_b == 0:
        return None

    inter = 0
    by_contig: dict[str, list[tuple[int, int]]] = {}
    for contig, start, end in mb:
        by_contig.setdefault(contig, []).append((start, end))

    for contig, start, end in ma:
        for other_start, other_end in by_contig.get(contig, ()):
            if other_start >= end:
                break
            if other_end > start:
                inter += min(end, other_end) - max(start, other_start)

    union = total_a + total_b - inter
    if union <= 0:
        return None
    return inter / union


def _merge(intervals: Sequence[tuple[str, int, int]]) -> list[tuple[str, int, int]]:
    """Sort and flatten intervals into a disjoint cover, per contig."""
    out: list[tuple[str, int, int]] = []
    for contig, start, end in sorted(intervals):
        if out and out[-1][0] == contig and start <= out[-1][2]:
            if end > out[-1][2]:
                out[-1] = (contig, out[-1][1], end)
        else:
            out.append((contig, start, end))
    return out


def ks_statistic(a: Sequence[float], b: Sequence[float]) -> float | None:
    """Two-sample Kolmogorov-Smirnov D: the max gap between two ECDFs.

    This is the fragment-size-distribution parity metric. D is reported rather
    than a p-value on purpose -- with millions of fragments any p-value is
    indistinguishable from zero and stops carrying information, whereas D stays
    interpretable as "the distributions differ by at most this fraction
    anywhere".
    """
    xs = sorted(v for v in a if math.isfinite(v))
    ys = sorted(v for v in b if math.isfinite(v))
    if not xs or not ys:
        return None

    # Both cursors must be advanced past *every* entry equal to the current
    # value before the ECDFs are compared. Advancing one side at a time is the
    # obvious implementation and it is wrong: on two identical samples it
    # reports D = 1/n, because it measures the gap in the middle of a tie that
    # neither distribution actually has. Fragment lengths are integers with
    # massive tie multiplicity, so this is the common case here, not a corner.
    i = j = 0
    d = 0.0
    n, m = len(xs), len(ys)
    while i < n or j < m:
        if j >= m or (i < n and xs[i] <= ys[j]):
            value = xs[i]
        else:
            value = ys[j]
        while i < n and xs[i] == value:
            i += 1
        while j < m and ys[j] == value:
            j += 1
        d = max(d, abs(i / n - j / m))
    return d


def histogram_ks(
    a: dict[int, float], b: dict[int, float]
) -> float | None:
    """KS distance between two *binned* distributions given as {value: weight}.

    Fragment size distributions arrive as histograms, not as raw samples, and
    expanding a 200-million-fragment histogram back into a sample list to reuse
    `ks_statistic` would be the one place this module actually cost real time.
    Walking the shared support and comparing running weight fractions gives the
    identical statistic without materialising anything.
    """
    total_a = sum(a.values())
    total_b = sum(b.values())
    if total_a <= 0 or total_b <= 0:
        return None

    cum_a = cum_b = 0.0
    d = 0.0
    for key in sorted(set(a) | set(b)):
        cum_a += a.get(key, 0.0)
        cum_b += b.get(key, 0.0)
        d = max(d, abs(cum_a / total_a - cum_b / total_b))
    return d


def relative_error(observed: float, expected: float) -> float | None:
    """|observed - expected| / |expected|, with an exact-zero special case.

    When the expected value is 0 the relative error is either 0 (both zero) or
    undefined, and reporting `inf` there would make the whole report unreadable
    for the common "no duplicates in a clean test file" case.
    """
    if not (math.isfinite(observed) and math.isfinite(expected)):
        return None
    if expected == 0.0:
        return 0.0 if observed == 0.0 else None
    return abs(observed - expected) / abs(expected)


def mean(values: Sequence[float]) -> float | None:
    finite = [v for v in values if math.isfinite(v)]
    if not finite:
        return None
    return sum(finite) / len(finite)
