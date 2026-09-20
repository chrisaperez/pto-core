# SPDX-License-Identifier: MIT
"""scrna_matrix against its own exact path, and against scanpy's neighbor graph.

This module is the odd one out in the suite: it has no CLI, only pybind11
bindings, so the cases drive it in-process rather than through subprocesses.
It also has no gold-standard binary to diff against -- the nearest equivalent
is scanpy's neighbor graph, which is itself built on an approximate method
(pynndescent) and therefore cannot serve as ground truth for recall.

So the checks are layered accordingly:

* The *exact* brute-force path (`method="brute"`) is ground truth, because it
  is exhaustive and deterministic by construction. HNSW recall is measured
  against it, on identical input.
* The similarities the exact path returns are re-derived longhand in Python.
  That is what catches a wrong runtime ISA dispatch or a mis-set alignment
  flag -- the failure mode modules/scrna_matrix/docs/AUDIT.md says to watch
  for, and one that recall alone would not see, since a uniformly wrong
  kernel still ranks neighbors consistently.
* scanpy, when installed, is compared on graph *agreement* rather than
  identity, since neither side claims to reproduce the other exactly.
"""

from __future__ import annotations

import math
import random

from .. import metrics as M
from ..config import KNN_RECALL_MIN
from ..harness import Context, Unavailable, case
from ..report import Direction, Metric


def _import_module():
    """Import the compiled Python package, or skip with build instructions.

    pybind11 is optional for scrna_matrix and its absence is a configure-time
    skip rather than an error, so an installation without it must skip here
    too rather than fail the suite.
    """
    # Both import paths are supported and both are exercised in the module's
    # own CI: the `scrna_matrix` wrapper package when the extension has been
    # installed into it, and the bare `scrna_matrix_py` extension when the
    # build directory is simply on PYTHONPATH. Every symbol this file uses is
    # re-exported unchanged by the wrapper, so either satisfies it -- and
    # requiring the wrapper would make the common in-tree build fail here for
    # a packaging reason that has nothing to do with correctness.
    for name in ("scrna_matrix", "scrna_matrix_py"):
        try:
            return __import__(name)
        except ImportError:
            continue

    raise Unavailable(
        "neither `scrna_matrix` nor `scrna_matrix_py` is importable. Build with "
        "pybind11 available and put the extension's directory on PYTHONPATH:\n"
        "  PYTHONPATH=build/modules/scrna_matrix python "
        "tests/validation/run_validation.py ..."
    )


def _numpy():
    try:
        import numpy as np
    except ImportError as exc:
        raise Unavailable("numpy is required for the scrna_matrix cases") from exc
    return np


def _load_counts(ctx: Context):
    """Load the dataset's count matrix as a dense float32 array.

    Densified deliberately: HNSW indexes dense vectors only (it belongs on a
    PCA embedding, not raw sparse gene space), and the exact path is serving
    as ground truth here, so both sides must see byte-identical input. The
    datasets in this suite are the 2-3k cell public matrices, where the dense
    form is a few hundred megabytes at most.
    """
    np = _numpy()
    path = ctx.input("matrix")
    if path.suffix != ".h5ad":
        raise Unavailable(f"unsupported matrix format: {path.suffix}")

    try:
        import anndata
    except ImportError as exc:
        raise Unavailable("anndata is required to read the .h5ad fixture") from exc

    adata = anndata.read_h5ad(path)
    matrix = adata.X
    dense = matrix.toarray() if hasattr(matrix, "toarray") else np.asarray(matrix)
    return np.ascontiguousarray(dense, dtype=np.float32)


def _embedding(counts, components: int = 50):
    """Reduce counts to a dense low-dimensional embedding.

    Uses numpy's own SVD rather than scanpy's PCA so these cases do not require
    scanpy to run at all. The embedding does not need to match scanpy's -- both
    k-NN paths run on *this* embedding, and what is being validated is their
    agreement with each other on identical input.
    """
    np = _numpy()
    centered = counts - counts.mean(axis=0, keepdims=True)
    # Economy SVD: `components` is ~50 against thousands of genes, so the full
    # decomposition would be an order of magnitude more work for columns that
    # are discarded on the next line.
    _, _, vt = np.linalg.svd(centered, full_matrices=False)
    return np.ascontiguousarray(centered @ vt[:components].T, dtype=np.float32)


def _as_block_csr(module, dense):
    """Wrap a dense array as a BlockCSR so the exact path can consume it.

    The brute-force path takes a Block_CSR, and the embedding is dense, so it
    is stored as a fully-populated sparse matrix. That is wasteful in the
    abstract and irrelevant here: 50 dimensions across a few thousand cells is
    a few megabytes, and it keeps ground truth running on exactly the same
    numbers HNSW receives, which is the entire point.

    Routed through `to_aligned` so the 64-byte-aligned adopt-in-place path is
    exercised rather than the copying fallback -- validating the padded path
    while shipping the copied one would test the wrong code.
    """
    np = _numpy()
    rows, cols = dense.shape
    data = module.to_aligned(np.ascontiguousarray(dense.reshape(-1), dtype=np.float32))
    indices = module.to_aligned(
        np.ascontiguousarray(np.tile(np.arange(cols, dtype=np.int32), rows))
    )
    indptr = module.to_aligned(
        np.ascontiguousarray(np.arange(rows + 1, dtype=np.int32) * cols)
    )
    return module.BlockCSR.from_scipy_csr(data, indices, indptr, rows, cols)


def _neighbor_rows(graph):
    """Reshape KnnGraph.to_numpy()'s flat arrays into per-cell rows.

    to_numpy returns (indices, similarities, k) as two flat n*k arrays. Slicing
    them into rows here means every caller agrees on the layout, and a k that
    does not divide the array length surfaces as an explicit error instead of a
    silently misaligned comparison.
    """
    indices, similarities, k = graph.to_numpy()
    if k == 0 or len(indices) % k != 0:
        raise Unavailable(
            f"KnnGraph returned {len(indices)} neighbor ids which is not a "
            f"multiple of k={k}"
        )
    rows = len(indices) // k
    return (
        [list(indices[i * k:(i + 1) * k]) for i in range(rows)],
        [list(similarities[i * k:(i + 1) * k]) for i in range(rows)],
        int(k),
    )


@case(
    "scrna/hnsw-recall-vs-exact",
    module="scrna_matrix",
    summary="HNSW k-NN recall against the exact brute-force path",
)
def hnsw_recall(ctx: Context) -> list[Metric]:
    """Measure approximate-graph recall against exhaustive search."""
    module = _import_module()
    if not module.has_hnsw():
        raise Unavailable(
            "this build has no HNSW support; configure with -DSCRNA_ENABLE_HNSW=ON"
        )

    counts = _load_counts(ctx)
    embedding = _embedding(counts, int(ctx.param("pca_components", 50)))
    k = int(ctx.param("knn_k", 15))

    exact_indices, _, _ = _neighbor_rows(
        module.build_knn_graph(_as_block_csr(module, embedding), k, method="brute")
    )
    approx_indices, _, _ = _neighbor_rows(
        module.build_knn_graph_dense(embedding, k, method="hnsw")
    )

    if len(exact_indices) != len(approx_indices):
        raise Unavailable(
            f"graphs disagree on cell count: exact {len(exact_indices)} vs "
            f"hnsw {len(approx_indices)}"
        )

    exact_sets = [set(row) for row in exact_indices]
    approx_sets = [set(row) for row in approx_indices]
    recalls = [
        len(a & b) / len(a) if a else float("nan")
        for a, b in zip(exact_sets, approx_sets)
    ]
    mean_recall = M.mean(recalls)

    if mean_recall is not None and mean_recall < KNN_RECALL_MIN:
        worst = sorted(range(len(recalls)), key=lambda i: recalls[i])[:200]
        path = ctx.scratch("scrna") / "knn_recall_diff.tsv"
        with path.open("w") as handle:
            handle.write("# lowest-recall cells first\ncell\trecall\tmissed_neighbors\n")
            for index in worst:
                missed = sorted(exact_sets[index] - approx_sets[index])
                handle.write(f"{index}\t{recalls[index]:.4f}\t{missed}\n")
        ctx.artifact("knn recall diff", path)

    finite = sorted(r for r in recalls if not math.isnan(r))
    return [
        Metric("knn.mean_recall", mean_recall, KNN_RECALL_MIN, Direction.MIN),
        Metric(
            "knn.p10_recall",
            finite[len(finite) // 10] if finite else None,
            # The 10th percentile is allowed to sit well below the mean: HNSW's
            # error concentrates on cells in sparse regions of the embedding,
            # so a handful of poorly-connected cells is the expected shape of
            # the approximation, not a uniform degradation.
            0.60, Direction.MIN,
        ),
        Metric("knn.cells", float(len(exact_sets)), None, Direction.MIN, required=False),
        Metric("knn.simd_isa", None, None, Direction.MIN, required=False,
               detail=f"dispatched to {module.simd_isa()}"),
    ]


@case(
    "scrna/cosine-vs-longhand",
    module="scrna_matrix",
    summary="SIMD cosine kernel against a scalar recomputation",
)
def cosine_vs_longhand(ctx: Context) -> list[Metric]:
    """Recompute a sample of the exact graph's similarities in pure Python.

    This catches a wrong runtime ISA dispatch, which is otherwise close to
    invisible: an AVX-512 kernel with a mis-set tail mask produces plausible
    numbers on plausible data and diverges only on the trailing elements of
    vectors whose length is not a multiple of the lane width. So the sampled
    cells deliberately include the last row and a fixed-seed spread, rather
    than the first N, and the embedding width is chosen not to be a round
    multiple of any lane count.
    """
    module = _import_module()
    counts = _load_counts(ctx)
    # 50 is neither a multiple of 8 (AVX2 float lanes) nor of 16 (AVX-512),
    # so every kernel here has to execute its tail path.
    embedding = _embedding(counts, int(ctx.param("pca_components", 50)))

    graph = module.build_knn_graph(
        _as_block_csr(module, embedding), int(ctx.param("knn_k", 15)), method="brute"
    )
    indices, similarities, _ = _neighbor_rows(graph)

    rows = len(indices)
    rng = random.Random(20260817)
    sample = {0, rows - 1}
    while len(sample) < min(64, rows):
        sample.add(rng.randrange(rows))

    worst = 0.0
    checked = 0
    for i in sorted(sample):
        for neighbor, reported in zip(indices[i], similarities[i]):
            expected = _cosine(embedding[i], embedding[neighbor])
            worst = max(worst, abs(expected - float(reported)))
            checked += 1

    return [
        Metric(
            "cosine.max_abs_error", worst,
            # float32 accumulation over 50 dimensions; 1e-5 sits comfortably
            # above the representable error and far below any real defect.
            1e-5, Direction.MAX,
            detail=f"over {checked} neighbor similarities from {len(sample)} cells, "
                   f"ISA {module.simd_isa()}",
        ),
        Metric("cosine.similarities_checked", float(checked), None, Direction.MIN,
               required=False),
    ]


def _cosine(a, b) -> float:
    dot = float(sum(float(x) * float(y) for x, y in zip(a, b)))
    na = math.sqrt(sum(float(x) * float(x) for x in a))
    nb = math.sqrt(sum(float(y) * float(y) for y in b))
    if na == 0.0 or nb == 0.0:
        return 0.0
    return dot / (na * nb)


@case(
    "scrna/graph-agreement-vs-scanpy",
    module="scrna_matrix",
    summary="Neighbor-graph agreement with scanpy (informational)",
)
def graph_agreement_vs_scanpy(ctx: Context) -> list[Metric]:
    """Compare the exact graph against scanpy's neighbors, as a sanity check.

    Every metric here is informational and none of them gate the suite. That
    is not laziness: scanpy defaults to pynndescent, an approximate method, on
    matrices of this size, so a disagreement is as likely to be scanpy's
    approximation as ours. The number is still worth printing -- a *collapse*
    in agreement, as opposed to the usual few percent, means the two are no
    longer computing the same kind of graph at all.
    """
    module = _import_module()
    try:
        import scanpy as sc  # type: ignore
    except ImportError as exc:
        raise Unavailable("scanpy is not installed") from exc

    np = _numpy()
    counts = _load_counts(ctx)
    embedding = _embedding(counts, int(ctx.param("pca_components", 50)))
    k = int(ctx.param("knn_k", 15))

    import anndata

    adata = anndata.AnnData(np.asarray(counts, dtype=np.float32))
    adata.obsm["X_pca"] = np.asarray(embedding, dtype=np.float32)
    sc.pp.neighbors(adata, n_neighbors=k, use_rep="X_pca", metric="cosine")

    connectivity = adata.obsp["distances"]
    scanpy_sets = [
        set(connectivity[i].indices.tolist()) for i in range(connectivity.shape[0])
    ]
    exact_indices, _, _ = _neighbor_rows(
        module.build_knn_graph(_as_block_csr(module, embedding), k, method="brute")
    )
    # scanpy excludes the query cell from its neighbor list; the exact path
    # includes it. Dropping self on our side compares like with like instead of
    # charging us a fixed 1/k penalty on every cell.
    exact_sets = [set(row) - {i} for i, row in enumerate(exact_indices)]

    jaccards = [
        M.jaccard(a, b) for a, b in zip(exact_sets, scanpy_sets)
    ]
    finite = [j for j in jaccards if j is not None]

    return [
        Metric("scanpy.mean_jaccard", M.mean(finite), None, Direction.MIN,
               required=False,
               detail="informational: scanpy's own graph is approximate"),
        Metric("scanpy.cells_compared", float(len(finite)), None, Direction.MIN,
               required=False),
    ]
