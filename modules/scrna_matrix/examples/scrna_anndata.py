"""AnnData integration helpers for scrna_matrix.

Provides a drop-in alternative to `scanpy.pp.neighbors` that computes the k-NN
graph with the Block_CSR engine and writes the result back into the AnnData
object using exactly the layout scanpy expects, so downstream tools
(`sc.tl.leiden`, `sc.tl.umap`, `sc.tl.paga`, ...) work unchanged:

    import scanpy as sc
    from scrna_anndata import neighbors

    adata = sc.datasets.pbmc3k_processed()
    neighbors(adata, n_neighbors=15)
    sc.tl.leiden(adata)

Connectivities are produced by scanpy's own fuzzy-simplicial-set routine when it
is importable, which makes them identical to what `sc.pp.neighbors` would have
written for the same neighbour sets. If that private helper is unavailable
(scanpy moved it), a clearly-labelled symmetrised fallback is used instead and
`adata.uns['neighbors']['params']['connectivity_impl']` records which ran.

Scale note: the underlying k-NN is exact brute force, O(n_obs^2). It is well
suited to tens of thousands of cells and is not a substitute for an approximate
index at 1M+ cells -- see docs/BENCHMARKS.md.
"""

from __future__ import annotations

import warnings
from typing import Any

import numpy as np
import scipy.sparse as sp

try:
    import scrna_matrix_py as _scrna
except ImportError as exc:  # pragma: no cover - import-time guidance
    raise ImportError(
        "scrna_matrix_py not found. Build with -DSCRNA_BUILD_PYTHON=ON and put the "
        "build directory on PYTHONPATH."
    ) from exc


def _as_csr(matrix: Any) -> sp.csr_matrix:
    """Return a CSR matrix with canonical (deduplicated, sorted) indices.

    Block_CSR requires strictly increasing column indices per row — the
    merge-join kernel walks two rows in lockstep and silently mis-computes
    otherwise. SciPy does not guarantee sortedness, so normalise here rather
    than let the C++ layer reject input the caller had no reason to suspect.
    Both calls are no-ops when the flags are already set.
    """
    if sp.isspmatrix_csr(matrix):
        out = matrix
    elif sp.issparse(matrix):
        out = matrix.tocsr()
    else:
        out = sp.csr_matrix(np.asarray(matrix))

    if not out.has_canonical_format:
        out = out.copy()          # never mutate the caller's adata.X in place
        out.sum_duplicates()      # implies sort_indices()
    return out


def blockcsr_from_anndata(adata, layer: str | None = None):
    """Build a BlockCSR from `adata.X` (or a named layer).

    The three CSR buffers are passed through `to_aligned` so BlockCSR adopts
    them in place instead of copying. The returned tuple's second element must be
    kept referenced for as long as the matrix is used -- it owns the memory the
    matrix borrows.
    """
    csr = _as_csr(adata.layers[layer] if layer is not None else adata.X)

    # Handed over at their own width, never pre-cast.
    #
    # This used to read `.astype("int32")`, and NumPy's astype WRAPS: a column
    # index of 2**32 + 5 became 5 before to_aligned() -- whose whole job is to
    # range-check that narrowing -- ever saw the real value. The check existed
    # and was bypassed by the caller. to_aligned() narrows int64->int32 and
    # float64->float32 itself and raises instead of wrapping.
    # (docs/AUDIT_2026-09-11_scrna_matrix.md, SM8.)
    #
    # float64 rather than float32 for the values: it is lossless for every
    # dtype AnnData stores here (float32, int32, and int64 counts below 2**53),
    # and leaves the one narrowing to the checked path.
    data = _scrna.to_aligned(csr.data.astype("float64", copy=False))
    indices = _scrna.to_aligned(csr.indices)
    indptr = _scrna.to_aligned(csr.indptr)

    mat = _scrna.BlockCSR.from_scipy_csr(
        data, indices, indptr, csr.shape[0], csr.shape[1]
    )
    return mat, (data, indices, indptr)


def knn_arrays(
    adata,
    n_neighbors: int = 15,
    layer: str | None = None,
    *,
    method: str = "brute",
    use_rep: str | None = None,
    **hnsw_kwargs,
):
    """Compute the k-NN graph and return (neighbor_idx, cosine_sim) as (n, k).

    `method="brute"` is exact and deterministic; `method="hnsw"` is approximate
    and scales to 1M+ cells.

    `use_rep` names a key in `adata.obsm` (e.g. `"X_pca"`) to search in place of
    `adata.X`. For `method="hnsw"` this matters a great deal: hnswlib indexes
    dense vectors, so searching in full gene space both costs
    n_obs x n_genes x 4 bytes and makes every distance an n_genes-wide dense
    operation — which is typically *slower* than the sparse exact path, not
    faster. A PCA embedding is the configuration where HNSW actually wins.
    """
    n = adata.shape[0]

    if use_rep is not None:
        rep = np.ascontiguousarray(adata.obsm[use_rep], dtype=np.float32)
        if method != "hnsw":
            raise ValueError("use_rep requires method='hnsw'")
        graph = _scrna.build_knn_graph_dense(rep, n_neighbors, method="hnsw", **hnsw_kwargs)
    else:
        mat, _keepalive = blockcsr_from_anndata(adata, layer=layer)
        graph = _scrna.build_knn_graph(mat, n_neighbors, method=method, **hnsw_kwargs)

    idx, sim, k = graph.to_numpy()
    return idx.reshape(n, k).copy(), sim.reshape(n, k).copy()


def _distances_matrix(idx: np.ndarray, sim: np.ndarray, n_obs: int) -> sp.csr_matrix:
    """Cosine *distance* CSR, self excluded, matching scanpy's obsp['distances']."""
    k = idx.shape[1]
    rows = np.repeat(np.arange(n_obs, dtype=np.int32), k)
    cols = idx.reshape(-1)
    dist = (1.0 - sim.reshape(-1)).astype(np.float64)

    # Drop unfilled slots (-1 sentinel, emitted when k exceeds the available
    # neighbour count) and self-edges, explicitly.
    #
    # Self-edges must be filtered here rather than via setdiag(0) +
    # eliminate_zeros(): a cosine distance of exactly 0.0 is a *legitimate* value
    # for two identical cells, and eliminate_zeros() cannot tell that apart from
    # the zeroed diagonal. Duplicated expression profiles are common in scRNA-seq
    # (doublets, low-complexity cells, technical replicates), and the previous
    # version silently deleted precisely those edges -- the strongest ones in the
    # graph -- leaving affected cells with fewer than k neighbours.
    keep = (cols >= 0) & (cols != rows)
    matrix = sp.csr_matrix(
        (dist[keep], (rows[keep], cols[keep])), shape=(n_obs, n_obs)
    )
    # Note: coo->csr sums duplicate (row, col) pairs. The top-k producers emit
    # distinct neighbours per row, so this cannot trigger; asserting it keeps a
    # future approximate backend from silently doubling an edge weight.
    assert matrix.nnz == int(keep.sum()), "duplicate neighbour entries were summed"
    return matrix


def _connectivities(idx: np.ndarray, sim: np.ndarray, n_obs: int) -> tuple[sp.csr_matrix, str]:
    """Fuzzy-simplicial-set connectivities, via scanpy's own routine when possible."""
    k = idx.shape[1]

    # scanpy/UMAP convention: column 0 is the point itself at distance 0, and
    # n_neighbors counts that self entry. Verified against sc.pp.neighbors output.
    knn_idx = np.empty((n_obs, k + 1), dtype=np.int32)
    knn_dist = np.empty((n_obs, k + 1), dtype=np.float64)
    knn_idx[:, 0] = np.arange(n_obs, dtype=np.int32)
    knn_dist[:, 0] = 0.0
    knn_idx[:, 1:] = idx
    knn_dist[:, 1:] = 1.0 - sim

    # Self-fill any unused slots so the array stays well formed. This is a
    # degenerate case (it means k >= n_obs - 1), and the repeated zero-distance
    # self entries distort UMAP's local-connectivity estimate, so say so rather
    # than papering over it.
    missing = knn_idx[:, 1:] < 0
    if missing.any():
        warnings.warn(
            f"{int(missing.any(axis=1).sum())} cell(s) have fewer than {k} neighbours; "
            "padding with self-references. Connectivities for those cells are "
            "unreliable -- reduce n_neighbors.",
            RuntimeWarning,
            stacklevel=3,
        )
        rows = np.repeat(np.arange(n_obs, dtype=np.int32)[:, None], k, axis=1)
        knn_idx[:, 1:][missing] = rows[missing]
        knn_dist[:, 1:][missing] = 0.0

    # Only the *import* is optional. Wrapping the call in the same `except`
    # would swallow genuine computation errors (bad shapes, non-finite
    # distances) and silently substitute a different algorithm, so the failure
    # would surface as subtly wrong connectivities rather than an error.
    try:
        from scanpy.neighbors._connectivity import umap as _umap_conn
    except (ImportError, AttributeError):
        _umap_conn = None

    if _umap_conn is not None:
        conn = _umap_conn(knn_idx, knn_dist, n_obs=n_obs, n_neighbors=k + 1)
        return sp.csr_matrix(conn), "scanpy.umap"

    warnings.warn(
        "scanpy's fuzzy-simplicial-set routine is unavailable; falling back to a "
        "symmetrised similarity graph. This is NOT the UMAP connectivity, so "
        "embeddings will differ from sc.pp.neighbors.",
        RuntimeWarning,
        stacklevel=3,
    )
    rows = np.repeat(np.arange(n_obs, dtype=np.int32), k)
    cols = idx.reshape(-1)
    vals = sim.reshape(-1).astype(np.float32)
    keep = (cols >= 0) & (cols != rows)
    adj = sp.csr_matrix((vals[keep], (rows[keep], cols[keep])), shape=(n_obs, n_obs))
    adj = adj.maximum(adj.T)
    return adj, "symmetrised-fallback"


def neighbors(
    adata,
    n_neighbors: int = 15,
    *,
    method: str = "brute",
    use_rep: str | None = None,
    layer: str | None = None,
    key_added: str | None = None,
    copy: bool = False,
    **hnsw_kwargs,
):
    """scanpy-compatible neighbours, computed with the Block_CSR engine.

    Writes `obsp['distances']`, `obsp['connectivities']` and `uns['neighbors']`
    (or the `key_added`-prefixed equivalents), so downstream scanpy tools work
    against the result exactly as they would after `sc.pp.neighbors`.

    Choosing a method:

    * `method="brute"` (default) — exact, deterministic, O(n_obs^2). Operates
      directly on the sparse matrix. Best up to roughly 20-50k cells.
    * `method="hnsw", use_rep="X_pca"` — approximate, near-linear, the
      configuration for 1M+ cells. Measured at 37 s / 1.0 GB for 1.3M cells on a
      50-dimensional embedding.

    `method="hnsw"` without `use_rep` searches full gene space, which is
    supported but usually a poor trade: hnswlib is dense, so it loses the
    sparsity the exact path exploits. See docs/BENCHMARKS.md.
    """
    adata = adata.copy() if copy else adata
    n_obs = adata.shape[0]

    idx, sim = knn_arrays(
        adata, n_neighbors=n_neighbors, layer=layer, method=method,
        use_rep=use_rep, **hnsw_kwargs,
    )
    dists = _distances_matrix(idx, sim, n_obs)
    conn, conn_impl = _connectivities(idx, sim, n_obs)

    if key_added is None:
        conn_key, dist_key, uns_key = "connectivities", "distances", "neighbors"
    else:
        conn_key = f"{key_added}_connectivities"
        dist_key = f"{key_added}_distances"
        uns_key = key_added

    adata.obsp[dist_key] = dists
    adata.obsp[conn_key] = conn
    adata.uns[uns_key] = {
        "connectivities_key": conn_key,
        "distances_key": dist_key,
        "params": {
            # scanpy counts the self entry, so k stored neighbours -> k + 1.
            "n_neighbors": n_neighbors + 1,
            "metric": "cosine",
            # Named for what actually ran, rather than borrowing scanpy's
            # 'umap'/'gauss' labels for a different neighbour search.
            "method": f"scrna_matrix_{method}",
            "use_rep": use_rep or "X",
            "connectivity_impl": conn_impl,
            # Brute force is exact and deterministic. HNSW is neither: its graph
            # depends on insertion order, so multi-threaded builds vary run to run.
            "exact": method == "brute",
            "random_state": 0,
            "simd": _scrna.simd_isa(),
        },
    }
    return adata if copy else None
