#!/usr/bin/env python3
"""Benchmark harness: Block_CSR/knn_graph vs scanpy.pp.neighbors / scipy.sparse.

Usage:
    python benchmark_1m_cells.py --input path/to/matrix.h5ad --backend scrna --k 15
    python benchmark_1m_cells.py --input path/to/matrix.h5ad --backend scanpy --k 15

Both backends accept the same .h5ad (or .npz) input, so a head-to-head run is a
matter of changing --backend only.

Reports wall-clock time and peak RSS for k-NN graph construction, per
docs/BENCHMARKS.md's protocol. `perf stat -e cache-misses` should wrap this
script externally to capture cache-miss counts, e.g.:

    perf stat -e cache-misses,cache-references \
        python benchmark_1m_cells.py --input 1M_neurons.h5ad --backend scrna

Scale warning: the scrna backend's k-NN is exact brute force, i.e. O(n_rows^2).
It is a correctness baseline, not a 1M-cell solution -- at 1.3M cells that is
~1e12 similarity evaluations. Use --max-cells to subsample to a size that
actually terminates, and see docs/BENCHMARKS.md for what would be required to
make the 1M+ claim real (an approximate index, which does not exist yet).

To generate a synthetic dataset for scale testing without staging real data:
    python benchmark_1m_cells.py --make-synthetic 50000 --output synth.h5ad
"""

from __future__ import annotations

import argparse
import resource
import sys
import time


def peak_rss_gb() -> float:
    """Peak resident set size in GB. ru_maxrss is KB on Linux, bytes on macOS."""
    raw = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform == "darwin":
        return raw / (1024.0**3)
    return raw / (1024.0**2)


def load_csr(input_path: str, max_cells: int | None = None):
    """Load a .h5ad or .npz into a scipy CSR matrix, optionally subsampled."""
    import scipy.sparse as sp

    if input_path.endswith(".npz"):
        csr = sp.load_npz(input_path).tocsr()
    elif input_path.endswith(".h5ad"):
        try:
            import anndata
        except ImportError as exc:
            raise SystemExit("anndata not installed; `pip install anndata` to read .h5ad.") from exc
        adata = anndata.read_h5ad(input_path)
        csr = sp.csr_matrix(adata.X)
    else:
        raise SystemExit(f"Unsupported input format for {input_path!r}; expected .h5ad or .npz.")

    if max_cells is not None and csr.shape[0] > max_cells:
        csr = csr[:max_cells]
    return csr


def run_scrna_backend(
    input_path: str, k: int, max_cells: int | None,
    method: str = "brute", n_pcs: int = 0,
) -> tuple[float, float]:
    try:
        import scrna_matrix_py as scrna
    except ImportError as exc:
        raise SystemExit(
            "scrna_matrix_py not built/installed. Build the CMake project with "
            "SCRNA_BUILD_PYTHON=ON and add the build dir to PYTHONPATH."
        ) from exc

    if method == "hnsw" and not scrna.has_hnsw():
        raise SystemExit("built without -DSCRNA_ENABLE_HNSW=ON; method='hnsw' unavailable.")

    csr = load_csr(input_path, max_cells)

    # --n-pcs reduces to a dense embedding first. This is the configuration HNSW
    # is for: it indexes dense vectors, so searching raw gene space costs
    # n_cells * n_genes * 4 bytes and is usually slower than the exact path.
    if n_pcs > 0:
        import numpy as np
        from scipy.sparse.linalg import svds

        u, s, _ = svds(csr.asfptype(), k=n_pcs)
        rep = np.ascontiguousarray((u * s)[:, ::-1], dtype=np.float32)
        est = scrna.hnsw_index_bytes(rep.shape[0], rep.shape[1]) / 2**30
        print(
            f"  [scrna] cells={rep.shape[0]} dims={rep.shape[1]} method={method} "
            f"simd={scrna.simd_isa()} est_index={est:.2f}GB",
            file=sys.stderr,
        )
        t0 = time.perf_counter()
        graph = scrna.build_knn_graph_dense(rep, k, method="hnsw")
        elapsed = time.perf_counter() - t0
        _ = graph.to_numpy()
        return elapsed, peak_rss_gb()

    # to_aligned guarantees the 64-byte alignment that lets BlockCSR adopt these
    # buffers in place rather than copying them.
    data = scrna.to_aligned(csr.data.astype("float32"))
    indices = scrna.to_aligned(csr.indices.astype("int32"))
    indptr = scrna.to_aligned(csr.indptr.astype("int32"))

    mat = scrna.BlockCSR.from_scipy_csr(data, indices, indptr, csr.shape[0], csr.shape[1])
    if method == "hnsw":
        est = scrna.hnsw_index_bytes(csr.shape[0], csr.shape[1]) / 2**30
        print(
            f"  [scrna] WARNING: HNSW over {csr.shape[1]} gene dimensions needs "
            f"~{est:.2f}GB and is usually slower than method='brute'. "
            f"Consider --n-pcs 50.",
            file=sys.stderr,
        )
    print(
        f"  [scrna] cells={csr.shape[0]} genes={csr.shape[1]} nnz={csr.nnz} "
        f"method={method} simd={scrna.simd_isa()} zero_copy={mat.is_zero_copy}",
        file=sys.stderr,
    )

    t0 = time.perf_counter()
    graph = scrna.build_knn_graph(mat, k, method=method)
    elapsed = time.perf_counter() - t0
    _ = graph.to_numpy()  # exercise zero-copy export path
    return elapsed, peak_rss_gb()


def run_scanpy_baseline(input_path: str, k: int, max_cells: int | None) -> tuple[float, float]:
    try:
        import anndata
        import scanpy as sc
    except ImportError as exc:
        raise SystemExit("scanpy not installed; `pip install scanpy` to run the baseline.") from exc

    if input_path.endswith(".h5ad"):
        adata = sc.read_h5ad(input_path)
        if max_cells is not None and adata.n_obs > max_cells:
            adata = adata[:max_cells].copy()
    else:
        adata = anndata.AnnData(X=load_csr(input_path, max_cells))

    t0 = time.perf_counter()
    sc.pp.neighbors(adata, n_neighbors=k, metric="cosine")
    elapsed = time.perf_counter() - t0
    return elapsed, peak_rss_gb()


def make_synthetic(n_cells: int, output: str, n_genes: int = 2000, density: float = 0.05) -> None:
    """Write a synthetic scRNA-seq-shaped .h5ad for scale testing.

    Deliberately not committed to the repo -- generate locally, as with the real
    10x datasets listed in docs/BENCHMARKS.md.
    """
    import anndata
    import numpy as np
    import scipy.sparse as sp

    rng = np.random.default_rng(0)
    x = sp.random(n_cells, n_genes, density=density, format="csr", random_state=rng)
    x.data = np.abs(rng.poisson(3.0, size=x.data.shape).astype("float32")) + 1.0
    adata = anndata.AnnData(X=x)
    adata.write_h5ad(output)
    print(f"wrote {output}: {n_cells} cells x {n_genes} genes, nnz={x.nnz}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", help="Path to .h5ad or .npz (both backends accept either)")
    parser.add_argument("--k", type=int, default=15, help="Number of nearest neighbors")
    parser.add_argument("--backend", choices=["scrna", "scanpy"], default="scrna")
    parser.add_argument(
        "--max-cells", type=int, default=None,
        help="Subsample to at most this many cells (brute force is O(n^2); use this to keep runs finite)",
    )
    parser.add_argument(
        "--method", choices=["brute", "hnsw"], default="brute",
        help="scrna backend only: exact brute force, or approximate HNSW",
    )
    parser.add_argument(
        "--n-pcs", type=int, default=0,
        help="Reduce to this many components before searching (recommended with --method hnsw)",
    )
    parser.add_argument("--make-synthetic", type=int, metavar="N_CELLS",
                        help="Generate a synthetic dataset of N_CELLS instead of benchmarking")
    parser.add_argument("--output", help="Output path for --make-synthetic")
    args = parser.parse_args()

    if args.make_synthetic:
        if not args.output:
            raise SystemExit("--make-synthetic requires --output")
        make_synthetic(args.make_synthetic, args.output)
        return

    if not args.input:
        raise SystemExit("--input is required (or use --make-synthetic)")

    if args.backend == "scrna":
        elapsed, rss_gb = run_scrna_backend(
            args.input, args.k, args.max_cells, args.method, args.n_pcs
        )
        label = f"backend=scrna method={args.method}" + (f" n_pcs={args.n_pcs}" if args.n_pcs else "")
    else:
        elapsed, rss_gb = run_scanpy_baseline(args.input, args.k, args.max_cells)
        label = "backend=scanpy"

    print(f"{label} k={args.k} wall_clock_s={elapsed:.3f} peak_rss_gb={rss_gb:.3f}")


if __name__ == "__main__":
    main()
