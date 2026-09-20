"""Cache-aligned sparse matrix engine and SIMD k-NN for single-cell data.

Thin re-export of the compiled extension. The extension is built by CMake and
installed into this package directory (see the install() rule in
``CMakeLists.txt``), so the import below is a sibling module, not a separate
top-level one.

    import scrna_matrix

    m = scrna_matrix.BlockCSR.from_scipy_csr(...)
    g = scrna_matrix.build_knn_graph(m, k=15, method="hnsw")

    # Build a reference atlas index once, query it for the rest of time.
    idx = scrna_matrix.HnswIndex(dim=50, capacity=1_300_000)
    idx.add_dense(atlas_pca)
    idx.save_index("pbmc1M.hnsw")
    ...
    idx = scrna_matrix.HnswIndex.load_index("pbmc1M.hnsw")
    nbr, sim, k = idx.query_dense(new_cells_pca, k=15).to_numpy()

The SIMD kernel is chosen at runtime, not at build time: a wheel built on a
baseline x86-64 host still uses AVX-512 on a host that has it. ``simd_isa()``
reports the kernel actually selected, and ``SCRNA_FORCE_ISA=scalar|neon|avx2|
avx512`` in the environment pins it, which is the supported way to check that a
fallback path is correct on hardware you do have.
"""

from __future__ import annotations

from .scrna_matrix_py import (  # noqa: F401
    BlockCSR,
    KnnGraph,
    build_knn_graph,
    build_knn_graph_dense,
    has_hnsw,
    hnsw_index_bytes,
    simd_isa,
    to_aligned,
)

__all__ = [
    "BlockCSR",
    "KnnGraph",
    "build_knn_graph",
    "build_knn_graph_dense",
    "has_hnsw",
    "hnsw_index_bytes",
    "simd_isa",
    "to_aligned",
]

# HnswIndex only exists when the engine was compiled with -DSCRNA_ENABLE_HNSW=ON.
# Published wheels always are, but a local build may not be, and an ImportError
# from this module would then be blamed on packaging rather than on the flag.
try:
    from .scrna_matrix_py import HnswIndex  # noqa: F401

    __all__.append("HnswIndex")
except ImportError:  # pragma: no cover - depends on build configuration
    pass

__version__ = "0.1.0"
