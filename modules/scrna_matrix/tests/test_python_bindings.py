#!/usr/bin/env python3
"""Tests for the pybind11 bindings: adoption vs copy, buffer lifetime, GIL release.

Run with the built extension on PYTHONPATH:

    PYTHONPATH=build python -m pytest tests/test_python_bindings.py -v
"""

from __future__ import annotations

import gc
import sys
import threading
import time

import numpy as np
import pytest
import scipy.sparse as sp

# Two import paths, both supported and both exercised in CI:
#   scrna_matrix_py  -- the bare extension, as an in-tree CMake build leaves it
#                       on PYTHONPATH.
#   scrna_matrix     -- the installed wheel, where the extension sits inside the
#                       package. cibuildwheel runs this file against that one,
#                       so the suite must not assume the in-tree layout.
try:
    import scrna_matrix_py as scrna
except ImportError:
    scrna = pytest.importorskip(
        "scrna_matrix",
        reason="build the CMake project with SCRNA_BUILD_PYTHON=ON and put the build "
        "dir on PYTHONPATH, or pip install the package",
    )

ALIGN = 64


def is_aligned(arr: np.ndarray) -> bool:
    return arr.ctypes.data % ALIGN == 0


def random_csr(n_rows: int = 60, n_cols: int = 80, density: float = 0.15, seed: int = 0):
    rng = np.random.default_rng(seed)
    m = sp.random(n_rows, n_cols, density=density, format="csr", random_state=rng)
    m.data = np.abs(m.data) + 0.1  # keep values positive, as counts are
    return m


def make_blockcsr(csr, aligned: bool):
    data = csr.data.astype("float32")
    indices = csr.indices.astype("int32")
    indptr = csr.indptr.astype("int32")
    if aligned:
        data = scrna.to_aligned(data)
        indices = scrna.to_aligned(indices)
        indptr = scrna.to_aligned(indptr)
    return scrna.BlockCSR.from_scipy_csr(data, indices, indptr, csr.shape[0], csr.shape[1]), (
        data,
        indices,
        indptr,
    )


def reference_knn(csr, k: int):
    """Dense cosine top-k, computed independently with numpy."""
    dense = np.asarray(csr.todense(), dtype=np.float64)
    norms = np.linalg.norm(dense, axis=1)
    norms[norms == 0] = 1.0
    sim = (dense @ dense.T) / np.outer(norms, norms)
    np.fill_diagonal(sim, -np.inf)  # exclude self
    order = np.argsort(-sim, axis=1, kind="stable")[:, :k]
    return order, np.take_along_axis(sim, order, axis=1)


def test_module_smoke():
    assert scrna.simd_isa() in {"scalar", "neon", "avx2", "avx512"}


def test_to_aligned_produces_aligned_buffers():
    raw = np.arange(100, dtype="float32")
    # An arbitrary numpy array is not guaranteed aligned; the aligned copy must be.
    aligned = scrna.to_aligned(raw)
    assert is_aligned(aligned)
    assert aligned.dtype == np.float32
    np.testing.assert_allclose(raw, aligned)


def misaligned_copy(arr: np.ndarray) -> np.ndarray:
    """Return an equal-valued array that is deliberately not 64-byte aligned."""
    pad = np.empty(arr.size + 16, dtype=arr.dtype)
    for off in range(1, 17):
        view = pad[off : off + arr.size]
        if not is_aligned(view):
            view[:] = arr
            return view
    raise AssertionError("could not construct a misaligned view")


def test_misaligned_arrays_fall_back_to_copying():
    """Whatever NumPy happens to hand us, an unaligned buffer must be copied."""
    csr = random_csr()
    data = misaligned_copy(csr.data.astype("float32"))
    indices = misaligned_copy(csr.indices.astype("int32"))
    indptr = misaligned_copy(csr.indptr.astype("int32"))
    assert not is_aligned(data)

    mat = scrna.BlockCSR.from_scipy_csr(data, indices, indptr, csr.shape[0], csr.shape[1])
    assert mat.shape == csr.shape
    assert mat.nnz == csr.nnz
    assert mat.is_zero_copy is False


def test_wrong_dtype_falls_back_to_copying():
    csr = random_csr()
    # float64 data / int64 indices: correct values, wrong dtype for adoption.
    mat = scrna.BlockCSR.from_scipy_csr(
        csr.data.astype("float64"), csr.indices.astype("int64"),
        csr.indptr.astype("int64"), csr.shape[0], csr.shape[1],
    )
    assert mat.is_zero_copy is False
    assert mat.nnz == csr.nnz


def test_writeable_numpy_buffers_are_never_adopted():
    """
    A stock, writeable NumPy array is COPIED even when it is perfectly aligned.

    This test previously asserted the opposite -- that realistic buffers are
    adopted in place without the caller doing anything -- and that behaviour was
    memory-unsafe. Borrowing a buffer Python can still mutate or reallocate
    means Block_CSR::validate() is a check-once on bytes that keep changing,
    while the AVX gather kernels index a dense scratch buffer with those
    (re-)validated indices using unchecked hardware gathers. Three escapes were
    confirmed under ASan; see the regression tests below.
    """
    csr = random_csr(n_rows=500, n_cols=400, density=0.1, seed=21)
    data = csr.data.astype("float32")
    indices = csr.indices.astype("int32")
    indptr = csr.indptr.astype("int32")
    assert data.flags.writeable
    mat = scrna.BlockCSR.from_scipy_csr(data, indices, indptr, *csr.shape)
    assert mat.is_zero_copy is False
    assert mat.nnz == csr.nnz

    # The in-place path is still available -- through to_aligned(), whose output
    # is frozen and capsule-backed.
    frozen = tuple(scrna.to_aligned(a) for a in (data, indices, indptr))
    assert all(not a.flags.writeable for a in frozen)
    adopted = scrna.BlockCSR.from_scipy_csr(*frozen, *csr.shape)
    assert adopted.is_zero_copy is True


# --------------------------------------------------------------------------
# Memory-safety regressions for the adopt path.
#
# Each of these reproduced a distinct ASan-confirmed memory-corruption bug
# against the pre-fix bindings, all reachable from pure Python with no ctypes:
#   * mutating `indices`  -> heap-buffer-overflow WRITE (arbitrary offset/value)
#   * mutating `indptr`   -> heap-buffer-overflow READ
#   * resize(refcheck=False) -> heap-use-after-free
# --------------------------------------------------------------------------

def _adopted_matrix(n_rows=4, n_cols=8):
    data = scrna.to_aligned(np.arange(1, n_rows * n_cols + 1, dtype=np.float32))
    indices = scrna.to_aligned(np.tile(np.arange(n_cols, dtype=np.int32), n_rows))
    indptr = scrna.to_aligned(
        np.arange(0, n_rows * n_cols + 1, n_cols, dtype=np.int32)
    )
    mat = scrna.BlockCSR.from_scipy_csr(data, indices, indptr, n_rows, n_cols)
    assert mat.is_zero_copy
    return mat, data, indices, indptr


def test_adopted_buffers_are_immutable():
    """Post-adoption writes must be refused, not silently corrupt the heap."""
    mat, data, indices, indptr = _adopted_matrix()
    for name, arr in (("data", data), ("indices", indices), ("indptr", indptr)):
        assert not arr.flags.writeable, name
        with pytest.raises(ValueError):
            arr[0] = 1 << 30
    # Still usable and still correct after the refusals.
    idx, sim, k = scrna.build_knn_graph(mat, 2).to_numpy()
    assert k == 2 and np.isfinite(sim).all()


def test_adopted_buffers_cannot_be_reallocated():
    """
    ndarray.resize(refcheck=False) reallocates and would leave the adopted
    pointer dangling. The keep-alive pins the *object*, not the buffer, and a
    live buffer export does not stop NumPy either -- both were measured. The
    defence is that to_aligned() output does not own its data, so resize fails.
    """
    _mat, data, indices, indptr = _adopted_matrix()
    for arr in (data, indices, indptr):
        with pytest.raises(ValueError):
            arr.resize(4, refcheck=False)


def test_adopted_buffers_cannot_be_unfrozen():
    """
    Re-enabling the writeable flag would reopen every escape above. NumPy allows
    it for arrays that own their data, which is exactly why adoption requires a
    capsule-backed array whose base is not writeable.
    """
    _mat, data, indices, indptr = _adopted_matrix()
    for arr in (data, indices, indptr):
        with pytest.raises(ValueError):
            arr.flags.writeable = True


def test_readonly_view_over_writeable_base_is_not_adopted():
    """
    The aliasing hole that defeated a flags-based adoption test.

    `v = w[:]; v.flags.writeable = False` produces a read-only *view* whose base
    `w` is still writeable. Under an `is_adoptable` that only asked
    `not arr.flags.writeable`, this was adopted, and mutating `w` afterwards
    restored the full heap-buffer-overflow WRITE primitive (re-confirmed under
    ASan). Flags describe one handle; they say nothing about other handles onto
    the same bytes. Adoption is gated on provenance instead.
    """
    n_rows, n_cols = 64, 64
    writeable = [
        np.ones(n_rows * n_cols, dtype=np.float32),
        np.tile(np.arange(n_cols, dtype=np.int32), n_rows),
        np.arange(0, n_rows * n_cols + 1, n_cols, dtype=np.int32),
    ]

    def readonly_view(a):
        v = a[:]
        v.flags.writeable = False
        return v

    views = [readonly_view(a) for a in writeable]
    for v, w in zip(views, writeable):
        assert not v.flags.writeable      # the view lies about immutability
        assert w.flags.writeable          # the base is the real story

    mat = scrna.BlockCSR.from_scipy_csr(*views, n_rows, n_cols)
    assert mat.is_zero_copy is False

    # Mutating the base is now harmless: the matrix owns a private copy.
    writeable[1][0] = n_cols + 4
    idx, sim, _ = scrna.build_knn_graph(mat, 2, method="brute").to_numpy()
    assert np.isfinite(sim).all()
    assert idx.min() >= -1 and idx.max() < n_rows


def test_only_our_own_capsule_backed_buffers_are_adopted():
    """Provenance, not flags: even a view of an adoptable array is copied."""
    _mat, data, indices, indptr = _adopted_matrix()
    assert type(data.base).__name__ == "PyCapsule"

    # A whole-array view of a frozen array is itself frozen, aligned and
    # contiguous -- but its base is the ndarray, not our capsule.
    view = data[:]
    assert not view.flags.writeable
    copied = scrna.BlockCSR.from_scipy_csr(view, indices, indptr, 4, 8)
    assert copied.is_zero_copy is False


def test_adopted_matrix_outlives_every_python_reference():
    """
    Lifetime anchoring: the capsule chain (array -> capsule -> AlignedVector)
    must keep the bytes alive with no Python name left for any of them.
    """
    n_rows, n_cols = 8, 8

    def buffers():
        return (
            scrna.to_aligned(np.arange(1, n_rows * n_cols + 1, dtype=np.float32)),
            scrna.to_aligned(np.tile(np.arange(n_cols, dtype=np.int32), n_rows)),
            scrna.to_aligned(np.arange(0, n_rows * n_cols + 1, n_cols, dtype=np.int32)),
        )

    mat = scrna.BlockCSR.from_scipy_csr(*buffers(), n_rows, n_cols)
    assert mat.is_zero_copy
    for _ in range(3):
        gc.collect()

    idx, sim, k = scrna.build_knn_graph(mat, 3, method="brute").to_numpy()
    assert k == 3
    assert np.isfinite(sim).all()
    assert mat.nnz == n_rows * n_cols


def test_adoption_does_not_leak_references():
    """
    NumpyOwner incref's three objects per adopt and must release all of them.
    A leak here is not just memory: it pins the source buffers forever, which
    silently defeats the whole point of the copy/adopt distinction.
    """
    n_rows, n_cols = 24, 24
    data = scrna.to_aligned(np.ones(n_rows * n_cols, dtype=np.float32))
    indices = scrna.to_aligned(np.tile(np.arange(n_cols, dtype=np.int32), n_rows))
    indptr = scrna.to_aligned(
        np.arange(0, n_rows * n_cols + 1, n_cols, dtype=np.int32)
    )

    before = tuple(sys.getrefcount(a) for a in (data, indices, indptr))
    for _ in range(50):
        mat = scrna.BlockCSR.from_scipy_csr(data, indices, indptr, n_rows, n_cols)
        assert mat.is_zero_copy
        del mat
    gc.collect()
    after = tuple(sys.getrefcount(a) for a in (data, indices, indptr))
    assert before == after, f"reference leak: {before} -> {after}"


def test_concurrent_adopt_and_knn_from_python_threads():
    """
    build_knn_graph releases the GIL, so NumpyOwner's incref/decref can run
    concurrently with other threads' Python work. Those must be GIL-guarded --
    an unguarded decref is silent heap corruption inside CPython rather than a
    clean failure, so this exercises the path rather than trusting the comment.
    """
    n_rows, n_cols = 24, 24

    def buffers():
        return (
            scrna.to_aligned(np.ones(n_rows * n_cols, dtype=np.float32)),
            scrna.to_aligned(np.tile(np.arange(n_cols, dtype=np.int32), n_rows)),
            scrna.to_aligned(
                np.arange(0, n_rows * n_cols + 1, n_cols, dtype=np.int32)
            ),
        )

    errors = []

    def worker():
        try:
            for _ in range(25):
                mat = scrna.BlockCSR.from_scipy_csr(*buffers(), n_rows, n_cols)
                assert mat.is_zero_copy
                scrna.build_knn_graph(mat, 3, method="brute")
                del mat
        except BaseException as exc:  # noqa: BLE001 - surfaced via `errors`
            errors.append(exc)

    threads = [threading.Thread(target=worker) for _ in range(8)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    assert not errors, errors
    gc.collect()


def test_knn_arrays_outlive_their_graph():
    """to_numpy() exports capsule-backed views; dropping the graph must not
    invalidate them."""
    csr = random_csr(n_rows=40, n_cols=50, seed=11)
    mat, _keep = make_blockcsr(csr, aligned=True)
    graph = scrna.build_knn_graph(mat, 3, method="brute")
    idx, sim, k = graph.to_numpy()
    del graph
    for _ in range(3):
        gc.collect()
    assert k == 3
    assert idx.shape == (40 * 3,)
    assert np.isfinite(sim).all()


def test_rows_above_the_magnitude_ceiling_are_rejected():
    """
    The SIMD dot kernels accumulate in float, so a large enough row saturates
    the dot product to +inf and the `!isfinite -> 0` guard reports 0 -- two
    identical rows scoring "completely dissimilar", silently.

    simd_math.hpp's dispatchers now catch that and recompute in double, so this
    ceiling is no longer the only thing standing between a user and a wrong
    neighbour graph. It is kept because it keeps the fallback unreachable for
    anything built through this API -- the fast path stays fast -- and because
    a squared norm above 1e34 is not expression data under any normalisation.
    """
    big = np.float32(1e30)
    data = np.array([big, big], dtype=np.float32)
    indices = np.array([0, 2], dtype=np.int32)
    indptr = np.array([0, 2], dtype=np.int32)
    with pytest.raises(ValueError, match="ceiling"):
        scrna.BlockCSR.from_scipy_csr(data, indices, indptr, 1, 4)

    # Large but safe magnitudes still build.
    ok = np.array([np.float32(1e15)] * 2, dtype=np.float32)
    mat = scrna.BlockCSR.from_scipy_csr(ok, indices, indptr, 1, 4)
    assert mat.nnz == 2


def test_aligned_arrays_are_adopted_zero_copy():
    csr = random_csr()
    mat, arrays = make_blockcsr(csr, aligned=True)
    assert all(is_aligned(a) for a in arrays)
    assert mat.is_zero_copy is True
    assert mat.nnz == csr.nnz
    assert mat.shape == csr.shape


def test_adopted_buffers_survive_source_deletion():
    """The owner token must pin the NumPy buffers past the caller's last reference."""
    csr = random_csr(n_rows=40, n_cols=50, seed=3)
    mat, arrays = make_blockcsr(csr, aligned=True)
    assert mat.is_zero_copy

    expected_nnz = mat.nnz
    del arrays
    gc.collect()

    # If the keep-alive were missing this would read freed memory.
    assert mat.nnz == expected_nnz
    graph = scrna.build_knn_graph(mat, 5)
    idx, sim, k = graph.to_numpy()
    assert k == 5
    assert idx.shape == (40 * 5,)
    assert np.isfinite(sim).all()


def test_zero_copy_and_copy_paths_agree():
    csr = random_csr(seed=7)
    mat_copy = scrna.BlockCSR.from_scipy_csr(
        misaligned_copy(csr.data.astype("float32")),
        misaligned_copy(csr.indices.astype("int32")),
        misaligned_copy(csr.indptr.astype("int32")),
        *csr.shape,
    )
    mat_adopt, _keep2 = make_blockcsr(csr, aligned=True)
    assert mat_copy.is_zero_copy is False
    assert mat_adopt.is_zero_copy is True

    idx_a, sim_a, _ = scrna.build_knn_graph(mat_copy, 6).to_numpy()
    idx_b, sim_b, _ = scrna.build_knn_graph(mat_adopt, 6).to_numpy()
    np.testing.assert_array_equal(idx_a, idx_b)
    np.testing.assert_allclose(sim_a, sim_b, atol=1e-6)


def test_knn_matches_numpy_reference():
    csr = random_csr(n_rows=50, n_cols=64, density=0.2, seed=11)
    mat, _keep = make_blockcsr(csr, aligned=True)
    k = 5
    idx, sim, _ = scrna.build_knn_graph(mat, k).to_numpy()
    idx = idx.reshape(csr.shape[0], k)
    sim = sim.reshape(csr.shape[0], k)

    ref_idx, ref_sim = reference_knn(csr, k)
    # Compare similarity values rather than raw index order, so exact ties do not
    # make this brittle.
    np.testing.assert_allclose(sim, ref_sim, atol=1e-4)
    # The top neighbour should agree outright.
    assert (idx[:, 0] == ref_idx[:, 0]).mean() > 0.95


def test_to_numpy_export_survives_graph_deletion():
    csr = random_csr(n_rows=30, n_cols=40, seed=13)
    mat, _keep = make_blockcsr(csr, aligned=True)
    graph = scrna.build_knn_graph(mat, 4)
    idx, sim, _ = graph.to_numpy()
    before = idx.copy()

    del graph
    gc.collect()
    # The capsule keeps the KnnGraph alive behind the arrays.
    np.testing.assert_array_equal(idx, before)
    assert np.isfinite(sim).all()


def test_build_knn_graph_releases_the_gil():
    """Another Python thread must make progress while the C++ sweep runs."""
    # Big enough that the O(n^2) sweep takes clearly longer than a scheduler
    # quantum -- see below for why size, not just duration, is what matters.
    csr = random_csr(n_rows=2400, n_cols=400, density=0.2, seed=17)
    mat, _keep = make_blockcsr(csr, aligned=True)

    ticks = 0
    stop = False

    # A tight loop, not sleep(0.001)-paced: the previous version measured
    # ticks accumulated only during a 0.02s window straddling the call, so on
    # a fast CI runner the whole build_knn_graph call could complete inside
    # one OS scheduler quantum before the spin thread was ever rescheduled,
    # failing the assertion despite the GIL actually being released
    # (measured flaky on a macOS runner: 0 ticks in that window). Counting
    # from thread start with no sleep removes both the baseline-window race
    # and the dependency on scheduler tick granularity -- it only needs the
    # spin thread to get ANY CPU time during the call, which a multi-core
    # runner gives it as soon as the GIL is actually released.
    def spin():
        nonlocal ticks
        while not stop:
            ticks += 1

    t = threading.Thread(target=spin, daemon=True)
    t.start()

    scrna.build_knn_graph(mat, 10)

    during = ticks
    stop = True
    t.join(timeout=1.0)

    # With the GIL held across the call this counter cannot advance at all.
    assert during > 0, "GIL was not released during build_knn_graph"


requires_hnsw = pytest.mark.skipif(
    not scrna.has_hnsw(), reason="built without -DSCRNA_ENABLE_HNSW=ON"
)


@requires_hnsw
def test_hnsw_recall_against_exact():
    csr = random_csr(n_rows=400, n_cols=64, density=0.25, seed=41)
    mat, _keep = make_blockcsr(csr, aligned=True)
    k = 10

    exact, _, _ = scrna.build_knn_graph(mat, k, method="brute").to_numpy()
    approx, _, _ = scrna.build_knn_graph(
        mat, k, method="hnsw", ef_construction=400, ef_search=200
    ).to_numpy()

    exact = exact.reshape(-1, k)
    approx = approx.reshape(-1, k)

    hits = sum(len(set(a) & set(e)) for a, e in zip(approx, exact))
    recall = hits / exact.size
    assert recall > 0.90, f"HNSW recall {recall:.3f} too low"


@requires_hnsw
def test_hnsw_never_returns_self():
    csr = random_csr(n_rows=200, n_cols=50, seed=43)
    mat, _keep = make_blockcsr(csr, aligned=True)
    k = 8
    idx, _, _ = scrna.build_knn_graph(mat, k, method="hnsw").to_numpy()
    idx = idx.reshape(-1, k)
    for i, row in enumerate(idx):
        assert i not in row, f"row {i} listed itself as its own neighbour"


@requires_hnsw
def test_hnsw_memory_guard_refuses_oversized_index():
    csr = random_csr(n_rows=200, n_cols=50, seed=45)
    mat, _keep = make_blockcsr(csr, aligned=True)
    with pytest.raises(Exception, match="(?i)mb|limit|index"):
        scrna.build_knn_graph(mat, 5, method="hnsw", max_index_bytes=1024)


@requires_hnsw
def test_hnsw_index_bytes_is_dense_in_dimension():
    """Sparsity does not reduce hnswlib's footprint; dimensionality drives it."""
    assert scrna.hnsw_index_bytes(1000, 2000) > 10 * scrna.hnsw_index_bytes(1000, 50)
    # Sanity-check the documented 1.3M-cell figures.
    assert scrna.hnsw_index_bytes(1_300_000, 50) < 1 * 2**30
    assert scrna.hnsw_index_bytes(1_300_000, 2000) > 8 * 2**30


@requires_hnsw
def test_hnsw_dense_path_matches_sparse_path():
    csr = random_csr(n_rows=200, n_cols=32, density=0.4, seed=47)
    mat, _keep = make_blockcsr(csr, aligned=True)
    dense = np.asarray(csr.todense(), dtype=np.float32)
    k = 6
    kw = dict(method="hnsw", ef_construction=300, ef_search=150, num_threads=1)

    a, _, _ = scrna.build_knn_graph(mat, k, **kw).to_numpy()
    b, _, _ = scrna.build_knn_graph_dense(dense, k, **kw).to_numpy()
    np.testing.assert_array_equal(a, b)


@requires_hnsw
def test_hnsw_rejects_bad_method_and_shape():
    csr = random_csr(n_rows=50, n_cols=20, seed=49)
    mat, _keep = make_blockcsr(csr, aligned=True)
    with pytest.raises(Exception):
        scrna.build_knn_graph(mat, 5, method="nonsense")
    with pytest.raises(Exception):
        scrna.build_knn_graph_dense(np.zeros(10, dtype=np.float32), 5)  # 1-D


@requires_hnsw
def test_hnsw_releases_the_gil():
    csr = random_csr(n_rows=3000, n_cols=200, density=0.2, seed=51)
    mat, _keep = make_blockcsr(csr, aligned=True)

    ticks = 0
    stop = False

    def spin():
        nonlocal ticks
        while not stop:
            ticks += 1
            time.sleep(0.001)

    t = threading.Thread(target=spin, daemon=True)
    t.start()
    time.sleep(0.02)
    baseline = ticks
    scrna.build_knn_graph(mat, 10, method="hnsw", ef_construction=200)
    during = ticks - baseline
    stop = True
    t.join(timeout=1.0)
    assert during > 0, "GIL was not released during HNSW build"


def test_rejects_mismatched_buffer_lengths():
    """Regression: short indices/indptr caused an out-of-bounds read (ASan-confirmed)."""
    data = scrna.to_aligned(np.ones(2000, dtype="float32"))
    short_idx = scrna.to_aligned(np.zeros(4, dtype="int32"))
    indptr = scrna.to_aligned(np.array([0, 2000], dtype="int32"))
    with pytest.raises(Exception, match="(?i)indices"):
        scrna.BlockCSR.from_scipy_csr(data, short_idx, indptr, 1, 10)

    good_idx = scrna.to_aligned(np.zeros(2000, dtype="int32"))
    short_indptr = scrna.to_aligned(np.array([0, 2000], dtype="int32"))
    with pytest.raises(Exception, match="(?i)indptr"):
        # n_rows=5 needs 6 indptr entries, only 2 supplied.
        scrna.BlockCSR.from_scipy_csr(data, good_idx, short_indptr, 5, 10)


def test_rejects_non_finite_values():
    """NaN/Inf reaching the top-k comparator is UB, not just a bad answer.

    Updated 2026-09-11 (SM11): this used to reach BlockCSR construction through
    to_aligned, because to_aligned checked float64 for finiteness and not
    float32. It now refuses one step earlier, so the rejection is asserted at
    both boundaries -- the aligned path and the copy path -- rather than only at
    whichever one happened to catch it first.
    """
    for bad in (np.nan, np.inf, -np.inf):
        vals = np.ones(8, dtype="float32")
        vals[3] = bad

        # The aligned path: refused when the buffer is made.
        with pytest.raises(Exception, match="(?i)finite|nan|inf"):
            scrna.to_aligned(vals)

        # The copy path: plain NumPy arrays, refused at construction.
        idx = np.arange(8, dtype="int32")
        indptr = np.array([0, 8], dtype="int32")
        with pytest.raises(Exception, match="(?i)finite|nan|inf"):
            scrna.BlockCSR.from_scipy_csr(vals, idx, indptr, 1, 8)


def test_to_aligned_refuses_lossy_narrowing():
    """int64->int32 used to wrap silently; float64->float32 used to become inf."""
    with pytest.raises(Exception, match="(?i)int32|64-bit"):
        scrna.to_aligned(np.array([0, 3_000_000_000], dtype="int64"))
    with pytest.raises(Exception, match="(?i)overflow|float32"):
        scrna.to_aligned(np.array([1e300], dtype="float64"))
    with pytest.raises(Exception, match="(?i)finite"):
        scrna.to_aligned(np.array([np.nan], dtype="float64"))
    # In-range narrowing still works.
    ok = scrna.to_aligned(np.array([1, 2, 3], dtype="int64"))
    assert ok.dtype == np.int32 and list(ok) == [1, 2, 3]


@requires_hnsw
def test_hnsw_index_bytes_saturates_instead_of_wrapping():
    """The estimate gates an allocation; wrapping would defeat the guard."""
    huge = scrna.hnsw_index_bytes(2**60, 1024)
    assert huge == 2**64 - 1, "overflow must saturate, not wrap to a small value"
    # And a saturated estimate must actually be refused.
    csr = random_csr(n_rows=50, n_cols=20, seed=61)
    mat, _keep = make_blockcsr(csr, aligned=True)
    with pytest.raises(Exception):
        scrna.build_knn_graph(mat, 5, method="hnsw", max_index_bytes=1)


@requires_hnsw
def test_hnsw_dense_rejects_non_finite():
    bad = np.ones((50, 8), dtype="float32")
    bad[7, 2] = np.nan
    with pytest.raises(Exception, match="(?i)finite"):
        scrna.build_knn_graph_dense(bad, 5)


def test_rejects_out_of_range_column_index():
    """Structural validation must reject indices the gather kernels would read OOB."""
    data = scrna.to_aligned(np.array([1.0, 2.0], dtype="float32"))
    indices = scrna.to_aligned(np.array([0, 999], dtype="int32"))  # 999 >= n_cols
    indptr = scrna.to_aligned(np.array([0, 2], dtype="int32"))
    with pytest.raises(Exception):
        scrna.BlockCSR.from_scipy_csr(data, indices, indptr, 1, 10)


# ---------------------------------------------------------------------------
# HnswIndex: build once, persist, query later
# ---------------------------------------------------------------------------


def _atlas(n=300, dim=16, seed=7):
    rng = np.random.default_rng(seed)
    return rng.standard_normal((n, dim)).astype(np.float32)


@requires_hnsw
def test_from_scipy_csr_refuses_lossy_narrowing():
    """SM7: the copy path forcecast, which converts without reporting.

    Measured on the unpatched module: an int64 column index of 2**32 + 5 -- a
    legal index in a matrix 2**33 columns wide -- wrapped to 5. That is an
    IN-RANGE column, so the matrix was accepted and every similarity computed
    from it used the wrong gene. to_aligned() refused the identical value.
    """
    with pytest.raises(Exception, match="(?i)int32|64-bit"):
        scrna.BlockCSR.from_scipy_csr(
            np.array([1.0], dtype="float32"),
            np.array([2**32 + 5], dtype="int64"),
            np.array([0, 1], dtype="int64"),
            1, 2**33,
        )
    # A float index array truncated 1.9 to column 1 the same way.
    with pytest.raises(Exception, match="(?i)int32 or int64"):
        scrna.BlockCSR.from_scipy_csr(
            np.array([1.0, 2.0], dtype="float32"),
            np.array([0.9, 1.9], dtype="float64"),
            np.array([0, 2], dtype="int32"),
            1, 4,
        )
    with pytest.raises(Exception, match="(?i)overflow|float32"):
        scrna.BlockCSR.from_scipy_csr(
            np.array([1e300, 1.0], dtype="float64"),
            np.array([0, 1], dtype="int32"),
            np.array([0, 2], dtype="int32"),
            1, 4,
        )
    # In-range float64/int64 still works -- this is the common SciPy shape.
    csr = random_csr(seed=77)
    ok = scrna.BlockCSR.from_scipy_csr(
        csr.data.astype("float64"), csr.indices.astype("int64"),
        csr.indptr.astype("int64"), csr.shape[0], csr.shape[1],
    )
    assert ok.nnz == csr.nnz
    # So does an integer count matrix, which cannot overflow float32.
    counts = random_csr(seed=78)
    counts.data = np.ceil(counts.data * 100).astype("int64")
    assert scrna.BlockCSR.from_scipy_csr(
        counts.data, counts.indices.astype("int32"), counts.indptr.astype("int32"),
        counts.shape[0], counts.shape[1],
    ).nnz == counts.nnz


@requires_hnsw
def test_a_k_that_overflows_the_result_is_refused():
    """SM1: n * k sizes the neighbour table and offsets every write into it.

    Brute force has checked that product since REVIEW_2026-08-15 finding 4. The
    HNSW producers did not, and k = 2**62 over 4 rows wrapped to 0, then stored
    through a null pointer (UBSan).
    """
    huge_k = 2**62
    dense = np.ascontiguousarray(np.random.default_rng(5).normal(size=(4, 4)), dtype="float32")
    csr = random_csr(n_rows=6, n_cols=8, seed=5)
    mat, _keep = make_blockcsr(csr, aligned=True)

    with pytest.raises(Exception, match="(?i)overflow"):
        scrna.build_knn_graph_dense(dense, huge_k)
    with pytest.raises(Exception, match="(?i)overflow"):
        scrna.build_knn_graph(mat, huge_k, method="hnsw")
    with pytest.raises(Exception, match="(?i)overflow"):
        scrna.build_knn_graph(mat, huge_k, method="brute")

    idx = scrna.HnswIndex(dim=4, capacity=4)
    idx.add_dense(dense)
    with pytest.raises(Exception, match="(?i)overflow"):
        idx.query_dense(dense, huge_k)


@requires_hnsw
def test_concurrent_add_dense_keeps_every_row():
    """SM2: add_dense reads size() as its first label and the binding releases
    the GIL, so two concurrent calls took the same base -- and hnswlib's
    addPoint UPDATES an existing label rather than failing. 400 + 400 rows left
    an index of 400, silently."""
    rng = np.random.default_rng(11)
    a = np.ascontiguousarray(rng.normal(size=(400, 8)), dtype="float32")
    b = np.ascontiguousarray(rng.normal(size=(400, 8)), dtype="float32")
    idx = scrna.HnswIndex(dim=8, capacity=800)

    errors = []

    def add(block):
        try:
            idx.add_dense(block)
        except Exception as exc:  # noqa: BLE001
            errors.append(exc)

    threads = [threading.Thread(target=add, args=(x,)) for x in (a, b)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()

    assert not errors
    assert len(idx) == 800


@requires_hnsw
def test_hnsw_index_refuses_M_outside_the_supported_range():
    """SM3: hnswlib clamps M > 10000 silently and does not validate 0 or 1.

    M=1 makes (int)inf undefined in getRandomLevel; M=0 walks a vector at
    offset -4; M=20000 built and saved an index that then refused to load,
    because our header kept 20000 and the payload recorded the clamped 10000.
    """
    for bad in (0, 1, 20000):
        with pytest.raises(Exception, match="(?i)outside"):
            scrna.HnswIndex(dim=4, capacity=10, M=bad, max_index_bytes=2**40)
    # The one-shot entry points take the same path.
    dense = np.ascontiguousarray(np.random.default_rng(6).normal(size=(20, 4)), dtype="float32")
    with pytest.raises(Exception, match="(?i)outside"):
        scrna.build_knn_graph_dense(dense, 3, M=1)
    assert scrna.HnswIndex(dim=4, capacity=10, M=16) is not None


@requires_hnsw
def test_save_refuses_labels_that_no_longer_match(tmp_path):
    """SM5: labels set, then more elements added, produced a meta describing a
    prefix of the index -- save succeeded and load_index refused the pair, long
    after the build that produced it was gone."""
    rng = np.random.default_rng(12)
    data = np.ascontiguousarray(rng.normal(size=(20, 8)), dtype="float32")
    idx = scrna.HnswIndex(dim=8, capacity=20)
    idx.add_dense(data[:10])
    idx.set_labels([f"cell{i}" for i in range(10)])
    idx.add_dense(data[10:])

    path = str(tmp_path / "stale.hnsw")
    with pytest.raises(Exception, match="(?i)label"):
        idx.save_index(path)
    assert not (tmp_path / "stale.hnsw").exists()

    idx.set_labels([f"cell{i}" for i in range(20)])
    idx.save_index(path)
    assert scrna.HnswIndex.load_index(path).labels == [f"cell{i}" for i in range(20)]


def test_to_aligned_refuses_non_finite_float32():
    """SM11: float64 input was checked for NaN/Inf and float32 was not.

    Block_CSR::validate() caught them either way, so nothing wrong was ever
    computed -- but to_aligned is documented as the checked path, and the
    rejection belongs where the reported index is still the caller's own.
    """
    for bad in (np.nan, np.inf, -np.inf):
        with pytest.raises(Exception, match="(?i)finite"):
            scrna.to_aligned(np.array([1.0, bad], dtype="float32"))
        with pytest.raises(Exception, match="(?i)finite"):
            scrna.to_aligned(np.array([1.0, bad], dtype="float64"))
    assert list(scrna.to_aligned(np.array([1.0, 2.0], dtype="float32"))) == [1.0, 2.0]


@pytest.mark.parametrize(
    "name, data, indices, indptr, rows, cols, match",
    [
        ("indptr decreasing",   [1.0, 2.0], [0, 1], [0, 2, 1], 2, 4, "(?i)non-decreasing"),
        ("indptr tail != nnz",  [1.0, 2.0], [0, 1], [0, 1],    1, 4, "(?i)must equal nnz"),
        ("indptr[0] != 0",      [1.0],      [0],    [1, 1],    1, 4, "(?i)must be 0"),
        ("index >= n_cols",     [1.0],      [7],    [0, 1],    1, 2, "(?i)out of range"),
        ("negative index",      [1.0],      [-1],   [0, 1],    1, 4, "(?i)out of range"),
        ("unsorted indices",    [1.0, 2.0], [2, 1], [0, 2],    1, 4, "(?i)increasing"),
        ("duplicate index",     [1.0, 2.0], [1, 1], [0, 2],    1, 4, "(?i)duplicate"),
        ("nan value",           [1.0, np.nan], [0, 1], [0, 2], 1, 4, "(?i)finite"),
        ("inf value",           [1.0, np.inf], [0, 1], [0, 2], 1, 4, "(?i)finite"),
    ],
)
def test_malformed_csr_structures_fail_closed(name, data, indices, indptr, rows, cols, match):
    """Every structural invariant, exercised through the boundary a user calls.

    The invariants live in Block_CSR::validate() and are covered in C++ by
    test_block_csr_validation. Nothing proved they survive the pybind11 layer,
    which is where real input arrives -- and the gather kernels index a dense
    scratch buffer with these values using unchecked hardware gathers, so a
    structural violation that got through would be an out-of-bounds read at
    vector speed rather than a wrong number.
    """
    with pytest.raises(Exception, match=match):
        scrna.BlockCSR.from_scipy_csr(
            np.array(data, dtype="float32"),
            np.array(indices, dtype="int32"),
            np.array(indptr, dtype="int32"),
            rows, cols,
        )


def test_an_empty_cell_scores_zero_not_nan():
    """A cell with no counts has norm 0. The result must be a defined zero, not
    a NaN reaching the top-k comparator (where it would be UB, not a bad score).
    """
    dense = np.array([[1.0, 2.0, 0.0], [0.0, 0.0, 0.0], [2.0, 1.0, 1.0]], dtype="float32")
    m = sp.csr_matrix(dense)
    mat = scrna.BlockCSR.from_scipy_csr(
        m.data, m.indices.astype("int32"), m.indptr.astype("int32"), 3, 3)
    _idx, sim, k = scrna.build_knn_graph(mat, 2, method="brute").to_numpy()
    assert np.isfinite(sim).all()
    assert (sim.reshape(3, k)[1] == 0.0).all()   # the empty cell matches nothing
    assert sim.reshape(3, k)[0][0] > 0.0         # the real cells still match


@requires_hnsw
def test_hnsw_index_refuses_a_capacity_past_int32(tmp_path):
    """SM10: neighbour ids are int32, so a larger index could answer with a
    wrapped, negative id -- the same value as the 'no neighbour' sentinel."""
    with pytest.raises(Exception, match="(?i)capacity"):
        scrna.HnswIndex(dim=1, capacity=2**31, max_index_bytes=2**62)


@requires_hnsw
def test_hnsw_k_at_the_edges_of_the_index():
    """k is a row stride, not a promise of k answers: k=0 is an empty graph and
    k past the index fills the remainder with the -1 sentinel."""
    rng = np.random.default_rng(21)
    data = np.ascontiguousarray(rng.normal(size=(8, 4)), dtype="float32")
    idx = scrna.HnswIndex(dim=4, capacity=8)
    idx.add_dense(data)

    empty_idx, _empty_sim, k0 = idx.query_dense(data, 0).to_numpy()
    assert k0 == 0 and empty_idx.size == 0

    wide, sims, k = idx.query_dense(data[:1], 50).to_numpy()
    assert k == 50 and wide.size == 50
    assert int((wide >= 0).sum()) == 8      # foreign query: every element is eligible
    assert wide[-1] == -1 and sims[-1] == 0.0


@requires_hnsw
def test_query_and_add_require_the_index_dimension():
    """Feature width must equal the index's dim at every dense entry point.

    The index stores fixed-width vectors, so a row of the wrong length would be
    read at the index's stride -- past the end of the caller's row.
    """
    rng = np.random.default_rng(31)
    idx = scrna.HnswIndex(dim=8, capacity=20)
    idx.add_dense(np.ascontiguousarray(rng.normal(size=(5, 8)), dtype="float32"))

    for bad_dim in (4, 16):
        wrong = np.ascontiguousarray(rng.normal(size=(3, bad_dim)), dtype="float32")
        with pytest.raises(Exception, match="(?i)columns but the index has dim"):
            idx.add_dense(wrong)
        with pytest.raises(Exception, match="(?i)columns but the index has dim"):
            idx.query_dense(wrong, 2)

    # A 1-D array is not an (n_obs, n_dims) matrix at all.
    with pytest.raises(Exception, match="(?i)2-D"):
        idx.query_dense(np.zeros(8, dtype="float32"), 2)

    # The CSR entry points hold the same rule against n_cols.
    csr = random_csr(n_rows=4, n_cols=4, seed=31)
    mat, _keep = make_blockcsr(csr, aligned=True)
    with pytest.raises(Exception, match="(?i)columns but the index has dim"):
        idx.add_csr(mat)
    with pytest.raises(Exception, match="(?i)columns but the index has dim"):
        idx.query_csr(mat, 2)


def test_empty_and_zero_dimension_matrices_are_defined():
    """Zero rows, zero columns and zero nonzeros are all legal shapes.

    They must build and search without allocating from a null pointer or
    dividing by a zero dimension -- defined and empty, not refused and not UB.
    """
    empty_data = np.array([], dtype="float32")
    empty_idx = np.array([], dtype="int32")

    zero_by_zero = scrna.BlockCSR.from_scipy_csr(
        empty_data, empty_idx, np.array([0], dtype="int32"), 0, 0)
    assert zero_by_zero.shape == (0, 0) and zero_by_zero.nnz == 0
    assert scrna.build_knn_graph(zero_by_zero, 5).to_numpy()[0].size == 0

    no_rows = scrna.BlockCSR.from_scipy_csr(
        empty_data, empty_idx, np.array([0], dtype="int32"), 0, 8)
    assert no_rows.shape == (0, 8)

    # Three cells, no counts at all: every row empty, so nothing matches
    # anything and every similarity is a defined zero.
    no_counts = scrna.BlockCSR.from_scipy_csr(
        empty_data, empty_idx, np.array([0, 0, 0, 0], dtype="int32"), 3, 4)
    nbr, sim, k = scrna.build_knn_graph(no_counts, 2).to_numpy()
    assert k == 2 and nbr.size == 6
    assert np.isfinite(sim).all() and (sim == 0.0).all()


@requires_hnsw
def test_hnsw_refuses_zero_dim_and_zero_capacity():
    with pytest.raises(Exception, match="(?i)dim must be non-zero"):
        scrna.HnswIndex(dim=0, capacity=10)
    with pytest.raises(Exception, match="(?i)capacity must be non-zero"):
        scrna.HnswIndex(dim=8, capacity=0)
    # An empty add is a no-op rather than an error.
    idx = scrna.HnswIndex(dim=8, capacity=10)
    idx.add_dense(np.zeros((0, 8), dtype="float32"))
    assert len(idx) == 0


def test_hnsw_index_round_trip_is_lossless(tmp_path):
    ref, qry = _atlas(), _atlas(n=8, seed=99)
    idx = scrna.HnswIndex(dim=16, capacity=300, num_threads=1)
    idx.add_dense(ref, num_threads=1)
    assert len(idx) == 300
    assert idx.dim == 16

    before = idx.query_dense(qry, k=5, num_threads=1).to_numpy()
    path = str(tmp_path / "atlas.hnsw")
    idx.save_index(path)

    loaded = scrna.HnswIndex.load_index(path)
    after = loaded.query_dense(qry, k=5, num_threads=1).to_numpy()

    # Exact equality on both arrays. A round trip through disk must not move a
    # single result; a tolerance here would hide a truncated payload.
    assert np.array_equal(before[0], after[0])
    assert np.array_equal(before[1], after[1])


@requires_hnsw
def test_hnsw_index_writes_both_files(tmp_path):
    idx = scrna.HnswIndex(dim=8, capacity=20, num_threads=1)
    idx.add_dense(_atlas(20, 8), num_threads=1)
    path = tmp_path / "a.hnsw"
    idx.save_index(str(path))
    # The sidecar is not optional: it carries the dimension and metric that make
    # the payload safe to load.
    assert path.exists()
    assert (tmp_path / "a.hnsw.meta").exists()


@requires_hnsw
def test_hnsw_index_query_does_not_exclude_self(tmp_path):
    """A foreign query identical to an indexed vector must find that vector.

    build_knn_graph_dense() drops self-matches because its query rows *are* the
    index. For a reference atlas that behaviour would silently delete a real
    neighbour whose label happened to equal the query row's position.
    """
    ref = _atlas(100, 12)
    idx = scrna.HnswIndex(dim=12, capacity=100, num_threads=1)
    idx.add_dense(ref, num_threads=1)

    # to_numpy() returns flat (neighbour_idx, neighbour_sim, k); reshape to rows.
    nbr, sim, k = idx.query_dense(ref[:3], k=3, num_threads=1).to_numpy()
    nbr = nbr.reshape(-1, k)
    sim = sim.reshape(-1, k)
    assert list(nbr[:, 0]) == [0, 1, 2]
    assert np.all(sim[:, 0] > 0.999)


@requires_hnsw
def test_hnsw_index_labels_round_trip(tmp_path):
    names = [f"AAACCTGAGC-{i}" for i in range(20)]
    idx = scrna.HnswIndex(dim=8, capacity=20, num_threads=1)
    idx.add_dense(_atlas(20, 8), num_threads=1)
    idx.set_labels(names)
    path = str(tmp_path / "labelled.hnsw")
    idx.save_index(path)

    assert scrna.HnswIndex.load_index(path).labels == names


@requires_hnsw
def test_hnsw_index_label_count_must_match(tmp_path):
    idx = scrna.HnswIndex(dim=8, capacity=20, num_threads=1)
    idx.add_dense(_atlas(20, 8), num_threads=1)
    with pytest.raises(ValueError):
        idx.set_labels(["only-one"])


@requires_hnsw
def test_hnsw_index_rejects_dim_mismatch():
    idx = scrna.HnswIndex(dim=16, capacity=50, num_threads=1)
    with pytest.raises(ValueError):
        idx.add_dense(_atlas(10, 8), num_threads=1)
    with pytest.raises(ValueError):
        idx.query_dense(_atlas(10, 8), k=3, num_threads=1)


@requires_hnsw
def test_hnsw_index_rejects_overfull_add():
    idx = scrna.HnswIndex(dim=8, capacity=10, num_threads=1)
    with pytest.raises(ValueError):
        idx.add_dense(_atlas(11, 8), num_threads=1)
    # Refused up front, so the index is untouched rather than half-populated.
    assert len(idx) == 0


@requires_hnsw
def test_hnsw_index_refuses_missing_or_corrupt_files(tmp_path):
    idx = scrna.HnswIndex(dim=8, capacity=30, num_threads=1)
    idx.add_dense(_atlas(30, 8), num_threads=1)
    path = tmp_path / "x.hnsw"
    idx.save_index(str(path))

    with pytest.raises(RuntimeError):
        scrna.HnswIndex.load_index(str(tmp_path / "does-not-exist.hnsw"))

    # Payload without its sidecar.
    (tmp_path / "x.hnsw.meta").unlink()
    with pytest.raises(RuntimeError):
        scrna.HnswIndex.load_index(str(path))


@requires_hnsw
def test_hnsw_index_refuses_truncated_payload(tmp_path):
    idx = scrna.HnswIndex(dim=8, capacity=30, num_threads=1)
    idx.add_dense(_atlas(30, 8), num_threads=1)
    path = tmp_path / "t.hnsw"
    idx.save_index(str(path))

    # What an interrupted save on a full disk actually leaves behind.
    body = path.read_bytes()
    path.write_bytes(body[: len(body) // 2])
    with pytest.raises(RuntimeError):
        scrna.HnswIndex.load_index(str(path))


@requires_hnsw
def test_hnsw_index_load_releases_the_gil(tmp_path):
    idx = scrna.HnswIndex(dim=32, capacity=4000, num_threads=1)
    idx.add_dense(_atlas(4000, 32), num_threads=1)
    path = str(tmp_path / "gil.hnsw")
    idx.save_index(path)

    ticks = [0]
    stop = threading.Event()

    def spin():
        while not stop.is_set():
            ticks[0] += 1
            time.sleep(0.001)

    t = threading.Thread(target=spin)
    t.start()
    try:
        scrna.HnswIndex.load_index(path)
    finally:
        stop.set()
        t.join()
    # The checksum pass is a full sequential read; if it held the GIL the
    # counter thread would not have advanced.
    assert ticks[0] > 0
