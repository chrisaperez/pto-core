#!/usr/bin/env bash
#
# pto-core end-to-end demo.
#
# Drives all three modules in one pass:
#
#   1. fastq_stream     synthetic FASTQ -> QC + adapter/quality trimming
#   2. scrna_matrix     sparse expression matrix -> cosine k-NN graph
#   3. cuttag_profiler  indexed BAM + promoter BED -> signal matrix + dashboard
#
# DESIGN NOTES
#
# Offline by default. Every fixture is either checked into the repo or
# synthesised locally, so the demo runs on an air-gapped machine -- which is the
# deployment these tools are built for. The one optional network step (scanpy's
# pbmc3k, ~5 MB) is used when available and cleanly substituted when not; pass
# --offline to skip the attempt entirely.
#
# Idempotent. Re-running reuses an existing Release build and existing fixtures
# unless --rebuild / --refresh-fixtures is given.
#
# Usage:
#   scripts/run_demo.sh                 # full demo, then serve the dashboard
#   scripts/run_demo.sh --no-serve      # pipeline only, no HTTP server
#   scripts/run_demo.sh --port 9000     # dashboard on another port
#   scripts/run_demo.sh --help
set -Eeuo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
readonly REPO_ROOT
BUILD_DIR="${BUILD_DIR:-${REPO_ROOT}/build}"
WORK_DIR="${WORK_DIR:-${REPO_ROOT}/demo-out}"
PORT=8080
DO_BUILD=1
DO_SERVE=1
REBUILD=0
REFRESH=0
OFFLINE=0
THREADS="$( (command -v nproc >/dev/null && nproc) || sysctl -n hw.ncpu 2>/dev/null || echo 4)"

c_blue=$'\033[1;34m'; c_green=$'\033[1;32m'; c_yellow=$'\033[1;33m'
c_red=$'\033[1;31m';  c_dim=$'\033[2m';      c_off=$'\033[0m'

step()  { printf '\n%s==>%s %s\n' "$c_blue"  "$c_off" "$*"; }
ok()    { printf '%s  ok%s  %s\n'  "$c_green" "$c_off" "$*"; }
note()  { printf '%s      %s%s\n'  "$c_dim"   "$*" "$c_off"; }
warn()  { printf '%s[warn]%s %s\n' "$c_yellow" "$c_off" "$*" >&2; }
die()   { printf '%s[fail]%s %s\n' "$c_red"   "$c_off" "$*" >&2; exit 1; }

trap 'die "aborted at line ${LINENO}"' ERR

usage() {
    sed -n '2,30p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
    cat <<'EOF'

OPTIONS
  --port N              dashboard port (default 8080)
  --build-dir DIR       CMake build directory (default ./build)
  --work-dir DIR        demo outputs (default ./demo-out)
  --threads N           worker threads (default: all cores)
  --skip-build          assume the build directory is already usable
  --rebuild             reconfigure and rebuild from scratch
  --refresh-fixtures    regenerate synthetic fixtures even if present
  --no-serve            run the pipeline, then stop (no HTTP server)
  --offline             never attempt a network download
  -h, --help            this message
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --port)             PORT="$2"; shift 2 ;;
        --build-dir)        BUILD_DIR="$2"; shift 2 ;;
        --work-dir)         WORK_DIR="$2"; shift 2 ;;
        --threads)          THREADS="$2"; shift 2 ;;
        --skip-build)       DO_BUILD=0; shift ;;
        --rebuild)          REBUILD=1; shift ;;
        --refresh-fixtures) REFRESH=1; shift ;;
        --no-serve)         DO_SERVE=0; shift ;;
        --offline)          OFFLINE=1; shift ;;
        -h|--help)          usage; exit 0 ;;
        *)                  die "unknown option: $1 (try --help)" ;;
    esac
done

FIXTURES="${WORK_DIR}/fixtures"
mkdir -p "${WORK_DIR}" "${FIXTURES}"

PYTHON="${PYTHON:-python3}"
command -v "${PYTHON}" >/dev/null 2>&1 || die "python3 not found; needed to synthesise fixtures"

# ---------------------------------------------------------------------------
# 0. Build
# ---------------------------------------------------------------------------
fastq_bin="${BUILD_DIR}/modules/fastq_stream/fastq_stream"
cuttag_bin="${BUILD_DIR}/modules/cuttag_profiler/cuttag_profiler"

if [[ ${REBUILD} -eq 1 ]]; then
    step "Rebuilding from scratch"
    rm -rf "${BUILD_DIR}"
elif [[ ${DO_BUILD} -eq 1 && -x "${fastq_bin}" && -x "${cuttag_bin}" ]]; then
    step "Build artifacts present"
    note "${BUILD_DIR}"
    note "pass --rebuild to force a clean reconfigure"
    DO_BUILD=0
fi

if [[ ${DO_BUILD} -eq 1 ]]; then
    step "Configuring and building (Release)"
    command -v cmake >/dev/null 2>&1 || die "cmake not found"

    # OpenMP is mandatory for scrna_matrix and is the single most common reason
    # a first build fails on macOS. Surface the fix here rather than letting the
    # user read it out of a CMake stack trace.
    if ! cmake -S "${REPO_ROOT}" -B "${BUILD_DIR}" \
            -DCMAKE_BUILD_TYPE=Release \
            >"${WORK_DIR}/configure.log" 2>&1; then
        tail -n 40 "${WORK_DIR}/configure.log" >&2
        echo >&2
        die "configure failed (full log: ${WORK_DIR}/configure.log).
    If the error mentions OpenMP:
        macOS  : brew install libomp
        Ubuntu : sudo apt-get install libomp-dev
    To build only the modules that do not need it:
        cmake -S . -B build -DPTO_BUILD_SCRNA_MATRIX=OFF"
    fi

    cmake --build "${BUILD_DIR}" -j "${THREADS}" >"${WORK_DIR}/build.log" 2>&1 \
        || { tail -n 40 "${WORK_DIR}/build.log" >&2; die "build failed (log: ${WORK_DIR}/build.log)"; }
    ok "built into ${BUILD_DIR}"
fi

[[ -x "${fastq_bin}" ]]  || die "fastq_stream binary not found at ${fastq_bin}"
[[ -x "${cuttag_bin}" ]] || die "cuttag_profiler binary not found at ${cuttag_bin}"

# ---------------------------------------------------------------------------
# 1. fastq_stream — QC and trimming
# ---------------------------------------------------------------------------
step "Stage 1/3  fastq_stream — streaming QC and trimming"

reads_fq="${FIXTURES}/demo_reads.fq.gz"
if [[ ${REFRESH} -eq 1 || ! -s "${reads_fq}" ]]; then
    note "synthesising 200k reads with adapter read-through and a quality ramp"
    # Synthetic rather than GIAB: scripts/download_giab_data.sh pulls tens of GB,
    # which is right for benchmarking and wrong for a demo. This produces the
    # same *shapes* the pipeline cares about -- adapter read-through on short
    # inserts, 3' quality decay, and some N-rich reads to exercise the filters.
    "${PYTHON}" - "$reads_fq" <<'PY'
import gzip, random, sys

out_path = sys.argv[1]
rng = random.Random(20260815)
TRUSEQ = "AGATCGGAAGAGC"
N_READS, READ_LEN = 200_000, 150

with gzip.open(out_path, "wt") as fh:
    for i in range(N_READS):
        insert = rng.randrange(40, READ_LEN + 40)
        seq = "".join(rng.choice("ACGT") for _ in range(min(insert, READ_LEN)))
        # Short inserts read through into the adapter, which is what the
        # seed-and-extend trimmer exists to find.
        if insert < READ_LEN:
            seq = (seq + TRUSEQ + "".join(rng.choice("ACGT") for _ in range(READ_LEN)))[:READ_LEN]
        # 1 read in 200 is N-heavy, to exercise --max-n-rate.
        if i % 200 == 0:
            seq = "N" * (READ_LEN // 2) + seq[READ_LEN // 2:]
        # Phred decays toward the 3' end, as on a real instrument.
        quals = []
        for pos in range(len(seq)):
            mean = 38 - 22 * (pos / len(seq)) ** 2
            quals.append(chr(33 + max(2, min(41, int(rng.gauss(mean, 3))))))
        fh.write(f"@read{i}\n{seq}\n+\n{''.join(quals)}\n")
PY
fi
ok "input: ${reads_fq} ($(du -h "${reads_fq}" | cut -f1))"

"${fastq_bin}" \
    --in "${reads_fq}" \
    --out "${WORK_DIR}/trimmed.fq" \
    --json "${WORK_DIR}/fastq_qc.json" \
    --threads "${THREADS}" \
    --trim-tail 20 --min-len 30 --min-mean-q 20 \
    2>"${WORK_DIR}/fastq_stream.log" \
    || { tail -n 20 "${WORK_DIR}/fastq_stream.log" >&2; die "fastq_stream failed"; }

"${PYTHON}" - "${WORK_DIR}/fastq_qc.json" <<'PY'
import json, sys
r = json.load(open(sys.argv[1]))
def find(*names):
    for n in names:
        if n in r: return r[n]
    for v in r.values():
        if isinstance(v, dict):
            for n in names:
                if n in v: return v[n]
    return None
kept, total = find("reads_out", "reads_kept"), find("reads_in", "reads_total")
if total: print(f"      reads in {total:,} -> kept {kept:,} ({100*kept/total:.1f}%)")
q = find("mean_q", "mean_quality")
if q is not None: print(f"      mean Phred {q:.1f}")
PY
ok "QC report: ${WORK_DIR}/fastq_qc.json"
ok "trimmed reads: ${WORK_DIR}/trimmed.fq ($(du -h "${WORK_DIR}/trimmed.fq" | cut -f1))"

# ---------------------------------------------------------------------------
# 2. scrna_matrix — k-NN graph
# ---------------------------------------------------------------------------
step "Stage 2/3  scrna_matrix — cosine k-NN graph"

pymod_dir="${BUILD_DIR}/modules/scrna_matrix"
if ! compgen -G "${pymod_dir}/scrna_matrix_py*.so" >/dev/null; then
    warn "Python bindings were not built (pybind11 missing at configure time)."
    warn "Skipping stage 2. Install pybind11 and re-run with --rebuild:"
    warn "    ${PYTHON} -m pip install pybind11 numpy scipy"
else
    KNN_OFFLINE="${OFFLINE}" PYTHONPATH="${pymod_dir}:${REPO_ROOT}/modules/scrna_matrix/examples" \
    "${PYTHON}" - "${WORK_DIR}" <<'PY'
import os, sys, time
import numpy as np
import scrna_matrix_py as scrna

work = sys.argv[1]
offline = os.environ.get("KNN_OFFLINE") == "1"

adata = None
if not offline:
    try:
        import scanpy as sc
        # pbmc3k is ~5 MB and cached by scanpy after the first run.
        adata = sc.datasets.pbmc3k()
        sc.pp.filter_genes(adata, min_cells=3)
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, n_top_genes=1000)
        adata = adata[:, adata.var.highly_variable].copy()
        print(f"      dataset: 10x pbmc3k  {adata.shape[0]:,} cells x {adata.shape[1]:,} genes")
    except Exception as exc:                      # noqa: BLE001
        print(f"      scanpy/pbmc3k unavailable ({type(exc).__name__}); using a synthetic matrix")
        adata = None

if adata is not None:
    from scrna_anndata import blockcsr_from_anndata
    mat, _keepalive = blockcsr_from_anndata(adata)
else:
    import scipy.sparse as sp
    rng = np.random.default_rng(0)
    n_cells, n_genes = 3000, 1000
    X = sp.random(n_cells, n_genes, density=0.06, format="csr",
                  dtype=np.float32, random_state=rng)
    X.data = np.abs(X.data).astype(np.float32) + 0.1
    X.sort_indices()
    print(f"      dataset: synthetic  {n_cells:,} cells x {n_genes:,} genes")
    # to_aligned() returns frozen, capsule-backed buffers -- the only shape the
    # engine will adopt without copying. A writeable array is always copied,
    # because borrowing memory Python can still mutate is not memory-safe.
    data    = scrna.to_aligned(X.data.astype(np.float32))
    indices = scrna.to_aligned(X.indices.astype(np.int32))
    indptr  = scrna.to_aligned(X.indptr.astype(np.int32))
    mat = scrna.BlockCSR.from_scipy_csr(data, indices, indptr, *X.shape)

print(f"      SIMD kernel: {scrna.simd_isa()}   zero-copy: {mat.is_zero_copy}   nnz: {mat.nnz:,}")

K = 15
method = "hnsw" if scrna.has_hnsw() else "brute"
t0 = time.perf_counter()
graph = scrna.build_knn_graph(mat, K, method=method)
elapsed = time.perf_counter() - t0

idx, sim, k = graph.to_numpy()
idx = idx.reshape(-1, k)
sim = sim.reshape(-1, k)
print(f"      method={method}  k={k}  {elapsed:.2f}s  "
      f"({mat.shape[0]/max(elapsed,1e-9):,.0f} cells/s)")
print(f"      mean similarity to nearest neighbour: {sim[:, 0].mean():.4f}")

np.savez_compressed(os.path.join(work, "knn_graph.npz"),
                    neighbor_idx=idx, neighbor_sim=sim, k=k)
PY
    ok "k-NN graph: ${WORK_DIR}/knn_graph.npz"
fi

# ---------------------------------------------------------------------------
# 3. cuttag_profiler — signal matrix
# ---------------------------------------------------------------------------
step "Stage 3/3  cuttag_profiler — CUT&Tag signal matrix"

data_root="${REPO_ROOT}/modules/cuttag_profiler/data"
bam="${data_root}/demo.bam"
bed="${data_root}/demo_promoters.bed"

if [[ ! -s "${bam}" || ! -s "${bed}" ]]; then
    note "checked-in fixtures missing; generating a synthetic dataset"
    command -v samtools >/dev/null 2>&1 \
        || die "demo BAM missing and samtools not on PATH.
    Install samtools, or restore modules/cuttag_profiler/data/demo.bam"
    data_root="${FIXTURES}/cuttag"
    mkdir -p "${data_root}"
    "${PYTHON}" "${REPO_ROOT}/modules/cuttag_profiler/scripts/make_demo_data.py" \
        --out-dir "${data_root}" --regions 2000 --reads 500000 >/dev/null
    bam="${data_root}/demo.bam"
    bed="${data_root}/demo_promoters.bed"
fi
ok "BAM: ${bam} ($(du -h "${bam}" | cut -f1))"

"${cuttag_bin}" profile \
    --bam "${bam}" --regions "${bed}" \
    --upstream 2000 --downstream 2000 --bin-size 50 \
    --normalization CPM --threads "${THREADS}" \
    --out-matrix "${WORK_DIR}/cuttag_matrix.tsv" \
    --out-profile "${WORK_DIR}/cuttag_profile.tsv" \
    2>&1 | sed 's/^/      /'

ok "signal matrix:  ${WORK_DIR}/cuttag_matrix.tsv ($(du -h "${WORK_DIR}/cuttag_matrix.tsv" | cut -f1))"
ok "meta-profile:   ${WORK_DIR}/cuttag_profile.tsv"

# ---------------------------------------------------------------------------
# Summary / dashboard
# ---------------------------------------------------------------------------
step "Demo complete"
printf '      outputs in %s\n' "${WORK_DIR}"
find "${WORK_DIR}" -maxdepth 1 -type f -exec du -h {} + | sed 's|^|      |'

if [[ ${DO_SERVE} -eq 0 ]]; then
    note "--no-serve given; skipping the dashboard"
    exit 0
fi

step "Launching the dashboard on http://127.0.0.1:${PORT}/"
# --data-root is passed explicitly and is NOT optional here.
#
# The binary now requires it (SECURITY_HTTP_2026-08-15 H2 is closed), so this is
# no longer load-bearing -- but it stays explicit because the demo should show
# the confined invocation, not rely on a default.
#
# There is still no authentication of any kind: the Host/Origin checks stop a
# hostile *web page*, but anything that can open a socket sets whatever headers
# it likes. Keep this bound to loopback. Never pass --allow-remote on a shared
# or internet-facing host -- see docs/DEMO.md for the tunnel pattern.
note "confined to --data-root ${data_root}"
note "press Ctrl-C to stop"
exec "${cuttag_bin}" serve \
    --port "${PORT}" \
    --host 127.0.0.1 \
    --data-root "${data_root}"
