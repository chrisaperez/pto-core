#!/usr/bin/env bash
#
# Provision an Ubuntu 24.04 EC2 instance to run cuttag_profiler.
#
# Run ON the instance, as a sudo-capable user:
#
#   git clone <repo> pto-core && cd pto-core
#   sudo ./deploy/provision_aws.sh
#
# or pass it as EC2 user-data (it is idempotent and non-interactive).
#
# WHAT THIS DOES NOT DO, DELIBERATELY
#
# It does not open port 8080 to the internet, and it does not print a security
# group rule that would. cuttag_profiler has NO authentication: the Host/Origin
# checks bind a browser, not anything that can open a socket, and the session
# token (SECURITY_HTTP H1c) is unimplemented.
#
# H2 (arbitrary file read) is closed -- the binary now requires --data-root --
# so an exposed port no longer reads the whole filesystem. It would still expose
# every BAM/BED under the data root, and the profiling API itself, to anyone who
# can reach the port.
#
# The supported access path is an SSH tunnel. The script prints the exact
# command at the end.
set -Eeuo pipefail

readonly APP_USER="${APP_USER:-ubuntu}"
readonly DATA_DIR="${DATA_DIR:-/srv/pto/data}"
readonly REPO_DIR="${REPO_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
readonly PORT="${PORT:-8080}"

log()  { printf '\033[1;34m==>\033[0m %s\n' "$*"; }
warn() { printf '\033[1;33m[warn]\033[0m %s\n' "$*" >&2; }
die()  { printf '\033[1;31m[fail]\033[0m %s\n' "$*" >&2; exit 1; }

trap 'die "aborted at line ${LINENO}"' ERR

[[ ${EUID} -eq 0 ]] || die "run with sudo: sudo $0"

if [[ -r /etc/os-release ]]; then
    # shellcheck disable=SC1091  # provided by the OS, not the repo
    . /etc/os-release
    [[ "${ID:-}" == "ubuntu" ]] || warn "expected Ubuntu; found '${ID:-unknown}' -- continuing"
    [[ "${VERSION_ID:-}" == "24.04" ]] || warn "expected 24.04; found '${VERSION_ID:-unknown}'"
fi

# ---------------------------------------------------------------------------
log "Installing base packages"
export DEBIAN_FRONTEND=noninteractive
apt-get update -qq
apt-get install -y -qq --no-install-recommends \
    ca-certificates curl gnupg git jq ufw

# ---------------------------------------------------------------------------
log "Installing Docker Engine and the compose plugin"
if ! command -v docker >/dev/null 2>&1; then
    install -m 0755 -d /etc/apt/keyrings
    curl -fsSL https://download.docker.com/linux/ubuntu/gpg \
        | gpg --dearmor -o /etc/apt/keyrings/docker.gpg
    chmod a+r /etc/apt/keyrings/docker.gpg
    # shellcheck disable=SC1091  # /etc/os-release is provided by the OS
    codename="$(. /etc/os-release && echo "${VERSION_CODENAME}")"
    cat >/etc/apt/sources.list.d/docker.list <<EOF
deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.gpg] \
https://download.docker.com/linux/ubuntu ${codename} stable
EOF
    apt-get update -qq
    apt-get install -y -qq --no-install-recommends \
        docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin
fi
systemctl enable --now docker
id -nG "${APP_USER}" 2>/dev/null | grep -qw docker || usermod -aG docker "${APP_USER}" || true

# ---------------------------------------------------------------------------
# Toolchain for building outside a container (the demo harness, the other two
# modules). The service image builds its own; this is for `scripts/run_demo.sh`.
log "Installing the build toolchain and Python dependencies"
apt-get install -y -qq --no-install-recommends \
    build-essential cmake ninja-build pkg-config \
    libomp-dev \
    libhts-dev zlib1g-dev libdeflate-dev samtools \
    python3 python3-dev python3-pip python3-venv \
    python3-numpy python3-scipy pybind11-dev

# libomp-dev, specifically: scrna_matrix requires OpenMP and the build hard-fails
# without it. That is deliberate -- without it every `#pragma omp` compiles to
# serial code and the test suite passes while never running a parallel iteration.

# ---------------------------------------------------------------------------
log "Creating the data directory"
mkdir -p "${DATA_DIR}"
chown -R "${APP_USER}:${APP_USER}" "${DATA_DIR}"
chmod 0755 "${DATA_DIR}"

if [[ -d "${REPO_DIR}/modules/cuttag_profiler/data" ]]; then
    cp -n "${REPO_DIR}/modules/cuttag_profiler/data/"demo* "${DATA_DIR}/" 2>/dev/null || true
    log "Seeded ${DATA_DIR} with the checked-in demo fixtures"
fi

# ---------------------------------------------------------------------------
log "Configuring the host firewall"
# Belt to the security group's braces. Note that Docker's DNAT rules sit ahead
# of ufw's INPUT chain, so ufw does NOT protect a port published on 0.0.0.0 --
# which is exactly why docker-compose.yml publishes to 127.0.0.1 instead. This
# is defence in depth for everything else on the box.
ufw --force reset >/dev/null
ufw default deny incoming >/dev/null
ufw default allow outgoing >/dev/null
ufw allow 22/tcp comment 'ssh' >/dev/null
# Note the absence of `ufw allow 8080`. That is intentional.
ufw --force enable >/dev/null
ufw status verbose | sed 's/^/    /'

# ---------------------------------------------------------------------------
log "Building and starting the service"
cd "${REPO_DIR}/deploy"
PTO_DATA_DIR="${DATA_DIR}" docker compose up -d --build

log "Waiting for the health check"
deadline=$(( $(date +%s) + 120 ))
until curl -fsS "http://127.0.0.1:${PORT}/api/health" >/dev/null 2>&1; do
    [[ $(date +%s) -lt ${deadline} ]] || {
        docker compose logs --tail 50
        die "service did not become healthy within 120s"
    }
    sleep 3
done
curl -fsS "http://127.0.0.1:${PORT}/api/health" | sed 's/^/    /'

# ---------------------------------------------------------------------------
public_ip="$(curl -fsS --max-time 3 http://169.254.169.254/latest/meta-data/public-ipv4 2>/dev/null || echo '<instance-ip>')"

cat <<EOF

$(printf '\033[1;32m')Provisioning complete.$(printf '\033[0m')

  Service   : cuttag-profiler, listening on 127.0.0.1:${PORT} of this instance
  Data root : ${DATA_DIR}  (mounted read-only into the container)

$(printf '\033[1;33m')ACCESS -- from your laptop, not from a browser on the instance:$(printf '\033[0m')

  ssh -N -L ${PORT}:127.0.0.1:${PORT} ${APP_USER}@${public_ip}

  then open  http://127.0.0.1:${PORT}/

$(printf '\033[1;33m')SECURITY GROUP$(printf '\033[0m')

  Inbound rules should be exactly one:

    Type: SSH   Protocol: TCP   Port: 22   Source: <your.ip.address>/32

  Do NOT add an inbound rule for ${PORT}. cuttag_profiler has no
  authentication of any kind (SECURITY_HTTP H1c is unimplemented), so an
  internet-reachable ${PORT} lets anyone profile, download and enumerate every
  file under the data root. H2 is closed, so that is now bounded by the mount
  rather than being the whole filesystem -- a smaller blast radius, not a safe
  one. The SSH tunnel above gives the same dashboard with the instance's
  existing key-based auth in front of it.

  If you genuinely need direct exposure, put an authenticating reverse proxy
  in front (ALB + Cognito, or nginx with mTLS) and restrict the source range.

MANAGEMENT

  cd ${REPO_DIR}/deploy
  docker compose logs -f
  docker compose restart
  docker compose down

EOF
