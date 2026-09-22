#!/bin/sh
#
# Container entrypoint for cuttag_profiler.
#
# A script rather than an inline `sh -c` in the Dockerfile, so the bind address
# can be validated before it reaches the binary.
#
# Why validate it at all: SECURITY_HTTP_2026-08-15 L8 (a std::system() shell
# injection via --host in open_in_browser()) is closed upstream -- that path
# now uses posix_spawnp with an explicit argv, so nothing in --host reaches a
# shell regardless of this script. The allowlist here is defence in depth on
# top of that fix, not a substitute for it: it means a mistyped or hostile
# PROFILER_BIND is rejected before it reaches the binary at all.
set -eu

PROFILER_BIND="${PROFILER_BIND:-0.0.0.0}"
PROFILER_PORT="${PROFILER_PORT:-8080}"
PROFILER_DATA_ROOT="${PROFILER_DATA_ROOT:-/data}"

case "${PROFILER_BIND}" in
    0.0.0.0|127.0.0.1|localhost|::1) ;;
    *)
        echo "entrypoint: PROFILER_BIND must be one of 0.0.0.0, 127.0.0.1, localhost, ::1" >&2
        echo "entrypoint: got '${PROFILER_BIND}'" >&2
        exit 64
        ;;
esac

case "${PROFILER_PORT}" in
    ''|*[!0-9]*)
        echo "entrypoint: PROFILER_PORT must be numeric, got '${PROFILER_PORT}'" >&2
        exit 64
        ;;
esac

if [ ! -d "${PROFILER_DATA_ROOT}" ]; then
    echo "entrypoint: data root '${PROFILER_DATA_ROOT}' is not a directory." >&2
    echo "entrypoint: mount your BAM/BED files there, e.g." >&2
    echo "            docker run -v /srv/experiments:/data:ro ..." >&2
    exit 66
fi

# Refuse to start unconfined. The binary now enforces this itself (H2 closed:
# --data-root is mandatory), but the check stays here so the image fails with a
# message about the *mount* the operator forgot rather than about a flag the
# entrypoint controls. There is no supported way to run this image unconfined.
if [ -z "${PROFILER_DATA_ROOT}" ]; then
    echo "entrypoint: refusing to start without a data root" >&2
    exit 64
fi

echo "cuttag_profiler: bind=${PROFILER_BIND}:${PROFILER_PORT} data-root=${PROFILER_DATA_ROOT} (read-only)"
echo "cuttag_profiler: session-token auth only (single shared secret, no login system) -- do not publish this port beyond loopback. See deploy/README.md."

# --allow-remote is required because 0.0.0.0 is a non-loopback bind. Inside a
# container that is unavoidable if the port is to be published at all; the real
# exposure boundary is the host-side publish address, which docker-compose pins
# to 127.0.0.1.
exec cuttag_profiler serve \
    --host "${PROFILER_BIND}" \
    --port "${PROFILER_PORT}" \
    --data-root "${PROFILER_DATA_ROOT}" \
    --no-browser \
    --allow-remote
