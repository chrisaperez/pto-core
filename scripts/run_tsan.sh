#!/usr/bin/env bash
#
# ThreadSanitizer verification for pto-core.
#
# Builds the whole tree with `-fsanitize=thread` (plus real OpenMP where the
# toolchain supports the combination) and runs every C++ suite under it.
#
# Designed to be run inside docker/Dockerfile.tsan, where clang + libomp-dev +
# libarcher make `-fopenmp -fsanitize=thread` a supported combination. It also
# runs directly on a Linux host with the same packages.
#
#   docker build -f docker/Dockerfile.tsan -t pto-tsan .
#   docker run --rm pto-tsan
#
# EXIT STATUS
#   0  every suite ran under TSan and reported no races
#   1  a race was detected, or a suite failed
#   2  the toolchain cannot support the requested configuration
#
# This script never degrades silently. If OpenMP is unavailable it refuses to
# "verify" scrna_matrix, because a serial build of an OpenMP library proves
# nothing about its parallel paths -- that exact false sense of coverage is what
# the hardening pass was cleaning up.

set -Eeuo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
readonly REPO_ROOT
BUILD_DIR="${BUILD_DIR:-${REPO_ROOT}/build-tsan}"
readonly BUILD_DIR
# Thread counts to sweep. Races are scheduling-dependent, so one run at one
# width is not evidence; oversubscribing past the core count perturbs the
# scheduler and is where the interesting interleavings show up.
read -r -a THREAD_COUNTS <<<"${TSAN_THREAD_COUNTS:-2 8 32}"
readonly THREAD_COUNTS
readonly REPEATS="${TSAN_REPEATS:-3}"

# Which modules to build. Defaults to the whole tree; `TSAN_MODULES=scrna_matrix`
# narrows it, which is what you want when the target is the OpenMP engine
# specifically -- fastq_stream and cuttag_profiler are std::thread-based and can
# be TSan'd on any host without this container.
readonly TSAN_MODULES="${TSAN_MODULES:-all}"
# Optional ctest -R regex, for iterating on one suite.
readonly TSAN_CTEST_FILTER="${TSAN_CTEST_FILTER:-}"

log()  { printf '\033[1;34m==>\033[0m %s\n' "$*"; }
warn() { printf '\033[1;33m[warn]\033[0m %s\n' "$*" >&2; }
die()  { printf '\033[1;31m[fail]\033[0m %s\n' "$*" >&2; exit "${2:-1}"; }

trap 'die "aborted at line ${LINENO}"' ERR

# --------------------------------------------------------------------------
# Toolchain discovery and validation
# --------------------------------------------------------------------------
CXX="${CXX:-clang++}"
CC="${CC:-clang}"
command -v "${CXX}" >/dev/null 2>&1 || die "compiler '${CXX}' not found" 2

log "toolchain: $("${CXX}" --version | head -n1)"

scratch="$(mktemp -d)"
trap 'rm -rf "${scratch}"' EXIT

# 1. TSan must actually link, not merely be accepted as a flag. Homebrew GCC on
#    macOS accepts -fsanitize=thread and then fails at link with "library 'tsan'
#    not found"; catching that here gives a usable message instead of a wall of
#    linker output 200 files later.
cat >"${scratch}/tsan_probe.cpp" <<'EOF'
#include <thread>
int main() { std::thread t([]{}); t.join(); return 0; }
EOF
if ! "${CXX}" -fsanitize=thread -o "${scratch}/tsan_probe" "${scratch}/tsan_probe.cpp" \
        >"${scratch}/tsan_probe.log" 2>&1; then
    sed 's/^/    /' "${scratch}/tsan_probe.log" >&2
    die "'${CXX} -fsanitize=thread' does not link on this host.
    On macOS this is expected: Apple Clang has no OpenMP and Homebrew GCC has no
    sanitizer runtimes. Use the container:
        docker build -f docker/Dockerfile.tsan -t pto-tsan .
        docker run --rm pto-tsan" 2
fi

# 2. OpenMP must combine with TSan, since scrna_matrix's parallel regions are
#    the primary target of this run.
cat >"${scratch}/omp_probe.cpp" <<'EOF'
#include <omp.h>
#ifndef _OPENMP
#error _OPENMP undefined
#endif
int main() { int n = 0;
#pragma omp parallel reduction(+ : n)
  { n += 1; }
  return n > 0 ? 0 : 1; }
EOF
if "${CXX}" -fopenmp -fsanitize=thread -o "${scratch}/omp_probe" \
        "${scratch}/omp_probe.cpp" >"${scratch}/omp_probe.log" 2>&1; then
    log "OpenMP + TSan: supported"
else
    sed 's/^/    /' "${scratch}/omp_probe.log" >&2
    die "'${CXX} -fopenmp -fsanitize=thread' does not build.
    scrna_matrix's parallel regions cannot be verified without it, and building
    it serially would report a clean run that proves nothing.
    Install libomp-dev (Debian/Ubuntu) or use docker/Dockerfile.tsan." 2
fi

# 3. libarcher. Without OMPT annotations TSan cannot see libomp's barriers and
#    reports false races inside every parallel region -- which trains people to
#    ignore the output, the worst possible outcome for this harness.
ARCHER_PATH="${ARCHER_PATH:-}"
if [[ -z "${ARCHER_PATH}" && -r /etc/pto-archer-path ]]; then
    ARCHER_PATH="$(cat /etc/pto-archer-path)"
fi
if [[ -z "${ARCHER_PATH}" ]]; then
    ARCHER_PATH="$(find /usr/lib /usr/local/lib -name 'libarcher.so*' 2>/dev/null | head -n1 || true)"
fi
if [[ -n "${ARCHER_PATH}" && -r "${ARCHER_PATH}" ]]; then
    export OMP_TOOL=enabled
    export OMP_TOOL_LIBRARIES="${ARCHER_PATH}"
    log "libarcher: ${ARCHER_PATH}"
else
    warn "libarcher not found. TSan cannot see libomp's internal happens-before"
    warn "edges, so reports inside parallel regions may be FALSE POSITIVES."
    warn "Install libomp-dev, or treat any race inside libomp with suspicion."
fi

SUPP=""
for candidate in /etc/pto-tsan.supp "${REPO_ROOT}/docker/tsan.supp"; do
    if [[ -r "${candidate}" ]]; then SUPP="${candidate}"; break; fi
done
[[ -n "${SUPP}" ]] && log "suppressions: ${SUPP}"

# halt_on_error=0 so one race does not mask the rest; history_size buys deeper
# stacks for the second thread, which is usually the one you actually need.
#
# detect_deadlocks defaults OFF. This is a race detector run, and the deadlock
# detector is a separate analysis that is unusable here for a specific,
# measured reason: on this platform its reports carry NO stack frames at all
# ("#0 <null> <null> (0x000000000042)"), so they cannot be attributed and --
# critically -- cannot be matched by any suppression rule, since TSan
# suppressions match on frames. Left on, it emits ~500 unattributable warnings
# per run and drowns the output the harness exists to surface.
#
# The reports were attributed by experiment rather than by assumption: building
# with -DSCRNA_ENABLE_HNSW=OFF drops them from 7,560 across six runs to ZERO
# while leaving all first-party parallel code in place. They are inversions
# between hnswlib's per-element mutexes. Our own parallel regions take no
# std::mutex at all -- one named `#pragma omp critical` -- so an inversion
# cannot originate there. Recorded in docker/tsan.supp §2 and in
# third_party/hnswlib/VERSION.txt.
#
# Set TSAN_DETECT_DEADLOCKS=1 to re-enable, e.g. after an hnswlib bump.
TSAN_OPTIONS="halt_on_error=0:exitcode=66:history_size=7"
TSAN_OPTIONS="${TSAN_OPTIONS}:detect_deadlocks=${TSAN_DETECT_DEADLOCKS:-0}"
if [[ "${TSAN_DETECT_DEADLOCKS:-0}" != "0" ]]; then
    TSAN_OPTIONS="${TSAN_OPTIONS}:second_deadlock_stack=1"
fi
[[ -n "${SUPP}" ]] && TSAN_OPTIONS="${TSAN_OPTIONS}:suppressions=${SUPP}"
export TSAN_OPTIONS
# Archer is quiet about non-races only if told to be.
export ARCHER_OPTIONS="${ARCHER_OPTIONS:-print_max_rss=0}"

# --------------------------------------------------------------------------
# Configure and build
# --------------------------------------------------------------------------
log "configuring in ${BUILD_DIR}"
rm -rf "${BUILD_DIR}"

module_flags=()
case "${TSAN_MODULES}" in
    all) log "modules: all" ;;
    scrna_matrix)
        module_flags=(-DPTO_BUILD_FASTQ_STREAM=OFF -DPTO_BUILD_CUTTAG_PROFILER=OFF)
        log "modules: scrna_matrix only" ;;
    *)
        die "TSAN_MODULES must be 'all' or 'scrna_matrix', got '${TSAN_MODULES}'" 2 ;;
esac

# RelWithDebInfo, not Debug: TSan needs symbols, and -O1+ keeps the suites fast
# enough to sweep thread counts. -fno-omit-frame-pointer keeps stacks readable.
cmake -S "${REPO_ROOT}" -B "${BUILD_DIR}" -G Ninja \
    -DCMAKE_BUILD_TYPE=RelWithDebInfo \
    -DCMAKE_C_COMPILER="${CC}" \
    -DCMAKE_CXX_COMPILER="${CXX}" \
    -DCMAKE_C_FLAGS="-fsanitize=thread -fno-omit-frame-pointer -g" \
    -DCMAKE_CXX_FLAGS="-fsanitize=thread -fno-omit-frame-pointer -g" \
    -DCMAKE_EXE_LINKER_FLAGS="-fsanitize=thread" \
    -DCMAKE_SHARED_LINKER_FLAGS="-fsanitize=thread" \
    -DFQ_NATIVE=OFF \
    -DSCRNA_HNSW_NATIVE_ISA=OFF \
    -DPROFILER_STATIC_HTSLIB=OFF \
    "${module_flags[@]}"

log "building"
cmake --build "${BUILD_DIR}" -j "$(nproc 2>/dev/null || echo 4)"

# --------------------------------------------------------------------------
# Run
# --------------------------------------------------------------------------
# A TSan-instrumented extension module dlopened into a stock CPython aborts with
# "interceptors are not working ... loaded too late", because the runtime must be
# present from process start. Preloading it lets the Python binding suite run
# under TSan rather than be silently skipped.
#
# Scoped to that suite ONLY. The C++ test binaries already link the same runtime,
# and LD_PRELOADing it into them a second time makes them SEGFAULT instantly at
# startup -- which is how this was found: test_matrix_ops died in 0.01s while the
# Python suite passed 34/34 in the same run.
TSAN_RUNTIME="${TSAN_RUNTIME:-}"
if [[ -z "${TSAN_RUNTIME}" && -r /etc/pto-tsan-runtime ]]; then
    TSAN_RUNTIME="$(cat /etc/pto-tsan-runtime)"
fi
if [[ -n "${TSAN_RUNTIME}" && -r "${TSAN_RUNTIME}" ]]; then
    log "TSan runtime for the Python suite: ${TSAN_RUNTIME}"
else
    warn "TSan runtime not located; the Python binding suite may fail to load"
    warn "the instrumented module rather than running under TSan."
fi

# Runs ctest with the right preload for the selected subset and echoes the exit
# status. `native` = binaries linked against TSan (no preload); `python` = the
# pytest suite that dlopens the instrumented module (preload required).
run_ctest_subset() {
    local kind="$1" threads="$2" outfile="$3"
    shift 3
    local -a args=(--test-dir "${BUILD_DIR}" --output-on-failure --timeout 900 "$@")
    [[ -n "${TSAN_CTEST_FILTER}" ]] && args+=(-R "${TSAN_CTEST_FILTER}")
    local rc=0
    if [[ "${kind}" == "python" && -n "${TSAN_RUNTIME}" ]]; then
        LD_PRELOAD="${TSAN_RUNTIME}" OMP_NUM_THREADS="${threads}" \
            ctest "${args[@]}" >"${outfile}" 2>&1 || rc=$?
    else
        OMP_NUM_THREADS="${threads}" ctest "${args[@]}" >"${outfile}" 2>&1 || rc=$?
    fi
    return "${rc}"
}

races=0
failures=0
ran=0
summary=()

for threads in "${THREAD_COUNTS[@]}"; do
    for ((rep = 1; rep <= REPEATS; rep++)); do
        log "ctest: OMP_NUM_THREADS=${threads} repeat ${rep}/${REPEATS}"
        out="${scratch}/ctest_${threads}_${rep}.log"
        # `|| rc=$?` inside run_ctest_subset rather than `set +e`: a bare
        # `set +e` does NOT disable an ERR trap, so an earlier version aborted
        # the whole harness on the first nonzero ctest instead of collecting the
        # result -- exactly the case this script exists to report on
        # (TSAN_OPTIONS sets exitcode=66, so a detected race *is* a nonzero
        # exit). A command on the left of `||` is exempt from both.
        rc=0
        pyout="${out%.log}.python.log"
        run_ctest_subset native "${threads}" "${out}" -E test_python_bindings || rc=$?
        rc_py=0
        run_ctest_subset python "${threads}" "${pyout}" -R test_python_bindings || rc_py=$?
        cat "${pyout}" >>"${out}"
        [[ ${rc} -eq 0 ]] && rc=${rc_py}
        ran=$((ran + 1))

        # Classify by report type. A data race and a lock-order inversion are
        # different analyses with different severities, and collapsing them into
        # one "RACE" verdict (an earlier version of this script did) makes a
        # clean race result indistinguishable from a deadlock-detector flood.
        n_race=$(grep -c "ThreadSanitizer: data race" "${out}" || true)
        n_other=$(grep -c "WARNING: ThreadSanitizer" "${out}" || true)
        n_other=$((n_other - n_race))
        if [[ ${n_race} -gt 0 ]]; then
            races=$((races + 1))
            printf '\033[1;31m--- ThreadSanitizer DATA RACE (threads=%s rep=%s) ---\033[0m\n' \
                "${threads}" "${rep}" >&2
            awk '/ThreadSanitizer: data race/,/^SUMMARY: ThreadSanitizer/' "${out}" >&2
            summary+=("threads=${threads} rep=${rep}: DATA RACE x${n_race}")
        elif [[ ${n_other} -gt 0 ]]; then
            # Non-race warnings that survived docker/tsan.supp. Not a pass:
            # every suppression in that file is justified in writing, so
            # anything getting through is unaccounted for and must be triaged.
            failures=$((failures + 1))
            printf '\033[1;33m--- %s unsuppressed non-race TSan warning(s) (threads=%s rep=%s) ---\033[0m\n' \
                "${n_other}" "${threads}" "${rep}" >&2
            grep -m5 "WARNING: ThreadSanitizer" "${out}" >&2 || true
            summary+=("threads=${threads} rep=${rep}: ${n_other} unsuppressed warning(s)")
        elif [[ ${rc} -ne 0 ]]; then
            failures=$((failures + 1))
            tail -n 60 "${out}" >&2
            summary+=("threads=${threads} rep=${rep}: TEST FAILURE (rc=${rc})")
        else
            summary+=("threads=${threads} rep=${rep}: clean")
        fi
    done
done

# --------------------------------------------------------------------------
# Report
# --------------------------------------------------------------------------
echo
if [[ -n "${TSAN_LOG_DIR:-}" ]]; then
    mkdir -p "${TSAN_LOG_DIR}"
    cp "${scratch}"/ctest_*.log "${TSAN_LOG_DIR}/" 2>/dev/null || true
    log "logs copied to ${TSAN_LOG_DIR}"
fi

log "summary"
printf '  %s\n' "${summary[@]}"
echo

if [[ ${ran} -eq 0 ]]; then
    die "no test runs executed; refusing to report success"
fi
if [[ ${races} -gt 0 ]]; then
    die "ThreadSanitizer reported races in ${races}/${ran} run(s)"
fi
if [[ ${failures} -gt 0 ]]; then
    die "${failures}/${ran} run(s) failed for reasons other than a data race"
fi

log "clean: ${ran} run(s) across thread counts [${THREAD_COUNTS[*]}], no races"
