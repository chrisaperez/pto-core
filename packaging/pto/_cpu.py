"""How many CPUs this process may actually use.

Not `os.cpu_count()`. That reads the machine, and every environment a
sequencing core runs jobs in gives a process less than the machine:

  * SLURM/LSF/SGE confine a task with cgroups and CPU affinity. A node with
    128 cores allocated `--cpus-per-task=4` still reports 128 from
    `hardware_concurrency()`.
  * Docker/Fargate/Kubernetes apply a cgroup CPU quota that no libc API
    reports as a core count.

pto-cloud already carries this rule for the cloud deployment -- see
`api/app/job_specs._threads_arg`, whose docstring records a 2-vCPU task
spawning 64 workers and spending its life in scheduler contention. The
binaries all default their worker count to
`std::thread::hardware_concurrency()`, so the fix is the same on a cluster
node as it is on Fargate: pass the thread count explicitly, from something
that knows what the scheduler actually granted.

This module is that something. It is stdlib-only and never raises.
"""

from __future__ import annotations

import os
from pathlib import Path

__all__ = ["available_cpus"]

_CGROUP_V2_MAX = Path("/sys/fs/cgroup/cpu.max")
_CGROUP_V1_QUOTA = Path("/sys/fs/cgroup/cpu/cpu.cfs_quota_us")
_CGROUP_V1_PERIOD = Path("/sys/fs/cgroup/cpu/cpu.cfs_period_us")


def _affinity() -> int | None:
    """CPUs this process is pinned to. Linux only; absent on macOS."""
    getaffinity = getattr(os, "sched_getaffinity", None)
    if getaffinity is None:
        return None
    try:
        return len(getaffinity(0))
    except OSError:
        return None


def _cgroup_quota() -> int | None:
    """The cgroup CPU quota, rounded DOWN, as a whole number of CPUs.

    Rounded down rather than up, and floored at 1: a 1.5-CPU quota gets one
    worker. Rounding up hands the scheduler more runnable threads than the
    quota can retire, which is throttling rather than parallelism.
    """
    try:
        if _CGROUP_V2_MAX.is_file():
            quota_s, period_s = _CGROUP_V2_MAX.read_text().split()
            if quota_s == "max":
                return None
            quota, period = int(quota_s), int(period_s)
        elif _CGROUP_V1_QUOTA.is_file() and _CGROUP_V1_PERIOD.is_file():
            quota = int(_CGROUP_V1_QUOTA.read_text().strip())
            period = int(_CGROUP_V1_PERIOD.read_text().strip())
            if quota <= 0:       # -1 == unlimited
                return None
        else:
            return None
        if period <= 0:
            return None
        return max(1, quota // period)
    except (OSError, ValueError):
        # An unreadable or unparseable cgroup file is "no information", not an
        # error worth propagating into a QC run.
        return None


def available_cpus() -> int:
    """Threads this process should ask a pto-core binary to start.

    The MINIMUM of every signal that says "fewer", never the maximum, and
    never less than 1. Being wrong low costs some throughput; being wrong
    high costs the contention this function exists to prevent.
    """
    candidates = [c for c in (_affinity(), _cgroup_quota(), os.cpu_count()) if c]
    return max(1, min(candidates)) if candidates else 1
