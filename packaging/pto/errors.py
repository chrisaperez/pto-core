"""Exceptions, arranged around what the caller can do about them.

The split between `ToolError` and `InvalidInput` is not cosmetic -- it is
`genomic_toolkit`'s exit-code contract, which distinguishes a system failure
from a deterministic property of the customer's data. pto-cloud encodes the
same distinction in `api/app/job_specs.TOOL_EXIT_CODES` because it decides
whether to retry a billed Fargate task; a workstation caller has the same
decision to make in a `for` loop over 300 samples, and the same answer:
retrying an `InvalidInput` reaches the same answer, more slowly.
"""

from __future__ import annotations

__all__ = ["BinaryNotFound", "InvalidInput", "PtoError", "ToolError"]


class PtoError(Exception):
    """Base class for everything this package raises."""


class BinaryNotFound(PtoError):
    """A pto-core binary could not be located. See `pto doctor`."""


class ToolError(PtoError):
    """A binary exited non-zero for a reason that is not the input's shape.

    Exit 1: unreadable input, a rejected option, an exception in the tool.
    Attributes let a caller log the failure without re-deriving it.
    """

    def __init__(self, tool: str, argv: list[str], returncode: int, stderr: str):
        self.tool = tool
        self.argv = list(argv)
        self.returncode = returncode
        self.stderr = stderr
        detail = stderr.strip().splitlines()
        tail = detail[-1] if detail else "(no stderr)"
        super().__init__(f"{tool} exited {returncode}: {tail}")


class InvalidInput(ToolError):
    """Exit 2 -- deterministic, user-actionable, and NOT worth retrying.

    `genomic_toolkit` returns 2 for exactly two conditions, both properties of
    the data rather than of the run:

      * `frip`: no fragment overlapped a peak while the stream's contigs are
        absent from the peak file -- almost always `chr1` vs `1` naming.
      * `markdup`: records arrived out of coordinate order, so the duplicate
        rate is an undercount. Fix with `unsorted=True` or a coordinate sort.

    Both produce a VALID report before exiting 2, so `result` carries it. That
    is the reason this is a distinct exception with a payload rather than a
    message: the number is still worth showing, next to the reason it is
    suspect.
    """

    def __init__(self, tool: str, argv: list[str], stderr: str, result: object = None):
        super().__init__(tool, argv, 2, stderr)
        self.result = result
