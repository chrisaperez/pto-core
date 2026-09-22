# SPDX-License-Identifier: MIT
"""Public-data end-to-end validation for pto-core.

Entry point is `tests/validation/run_validation.py`; this package holds the
machinery it drives. Standard library only, by the same reasoning that keeps
the C++ test suites framework-free -- these tools run on clusters that have a
compiler and an interpreter and nothing else, and the harness that proves they
are correct has to run in the same place.

The one exception is the scrna_matrix cases, which need numpy and anndata to
read a count matrix at all. They skip cleanly when those are absent rather
than making them a suite-wide requirement.
"""

__version__ = "1.0.0"
