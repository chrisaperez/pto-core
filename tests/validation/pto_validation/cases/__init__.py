# SPDX-License-Identifier: MIT
"""Validation cases, one module per pto module.

Importing any submodule here registers its cases with
`pto_validation.harness.REGISTRY` as a side effect of the `@case`
decorator. `harness.load_cases()` performs those imports; nothing else should
need to import these modules directly.
"""
