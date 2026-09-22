## Summary

<!-- What does this change do, and why? Link any related issue. -->

## Which module(s) does this touch?

- [ ] `modules/scrna_matrix`
- [ ] `modules/fastq_stream`
- [ ] `modules/cuttag_profiler`
- [ ] `modules/genomic_toolkit`
- [ ] `modules/peaks`
- [ ] `packaging/` (the `pto-core` Python wheel/CLI)
- [ ] Repo-wide (CI, docs, build scripts)

## Testing

- [ ] Ran the affected module's C++ test suite locally and it passes
      (e.g. `ctest --test-dir build -R <suite> --output-on-failure`, or the
      test binary directly)
- [ ] Ran the full top-level `ctest --test-dir build --output-on-failure` if
      this could plausibly affect more than one module
- [ ] Ran `python -m pytest packaging/tests -q` if this touches
      `packaging/pto/` or the `pto-core` CLI
- [ ] Added or updated a test that would have caught this bug/regression, or
      explained below why one isn't practical
- [ ] Ran the relevant sanitizer build (ASan/UBSan/TSan) if this touches
      concurrency, pointer arithmetic, or byte-level parsing

## Dependencies

- [ ] This change introduces **no new external runtime dependency**
      (Python package or system library) — or, if it does, it's named and
      justified below
- [ ] If a new dependency was unavoidable, it's optional / configure-time
      detected rather than hard-required, consistent with the "zero runtime
      dependencies" goal (see the README's `Architecture & System Design`
      section)
- [ ] The affected module still configures and builds standalone with the
      other four deleted (`cmake -S modules/<name> -B build/<name>`)

## Checklist

- [ ] This does not introduce a shared top-level `include/`/`python/`/lib
      layer between modules — that was tried, reverted, and is a deliberate
      non-goal (README, "Repository layout")
- [ ] `shellcheck -S style scripts/*.sh` passes, if a shell script changed
- [ ] Docs updated (module README, `CLAUDE.md`, or a `docs/AUDIT_*.md` /
      `docs/TORTURE_*.md` entry) if behavior, a security-relevant invariant,
      or a public flag changed

<!--
If this fixes a defect found by inspection, fuzzing, or a sanitizer rather
than a user bug report, consider whether it belongs in one of the
docs/AUDIT_*.md or docs/TORTURE_*.md living records rather than only in the
commit message.
-->
