# Security Policy

## Reporting a vulnerability

Please report suspected security issues privately, not as a public GitHub
issue — even a low-severity one, since a public report is also a disclosure
to anyone who wants to exploit it before a fix ships.

Preferred channel: **GitHub Security Advisories**
(Security tab → "Report a vulnerability" on this repository). It reaches the
maintainer directly, supports a private discussion thread, and can mint a
CVE once a fix is ready.

If that option isn't available to you for some reason, open a regular issue
asking for a private contact channel and describe severity only in broad
terms until one is established — don't post exploit details there.

Please include, where relevant:

- Which module (`fastq_stream`, `cuttag_profiler`, `genomic_toolkit`,
  `peaks`, `scrna_matrix`) and which entry point.
- A minimal reproduction — a crafted input file or command line is far more
  useful than a description of the class of bug.
- What you observed versus what you expected (a crash, a wrong answer with
  no error, a memory-safety report from a sanitizer, etc.).
- Whether you've checked it against the audit documents below — if it's
  already a recorded, closed finding, saying so saves a round trip.

## What's in scope

This is a set of command-line tools and libraries that run locally on data
you already control, not a hosted service — most of the usual web
vulnerability classes don't apply. What's realistically in scope:

- **`cuttag_profiler`'s embedded HTTP server** (`src/http_server.cpp`,
  `include/profiler/http_security.hpp`, `include/profiler/secrets.hpp`) —
  the one component in this repo that binds a socket. Path confinement,
  session-token handling, `Host`/`Origin` validation.
- **Parsers over untrusted or malformed input**: BAM/BGZF decoding (`peaks`,
  `cuttag_profiler`, `fastq_stream`, `genomic_toolkit`), BED/BEDPE/GTF
  readers, FASTQ record framing. A crafted file that causes a crash, an
  out-of-bounds read/write, or — worse — a plausible but silently wrong
  answer with exit code 0.
- **Memory safety and concurrency** anywhere: the HNSW index in
  `scrna_matrix`, the lock-free queue in `fastq_stream`, anything a sanitizer
  build would flag.
- **Path handling**: an output path that resolves to an input path, symlink
  or FIFO handling, anything under `file_identity.hpp`/`safe_open.hpp` in
  each module.

Out of scope: the demo dashboard's authentication model is intentionally a
single shared per-run token, not a login system — see
[`deploy/README.md`](deploy/README.md) for what it's actually meant to
protect against (a hostile web page on the same tunnel, not a hostile
co-tenant) before reporting it as a missing-auth finding.

## Track record

This project has an unusually explicit adversarial-review history for its
size — worth reading before assuming a class of bug hasn't already been
looked for. Findings are logged before they're fixed, and residuals are
stated rather than hidden:

| Document | Scope |
|:--|:--|
| [`docs/REVIEW_2026-08-15.md`](docs/REVIEW_2026-08-15.md) | Cross-module architecture and security review |
| [`docs/SECURITY_HTTP_2026-08-15.md`](docs/SECURITY_HTTP_2026-08-15.md) | `cuttag_profiler`'s HTTP server and path confinement (9 findings) |
| [`docs/TORTURE_2026-08-17.md`](docs/TORTURE_2026-08-17.md) | Adversarial fuzzing, sanitizer, and concurrency sweep |
| [`docs/TORTURE_2026-09-10.md`](docs/TORTURE_2026-09-10.md) | `fastq_stream`, then `pto-peaks` |
| [`docs/AUDIT_2026-09-11_fastq_stream.md`](docs/AUDIT_2026-09-11_fastq_stream.md) | Phred offsets, record framing, truncation |
| [`docs/AUDIT_2026-09-11_genomic_toolkit.md`](docs/AUDIT_2026-09-11_genomic_toolkit.md) | `sizes`, `markdup`, `frip` |
| [`docs/AUDIT_2026-09-11_cuttag_profiler.md`](docs/AUDIT_2026-09-11_cuttag_profiler.md) | `profile`, bin geometry, skipped regions |
| [`docs/AUDIT_2026-09-11_scrna_matrix.md`](docs/AUDIT_2026-09-11_scrna_matrix.md) | HNSW producers, index concurrency, NumPy boundary |

74 confirmed defects across these documents, 24 rated High or Critical, every
one with a regression test. If you find something these missed, it's a
genuinely useful report, not a sign the process doesn't work — four of the
five findings in the cross-module review were the same defect recurring in a
sibling function a previous pass hadn't walked, which is exactly the kind of
thing an external report is good at catching.

## Supported versions

Pre-1.0 (`0.x`): only the latest released version is supported. There is no
LTS branch and no backport policy at this stage.

## Disclosure

No fixed disclosure timeline is promised, but the intent — matching the
audit documents above — is to fix before disclosing, credit the reporter
(unless you'd rather stay anonymous), and add a regression test in the
commit that closes it.
