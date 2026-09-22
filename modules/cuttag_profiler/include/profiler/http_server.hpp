// SPDX-License-Identifier: MIT
//
// Embedded HTTP server exposing the profiling core to a local browser.
//
// Air-gap guarantees enforced here:
//   * every asset is served from a byte array compiled into the binary, so no
//     CDN, font host or analytics endpoint is ever contacted;
//   * a Content-Security-Policy of `default-src 'self'` is attached to the
//     document response, so even a modified asset cannot exfiltrate data.
//
// Inbound access is a separate question, and the loopback bind is NOT the
// answer to it. Binding 127.0.0.1 restricts which network interfaces reach the
// socket; it does not restrict which web origins can make the operator's own
// browser issue requests to it, and a DNS-rebinding page reads the replies
// (SECURITY_HTTP_2026-08-15 H1). Three controls do that work, all in a
// pre-routing handler covering /api/:
//   * a loopback-only `Host` allowlist, which is what actually defeats
//     rebinding -- a rebound page still sends the attacker's hostname;
//   * an `Origin` check plus a Content-Type gate on POST, for CSRF; and
//   * a session token, minted per run from the platform CSPRNG, which is the
//     only one of the three that stops a hostile process on the same machine.
//     A local process sets whatever headers it likes; it cannot guess this.
//
// This whole file is compiled out under PTO_CLOUD_BUILD. The token makes the
// dashboard safe for one operator on one workstation; it does not make a
// listening port acceptable inside a multi-tenant worker, where per-job task
// isolation is the boundary. See the option's comment in CMakeLists.txt.
#pragma once

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "profiler/bam_reader.hpp"
#include "profiler/types.hpp"

namespace profiler {

struct ServerOptions {
    std::string host = "127.0.0.1";
    int port = 8080;

    // Root directory every request is confined to. Requests that escape it
    // after canonicalisation are rejected.
    //
    // REQUIRED. An empty value is no longer "no confinement" -- that was
    // SECURITY_HTTP_2026-08-15 finding H2, a Critical unauthenticated arbitrary
    // file read that shipped as the DEFAULT configuration. The server refuses
    // to start without either this or an explicit `confinement_disabled`.
    std::string data_root;

    // Explicit, deliberate opt-out of confinement.
    //
    // A separate flag rather than "data_root is empty" on purpose. The original
    // defect was that the unconfined state was reachable by *omission* -- you
    // got it by not thinking about it. Making it a named boolean means the
    // unsafe configuration cannot be entered by accident, shows up in `ps` and
    // in shell history, and is greppable in any wrapper that sets it.
    bool confinement_disabled = false;

    // Open the dashboard in the user's default browser once bound.
    bool open_browser = true;
    // Refuse to bind to anything other than a loopback address.
    bool loopback_only = true;

    // Mint a session token at startup and require it on every /api/ request
    // (SECURITY_HTTP H1c). The dashboard URL printed at startup carries it, so
    // the operator never types or copies it.
    //
    // ON by default, and off only via --no-auth, for the same reason
    // confinement works that way after H2: the unsafe configuration must not
    // be reachable by omission. Turning it off is legitimate for a scripted
    // local client that would otherwise have to scrape the token out of the
    // banner -- but validate() refuses to combine it with --allow-remote.
    bool require_token = true;

    // Throws std::invalid_argument unless the confinement configuration is
    // coherent: a data_root that exists and is a directory, or an explicit
    // opt-out. Also canonicalises data_root so the confinement check compares
    // like with like.
    //
    // Defined on the options type, not in main.cpp's argument parser, so that
    // every construction path is covered. HttpServer's constructor calls it.
    // The CLI validated S4/S5/S6 while the HTTP handler did not, and that gap
    // was REVIEW_2026-08-15 finding 2; validation belongs on the type.
    void validate();
};

class HttpServer {
public:
    explicit HttpServer(ServerOptions options);
    ~HttpServer();

    HttpServer(const HttpServer&) = delete;
    HttpServer& operator=(const HttpServer&) = delete;

    // Binds and serves until `stop()` is called. Throws std::runtime_error if
    // the port is unavailable or a non-loopback host was requested while
    // `loopback_only` is set.
    void run();
    void stop();

    // The bound port; differs from the requested one when port 0 was asked for.
    [[nodiscard]] int bound_port() const noexcept { return bound_port_; }

    // The session token this server will require, minted in the constructor so
    // it is available before run() blocks. Empty when require_token is off.
    //
    // Exposed for the dashboard URL and for tests. It is a secret: never log
    // it, and compare it with profiler::constant_time_equals.
    [[nodiscard]] const std::string& session_token() const noexcept {
        return session_token_;
    }

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
    ServerOptions options_;
    int bound_port_ = 0;
    std::string session_token_;
};

// Opens `url` in the platform's default browser. Best-effort and non-fatal.
void open_in_browser(const std::string& url);

}  // namespace profiler
