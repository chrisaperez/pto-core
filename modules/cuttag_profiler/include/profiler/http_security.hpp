// SPDX-License-Identifier: MIT
//
// Header-validation predicates for the embedded HTTP server.
//
// These live in a header rather than in http_server.cpp's anonymous namespace
// for one reason: they are security controls, and a security control with no
// test is a control that rots. Parsing a Host or Origin header is exactly the
// kind of string handling that acquires bugs quietly -- an IPv6 literal, a
// missing port, a suffix that merely *contains* "127.0.0.1" -- so the cases are
// pinned in tests/test_http_security.cpp.
//
// Threat model, stated so callers do not over-trust these:
//   * They stop a remote web page (DNS rebinding, cross-origin CSRF), because a
//     browser controls Host and Origin and will not let script forge them.
//   * They do NOT stop a hostile process on the same machine, which can send
//     any header it likes. That needs the session token from
//     docs/SECURITY_HTTP_2026-08-15.md H1(c), which is not implemented.
#pragma once

#include <cstddef>
#include <string>
#include <string_view>

namespace profiler {

// True when `host` (a raw Host header value) names this loopback service on
// `port`.
//
// Host is the field that exposes DNS rebinding: after the attacker's DNS record
// flips to 127.0.0.1, the browser still sends the ATTACKER's hostname here,
// because Host reflects the URL the page was loaded from rather than the
// address it resolved to. So anything that is not a loopback name is a rebound
// request, and an absent Host (illegal in HTTP/1.1) is not trusted either.
[[nodiscard]] inline bool host_header_is_local(std::string_view host, int port) {
    if (host.empty()) return false;

    std::string_view name = host;
    const std::string port_digits = std::to_string(port);

    if (name.front() == '[') {
        // Bracketed IPv6 literal: "[::1]" or "[::1]:8080".
        const auto close = name.find(']');
        if (close == std::string_view::npos) return false;
        const std::string_view rest = name.substr(close + 1);
        if (!rest.empty()) {
            if (rest.size() < 2 || rest.front() != ':') return false;
            if (rest.substr(1) != port_digits) return false;
        }
        name = name.substr(1, close - 1);
    } else if (const auto colon = name.rfind(':'); colon != std::string_view::npos) {
        // A bare IPv6 literal without brackets has several colons and is not a
        // legal Host value; rfind would mis-split it, so require exactly one.
        if (name.find(':') != colon) return false;
        if (name.substr(colon + 1) != port_digits) return false;
        name = name.substr(0, colon);
    }

    // Exact match only. A substring test would accept "127.0.0.1.evil.com".
    return name == "127.0.0.1" || name == "localhost" || name == "::1";
}

// True when `origin` (a raw Origin header value) identifies this same server.
//
// Deliberately NOT "reject any request carrying Origin". Browsers send Origin
// on same-origin POSTs too, and the shipped dashboard's own
// POST /api/profile carries `Origin: http://127.0.0.1:<port>`. Rejecting every
// origin would break the UI while still looking like a working control.
//
// An absent Origin is handled by the caller (it is normal for same-origin GETs
// and for curl). The literal "null" -- sent by sandboxed iframes and some
// privacy modes -- does not parse as an http:// origin and is therefore
// refused, which is the intended outcome.
[[nodiscard]] inline bool origin_is_self(std::string_view origin, int port) {
    constexpr std::string_view kScheme = "http://";
    if (origin.rfind(kScheme, 0) != 0) return false;

    std::string_view authority = origin.substr(kScheme.size());
    // An Origin has no path; anything after the authority makes it malformed.
    if (authority.find('/') != std::string_view::npos) return false;

    // Origin omits the port when it is the scheme default (80). This server
    // never binds 80 by default, but be exact rather than lenient.
    return host_header_is_local(authority, port);
}

}  // namespace profiler
