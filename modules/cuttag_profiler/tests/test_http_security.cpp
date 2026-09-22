// SPDX-License-Identifier: MIT
//
// Unit tests for the Host/Origin validation that backs SECURITY_HTTP H1(a),
// H1(b) and H3.
//
// These are the checks standing between the operator's browser and an
// unauthenticated attacker's web page, and they are pure string parsing -- the
// category that acquires quiet bugs. Every case below is one an attacker would
// actually try.
#include <cstdio>

#include <cstddef>
#include <filesystem>
#include <fstream>
#include <string>

#include "profiler/http_security.hpp"
#include "profiler/http_server.hpp"
#include "test_util.hpp"

using profiler::host_header_is_local;
using profiler::HttpServer;
using profiler::origin_is_self;
using profiler::ServerOptions;

namespace {

// Returns true if validate() rejected the configuration.
bool rejects(ServerOptions opts) {
    try {
        opts.validate();
        return false;
    } catch (const std::invalid_argument&) {
        return true;
    }
}

}  // namespace

int main() {
    testing::Suite suite{"http_security", {}};
    constexpr int kPort = 8080;

    suite.add("loopback Host values on the bound port are accepted", [] {
        CHECK(host_header_is_local("127.0.0.1:8080", kPort));
        CHECK(host_header_is_local("localhost:8080", kPort));
        CHECK(host_header_is_local("[::1]:8080", kPort));
        // Port omitted: some clients drop it. The name is what matters.
        CHECK(host_header_is_local("127.0.0.1", kPort));
        CHECK(host_header_is_local("localhost", kPort));
        CHECK(host_header_is_local("[::1]", kPort));
    });

    suite.add("a rebound Host is rejected", [] {
        // The DNS-rebinding case: the address resolves to 127.0.0.1 but the
        // browser still sends the attacker's name here.
        CHECK(!host_header_is_local("evil.example.com:8080", kPort));
        CHECK(!host_header_is_local("evil.example.com", kPort));
        // Absent Host is illegal in HTTP/1.1 and must not be trusted.
        CHECK(!host_header_is_local("", kPort));
    });

    suite.add("substring and suffix tricks are rejected", [] {
        // These are why the comparison is exact rather than a `find`.
        CHECK(!host_header_is_local("127.0.0.1.evil.com:8080", kPort));
        CHECK(!host_header_is_local("evil.com/127.0.0.1", kPort));
        CHECK(!host_header_is_local("notlocalhost:8080", kPort));
        CHECK(!host_header_is_local("localhost.evil.com:8080", kPort));
        CHECK(!host_header_is_local("xlocalhost", kPort));
    });

    suite.add("a mismatched port is rejected", [] {
        // Another service on the same loopback interface is a different origin.
        CHECK(!host_header_is_local("127.0.0.1:9999", kPort));
        CHECK(!host_header_is_local("localhost:80", kPort));
        CHECK(!host_header_is_local("[::1]:9999", kPort));
        // Prefix of the right port must not pass.
        CHECK(!host_header_is_local("127.0.0.1:808", kPort));
        CHECK(!host_header_is_local("127.0.0.1:80800", kPort));
    });

    suite.add("malformed Host values are rejected, not mis-parsed", [] {
        CHECK(!host_header_is_local("[::1:8080", kPort));   // unclosed bracket
        CHECK(!host_header_is_local("::1:8080", kPort));    // unbracketed IPv6
        CHECK(!host_header_is_local("127.0.0.1:", kPort));  // empty port
        CHECK(!host_header_is_local(":8080", kPort));       // empty name
        CHECK(!host_header_is_local("[]:8080", kPort));     // empty literal
    });

    suite.add("Origin must name this same server", [] {
        CHECK(origin_is_self("http://127.0.0.1:8080", kPort));
        CHECK(origin_is_self("http://localhost:8080", kPort));
        CHECK(origin_is_self("http://[::1]:8080", kPort));
    });

    suite.add("cross-origin and opaque Origins are refused", [] {
        CHECK(!origin_is_self("http://evil.example.com", kPort));
        CHECK(!origin_is_self("http://127.0.0.1:9999", kPort));
        // Sandboxed iframes and some privacy modes send the literal "null".
        CHECK(!origin_is_self("null", kPort));
        CHECK(!origin_is_self("", kPort));
        // https to a loopback http listener is still a different origin.
        CHECK(!origin_is_self("https://127.0.0.1:8080", kPort));
        // An Origin carries no path; anything with one is malformed.
        CHECK(!origin_is_self("http://127.0.0.1:8080/", kPort));
        CHECK(!origin_is_self("http://evil.com/http://127.0.0.1:8080", kPort));
    });

    // ---------------------------------------------------------- H2 regressions
    //
    // SECURITY_HTTP H2: an empty data_root used to mean "no confinement", so
    // the SHIPPED DEFAULT was an unauthenticated arbitrary file read. These
    // pin the property that the unconfined state cannot be reached by omission.
    suite.add("H2: serve refuses to start without a data root", [] {
        ServerOptions bare;                       // exactly the old default
        CHECK(bare.data_root.empty());
        CHECK(!bare.confinement_disabled);
        CHECK(rejects(bare));
    });

    suite.add("H2: unconfined requires saying so explicitly", [] {
        ServerOptions opted_out;
        opted_out.confinement_disabled = true;
        CHECK(!rejects(opted_out));               // allowed, but only on request

        // And the two ways of expressing intent cannot be combined, so a
        // wrapper that sets both does not get a silent winner.
        ServerOptions both;
        both.data_root = std::filesystem::temp_directory_path().string();
        both.confinement_disabled = true;
        CHECK(rejects(both));
    });

    suite.add("H2: a data root must exist and be a directory", [] {
        ServerOptions missing;
        missing.data_root = "/nonexistent/path/that/should/not/be/here";
        CHECK(rejects(missing));

        // A regular file is not a jail. Accepting one would confine every
        // request to a single inode and reject everything -- fail-closed, but
        // baffling to debug.
        const auto tmp = std::filesystem::temp_directory_path() /
                         "pto_h2_probe.txt";
        { std::ofstream out(tmp); out << "x"; }
        ServerOptions file_root;
        file_root.data_root = tmp.string();
        CHECK(rejects(file_root));
        std::filesystem::remove(tmp);
    });

    // ------------------------------------------------------- H1(c) regressions
    //
    // The session token. Host and Origin validation bind the BROWSER threat
    // model -- they work because a browser controls those fields and will not
    // let script forge them. A hostile process on the same machine sends
    // whatever headers it likes, and until the token existed nothing stopped
    // it. These pin the configuration half; the comparison itself is
    // tests/test_secrets.cpp's subject.
    suite.add("H1c: the token is on by default and is a real secret", [] {
        ServerOptions opts;
        CHECK(opts.require_token);          // reachable by omission: the SAFE state

        opts.data_root = std::filesystem::temp_directory_path().string();
        HttpServer server(std::move(opts));
        // 32 bytes of CSPRNG rendered as fixed-width hex.
        CHECK_EQ(server.session_token().size(), std::size_t{64});
    });

    suite.add("H1c: each run mints its own token", [] {
        // A token that survived a restart would keep authenticating a URL the
        // operator pasted into a chat months ago.
        const auto root = std::filesystem::temp_directory_path().string();
        ServerOptions a;
        a.data_root = root;
        ServerOptions b;
        b.data_root = root;
        HttpServer first(std::move(a));
        HttpServer second(std::move(b));
        CHECK(first.session_token() != second.session_token());
    });

    suite.add("H1c: --no-auth mints nothing", [] {
        ServerOptions opts;
        opts.data_root = std::filesystem::temp_directory_path().string();
        opts.require_token = false;
        HttpServer server(std::move(opts));
        // Empty, and constant_time_equals("", "") is true -- which is exactly
        // why the pre-routing handler skips the check outright when the token
        // is disabled rather than comparing against an empty string.
        CHECK(server.session_token().empty());
    });

    suite.add("H1c: --allow-remote refuses to run unauthenticated", [] {
        // Binding a file reader to a public interface behind nothing but a
        // printf warning is not a combination worth supporting. On validate()
        // rather than in main.cpp's parser, so a library embedding or a second
        // entry point cannot reach it either -- the lesson of REVIEW finding 2.
        ServerOptions remote_no_auth;
        remote_no_auth.data_root = std::filesystem::temp_directory_path().string();
        remote_no_auth.loopback_only = false;
        remote_no_auth.require_token = false;
        CHECK(rejects(remote_no_auth));

        // With the token it is permitted -- deliberate, warned about, and
        // authenticated.
        ServerOptions remote_with_auth;
        remote_with_auth.data_root = std::filesystem::temp_directory_path().string();
        remote_with_auth.loopback_only = false;
        CHECK(!rejects(remote_with_auth));
    });

    suite.add("H1c: --allow-remote refuses to run unconfined", [] {
        // Previously enforced in main.cpp only, so nothing stopped a second
        // caller constructing it.
        ServerOptions opts;
        opts.confinement_disabled = true;
        opts.loopback_only = false;
        CHECK(rejects(opts));
    });

    suite.add("H2: a valid data root is accepted and canonicalised", [] {
        const auto dir = std::filesystem::temp_directory_path();
        ServerOptions opts;
        // Deliberately non-canonical: trailing "/." must not defeat the
        // confinement comparison, which canonicalises the request path.
        opts.data_root = (dir / ".").string();
        CHECK(!rejects(opts));
        opts.validate();
        CHECK(opts.data_root == std::filesystem::weakly_canonical(dir).string());
        CHECK(!opts.data_root.empty());
    });

    return suite.run();
}
