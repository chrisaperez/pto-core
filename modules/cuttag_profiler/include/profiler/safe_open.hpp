// SPDX-License-Identifier: MIT
//
// Opening an untrusted path exactly once.
//
// SECURITY_HTTP_2026-08-15 finding M6: the HTTP handlers used to check a path
// with `fs::exists()` and then hand the *name* to htslib, which looked it up
// again. Two lookups of a mutable namespace with a gap in between, so anything
// that could write inside the data root could pass the check with a real file
// and be read through a symlink -- the canonicalisation that would have caught
// it having already run. The data root is exactly the directory a pipeline
// drops files into, and those are often group-writable.
//
// The fix is to stop checking a name. `open_regular_file` performs one
// `open()` and asks every subsequent question of the returned descriptor, so
// there is no second lookup to race. Callers that can consume a descriptor
// (see `read_regions_from_fd`) then never touch the name again.
#pragma once

#include <string>
#include <utility>

namespace profiler {

// Owning file descriptor. Small enough to hand-roll; a unique_ptr with a
// custom deleter would have to represent -1 as a valid state.
class UniqueFd {
public:
    UniqueFd() noexcept = default;
    explicit UniqueFd(int fd) noexcept : fd_(fd) {}
    ~UniqueFd() { reset(); }

    UniqueFd(const UniqueFd&) = delete;
    UniqueFd& operator=(const UniqueFd&) = delete;

    UniqueFd(UniqueFd&& other) noexcept : fd_(std::exchange(other.fd_, -1)) {}
    UniqueFd& operator=(UniqueFd&& other) noexcept {
        if (this != &other) reset(std::exchange(other.fd_, -1));
        return *this;
    }

    [[nodiscard]] int get() const noexcept { return fd_; }
    [[nodiscard]] explicit operator bool() const noexcept { return fd_ >= 0; }

    // Gives up ownership. Used when handing the descriptor to a library that
    // closes it itself -- htslib's bgzf_hdopen does.
    [[nodiscard]] int release() noexcept { return std::exchange(fd_, -1); }

    void reset(int fd = -1) noexcept;

private:
    int fd_ = -1;
};

// Opens `path` read-only and returns the descriptor, having verified through
// that descriptor -- never through the name -- that it refers to a regular
// file.
//
// Three properties, each of which was a real hole:
//
//   * ONE lookup. Nothing re-resolves `path` afterwards, so the check/open
//     race of M6 has no window to open in.
//   * `O_NOFOLLOW`, so the final component cannot be a symlink. Confinement
//     canonicalises before this point; O_NOFOLLOW is what stops the name being
//     repointed at a symlink *after* that canonicalisation ran.
//   * `S_ISREG`, from `fstat` on the descriptor. A FIFO inside the data root
//     passes every path check and then blocks the calling thread in
//     `bgzf_open` until a writer appears -- a handful of such requests
//     exhausts httplib's fixed thread pool. A device node is never valid input
//     and reading one can block or hang the process just as thoroughly.
//
// `O_CLOEXEC` is set so a descriptor cannot leak into the browser process
// spawned by `open_in_browser`.
//
// Throws std::runtime_error on any failure. The message names the path and is
// therefore only ever safe to show the operator: HTTP callers must funnel it
// through the sanitiser (M4), which they do by catching it as a PathError.
[[nodiscard]] UniqueFd open_regular_file(const std::string& path);

}  // namespace profiler
