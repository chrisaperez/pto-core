// SPDX-License-Identifier: MIT
#include "profiler/safe_open.hpp"

#include <cerrno>
#include <cstring>
#include <stdexcept>

#if defined(_WIN32)
#include <io.h>
#include <sys/stat.h>
#include <sys/types.h>
#else
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#endif

namespace profiler {
namespace {

std::string describe_errno(int saved) {
    // strerror is not thread-safe and strerror_r has two incompatible
    // signatures across platforms; the numeric code plus the caller's context
    // is enough for an operator log line and cannot be misformatted.
    return "errno " + std::to_string(saved);
}

}  // namespace

void UniqueFd::reset(int fd) noexcept {
    if (fd_ >= 0) {
#if defined(_WIN32)
        ::_close(fd_);
#else
        // close() can report EINTR, but on every platform this targets the
        // descriptor is released regardless, and retrying risks closing a
        // descriptor another thread has since been handed. Do not loop.
        ::close(fd_);
#endif
    }
    fd_ = fd;
}

UniqueFd open_regular_file(const std::string& path) {
#if defined(_WIN32)
    // No O_NOFOLLOW: Windows reparse points are not followed by _open in the
    // first place for the cases that matter here, and the symlink half of M6
    // is a POSIX-shared-host concern. The S_ISREG check below still applies.
    const int fd = ::_open(path.c_str(), _O_RDONLY | _O_BINARY);
#else
    // O_NONBLOCK is load-bearing, not a tuning choice. Opening a FIFO
    // O_RDONLY *blocks in open() itself* until a writer appears -- so without
    // it this function hangs on the exact input the S_ISREG check below exists
    // to reject, and the check never runs. With it, the open returns
    // immediately, fstat says S_ISFIFO, and the request is refused. It is
    // cleared again once the descriptor is known to be a regular file.
    const int fd =
        ::open(path.c_str(), O_RDONLY | O_CLOEXEC | O_NOFOLLOW | O_NONBLOCK);
#endif
    if (fd < 0) {
        const int saved = errno;
        // ELOOP here means the final component is a symlink, which O_NOFOLLOW
        // refused. That is a rejection, not a missing file, and it is the case
        // this function exists for -- so say which it was in the operator log.
        throw std::runtime_error(
            (saved == ELOOP ? "refusing to open a symlink: "
                            : "cannot open file: ") +
            path + " (" + describe_errno(saved) + ")");
    }

    UniqueFd owned(fd);

#if defined(_WIN32)
    struct _stat64 st {};
    const int rc = ::_fstat64(owned.get(), &st);
    const bool regular = (st.st_mode & _S_IFMT) == _S_IFREG;
#else
    struct stat st {};
    const int rc = ::fstat(owned.get(), &st);
    const bool regular = S_ISREG(st.st_mode);
#endif
    if (rc != 0) {
        throw std::runtime_error("cannot stat open file: " + path + " (" +
                                 describe_errno(errno) + ")");
    }
    // Asked of the DESCRIPTOR, not the name. fs::is_regular_file(path) would
    // be a second lookup and could disagree with what was actually opened --
    // which is the whole of M6.
    if (!regular) {
        throw std::runtime_error("not a regular file: " + path);
    }

#if !defined(_WIN32)
    // Drop O_NONBLOCK now that the descriptor is known to be a regular file.
    // POSIX says the flag has no effect on regular-file reads, but htslib is
    // handed this descriptor and nothing here needs to depend on that: leaving
    // a flag set that the consumer did not ask for is how a short read starts
    // being interpreted as end of file.
    const int flags = ::fcntl(owned.get(), F_GETFL);
    if (flags >= 0) {
        (void)::fcntl(owned.get(), F_SETFL, flags & ~O_NONBLOCK);
    }
#endif
    return owned;
}

}  // namespace profiler
