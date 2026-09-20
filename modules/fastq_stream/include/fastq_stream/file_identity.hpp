// SPDX-License-Identifier: MIT
// fastq_stream — would writing an output path destroy an input?
//
// `fastq_stream -i reads.fq -o reads.fq` opened the output with O_TRUNC while
// the reader had staged 1 MiB, and left a zero-byte file where 42 MB of reads
// had been (docs/TORTURE_2026-09-10.md F1). The same defect, one module over,
// overwrote pto-peaks' input BED and --chrom-sizes file with narrowPeak rows.
//
// ONE HEADER PER MODULE, deliberately. An identical copy of this file lives in
// fastq_stream, peaks, genomic_toolkit and cuttag_profiler, each in that
// module's own namespace and exercised by that module's own tests.
// pto-core/CLAUDE.md forbids a shared cross-module layer (modules share no
// code, and each must stay independently buildable); a guard this size is not
// the thing to break that for. If you change the semantics here, change all
// four -- grep for `output_overwrites_input`.
//
// Semantics
// ---------
// Named paths are compared with std::filesystem::equivalent, i.e. (device,
// inode) after following links: a hard link, a symlink, `./x` vs `x` and a path
// through `..` are all the same file. Only paths that EXIST and are REGULAR
// FILES are compared. A missing output cannot clobber anything; and two names
// for one terminal, pipe or /dev/null are not a hazard -- refusing them would
// break `tool < /dev/tty` and writing to /dev/null.
//
// "-" is a standard stream, and it is asked of the DESCRIPTOR (fstat). That is
// what catches `tool -o x < x`, which the shell does not truncate, and
// `tool --stdout >> x` with x an input. std::filesystem has no descriptor API,
// and the path-shaped substitute was tried first: /dev/stdin does not resolve
// to the redirected file on macOS, so `pto-peaks -o x < x` passed that check
// and overwrote x. The peaks CLI test is what caught it.
#pragma once

#include <sys/stat.h>
#include <unistd.h>

#include <filesystem>
#include <optional>
#include <string>
#include <system_error>
#include <utility>
#include <vector>

namespace fq {

// True when both paths exist, are regular files after following links, and are
// the same file. Never throws for a path it cannot stat: an unreadable path is
// not provably the same file, and the open that follows reports it.
[[nodiscard]] inline bool same_regular_file(const std::filesystem::path& a,
                                            const std::filesystem::path& b) {
  std::error_code ec;
  if (!std::filesystem::is_regular_file(a, ec) || ec) return false;
  if (!std::filesystem::is_regular_file(b, ec) || ec) return false;
  const bool same = std::filesystem::equivalent(a, b, ec);
  return !ec && same;
}

// (device, inode) of the regular file a CLI argument names, or nullopt. "-" is
// stdin for an input and stdout for an output, asked of the descriptor.
[[nodiscard]] inline std::optional<std::pair<dev_t, ino_t>> regular_file_id(
    const std::string& arg, bool output) {
  struct stat st {};
  const int rc = arg == "-" ? ::fstat(output ? STDOUT_FILENO : STDIN_FILENO, &st)
                            : ::stat(arg.c_str(), &st);
  if (rc != 0 || !S_ISREG(st.st_mode)) return std::nullopt;
  return std::make_pair(st.st_dev, st.st_ino);
}

// Whether writing `output` would overwrite `input`.
[[nodiscard]] inline bool overwrites(const std::string& output, const std::string& input) {
  if (output == "-" || input == "-") {
    const auto out = regular_file_id(output, /*output=*/true);
    const auto in = regular_file_id(input, /*output=*/false);
    return out.has_value() && in.has_value() && *out == *in;
  }
  return same_regular_file(output, input);
}

// The first of `inputs` that writing `output` would overwrite, or nullopt.
// Empty strings stand for options the user did not pass and are skipped.
[[nodiscard]] inline std::optional<std::string> output_overwrites_input(
    const std::string& output, const std::vector<std::string>& inputs) {
  if (output.empty()) return std::nullopt;
  for (const std::string& in : inputs) {
    if (in.empty()) continue;
    if (overwrites(output, in)) return in;
  }
  return std::nullopt;
}

}  // namespace fq
