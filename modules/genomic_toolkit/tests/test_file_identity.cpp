// SPDX-License-Identifier: MIT
//
// genomic_toolkit: the output-overwrites-input guard (include/toolkit/file_identity.hpp).
//
// Each module that writes files has its own copy of that header and this suite
// (the header says why). The CLI-level shapes -- `-o` naming an input, an
// input redirected into stdin -- are driven through the binary elsewhere; this
// pins the semantics they rest on.
#include <fcntl.h>
#include <unistd.h>

#include <filesystem>
#include <fstream>
#include <random>
#include <string>
#include <system_error>
#include <vector>

#include "test_util.hpp"

#include "toolkit/file_identity.hpp"

namespace {

namespace fs = std::filesystem;
using toolkit::output_overwrites_input;
using toolkit::same_regular_file;

struct TempDir {
    fs::path path;
    TempDir() {
        std::random_device rd;
        path = fs::temp_directory_path() /
               ("fid_" + std::to_string(rd()) + "_" + std::to_string(rd()));
        fs::create_directories(path);
    }
    ~TempDir() {
        std::error_code ec;
        fs::remove_all(path, ec);
    }
};

void write_file(const fs::path& p, const std::string& text) {
    std::ofstream out(p);
    out << text;
}

void test_one_file_under_every_spelling_is_one_file() {
    TempDir d;
    const fs::path f = d.path / "in.dat";
    write_file(f, "x\n");
    fs::create_hard_link(f, d.path / "hard.dat");
    fs::create_symlink(f, d.path / "soft.dat");
    fs::create_directories(d.path / "sub");
    const std::vector<fs::path> aliases = {f, d.path / "hard.dat", d.path / "soft.dat",
                                           d.path / "sub" / ".." / "in.dat",
                                           d.path / "." / "in.dat"};
    for (const fs::path& alias : aliases) {
        CHECK(same_regular_file(f, alias));
        CHECK(output_overwrites_input(alias.string(), {f.string()}).has_value());
    }
}

void test_different_and_missing_files_never_clash() {
    TempDir d;
    const fs::path a = d.path / "a.dat";
    const fs::path b = d.path / "b.dat";
    write_file(a, "same bytes\n");
    write_file(b, "same bytes\n");  // equal content is not identity
    CHECK(!same_regular_file(a, b));
    CHECK(!output_overwrites_input(b.string(), {a.string()}).has_value());
    // An output that does not exist yet cannot overwrite anything.
    CHECK(!output_overwrites_input((d.path / "new.out").string(), {a.string()}).has_value());
    // Nor can a missing input be overwritten; the open reports that one.
    CHECK(!output_overwrites_input(a.string(), {(d.path / "gone.dat").string()}).has_value());
    fs::create_symlink(d.path / "gone.dat", d.path / "dangling.dat");
    CHECK(!same_regular_file(d.path / "dangling.dat", d.path / "dangling.dat"));
}

void test_non_regular_files_never_clash() {
    TempDir d;
    CHECK(!same_regular_file("/dev/null", "/dev/null"));
    CHECK(!output_overwrites_input("/dev/null", {"/dev/null"}).has_value());
    CHECK(!same_regular_file(d.path, d.path));
}

void test_unset_inputs_are_skipped_and_the_clash_is_named() {
    TempDir d;
    const fs::path other = d.path / "other.dat";
    const fs::path f = d.path / "in.dat";
    write_file(other, "x\n");
    write_file(f, "y\n");
    const auto clash = output_overwrites_input(f.string(), {other.string(), "", f.string()});
    CHECK(clash.has_value());
    CHECK(clash.value_or("") == f.string());
    CHECK(!output_overwrites_input("", {f.string()}).has_value());
}

// `tool -o x < x`: "-" must resolve through the descriptor. A /dev/stdin path
// does not reach the redirected file on macOS, which is how this was first
// implemented and how it first failed.
void test_dash_is_asked_of_the_descriptor() {
    TempDir d;
    const fs::path f = d.path / "in.dat";
    const fs::path other = d.path / "other.dat";
    write_file(f, "x\n");
    write_file(other, "y\n");

    const int saved = ::dup(STDIN_FILENO);
    const int fd = ::open(f.c_str(), O_RDONLY);
    CHECK(saved >= 0);
    CHECK(fd >= 0);
    if (saved < 0 || fd < 0) return;
    ::dup2(fd, STDIN_FILENO);
    ::close(fd);
    const bool caught = output_overwrites_input(f.string(), {"-"}).has_value();
    const bool other_clear = !output_overwrites_input(other.string(), {"-"}).has_value();
    ::dup2(saved, STDIN_FILENO);
    ::close(saved);

    CHECK(caught);
    CHECK(other_clear);
}

}  // namespace

int main() {
    testing::Suite suite{"gtk/file_identity", {}};
    suite.add("one file under every spelling is one file",
              test_one_file_under_every_spelling_is_one_file);
    suite.add("different and missing files never clash",
              test_different_and_missing_files_never_clash);
    suite.add("non-regular files never clash", test_non_regular_files_never_clash);
    suite.add("unset inputs are skipped and the clash is named",
              test_unset_inputs_are_skipped_and_the_clash_is_named);
    suite.add("dash is asked of the descriptor", test_dash_is_asked_of_the_descriptor);
    return suite.run();
}
