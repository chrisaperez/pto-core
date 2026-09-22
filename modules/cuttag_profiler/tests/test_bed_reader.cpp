#include "profiler/bed_reader.hpp"
#include "test_util.hpp"
#include <fstream>
#include <string>
#include <stdexcept>

void test_bounded_line_growth() {
    const std::string path = "huge_line.bed";
    std::ofstream out(path, std::ios::binary);
    // Write 65 MiB of 'A' without a newline
    std::string chunk(1024 * 1024, 'A');
    for (int i = 0; i < 65; ++i) {
        out.write(chunk.data(), chunk.size());
    }
    out.close();

    bool caught = false;
    try {
        profiler::read_bed(path);
    } catch (const std::runtime_error& e) {
        std::string msg = e.what();
        if (msg.find("limit") != std::string::npos) {
            caught = true;
        }
    }
    std::remove(path.c_str());
    CHECK(caught);
}

int main() {
    test_bounded_line_growth();
    return testing::g_failures > 0 ? 1 : 0;
}
