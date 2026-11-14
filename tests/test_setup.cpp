// Tests of the project setup itself (step 0.1): the library links, and the
// tracked test scans are present and complete.

#include <cstdint>
#include <fstream>
#include <string>

#include "doctest/doctest.h"
#include "version.hpp"

namespace {

/**
 * @brief Size of a file in bytes.
 * @param path File path, relative to the repository root.
 * @return The size, or -1 if the file cannot be opened.
 */
std::int64_t file_size(const std::string& path) {
    std::ifstream file(path, std::ios::binary | std::ios::ate);
    if (!file) {
        return -1;
    }
    return static_cast<std::int64_t>(file.tellg());
}

}  // namespace

TEST_CASE("version is major.minor.patch") {
    const std::string v = denoiser::version();
    int dots = 0;
    for (char c : v) {
        const bool digit = (c >= '0' && c <= '9');
        CHECK((digit || c == '.'));
        if (c == '.') {
            ++dots;
        }
    }
    CHECK(dots == 2);
    CHECK(v.front() != '.');
    CHECK(v.back() != '.');
}

TEST_CASE("tracked test scans are present and complete") {
    // 64 x 64 x 64 voxels, 2 bytes each.
    const std::int64_t expected = 64 * 64 * 64 * 2;

    for (const std::string scan : {"Sample1_18p15_5mu_c64", "Sample2_5mu_c64"}) {
        CAPTURE(scan);
        const std::string base = "tests/data/" + scan + "/" + scan;
        CHECK(file_size(base + ".vol") == expected);
        CHECK(file_size(base + ".vgi") > 0);
    }
}
