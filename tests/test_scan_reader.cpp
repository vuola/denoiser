// Tests of the scan reader (include/scan_reader.hpp).
//
// Small voxel files with known content are written into build/test_tmp/; the
// tracked test scans in tests/data/ are read as well. The tests run from the
// repository root.

#include "doctest/doctest.h"

#include "scan_reader.hpp"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using denoiser::ByteOrder;
using denoiser::Region;
using denoiser::SampleType;
using denoiser::ScanInfo;
using denoiser::ScanStatistics;
using denoiser::VolumeSize;

namespace {

const std::string kTmp = "build/test_tmp";

/// Writes bytes to a file below build/test_tmp/ and returns its path.
std::string writeFile(const std::string& name, const std::vector<unsigned char>& bytes) {
    std::filesystem::create_directories(std::filesystem::path(kTmp + "/" + name).parent_path());
    const std::string path = kTmp + "/" + name;
    std::ofstream file(path, std::ios::binary);
    file.write(reinterpret_cast<const char*>(bytes.data()), static_cast<std::streamsize>(bytes.size()));
    return path;
}

/// Writes a text file below build/test_tmp/ and returns its path.
std::string writeText(const std::string& name, const std::string& text) {
    return writeFile(name, std::vector<unsigned char>(text.begin(), text.end()));
}

/// The value the test volumes hold at voxel (x, y, z): different for every
/// voxel of a volume up to 10 x 10 x 10, and above 255 so that both bytes
/// of a 16-bit sample matter.
unsigned valueAt(int x, int y, int z) {
    return 1000u + 100u * static_cast<unsigned>(z) + 10u * static_cast<unsigned>(y) +
           static_cast<unsigned>(x);
}

/// Returns the bytes of a 16-bit volume holding valueAt(), after `header`
/// bytes of the value 0xEE.
std::vector<unsigned char> volume16(const VolumeSize& size, bool little, int header = 0) {
    std::vector<unsigned char> bytes(static_cast<std::size_t>(header), 0xEE);
    for (int z = 0; z < size.depth; ++z) {
        for (int y = 0; y < size.height; ++y) {
            for (int x = 0; x < size.width; ++x) {
                const unsigned v = valueAt(x, y, z);
                const unsigned char low = static_cast<unsigned char>(v & 0xFF);
                const unsigned char high = static_cast<unsigned char>(v >> 8);
                bytes.push_back(little ? low : high);
                bytes.push_back(little ? high : low);
            }
        }
    }
    return bytes;
}

/// The .vgi text of a lab scan, shortened: Windows line ends, indented keys,
/// a second "Size" in mixed case, and keys that the reader ignores.
const char* const kLabVgi =
    "{volume1}\r\n"
    "[representation]\r\n"
    "size = 5 4 3\r\n"
    "bitsperelement = 16\r\n"
    "datatype = unsigned integer\r\n"
    "datarange = 0 65535\r\n"
    "[file1]\r\n"
    "Mirror = 0 1 1 \r\n"
    "  RegionOfInterestStart = 0 0 0 \r\n"
    " FileFormat = raw\r\n"
    "Size = 99 99 99\r\n"
    "Name = X:\\somewhere\\else\\scan.vol\r\n"
    "{volumeprimitive9}\r\n"
    "[section1]\r\n"
    " offset = 0 \r\n"
    "[geometry]\r\n"
    "resolution = 0.005000 0.005000 0.005000\r\n"
    "unit = mm\r\n";

}  // namespace

TEST_CASE("parseVgi reads a lab-style file") {
    const ScanInfo info = denoiser::parseVgi(kLabVgi, "test.vgi");
    CHECK(info.size.width == 5);
    CHECK(info.size.height == 4);
    CHECK(info.size.depth == 3);
    CHECK(info.sample_type == SampleType::kUInt16);
    CHECK(info.byte_order == ByteOrder::kLittleEndian);
    CHECK(info.header_bytes == 0);
    CHECK(info.voxel_size[0] == doctest::Approx(0.005));
    CHECK(info.voxel_size[1] == doctest::Approx(0.005));
    CHECK(info.voxel_size[2] == doctest::Approx(0.005));
    CHECK(info.unit == "mm");
    CHECK(info.volume_path.empty());
}

TEST_CASE("parseVgi: sample types, header and optional keys") {
    SUBCASE("8-bit, no data type, no resolution") {
        const ScanInfo info = denoiser::parseVgi("size = 2 3 4\nbitsperelement = 8\n", "t");
        CHECK(info.sample_type == SampleType::kUInt8);
        CHECK(info.voxel_size[0] == 0.0);
        CHECK(info.unit.empty());
    }
    SUBCASE("float with a header") {
        const ScanInfo info = denoiser::parseVgi(
            "size = 2 3 4\nbitsperelement = 32\ndatatype = Float\nSkipHeader = 128\n", "t");
        CHECK(info.sample_type == SampleType::kFloat32);
        CHECK(info.header_bytes == 128);
    }
}

TEST_CASE("parseVgi rejects what it cannot read") {
    using denoiser::parseVgi;
    CHECK_THROWS_AS(parseVgi("bitsperelement = 16\n", "t"), std::runtime_error);
    CHECK_THROWS_AS(parseVgi("size = 2 3 4\n", "t"), std::runtime_error);
    CHECK_THROWS_AS(parseVgi("size = 2 3\nbitsperelement = 16\n", "t"), std::runtime_error);
    CHECK_THROWS_AS(parseVgi("size = 2 3 4 5\nbitsperelement = 16\n", "t"), std::runtime_error);
    CHECK_THROWS_AS(parseVgi("size = 2 0 4\nbitsperelement = 16\n", "t"), std::runtime_error);
    CHECK_THROWS_AS(parseVgi("size = 2 x 4\nbitsperelement = 16\n", "t"), std::runtime_error);
    CHECK_THROWS_AS(parseVgi("size = 2 3 4\nbitsperelement = 12\n", "t"), std::runtime_error);
    CHECK_THROWS_AS(parseVgi("size = 2 3 4\nbitsperelement = 16\ndatatype = signed integer\n", "t"),
                    std::runtime_error);
    CHECK_THROWS_AS(parseVgi("size = 2 3 4\nbitsperelement = 16\ndatatype = float\n", "t"),
                    std::runtime_error);
    CHECK_THROWS_AS(parseVgi("size = 2 3 4\nbitsperelement = 16\nskipheader = -1\n", "t"),
                    std::runtime_error);
    CHECK_THROWS_AS(parseVgi("size = 2 3 4\nbitsperelement = 16\nresolution = 1 1\n", "t"),
                    std::runtime_error);
    // The message names the file.
    try {
        parseVgi("", "my_scan.vgi");
        FAIL("no exception");
    } catch (const std::runtime_error& error) {
        CHECK(std::string(error.what()).find("my_scan.vgi") != std::string::npos);
    }
}

TEST_CASE("openScan finds the scan from a folder, a .vgi file or a voxel file") {
    const VolumeSize size{5, 4, 3};
    writeFile("scan_a/scan_a.vol", volume16(size, true));
    writeText("scan_a/scan_a.vgi", kLabVgi);

    for (const std::string& path : {kTmp + "/scan_a", kTmp + "/scan_a/scan_a.vgi",
                                   kTmp + "/scan_a/scan_a.vol"}) {
        CAPTURE(path);
        const ScanInfo info = denoiser::openScan(path);
        CHECK(info.volume_path == kTmp + "/scan_a/scan_a.vol");
        CHECK(info.size.width == 5);
        const std::vector<float> values = denoiser::readVolume(info);
        REQUIRE(values.size() == 60);
        CHECK(values[0] == 1000.0f);
        CHECK(values[(2 * 4 + 3) * 5 + 4] == static_cast<float>(valueAt(4, 3, 2)));
    }
}

TEST_CASE("openScan reports missing and mismatching files") {
    const VolumeSize size{5, 4, 3};

    CHECK_THROWS_AS(denoiser::openScan(kTmp + "/does_not_exist"), std::runtime_error);

    std::filesystem::create_directories(kTmp + "/empty_folder");
    CHECK_THROWS_AS(denoiser::openScan(kTmp + "/empty_folder"), std::runtime_error);

    writeText("two_vgi/a.vgi", kLabVgi);
    writeText("two_vgi/b.vgi", kLabVgi);
    CHECK_THROWS_AS(denoiser::openScan(kTmp + "/two_vgi"), std::runtime_error);

    writeText("no_vol/scan.vgi", kLabVgi);
    CHECK_THROWS_AS(denoiser::openScan(kTmp + "/no_vol"), std::runtime_error);

    // One byte too few.
    std::vector<unsigned char> bytes = volume16(size, true);
    bytes.pop_back();
    writeFile("short/scan.vol", bytes);
    writeText("short/scan.vgi", kLabVgi);
    try {
        denoiser::openScan(kTmp + "/short");
        FAIL("no exception");
    } catch (const std::runtime_error& error) {
        const std::string message = error.what();
        CHECK(message.find("119 bytes") != std::string::npos);
        CHECK(message.find("need 120") != std::string::npos);
    }
}

TEST_CASE("readRegion returns the right voxels for every sample type and byte order") {
    const VolumeSize size{5, 4, 3};
    const Region region{1, 2, 1, 3, 2, 2};

    SUBCASE("16-bit little-endian") {
        const std::string path = writeFile("u16le.raw", volume16(size, true));
        const ScanInfo info = denoiser::describeRaw(path, size, SampleType::kUInt16,
                                                    ByteOrder::kLittleEndian);
        const std::vector<float> values = denoiser::readRegion(info, region);
        REQUIRE(values.size() == 12);
        std::size_t i = 0;
        for (int z = 1; z < 3; ++z) {
            for (int y = 2; y < 4; ++y) {
                for (int x = 1; x < 4; ++x) {
                    CHECK(values[i++] == static_cast<float>(valueAt(x, y, z)));
                }
            }
        }
    }
    SUBCASE("16-bit big-endian with a header") {
        const std::string path = writeFile("u16be.raw", volume16(size, false, 7));
        const ScanInfo info = denoiser::describeRaw(path, size, SampleType::kUInt16,
                                                    ByteOrder::kBigEndian, 7);
        const std::vector<float> values = denoiser::readRegion(info, region);
        REQUIRE(values.size() == 12);
        CHECK(values.front() == static_cast<float>(valueAt(1, 2, 1)));
        CHECK(values.back() == static_cast<float>(valueAt(3, 3, 2)));
    }
    SUBCASE("8-bit") {
        std::vector<unsigned char> bytes;
        for (int i = 0; i < 60; ++i) {
            bytes.push_back(static_cast<unsigned char>(100 + i));
        }
        const std::string path = writeFile("u8.raw", bytes);
        const ScanInfo info = denoiser::describeRaw(path, size, SampleType::kUInt8,
                                                    ByteOrder::kLittleEndian);
        const std::vector<float> values = denoiser::readRegion(info, region);
        CHECK(values.front() == static_cast<float>(100 + (1 * 4 + 2) * 5 + 1));
        CHECK(values.back() == static_cast<float>(100 + (2 * 4 + 3) * 5 + 3));
    }
    SUBCASE("32-bit float, both byte orders") {
        for (const bool little : {true, false}) {
            CAPTURE(little);
            std::vector<unsigned char> bytes;
            for (int i = 0; i < 60; ++i) {
                const float value = -1.5f + 0.25f * static_cast<float>(i);
                std::uint32_t bits = 0;
                std::memcpy(&bits, &value, sizeof(bits));
                for (int b = 0; b < 4; ++b) {
                    const int shift = little ? 8 * b : 8 * (3 - b);
                    bytes.push_back(static_cast<unsigned char>((bits >> shift) & 0xFF));
                }
            }
            const std::string path = writeFile(little ? "f32le.raw" : "f32be.raw", bytes);
            const ScanInfo info = denoiser::describeRaw(
                path, size, SampleType::kFloat32,
                little ? ByteOrder::kLittleEndian : ByteOrder::kBigEndian);
            const std::vector<float> values = denoiser::readVolume(info);
            REQUIRE(values.size() == 60);
            CHECK(values[0] == -1.5f);
            CHECK(values[59] == -1.5f + 0.25f * 59.0f);
        }
    }
}

TEST_CASE("readRegion: full-width regions and regions that do not fit") {
    const VolumeSize size{5, 4, 3};
    const std::string path = writeFile("fit.raw", volume16(size, true));
    const ScanInfo info = denoiser::describeRaw(path, size, SampleType::kUInt16,
                                                ByteOrder::kLittleEndian);

    // Full width, part of the height: the path that reads a slice at once.
    const std::vector<float> rows = denoiser::readRegion(info, Region{0, 1, 2, 5, 2, 1});
    REQUIRE(rows.size() == 10);
    CHECK(rows.front() == static_cast<float>(valueAt(0, 1, 2)));
    CHECK(rows.back() == static_cast<float>(valueAt(4, 2, 2)));

    // The last voxel alone.
    const std::vector<float> last = denoiser::readRegion(info, Region{4, 3, 2, 1, 1, 1});
    REQUIRE(last.size() == 1);
    CHECK(last[0] == static_cast<float>(valueAt(4, 3, 2)));

    CHECK_THROWS_AS(denoiser::readRegion(info, Region{0, 0, 0, 6, 4, 3}), std::runtime_error);
    CHECK_THROWS_AS(denoiser::readRegion(info, Region{1, 0, 0, 5, 4, 3}), std::runtime_error);
    CHECK_THROWS_AS(denoiser::readRegion(info, Region{0, 0, 3, 1, 1, 1}), std::runtime_error);
    CHECK_THROWS_AS(denoiser::readRegion(info, Region{-1, 0, 0, 2, 2, 2}), std::runtime_error);
    CHECK_THROWS_AS(denoiser::readRegion(info, Region{0, 0, 0, 0, 4, 3}), std::runtime_error);
}

TEST_CASE("describeRaw checks the description against the file") {
    const VolumeSize size{5, 4, 3};
    const std::string path = writeFile("check.raw", volume16(size, true));

    CHECK_NOTHROW(denoiser::describeRaw(path, size, SampleType::kUInt16, ByteOrder::kLittleEndian));
    // The same bytes also fit 10 x 4 x 3 voxels of 8 bits.
    CHECK_NOTHROW(denoiser::describeRaw(path, VolumeSize{10, 4, 3}, SampleType::kUInt8,
                                        ByteOrder::kLittleEndian));
    CHECK_THROWS_AS(denoiser::describeRaw(path, VolumeSize{5, 4, 4}, SampleType::kUInt16,
                                          ByteOrder::kLittleEndian),
                    std::runtime_error);
    CHECK_THROWS_AS(denoiser::describeRaw(path, size, SampleType::kUInt16,
                                          ByteOrder::kLittleEndian, 2),
                    std::runtime_error);
    CHECK_THROWS_AS(denoiser::describeRaw(path, size, SampleType::kUInt16,
                                          ByteOrder::kLittleEndian, -1),
                    std::runtime_error);
    CHECK_THROWS_AS(denoiser::describeRaw(path, VolumeSize{0, 4, 3}, SampleType::kUInt16,
                                          ByteOrder::kLittleEndian),
                    std::runtime_error);
    CHECK_THROWS_AS(denoiser::describeRaw(kTmp + "/missing.raw", size, SampleType::kUInt16,
                                          ByteOrder::kLittleEndian),
                    std::runtime_error);
}

TEST_CASE("scanStatistics") {
    SUBCASE("known values") {
        // 2 x 2 x 2 voxels: 1, 2, ..., 8. Mean 4.5, variance 5.25.
        std::vector<unsigned char> bytes = {1, 2, 3, 4, 5, 6, 7, 8};
        const std::string path = writeFile("stats.raw", bytes);
        const ScanStatistics stats = denoiser::scanStatistics(denoiser::describeRaw(
            path, VolumeSize{2, 2, 2}, SampleType::kUInt8, ByteOrder::kLittleEndian));
        CHECK(stats.voxels == 8);
        CHECK(stats.non_finite == 0);
        CHECK(stats.minimum == 1.0);
        CHECK(stats.maximum == 8.0);
        CHECK(stats.mean == doctest::Approx(4.5));
        CHECK(stats.std_dev == doctest::Approx(std::sqrt(5.25)));
    }
    SUBCASE("NaN and infinity are counted and left out") {
        const float values[4] = {2.0f, std::numeric_limits<float>::quiet_NaN(), 4.0f,
                                 std::numeric_limits<float>::infinity()};
        std::vector<unsigned char> bytes(sizeof(values));
        std::memcpy(bytes.data(), values, sizeof(values));
        // The bytes are in the order of this machine; describe them so.
        const std::uint16_t probe = 1;
        unsigned char first_byte = 0;
        std::memcpy(&first_byte, &probe, 1);
        const std::string path = writeFile("stats_nan.raw", bytes);
        const ScanStatistics stats = denoiser::scanStatistics(denoiser::describeRaw(
            path, VolumeSize{4, 1, 1}, SampleType::kFloat32,
            first_byte == 1 ? ByteOrder::kLittleEndian : ByteOrder::kBigEndian));
        CHECK(stats.voxels == 4);
        CHECK(stats.non_finite == 2);
        CHECK(stats.minimum == 2.0);
        CHECK(stats.maximum == 4.0);
        CHECK(stats.mean == doctest::Approx(3.0));
        CHECK(stats.std_dev == doctest::Approx(1.0));
    }
}

TEST_CASE("the tracked test scans are read correctly") {
    // Sum, extremes and single voxels were computed independently with
    // Python's array module from the .vol files.
    struct Expected {
        const char* scan;
        double sum;
        double minimum;
        double maximum;
        float first;
        float at_3_7_5;
        float last;
    };
    const Expected expected[] = {
        {"Sample1_18p15_5mu_c64", 8619772928.0, 11048.0, 53022.0, 29333.0f, 33520.0f, 39858.0f},
        {"Sample2_5mu_c64", 7245312265.0, 9442.0, 46217.0, 31139.0f, 20215.0f, 30361.0f},
    };
    for (const Expected& e : expected) {
        CAPTURE(e.scan);
        const ScanInfo info = denoiser::openScan(std::string("tests/data/") + e.scan);
        CHECK(info.size.width == 64);
        CHECK(info.size.height == 64);
        CHECK(info.size.depth == 64);
        CHECK(info.sample_type == SampleType::kUInt16);
        CHECK(info.voxel_size[2] == doctest::Approx(0.005));
        CHECK(info.unit == "mm");

        const std::vector<float> values = denoiser::readVolume(info);
        REQUIRE(values.size() == 64u * 64u * 64u);
        CHECK(values.front() == e.first);
        CHECK(values[(5 * 64 + 7) * 64 + 3] == e.at_3_7_5);
        CHECK(values.back() == e.last);

        // A region gives the same voxels as the whole volume.
        const std::vector<float> one = denoiser::readRegion(info, Region{3, 7, 5, 1, 1, 1});
        CHECK(one[0] == e.at_3_7_5);

        const ScanStatistics stats = denoiser::scanStatistics(info);
        CHECK(stats.minimum == e.minimum);
        CHECK(stats.maximum == e.maximum);
        CHECK(stats.mean == doctest::Approx(e.sum / (64.0 * 64.0 * 64.0)).epsilon(1e-12));
    }
}
