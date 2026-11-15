// Tests of the scan writer (include/scan_writer.hpp).
//
// Files are written into build/test_tmp/writer/. The tests run from the
// repository root.

#include "doctest/doctest.h"

#include "scan_reader.hpp"
#include "scan_writer.hpp"

#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using denoiser::ByteOrder;
using denoiser::Region;
using denoiser::SampleType;
using denoiser::ScanInfo;
using denoiser::VolumeSize;
using denoiser::WriteReport;

namespace {

const std::string kOut = "build/test_tmp/writer";

/// Returns a description of a scan to be written below build/test_tmp/writer/.
ScanInfo describe(const std::string& name, const VolumeSize& size, SampleType type,
                  ByteOrder order = ByteOrder::kLittleEndian) {
    std::filesystem::create_directories(kOut);
    ScanInfo info;
    info.volume_path = kOut + "/" + name;
    info.size = size;
    info.sample_type = type;
    info.byte_order = order;
    return info;
}

/// Returns all bytes of a file.
std::vector<unsigned char> fileBytes(const std::string& path) {
    std::ifstream file(path, std::ios::binary);
    REQUIRE(file.good());
    return std::vector<unsigned char>((std::istreambuf_iterator<char>(file)),
                                      std::istreambuf_iterator<char>());
}

/// Returns width * height * depth values: 1000, 1001, 1002, ...
std::vector<float> counting(const VolumeSize& size) {
    std::vector<float> values(static_cast<std::size_t>(size.width) * size.height * size.depth);
    for (std::size_t i = 0; i < values.size(); ++i) {
        values[i] = 1000.0f + static_cast<float>(i);
    }
    return values;
}

}  // namespace

TEST_CASE("vgiText") {
    ScanInfo info;
    info.volume_path = "some/folder/result.vol";
    info.size = VolumeSize{5, 4, 3};

    SUBCASE("16-bit with voxel size") {
        info.voxel_size[0] = info.voxel_size[1] = info.voxel_size[2] = 0.005;
        info.unit = "mm";
        CHECK(denoiser::vgiText(info) ==
              "{volume1}\n[representation]\nsize = 5 4 3\nbitsperelement = 16\n"
              "datatype = unsigned integer\ndatarange = 0 65535\n"
              "[file1]\nFileFormat = raw\nSize = 5 4 3\nName = result.vol\nbitsperelement = 16\n"
              "datatype = unsigned integer\ndatarange = 0 65535\n"
              "{volumeprimitive9}\n[geometry]\nresolution = 0.005000 0.005000 0.005000\nunit = mm\n");
    }
    SUBCASE("float without voxel size") {
        info.sample_type = SampleType::kFloat32;
        const std::string text = denoiser::vgiText(info);
        CHECK(text.find("bitsperelement = 32\ndatatype = float\n") != std::string::npos);
        CHECK(text.find("datarange") == std::string::npos);
        CHECK(text.find("resolution") == std::string::npos);
        // No folder gets into the text.
        CHECK(text.find("some") == std::string::npos);
    }
    SUBCASE("the reader understands the text") {
        info.sample_type = SampleType::kUInt8;
        info.voxel_size[0] = 0.001;
        info.voxel_size[1] = 0.002;
        info.voxel_size[2] = 0.003;
        info.unit = "mm";
        const ScanInfo back = denoiser::parseVgi(denoiser::vgiText(info), "t");
        CHECK(back.size.width == 5);
        CHECK(back.size.height == 4);
        CHECK(back.size.depth == 3);
        CHECK(back.sample_type == SampleType::kUInt8);
        CHECK(back.voxel_size[0] == doctest::Approx(0.001));
        CHECK(back.voxel_size[1] == doctest::Approx(0.002));
        CHECK(back.voxel_size[2] == doctest::Approx(0.003));
        CHECK(back.unit == "mm");
    }
    SUBCASE("bad size") {
        info.size.depth = 0;
        CHECK_THROWS_AS(denoiser::vgiText(info), std::runtime_error);
    }
}

TEST_CASE("createScan makes both files, the voxel file full of zeros") {
    const ScanInfo info = describe("created.vol", VolumeSize{5, 4, 3}, SampleType::kUInt16);
    denoiser::createScan(info);

    const std::vector<unsigned char> bytes = fileBytes(info.volume_path);
    REQUIRE(bytes.size() == 120);
    for (const unsigned char byte : bytes) {
        CHECK(byte == 0);
    }
    CHECK(denoiser::vgiPathFor(info.volume_path) == kOut + "/created.vgi");

    // The reader opens the result from its folder name alone.
    const ScanInfo opened = denoiser::openScan(kOut + "/created.vgi");
    CHECK(opened.volume_path == info.volume_path);
    CHECK(opened.size.depth == 3);

    // Creating again replaces what was there.
    denoiser::writeRegion(info, Region{0, 0, 0, 1, 1, 1}, {77.0f});
    denoiser::createScan(info);
    CHECK(denoiser::readVolume(info)[0] == 0.0f);
}

TEST_CASE("writeVolume and readVolume agree for every sample type and byte order") {
    const VolumeSize size{5, 4, 3};
    const std::vector<float> values = counting(size);

    SUBCASE("16-bit little-endian: bytes on disk") {
        const ScanInfo info = describe("u16le.vol", size, SampleType::kUInt16);
        const WriteReport report = denoiser::writeVolume(info, values);
        CHECK(report.changed() == 0);
        const std::vector<unsigned char> bytes = fileBytes(info.volume_path);
        REQUIRE(bytes.size() == 120);
        CHECK(bytes[0] == 0xE8);  // 1000 = 0x03E8
        CHECK(bytes[1] == 0x03);
        CHECK(bytes[118] == 0x23);  // 1059 = 0x0423
        CHECK(bytes[119] == 0x04);
        CHECK(denoiser::readVolume(info) == values);
    }
    SUBCASE("16-bit big-endian: bytes on disk") {
        const ScanInfo info = describe("u16be.vol", size, SampleType::kUInt16, ByteOrder::kBigEndian);
        denoiser::writeVolume(info, values);
        const std::vector<unsigned char> bytes = fileBytes(info.volume_path);
        CHECK(bytes[0] == 0x03);
        CHECK(bytes[1] == 0xE8);
        CHECK(denoiser::readVolume(info) == values);
    }
    SUBCASE("8-bit") {
        std::vector<float> small(values.size());
        for (std::size_t i = 0; i < small.size(); ++i) {
            small[i] = static_cast<float>(100 + i);
        }
        const ScanInfo info = describe("u8.vol", size, SampleType::kUInt8);
        CHECK(denoiser::writeVolume(info, small).changed() == 0);
        CHECK(fileBytes(info.volume_path).size() == 60);
        CHECK(denoiser::readVolume(info) == small);
    }
    SUBCASE("float, both byte orders") {
        std::vector<float> reals(values.size());
        for (std::size_t i = 0; i < reals.size(); ++i) {
            reals[i] = -3.25f + 0.1f * static_cast<float>(i);
        }
        for (const ByteOrder order : {ByteOrder::kLittleEndian, ByteOrder::kBigEndian}) {
            const ScanInfo info = describe("f32.vol", size, SampleType::kFloat32, order);
            CHECK(denoiser::writeVolume(info, reals).changed() == 0);
            CHECK(fileBytes(info.volume_path).size() == 240);
            CHECK(denoiser::readVolume(info) == reals);
        }
    }
}

TEST_CASE("values that do not fit the sample type are limited and counted") {
    const float nan = std::numeric_limits<float>::quiet_NaN();
    const float inf = std::numeric_limits<float>::infinity();

    SUBCASE("16-bit") {
        const std::vector<float> in = {-0.4f, -0.6f, 0.5f,  1.49f,     2.5f,
                                       65535.4f, 65535.6f, 70000.0f, nan, inf, -inf, 12345.0f};
        const std::vector<float> expected = {0.0f, 0.0f, 1.0f, 1.0f, 3.0f,
                                             65535.0f, 65535.0f, 65535.0f, 0.0f, 65535.0f, 0.0f, 12345.0f};
        const ScanInfo info = describe("limits16.vol", VolumeSize{12, 1, 1}, SampleType::kUInt16);
        const WriteReport report = denoiser::writeVolume(info, in);
        CHECK(report.below == 2);       // -0.6 and -inf
        CHECK(report.above == 3);       // 65535.6, 70000 and inf
        CHECK(report.non_finite == 1);
        CHECK(report.changed() == 6);
        CHECK(denoiser::readVolume(info) == expected);
    }
    SUBCASE("8-bit") {
        const ScanInfo info = describe("limits8.vol", VolumeSize{4, 1, 1}, SampleType::kUInt8);
        const WriteReport report = denoiser::writeVolume(info, {-1.0f, 255.4f, 256.0f, 7.0f});
        CHECK(report.below == 1);
        CHECK(report.above == 1);
        CHECK(denoiser::readVolume(info) == std::vector<float>{0.0f, 255.0f, 255.0f, 7.0f});
    }
    SUBCASE("float keeps everything") {
        const ScanInfo info = describe("limitsf.vol", VolumeSize{4, 1, 1}, SampleType::kFloat32);
        const WriteReport report = denoiser::writeVolume(info, {-1.0e30f, nan, inf, 0.1f});
        CHECK(report.changed() == 0);
        const std::vector<float> back = denoiser::readVolume(info);
        CHECK(back[0] == -1.0e30f);
        CHECK(std::isnan(back[1]));
        CHECK(back[2] == inf);
        CHECK(back[3] == 0.1f);
    }
}

TEST_CASE("writeRegion writes the box and nothing else") {
    const VolumeSize size{5, 4, 3};
    const ScanInfo info = describe("region.vol", size, SampleType::kUInt16);
    denoiser::createScan(info);

    const Region box{1, 2, 1, 3, 2, 2};
    std::vector<float> values(12);
    for (std::size_t i = 0; i < values.size(); ++i) {
        values[i] = 500.0f + static_cast<float>(i);
    }
    CHECK(denoiser::writeRegion(info, box, values).changed() == 0);

    const std::vector<float> all = denoiser::readVolume(info);
    std::size_t next = 0;
    for (int z = 0; z < 3; ++z) {
        for (int y = 0; y < 4; ++y) {
            for (int x = 0; x < 5; ++x) {
                const bool inside = x >= 1 && x < 4 && y >= 2 && z >= 1;
                const float value = all[static_cast<std::size_t>((z * 4 + y) * 5 + x)];
                if (inside) {
                    CHECK(value == values[next++]);
                } else {
                    CHECK(value == 0.0f);
                }
            }
        }
    }
    CHECK(next == 12);
}

TEST_CASE("a volume written in blocks equals the volume written at once") {
    const VolumeSize size{7, 6, 5};
    const std::vector<float> values = counting(size);

    const ScanInfo whole = describe("whole.vol", size, SampleType::kUInt16);
    denoiser::writeVolume(whole, values);

    // Blocks of up to 4 x 4 x 2 voxels, written last block first. Blocks
    // along the full width use the slice-at-once path, the others the
    // row-by-row path.
    const ScanInfo blocks = describe("blocks.vol", size, SampleType::kUInt16);
    denoiser::createScan(blocks);
    for (int z0 = 4; z0 >= 0; z0 -= 2) {
        for (int y0 = 4; y0 >= 0; y0 -= 4) {
            for (int x0 = 4; x0 >= 0; x0 -= 4) {
                const Region box{x0, y0, z0, x0 == 4 ? 3 : 4, y0 == 4 ? 2 : 4, z0 == 4 ? 1 : 2};
                denoiser::writeRegion(blocks, box, denoiser::readRegion(whole, box));
            }
        }
    }
    CHECK(fileBytes(blocks.volume_path) == fileBytes(whole.volume_path));

    const ScanInfo slabs = describe("slabs.vol", size, SampleType::kUInt16);
    denoiser::createScan(slabs);
    for (int z0 = 0; z0 < 5; ++z0) {
        const Region slab{0, 0, z0, 7, 6, 1};
        denoiser::writeRegion(slabs, slab, denoiser::readRegion(whole, slab));
    }
    CHECK(fileBytes(slabs.volume_path) == fileBytes(whole.volume_path));
}

TEST_CASE("the writer refuses what it cannot do") {
    const VolumeSize size{5, 4, 3};
    ScanInfo info = describe("refuse.vol", size, SampleType::kUInt16);
    denoiser::createScan(info);
    const std::vector<float> twelve(12, 1.0f);

    // Regions that do not fit, and a wrong number of values.
    CHECK_THROWS_AS(denoiser::writeRegion(info, Region{3, 2, 1, 3, 2, 2}, twelve), std::runtime_error);
    CHECK_THROWS_AS(denoiser::writeRegion(info, Region{1, 2, 1, 3, 2, 3}, twelve), std::runtime_error);
    CHECK_THROWS_AS(denoiser::writeRegion(info, Region{1, 2, 1, 3, 2, 1}, twelve), std::runtime_error);
    CHECK_THROWS_AS(denoiser::writeVolume(info, twelve), std::runtime_error);

    // A file that createScan() did not make for this description.
    ScanInfo other = info;
    other.size.depth = 4;
    CHECK_THROWS_AS(denoiser::writeRegion(other, Region{0, 0, 0, 3, 2, 2}, twelve), std::runtime_error);
    ScanInfo missing = describe("never_created.vol", size, SampleType::kUInt16);
    CHECK_THROWS_AS(denoiser::writeRegion(missing, Region{0, 0, 0, 3, 2, 2}, twelve), std::runtime_error);

    // Descriptions the writer does not accept.
    ScanInfo header = info;
    header.header_bytes = 16;
    CHECK_THROWS_AS(denoiser::createScan(header), std::runtime_error);
    ScanInfo no_path = info;
    no_path.volume_path.clear();
    CHECK_THROWS_AS(denoiser::createScan(no_path), std::runtime_error);
    ScanInfo vgi_name = info;
    vgi_name.volume_path = kOut + "/wrong.vgi";
    CHECK_THROWS_AS(denoiser::createScan(vgi_name), std::runtime_error);
    ScanInfo no_folder = info;
    no_folder.volume_path = kOut + "/no_such_folder/x.vol";
    CHECK_THROWS_AS(denoiser::createScan(no_folder), std::runtime_error);

    // The refused calls left the file as it was.
    CHECK(denoiser::readVolume(info) == std::vector<float>(60, 0.0f));
}

TEST_CASE("a tracked test scan is copied byte for byte") {
    for (const std::string& scan : {std::string("Sample1_18p15_5mu_c64"), std::string("Sample2_5mu_c64")}) {
        CAPTURE(scan);
        const ScanInfo source = denoiser::openScan("tests/data/" + scan);

        std::filesystem::create_directories(kOut + "/" + scan);
        ScanInfo copy = source;
        copy.volume_path = kOut + "/" + scan + "/" + scan + ".vol";
        CHECK(denoiser::writeVolume(copy, denoiser::readVolume(source)).changed() == 0);

        const std::string base = "tests/data/" + scan + "/" + scan;
        CHECK(fileBytes(copy.volume_path) == fileBytes(base + ".vol"));
        CHECK(fileBytes(denoiser::vgiPathFor(copy.volume_path)) == fileBytes(base + ".vgi"));
    }
}
