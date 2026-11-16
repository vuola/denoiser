// Tests of the pre-processing (include/preprocess.hpp) and of the brightness
// profile (include/radial_profile.hpp).
//
// Files are written into build/test_tmp/preprocess/. The tests run from the
// repository root.

#include "doctest/doctest.h"

#include "preprocess.hpp"
#include "radial_profile.hpp"
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

using denoiser::IntensityMap;
using denoiser::ProfileSettings;
using denoiser::RadialRing;
using denoiser::Region;
using denoiser::SampleType;
using denoiser::ScanInfo;
using denoiser::VolumeSize;

namespace {

const std::string kOut = "build/test_tmp/preprocess";

/// Returns a description of a 16-bit scan to be written below kOut.
ScanInfo describe(const std::string& name, const VolumeSize& size,
                  SampleType type = SampleType::kUInt16) {
    std::filesystem::create_directories(kOut);
    ScanInfo info;
    info.volume_path = kOut + "/" + name;
    info.size = size;
    info.sample_type = type;
    return info;
}

/// Returns all bytes of a file.
std::vector<unsigned char> fileBytes(const std::string& path) {
    std::ifstream file(path, std::ios::binary);
    REQUIRE(file.good());
    return std::vector<unsigned char>((std::istreambuf_iterator<char>(file)),
                                      std::istreambuf_iterator<char>());
}

}  // namespace

TEST_CASE("IntensityMap: the default mapping is the one of the old workflow") {
    const IntensityMap map = denoiser::defaultMap(SampleType::kUInt16);
    CHECK(map.low == 0.0);
    CHECK(map.high == 65535.0);
    CHECK(map.toWork(0.0f) == 0.0f);
    CHECK(map.toWork(65535.0f) == 255.0f);
    // count * 255 / 65535, computed in double and rounded to float once.
    CHECK(map.toWork(31139.0f) == static_cast<float>(31139.0 * 255.0 / 65535.0));

    CHECK(denoiser::defaultMap(SampleType::kUInt8).high == 255.0);
    // Float scans are taken as they are.
    const IntensityMap same = denoiser::defaultMap(SampleType::kFloat32);
    CHECK(same.toWork(123.25f) == 123.25f);
    CHECK(same.toStored(123.25f) == 123.25f);
}

TEST_CASE("IntensityMap: every 16-bit count survives the way there and back") {
    for (const IntensityMap map : {IntensityMap{0.0, 65535.0}, IntensityMap{17900.0, 48000.0},
                                   IntensityMap{-1000.0, 70000.0}}) {
        CAPTURE(map.low);
        int wrong = 0;
        for (int count = 0; count <= 65535; ++count) {
            const float back = map.toStored(map.toWork(static_cast<float>(count)));
            if (std::lround(back) != count) {
                ++wrong;
            }
        }
        CHECK(wrong == 0);
    }
}

TEST_CASE("IntensityMap: linear, and nothing is clipped") {
    const IntensityMap map{10000.0, 20000.0};
    CHECK(map.toWork(10000.0f) == 0.0f);
    CHECK(map.toWork(20000.0f) == 255.0f);
    CHECK(map.toWork(15000.0f) == doctest::Approx(127.5));
    // Outside the mapped range the line simply continues.
    CHECK(map.toWork(0.0f) == doctest::Approx(-255.0));
    CHECK(map.toWork(30000.0f) == doctest::Approx(510.0));
    CHECK(map.toStored(-255.0f) == doctest::Approx(0.0));
    CHECK(map.toStored(510.0f) == doctest::Approx(30000.0));

    std::vector<float> values = {0.0f, 10000.0f, 15000.0f, 30000.0f};
    denoiser::mapToWork(values, map);
    CHECK(values[0] == doctest::Approx(-255.0));
    CHECK(values[2] == doctest::Approx(127.5));
    denoiser::mapToStored(values, map);
    CHECK(values[0] == doctest::Approx(0.0));
    CHECK(values[3] == doctest::Approx(30000.0));
}

TEST_CASE("IntensityMap: bad ranges are refused") {
    std::vector<float> values = {1.0f};
    CHECK_THROWS_AS(denoiser::mapToWork(values, IntensityMap{5.0, 5.0}), std::invalid_argument);
    CHECK_THROWS_AS(denoiser::mapToWork(values, IntensityMap{6.0, 5.0}), std::invalid_argument);
    CHECK_THROWS_AS(denoiser::mapToStored(values, IntensityMap{0.0, std::nan("")}),
                    std::invalid_argument);
    CHECK_THROWS_AS(denoiser::checkMap(IntensityMap{-std::numeric_limits<double>::infinity(), 1.0}),
                    std::invalid_argument);
    CHECK(values[0] == 1.0f);
}

TEST_CASE("centredCube") {
    // The 250-voxel cube of the old workflow in a 1142 x 1142 x 1000 scan
    // starts at (446, 446, 375). The old macro started z at 374 because
    // ImageJ counts slices from 1 (README section 2.3).
    const Region old = denoiser::centredCube(VolumeSize{1142, 1142, 1000}, 250);
    CHECK(old.x0 == 446);
    CHECK(old.y0 == 446);
    CHECK(old.z0 == 375);
    CHECK(old.width == 250);
    CHECK(old.height == 250);
    CHECK(old.depth == 250);

    // Odd space left over: the cube starts half a voxel towards the origin.
    const Region odd = denoiser::centredCube(VolumeSize{10, 11, 12}, 5);
    CHECK(odd.x0 == 2);
    CHECK(odd.y0 == 3);
    CHECK(odd.z0 == 3);

    // The whole volume, and cubes that do not fit.
    CHECK(denoiser::centredCube(VolumeSize{7, 7, 7}, 7).x0 == 0);
    CHECK_THROWS_AS(denoiser::centredCube(VolumeSize{7, 7, 6}, 7), std::invalid_argument);
    CHECK_THROWS_AS(denoiser::centredCube(VolumeSize{7, 7, 7}, 0), std::invalid_argument);
}

TEST_CASE("replaceNonFinite") {
    const float nan = std::numeric_limits<float>::quiet_NaN();
    const float inf = std::numeric_limits<float>::infinity();
    std::vector<float> values = {1.0f, nan, 3.0f, inf, -inf, 6.0f};
    CHECK(denoiser::replaceNonFinite(values, -1.0f) == 3);
    CHECK(values == std::vector<float>{1.0f, -1.0f, 3.0f, -1.0f, -1.0f, 6.0f});
    CHECK(denoiser::replaceNonFinite(values, 0.0f) == 0);
}

TEST_CASE("cropScan cuts the right voxels and keeps the description") {
    const ScanInfo source = denoiser::openScan("tests/data/Sample2_5mu_c64");
    std::filesystem::create_directories(kOut);
    const Region region{3, 7, 5, 20, 10, 30};

    const ScanInfo crop = denoiser::cropScan(source, region, kOut + "/crop.vol");
    CHECK(crop.size.width == 20);
    CHECK(crop.size.height == 10);
    CHECK(crop.size.depth == 30);
    CHECK(crop.sample_type == SampleType::kUInt16);

    // The files on disk say the same, including the voxel size.
    const ScanInfo opened = denoiser::openScan(kOut + "/crop.vgi");
    CHECK(opened.size.depth == 30);
    CHECK(opened.voxel_size[0] == doctest::Approx(0.005));
    CHECK(opened.unit == "mm");

    CHECK(denoiser::readVolume(opened) == denoiser::readRegion(source, region));
    // Voxel (3, 7, 5) of the source, a value known from step 1.1.
    CHECK(denoiser::readVolume(opened).front() == 20215.0f);

    // The whole scan as its own "crop" is a byte-for-byte copy.
    denoiser::cropScan(source, denoiser::wholeVolume(source.size), kOut + "/all.vol");
    CHECK(fileBytes(kOut + "/all.vol") ==
          fileBytes("tests/data/Sample2_5mu_c64/Sample2_5mu_c64.vol"));

    // A region that does not fit, and writing over the source.
    CHECK_THROWS_AS(denoiser::cropScan(source, Region{60, 0, 0, 5, 5, 5}, kOut + "/bad.vol"),
                    std::runtime_error);
    CHECK_THROWS_AS(denoiser::cropScan(opened, Region{0, 0, 0, 5, 5, 5}, kOut + "/crop.vol"),
                    std::runtime_error);
    CHECK(denoiser::readVolume(opened).front() == 20215.0f);
}

TEST_CASE("radialProfile on a volume made of rings") {
    // 41 x 41 x 3 voxels around the axis (20, 20). The value depends on the
    // distance r from the axis only: 1000 * (whole part of r / 4).
    const VolumeSize size{41, 41, 3};
    std::vector<float> values;
    for (int z = 0; z < 3; ++z) {
        for (int y = 0; y < 41; ++y) {
            for (int x = 0; x < 41; ++x) {
                const double r = std::hypot(x - 20.0, y - 20.0);
                values.push_back(1000.0f * static_cast<float>(static_cast<int>(r / 4.0)));
            }
        }
    }
    const ScanInfo info = describe("rings.vol", size);
    denoiser::writeVolume(info, values);

    ProfileSettings settings = denoiser::defaultProfileSettings(info);
    CHECK(settings.centre_x == 20.0);
    CHECK(settings.centre_y == 20.0);
    CHECK(settings.ring_width == 8.0);
    CHECK(settings.value_high == 65535.0);
    settings.ring_width = 4.0;

    const std::vector<RadialRing> rings = denoiser::radialProfile(info, settings);
    // The farthest corner is at r = 28.3, in ring 7.
    REQUIRE(rings.size() == 8);
    std::int64_t total = 0;
    for (std::size_t k = 0; k < rings.size(); ++k) {
        CAPTURE(k);
        CHECK(rings[k].r_inner == 4.0 * static_cast<double>(k));
        CHECK(rings[k].r_outer == 4.0 * static_cast<double>(k + 1));
        CHECK(rings[k].voxels > 0);
        CHECK(rings[k].mean == doctest::Approx(1000.0 * static_cast<double>(k)));
        // All voxels of a ring have one value; the percentiles hit it to
        // within one histogram class (64 counts).
        CHECK(std::fabs(rings[k].p10 - 1000.0 * static_cast<double>(k)) <= 64.0);
        CHECK(std::fabs(rings[k].p50 - 1000.0 * static_cast<double>(k)) <= 64.0);
        CHECK(std::fabs(rings[k].p90 - 1000.0 * static_cast<double>(k)) <= 64.0);
        CHECK(rings[k].at_low == (k == 0 ? rings[k].voxels : 0));
        CHECK(rings[k].at_high == 0);
        total += rings[k].voxels;
    }
    CHECK(total == 41 * 41 * 3);
    // Ring 0 holds the voxels with r < 4: 45 per slice.
    CHECK(rings[0].voxels == 45 * 3);

    SUBCASE("a part of the slices") {
        settings.z0 = 1;
        settings.depth = 1;
        const std::vector<RadialRing> one = denoiser::radialProfile(info, settings);
        CHECK(one[0].voxels == 45);
    }
    SUBCASE("an axis off the middle changes the rings") {
        settings.centre_x = 0.0;
        settings.centre_y = 0.0;
        const std::vector<RadialRing> corner = denoiser::radialProfile(info, settings);
        // The farthest corner is now at r = 56.6, in ring 14.
        CHECK(corner.size() == 15);
        // Quarter disc of radius 4 around the corner voxel: 15 voxels per slice.
        CHECK(corner[0].voxels == 15 * 3);
    }
    SUBCASE("settings out of range") {
        ProfileSettings bad = settings;
        bad.ring_width = 0.0;
        CHECK_THROWS_AS(denoiser::radialProfile(info, bad), std::invalid_argument);
        bad = settings;
        bad.z0 = 3;
        CHECK_THROWS_AS(denoiser::radialProfile(info, bad), std::invalid_argument);
        bad = settings;
        bad.z0 = 1;
        bad.depth = 3;
        CHECK_THROWS_AS(denoiser::radialProfile(info, bad), std::invalid_argument);
        bad = settings;
        bad.value_high = bad.value_low;
        CHECK_THROWS_AS(denoiser::radialProfile(info, bad), std::invalid_argument);
    }
}

TEST_CASE("radialProfile: percentiles and voxels at the limits") {
    // One ring that holds everything: 100 x 1 x 1 voxels, axis far away is
    // not needed; a ring width of 1000 puts all voxels into ring 0.
    // Values: 0, 640, 1280, ..., 63360, so that every value sits at the
    // lower edge of its own histogram class (64 counts wide), except that
    // the last two voxels are set to 65535.
    std::vector<float> values;
    for (int i = 0; i < 100; ++i) {
        values.push_back(i < 98 ? 640.0f * static_cast<float>(i) : 65535.0f);
    }
    const ScanInfo info = describe("line.vol", VolumeSize{100, 1, 1});
    denoiser::writeVolume(info, values);

    ProfileSettings settings = denoiser::defaultProfileSettings(info);
    settings.ring_width = 1000.0;
    const std::vector<RadialRing> rings = denoiser::radialProfile(info, settings);
    REQUIRE(rings.size() == 1);
    CHECK(rings[0].voxels == 100);
    CHECK(rings[0].at_low == 1);
    CHECK(rings[0].at_high == 2);
    // The 10th value from below is 5760, the 50th 31360, the 90th 56960;
    // the percentile lies within one class above each.
    CHECK(rings[0].p10 >= 5760.0);
    CHECK(rings[0].p10 <= 5760.0 + 64.0);
    CHECK(rings[0].p50 >= 31360.0);
    CHECK(rings[0].p50 <= 31360.0 + 64.0);
    CHECK(rings[0].p90 >= 56960.0);
    CHECK(rings[0].p90 <= 56960.0 + 64.0);
}

TEST_CASE("scanStatistics counts the voxels at the limits of the range") {
    const std::vector<float> values = {0.0f, 0.0f, 1.0f, 65535.0f, 65534.0f, 300.0f};
    const ScanInfo info = describe("limits.vol", VolumeSize{6, 1, 1});
    denoiser::writeVolume(info, values);
    const denoiser::ScanStatistics stats = denoiser::scanStatistics(info);
    CHECK(stats.at_low_limit == 2);
    CHECK(stats.at_high_limit == 1);

    // The tracked test scans have none.
    const denoiser::ScanStatistics tracked =
        denoiser::scanStatistics(denoiser::openScan("tests/data/Sample1_18p15_5mu_c64"));
    CHECK(tracked.at_low_limit == 0);
    CHECK(tracked.at_high_limit == 0);

    // For float scans the limits do not apply.
    ScanInfo real = describe("limits_f.vol", VolumeSize{6, 1, 1}, SampleType::kFloat32);
    denoiser::writeVolume(real, values);
    const denoiser::ScanStatistics real_stats = denoiser::scanStatistics(real);
    CHECK(real_stats.at_low_limit == 0);
    CHECK(real_stats.at_high_limit == 0);
}
