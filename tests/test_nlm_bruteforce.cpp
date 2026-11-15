// Tests of the brute-force Non-Local Means filter (include/nlm_bruteforce.hpp).

#include "doctest/doctest.h"

#include "nlm_bruteforce.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

using denoiser::NlmParams;
using denoiser::VolumeSize;
using denoiser::nlmBruteForce;
using denoiser::originalParams;

namespace {

/// Returns a volume of repeatable pseudo-random values in 0...100.
/// A small integer generator is used so that every platform gets the same
/// numbers.
std::vector<float> randomVolume(const VolumeSize& size, std::uint32_t seed) {
    std::vector<float> volume(static_cast<std::size_t>(size.width) * size.height * size.depth);
    std::uint32_t state = seed;
    for (float& value : volume) {
        state = state * 1664525u + 1013904223u;
        value = static_cast<float>(state >> 16) * (100.0f / 65535.0f);
    }
    return volume;
}

std::size_t indexOf(const VolumeSize& size, int x, int y, int z) {
    return (static_cast<std::size_t>(z) * size.height + y) * size.width + x;
}

NlmParams smallParams() {
    NlmParams params;
    params.sigma = 5.0;
    params.h = 4.0;
    for (int axis = 0; axis < 3; ++axis) {
        params.patch_radius[axis] = 1;
        params.search_radius[axis] = 2;
    }
    return params;
}

}  // namespace

TEST_CASE("originalParams follows the table of the original program") {
    struct Row { double sigma; int patch; int search; double h; };
    const Row rows[] = {
        {10.0, 1, 10, 4.0}, {15.0, 1, 10, 6.0}, {17.0, 2, 10, 6.8}, {30.0, 2, 10, 12.0},
        {40.0, 3, 17, 14.0}, {60.0, 4, 17, 21.0}, {100.0, 5, 17, 30.0},
    };
    for (const Row& row : rows) {
        CAPTURE(row.sigma);
        const NlmParams params = originalParams(row.sigma);
        CHECK(params.sigma == doctest::Approx(row.sigma));
        CHECK(params.h == doctest::Approx(row.h).epsilon(1e-6));
        for (int axis = 0; axis < 3; ++axis) {
            CHECK(params.patch_radius[axis] == row.patch);
            CHECK(params.search_radius[axis] == row.search);
        }
    }
    CHECK_THROWS_AS(originalParams(0.0), std::invalid_argument);
    CHECK_THROWS_AS(originalParams(-1.0), std::invalid_argument);
}

TEST_CASE("nlmBruteForce leaves a constant volume unchanged") {
    const VolumeSize size{5, 4, 3};
    const std::vector<float> volume(5 * 4 * 3, 42.5f);
    const std::vector<float> result = nlmBruteForce(volume, size, smallParams());
    REQUIRE(result.size() == volume.size());
    for (float value : result) {
        CHECK(value == 42.5f);
    }
}

TEST_CASE("nlmBruteForce matches a case worked out by hand") {
    // Three voxels in a row: 0, 0, 9. Patch of one voxel, search radius 1
    // along x, sigma = 1, h = 1, so the threshold is 2 and h^2 is 1.
    //
    // Voxel 0: both neighbours (the left one is itself by edge replication)
    //          have value 0, distance 0, weight 1.            -> 0
    // Voxel 1: left 0 with weight 1; right 9 with distance 81, excess 79,
    //          weight 1/80; itself 0 with the largest weight 1.
    //          (9/80) / (2 + 1/80)                             -> 9/161
    // Voxel 2: left 0 with weight 1/80; right is itself, 9 with weight 1;
    //          itself 9 with weight 1.  18 / (2 + 1/80)        -> 1440/161
    const VolumeSize size{3, 1, 1};
    NlmParams params;
    params.sigma = 1.0;
    params.h = 1.0;
    for (int axis = 0; axis < 3; ++axis) {
        params.patch_radius[axis] = 0;
        params.search_radius[axis] = 0;
    }
    params.search_radius[0] = 1;

    const std::vector<float> result = nlmBruteForce({0.0f, 0.0f, 9.0f}, size, params);
    CHECK(result[0] == doctest::Approx(0.0));
    CHECK(result[1] == doctest::Approx(9.0 / 161.0));
    CHECK(result[2] == doctest::Approx(1440.0 / 161.0));
}

TEST_CASE("nlmBruteForce with a search window of one voxel returns the input") {
    const VolumeSize size{4, 3, 2};
    const std::vector<float> volume = randomVolume(size, 7);
    NlmParams params = smallParams();
    for (int axis = 0; axis < 3; ++axis) {
        params.search_radius[axis] = 0;
    }
    CHECK(nlmBruteForce(volume, size, params) == volume);
}

TEST_CASE("nlmBruteForce results stay within the range of the input") {
    const VolumeSize size{7, 6, 5};
    const std::vector<float> volume = randomVolume(size, 1);
    const auto [lowest, highest] = std::minmax_element(volume.begin(), volume.end());
    for (float value : nlmBruteForce(volume, size, smallParams())) {
        CHECK(value >= *lowest);
        CHECK(value <= *highest);
    }
}

TEST_CASE("nlmBruteForce commutes with mirroring the volume along x") {
    const VolumeSize size{7, 5, 4};
    const std::vector<float> volume = randomVolume(size, 2);
    std::vector<float> mirrored(volume.size());
    for (int z = 0; z < size.depth; ++z)
        for (int y = 0; y < size.height; ++y)
            for (int x = 0; x < size.width; ++x)
                mirrored[indexOf(size, size.width - 1 - x, y, z)] = volume[indexOf(size, x, y, z)];

    const std::vector<float> result = nlmBruteForce(volume, size, smallParams());
    const std::vector<float> result_mirrored = nlmBruteForce(mirrored, size, smallParams());
    for (int z = 0; z < size.depth; ++z)
        for (int y = 0; y < size.height; ++y)
            for (int x = 0; x < size.width; ++x)
                CHECK(result_mirrored[indexOf(size, size.width - 1 - x, y, z)]
                      == doctest::Approx(result[indexOf(size, x, y, z)]).epsilon(1e-6));
}

TEST_CASE("nlmBruteForce applies per-axis radii to the right axes") {
    // Swapping the x and z axes of the volume and of the radii must give the
    // swapped result.
    const VolumeSize size{6, 4, 3};
    const VolumeSize swapped_size{3, 4, 6};
    const std::vector<float> volume = randomVolume(size, 3);
    std::vector<float> swapped(volume.size());
    for (int z = 0; z < size.depth; ++z)
        for (int y = 0; y < size.height; ++y)
            for (int x = 0; x < size.width; ++x)
                swapped[indexOf(swapped_size, z, y, x)] = volume[indexOf(size, x, y, z)];

    NlmParams params = smallParams();
    params.patch_radius[0] = 2; params.patch_radius[1] = 1; params.patch_radius[2] = 0;
    params.search_radius[0] = 3; params.search_radius[1] = 1; params.search_radius[2] = 2;
    NlmParams swapped_params = params;
    std::swap(swapped_params.patch_radius[0], swapped_params.patch_radius[2]);
    std::swap(swapped_params.search_radius[0], swapped_params.search_radius[2]);

    const std::vector<float> result = nlmBruteForce(volume, size, params);
    const std::vector<float> result_swapped = nlmBruteForce(swapped, swapped_size, swapped_params);
    for (int z = 0; z < size.depth; ++z)
        for (int y = 0; y < size.height; ++y)
            for (int x = 0; x < size.width; ++x)
                CHECK(result_swapped[indexOf(swapped_size, z, y, x)]
                      == doctest::Approx(result[indexOf(size, x, y, z)]).epsilon(1e-6));
}

TEST_CASE("nlmBruteForce rejects bad arguments") {
    const VolumeSize size{2, 2, 2};
    const std::vector<float> volume(8, 1.0f);
    NlmParams params = smallParams();

    CHECK_THROWS_AS(nlmBruteForce(volume, VolumeSize{2, 2, 3}, params), std::invalid_argument);
    CHECK_THROWS_AS(nlmBruteForce(volume, VolumeSize{0, 2, 2}, params), std::invalid_argument);

    params.sigma = 0.0;
    CHECK_THROWS_AS(nlmBruteForce(volume, size, params), std::invalid_argument);
    params = smallParams();
    params.h = -1.0;
    CHECK_THROWS_AS(nlmBruteForce(volume, size, params), std::invalid_argument);
    params = smallParams();
    params.patch_radius[1] = -1;
    CHECK_THROWS_AS(nlmBruteForce(volume, size, params), std::invalid_argument);
}
