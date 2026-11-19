// Tests of the fast Non-Local Means filter (include/nlm_fast.hpp) against
// the brute-force filter (include/nlm_bruteforce.hpp), which is the
// reference: both must give the same result to within rounding for every
// combination of settings.

#include "doctest/doctest.h"

#include "nlm_bruteforce.hpp"
#include "nlm_fast.hpp"
#include "scan_reader.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

using denoiser::NlmFastOptions;
using denoiser::NlmKernel;
using denoiser::NlmParams;
using denoiser::VolumeSize;

namespace {

/// Repeatable pseudo-random values in 0...100 with some structure: a
/// smooth part along z plus noise, so that patch distances vary.
std::vector<float> testVolume(const VolumeSize& size, std::uint32_t seed) {
    std::vector<float> volume;
    std::uint32_t state = seed;
    for (int z = 0; z < size.depth; ++z) {
        for (int y = 0; y < size.height; ++y) {
            for (int x = 0; x < size.width; ++x) {
                state = state * 1664525u + 1013904223u;
                const float noise = static_cast<float>(state >> 16) * (30.0f / 65535.0f);
                const float structure = ((x / 3 + y / 2) % 2 == 0) ? 20.0f : 60.0f;
                volume.push_back(structure + noise);
            }
        }
    }
    return volume;
}

NlmParams makeParams(double sigma, double h, int p0, int p1, int p2, int s0, int s1, int s2) {
    NlmParams params;
    params.sigma = sigma;
    params.h = h;
    params.patch_radius[0] = p0;
    params.patch_radius[1] = p1;
    params.patch_radius[2] = p2;
    params.search_radius[0] = s0;
    params.search_radius[1] = s1;
    params.search_radius[2] = s2;
    return params;
}

/// Largest absolute difference between two volumes.
double largestDifference(const std::vector<float>& a, const std::vector<float>& b) {
    REQUIRE(a.size() == b.size());
    double largest = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        largest = std::max(largest, std::fabs(static_cast<double>(a[i]) - b[i]));
    }
    return largest;
}

}  // namespace

TEST_CASE("nlmFast equals nlmBruteForce for many settings") {
    // Blocks of 4 voxels, so that every volume is cut into several blocks
    // and the block edges are tested.
    const VolumeSize size{11, 9, 13};
    const std::vector<float> volume = testVolume(size, 3);
    NlmFastOptions options;
    options.block_size = 4;
    options.single_precision = false;
    NlmFastOptions single = options;
    single.single_precision = true;

    struct Case {
        const char* name;
        NlmParams params;
    };
    std::vector<Case> cases = {
        {"cube patch", makeParams(5, 4, 1, 1, 1, 2, 2, 2)},
        {"column patch", makeParams(5, 4, 0, 0, 3, 1, 1, 4)},
        {"patch along x only", makeParams(5, 4, 2, 0, 0, 3, 1, 1)},
        {"patch along y only", makeParams(5, 4, 0, 2, 0, 1, 3, 1)},
        {"patch of one voxel", makeParams(5, 4, 0, 0, 0, 2, 2, 2)},
        {"uneven patch", makeParams(5, 4, 1, 2, 3, 1, 2, 3)},
        {"search window larger than the volume", makeParams(5, 4, 1, 1, 1, 6, 5, 7)},
        {"search along z only", makeParams(5, 4, 1, 1, 1, 0, 0, 5)},
        {"search window of one voxel", makeParams(5, 4, 1, 1, 1, 0, 0, 0)},
        {"patch larger than the block", makeParams(5, 4, 5, 5, 6, 1, 1, 1)},
    };
    NlmParams normalised = makeParams(5, 4, 1, 1, 2, 2, 2, 2);
    normalised.normalise_distance = true;
    cases.push_back({"mean distance", normalised});
    NlmParams exponential = makeParams(5, 20, 1, 1, 2, 2, 2, 2);
    exponential.kernel = NlmKernel::kExponential;
    cases.push_back({"exponential weight", exponential});
    NlmParams ipol = exponential;
    ipol.normalise_distance = true;
    ipol.h = 4;
    cases.push_back({"IPOL variant", ipol});

    for (const Case& c : cases) {
        CAPTURE(std::string(c.name));
        const std::vector<float> reference = denoiser::nlmBruteForce(volume, size, c.params);
        const std::vector<float> fast = denoiser::nlmFast(volume, size, c.params, options);
        CHECK(largestDifference(fast, reference) < 1e-4);
        CHECK(largestDifference(denoiser::nlmFast(volume, size, c.params, single), reference) < 0.01);
    }
}

TEST_CASE("nlmFast with the original table on a tracked test scan") {
    // The 64^3 scan on the 0...255 scale, sigma 17 (patch 2, search 10),
    // compared on a 24^3 part to keep the brute-force run short.
    const denoiser::ScanInfo info = denoiser::openScan("tests/data/Sample2_5mu_c64");
    const VolumeSize size{24, 24, 24};
    std::vector<float> volume = denoiser::readRegion(info, denoiser::Region{20, 20, 20, 24, 24, 24});
    for (float& v : volume) {
        v *= 255.0f / 65535.0f;
    }
    const NlmParams params = denoiser::originalParams(17.0);
    const std::vector<float> reference = denoiser::nlmBruteForce(volume, size, params);
    NlmFastOptions in_double;
    in_double.single_precision = false;
    CHECK(largestDifference(denoiser::nlmFast(volume, size, params, in_double), reference) < 1e-4);

    // The default is single precision.
    CHECK(NlmFastOptions{}.single_precision);
    CHECK(largestDifference(denoiser::nlmFast(volume, size, params), reference) < 0.01);
}

TEST_CASE("nlmFast: block size and threads change nothing but rounding") {
    const VolumeSize size{20, 17, 23};
    const std::vector<float> volume = testVolume(size, 9);
    const NlmParams params = makeParams(6, 5, 0, 0, 4, 2, 2, 6);

    NlmFastOptions small;
    small.block_size = 5;
    small.single_precision = false;
    NlmFastOptions one_block = small;
    one_block.block_size = 64;
    const std::vector<float> a = denoiser::nlmFast(volume, size, params, small);
    const std::vector<float> b = denoiser::nlmFast(volume, size, params, one_block);
    CHECK(largestDifference(a, b) < 1e-4);

    // The same call twice gives the same bits, whatever the thread timing.
    CHECK(denoiser::nlmFast(volume, size, params, small) == a);
}

TEST_CASE("nlmFast keeps a constant volume and the range of the input") {
    const VolumeSize size{6, 5, 7};
    const NlmParams params = makeParams(5, 4, 1, 1, 2, 2, 2, 3);
    const std::vector<float> constant(6 * 5 * 7, 42.5f);
    for (float value : denoiser::nlmFast(constant, size, params)) {
        CHECK(value == 42.5f);
    }
    const std::vector<float> volume = testVolume(size, 4);
    const auto [lowest, highest] = std::minmax_element(volume.begin(), volume.end());
    for (float value : denoiser::nlmFast(volume, size, params)) {
        CHECK(value >= *lowest);
        CHECK(value <= *highest);
    }
}

TEST_CASE("nlmFast rejects bad arguments") {
    const VolumeSize size{2, 2, 2};
    const std::vector<float> volume(8, 1.0f);
    const NlmParams params = makeParams(5, 4, 1, 1, 1, 1, 1, 1);
    CHECK_THROWS_AS(denoiser::nlmFast(volume, VolumeSize{2, 2, 3}, params), std::invalid_argument);
    NlmParams bad = params;
    bad.h = 0.0;
    CHECK_THROWS_AS(denoiser::nlmFast(volume, size, bad), std::invalid_argument);
    bad = params;
    bad.search_radius[2] = -1;
    CHECK_THROWS_AS(denoiser::nlmFast(volume, size, bad), std::invalid_argument);
    NlmFastOptions options;
    options.block_size = 0;
    CHECK_THROWS_AS(denoiser::nlmFast(volume, size, params, options), std::invalid_argument);
}
