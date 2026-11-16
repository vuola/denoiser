// Tests of the noise estimate (include/noise_estimate.hpp).
//
// Files are written into build/test_tmp/noise/. The tests run from the
// repository root.

#include "doctest/doctest.h"

#include "noise_estimate.hpp"
#include "scan_reader.hpp"
#include "scan_writer.hpp"

#include <cmath>
#include <cstddef>
#include <filesystem>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

using denoiser::NoiseEstimate;
using denoiser::NoiseSummary;
using denoiser::Region;
using denoiser::SampleType;
using denoiser::ScanInfo;
using denoiser::SliceNoise;
using denoiser::VolumeSize;

namespace {

const std::string kOut = "build/test_tmp/noise";

/// Returns a description of a float scan to be written below kOut.
ScanInfo describe(const std::string& name, const VolumeSize& size) {
    std::filesystem::create_directories(kOut);
    ScanInfo info;
    info.volume_path = kOut + "/" + name;
    info.size = size;
    info.sample_type = SampleType::kFloat32;
    return info;
}

/// Noise with mean 0 and standard deviation 1, close to normal: the sum of
/// twelve uniform numbers minus 6. std::mt19937 gives the same numbers on
/// every platform, which std::normal_distribution does not.
class Noise {
public:
    explicit Noise(unsigned seed) : generator_(seed) {}
    double next() {
        double sum = 0.0;
        for (int i = 0; i < 12; ++i) {
            sum += static_cast<double>(generator_()) / 4294967296.0;
        }
        return sum - 6.0;
    }

private:
    std::mt19937 generator_;
};

/// A slice-independent pattern standing in for material: stripes 4 voxels
/// wide in x with values 0 and `contrast`.
float pattern(int x, float contrast) {
    return (x / 4) % 2 == 0 ? 0.0f : contrast;
}

}  // namespace

TEST_CASE("sampleSlices spreads the slices from first to last") {
    const std::vector<int> twenty = denoiser::sampleSlices(1000, 20);
    REQUIRE(twenty.size() == 20);
    CHECK(twenty.front() == 0);
    CHECK(twenty[1] == 52);   // 999 / 19 = 52.58, rounded down
    CHECK(twenty[10] == 525);  // 9990 / 19 = 525.79
    CHECK(twenty.back() == 999);

    // More slices asked for than there are: every slice once.
    CHECK(denoiser::sampleSlices(5, 20) == std::vector<int>{0, 1, 2, 3, 4});
    CHECK(denoiser::sampleSlices(10, 1) == std::vector<int>{0});
    CHECK(denoiser::sampleSlices(1, 1) == std::vector<int>{0});

    CHECK_THROWS_AS(denoiser::sampleSlices(0, 5), std::invalid_argument);
    CHECK_THROWS_AS(denoiser::sampleSlices(5, 0), std::invalid_argument);
}

TEST_CASE("summarise: median, mean and quartiles as numpy gives them") {
    const NoiseSummary s = denoiser::summarise({4.0, 1.0, 3.0, 2.0});
    CHECK(s.count == 4);
    CHECK(s.median == doctest::Approx(2.5));
    CHECK(s.mean == doctest::Approx(2.5));
    CHECK(s.q1 == doctest::Approx(1.75));
    CHECK(s.q3 == doctest::Approx(3.25));

    // NaN stands for "no value from this slice" and is left out.
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const NoiseSummary with_nan = denoiser::summarise({5.0, nan, 1.0, 3.0});
    CHECK(with_nan.count == 3);
    CHECK(with_nan.median == doctest::Approx(3.0));
    CHECK(with_nan.q1 == doctest::Approx(2.0));

    const NoiseSummary none = denoiser::summarise({nan});
    CHECK(none.count == 0);
    CHECK(std::isnan(none.median));
    CHECK(std::isnan(none.mean));
}

TEST_CASE("sliceNoise on four voxels, computed by hand") {
    // Slice:  0 2      next slice:  1 3
    //         4 6                   5 9
    const std::vector<float> slice = {0.0f, 2.0f, 4.0f, 6.0f};
    const std::vector<float> next = {1.0f, 3.0f, 5.0f, 9.0f};
    const SliceNoise n = denoiser::sliceNoise(slice.data(), next.data(), 2, 2);

    CHECK(n.std_dev == doctest::Approx(std::sqrt(5.0)));
    // Median 3, distances from it 3 1 1 3, their median 2.
    CHECK(n.mad == doctest::Approx(1.4826 * 2.0));
    // Differences 2 2 (along x) and 4 4 (along y): std 1.
    CHECK(n.diff == doctest::Approx(1.0 / std::sqrt(2.0)));
    // 3 x 3 means with repeated edge voxels: 2, 8/3, 10/3, 4. Residuals
    // -2, -2/3, 2/3, 2; median 0; distances 2, 2/3, 2/3, 2; median 4/3.
    CHECK(n.resid == doctest::Approx(1.4826 * 4.0 / 3.0));
    // Differences to the next slice 1 1 1 3: std sqrt(0.75); MAD 0.
    CHECK(n.zdiff == doctest::Approx(std::sqrt(0.75) / std::sqrt(2.0)));
    CHECK(n.zdiff_robust == doctest::Approx(0.0));

    // Without a next slice there are no z estimates.
    const SliceNoise alone = denoiser::sliceNoise(slice.data(), nullptr, 2, 2);
    CHECK(std::isnan(alone.zdiff));
    CHECK(std::isnan(alone.zdiff_robust));
    CHECK(alone.diff == doctest::Approx(n.diff));

    // A single voxel has no neighbours: no spread, no differences.
    const SliceNoise single = denoiser::sliceNoise(slice.data(), nullptr, 1, 1);
    CHECK(single.std_dev == 0.0);
    CHECK(single.diff == 0.0);
    CHECK(single.resid == 0.0);

    CHECK_THROWS_AS(denoiser::sliceNoise(slice.data(), nullptr, 0, 2), std::invalid_argument);
}

TEST_CASE("sliceNoise: independent noise without material") {
    // 200 x 200 voxels of noise with sigma 10 around the value 100.
    const int n = 200;
    Noise noise(1);
    std::vector<float> slice(static_cast<std::size_t>(n) * n);
    std::vector<float> next(slice.size());
    for (float& v : slice) {
        v = 100.0f + 10.0f * static_cast<float>(noise.next());
    }
    for (float& v : next) {
        v = 100.0f + 10.0f * static_cast<float>(noise.next());
    }
    const SliceNoise r = denoiser::sliceNoise(slice.data(), next.data(), n, n);
    CHECK(r.std_dev == doctest::Approx(10.0).epsilon(0.02));
    CHECK(r.mad == doctest::Approx(10.0).epsilon(0.03));
    CHECK(r.diff == doctest::Approx(10.0).epsilon(0.02));
    CHECK(r.zdiff == doctest::Approx(10.0).epsilon(0.02));
    CHECK(r.zdiff_robust == doctest::Approx(10.0).epsilon(0.03));
    // The 3 x 3 mean contains the voxel itself, so the residual keeps only
    // 8/9 of the noise variance: this estimator reads 5.7 % low.
    CHECK(r.resid == doctest::Approx(10.0 * std::sqrt(8.0 / 9.0)).epsilon(0.03));
}

TEST_CASE("sliceNoise: material that is the same in both slices") {
    // Stripes with a contrast of 60 plus noise with sigma 10. The estimators
    // within the slice see the stripes; those along z do not.
    const int n = 200;
    Noise noise(2);
    std::vector<float> slice(static_cast<std::size_t>(n) * n);
    std::vector<float> next(slice.size());
    for (int y = 0; y < n; ++y) {
        for (int x = 0; x < n; ++x) {
            const std::size_t i = static_cast<std::size_t>(y) * n + static_cast<std::size_t>(x);
            slice[i] = pattern(x, 60.0f) + 10.0f * static_cast<float>(noise.next());
            next[i] = pattern(x, 60.0f) + 10.0f * static_cast<float>(noise.next());
        }
    }
    const SliceNoise r = denoiser::sliceNoise(slice.data(), next.data(), n, n);
    CHECK(r.std_dev > 30.0);
    CHECK(r.diff > 14.0);
    CHECK(r.resid > 12.0);
    CHECK(r.zdiff == doctest::Approx(10.0).epsilon(0.02));
    CHECK(r.zdiff_robust == doctest::Approx(10.0).epsilon(0.03));
}

TEST_CASE("estimateNoise: slices, summaries and the earlier recommendation") {
    // 48 x 40 x 30 voxels: stripes plus independent noise with sigma 5.
    const VolumeSize size{48, 40, 30};
    Noise noise(3);
    std::vector<float> values;
    for (int z = 0; z < size.depth; ++z) {
        for (int y = 0; y < size.height; ++y) {
            for (int x = 0; x < size.width; ++x) {
                values.push_back(pattern(x, 20.0f) + 5.0f * static_cast<float>(noise.next()));
            }
        }
    }
    const ScanInfo info = describe("stripes.vol", size);
    denoiser::writeVolume(info, values);

    const NoiseEstimate e = denoiser::estimateNoise(info, denoiser::wholeVolume(size), 20);
    REQUIRE(e.slices.size() == 20);
    REQUIRE(e.per_slice.size() == 20);
    CHECK(e.slices.front() == 0);
    CHECK(e.slices.back() == 29);
    // The last slice has no next slice.
    CHECK(std::isnan(e.per_slice.back().zdiff));
    CHECK(e.diff.count == 20);
    CHECK(e.zdiff.count == 19);
    CHECK(e.zdiff_robust.count == 19);
    CHECK(e.recommended == doctest::Approx(0.5 * (e.diff.median + e.resid.median)));
    CHECK(e.zdiff_robust.median == doctest::Approx(5.0).epsilon(0.04));
    CHECK(e.recommended > 6.0);  // the stripes are taken for noise

    // A region: slice numbers are those of the scan; the first slice of the
    // region gives the same values as the same slice of the whole scan.
    const Region part{8, 0, 10, 32, 40, 12};
    const NoiseEstimate p = denoiser::estimateNoise(info, part, 3);
    CHECK(p.slices == std::vector<int>{10, 15, 21});
    const std::vector<float> slice10 = denoiser::readRegion(info, Region{8, 0, 10, 32, 40, 1});
    const std::vector<float> slice11 = denoiser::readRegion(info, Region{8, 0, 11, 32, 40, 1});
    const SliceNoise direct = denoiser::sliceNoise(slice10.data(), slice11.data(), 32, 40);
    CHECK(p.per_slice[0].diff == direct.diff);
    CHECK(p.per_slice[0].zdiff_robust == direct.zdiff_robust);

    CHECK_THROWS_AS(denoiser::estimateNoise(info, Region{0, 0, 0, 49, 40, 30}, 20),
                    std::invalid_argument);
    CHECK_THROWS_AS(denoiser::estimateNoise(info, part, 0), std::invalid_argument);
}

TEST_CASE("estimateNoise refuses NaN") {
    const VolumeSize size{4, 4, 2};
    std::vector<float> values(32, 1.0f);
    values[5] = std::numeric_limits<float>::quiet_NaN();
    const ScanInfo info = describe("nan.vol", size);
    denoiser::writeVolume(info, values);
    CHECK_THROWS_AS(denoiser::estimateNoise(info, denoiser::wholeVolume(size), 2),
                    std::runtime_error);
}

TEST_CASE("estimateNoise on a tracked test scan") {
    // Values of step 1.4, in counts; they guard against changes of the
    // method and must be the same on every platform.
    const ScanInfo info = denoiser::openScan("tests/data/Sample2_5mu_c64");
    const NoiseEstimate e = denoiser::estimateNoise(info, denoiser::wholeVolume(info.size), 20);
    CHECK(e.recommended == doctest::Approx(2999.880).epsilon(1e-6));
    CHECK(e.zdiff_robust.median == doctest::Approx(3092.128).epsilon(1e-6));
}

TEST_CASE("differenceByLag shows noise that is correlated between neighbours") {
    // Stripes that do not change along z, plus noise with sigma 10 that is
    // correlated along z only: independent noise w is replaced by
    // w[z] + a * (w[z-1] + w[z+1]) with a = 0.25. Its correlation is
    // 2a / (1 + 2a^2) = 0.444 at distance 1, a^2 / (1 + 2a^2) = 0.056 at
    // distance 2 and 0 beyond.
    const VolumeSize size{64, 64, 40};
    const double a = 0.25;
    const double scale = 10.0 / std::sqrt(1.0 + 2.0 * a * a);
    const std::size_t slice_voxels = 64 * 64;
    Noise noise(4);
    std::vector<double> white(slice_voxels * static_cast<std::size_t>(size.depth + 2));
    for (double& v : white) {
        v = noise.next();
    }
    std::vector<float> values;
    for (int z = 0; z < size.depth; ++z) {
        for (std::size_t i = 0; i < slice_voxels; ++i) {
            const std::size_t centre = static_cast<std::size_t>(z + 1) * slice_voxels + i;
            const double n =
                white[centre] + a * (white[centre - slice_voxels] + white[centre + slice_voxels]);
            values.push_back(pattern(static_cast<int>(i % 64), 60.0f) +
                             static_cast<float>(scale * n));
        }
    }
    const ScanInfo info = describe("correlated.vol", size);
    denoiser::writeVolume(info, values);
    const Region whole = denoiser::wholeVolume(size);

    const std::vector<denoiser::LagSigma> lags = denoiser::differenceByLag(info, whole, 20, 5);
    REQUIRE(lags.size() == 5);
    CHECK(lags[0].lag == 1);
    CHECK(lags[4].lag == 5);
    // Along z: too low at distance 1, nearly right at 2, right from 3 on.
    CHECK(lags[0].z == doctest::Approx(10.0 * std::sqrt(1.0 - 0.4444)).epsilon(0.04));
    CHECK(lags[1].z == doctest::Approx(10.0 * std::sqrt(1.0 - 0.0556)).epsilon(0.04));
    CHECK(lags[2].z == doctest::Approx(10.0).epsilon(0.04));
    CHECK(lags[3].z == doctest::Approx(10.0).epsilon(0.04));
    CHECK(lags[4].z == doctest::Approx(10.0).epsilon(0.04));
    // Along y nothing but noise changes, and that noise is independent.
    CHECK(lags[0].y == doctest::Approx(10.0).epsilon(0.04));
    // Along x the stripes are in the differences.
    CHECK(lags[1].x > 20.0);

    // The neighbour estimate of estimateNoise is the entry at distance 1.
    const NoiseEstimate e = denoiser::estimateNoise(info, whole, 20);
    CHECK(e.zdiff_robust.median == doctest::Approx(lags[0].z).epsilon(0.04));

    // A region narrower than the distance gives no value along that axis.
    const std::vector<denoiser::LagSigma> narrow =
        denoiser::differenceByLag(info, Region{0, 0, 0, 2, 64, 40}, 5, 3);
    CHECK(!std::isnan(narrow[0].x));
    CHECK(std::isnan(narrow[1].x));
    CHECK(!std::isnan(narrow[2].y));

    CHECK_THROWS_AS(denoiser::differenceByLag(info, whole, 20, 0), std::invalid_argument);
    CHECK_THROWS_AS(denoiser::differenceByLag(info, whole, 20, 40), std::invalid_argument);
    CHECK_THROWS_AS(denoiser::differenceByLag(info, whole, 0, 3), std::invalid_argument);
}
