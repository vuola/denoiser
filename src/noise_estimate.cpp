// denoiser - 3D Non-Local Means filter for CT volumes.
// Derived from the C++ filter by Samuel Väisänen (2024), itself a port of the
// ImageJ plugin by Pascal Behnel and Thorsten Wagner. GPL-3.0, see LICENSE.

#include "noise_estimate.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>

namespace denoiser {

namespace {

/// 1.4826 * MAD is the standard deviation of normally distributed values.
const double kMadToSigma = 1.4826;

const double kNaN = std::numeric_limits<double>::quiet_NaN();

/// Value at a fractional position of the sorted values, with linear
/// interpolation between the two nearest (numpy's default percentile rule).
/// Reorders `values`. `values` must not be empty.
double sortedValueAt(std::vector<double>& values, double position) {
    const std::size_t below = static_cast<std::size_t>(std::floor(position));
    const double fraction = position - static_cast<double>(below);
    std::nth_element(values.begin(), values.begin() + static_cast<std::ptrdiff_t>(below),
                     values.end());
    const double low = values[below];
    if (fraction == 0.0 || below + 1 >= values.size()) {
        return low;
    }
    // After nth_element everything behind `below` is >= low; the next sorted
    // value is the smallest of those.
    const double high =
        *std::min_element(values.begin() + static_cast<std::ptrdiff_t>(below + 1), values.end());
    return low + fraction * (high - low);
}

/// Median; the mean of the two middle values for an even count. Reorders
/// `values`. `values` must not be empty.
double median(std::vector<double>& values) {
    return sortedValueAt(values, 0.5 * static_cast<double>(values.size() - 1));
}

/// Median of |v - median(v)|. Changes `values`. `values` must not be empty.
double medianAbsoluteDeviation(std::vector<double>& values) {
    const double centre = median(values);
    for (double& v : values) {
        v = std::fabs(v - centre);
    }
    return median(values);
}

/// Standard deviation (dividing by the count). 0 for an empty vector.
double standardDeviation(const std::vector<double>& values) {
    if (values.empty()) {
        return 0.0;
    }
    double sum = 0.0;
    for (double v : values) {
        sum += v;
    }
    const double mean = sum / static_cast<double>(values.size());
    double squares = 0.0;
    for (double v : values) {
        squares += (v - mean) * (v - mean);
    }
    return std::sqrt(squares / static_cast<double>(values.size()));
}

/// Differences between voxels `lag` apart along x (step 1) or y (step
/// width) within one slice.
std::vector<double> differencesInSlice(const float* slice, int width, int height, int lag,
                                       bool along_x) {
    std::vector<double> diffs;
    const int nx = along_x ? width - lag : width;
    const int ny = along_x ? height : height - lag;
    if (nx <= 0 || ny <= 0) {
        return diffs;
    }
    diffs.reserve(static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny));
    const std::ptrdiff_t step = along_x ? lag : static_cast<std::ptrdiff_t>(lag) * width;
    for (int y = 0; y < ny; ++y) {
        const float* row = slice + static_cast<std::ptrdiff_t>(y) * width;
        for (int x = 0; x < nx; ++x) {
            diffs.push_back(static_cast<double>(row[x + step]) - static_cast<double>(row[x]));
        }
    }
    return diffs;
}

/// Differences between two slices, voxel by voxel.
std::vector<double> differencesBetween(const float* slice, const float* next, std::size_t count) {
    std::vector<double> diffs(count);
    for (std::size_t i = 0; i < count; ++i) {
        diffs[i] = static_cast<double>(next[i]) - static_cast<double>(slice[i]);
    }
    return diffs;
}

/// Robust sigma from differences: 1.4826 * MAD / sqrt(2). NaN if empty.
double robustSigma(std::vector<double> diffs) {
    if (diffs.empty()) {
        return kNaN;
    }
    return kMadToSigma * medianAbsoluteDeviation(diffs) / std::sqrt(2.0);
}

/// Reads one slice of the region and refuses NaN and infinite values.
std::vector<float> readSlice(const ScanInfo& info, const Region& region, int z) {
    const Region one{region.x0, region.y0, region.z0 + z, region.width, region.height, 1};
    std::vector<float> slice = readRegion(info, one);
    for (float v : slice) {
        if (!std::isfinite(v)) {
            throw std::runtime_error(info.volume_path + ": slice " +
                                     std::to_string(region.z0 + z) +
                                     " contains NaN or infinite values");
        }
    }
    return slice;
}

void checkRegion(const ScanInfo& info, const Region& region, int slices) {
    if (!regionFits(region, info.size)) {
        throw std::invalid_argument("the region does not fit the scan");
    }
    if (slices < 1) {
        throw std::invalid_argument("the number of slices must be at least 1");
    }
}

}  // namespace

SliceNoise sliceNoise(const float* slice, const float* next, int width, int height) {
    if (width < 1 || height < 1) {
        throw std::invalid_argument("sliceNoise: width and height must be at least 1");
    }
    const std::size_t count = static_cast<std::size_t>(width) * static_cast<std::size_t>(height);
    SliceNoise result;

    std::vector<double> values(slice, slice + count);
    result.std_dev = standardDeviation(values);
    result.mad = kMadToSigma * medianAbsoluteDeviation(values);

    // Neighbour differences along x and along y, taken together.
    std::vector<double> diffs = differencesInSlice(slice, width, height, 1, true);
    const std::vector<double> diffs_y = differencesInSlice(slice, width, height, 1, false);
    diffs.insert(diffs.end(), diffs_y.begin(), diffs_y.end());
    result.diff = standardDeviation(diffs) / std::sqrt(2.0);

    // Residual against the mean of the 3 x 3 neighbourhood, centre included.
    // The nine values are added in the order (dy, dx) = (-1,-1), (-1,0), ...
    std::vector<double> residual(count);
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            double sum = 0.0;
            for (int dy = -1; dy <= 1; ++dy) {
                const int yy = std::min(std::max(y + dy, 0), height - 1);
                const float* row = slice + static_cast<std::ptrdiff_t>(yy) * width;
                for (int dx = -1; dx <= 1; ++dx) {
                    const int xx = std::min(std::max(x + dx, 0), width - 1);
                    sum += static_cast<double>(row[xx]);
                }
            }
            const std::size_t index =
                static_cast<std::size_t>(y) * static_cast<std::size_t>(width) +
                static_cast<std::size_t>(x);
            residual[index] = static_cast<double>(slice[index]) - sum / 9.0;
        }
    }
    result.resid = kMadToSigma * medianAbsoluteDeviation(residual);

    if (next != nullptr) {
        std::vector<double> dz = differencesBetween(slice, next, count);
        result.zdiff = standardDeviation(dz) / std::sqrt(2.0);
        result.zdiff_robust = robustSigma(std::move(dz));
    } else {
        result.zdiff = kNaN;
        result.zdiff_robust = kNaN;
    }
    return result;
}

NoiseSummary summarise(const std::vector<double>& values) {
    std::vector<double> finite;
    finite.reserve(values.size());
    double sum = 0.0;
    for (double v : values) {
        if (!std::isnan(v)) {
            finite.push_back(v);
            sum += v;
        }
    }
    NoiseSummary summary;
    summary.count = static_cast<int>(finite.size());
    if (finite.empty()) {
        summary.median = summary.mean = summary.q1 = summary.q3 = kNaN;
        return summary;
    }
    const double last = static_cast<double>(finite.size() - 1);
    summary.mean = sum / static_cast<double>(finite.size());
    summary.median = sortedValueAt(finite, 0.5 * last);
    summary.q1 = sortedValueAt(finite, 0.25 * last);
    summary.q3 = sortedValueAt(finite, 0.75 * last);
    return summary;
}

std::vector<int> sampleSlices(int depth, int count) {
    if (depth < 1 || count < 1) {
        throw std::invalid_argument("sampleSlices: depth and count must be at least 1");
    }
    count = std::min(count, depth);
    std::vector<int> slices(static_cast<std::size_t>(count), 0);
    for (int i = 1; i < count; ++i) {
        slices[static_cast<std::size_t>(i)] = static_cast<int>(
            static_cast<std::int64_t>(i) * (depth - 1) / (count - 1));
    }
    return slices;
}

NoiseEstimate estimateNoise(const ScanInfo& info, const Region& region, int slices) {
    checkRegion(info, region, slices);

    NoiseEstimate estimate;
    const std::vector<int> picked = sampleSlices(region.depth, slices);
    for (int z : picked) {
        const std::vector<float> slice = readSlice(info, region, z);
        std::vector<float> next;
        if (z + 1 < region.depth) {
            next = readSlice(info, region, z + 1);
        }
        estimate.slices.push_back(region.z0 + z);
        estimate.per_slice.push_back(sliceNoise(slice.data(), next.empty() ? nullptr : next.data(),
                                                region.width, region.height));
    }

    const auto summary_of = [&estimate](double SliceNoise::* member) {
        std::vector<double> values;
        for (const SliceNoise& s : estimate.per_slice) {
            values.push_back(s.*member);
        }
        return summarise(values);
    };
    estimate.std_dev = summary_of(&SliceNoise::std_dev);
    estimate.mad = summary_of(&SliceNoise::mad);
    estimate.diff = summary_of(&SliceNoise::diff);
    estimate.resid = summary_of(&SliceNoise::resid);
    estimate.zdiff = summary_of(&SliceNoise::zdiff);
    estimate.zdiff_robust = summary_of(&SliceNoise::zdiff_robust);
    estimate.recommended = 0.5 * (estimate.diff.median + estimate.resid.median);
    return estimate;
}

std::vector<LagSigma> differenceByLag(const ScanInfo& info, const Region& region, int slices,
                                      int max_lag) {
    checkRegion(info, region, slices);
    if (max_lag < 1) {
        throw std::invalid_argument("the largest lag must be at least 1");
    }
    if (region.depth <= max_lag) {
        throw std::invalid_argument("the region must have more slices than the largest lag");
    }

    // Slices are picked so that every one of them has max_lag slices behind it.
    const std::vector<int> picked = sampleSlices(region.depth - max_lag, slices);
    const std::size_t lags = static_cast<std::size_t>(max_lag);
    std::vector<std::vector<double>> x(lags), y(lags), z(lags);
    const std::size_t count =
        static_cast<std::size_t>(region.width) * static_cast<std::size_t>(region.height);

    for (int slice_z : picked) {
        const std::vector<float> slice = readSlice(info, region, slice_z);
        for (int lag = 1; lag <= max_lag; ++lag) {
            const std::size_t k = static_cast<std::size_t>(lag - 1);
            x[k].push_back(robustSigma(
                differencesInSlice(slice.data(), region.width, region.height, lag, true)));
            y[k].push_back(robustSigma(
                differencesInSlice(slice.data(), region.width, region.height, lag, false)));
            const std::vector<float> other = readSlice(info, region, slice_z + lag);
            z[k].push_back(robustSigma(differencesBetween(slice.data(), other.data(), count)));
        }
    }

    std::vector<LagSigma> result;
    for (int lag = 1; lag <= max_lag; ++lag) {
        const std::size_t k = static_cast<std::size_t>(lag - 1);
        result.push_back(LagSigma{lag, summarise(x[k]).median, summarise(y[k]).median,
                                  summarise(z[k]).median});
    }
    return result;
}

}  // namespace denoiser
