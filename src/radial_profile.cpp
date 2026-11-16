// denoiser - 3D Non-Local Means filter for CT volumes.
// Derived from the C++ filter by Samuel Väisänen (2024), itself a port of the
// ImageJ plugin by Pascal Behnel and Thorsten Wagner. GPL-3.0, see LICENSE.

#include "radial_profile.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>

#include "volume.hpp"

namespace denoiser {

namespace {

const int kClasses = 1024;  ///< histogram classes per ring

/// Returns the value below which `share` of the histogram's entries lie.
/// `counts` has kClasses entries for the range low...high and `total` is
/// their sum (> 0). The position inside a class is interpolated linearly.
double percentile(const std::int64_t* counts, std::int64_t total, double share, double low,
                  double high) {
    const double target = share * static_cast<double>(total);
    const double class_width = (high - low) / kClasses;
    double below = 0.0;
    for (int c = 0; c < kClasses; ++c) {
        const double here = static_cast<double>(counts[c]);
        if (here > 0.0 && below + here >= target) {
            return low + (c + (target - below) / here) * class_width;
        }
        below += here;
    }
    return high;
}

}  // namespace

ProfileSettings defaultProfileSettings(const ScanInfo& info) {
    ProfileSettings settings;
    settings.centre_x = (info.size.width - 1) / 2.0;
    settings.centre_y = (info.size.height - 1) / 2.0;
    settings.value_high = info.sample_type == SampleType::kUInt8 ? 255.0 : 65535.0;
    return settings;
}

std::vector<RadialRing> radialProfile(const ScanInfo& info, const ProfileSettings& settings) {
    if (!(settings.ring_width > 0.0) || !std::isfinite(settings.ring_width)) {
        throw std::invalid_argument("the ring width must be a positive number");
    }
    if (!std::isfinite(settings.centre_x) || !std::isfinite(settings.centre_y)) {
        throw std::invalid_argument("the centre must be two finite numbers");
    }
    if (!std::isfinite(settings.value_low) || !std::isfinite(settings.value_high) ||
        !(settings.value_high > settings.value_low)) {
        throw std::invalid_argument("the value range needs two finite numbers, the second larger than the first");
    }
    const int depth = settings.depth == 0 ? info.size.depth - settings.z0 : settings.depth;
    if (settings.z0 < 0 || depth < 1 || depth > info.size.depth - settings.z0) {
        throw std::invalid_argument("the slices of the profile do not fit the volume");
    }

    const int width = info.size.width;
    const int height = info.size.height;

    // The ring of every (x, y) position; the same for all slices.
    std::vector<int> ring_of(static_cast<std::size_t>(width) * height);
    int rings = 0;
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            const double r = std::hypot(x - settings.centre_x, y - settings.centre_y);
            const int ring = static_cast<int>(r / settings.ring_width);
            ring_of[static_cast<std::size_t>(y) * width + x] = ring;
            rings = std::max(rings, ring + 1);
        }
    }

    std::vector<std::int64_t> counts(static_cast<std::size_t>(rings) * kClasses, 0);
    std::vector<double> sums(static_cast<std::size_t>(rings), 0.0);
    std::vector<RadialRing> result(static_cast<std::size_t>(rings));

    const double scale = kClasses / (settings.value_high - settings.value_low);
    for (int z = settings.z0; z < settings.z0 + depth; ++z) {
        const std::vector<float> slice = readRegion(info, Region{0, 0, z, width, height, 1});
        for (std::size_t i = 0; i < slice.size(); ++i) {
            const double value = slice[i];
            if (!std::isfinite(value)) {
                continue;
            }
            const std::size_t ring = static_cast<std::size_t>(ring_of[i]);
            RadialRing& entry = result[ring];
            ++entry.voxels;
            sums[ring] += value;
            if (value <= settings.value_low) {
                ++entry.at_low;
            }
            if (value >= settings.value_high) {
                ++entry.at_high;
            }
            const int cls = std::clamp(static_cast<int>((value - settings.value_low) * scale), 0,
                                       kClasses - 1);
            ++counts[ring * kClasses + static_cast<std::size_t>(cls)];
        }
    }

    for (int ring = 0; ring < rings; ++ring) {
        RadialRing& entry = result[static_cast<std::size_t>(ring)];
        entry.r_inner = ring * settings.ring_width;
        entry.r_outer = (ring + 1) * settings.ring_width;
        if (entry.voxels == 0) {
            continue;
        }
        const std::int64_t* ring_counts = &counts[static_cast<std::size_t>(ring) * kClasses];
        entry.mean = sums[static_cast<std::size_t>(ring)] / static_cast<double>(entry.voxels);
        entry.p10 = percentile(ring_counts, entry.voxels, 0.10, settings.value_low, settings.value_high);
        entry.p50 = percentile(ring_counts, entry.voxels, 0.50, settings.value_low, settings.value_high);
        entry.p90 = percentile(ring_counts, entry.voxels, 0.90, settings.value_low, settings.value_high);
    }
    return result;
}

}  // namespace denoiser
