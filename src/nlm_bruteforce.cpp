// denoiser - 3D Non-Local Means filter for CT volumes.
// Derived from the C++ filter by Samuel Väisänen (2024), itself a port of the
// ImageJ plugin by Pascal Behnel and Thorsten Wagner. GPL-3.0, see LICENSE.

#include "nlm_bruteforce.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace denoiser {

NlmParams originalParams(double sigma) {
    if (!(sigma > 0.0) || !std::isfinite(sigma)) {
        throw std::invalid_argument("sigma must be a positive finite number");
    }

    // The original program holds sigma and h as 32-bit floats; the same
    // rounding is applied here so that both use identical settings.
    const float sigma_f = static_cast<float>(sigma);
    float h_factor = 0.0f;
    int patch = 0;
    int search = 0;
    if (sigma_f <= 15.0f) {
        h_factor = 0.4f;  patch = 1; search = 10;
    } else if (sigma_f <= 30.0f) {
        h_factor = 0.4f;  patch = 2; search = 10;
    } else if (sigma_f <= 45.0f) {
        h_factor = 0.35f; patch = 3; search = 17;
    } else if (sigma_f <= 75.0f) {
        h_factor = 0.35f; patch = 4; search = 17;
    } else {
        h_factor = 0.3f;  patch = 5; search = 17;
    }

    NlmParams params;
    params.sigma = sigma_f;
    params.h = h_factor * sigma_f;  // float product, as in the original
    for (int axis = 0; axis < 3; ++axis) {
        params.patch_radius[axis] = patch;
        params.search_radius[axis] = search;
    }
    return params;
}

double patchVoxels(const NlmParams& params) {
    return static_cast<double>(2 * params.patch_radius[0] + 1) *
           (2 * params.patch_radius[1] + 1) * (2 * params.patch_radius[2] + 1);
}

double nlmWeight(double distance, const NlmParams& params) {
    const double voxels = patchVoxels(params);
    const double noise = 2.0 * params.sigma * params.sigma;
    const double excess =
        params.normalise_distance ? distance / voxels - noise : distance - noise * voxels;
    if (excess <= 0.0) {
        return 1.0;
    }
    const double h2 = params.h * params.h;
    return params.kernel == NlmKernel::kRational ? h2 / (h2 + excess) : std::exp(-excess / h2);
}

void checkNlmInput(const std::vector<float>& volume, const VolumeSize& size,
                   const NlmParams& params) {
    if (size.width < 1 || size.height < 1 || size.depth < 1) {
        throw std::invalid_argument("volume dimensions must be at least 1");
    }
    const std::size_t voxels = static_cast<std::size_t>(size.width) * size.height * size.depth;
    if (volume.size() != voxels) {
        throw std::invalid_argument("volume length does not match its size");
    }
    if (!(params.sigma > 0.0) || !(params.h > 0.0)) {
        throw std::invalid_argument("sigma and h must be positive");
    }
    for (int axis = 0; axis < 3; ++axis) {
        if (params.patch_radius[axis] < 0 || params.search_radius[axis] < 0) {
            throw std::invalid_argument("radii must not be negative");
        }
    }
}

std::vector<float> nlmBruteForce(const std::vector<float>& volume, const VolumeSize& size,
                                 const NlmParams& params) {
    checkNlmInput(volume, size, params);

    const int px = params.patch_radius[0], py = params.patch_radius[1], pz = params.patch_radius[2];
    const int sx = params.search_radius[0], sy = params.search_radius[1], sz = params.search_radius[2];

    // Copy the volume into a larger array with replicated edges, wide enough
    // for the farthest voxel any patch of any shifted position can touch.
    const int bx = px + sx, by = py + sy, bz = pz + sz;
    const std::size_t pw = static_cast<std::size_t>(size.width) + 2 * bx;
    const std::size_t ph = static_cast<std::size_t>(size.height) + 2 * by;
    const std::size_t pd = static_cast<std::size_t>(size.depth) + 2 * bz;
    std::vector<float> padded(pw * ph * pd);
    for (std::size_t z = 0; z < pd; ++z) {
        const int src_z = std::clamp(static_cast<int>(z) - bz, 0, size.depth - 1);
        for (std::size_t y = 0; y < ph; ++y) {
            const int src_y = std::clamp(static_cast<int>(y) - by, 0, size.height - 1);
            for (std::size_t x = 0; x < pw; ++x) {
                const int src_x = std::clamp(static_cast<int>(x) - bx, 0, size.width - 1);
                padded[(z * ph + y) * pw + x] =
                    volume[(static_cast<std::size_t>(src_z) * size.height + src_y) * size.width + src_x];
            }
        }
    }

    // The weight in terms of the summed distance, for every combination of
    // settings: excess = distance * scale - offset.
    const double patch_voxels = patchVoxels(params);
    const double scale = params.normalise_distance ? 1.0 / patch_voxels : 1.0;
    const double offset =
        2.0 * params.sigma * params.sigma * (params.normalise_distance ? 1.0 : patch_voxels);
    const double h2 = params.h * params.h;
    const bool rational = params.kernel == NlmKernel::kRational;

    // Index steps in the padded array.
    const std::ptrdiff_t step_y = static_cast<std::ptrdiff_t>(pw);
    const std::ptrdiff_t step_z = static_cast<std::ptrdiff_t>(pw * ph);

    std::vector<float> result(volume.size());

    #pragma omp parallel for schedule(dynamic)
    for (int z = 0; z < size.depth; ++z) {
        for (int y = 0; y < size.height; ++y) {
            for (int x = 0; x < size.width; ++x) {
                const float* centre = padded.data() + (z + bz) * step_z + (y + by) * step_y + (x + bx);

                double weighted_sum = 0.0;
                double weight_sum = 0.0;
                double max_weight = -1.0;

                for (int dz = -sz; dz <= sz; ++dz) {
                    for (int dy = -sy; dy <= sy; ++dy) {
                        for (int dx = -sx; dx <= sx; ++dx) {
                            if (dx == 0 && dy == 0 && dz == 0) {
                                continue;  // the voxel itself is added after the loop
                            }
                            const float* other = centre + dz * step_z + dy * step_y + dx;

                            double distance = 0.0;
                            for (int qz = -pz; qz <= pz; ++qz) {
                                for (int qy = -py; qy <= py; ++qy) {
                                    const std::ptrdiff_t row = qz * step_z + qy * step_y;
                                    for (int qx = -px; qx <= px; ++qx) {
                                        const double difference =
                                            static_cast<double>(centre[row + qx]) - other[row + qx];
                                        distance += difference * difference;
                                    }
                                }
                            }

                            const double excess = distance * scale - offset;
                            const double weight =
                                excess <= 0.0 ? 1.0
                                : rational    ? h2 / (h2 + excess)
                                              : std::exp(-excess / h2);
                            weighted_sum += weight * *other;
                            weight_sum += weight;
                            max_weight = std::max(max_weight, weight);
                        }
                    }
                }

                if (max_weight < 0.0) {
                    // Search window of one voxel: nothing to average with.
                    max_weight = 1.0;
                }
                weighted_sum += max_weight * *centre;
                weight_sum += max_weight;
                result[(static_cast<std::size_t>(z) * size.height + y) * size.width + x] =
                    static_cast<float>(weighted_sum / weight_sum);
            }
        }
    }
    return result;
}

}  // namespace denoiser
