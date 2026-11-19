// denoiser - 3D Non-Local Means filter for CT volumes.
// Derived from the C++ filter by Samuel Väisänen (2024), itself a port of the
// ImageJ plugin by Pascal Behnel and Thorsten Wagner. GPL-3.0, see LICENSE.

#include "nlm_fast.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace denoiser {

namespace {

/// Running sums along x. `src` holds rows of out_w + 2r values; `dst`
/// receives rows of out_w values, each the sum of 2r + 1 neighbours.
template <typename Real>
void sumAlongX(const Real* src, std::size_t rows, int out_w, int r, Real* dst) {
    const std::size_t src_w = static_cast<std::size_t>(out_w) + 2 * static_cast<std::size_t>(r);
    for (std::size_t row = 0; row < rows; ++row) {
        const Real* s = src + row * src_w;
        Real* d = dst + row * static_cast<std::size_t>(out_w);
        Real sum = 0;
        for (int k = 0; k <= 2 * r; ++k) {
            sum += s[k];
        }
        d[0] = sum;
        for (int i = 1; i < out_w; ++i) {
            sum += s[i + 2 * r] - s[i - 1];
            d[i] = sum;
        }
    }
}

/// Running sums across whole rows: for y (a row is one line of x values,
/// `stacks` = number of slices) and for z (a row is one slice, stacks = 1).
/// Each stack of `src` holds out_rows + 2r rows of `row_len` values; each
/// output row is the sum of 2r + 1 neighbouring input rows.
template <typename Real>
void sumAcrossRows(const Real* src, std::size_t row_len, int out_rows, int r,
                   std::size_t stacks, Real* dst) {
    const std::size_t src_rows = static_cast<std::size_t>(out_rows) + 2 * static_cast<std::size_t>(r);
    for (std::size_t stack = 0; stack < stacks; ++stack) {
        const Real* s = src + stack * src_rows * row_len;
        Real* d = dst + stack * static_cast<std::size_t>(out_rows) * row_len;
        for (std::size_t x = 0; x < row_len; ++x) {
            d[x] = 0;
        }
        for (int k = 0; k <= 2 * r; ++k) {
            const Real* in = s + static_cast<std::size_t>(k) * row_len;
            for (std::size_t x = 0; x < row_len; ++x) {
                d[x] += in[x];
            }
        }
        for (int i = 1; i < out_rows; ++i) {
            const Real* previous = d + static_cast<std::size_t>(i - 1) * row_len;
            const Real* added = s + static_cast<std::size_t>(i + 2 * r) * row_len;
            const Real* removed = s + static_cast<std::size_t>(i - 1) * row_len;
            Real* out = d + static_cast<std::size_t>(i) * row_len;
            for (std::size_t x = 0; x < row_len; ++x) {
                out[x] = previous[x] + (added[x] - removed[x]);
            }
        }
    }
}

/// Buffers of one block, kept between blocks of the same thread.
template <typename Real>
struct BlockBuffers {
    std::vector<float> input;     // the block with a margin of patch + search radius
    std::vector<Real> squares;    // squared differences on the block plus patch radius
    std::vector<Real> along_x;    // after the sum along x
    std::vector<Real> along_y;    // after the sum along y
    std::vector<Real> distance;   // patch distance of every block voxel
    std::vector<double> weighted_sum;
    std::vector<double> weight_sum;
    std::vector<double> max_weight;
};

/// Filters the voxels of one block and writes them into `result`.
template <typename Real>
void filterBlock(const std::vector<float>& volume, const VolumeSize& size, const NlmParams& params,
                 const Region& block, BlockBuffers<Real>& buf, std::vector<float>& result) {
    const int px = params.patch_radius[0], py = params.patch_radius[1], pz = params.patch_radius[2];
    const int sx = params.search_radius[0], sy = params.search_radius[1], sz = params.search_radius[2];
    const int bw = block.width, bh = block.height, bd = block.depth;

    // Input with a margin of patch + search radius, edges replicated as in
    // nlmBruteForce().
    const int mx = px + sx, my = py + sy, mz = pz + sz;
    const std::size_t iw = static_cast<std::size_t>(bw + 2 * mx);
    const std::size_t ih = static_cast<std::size_t>(bh + 2 * my);
    const std::size_t id = static_cast<std::size_t>(bd + 2 * mz);
    buf.input.resize(iw * ih * id);
    for (std::size_t z = 0; z < id; ++z) {
        const int src_z = std::clamp(block.z0 + static_cast<int>(z) - mz, 0, size.depth - 1);
        for (std::size_t y = 0; y < ih; ++y) {
            const int src_y = std::clamp(block.y0 + static_cast<int>(y) - my, 0, size.height - 1);
            const float* row = volume.data() +
                               (static_cast<std::size_t>(src_z) * size.height + src_y) * size.width;
            float* out = buf.input.data() + (z * ih + y) * iw;
            for (std::size_t x = 0; x < iw; ++x) {
                const int src_x = std::clamp(block.x0 + static_cast<int>(x) - mx, 0, size.width - 1);
                out[x] = row[src_x];
            }
        }
    }

    // Squared differences are needed on the block plus the patch radius.
    const std::size_t ew = static_cast<std::size_t>(bw + 2 * px);
    const std::size_t eh = static_cast<std::size_t>(bh + 2 * py);
    const std::size_t ed = static_cast<std::size_t>(bd + 2 * pz);
    const std::size_t block_voxels = static_cast<std::size_t>(bw) * bh * bd;
    buf.squares.resize(ew * eh * ed);
    buf.along_x.resize(static_cast<std::size_t>(bw) * eh * ed);
    buf.along_y.resize(static_cast<std::size_t>(bw) * bh * ed);
    buf.distance.resize(block_voxels);
    buf.weighted_sum.assign(block_voxels, 0.0);
    buf.weight_sum.assign(block_voxels, 0.0);
    buf.max_weight.assign(block_voxels, -1.0);

    // The weight as a function of the summed distance, as in nlmWeight():
    // excess = distance * scale - offset.
    const double patch_voxels = patchVoxels(params);
    const double scale = params.normalise_distance ? 1.0 / patch_voxels : 1.0;
    const double offset =
        2.0 * params.sigma * params.sigma * (params.normalise_distance ? 1.0 : patch_voxels);
    const double h2 = params.h * params.h;
    const bool rational = params.kernel == NlmKernel::kRational;

    const std::ptrdiff_t step_y = static_cast<std::ptrdiff_t>(iw);
    const std::ptrdiff_t step_z = static_cast<std::ptrdiff_t>(iw * ih);

    for (int dz = -sz; dz <= sz; ++dz) {
        for (int dy = -sy; dy <= sy; ++dy) {
            for (int dx = -sx; dx <= sx; ++dx) {
                if (dx == 0 && dy == 0 && dz == 0) {
                    continue;  // the voxel itself is added after the loop
                }
                const std::ptrdiff_t shift = dz * step_z + dy * step_y + dx;

                // Squared differences. Position (0, 0, 0) of `squares` is the
                // block's first voxel minus the patch radius, which lies at
                // (sx, sy, sz) of `input`.
                for (std::size_t z = 0; z < ed; ++z) {
                    for (std::size_t y = 0; y < eh; ++y) {
                        const float* a = buf.input.data() + (z + sz) * iw * ih + (y + sy) * iw + sx;
                        const float* b = a + shift;
                        Real* out = buf.squares.data() + (z * eh + y) * ew;
                        for (std::size_t x = 0; x < ew; ++x) {
                            const Real difference = static_cast<Real>(a[x]) - static_cast<Real>(b[x]);
                            out[x] = difference * difference;
                        }
                    }
                }

                // Box sums over the patch, one axis at a time.
                const Real* current = buf.squares.data();
                if (px > 0) {
                    sumAlongX(current, eh * ed, bw, px, buf.along_x.data());
                    current = buf.along_x.data();
                }
                if (py > 0) {
                    sumAcrossRows(current, static_cast<std::size_t>(bw), bh, py, ed, buf.along_y.data());
                    current = buf.along_y.data();
                }
                if (pz > 0) {
                    sumAcrossRows(current, static_cast<std::size_t>(bw) * bh, bd, pz, 1,
                                  buf.distance.data());
                    current = buf.distance.data();
                }

                // Weights. `current` now holds bw * bh * bd distances.
                for (int z = 0; z < bd; ++z) {
                    for (int y = 0; y < bh; ++y) {
                        const std::size_t first = (static_cast<std::size_t>(z) * bh + y) * bw;
                        const float* other = buf.input.data() + (z + mz) * iw * ih + (y + my) * iw + mx + shift;
                        for (int x = 0; x < bw; ++x) {
                            const std::size_t i = first + static_cast<std::size_t>(x);
                            const double excess = static_cast<double>(current[i]) * scale - offset;
                            const double weight = excess <= 0.0 ? 1.0
                                                  : rational    ? h2 / (h2 + excess)
                                                                : std::exp(-excess / h2);
                            buf.weighted_sum[i] += weight * other[x];
                            buf.weight_sum[i] += weight;
                            buf.max_weight[i] = std::max(buf.max_weight[i], weight);
                        }
                    }
                }
            }
        }
    }

    // The voxel itself, at the largest weight found for it.
    for (int z = 0; z < bd; ++z) {
        for (int y = 0; y < bh; ++y) {
            const float* centre = buf.input.data() + (z + mz) * iw * ih + (y + my) * iw + mx;
            float* out = result.data() +
                         (static_cast<std::size_t>(block.z0 + z) * size.height + block.y0 + y) *
                             size.width + block.x0;
            for (int x = 0; x < bw; ++x) {
                const std::size_t i = (static_cast<std::size_t>(z) * bh + y) * bw + x;
                double max_weight = buf.max_weight[i];
                if (max_weight < 0.0) {
                    max_weight = 1.0;  // search window of one voxel
                }
                const double weighted = buf.weighted_sum[i] + max_weight * centre[x];
                const double total = buf.weight_sum[i] + max_weight;
                out[x] = static_cast<float>(weighted / total);
            }
        }
    }
}

template <typename Real>
std::vector<float> filterAll(const std::vector<float>& volume, const VolumeSize& size,
                             const NlmParams& params, int block_size) {
    std::vector<Region> blocks;
    for (int z = 0; z < size.depth; z += block_size) {
        for (int y = 0; y < size.height; y += block_size) {
            for (int x = 0; x < size.width; x += block_size) {
                blocks.push_back(Region{x, y, z, std::min(block_size, size.width - x),
                                        std::min(block_size, size.height - y),
                                        std::min(block_size, size.depth - z)});
            }
        }
    }

    std::vector<float> result(volume.size());
    #pragma omp parallel
    {
        BlockBuffers<Real> buffers;
        #pragma omp for schedule(dynamic)
        for (std::ptrdiff_t b = 0; b < static_cast<std::ptrdiff_t>(blocks.size()); ++b) {
            filterBlock(volume, size, params, blocks[static_cast<std::size_t>(b)], buffers, result);
        }
    }
    return result;
}

}  // namespace

std::vector<float> nlmFast(const std::vector<float>& volume, const VolumeSize& size,
                           const NlmParams& params, const NlmFastOptions& options) {
    checkNlmInput(volume, size, params);
    if (options.block_size < 1) {
        throw std::invalid_argument("the block size must be at least 1");
    }
    if (options.single_precision) {
        return filterAll<float>(volume, size, params, options.block_size);
    }
    return filterAll<double>(volume, size, params, options.block_size);
}

}  // namespace denoiser
