// denoiser - 3D Non-Local Means filter for CT volumes.
// Derived from the C++ filter by Samuel Väisänen (2024), itself a port of the
// ImageJ plugin by Pascal Behnel and Thorsten Wagner. GPL-3.0, see LICENSE.

/// @file
/// The Non-Local Means filter, computed shift by shift.
///
/// This computes the filter defined at nlmBruteForce() with the same
/// settings (NlmParams), but without summing every patch distance voxel by
/// voxel. For one shift d of the search window at a time, the squared
/// differences (v(u) - v(u + d))^2 are formed once for all positions u, and
/// the patch distance of every voxel is their sum over the patch box. That
/// box sum is taken one axis after the other with running sums: along x,
/// then y, then z, each costing two operations per voxel whatever the patch
/// radius. An axis with patch radius 0 costs nothing.
///
/// Unlike the original program, which uses one summed-area table per block
/// of 256^3 voxels and loses precision in the differences of its large
/// entries, the running sums here never span more than one block edge plus
/// the patch. In single precision (the default) the result differs from
/// nlmBruteForce() by at most one count of a 16-bit output; in double
/// precision it was identical in every case measured (README chapter 13).
///
/// The output is computed in cubic blocks, each independently and with the
/// shifts in the fixed order of nlmBruteForce(), so the result does not
/// depend on the number of threads. It does depend, at the level of
/// rounding, on the block size, because the running sums start anew at each
/// block.

#ifndef DENOISER_NLM_FAST_HPP
#define DENOISER_NLM_FAST_HPP

#include <vector>

#include "nlm_bruteforce.hpp"
#include "volume.hpp"

namespace denoiser {

/// Settings of nlmFast() that change the arithmetic, not the filter.
struct NlmFastOptions {
    /// Edge of the cubic blocks the output is computed in, >= 1. Blocks at
    /// the end of an axis are smaller.
    int block_size = 64;
    /// true (default): squared differences and their running sums in single
    /// precision (float), as a GPU version would use them; faster, and off
    /// from nlmBruteForce() by far less than the noise.
    /// false: in double precision, which reproduces nlmBruteForce() to
    /// rounding; used by the checks against it.
    bool single_precision = true;
};

/// Filters a volume with Non-Local Means.
///
/// The result agrees with nlmBruteForce() for the same settings to within
/// rounding (README chapter 13 gives the measured difference). Memory use
/// besides input and output is about 16 * (block + 2 * (patch + search))^3
/// bytes per thread.
///
/// @param volume   voxel values, x fastest; size must be width * height * depth
/// @param size     volume size, every dimension at least 1
/// @param params   filter settings; sigma > 0, h > 0, radii >= 0
/// @param options  block size and precision
/// @return the filtered volume, same size and order as the input
/// @throws std::invalid_argument if a size, radius, setting or option is
///         out of range or the vector length does not match the size
std::vector<float> nlmFast(const std::vector<float>& volume, const VolumeSize& size,
                           const NlmParams& params, const NlmFastOptions& options = {});

}  // namespace denoiser

#endif  // DENOISER_NLM_FAST_HPP
