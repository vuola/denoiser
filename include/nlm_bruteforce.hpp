// denoiser - 3D Non-Local Means filter for CT volumes.
// Derived from the C++ filter by Samuel Väisänen (2024), itself a port of the
// ImageJ plugin by Pascal Behnel and Thorsten Wagner. GPL-3.0, see LICENSE.

/// @file
/// The Non-Local Means filter written straight from its definition.
///
/// This version is slow (every patch distance is summed voxel by voxel) and
/// exists as ground truth: fast versions are checked against it on small
/// volumes. With the settings of originalParams() it computes what the
/// original program computes, with these differences in arithmetic only: all
/// sums are in double precision and the patch distance is summed directly
/// instead of through an integral image.
///
/// The settings (NlmParams) also hold the choices that the original program
/// has built in: one radius for all three axes, a patch distance that is not
/// divided by the patch size, and a rational weight function. Each can be
/// set differently here (README chapter 12).

#ifndef DENOISER_NLM_BRUTEFORCE_HPP
#define DENOISER_NLM_BRUTEFORCE_HPP

#include <vector>

#include "volume.hpp"

namespace denoiser {

/// How a patch distance above the noise level is turned into a weight.
/// `excess` is defined at nlmBruteForce().
enum class NlmKernel {
    kRational,     ///< weight = h^2 / (h^2 + excess); the original program and the ImageJ plugin
    kExponential,  ///< weight = exp(-excess / h^2); the IPOL paper (Buades, Coll, Morel 2011)
};

/// Settings of the filter. Radii are in voxels and given per axis in the
/// order x, y, z; radius r means 2r + 1 voxels along that axis.
struct NlmParams {
    double sigma = 0.0;              ///< noise standard deviation, in the unit of the voxel values
    double h = 0.0;                  ///< filtering strength, in the unit of the voxel values
    int patch_radius[3] = {1, 1, 1};     ///< half size of the compared patches
    int search_radius[3] = {10, 10, 10}; ///< half size of the search window
    /// false: the patch distance is the sum over the patch (original program).
    /// true: it is the mean over the patch (IPOL paper), so that h means the
    /// same for every patch size.
    bool normalise_distance = false;
    NlmKernel kernel = NlmKernel::kRational;  ///< weight function
};

/// Returns the settings the original program chooses for a noise level.
///
/// The original program expects voxel values on the scale 0...255 and picks
/// patch radius, search radius and h from sigma by this table:
///
/// | sigma      | patch radius | search radius | h            |
/// |------------|--------------|---------------|--------------|
/// | up to 15   | 1            | 10            | 0.40 * sigma |
/// | up to 30   | 2            | 10            | 0.40 * sigma |
/// | up to 45   | 3            | 17            | 0.35 * sigma |
/// | up to 75   | 4            | 17            | 0.35 * sigma |
/// | above 75   | 5            | 17            | 0.30 * sigma |
///
/// @param sigma  noise standard deviation on the 0...255 scale, > 0
/// @return settings with the same radius on all three axes, the distance
///         not normalised and the rational kernel
/// @throws std::invalid_argument if sigma is not a positive finite number
NlmParams originalParams(double sigma);

/// Returns the number of voxels in a patch: the product of 2r + 1 over the
/// three axes.
/// @param params  filter settings with radii >= 0
double patchVoxels(const NlmParams& params);

/// Returns the weight of one shifted position.
///
/// @param distance  sum over the patch of the squared differences, >= 0
/// @param params    filter settings; sigma > 0, h > 0, radii >= 0
/// @return a weight between 0 and 1; see nlmBruteForce() for the formula
double nlmWeight(double distance, const NlmParams& params);

/// Checks the arguments of the filters (nlmBruteForce(), nlmFast()).
/// @throws std::invalid_argument if a dimension is below 1, the vector
///         length does not match the size, sigma or h is not positive, or a
///         radius is negative
void checkNlmInput(const std::vector<float>& volume, const VolumeSize& size,
                   const NlmParams& params);

/// Filters a volume with Non-Local Means, computed from the definition.
///
/// For every voxel p and every shift d inside the search window (d != 0):
///
///     distance(p, d) = sum over the patch offsets q of
///                      (v(p + q) - v(p + q + d))^2
///     N              = number of voxels in a patch
///     excess         = distance - 2 * sigma^2 * N        (normalise_distance false)
///                      distance / N - 2 * sigma^2        (normalise_distance true)
///     weight(p, d)   = 1                                 if excess <= 0
///                      h^2 / (h^2 + excess)              otherwise, kRational
///                      exp(-excess / h^2)                otherwise, kExponential
///
/// The result at p is the weighted mean of the voxels v(p + d), with the voxel
/// p itself included at the largest weight found for p. Coordinates outside
/// the volume are moved to the nearest voxel inside it (edge replication).
/// The defaults (distance not divided by the patch size, rational kernel)
/// are the behaviour of the original program. Without normalisation the
/// same h acts sqrt(N) times more sharply than with it.
///
/// Voxels are computed independently of each other and each in a fixed order
/// of shifts (z outermost, x innermost), so the result does not depend on the
/// number of threads.
///
/// Cost: voxels * search window size * patch size operations. Meant for
/// volumes up to about 64^3.
///
/// @param volume  voxel values, x fastest; size must be width * height * depth
/// @param size    volume size, every dimension at least 1
/// @param params  filter settings; sigma > 0, h > 0, radii >= 0
/// @return the filtered volume, same size and order as the input
/// @throws std::invalid_argument if a size, radius or setting is out of range
///         or the vector length does not match the size
std::vector<float> nlmBruteForce(const std::vector<float>& volume, const VolumeSize& size,
                                 const NlmParams& params);

}  // namespace denoiser

#endif
