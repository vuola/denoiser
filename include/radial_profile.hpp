// denoiser - 3D Non-Local Means filter for CT volumes.
// Derived from the C++ filter by Samuel Väisänen (2024), itself a port of the
// ImageJ plugin by Pascal Behnel and Thorsten Wagner. GPL-3.0, see LICENSE.

/// @file
/// Brightness against distance from an axis parallel to z.
///
/// The rods are scanned with their axis along z, so "distance from the
/// axis" is the distance in the x-y plane from a centre point. The profile
/// shows whether the scan is equally bright at the centre and at the rim
/// (README requirement R4), where a bore hole or the air around the rod
/// lies, and where voxels sit at the limits of the value range.
///
/// The voxels are sorted into rings by their distance from the centre:
/// ring k holds the distances from k * ring_width up to (k + 1) * ring_width
/// voxels. For every ring the mean and three percentiles are computed. In a
/// rod, the 10th percentile follows the darker material (epoxy) and the
/// 90th the brighter one (glass), so their change with distance shows drift
/// separately from changes in the mix of materials.

#ifndef DENOISER_RADIAL_PROFILE_HPP
#define DENOISER_RADIAL_PROFILE_HPP

#include <cstdint>
#include <vector>

#include "scan_reader.hpp"

namespace denoiser {

/// Settings of radialProfile().
struct ProfileSettings {
    double centre_x = 0.0;     ///< x of the axis, in voxel coordinates (voxel centres are whole numbers)
    double centre_y = 0.0;     ///< y of the axis
    double ring_width = 8.0;   ///< width of a ring in voxels, > 0
    int z0 = 0;                ///< first slice used
    int depth = 0;             ///< number of slices used; 0 = all from z0 to the end
    double value_low = 0.0;    ///< lower end of the value range, see RadialRing
    double value_high = 65535.0;  ///< upper end of the value range, > value_low
};

/// Result for one ring.
///
/// Percentiles are read from a histogram of 1024 equal classes between
/// value_low and value_high, with linear interpolation inside a class, so
/// they are accurate to about one class width ((value_high - value_low) /
/// 1024; 64 counts for the full 16-bit range). The mean is exact. NaN and
/// infinite values are left out of everything.
struct RadialRing {
    double r_inner = 0.0;       ///< inner radius in voxels (included)
    double r_outer = 0.0;       ///< outer radius in voxels (not included)
    std::int64_t voxels = 0;    ///< finite voxels in the ring
    double mean = 0.0;          ///< mean value
    double p10 = 0.0;           ///< 10th percentile
    double p50 = 0.0;           ///< median
    double p90 = 0.0;           ///< 90th percentile
    std::int64_t at_low = 0;    ///< voxels with value <= value_low
    std::int64_t at_high = 0;   ///< voxels with value >= value_high
};

/// Returns settings for a scan: the axis through the middle of the x-y
/// plane, all slices, rings of 8 voxels, and the value range of the sample
/// type (0...255 for 8 bits, 0...65535 for 16 bits and for float).
/// @param info  the scan description
ProfileSettings defaultProfileSettings(const ScanInfo& info);

/// Computes the brightness profile of a scan around an axis parallel to z.
///
/// The scan is read slice by slice. Rings run from the axis to the corner
/// of the x-y plane that is farthest from it; a ring without voxels has
/// voxels == 0 and zeros elsewhere.
///
/// @param info      the scan description
/// @param settings  axis, ring width, slices and value range
/// @return one entry per ring, innermost first
/// @throws std::invalid_argument if a setting is out of range
/// @throws std::runtime_error if the scan cannot be read
std::vector<RadialRing> radialProfile(const ScanInfo& info, const ProfileSettings& settings);

}  // namespace denoiser

#endif  // DENOISER_RADIAL_PROFILE_HPP
