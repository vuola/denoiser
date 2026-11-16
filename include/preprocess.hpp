// denoiser - 3D Non-Local Means filter for CT volumes.
// Derived from the C++ filter by Samuel Väisänen (2024), itself a port of the
// ImageJ plugin by Pascal Behnel and Thorsten Wagner. GPL-3.0, see LICENSE.

/// @file
/// Pre-processing: what happens to voxel values and to the volume between
/// reading a scan and filtering it, and the way back afterwards.
///
///   - IntensityMap: the linear mapping between the stored values (for
///     example 16-bit counts) and the filter's working range 0...255.
///   - centredCube(): a cube-shaped region in the middle of a volume.
///   - replaceNonFinite(): NaN and infinite values, which the filter cannot
///     work with.
///   - cropScan(): a region of a scan written as a scan of its own.
///
/// The mapping is linear and never clips (README requirement R2): a value
/// outside the mapped range simply lands outside 0...255, and comes back
/// unchanged on the way out.

#ifndef DENOISER_PREPROCESS_HPP
#define DENOISER_PREPROCESS_HPP

#include <cstdint>
#include <string>
#include <vector>

#include "scan_reader.hpp"
#include "volume.hpp"

namespace denoiser {

/// Linear mapping between stored values and the working range.
///
///     working = (stored - low) * 255 / (high - low)
///     stored  = working * (high - low) / 255 + low
///
/// The stored value `low` becomes 0 and `high` becomes 255. Nothing is
/// clipped. The filter chooses its settings from the noise sigma on the
/// working scale, so the mapping decides what a given sigma means.
struct IntensityMap {
    double low = 0.0;       ///< stored value that becomes 0
    double high = 65535.0;  ///< stored value that becomes 255

    /// Maps one stored value to the working range. The arithmetic is in
    /// double precision; the result is rounded to float once.
    float toWork(float stored) const {
        return static_cast<float>((static_cast<double>(stored) - low) * 255.0 / (high - low));
    }

    /// Maps one working value back to the scale of the stored values.
    float toStored(float work) const {
        return static_cast<float>(static_cast<double>(work) * (high - low) / 255.0 + low);
    }
};

/// Returns the default mapping for a sample type: the full range of the
/// type for integers (0...255 for 8 bits, 0...65535 for 16 bits), which is
/// how the samples of the old workflow were scaled (README section 2.3).
/// Float scans are taken to be on the working scale already (0...255).
/// @param type  the sample type of the scan
IntensityMap defaultMap(SampleType type);

/// Checks a mapping.
/// @param map  the mapping
/// @throws std::invalid_argument unless low and high are finite and high > low
void checkMap(const IntensityMap& map);

/// Maps all values from the stored scale to the working range, in place.
/// @param values  values on the stored scale; replaced by working values
/// @param map     the mapping
/// @throws std::invalid_argument if the mapping is not valid
void mapToWork(std::vector<float>& values, const IntensityMap& map);

/// Maps all values from the working range back to the stored scale, in
/// place. The values are not rounded or limited; the writer does that and
/// reports it (scan_writer.hpp).
/// @param values  working values; replaced by values on the stored scale
/// @param map     the mapping
/// @throws std::invalid_argument if the mapping is not valid
void mapToStored(std::vector<float>& values, const IntensityMap& map);

/// Returns the cube of edge length `edge` in the middle of a volume. Where
/// the space left over is odd, the cube starts half a voxel towards the
/// origin: start = (dimension - edge) / 2, rounded down.
/// @param size  the volume size
/// @param edge  edge length of the cube in voxels, >= 1
/// @return the region
/// @throws std::invalid_argument if the cube does not fit the volume
Region centredCube(const VolumeSize& size, int edge);

/// Replaces NaN and infinite values, in place.
/// @param values       the values to check
/// @param replacement  the value written in place of each NaN or infinity
/// @return how many values were replaced
std::int64_t replaceNonFinite(std::vector<float>& values, float replacement);

/// Writes a region of a scan as a scan of its own: voxel file and .vgi
/// file, in the same sample type, byte order and voxel size as the source.
/// The values are copied unchanged. The source is read slab by slab, so a
/// large region does not need its size in memory.
///
/// @param source       the scan to cut from
/// @param region       the region to cut; it must lie inside the volume
/// @param volume_path  the voxel file to write; its folder must exist; the
///                     .vgi file is written next to it
/// @return the description of the new scan
/// @throws std::runtime_error if the region does not fit, the output is the
///         source file itself, or a file cannot be read or written
ScanInfo cropScan(const ScanInfo& source, const Region& region, const std::string& volume_path);

}  // namespace denoiser

#endif  // DENOISER_PREPROCESS_HPP
