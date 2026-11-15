// denoiser - 3D Non-Local Means filter for CT volumes.
// Derived from the C++ filter by Samuel Väisänen (2024), itself a port of the
// ImageJ plugin by Pascal Behnel and Thorsten Wagner. GPL-3.0, see LICENSE.

/// @file
/// Writing a volume in the CT lab's format: a voxel file (.vol, no header,
/// x fastest) and a text file (.vgi) that describes it, so that the result
/// opens wherever the input does.
///
/// A scan is written in two steps, which lets a large volume be written
/// block by block:
///
///     createScan(info);                      // both files; all voxels zero
///     writeRegion(info, region, values);     // as often as needed
///
/// writeVolume() does both for a volume held in memory.
///
/// The description is the same ScanInfo that the reader returns
/// (scan_reader.hpp). To write a result "in the input's sample type", copy
/// the input's ScanInfo and change volume_path (and the size, if a region
/// was cut).
///
/// Values are given as floats and stored in the scan's sample type:
///   - 8- and 16-bit integers: rounded to the nearest whole number (halves
///     away from zero); values below 0 or above the largest storable number
///     are stored as that limit, NaN as 0. How often this happened is
///     reported, so that it never goes unnoticed.
///   - 32-bit floats: stored unchanged, NaN and infinity included.
///
/// All functions report problems by throwing std::runtime_error with a
/// message that names the file.

#ifndef DENOISER_SCAN_WRITER_HPP
#define DENOISER_SCAN_WRITER_HPP

#include <cstdint>
#include <string>
#include <vector>

#include "scan_reader.hpp"
#include "volume.hpp"

namespace denoiser {

/// What happened to the values of one write call on their way into the
/// file's sample type. All counts are zero for float scans.
struct WriteReport {
    std::int64_t below = 0;       ///< values below 0, stored as 0
    std::int64_t above = 0;       ///< values above the largest storable number, stored as that number
    std::int64_t non_finite = 0;  ///< NaN values, stored as 0 (infinities count as below or above)

    /// Number of values that could not be stored as they were.
    std::int64_t changed() const { return below + above + non_finite; }
};

/// Returns the path of the .vgi file that belongs to a voxel file: the same
/// name with the extension .vgi.
/// @param volume_path  path of the voxel file
std::string vgiPathFor(const std::string& volume_path);

/// Returns the text of a .vgi file for a scan.
///
/// The text holds the volume size, the sample type, the voxel file's name
/// without any folder, and the voxel size if it is known. Lines end with a
/// line feed. For 16-bit scans the text is the same as that of the tracked
/// test scans. No path is written, so nothing about the computer that made
/// the file gets into it.
///
/// @param info  the scan description; volume_path gives the file name
/// @return the text
/// @throws std::runtime_error if a dimension is below 1
std::string vgiText(const ScanInfo& info);

/// Creates the two files of a scan: the .vgi file, and the voxel file at
/// its full size with every voxel zero. Existing files of the same names
/// are replaced. The folder must exist.
///
/// @param info  the scan description. volume_path is the voxel file to
///              create; header_bytes must be 0 (files are written without
///              a header)
/// @throws std::runtime_error if the description is out of range or a file
///         cannot be created
void createScan(const ScanInfo& info);

/// Writes a box of voxels into a scan made with createScan().
///
/// Voxels outside the box keep their content, so a volume can be written
/// in blocks in any order.
///
/// @param info    the scan description used with createScan()
/// @param region  the box to write; it must lie inside the volume
/// @param values  region.width * region.height * region.depth values,
///                x fastest
/// @return how many values had to be changed to fit the sample type
/// @throws std::runtime_error if the region does not fit the volume, the
///         number of values does not match the region, or the file cannot
///         be written
WriteReport writeRegion(const ScanInfo& info, const Region& region,
                        const std::vector<float>& values);

/// Writes a whole volume held in memory: createScan() followed by one
/// writeRegion() for all voxels.
///
/// @param info    the scan description, as for createScan()
/// @param values  width * height * depth values, x fastest
/// @return how many values had to be changed to fit the sample type
/// @throws std::runtime_error as createScan() and writeRegion()
WriteReport writeVolume(const ScanInfo& info, const std::vector<float>& values);

}  // namespace denoiser

#endif  // DENOISER_SCAN_WRITER_HPP
