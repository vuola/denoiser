// denoiser - 3D Non-Local Means filter for CT volumes.
// Derived from the C++ filter by Samuel Väisänen (2024), itself a port of the
// ImageJ plugin by Pascal Behnel and Thorsten Wagner. GPL-3.0, see LICENSE.

/// @file
/// Reading a CT scan as the lab delivers it.
///
/// A scan is a voxel file (.vol: samples without any header, x fastest, then
/// y, then z) and a text file (.vgi) that states the volume size, the sample
/// type and the voxel size. openScan() reads the .vgi file; describeRaw()
/// takes the same description by hand, for raw files that have no .vgi.
/// Either way the result is a ScanInfo, and readRegion() reads voxels with
/// it.
///
/// Voxel values are returned as 32-bit floats holding the stored numbers
/// unchanged (a 16-bit count of 31139 becomes 31139.0f; every 8- and 16-bit
/// value is exact in a float). Mapping to the working range is a separate
/// step.
///
/// All functions report problems by throwing std::runtime_error with a
/// message that names the file.

#ifndef DENOISER_SCAN_READER_HPP
#define DENOISER_SCAN_READER_HPP

#include <cstdint>
#include <string>
#include <vector>

#include "volume.hpp"

namespace denoiser {

/// How one voxel is stored in the file.
enum class SampleType {
    kUInt8,    ///< 8-bit unsigned integer
    kUInt16,   ///< 16-bit unsigned integer
    kFloat32,  ///< 32-bit IEEE float
};

/// Order of the bytes of one sample in the file.
enum class ByteOrder {
    kLittleEndian,  ///< lowest byte first (the CT lab's files)
    kBigEndian,     ///< highest byte first (files saved by ImageJ)
};

/// Everything needed to read voxels from a scan.
struct ScanInfo {
    std::string volume_path;                          ///< the voxel file
    VolumeSize size;                                  ///< volume size in voxels
    SampleType sample_type = SampleType::kUInt16;     ///< how a voxel is stored
    ByteOrder byte_order = ByteOrder::kLittleEndian;  ///< byte order of a sample
    std::int64_t header_bytes = 0;                    ///< bytes to skip before the first voxel
    double voxel_size[3] = {0.0, 0.0, 0.0};           ///< voxel edge lengths along x, y, z; 0 = not known
    std::string unit;                                 ///< unit of voxel_size, e.g. "mm"; empty = not known
};

/// Value statistics of a scan, in the unit of the stored numbers.
struct ScanStatistics {
    std::int64_t voxels = 0;      ///< number of voxels read
    std::int64_t non_finite = 0;  ///< NaN and infinite values (float scans only); left out of the other figures
    double minimum = 0.0;         ///< smallest value
    double maximum = 0.0;         ///< largest value
    double mean = 0.0;            ///< arithmetic mean
    double std_dev = 0.0;         ///< standard deviation (population formula)
};

/// Returns the number of bytes one sample takes in the file.
/// @param type  the sample type
int bytesPerSample(SampleType type);

/// Returns a short description of a sample type, e.g. "16-bit unsigned integer".
/// @param type  the sample type
std::string sampleTypeName(SampleType type);

/// Reads the volume description from the text of a .vgi file.
///
/// Keys are matched without regard to upper and lower case, and spaces
/// around keys and values are ignored. Where a key occurs several times the
/// first occurrence counts; in the lab's files that is the [representation]
/// block for size and type and the [geometry] block for the voxel size.
///
/// | Key              | Meaning                         | Needed |
/// |------------------|---------------------------------|--------|
/// | size             | width height depth in voxels    | yes    |
/// | bitsperelement   | 8, 16 or 32                     | yes    |
/// | datatype         | "unsigned integer" or "float"   | no, default "unsigned integer" |
/// | skipheader       | bytes before the first voxel    | no, default 0 |
/// | resolution       | voxel edge lengths x y z        | no     |
/// | unit             | unit of the resolution          | no     |
///
/// The .vgi format does not state the byte order; little-endian is assumed.
/// The voxel file's name is not taken from the text (volume_path stays
/// empty), because the lab's files hold a path of the lab's own computer
/// there.
///
/// @param text  the whole content of the .vgi file
/// @param name  name of the file, used in error messages only
/// @return the description, with volume_path empty
/// @throws std::runtime_error if a needed key is missing, a value cannot be
///         read, a dimension is below 1, or the sample type is not supported
ScanInfo parseVgi(const std::string& text, const std::string& name);

/// Opens a scan delivered by the CT lab.
///
/// `path` may be:
///   - a folder that holds exactly one .vgi file,
///   - a .vgi file,
///   - a voxel file (any other extension) with a .vgi file of the same name
///     next to it.
/// The voxel file is the .vgi file's name with the extension .vol, in the
/// same folder. Its size is checked with checkScan().
///
/// @param path  scan folder, .vgi file or voxel file
/// @return the description of the scan, ready for readRegion()
/// @throws std::runtime_error if a file is missing, the folder holds no or
///         several .vgi files, the .vgi file cannot be understood, or the
///         voxel file's size does not match the description
ScanInfo openScan(const std::string& path);

/// Describes a raw voxel file by hand, for files without a .vgi file.
///
/// @param volume_path   the voxel file
/// @param size          volume size in voxels, every dimension >= 1
/// @param sample_type   how a voxel is stored
/// @param byte_order    byte order of a sample
/// @param header_bytes  bytes to skip before the first voxel, >= 0
/// @return the description, checked with checkScan()
/// @throws std::runtime_error if a value is out of range, the file is
///         missing, or its size does not match the description
ScanInfo describeRaw(const std::string& volume_path, const VolumeSize& size,
                     SampleType sample_type, ByteOrder byte_order,
                     std::int64_t header_bytes = 0);

/// Returns the size in bytes the voxel file must have according to the
/// description: header plus voxels times bytes per sample.
/// @param info  the scan description
std::int64_t expectedFileSize(const ScanInfo& info);

/// Checks that the voxel file exists and has exactly the expected size.
/// @param info  the scan description
/// @throws std::runtime_error if the description is out of range, the file
///         cannot be opened, or its size differs from expectedFileSize()
void checkScan(const ScanInfo& info);

/// Reads a box of voxels from a scan.
///
/// Only the requested voxels are read from the file, so a small region of a
/// large scan is cheap.
///
/// @param info    the scan description, from openScan() or describeRaw()
/// @param region  the box to read; it must lie inside the volume
/// @return region.width * region.height * region.depth values, x fastest,
///         each holding the stored number unchanged
/// @throws std::runtime_error if the region does not fit the volume or the
///         file cannot be read
std::vector<float> readRegion(const ScanInfo& info, const Region& region);

/// Reads all voxels of a scan. Meant for small volumes: the result takes
/// four bytes per voxel.
/// @param info  the scan description
/// @return width * height * depth values, x fastest
/// @throws std::runtime_error if the file cannot be read
std::vector<float> readVolume(const ScanInfo& info);

/// Computes value statistics of a whole scan, reading it slice by slice so
/// that the memory needed is that of one slice.
/// @param info  the scan description
/// @return minimum, maximum, mean and standard deviation of the finite values
/// @throws std::runtime_error if the file cannot be read
ScanStatistics scanStatistics(const ScanInfo& info);

}  // namespace denoiser

#endif  // DENOISER_SCAN_READER_HPP
