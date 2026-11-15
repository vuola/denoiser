// denoiser - 3D Non-Local Means filter for CT volumes.
// Derived from the C++ filter by Samuel Väisänen (2024), itself a port of the
// ImageJ plugin by Pascal Behnel and Thorsten Wagner. GPL-3.0, see LICENSE.

/// @file
/// Basic descriptions of a volume that every part of the program shares:
/// its size and a box-shaped region inside it.

#ifndef DENOISER_VOLUME_HPP
#define DENOISER_VOLUME_HPP

#include <cstdint>

namespace denoiser {

/// Size of a volume in voxels. Voxel (x, y, z) is stored at index
/// (z * height + y) * width + x, that is, x runs fastest.
struct VolumeSize {
    int width = 0;   ///< number of voxels along x
    int height = 0;  ///< number of voxels along y
    int depth = 0;   ///< number of voxels along z
};

/// A box of voxels inside a volume. Coordinates are zero-based: the box
/// holds the voxels x0 ... x0 + width - 1 along x, and likewise along y and z.
struct Region {
    int x0 = 0;      ///< first voxel along x
    int y0 = 0;      ///< first voxel along y
    int z0 = 0;      ///< first voxel along z
    int width = 0;   ///< number of voxels along x
    int height = 0;  ///< number of voxels along y
    int depth = 0;   ///< number of voxels along z
};

/// Returns the number of voxels of a volume as a 64-bit number, so that
/// volumes of more than 2^31 voxels are counted correctly.
/// @param size  volume size; dimensions are expected to be >= 0
inline std::int64_t voxelCount(const VolumeSize& size) {
    return static_cast<std::int64_t>(size.width) * size.height * size.depth;
}

/// Returns the region that covers a whole volume.
/// @param size  volume size
inline Region wholeVolume(const VolumeSize& size) {
    return Region{0, 0, 0, size.width, size.height, size.depth};
}

/// Tells whether a region is non-empty and lies completely inside a volume.
/// @param region  the region to check
/// @param size    the volume size
inline bool regionFits(const Region& region, const VolumeSize& size) {
    return region.x0 >= 0 && region.y0 >= 0 && region.z0 >= 0 &&
           region.width >= 1 && region.height >= 1 && region.depth >= 1 &&
           region.width <= size.width - region.x0 &&
           region.height <= size.height - region.y0 &&
           region.depth <= size.depth - region.z0;
}

}  // namespace denoiser

#endif  // DENOISER_VOLUME_HPP
