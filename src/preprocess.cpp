// denoiser - 3D Non-Local Means filter for CT volumes.
// Derived from the C++ filter by Samuel Väisänen (2024), itself a port of the
// ImageJ plugin by Pascal Behnel and Thorsten Wagner. GPL-3.0, see LICENSE.

#include "preprocess.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <stdexcept>

#include "scan_writer.hpp"

namespace denoiser {

IntensityMap defaultMap(SampleType type) {
    IntensityMap map;
    map.low = 0.0;
    map.high = type == SampleType::kUInt16 ? 65535.0 : 255.0;
    return map;
}

void checkMap(const IntensityMap& map) {
    if (!std::isfinite(map.low) || !std::isfinite(map.high) || !(map.high > map.low)) {
        throw std::invalid_argument("the mapped range needs two finite values, the second larger than the first");
    }
}

void mapToWork(std::vector<float>& values, const IntensityMap& map) {
    checkMap(map);
    for (float& value : values) {
        value = map.toWork(value);
    }
}

void mapToStored(std::vector<float>& values, const IntensityMap& map) {
    checkMap(map);
    for (float& value : values) {
        value = map.toStored(value);
    }
}

Region centredCube(const VolumeSize& size, int edge) {
    if (edge < 1 || edge > size.width || edge > size.height || edge > size.depth) {
        throw std::invalid_argument("a centred cube of " + std::to_string(edge) +
                                    " voxels does not fit a volume of " + std::to_string(size.width) +
                                    " x " + std::to_string(size.height) + " x " +
                                    std::to_string(size.depth));
    }
    return Region{(size.width - edge) / 2, (size.height - edge) / 2, (size.depth - edge) / 2,
                  edge, edge, edge};
}

std::int64_t replaceNonFinite(std::vector<float>& values, float replacement) {
    std::int64_t replaced = 0;
    for (float& value : values) {
        if (!std::isfinite(value)) {
            value = replacement;
            ++replaced;
        }
    }
    return replaced;
}

ScanInfo cropScan(const ScanInfo& source, const Region& region, const std::string& volume_path) {
    if (!regionFits(region, source.size)) {
        throw std::runtime_error(source.volume_path + ": the region to cut does not fit the volume");
    }
    std::error_code error;
    if (std::filesystem::equivalent(source.volume_path, volume_path, error)) {
        throw std::runtime_error(volume_path + ": the output would overwrite the source scan");
    }

    ScanInfo result = source;
    result.volume_path = volume_path;
    result.size = VolumeSize{region.width, region.height, region.depth};
    result.header_bytes = 0;
    createScan(result);

    // Slabs of about 16 million voxels (64 MB as floats) at a time.
    const std::int64_t slice = static_cast<std::int64_t>(region.width) * region.height;
    const int slab_depth = static_cast<int>(std::max<std::int64_t>(1, 16000000 / slice));
    for (int z = 0; z < region.depth; z += slab_depth) {
        const int depth = std::min(slab_depth, region.depth - z);
        const Region from{region.x0, region.y0, region.z0 + z, region.width, region.height, depth};
        const Region to{0, 0, z, region.width, region.height, depth};
        writeRegion(result, to, readRegion(source, from));
    }
    return result;
}

}  // namespace denoiser
