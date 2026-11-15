// denoiser - 3D Non-Local Means filter for CT volumes.
// Derived from the C++ filter by Samuel Väisänen (2024), itself a port of the
// ImageJ plugin by Pascal Behnel and Thorsten Wagner. GPL-3.0, see LICENSE.

#include "scan_writer.hpp"

#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <stdexcept>

namespace denoiser {

namespace {

namespace fs = std::filesystem;

/// Checks the parts of a description that the writer depends on.
void checkDescription(const ScanInfo& info) {
    if (info.volume_path.empty()) {
        throw std::runtime_error("the voxel file's path is empty");
    }
    if (info.size.width < 1 || info.size.height < 1 || info.size.depth < 1) {
        throw std::runtime_error(info.volume_path + ": every volume dimension must be at least 1");
    }
    if (info.header_bytes != 0) {
        throw std::runtime_error(info.volume_path + ": files are written without a header");
    }
}

/// Converts one value to a whole number in 0 ... limit and counts in
/// `report` if it had to be changed.
unsigned toInteger(float value, float limit, WriteReport& report) {
    if (std::isnan(value)) {
        ++report.non_finite;
        return 0u;
    }
    const float rounded = std::round(value);
    if (rounded < 0.0f) {
        ++report.below;
        return 0u;
    }
    if (rounded > limit) {
        ++report.above;
        return static_cast<unsigned>(limit);
    }
    return static_cast<unsigned>(rounded);
}

/// Converts `count` values to file bytes of the scan's sample type.
void encode(const float* values, std::size_t count, const ScanInfo& info, unsigned char* bytes,
            WriteReport& report) {
    const bool little = info.byte_order == ByteOrder::kLittleEndian;
    switch (info.sample_type) {
        case SampleType::kUInt8:
            for (std::size_t i = 0; i < count; ++i) {
                bytes[i] = static_cast<unsigned char>(toInteger(values[i], 255.0f, report));
            }
            break;
        case SampleType::kUInt16:
            for (std::size_t i = 0; i < count; ++i) {
                const unsigned v = toInteger(values[i], 65535.0f, report);
                const unsigned char low = static_cast<unsigned char>(v & 0xFFu);
                const unsigned char high = static_cast<unsigned char>(v >> 8);
                bytes[2 * i] = little ? low : high;
                bytes[2 * i + 1] = little ? high : low;
            }
            break;
        case SampleType::kFloat32:
            for (std::size_t i = 0; i < count; ++i) {
                std::uint32_t bits = 0;
                std::memcpy(&bits, &values[i], sizeof(bits));
                for (int b = 0; b < 4; ++b) {
                    const int shift = little ? 8 * b : 8 * (3 - b);
                    bytes[4 * i + static_cast<std::size_t>(b)] =
                        static_cast<unsigned char>((bits >> shift) & 0xFFu);
                }
            }
            break;
    }
}

/// Formats a voxel edge length with six decimals, as the lab's files do.
std::string sixDecimals(double value) {
    char text[64];
    std::snprintf(text, sizeof(text), "%.6f", value);
    return text;
}

}  // namespace

std::string vgiPathFor(const std::string& volume_path) {
    fs::path path(volume_path);
    path.replace_extension(".vgi");
    return path.generic_string();
}

std::string vgiText(const ScanInfo& info) {
    if (info.size.width < 1 || info.size.height < 1 || info.size.depth < 1) {
        throw std::runtime_error(info.volume_path + ": every volume dimension must be at least 1");
    }
    const std::string size = std::to_string(info.size.width) + " " +
                             std::to_string(info.size.height) + " " +
                             std::to_string(info.size.depth);
    const int bits = 8 * bytesPerSample(info.sample_type);
    std::string type = "bitsperelement = " + std::to_string(bits) + "\n";
    if (info.sample_type == SampleType::kFloat32) {
        type += "datatype = float\n";
    } else {
        type += "datatype = unsigned integer\n";
        type += "datarange = 0 " + std::to_string((1 << bits) - 1) + "\n";
    }

    std::string text = "{volume1}\n[representation]\nsize = " + size + "\n" + type +
                       "[file1]\nFileFormat = raw\nSize = " + size + "\nName = " +
                       fs::path(info.volume_path).filename().string() + "\n" + type;
    if (info.voxel_size[0] > 0.0 && info.voxel_size[1] > 0.0 && info.voxel_size[2] > 0.0) {
        text += "{volumeprimitive9}\n[geometry]\nresolution = " + sixDecimals(info.voxel_size[0]) +
                " " + sixDecimals(info.voxel_size[1]) + " " + sixDecimals(info.voxel_size[2]) +
                "\nunit = " + (info.unit.empty() ? std::string("mm") : info.unit) + "\n";
    }
    return text;
}

void createScan(const ScanInfo& info) {
    checkDescription(info);

    const std::string vgi_path = vgiPathFor(info.volume_path);
    if (vgi_path == fs::path(info.volume_path).generic_string()) {
        throw std::runtime_error(info.volume_path + ": the voxel file must not have the extension .vgi");
    }
    {
        std::ofstream vgi(vgi_path, std::ios::binary | std::ios::trunc);
        vgi << vgiText(info);
        if (!vgi) {
            throw std::runtime_error(vgi_path + ": cannot write the .vgi file");
        }
    }
    {
        std::ofstream volume(info.volume_path, std::ios::binary | std::ios::trunc);
        if (!volume) {
            throw std::runtime_error(info.volume_path + ": cannot create the voxel file");
        }
    }
    // Growing an empty file fills it with zero bytes without writing them
    // one by one.
    std::error_code error;
    fs::resize_file(info.volume_path, static_cast<std::uintmax_t>(expectedFileSize(info)), error);
    if (error) {
        throw std::runtime_error(info.volume_path + ": cannot set the file size (" +
                                 error.message() + ")");
    }
}

WriteReport writeRegion(const ScanInfo& info, const Region& region,
                        const std::vector<float>& values) {
    checkDescription(info);
    if (!regionFits(region, info.size)) {
        throw std::runtime_error(info.volume_path + ": the region to write does not fit the volume");
    }
    const std::size_t row = static_cast<std::size_t>(region.width);
    const std::size_t rows = static_cast<std::size_t>(region.height);
    if (values.size() != row * rows * static_cast<std::size_t>(region.depth)) {
        throw std::runtime_error(info.volume_path +
                                 ": the number of values does not match the region to write");
    }
    // The file must be the one createScan() made: a wrong size would mean
    // writing into some other file.
    checkScan(info);

    std::fstream file(info.volume_path, std::ios::binary | std::ios::in | std::ios::out);
    if (!file) {
        throw std::runtime_error(info.volume_path + ": cannot open the voxel file for writing");
    }

    const int sample_bytes = bytesPerSample(info.sample_type);
    // Rows that span the full width follow each other in the file without a
    // gap, so a whole slice of the region is written at once.
    const bool full_width = region.x0 == 0 && region.width == info.size.width;
    const std::size_t chunk = full_width ? row * rows : row;
    const std::size_t chunks_per_slice = full_width ? 1 : rows;
    std::vector<unsigned char> buffer(chunk * static_cast<std::size_t>(sample_bytes));

    WriteReport report;
    const float* in = values.data();
    for (int z = region.z0; z < region.z0 + region.depth; ++z) {
        for (std::size_t c = 0; c < chunks_per_slice; ++c) {
            const std::int64_t first_voxel =
                (static_cast<std::int64_t>(z) * info.size.height + region.y0 +
                 static_cast<std::int64_t>(c)) * info.size.width + region.x0;
            encode(in, chunk, info, buffer.data(), report);
            file.seekp(static_cast<std::streamoff>(first_voxel * sample_bytes));
            file.write(reinterpret_cast<const char*>(buffer.data()),
                       static_cast<std::streamsize>(buffer.size()));
            in += chunk;
        }
    }
    file.flush();
    if (!file) {
        throw std::runtime_error(info.volume_path + ": writing voxels failed");
    }
    return report;
}

WriteReport writeVolume(const ScanInfo& info, const std::vector<float>& values) {
    checkDescription(info);
    if (static_cast<std::int64_t>(values.size()) != voxelCount(info.size)) {
        throw std::runtime_error(info.volume_path +
                                 ": the number of values does not match the volume size");
    }
    createScan(info);
    return writeRegion(info, wholeVolume(info.size), values);
}

}  // namespace denoiser
