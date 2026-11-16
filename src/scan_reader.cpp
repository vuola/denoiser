// denoiser - 3D Non-Local Means filter for CT volumes.
// Derived from the C++ filter by Samuel Väisänen (2024), itself a port of the
// ImageJ plugin by Pascal Behnel and Thorsten Wagner. GPL-3.0, see LICENSE.

#include "scan_reader.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>

namespace denoiser {

namespace {

namespace fs = std::filesystem;

/// Returns the text without spaces, tabs and line ends at both ends.
std::string trimmed(const std::string& text) {
    const char* blanks = " \t\r\n";
    const std::size_t first = text.find_first_not_of(blanks);
    if (first == std::string::npos) {
        return "";
    }
    return text.substr(first, text.find_last_not_of(blanks) - first + 1);
}

/// Returns the text in lower case (ASCII letters only).
std::string lowered(std::string text) {
    for (char& c : text) {
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    }
    return text;
}

/// Reads exactly `count` numbers from a text value such as "64 64 64".
/// Returns false if there are fewer or more, or something is not a number.
template <typename Number>
bool readNumbers(const std::string& value, int count, Number* numbers) {
    std::istringstream stream(value);
    for (int i = 0; i < count; ++i) {
        if (!(stream >> numbers[i])) {
            return false;
        }
    }
    std::string rest;
    return !(stream >> rest);
}

/// Returns the first value of every "key = value" line, keys in lower case.
std::map<std::string, std::string> firstValues(const std::string& text) {
    std::map<std::string, std::string> values;
    std::istringstream stream(text);
    std::string line;
    while (std::getline(stream, line)) {
        const std::size_t equals = line.find('=');
        if (equals == std::string::npos) {
            continue;
        }
        values.emplace(lowered(trimmed(line.substr(0, equals))), trimmed(line.substr(equals + 1)));
    }
    return values;
}

/// Checks the parts of a description that do not need the file.
void checkDescription(const ScanInfo& info) {
    if (info.size.width < 1 || info.size.height < 1 || info.size.depth < 1) {
        throw std::runtime_error(info.volume_path + ": every volume dimension must be at least 1");
    }
    if (info.header_bytes < 0) {
        throw std::runtime_error(info.volume_path + ": the header size must not be negative");
    }
}

/// Converts `count` samples from file bytes to floats.
void decode(const unsigned char* bytes, std::size_t count, const ScanInfo& info, float* out) {
    const bool little = info.byte_order == ByteOrder::kLittleEndian;
    switch (info.sample_type) {
        case SampleType::kUInt8:
            for (std::size_t i = 0; i < count; ++i) {
                out[i] = static_cast<float>(bytes[i]);
            }
            break;
        case SampleType::kUInt16:
            for (std::size_t i = 0; i < count; ++i) {
                const unsigned b0 = bytes[2 * i];
                const unsigned b1 = bytes[2 * i + 1];
                out[i] = static_cast<float>(little ? (b0 | (b1 << 8)) : (b1 | (b0 << 8)));
            }
            break;
        case SampleType::kFloat32:
            for (std::size_t i = 0; i < count; ++i) {
                const unsigned char* b = bytes + 4 * i;
                const std::uint32_t bits =
                    little ? (static_cast<std::uint32_t>(b[0]) | (static_cast<std::uint32_t>(b[1]) << 8) |
                              (static_cast<std::uint32_t>(b[2]) << 16) | (static_cast<std::uint32_t>(b[3]) << 24))
                           : (static_cast<std::uint32_t>(b[3]) | (static_cast<std::uint32_t>(b[2]) << 8) |
                              (static_cast<std::uint32_t>(b[1]) << 16) | (static_cast<std::uint32_t>(b[0]) << 24));
                std::memcpy(&out[i], &bits, sizeof(float));
            }
            break;
    }
}

/// Reads `count` samples that start at voxel index `first_voxel` of the file.
void readSamples(std::ifstream& file, const ScanInfo& info, std::int64_t first_voxel,
                 std::size_t count, std::vector<unsigned char>& buffer, float* out) {
    const int sample_bytes = bytesPerSample(info.sample_type);
    const std::size_t byte_count = count * static_cast<std::size_t>(sample_bytes);
    buffer.resize(byte_count);
    file.seekg(static_cast<std::streamoff>(info.header_bytes + first_voxel * sample_bytes));
    file.read(reinterpret_cast<char*>(buffer.data()), static_cast<std::streamsize>(byte_count));
    if (!file || file.gcount() != static_cast<std::streamsize>(byte_count)) {
        throw std::runtime_error(info.volume_path + ": reading voxels failed");
    }
    decode(buffer.data(), count, info, out);
}

/// Opens the voxel file for reading.
std::ifstream openVolumeFile(const ScanInfo& info) {
    std::ifstream file(info.volume_path, std::ios::binary);
    if (!file) {
        throw std::runtime_error(info.volume_path + ": cannot open the voxel file");
    }
    return file;
}

}  // namespace

int bytesPerSample(SampleType type) {
    switch (type) {
        case SampleType::kUInt8:   return 1;
        case SampleType::kUInt16:  return 2;
        case SampleType::kFloat32: return 4;
    }
    return 0;
}

std::string sampleTypeName(SampleType type) {
    switch (type) {
        case SampleType::kUInt8:   return "8-bit unsigned integer";
        case SampleType::kUInt16:  return "16-bit unsigned integer";
        case SampleType::kFloat32: return "32-bit float";
    }
    return "";
}

ScanInfo parseVgi(const std::string& text, const std::string& name) {
    const std::map<std::string, std::string> values = firstValues(text);
    const auto value_of = [&values](const char* key) {
        const auto found = values.find(key);
        return found == values.end() ? std::string() : found->second;
    };

    ScanInfo info;

    int size[3] = {0, 0, 0};
    if (values.count("size") == 0) {
        throw std::runtime_error(name + ": the key 'size' is missing");
    }
    if (!readNumbers(value_of("size"), 3, size)) {
        throw std::runtime_error(name + ": 'size' must be three whole numbers, found '" +
                                 value_of("size") + "'");
    }
    if (size[0] < 1 || size[1] < 1 || size[2] < 1) {
        throw std::runtime_error(name + ": every volume dimension must be at least 1");
    }
    info.size = VolumeSize{size[0], size[1], size[2]};

    int bits = 0;
    if (values.count("bitsperelement") == 0) {
        throw std::runtime_error(name + ": the key 'bitsperelement' is missing");
    }
    if (!readNumbers(value_of("bitsperelement"), 1, &bits)) {
        throw std::runtime_error(name + ": 'bitsperelement' must be a whole number, found '" +
                                 value_of("bitsperelement") + "'");
    }
    const std::string datatype =
        values.count("datatype") ? lowered(value_of("datatype")) : "unsigned integer";
    if (datatype == "unsigned integer" && bits == 8) {
        info.sample_type = SampleType::kUInt8;
    } else if (datatype == "unsigned integer" && bits == 16) {
        info.sample_type = SampleType::kUInt16;
    } else if (datatype == "float" && bits == 32) {
        info.sample_type = SampleType::kFloat32;
    } else {
        throw std::runtime_error(name + ": sample type '" + datatype + "' with " +
                                 std::to_string(bits) + " bits is not supported (supported: " +
                                 "unsigned integer with 8 or 16 bits, float with 32 bits)");
    }

    if (values.count("skipheader")) {
        long long header = 0;
        if (!readNumbers(value_of("skipheader"), 1, &header) || header < 0) {
            throw std::runtime_error(name + ": 'skipheader' must be a whole number >= 0, found '" +
                                     value_of("skipheader") + "'");
        }
        info.header_bytes = header;
    }

    if (values.count("resolution")) {
        if (!readNumbers(value_of("resolution"), 3, info.voxel_size) || info.voxel_size[0] <= 0.0 ||
            info.voxel_size[1] <= 0.0 || info.voxel_size[2] <= 0.0) {
            throw std::runtime_error(name + ": 'resolution' must be three positive numbers, found '" +
                                     value_of("resolution") + "'");
        }
        info.unit = value_of("unit");
    }
    return info;
}

ScanInfo openScan(const std::string& path) {
    std::error_code error;
    fs::path vgi_path;
    fs::path volume_path;

    if (fs::is_directory(path, error)) {
        std::vector<fs::path> found;
        for (const fs::directory_entry& entry : fs::directory_iterator(path)) {
            if (entry.is_regular_file() && lowered(entry.path().extension().string()) == ".vgi") {
                found.push_back(entry.path());
            }
        }
        if (found.size() != 1) {
            throw std::runtime_error(path + ": expected one .vgi file in the folder, found " +
                                     std::to_string(found.size()));
        }
        vgi_path = found.front();
        volume_path = vgi_path;
        volume_path.replace_extension(".vol");
    } else if (lowered(fs::path(path).extension().string()) == ".vgi") {
        vgi_path = path;
        volume_path = vgi_path;
        volume_path.replace_extension(".vol");
    } else {
        volume_path = path;
        vgi_path = volume_path;
        vgi_path.replace_extension(".vgi");
    }

    std::ifstream vgi(vgi_path, std::ios::binary);
    if (!vgi) {
        throw std::runtime_error(vgi_path.generic_string() + ": cannot open the .vgi file");
    }
    std::ostringstream text;
    text << vgi.rdbuf();

    ScanInfo info = parseVgi(text.str(), vgi_path.generic_string());
    info.volume_path = volume_path.generic_string();
    checkScan(info);
    return info;
}

ScanInfo describeRaw(const std::string& volume_path, const VolumeSize& size,
                     SampleType sample_type, ByteOrder byte_order, std::int64_t header_bytes) {
    ScanInfo info;
    info.volume_path = volume_path;
    info.size = size;
    info.sample_type = sample_type;
    info.byte_order = byte_order;
    info.header_bytes = header_bytes;
    checkScan(info);
    return info;
}

std::int64_t expectedFileSize(const ScanInfo& info) {
    return info.header_bytes + voxelCount(info.size) * bytesPerSample(info.sample_type);
}

void checkScan(const ScanInfo& info) {
    checkDescription(info);
    std::error_code error;
    const std::uintmax_t actual = fs::file_size(info.volume_path, error);
    if (error) {
        throw std::runtime_error(info.volume_path + ": cannot open the voxel file");
    }
    const std::int64_t expected = expectedFileSize(info);
    if (actual != static_cast<std::uintmax_t>(expected)) {
        throw std::runtime_error(
            info.volume_path + ": the file has " + std::to_string(actual) + " bytes, but " +
            std::to_string(info.size.width) + " x " + std::to_string(info.size.height) + " x " +
            std::to_string(info.size.depth) + " voxels of type " + sampleTypeName(info.sample_type) +
            " with a header of " + std::to_string(info.header_bytes) + " bytes need " +
            std::to_string(expected));
    }
}

std::vector<float> readRegion(const ScanInfo& info, const Region& region) {
    checkDescription(info);
    if (!regionFits(region, info.size)) {
        throw std::runtime_error(info.volume_path + ": the requested region does not fit the volume");
    }
    std::ifstream file = openVolumeFile(info);

    const std::size_t row = static_cast<std::size_t>(region.width);
    const std::size_t rows = static_cast<std::size_t>(region.height);
    std::vector<float> values(row * rows * static_cast<std::size_t>(region.depth));
    std::vector<unsigned char> buffer;

    // Rows that span the full width follow each other in the file without a
    // gap, so a whole slice of the region is read at once.
    const bool full_width = region.x0 == 0 && region.width == info.size.width;
    float* out = values.data();
    for (int z = region.z0; z < region.z0 + region.depth; ++z) {
        const std::int64_t slice_start =
            (static_cast<std::int64_t>(z) * info.size.height + region.y0) * info.size.width + region.x0;
        if (full_width) {
            readSamples(file, info, slice_start, row * rows, buffer, out);
            out += row * rows;
        } else {
            for (std::size_t y = 0; y < rows; ++y) {
                readSamples(file, info, slice_start + static_cast<std::int64_t>(y) * info.size.width,
                            row, buffer, out);
                out += row;
            }
        }
    }
    return values;
}

std::vector<float> readVolume(const ScanInfo& info) {
    return readRegion(info, wholeVolume(info.size));
}

ScanStatistics scanStatistics(const ScanInfo& info) {
    checkDescription(info);
    std::ifstream file = openVolumeFile(info);

    const std::size_t slice = static_cast<std::size_t>(info.size.width) * info.size.height;
    std::vector<float> values(slice);
    std::vector<unsigned char> buffer;

    ScanStatistics result;
    double minimum = std::numeric_limits<double>::infinity();
    double maximum = -std::numeric_limits<double>::infinity();
    // Sums are kept per slice and then added up, which keeps the rounding
    // error small on volumes of billions of voxels.
    long double sum = 0.0L;
    long double sum_squares = 0.0L;
    const bool integer_type = info.sample_type != SampleType::kFloat32;
    const float high_limit = info.sample_type == SampleType::kUInt8 ? 255.0f : 65535.0f;
    for (int z = 0; z < info.size.depth; ++z) {
        readSamples(file, info, static_cast<std::int64_t>(z) * static_cast<std::int64_t>(slice), slice,
                    buffer, values.data());
        double slice_sum = 0.0;
        double slice_squares = 0.0;
        for (const float value : values) {
            if (!std::isfinite(value)) {
                ++result.non_finite;
                continue;
            }
            if (integer_type && value == 0.0f) {
                ++result.at_low_limit;
            } else if (integer_type && value == high_limit) {
                ++result.at_high_limit;
            }
            const double v = value;
            minimum = std::min(minimum, v);
            maximum = std::max(maximum, v);
            slice_sum += v;
            slice_squares += v * v;
        }
        sum += slice_sum;
        sum_squares += slice_squares;
    }

    result.voxels = voxelCount(info.size);
    const std::int64_t finite = result.voxels - result.non_finite;
    if (finite > 0) {
        const long double mean = sum / finite;
        const long double variance = sum_squares / finite - mean * mean;
        result.minimum = minimum;
        result.maximum = maximum;
        result.mean = static_cast<double>(mean);
        result.std_dev = variance > 0.0L ? static_cast<double>(std::sqrt(variance)) : 0.0;
    }
    return result;
}

}  // namespace denoiser
