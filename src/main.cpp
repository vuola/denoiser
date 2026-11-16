// denoiser - 3D Non-Local Means denoising of CT scans.
// Derived from the 3D C++ filter by Samuel Väisänen (2024) and, through it,
// from the ImageJ plugin by Pascal Behnel and Thorsten Wagner.
// Licensed under the GNU General Public License v3; see LICENSE.
//
// Command-line entry point. The commands "info" (step 1.1), "crop" and
// "profile" (step 1.3) are implemented; "noise" and "run" are recognised and
// follow in later steps of phase 1 (README sections 1.6 and 2.4).

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "preprocess.hpp"
#include "radial_profile.hpp"
#include "scan_reader.hpp"
#include "version.hpp"

namespace {

/// Exit codes of the program.
enum ExitCode {
    kOk = 0,              ///< The command succeeded.
    kUsageError = 1,      ///< Unknown or missing command, or wrong options.
    kNotImplemented = 2,  ///< The command exists in the plan but not yet in code.
    kFailed = 3,          ///< The command was understood but could not be carried out.
};

/**
 * @brief Print the command summary.
 * @param out Stream to print to.
 */
void print_usage(std::ostream& out) {
    out << "denoiser " << denoiser::version()
        << " - 3D Non-Local Means denoising of CT scans\n"
           "\n"
           "Usage:\n"
           "  denoiser info    <scan>                    show what the scan contains\n"
           "  denoiser crop    <scan> <region> --out DIR write a region as a scan of its own\n"
           "  denoiser profile <scan> [options]          brightness against distance from the axis\n"
           "  denoiser noise   <scan>                    estimate the noise sigma\n"
           "  denoiser run     <scan> [options]          denoise the scan\n"
           "  denoiser --version\n"
           "  denoiser --help\n"
           "\n"
           "<scan> is a scan folder as delivered by the CT lab, its .vgi file or its\n"
           ".vol file. A raw voxel file without a .vgi file is described by hand:\n"
           "  --raw W H D        volume size in voxels (required for such a file)\n"
           "  --type T           u8, u16 (default) or f32\n"
           "  --big-endian       samples are stored highest byte first\n"
           "  --header N         bytes to skip before the first voxel (default 0)\n"
           "\n"
           "crop: <region> is one of\n"
           "  --roi X0 Y0 Z0 W H D   first voxel (counted from 0) and size\n"
           "  --centre-cube N        cube of N voxels in the middle of the scan\n"
           "The result is DIR/<name>.vol and .vgi, <name> being the last part of DIR,\n"
           "and DIR/crop_record.txt with the source and the region.\n"
           "\n"
           "profile: prints a table (CSV) with one line per ring around the axis\n"
           "  --centre X Y       the axis, in voxels (default: middle of the x-y plane)\n"
           "  --ring W           ring width in voxels (default 8)\n"
           "  --slices Z0 N      use N slices from slice Z0 (default: all)\n"
           "  --out FILE         write the table to FILE instead of the screen\n"
           "\n"
           "The commands noise and run are not implemented yet.\n";
}

/**
 * @brief Read a number from a command-line word.
 * @tparam Number std::int64_t or double.
 * @param word The word, e.g. "64".
 * @param what What the number is, for the error message.
 * @return The number.
 * @throws std::invalid_argument if the word is not such a number.
 */
template <typename Number>
Number to_number(const std::string& word, const std::string& what) {
    std::istringstream stream(word);
    Number number = 0;
    std::string rest;
    if (!(stream >> number) || (stream >> rest)) {
        throw std::invalid_argument(what + " must be a number, found '" + word + "'");
    }
    return number;
}

/**
 * @brief Read a voxel count or voxel coordinate from a command-line word.
 * @param word The word, e.g. "64".
 * @param what What the number is, for the error message.
 * @param lowest Smallest allowed value (0 for coordinates, 1 for sizes).
 * @return The number, between `lowest` and 1 000 000.
 * @throws std::invalid_argument if the word is not a whole number in range.
 */
int to_voxels(const std::string& word, const std::string& what, int lowest) {
    const std::int64_t number = to_number<std::int64_t>(word, what);
    if (number < lowest || number > 1000000) {
        throw std::invalid_argument(what + " must be between " + std::to_string(lowest) +
                                    " and 1000000, found '" + word + "'");
    }
    return static_cast<int>(number);
}

/**
 * @brief Take an option and its values out of the command-line words.
 *
 * @param words The command-line words; the option and its values are removed.
 * @param name The option, e.g. "--roi".
 * @param count Number of values that follow the option.
 * @param values Receives the values if the option is present.
 * @return true if the option was present.
 * @throws std::invalid_argument if the option is given twice or values are
 *         missing.
 */
bool take_option(std::vector<std::string>& words, const std::string& name, std::size_t count,
                 std::vector<std::string>& values) {
    bool found = false;
    for (std::size_t i = 0; i < words.size();) {
        if (words[i] != name) {
            ++i;
            continue;
        }
        if (found) {
            throw std::invalid_argument(name + " is given twice");
        }
        if (count > 0 && i + count >= words.size()) {
            throw std::invalid_argument(name + " needs " + std::to_string(count) +
                                        (count == 1 ? " value" : " values"));
        }
        values.assign(words.begin() + static_cast<std::ptrdiff_t>(i + 1),
                      words.begin() + static_cast<std::ptrdiff_t>(i + 1 + count));
        words.erase(words.begin() + static_cast<std::ptrdiff_t>(i),
                    words.begin() + static_cast<std::ptrdiff_t>(i + 1 + count));
        found = true;
    }
    return found;
}

/**
 * @brief Open the scan named on the command line.
 *
 * Takes the raw-file options and the scan path out of the words. Without
 * --raw the scan is opened through its .vgi file; with --raw the description
 * is taken from the options. Any word left over is an error, so the caller
 * takes its own options out first.
 *
 * @param words The command-line words after the command's own options were
 *        removed.
 * @return The scan description, checked against the voxel file.
 * @throws std::invalid_argument if the words are not a valid combination.
 * @throws std::runtime_error if the scan cannot be opened.
 */
denoiser::ScanInfo take_scan(std::vector<std::string> words) {
    std::vector<std::string> values;
    denoiser::VolumeSize size;
    const bool raw = take_option(words, "--raw", 3, values);
    if (raw) {
        size.width = to_voxels(values[0], "the width", 1);
        size.height = to_voxels(values[1], "the height", 1);
        size.depth = to_voxels(values[2], "the depth", 1);
    }

    bool raw_option_seen = false;
    denoiser::SampleType type = denoiser::SampleType::kUInt16;
    if (take_option(words, "--type", 1, values)) {
        raw_option_seen = true;
        if (values[0] == "u8") {
            type = denoiser::SampleType::kUInt8;
        } else if (values[0] == "u16") {
            type = denoiser::SampleType::kUInt16;
        } else if (values[0] == "f32") {
            type = denoiser::SampleType::kFloat32;
        } else {
            throw std::invalid_argument("--type must be u8, u16 or f32, found '" + values[0] + "'");
        }
    }
    denoiser::ByteOrder order = denoiser::ByteOrder::kLittleEndian;
    if (take_option(words, "--big-endian", 0, values)) {
        raw_option_seen = true;
        order = denoiser::ByteOrder::kBigEndian;
    }
    std::int64_t header = 0;
    if (take_option(words, "--header", 1, values)) {
        raw_option_seen = true;
        header = to_number<std::int64_t>(values[0], "the header size");
    }

    std::string path;
    for (const std::string& word : words) {
        if (!word.empty() && word[0] == '-') {
            throw std::invalid_argument("unknown option '" + word + "'");
        }
        if (!path.empty()) {
            throw std::invalid_argument("more than one scan given: '" + path + "' and '" + word + "'");
        }
        path = word;
    }
    if (path.empty()) {
        throw std::invalid_argument("no scan given");
    }
    if (!raw) {
        if (raw_option_seen) {
            throw std::invalid_argument("--type, --big-endian and --header need --raw W H D");
        }
        return denoiser::openScan(path);
    }
    return denoiser::describeRaw(path, size, type, order, header);
}

/**
 * @brief Print the description of a scan.
 * @param out Stream to print to.
 * @param info The scan description.
 */
void print_description(std::ostream& out, const denoiser::ScanInfo& info) {
    const bool little = info.byte_order == denoiser::ByteOrder::kLittleEndian;
    out << "Voxel file:   " << info.volume_path << "\n"
        << "Size:         " << info.size.width << " x " << info.size.height << " x "
        << info.size.depth << " voxels (" << denoiser::voxelCount(info.size) << ")\n"
        << "Sample type:  " << denoiser::sampleTypeName(info.sample_type);
    if (info.sample_type != denoiser::SampleType::kUInt8) {
        out << ", " << (little ? "little-endian" : "big-endian");
    }
    out << "\nHeader:       " << info.header_bytes << " bytes\n"
        << "Voxel size:   ";
    if (info.voxel_size[0] > 0.0) {
        out << info.voxel_size[0] << " x " << info.voxel_size[1] << " x " << info.voxel_size[2]
            << " " << info.unit << "\n";
    } else {
        out << "not stated\n";
    }
}

/**
 * @brief The command "info": print what a scan contains.
 *
 * Prints the description, confirms the file size and prints the value
 * statistics of the whole volume, including how many voxels sit at the ends
 * of the storable range.
 *
 * @param words The command-line words after "info".
 * @return kOk.
 * @throws std::invalid_argument for a wrong command line.
 * @throws std::runtime_error if the scan cannot be read.
 */
int run_info(const std::vector<std::string>& words) {
    const denoiser::ScanInfo info = take_scan(words);
    print_description(std::cout, info);
    std::cout << "File size:    " << denoiser::expectedFileSize(info) << " bytes, as expected\n";

    const denoiser::ScanStatistics stats = denoiser::scanStatistics(info);
    std::cout << std::setprecision(8)
              << "Minimum:      " << stats.minimum << "\n"
              << "Maximum:      " << stats.maximum << "\n"
              << "Mean:         " << stats.mean << "\n"
              << "Std. dev.:    " << stats.std_dev << "\n";
    if (info.sample_type == denoiser::SampleType::kFloat32) {
        std::cout << "NaN or Inf:   " << stats.non_finite << " voxels\n";
    } else {
        const double percent = 100.0 / static_cast<double>(stats.voxels);
        const int top = info.sample_type == denoiser::SampleType::kUInt8 ? 255 : 65535;
        std::cout << std::setprecision(3)
                  << "At 0:         " << stats.at_low_limit << " voxels ("
                  << static_cast<double>(stats.at_low_limit) * percent << " %)\n"
                  << "At " << std::left << std::setw(11) << (std::to_string(top) + ":")
                  << stats.at_high_limit << " voxels ("
                  << static_cast<double>(stats.at_high_limit) * percent << " %)\n";
    }
    return kOk;
}

/**
 * @brief The command "crop": write a region of a scan as a scan of its own.
 *
 * Writes DIR/<name>.vol, DIR/<name>.vgi and DIR/crop_record.txt, where
 * <name> is the last part of DIR. The record states the source and the
 * region, so that the crop can be placed back into the source's voxel grid.
 *
 * @param words The command-line words after "crop".
 * @return kOk.
 * @throws std::invalid_argument for a wrong command line.
 * @throws std::runtime_error if a file cannot be read or written.
 */
int run_crop(std::vector<std::string> words) {
    std::vector<std::string> roi;
    std::vector<std::string> cube;
    std::vector<std::string> out;
    const bool has_roi = take_option(words, "--roi", 6, roi);
    const bool has_cube = take_option(words, "--centre-cube", 1, cube);
    if (!take_option(words, "--out", 1, out)) {
        throw std::invalid_argument("--out DIR is required");
    }
    if (has_roi == has_cube) {
        throw std::invalid_argument("give the region with either --roi or --centre-cube");
    }
    const denoiser::ScanInfo source = take_scan(words);

    denoiser::Region region;
    if (has_cube) {
        region = denoiser::centredCube(source.size, to_voxels(cube[0], "the cube size", 1));
    } else {
        region = denoiser::Region{to_voxels(roi[0], "X0", 0), to_voxels(roi[1], "Y0", 0),
                                  to_voxels(roi[2], "Z0", 0), to_voxels(roi[3], "W", 1),
                                  to_voxels(roi[4], "H", 1),  to_voxels(roi[5], "D", 1)};
        if (!denoiser::regionFits(region, source.size)) {
            throw std::invalid_argument("the region does not fit the scan of " +
                                        std::to_string(source.size.width) + " x " +
                                        std::to_string(source.size.height) + " x " +
                                        std::to_string(source.size.depth) + " voxels");
        }
    }

    namespace fs = std::filesystem;
    const fs::path folder = fs::path(out[0]).lexically_normal();
    fs::path name = folder.filename();
    if (name.empty()) {
        name = folder.parent_path().filename();  // DIR was given with a trailing slash
    }
    if (name.empty() || name == "." || name == "..") {
        throw std::invalid_argument("--out needs a folder name, found '" + out[0] + "'");
    }
    std::error_code error;
    fs::create_directories(folder, error);
    if (error) {
        throw std::runtime_error(out[0] + ": cannot create the folder (" + error.message() + ")");
    }

    const std::string volume_path = (folder / name).generic_string() + ".vol";
    const denoiser::ScanInfo result = denoiser::cropScan(source, region, volume_path);

    const std::string record_path = (folder / "crop_record.txt").generic_string();
    std::ofstream record(record_path, std::ios::binary);
    record << "denoiser " << denoiser::version() << " crop\n"
           << "source = " << source.volume_path << "\n"
           << "source size = " << source.size.width << " " << source.size.height << " "
           << source.size.depth << "\n"
           << "first voxel (x y z, counted from 0) = " << region.x0 << " " << region.y0 << " "
           << region.z0 << "\n"
           << "size (x y z) = " << region.width << " " << region.height << " " << region.depth << "\n";
    if (!record) {
        throw std::runtime_error(record_path + ": cannot write the record");
    }

    print_description(std::cout, result);
    std::cout << "First voxel:  " << region.x0 << " " << region.y0 << " " << region.z0
              << " of the source\n";
    return kOk;
}

/**
 * @brief The command "profile": brightness against distance from the axis.
 *
 * Prints one CSV line per ring: radii, voxel count, mean, 10th, 50th and
 * 90th percentile, and the number of voxels at either end of the value
 * range.
 *
 * @param words The command-line words after "profile".
 * @return kOk.
 * @throws std::invalid_argument for a wrong command line.
 * @throws std::runtime_error if a file cannot be read or written.
 */
int run_profile(std::vector<std::string> words) {
    std::vector<std::string> centre;
    std::vector<std::string> ring;
    std::vector<std::string> slices;
    std::vector<std::string> out;
    const bool has_centre = take_option(words, "--centre", 2, centre);
    const bool has_ring = take_option(words, "--ring", 1, ring);
    const bool has_slices = take_option(words, "--slices", 2, slices);
    const bool has_out = take_option(words, "--out", 1, out);
    const denoiser::ScanInfo info = take_scan(words);

    denoiser::ProfileSettings settings = denoiser::defaultProfileSettings(info);
    if (has_centre) {
        settings.centre_x = to_number<double>(centre[0], "the centre's x");
        settings.centre_y = to_number<double>(centre[1], "the centre's y");
    }
    if (has_ring) {
        settings.ring_width = to_number<double>(ring[0], "the ring width");
    }
    if (has_slices) {
        settings.z0 = to_voxels(slices[0], "the first slice", 0);
        settings.depth = to_voxels(slices[1], "the number of slices", 1);
    }

    const std::vector<denoiser::RadialRing> rings = denoiser::radialProfile(info, settings);

    std::ofstream file;
    if (has_out) {
        file.open(out[0], std::ios::binary);
        if (!file) {
            throw std::runtime_error(out[0] + ": cannot write the table");
        }
    }
    std::ostream& table = has_out ? file : std::cout;
    table << "r_inner,r_outer,voxels,mean,p10,p50,p90,at_low,at_high\n" << std::setprecision(8);
    for (const denoiser::RadialRing& r : rings) {
        table << r.r_inner << "," << r.r_outer << "," << r.voxels << "," << r.mean << "," << r.p10
              << "," << r.p50 << "," << r.p90 << "," << r.at_low << "," << r.at_high << "\n";
    }
    if (has_out) {
        file.flush();
        if (!file) {
            throw std::runtime_error(out[0] + ": cannot write the table");
        }
        std::cout << rings.size() << " rings written to " << out[0] << "\n";
    }
    return kOk;
}

}  // namespace

/**
 * @brief Program entry point.
 * @param argc Number of command-line words, including the program name.
 * @param argv The command-line words.
 * @return An ExitCode value.
 */
int main(int argc, char* argv[]) {
    if (argc < 2) {
        print_usage(std::cerr);
        return kUsageError;
    }

    const std::string command = argv[1];
    const std::vector<std::string> words(argv + 2, argv + argc);

    if (command == "--help" || command == "-h") {
        print_usage(std::cout);
        return kOk;
    }
    if (command == "--version") {
        std::cout << denoiser::version() << "\n";
        return kOk;
    }
    if (command == "noise" || command == "run") {
        std::cerr << "denoiser: command '" << command << "' is not implemented yet\n";
        return kNotImplemented;
    }

    try {
        if (command == "info") {
            return run_info(words);
        }
        if (command == "crop") {
            return run_crop(words);
        }
        if (command == "profile") {
            return run_profile(words);
        }
    } catch (const std::invalid_argument& error) {
        std::cerr << "denoiser " << command << ": " << error.what() << "\n";
        return kUsageError;
    } catch (const std::runtime_error& error) {
        std::cerr << "denoiser " << command << ": " << error.what() << "\n";
        return kFailed;
    }

    std::cerr << "denoiser: unknown command '" << command << "'\n\n";
    print_usage(std::cerr);
    return kUsageError;
}
