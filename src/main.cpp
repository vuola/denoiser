// denoiser - 3D Non-Local Means denoising of CT scans.
// Derived from the 3D C++ filter by Samuel Väisänen (2024) and, through it,
// from the ImageJ plugin by Pascal Behnel and Thorsten Wagner.
// Licensed under the GNU General Public License v3; see LICENSE.
//
// Command-line entry point. The command "info" is implemented (step 1.1);
// "noise" and "run" are recognised and follow in later steps of phase 1
// (README sections 1.6 and 2.4).

#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

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
           "  denoiser info  <scan> [raw options]      show what the scan contains\n"
           "  denoiser noise <scan folder>             estimate the noise sigma\n"
           "  denoiser run   <scan folder> [options]   denoise the scan\n"
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
           "The commands noise and run are not implemented yet.\n";
}

/**
 * @brief Read a whole number from a command-line word.
 * @param word The word, e.g. "64".
 * @param what What the number is, for the error message.
 * @return The number.
 * @throws std::invalid_argument if the word is not a whole number.
 */
std::int64_t to_integer(const std::string& word, const std::string& what) {
    std::istringstream stream(word);
    std::int64_t number = 0;
    std::string rest;
    if (!(stream >> number) || (stream >> rest)) {
        throw std::invalid_argument(what + " must be a whole number, found '" + word + "'");
    }
    return number;
}

/**
 * @brief Read a volume dimension from a command-line word.
 * @param word The word, e.g. "64".
 * @param what Which dimension it is, for the error message.
 * @return The dimension, between 1 and 1 000 000.
 * @throws std::invalid_argument if the word is not a whole number in range.
 */
int to_dimension(const std::string& word, const std::string& what) {
    const std::int64_t number = to_integer(word, what);
    if (number < 1 || number > 1000000) {
        throw std::invalid_argument(what + " must be between 1 and 1000000, found '" + word + "'");
    }
    return static_cast<int>(number);
}

/**
 * @brief Build the scan description from the words after "denoiser info".
 *
 * Without --raw the scan is opened through its .vgi file; with --raw the
 * description is taken from the options.
 *
 * @param words The command-line words after the command name.
 * @return The scan description, checked against the voxel file.
 * @throws std::invalid_argument if the words are not a valid combination.
 * @throws std::runtime_error if the scan cannot be opened.
 */
denoiser::ScanInfo scan_from_words(const std::vector<std::string>& words) {
    std::string path;
    bool raw = false;
    bool raw_option_seen = false;
    denoiser::VolumeSize size;
    denoiser::SampleType type = denoiser::SampleType::kUInt16;
    denoiser::ByteOrder order = denoiser::ByteOrder::kLittleEndian;
    std::int64_t header = 0;

    for (std::size_t i = 0; i < words.size(); ++i) {
        const std::string word = words[i];
        const auto next = [&](const std::string& what) -> const std::string& {
            if (i + 1 >= words.size()) {
                throw std::invalid_argument(word + " needs " + what);
            }
            return words[++i];
        };
        if (word == "--raw") {
            raw = true;
            size.width = to_dimension(next("three numbers"), "the width");
            size.height = to_dimension(next("three numbers"), "the height");
            size.depth = to_dimension(next("three numbers"), "the depth");
        } else if (word == "--type") {
            raw_option_seen = true;
            const std::string& name = next("u8, u16 or f32");
            if (name == "u8") {
                type = denoiser::SampleType::kUInt8;
            } else if (name == "u16") {
                type = denoiser::SampleType::kUInt16;
            } else if (name == "f32") {
                type = denoiser::SampleType::kFloat32;
            } else {
                throw std::invalid_argument("--type must be u8, u16 or f32, found '" + name + "'");
            }
        } else if (word == "--big-endian") {
            raw_option_seen = true;
            order = denoiser::ByteOrder::kBigEndian;
        } else if (word == "--header") {
            raw_option_seen = true;
            header = to_integer(next("a number of bytes"), "the header size");
        } else if (!word.empty() && word[0] == '-') {
            throw std::invalid_argument("unknown option '" + word + "'");
        } else if (path.empty()) {
            path = word;
        } else {
            throw std::invalid_argument("more than one scan given: '" + path + "' and '" + word + "'");
        }
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
 * @brief The command "info": print what a scan contains.
 *
 * Prints the description, confirms the file size and prints the value
 * statistics of the whole volume.
 *
 * @param words The command-line words after "info".
 * @return kOk, kUsageError for a wrong command line, kFailed if the scan
 *         cannot be read.
 */
int run_info(const std::vector<std::string>& words) {
    denoiser::ScanInfo info;
    try {
        info = scan_from_words(words);
    } catch (const std::invalid_argument& error) {
        std::cerr << "denoiser info: " << error.what() << "\n";
        return kUsageError;
    } catch (const std::runtime_error& error) {
        std::cerr << "denoiser info: " << error.what() << "\n";
        return kFailed;
    }

    const bool little = info.byte_order == denoiser::ByteOrder::kLittleEndian;
    std::cout << "Voxel file:   " << info.volume_path << "\n"
              << "Size:         " << info.size.width << " x " << info.size.height << " x "
              << info.size.depth << " voxels (" << denoiser::voxelCount(info.size) << ")\n"
              << "Sample type:  " << denoiser::sampleTypeName(info.sample_type);
    if (info.sample_type != denoiser::SampleType::kUInt8) {
        std::cout << ", " << (little ? "little-endian" : "big-endian");
    }
    std::cout << "\nHeader:       " << info.header_bytes << " bytes\n"
              << "Voxel size:   ";
    if (info.voxel_size[0] > 0.0) {
        std::cout << info.voxel_size[0] << " x " << info.voxel_size[1] << " x "
                  << info.voxel_size[2] << " " << info.unit << "\n";
    } else {
        std::cout << "not stated\n";
    }
    std::cout << "File size:    " << denoiser::expectedFileSize(info) << " bytes, as expected\n";

    try {
        const denoiser::ScanStatistics stats = denoiser::scanStatistics(info);
        std::cout << std::setprecision(8)
                  << "Minimum:      " << stats.minimum << "\n"
                  << "Maximum:      " << stats.maximum << "\n"
                  << "Mean:         " << stats.mean << "\n"
                  << "Std. dev.:    " << stats.std_dev << "\n";
        if (info.sample_type == denoiser::SampleType::kFloat32) {
            std::cout << "NaN or Inf:   " << stats.non_finite << " voxels\n";
        }
    } catch (const std::runtime_error& error) {
        std::cerr << "denoiser info: " << error.what() << "\n";
        return kFailed;
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

    if (command == "--help" || command == "-h") {
        print_usage(std::cout);
        return kOk;
    }
    if (command == "--version") {
        std::cout << denoiser::version() << "\n";
        return kOk;
    }
    if (command == "info") {
        return run_info(std::vector<std::string>(argv + 2, argv + argc));
    }
    if (command == "noise" || command == "run") {
        std::cerr << "denoiser: command '" << command << "' is not implemented yet\n";
        return kNotImplemented;
    }

    std::cerr << "denoiser: unknown command '" << command << "'\n\n";
    print_usage(std::cerr);
    return kUsageError;
}
