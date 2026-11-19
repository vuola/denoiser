// Command-line program around the brute-force filter, with the same arguments
// and file format as the original program, so that
// tests/reference/run_reference.py can run either one:
//
//     denoiser_bruteforce input output sigma [options]
//
//   input   headerless 32-bit big-endian floats, x fastest; the volume size is
//           read from the end of the file name: "<anything>WxHxD.<ext>"
//   output  the filtered volume in the same format
//   sigma   noise standard deviation on the 0...255 scale; patch radius,
//           search radius and h follow the table of the original program
//
// Options, each replacing one choice of that table (README chapter 12):
//   --patch PX PY PZ     patch radius along x, y and z
//   --search SX SY SZ    search radius along x, y and z
//   --h-factor F         h = F * sigma
//   --normalise          patch distance as mean over the patch, not as sum
//   --kernel K           rational (default) or exponential
//
// and options that choose the computation instead of the definition
// (README chapter 13):
//   --fast               the fast filter (src/nlm_fast.cpp)
//   --single             with --fast: distances in single precision (default)
//   --double             with --fast: distances in double precision
//   --block N            with --fast: block edge in voxels (default 64)
//
// Exit code 0 on success, 1 on failure with the reason on standard error.

#include "nlm_bruteforce.hpp"
#include "nlm_fast.hpp"

#include <cctype>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <exception>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

/// Reads the volume size from a file name ending in "WxHxD.<ext>".
/// Throws std::runtime_error if the name does not end that way.
denoiser::VolumeSize sizeFromName(const std::string& name) {
    const std::string stem = name.substr(0, name.rfind('.'));
    int numbers[3] = {0, 0, 0};
    std::size_t end = stem.size();
    for (int i = 2; i >= 0; --i) {
        std::size_t begin = end;
        while (begin > 0 && std::isdigit(static_cast<unsigned char>(stem[begin - 1]))) {
            --begin;
        }
        if (begin == end || end - begin > 6) {
            throw std::runtime_error("cannot read the volume size from the file name");
        }
        numbers[i] = std::stoi(stem.substr(begin, end - begin));
        if (i > 0) {
            if (begin == 0 || stem[begin - 1] != 'x') {
                throw std::runtime_error("cannot read the volume size from the file name");
            }
            end = begin - 1;
        }
    }
    return denoiser::VolumeSize{numbers[0], numbers[1], numbers[2]};
}

/// Reads `count` big-endian 32-bit floats; throws if the file has another size.
std::vector<float> readBigEndianFloats(const std::string& path, std::size_t count) {
    std::ifstream file(path, std::ios::binary);
    if (!file) {
        throw std::runtime_error("cannot open the input file");
    }
    std::vector<unsigned char> bytes(count * 4);
    file.read(reinterpret_cast<char*>(bytes.data()), static_cast<std::streamsize>(bytes.size()));
    if (static_cast<std::size_t>(file.gcount()) != bytes.size() || file.peek() != EOF) {
        throw std::runtime_error("the file size does not match the size in its name");
    }
    std::vector<float> values(count);
    for (std::size_t i = 0; i < count; ++i) {
        const std::uint32_t bits = (std::uint32_t{bytes[4 * i]} << 24) | (std::uint32_t{bytes[4 * i + 1]} << 16)
                                 | (std::uint32_t{bytes[4 * i + 2]} << 8) | std::uint32_t{bytes[4 * i + 3]};
        std::memcpy(&values[i], &bits, 4);
    }
    return values;
}

/// Writes floats as big-endian 32-bit values; throws on failure.
void writeBigEndianFloats(const std::string& path, const std::vector<float>& values) {
    std::vector<unsigned char> bytes(values.size() * 4);
    for (std::size_t i = 0; i < values.size(); ++i) {
        std::uint32_t bits = 0;
        std::memcpy(&bits, &values[i], 4);
        bytes[4 * i] = static_cast<unsigned char>(bits >> 24);
        bytes[4 * i + 1] = static_cast<unsigned char>(bits >> 16);
        bytes[4 * i + 2] = static_cast<unsigned char>(bits >> 8);
        bytes[4 * i + 3] = static_cast<unsigned char>(bits);
    }
    std::ofstream file(path, std::ios::binary | std::ios::trunc);
    file.write(reinterpret_cast<const char*>(bytes.data()), static_cast<std::streamsize>(bytes.size()));
    if (!file) {
        throw std::runtime_error("cannot write the output file");
    }
}

}  // namespace

/// Applies the options that follow the three fixed arguments to `params`.
/// Throws std::runtime_error for an unknown option or a missing or bad value.
void applyOptions(int argc, char** argv, denoiser::NlmParams& params, bool& fast,
                  denoiser::NlmFastOptions& fast_options) {
    int i = 4;
    const auto value = [&]() -> std::string {
        if (i + 1 >= argc) {
            throw std::runtime_error(std::string(argv[i]) + " needs more values");
        }
        return argv[++i];
    };
    for (; i < argc; ++i) {
        const std::string option = argv[i];
        if (option == "--patch" || option == "--search") {
            int* radius = option == "--patch" ? params.patch_radius : params.search_radius;
            for (int axis = 0; axis < 3; ++axis) {
                radius[axis] = std::stoi(value());
            }
        } else if (option == "--h-factor") {
            params.h = std::stod(value()) * params.sigma;
        } else if (option == "--normalise") {
            params.normalise_distance = true;
        } else if (option == "--kernel") {
            const std::string kernel = value();
            if (kernel == "rational") {
                params.kernel = denoiser::NlmKernel::kRational;
            } else if (kernel == "exponential") {
                params.kernel = denoiser::NlmKernel::kExponential;
            } else {
                throw std::runtime_error("--kernel must be rational or exponential");
            }
        } else if (option == "--fast") {
            fast = true;
        } else if (option == "--single") {
            fast_options.single_precision = true;
        } else if (option == "--double") {
            fast_options.single_precision = false;
        } else if (option == "--block") {
            fast_options.block_size = std::stoi(value());
        } else {
            throw std::runtime_error("unknown option " + option);
        }
    }
}

int main(int argc, char** argv) {
    if (argc < 4) {
        std::fprintf(stderr, "Usage: denoiser_bruteforce input output sigma [options]\n"
                             "  input, output: 32-bit big-endian float files; the input name ends in WxHxD.ext\n"
                             "  sigma: noise standard deviation on the 0...255 scale\n"
                             "  options: --patch PX PY PZ, --search SX SY SZ, --h-factor F,\n"
                             "           --normalise, --kernel rational|exponential,\n"
                             "           --fast, --single, --double, --block N\n");
        return 1;
    }
    try {
        const denoiser::VolumeSize size = sizeFromName(argv[1]);
        denoiser::NlmParams params = denoiser::originalParams(std::stod(argv[3]));
        bool fast = false;
        denoiser::NlmFastOptions fast_options;
        applyOptions(argc, argv, params, fast, fast_options);
        const std::size_t voxels = static_cast<std::size_t>(size.width) * size.height * size.depth;
        const std::vector<float> volume = readBigEndianFloats(argv[1], voxels);
        writeBigEndianFloats(argv[2], fast ? denoiser::nlmFast(volume, size, params, fast_options)
                                           : denoiser::nlmBruteForce(volume, size, params));
    } catch (const std::exception& error) {
        std::fprintf(stderr, "Failed: %s\n", error.what());
        return 1;
    }
    return 0;
}
