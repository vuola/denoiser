// denoiser - 3D Non-Local Means denoising of CT scans.
// Derived from the 3D C++ filter by Samuel Väisänen (2024) and, through it,
// from the ImageJ plugin by Pascal Behnel and Thorsten Wagner.
// Licensed under the GNU General Public License v3; see LICENSE.
//
// Command-line entry point. In step 0.1 the commands are only recognised;
// they are implemented in phase 1 (README sections 1.6 and 2.4).

#include <iostream>
#include <string>

#include "version.hpp"

namespace {

/// Exit codes of the program.
enum ExitCode {
    kOk = 0,              ///< The command succeeded.
    kUsageError = 1,      ///< Unknown or missing command.
    kNotImplemented = 2,  ///< The command exists in the plan but not yet in code.
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
           "  denoiser info  <scan folder>             show what the scan contains\n"
           "  denoiser noise <scan folder>             estimate the noise sigma\n"
           "  denoiser run   <scan folder> [options]   denoise the scan\n"
           "  denoiser --version\n"
           "  denoiser --help\n"
           "\n"
           "The commands info, noise and run are not implemented yet.\n";
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
    if (command == "info" || command == "noise" || command == "run") {
        std::cerr << "denoiser: command '" << command << "' is not implemented yet\n";
        return kNotImplemented;
    }

    std::cerr << "denoiser: unknown command '" << command << "'\n\n";
    print_usage(std::cerr);
    return kUsageError;
}
