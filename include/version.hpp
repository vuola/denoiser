// denoiser - 3D Non-Local Means denoising of CT scans.
// Derived from the 3D C++ filter by Samuel Väisänen (2024) and, through it,
// from the ImageJ plugin by Pascal Behnel and Thorsten Wagner.
// Licensed under the GNU General Public License v3; see LICENSE.

#ifndef DENOISER_VERSION_HPP
#define DENOISER_VERSION_HPP

#include <string>

namespace denoiser {

/**
 * @brief Version of the program.
 *
 * The number is "major.minor.patch". Major and minor follow the project
 * step that the code implements (README section 1.6); for example "0.1.0"
 * is step 0.1.
 *
 * @return The version as text, e.g. "0.1.0". Never empty.
 */
std::string version();

}  // namespace denoiser

#endif  // DENOISER_VERSION_HPP
