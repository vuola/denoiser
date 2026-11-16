// denoiser - 3D Non-Local Means filter for CT volumes.
// Derived from the C++ filter by Samuel Väisänen (2024), itself a port of the
// ImageJ plugin by Pascal Behnel and Thorsten Wagner. GPL-3.0, see LICENSE.

/// @file
/// Estimate of the noise standard deviation (sigma) of a scan.
///
/// The filter needs sigma: it sets the patch size, the search size and the
/// smoothing strength. This file reproduces the method of the earlier script
/// `estimate_noise.py` (README section 2.2), so that old and new runs use the
/// same number, and adds two estimators that use differences along z.
///
/// A number of slices, evenly spread over the depth, is examined. Each slice
/// gives one value per estimator; the values of all slices are then
/// summarised by median, mean and quartiles.
///
/// Estimators per slice (s = the slice, n = its next slice in z):
///
/// | Name         | Formula                                              |
/// |--------------|------------------------------------------------------|
/// | std_dev      | standard deviation of s                              |
/// | mad          | 1.4826 * MAD(s)                                      |
/// | diff         | std of all x and y neighbour differences / sqrt(2)   |
/// | resid        | 1.4826 * MAD(s - mean of the 3 x 3 neighbourhood)    |
/// | zdiff        | std of (n - s) / sqrt(2)                             |
/// | zdiff_robust | 1.4826 * MAD(n - s) / sqrt(2)                        |
///
/// MAD(v) is the median of |v - median(v)|; 1.4826 * MAD equals the standard
/// deviation for normally distributed values. The first five are the
/// estimators of the script; its recommendation is the mean of the medians
/// of `diff` and `resid`. `zdiff_robust` is new.
///
/// What each estimator measures besides noise:
/// - `std_dev` and `mad` measure the whole spread of the slice, materials
///   included. They are upper limits, not noise estimates.
/// - `diff` and `resid` contain the material edges within the slice.
/// - `zdiff` and `zdiff_robust` contain only what changes from one slice to
///   the next. In a rod with fibres along z that is little.
/// - All difference estimators assume that the noise of neighbouring voxels
///   is independent. If it is not, they read low. differenceByLag() shows
///   whether that is the case.
///
/// All results are in the units of the stored values (counts for a 16-bit
/// scan). Multiply by 255 / (high - low) of the intensity mapping to get
/// sigma on the working scale.

#ifndef DENOISER_NOISE_ESTIMATE_HPP
#define DENOISER_NOISE_ESTIMATE_HPP

#include <vector>

#include "scan_reader.hpp"
#include "volume.hpp"

namespace denoiser {

/// The estimates from one slice. See the file comment for the formulas.
struct SliceNoise {
    double std_dev = 0.0;       ///< standard deviation of the slice
    double mad = 0.0;           ///< 1.4826 * MAD of the slice
    double diff = 0.0;          ///< from x and y neighbour differences; 0 for a single voxel
    double resid = 0.0;         ///< from the residual against the 3 x 3 mean
    double zdiff = 0.0;         ///< from the difference to the next slice; NaN without a next slice
    double zdiff_robust = 0.0;  ///< robust version of zdiff; NaN without a next slice
};

/// Median, mean and quartiles of one estimator over the examined slices.
/// Quartiles are interpolated linearly between sorted values, as numpy does.
struct NoiseSummary {
    int count = 0;        ///< number of slices that gave a value
    double median = 0.0;  ///< NaN if count == 0, as are the other three
    double mean = 0.0;
    double q1 = 0.0;      ///< 25th percentile
    double q3 = 0.0;      ///< 75th percentile
};

/// Result of estimateNoise().
struct NoiseEstimate {
    std::vector<int> slices;            ///< examined slices, counted from 0 in the scan
    std::vector<SliceNoise> per_slice;  ///< the estimates of each examined slice
    NoiseSummary std_dev;
    NoiseSummary mad;
    NoiseSummary diff;
    NoiseSummary resid;
    NoiseSummary zdiff;
    NoiseSummary zdiff_robust;
    /// Sigma as the earlier script recommends it: the mean of the medians of
    /// `diff` and `resid`.
    double recommended = 0.0;
};

/// One line of differenceByLag(): the noise estimate from differences
/// between voxels `lag` voxels apart, along each axis.
struct LagSigma {
    int lag = 0;     ///< distance between the two voxels
    double x = 0.0;  ///< 1.4826 * MAD(difference along x) / sqrt(2), median over slices; NaN if the region is narrower than lag + 1
    double y = 0.0;  ///< the same along y
    double z = 0.0;  ///< the same along z
};

/// Computes the estimates of one slice.
///
/// Edge voxels of the 3 x 3 mean use the nearest voxel inside the slice in
/// place of a missing neighbour.
///
/// @param slice   width * height values, x fastest
/// @param next    the next slice in z, same layout, or nullptr if there is none
/// @param width   slice width in voxels, >= 1
/// @param height  slice height in voxels, >= 1
/// @return the six estimates
/// @throws std::invalid_argument if width or height is below 1
SliceNoise sliceNoise(const float* slice, const float* next, int width, int height);

/// Summarises one estimator over the slices. NaN values are left out.
/// @param values  one value per slice
NoiseSummary summarise(const std::vector<double>& values);

/// Picks `count` slice numbers spread evenly from the first to the last
/// slice: slice i is i * (depth - 1) / (count - 1), rounded down. A count
/// above depth is reduced to depth; a count of 1 gives slice 0.
/// @param depth  number of slices, >= 1
/// @param count  number of slices wanted, >= 1
/// @return slice numbers in rising order, counted from 0
/// @throws std::invalid_argument if depth or count is below 1
std::vector<int> sampleSlices(int depth, int count);

/// Estimates the noise of a region of a scan.
///
/// The region is read one slice at a time, so memory use is two slices.
/// The last slice of the region has no next slice and gives no z estimates.
///
/// @param info    the scan description
/// @param region  the part of the scan to use
/// @param slices  number of slices to examine (the script used 20)
/// @return the estimates, in units of the stored values
/// @throws std::invalid_argument if the region does not fit or slices < 1
/// @throws std::runtime_error if the scan cannot be read or the examined
///         slices contain NaN or infinite values
NoiseEstimate estimateNoise(const ScanInfo& info, const Region& region, int slices = 20);

/// Estimates the noise from differences between voxels 1, 2, ... max_lag
/// voxels apart, along each axis.
///
/// If the noise of neighbouring voxels were independent and the material did
/// not change over these distances, every entry would equal sigma. Entries
/// that rise with the lag and then level off show noise that is correlated
/// between neighbours: the level is the better sigma, and the entry at lag 1
/// is too low by the factor sqrt(1 - correlation). Entries that keep rising
/// show material structure.
///
/// @param info     the scan description
/// @param region   the part of the scan to use; its depth must exceed max_lag
/// @param slices   number of slices to examine
/// @param max_lag  largest distance, >= 1
/// @return one entry per lag, in units of the stored values
/// @throws std::invalid_argument if an argument is out of range
/// @throws std::runtime_error if the scan cannot be read or the examined
///         slices contain NaN or infinite values
std::vector<LagSigma> differenceByLag(const ScanInfo& info, const Region& region, int slices,
                                      int max_lag);

}  // namespace denoiser

#endif  // DENOISER_NOISE_ESTIMATE_HPP
