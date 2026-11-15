#!/usr/bin/env python3
"""Filter a volume along z only: moving median or moving mean per voxel column.

This is the simple filter used on the rods before this project (ImageJ
"Median 3D" with radii x = 0, y = 0, z = 16; README section 5.5). It serves as
a competitor in the baseline: a Non-Local Means result should beat it.

For every voxel the window is the 2 * radius + 1 voxels above and below it in
the same (x, y) column. Near the first and last slice the window is cut at the
volume's end, so it holds fewer voxels there. The median of an even number of
values is the mean of the two middle ones. Results are rounded to whole
counts.

Output, in the folder given with --out (default build/zfilter/):

    <scan>_z<kind>_r<radius>.vol + .vgi     16-bit unsigned, same size as the input

Usage (from the repository root):
    python tests/phantom/z_filter.py build/phantom/phantom_coarse --kind median --radius 16

The volume is held in memory, so the script is meant for phantoms and small
crops. Needs only the Python standard library.
"""

import argparse
import array
import os
import sys

REPO_ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(REPO_ROOT, "tests", "data"))
sys.path.insert(0, os.path.join(REPO_ROOT, "tests", "reference"))
from make_crops import read_vgi, write_vgi  # noqa: E402
from slice_images import find_volume        # noqa: E402


def filter_column(column, radius, kind):
    """Filter one column of values along its length.

    Parameters:
        column  sequence of integer values, one per slice
        radius  half window size in slices, >= 0
        kind    "median" or "mean"

    Returns a list of integers of the same length.
    """
    depth = len(column)
    result = [0] * depth
    if kind == "mean":
        sums = [0] * (depth + 1)
        for z, value in enumerate(column):
            sums[z + 1] = sums[z] + value
        for z in range(depth):
            first = max(0, z - radius)
            last = min(depth, z + radius + 1)
            count = last - first
            result[z] = (2 * (sums[last] - sums[first]) + count) // (2 * count)
    else:
        for z in range(depth):
            window = sorted(column[max(0, z - radius):z + radius + 1])
            middle = len(window) // 2
            if len(window) % 2:
                result[z] = window[middle]
            else:
                result[z] = (window[middle - 1] + window[middle] + 1) // 2
    return result


def z_filter(scan_path, kind, radius, out_dir):
    """Filter the scan at `scan_path` along z and write the result.

    Parameters:
        scan_path  scan folder or .vol file with its .vgi (16-bit unsigned)
        kind       "median" or "mean"
        radius     half window size in slices, >= 0
        out_dir    folder for the result (created if missing)

    Returns the path of the written .vol file.
    Raises ValueError if the scan cannot be read.
    """
    vol_path, (width, height, depth) = find_volume(scan_path)
    info = read_vgi(os.path.splitext(vol_path)[0] + ".vgi")
    voxels = width * height * depth
    values = array.array("H")
    with open(vol_path, "rb") as handle:
        values.fromfile(handle, voxels)
    if sys.byteorder == "big":
        values.byteswap()

    plane = width * height
    result = array.array("H", bytes(2 * voxels))
    for start in range(plane):
        result[start::plane] = array.array("H", filter_column(values[start::plane], radius, kind))

    name = f"{os.path.splitext(os.path.basename(vol_path))[0]}_z{kind}_r{radius}"
    os.makedirs(out_dir, exist_ok=True)
    out_vol = os.path.join(out_dir, name + ".vol")
    if sys.byteorder == "big":
        result.byteswap()
    with open(out_vol, "wb") as handle:
        result.tofile(handle)
    write_vgi(os.path.join(out_dir, name + ".vgi"), name + ".vol", (width, height, depth), info)
    return out_vol


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("scan", help="scan folder or .vol file")
    parser.add_argument("--kind", choices=["median", "mean"], default="median")
    parser.add_argument("--radius", type=int, default=16,
                        help="half window size in slices (default: 16, a window of 33)")
    parser.add_argument("--out", default=os.path.join(REPO_ROOT, "build", "zfilter"),
                        help="output folder (default: build/zfilter)")
    args = parser.parse_args()
    if args.radius < 0:
        parser.error("--radius must not be negative")

    try:
        out_vol = z_filter(args.scan, args.kind, args.radius, args.out)
    except (ValueError, OSError) as error:
        print(f"Failed: {error}", file=sys.stderr)
        return 1
    print(os.path.relpath(out_vol))
    return 0


if __name__ == "__main__":
    sys.exit(main())
