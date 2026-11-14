#!/usr/bin/env python3
"""Compare two volume files voxel by voxel.

Both files must have the same type, chosen by the file extension:

    .vol   16-bit unsigned little-endian; differences are in counts
    .f32   32-bit little-endian float; differences are in the file's own unit

Prints the number of voxels, how many differ, the largest absolute difference
and the root-mean-square difference. With --tolerance the exit code tells
whether the largest difference is within the tolerance (0) or not (1).

The files are held in memory, so the script is meant for the small test scans.

Usage (from the repository root):
    python tests/reference/compare_volumes.py a.vol b.vol
    python tests/reference/compare_volumes.py a.f32 b.f32 --tolerance 0.001

Needs only the Python standard library.
"""

import argparse
import array
import math
import os
import sys

# File extension -> (array type code, bytes per voxel)
FORMATS = {".vol": ("H", 2), ".f32": ("f", 4)}


def read_volume(path):
    """Read a .vol or .f32 file into an array.

    Raises ValueError for another extension or a size that is not a whole
    number of voxels.
    """
    extension = os.path.splitext(path)[1].lower()
    if extension not in FORMATS:
        raise ValueError(f"{path}: extension must be .vol or .f32")
    type_code, voxel_bytes = FORMATS[extension]
    size = os.path.getsize(path)
    if size % voxel_bytes:
        raise ValueError(f"{path}: size is not a multiple of {voxel_bytes} bytes")
    values = array.array(type_code)
    with open(path, "rb") as volume:
        values.fromfile(volume, size // voxel_bytes)
    if sys.byteorder == "big":
        values.byteswap()
    return values


def compare(path_a, path_b):
    """Compare two volume files of the same type and size.

    Returns a dict with:
        voxels     number of voxels
        differing  number of voxels whose values differ
        max_abs    largest absolute difference
        rms        root-mean-square difference
    Raises ValueError if the types or sizes differ.
    """
    if os.path.splitext(path_a)[1].lower() != os.path.splitext(path_b)[1].lower():
        raise ValueError("the two files have different types")
    a = read_volume(path_a)
    b = read_volume(path_b)
    if len(a) != len(b):
        raise ValueError(f"different sizes: {len(a)} and {len(b)} voxels")

    differing = 0
    max_abs = 0.0
    sum_squares = 0.0
    for value_a, value_b in zip(a, b):
        difference = abs(value_a - value_b)
        if difference:
            differing += 1
            sum_squares += difference * difference
            if difference > max_abs:
                max_abs = difference
    rms = math.sqrt(sum_squares / len(a)) if len(a) else 0.0
    return {"voxels": len(a), "differing": differing, "max_abs": max_abs, "rms": rms}


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("a", help="first file (.vol or .f32)")
    parser.add_argument("b", help="second file, same type")
    parser.add_argument("--tolerance", type=float,
                        help="largest allowed absolute difference; sets the exit code")
    args = parser.parse_args()

    try:
        result = compare(args.a, args.b)
    except (ValueError, OSError) as error:
        print(f"Failed: {error}", file=sys.stderr)
        return 2

    print(f"voxels {result['voxels']}, differing {result['differing']}, "
          f"max abs {result['max_abs']:g}, rms {result['rms']:g}")
    if args.tolerance is not None and result["max_abs"] > args.tolerance:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
