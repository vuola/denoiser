#!/usr/bin/env python3
"""Cut small test crops from the full CT volumes in data/.

The full volumes (data/*_WxHxD.raw) are not tracked by git. This script cuts a
centred cube from each of them and writes it next to this script, where it IS
tracked, so that tests can run from a fresh clone.

Format of both input and output: headerless raw, 32-bit big-endian float,
x fastest, then y, then z. Dimensions are taken from the end of the file name
(name_WxHxD.raw), the same convention as the original denoiser.

The crop is a plain byte copy: no value is converted, so the crop is
bit-identical to the corresponding region of the full volume.

Usage (from the repository root):
    python tests/data/make_crops.py            # 64-voxel cube
    python tests/data/make_crops.py --size 32  # other cube edge length

Needs only the Python standard library.
"""

import argparse
import glob
import os
import re
import sys

BYTES_PER_VOXEL = 4  # 32-bit float
NAME_PATTERN = re.compile(r"^(?P<stem>.*?)(?P<w>\d+)x(?P<h>\d+)x(?P<d>\d+)\.raw$")


def crop_centre(src_path, dst_dir, size):
    """Write the centred size^3 cube of src_path into dst_dir.

    Returns the output path. Raises ValueError if the file name carries no
    dimensions, the file size does not match them, or the volume is smaller
    than the requested cube.
    """
    name = os.path.basename(src_path)
    match = NAME_PATTERN.match(name)
    if not match:
        raise ValueError(f"{name}: no WxHxD dimensions at the end of the name")
    width, height, depth = (int(match[k]) for k in ("w", "h", "d"))

    expected = width * height * depth * BYTES_PER_VOXEL
    actual = os.path.getsize(src_path)
    if actual != expected:
        raise ValueError(f"{name}: size {actual} bytes, expected {expected}")
    if min(width, height, depth) < size:
        raise ValueError(f"{name}: volume smaller than {size}^3")

    x0 = (width - size) // 2
    y0 = (height - size) // 2
    z0 = (depth - size) // 2

    dst_path = os.path.join(dst_dir, f"{match['stem']}{size}x{size}x{size}.raw")
    row_bytes = size * BYTES_PER_VOXEL
    with open(src_path, "rb") as src, open(dst_path, "wb") as dst:
        for z in range(z0, z0 + size):
            for y in range(y0, y0 + size):
                src.seek(((z * height + y) * width + x0) * BYTES_PER_VOXEL)
                dst.write(src.read(row_bytes))
    return dst_path


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    default_src = os.path.normpath(os.path.join(here, "..", "..", "data"))

    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--src", default=default_src, help="folder with the full volumes")
    parser.add_argument("--size", type=int, default=64, help="cube edge length in voxels")
    args = parser.parse_args()

    sources = sorted(glob.glob(os.path.join(args.src, "*.raw")))
    if not sources:
        print(f"No .raw files in {args.src}", file=sys.stderr)
        return 1

    for src_path in sources:
        try:
            dst_path = crop_centre(src_path, here, args.size)
        except ValueError as error:
            print(f"Skipped {error}", file=sys.stderr)
            return 1
        print(f"{os.path.basename(src_path)} -> {os.path.relpath(dst_path)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
