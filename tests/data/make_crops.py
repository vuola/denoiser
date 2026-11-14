#!/usr/bin/env python3
"""Cut small test scans from the full CT scans in data/RAW/.

The full scans are not tracked by git. This script cuts a small cube from each
of them and writes it below this script's folder as a miniature scan folder in
the CT lab's format, which IS tracked, so that tests can run from a fresh
clone:

    tests/data/<scan>_c<size>/<scan>_c<size>.vol   voxels
    tests/data/<scan>_c<size>/<scan>_c<size>.vgi   metadata

Input and output voxel format: headerless, unsigned integer, little-endian,
x fastest, then y, then z. Volume size and bits per voxel are read from the
scan's .vgi file. The crop is a plain byte copy: no value is converted.

The default region, origin (539, 539, 467) and size 64 (zero-based voxel
coordinates), is the centre of the 250^3 cube at (446, 446, 374) that the old
manual workflow cut from these scans.

Usage (from the repository root):
    python tests/data/make_crops.py
    python tests/data/make_crops.py --origin 100 100 100 --size 32

Needs only the Python standard library.
"""

import argparse
import glob
import os
import sys

DEFAULT_ORIGIN = (539, 539, 467)
DEFAULT_SIZE = 64


def read_vgi(vgi_path):
    """Return the volume description stored in a .vgi file.

    Only the first occurrence of each key is used, which is the
    [representation] block at the top of the file.

    Returns a dict with:
        size        (width, height, depth) in voxels
        bits        bits per voxel
        datatype    text, e.g. "unsigned integer"
        resolution  (rx, ry, rz) voxel edge lengths, or None if absent
        unit        text unit of the resolution, or None if absent
    Raises ValueError if size or bits per voxel are missing.
    """
    first = {}
    with open(vgi_path, "r", encoding="latin-1") as vgi:
        for line in vgi:
            if "=" not in line:
                continue
            key, value = (part.strip() for part in line.split("=", 1))
            first.setdefault(key.lower(), value)

    if "size" not in first or "bitsperelement" not in first:
        raise ValueError(f"{vgi_path}: size or bitsperelement missing")

    resolution = first.get("resolution")
    return {
        "size": tuple(int(n) for n in first["size"].split()),
        "bits": int(first["bitsperelement"]),
        "datatype": first.get("datatype", "unsigned integer"),
        "resolution": tuple(float(r) for r in resolution.split()) if resolution else None,
        "unit": first.get("unit"),
    }


def write_vgi(vgi_path, vol_name, size, info):
    """Write a minimal .vgi file for a cube of edge `size` voxels.

    `info` is the dict returned by read_vgi() for the source scan; data type,
    bits per voxel and voxel size are copied from it. No path or name from the
    source file is copied.
    """
    max_value = 2 ** info["bits"] - 1
    lines = [
        "{volume1}",
        "[representation]",
        f"size = {size} {size} {size}",
        f"bitsperelement = {info['bits']}",
        f"datatype = {info['datatype']}",
        f"datarange = 0 {max_value}",
        "[file1]",
        "FileFormat = raw",
        f"Size = {size} {size} {size}",
        f"Name = {vol_name}",
        f"bitsperelement = {info['bits']}",
        f"datatype = {info['datatype']}",
        f"datarange = 0 {max_value}",
    ]
    if info["resolution"]:
        lines += [
            "{volumeprimitive9}",
            "[geometry]",
            "resolution = " + " ".join(f"{r:.6f}" for r in info["resolution"]),
            f"unit = {info['unit'] or 'mm'}",
        ]
    with open(vgi_path, "w", encoding="ascii", newline="\n") as vgi:
        vgi.write("\n".join(lines) + "\n")


def crop_scan(scan_dir, dst_root, origin, size):
    """Cut a size^3 cube starting at `origin` from the scan in scan_dir.

    Returns the output folder. Raises ValueError if the scan folder does not
    hold exactly one .vgi with a matching .vol, the .vol size disagrees with
    the .vgi, or the region does not fit inside the volume.
    """
    vgi_files = glob.glob(os.path.join(scan_dir, "*.vgi"))
    if len(vgi_files) != 1:
        raise ValueError(f"{scan_dir}: expected one .vgi file, found {len(vgi_files)}")
    vol_path = os.path.splitext(vgi_files[0])[0] + ".vol"
    if not os.path.isfile(vol_path):
        raise ValueError(f"{vol_path}: not found")

    info = read_vgi(vgi_files[0])
    width, height, depth = info["size"]
    bytes_per_voxel = info["bits"] // 8

    expected = width * height * depth * bytes_per_voxel
    actual = os.path.getsize(vol_path)
    if actual != expected:
        raise ValueError(f"{vol_path}: size {actual} bytes, expected {expected}")

    x0, y0, z0 = origin
    if min(origin) < 0 or x0 + size > width or y0 + size > height or z0 + size > depth:
        raise ValueError(f"{scan_dir}: region does not fit in {width}x{height}x{depth}")

    name = f"{os.path.basename(os.path.normpath(scan_dir))}_c{size}"
    dst_dir = os.path.join(dst_root, name)
    os.makedirs(dst_dir, exist_ok=True)

    row_bytes = size * bytes_per_voxel
    with open(vol_path, "rb") as src, open(os.path.join(dst_dir, name + ".vol"), "wb") as dst:
        for z in range(z0, z0 + size):
            for y in range(y0, y0 + size):
                src.seek(((z * height + y) * width + x0) * bytes_per_voxel)
                dst.write(src.read(row_bytes))

    write_vgi(os.path.join(dst_dir, name + ".vgi"), name + ".vol", size, info)
    return dst_dir


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    default_src = os.path.normpath(os.path.join(here, "..", "..", "data", "RAW"))

    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--src", default=default_src, help="folder holding the scan folders")
    parser.add_argument("--origin", type=int, nargs=3, default=DEFAULT_ORIGIN,
                        metavar=("X", "Y", "Z"), help="first voxel of the cube, zero-based")
    parser.add_argument("--size", type=int, default=DEFAULT_SIZE, help="cube edge in voxels")
    args = parser.parse_args()

    scans = sorted(d for d in glob.glob(os.path.join(args.src, "*")) if os.path.isdir(d))
    if not scans:
        print(f"No scan folders in {args.src}", file=sys.stderr)
        return 1

    for scan_dir in scans:
        try:
            dst_dir = crop_scan(scan_dir, here, tuple(args.origin), args.size)
        except ValueError as error:
            print(f"Failed: {error}", file=sys.stderr)
            return 1
        print(f"{scan_dir} -> {os.path.relpath(dst_dir)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
