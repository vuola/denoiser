#!/usr/bin/env python3
"""Write PNG images of slices for visual inspection of a filtering result.

Each image has three panels of the same slice, left to right:

    input | filtered | difference (filtered - input)

The input and filtered panels share one grey window, so they can be compared
directly. In the difference panel mid-grey means "no change", brighter means
the filter raised the value and darker means it lowered it; the difference is
shown on the same grey scale as the other panels, times --gain. A difference panel that looks like pure noise means the
filter removed only noise; visible outlines of structures mean it also changed
the structures.

Both volumes are given as a scan folder (holding one .vgi and its .vol) or as
a .vol file with its .vgi next to it. They must have the same size and be
16-bit unsigned. Only the requested slices are read, so full scans work too.

Output, in the folder given with --out (default build/slices/):

    <filtered name>_<axis><slice>.png      e.g. Sample2_5mu_c64_ref_s10_z32.png

Slice planes:  z-slice: x to the right, y down
               y-slice: x to the right, z down
               x-slice: y to the right, z down

Usage (from the repository root):
    python tests/reference/slice_images.py tests/data/Sample2_5mu_c64 \n        build/reference/Sample2_5mu_c64_ref_s10.vol
    python tests/reference/slice_images.py before.vol after.vol --axis z --slice 500

Needs only the Python standard library.
"""

import argparse
import array
import glob
import os
import struct
import sys
import zlib

REPO_ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(REPO_ROOT, "tests", "data"))
from make_crops import read_vgi  # noqa: E402  (shared .vgi reader)

GAP_PIXELS = 4      # white gap between panels
MIN_PANEL = 256     # small slices are enlarged to at least this many pixels


def find_volume(path):
    """Return (vol_path, (width, height, depth)) for a scan folder or .vol file.

    Raises ValueError if the .vgi or .vol is missing, a folder holds more than
    one .vgi, the volume is not 16-bit unsigned, or the .vol size disagrees
    with the .vgi.
    """
    if os.path.isdir(path):
        vgi_files = glob.glob(os.path.join(path, "*.vgi"))
        if len(vgi_files) != 1:
            raise ValueError(f"{path}: expected one .vgi file, found {len(vgi_files)}; "
                             "give the .vol file instead")
        vgi_path = vgi_files[0]
    else:
        vgi_path = os.path.splitext(path)[0] + ".vgi"
    vol_path = os.path.splitext(vgi_path)[0] + ".vol"
    for needed in (vgi_path, vol_path):
        if not os.path.isfile(needed):
            raise ValueError(f"{needed}: not found")

    info = read_vgi(vgi_path)
    if info["bits"] != 16 or "unsigned" not in info["datatype"].lower():
        raise ValueError(f"{vgi_path}: only 16-bit unsigned volumes are supported")
    width, height, depth = info["size"]
    if os.path.getsize(vol_path) != width * height * depth * 2:
        raise ValueError(f"{vol_path}: size does not match the .vgi file")
    return vol_path, (width, height, depth)


def read_slice(vol_path, size, axis, index):
    """Read one slice of a 16-bit unsigned little-endian volume.

    Parameters:
        vol_path  the .vol file
        size      (width, height, depth) of the volume
        axis      "x", "y" or "z": the axis the slice is perpendicular to
        index     zero-based position of the slice along that axis

    Returns (values, columns, rows): values is an array('H') of
    columns * rows counts, row by row, in the plane orientation described in
    the module text.
    """
    width, height, depth = size
    values = array.array("H")
    with open(vol_path, "rb") as vol:
        if axis == "z":
            vol.seek(index * height * width * 2)
            values.fromfile(vol, width * height)
            columns, rows = width, height
        elif axis == "y":
            for z in range(depth):
                vol.seek(((z * height + index) * width) * 2)
                values.fromfile(vol, width)
            columns, rows = width, depth
        else:
            for z in range(depth):
                vol.seek(z * height * width * 2)
                plane = array.array("H")
                plane.fromfile(vol, width * height)
                values.extend(plane[index::width])
            columns, rows = height, depth
    if sys.byteorder == "big":
        values.byteswap()
    return values, columns, rows


def percentile(sorted_values, percent):
    """Return the value at `percent` (0...100) of an already sorted sequence."""
    position = int(round((len(sorted_values) - 1) * percent / 100.0))
    return sorted_values[position]


def to_grey(values, low, high):
    """Map counts to display grey 0...255: `low` -> 0, `high` -> 255, limited."""
    scale = 255.0 / (high - low)
    return bytes(min(255, max(0, int((value - low) * scale + 0.5))) for value in values)


def difference_grey(before, after, low, high, gain):
    """Map (after - before) to display grey with 128 = no change.

    The difference uses the same counts-to-grey scale as the window
    low...high, multiplied by `gain`.
    """
    scale = gain * 255.0 / (high - low)
    return bytes(min(255, max(0, int(128.5 + (a - b) * scale)))
                 for b, a in zip(before, after))


def write_png(path, pixels, columns, rows):
    """Write 8-bit grey pixels (bytes, row by row) as a PNG file."""
    def chunk(kind, data):
        body = kind + data
        return struct.pack(">I", len(data)) + body + struct.pack(">I", zlib.crc32(body))

    raw = b"".join(b"\x00" + pixels[row * columns:(row + 1) * columns] for row in range(rows))
    with open(path, "wb") as png:
        png.write(b"\x89PNG\r\n\x1a\n")
        png.write(chunk(b"IHDR", struct.pack(">IIBBBBB", columns, rows, 8, 0, 0, 0, 0)))
        png.write(chunk(b"IDAT", zlib.compress(raw, 6)))
        png.write(chunk(b"IEND", b""))


def join_panels(panels, columns, rows, zoom):
    """Place equally sized grey panels side by side, enlarged `zoom` times.

    Returns (pixels, total_columns, total_rows) with a white gap between
    panels. Enlargement repeats pixels, so no values are invented.
    """
    gap = b"\xff" * GAP_PIXELS
    lines = []
    for row in range(rows):
        parts = []
        for panel in panels:
            line = panel[row * columns:(row + 1) * columns]
            if zoom > 1:
                line = b"".join(bytes([value]) * zoom for value in line)
            parts.append(line)
        lines.extend([gap.join(parts)] * zoom)
    total_columns = len(panels) * columns * zoom + (len(panels) - 1) * GAP_PIXELS
    return b"".join(lines), total_columns, rows * zoom


def make_image(before_vol, after_vol, size, axis, index, out_path, window, gain):
    """Write the three-panel image of one slice.

    Parameters:
        before_vol, after_vol  .vol files of the input and the filtered volume
        size      (width, height, depth) of both volumes
        axis      "x", "y" or "z"
        index     zero-based slice position along `axis`
        out_path  PNG file to write
        window    (low, high) counts shown as black and white, or None to use
                  the 0.5 and 99.5 percentiles of the input slice
        gain      amplification of the difference panel

    Returns a dict with the window used ("window") and the root-mean-square
    difference of the slice in counts ("rms").
    """
    before, columns, rows = read_slice(before_vol, size, axis, index)
    after, _, _ = read_slice(after_vol, size, axis, index)

    if window is None:
        ordered = sorted(before)
        window = (percentile(ordered, 0.5), percentile(ordered, 99.5))
    low, high = window
    if high <= low:
        high = low + 1

    panels = [to_grey(before, low, high), to_grey(after, low, high),
              difference_grey(before, after, low, high, gain)]
    zoom = max(1, MIN_PANEL // max(columns, rows))
    pixels, total_columns, total_rows = join_panels(panels, columns, rows, zoom)
    write_png(out_path, pixels, total_columns, total_rows)

    rms = (sum((a - b) ** 2 for b, a in zip(before, after)) / len(before)) ** 0.5
    return {"window": (low, high), "rms": rms}


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("before", help="input volume: scan folder or .vol file")
    parser.add_argument("after", help="filtered volume: scan folder or .vol file")
    parser.add_argument("--axis", choices=["x", "y", "z", "all"], default="all",
                        help="axis the slice is perpendicular to (default: all three)")
    parser.add_argument("--slice", type=int,
                        help="zero-based slice position (default: the middle)")
    parser.add_argument("--window", type=int, nargs=2, metavar=("LOW", "HIGH"),
                        help="counts shown as black and white (default: 0.5 and 99.5 "
                             "percentiles of the input slice)")
    parser.add_argument("--gain", type=float, default=1.0,
                        help="amplification of the difference panel (default: 1)")
    parser.add_argument("--out", default=os.path.join(REPO_ROOT, "build", "slices"),
                        help="output folder (default: build/slices)")
    args = parser.parse_args()

    try:
        before_vol, size = find_volume(args.before)
        after_vol, after_size = find_volume(args.after)
        if size != after_size:
            raise ValueError(f"the volumes have different sizes: {size} and {after_size}")

        os.makedirs(args.out, exist_ok=True)
        name = os.path.splitext(os.path.basename(after_vol))[0]
        for axis in ("x", "y", "z") if args.axis == "all" else (args.axis,):
            length = size["xyz".index(axis)]
            index = length // 2 if args.slice is None else args.slice
            if not 0 <= index < length:
                raise ValueError(f"slice {index} is outside 0...{length - 1} on axis {axis}")
            out_path = os.path.join(args.out, f"{name}_{axis}{index}.png")
            result = make_image(before_vol, after_vol, size, axis, index, out_path,
                                args.window, args.gain)
            low, high = result["window"]
            print(f"{os.path.relpath(out_path)}: window {low}...{high} counts, "
                  f"rms difference {result['rms']:.0f} counts")
    except (ValueError, OSError) as error:
        print(f"Failed: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
