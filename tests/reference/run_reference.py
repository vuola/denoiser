#!/usr/bin/env python3
"""Run the original, unchanged filter on a scan in the CT lab's format.

The original program (ORIGINAL_CODE/denoiserVersionA.cpp, built by
"make reference") reads only headerless 32-bit big-endian float files whose
size is given in the file name. This script is the adapter around it:

    1. read the scan's .vgi and .vol (16-bit unsigned, little-endian)
    2. write a temporary float file, value = count * 255 / 65535
    3. run the original program on it
    4. convert the result back to 16-bit and write it as .vol + .vgi
    5. delete the temporary float files

Output, in the folder given with --out (default build/reference/):

    <scan>_ref_s<sigma>.vol   filtered voxels, 16-bit unsigned little-endian
    <scan>_ref_s<sigma>.vgi   metadata
    <scan>_ref_s<sigma>.f32   only with --keep-float: the program's float
                              result before rounding, 32-bit little-endian
                              float on the 0...255 scale

With --bruteforce the brute-force filter (built by "make bruteforce") is run
instead of the original program, and the files are named <scan>_bf_s<sigma>.

The whole volume is held in memory, so the script is meant for the small test
scans, not for full scans.

Usage (from the repository root):
    python tests/reference/run_reference.py tests/data/Sample2_5mu_c64 --sigma 10

Needs only the Python standard library.
"""

import argparse
import array
import glob
import os
import shutil
import subprocess
import sys
import tempfile
import time

REPO_ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(REPO_ROOT, "tests", "data"))
from make_crops import read_vgi, write_vgi  # noqa: E402  (shared .vgi reader and writer)

# The original program expects values 0...255; 16-bit counts are 0...65535.
COUNTS_TO_WORK = 255.0 / 65535.0
WORK_TO_COUNTS = 65535.0 / 255.0


def default_program(one_worker, bruteforce=False):
    """Return the path of a program built by "make reference" or "make bruteforce".

    one_worker selects the repeatable one-worker reference program instead of
    the multi-threaded one; bruteforce selects the brute-force filter.
    """
    if bruteforce:
        name = "denoiser_bruteforce"
    else:
        name = "denoiser_reference_1w" if one_worker else "denoiser_reference"
    suffix = ".exe" if os.name == "nt" else ""
    return os.path.join(REPO_ROOT, "build", name + suffix)


def find_scan(scan_dir):
    """Return (vgi_path, vol_path) of the scan in scan_dir.

    Raises ValueError unless the folder holds exactly one .vgi file with a
    .vol file of the same name.
    """
    vgi_files = glob.glob(os.path.join(scan_dir, "*.vgi"))
    if len(vgi_files) != 1:
        raise ValueError(f"{scan_dir}: expected one .vgi file, found {len(vgi_files)}")
    vol_path = os.path.splitext(vgi_files[0])[0] + ".vol"
    if not os.path.isfile(vol_path):
        raise ValueError(f"{vol_path}: not found")
    return vgi_files[0], vol_path


def read_counts(vol_path, voxel_count):
    """Read a 16-bit unsigned little-endian .vol file.

    Returns an array('H') of voxel_count values in file order (x fastest).
    Raises ValueError if the file size is not voxel_count * 2 bytes.
    """
    if os.path.getsize(vol_path) != voxel_count * 2:
        raise ValueError(f"{vol_path}: size does not match the .vgi file")
    counts = array.array("H")
    with open(vol_path, "rb") as vol:
        counts.fromfile(vol, voxel_count)
    if sys.byteorder == "big":
        counts.byteswap()
    return counts


def write_big_endian_floats(path, values):
    """Write an array('f') to `path` as 32-bit big-endian floats."""
    out = array.array("f", values)
    if sys.byteorder == "little":
        out.byteswap()
    with open(path, "wb") as raw:
        out.tofile(raw)


def read_big_endian_floats(path, voxel_count):
    """Read voxel_count 32-bit big-endian floats from `path` as array('f').

    Raises ValueError if the file size is not voxel_count * 4 bytes.
    """
    if os.path.getsize(path) != voxel_count * 4:
        raise ValueError(f"{path}: the reference program wrote an incomplete file")
    values = array.array("f")
    with open(path, "rb") as raw:
        values.fromfile(raw, voxel_count)
    if sys.byteorder == "little":
        values.byteswap()
    return values


def to_counts(work_values):
    """Convert working-scale floats (0...255) to 16-bit counts.

    Values are rounded to the nearest count and limited to 0...65535.
    Returns (counts as array('H'), number of values that had to be limited).
    """
    counts = array.array("H", bytes(2 * len(work_values)))
    limited = 0
    for i, value in enumerate(work_values):
        count = int(value * WORK_TO_COUNTS + 0.5)
        if count < 0 or count > 65535:
            limited += 1
            count = min(max(count, 0), 65535)
        counts[i] = count
    return counts, limited


def run_reference(scan_dir, sigma, out_dir, program, keep_float, tag="ref"):
    """Filter the scan in scan_dir with the original program.

    Parameters:
        scan_dir    folder holding one .vgi and its .vol
        sigma       noise standard deviation on the 0...255 scale, > 0
        out_dir     folder for the result (created if missing)
        program     path of the reference program
        keep_float  also keep the unrounded float result as .f32
        tag         short text in the output name: <scan>_<tag>_s<sigma>

    Returns a dict with the result .vol path ("vol"), the run time of the
    program in seconds ("seconds"), the volume size ("size") and the number of
    result values outside 0...65535 ("limited").
    Raises ValueError on bad input and RuntimeError if the program fails.
    """
    if not os.path.isfile(program):
        raise ValueError(f"{program}: not found; build it with 'make reference' "
                         "or 'make bruteforce'")
    # Windows does not start a program given as a relative path with "/".
    program = os.path.abspath(program)

    vgi_path, vol_path = find_scan(scan_dir)
    info = read_vgi(vgi_path)
    if info["bits"] != 16 or "unsigned" not in info["datatype"].lower():
        raise ValueError(f"{vgi_path}: only 16-bit unsigned scans are supported")
    width, height, depth = info["size"]
    voxel_count = width * height * depth

    counts = read_counts(vol_path, voxel_count)
    work = array.array("f", (count * COUNTS_TO_WORK for count in counts))

    name = f"{os.path.basename(os.path.normpath(scan_dir))}_{tag}_s{sigma:g}"
    os.makedirs(out_dir, exist_ok=True)

    temp_dir = tempfile.mkdtemp(prefix="denoiser_ref_")
    try:
        # The original program reads the volume size from the file name.
        float_in = os.path.join(temp_dir, f"in_{width}x{height}x{depth}.raw")
        float_out = os.path.join(temp_dir, "out.raw")
        write_big_endian_floats(float_in, work)

        start = time.perf_counter()
        done = subprocess.run([program, float_in, float_out, f"{sigma:g}"],
                              stdout=subprocess.DEVNULL)
        seconds = time.perf_counter() - start
        if done.returncode != 0:
            raise RuntimeError(f"reference program failed with exit code {done.returncode}")

        result = read_big_endian_floats(float_out, voxel_count)
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)

    result_counts, limited = to_counts(result)
    if sys.byteorder == "big":
        result_counts.byteswap()
    out_vol = os.path.join(out_dir, name + ".vol")
    with open(out_vol, "wb") as vol:
        result_counts.tofile(vol)

    write_vgi(os.path.join(out_dir, name + ".vgi"), name + ".vol",
              (width, height, depth), info)

    if keep_float:
        if sys.byteorder == "big":
            result.byteswap()
        with open(os.path.join(out_dir, name + ".f32"), "wb") as raw:
            result.tofile(raw)

    return {"vol": out_vol, "seconds": seconds, "size": (width, height, depth),
            "limited": limited}


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("scan", help="scan folder holding a .vgi and a .vol file")
    parser.add_argument("--sigma", type=float, required=True,
                        help="noise standard deviation on the 0...255 scale")
    parser.add_argument("--out", default=os.path.join(REPO_ROOT, "build", "reference"),
                        help="output folder (default: build/reference)")
    parser.add_argument("--one-worker", action="store_true",
                        help="use the repeatable one-worker program "
                             "(build/denoiser_reference_1w)")
    parser.add_argument("--bruteforce", action="store_true",
                        help="use the brute-force filter (build/denoiser_bruteforce); "
                             "the output is named <scan>_bf_s<sigma>")
    parser.add_argument("--program",
                        help="reference program, if not one of the two in build/")
    parser.add_argument("--keep-float", action="store_true",
                        help="also keep the unrounded float result as .f32")
    args = parser.parse_args()

    if not args.sigma > 0:
        parser.error("--sigma must be greater than 0")

    try:
        program = args.program or default_program(args.one_worker, args.bruteforce)
        result = run_reference(args.scan, args.sigma, args.out, program, args.keep_float,
                               tag="bf" if args.bruteforce else "ref")
    except (ValueError, RuntimeError) as error:
        print(f"Failed: {error}", file=sys.stderr)
        return 1

    width, height, depth = result["size"]
    print(f"{os.path.relpath(result['vol'])}: {width}x{height}x{depth}, "
          f"sigma {args.sigma:g}, {result['seconds']:.1f} s, "
          f"{result['limited']} values outside 0...65535")
    return 0


if __name__ == "__main__":
    sys.exit(main())
