#!/usr/bin/env python3
"""Score a volume against a phantom: noise left, and how well materials separate.

The volume is usually a filter result for the phantom's noisy scan; giving the
noisy scan itself shows the starting point. The phantom's noise-free volume,
labels and grey levels (written by tests/phantom/make_phantom.py) are the
truth.

What is measured (README requirement R1):

  Error against the noise-free volume
      rmse      root-mean-square error, on the 0...255 scale and in counts
      psnr      20 log10(255 / rmse), in dB, rmse on the 0...255 scale

  Classification
      Every voxel is classified from its value alone, with two fixed
      thresholds halfway between the phantom's grey levels:
          value < (air + epoxy) / 2      -> air
          value < (epoxy + glass) / 2    -> epoxy
          otherwise                      -> glass
      and compared with the true material. Void voxels count as air.
      wrong_all       share of all voxels classified wrongly
      wrong_interior  the same among voxels that hold one material only
      wrong_edge      the same among edge voxels (a mixture of two materials)
      Shares per true material are printed as a table.

  Voids
      void_contrast   how much of the step from epoxy down to air is left in
                      the void voxels: (epoxy - mean void value) / (epoxy - air).
                      1 = voids as dark as air, 0 = voids filled in to epoxy.

If the phantom was made with a brightness drift, the fixed thresholds are part
of what is tested: drift moves voxels across them.

Usage (from the repository root):
    python tests/phantom/score.py build/phantom/phantom_fine build/phantom/phantom_fine
    python tests/phantom/score.py build/phantom/phantom_fine build/reference/phantom_fine_ref_s14.vol

The first argument is the phantom's noisy scan folder, the second the volume
to score (scan folder or .vol file with its .vgi). Volumes are held in
memory, so the script is meant for phantoms, not for full scans.

Needs only the Python standard library.
"""

import argparse
import array
import math
import os
import sys

REPO_ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(REPO_ROOT, "tests", "data"))
sys.path.insert(0, os.path.join(REPO_ROOT, "tests", "reference"))
from slice_images import find_volume  # noqa: E402  (shared scan-folder check)

COUNTS_TO_WORK = 255.0 / 65535.0
CLASS_NAMES = ("air", "epoxy", "glass")
LABEL_VOID = 3
LABEL_EDGE = 16


def read_record(phantom_dir):
    """Return the key = value pairs of the phantom's phantom.txt as a dict of text."""
    record = {}
    with open(os.path.join(phantom_dir, "phantom.txt"), "r", encoding="ascii") as handle:
        for line in handle:
            if "=" in line:
                key, value = (part.strip() for part in line.split("=", 1))
                record[key] = value
    return record


def read_counts(vol_path, voxel_count):
    """Read a 16-bit unsigned little-endian .vol file into an array('H')."""
    values = array.array("H")
    with open(vol_path, "rb") as handle:
        values.fromfile(handle, voxel_count)
    if sys.byteorder == "big":
        values.byteswap()
    return values


def score(phantom_dir, volume_path):
    """Score the volume at `volume_path` against the phantom in `phantom_dir`.

    Parameters:
        phantom_dir  the phantom's noisy scan folder; the noise-free volume
                     and the labels are expected next to it as <name>_clean
                     and <name>_labels
        volume_path  scan folder or .vol file of the volume to score

    Returns a dict with rmse_work, rmse_counts, psnr, wrong_all,
    wrong_interior, wrong_edge, void_contrast (None without voids), voxels,
    and "confusion": a 3 x 3 list, confusion[true][found] in voxels.
    Raises ValueError if files are missing or sizes disagree.
    """
    phantom_dir = os.path.normpath(phantom_dir)
    record = read_record(phantom_dir)
    air, epoxy, glass = (int(record[f"level_{name}"]) for name in CLASS_NAMES)
    air_epoxy = (air + epoxy) / 2.0
    epoxy_glass = (epoxy + glass) / 2.0

    clean_path, size = find_volume(phantom_dir + "_clean")
    volume_vol, volume_size = find_volume(volume_path)
    if size != volume_size:
        raise ValueError(f"the volume is {volume_size}, the phantom {size}")
    voxel_count = size[0] * size[1] * size[2]

    label_dir = phantom_dir + "_labels"
    label_path = os.path.join(label_dir, os.path.basename(label_dir) + ".vol")
    if not os.path.isfile(label_path) or os.path.getsize(label_path) != voxel_count:
        raise ValueError(f"{label_path}: missing or of the wrong size")
    with open(label_path, "rb") as handle:
        labels = handle.read()

    clean = read_counts(clean_path, voxel_count)
    values = read_counts(volume_vol, voxel_count)

    sum_squares = 0
    confusion = [[0] * 3 for _ in range(3)]
    wrong = {"interior": 0, "edge": 0}
    total = {"interior": 0, "edge": 0}
    void_sum = 0
    void_voxels = 0
    for value, truth, label in zip(values, clean, labels):
        difference = value - truth
        sum_squares += difference * difference

        material = label & 15
        if material == LABEL_VOID:
            void_sum += value
            void_voxels += 1
            material = 0
        found = 0 if value < air_epoxy else (1 if value < epoxy_glass else 2)
        confusion[material][found] += 1
        where = "edge" if label & LABEL_EDGE else "interior"
        total[where] += 1
        if found != material:
            wrong[where] += 1

    rmse_counts = math.sqrt(sum_squares / voxel_count)
    rmse_work = rmse_counts * COUNTS_TO_WORK
    return {
        "voxels": voxel_count,
        "rmse_counts": rmse_counts,
        "rmse_work": rmse_work,
        "psnr": 20.0 * math.log10(255.0 / rmse_work) if rmse_work > 0 else math.inf,
        "wrong_all": (wrong["interior"] + wrong["edge"]) / voxel_count,
        "wrong_interior": wrong["interior"] / max(1, total["interior"]),
        "wrong_edge": wrong["edge"] / max(1, total["edge"]),
        "void_contrast": ((epoxy - void_sum / void_voxels) / (epoxy - air)
                          if void_voxels and epoxy != air else None),
        "void_voxels": void_voxels,
        "confusion": confusion,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("phantom", help="the phantom's noisy scan folder")
    parser.add_argument("volume", help="volume to score: scan folder or .vol file")
    args = parser.parse_args()

    try:
        result = score(args.phantom, args.volume)
    except (ValueError, OSError, KeyError) as error:
        print(f"Failed: {error}", file=sys.stderr)
        return 1

    print(f"rmse            {result['rmse_work']:.3f} on the 0...255 scale "
          f"({result['rmse_counts']:.0f} counts)")
    print(f"psnr            {result['psnr']:.2f} dB")
    print(f"wrong_all       {100 * result['wrong_all']:.3f} %")
    print(f"wrong_interior  {100 * result['wrong_interior']:.3f} %")
    print(f"wrong_edge      {100 * result['wrong_edge']:.3f} %")
    if result["void_contrast"] is not None:
        print(f"void_contrast   {result['void_contrast']:.3f} "
              f"({result['void_voxels']} void voxels)")
    print("true material -> classified as (% of the true material)")
    print("            " + "".join(f"{name:>9}" for name in CLASS_NAMES))
    for name, row in zip(CLASS_NAMES, result["confusion"]):
        row_total = max(1, sum(row))
        print(f"    {name:<8}" + "".join(f"{100 * n / row_total:9.2f}" for n in row))
    return 0


if __name__ == "__main__":
    sys.exit(main())
