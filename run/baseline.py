#!/usr/bin/env python3
"""Baseline of step 0.4: numbers that later steps are compared with.

Three parts, each written as a CSV file into run/baseline/<tag>/ :

  scores.csv    Phantom scores (tests/phantom/score.py) of the unfiltered
                scan, the original program over a range of sigma values, and
                the median and mean filters along z over a range of windows.

  accuracy.csv  The original program (one worker) against the brute-force
                filter on the same volume: largest and RMS difference of the
                float results on the 0...255 scale, and of the 16-bit results
                in counts. Volumes: the tracked 64^3 test scans, and larger
                local crops in data/crops/ if present.

  timing.csv    Run time (best of --repeats runs) and peak memory of the
                original program, multi-worker and one-worker, and of the
                brute-force filter, on the tracked 64^3 test scans. Peak
                memory is measured on Linux only.

<tag> names the machine in the output folder; it defaults to the operating
system ("windows" or "linux"). Rows are written as soon as they are ready, so
an interrupted run leaves usable files.

Needs "make reference" and "make bruteforce" to have been run. Intermediate
volumes go to build/baseline/ and are not kept in git.

Usage (from the repository root):
    python run/baseline.py                       all three parts
    python run/baseline.py --parts timing
    python run/baseline.py --parts accuracy --sizes 64 128

Needs only the Python standard library.
"""

import argparse
import csv
import glob
import os
import sys
import time

REPO_ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
for folder in ("tests/data", "tests/reference", "tests/phantom"):
    sys.path.insert(0, os.path.join(REPO_ROOT, *folder.split("/")))
import compare_volumes   # noqa: E402
import make_phantom      # noqa: E402
import run_reference     # noqa: E402
import score             # noqa: E402
import z_filter          # noqa: E402
from slice_images import find_volume  # noqa: E402

WORK_DIR = os.path.join(REPO_ROOT, "build", "baseline")
PHANTOM_DIR = os.path.join(REPO_ROOT, "build", "phantom")
TEST_SCANS = ("Sample1_18p15_5mu_c64", "Sample2_5mu_c64")

# Sigma values for the phantom scores: both sides of the limits 15 and 30 of
# the original parameter table, around the phantom's true sigma of 14.
SCORE_SIGMAS = (8, 10, 12, 14, 15, 16, 17, 20, 25, 30, 35)
Z_RADII = (1, 2, 4, 8, 16, 32)
# Sigma 10 uses patch radius 1, sigma 17 (the value used on the real samples)
# patch radius 2.
ACCURACY_SIGMAS = (10, 17)


class CsvLog:
    """A CSV file that is written and flushed row by row."""

    def __init__(self, path, columns):
        os.makedirs(os.path.dirname(path), exist_ok=True)
        self.handle = open(path, "w", newline="", encoding="ascii")
        self.writer = csv.writer(self.handle, lineterminator="\n")
        self.columns = columns
        self.writer.writerow(columns)
        self.handle.flush()

    def add(self, **values):
        """Append one row; missing columns stay empty, floats get 6 digits."""
        row = []
        for column in self.columns:
            value = values.get(column, "")
            row.append(f"{value:.6g}" if isinstance(value, float) else value)
        self.writer.writerow(row)
        self.handle.flush()
        print("  " + ", ".join(f"{c}={v}" for c, v in zip(self.columns, row) if v != ""), flush=True)


def ensure_phantom(preset):
    """Make the phantom of the given preset unless it exists; return its folder."""
    folder = os.path.join(PHANTOM_DIR, f"phantom_{preset}")
    if not os.path.isfile(os.path.join(folder, "phantom.txt")):
        args = argparse.Namespace(preset=preset, name=None, **{
            key: None for key in ("size", "fibre_diameter", "bore_diameter", "twist", "tilt",
                                  "glass_share", "void_diameters", "voids_per_diameter",
                                  "air", "epoxy", "glass", "sigma", "drift", "seed")})
        make_phantom.make_phantom(make_phantom.resolve_settings(args), PHANTOM_DIR)
    return folder


def score_columns(result):
    """Return the score values that go into scores.csv."""
    return {key: result[key] for key in ("psnr", "rmse_work", "wrong_all", "wrong_interior",
                                         "wrong_edge", "void_contrast")
            if result[key] is not None}


def part_scores(out_dir, presets):
    """Write scores.csv: phantom scores of the original program and the z filters."""
    log = CsvLog(os.path.join(out_dir, "scores.csv"),
                 ["phantom", "method", "parameter", "psnr", "rmse_work", "wrong_all",
                  "wrong_interior", "wrong_edge", "void_contrast", "seconds"])
    program = run_reference.default_program(one_worker=False)
    for preset in presets:
        phantom = ensure_phantom(preset)
        name = os.path.basename(phantom)
        log.add(phantom=name, method="unfiltered", **score_columns(score.score(phantom, phantom)))

        for sigma in SCORE_SIGMAS:
            run = run_reference.run_reference(phantom, sigma, WORK_DIR, program, keep_float=False)
            log.add(phantom=name, method="original", parameter=f"sigma {sigma}",
                    seconds=run["seconds"], **score_columns(score.score(phantom, run["vol"])))

        depth = find_volume(phantom)[1][2]
        for kind in ("median", "mean"):
            for radius in Z_RADII:
                if 2 * radius + 1 > depth:
                    continue
                start = time.perf_counter()
                out_vol = z_filter.z_filter(phantom, kind, radius, WORK_DIR)
                seconds = time.perf_counter() - start
                log.add(phantom=name, method=f"z {kind}", parameter=f"radius {radius}",
                        seconds=seconds, **score_columns(score.score(phantom, out_vol)))


def accuracy_volumes(sizes):
    """Return the scan folders for the accuracy part, smallest first.

    The tracked 64^3 test scans, then every crop in data/crops/. `sizes`, if
    given, keeps only cubes with one of these edge lengths.
    """
    folders = [os.path.join(REPO_ROOT, "tests", "data", scan) for scan in TEST_SCANS]
    folders += sorted(glob.glob(os.path.join(REPO_ROOT, "data", "crops", "*")))
    volumes = []
    for folder in folders:
        try:
            _, size = find_volume(folder)
        except ValueError:
            continue
        if sizes is None or size[0] in sizes:
            volumes.append((size[0] * size[1] * size[2], folder, size))
    return [(folder, size) for _, folder, size in sorted(volumes)]


def part_accuracy(out_dir, sizes):
    """Write accuracy.csv: original program (one worker) against brute force."""
    log = CsvLog(os.path.join(out_dir, "accuracy.csv"),
                 ["volume", "size", "sigma", "float_max_abs", "float_rms",
                  "counts_max_abs", "counts_rms", "counts_differing_share",
                  "original_seconds", "bruteforce_seconds"])
    original = run_reference.default_program(one_worker=True)
    bruteforce = run_reference.default_program(False, bruteforce=True)
    for folder, size in accuracy_volumes(sizes):
        name = os.path.basename(folder)
        for sigma in ACCURACY_SIGMAS:
            print(f"accuracy: {name}, sigma {sigma}", flush=True)
            ref = run_reference.run_reference(folder, sigma, WORK_DIR, original, keep_float=True)
            bf = run_reference.run_reference(folder, sigma, WORK_DIR, bruteforce,
                                             keep_float=True, tag="bf")
            as_float = compare_volumes.compare(os.path.splitext(ref["vol"])[0] + ".f32",
                                               os.path.splitext(bf["vol"])[0] + ".f32")
            as_counts = compare_volumes.compare(ref["vol"], bf["vol"])
            log.add(volume=name, size="x".join(str(n) for n in size), sigma=sigma,
                    float_max_abs=float(as_float["max_abs"]), float_rms=float(as_float["rms"]),
                    counts_max_abs=int(as_counts["max_abs"]), counts_rms=float(as_counts["rms"]),
                    counts_differing_share=as_counts["differing"] / as_counts["voxels"],
                    original_seconds=ref["seconds"], bruteforce_seconds=bf["seconds"])


def part_timing(out_dir, repeats):
    """Write timing.csv: run time and peak memory on the tracked test scans."""
    log = CsvLog(os.path.join(out_dir, "timing.csv"),
                 ["volume", "program", "sigma", "best_seconds", "all_seconds", "peak_mb"])
    programs = (("original", run_reference.default_program(one_worker=False), "ref"),
                ("original one worker", run_reference.default_program(one_worker=True), "ref"),
                ("brute force", run_reference.default_program(False, bruteforce=True), "bf"))
    for scan in TEST_SCANS:
        folder = os.path.join(REPO_ROOT, "tests", "data", scan)
        for sigma in ACCURACY_SIGMAS:
            for label, program, tag in programs:
                runs = [run_reference.run_reference(folder, sigma, WORK_DIR, program,
                                                    keep_float=False, tag=tag)
                        for _ in range(repeats)]
                seconds = [run["seconds"] for run in runs]
                peaks = [run["peak_mb"] for run in runs if run["peak_mb"] is not None]
                log.add(volume=scan, program=label, sigma=sigma, best_seconds=min(seconds),
                        all_seconds=" ".join(f"{s:.1f}" for s in seconds),
                        peak_mb=max(peaks) if peaks else "")


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--parts", nargs="+", choices=["scores", "accuracy", "timing"],
                        default=["timing", "scores", "accuracy"])
    parser.add_argument("--tag", default="windows" if os.name == "nt" else "linux",
                        help="name of the machine in the output folder")
    parser.add_argument("--presets", nargs="+", choices=sorted(make_phantom.PRESETS),
                        default=["coarse", "fine"], help="phantoms for the scores part")
    parser.add_argument("--sizes", type=int, nargs="+",
                        help="accuracy part: only cubes with these edge lengths")
    parser.add_argument("--repeats", type=int, default=3, help="timing part: runs per case")
    args = parser.parse_args()

    out_dir = os.path.join(REPO_ROOT, "run", "baseline", args.tag)
    try:
        if "timing" in args.parts:
            part_timing(out_dir, args.repeats)
        if "scores" in args.parts:
            part_scores(out_dir, args.presets)
        if "accuracy" in args.parts:
            part_accuracy(out_dir, args.sizes)
    except (ValueError, RuntimeError, OSError) as error:
        print(f"Failed: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
