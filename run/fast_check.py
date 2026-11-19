#!/usr/bin/env python3
"""Check of step 1.6: the fast filter against the brute-force filter.

For every volume and every setting, the brute-force filter (the reference),
the fast filter in double precision and the fast filter in single precision
are run, and the two fast results are compared with the reference. One row
per volume, setting and computation goes into run/fast/<tag>/check.csv, as
soon as it is ready:

    float_max_abs, float_rms   difference of the float results on the 0...255
                               scale (empty for the reference itself)
    counts_max_abs             largest difference after rounding to 16 bits
    seconds, peak_mb           run time and peak memory (Linux only)

Volumes: the tracked 64^3 test scans, the 232^3 crop in data/crops/ if
present, and the phantoms fine_cn and coarse_cn (made if missing).

Settings:
    table 10   the original table at sigma 10 (patch 1, search 10)
    table 17   the original table at sigma 17 (patch 2, search 10)
    rods       the defaults for rods of README section 12.7: sigma 14,
               patch 0 0 8, search 3 3 30, h = 0.4 * sigma

Needs "make bruteforce". Intermediate volumes go to build/fast/.

Usage (from the repository root):
    python run/fast_check.py
    python run/fast_check.py --volumes Sample2_5mu_c64 --settings rods

Needs only the Python standard library.
"""

import argparse
import glob
import os
import sys

REPO_ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
sys.path.insert(0, os.path.join(REPO_ROOT, "run"))
import baseline          # noqa: E402  (also sets the path to the other scripts)
import compare_volumes   # noqa: E402
import params            # noqa: E402
import run_reference     # noqa: E402

WORK_DIR = os.path.join(REPO_ROOT, "build", "fast")

SETTINGS = {
    "table 10": (10.0, []),
    "table 17": (17.0, []),
    "rods": (14.0, ["--patch", "0", "0", "8", "--search", "3", "3", "30", "--h-factor", "0.4"]),
}
COMPUTATIONS = (("brute force", "bf", []),
                ("fast double", "fd", ["--fast", "--double"]),
                ("fast single", "fs", ["--fast", "--single"]))
COLUMNS = ["volume", "size", "setting", "computation", "float_max_abs", "float_rms",
           "counts_max_abs", "seconds", "peak_mb"]


def volumes():
    """Return the volumes as (name, folder), smallest first."""
    found = [(scan, os.path.join(REPO_ROOT, "tests", "data", scan)) for scan in baseline.TEST_SCANS]
    found += [(os.path.basename(f), f)
              for f in sorted(glob.glob(os.path.join(REPO_ROOT, "data", "crops", "*_c232")))]
    found += [(f"phantom_{key}", None) for key in ("coarse_cn", "fine_cn")]
    return found


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--volumes", nargs="+", help="names of the volumes to use (default: all)")
    parser.add_argument("--settings", nargs="+", choices=sorted(SETTINGS), default=list(SETTINGS))
    parser.add_argument("--tag", default="windows" if os.name == "nt" else "linux",
                        help="name of the machine in the output folder")
    parser.add_argument("--name", default="check", help="name of the CSV file, without .csv")
    args = parser.parse_args()

    log = baseline.CsvLog(os.path.join(REPO_ROOT, "run", "fast", args.tag, args.name + ".csv"),
                          COLUMNS)
    program = run_reference.default_program(False, bruteforce=True)
    try:
        for name, folder in volumes():
            if args.volumes and name not in args.volumes:
                continue
            if folder is None:
                folder = params.ensure_phantom(name[len("phantom_"):])
            for setting in args.settings:
                sigma, options = SETTINGS[setting]
                runs = {}
                for label, tag, extra in COMPUTATIONS:
                    run = run_reference.run_reference(
                        folder, sigma, WORK_DIR, program, keep_float=True,
                        tag=f"{tag}_{setting.replace(' ', '')}", options=options + extra)
                    runs[label] = run
                    row = {"volume": name, "size": "x".join(str(n) for n in run["size"]),
                           "setting": setting, "computation": label, "seconds": run["seconds"],
                           "peak_mb": run["peak_mb"] if run["peak_mb"] is not None else ""}
                    if label != "brute force":
                        reference = runs["brute force"]["vol"]
                        as_float = compare_volumes.compare(
                            os.path.splitext(run["vol"])[0] + ".f32",
                            os.path.splitext(reference)[0] + ".f32")
                        as_counts = compare_volumes.compare(run["vol"], reference)
                        row.update(float_max_abs=float(as_float["max_abs"]),
                                   float_rms=float(as_float["rms"]),
                                   counts_max_abs=int(as_counts["max_abs"]))
                    log.add(**row)
                for run in runs.values():
                    stem = os.path.splitext(run["vol"])[0]
                    for extension in (".vol", ".vgi", ".f32"):
                        os.remove(stem + extension)
    except (ValueError, RuntimeError, OSError) as error:
        print(f"Failed: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
