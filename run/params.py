#!/usr/bin/env python3
"""Experiments of step 1.5: which filter settings classify the phantom best.

The brute-force filter is run on the phantoms with different settings, and
every result is scored (tests/phantom/score.py). Rows are written into
run/params/<tag>/<name>.csv as soon as they are ready.

Phantoms (made in build/phantom/ if missing), true noise sigma 14:

    coarse, fine          independent noise (as in chapter 7 of the README)
    coarse_cn, fine_cn    noise shared between neighbours, as in the real
                          scans (README chapter 11)

Sets of settings (--sets); unless stated, sigma 14, h = 0.4 * sigma, patch
radius 1 1 1, search radius 10 10 10, summed distance, rational kernel:

    original  the original program over a range of sigma, and the mean and
              median along z, as in the baseline of step 0.4
    A         patch radius 1 and 2 at sigma 14 and 16
    A2        patch radius 1 and 2 at sigma 12 and 18
    B         search windows that are longer along z than across
    C         patches that are longer along z than across, search 5 5 20
    D         distance as mean over the patch, exponential kernel
    E, F      larger h for the patches 1 1 4 and 0 0 4, search 5 5 20
    G         longer patches and narrower, longer search windows
    H         the patch 0 0 4 at other sigma values; the IPOL variant with
              long patches
    I         column patches of 25 and 33 voxels; other windows and sigma
              values for the column of 17

Needs "make bruteforce" (and "make reference" for the set "original").
Intermediate volumes go to build/params/ and are not kept in git.

Usage (from the repository root):
    python run/params.py --phantoms coarse_cn --sets A B --name round1
    python run/params.py --phantoms fine_cn coarse_cn --sets original --name original

Needs only the Python standard library.
"""

import argparse
import os
import sys
import time

REPO_ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
sys.path.insert(0, os.path.join(REPO_ROOT, "run"))
import baseline          # noqa: E402  (also sets the path to the other scripts)
import make_phantom      # noqa: E402
import run_reference     # noqa: E402
import score             # noqa: E402
import z_filter          # noqa: E402
from slice_images import find_volume  # noqa: E402

WORK_DIR = os.path.join(REPO_ROOT, "build", "params")
PHANTOM_DIR = os.path.join(REPO_ROOT, "build", "phantom")
PHANTOMS = {"coarse": ("coarse", False), "coarse_cn": ("coarse", True),
            "fine": ("fine", False), "fine_cn": ("fine", True)}

DEFAULT = {"sigma": 14.0, "patch": (1, 1, 1), "search": (10, 10, 10), "h_factor": 0.4,
           "normalise": False, "kernel": "rational"}


def variant(**changes):
    """Return the default settings with the given changes."""
    settings = dict(DEFAULT)
    settings.update(changes)
    return settings


P1, P2 = (1, 1, 1), (2, 2, 2)
LONG = (5, 5, 20)
COLUMN = (0, 0, 4)   # a patch of one column of nine voxels along z
SETS = {
    "A": [variant(patch=p, sigma=s) for s in (14.0, 16.0) for p in (P1, P2)],
    "A2": [variant(patch=p, sigma=s) for s in (12.0, 18.0) for p in (P1, P2)],
    "B": [variant(search=LONG), variant(search=(3, 3, 20)), variant(search=(5, 5, 10)),
          variant(patch=P2, search=LONG), variant(patch=P2, search=(3, 3, 20))],
    "C": [variant(patch=(1, 1, 2), search=LONG), variant(patch=(1, 1, 4), search=LONG),
          variant(patch=(2, 2, 4), search=LONG), variant(patch=(0, 0, 4), search=LONG)],
    "D": [variant(normalise=True, kernel="exponential"),
          variant(normalise=True, kernel="exponential", patch=P2),
          variant(normalise=True),
          variant(kernel="exponential")],
    # Second round, chosen after the first: see README chapter 12.
    "E": [variant(patch=(1, 1, 4), search=LONG, h_factor=f) for f in (0.8, 1.6, 3.2)],
    "F": [variant(patch=COLUMN, search=LONG, h_factor=f) for f in (0.8, 1.6, 3.2)],
    "G": [variant(patch=(0, 0, 8), search=LONG), variant(patch=(1, 1, 8), search=LONG),
          variant(patch=(0, 0, 8), search=(3, 3, 30)), variant(patch=COLUMN, search=(3, 3, 20)),
          variant(patch=COLUMN, search=(1, 1, 20)), variant(patch=(0, 0, 8), search=(1, 1, 30))],
    "H": [variant(patch=COLUMN, search=LONG, sigma=s) for s in (12.0, 16.0, 18.0)]
         + [variant(patch=p, search=LONG, normalise=True, kernel="exponential")
            for p in ((1, 1, 4), COLUMN)],
    "I": [variant(patch=(0, 0, 12), search=(3, 3, 30)), variant(patch=(0, 0, 16), search=(3, 3, 30)),
          variant(patch=(0, 0, 8), search=(3, 3, 40)), variant(patch=(0, 0, 8), search=(2, 2, 30)),
          variant(patch=(0, 0, 8), search=(3, 3, 30), sigma=16.0),
          variant(patch=(0, 0, 8), search=(3, 3, 30), sigma=12.0)],
}

COLUMNS = ["phantom", "set", "method", "sigma", "patch", "search", "h_factor", "distance",
           "kernel", "psnr", "rmse_work", "wrong_all", "wrong_interior", "wrong_edge",
           "void_contrast", "seconds"]


def ensure_phantom(key):
    """Make the phantom named by `key` unless it exists; return its folder."""
    preset, correlated = PHANTOMS[key]
    folder = os.path.join(PHANTOM_DIR, f"phantom_{key}")
    if not os.path.isfile(os.path.join(folder, "phantom.txt")):
        args = argparse.Namespace(preset=preset, name=None,
                                  correlated_noise=True if correlated else None)
        make_phantom.make_phantom(make_phantom.resolve_settings(args), PHANTOM_DIR)
    return folder


def options_of(settings):
    """Return the command-line words of the brute-force filter for `settings`."""
    words = ["--patch", *(str(r) for r in settings["patch"]),
             "--search", *(str(r) for r in settings["search"]),
             "--h-factor", f"{settings['h_factor']:g}", "--kernel", settings["kernel"]]
    if settings["normalise"]:
        words.append("--normalise")
    return words


def run_original(log, phantom, name):
    """The rows of the baseline of step 0.4: original program and z filters."""
    program = run_reference.default_program(one_worker=False)
    log.add(phantom=name, set="original", method="unfiltered",
            **baseline.score_columns(score.score(phantom, phantom)))
    for sigma in baseline.SCORE_SIGMAS:
        run = run_reference.run_reference(phantom, sigma, WORK_DIR, program, keep_float=False)
        log.add(phantom=name, set="original", method="original", sigma=float(sigma),
                seconds=run["seconds"], **baseline.score_columns(score.score(phantom, run["vol"])))
    depth = find_volume(phantom)[1][2]
    for kind in ("median", "mean"):
        for radius in baseline.Z_RADII:
            if 2 * radius + 1 > depth:
                continue
            start = time.perf_counter()
            out_vol = z_filter.z_filter(phantom, kind, radius, WORK_DIR)
            log.add(phantom=name, set="original", method=f"z {kind}", search=f"0 0 {radius}",
                    seconds=time.perf_counter() - start,
                    **baseline.score_columns(score.score(phantom, out_vol)))


def run_set(log, phantom, name, set_name):
    """Run the brute-force filter with every setting of one set and score it."""
    program = run_reference.default_program(False, bruteforce=True)
    for index, settings in enumerate(SETS[set_name]):
        run = run_reference.run_reference(phantom, settings["sigma"], WORK_DIR, program,
                                          keep_float=False, tag=f"bf_{set_name}{index}",
                                          options=options_of(settings))
        log.add(phantom=name, set=set_name, method="brute force", sigma=settings["sigma"],
                patch=" ".join(str(r) for r in settings["patch"]),
                search=" ".join(str(r) for r in settings["search"]),
                h_factor=settings["h_factor"],
                distance="mean" if settings["normalise"] else "sum",
                kernel=settings["kernel"], seconds=run["seconds"],
                **baseline.score_columns(score.score(phantom, run["vol"])))
        os.remove(run["vol"])
        os.remove(os.path.splitext(run["vol"])[0] + ".vgi")


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--phantoms", nargs="+", choices=sorted(PHANTOMS), required=True)
    parser.add_argument("--sets", nargs="+", choices=["original"] + sorted(SETS), required=True)
    parser.add_argument("--name", required=True, help="name of the CSV file, without .csv")
    parser.add_argument("--tag", default="windows" if os.name == "nt" else "linux",
                        help="name of the machine in the output folder")
    args = parser.parse_args()

    out_path = os.path.join(REPO_ROOT, "run", "params", args.tag, args.name + ".csv")
    log = baseline.CsvLog(out_path, COLUMNS)
    try:
        for key in args.phantoms:
            phantom = ensure_phantom(key)
            for set_name in args.sets:
                if set_name == "original":
                    run_original(log, phantom, key)
                else:
                    run_set(log, phantom, key, set_name)
    except (ValueError, RuntimeError, OSError) as error:
        print(f"Failed: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
