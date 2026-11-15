#!/usr/bin/env python3
"""Make a synthetic glass-fibre rod scan with known noise and known materials.

The phantom imitates a CT scan taken inside a pultruded glass-fibre rod
(README chapter 5): glass fibres along z in epoxy, an air-filled bore hole
along the rod axis, and small air voids in the epoxy. Unlike a real scan it
comes with the noise-free volume and the true material of every voxel, so a
filter result can be scored (tests/phantom/score.py).

Output, below the folder given with --out (default build/phantom/):

    <name>/<name>.vol + .vgi                noisy scan, 16-bit unsigned
    <name>/phantom.txt                      every setting, the grey levels,
                                            voxel counts and file checksums
    <name>/voids.csv                        the voids: centre, diameter, voxels
    <name>_clean/<name>_clean.vol + .vgi    the same volume without noise
    <name>_labels/<name>_labels.vol + .vgi  material of every voxel, 8-bit

All .vol files are headerless, little-endian, x fastest, then y, then z.

Label values:  0 air in the bore hole   1 epoxy   2 glass   3 air in a void
               +16 if the voxel is a mixture of two materials (an edge voxel)

Geometry, all lengths in voxels. The rod axis runs along z through the centre
of the x-y plane. Fibres are circles in each slice. Their centres lie on a
hexagonal grid whose spacing gives the requested glass share, each moved by a
random offset small enough that fibres never overlap. From slice to slice all
fibres move together: a rotation about the rod axis (--twist) plus a sideways
shift (--tilt). The bore hole stays on the axis. A voxel at a fibre edge or at
the bore wall gets a mixed value, by a linear ramp one voxel wide.

Voids are balls placed at random where the centre voxel is pure epoxy. Only
epoxy voxels inside the ball become air, so a void larger than the gap between
fibres is cut by the fibres around it; voids.csv gives the voxels it kept.

Noise is added to every voxel independently. It is the sum of 48 random bytes,
scaled to the requested standard deviation, which is very close to a normal
distribution. Only integer random numbers and the operations + - * / sqrt are
used, so the same seed gives the same files on every platform.

Usage (from the repository root):
    python tests/phantom/make_phantom.py --preset fine
    python tests/phantom/make_phantom.py --preset coarse --sigma 17
    python tests/phantom/make_phantom.py --preset fine --size 48 48 24 --name tiny

Needs only the Python standard library (3.9 or newer).
"""

import argparse
import array
import csv
import hashlib
import math
import os
import random
import sys

REPO_ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(REPO_ROOT, "tests", "data"))
from make_crops import write_vgi  # noqa: E402  (shared .vgi writer)

LABEL_AIR, LABEL_EPOXY, LABEL_GLASS, LABEL_VOID = 0, 1, 2, 3
LABEL_EDGE = 16     # added to the label of a voxel that mixes two materials

# Noise: sum of NOISE_BYTES uniform bytes (0...255).
NOISE_BYTES = 48
NOISE_MEAN = NOISE_BYTES * 127.5
NOISE_STD = math.sqrt(NOISE_BYTES * (256 ** 2 - 1) / 12.0)

# Sigma is given on the 0...255 working scale; the files hold 16-bit counts.
WORK_TO_COUNTS = 65535.0 / 255.0

# Presets. "fine" imitates the coming 1 um scans, "coarse" the present 5 um
# scans (the same rod seen with five times larger voxels). Grey levels are in
# 16-bit counts; README section 6.2 explains where they come from.
PRESETS = {
    "fine": {
        "size": (160, 160, 96), "fibre_diameter": 20.0, "bore_diameter": 40.0,
        "twist": 0.0005, "void_diameters": (1, 3, 5, 7), "voxel_mm": 0.001,
    },
    "coarse": {
        "size": (64, 64, 64), "fibre_diameter": 4.0, "bore_diameter": 16.0,
        "twist": 0.0025, "void_diameters": (1, 3), "voxel_mm": 0.005,
    },
}
COMMON_DEFAULTS = {
    "glass_share": 0.60, "tilt": (0.02, 0.01), "voids_per_diameter": 6,
    "air": 18000, "epoxy": 26000, "glass": 38000,
    "sigma": 14.0, "drift": 0.0, "seed": 1,
}


def fibre_grid(spacing, jitter, reach, rng):
    """Return fibre centres (x, y) relative to the rod axis, for slice z = mid.

    Parameters:
        spacing  distance between neighbouring grid points
        jitter   largest random offset of a centre from its grid point
        reach    centres are made for all grid points within this distance
                 of the axis
        rng      random.Random used for the offsets
    """
    row_step = spacing * math.sqrt(3.0) / 2.0
    rows = int(reach / row_step) + 2
    columns = int(reach / spacing) + 2
    centres = []
    for row in range(-rows, rows + 1):
        for column in range(-columns, columns + 1):
            x = (column + (0.5 if row % 2 else 0.0)) * spacing
            y = row * row_step
            if x * x + y * y > reach * reach:
                continue
            # Uniform offset inside a circle, by rejection from a square.
            while True:
                ox = (2.0 * rng.random() - 1.0) * jitter
                oy = (2.0 * rng.random() - 1.0) * jitter
                if ox * ox + oy * oy <= jitter * jitter:
                    break
            centres.append((x + ox, y + oy))
    return centres


def small_angle_cos_sin(angle):
    """Return (cos, sin) of a small angle from their power series.

    Written out with + and * only, so that the result is the same on every
    platform. Accurate to better than 1e-9 for |angle| < 0.5 rad.
    """
    a2 = angle * angle
    cos = 1.0 - a2 / 2.0 * (1.0 - a2 / 12.0 * (1.0 - a2 / 30.0 * (1.0 - a2 / 56.0)))
    sin = angle * (1.0 - a2 / 6.0 * (1.0 - a2 / 20.0 * (1.0 - a2 / 42.0)))
    return cos, sin


class Phantom:
    """Geometry of one phantom; answers "what is at this voxel?" slice by slice."""

    def __init__(self, settings):
        """`settings` is the dict built by resolve_settings()."""
        self.s = settings
        self.width, self.height, self.depth = settings["size"]
        self.axis_x = (self.width - 1) / 2.0
        self.axis_y = (self.height - 1) / 2.0
        self.mid_z = (self.depth - 1) / 2.0
        self.fibre_radius = settings["fibre_diameter"] / 2.0
        self.bore_radius = settings["bore_diameter"] / 2.0

        diameter = settings["fibre_diameter"]
        # Hexagonal grid: glass share = pi d^2 / (2 sqrt(3) spacing^2).
        self.spacing = diameter * math.sqrt(math.pi / (2.0 * math.sqrt(3.0) * settings["glass_share"]))
        # Offsets up to 45 % of the gap keep a tenth of the gap between fibres.
        self.jitter = 0.45 * (self.spacing - diameter)

        half_z = self.depth / 2.0
        tilt_x, tilt_y = settings["tilt"]
        largest_shift = half_z * math.sqrt(tilt_x ** 2 + tilt_y ** 2)
        half_diagonal = math.sqrt(self.width ** 2 + self.height ** 2) / 2.0
        reach = half_diagonal + largest_shift + 2.0 * self.spacing
        if abs(settings["twist"] * half_z) >= 0.5:
            raise ValueError("twist times half the depth must stay below 0.5 rad")

        self.rng = random.Random(settings["seed"])
        self.base_centres = fibre_grid(self.spacing, self.jitter, reach, self.rng)

    def fibre_centres(self, z):
        """Return the fibre centres (x, y) in voxel coordinates for slice z."""
        dz = z - self.mid_z
        cos, sin = small_angle_cos_sin(self.s["twist"] * dz)
        shift_x = self.axis_x + self.s["tilt"][0] * dz
        shift_y = self.axis_y + self.s["tilt"][1] * dz
        return [(bx * cos - by * sin + shift_x, bx * sin + by * cos + shift_y)
                for bx, by in self.base_centres]

    def glass_cover(self, z):
        """Return the glass share (0...1) of every voxel of slice z, row by row."""
        width, height, radius = self.width, self.height, self.fibre_radius
        cover = [0.0] * (width * height)
        reach = radius + 0.5
        for fx, fy in self.fibre_centres(z):
            x_first = max(0, int(math.ceil(fx - reach)))
            x_last = min(width - 1, int(math.floor(fx + reach)))
            y_first = max(0, int(math.ceil(fy - reach)))
            y_last = min(height - 1, int(math.floor(fy + reach)))
            for y in range(y_first, y_last + 1):
                dy2 = (y - fy) ** 2
                row = y * width
                for x in range(x_first, x_last + 1):
                    part = reach - math.sqrt((x - fx) ** 2 + dy2)
                    if part > 0.0:
                        cover[row + x] = min(1.0, cover[row + x] + min(1.0, part))
        return cover

    def air_cover(self):
        """Return the bore-hole air share (0...1) of every voxel of a slice.

        The bore hole does not move with z, so one slice serves all.
        """
        cover = [0.0] * (self.width * self.height)
        if self.bore_radius <= 0.0:
            return cover
        reach = self.bore_radius + 0.5
        for y in range(self.height):
            dy2 = (y - self.axis_y) ** 2
            for x in range(self.width):
                part = reach - math.sqrt((x - self.axis_x) ** 2 + dy2)
                if part > 0.0:
                    cover[y * self.width + x] = min(1.0, part)
        return cover

    def is_pure_epoxy(self, x, y, z):
        """Return True if voxel (x, y, z) holds neither glass nor bore-hole air."""
        if math.sqrt((x - self.axis_x) ** 2 + (y - self.axis_y) ** 2) < self.bore_radius + 0.5:
            return False
        reach = self.fibre_radius + 0.5
        return all((x - fx) ** 2 + (y - fy) ** 2 >= reach * reach
                   for fx, fy in self.fibre_centres(z))

    def place_voids(self):
        """Choose void centres; return a list of dicts (id, x, y, z, diameter).

        Each void is a ball of the given diameter centred on a voxel. The
        centre voxel is pure epoxy, the ball lies inside the volume and balls
        of different voids do not touch. Raises ValueError if no place is
        found for a void.
        """
        voids = []
        for diameter in self.s["void_diameters"]:
            radius = diameter / 2.0
            margin = int(math.ceil(radius))
            if min(self.width, self.height, self.depth) <= 2 * margin:
                raise ValueError(f"the volume is too small for a void of diameter {diameter}")
            for _ in range(self.s["voids_per_diameter"]):
                for _attempt in range(20000):
                    x = self.rng.randrange(margin, self.width - margin)
                    y = self.rng.randrange(margin, self.height - margin)
                    z = self.rng.randrange(margin, self.depth - margin)
                    apart = all((x - v["x"]) ** 2 + (y - v["y"]) ** 2 + (z - v["z"]) ** 2
                                > (radius + v["diameter"] / 2.0 + 2.0) ** 2 for v in voids)
                    if apart and self.is_pure_epoxy(x, y, z):
                        voids.append({"id": len(voids) + 1, "x": x, "y": y, "z": z,
                                      "diameter": diameter, "voxels": 0})
                        break
                else:
                    raise ValueError(f"no free place found for a void of diameter {diameter}")
        return voids


def resolve_settings(args):
    """Combine the preset with the command-line values into one settings dict.

    Raises ValueError for values that make no sense.
    """
    settings = dict(COMMON_DEFAULTS)
    settings.update(PRESETS[args.preset])
    settings["preset"] = args.preset
    for key in ("size", "fibre_diameter", "bore_diameter", "twist", "tilt", "glass_share",
                "void_diameters", "voids_per_diameter", "air", "epoxy", "glass",
                "sigma", "drift", "seed"):
        value = getattr(args, key)
        if value is not None:
            settings[key] = tuple(value) if isinstance(value, list) else value
    settings["name"] = args.name or f"phantom_{args.preset}"

    if min(settings["size"]) < 1:
        raise ValueError("size must be positive")
    if settings["fibre_diameter"] <= 0:
        raise ValueError("fibre diameter must be positive")
    if not 0.0 < settings["glass_share"] < 0.9:
        raise ValueError("glass share must be between 0 and 0.9 (touching circles reach 0.907)")
    if settings["sigma"] < 0:
        raise ValueError("sigma must not be negative")
    for level in ("air", "epoxy", "glass"):
        if not 0 <= settings[level] <= 65535:
            raise ValueError(f"{level} level must be within 0...65535")
    return settings


def sha256_of(path):
    """Return the SHA-256 checksum of a file as hexadecimal text."""
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def make_phantom(settings, out_root):
    """Write the phantom described by `settings` below `out_root`.

    Returns a dict with the output folder of the noisy scan ("folder"), the
    number of voxels per label ("counts", edge voxels counted with their
    material), the number of noisy values that had to be limited to 0...65535
    ("limited") and the list of voids ("voids").
    """
    phantom = Phantom(settings)
    width, height, depth = settings["size"]
    name = settings["name"]
    air, epoxy, glass = settings["air"], settings["epoxy"], settings["glass"]
    sigma_counts = settings["sigma"] * WORK_TO_COUNTS
    noise_scale = sigma_counts / NOISE_STD
    drift = settings["drift"]
    drift_reference = max(width, height) / 2.0

    voids = phantom.place_voids()
    bore = phantom.air_cover()

    folders = {kind: os.path.join(out_root, name + suffix)
               for kind, suffix in (("noisy", ""), ("clean", "_clean"), ("labels", "_labels"))}
    for folder in folders.values():
        os.makedirs(folder, exist_ok=True)
    paths = {kind: os.path.join(folder, os.path.basename(folder) + ".vol")
             for kind, folder in folders.items()}

    counts = {LABEL_AIR: 0, LABEL_EPOXY: 0, LABEL_GLASS: 0, LABEL_VOID: 0}
    edge_voxels = 0
    limited = 0
    pixels = width * height
    noise_rng = random.Random(settings["seed"] + 1000003)

    with open(paths["noisy"], "wb") as noisy_file, open(paths["clean"], "wb") as clean_file, \
            open(paths["labels"], "wb") as label_file:
        for z in range(depth):
            cover = phantom.glass_cover(z)
            values = [0.0] * pixels
            labels = bytearray(pixels)
            for i in range(pixels):
                glass_part = cover[i]
                air_part = bore[i]
                value = epoxy + glass_part * (glass - epoxy)
                if air_part > 0.0:
                    value = value * (1.0 - air_part) + air * air_part
                if air_part >= 0.5:
                    label = LABEL_AIR
                elif glass_part >= 0.5:
                    label = LABEL_GLASS
                else:
                    label = LABEL_EPOXY
                if 0.0 < air_part < 1.0 or (air_part == 0.0 and 0.0 < glass_part < 1.0):
                    label += LABEL_EDGE
                values[i] = value
                labels[i] = label

            # Voids: epoxy voxels inside a ball become air.
            for void in voids:
                radius = void["diameter"] / 2.0
                dz2 = (z - void["z"]) ** 2
                if dz2 >= radius * radius:
                    continue
                reach = int(math.ceil(radius))
                for y in range(void["y"] - reach, void["y"] + reach + 1):
                    for x in range(void["x"] - reach, void["x"] + reach + 1):
                        inside = (x - void["x"]) ** 2 + (y - void["y"]) ** 2 + dz2 < radius * radius
                        i = y * width + x
                        if inside and labels[i] & 15 == LABEL_EPOXY:
                            # An epoxy voxel that was mixed with glass keeps
                            # its edge mark; its glass part is replaced too.
                            labels[i] = LABEL_VOID + (labels[i] & LABEL_EDGE)
                            values[i] = float(air)
                            void["voxels"] += 1

            if drift != 0.0:
                for y in range(height):
                    dy2 = (y - phantom.axis_y) ** 2
                    for x in range(width):
                        r2 = ((x - phantom.axis_x) ** 2 + dy2) / (drift_reference ** 2)
                        i = y * width + x
                        values[i] = air + (values[i] - air) * (1.0 + drift * r2)

            clean = array.array("H", bytes(2 * pixels))
            noisy = array.array("H", bytes(2 * pixels))
            noise_bytes = noise_rng.randbytes(NOISE_BYTES * pixels)
            for i in range(pixels):
                value = values[i]
                clean[i] = min(65535, max(0, int(value + 0.5)))
                noise = (sum(noise_bytes[i * NOISE_BYTES:(i + 1) * NOISE_BYTES]) - NOISE_MEAN) * noise_scale
                count = int(math.floor(value + noise + 0.5))
                if count < 0 or count > 65535:
                    limited += 1
                    count = min(65535, max(0, count))
                noisy[i] = count
                counts[labels[i] & 15] += 1
                edge_voxels += labels[i] >> 4
            if sys.byteorder == "big":
                clean.byteswap()
                noisy.byteswap()
            clean.tofile(clean_file)
            noisy.tofile(noisy_file)
            label_file.write(labels)

    voxel_mm = settings["voxel_mm"]
    for kind, bits in (("noisy", 16), ("clean", 16), ("labels", 8)):
        info = {"bits": bits, "datatype": "unsigned integer",
                "resolution": (voxel_mm,) * 3, "unit": "mm"}
        vol_name = os.path.basename(paths[kind])
        write_vgi(os.path.splitext(paths[kind])[0] + ".vgi", vol_name, settings["size"], info)

    with open(os.path.join(folders["noisy"], "voids.csv"), "w", newline="", encoding="ascii") as out:
        writer = csv.writer(out, lineterminator="\n")
        writer.writerow(["id", "x", "y", "z", "diameter", "voxels"])
        for void in voids:
            writer.writerow([void[key] for key in ("id", "x", "y", "z", "diameter", "voxels")])

    total = pixels * depth
    record = [
        ("name", name), ("preset", settings["preset"]),
        ("size", " ".join(str(n) for n in settings["size"])),
        ("voxel_mm", voxel_mm),
        ("fibre_diameter", settings["fibre_diameter"]),
        ("glass_share_requested", settings["glass_share"]),
        ("fibre_spacing", f"{phantom.spacing:.4f}"),
        ("bore_diameter", settings["bore_diameter"]),
        ("twist_rad_per_voxel", settings["twist"]),
        ("tilt", " ".join(str(t) for t in settings["tilt"])),
        ("void_diameters", " ".join(str(d) for d in settings["void_diameters"])),
        ("voids_per_diameter", settings["voids_per_diameter"]),
        ("level_air", air), ("level_epoxy", epoxy), ("level_glass", glass),
        ("sigma_work", settings["sigma"]), ("sigma_counts", f"{sigma_counts:.2f}"),
        ("drift", drift), ("seed", settings["seed"]),
        ("voxels", total),
        ("voxels_air", counts[LABEL_AIR]), ("voxels_epoxy", counts[LABEL_EPOXY]),
        ("voxels_glass", counts[LABEL_GLASS]), ("voxels_void", counts[LABEL_VOID]),
        ("voxels_edge", edge_voxels),
        ("glass_share_outside_bore",
         f"{counts[LABEL_GLASS] / max(1, total - counts[LABEL_AIR]):.4f}"),
        ("noisy_values_limited", limited),
    ]
    record += [(f"sha256_{kind}", sha256_of(paths[kind])) for kind in ("noisy", "clean", "labels")]
    with open(os.path.join(folders["noisy"], "phantom.txt"), "w", encoding="ascii", newline="\n") as out:
        out.write("".join(f"{key} = {value}\n" for key, value in record))

    return {"folder": folders["noisy"], "counts": counts, "limited": limited, "voids": voids,
            "record": dict(record)}


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--preset", choices=sorted(PRESETS), default="fine",
                        help="fine: 1 um scan (default); coarse: 5 um scan")
    parser.add_argument("--name", help="name of the phantom (default: phantom_<preset>)")
    parser.add_argument("--out", default=os.path.join(REPO_ROOT, "build", "phantom"),
                        help="output folder (default: build/phantom)")
    parser.add_argument("--size", type=int, nargs=3, metavar=("W", "H", "D"),
                        help="volume size in voxels")
    parser.add_argument("--fibre-diameter", type=float, help="fibre diameter in voxels")
    parser.add_argument("--glass-share", type=float, help="glass share of the area, 0...0.9")
    parser.add_argument("--bore-diameter", type=float, help="bore-hole diameter in voxels; 0 = none")
    parser.add_argument("--twist", type=float, help="rotation about the axis, rad per voxel of z")
    parser.add_argument("--tilt", type=float, nargs=2, metavar=("TX", "TY"),
                        help="sideways shift in voxels per voxel of z")
    parser.add_argument("--void-diameters", type=int, nargs="*", metavar="D",
                        help="void diameters in voxels; give none for no voids")
    parser.add_argument("--voids-per-diameter", type=int, help="number of voids of each diameter")
    parser.add_argument("--air", type=int, help="grey level of air, counts")
    parser.add_argument("--epoxy", type=int, help="grey level of epoxy, counts")
    parser.add_argument("--glass", type=int, help="grey level of glass, counts")
    parser.add_argument("--sigma", type=float,
                        help="noise standard deviation on the 0...255 scale (1 = 257 counts)")
    parser.add_argument("--drift", type=float,
                        help="brightness drift: relative gain at the edge of the slice "
                             "against the axis, e.g. 0.1")
    parser.add_argument("--seed", type=int, help="random seed")
    args = parser.parse_args()

    try:
        settings = resolve_settings(args)
        result = make_phantom(settings, args.out)
    except ValueError as error:
        print(f"Failed: {error}", file=sys.stderr)
        return 1

    record = result["record"]
    print(f"{os.path.relpath(result['folder'])}: {record['size'].replace(' ', 'x')}, "
          f"glass share outside the bore {record['glass_share_outside_bore']}, "
          f"{len(result['voids'])} voids, sigma {settings['sigma']:g}, "
          f"{result['limited']} noisy values limited")
    return 0


if __name__ == "__main__":
    sys.exit(main())
