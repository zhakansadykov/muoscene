"""slab — films under a uniform slab of soil: a scene whose answer is known.

Run from any folder:

    python run.py

Three horizontal films lie under 0.5, 2 and 20 m of flat ground. Two answers are known. The
opacity along every line of sight is the density of the soil times the depth over the cosine of
the zenith angle, and the map holds exactly that. The flux under the slab is what Mulder gives:
every run compares the fast path with Mulder on random lines of sight, and this example also
computes the flux along the vertical and along 45 degrees, the directions of the numbers frozen
for Mulder in the tests of the program. It prints "checks passed: N of N" and exits with 0 when
everything holds, with 1 otherwise. The first run builds three flux tables and compiles the
physics tables of Mulder, up to a minute; later runs take seconds. The figures — counts,
opacity and transmission of each film above each threshold — go to the folder `figures` under
the folder the example is run from (FIGURES below; None writes none).
"""

import json
import math
import pathlib
import sys
import time

import numpy as np

from muoscene import geometry as G
from muoscene import materials as M
from muoscene.detectors import DirectionGrid, EmulsionStack, Pose
from muoscene.flux import EnergyResponse, MulderTable
from muoscene.model import Hypothesis, Station, expected_map
from muoscene.site import Site
from muoscene.transport import DEFINED, NOT_ACCEPTED
from muoscene.viz import maps as figure_maps
from muoscene.viz import style

# ------------------------------------------------------------------ parameters

SITE = {"latitude_deg": 43.35, "longitude_deg": 79.08, "altitude_m": 1150.0}
SOIL = {"dry_density_g_cm3": 1.60, "water_percent": 15.0}  # loam
GROUND = {"half_size_m": 200.0, "depth_m": 60.0}  # flat, its surface at z = 0
DETECTORS = [
    # name    x_m    y_m  depth_m  exposure_days   (ENU; horizontal films inside the ground)
    ("0.5 m", -10.0, 0.0, 0.5, 60),
    ("2 m", 0.0, 0.0, 2.0, 60),
    ("20 m", 10.0, 0.0, 20.0, 60),
]
FILM = {"area_cm2": 120.0, "n_films": 6, "tan_max": 1.0, "efficiency": 0.9}
THRESHOLDS_GEV = (1.0, 0.2)
SEED = 10
CHECK_DIRECTIONS = 256
FIGURES = "figures"  # folder for the figures, under the folder the example is run from; None for none

# ------------------------------------------------------------------ nothing to edit below

M2_PER_CM2 = 1.0e-4
DAY_S = 86400.0
PENCIL = DirectionGrid([-1.0e-3, 1.0e-3], [-1.0e-3, 1.0e-3])
"""One bin 0.11 degrees wide around the normal of a film: its counts give the flux along the normal."""
ZENITH_DEG = (0.0, 45.0)
"""Directions of the flux along a line: those of the frozen `threshold_slab` cases of the tests."""


def film(threshold_gev: float) -> EmulsionStack:
    return EmulsionStack(**FILM, response=EnergyResponse.threshold(threshold_gev))


def run(check_directions=CHECK_DIRECTIONS, soil=None, density_g_cm3=None) -> dict:
    """The maps on the default grid and along two lines, for every threshold.

    soil, density_g_cm3   another material of the slab and its density (the tests give the loam of
                          the frozen cases); by default the loam of SOIL
    """
    site = Site(**SITE)
    soil = M.loam(SOIL["dry_density_g_cm3"], SOIL["water_percent"]) if soil is None else soil
    half = GROUND["half_size_m"]
    ground = G.Ground(
        "ground", [(-half, 0.0), (half, 0.0)], -GROUND["depth_m"], -half, half, soil, density_g_cm3
    )
    scene = site.scene([ground])
    material = soil.material if isinstance(soil, M.Soil) else soil
    density = ground.bulk_density_g_cm3
    tables = {depth: MulderTable(site.altitude_m - depth, material) for _, _, _, depth, _ in DETECTORS}
    out = {"site": site, "scene": scene, "density_g_cm3": density, "tables": tables, "maps": {}, "lines": {}}
    for threshold in THRESHOLDS_GEV:
        stations = [
            Station(name, film(threshold), Pose((x, y, -depth)), tables[depth], days)
            for name, x, y, depth, days in DETECTORS
        ]
        out["maps"][threshold] = expected_map(
            site,
            scene,
            stations,
            DirectionGrid.square(),
            Hypothesis.as_built("slab"),
            rng=np.random.default_rng(SEED),
            undefined="mask",
            check_directions=check_directions,
        )
        # along a line: a film whose normal is that line, and one bin around it; every sub-sample
        # of the bin is compared with Mulder
        pencils = [
            Station(
                f"{name} at {zenith:g}",
                film(threshold),
                Pose((x, y, -depth), 0.0, zenith),
                tables[depth],
                days,
            )
            for name, x, y, depth, days in DETECTORS
            for zenith in ZENITH_DEG
        ]
        out["lines"][threshold] = expected_map(
            site,
            scene,
            pencils,
            PENCIL,
            Hypothesis.as_built("slab"),
            rng=np.random.default_rng(SEED),
            undefined="refuse",
            check_directions=16,
        )
    return out


def flux_along(expected, name: str, threshold_gev: float, days: float, tilt_deg: float) -> float:
    """Flux along the normal of a pencil station, m-2 s-1 sr-1: its counts over area, angle and time."""
    detector = film(threshold_gev)
    lines = detector.lines_of_sight(Pose((0.0, 0.0, 0.0), 0.0, tilt_deg), PENCIL)
    exposure = float(
        np.sum(detector.projected_area_cm2(lines.tx, lines.ty) * M2_PER_CM2 * lines.solid_angle_sr)
        * detector.efficiency
        * days
        * DAY_S
    )
    return float(expected["counts"].sel(detector=name).values.sum()) / exposure


def along(result) -> dict:
    """Flux along the vertical and along 45 degrees under every depth, with the deviation from Mulder."""
    out = {}
    for threshold, expected in result["lines"].items():
        checks = json.loads(expected.attrs["fast_vs_direct"])
        for name, _, _, depth, days in DETECTORS:
            for zenith in ZENITH_DEG:
                station = f"{name} at {zenith:g}"
                out[(threshold, depth, zenith)] = (
                    flux_along(expected, station, threshold, days, zenith),
                    checks[station]["mean"],
                )
    return out


def checks(result) -> list[tuple[bool, str]]:
    """The known answers of the slab, each with whether it holds and what was found."""
    out, worst, defined, total = [], 0.0, 0, 0
    density = result["density_g_cm3"]
    for threshold, expected in result["maps"].items():
        for name, _, _, depth, _ in DETECTORS:
            lines = film(threshold).lines_of_sight(Pose((0.0, 0.0, -depth)), DirectionGrid.square())
            n = DirectionGrid.square().shape[0] * DirectionGrid.square().shape[1]
            omega = np.bincount(lines.bin, lines.solid_angle_sr, n)
            analytic = np.bincount(
                lines.bin, lines.solid_angle_sr * density * depth * 100.0 / lines.direction[:, 2], n
            )
            status = expected["status"].sel(detector=name).values[0].ravel()
            got = expected["opacity_gcm2"].sel(detector=name).values[0].ravel()
            inside = status != NOT_ACCEPTED
            worst = max(worst, float(np.max(np.abs(got[inside] / (analytic[inside] / omega[inside]) - 1.0))))
            defined += int(np.count_nonzero(status == DEFINED))
            total += int(np.count_nonzero(inside))
    out.append((worst < 1e-9, f"opacity is density x depth / cos(zenith) in every bin: within {worst:.1e}"))
    out.append(
        (defined == total, f"under 0.5 m of soil and more every bin has a number: {defined} of {total}")
    )
    deviations = [
        check["max"]
        for expected in result["maps"].values()
        for check in json.loads(expected.attrs["fast_vs_direct"]).values()
        if check["compared"]
    ]
    largest = max(deviations) if deviations else math.nan
    out.append(
        (
            bool(deviations) and largest < 5e-3,
            f"fast path against Mulder in the maps: largest {largest * 100:.3f} %",
        )
    )
    lines = [abs(deviation) for _, deviation in along(result).values()]
    out.append(
        (
            max(lines) < 5e-3,
            f"along the vertical and 45 degrees against Mulder: largest {max(lines) * 100:.3f} %",
        )
    )
    return out


def report(result) -> list[str]:
    out = [
        f"Uniform slab of {result['density_g_cm3']:.3f} g/cm3 at {SITE['altitude_m']:g} m;"
        " sea-level reference flux (GCCLY15) carried up by Mulder; flux in muons m-2 s-1 sr-1."
    ]
    lines = along(result)
    for threshold, expected in result["maps"].items():
        out.append(f"Above {threshold:g} GeV:")
        for name, _, _, depth, _ in DETECTORS:
            counts = float(np.nansum(expected["counts"].sel(detector=name).values))
            vertical, v_dev = lines[(threshold, depth, 0.0)]
            slant, s_dev = lines[(threshold, depth, 45.0)]
            out.append(
                f"  under {name:5}: {counts:.4g} tracks; vertical {vertical:.4g}"
                f" ({v_dev * 100:+.3f} % from Mulder), 45 degrees {slant:.4g} ({s_dev * 100:+.3f} %)"
            )
    return out


def figures(result, folder=FIGURES, names=None) -> list[pathlib.Path]:
    """Write a sheet of maps for each film above each threshold into `folder`: counts, opacity and
    transmission. `names` limits the films. Returns the paths."""
    folder = pathlib.Path(folder)
    names = [row[0] for row in DETECTORS] if names is None else list(names)
    out = []
    for threshold, expected in result["maps"].items():
        for name in names:
            figure = figure_maps.sheet_figure(expected, name)
            out.append(style.save(figure, folder / f"slab_{name.replace(' ', '')}_{threshold:g}GeV.png"))
    return out


def main() -> int:
    """Run, print the numbers and the checks; 0 when every check holds, 1 otherwise."""
    start = time.perf_counter()
    try:
        result = run()
    except (ValueError, RuntimeError) as error:
        print(f"slab: {error}", file=sys.stderr)
        return 1
    print("\n".join(report(result)))
    found = checks(result)
    for holds, text in found:
        print(f"  [{'ok' if holds else 'FAILED'}] {text}")
    if FIGURES is not None:
        written = figures(result, FIGURES)
        print(f"figures: {len(written)} files in {pathlib.Path(FIGURES).resolve()}")
    passed = sum(holds for holds, _ in found)
    print(f"checks passed: {passed} of {len(found)} ({time.perf_counter() - start:.1f} s)")
    return 0 if passed == len(found) else 1


if __name__ == "__main__":
    raise SystemExit(main())
