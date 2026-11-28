"""open_sky — two films under the open sky: a scene whose answer is known.

Run from any folder:

    python run.py

Nothing stands between the films and the sky, so the expected map is the flux table itself:
in every bin the counts are the flux of the table times the projected area, the solid angle,
the efficiency and the exposure. The run checks that, checks the fast path against Mulder, and
shows what the sea-level reference flux cannot give at 1150 m: no number above 0.2 GeV in the
open sky, and none below 15.7 degrees of elevation above 1 GeV. It is also the first check of
an installation: it prints "checks passed: N of N" and exits with 0 when everything holds,
with 1 otherwise. The first run builds a flux table and compiles the physics tables of Mulder,
up to a minute; later runs take seconds. The figures — the map of the counts of each film
above each threshold and the flux of the table against the zenith angle — go to the folder
`figures` under the folder the example is run from (FIGURES below; None writes none).
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
from muoscene.transport import BELOW_CUTOFF, DEFINED, LOW_ELEVATION, NOT_ACCEPTED, flux_domain
from muoscene.viz import curves, style
from muoscene.viz import maps as figure_maps

# ------------------------------------------------------------------ parameters

SITE = {"latitude_deg": 43.35, "longitude_deg": 79.08, "altitude_m": 1150.0}
GROUND = {
    "dry_density_g_cm3": 1.60,
    "water_percent": 15.0,
    "half_size_m": 100.0,
    "depth_m": 20.0,
}  # flat, z = 0
DETECTORS = [
    # name     x_m  y_m  height_m  azimuth_deg  tilt_deg  exposure_days   (ENU; height above the ground)
    ("flat", 0.0, 0.0, 0.5, 0.0, 0.0, 60),
    ("tilted", 5.0, 0.0, 0.5, 180.0, 60.0, 60),  # looks south, down to 15 degrees below the horizon
]
FILM = {"area_cm2": 120.0, "n_films": 6, "tan_max": 1.0, "efficiency": 0.9}
THRESHOLDS_GEV = (1.0, 0.2)
SEED = 10
CHECK_DIRECTIONS = 256
FIGURES = "figures"  # folder for the figures, under the folder the example is run from; None for none

# ------------------------------------------------------------------ nothing to edit below

M2_PER_CM2 = 1.0e-4
DAY_S = 86400.0


def film(threshold_gev: float) -> EmulsionStack:
    return EmulsionStack(**FILM, response=EnergyResponse.threshold(threshold_gev))


def run(check_directions=CHECK_DIRECTIONS):
    """The scene, one flux table and an ExpectedMap per threshold: (site, scene, poses, table, maps)."""
    site = Site(**SITE)
    soil = M.loam(GROUND["dry_density_g_cm3"], GROUND["water_percent"])
    half = GROUND["half_size_m"]
    ground = G.Ground("ground", [(-half, 0.0), (half, 0.0)], -GROUND["depth_m"], -half, half, soil)
    scene = site.scene([ground])
    poses = {name: Pose((x, y, h), azimuth, tilt) for name, x, y, h, azimuth, tilt, _ in DETECTORS}
    days = {row[0]: row[6] for row in DETECTORS}
    tables = {}
    for pose in poses.values():
        altitude = site.altitude_m + pose.position_m[2]
        if altitude not in tables:
            tables[altitude] = MulderTable(altitude, soil.material)
    maps = {}
    for threshold in THRESHOLDS_GEV:
        stations = [
            Station(name, film(threshold), pose, tables[site.altitude_m + pose.position_m[2]], days[name])
            for name, pose in poses.items()
        ]
        maps[threshold] = expected_map(
            site,
            scene,
            stations,
            DirectionGrid.square(),
            Hypothesis.as_built("open sky"),
            rng=np.random.default_rng(SEED),
            undefined="mask",
            check_directions=check_directions,
        )
    return site, scene, poses, tables, maps


def by_hand(table, threshold_gev: float, pose: Pose, days: float) -> np.ndarray:
    """Counts per bin from the table alone: flux × projected area × solid angle × efficiency × time."""
    detector = film(threshold_gev)
    lines = detector.lines_of_sight(pose, DirectionGrid.square())
    shape = DirectionGrid.square().shape
    out = np.zeros(shape[0] * shape[1])
    sine = lines.direction[:, 2]
    good = flux_domain(table, detector.response, sine, np.zeros(len(lines))) == DEFINED
    flux = np.zeros(len(lines))
    flux[good] = table.integral(detector.response, sine[good], 0.0)
    area = detector.projected_area_cm2(lines.tx, lines.ty) * M2_PER_CM2
    np.add.at(out, lines.bin, flux * area * lines.solid_angle_sr * detector.efficiency * days * DAY_S)
    return out.reshape(shape)


def checks(site, scene, poses, tables, maps) -> list[tuple[bool, str]]:
    """The known answers of the scene, each with whether it holds and what was found."""
    out = []
    days = {row[0]: row[6] for row in DETECTORS}
    worst = 0.0
    for threshold, expected in maps.items():
        for name, pose in poses.items():
            status = expected["status"].sel(detector=name).values[0]
            counts = expected["counts"].sel(detector=name).values[0]
            if np.any(status == DEFINED):
                hand = by_hand(tables[site.altitude_m + pose.position_m[2]], threshold, pose, days[name])
                inside = status == DEFINED
                worst = max(worst, float(np.max(np.abs(counts[inside] / hand[inside] - 1.0))))
    out.append(
        (worst < 1e-12, f"counts are those of the flux table in every bin with a number: within {worst:.1e}")
    )
    sky = all(
        bool(np.all(expected["reference"].values[(expected["status"].values[:, 0] != NOT_ACCEPTED)]))
        for expected in maps.values()
    )
    out.append((sky, "every accepted bin sees the open sky"))
    deviations = []
    for expected in maps.values():
        for check in json.loads(expected.attrs["fast_vs_direct"]).values():
            if check["compared"]:
                deviations.append(check["max"])
    largest = max(deviations) if deviations else math.nan
    out.append(
        (
            bool(deviations) and largest < 5e-3,
            f"fast path against Mulder: largest deviation {largest * 100:.3f} %",
        )
    )
    if 0.2 in maps:
        n = int(np.count_nonzero(maps[0.2]["status"].values == DEFINED))
        out.append(
            (n == 0, f"above 0.2 GeV no bin has a number, as the cut-off of the reference says: {n} have")
        )
    if 1.0 in maps and "tilted" in poses:
        expected = maps[1.0]
        status = expected["status"].sel(detector="tilted").values[0]
        zenith = expected["zenith_deg"].sel(detector="tilted").values
        below = zenith > 92.0
        low = int(np.count_nonzero(status == LOW_ELEVATION))
        cut = int(np.count_nonzero(status == BELOW_CUTOFF))
        holds = bool(np.all(status[below] == NOT_ACCEPTED)) and low > 0 and cut > 0
        out.append(
            (
                holds,
                f"tilted film above 1 GeV: nothing counted below the horizon; {low} bins below 14.5"
                f" degrees and {cut} below the cut-off have no number",
            )
        )
    return out


def report(site, scene, poses, tables, maps) -> list[str]:
    out = [f"Open sky at {site.altitude_m:g} m; sea-level reference flux (GCCLY15) carried up by Mulder."]
    for altitude, table in tables.items():
        vertical = float(table.integral(EnergyResponse.threshold(1.0), 1.0, 0.0))
        out.append(
            f"  flux table at {altitude:g} m: cut-off {float(table.cutoff(1.0, 0.0)):.3f} GeV at the"
            f" vertical; {vertical:.2f} muons m-2 s-1 sr-1 above 1 GeV at the vertical"
        )
    for threshold, expected in maps.items():
        out.append(f"Above {threshold:g} GeV:")
        for name in poses:
            status = expected["status"].sel(detector=name).values[0]
            counts = float(np.nansum(expected["counts"].sel(detector=name).values[0]))
            defined = int(np.count_nonzero(status == DEFINED))
            tracks = f"tracks {counts:.4g}" if defined else "no number of tracks"
            out.append(
                f"  {name:7} bins: {defined} with a number,"
                f" {int(np.count_nonzero(status == BELOW_CUTOFF))} below the cut-off,"
                f" {int(np.count_nonzero(status == LOW_ELEVATION))} below 14.5 degrees,"
                f" {int(np.count_nonzero(status == NOT_ACCEPTED))} not accepted; {tracks}"
            )
    return out


def figures(result, folder=FIGURES) -> list[pathlib.Path]:
    """Write the figures of a run into `folder`: the counts of each film above each threshold, and
    the flux of the table against the zenith angle above the thresholds. Returns the paths."""
    site, scene, poses, tables, maps = result
    folder = pathlib.Path(folder)
    out = []
    for threshold, expected in maps.items():
        for name in poses:
            figure = figure_maps.map_figure(expected, name, "counts")
            out.append(style.save(figure, folder / f"open_sky_{name}_{threshold:g}GeV.png"))
    responses = [EnergyResponse.threshold(t) for t in THRESHOLDS_GEV]
    for n, table in enumerate(tables.values()):
        name = "open_sky_flux.png" if n == 0 else f"open_sky_flux_{n}.png"
        out.append(style.save(curves.flux_figure(table, responses), folder / name))
    return out


def main() -> int:
    """Run, print the numbers and the checks; 0 when every check holds, 1 otherwise."""
    start = time.perf_counter()
    try:
        result = run()
    except (ValueError, RuntimeError) as error:
        print(f"open_sky: {error}", file=sys.stderr)
        return 1
    print("\n".join(report(*result)))
    found = checks(*result)
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
