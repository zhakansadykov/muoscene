"""canyon_demo — a synthetic canyon with a pier on each bank and detectors on its walls.

Run from any folder:

    python run.py

For detector thresholds of 1 and 0.2 GeV it prints what every detector is expected to see in
its exposure, how much the void and the decompaction zone change that, and how a detector in a
niche compares with one on pins in front of the wall. The sizes of the scene are invented and
differ from any real site (README.md). Change the parameters below; nothing under the line
"nothing to edit below" needs changing. The first run builds flux tables and compiles the
physics tables of Mulder, up to a minute or two; later runs take seconds. The figures — for
each detector above each threshold the counts, opacity and transmission without the anomalies
and how much each anomaly changes the counts — go to the folder `figures` under the folder the
example is run from (FIGURES below; None writes none).
"""

import json
import math
import pathlib
import sys
import time
from dataclasses import dataclass

import numpy as np

from muoscene import geometry as G
from muoscene import materials as M
from muoscene.detectors import DirectionGrid, EmulsionStack, Pose
from muoscene.flux import EnergyResponse, MulderTable
from muoscene.model import Hypothesis, Station, expected_counts, expected_map
from muoscene.site import Site
from muoscene.transport import BELOW_CUTOFF, DEFINED, LOW_ELEVATION, NOT_ACCEPTED, TOO_OPAQUE
from muoscene.viz import maps as figure_maps
from muoscene.viz import style

# ------------------------------------------------------------------ parameters

SITE = {"latitude_deg": 43.35, "longitude_deg": 79.08, "altitude_m": 1150.0}
CANYON_AZIMUTH_DEG = 30.0  # the canyon runs towards this azimuth from true north: +y of the scene

# Axes of the scene: x across the canyon from the left wall (x = 0) to the right one, y along it,
# z up from the edge of the left brow. "Left" and "right" as seen looking towards CANYON_AZIMUTH_DEG.
CANYON = {
    "floor_z_m": -40.0,  # floor of the canyon
    "width_m": 30.0,  # between the two vertical walls
    "behind_m": 300.0,  # the model goes this far behind each wall,
    "half_length_m": 300.0,  # this far along the canyon either way,
    "z_bottom_m": -150.0,  # and down to here
    "niche_m": (0.5, 0.3, 0.3),  # a niche: depth into the wall, width, height
}
BROWS = {  # the edge is the top of the wall; the surface rises by rise_m over rise_over_m, then is level
    "left": {"edge_z_m": 0.0, "rise_m": 0.0, "rise_over_m": 0.0},
    "right": {"edge_z_m": 5.0, "rise_m": 1.5, "rise_over_m": 8.0},
}
SOIL = {"dry_density_g_cm3": 1.60, "water_percent": 15.0}  # loam, from the surface down to layer 1
LAYERS = [  # loam of the same water content, denser with depth, each down to the next
    # name      top_z_m  dry_density_g_cm3  water_percent
    ("layer 1", -8.0, 1.70, 15.0),
    ("layer 2", -18.0, 1.80, 15.0),
]
PIERS = [  # a pile cap level with the surface and piles under it, of concrete of 2.30 g/cm3
    {
        "name": "left pier",
        "bank": "left",
        "along_m": 0.0,  # axis of the pier along the canyon
        "distance_m": 8.0,  # from the wall to the near face of the cap
        "cap_m": (6.0, 6.0, 1.5),  # across the canyon, along it, thickness
        "piles": (3, 3),
        "pile_diameter_m": 0.8,
        "pile_length_m": 8.0,
        "pile_pitch_m": 2.2,
    },
    {
        "name": "right pier",
        "bank": "right",
        "along_m": 0.0,
        "distance_m": 9.0,
        "cap_m": (6.0, 6.0, 2.0),
        "piles": (3, 3),
        "pile_diameter_m": 1.0,
        "pile_length_m": 10.0,
        "pile_pitch_m": 2.0,
    },
]
ANOMALIES = [  # spheres; drho_rel is the relative change of density in the scene as built: -1 a void
    # name          bank     along_m  behind_m  depth_m  diameter_m  drho_rel
    ("void", "left", 0.0, 6.0, 6.0, 2.0, -1.0),
    ("decompaction", "right", 0.0, 6.0, 5.0, 3.0, -0.2),
]
HYPOTHESES = {  # what the maps are computed for: drho_rel of every anomaly; the first is the baseline
    "none": {"void": 0.0, "decompaction": 0.0},
    "void": {"void": -1.0, "decompaction": 0.0},
    "decompaction": {"void": 0.0, "decompaction": -0.2},
}

# Detectors, one line each. To move a detector, change its line; to put it on pins instead of a
# niche, change "niche" to "wall".
#   bank           "left" or "right"
#   along_m        along the canyon, m, positive towards CANYON_AZIMUTH_DEG; the piers stand at 0
#   depth_m        below the edge of the brow of that bank, m
#   mount          "niche": a niche is cut into the wall there; "wall": on pins, STANDOFF_M in front of it
#   azimuth_deg    of the normal of the film, from true north, clockwise: 300 points into the left
#                  massif (and from the right wall across the canyon), 120 into the right massif
#   tilt_deg       of the normal from the vertical: 0 is a horizontal film; up to 30 is computed in full
#   exposure_days  length of the exposure
DETECTORS = [
    # name   bank     along_m  depth_m  mount    azimuth_deg  tilt_deg  exposure_days
    ("L12n", "left", -1.0, 12.0, "niche", 300.0, 20.0, 60),
    ("L12w", "left", 1.0, 12.0, "wall", 300.0, 20.0, 60),
    ("L22n", "left", -1.0, 22.0, "niche", 300.0, 0.0, 60),
    ("L22w", "left", 1.0, 22.0, "wall", 300.0, 0.0, 60),
    ("R12", "right", 0.0, 12.0, "niche", 120.0, 20.0, 60),
    ("R22", "right", 0.0, 22.0, "wall", 300.0, 30.0, 90),
]
STANDOFF_M = 0.15  # a detector on the wall: its centre this far in front of the wall
FILM = {"area_cm2": 120.0, "n_films": 6, "tan_max": 1.0, "efficiency": 0.9}  # every detector
THRESHOLDS_GEV = (1.0, 0.2)  # kinetic energy a muon needs at the detector to leave a track
SEED = 10  # of the lines of sight on which every run checks the fast path against Mulder
CHECK_DIRECTIONS = 256  # lines of sight of that check, per detector and threshold
VISIBLE = 0.01  # a bin "sees" an anomaly when its counts change by more than this share
FIGURES = "figures"  # folder for the figures, under the folder the example is run from; None for none

# ------------------------------------------------------------------ nothing to edit below

BANKS = ("left", "right")
MOUNTS = ("niche", "wall")
COLUMNS = ("name", "bank", "along_m", "depth_m", "mount", "azimuth_deg", "tilt_deg", "exposure_days")
REASONS = {  # why a bin has no number: what is printed, and the explanation printed once
    BELOW_CUTOFF: (
        "below the cut-off of the reference flux",
        "the sea-level reference flux, carried up to the site, has no muons below its cut-off in the"
        " open sky and under thin soil; a reference at the site altitude is roadmap 0.1-05b",
    ),
    LOW_ELEVATION: ("below 14.5 degrees of elevation", "the flux table starts at 14.5 degrees of elevation"),
    TOO_OPAQUE: ("behind too much matter", "the flux table ends at 3e4 g/cm2 of matter"),
}


class InputError(ValueError):
    """A mistake in the parameters, reported in one line without a traceback."""


@dataclass(frozen=True)
class Row:
    """One line of DETECTORS, checked."""

    name: str
    bank: str
    along_m: float
    depth_m: float
    mount: str
    azimuth_deg: float
    tilt_deg: float
    exposure_days: float


@dataclass(frozen=True)
class Setup:
    """The scene and where every detector is: in the axes of the scene and in ENU, m."""

    site: Site
    canyon: G.Canyon
    scene: G.Scene
    soil: M.Soil
    rows: tuple
    scene_points: dict
    enu: dict


@dataclass(frozen=True)
class Result:
    """What a run gives: the setup, the flux tables by altitude and an ExpectedMap per threshold."""

    setup: Setup
    tables: dict
    maps: dict
    check_directions: int
    seconds: float


def _number(row, column, value):
    if isinstance(value, bool) or not isinstance(value, int | float) or not math.isfinite(value):
        raise InputError(f"detector {row}: {column} must be a number, not {value!r}")
    return float(value)


def check_rows(detectors) -> tuple:
    """The lines of DETECTORS as `Row`s; InputError with the name of the detector for a wrong one."""
    rows, names = [], set()
    for i, line in enumerate(detectors, start=1):
        if not isinstance(line, tuple | list) or len(line) != len(COLUMNS):
            raise InputError(f"line {i} of DETECTORS: {len(COLUMNS)} values are needed: {', '.join(COLUMNS)}")
        name, bank, along, depth, mount, azimuth, tilt, days = line
        if not isinstance(name, str) or not name.strip():
            raise InputError(f"line {i} of DETECTORS: the first value is the name of the detector")
        if name in names:
            raise InputError(f"detector {name}: two detectors have this name")
        names.add(name)
        if bank not in BANKS:
            raise InputError(f"detector {name}: bank is 'left' or 'right', not {bank!r}")
        if mount not in MOUNTS:
            raise InputError(f"detector {name}: mount is 'niche' or 'wall', not {mount!r}")
        tilt = _number(name, "tilt_deg", tilt)
        if not 0 <= tilt <= 90:
            raise InputError(f"detector {name}: tilt_deg is between 0 (horizontal) and 90 (vertical)")
        days = _number(name, "exposure_days", days)
        if not days > 0:
            raise InputError(f"detector {name}: exposure_days must be positive")
        values = (
            _number(name, c, v) for c, v in (("along_m", along), ("depth_m", depth), ("azimuth_deg", azimuth))
        )
        along, depth, azimuth = values
        rows.append(Row(name, bank, along, depth, mount, azimuth, tilt, days))
    if not rows:
        raise InputError("DETECTORS is empty")
    return tuple(rows)


def build(detectors=DETECTORS) -> Setup:
    """The site, the canyon with a niche for every detector that has one, the scene and the points."""
    rows = check_rows(detectors)
    site = Site(**SITE, scene_azimuth_deg=CANYON_AZIMUTH_DEG)
    soil = M.loam(SOIL["dry_density_g_cm3"], SOIL["water_percent"])
    layers = tuple(G.Layer(name, top, M.loam(dry, water)) for name, top, dry, water in LAYERS)
    piers = tuple(G.Pier(**pier, material=M.CONCRETE) for pier in PIERS)
    canyon = G.Canyon(
        left=G.Bank(**BROWS["left"]),
        right=G.Bank(**BROWS["right"]),
        soil=soil,
        layers=layers,
        piers=piers,
        **CANYON,
    )
    extra, points = [], {}
    for row in rows:
        try:
            if row.mount == "niche":
                cut, points[row.name] = canyon.niche(f"niche {row.name}", row.bank, row.along_m, row.depth_m)
                extra.append(cut)
            else:
                points[row.name] = canyon.mount(row.bank, row.along_m, row.depth_m, STANDOFF_M)
        except ValueError as error:
            raise InputError(f"detector {row.name}: {error}") from error
    for name, bank, along, behind, depth, diameter, drho in ANOMALIES:
        centre = canyon.point(bank, along, depth, behind)
        extra.append(G.Anomaly(name, G.Sphere(centre, diameter / 2), drho))
    scene = site.scene(canyon.bodies(extra))
    enu = {name: tuple(float(v) for v in site.enu_from_scene(point)) for name, point in points.items()}
    return Setup(site, canyon, scene, soil, rows, points, enu)


def film(threshold_gev: float) -> EmulsionStack:
    return EmulsionStack(**FILM, response=EnergyResponse.threshold(threshold_gev))


def _guilty(setup: Setup, threshold_gev: float, grid: DirectionGrid):
    """The first detector whose lines of sight the scan refuses, with the reason; None if none."""
    for row in setup.rows:
        pose = Pose(setup.enu[row.name], row.azimuth_deg, row.tilt_deg)
        lines = film(threshold_gev).lines_of_sight(pose, grid)
        try:
            setup.site.scan(setup.scene, pose.position_m, lines.direction)
        except ValueError as error:
            return (
                f"detector {row.name}: {error}. In this example: lower tilt_deg of this detector (into the"
                " massif up to 40 degrees passes), or give the model more room (behind_m, half_length_m)"
            )
    return None


def run(detectors=DETECTORS, check_directions=CHECK_DIRECTIONS) -> Result:
    """Build the scene and compute an ExpectedMap for every threshold of THRESHOLDS_GEV."""
    start = time.perf_counter()
    setup = build(detectors)
    grid = DirectionGrid.square()
    hypotheses = [Hypothesis(name, drho_rel=values) for name, values in HYPOTHESES.items()]
    tables, maps = {}, {}
    for threshold in THRESHOLDS_GEV:
        stations = []
        for row in setup.rows:
            altitude = setup.site.altitude_m + setup.enu[row.name][2]
            if altitude not in tables:  # the flux table of the altitude of the detector
                tables[altitude] = MulderTable(altitude, setup.soil.material)
            pose = Pose(setup.enu[row.name], row.azimuth_deg, row.tilt_deg)
            stations.append(Station(row.name, film(threshold), pose, tables[altitude], row.exposure_days))
        try:
            maps[threshold] = expected_map(
                setup.site,
                setup.scene,
                stations,
                grid,
                hypotheses,
                rng=np.random.default_rng(SEED),
                undefined="mask",
                check_directions=check_directions,
            )
        except ValueError as error:
            raise InputError(_guilty(setup, threshold, grid) or str(error)) from error
    return Result(setup, tables, maps, check_directions, time.perf_counter() - start)


# ------------------------------------------------------------------ what is printed


def echo(setup: Setup) -> list[str]:
    """The detectors as the program has read them: where they are in ENU and where they look."""
    lines = [
        "Detectors as read (ENU: x east, y north, z up from the site, m):",
        f"  {'name':6} {'mount':6} {'x':>8} {'y':>8} {'z':>7} {'altitude':>9}"
        "  normal (azimuth, tilt)  lowest elevation",
    ]
    for row in setup.rows:
        x, y, z = setup.enu[row.name]
        lowest = 90.0 - row.tilt_deg - math.degrees(math.atan(FILM["tan_max"]))
        low = f"{lowest:.0f} deg" if lowest > 0 else "below the horizon (not counted)"
        lines.append(
            f"  {row.name:6} {row.mount:6} {x:8.2f} {y:8.2f} {z:7.2f} {setup.site.altitude_m + z:8.1f} m"
            f"  {row.azimuth_deg:5.0f}, {row.tilt_deg:2.0f} deg          {low}"
        )
    return lines


def _effect(expected, detector: str, hypothesis: str, baseline: str):
    base = expected_counts(expected, detector, baseline)
    other = expected_counts(expected, detector, hypothesis)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.ma.masked_invalid(other / np.ma.masked_equal(base, 0.0) - 1.0)
    return float(other.sum() - base.sum()), float(base.sum()), ratio


def report(result: Result) -> list[str]:
    """The summary of a run, one list of lines per threshold."""
    setup = result.setup
    baseline, *others = HYPOTHESES
    out = [
        f"Scene: synthetic canyon at {SITE['altitude_m']:g} m, {len(setup.rows)} detectors,"
        f" hypotheses {', '.join(HYPOTHESES)}; continuous transport, sea-level reference flux (GCCLY15).",
        f"Every run checks the fast path against Mulder on {result.check_directions} random lines of sight"
        f" per detector (seed {SEED}); the run is refused above 0.5 %.",
    ]
    for threshold, expected in result.maps.items():
        checks = json.loads(expected.attrs["fast_vs_direct"])
        out += ["", f"Above {threshold:g} GeV at the detector"]
        out.append(
            f"  {'detector':8} {'days':>4}  {'bins with a number':18}  {'tracks in them':>14}"
            "  fast path against Mulder"
        )
        undefined, reasons = [], []
        for row in setup.rows:
            status = expected["status"].sel(detector=row.name, hypothesis=baseline).values
            accepted = int(np.count_nonzero(status != NOT_ACCEPTED))
            defined = int(np.count_nonzero(status == DEFINED))
            tracks = float(expected_counts(expected, row.name, baseline).sum())
            check = checks[row.name]
            if check["compared"]:
                text = (
                    f"mean {check['mean'] * 100:+.3f} %, rms {check['rms'] * 100:.3f} %,"
                    f" largest {check['max'] * 100:.3f} % ({check['n_directions']} lines)"
                )
            else:
                text = f"not compared: {check['why']}"
            out.append(
                f"  {row.name:8} {row.exposure_days:4.0f}  {defined:4d} of {accepted:<10d}"
                f"  {tracks:14.4g}  {text}"
            )
            for code, (what, why) in REASONS.items():
                n = int(np.count_nonzero(status == code))
                if n:
                    undefined.append(f"    {row.name}: {n} {what}")
                    reasons += [] if why in reasons else [why]
        if undefined:
            out += ["  Bins without a number:", *undefined, *(f"    ({why})" for why in reasons)]
        else:
            out.append("  Bins without a number: none.")
        for hypothesis in others:
            out.append(f"  {hypothesis} against {baseline}:")
            for row in setup.rows:
                more, total, ratio = _effect(expected, row.name, hypothesis, baseline)
                seen = int(np.ma.sum(np.abs(ratio) > VISIBLE))
                if seen == 0:
                    out.append(f"    {row.name:6} not seen: no bin changes by more than {VISIBLE * 100:g} %")
                    continue
                j, i = np.unravel_index(np.ma.argmax(np.abs(ratio)), ratio.shape)
                zenith = float(expected["zenith_deg"].sel(detector=row.name).values[j, i])
                azimuth = float(expected["azimuth_deg"].sel(detector=row.name).values[j, i])
                out.append(
                    f"    {row.name:6} {more:+8.0f} tracks ({more / total * 100:+.2f} %), {seen} bins beyond"
                    f" {VISIBLE * 100:g} %; largest {ratio[j, i] * 100:+.1f} % at azimuth {azimuth:.0f},"
                    f" elevation {90.0 - zenith:.0f} deg"
                )
        pairs = _pairs(setup.rows)
        if pairs:
            out.append("  Niche against wall (same bank, depth, azimuth and tilt):")
        for niche, wall in pairs:
            numbers = []
            for row in (niche, wall):
                status = expected["status"].sel(detector=row.name, hypothesis=baseline).values
                sky = int(expected["reference"].sel(detector=row.name).values.sum())
                effects = [_effect(expected, row.name, h, baseline)[0] for h in others]
                numbers.append(
                    (
                        float(expected_counts(expected, row.name, baseline).sum()),
                        effects,
                        int(np.count_nonzero(status == DEFINED)),
                        sky,
                    )
                )
            (n0, e0, d0, s0), (n1, e1, d1, s1) = numbers
            effects = "; ".join(
                f"{h} {a:+.0f} / {b:+.0f}"
                for h, a, b in zip(others, e0, e1, strict=True)
                if abs(a) + abs(b) >= 0.5
            )
            out.append(
                f"    {niche.name} / {wall.name}: tracks {n0:.4g} / {n1:.4g}; {effects or 'no anomaly seen'};"
                f" bins with a number {d0} / {d1}; open-sky bins {s0} / {s1}"
            )
    out += ["", f"Done in {result.seconds:.1f} s."]
    return out


def _pairs(rows) -> list:
    """Detectors in a niche and on the wall at the same bank, depth, azimuth and tilt."""
    out = []
    for niche in rows:
        if niche.mount != "niche":
            continue
        for wall in rows:
            same = (wall.bank, wall.depth_m, wall.azimuth_deg, wall.tilt_deg) == (
                niche.bank,
                niche.depth_m,
                niche.azimuth_deg,
                niche.tilt_deg,
            )
            if wall.mount == "wall" and same:
                out.append((niche, wall))
                break
    return out


def figures(result: Result, folder=FIGURES, names=None) -> list[pathlib.Path]:
    """Write a sheet of maps for every detector above every threshold into `folder`: counts,
    opacity and transmission under the first hypothesis, the ratio of the counts under each other
    hypothesis to it. `names` limits the detectors. Returns the paths."""
    folder = pathlib.Path(folder)
    names = [row.name for row in result.setup.rows] if names is None else list(names)
    out = []
    for threshold, expected in result.maps.items():
        for name in names:
            figure = figure_maps.sheet_figure(expected, name)
            out.append(style.save(figure, folder / f"canyon_demo_{name}_{threshold:g}GeV.png"))
    return out


def main(detectors=None) -> int:
    """Print the input, run, print the summary. 0 on success, 2 for a mistake in the parameters,
    1 when the check against Mulder refuses the run."""
    detectors = DETECTORS if detectors is None else detectors
    try:
        print("\n".join(echo(build(detectors))))
        print(
            "\nComputing; the first run builds flux tables and compiles Mulder's tables,"
            " up to a minute or two."
        )
        result = run(detectors)
        print("\n".join(report(result)))
        if FIGURES is not None:
            written = figures(result, FIGURES)
            print(f"Figures: {len(written)} files in {pathlib.Path(FIGURES).resolve()}")
    except ValueError as error:
        print(f"canyon_demo: {error}", file=sys.stderr)
        return 2
    except RuntimeError as error:
        print(f"canyon_demo: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
