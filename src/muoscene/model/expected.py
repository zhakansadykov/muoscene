"""model.expected — `expected_map`: scene × flux table × detector × exposure → `ExpectedMap`.

The fast path of the conventions ("Fast and direct paths"), for every detector of a run:

    lines of sight of the detector through the grid  →  `Site.scan`  →  `PathTable`
    →  opacity by material under every hypothesis    →  equivalent opacity (`transport.table`)
    →  flux table, called by the detector with its own energy response  →  counts per bin

and, in every run, the same flux from Mulder itself through the scene on a random sample of
those lines, with the deviation stored in the result (`transport.direct`).

An `ExpectedMap` is an `xarray.Dataset` with the dimensions `detector`, `hypothesis`, `ty`, `tx`
(conventions, "Expected map"). Decisions of 21 November 2026 (`docs/notes/forward.md`):

- `solid_angle_sr` is the whole solid angle of the rectangle of a bin, a property of the grid.
  What a detector accepts of it — within its circle, in front of its plane, above the horizon —
  is `accepted_solid_angle_sr`, and that is what its counts are summed over.
- Bins without a flux. With a sea-level reference above sea level the table has no flux below
  its cut-off: at 1150 m the whole open sky for a response from 0.2 GeV, and elevations below
  15.7 degrees for one from 1 GeV. `undefined` is always given: "refuse" refuses the run,
  "mask" computes the rest and marks the bins. `status` holds the reason for every bin; a bin
  with one line of sight without a flux is undefined as a whole, and its `counts` are nan, the
  mark of a missing number, never a zero. A bin the detector does not accept has the status
  "not accepted" and zero counts: nothing is expected there.
- `reference` is geometry: no accepted line of sight of the bin meets matter under any
  hypothesis. It is defined whether or not the table has a flux for the open sky.
- `transmission` is the counts over the counts the same lines would give with no matter. Where
  the open sky has no flux (`open_sky_status`) it is nan, though the counts under soil exist.
- `zenith_deg` and `azimuth_deg` are those of the centre of a bin on the side of the normal of
  the detector (`Pose.direction`). For a tilted detector a zenith angle above 90 degrees says
  that the centre looks below the horizon; the counts of such a bin come from its lines of
  sight above the horizon alone.
- `opacity_gcm2` is the equivalent opacity of the matter; the table is entered with it less the
  air that matter displaces (`transport.table`).

Not here: significance (roadmap 0.1-09), the configuration hash, the run folder and NetCDF
(0.1-11), scattering and systematic variants (0.3).

Units: counts are tracks, opacity g/cm², solid angles sr, angles degrees, exposure days.

Layer 3.
"""

from __future__ import annotations

import datetime
import json
import math
from collections.abc import Iterable
from dataclasses import dataclass
from importlib import metadata
from typing import Literal

import numpy as np
import xarray as xr

from .. import __version__
from ..detectors import DEFAULT_SUPERSAMPLE, Detector, DirectionGrid, LinesOfSight, Pose
from ..flux import LEGACY_RESPONSE, FluxTable, LegacyTable
from ..geometry import PathTable, Scene
from ..materials import register
from ..site import Site
from ..transport.direct import K_MAX_GEV, LIMIT, NODES_PER_DECADE, DirectCheck, compare_with_direct, seed_of
from ..transport.table import (
    BELOW_CUTOFF,
    DEFINED,
    LOW_ELEVATION,
    NOT_ACCEPTED,
    STATUS_MEANINGS,
    TOO_OPAQUE,
    equivalent_opacity,
    flux_domain,
)
from .hypotheses import Hypothesis, HypothesisSet

CHECK_DIRECTIONS = 256
"""Lines of sight per detector on which a run compares the fast path with Mulder: 0.3 s on a Linux VM"""
TABLE_ALTITUDE_M = 1.0
"""A flux table this close to the altitude of its detector is its own: the open sky moves by 0.02 % a metre"""
SOLID = ("ground", "layer", "structure", "anomaly")
REASONS = (BELOW_CUTOFF, LOW_ELEVATION, TOO_OPAQUE)
VARIABLES = (
    "solid_angle_sr",
    "accepted_solid_angle_sr",
    "opacity_gcm2",
    "counts",
    "transmission",
    "status",
    "open_sky_status",
    "reference",
    "zenith_deg",
    "azimuth_deg",
    "exposure_days",
)
"""Variables of an `ExpectedMap`."""


@dataclass(frozen=True)
class Station:
    """A detector placed in the scene for an exposure, with the flux table of its altitude.

    name            label along the `detector` dimension of a map
    detector        the detector: acceptance, energy response, efficiency, area
    pose            where it is and how it is turned, in ENU
    table           flux table at the altitude of the detector: a `FluxTable` made for the
                    atmosphere of the scene and without a geomagnetic field, or a `LegacyTable`
    exposure_days   length of the exposure, days
    """

    name: str
    detector: Detector
    pose: Pose
    table: FluxTable | LegacyTable
    exposure_days: float

    def __post_init__(self) -> None:
        if not isinstance(self.name, str) or not self.name.strip():
            raise ValueError("a station needs a name")
        if not isinstance(self.detector, Detector):
            raise TypeError(f"{self.name}: detector must be a Detector")
        if not isinstance(self.pose, Pose):
            raise TypeError(f"{self.name}: pose must be a Pose")
        if not isinstance(self.table, FluxTable | LegacyTable):
            raise TypeError(f"{self.name}: table must be a FluxTable or a LegacyTable")
        days = self.exposure_days
        if (
            isinstance(days, bool)
            or not isinstance(days, int | float)
            or not (math.isfinite(days) and days > 0)
        ):
            raise ValueError(f"{self.name}: exposure_days must be a positive number")
        object.__setattr__(self, "exposure_days", float(days))
        legacy = self.detector.response is LEGACY_RESPONSE
        if legacy != isinstance(self.table, LegacyTable):
            raise ValueError(
                f"{self.name}: LEGACY_RESPONSE goes with a LegacyTable, an EnergyResponse with a FluxTable"
            )

    def to_dict(self) -> dict:
        """Plain values for a run record: the detector with its energy response, the pose, the table."""
        detector, pose, table = self.detector, self.pose, self.table
        response = (
            "legacy"
            if detector.response is LEGACY_RESPONSE
            else {
                "kinetic_energy_gev": list(detector.response.kinetic_energy_gev),
                "efficiency": list(detector.response.efficiency),
            }
        )
        described = {"class": type(table).__name__}
        if isinstance(table, FluxTable):
            described |= {
                "altitude_m": table.altitude_m,
                "reference_material": table.meta.get("material_name"),
                "key": getattr(table, "key", None),
            }
        return {
            "detector": {
                "class": type(detector).__name__,
                "tan_max": detector.tan_max,
                "efficiency": detector.efficiency,
                "area_cm2": getattr(detector, "area_cm2", None),
                "energy_response": response,
            },
            "pose": {
                "position_m": list(pose.position_m),
                "azimuth_deg": pose.azimuth_deg,
                "tilt_deg": pose.tilt_deg,
                "roll_deg": pose.roll_deg,
            },
            "exposure_days": self.exposure_days,
            "table": described,
        }


def _materials(scene: Scene) -> dict:
    """Materials of the solid bodies of a scene by their names in Mulder."""
    out = {}
    for body in scene.bodies:
        for part in (body, *getattr(body, "layers", ())):
            composition = getattr(part, "composition", None)
            if composition is not None:
                out[register(composition)] = composition
    return out


def _check_table(station: Station, site: Site, scene: Scene) -> None:
    table = station.table
    if not isinstance(table, FluxTable):
        return
    altitude = site.altitude_m + station.pose.position_m[2]
    if abs(table.altitude_m - altitude) > TABLE_ALTITUDE_M:
        raise ValueError(
            f"{station.name}: the flux table is at {table.altitude_m:g} m and the detector at {altitude:g} m:"
            " a table holds the open sky at the altitude of its detector"
        )
    inputs = table.meta.get("inputs")
    if not inputs:
        return
    if inputs.get("atmosphere") != scene.atmosphere:
        raise ValueError(
            f"{station.name}: the flux table was made for the atmosphere {inputs.get('atmosphere')!r},"
            f" the scene for {scene.atmosphere!r}"
        )
    if inputs.get("date") is not None:
        raise ValueError(
            f"{station.name}: a table with the geomagnetic field of a date holds one azimuth; maps with the"
            " field come with version 0.3"
        )


def _opacity(station: Station, hypothesis: Hypothesis, paths: PathTable, materials, length, air_density):
    """Equivalent opacity of the matter and the opacity the table is entered with, g/cm²."""
    if isinstance(station.table, LegacyTable):
        plain = np.asarray(hypothesis.opacity_g_cm2(paths), dtype=float)
        return plain, plain  # the preliminary calculation: a plain sum of density times length
    equivalent = equivalent_opacity(
        station.table,
        hypothesis.opacity_by_material_g_cm2(paths),
        materials,
        station.detector.response.k_min_gev,
        solid_length_m=length,
        air_density_g_cm3=air_density,
    )
    return equivalent.matter_g_cm2, equivalent.lookup_g_cm2


def _per_bin(lines: LinesOfSight, weights, n_bins: int) -> np.ndarray:
    return np.bincount(lines.bin, weights=weights, minlength=n_bins)


def _counts(station: Station, lines: LinesOfSight, lookup, code, accepted, shape):
    """Counts per bin and their status: nan where a line of sight of the bin has no flux."""
    n_bins = shape[0] * shape[1]
    undefined = np.array([_per_bin(lines, lines.solid_angle_sr * (code == r), n_bins) for r in REASONS])
    partly = undefined.sum(axis=0) > 0
    status = np.where(partly, np.array(REASONS)[undefined.argmax(axis=0)], DEFINED)
    status = np.where(accepted, status, NOT_ACCEPTED).astype(np.int8).reshape(shape)
    good = code == DEFINED
    counts = station.detector.counts(station.table, lines.select(good), lookup[good], station.exposure_days)
    return np.where(status > NOT_ACCEPTED, np.nan, counts), status


def _refusal(station: Station, hypothesis: str, lines: LinesOfSight, code) -> str:
    sine = lines.direction[:, 2]
    parts = []
    for reason in REASONS:
        hit = code == reason
        if np.any(hit):
            low, high = np.degrees(np.arcsin(sine[hit].min())), np.degrees(np.arcsin(sine[hit].max()))
            parts.append(
                f"{int(hit.sum())} {STATUS_MEANINGS[reason].replace('_', ' ')}"
                f" (elevations {low:.1f} to {high:.1f} degrees)"
            )
    return (
        f"{station.name}, hypothesis {hypothesis!r}: the flux table has no flux for the energy response"
        f" {station.detector.response!r} on {int(np.sum(code != DEFINED))} of {len(lines)} accepted lines of"
        f" sight: {'; '.join(parts)}. Nothing is put in its place. Raise the threshold, narrow the"
        ' acceptance, or ask for undefined="mask" to compute the rest and mark these bins; the'
        " reference at the site altitude is roadmap 0.1-05b"
    )


def expected_map(
    site: Site,
    scene: Scene,
    stations: Station | Iterable[Station],
    grid: DirectionGrid,
    hypotheses: Hypothesis | HypothesisSet | Iterable[Hypothesis],
    *,
    rng: np.random.Generator,
    undefined: Literal["refuse", "mask"],
    supersample: int = DEFAULT_SUPERSAMPLE,
    check_directions: int = CHECK_DIRECTIONS,
    check_limit: float | None = LIMIT,
    reference_opacity_g_cm2: float = 0.0,
) -> xr.Dataset:
    """Expected counts of every detector in every bin under every hypothesis: an `ExpectedMap`.

    site, scene        the site and a scene of it (`Site.scene`)
    stations           the detectors, each with its pose, flux table and exposure; different names
    grid               bins in the slopes, the same for every detector
    hypotheses         densities the map is computed for; none of them needs a new scan
    rng                generator the sample of the check is drawn from, made with a seed
                       (`numpy.random.default_rng(seed)`); the seed is written to the map
    undefined          what to do with lines of sight the table has no flux for: "refuse" the
                       run, or "mask" the bins. Always given
    supersample        sub-samples per side of a bin
    check_directions   lines of sight per detector on which the fast path is compared with
                       Mulder through the scene as built; 0 leaves the comparison out and says so
    check_limit        largest deviation accepted in that comparison (RuntimeError above it);
                       None only records
    reference_opacity_g_cm2
                       a bin is a `reference` bin when no accepted line of sight has more
                       equivalent opacity than this under any hypothesis, g/cm²

    The variables and what they mean are in the docstring of the module; the attributes carry
    the versions of the package and of the backends, the stations with their energy responses
    and exposures, the hypotheses, and the deviation of the fast path from the direct one with
    the seed of its sample. Nothing is written to disk (roadmap 0.1-11).
    """
    if not isinstance(site, Site) or not isinstance(scene, Scene):
        raise TypeError("expected_map takes a Site and a Scene of it")
    if not isinstance(grid, DirectionGrid):
        raise TypeError("grid must be a DirectionGrid")
    if undefined not in ("refuse", "mask"):
        raise ValueError('undefined is always given: "refuse" or "mask"')
    stations = (stations,) if isinstance(stations, Station) else tuple(stations)
    if not stations or not all(isinstance(s, Station) for s in stations):
        raise TypeError("stations: at least one Station")
    if len({s.name for s in stations}) != len(stations):
        raise ValueError("the stations of a map have different names")
    hypotheses = hypotheses if isinstance(hypotheses, HypothesisSet) else HypothesisSet(hypotheses)
    if isinstance(check_directions, bool) or not isinstance(check_directions, int) or check_directions < 0:
        raise ValueError("check_directions must be a non-negative integer")
    if not (math.isfinite(reference_opacity_g_cm2) and reference_opacity_g_cm2 >= 0):
        raise ValueError("reference_opacity_g_cm2 must be finite and non-negative")
    seed = seed_of(rng)
    for station in stations:
        _check_table(station, site, scene)

    materials = _materials(scene)
    shape = grid.shape
    n_bins = shape[0] * shape[1]
    solid = [i for i, body in enumerate(scene.table_bodies) if body.kind in SOLID]
    out = {name: [] for name in VARIABLES[1:]}
    checks = {}
    for station in stations:
        detector, pose, table = station.detector, station.pose, station.table
        lines = detector.lines_of_sight(pose, grid, supersample)
        if len(lines) == 0:
            raise ValueError(
                f"{station.name}: no line of sight of the grid is accepted and above the horizon"
            )
        paths = site.scan(scene, pose.position_m, lines.direction)
        hypotheses.check(paths)
        length = paths.length_m[:, solid].sum(axis=1)
        sine = lines.direction[:, 2]
        omega = _per_bin(lines, lines.solid_angle_sr, n_bins)
        accepted = omega > 0
        safe = np.where(accepted, omega, 1.0)

        opacity, counts, status = [], [], []
        heaviest = np.zeros(len(lines))
        for hypothesis in hypotheses:
            matter, lookup = _opacity(station, hypothesis, paths, materials, length, scene.air_density_g_cm3)
            code = flux_domain(table, detector.response, sine, lookup)
            if undefined == "refuse" and np.any(code != DEFINED):
                raise ValueError(_refusal(station, hypothesis.name, lines, code))
            n, s = _counts(station, lines, lookup, code, accepted, shape)
            mean = _per_bin(lines, lines.solid_angle_sr * matter, n_bins) / safe
            opacity.append(np.where(accepted, mean, np.nan).reshape(shape))
            counts.append(n)
            status.append(s)
            heaviest = np.maximum(heaviest, matter)

        sky_code = flux_domain(table, detector.response, sine, np.zeros(len(lines)))
        sky, sky_status = _counts(station, lines, np.zeros(len(lines)), sky_code, accepted, shape)
        counts, status = np.array(counts), np.array(status)
        with np.errstate(invalid="ignore", divide="ignore"):
            ratio = counts / sky
        both = (status == DEFINED) & (sky_status == DEFINED) & (sky > 0)
        centre = grid.directions(pose)
        azimuth = np.degrees(np.arctan2(centre[..., 0], centre[..., 1])) % 360.0
        level = np.hypot(centre[..., 0], centre[..., 1]) < 1e-12
        seen = _per_bin(lines, lines.solid_angle_sr * (heaviest > reference_opacity_g_cm2), n_bins)

        out["accepted_solid_angle_sr"].append(omega.reshape(shape))
        out["opacity_gcm2"].append(np.array(opacity))
        out["counts"].append(counts)
        out["transmission"].append(np.where(both, ratio, np.nan))
        out["status"].append(status)
        out["open_sky_status"].append(sky_status)
        out["reference"].append((accepted & (seen == 0)).reshape(shape))
        out["zenith_deg"].append(np.degrees(np.arccos(np.clip(centre[..., 2], -1.0, 1.0))))
        out["azimuth_deg"].append(np.where(level, 0.0, azimuth))
        out["exposure_days"].append(station.exposure_days)
        checks[station.name] = _check(
            station, site, scene, lines, paths, materials, length, rng, check_directions, check_limit
        )

    cell, per_detector, per_map = (
        ("ty", "tx"),
        ("detector", "ty", "tx"),
        ("detector", "hypothesis", "ty", "tx"),
    )
    flags = {
        "flag_values": np.arange(len(STATUS_MEANINGS), dtype=np.int8),
        "flag_meanings": " ".join(STATUS_MEANINGS),
    }
    dims = {
        "accepted_solid_angle_sr": per_detector,
        "opacity_gcm2": per_map,
        "counts": per_map,
        "transmission": per_map,
        "status": per_map,
        "open_sky_status": per_detector,
        "reference": per_detector,
        "zenith_deg": per_detector,
        "azimuth_deg": per_detector,
        "exposure_days": ("detector",),
    }
    meaning = {
        "solid_angle_sr": "solid angle of the bin, the whole rectangle in the slopes",
        "accepted_solid_angle_sr": "solid angle of the accepted lines of sight of the bin",
        "opacity_gcm2": "mean equivalent opacity of the matter over the accepted lines of sight, g/cm2",
        "counts": "expected number of tracks; nan where status is above 1",
        "transmission": "counts over the counts of the same lines of sight without matter",
        "status": "whether the bin has counts, and if not, why",
        "open_sky_status": "whether the bin has counts without matter, and if not, why",
        "reference": "no accepted line of sight of the bin meets matter under any hypothesis",
        "zenith_deg": "zenith angle of the bin centre in ENU, on the side of the detector normal",
        "azimuth_deg": "azimuth of the bin centre from true north, clockwise",
        "exposure_days": "length of the exposure, days",
    }
    data = {"solid_angle_sr": (cell, np.array(grid.solid_angle_sr), {"long_name": meaning["solid_angle_sr"]})}
    for name, dimensions in dims.items():
        attributes = {"long_name": meaning[name]} | (flags if name.endswith("status") else {})
        data[name] = (dimensions, np.array(out[name]), attributes)
    compared = [c["max"] for c in checks.values() if c["compared"]]
    return xr.Dataset(
        data,
        coords={
            "detector": [s.name for s in stations],
            "hypothesis": list(hypotheses.names),
            "ty": grid.ty_centers,
            "tx": grid.tx_centers,
            "ty_edges": ("ty_edge", np.array(grid.ty_edges)),
            "tx_edges": ("tx_edge", np.array(grid.tx_edges)),
        },
        attrs={
            "title": "muoscene expected map",
            "conventions": "muoscene docs/conventions.md, revision 4",
            "muoscene_version": __version__,
            "mulder_version": metadata.version("mulder"),
            "calzone_version": metadata.version("calzone"),
            "created": datetime.datetime.now(datetime.UTC).isoformat(timespec="seconds"),
            "transport": "continuous, straight lines of sight",
            "undefined": undefined,
            "supersample": supersample,
            "reference_opacity_g_cm2": float(reference_opacity_g_cm2),
            "site": json.dumps(
                {
                    "latitude_deg": site.latitude_deg,
                    "longitude_deg": site.longitude_deg,
                    "altitude_m": site.altitude_m,
                    "scene_azimuth_deg": site.scene_azimuth_deg,
                }
            ),
            "stations": json.dumps({s.name: s.to_dict() for s in stations}),
            "hypotheses": json.dumps({h.name: h.to_dict() for h in hypotheses}),
            "check_seed": seed,
            "check_directions": check_directions,
            "fast_vs_direct": json.dumps(checks),
            "fast_vs_direct_max": max(compared) if compared else float("nan"),
        },
    )


def _check(station, site, scene, lines, paths, materials, length, rng, size, limit) -> dict:
    """The fast path against Mulder through the scene as built, on a random sample of the lines of sight."""
    if isinstance(station.table, LegacyTable):
        return {"compared": False, "why": "the preliminary calculation is not the physics of Mulder"}
    if size == 0:
        return {"compared": False, "why": "check_directions is 0"}
    detector = station.detector
    _, lookup = _opacity(station, Hypothesis.as_built(), paths, materials, length, scene.air_density_g_cm3)
    defined = np.flatnonzero(
        flux_domain(station.table, detector.response, lines.direction[:, 2], lookup) == DEFINED
    )
    if defined.size == 0:
        return {"compared": False, "why": "no line of sight of the scene as built has a flux"}
    chosen = np.zeros(len(lines), dtype=bool)
    chosen[rng.choice(defined, size=min(size, defined.size), replace=False)] = True
    sample = lines.select(chosen)
    fast = detector.flux_m2_s_sr(station.table, sample, lookup[chosen])
    check: DirectCheck = compare_with_direct(
        site, scene, station.pose.position_m, sample.direction, detector.response, fast, limit=limit
    )
    return {
        "compared": True,
        "hypothesis": "the scene as built",
        "n_directions": check.n_directions,
        "mean": check.mean,
        "rms": check.rms,
        "max": check.max,
        "max_at": {"tx": float(sample.tx[check.worst]), "ty": float(sample.ty[check.worst])},
        "limit": limit,
        "k_max_gev": K_MAX_GEV,
        "nodes_per_decade": NODES_PER_DECADE,
    }


def grid_of(expected: xr.Dataset) -> DirectionGrid:
    """The `DirectionGrid` an `ExpectedMap` was computed on, from its bin edges."""
    return DirectionGrid(expected["tx_edges"].values, expected["ty_edges"].values)


def expected_counts(expected: xr.Dataset, detector: str, hypothesis: str) -> np.ma.MaskedArray:
    """Expected counts of one detector under one hypothesis; shape (bins in ty, bins in tx).

    A masked array: the bins the map has no number for (`status` above 1) are masked, so that a
    sum or a Poisson draw over them cannot take a missing number for a zero. Bins the detector
    does not accept hold zero and are not masked.
    """
    chosen = expected.sel(detector=detector, hypothesis=hypothesis)
    status = chosen["status"].values
    return np.ma.MaskedArray(
        np.where(status > NOT_ACCEPTED, 0.0, chosen["counts"].values), mask=status > NOT_ACCEPTED
    )
