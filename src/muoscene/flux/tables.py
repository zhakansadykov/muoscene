"""flux.tables — `MulderTable`: table built with Mulder at the detector altitude; disk cache.

Mulder computes the open-sky differential flux at the altitude of the detector on a grid of
kinetic energy and sine of elevation, for an atmosphere model and, optionally, the geomagnetic
field of a date; the CSDA range and stopping power of the reference material carry it under an
opacity (`base`). Every table is then compared with Mulder itself under a dense slab of the
reference material at random energies, elevations and opacities, and refused if they disagree.

Reference flux: version 0.1.0 has GCCLY15 at sea level only. Above sea level it has a cut-off,
which is found for every elevation, stored in the table and enforced by it. A reference at the
site altitude (roadmap 0.1-05b) will be passed here as a table.

Units at this boundary: g/cm³ and g/cm² to kg/m³ and m of Mulder; energies are kinetic, in GeV.
"""

from __future__ import annotations

import hashlib
import json
import os
import pathlib
import zipfile

import mulder
import numpy as np

from ..materials.backend import compiled, definition, mulder_density_kg_m3, register
from ..materials.composition import Material, Mixture
from .atmosphere import REFERENCE, check_atmosphere, fluxmeter, iso_date, open_sky_cutoff, states
from .base import OPACITY_MAX_G_CM2, FluxTable

K_MIN_GEV = 1.0e-3
"""Lowest kinetic energy of a table: the first energy of the physics tables of Mulder."""
K_MAX_GEV = 1.0e4
"""Highest kinetic energy at the detector: 3e-6 of the muons under 3e4 g/cm² are above it."""
CUTOFF_MARGIN_GEV = 0.01
"""With a cut-off the table starts this much above it (`docs/notes/flux.md`)."""
NODES_PER_DECADE = 100
"""Energy nodes of the flux per decade."""
RANGE_NODES_PER_DECADE = 200
"""Energy nodes of the range and stopping power per decade."""
SIN_ELEVATION = tuple(np.linspace(0.25, 1.0, 151).tolist())
"""Default nodes of the sine of elevation: from 14.5 degrees up, in steps of 0.005."""
DENSE_SLAB_KG_M3 = 1.0e5
"""Density of the slab of the check: 3e4 g/cm² are 3 m thick, so the slab stays at the detector altitude."""
CHECK_POINTS = 64
"""Random points at which a new table is compared with Mulder under the slab."""
CHECK_SEED = 5
CHECK_K_MAX_GEV = 1.0e3
"""Highest energy of the comparison: above 1 TeV the flux of Mulder under a thick slab ripples by 0.3 %."""
CHECK_LIMIT = 2.0e-3
"""Largest relative difference at which a table is still accepted; 20 000 points stay within 4.5e-4."""
BUILD = 1
"""Version of the algorithm: part of the cache key."""


def default_cache_dir() -> pathlib.Path:
    """Where tables are cached: $MUOSCENE_CACHE, else $XDG_CACHE_HOME/muoscene, else ~/.cache/muoscene."""
    explicit = os.environ.get("MUOSCENE_CACHE")
    if explicit:
        return pathlib.Path(explicit).expanduser()
    base = os.environ.get("XDG_CACHE_HOME")
    return (pathlib.Path(base).expanduser() if base else pathlib.Path.home() / ".cache") / "muoscene"


def _reference(reference) -> str:
    """The reference flux chosen by the conventions; only the sea-level one exists in 0.1.0."""
    if isinstance(reference, str) and reference in ("auto", "sea-level"):
        return REFERENCE
    if isinstance(reference, str) and reference != "site-altitude":
        raise ValueError(f'reference must be "auto", "sea-level" or "site-altitude", not {reference!r}')
    raise NotImplementedError(
        "a reference flux at the site altitude, by name or as a table, comes with roadmap 0.1-05b"
        " (version 0.1.1); version 0.1.0 has the sea-level reference only"
    )


class MulderTable(FluxTable):
    """Flux table of a detector at `altitude_m`, built by Mulder and cached on disk.

    altitude_m        altitude of the detector above sea level, m (negative below it)
    material          reference material of the opacity: a `materials.Material` or `Mixture`
    latitude_deg, longitude_deg
                      of the site, degrees; needed with `date`, otherwise those of Mulder (0, 0)
    atmosphere        density profile of the air: a name from `atmosphere_models()`
    date              "YYYY-MM-DD" for the geomagnetic field (IGRF) of that date; None for no field
    azimuth_deg       azimuth of the lines of sight, from north, clockwise; needed with `date`,
                      because the field makes the flux depend on it, and refused without
    reference         "auto" follows the conventions: GCCLY15 at sea level in version 0.1.0
    cache, cache_dir  whether to read and write the table on disk, and in which directory;
                      by default `default_cache_dir() / "flux"`, outside any repository
    sin_elevation     nodes of the sine of elevation
    nodes_per_decade  energy nodes per decade
    check_points      random points of the comparison with Mulder under a dense slab

    The cache key is a hash of everything above, of the definition of the material, of the
    version of Mulder and of the algorithm. `from_cache` tells whether the table was read;
    `meta` holds the inputs and the result of the comparison.
    """

    def __init__(
        self,
        altitude_m: float,
        material: Material | Mixture,
        *,
        latitude_deg: float | None = None,
        longitude_deg: float | None = None,
        atmosphere: str = "us-standard",
        date=None,
        azimuth_deg: float | None = None,
        reference="auto",
        cache: bool = True,
        cache_dir=None,
        sin_elevation=SIN_ELEVATION,
        nodes_per_decade: int = NODES_PER_DECADE,
        check_points: int = CHECK_POINTS,
    ) -> None:
        if not isinstance(material, Material | Mixture):
            raise TypeError(
                "the reference material must be a Material or a Mixture (for a soil: soil.material)"
            )
        date = iso_date(date)
        if date is not None and (latitude_deg is None or longitude_deg is None or azimuth_deg is None):
            raise ValueError(
                "with a date the geomagnetic field needs latitude_deg, longitude_deg and azimuth_deg"
            )
        if date is None and azimuth_deg is not None:
            raise ValueError(
                "azimuth_deg has a meaning only with a date: without a field the table has no azimuth"
            )
        if not np.isfinite(altitude_m) or abs(altitude_m) > 9000:
            raise ValueError("altitude_m must lie within 9000 m of sea level")
        for value, limit, what in ((latitude_deg, 90, "latitude_deg"), (longitude_deg, 360, "longitude_deg")):
            if value is not None and not abs(value) <= limit:
                raise ValueError(f"{what} is out of range")
        sines = np.array(sin_elevation, dtype=float)
        if (
            sines.ndim != 1
            or sines.size < 2
            or not np.all(np.diff(sines) > 0)
            or sines[0] <= 0
            or sines[-1] > 1
        ):
            raise ValueError("sin_elevation: at least two increasing values in (0, 1]")
        if not (isinstance(nodes_per_decade, int) and nodes_per_decade >= 10):
            raise ValueError("nodes_per_decade must be an integer of at least 10")
        if not (isinstance(check_points, int) and check_points >= 0):
            raise ValueError("check_points must be a non-negative integer")
        inputs = {
            "build": BUILD,
            "mulder": mulder.config.VERSION,
            "reference": _reference(reference),
            "altitude_m": float(altitude_m),
            "latitude_deg": None if latitude_deg is None else float(latitude_deg),
            "longitude_deg": None if longitude_deg is None else float(longitude_deg),
            "azimuth_deg": None if azimuth_deg is None else float(azimuth_deg),
            "atmosphere": check_atmosphere(atmosphere),
            "date": date,
            "material": definition(material),
            "sin_elevation": sines.tolist(),
            "nodes_per_decade": nodes_per_decade,
            "k_min_gev": K_MIN_GEV,
            "k_max_gev": K_MAX_GEV,
            "cutoff_margin_gev": CUTOFF_MARGIN_GEV,
            "range_nodes_per_decade": RANGE_NODES_PER_DECADE,
            "opacity_max_g_cm2": OPACITY_MAX_G_CM2,
            "check": {
                "points": check_points,
                "seed": CHECK_SEED,
                "slab_kg_m3": DENSE_SLAB_KG_M3,
                "k_max_gev": CHECK_K_MAX_GEV,
                "limit": CHECK_LIMIT,
            },
        }
        self.key = hashlib.sha256(json.dumps(inputs, sort_keys=True).encode()).hexdigest()
        """Hash of the inputs: the name of the table in the cache."""
        self.cache_path = None
        """File of the table in the cache; None without a cache."""
        data = None
        if cache:
            directory = default_cache_dir() / "flux" if cache_dir is None else pathlib.Path(cache_dir)
            self.cache_path = directory / f"{self.key[:24]}.npz"
            data = _cached(self.cache_path, inputs)
        self.from_cache = data is not None
        """Whether the table was read from the cache instead of being built."""
        if data is None:
            data = _build(inputs, material)
        super().__init__(**data)
        if not self.from_cache:
            self.meta["check"] = _check(self, inputs, material)
            if self.cache_path is not None:
                self.cache_path.parent.mkdir(parents=True, exist_ok=True)
                self.save(self.cache_path)


def _cached(path: pathlib.Path, inputs: dict) -> dict | None:
    """Arguments of the table stored at `path`, if it is there, readable and made from `inputs`."""
    if not path.is_file():
        return None
    try:
        data = FluxTable.read(path)
        FluxTable(**data)  # a file that does not make a table is not a table
    except (OSError, TypeError, ValueError, KeyError, zipfile.BadZipFile):
        return None
    return data if data["meta"].get("inputs") == inputs else None


def _where(inputs: dict) -> dict:
    return {key: inputs[key] for key in ("latitude_deg", "longitude_deg", "azimuth_deg")}


def _build(inputs: dict, material: Material | Mixture) -> dict:
    """Open-sky flux from Mulder, and range and stopping power of the reference material."""
    physics = compiled(material)
    per_metre = mulder_density_kg_m3(material) / 10  # g/cm² in a metre of the material at its Mulder density
    altitude, sines = inputs["altitude_m"], np.array(inputs["sin_elevation"])
    sky = fluxmeter(atmosphere=inputs["atmosphere"], date=inputs["date"])
    cutoff = open_sky_cutoff(sky, altitude, sines, **_where(inputs))
    if cutoff.max() < K_MIN_GEV / 2:
        cutoff = np.zeros(sines.shape)  # below the first energy of the table: no cut-off
    elif cutoff.min() <= 0:
        raise RuntimeError("the reference flux has a cut-off at some elevations and none at others")
    first = CUTOFF_MARGIN_GEV if cutoff[0] > 0 else K_MIN_GEV
    reach = float(physics.range(K_MAX_GEV)) * per_metre + OPACITY_MAX_G_CM2
    top = 1.001 * float(
        physics.inverse_range(reach / per_metre)
    )  # open-sky energy of K_MAX_GEV under the largest opacity
    d = np.logspace(
        np.log10(first), np.log10(top), int(np.ceil(np.log10(top / first) * inputs["nodes_per_decade"])) + 1
    )
    energy = cutoff[:, None] + d
    state = states(
        np.full(energy.shape, altitude), np.broadcast_to(sines[:, None], energy.shape), **_where(inputs)
    )
    flux = np.asarray(sky.flux(energy=energy.ravel(), **{k: v.ravel() for k, v in state.items()})).reshape(
        energy.shape
    )
    if not np.all(flux > 0):
        raise RuntimeError("Mulder gave no flux at an energy above the cut-off")
    highest = 1.001 * (cutoff.max() + top)
    k = np.logspace(
        np.log10(K_MIN_GEV),
        np.log10(highest),
        int(np.ceil(np.log10(highest / K_MIN_GEV) * RANGE_NODES_PER_DECADE)) + 1,
    )
    return {
        "altitude_m": altitude,
        "sin_elevation": sines,
        "cutoff_gev": cutoff,
        "energy_above_cutoff_gev": d,
        "flux_m2_s_sr_gev": flux,
        "range_energy_gev": k,
        "range_g_cm2": np.asarray(physics.range(k)) * per_metre,
        "stopping_power_gev_cm2_g": np.asarray(physics.stopping_power(k)) / per_metre,
        "opacity_max_g_cm2": OPACITY_MAX_G_CM2,
        "meta": {"inputs": inputs, "material_name": material.name, "mulder_material": register(material)},
    }


def check_sample(table: FluxTable, points: int, seed: int, k_max_gev: float = CHECK_K_MAX_GEV):
    """Random energies up to `k_max_gev`, sines of elevation and opacities over the table.

    Sines are uniform, opacities log-uniform from 1 g/cm², energies log-uniform from the lowest
    one the table holds in each direction.
    """
    rng = np.random.default_rng(seed)
    sine = rng.uniform(table.sin_elevation[0], table.sin_elevation[-1], points)
    opacity = np.exp(rng.uniform(0.0, np.log(table.opacity_max_g_cm2), points))
    low = table.energy_min(sine, opacity)
    energy = low * np.exp(rng.uniform(0.0, 1.0, points) * np.log(k_max_gev / low))
    return energy, sine, opacity


def slab_flux(inputs: dict, material: Material | Mixture, energy, sine, opacity) -> np.ndarray:
    """Differential flux from Mulder under a dense slab of the material whose top is at the table altitude."""
    slab = mulder.Layer(inputs["altitude_m"], material=register(material), density=DENSE_SLAB_KG_M3)
    meter = fluxmeter(slab, atmosphere=inputs["atmosphere"], date=inputs["date"])
    depth = np.asarray(opacity) * 10 / DENSE_SLAB_KG_M3 * np.asarray(sine)  # m below the top
    return np.asarray(
        meter.flux(energy=energy, **states(inputs["altitude_m"] - depth, sine, **_where(inputs)))
    )


def _check(table: FluxTable, inputs: dict, material: Material | Mixture) -> dict:
    """Compare the table with Mulder under a dense slab; refuse a table that is off."""
    points = inputs["check"]["points"]
    if points == 0:
        return {"points": 0}
    energy, sine, opacity = check_sample(table, points, inputs["check"]["seed"], inputs["check"]["k_max_gev"])
    ratio = table.differential(energy, sine, opacity) / slab_flux(inputs, material, energy, sine, opacity) - 1
    worst = int(np.argmax(np.abs(ratio)))
    result = {
        "points": points,
        "rms": float(np.sqrt(np.mean(ratio**2))),
        "max": float(np.abs(ratio).max()),
        "max_at": {
            "kinetic_energy_gev": float(energy[worst]),
            "sin_elevation": float(sine[worst]),
            "opacity_g_cm2": float(opacity[worst]),
        },
    }
    if not result["max"] <= inputs["check"]["limit"]:
        raise RuntimeError(
            f"the table differs from Mulder under a slab by {100 * result['max']:.2f} %"
            f" at {result['max_at']}: continuous transport does not carry the open-sky flux"
            " through this material as the table assumes"
        )
    return result
