"""flux.atmosphere — atmosphere model, date and geomagnetic field passed to Mulder; altitude cut-off.

The reference flux of version 0.1.0 is GCCLY15 (Guan et al. 2015) at sea level. Mulder carries
it through the air between sea level and the site: a model of the density profile, and, when a
date is given, the geomagnetic field (IGRF) of that date, which bends the muons in that air.
Above sea level this defines the flux only above the altitude cut-off, the energy lost in that
air (`open_sky_cutoff`); below sea level there is no cut-off.

Units at this boundary: sine of elevation to degrees of elevation; energies are kinetic, in GeV,
on both sides; fluxes are m⁻² s⁻¹ sr⁻¹ GeV⁻¹ on both sides.
"""

from __future__ import annotations

import datetime

import mulder
import numpy as np

REFERENCE = "GCCLY15"
"""Reference flux of Mulder used by version 0.1.0: Guan et al. (2015), at sea level."""
CUTOFF_SEARCH_GEV = (1.0e-7, 100.0)
"""Kinetic energies between which the cut-off is looked for."""
CUTOFF_STEPS = 70
"""Bisection steps in the logarithm of energy: the cut-off is found to 1e-15 of its value."""


def atmosphere_models() -> tuple[str, ...]:
    """Names of the density profiles Mulder knows, e.g. "us-standard" and "midlatitude-winter"."""
    return tuple(mulder.Atmosphere.models)


def check_atmosphere(name: str) -> str:
    """The name of an atmosphere model, or ValueError with the list of those available."""
    if name not in atmosphere_models():
        raise ValueError(f"unknown atmosphere {name!r}; Mulder has {', '.join(atmosphere_models())}")
    return name


def iso_date(date) -> str | None:
    """A date for the geomagnetic field as "YYYY-MM-DD"; None stays None (no field)."""
    if date is None:
        return None
    if isinstance(date, datetime.date):
        return date.isoformat()[:10]
    try:
        return datetime.date.fromisoformat(str(date)).isoformat()
    except ValueError as error:
        raise ValueError(f"the date of the geomagnetic field must be YYYY-MM-DD, not {date!r}") from error


def fluxmeter(*layers, atmosphere: str, date: str | None = None):
    """A Mulder fluxmeter in continuous mode over the given layers, with the GCCLY15 reference."""
    options = {"atmosphere": check_atmosphere(atmosphere)}
    if date is not None:
        options["date"] = iso_date(date)
    return mulder.Fluxmeter(*layers, **options)


def states(
    altitude_m,
    sin_elevation,
    *,
    latitude_deg: float | None = None,
    longitude_deg: float | None = None,
    azimuth_deg: float | None = None,
) -> dict:
    """Arguments of `Fluxmeter.flux` for lines of sight at the given altitudes (m) and sines of elevation.

    The azimuth is measured from north, clockwise, as in the conventions. Latitude and longitude
    default to those of Mulder (0, 0); without a geomagnetic field they change the flux by
    2e-6 (`docs/notes/flux.md`).
    """
    altitude, sine = np.broadcast_arrays(np.asarray(altitude_m, float), np.asarray(sin_elevation, float))
    out = {"altitude": altitude, "elevation": np.degrees(np.arcsin(sine))}
    for key, value in (("latitude", latitude_deg), ("longitude", longitude_deg), ("azimuth", azimuth_deg)):
        if value is not None:
            out[key] = np.full(altitude.shape, float(value))
    return out


def open_sky_cutoff(meter, altitude_m: float, sin_elevation, **where) -> np.ndarray:
    """Lowest kinetic energy (GeV) at which the fluxmeter gives a flux, per sine of elevation.

    0 where the flux is positive down to 1e-7 GeV: at and below the altitude of the reference.
    `where` are the latitude, longitude and azimuth of `states`.
    """
    sine = np.atleast_1d(np.asarray(sin_elevation, dtype=float))
    state = states(np.full(sine.shape, float(altitude_m)), sine, **where)

    def positive(ln_energy):
        return np.asarray(meter.flux(energy=np.exp(ln_energy), **state)) > 0

    low = np.full(sine.shape, np.log(CUTOFF_SEARCH_GEV[0]))
    high = np.full(sine.shape, np.log(CUTOFF_SEARCH_GEV[1]))
    none = positive(low)
    if not np.all(positive(high)):
        raise RuntimeError(f"no flux at {CUTOFF_SEARCH_GEV[1]:g} GeV: the site is outside what Mulder models")
    for _ in range(CUTOFF_STEPS):
        middle = 0.5 * (low + high)
        above = positive(middle)
        high = np.where(above, middle, high)
        low = np.where(above, low, middle)
    return np.where(none, 0.0, np.exp(high))
