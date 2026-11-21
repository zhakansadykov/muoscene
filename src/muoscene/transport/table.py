"""transport.table — fast path: equivalent opacity of a path and where a flux table has a flux.

A flux table takes the opacity of its reference material (`muoscene.flux`). A line of sight
through a scene crosses several materials, and a `PathTable` knows how much of each but not in
which order. Two decisions of 21 November 2026 turn those opacities into the one number the
table is entered with (`docs/notes/forward.md`):

Weights. The opacity X_m of material m counts as w_m X_m of the reference material, with

    w_m = (R_ref(K_top) − R_ref(K_min)) / (R_m(K_top) − R_m(K_min)),    R_ref(K_top) = R_ref(K_min) + X_eq,

R the CSDA range, K_min the first energy of the detector response and K_top the energy a muon
needs above the path to arrive with K_min. X_eq = Σ w_m X_m is found by iteration. It is the
weight of the conventions, R_ref(K) / R_m(K), for a detector with a threshold: exact for a path
through one material, and the same whichever material the table has as its reference (2e-7).
On 3000 g/cm² of loam under 1840 g/cm² of concrete it is +0.023 % in flux from the chain of
ranges in the true order; that chain in the wrong order is +0.048 %, a weight taken at one
energy +0.17, +0.06 and −0.08 % at 1, 7 and 100 GeV, and no weight at all −0.17 %
(`tests/unit/test_transport_table.py`).

Displaced air. A table holds the open sky at the altitude of the detector, so along a line of
sight under soil it has air where the soil is. That air, ρ_air times the length of the path in
solid bodies, weighted like any other material, is taken out of the opacity the table is
entered with. Without it the flux under 20 m of soil is 0.074–0.078 % low; with it 0.012–0.021 %
high, the decay of muons inside the soil, which a table does not have. A table at the altitude
of the surface would be 0.03–0.05 % high under soil (`docs/notes/flux.md`) and 0.31–0.39 %
high in the open sky seen from a niche 20 m below the brow.

`flux_domain` tells, line by line, whether a table has a flux for a detector response, and if
not, why: the forward model uses it to refuse a run or to mark the bins it has no number for.

Units: opacity in g/cm², lengths in m, densities in g/cm³, kinetic energy in GeV.

Layer 2.
"""

from __future__ import annotations

import functools
from collections.abc import Mapping
from dataclasses import dataclass

import numpy as np

from ..flux import LEGACY_RESPONSE, EnergyResponse, FluxTable, LegacyTable
from ..materials import AIR, Material, Mixture, csda_range

DEFINED = 0
"""The bin or line of sight has a number."""
NOT_ACCEPTED = 1
"""No accepted line of sight: outside the acceptance circle, behind the detector or below the horizon."""
BELOW_CUTOFF = 2
"""The response starts below the lowest energy the table holds there: the cut-off of the reference flux."""
LOW_ELEVATION = 3
"""The elevation is outside the flux table."""
TOO_OPAQUE = 4
"""The opacity is beyond the flux table: continuous transport is not valid there."""
STATUS_MEANINGS = (
    "defined",
    "not_accepted",
    "below_cutoff_of_reference_flux",
    "elevation_outside_flux_table",
    "opacity_beyond_flux_table",
)
"""What the status codes 0 to 4 mean, in the order of the codes."""

EDGE = 1e-9
"""Relative slack at the edges of a table, as in `muoscene.flux`."""
LEGACY_SIN_ELEVATION = (0.25, 1.0)
"""Sines of elevation a `LegacyTable` serves."""
SPAN_MIN_G_CM2 = 1.0e-3
"""Weights are taken over at least this much opacity: below it the difference of two ranges is round-off."""
TOLERANCE = 1.0e-9
"""The iteration stops when the equivalent opacity moves by less than this, relatively.

Each pass gains three to four digits, so the result is closer than that; below 1e-11 the
difference of two ranges over a thin path is round-off and the passes would alternate."""
ITERATIONS_MAX = 60


@dataclass(frozen=True)
class EquivalentOpacity:
    """Equivalent opacity of lines of sight in the reference material of a flux table.

    matter_g_cm2          Σ w_m X_m over the solid materials, g/cm²: the equivalent opacity of the
                          conventions, what a map reports
    displaced_air_g_cm2   opacity of the air those bodies displace, in the reference material, g/cm²
    lookup_g_cm2          what the flux table is entered with: the matter less the displaced air,
                          not below zero, g/cm²
    iterations            passes the weights took to settle
    """

    matter_g_cm2: np.ndarray
    displaced_air_g_cm2: np.ndarray
    lookup_g_cm2: np.ndarray
    iterations: int


@functools.lru_cache(maxsize=64)
def _ln_range(material: Material | Mixture, nodes: bytes) -> np.ndarray:
    """Logarithm of the CSDA range of a material at the energy nodes of a flux table."""
    return np.log(np.asarray(csda_range(material, np.frombuffer(nodes, dtype=float)), dtype=float))


def equivalent_opacity(
    table: FluxTable,
    opacity_by_material_g_cm2: Mapping[str, np.ndarray],
    materials: Mapping[str, Material | Mixture],
    k_min_gev: float,
    *,
    solid_length_m=None,
    air_density_g_cm3: float | None = None,
) -> EquivalentOpacity:
    """Opacity of the reference material of `table` that stands for a path through several materials.

    table                       the flux table; its range table is the reference material
    opacity_by_material_g_cm2   opacity in every solid material along every line of sight, g/cm²,
                                by the name of the material in Mulder
                                (`PathTable.opacity_by_material_g_cm2`)
    materials                   the materials by the same names
    k_min_gev                   first kinetic energy of the detector response, GeV
    solid_length_m              length of the path in solid bodies along every line, m, and
    air_density_g_cm3           density of the air of the scene, g/cm³: with both, the air the
                                bodies displace is taken out of `lookup_g_cm2`; with neither,
                                `lookup_g_cm2` is the matter

    The material whose name is that of the reference of the table has the weight one. A name
    without a material, a negative opacity and a threshold outside the range table are refused.
    """
    if not isinstance(table, FluxTable):
        raise TypeError("the equivalent opacity is that of the reference material of a FluxTable")
    nodes = table.range_energy_gev
    if not nodes[0] <= k_min_gev <= nodes[-1]:
        raise ValueError(f"k_min_gev must lie within the range table, {nodes[0]:g} to {nodes[-1]:g} GeV")
    if (solid_length_m is None) != (air_density_g_cm3 is None):
        raise ValueError("give both solid_length_m and air_density_g_cm3, or neither")
    ln_e, ln_ref = np.log(nodes), np.log(table.range_g_cm2)
    reference = table.meta.get("mulder_material")
    names = list(opacity_by_material_g_cm2)
    curves = []
    for name in names:
        if name == reference:
            curves.append(ln_ref)
        elif name in materials:
            curves.append(_ln_range(materials[name], nodes.tobytes()))
        else:
            raise ValueError(f"no material given for {name!r}")
    x = [np.asarray(opacity_by_material_g_cm2[name], dtype=float) for name in names]
    length = None if solid_length_m is None else np.asarray(solid_length_m, dtype=float)
    shape = np.broadcast_shapes(*(a.shape for a in x), *(() if length is None else (length.shape,)))
    x = [np.broadcast_to(a, shape) for a in x]
    if any(not np.all(np.isfinite(a) & (a >= 0)) for a in x):
        raise ValueError("opacities must be finite and non-negative")

    ln_k = np.log(k_min_gev)
    r0_ref = np.exp(np.interp(ln_k, ln_e, ln_ref))
    r_top_max = table.range_g_cm2[-1]

    def weights(equivalent, curve_list):
        """Range-difference weights of the materials for the equivalent opacity reached so far."""
        span = np.minimum(np.maximum(equivalent, SPAN_MIN_G_CM2), r_top_max - r0_ref)
        ln_top = np.interp(np.log(r0_ref + span), ln_ref, ln_e)
        return [
            1.0
            if curve is ln_ref
            else span / (np.exp(np.interp(ln_top, ln_e, curve)) - np.exp(np.interp(ln_k, ln_e, curve)))
            for curve in curve_list
        ]

    matter = np.zeros(shape)
    for a in x:
        matter = matter + a
    iterations = 0
    if any(curve is not ln_ref for curve in curves):
        for iterations in range(1, ITERATIONS_MAX + 1):  # noqa: B007
            new = np.zeros(shape)
            for w, a in zip(weights(matter, curves), x, strict=True):
                new = new + w * a
            moved = float(np.max(np.abs(new - matter) / np.maximum(new, SPAN_MIN_G_CM2), initial=0.0))
            matter = new
            if moved < TOLERANCE:
                break
        else:
            raise RuntimeError("the weights of the equivalent opacity did not settle")

    if length is None:
        air = np.zeros(shape)
    else:
        if not np.all(np.isfinite(length) & (length >= 0)):
            raise ValueError("solid_length_m must be finite and non-negative")
        if not (np.isfinite(air_density_g_cm3) and air_density_g_cm3 >= 0):
            raise ValueError("air_density_g_cm3 must be finite and non-negative")
        (w_air,) = weights(matter, [_ln_range(AIR, nodes.tobytes())])
        air = w_air * air_density_g_cm3 * np.broadcast_to(length, shape) * 100.0
    arrays = [matter, air, np.maximum(matter - air, 0.0)]
    for array in arrays:
        array.flags.writeable = False
    return EquivalentOpacity(*arrays, iterations=iterations)


def flux_domain(table, response, sin_elevation, opacity_g_cm2) -> np.ndarray:
    """Whether a flux table has a flux for a detector response along lines of sight, and if not, why.

    table            a `FluxTable`, or a `LegacyTable` with `LEGACY_RESPONSE`
    response         the energy response of the detector
    sin_elevation    sine of the elevation of every line of sight
    opacity_g_cm2    opacity the table is entered with along each, g/cm²

    Returns one code per line, an int8 array: `DEFINED`, `LOW_ELEVATION` (outside the elevations
    of the table), `TOO_OPAQUE` (beyond its largest opacity) or `BELOW_CUTOFF` (the response
    starts below the lowest energy the table holds there, which with a sea-level reference above
    sea level is its cut-off). `table.integral` succeeds on the lines that are `DEFINED` and on
    no others. A negative opacity and a response the table cannot take at all are refused.
    """
    s, x = np.broadcast_arrays(np.asarray(sin_elevation, dtype=float), np.asarray(opacity_g_cm2, dtype=float))
    if not (np.all(np.isfinite(s)) and np.all(np.isfinite(x)) and np.all(x >= 0)):
        raise ValueError("sines of elevation must be finite and opacities finite and non-negative")
    if isinstance(table, FluxTable):
        if not isinstance(response, EnergyResponse):
            raise TypeError("a FluxTable takes an EnergyResponse")
        if response.kinetic_energy_gev[-1] > table.k_max_gev * (1 + EDGE):
            raise ValueError(
                f"the energy response reaches {response.kinetic_energy_gev[-1]:g} GeV, above"
                f" {table.k_max_gev:.4g} GeV, the highest energy of the table"
            )
        low, high = float(table.sin_elevation[0]), float(table.sin_elevation[-1])
    elif isinstance(table, LegacyTable):
        if response is not LEGACY_RESPONSE:
            raise TypeError("a LegacyTable takes LEGACY_RESPONSE")
        low, high = LEGACY_SIN_ELEVATION
    else:
        raise TypeError("the domain is known for a FluxTable and for a LegacyTable")
    code = np.full(s.shape, DEFINED, dtype=np.int8)
    code[(s < low * (1 - EDGE)) | (s > high * (1 + EDGE))] = LOW_ELEVATION
    code[(code == DEFINED) & (x > table.opacity_max_g_cm2 * (1 + EDGE))] = TOO_OPAQUE
    if isinstance(table, FluxTable):
        inside = code == DEFINED
        if np.any(inside):
            lowest = np.atleast_1d(table.energy_min(s[inside], x[inside]))
            below = np.zeros(s.shape, dtype=bool)
            below[inside] = response.k_min_gev < lowest * (1 - EDGE)
            code[below] = BELOW_CUTOFF
    return code
