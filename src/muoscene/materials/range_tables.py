"""materials.range_tables — muon ranges from Mulder and equivalent opacity.

The CSDA range R_m(K) of a muon of kinetic energy K in material m, in g/cm², is what decides the
energy needed to cross a given opacity. Two materials compare through their ranges:

- weight w_m(K) = R_ref(K) / R_m(K): one g/cm² of m counts as w_m g/cm² of the reference for a
  muon that stops in it from energy K (0.995 for concrete against loam at 7 GeV);
- equivalent opacity X_eq = R_ref(R_m⁻¹(X)): the opacity of the reference that needs the same
  energy as X g/cm² of m. It is exact for a path through one material; how paths through several
  materials are combined is the business of the forward model (roadmap 0.1-08).

Energies are kinetic, in GeV; opacities in g/cm². Energies from 0 to `K_MAX_GEV` are served;
a negative, larger or undefined value is refused, not extrapolated.
"""

from __future__ import annotations

import numpy as np

from .backend import compiled, mulder_density_kg_m3
from .composition import Material, Mixture

AnyMaterial = Material | Mixture

K_MAX_GEV = 1e6
"""Upper end of the energies served: 1.7e6 g/cm² of standard rock, far beyond any scene."""


def _energies(kinetic_energy_gev):
    k = np.asarray(kinetic_energy_gev, dtype=float)
    if not np.all((k >= 0) & (k <= K_MAX_GEV)):
        raise ValueError(f"kinetic energy must be within 0 to {K_MAX_GEV:g} GeV")
    return k


def _out(value, like):
    value = np.asarray(value, dtype=float)
    return float(value) if np.ndim(like) == 0 else value


def csda_range(material: AnyMaterial, kinetic_energy_gev):
    """CSDA range in g/cm² at kinetic energy K (GeV)."""
    k = _energies(kinetic_energy_gev)
    metres = np.asarray(compiled(material).range(k), dtype=float)
    return _out(metres * mulder_density_kg_m3(material) / 10, kinetic_energy_gev)


def energy_to_cross(material: AnyMaterial, opacity_g_cm2):
    """Kinetic energy (GeV) a muon needs to cross `opacity_g_cm2` of the material (inverse range)."""
    x = np.asarray(opacity_g_cm2, dtype=float)
    if not np.all(x >= 0):
        raise ValueError("opacity must be non-negative")
    if np.any(x > csda_range(material, K_MAX_GEV)):
        raise ValueError(f"opacity beyond the range of a {K_MAX_GEV:g} GeV muon")
    metres = x * 10 / mulder_density_kg_m3(material)
    return _out(compiled(material).inverse_range(metres), opacity_g_cm2)


def weight(material: AnyMaterial, reference: AnyMaterial, kinetic_energy_gev):
    """w_m(K) = R_ref(K) / R_m(K), for K > 0."""
    k = _energies(kinetic_energy_gev)
    if np.any(k == 0):
        raise ValueError("the weight is not defined at zero energy")
    ratio = np.asarray(csda_range(reference, k)) / np.asarray(csda_range(material, k))
    return _out(ratio, kinetic_energy_gev)


def equivalent_opacity(material: AnyMaterial, reference: AnyMaterial, opacity_g_cm2):
    """Opacity of the reference (g/cm²) that needs the same energy as `opacity_g_cm2` of the material."""
    k = np.asarray(energy_to_cross(material, opacity_g_cm2))
    return _out(csda_range(reference, k), opacity_g_cm2)
