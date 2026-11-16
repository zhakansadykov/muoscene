"""materials.backend — materials registered with Mulder under names derived from their definition.

Mulder keeps one process-wide registry of materials and refuses a second definition under a
name already taken. muoscene therefore names every material by a hash of its full definition:
"loam, w = 12 %" and "loam, w = 18 %" never collide, and defining the same material twice is a
no-op. Mulder caches compiled physics tables on disk by the same kind of hash
(`mulder.config.CACHE`), so a material compiles once per machine: about 3 s for a new one.

Units at this boundary: g/cm³ → kg/m³, eV → GeV.
"""

from __future__ import annotations

import hashlib
import re
from functools import cache

import mulder
import mulder.materials as mm

from .composition import Material, Mixture

PREFIX = "ms"


def _number(x: float) -> str:
    """12 significant digits: definitions that differ only by round-off share one Mulder material."""
    return f"{x:.12g}"


@cache
def definition(material: Material | Mixture) -> str:
    """Canonical text of what Mulder will be given; equal definitions give equal text."""
    if isinstance(material, Material):
        fractions = ",".join(f"{s}:{_number(w)}" for s, w in material.mass_fractions)
        i = "mulder" if material.i_ev is None else _number(material.i_ev)
        return f"material[{fractions}|rho={_number(material.density_g_cm3)}|I={i}]"
    parts = ",".join(sorted(f"{mulder_name(m)}:{_number(w)}" for m, w in material.parts))
    return f"mixture[{parts}]"


@cache
def mulder_name(material: Material | Mixture) -> str:
    """Name under which Mulder knows this material: prefix, readable label and a hash.

    Only letters, digits and hyphens: Mulder writes the names as bare TOML keys into its disk
    cache, and a name with a dot or a space never matches its cache entry again (each process
    then recompiles the tables, 2-11 s per material).
    """
    digest = hashlib.sha256(definition(material).encode()).hexdigest()[:12]
    label = re.sub(r"[^0-9a-z]+", "-", material.name.lower()).strip("-")[:32] or "material"
    return f"{PREFIX}-{label}-{digest}"


@cache
def register(material: Material | Mixture) -> str:
    """Define the material in Mulder if it is not there yet; returns its Mulder name."""
    name = mulder_name(material)
    if isinstance(material, Material):
        if name not in mm.Material.all():
            mm.Material.define(
                name,
                composition=dict(material.mass_fractions),
                density=1e3 * material.density_g_cm3,
                I=None if material.i_ev is None else 1e-9 * material.i_ev,
            )
    else:
        parts = {register(m): w for m, w in material.parts}
        if name not in mm.Composite.all():
            mm.Composite.define(name, composition=parts)
    return name


def mulder_density_kg_m3(material: Material | Mixture) -> float:
    """Density Mulder attaches to the material, kg/m³ (a composite's is that of its parts, no pores)."""
    name = register(material)
    if isinstance(material, Material):
        return float(mm.Material(name).density)
    return float(mm.Composite(name).density)


@cache
def _compile(name: str):
    return mulder.Physics().compile(name)


def compiled(material: Material | Mixture):
    """Mulder physics tables of the material (`mulder.CompiledMaterial`)."""
    return _compile(register(material))
