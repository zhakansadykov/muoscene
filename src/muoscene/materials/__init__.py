"""Materials.

Compositions by formula, element mass fractions or oxide analysis, mixtures, and soils given by
dry density and water content (`composition`); a library of soils, concrete, standard rock,
water and air with sources (`library`); registration with Mulder under names derived from the
definition (`backend`); muon ranges, equivalent-opacity weights and equivalent opacity
(`range_tables`). Choices and their size: `docs/notes/materials.md`.

Layer 1: may import only `muoscene` modules of layer <= 1 (see `muoscene._layers`).
"""

from .backend import mulder_name, register
from .composition import Material, Mixture, Soil, berger_seltzer_mean_excitation
from .library import (
    AIR,
    CALCITE,
    CONCRETE,
    LOESS_SKELETON,
    QUARTZ,
    SILICATE_SKELETON,
    STANDARD_ROCK,
    WATER,
    loam,
    loess,
    sand,
    soil,
)
from .range_tables import csda_range, energy_to_cross, equivalent_opacity, weight

__all__ = [
    "AIR",
    "CALCITE",
    "CONCRETE",
    "LOESS_SKELETON",
    "QUARTZ",
    "SILICATE_SKELETON",
    "STANDARD_ROCK",
    "WATER",
    "Material",
    "Mixture",
    "Soil",
    "berger_seltzer_mean_excitation",
    "csda_range",
    "energy_to_cross",
    "equivalent_opacity",
    "loam",
    "loess",
    "mulder_name",
    "register",
    "sand",
    "soil",
    "weight",
]
