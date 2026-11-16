"""Muon flux at the detector as a function of what lies along the line of sight.

`FluxTable` holds the open-sky differential flux at the detector altitude on a grid of kinetic
energy and sine of elevation, carries it under an opacity of its reference material by the CSDA
range of that material, and integrates it over a detector energy response (`EnergyResponse`).
`MulderTable` builds it with Mulder for an atmosphere model and optionally the geomagnetic field
of a date, checks it against Mulder under a dense slab and caches it on disk; the reference is
GCCLY15 (Guan et al. 2015) at sea level, whose altitude cut-off the table records and enforces;
a tabulated reference at the site altitude, needed for muons below the cut-off, comes with
roadmap 0.1-05b. `LegacyTable` reproduces the preliminary calculation (sea level, a + bE) for
regression tests. Design and numbers: `docs/notes/flux.md`.

Layer 1: may import only `muoscene` modules of layer <= 1 (see `muoscene._layers`).
"""

from .atmosphere import atmosphere_models
from .base import OPACITY_MAX_G_CM2, EnergyResponse, FluxTable
from .parametrizations import LEGACY_RESPONSE, LegacyTable
from .tables import MulderTable, default_cache_dir

__all__ = [
    "LEGACY_RESPONSE",
    "OPACITY_MAX_G_CM2",
    "EnergyResponse",
    "FluxTable",
    "LegacyTable",
    "MulderTable",
    "atmosphere_models",
    "default_cache_dir",
]
