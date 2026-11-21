"""From paths to flux.

Fast path: the equivalent opacity of a `PathTable` under a density hypothesis and a `FluxTable`
lookup. Direct path: Mulder continuous flux on the built scene, used to check the fast path on a
sample of directions in every run, and Mulder discrete flux, which includes scattering and
energy-loss fluctuations. A scattering kernel calibrated on discrete runs serves planning loops.

Layer 2: may import only `muoscene` modules of layer <= 2 (see `muoscene._layers`).
"""

from .direct import DirectCheck, compare_with_direct, direct_flux, energy_nodes, seed_of
from .table import (
    BELOW_CUTOFF,
    DEFINED,
    LOW_ELEVATION,
    NOT_ACCEPTED,
    STATUS_MEANINGS,
    TOO_OPAQUE,
    EquivalentOpacity,
    equivalent_opacity,
    flux_domain,
)

__all__ = [
    "BELOW_CUTOFF",
    "DEFINED",
    "LOW_ELEVATION",
    "NOT_ACCEPTED",
    "STATUS_MEANINGS",
    "TOO_OPAQUE",
    "DirectCheck",
    "EquivalentOpacity",
    "compare_with_direct",
    "direct_flux",
    "energy_nodes",
    "equivalent_opacity",
    "flux_domain",
    "seed_of",
]
