"""Comparison with measured data.

Tracks -> counts in the same DirectionGrid with the same cuts; data/model ratio, pulls, chi2;
fit of densities by region with nuisance parameters. File reading lives in `io`.

Version 0.1 has `bin_tracks` and `BinnedTracks` (`docs/notes/tracks.md`); the rest comes later.

Layer 4: may import only `muoscene` modules of layer <= 4 (see `muoscene._layers`).
"""

from .binning import BinnedTracks, bin_tracks

__all__ = ["BinnedTracks", "bin_tracks"]
