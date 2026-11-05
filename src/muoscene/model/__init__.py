"""Forward model.

Scene x flux x detector x exposure -> ExpectedMap (xarray Dataset: expected counts per bin,
through-mass and reference components, hypotheses, systematic variants, run metadata).
Hypothesis sets (void, decompaction, none) and nuisance parameters. The only place where
geometry, transport, flux and detectors meet.

Layer 3: may import only `muoscene` modules of layer <= 3 (see `muoscene._layers`).
"""
