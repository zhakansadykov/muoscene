"""Atmospheric muon flux at the surface.

FluxModel interface I(>E, cos theta); parametrisations (Guan et al. 2015, Reyna 2006,
Bugaev et al. 1998, Gaisser); altitude correction; seasonal correction with optional MCEq backend;
precomputed integral tables with interpolation and caching. Model spread is a systematic.

Layer 1: may import only `muoscene` modules of layer <= 1 (see `muoscene._layers`).
"""
