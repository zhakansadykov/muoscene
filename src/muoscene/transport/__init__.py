"""From paths to flux.

Fast path: the equivalent opacity of a `PathTable` under a density hypothesis and a `FluxTable`
lookup. Direct path: Mulder continuous flux on the built scene, used to check the fast path on a
sample of directions in every run, and Mulder discrete flux, which includes scattering and
energy-loss fluctuations. A scattering kernel calibrated on discrete runs serves planning loops.

Layer 2: may import only `muoscene` modules of layer <= 2 (see `muoscene._layers`).
"""
