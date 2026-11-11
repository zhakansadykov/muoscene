"""Orchestration.

Builds the site, materials, scene, flux model and detectors from a resolved configuration and
runs the stages (check, scene, plan, synth, compare, invert). The command-line interface and
Python users call the same functions here; `muoscene.cli` only parses arguments. Run folders are
written through `muoscene.provenance`.

Layer 5: may import only `muoscene` modules of layer <= 5 (see `muoscene._layers`).
"""
