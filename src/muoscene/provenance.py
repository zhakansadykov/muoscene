"""Run folders.

`runs/<date>_<name>_<hash>/` with the original and resolved configuration, hashes of input
files, package version and git commit, environment and log. Nothing without a run folder goes
into a paper.

Layer 5: may import only `muoscene` modules of layer <= 5 (see `muoscene._layers`).
"""
