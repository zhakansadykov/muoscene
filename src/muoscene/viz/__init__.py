"""Visualisation.

Angular maps of a detector — an `ExpectedMap` or any array on its grid — and the flux of a table
against the zenith angle, in one style for papers and reports, with matplotlib (`maps`,
`curves`, `style`; roadmap 0.1-12a, `docs/notes/figures.md`). Later: the scene in plan and in
section (0.1-12) and a 3D scene with PyVista exported to HTML (`scene3d`, 0.1-12b). The core never
imports this package, and this package imports matplotlib inside its functions.

Layer 5: may import only `muoscene` modules of layer <= 5 (see `muoscene._layers`).
"""
