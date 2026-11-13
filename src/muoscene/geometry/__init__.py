"""Scene geometry.

A scene is described by bodies: the ground as an extruded cross-section or a closed mesh, layers,
boxes, cylinders, pile grids, niches, spherical anomalies. The builder turns it into a Calzone
geometry and Mulder scans it. `Scene.scan` returns a `PathTable`, the path length in every body
for every direction; opacity follows from it for any set of densities, so density hypotheses
need no new scan. The builder owns the unit conversion (metres to Calzone centimetres), the
placement rules (positions relative to mother volumes, subtraction of cuts, splitting of
structures at layer boundaries) and the checks (an anomaly lies within one host, no overlaps).
A voxel overlay with exact path lengths serves coverage maps and inversion.

Knows materials by name only; no energies, no fluxes. Calzone is imported only in this package.

Layer 1: may import only `muoscene` modules of layer <= 1 (see `muoscene._layers`).
"""
