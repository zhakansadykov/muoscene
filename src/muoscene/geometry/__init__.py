"""Scene geometry.

A scene is described by bodies: the ground as an extruded cross-section or a closed mesh, layers,
boxes, cylinders, pile grids, niches, spherical anomalies. The builder turns it into a Calzone
geometry and Mulder scans it. `Scene.scan` returns a `PathTable`, the path length in every body
for every direction; opacity follows from it for any set of densities, so density hypotheses
need no new scan. The builder owns the unit conversion (metres to Calzone centimetres), the
placement rules (positions relative to mother volumes, subtraction of cuts, splitting of
structures at layer boundaries) and the checks (an anomaly lies within one host, no overlaps).
The air of a scene has the density of the atmosphere at the altitude of the site. Closed
meshes, terrain models and a voxel overlay with exact path lengths (coverage maps, inversion)
come with version 0.2. Design and numbers: `docs/notes/geometry.md`.

Knows materials by name only; no energies, no fluxes. Calzone is imported only in this package.

Layer 1: may import only `muoscene` modules of layer <= 1 (see `muoscene._layers`).
"""

from .pathtable import BodyInfo, PathTable
from .primitives import Anomaly, Box, Cut, Cylinder, Ground, Layer, Sphere, Structure, pile_grid
from .scene import Scene, air_density_g_cm3

__all__ = [
    "Anomaly",
    "BodyInfo",
    "Box",
    "Cut",
    "Cylinder",
    "Ground",
    "Layer",
    "PathTable",
    "Scene",
    "Sphere",
    "Structure",
    "air_density_g_cm3",
    "pile_grid",
]
