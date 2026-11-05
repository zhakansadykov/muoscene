"""Scene geometry and ray tracing.

Ray batches and per-ray intervals; analytic primitives (box, cylinder, sphere, prism, pile grid);
closed triangle meshes (terrain, canyon wall, structures); height-field surfaces for layer
boundaries; voxel overlay with exact path lengths (Amanatides-Woo). The Scene resolves overlapping
bodies by priority (cuts 100 > structures 80 > replacing anomalies 60 > layers 40 > ground 20,
air elsewhere) and returns a PathTable: path length per material and body, opacity X in g/cm^2,
and optionally a sparse ray-by-voxel length matrix. Density-modifier anomalies and voxels scale
density instead of replacing material.

Knows materials by name only; no energies, no fluxes.

Layer 1: may import only `muoscene` modules of layer <= 1 (see `muoscene._layers`).
"""
