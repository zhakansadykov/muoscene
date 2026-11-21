"""Detectors.

Detector interface: pose in ENU; direction grid in the detector frame (tan theta_x, tan theta_y,
solid angles) and its transform to ENU; acceptance, efficiency eps(direction, time), blur kernel,
background. Emulsion stack: cone, angular efficiency, latent-image fading from temperature and
humidity logs, track-density budget (exposure plus transport and storage), number of films and
scanning cost. Straw telescope (level-1 response). Generic two-plane telescope. `TrackSet`, the
container of measured and synthetic tracks, is defined here so that `synth`, `data` and `io` share it.

Version 0.1 has `Pose`, `DirectionGrid`, `LinesOfSight`, the `Detector` interface and
`EmulsionStack` with its track density (`docs/notes/detectors.md`); the rest comes later.

Layer 2: may import only `muoscene` modules of layer <= 2 (see `muoscene._layers`).
"""

from .base import Detector, LinesOfSight, Pose, lines_of_sight
from .emulsion import EmulsionStack, TrackDensity
from .grid import (
    DEFAULT_BIN_WIDTH,
    DEFAULT_SUPERSAMPLE,
    DEFAULT_TAN_MAX,
    DirectionGrid,
    Subsamples,
    clipped_solid_angle_sr,
    rectangle_solid_angle_sr,
)

__all__ = [
    "DEFAULT_BIN_WIDTH",
    "DEFAULT_SUPERSAMPLE",
    "DEFAULT_TAN_MAX",
    "Detector",
    "DirectionGrid",
    "EmulsionStack",
    "LinesOfSight",
    "Pose",
    "Subsamples",
    "TrackDensity",
    "clipped_solid_angle_sr",
    "lines_of_sight",
    "rectangle_solid_angle_sr",
]
