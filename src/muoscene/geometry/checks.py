"""geometry.checks — what a scene must satisfy before it is built, and what a scan must not cross.

Analytic checks on the description, made when a scene is built:

- names are unique, and "air" is the air of the scene;
- shapes of structures, cuts and anomalies may touch but not overlap; in version 0.1 a cut that
  crosses a structure is refused, not resolved;
- a structure lies inside one ground or in the air (`host_ground`); the builder then splits it
  at the tops of the layers;
- an anomaly lies inside exactly one piece of one body (`builder`).

Checks of a scan (`Scene.scan`): the detector is not on a boundary between bodies
(`BOUNDARY_M`), and no line of sight is cut by the limits of a ground (`cut_by_limits`).

`check_overlaps` is the check of Calzone itself (Geant4 samples the surfaces of all volumes):
slower, independent of the above, and available as `Scene.check`.

Lengths in metres.
"""

from __future__ import annotations

import math

import calzone
import numpy as np

from .primitives import Anomaly, Box, Cut, Cylinder, Ground, Structure

TOUCH_M = 1e-9
"""Shapes closer than this touch; a deeper intersection is an overlap. Round-off of coordinates
of hundreds of metres is 1e-13 m."""
BOUNDARY_M = 1e-6
"""A scan refuses a point closer than this to a boundary between bodies: there the medium is
ambiguous. The frozen cases of the spike keep the detector 1 mm above the niche floor."""
AIR = "air"
"""Name of the air of the scene in a `PathTable`."""


def check_names(bodies) -> None:
    """Refuse bodies of an unknown kind, repeated names and the reserved name of the air."""
    seen: set[str] = set()
    for body in bodies:
        if not isinstance(body, Ground | Structure | Cut | Anomaly):
            raise TypeError(f"a scene is made of Ground, Structure, Cut and Anomaly bodies, not {body!r}")
        names = [body.name, *(layer.name for layer in getattr(body, "layers", ()))]
        for name in names:
            if name == AIR:
                raise ValueError(f'the name "{AIR}" is that of the air of the scene')
            if name in seen:
                raise ValueError(f"two bodies are named {name!r}")
            seen.add(name)


def overlap(a: Box | Cylinder, b: Box | Cylinder) -> bool:
    """Whether two shapes share a volume: touching faces do not count."""
    (a_lo, a_hi), (b_lo, b_hi) = a.aabb(), b.aabb()
    if min(a_hi[2], b_hi[2]) - max(a_lo[2], b_lo[2]) <= TOUCH_M:
        return False
    if isinstance(a, Box) and isinstance(b, Box):
        return all(min(a_hi[i], b_hi[i]) - max(a_lo[i], b_lo[i]) > TOUCH_M for i in (0, 1))
    if isinstance(a, Cylinder) and isinstance(b, Cylinder):
        (ax, ay), (bx, by) = a.centre_xy_m, b.centre_xy_m
        return a.radius_m + b.radius_m - math.hypot(ax - bx, ay - by) > TOUCH_M
    box, cylinder = (a, b) if isinstance(a, Box) else (b, a)
    (lo, hi), (cx, cy) = box.aabb(), cylinder.centre_xy_m
    dx = max(lo[0] - cx, 0.0, cx - hi[0])
    dy = max(lo[1] - cy, 0.0, cy - hi[1])
    return cylinder.radius_m - math.hypot(dx, dy) > TOUCH_M


def check_overlaps_of_shapes(structures, cuts) -> None:
    """Refuse structures that overlap each other, cuts that overlap, and a cut that crosses a structure."""
    solid = [(s.name, shape) for s in structures for shape in s.shapes]
    for i, (name, shape) in enumerate(solid):
        for other_name, other in solid[i + 1 :]:
            if overlap(shape, other):
                both = f"{name!r}" if name == other_name else f"{name!r} and {other_name!r}"
                raise ValueError(f"shapes of {both} overlap: structures may touch but not overlap")
    for i, cut in enumerate(cuts):
        for other in cuts[i + 1 :]:
            if overlap(cut.shape, other.shape):
                raise ValueError(
                    f"cuts {cut.name!r} and {other.name!r} overlap: give them as one cut or apart"
                )
        for name, shape in solid:
            if overlap(cut.shape, shape):
                raise ValueError(
                    f"cut {cut.name!r} crosses structure {name!r}: version 0.1 does not cut structures;"
                    " shape the structure around the cut"
                )


def check_grounds(grounds) -> None:
    """Refuse grounds whose footprints overlap: two grounds stand side by side."""
    for i, a in enumerate(grounds):
        for b in grounds[i + 1 :]:
            along_x = min(a.x_max_m, b.x_max_m) - max(a.x_min_m, b.x_min_m)
            along_y = min(a.y_max_m, b.y_max_m) - max(a.y_min_m, b.y_min_m)
            if along_x > TOUCH_M and along_y > TOUCH_M:
                raise ValueError(f"grounds {a.name!r} and {b.name!r} overlap in plan")


def host_ground(shape: Box | Cylinder, grounds, name: str) -> Ground | None:
    """The ground a shape of a structure lies in; None if it stands in the air.

    The ground is the same along y, so a shape is inside it when its outline in the x-z plane
    is under the profile and its extent along y within the limits. A shape that crosses the
    surface or a limit of a ground is refused.
    """
    lo, hi = shape.aabb()
    host = None
    for ground in grounds:
        if hi[1] <= ground.y_min_m + TOUCH_M or lo[1] >= ground.y_max_m - TOUCH_M:
            continue
        if hi[2] <= ground.z_bottom_m + TOUCH_M:
            continue
        low, high = ground.surface_range_m(lo[0], hi[0])
        if math.isnan(low) or lo[2] >= high - TOUCH_M:
            continue  # beside the ground or above its surface
        inside = (
            lo[0] >= ground.x_min_m - TOUCH_M
            and hi[0] <= ground.x_max_m + TOUCH_M
            and lo[1] >= ground.y_min_m - TOUCH_M
            and hi[1] <= ground.y_max_m + TOUCH_M
            and lo[2] >= ground.z_bottom_m - TOUCH_M
            and hi[2] <= low + TOUCH_M
        )
        if not inside:
            raise ValueError(
                f"a shape of structure {name!r} crosses the surface or a limit of ground {ground.name!r}"
                f" (it spans z from {lo[2]:g} to {hi[2]:g} m where the surface lies between {low:g} and"
                f" {high:g} m): give the part in the ground and the part in the air as two structures"
            )
        host = ground
    return host


def cut_by_limits(ground: Ground, position_m: np.ndarray, direction: np.ndarray) -> np.ndarray:
    """Which lines of sight cross a limit of the model inside the ground.

    position_m, direction   arrays of shape (n, 3): points in metres and unit vectors towards the sky
    Returns a boolean array of shape (n,). The limits are the bottom of the ground, its two faces
    across y and the vertical faces under the first and the last point of its profile; beyond
    them the model has nothing, so a path through them is not the path in the real ground.
    """
    o, u = position_m, direction
    out = np.zeros(len(o), dtype=bool)
    for axis, level in ((1, ground.y_min_m), (1, ground.y_max_m), (0, ground.x_min_m), (0, ground.x_max_m)):
        with np.errstate(divide="ignore", invalid="ignore"):
            t = (level - o[:, axis]) / u[:, axis]
        ahead = np.isfinite(t) & (t > 0)
        t = np.where(ahead, t, 0.0)
        x, y, z = (o[:, i] + t * u[:, i] for i in range(3))
        if axis == 1:
            top = ground.surface_m(x)
            hit = (x > ground.x_min_m + TOUCH_M) & (x < ground.x_max_m - TOUCH_M)
        else:
            end = ground.profile_xz_m[0][1] if level == ground.x_min_m else ground.profile_xz_m[-1][1]
            top = np.full(len(o), end)
            hit = (y > ground.y_min_m + TOUCH_M) & (y < ground.y_max_m - TOUCH_M)
        with np.errstate(invalid="ignore"):
            hit &= (z > ground.z_bottom_m + TOUCH_M) & (z < top - TOUCH_M)
        out |= ahead & hit
    # the bottom, for a point below it
    with np.errstate(divide="ignore", invalid="ignore"):
        t = (ground.z_bottom_m - o[:, 2]) / u[:, 2]
    ahead = np.isfinite(t) & (t > 0)
    t = np.where(ahead, t, 0.0)
    x, y = o[:, 0] + t * u[:, 0], o[:, 1] + t * u[:, 1]
    with np.errstate(invalid="ignore"):
        under = ground.surface_m(x) > ground.z_bottom_m + TOUCH_M
    inside = (x > ground.x_min_m) & (x < ground.x_max_m) & (y > ground.y_min_m) & (y < ground.y_max_m)
    return out | (ahead & inside & under)


def check_overlaps(geometry: calzone.Geometry, resolution: int | None = None) -> None:
    """The overlap check of Calzone: ValueError with the report of Geant4 if two volumes overlap.

    `resolution` is the number of points Geant4 samples on the surface of each volume (the
    default of Calzone when None). Volumes that only touch pass.
    """
    try:
        if resolution is None:
            geometry.check()
        else:
            geometry.check(resolution)
    except calzone.Geant4Exception as error:
        raise ValueError(f"the built geometry has overlapping volumes: {error}") from error
