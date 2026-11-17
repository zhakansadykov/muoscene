"""geometry.primitives — body descriptions: extruded cross-section, box, cylinder, pile grid, sphere, niche.

A scene is a list of bodies in the world frame (local ENU, metres; `docs/conventions.md`):

- `Ground` — the ground as an extruded cross-section: a profile of the surface in the x-z plane,
  a flat bottom and two limits along y. Horizontal `Layer`s divide it by elevation.
- `Structure` — a box, a vertical cylinder or several of them (a pile cap, a pile grid).
- `Cut` — a niche or an excavation: a box of air that takes precedence over everything else.
- `Anomaly` — a sphere that scales the density of the body it lies in by (1 + `drho_rel`).

A body has a name, a shape, a material and a bulk density. The material is a `Material`, a
`Mixture` or a `Soil` of `muoscene.materials`; a soil brings its own bulk density ρ_d (1 + w), for
the others the density of the material is the default. Nothing here knows Calzone or Mulder:
shapes only answer geometric questions (extent, signed distance of a point), which the builder
and the checks use.

Version 0.1 has no rotations: boxes are aligned with the axes, cylinders are vertical and the
ground is extruded along y. Meshes, terrain models and non-planar layers come with version 0.2.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from ..materials.composition import Material, Mixture, Soil

AnyMaterial = Material | Mixture


def _vector(values, n: int, what: str) -> tuple[float, ...]:
    out = tuple(float(v) for v in np.asarray(values, dtype=float).ravel())
    if len(out) != n or not all(math.isfinite(v) for v in out):
        raise ValueError(f"{what}: {n} finite numbers are needed")
    return out


def _positive(value, what: str) -> float:
    value = float(value)
    if not (math.isfinite(value) and value > 0):
        raise ValueError(f"{what} must be positive")
    return value


def _name(name, what: str) -> str:
    if not isinstance(name, str) or not name.strip():
        raise ValueError(f"{what} needs a name")
    return name


def resolve_material(material, density_g_cm3, what: str) -> tuple[AnyMaterial, float]:
    """Material and bulk density (g/cm³) of a body.

    A `Soil` gives its mixture of skeleton and water and its bulk density ρ_d (1 + w), and takes
    no density of its own; a `Material` or a `Mixture` defaults to its own density.
    """
    if isinstance(material, Soil):
        if density_g_cm3 is not None:
            raise ValueError(
                f"{what}: a soil brings its bulk density (dry density and water content);"
                " change those instead of giving density_g_cm3"
            )
        return material.material, float(material.bulk_density_g_cm3)
    if not isinstance(material, Material | Mixture):
        raise TypeError(f"{what}: the material must be a Material, a Mixture or a Soil of muoscene.materials")
    if density_g_cm3 is None:
        return material, float(material.density_g_cm3)
    return material, _positive(density_g_cm3, f"{what}: density_g_cm3")


# ------------------------------------------------------------------ shapes


@dataclass(frozen=True)
class Box:
    """A box aligned with the axes, given by two opposite corners in metres."""

    lo_m: tuple[float, float, float]
    hi_m: tuple[float, float, float]

    def __post_init__(self) -> None:
        object.__setattr__(self, "lo_m", _vector(self.lo_m, 3, "lo_m"))
        object.__setattr__(self, "hi_m", _vector(self.hi_m, 3, "hi_m"))
        if not all(h > low for low, h in zip(self.lo_m, self.hi_m, strict=True)):
            raise ValueError("a box needs hi_m above lo_m along every axis")

    @classmethod
    def centred(cls, centre_m, size_m) -> Box:
        """From the centre and the full sizes along x, y and z, in metres."""
        c, s = np.array(_vector(centre_m, 3, "centre_m")), np.array(_vector(size_m, 3, "size_m"))
        return cls(tuple(c - s / 2), tuple(c + s / 2))

    @property
    def centre_m(self) -> tuple[float, float, float]:
        """Centre, m."""
        return tuple((low + h) / 2 for low, h in zip(self.lo_m, self.hi_m, strict=True))

    @property
    def size_m(self) -> tuple[float, float, float]:
        """Full sizes along x, y and z, m."""
        return tuple(h - low for low, h in zip(self.lo_m, self.hi_m, strict=True))

    def aabb(self) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
        """Lowest and highest corner, m."""
        return self.lo_m, self.hi_m

    def clip_z(self, z_lo_m: float, z_hi_m: float) -> Box:
        """The part of the box between two elevations, m."""
        lo, hi = max(self.lo_m[2], z_lo_m), min(self.hi_m[2], z_hi_m)
        return Box((*self.lo_m[:2], lo), (*self.hi_m[:2], hi))

    def signed_distance_m(self, point_m):
        """Distance of points (shape (..., 3), m) to the surface, m: negative inside, positive outside."""
        p = np.asarray(point_m, dtype=float)
        q = np.abs(p - np.array(self.centre_m)) - np.array(self.size_m) / 2
        return _signed(np.linalg.norm(np.maximum(q, 0.0), axis=-1), q.max(axis=-1))


@dataclass(frozen=True)
class Cylinder:
    """A vertical cylinder: axis at (x, y), from `z_bottom_m` to `z_top_m`, radius in metres."""

    centre_xy_m: tuple[float, float]
    z_bottom_m: float
    z_top_m: float
    radius_m: float

    def __post_init__(self) -> None:
        object.__setattr__(self, "centre_xy_m", _vector(self.centre_xy_m, 2, "centre_xy_m"))
        for field in ("z_bottom_m", "z_top_m"):
            object.__setattr__(self, field, _vector([getattr(self, field)], 1, field)[0])
        object.__setattr__(self, "radius_m", _positive(self.radius_m, "radius_m"))
        if not self.z_top_m > self.z_bottom_m:
            raise ValueError("a cylinder needs z_top_m above z_bottom_m")

    @property
    def centre_m(self) -> tuple[float, float, float]:
        """Centre, m."""
        return (*self.centre_xy_m, (self.z_bottom_m + self.z_top_m) / 2)

    @property
    def length_m(self) -> float:
        """Length along the axis, m."""
        return self.z_top_m - self.z_bottom_m

    def aabb(self) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
        """Lowest and highest corner of the bounding box, m."""
        (x, y), r = self.centre_xy_m, self.radius_m
        return (x - r, y - r, self.z_bottom_m), (x + r, y + r, self.z_top_m)

    def clip_z(self, z_lo_m: float, z_hi_m: float) -> Cylinder:
        """The part of the cylinder between two elevations, m."""
        return Cylinder(
            self.centre_xy_m, max(self.z_bottom_m, z_lo_m), min(self.z_top_m, z_hi_m), self.radius_m
        )

    def signed_distance_m(self, point_m):
        """Distance of points (shape (..., 3), m) to the surface, m: negative inside, positive outside."""
        p = np.asarray(point_m, dtype=float)
        radial = np.hypot(p[..., 0] - self.centre_xy_m[0], p[..., 1] - self.centre_xy_m[1]) - self.radius_m
        axial = np.abs(p[..., 2] - self.centre_m[2]) - self.length_m / 2
        return _signed(np.hypot(np.maximum(radial, 0.0), np.maximum(axial, 0.0)), np.maximum(radial, axial))


@dataclass(frozen=True)
class Sphere:
    """A sphere: centre and radius in metres."""

    centre_m: tuple[float, float, float]
    radius_m: float

    def __post_init__(self) -> None:
        object.__setattr__(self, "centre_m", _vector(self.centre_m, 3, "centre_m"))
        object.__setattr__(self, "radius_m", _positive(self.radius_m, "radius_m"))

    def aabb(self) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
        """Lowest and highest corner of the bounding box, m."""
        c, r = self.centre_m, self.radius_m
        return (c[0] - r, c[1] - r, c[2] - r), (c[0] + r, c[1] + r, c[2] + r)

    def signed_distance_m(self, point_m):
        """Distance of points (shape (..., 3), m) to the surface, m: negative inside, positive outside."""
        distance = np.linalg.norm(np.asarray(point_m, dtype=float) - np.array(self.centre_m), axis=-1)
        return _signed(distance - self.radius_m, distance - self.radius_m)


@dataclass(frozen=True)
class Prism:
    """A polygon in the x-z plane extruded along y: one piece of a ground or of a layer.

    `polygon_xz_m` is a simple polygon, counter-clockwise for x to the right and z up. Made by
    `Ground.pieces`; not a body of its own.
    """

    polygon_xz_m: tuple[tuple[float, float], ...]
    y_min_m: float
    y_max_m: float

    def aabb(self) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
        """Lowest and highest corner of the bounding box, m."""
        x, z = zip(*self.polygon_xz_m, strict=True)
        return (min(x), self.y_min_m, min(z)), (max(x), self.y_max_m, max(z))

    def signed_distance_m(self, point_m):
        """Distance of points (shape (..., 3), m) to the surface, m: negative inside, positive outside."""
        p = np.asarray(point_m, dtype=float)
        in_plane = polygon_signed_distance(self.polygon_xz_m, p[..., 0], p[..., 2])
        along = np.abs(p[..., 1] - (self.y_min_m + self.y_max_m) / 2) - (self.y_max_m - self.y_min_m) / 2
        return _signed(
            np.hypot(np.maximum(in_plane, 0.0), np.maximum(along, 0.0)), np.maximum(in_plane, along)
        )


def _signed(outside, inside):
    """The distance outside where it is positive, the negative one inside elsewhere; a float for one point."""
    value = np.where(outside > 0, outside, inside)
    return float(value) if value.ndim == 0 else value


def polygon_signed_distance(polygon, x_m, z_m) -> np.ndarray:
    """Distance of the points (x, z) to the outline of a simple polygon, m: negative inside."""
    a = np.asarray(polygon, dtype=float)
    b = np.roll(a, -1, axis=0)
    edge = b - a
    point = np.stack(np.broadcast_arrays(np.asarray(x_m, dtype=float), np.asarray(z_m, dtype=float)), axis=-1)
    to_point = point[..., None, :] - a
    along = np.clip((to_point * edge).sum(axis=-1) / (edge * edge).sum(axis=-1), 0.0, 1.0)
    distance = np.linalg.norm(to_point - along[..., None] * edge, axis=-1).min(axis=-1)
    # even-odd rule along a ray towards +x
    x, z = point[..., None, 0], point[..., None, 1]
    crosses = (a[:, 1] > z) != (b[:, 1] > z)
    with np.errstate(divide="ignore", invalid="ignore"):
        x_at = a[:, 0] + (z - a[:, 1]) * edge[:, 0] / edge[:, 1]
    inside = np.count_nonzero(crosses & (x_at > x), axis=-1) % 2 == 1
    return np.where(inside, -distance, distance)


def pile_grid(
    centre_xy_m,
    z_top_m: float,
    length_m: float,
    radius_m: float,
    n_x: int,
    n_y: int,
    pitch_x_m: float,
    pitch_y_m: float | None = None,
) -> tuple[Cylinder, ...]:
    """Vertical piles on a rectangular grid: `n_x` by `n_y` cylinders around `centre_xy_m`.

    The piles hang from `z_top_m` (the underside of a pile cap) down by `length_m`; the pitch is
    the distance between axes, in metres, the same along y unless `pitch_y_m` is given. The
    order is x first: pile i * n_y + j stands at column i along x and row j along y.
    """
    cx, cy = _vector(centre_xy_m, 2, "centre_xy_m")
    if not (isinstance(n_x, int) and isinstance(n_y, int) and n_x >= 1 and n_y >= 1):
        raise ValueError("n_x and n_y must be positive integers")
    length, radius = _positive(length_m, "length_m"), _positive(radius_m, "radius_m")
    px = _positive(pitch_x_m, "pitch_x_m")
    py = px if pitch_y_m is None else _positive(pitch_y_m, "pitch_y_m")
    if (n_x > 1 and px < 2 * radius) or (n_y > 1 and py < 2 * radius):
        raise ValueError("the pitch of a pile grid is below the diameter of a pile")
    z_top = _vector([z_top_m], 1, "z_top_m")[0]
    return tuple(
        Cylinder(
            (cx + (i - (n_x - 1) / 2) * px, cy + (j - (n_y - 1) / 2) * py), z_top - length, z_top, radius
        )
        for i in range(n_x)
        for j in range(n_y)
    )


# ------------------------------------------------------------------ bodies


@dataclass(frozen=True)
class Layer:
    """A horizontal layer of a ground: everything below `z_top_m` down to the next layer.

    name            name of the body in a `PathTable`
    z_top_m         elevation of its top in the world frame, m
    material        `Material`, `Mixture` or `Soil`
    density_g_cm3   bulk density; by default that of the soil or of the material
    """

    name: str
    z_top_m: float
    material: AnyMaterial | Soil
    density_g_cm3: float | None = None

    def __post_init__(self) -> None:
        _name(self.name, "a layer")
        object.__setattr__(self, "z_top_m", _vector([self.z_top_m], 1, f"{self.name}: z_top_m")[0])
        resolve_material(self.material, self.density_g_cm3, self.name)

    @property
    def composition(self) -> AnyMaterial:
        """Material that Mulder is given."""
        return resolve_material(self.material, self.density_g_cm3, self.name)[0]

    @property
    def bulk_density_g_cm3(self) -> float:
        """Bulk density, g/cm³."""
        return resolve_material(self.material, self.density_g_cm3, self.name)[1]


@dataclass(frozen=True)
class GroundPiece:
    """One connected piece of a ground or of a layer: the body it belongs to and its shape."""

    owner: str
    z_lo_m: float
    z_hi_m: float
    shape: Prism


@dataclass(frozen=True)
class Ground:
    """The ground as an extruded cross-section.

    name            name of the body in a `PathTable`
    profile_xz_m    the surface in the x-z plane, m: points (x, z) with x not decreasing. Two
                    points with the same x make a vertical wall; overhangs cannot be given (a
                    niche is a `Cut`). The ground is everything below the profile.
    z_bottom_m      flat bottom of the model, m
    y_min_m, y_max_m  limits of the extrusion along y, m
    material        `Material`, `Mixture` or `Soil` of the part above the first layer
    density_g_cm3   bulk density; by default that of the soil or of the material
    layers          `Layer`s from the top down

    The bottom, the two faces across y and the vertical faces under the first and the last point
    of the profile are limits of the model, not surfaces of the object: a line of sight that
    crosses one of them inside the ground is refused by `Scene.scan`. A canyon wall is therefore
    a vertical segment of the profile, and a profile may start or end at `z_bottom_m`.
    """

    name: str
    profile_xz_m: tuple[tuple[float, float], ...]
    z_bottom_m: float
    y_min_m: float
    y_max_m: float
    material: AnyMaterial | Soil
    density_g_cm3: float | None = None
    layers: tuple[Layer, ...] = ()

    def __post_init__(self) -> None:
        _name(self.name, "a ground")
        profile = tuple(
            _vector(point, 2, f"{self.name}: a point of the profile") for point in self.profile_xz_m
        )
        object.__setattr__(self, "profile_xz_m", profile)
        for field in ("z_bottom_m", "y_min_m", "y_max_m"):
            object.__setattr__(self, field, _vector([getattr(self, field)], 1, f"{self.name}: {field}")[0])
        object.__setattr__(self, "layers", tuple(self.layers))
        resolve_material(self.material, self.density_g_cm3, self.name)
        if not self.y_max_m > self.y_min_m:
            raise ValueError(f"{self.name}: y_max_m must be above y_min_m")
        x, z = [p[0] for p in profile], [p[1] for p in profile]
        if len(profile) < 2 or not x[-1] > x[0]:
            raise ValueError(f"{self.name}: the profile needs at least two points and a width along x")
        if any(b < a for a, b in zip(x[:-1], x[1:], strict=True)):
            raise ValueError(f"{self.name}: x must not decrease along the profile (no overhangs)")
        if any(p == q for p, q in zip(profile[:-1], profile[1:], strict=True)):
            raise ValueError(f"{self.name}: the profile repeats a point")
        if any(a == b == c for a, b, c in zip(x[:-2], x[1:-1], x[2:], strict=True)):
            raise ValueError(f"{self.name}: three points of the profile share one x")
        if min(z) < self.z_bottom_m or not max(z) > self.z_bottom_m:
            raise ValueError(f"{self.name}: the profile must lie above z_bottom_m")
        if not all(isinstance(layer, Layer) for layer in self.layers):
            raise TypeError(f"{self.name}: layers must be Layer objects")
        tops = [layer.z_top_m for layer in self.layers]
        if any(b >= a for a, b in zip(tops[:-1], tops[1:], strict=True)):
            raise ValueError(f"{self.name}: layers go from the top down, each top below the one before")
        if tops and not (tops[0] < max(z) and tops[-1] > self.z_bottom_m):
            raise ValueError(
                f"{self.name}: the tops of the layers must lie below the highest point of the profile"
                " and above z_bottom_m (the material of the ground is that of its top part)"
            )

    @property
    def composition(self) -> AnyMaterial:
        """Material that Mulder is given for the part above the first layer."""
        return resolve_material(self.material, self.density_g_cm3, self.name)[0]

    @property
    def bulk_density_g_cm3(self) -> float:
        """Bulk density of the part above the first layer, g/cm³."""
        return resolve_material(self.material, self.density_g_cm3, self.name)[1]

    @property
    def x_min_m(self) -> float:
        """x of the first point of the profile, m."""
        return self.profile_xz_m[0][0]

    @property
    def x_max_m(self) -> float:
        """x of the last point of the profile, m."""
        return self.profile_xz_m[-1][0]

    @property
    def z_max_m(self) -> float:
        """Highest point of the profile, m."""
        return max(z for _, z in self.profile_xz_m)

    def surface_m(self, x_m, upper: bool = True) -> np.ndarray:
        """Elevation of the surface (m) above the points `x_m`; nan outside the profile.

        At a vertical wall the surface has two values: the upper one, or the lower one with
        `upper=False`.
        """
        x = np.asarray(x_m, dtype=float)
        px = np.array([p[0] for p in self.profile_xz_m])
        pz = np.array([p[1] for p in self.profile_xz_m])
        out = np.full(x.shape, -np.inf if upper else np.inf)
        pick = np.maximum if upper else np.minimum
        for i in range(len(px) - 1):
            if px[i + 1] > px[i]:
                on = (x >= px[i]) & (x <= px[i + 1])
                height = pz[i] + (pz[i + 1] - pz[i]) * (x - px[i]) / (px[i + 1] - px[i])
            else:  # a vertical wall: its lower or its upper end
                on = x == px[i]
                height = pick(pz[i], pz[i + 1])
            out = np.where(on, pick(out, height), out)
        return np.where(np.isfinite(out), out, np.nan)

    def surface_range_m(self, x_lo_m: float, x_hi_m: float) -> tuple[float, float]:
        """Lowest and highest elevation of the surface (m) between two x, the ends excluded.

        A vertical wall exactly at an end does not count: a body that stands flush against a
        wall sees only the surface on its own side. nan when the interval misses the profile.
        """
        px = [p[0] for p in self.profile_xz_m]
        pz = [p[1] for p in self.profile_xz_m]
        lo, hi = max(x_lo_m, px[0]), min(x_hi_m, px[-1])
        if not hi > lo:
            return math.nan, math.nan
        heights = [z for x, z in self.profile_xz_m if lo < x < hi]
        for i in range(len(px) - 1):
            if px[i + 1] > px[i]:
                slope = (pz[i + 1] - pz[i]) / (px[i + 1] - px[i])
                if px[i] <= lo < px[i + 1]:
                    heights.append(pz[i] + slope * (lo - px[i]))
                if px[i] < hi <= px[i + 1]:
                    heights.append(pz[i] + slope * (hi - px[i]))
        return min(heights), max(heights)

    def pieces(self) -> tuple[GroundPiece, ...]:
        """The ground cut at the tops of its layers: connected pieces, each with its owner.

        A band between two elevations is one polygon for every stretch of the profile that rises
        above its lower edge, so the two banks of a canyon give two pieces of the upper layers.
        Neighbouring pieces share their vertices exactly.
        """
        owners = [self.name, *(layer.name for layer in self.layers)]
        tops = [None, *(layer.z_top_m for layer in self.layers)]
        bottoms = [*(layer.z_top_m for layer in self.layers), self.z_bottom_m]
        out = []
        for owner, z_hi, z_lo in zip(owners, tops, bottoms, strict=True):
            for polygon in _bands(self.profile_xz_m, z_lo, z_hi):
                top = max(z for _, z in polygon)
                out.append(GroundPiece(owner, z_lo, top, Prism(polygon, self.y_min_m, self.y_max_m)))
        return tuple(out)


def _bands(profile, z_lo: float, z_hi: float | None) -> list[tuple[tuple[float, float], ...]]:
    """Polygons of the region under the profile between `z_lo` and `z_hi` (None: no upper limit)."""
    fine = [profile[0]]
    for (x0, z0), (x1, z1) in zip(profile[:-1], profile[1:], strict=True):
        crossings = []
        for level in (z_lo, z_hi):
            if level is not None and (z0 - level) * (z1 - level) < 0:
                t = (level - z0) / (z1 - z0)
                crossings.append((t, (x0 + t * (x1 - x0), level)))
        fine += [point for _, point in sorted(crossings)] + [(x1, z1)]
    runs, run = [], []
    for (x0, z0), (x1, z1) in zip(fine[:-1], fine[1:], strict=True):
        if max(z0, z1) > z_lo:
            run = run or [(x0, z0)]
            run.append((x1, z1))
            if z1 <= z_lo:  # the profile comes down to the lower edge: the piece ends here
                runs.append(run)
                run = []
        elif run:
            runs.append(run)
            run = []
    if run:
        runs.append(run)
    out = []
    for run in runs:
        top = [(x, z if z_hi is None else min(z, z_hi)) for x, z in run]
        polygon = _without_spare_vertices([(top[0][0], z_lo), (top[-1][0], z_lo), *top[::-1]])
        if len(polygon) >= 3 and polygon_area_m2(polygon) > SLIVER_M2:
            out.append(tuple(polygon))
    return out


SLIVER_M2 = 1e-9
"""Pieces of a cross-section below this area, m², are dropped: a mesh cannot be made of them."""


def polygon_area_m2(polygon) -> float:
    """Area of a polygon in the x-z plane, m²: positive when it is counter-clockwise."""
    a = np.asarray(polygon, dtype=float)
    b = np.roll(a, -1, axis=0)
    return float(np.sum(a[:, 0] * b[:, 1] - b[:, 0] * a[:, 1])) / 2


def _without_spare_vertices(polygon: list[tuple[float, float]]) -> list[tuple[float, float]]:
    """A polygon without repeated vertices and without vertices on a straight line."""
    points = [p for i, p in enumerate(polygon) if p != polygon[i - 1]]
    changed = True
    while changed and len(points) >= 3:
        changed = False
        for i in range(len(points)):
            (ax, az), (bx, bz), (cx, cz) = points[i - 1], points[i], points[(i + 1) % len(points)]
            if (bx - ax) * (cz - az) - (bz - az) * (cx - ax) == 0:
                points.pop(i)
                changed = True
                break
    return points


@dataclass(frozen=True)
class Structure:
    """A structure of one material: a box, a vertical cylinder or several of them.

    name            name of the body in a `PathTable`; all its shapes share one column
    shapes          a `Box`, a `Cylinder` or a sequence of them (e.g. `pile_grid(...)`)
    material        `Material`, `Mixture` or `Soil`
    density_g_cm3   bulk density; by default that of the soil or of the material

    A shape lies either inside one ground or in the air; where it crosses the top of a layer the
    builder splits it. Shapes of structures may touch but not overlap.
    """

    name: str
    shapes: Box | Cylinder | tuple[Box | Cylinder, ...]
    material: AnyMaterial | Soil
    density_g_cm3: float | None = None

    def __post_init__(self) -> None:
        _name(self.name, "a structure")
        shapes = (self.shapes,) if isinstance(self.shapes, Box | Cylinder) else tuple(self.shapes)
        if not shapes or not all(isinstance(shape, Box | Cylinder) for shape in shapes):
            raise TypeError(f"{self.name}: a structure is made of Box and Cylinder shapes")
        object.__setattr__(self, "shapes", shapes)
        resolve_material(self.material, self.density_g_cm3, self.name)

    @property
    def composition(self) -> AnyMaterial:
        """Material that Mulder is given."""
        return resolve_material(self.material, self.density_g_cm3, self.name)[0]

    @property
    def bulk_density_g_cm3(self) -> float:
        """Bulk density, g/cm³."""
        return resolve_material(self.material, self.density_g_cm3, self.name)[1]


@dataclass(frozen=True)
class Cut:
    """A niche or an excavation: a box of air. Cuts take precedence over every other body.

    name    name of the body in a `PathTable`: the path through the air of the cut
    shape   a `Box`; it may reach out of the ground into the open air
    """

    name: str
    shape: Box

    def __post_init__(self) -> None:
        _name(self.name, "a cut")
        if not isinstance(self.shape, Box):
            raise TypeError(f"{self.name}: a cut is a Box")


@dataclass(frozen=True)
class Anomaly:
    """A sphere that scales the density of its host by (1 + `drho_rel`).

    name       name of the body in a `PathTable`
    shape      a `Sphere`; it must lie inside exactly one piece of one body, its host, which the
               builder finds: a layer, a ground or a structure
    drho_rel   relative change of density: -1 is a void, -0.2 a 20 % decompaction
    """

    name: str
    shape: Sphere
    drho_rel: float

    def __post_init__(self) -> None:
        _name(self.name, "an anomaly")
        if not isinstance(self.shape, Sphere):
            raise TypeError(f"{self.name}: an anomaly is a Sphere")
        value = float(self.drho_rel)
        if not (math.isfinite(value) and value >= -1):
            raise ValueError(f"{self.name}: drho_rel must be at least -1 (a void)")
        object.__setattr__(self, "drho_rel", value)


Body = Ground | Structure | Cut | Anomaly
