"""geometry.builder — bodies to a Calzone geometry tree: units, placement, subtraction, splitting.

The precedence of the conventions — cuts, then structures, then layers, then the ground — is
realised with what Calzone offers:

- the pieces of every ground and of its layers are closed meshes, one for each connected piece
  between two layer tops (`Ground.pieces`), placed side by side in the world; they touch along
  the layer tops and do not overlap;
- a cut is a box of air in the world that every piece of ground it reaches subtracts (Calzone
  subtracts only volumes that have the same mother, which is why layers are not nested);
- a shape of a structure is split at the tops of the layers it crosses and each part becomes a
  daughter of the piece of ground around it; a structure in the air is a daughter of the world;
- an anomaly is a daughter of the one volume it lies in, its host.

Units and frames at this boundary: muoscene metres become Calzone centimetres; a daughter is
placed relative to the origin of its mother — the centre of a box or of a cylinder, the origin
of the scene for a mesh and for the world, whose own position only moves its box.

Calzone takes a mesh only from a file: the meshes are written as OBJ files with face normals
into a temporary directory of the system, which is removed as soon as the geometry is built.
The materials of the Calzone volumes are placeholders (air and water of Geant4); the materials
and densities of the scene are set on the media of Mulder (`scene`).
"""

from __future__ import annotations

import pathlib
import tempfile
from dataclasses import dataclass

import calzone
import numpy as np

from . import checks
from .primitives import (
    Anomaly,
    Box,
    Cut,
    Cylinder,
    Ground,
    Prism,
    Sphere,
    Structure,
    polygon_area_m2,
)

CM_PER_M = 100.0
"""Calzone lengths are centimetres."""
MESH_ALGORITHM = "voxels"
"""Mesh traversal of Calzone. With "bvh" lines of sight that cross the face shared by two touching
meshes lose path: 2 % of them, by up to 6 m, in a ground with one layer top (`docs/notes/geometry.md`)."""
WORLD = "Scene"
"""Calzone name of the world volume, the air of the scene."""
PLACEHOLDER_AIR, PLACEHOLDER_SOLID = "G4_AIR", "G4_WATER"
MIN_PART_M = 1e-9
"""A part of a split shape thinner than this is dropped."""


@dataclass(frozen=True)
class Volume:
    """One Calzone volume and the body it belongs to.

    name     Calzone name, unique in the geometry (letters and digits)
    body     name of the body in a `PathTable`; "air" for the world
    kind     "air", "ground", "layer", "structure", "cut" or "anomaly"
    shape    its shape in the world frame, metres
    mother   Calzone name of the volume it is a daughter of
    """

    name: str
    body: str
    kind: str
    shape: Box | Cylinder | Sphere | Prism
    mother: str | None


@dataclass(frozen=True)
class Assembly:
    """A built scene: the Calzone geometry, its volumes and the box of the world.

    geometry      `calzone.Geometry`
    volumes       every volume, the world first
    lo_m, hi_m    corners of the world box in the world frame, m
    hosts         name of the host body of every anomaly
    """

    geometry: calzone.Geometry
    volumes: tuple[Volume, ...]
    lo_m: tuple[float, float, float]
    hi_m: tuple[float, float, float]
    hosts: dict[str, str]


def assemble(bodies, margin_m: float) -> Assembly:
    """Check the bodies and build the Calzone geometry; the world box is `margin_m` larger than them."""
    checks.check_names(bodies)
    grounds = [b for b in bodies if isinstance(b, Ground)]
    structures = [b for b in bodies if isinstance(b, Structure)]
    cuts = [b for b in bodies if isinstance(b, Cut)]
    anomalies = [b for b in bodies if isinstance(b, Anomaly)]
    if not (grounds or structures):
        raise ValueError("a scene needs at least one ground or structure")
    checks.check_grounds(grounds)
    checks.check_overlaps_of_shapes(structures, cuts)

    volumes: list[Volume] = []

    def add(body: str, kind: str, shape, mother: str) -> Volume:
        volume = Volume(f"V{len(volumes)}", body, kind, shape, mother)
        volumes.append(volume)
        return volume

    pieces: dict[str, list[tuple[Volume, float, float]]] = {}
    for ground in grounds:
        for piece in ground.pieces():
            kind = "ground" if piece.owner == ground.name else "layer"
            volume = add(piece.owner, kind, piece.shape, WORLD)
            pieces.setdefault(ground.name, []).append((volume, piece.z_lo_m, piece.z_hi_m))
    for structure in structures:
        for shape in structure.shapes:
            ground = checks.host_ground(shape, grounds, structure.name)
            if ground is None:
                add(structure.name, "structure", shape, WORLD)
                continue
            for volume, z_lo, z_hi in pieces[ground.name]:
                part = _part_inside(shape, volume.shape, z_lo, z_hi)
                if part is not None:
                    add(structure.name, "structure", part, volume.name)
    for cut in cuts:
        add(cut.name, "cut", cut.shape, WORLD)
    hosts = {}
    for i, anomaly in enumerate(anomalies):
        for other in anomalies[:i]:
            gap = np.linalg.norm(np.subtract(anomaly.shape.centre_m, other.shape.centre_m))
            if gap < anomaly.shape.radius_m + other.shape.radius_m - checks.TOUCH_M:
                raise ValueError(f"anomalies {anomaly.name!r} and {other.name!r} overlap")
        host = _host(anomaly, volumes, cuts)
        hosts[anomaly.name] = host.body
        add(anomaly.name, "anomaly", anomaly.shape, host.name)

    margin = float(margin_m)
    if not (np.isfinite(margin) and margin > 0):
        raise ValueError("margin_m must be positive")
    corners = np.array([v.shape.aabb() for v in volumes if v.mother == WORLD])
    lo, hi = corners[:, 0].min(axis=0) - margin, corners[:, 1].max(axis=0) + margin
    world = Volume(WORLD, checks.AIR, "air", Box(tuple(lo), tuple(hi)), None)
    with tempfile.TemporaryDirectory(prefix="muoscene-") as directory:
        geometry = calzone.Geometry({WORLD: _tree(world, volumes, pathlib.Path(directory))})
    return Assembly(geometry, (world, *volumes), tuple(lo.tolist()), tuple(hi.tolist()), hosts)


def _part_inside(shape: Box | Cylinder, piece: Prism, z_lo: float, z_hi: float):
    """The part of a shape of a structure inside one piece of its ground; None if there is none."""
    lo, hi = shape.aabb()
    (piece_lo, piece_hi) = piece.aabb()
    if min(hi[2], z_hi) - max(lo[2], z_lo) <= MIN_PART_M:
        return None
    if lo[0] < piece_lo[0] - checks.TOUCH_M or hi[0] > piece_hi[0] + checks.TOUCH_M:
        return None  # another piece of the same band: the other bank of a canyon
    return shape.clip_z(z_lo, z_hi)


def _host(anomaly: Anomaly, volumes: list[Volume], cuts: list[Cut]) -> Volume:
    """The volume an anomaly lies in; ValueError if it is not inside exactly one."""
    centre, radius = anomaly.shape.centre_m, anomaly.shape.radius_m
    host, crossed = None, []
    for kinds in (("structure",), ("ground", "layer")):
        for volume in volumes:
            if volume.kind not in kinds:
                continue
            distance = volume.shape.signed_distance_m(centre)
            if distance <= -radius + checks.TOUCH_M:
                host = volume
                break
            if distance < radius - checks.TOUCH_M:
                crossed.append(volume.body)
        if host is not None or crossed:
            break
    if host is None:
        if not crossed:
            raise ValueError(f"anomaly {anomaly.name!r} lies in the air: it needs a host body")
        names = ", ".join(repr(name) for name in dict.fromkeys(crossed))
        raise ValueError(
            f"anomaly {anomaly.name!r} crosses a boundary of {names}: an anomaly lies inside one body,"
            " clear of layer tops, surfaces, structures and the limits of the model"
            " (a structure is split where it crosses the top of a layer)"
        )
    for cut in cuts:
        if cut.shape.signed_distance_m(centre) < radius - checks.TOUCH_M:
            raise ValueError(f"anomaly {anomaly.name!r} reaches into cut {cut.name!r}")
    return host


# ------------------------------------------------------------------ the Calzone tree


def _tree(world: Volume, volumes: list[Volume], directory: pathlib.Path) -> dict:
    """Definition of the world volume with every other volume inside, as Calzone reads it."""
    by_name = {volume.name: volume for volume in (world, *volumes)}
    cuts = [volume for volume in volumes if volume.kind == "cut"]
    centre = np.array(world.shape.centre_m)
    nodes = {
        WORLD: {
            "box": (np.array(world.shape.size_m) * CM_PER_M).tolist(),
            "position": (centre * CM_PER_M).tolist(),
            "material": PLACEHOLDER_AIR,
        }
    }
    for volume in volumes:
        mother = by_name[volume.mother]
        origin = (
            np.zeros(3)
            if isinstance(mother.shape, Prism) or mother is world
            else np.array(mother.shape.centre_m)
        )
        node: dict = {"material": PLACEHOLDER_AIR if volume.kind == "cut" else PLACEHOLDER_SOLID}
        shape = volume.shape
        if isinstance(shape, Prism):
            path = directory / f"{volume.name}.obj"
            write_prism(path, shape)
            node["mesh"] = {"path": str(path), "units": "m", "algorithm": MESH_ALGORITHM}
            reached = [cut.name for cut in cuts if _boxes_overlap(shape.aabb(), cut.shape.aabb())]
            if reached:
                node["subtract"] = reached
        else:
            if isinstance(shape, Box):
                node["box"] = (np.array(shape.size_m) * CM_PER_M).tolist()
            elif isinstance(shape, Cylinder):
                node["cylinder"] = {"length": shape.length_m * CM_PER_M, "radius": shape.radius_m * CM_PER_M}
            else:
                node["sphere"] = {"radius": shape.radius_m * CM_PER_M}
            node["position"] = ((np.array(shape.centre_m) - origin) * CM_PER_M).tolist()
        nodes[volume.name] = node
        nodes[volume.mother][volume.name] = node
    return nodes[WORLD]


def _boxes_overlap(a, b) -> bool:
    return all(min(a[1][i], b[1][i]) - max(a[0][i], b[0][i]) > checks.TOUCH_M for i in range(3))


def triangulate(polygon) -> list[tuple[int, int, int]]:
    """Triangles of a simple counter-clockwise polygon, as triples of vertex indices (ear clipping)."""
    p = np.asarray(polygon, dtype=float)

    def turn(a: int, b: int, c: int) -> float:
        return (p[b, 0] - p[a, 0]) * (p[c, 1] - p[a, 1]) - (p[b, 1] - p[a, 1]) * (p[c, 0] - p[a, 0])

    left = list(range(len(p)))
    triangles = []
    while len(left) > 3:
        for k in range(len(left)):
            a, b, c = left[k - 1], left[k], left[(k + 1) % len(left)]
            if turn(a, b, c) <= 0:
                continue
            others = (q for q in left if q not in (a, b, c))
            if not any(turn(a, b, q) >= 0 and turn(b, c, q) >= 0 and turn(c, a, q) >= 0 for q in others):
                triangles.append((a, b, c))
                left.pop(k)
                break
        else:
            raise ValueError("the cross-section cannot be triangulated: its outline is not a simple polygon")
    triangles.append((left[0], left[1], left[2]))
    return triangles


def prism_mesh(prism: Prism) -> tuple[np.ndarray, list[tuple[int, int, int]], np.ndarray]:
    """The closed mesh of a prism, as `write_prism` writes it and Calzone reads it.

    Returns the vertices (m, shape (2n, 3): the polygon at y_min, then at y_max), the triangles as
    triples of vertex indices, counter-clockwise seen from outside, and their outward unit normals
    (shape (number of triangles, 3)). Every side face is two triangles that share a diagonal; the
    caps are triangulated by `triangulate`. `muoscene.geometry.edges` takes the edges of the
    meshes of a scene from here.
    """
    polygon = list(prism.polygon_xz_m)
    if polygon_area_m2(polygon) < 0:
        polygon.reverse()
    n = len(polygon)
    vertices = [(x, prism.y_min_m, z) for x, z in polygon] + [(x, prism.y_max_m, z) for x, z in polygon]
    faces = []
    for i in range(n):  # the sides: the outward normal of the edge a -> b is (dz, 0, -dx)
        j = (i + 1) % n
        faces += [(i, n + i, n + j), (i, n + j, j)]
    for a, b, c in triangulate(polygon):  # the caps at y_min (normal -y) and at y_max (normal +y)
        faces += [(a, b, c), (n + a, n + c, n + b)]
    v = np.array(vertices)
    normals = []
    for a, b, c in faces:
        normal = np.cross(v[b] - v[a], v[c] - v[a])
        normal /= np.linalg.norm(normal)
        normals.append(normal)
    return v, faces, np.array(normals)


def write_prism(path: pathlib.Path, prism: Prism) -> None:
    """Write a prism as a closed OBJ mesh in metres with one outward normal per face (`prism_mesh`)."""
    vertices, faces, normals = prism_mesh(prism)
    lines = [f"v {x!r} {y!r} {z!r}" for x, y, z in vertices.tolist()]
    rows = []
    for k, ((a, b, c), normal) in enumerate(zip(faces, normals, strict=True), start=1):
        lines.append("vn {!r} {!r} {!r}".format(*normal.tolist()))
        rows.append(f"f {a + 1}//{k} {b + 1}//{k} {c + 1}//{k}")
    path.write_text("\n".join(lines + rows) + "\n", encoding="utf-8")
