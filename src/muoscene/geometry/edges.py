"""geometry.edges — lines of sight that meet an edge or a vertex of a mesh, moved off it.

The pieces of a ground and of its layers are closed meshes of triangles (`builder.prism_mesh`):
every side face of a prism is two triangles that share a diagonal, and the two caps are
triangulated. Where a line of sight meets an edge or a vertex of these triangles exactly,
Calzone 1.2.2 and Mulder 0.3.7 fail without an error (roadmap 0.1-10b, `docs/notes/geometry.md`):

- on macOS the flux of Mulder along a line that leaves the ground through the diagonal of a face
  never returns (the navigator of Geant4 turns in its check of overlaps), and the scan gives
  no path in the ground along it;
- on every system a line that touches a concave edge of a mesh — the foot of a canyon wall,
  where the line passes from the ground under the floor into the ground of the bank — gets the
  path beyond the edge as air: 29.48 m of solid instead of 51.84 m in `tests/unit/test_geometry_edges.py`.

Such lines do not come by chance but from symmetric scenes: a detector at the centre of a square
ground looks through the diagonal of its top along every direction with tx = −ty. The trap is
narrow: on macOS the flux along a line of `examples/slab` 1e-14 m off the diagonal never returns
and 1e-13 m off it is right; on Linux a line 1e-13 m off a concave edge gets the wrong path and
1e-12 m off it the right one — a thousand times below `TOLERANCE_M`.

`Scene.scan` and `muoscene.site.Site.differential_flux` therefore give Mulder every line of
sight through `move_off_edges` (decided on 27 November 2026): a line that passes within
`TOLERANCE_M` (1 nm) of an edge or a vertex of a mesh has its point moved by 3 nm across the line
and the edge; when that lands near another edge — at a vertex, or where the diagonals of two
touching faces cross — the next of 16 fixed moves of 3 and 9 nm across the line is taken. The
direction is kept, and the lines and the order of the moves are the same on every system. A line
that runs at less than `GRAZING` (0.57 degrees) to a face at such an edge, or along the edge, is
refused: there a move of a few nanometres changes the path by more than a micrometre, by
millimetres for a line nearly in the face, and no path is the right one.

What the move costs: where the moved line crosses a face, its path changes by at most the size
of the move over |u · n| of that face — 9 nm / 0.01 = 0.9 µm at the edge it was moved off, and
3 nm over the sine of the elevation at a flat top. A point stays at least 1 µm − 9 nm from a
boundary (`checks.BOUNDARY_M`). In `examples/slab` 114 of the 20 388 lines of the map of the film
under 2 m and 106 of the film under 20 m are moved, all by 3 nm; the counts of the film under 2 m
are bit for bit those of before (its moves are horizontal) and those of 28 bins of the film under
20 m change by at most 3.4e-11. The lines of `examples/canyon_demo` pass no closer than 0.3 mm
to an edge, and its maps are bit for bit those of before.

The edges are those of every mesh as it is written for Calzone, the shared faces of touching
pieces and the parts of meshes inside a cut included: the top of a ground and of each layer, the
walls and the slopes of its profile, the faces where pieces meet at the top of a layer (the two
pieces split such a face along crossing diagonals), and the bottom, the faces across y and the
faces under the ends of the profile, which are limits of the model. Boxes, cylinders and spheres
are solids of Geant4, not meshes: lines through the edges and a vertex of a box structure and
through the edges of a niche answered right on macOS; cylinders and spheres were not tried.

Units: metres. Layer 1.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

TOLERANCE_M = 1e-9
"""A line of sight closer than this to an edge or a vertex of a mesh is moved off it."""
SHIFTS_M = (3e-9, 9e-9)
"""Sizes of the moves, tried in this order: well above the tolerance and below any length that matters."""
ANGLES = 8
"""Directions of a move across a line: first across the line and the edge, both ways, then the others."""
GRAZING = 0.01
"""A line within the tolerance of an edge whose |u · n| for a face at that edge is below this is
refused: 0.57 degrees to the face; a move of 9 nm then shifts the path by up to 0.9 µm."""
CHUNK = 65536
"""Lines of sight examined at a time: bounds the memory of the screening."""
SCREEN = 10.0
"""Candidates are screened against SCREEN × TOLERANCE_M with the direction of their group of edges,
then checked against TOLERANCE_M with their own: covers the rounding of the grouping, 3e-9 m at
3 km."""


@dataclass(frozen=True)
class EdgeGroup:
    """Edges of one direction.

    direction   unit vector, its first nonzero component positive
    start_m     one end of each edge, m; shape (m, 3)
    unit        the direction of each edge itself; shape (m, 3)
    length_m    shape (m,)
    normals     outward unit normals of the two triangles at each edge; shape (m, 2, 3)
    """

    direction: np.ndarray
    start_m: np.ndarray
    unit: np.ndarray
    length_m: np.ndarray
    normals: np.ndarray


@dataclass(frozen=True)
class MeshEdges:
    """The edges of the meshes of a scene, in its frame (m), grouped by direction."""

    groups: tuple[EdgeGroup, ...]

    @property
    def count(self) -> int:
        """Number of edges."""
        return sum(len(group.length_m) for group in self.groups)


def mesh_edges(meshes) -> MeshEdges:
    """Edges of closed meshes given as `builder.prism_mesh` returns them: (vertices, triangles, normals)."""
    rows: dict[tuple, list] = {}
    for vertices, faces, normals in meshes:
        vertices = np.asarray(vertices, dtype=float)
        at: dict[tuple[int, int], list[np.ndarray]] = {}
        for (a, b, c), normal in zip(faces, normals, strict=True):
            for i, j in ((a, b), (b, c), (c, a)):
                at.setdefault((min(i, j), max(i, j)), []).append(np.asarray(normal, dtype=float))
        for (i, j), found in at.items():
            p, q = vertices[i], vertices[j]
            e = q - p
            first = e[np.flatnonzero(e)[0]]
            if first < 0:
                p, q, e = q, p, -e
            length = float(np.linalg.norm(e))
            e = e / length
            pair = (found + found)[:2]  # a closed mesh has two triangles at every edge
            key = tuple(np.round(e, 12).tolist())
            rows.setdefault(key, []).append((p, e, length, np.array(pair)))
    groups = []
    for key in sorted(rows):
        found = rows[key]
        groups.append(
            EdgeGroup(
                direction=np.array(key) / np.linalg.norm(key),
                start_m=np.array([r[0] for r in found]),
                unit=np.array([r[1] for r in found]),
                length_m=np.array([r[2] for r in found]),
                normals=np.array([r[3] for r in found]),
            )
        )
    return MeshEdges(tuple(groups))


def _distance(o, u, a, e, length):
    """Distance (m) of rays o + t u (t >= 0) from segments a + s e (0 <= s <= length); arrays (k, 3), (k,)."""
    w = o - a
    ue = np.einsum("ij,ij->i", u, e)
    wu = np.einsum("ij,ij->i", w, u)
    we = np.einsum("ij,ij->i", w, e)
    den = 1.0 - ue * ue
    skew = den > 1e-24
    safe = np.where(skew, den, 1.0)
    # the closest points of the two infinite lines, then clamped to the ray and the segment
    s = np.where(skew, (we - ue * wu) / safe, np.clip(we, 0.0, length))
    s = np.clip(s, 0.0, length)
    t = np.maximum(s * ue - wu, 0.0)
    s = np.clip(we + t * ue, 0.0, length)
    t = np.maximum(s * ue - wu, 0.0)
    gap = w + t[:, None] * u - s[:, None] * e
    return np.linalg.norm(gap, axis=1)


def near_edges(edges: MeshEdges, position_m, direction, tolerance_m: float = TOLERANCE_M):
    """Pairs (line, edge) closer than `tolerance_m`: indices of the lines, groups and edges, distances.

    position_m, direction   arrays (n, 3): points in metres, unit vectors towards the sky
    """
    o = np.asarray(position_m, dtype=float)
    u = np.asarray(direction, dtype=float)
    lines, groups, members, distances = [], [], [], []
    for g, group in enumerate(edges.groups):
        n = np.cross(u, group.direction)
        nn = np.linalg.norm(n, axis=1)
        # distance of the infinite line of sight from the infinite line of each edge, screened
        reach = np.abs(n @ group.start_m.T - np.einsum("ij,ij->i", n, o)[:, None])
        candidate = reach <= SCREEN * tolerance_m * nn[:, None]
        parallel = nn <= 1e-12
        if np.any(parallel):
            w = group.start_m[None, :, :] - o[parallel, None, :]
            off = np.linalg.norm(np.cross(w, group.direction), axis=2)
            candidate[parallel] = off <= SCREEN * tolerance_m
        i, k = np.nonzero(candidate)
        if i.size == 0:
            continue
        d = _distance(o[i], u[i], group.start_m[k], group.unit[k], group.length_m[k])
        hit = d < tolerance_m
        lines.append(i[hit])
        groups.append(np.full(int(hit.sum()), g))
        members.append(k[hit])
        distances.append(d[hit])
    if not lines:
        empty = np.zeros(0, dtype=int)
        return empty, empty, empty, np.zeros(0)
    return tuple(np.concatenate(x) for x in (lines, groups, members, distances))


def move_off_edges(edges: MeshEdges, position_m, direction) -> tuple[np.ndarray, np.ndarray]:
    """Points of lines of sight moved off the edges of the meshes, and which were moved.

    position_m, direction   arrays (n, 3) in the frame of the scene: points in metres, unit
                            vectors towards the sky
    Returns the points (a new array; the others are bit for bit those given) and a boolean array
    of shape (n,). ValueError for lines that run at less than 0.57 degrees to a face at the edge
    they meet, or that no move takes away from every edge.
    """
    o = np.array(position_m, dtype=float)
    u = np.asarray(direction, dtype=float)
    moved = np.zeros(len(o), dtype=bool)
    if not edges.groups or len(o) == 0:
        return o, moved
    grazing, stuck = [], []
    for start in range(0, len(o), CHUNK):
        part = slice(start, start + CHUNK)
        lines, groups, members, distances = near_edges(edges, o[part], u[part])
        for line in np.unique(lines):
            index = start + int(line)
            mine = lines == line
            near = list(zip(groups[mine], members[mine], distances[mine], strict=True))
            normals = np.concatenate([edges.groups[g].normals[k] for g, k, _ in near])
            if np.min(np.abs(normals @ u[index])) < GRAZING:
                grazing.append(index)
                continue
            g, k, _ = min(near, key=lambda item: item[2])
            moved_to = _move(edges, o[index], u[index], edges.groups[g].unit[k])
            if moved_to is None:
                stuck.append(index)
                continue
            o[index] = moved_to
            moved[index] = True
    if grazing or stuck:
        raise ValueError(_refusal(grazing, stuck, o, u))
    return o, moved


def _move(edges: MeshEdges, o, u, edge):
    """The first of the fixed moves of `o` across `u` that leaves every edge farther than the tolerance."""
    across = np.cross(u, edge)
    across /= np.linalg.norm(across)
    other = np.cross(u, across)
    order = (0, 4, 2, 6, 1, 5, 3, 7)  # across the line and the edge, both ways, then the rest
    for size in SHIFTS_M:
        for step in order:
            angle = 2.0 * math.pi * step / ANGLES
            point = o + size * (math.cos(angle) * across + math.sin(angle) * other)
            lines, *_ = near_edges(edges, point[None, :], u[None, :])
            if lines.size == 0:
                return point
    return None


def _refusal(grazing, stuck, o, u) -> str:
    first = (grazing or stuck)[0]
    elevation = math.degrees(math.asin(min(1.0, float(u[first, 2]))))
    parts = []
    if grazing:
        parts.append(
            f"{len(grazing)} pass within {TOLERANCE_M:g} m of an edge of a mesh at less than"
            f" {math.degrees(math.asin(GRAZING)):.2f} degrees to a face there, or along the edge"
        )
    if stuck:
        parts.append(f"{len(stuck)} pass within {TOLERANCE_M:g} m of edges that no move of 9 nm leaves")
    return (
        f"{len(grazing) + len(stuck)} lines of sight meet an edge of a mesh of the scene: "
        + "; ".join(parts)
        + f". The first starts at {tuple(round(v, 9) for v in o[first].tolist())} m at an elevation of"
        f" {elevation:.3f} degrees. Along such a line the path in the ground is not defined:"
        " move the detector by a millimetre or leave these directions out (roadmap 0.1-10b)"
    )
