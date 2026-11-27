"""`muoscene.geometry.edges`: lines of sight that meet an edge or a vertex of a mesh are moved off it.

Roadmap 0.1-10b. Geometry and scans only: no flux of Mulder is computed here, because along a
line through an edge it may never return on macOS; `test_example_slab.py` gives Mulder such
lines in a child process with a timeout. The slab and the pencil are those of `examples/slab`;
the canyon is an invented one of one ground and one layer, with a wall foot on each side.
"""

import importlib.util
import math
import pathlib
import tempfile

import numpy as np
import pytest

from muoscene import geometry as G
from muoscene import materials as M
from muoscene.detectors import DirectionGrid, EmulsionStack, Pose
from muoscene.flux import EnergyResponse
from muoscene.geometry import builder, edges
from muoscene.geometry.primitives import Prism

SOIL = M.loam(1.60, 15.0)
FILM = EmulsionStack(
    area_cm2=120.0, n_films=6, tan_max=1.0, efficiency=0.9, response=EnergyResponse.threshold(1.0)
)
PENCIL = DirectionGrid([-1.0e-3, 1.0e-3], [-1.0e-3, 1.0e-3])  # that of examples/slab
PROFILE = ((-100.0, 0.0), (0.0, 0.0), (0.0, -30.0), (40.0, -30.0), (40.0, 5.0), (60.0, 6.5), (140.0, 6.5))
BOTTOM, HALF_Y, LAYER_TOP = -80.0, 100.0, -10.0
EXAMPLES = pathlib.Path(__file__).resolve().parents[2] / "examples"


def slab():
    """The ground of examples/slab: 400 x 400 m, 60 m deep, its top at z = 0."""
    return G.Scene(
        [G.Ground("ground", [(-200.0, 0.0), (200.0, 0.0)], -60.0, -200.0, 200.0, SOIL)], altitude_m=1150.0
    )


def canyon(extra=()):
    layer = G.Layer("layer", LAYER_TOP, SOIL.material, 1.98)
    ground = G.Ground("ground", PROFILE, BOTTOM, -HALF_Y, HALF_Y, SOIL, layers=(layer,))
    return G.Scene([ground, *extra], altitude_m=1150.0)


def unit(v):
    v = np.asarray(v, dtype=float)
    return v / np.linalg.norm(v)


def distance_to_segment(o, u, a, b):
    """Distance of the ray o + t u (t >= 0) from the segment a-b: the minimum over a fine search.

    Independent of `edges`: the distance along the ray is minimised over the segment by golden
    section, the function being convex in the parameter of the segment.
    """

    def ray(p):
        return np.linalg.norm(p - o - max(float((p - o) @ u), 0.0) * u)

    lo, hi = 0.0, 1.0
    golden = (math.sqrt(5.0) - 1.0) / 2.0
    for _ in range(200):
        c, d = hi - golden * (hi - lo), lo + golden * (hi - lo)
        if ray(a + c * (b - a)) < ray(a + d * (b - a)):
            hi = d
        else:
            lo = c
    return min(ray(a + lo * (b - a)), ray(a), ray(b))


def obj_edges(scene):
    """Edges of the meshes of a scene as Calzone reads them: from the OBJ files the builder writes."""
    found = set()
    with tempfile.TemporaryDirectory() as folder:
        for i, volume in enumerate(scene._volumes):
            if not isinstance(volume.shape, Prism):
                continue
            path = pathlib.Path(folder) / f"{i}.obj"
            builder.write_prism(path, volume.shape)
            vertices = []
            for line in path.read_text().splitlines():
                if line.startswith("v "):
                    vertices.append(tuple(float(x) for x in line.split()[1:]))
                elif line.startswith("f "):
                    k = [int(item.split("/")[0]) - 1 for item in line.split()[1:]]
                    for p, q in ((k[0], k[1]), (k[1], k[2]), (k[2], k[0])):
                        found.add((i, *sorted((vertices[p], vertices[q]))))
    return found


def all_edges(scene):
    """Ends of every edge of `scene._edges`: arrays (n, 3)."""
    starts = np.concatenate([g.start_m for g in scene._edges.groups])
    ends = np.concatenate([g.start_m + g.unit * g.length_m[:, None] for g in scene._edges.groups])
    return starts, ends


def solid_path(o, u):
    """Length of the ray o + t u in the ground of `canyon` (layer included), m: analytic."""
    poly = [*PROFILE, (PROFILE[-1][0], BOTTOM), (PROFILE[0][0], BOTTOM)]
    ts = {0.0, 2.0e3}
    for y in (-HALF_Y, HALF_Y):
        if u[1]:
            ts.add((y - o[1]) / u[1])
    for (ax, az), (bx, bz) in zip(poly, poly[1:] + poly[:1], strict=True):
        ex, ez = bx - ax, bz - az
        den = u[0] * ez - u[2] * ex
        if den:
            wx, wz = ax - o[0], az - o[2]
            if -1e-12 <= (wx * u[2] - wz * u[0]) / den <= 1 + 1e-12:
                ts.add((wx * ez - wz * ex) / den)
    ts = sorted(t for t in ts if 0.0 <= t <= 2.0e3)
    total = 0.0
    for t0, t1 in zip(ts[:-1], ts[1:], strict=True):
        x, y, z = o + 0.5 * (t0 + t1) * u
        inside = False
        for (ax, az), (bx, bz) in zip(poly, poly[1:] + poly[:1], strict=True):
            if (az > z) != (bz > z) and x < ax + (z - az) * (bx - ax) / (bz - az):
                inside = not inside
        total += (t1 - t0) if inside and -HALF_Y < y < HALF_Y else 0.0
    return total


def test_the_edges_are_those_of_the_meshes_calzone_reads():
    scene = canyon()
    written = {edge[1:] for edge in obj_edges(scene)}
    starts, ends = all_edges(scene)
    ours = {tuple(sorted((tuple(a), tuple(b)))) for a, b in zip(starts.tolist(), ends.tolist(), strict=True)}
    assert len(written) == len({tuple(np.round(np.array(e).ravel(), 9)) for e in written})
    rounded = {tuple(np.round(np.array(e).ravel(), 9)) for e in ours}
    assert rounded == {tuple(np.round(np.array(e).ravel(), 9)) for e in written}
    assert scene._edges.count == len(obj_edges(scene))  # every mesh keeps its own edges


def test_the_lines_of_the_pencil_through_the_diagonal_are_moved():
    """The film "2 m" of examples/slab: 4 of its 16 pencil lines leave through the diagonal x = -y."""
    lines = FILM.lines_of_sight(Pose((0.0, 0.0, -2.0)), PENCIL)
    o, u, moved = slab().traced_lines((0.0, 0.0, -2.0), lines.direction)
    assert np.flatnonzero(moved).tolist() == [3, 6, 9, 12]
    assert np.all(lines.tx[moved] == -lines.ty[moved])
    step = o - np.array([0.0, 0.0, -2.0])
    assert np.allclose(np.linalg.norm(step[moved], axis=1), 3e-9, rtol=1e-6)
    assert np.allclose(np.einsum("ij,ij->i", step[moved], u[moved]), 0.0, atol=1e-20)
    assert np.all(step[~moved] == 0.0)
    assert np.array_equal(u, lines.direction / np.linalg.norm(lines.direction, axis=1, keepdims=True))


def test_after_the_move_every_line_is_off_every_edge():
    """The map of the film "2 m": 114 of 20 388 lines meet the diagonal and none is left near an edge."""
    scene = slab()
    u = FILM.lines_of_sight(Pose((0.0, 0.0, -2.0)), DirectionGrid.square()).direction
    o, u, moved = scene.traced_lines((0.0, 0.0, -2.0), u)
    assert int(moved.sum()) == 114 and np.all(np.abs(u[moved, 0] + u[moved, 1]) < 1e-9)
    starts, ends = all_edges(scene)
    for i in np.flatnonzero(moved):
        nearest = min(distance_to_segment(o[i], u[i], a, b) for a, b in zip(starts, ends, strict=True))
        assert nearest > edges.TOLERANCE_M, (i, nearest)
    again_o, _, again = scene.traced_lines(o, u)
    assert not again.any() and np.array_equal(again_o, o)  # a second pass moves nothing


@pytest.mark.parametrize(
    "what, origin, target",
    [
        ("left brow", (-6.0, 0.0, -9.0), (0.0, 0.0, 0.0)),
        ("layer top on the left wall", (-6.0, 0.0, -14.0), (0.0, 0.0, -10.0)),
        (
            "crossing diagonals of the layer top and the diagonal of the top",
            (-50.0, 0.0, -20.0),
            (-50.0, 0.0, -10.0),
        ),
        ("crossing diagonals, obliquely", (-55.0, -3.0, -16.0), (-50.0, 0.0, -10.0)),
        ("diagonal of a slope", (50.0, 0.0, -5.0), (50.0, 0.0, 5.75)),
        ("foot of the right wall, concave", (35.0, 0.0, -40.0), (40.0, 0.0, -30.0)),
        ("foot of the left wall, concave", (5.0, 0.0, -40.0), (0.0, 0.0, -30.0)),
        ("right brow", (45.0, 0.0, -2.0), (40.0, 0.0, 5.0)),
        ("top and the face y = y_max", (-50.0, 95.0, -5.0), (-50.0, 100.0, 0.0)),
        ("vertex of the top, x = x_min and y = y_max", (-95.0, 95.0, -5.0), (-100.0, 100.0, 0.0)),
    ],
)
def test_lines_through_edges_get_their_analytic_path(what, origin, target):
    """Moved by 3 or 9 nm, the line keeps its path within 1e-7 m: 9 nm over |u · n| = 0.24 at most.

    Without the move the feet of the walls give 29.48 m instead of 51.84 m and 22.36 m instead of
    44.72 m on Linux (Calzone 1.2.2, Mulder 0.3.7): the path in the bank beyond the edge is taken
    for air; macOS gave the same for another line through a foot, 41.2 m instead of 123.7 m.
    """
    scene = canyon()
    o = np.array(origin)
    u = unit(np.array(target) - o)
    table = scene.scan(o, u)
    solid = float(table.length_in_m("ground") + table.length_in_m("layer"))
    assert solid == pytest.approx(solid_path(o, u), abs=1e-7), what
    assert np.array_equal(table.position_m, o)  # the table keeps the point as given
    assert bool(scene.traced_lines(o, u)[2]), what


def test_a_line_that_grazes_a_face_at_an_edge_is_refused():
    scene = canyon()
    o = np.array([-0.001, 0.0, -9.5])  # up the left wall to its brow, 0.006 degrees from the wall
    with pytest.raises(ValueError, match="meet an edge of a mesh of the scene: 1 pass within 1e-09 m"):
        scene.scan(o, unit(-o))
    top = slab()
    for elevation_deg, refused in ((0.3, True), (1.0, False)):  # through the diagonal of the top
        u = unit(
            [
                math.cos(math.radians(elevation_deg)) / math.sqrt(2.0),
                -math.cos(math.radians(elevation_deg)) / math.sqrt(2.0),
                math.sin(math.radians(elevation_deg)),
            ]
        )
        o = np.array([0.0, 0.0, 0.0]) - 0.01 / u[2] * u
        if refused:
            with pytest.raises(ValueError, match="at less than 0.57 degrees to a face"):
                top.traced_lines(o, u)
        else:
            assert top.traced_lines(o, u)[2]


def test_traced_lines_keep_the_shapes_of_scan():
    scene = slab()
    u = FILM.lines_of_sight(Pose((0.0, 0.0, -2.0)), PENCIL).direction.reshape(4, 4, 3)
    points = np.array([[0.0, 0.0, -2.0], [0.37, 0.21, -2.0]])[:, None, None, :]
    o, d, moved = scene.traced_lines(points, u)
    assert o.shape == d.shape == (2, 4, 4, 3) and moved.shape == (2, 4, 4)
    assert moved[0].sum() == 4 and not moved[1].any()
    assert scene.scan(points, u).shape == (2, 4, 4)
    with pytest.raises(ValueError, match="outside the box"):
        scene.traced_lines((0.0, 0.0, 500.0), u)


def test_the_lines_of_the_canyon_demo_pass_no_closer_than_0_3_mm():
    """None of its 1.2e5 lines of sight is moved; the nearest passes 0.32 mm from an edge."""
    spec = importlib.util.spec_from_file_location("canyon_demo_run", EXAMPLES / "canyon_demo" / "run.py")
    example = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(example)
    setup = example.build()
    for row in setup.rows:
        pose = Pose(setup.enu[row.name], row.azimuth_deg, row.tilt_deg)
        u = setup.site.scene_from_enu(
            example.film(1.0).lines_of_sight(pose, DirectionGrid.square()).direction
        )
        o = np.broadcast_to(setup.site.scene_from_enu(np.array(pose.position_m)), u.shape)
        lines, *_ = edges.near_edges(setup.scene._edges, o, u, 3e-4)
        assert lines.size == 0, row.name
        assert not setup.scene.traced_lines(o, u)[2].any()
