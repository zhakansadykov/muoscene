"""What a scene and a scan refuse: descriptions that cannot be built and lines of sight that mean nothing.

Sizes are those of `tests/regression/legacy_cases.py` (wall, banks, niches, pile cap) and of
`mulder_cases.py` (reference void and the depth of its niche, which is the layer top here); piles
have the radius and the length of the invented pier of the spike. Bodies that must cross, touch
or overlap are these same boxes, moved. Nothing here needs physics tables.
"""

import calzone
import numpy as np
import pytest

from muoscene import geometry as G
from muoscene import materials as M
from muoscene.geometry import builder, checks
from muoscene.geometry.primitives import Prism, polygon_area_m2
from tests.regression import legacy, legacy_cases
from tests.regression import mulder_cases as mc

NICHE = legacy.Niche()
DEPTH = legacy_cases.LOWER
CAP_LO, CAP_HI = legacy_cases.pile_cap(0.0)
CAP = G.Box(CAP_LO, CAP_HI)
VOID = G.Sphere(mc.REFERENCE_VOID[:3], mc.REFERENCE_VOID[3])  # 3 m from the wall, 3 m deep, radius 0.5 m
SOIL = M.loam(1.60, 15.0)
RISING = legacy_cases.BANKS["rising"]
PROFILE = [(0.0, -mc.GROUND[2]), (0.0, 0.0), (mc.GROUND[0], 0.0)]
UP = np.array([0.0, 0.0, 1.0])
CLAY_TOP = -mc.VOID_DEPTH  # m: a layer top 10 m below the brow


def ground(profile=PROFILE, layers=(), name="ground", y=(-mc.GROUND[1] / 2, mc.GROUND[1] / 2), **options):
    options.setdefault("material", SOIL)
    return G.Ground(name, profile, -mc.GROUND[2], y[0], y[1], layers=layers, **options)


def niche(depth=DEPTH, name="niche", x=0.0):
    """The niche of the preliminary calculation, `depth` below the brow, its mouth at `x`."""
    lo, hi = (x, -NICHE.w / 2, -depth - mc.EPS), (x + NICHE.xn, NICHE.w / 2, -depth + NICHE.hn)
    return G.Cut(name, G.Box(lo, hi))


def scene(*extra, layers=(), **options):
    options.setdefault("altitude_m", mc.ALTITUDE)
    return G.Scene([ground(layers=layers), niche(), *extra], **options)


# ------------------------------------------------ shapes


def test_signed_distances_of_shapes():
    """Negative inside, positive outside, zero on the surface; exact for box, cylinder, sphere and prism."""
    assert CAP.signed_distance_m((12.0, 0.0, -1.0)) == pytest.approx(-1.0)  # 1 m from top and bottom
    assert CAP.signed_distance_m((16.0, 0.0, -1.0)) == pytest.approx(1.0)
    assert CAP.signed_distance_m((16.0, 4.0, 1.0)) == pytest.approx(np.sqrt(3.0))  # off a corner
    assert CAP.signed_distance_m((9.0, 0.0, -1.0)) == 0.0
    pile = G.Cylinder((12.0, 0.0), -10.0, -2.0, 0.4)
    assert pile.signed_distance_m((12.0, 0.0, -6.0)) == pytest.approx(-0.4)
    assert pile.signed_distance_m((12.0, 0.0, -2.1)) == pytest.approx(-0.1)
    assert pile.signed_distance_m((13.0, 0.0, -6.0)) == pytest.approx(0.6)
    assert pile.signed_distance_m((12.7, 0.0, -1.6)) == pytest.approx(0.5)  # 0.3 out, 0.4 up
    assert VOID.signed_distance_m((3.0, 0.0, -3.0)) == pytest.approx(-0.5)
    assert VOID.signed_distance_m((3.0, 2.0, -3.0)) == pytest.approx(1.5)
    bank = ground([(0.0, -400.0), (0.0, 0.0), (RISING["xr"], RISING["rise"]), (1500.0, RISING["rise"])])
    (piece,) = bank.pieces()
    assert isinstance(piece.shape, Prism) and polygon_area_m2(piece.shape.polygon_xz_m) > 0
    assert piece.shape.signed_distance_m((100.0, 0.0, -1.5)) == pytest.approx(-3.0)  # the brow is at 1.5 m
    assert piece.shape.signed_distance_m((-2.0, 0.0, -10.0)) == pytest.approx(2.0)  # in front of the wall
    assert piece.shape.signed_distance_m((100.0, 1502.0, -10.0)) == pytest.approx(
        2.0
    )  # beyond the face across y
    assert piece.shape.signed_distance_m((4.0, 0.0, 0.75)) == pytest.approx(0.0, abs=1e-12)  # on the slope
    assert G.Box.centred((12.0, 0.0, -1.0), (6.0, 6.0, 2.0)) == CAP
    assert CAP.clip_z(-1.0, 5.0) == G.Box((9.0, -3.0, -1.0), (15.0, 3.0, 0.0))
    assert pile.clip_z(-6.0, 0.0) == G.Cylinder((12.0, 0.0), -6.0, -2.0, 0.4)


@pytest.mark.parametrize(
    "make",
    [
        lambda: G.Box((0, 0, 0), (1, 1, 0)),  # no height
        lambda: G.Box((0, 0, 0), (1, float("nan"), 1)),
        lambda: G.Box((0, 0), (1, 1)),
        lambda: G.Cylinder((0, 0), -2.0, -8.0, 0.4),  # top below bottom
        lambda: G.Cylinder((0, 0), -8.0, -2.0, 0.0),
        lambda: G.Sphere((0, 0, 0), -1.0),
        lambda: G.Anomaly("void", VOID, -1.5),  # less than a void
        lambda: G.Anomaly("", VOID, -1.0),
        lambda: G.pile_grid((0, 0), -2.0, 8.0, 0.4, 3, 3, 0.7),  # piles 0.8 m wide, 0.7 m apart
        lambda: G.pile_grid((0, 0), -2.0, 8.0, 0.4, 0, 3, 2.2),
        lambda: G.pile_grid((0, 0), -2.0, -8.0, 0.4, 3, 3, 2.2),
    ],
)
def test_shapes_are_checked(make):
    with pytest.raises(ValueError):
        make()


@pytest.mark.parametrize(
    "make",
    [
        lambda: G.Structure("cap", VOID, M.CONCRETE),  # a sphere is an anomaly
        lambda: G.Structure("cap", (), M.CONCRETE),
        lambda: G.Structure("cap", CAP, "concrete"),  # a name is not a material
        lambda: G.Cut("niche", G.Cylinder((0, 0), -8.0, -2.0, 0.4)),
        lambda: G.Anomaly("void", CAP, -1.0),
        lambda: ground(layers=("clay",)),
    ],
)
def test_kinds_of_shapes_and_materials_are_checked(make):
    with pytest.raises(TypeError):
        make()


def test_pile_grid():
    """Three by two piles under the cap, x first: a pitch apart along x, two pitches along y."""
    piles = G.pile_grid((12.0, 0.0), -2.0, 8.0, 0.4, 3, 2, 2.2, 4.4)
    assert len(piles) == 6
    assert [pile.centre_xy_m for pile in piles[:3]] == [(9.8, -2.2), (9.8, 2.2), (12.0, -2.2)]
    assert all((pile.z_bottom_m, pile.z_top_m, pile.radius_m) == (-10.0, -2.0, 0.4) for pile in piles)
    (single,) = G.pile_grid((12.0, 0.0), -2.0, 8.0, 0.4, 1, 1, 0.1)  # one pile needs no pitch
    assert single.centre_xy_m == (12.0, 0.0) and single.length_m == 8.0


# ------------------------------------------------ the ground


@pytest.mark.parametrize(
    "arguments",
    [
        {"profile": [(0.0, 0.0)]},
        {"profile": [(0.0, 0.0), (0.0, -5.0)]},  # no width
        {"profile": [(0.0, 0.0), (10.0, 0.0), (5.0, 1.0)]},  # x goes back: an overhang
        {"profile": [(0.0, 0.0), (10.0, 0.0), (10.0, 0.0), (20.0, 0.0)]},  # a point twice
        {
            "profile": [(0.0, 0.0), (10.0, 0.0), (10.0, -5.0), (10.0, 1.0), (20.0, 1.0)]
        },  # three points on one x
        {"profile": [(0.0, 0.0), (10.0, -500.0), (20.0, 0.0)]},  # below the bottom of 400 m
        {"profile": [(0.0, -400.0), (20.0, -400.0)]},  # nothing above the bottom
        {"profile": [(0.0, 0.0), (10.0, float("inf"))]},
        {"y": (5.0, 5.0)},
        {"layers": (G.Layer("a", -22.0, SOIL), G.Layer("b", -12.0, SOIL))},  # from the bottom up
        {"layers": (G.Layer("a", -12.0, SOIL), G.Layer("b", -12.0, SOIL))},
        {"layers": (G.Layer("a", 0.0, SOIL),)},  # starts at the surface: that is the ground itself
        {"layers": (G.Layer("a", -400.0, SOIL),)},  # at the bottom: nothing of it
        {"density_g_cm3": 1.9},  # a soil brings its own density
        {"material": M.CONCRETE, "density_g_cm3": -2.4},
        {"name": " "},
    ],
)
def test_a_ground_is_checked(arguments):
    ground()
    with pytest.raises(ValueError):
        ground(**arguments)


def test_pieces_of_a_ground_share_their_vertices():
    """Layer tops at 3 m, through the slope of the rising brow and below it: touching pieces, no gaps.

    The areas of the pieces add up to the area under the profile to 1e-9 m2 of 6e5 m2, and the
    point where a layer top cuts the slope is the same number in both pieces.
    """
    profile = [(0.0, -400.0), (0.0, 0.0), (RISING["xr"], RISING["rise"]), (1500.0, RISING["rise"])]
    tops = (RISING["rise"] / 3, -mc.REFERENCE_VOID[0], -DEPTH)
    bank = ground(profile, layers=tuple(G.Layer(f"layer {i}", z, SOIL) for i, z in enumerate(tops)))
    pieces = bank.pieces()
    assert [piece.owner for piece in pieces] == ["ground", "layer 0", "layer 1", "layer 2"]
    assert [(piece.z_lo_m, piece.z_hi_m) for piece in pieces] == [
        (0.5, 1.5),
        (-3.0, 0.5),
        (-22.0, -3.0),
        (-400.0, -22.0),
    ]
    whole = polygon_area_m2(ground(profile).pieces()[0].shape.polygon_xz_m)
    areas = [polygon_area_m2(piece.shape.polygon_xz_m) for piece in pieces]
    assert min(areas) > 0 and sum(areas) == pytest.approx(whole, abs=1e-9)
    cut = (RISING["xr"] / 3, RISING["rise"] / 3)  # where the first layer top meets the slope
    assert cut in pieces[0].shape.polygon_xz_m and cut in pieces[1].shape.polygon_xz_m
    np.testing.assert_allclose(bank.surface_m([0.0, 4.0, 8.0, 700.0]), [0.0, 0.75, 1.5, 1.5], atol=1e-15)
    assert bank.surface_m([0.0], upper=False)[0] == -400.0 and np.isnan(bank.surface_m([-1.0, 1501.0])).all()
    assert bank.surface_range_m(0.0, 4.0) == (0.0, 0.75)  # the wall at its end does not count
    assert bank.surface_range_m(-5.0, 4.0) == (0.0, 0.75) and bank.surface_range_m(-5.0, -1.0)[0] != -5.0


def test_triangulation_covers_a_polygon_that_is_not_convex():
    """A trench between two banks: 8 vertices, 6 triangles, all counter-clockwise, the same area."""
    polygon = [(-20, -30), (20, -30), (20, 1.5), (0, 0), (0, -12), (-6, -12), (-6, 0), (-20, 1.5)]
    triangles = builder.triangulate(polygon)
    areas = [polygon_area_m2([polygon[i] for i in triangle]) for triangle in triangles]
    assert len(triangles) == len(polygon) - 2 and min(areas) > 0
    assert sum(areas) == pytest.approx(polygon_area_m2(polygon), rel=1e-14)
    with pytest.raises(ValueError, match="not a simple polygon"):
        builder.triangulate(polygon[::-1])  # clockwise


# ------------------------------------------------ bodies against each other


def test_names_are_unique_and_air_is_taken():
    cap = G.Structure("cap", CAP, M.CONCRETE)
    for bodies, message in (
        ([ground(), cap, cap], "two bodies are named 'cap'"),
        ([ground(layers=(G.Layer("cap", CLAY_TOP, SOIL),)), cap], "two bodies are named 'cap'"),
        ([ground(name="air")], "air of the scene"),
        ([ground(), G.Cut("ground", CAP)], "two bodies are named 'ground'"),
    ):
        with pytest.raises(ValueError, match=message):
            G.Scene(bodies, altitude_m=0.0)
    with pytest.raises(TypeError, match="Ground, Structure, Cut and Anomaly"):
        G.Scene([ground(), CAP], altitude_m=0.0)  # a shape is not a body
    with pytest.raises(ValueError, match="at least one ground or structure"):
        G.Scene([niche()], altitude_m=0.0)
    with pytest.raises(ValueError, match="margin_m"):
        scene(margin_m=0.0)
    with pytest.raises(ValueError, match="both latitude_deg and longitude_deg"):
        scene(latitude_deg=43.35)
    with pytest.raises(TypeError, match="air of a scene"):
        scene(air=SOIL)


def test_structures_touch_but_do_not_overlap():
    """A pile under the cap touches it; a pile that reaches 0.5 m into it overlaps; so do two caps."""
    under = G.Structure("pile", G.Cylinder((12.0, 0.0), -10.0, CAP_LO[2], 0.4), M.CONCRETE)
    beside = G.Structure("second cap", G.Box((15.0, -3.0, -2.0), (21.0, 3.0, 0.0)), M.CONCRETE)
    built = scene(G.Structure("cap", CAP, M.CONCRETE), under, beside)
    built.check()
    for other in (
        G.Structure("pile", G.Cylinder((12.0, 0.0), -10.0, CAP_LO[2] + 0.5, 0.4), M.CONCRETE),
        G.Structure("second cap", G.Box((14.0, -3.0, -2.0), (20.0, 3.0, 0.0)), M.CONCRETE),
        G.Structure("piles", G.pile_grid((30.0, 0.0), -2.0, 8.0, 0.4, 2, 1, 0.8)[:1] * 2, M.CONCRETE),
    ):
        with pytest.raises(ValueError, match="overlap"):
            scene(G.Structure("cap", CAP, M.CONCRETE), other)
    assert checks.overlap(CAP, G.Cylinder((15.3, 0.0), -2.0, 0.0, 0.4))  # 0.1 m into the side of the cap
    assert not checks.overlap(CAP, G.Cylinder((15.4, 0.0), -2.0, 0.0, 0.4))
    assert not checks.overlap(CAP, G.Cylinder((15.3, 3.3), -2.0, 0.0, 0.4))  # off the corner: 0.42 m away


def test_a_structure_lies_in_the_ground_or_in_the_air():
    """Flush with the brow or standing on it; one that crosses the brow is to be given in two parts."""
    block = G.Box((CAP_LO[0], CAP_LO[1], 0.0), (CAP_HI[0], CAP_HI[1], 2.0))
    scene(G.Structure("cap", CAP, M.CONCRETE), G.Structure("block", block, M.CONCRETE))
    for shape in (
        G.Box((9.0, -3.0, -1.0), (15.0, 3.0, 1.0)),  # the cap raised by 1 m: through the brow
        G.Box((-3.0, -3.0, -2.0), (3.0, 3.0, 0.0)),  # moved to the wall: through it
        G.Box((9.0, 1497.0, -2.0), (15.0, 1503.0, 0.0)),  # moved along it: through a limit of the model
    ):
        with pytest.raises(ValueError, match="crosses the surface or a limit of ground 'ground'"):
            scene(G.Structure("cap", shape, M.CONCRETE))
    # under a rising brow a cap whose top is level stands half in the slope
    sloping = ground([(0.0, -400.0), (0.0, 0.0), (RISING["xr"], RISING["rise"]), (1500.0, RISING["rise"])])
    with pytest.raises(ValueError, match="crosses the surface"):
        G.Scene(
            [sloping, G.Structure("cap", G.Box((2.0, -3.0, -2.0), (8.0, 3.0, 0.5)), M.CONCRETE)],
            altitude_m=0.0,
        )
    # the cap against the wall, from the inside
    assert scene(G.Structure("cap", G.Box((0.0, -3.0, -2.0), (6.0, 3.0, 0.0)), M.CONCRETE))


def test_cuts_do_not_cross_structures_or_each_other():
    """A niche half in the face of the cap is refused; one that ends at that face is cut from the ground."""
    cap = G.Structure("cap", CAP, M.CONCRETE)
    level = -CAP_LO[2] / 2  # 1 m below the brow, the middle of the cap
    with pytest.raises(ValueError, match="cut 'pit' crosses structure 'cap'"):
        scene(cap, niche(level, "pit", x=CAP_LO[0] - NICHE.xn / 2))
    with pytest.raises(ValueError, match="cuts 'niche' and 'second niche' overlap"):
        scene(niche(DEPTH, "second niche", x=NICHE.xn / 2))
    beside = niche(level, "pit", x=CAP_LO[0] - NICHE.xn)  # touches the cap
    built = scene(cap, beside, niche(legacy_cases.UPPER, "upper niche"))
    built.check()
    inside = [(CAP_LO[0] - NICHE.xd, 0.0, -level + 0.1), (NICHE.xd, 0.0, -legacy_cases.UPPER + 0.1)]
    assert list(built.locate(inside)) == ["pit", "upper niche"]


@pytest.mark.parametrize(
    ("centre", "radius", "message"),
    [
        ((3.0, 0.0, -0.3), 0.5, "crosses a boundary of 'ground'"),  # through the brow
        ((0.3, 0.0, -3.0), 0.5, "crosses a boundary of 'ground'"),  # through the wall
        ((3.0, 0.0, CLAY_TOP - 0.2), 0.5, "crosses a boundary of 'ground', 'clay'"),  # through a layer top
        ((3.0, 0.0, 3.0), 0.5, "lies in the air"),
        ((9.2, 0.0, -1.0), 0.5, "crosses a boundary of 'cap'"),  # half in the cap
        ((12.0, 0.0, -1.0), 1.5, "crosses a boundary of 'cap'"),  # larger than the cap is thick
        ((0.9, 0.0, -DEPTH), 0.5, "reaches into cut 'niche'"),
        ((3.0, 0.7, -3.0), 0.5, "anomalies 'second' and 'void' overlap"),
    ],
)
def test_an_anomaly_lies_inside_one_host(centre, radius, message):
    """An anomaly that crosses a surface, a layer top, a structure, a cut or another anomaly is refused."""
    cap = G.Structure("cap", CAP, M.CONCRETE)
    layers = (G.Layer("clay", CLAY_TOP, SOIL),)
    void = G.Anomaly("void", VOID, -1.0)
    with pytest.raises(ValueError, match=message):
        scene(cap, void, G.Anomaly("second", G.Sphere(centre, radius), -0.2), layers=layers)


def test_an_anomaly_may_touch_the_boundary_of_its_host():
    """Half a metre under the top of a layer with a radius of half a metre: inside, and Geant4 agrees."""
    layers = (G.Layer("clay", CLAY_TOP, SOIL),)
    touching = G.Anomaly("zone", G.Sphere((3.0, 0.0, CLAY_TOP - 0.5), 0.5), -0.2)
    beside = G.Anomaly("void", G.Sphere((3.0, 1.0, CLAY_TOP - 0.5), 0.5), -1.0)  # and touches the first
    in_cap = G.Anomaly("bubble", G.Sphere((12.0, 0.0, -1.0), 1.0), -1.0)  # as thick as the cap
    built = scene(G.Structure("cap", CAP, M.CONCRETE), touching, beside, in_cap, layers=layers)
    assert built.hosts == {"zone": "clay", "void": "clay", "bubble": "cap"}
    built.check()


def test_grounds_stand_side_by_side():
    """Two grounds with a gap of 6 m between their walls; a second ground on the first is refused."""
    left = ground([(-mc.GROUND[0], 0.0), (-6.0, 0.0), (-6.0, -400.0)], name="far bank")
    built = G.Scene([left, ground()], altitude_m=0.0)
    assert list(built.locate([(-10.0, 0.0, -3.0), (10.0, 0.0, -3.0), (-3.0, 0.0, -3.0)])) == [
        "far bank",
        "ground",
        "air",
    ]
    with pytest.raises(ValueError, match="grounds 'ground' and 'second' overlap in plan"):
        G.Scene([ground(), ground([(1000.0, 0.0), (2000.0, 0.0)], name="second")], altitude_m=0.0)


# ------------------------------------------------ what a scan refuses


@pytest.mark.parametrize(
    ("position", "direction", "message"),
    [
        ((NICHE.xd, 0.0, 2.0), UP, "outside the box of the scene"),  # 1 m of air above the brow
        ((-2.0, 0.0, -DEPTH), UP, "outside the box"),
        ((NICHE.xd, 0.0, float("nan")), UP, "position_m must be finite"),
        ((NICHE.xd, 0.0), UP, "three coordinates"),
        ((NICHE.xd, 0.0, -DEPTH), (0.0, float("nan"), 1.0), "direction must be finite"),
        ((NICHE.xd, 0.0, -DEPTH), (0.0, 0.0, 2.0), "unit vectors"),  # Mulder would halve the lengths
        ((NICHE.xd, 0.0, -DEPTH), (0.0, 0.0, 0.0), "unit vectors"),
        ((NICHE.xd, 0.0, -DEPTH), (0.6, 0.8), "three components"),
        ((NICHE.xd, 0.0, -DEPTH), (0.0, 0.0, -1.0), "towards the sky"),
        ((NICHE.xd, 0.0, -DEPTH), (1.0, 0.0, 0.0), "towards the sky"),
        ((NICHE.xd, 0.0, -DEPTH - mc.EPS), UP, "on a boundary of 'niche'"),  # on the floor of the niche
        ((NICHE.xd, 0.0, -DEPTH - mc.EPS + 1e-7), UP, "on a boundary of 'niche'"),
        ((0.0, 0.0, -DEPTH + 0.1), UP, "on a boundary of"),  # in the mouth of the niche
        ((5.0, 0.0, CLAY_TOP), UP, "on a boundary of"),  # on the top of the layer
        ((5.0, 0.0, 0.0), UP, "on a boundary of 'ground'"),  # on the brow
    ],
)
def test_a_scan_refuses_points_and_directions_it_cannot_answer(position, direction, message):
    """Mulder crashes on a point outside its geometry, hangs on a nan and scales with a long vector."""
    built = scene(layers=(G.Layer("clay", CLAY_TOP, SOIL),))
    with pytest.raises(ValueError, match=message):
        built.scan(position, direction)


def test_the_detector_keeps_off_boundaries_by_a_micron():
    """1 mm above the floor of its niche, as in the frozen cases, or 10 um: fine; 0.1 um is refused."""
    built = scene()
    distance, body = built.distance_to_boundary_m(mc.detector(DEPTH))
    assert (round(distance, 9), body) == (mc.EPS, "niche")
    assert built.distance_to_boundary_m((5.0, 0.0, -3.0)) == (3.0, "ground")
    close = built.scan((NICHE.xd, 0.0, -DEPTH - mc.EPS + 1e-5), UP)
    assert close.length_in_m("ground") == pytest.approx(DEPTH - NICHE.hn, abs=1e-9)
    assert checks.BOUNDARY_M == 1e-6
    # a cut that reaches out of the wall: inside it the plane of the wall is no boundary
    porch = niche(
        legacy_cases.UPPER, "porch", x=-NICHE.xn / 2
    )  # the upper niche, half of it in front of the wall
    middle = (0.0, 0.0, -legacy_cases.UPPER + NICHE.hn / 2)
    table = scene(porch).scan(middle, (-0.6, 0.0, 0.8))  # out through its roof, over the canyon
    assert table.length_in_m("porch") == pytest.approx(NICHE.hn / 2 / 0.8, abs=1e-9)
    assert table.length_in_m("ground") == 0.0


def test_lines_of_sight_cut_by_the_limits_of_the_model_are_refused():
    """A ground that ends 15 m behind the wall and 3 m to each side: what leaves through its ends is no path.

    Vertical and steep lines of sight leave through the brow and pass; a line that reaches an
    end, a face across y or the bottom inside the ground is refused with its elevation.
    """
    narrow = ground([(0.0, -400.0), (0.0, 0.0), (CAP_HI[0], 0.0)], y=(CAP_LO[1], CAP_HI[1]))
    built = G.Scene([narrow, niche()], altitude_m=mc.ALTITUDE)
    origin = mc.detector(DEPTH)
    steep = np.array([[0.0, 0.0, 1.0], [0.5, 0.0, 1.0], [0.0, 0.1, 1.0], [-1.0, 0.0, 1.0], [-3.0, 0.2, 1.0]])
    steep /= np.linalg.norm(steep, axis=1)[:, None]
    table = built.scan(origin, steep)
    assert table.length_in_m("ground")[0] == pytest.approx(DEPTH - NICHE.hn, abs=1e-9)
    np.testing.assert_array_equal(table.length_in_m("ground")[3:], 0.0)  # out through the mouth of the niche
    for direction, elevation in (
        ((1.0, 0.0, 1.0), "45.00"),
        ((0.0, 1.0, 1.0), "45.00"),
        ((3.0, 0.0, 1.0), "18.43"),
    ):
        u = np.array(direction) / np.linalg.norm(direction)
        with pytest.raises(
            ValueError, match=f"1 of 1 lines of sight cross a limit .* {elevation} to {elevation}"
        ):
            built.scan(origin, u)
    with pytest.raises(ValueError, match="2 of 7 lines of sight"):
        built.scan(origin, np.concatenate([steep, [[0.6, 0.0, 0.8], [0.0, -0.6, 0.8]]]))
    below = (5.0, 0.0, -mc.GROUND[2] - 0.5)  # under the bottom, inside the box of the scene
    with pytest.raises(ValueError, match="cross a limit of the model inside ground 'ground'"):
        built.scan(below, UP)
    # the same lines in the wide ground of the frozen cases are fine: the wall is part of the profile
    inward = np.array([[1.0, 0.0, 1.0]]) / np.sqrt(2.0)
    wide = scene().scan(origin, inward)
    np.testing.assert_allclose(
        wide.length_in_m("ground"), legacy.soil_length(legacy.Model(), inward, DEPTH), rtol=0, atol=1e-9
    )
    assert not checks.cut_by_limits(ground(), origin[None, :], legacy_cases.direction_set()[:1]).any()


# ------------------------------------------------ the check of Calzone


def test_the_overlap_check_of_calzone():
    """Touching volumes pass; two boxes that share half their volume are reported as a ValueError."""
    scene(G.Structure("cap", CAP, M.CONCRETE), G.Anomaly("void", VOID, -1.0)).check(resolution=1000)
    box = {"box": [600.0, 600.0, 200.0], "material": "G4_WATER"}
    world = {"box": [4000.0] * 3, "material": "G4_AIR", "A": box, "B": box | {"position": [300.0, 0.0, 0.0]}}
    with pytest.raises(ValueError, match="overlapping volumes"):
        checks.check_overlaps(calzone.Geometry({"World": world}))
    checks.check_overlaps(calzone.Geometry({"World": world | {"B": box | {"position": [600.0, 0.0, 0.0]}}}))


def test_meshes_are_traversed_by_voxels():
    """The other algorithm of Calzone, "bvh", loses path where two meshes touch (`docs/notes/geometry.md`)."""
    assert builder.MESH_ALGORITHM == "voxels"
