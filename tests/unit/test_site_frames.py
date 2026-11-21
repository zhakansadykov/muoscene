"""`Site`: the azimuth of the axes of a scene, and the one place where it is applied.

A round trip does not catch a sign that is wrong both ways, so the direction of the azimuth is
checked on paths worked out by hand in the wall with a niche of `tests/regression/legacy_cases.py`,
and against the local frame of Mulder, which has the same angle under the name `azimuth`.
Sizes are those the repository already has: the bank and the niche of the preliminary
calculation, the ground of `mulder_cases.GROUND`.
"""

import math

import mulder
import numpy as np
import pytest

from muoscene import geometry as G
from muoscene import materials as M
from muoscene.site import Site
from muoscene.site.frames import cos_sin_deg
from tests.regression import legacy, legacy_cases
from tests.regression import mulder_cases as mc

NICHE = legacy.Niche()
DEPTH = legacy_cases.UPPER  # 12 m below the brow
WHERE = (mc.SITE["latitude"], mc.SITE["longitude"], mc.ALTITUDE)
ATOL_PATH = 1e-9  # m, as for the scenes of roadmap 0.1-06; met here: 2e-15 m


def bank(site, depth=DEPTH):
    """The flat bank of legacy_cases behind a wall at x = 0 of the scene, with the niche."""
    ground = G.Ground(
        "ground",
        [(0.0, -mc.GROUND[2]), (0.0, 0.0), (mc.GROUND[0], 0.0)],
        z_bottom_m=-mc.GROUND[2],
        y_min_m=-mc.GROUND[1] / 2,
        y_max_m=mc.GROUND[1] / 2,
        material=M.STANDARD_ROCK,
    )
    niche = G.Cut(
        "niche", G.Box((0.0, -NICHE.w / 2, -depth - mc.EPS), (NICHE.xn, NICHE.w / 2, -depth + NICHE.hn))
    )
    return site.scene([ground, niche])


def towards(azimuth_deg, elevation_deg):
    """Unit vector in ENU by the formula of the conventions: azimuth from north, clockwise."""
    a, e = math.radians(azimuth_deg), math.radians(elevation_deg)
    return np.array([math.cos(e) * math.sin(a), math.cos(e) * math.cos(a), math.sin(e)])


def test_axes_of_a_scene_by_hand():
    """At 0 the scene is ENU; at 90 its +y points east and its +x south; at 30 by the formulas."""
    np.testing.assert_array_equal(Site(*WHERE).scene_axes, np.eye(3))
    np.testing.assert_array_equal(Site(*WHERE, 90.0).scene_axes, [[0, -1, 0], [1, 0, 0], [0, 0, 1]])
    np.testing.assert_array_equal(Site(*WHERE, 180.0).scene_axes, [[-1, 0, 0], [0, -1, 0], [0, 0, 1]])
    half, root = 0.5, math.sqrt(3) / 2
    np.testing.assert_allclose(
        Site(*WHERE, 30.0).scene_axes, [[root, -half, 0], [half, root, 0], [0, 0, 1]], atol=1e-15
    )
    site = Site(*WHERE, 30.0)
    # the +y axis of the scene lies at the azimuth of the site, its +x axis 90 degrees further
    np.testing.assert_allclose(site.enu_from_scene([0.0, 1.0, 0.0]), towards(30.0, 0.0), atol=1e-15)
    np.testing.assert_allclose(site.enu_from_scene([1.0, 0.0, 0.0]), towards(120.0, 0.0), atol=1e-15)
    # a point 10 m east of the origin is 10 cos 30 along x of the scene and 10 sin 30 along its y
    np.testing.assert_allclose(site.scene_from_enu([10.0, 0.0, 2.0]), [10 * root, 5.0, 2.0], atol=1e-14)
    assert Site(*WHERE, -90.0).scene_azimuth_deg == 270.0
    assert cos_sin_deg(270.0) == (0.0, -1.0)


def test_round_trip_and_shapes():
    site = Site(*WHERE, 217.3)
    v = np.random.default_rng(3).normal(size=(4, 5, 3))
    np.testing.assert_allclose(site.enu_from_scene(site.scene_from_enu(v)), v, atol=5e-16, rtol=0)
    np.testing.assert_allclose(
        site.scene_from_enu(v)[..., 2], v[..., 2], atol=0, rtol=0
    )  # the vertical stays
    np.testing.assert_allclose(
        np.linalg.norm(site.scene_from_enu(v), axis=-1), np.linalg.norm(v, axis=-1), rtol=4e-16
    )
    assert site.scene_from_enu([1.0, 2.0, 3.0]).shape == (3,)


@pytest.mark.parametrize("azimuth", [0.0, 30.0, 90.0, 217.3])
def test_the_azimuth_is_that_of_the_local_frame_of_mulder(azimuth):
    """`LocalFrame(azimuth=ψ)` of Mulder turns vectors into ENU as the site does."""
    site = Site(*WHERE, azimuth)
    frame, enu = site.mulder_frame(), Site(*WHERE).mulder_frame()
    assert frame.azimuth == azimuth and enu.azimuth == 0.0
    assert (frame.latitude, frame.longitude, frame.altitude) == WHERE
    v = np.random.default_rng(4).normal(size=(6, 3))
    by_mulder = np.asarray(frame.transform(v, destination=enu, mode="vector"))
    np.testing.assert_allclose(site.enu_from_scene(v), by_mulder, atol=3e-15, rtol=0)  # round-off, met: 9e-16


def test_paths_worked_out_by_hand_in_a_turned_scene():
    """Azimuth of the scene 30 degrees: the massif lies at the azimuth 120, the canyon at 300.

    Detector on the floor of the niche (0.5 m deep, 0.3 m wide and high, detector 0.25 m from
    the face), 12 m below a flat brow. Along a line of sight the path in the ground is the
    distance to the brow less the distance to where the line leaves the niche:
    - azimuth 120, elevation 60: through the roof, (12 − 0.3) / sin 60;
    - azimuth 120, elevation 30: through the back wall, 12 / sin 30 − 0.25 / cos 30;
    - azimuth 300, elevation 60: through the roof and the lip above the niche, out of the face
      of the wall, 0.25 / cos 60 − 0.3 / sin 60;
    - azimuth 30 and 210, elevation 45, along the wall: through a side, 12 / sin 45 − 0.15 / cos 45;
    - azimuth 75, elevation 45, at 45 degrees to the wall: the side comes first, at 0.15 / 0.5,
      so 12 / sin 45 − 0.3. With the azimuth counted the other way the line would leave
      through the back wall and the path would be 0.066 m shorter.
    """
    site = Site(*WHERE, 30.0)
    scene = bank(site)
    half, root = 0.5, math.sqrt(3) / 2
    detector = np.array([NICHE.xd * root, -NICHE.xd * half, -DEPTH])  # 0.25 m along x of the scene
    np.testing.assert_allclose(site.scene_from_enu(detector), [NICHE.xd, 0.0, -DEPTH], atol=1e-15)
    s45 = math.sin(math.radians(45.0))
    cases = [
        (120.0, 60.0, (DEPTH - NICHE.hn) / root),
        (120.0, 30.0, DEPTH / half - (NICHE.xn - NICHE.xd) / root),
        (300.0, 60.0, NICHE.xd / half - NICHE.hn / root),
        (30.0, 45.0, DEPTH / s45 - NICHE.w / 2 / s45),
        (210.0, 45.0, DEPTH / s45 - NICHE.w / 2 / s45),
        (75.0, 45.0, DEPTH / s45 - NICHE.w / 2 / 0.5),
    ]
    u = np.array([towards(azimuth, elevation) for azimuth, elevation, _ in cases])
    table = site.scan(scene, detector, u)
    np.testing.assert_allclose(table.length_in_m("ground"), [c[2] for c in cases], atol=ATOL_PATH, rtol=0)
    assert [round(c[2], 4) for c in cases] == [13.51, 23.7113, 0.1536, 16.7584, 16.7584, 16.6706]
    # the table speaks ENU: it holds what was given, not what the scene was scanned with
    np.testing.assert_array_equal(table.direction, u)
    np.testing.assert_array_equal(table.position_m, np.broadcast_to(detector, u.shape))
    # the same lines in a scene whose azimuth is counted the other way give other paths
    mirrored = Site(*WHERE, -30.0)
    wrong = mirrored.scan(bank(mirrored), mirrored.enu_from_scene([NICHE.xd, 0.0, -DEPTH]), u)
    assert np.abs(wrong.length_in_m("ground") - table.length_in_m("ground")).max() > 0.5


@pytest.mark.parametrize("azimuth", [0.0, 90.0, 143.0])
def test_scan_is_the_scan_of_the_scene_along_its_own_axes(azimuth):
    site = Site(*WHERE, azimuth)
    scene = bank(site)
    rng = np.random.default_rng(8)
    t = rng.uniform(-1.2, 1.2, (300, 2))
    u_scene = np.column_stack([t, np.ones(len(t))])
    u_scene /= np.linalg.norm(u_scene, axis=1)[:, None]
    o_scene = np.array([NICHE.xd, 0.0, -DEPTH])
    own = scene.scan(o_scene, u_scene)
    ours = site.scan(scene, site.enu_from_scene(o_scene), site.enu_from_scene(u_scene))
    # 2e-13 m: the turn there and back moves a direction by 2e-16, a path of 20 m by as much times its length
    np.testing.assert_allclose(ours.length_m, own.length_m, atol=1e-11, rtol=0)
    assert ours.names == own.names and ours.shape == own.shape == (300,)
    if azimuth == 0.0:
        np.testing.assert_array_equal(ours.length_m, own.length_m)  # nothing is turned at all


def test_a_scene_of_a_site_tells_mulder_where_north_is():
    """With the geomagnetic field the flux depends on the azimuth: east against west 0.25 %.

    Scene turned by 90 degrees: its +y axis points east. The open-sky flux above 1 GeV along +y
    and along −y of the scene, through the local geometry of the scene, equals that of Mulder
    without a local geometry towards the geographic east and west. A scene that kept its y
    axis for north, as before `azimuth_deg`, would give north and south.
    """
    site = Site(*WHERE, 90.0)
    scene = bank(site)
    assert scene.azimuth_deg == 90.0 and scene.mulder.frame.azimuth == 90.0
    k = np.logspace(0.0, 3.0, 61)

    def above_1_gev(meter, **state):
        flux = np.asarray(meter.flux(energy=k, **{name: np.full(k.shape, v) for name, v in state.items()}))
        return float(np.trapezoid(flux * k, np.log(k)))

    local = mulder.Fluxmeter(geometry=scene.mulder, date=mc.DATE)
    sky = mulder.Fluxmeter(date=mc.DATE)
    point = np.array([-0.5, 0.0, 0.5])  # in the air of the scene, in front of the wall
    s = math.sin(math.radians(45.0))
    out = {}
    for name, sign, azimuth in (("east", 1.0, 90.0), ("west", -1.0, 270.0)):
        through_scene = float(
            np.trapezoid(
                np.asarray(
                    local.flux(
                        energy=k,
                        position=np.broadcast_to(point, (k.size, 3)),
                        direction=np.broadcast_to([0.0, sign * s, s], (k.size, 3)),
                    )
                )
                * k,
                np.log(k),
            )
        )
        geographic = above_1_gev(
            sky,
            latitude=WHERE[0],
            longitude=WHERE[1],
            altitude=WHERE[2] + point[2],
            azimuth=azimuth,
            elevation=45.0,
        )
        # 1e-5: the air of the scene against the atmosphere of Mulder is 2.6e-6 (docs/numbers.md)
        assert through_scene == pytest.approx(geographic, rel=1e-5)
        out[name] = through_scene
    assert out["east"] / out["west"] - 1 == pytest.approx(0.0025, abs=0.0005)


def test_scene_azimuth_changes_no_scan():
    """`Scene(azimuth_deg=...)` only tells Mulder where north is: no path changes."""
    u = np.array([towards(a, e) for a in (10.0, 100.0, 200.0, 300.0) for e in (20.0, 50.0, 80.0)])
    o = (NICHE.xd, 0.0, -DEPTH)
    paths = [bank(Site(*WHERE, azimuth)).scan(o, u).length_m for azimuth in (0.0, 40.0, 180.0)]
    np.testing.assert_array_equal(paths[0], paths[1])
    np.testing.assert_array_equal(paths[0], paths[2])
    plain = G.Scene([G.Structure("box", G.Box((0, 0, 0), (1, 1, 1)), M.CONCRETE)], altitude_m=0.0)
    assert plain.azimuth_deg == 0.0 and plain.latitude_deg is None and plain.mulder.frame.azimuth == 0.0
    with pytest.raises(ValueError, match="azimuth_deg"):
        G.Scene(
            [G.Structure("box", G.Box((0, 0, 0), (1, 1, 1)), M.CONCRETE)],
            altitude_m=0.0,
            azimuth_deg=math.nan,
        )


def test_a_scene_of_another_site_is_refused():
    site = Site(*WHERE, 30.0)
    for other in (
        Site(*WHERE, 31.0),
        Site(WHERE[0], WHERE[1], 1000.0, 30.0),
        Site(40.0, WHERE[1], WHERE[2], 30.0),
    ):
        with pytest.raises(ValueError, match="another site"):
            site.scan(bank(other), (0.2, 0.0, -DEPTH), towards(0.0, 60.0))
    ground = G.Ground("ground", [(0, -50), (0, 0), (100, 0)], -50, -100, 100, M.STANDARD_ROCK)
    with pytest.raises(ValueError, match="another site"):  # no latitude: not a scene of this site
        Site(*WHERE).scan(G.Scene([ground], altitude_m=mc.ALTITUDE), (-0.5, 0.0, 0.5), towards(0.0, 60.0))
    with pytest.raises(TypeError):
        site.scan("scene", (0.0, 0.0, 0.0), towards(0.0, 60.0))
    with pytest.raises(ValueError, match="that of the site"):
        site.scene([ground], altitude_m=0.0)
    # what the scene refuses stays refused: a direction below the horizon, a point outside
    scene = bank(site)
    with pytest.raises(ValueError, match="towards the sky"):
        site.scan(scene, site.enu_from_scene([NICHE.xd, 0.0, -DEPTH]), towards(120.0, -10.0))
    with pytest.raises(ValueError, match="outside the box"):
        site.scan(scene, (0.0, 0.0, 500.0), towards(120.0, 60.0))
    with pytest.raises(ValueError, match="three components"):
        site.scan(scene, (0.0, 0.0), towards(120.0, 60.0))


@pytest.mark.parametrize(
    "arguments",
    [
        (91.0, 79.08, 1150.0),
        (43.35, 400.0, 1150.0),
        (43.35, 79.08, 1.0e4),
        (43.35, 79.08, math.nan),
        (43.35, 79.08, 1150.0, math.inf),
        ("43", 79.08, 1150.0),
        (True, 79.08, 1150.0),
    ],
)
def test_a_site_checks_its_numbers(arguments):
    with pytest.raises(ValueError):
        Site(*arguments)


def test_a_site_is_immutable_and_passes_options_to_its_scene():
    site = Site(*WHERE, 30.0)
    with pytest.raises(AttributeError):
        site.altitude_m = 0.0
    assert site == Site(*WHERE, 390.0)
    ground = G.Ground("ground", [(0, -50), (0, 0), (100, 0)], -50, -100, 100, M.STANDARD_ROCK)
    scene = site.scene([ground], margin_m=5.0, atmosphere="midlatitude-winter")
    assert scene.hi_m[2] == 5.0 and scene.atmosphere == "midlatitude-winter"
    assert (scene.altitude_m, scene.latitude_deg, scene.longitude_deg, scene.azimuth_deg) == (
        *WHERE[2:],
        *WHERE[:2],
        30.0,
    )
