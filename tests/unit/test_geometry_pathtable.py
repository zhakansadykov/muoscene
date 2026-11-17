"""`PathTable`: opacity under density hypotheses, without a new scan.

The scene is the flat bank of `tests/regression/legacy_cases.py` behind its wall, with the lower
niche, the pile cap, the reference void of `mulder_cases.py` and a layer below the upper niche
level. Opacities are compared with what Mulder itself sums along the same lines of sight
(`scan(output="grammage")`) and with scenes rebuilt at other densities.
"""

import numpy as np
import pytest

from muoscene import geometry as G
from muoscene import materials as M
from muoscene.geometry import BodyInfo, PathTable
from muoscene.materials.backend import register
from tests.regression import legacy, legacy_cases
from tests.regression import mulder_cases as mc

NICHE = legacy.Niche()
DEPTH = legacy_cases.LOWER
CAP_LO, CAP_HI = legacy_cases.pile_cap(0.0)
VOID_CENTRE, VOID_RADIUS = np.array(mc.REFERENCE_VOID[:3]), mc.REFERENCE_VOID[3]
# Opacities of a few thousand g/cm2 are sums of a dozen products: round-off is 1e-12 g/cm2.
RTOL = 1e-12


def scene(upper=(1.60, 15.0), lower=(1.70, 10.0), concrete_g_cm3=None, drho_rel=-1.0):
    """Loam above 12 m, sand below; dry density (g/cm3) and water content (%) of each, the cap, the void."""
    profile = [(0.0, -mc.GROUND[2]), (0.0, 0.0), (mc.GROUND[0], 0.0)]
    ground = G.Ground(
        "loam",
        profile,
        -mc.GROUND[2],
        -mc.GROUND[1] / 2,
        mc.GROUND[1] / 2,
        M.loam(*upper),
        layers=(G.Layer("sand", -legacy_cases.UPPER, M.sand(*lower)),),
    )
    niche = G.Box((0.0, -NICHE.w / 2, -DEPTH - mc.EPS), (NICHE.xn, NICHE.w / 2, -DEPTH + NICHE.hn))
    bodies = [
        ground,
        G.Structure("cap", G.Box(CAP_LO, CAP_HI), M.CONCRETE, concrete_g_cm3),
        G.Cut("niche", niche),
        G.Anomaly("void", G.Sphere(VOID_CENTRE, VOID_RADIUS), drho_rel),
    ]
    return G.Scene(bodies, altitude_m=mc.ALTITUDE)


def lines():
    """The directions of legacy_cases, a cone over the cap and a cone on the void."""
    origin = mc.detector(DEPTH)
    towards_cap = (CAP_LO + CAP_HI) / 2 - origin
    towards_void = VOID_CENTRE - origin
    cones = [
        mc.cone(axis / np.linalg.norm(axis), alpha)[0]
        for axis, alpha in ((towards_cap, 0.3), (towards_void, 0.04))
    ]
    return origin, np.concatenate([legacy_cases.direction_set(), *cones])


@pytest.fixture(scope="module")
def scanned():
    built = scene()
    origin, u = lines()
    return built, built.scan(origin, u)


def grammage_g_cm2(built, media=lambda medium: True):
    """Opacity along `lines()` as Mulder sums it over the chosen media of a scene, g/cm2."""
    origin, u = lines()
    position = np.broadcast_to(origin, u.shape)
    per_medium = np.asarray(built.mulder.scan(position=position, direction=u, output="grammage"))
    chosen = [i for i, medium in enumerate(built.mulder.media) if media(medium)]
    return per_medium[:, chosen].sum(axis=1) / 10  # kg/m2 to g/cm2


def test_bodies_of_the_table(scanned):
    built, table = scanned
    assert table.names == ("loam", "sand", "cap", "niche", "void", "air")
    assert [body.kind for body in table.bodies] == ["ground", "layer", "structure", "cut", "anomaly", "air"]
    assert table.bodies == built.table_bodies
    assert table.shape == (len(lines()[1]),)
    void = table.body("void")
    assert (void.host, void.drho_rel, void.density_g_cm3) == ("loam", -1.0, 0.0)
    assert void.material == table.body("loam").material == register(M.loam(1.60, 15.0).material)
    np.testing.assert_array_equal(table.length_in_m("cap"), table.length_m[:, 2])
    with pytest.raises(ValueError, match="no body named 'clay'"):
        table.length_in_m("clay")


def test_opacity_is_density_times_length_as_mulder_sums_it(scanned):
    """The sum over bodies equals the grammage Mulder integrates along the same lines, to 1e-12."""
    built, table = scanned
    air = register(M.AIR)
    np.testing.assert_allclose(table.opacity_g_cm2(air=True), grammage_g_cm2(built), rtol=RTOL)
    solid = grammage_g_cm2(built, lambda medium: medium.material != air)
    np.testing.assert_allclose(table.opacity_g_cm2(), solid, rtol=RTOL, atol=1e-12)
    # the air of the box and of the niche: at most 30 m of it, 3 g/cm2
    in_air = table.opacity_g_cm2(air=True) - table.opacity_g_cm2()
    assert 0 < in_air.min() and in_air.max() < 3.5
    by_body = table.opacity_by_body_g_cm2()
    assert by_body.shape == (*table.shape, 6)
    np.testing.assert_allclose(by_body[:, 0], 100 * 1.60 * 1.15 * table.length_in_m("loam"), rtol=RTOL)
    np.testing.assert_array_equal(by_body[:, 4], 0.0)  # a void


def test_a_hypothesis_without_a_new_scan_equals_a_rebuilt_scene(scanned):
    """Other soils, a heavier cap and a decompaction zone in place of the void: one table, no new scan.

    The table under the hypothesis equals the scan of a scene rebuilt with those densities, and
    the grammage Mulder sums in it, to 1e-12.
    """
    _, table = scanned
    upper, lower = (1.50, 10.0), (1.75, 12.0)
    hypothesis = {"loam": 1.50 * 1.10, "sand": 1.75 * 1.12, "cap": 2.4}
    rebuilt = scene(upper, lower, concrete_g_cm3=2.4, drho_rel=-0.2)
    again = rebuilt.scan(*lines())
    np.testing.assert_array_equal(again.length_m, table.length_m)  # the geometry is the same
    under = table.opacity_g_cm2(hypothesis, {"void": -0.2})
    np.testing.assert_allclose(under, again.opacity_g_cm2(), rtol=RTOL)
    air = register(M.AIR)
    np.testing.assert_allclose(under, grammage_g_cm2(rebuilt, lambda m: m.material != air), rtol=RTOL)
    np.testing.assert_allclose(
        table.densities_g_cm3(hypothesis, {"void": -0.2}), again.densities_g_cm3(), rtol=1e-15
    )
    assert np.abs(under - table.opacity_g_cm2()).max() > 100  # g/cm2: and it is another opacity


def test_a_void_takes_its_chord_out_of_the_opacity(scanned):
    """Void against no anomaly: the opacity differs by the density of the host times the chord."""
    _, table = scanned
    origin, u = lines()
    chord = legacy.sphere_chord(u, origin, VOID_CENTRE, VOID_RADIUS)
    assert np.count_nonzero(chord) > 400
    difference = table.opacity_g_cm2(drho_rel={"void": 0.0}) - table.opacity_g_cm2()
    np.testing.assert_allclose(difference, 100 * 1.60 * 1.15 * chord, rtol=0, atol=1e-7)  # 1e-9 m of soil
    half = table.opacity_g_cm2(drho_rel={"void": -0.5}) - table.opacity_g_cm2()
    np.testing.assert_allclose(half, difference / 2, rtol=0, atol=1e-7)


def test_an_anomaly_follows_the_density_of_its_host(scanned):
    """Under a hypothesis the anomaly is its host of that hypothesis times (1 + drho_rel)."""
    _, table = scanned
    densities = dict(zip(table.names, table.densities_g_cm3({"loam": 2.0}, {"void": -0.2}), strict=True))
    assert densities["void"] == pytest.approx(2.0 * 0.8, rel=1e-15)
    assert densities["sand"] == pytest.approx(1.70 * 1.10, rel=1e-15)  # not named: as in the scene
    assert densities["air"] == densities["niche"] == pytest.approx(1.0997e-3, rel=1e-4)
    plain = dict(zip(table.names, table.densities_g_cm3(), strict=True))
    assert plain["void"] == 0.0 and plain["cap"] == 2.300
    denser = dict(zip(table.names, table.densities_g_cm3(drho_rel={"void": 0.5}), strict=True))
    assert denser["void"] == pytest.approx(1.5 * 1.60 * 1.15, rel=1e-15)


def test_opacity_by_material(scanned):
    """What the equivalent opacity of the conventions is summed over: one entry per solid material."""
    _, table = scanned
    loam, sand = register(M.loam(1.60, 15.0).material), register(M.sand(1.70, 10.0).material)
    by_material = table.opacity_by_material_g_cm2(drho_rel={"void": -0.2})
    assert set(by_material) == {loam, sand, register(M.CONCRETE)}
    total = sum(by_material.values())
    np.testing.assert_allclose(total, table.opacity_g_cm2(drho_rel={"void": -0.2}), rtol=RTOL)
    by_body = table.opacity_by_body_g_cm2(drho_rel={"void": -0.2})
    np.testing.assert_allclose(
        by_material[loam], by_body[:, 0] + by_body[:, 4], rtol=RTOL
    )  # with its anomaly
    np.testing.assert_array_equal(by_material[register(M.CONCRETE)], by_body[:, 2])


@pytest.mark.parametrize(
    ("densities", "changes"),
    [
        ({"clay": 1.9}, None),  # no such body
        ({"loam": -1.0}, None),
        ({"loam": float("nan")}, None),
        ({"void": 0.5}, None),  # an anomaly takes drho_rel
        ({"air": 1.2e-3}, None),  # the air is that of the scene
        ({"niche": 1.0}, None),
        (None, {"loam": -0.2}),  # not an anomaly
        (None, {"void": -1.5}),  # less than nothing
        (None, {"void": float("inf")}),
        (None, {"bubble": -1.0}),
    ],
)
def test_hypotheses_are_checked(scanned, densities, changes):
    _, table = scanned
    for method in (table.densities_g_cm3, table.opacity_g_cm2, table.opacity_by_material_g_cm2):
        with pytest.raises(ValueError):
            method(densities, changes)


def test_a_table_by_hand():
    """Two lines of sight through 10 m and 20 m of soil of 1.8 g/cm3 with 1 m of a void in the first."""
    bodies = [
        BodyInfo("soil", "ground", "Soil", 1.8),
        BodyInfo("void", "anomaly", "Soil", 0.0, host="soil", drho_rel=-1.0),
        BodyInfo("air", "air", "Air", 1.1e-3),
    ]
    up = [[0.0, 0.0, 1.0]] * 2
    table = PathTable(bodies, [[0.0, 0.0, -10.0], [0.0, 0.0, -20.0]], up, [[9.0, 1.0, 1.0], [20.0, 0.0, 1.0]])
    np.testing.assert_allclose(table.opacity_g_cm2(), [1620.0, 3600.0], rtol=1e-15)
    np.testing.assert_allclose(table.opacity_g_cm2(drho_rel={"void": 0.0}), [1800.0, 3600.0], rtol=1e-15)
    np.testing.assert_allclose(table.opacity_g_cm2(air=True), [1620.11, 3600.11], rtol=1e-12)
    np.testing.assert_allclose(
        table.opacity_g_cm2({"soil": 2.0}, {"void": -0.2}), [1960.0, 4000.0], rtol=1e-15
    )
    assert set(table.opacity_by_material_g_cm2()) == {"Soil"}
    assert table.shape == (2,) and table.names == ("soil", "void", "air")
    with pytest.raises(ValueError, match="read-only"):
        table.length_m[0, 0] = 0.0


@pytest.mark.parametrize(
    "arguments",
    [
        {"bodies": []},
        {"bodies": [BodyInfo("soil", "ground", "Soil", 1.8)] * 2},  # one name twice
        {"bodies": [BodyInfo("soil", "rock", "Soil", 1.8), BodyInfo("air", "air", "Air", 1.1e-3)]},
        {"bodies": [BodyInfo("soil", "ground", "Soil", 1.8), BodyInfo("void", "anomaly", "Soil", 0.0)]},
        {
            "bodies": [
                BodyInfo("soil", "ground", "Soil", 1.8, host="soil", drho_rel=0.0),
                BodyInfo("a", "air", "A", 1),
            ]
        },
        {
            "bodies": [
                BodyInfo("a", "air", "A", 1),
                BodyInfo("void", "anomaly", "A", 0.0, host="a", drho_rel=-1.0),
            ]
        },
        {"length_m": [[9.0, 1.0, 1.0]]},  # one line of sight, two points
        {"length_m": [[9.0, 1.0], [20.0, 0.0]]},  # two columns, three bodies
        {"length_m": [[9.0, -1.0, 1.0], [20.0, 0.0, 1.0]]},
        {"length_m": [[9.0, float("nan"), 1.0], [20.0, 0.0, 1.0]]},
        {"direction": [[0.0, 0.0, 1.0]]},
    ],
)
def test_a_table_refuses_what_does_not_fit(arguments):
    good = {
        "bodies": [
            BodyInfo("soil", "ground", "Soil", 1.8),
            BodyInfo("void", "anomaly", "Soil", 0.0, host="soil", drho_rel=-1.0),
            BodyInfo("air", "air", "Air", 1.1e-3),
        ],
        "position_m": [[0.0, 0.0, -10.0], [0.0, 0.0, -20.0]],
        "direction": [[0.0, 0.0, 1.0]] * 2,
        "length_m": [[9.0, 1.0, 1.0], [20.0, 0.0, 1.0]],
    }
    PathTable(**good)
    with pytest.raises(ValueError):
        PathTable(**(good | arguments))
