"""`Scene.scan`: the criteria of roadmap 0.1-06 against the frozen paths and analytic chords.

Every size is one the repository already has: the wall, the banks, the niche and the pile cap of
`tests/regression/legacy_cases.py`, the reference void of `mulder_cases.py`, the invented pier of
`docs/notes/spike/mulder_spike.py` (read from the script itself) and the soil of the example of
`docs/notes/flux.md`. Layer tops are made of the same numbers.

Path lengths do not depend on materials. The tests that put a flux through a scene use the loam
of the spike, whose physics tables the flux tests compile anyway, and the air of the library.
"""

import importlib.util
import json
import pathlib
import time

import mulder
import numpy as np
import pytest

from muoscene import geometry as G
from muoscene import materials as M
from muoscene.materials.backend import register
from muoscene.materials.composition import Material, Mixture
from tests.regression import legacy, legacy_cases
from tests.regression import mulder_cases as mc

ROOT = pathlib.Path(__file__).resolve().parents[2]
FROZEN = json.loads(
    (pathlib.Path(mc.__file__).parent / "data" / "mulder_reference.json").read_text(encoding="utf-8")
)["cases"]
ATOL_PATH = 1e-9  # m, the criterion of the roadmap; the largest difference met is 9e-12 m (piles, 401 x 401)
NICHE = legacy.Niche()
CAP_LO, CAP_HI = legacy_cases.pile_cap(
    0.0
)  # the pile cap of legacy_cases: 6 m x 6 m x 2 m, 9 m from the wall
SOIL = M.loam(1.60, 15.0)  # the soil of the example of docs/notes/flux.md: bulk density 1.84 g/cm3


def _spike():
    """The spike script as a module: the invented pier is defined there and nowhere else."""
    path = ROOT / "docs" / "notes" / "spike" / "mulder_spike.py"
    spec = importlib.util.spec_from_file_location("mulder_spike", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


SPIKE = _spike()
PIER_X0, PIER_SIZE, PIER_CAP, PILE_RADIUS, PILE_LENGTH, PILES, PILE_PITCH = SPIKE.PIER
PIER_VOID = SPIKE.PIER_VOID  # x, y, z, radius of the void in front of the pier
# Layer tops for the pier: the depth of the reference void (3 m), which cuts the piles, and the
# pile length (8 m), which cuts them again; the void in front of the pier lies between the two.
PIER_LAYERS = (-mc.REFERENCE_VOID[0], -PILE_LENGTH)


def spike_loam():
    """The loam of the spike: the `RegLoam` of the frozen cases as a muoscene material."""
    quartz = Material.from_formula("spike quartz", "SiO2", 2.65, i_ev="mulder")
    calcite = Material.from_formula("spike calcite", "CaCO3", 2.71, i_ev="mulder")
    alumina = Material.from_formula("spike alumina", "Al2O3", 3.95, i_ev="mulder")
    water = Material.from_formula("spike water", "H2O", 1.02, i_ev="mulder")  # Mulder's Water
    return Mixture.of("spike loam", {quartz: 0.60, calcite: 0.15, alumina: 0.12, water: 0.13})


def wall_profile(bank):
    """The bank of legacy_cases behind a vertical wall at x = 0, as a profile of `mulder_cases.GROUND`."""
    rise, xr = legacy_cases.BANKS[bank]["rise"], legacy_cases.BANKS[bank]["xr"]
    brow = [(xr, rise)] if rise > 0 else []
    return [(0.0, -mc.GROUND[2]), (0.0, 0.0), *brow, (mc.GROUND[0], rise)]


def niche_cut(depth):
    """The niche of the preliminary calculation, its floor `mc.EPS` below the detector."""
    return G.Cut(
        "niche", G.Box((0.0, -NICHE.w / 2, -depth - mc.EPS), (NICHE.xn, NICHE.w / 2, -depth + NICHE.hn))
    )


def wall_scene(depth, bank="flat", layer_tops=(), extra=(), material=SOIL, density_g_cm3=None, **options):
    """The scene of `mulder_cases.build_scene` through muoscene, optionally with layers and more bodies."""
    layers = tuple(G.Layer(f"layer {i + 1}", z, M.STANDARD_ROCK) for i, z in enumerate(layer_tops))
    ground = G.Ground(
        "ground",
        wall_profile(bank),
        z_bottom_m=-mc.GROUND[2],
        y_min_m=-mc.GROUND[1] / 2,
        y_max_m=mc.GROUND[1] / 2,
        material=material,
        density_g_cm3=density_g_cm3,
        layers=layers,
    )
    options.setdefault("altitude_m", mc.ALTITUDE)
    return G.Scene([ground, niche_cut(depth), *extra], **options)


def pier_scene(layer_tops=PIER_LAYERS):
    """The invented pier of the spike: pile cap, 3 x 3 piles, a void in front, a niche 20 m below the brow."""
    layers = tuple(G.Layer(f"layer {i + 1}", z, M.sand(1.70, 10.0)) for i, z in enumerate(layer_tops))
    ground = G.Ground(
        "soil", [(0.0, -200.0), (0.0, 0.0), (400.0, 0.0)], -200.0, -400.0, 400.0, SOIL, layers=layers
    )
    centre = (PIER_X0 + PIER_SIZE / 2, 0.0)
    cap = G.Box((PIER_X0, -PIER_SIZE / 2, -PIER_CAP), (PIER_X0 + PIER_SIZE, PIER_SIZE / 2, 0.0))
    piles = G.pile_grid(centre, -PIER_CAP, PILE_LENGTH, PILE_RADIUS, PILES, PILES, PILE_PITCH)
    bodies = [
        ground,
        G.Structure("cap", cap, M.CONCRETE),
        G.Structure("piles", piles, M.CONCRETE),
        G.Cut("niche", G.Box((0.0, -0.15, -SPIKE.PIER_NICHE), (0.5, 0.15, -SPIKE.PIER_NICHE + 0.3))),
        G.Anomaly("void", G.Sphere(PIER_VOID[:3], PIER_VOID[3]), drho_rel=-1.0),
    ]
    return G.Scene(bodies, altitude_m=mc.ALTITUDE), cap, piles


def slope_grid(n, t_max=1.0):
    """n x n directions on a grid of slopes up to `t_max`, as the spike scans."""
    tx, ty = np.meshgrid(np.linspace(-t_max, t_max, n), np.linspace(-t_max, t_max, n))
    norm = np.sqrt(1 + tx**2 + ty**2)
    return np.stack([tx / norm, ty / norm, 1 / norm], axis=-1)


def cylinder_chord(u, origin, cylinder):
    """Chord (m) of rays from `origin` along u through a vertical cylinder; u_z > 0."""
    ox, oy = origin[0] - cylinder.centre_xy_m[0], origin[1] - cylinder.centre_xy_m[1]
    a = u[:, 0] ** 2 + u[:, 1] ** 2
    b = ox * u[:, 0] + oy * u[:, 1]
    disc = b * b - a * (ox * ox + oy * oy - cylinder.radius_m**2)
    root = np.sqrt(np.maximum(disc, 0.0))
    with np.errstate(divide="ignore", invalid="ignore"):
        t_in = np.where(disc > 0, (-b - root) / a, np.inf)
        t_out = np.where(disc > 0, (-b + root) / a, -np.inf)
    t_bottom, t_top = (cylinder.z_bottom_m - origin[2]) / u[:, 2], (cylinder.z_top_m - origin[2]) / u[:, 2]
    return np.clip(np.minimum(t_out, t_top) - np.maximum.reduce([t_in, t_bottom, np.zeros(len(u))]), 0, None)


def path_under_profile(profile, origin, u, z_lo=-np.inf, z_hi=np.inf):
    """Length (m) of the ray origin + t u under a profile, between two elevations; u_z > 0.

    The ray is cut at its crossings with the segments of the profile and is alternately under
    and above it; it must not leave through the ends of the profile.
    """
    (ox, _, oz), (ux, _, uz) = origin, u  # u is a unit vector in space, so t is a length
    crossings = []
    for (ax, az), (bx, bz) in zip(profile[:-1], profile[1:], strict=True):
        dx, dz = bx - ax, bz - az
        det = dx * uz - ux * dz
        if det == 0:
            continue
        t = (dx * (az - oz) - dz * (ax - ox)) / det
        s = (ux * (az - oz) - uz * (ax - ox)) / det
        if t > 0 and 0 <= s <= 1:
            crossings.append(t)
    px, pz = zip(*profile, strict=True)
    under = oz < np.interp(ox, px, pz)
    total, last = 0.0, 0.0
    for t in sorted(crossings):
        if under:
            lo, hi = max(last, (z_lo - oz) / uz), min(t, (z_hi - oz) / uz)
            total += max(hi - lo, 0.0)
        under, last = not under, t
    assert not under, "the ray ends under the profile"
    return total


# ------------------------------------------------ the wall with a niche: parity with the reference


@pytest.mark.parametrize("bank", sorted(legacy_cases.BANKS))
@pytest.mark.parametrize("depth", [legacy_cases.UPPER, legacy_cases.LOWER])
def test_wall_paths_equal_the_frozen_ones(bank, depth):
    """The path in the ground behind the wall with a niche is that of `mulder_cases.wall` to 1e-9 m.

    The frozen scene is a box or a five-sided mesh written with six decimals; here it is a
    profile turned into a mesh by the builder. They agree to 1.1e-13 m.
    """
    table = wall_scene(depth, bank).scan(mc.detector(depth), legacy_cases.direction_set())
    np.testing.assert_allclose(
        table.length_in_m("ground"), FROZEN[f"wall/{bank}/{depth:.0f}m"], rtol=0, atol=ATOL_PATH
    )


@pytest.mark.parametrize("depth", [legacy_cases.UPPER, legacy_cases.LOWER])
def test_wall_paths_equal_the_preliminary_calculation(depth):
    """Under a flat brow the path is that of `legacy.soil_length` for every direction of its set."""
    u = legacy_cases.direction_set()
    table = wall_scene(depth).scan(mc.detector(depth), u)
    prelim = legacy.soil_length(legacy.Model(), u, depth)
    np.testing.assert_allclose(table.length_in_m("ground"), prelim, rtol=0, atol=ATOL_PATH)
    # what is not ground is the air of the niche and of the scene box
    assert np.all(table.length_in_m("niche") > 0)
    assert table.names == ("ground", "niche", "air")


# ------------------------------------------------ layers add up to the ground


@pytest.mark.parametrize("bank", sorted(legacy_cases.BANKS))
@pytest.mark.parametrize(
    "layer_tops",
    [
        (-mc.REFERENCE_VOID[0],),  # 3 m: above both niches
        (
            -mc.REFERENCE_VOID[0],
            -legacy_cases.UPPER - 0.1,
            -legacy_cases.LOWER - 0.1,
        ),  # between and below them
        (-legacy_cases.UPPER + NICHE.hn / 2, -legacy_cases.LOWER + NICHE.hn / 2),  # through the niches
    ],
)
def test_layers_add_up_to_the_ground(bank, layer_tops):
    """Cutting the ground into layers changes no path: the layers sum to the frozen ground to 1e-9 m."""
    for depth in (legacy_cases.UPPER, legacy_cases.LOWER):
        scene = wall_scene(depth, bank, layer_tops)
        table = scene.scan(mc.detector(depth), legacy_cases.direction_set())
        assert table.names == ("ground", *(f"layer {i + 1}" for i in range(len(layer_tops))), "niche", "air")
        total = table.length_m[:, : len(layer_tops) + 1].sum(axis=1)
        np.testing.assert_allclose(total, FROZEN[f"wall/{bank}/{depth:.0f}m"], rtol=0, atol=ATOL_PATH)
        scene.check()  # touching pieces are not overlaps for Geant4


def test_a_layer_holds_its_thickness_over_the_sine_of_elevation():
    """A layer between 3 m and 8 m of depth holds 5 m / u_z of every line of sight that leaves by the brow."""
    depth = legacy_cases.LOWER
    u = legacy_cases.direction_set()
    inward = u[:, 0] >= 0  # towards the canyon a line of sight leaves through the wall instead
    table = wall_scene(depth, "flat", PIER_LAYERS).scan(mc.detector(depth), u[inward])
    thickness = PIER_LAYERS[0] - PIER_LAYERS[1]
    np.testing.assert_allclose(table.length_in_m("layer 1"), thickness / u[inward, 2], rtol=0, atol=ATOL_PATH)
    np.testing.assert_allclose(
        table.length_in_m("ground"), -PIER_LAYERS[0] / u[inward, 2], rtol=0, atol=ATOL_PATH
    )


def test_two_banks_and_layers_against_the_analytic_path():
    """A trench between two banks: every layer has the analytic path, also where it is in two pieces.

    The trench is 6 m wide (the side of the pile cap) and 12 m deep (the upper niche level), the
    banks rise like the rising bank of legacy_cases. Layer tops at 3 m (two pieces, one a bank),
    at the trench floor itself and at 22 m. Lines of sight start in the far bank, below the
    floor of the trench, and cross it. Largest difference met: 2e-14 m.
    """
    width, floor = CAP_HI[0] - CAP_LO[0], -legacy_cases.UPPER
    rise, xr = legacy_cases.BANKS["rising"]["rise"], legacy_cases.BANKS["rising"]["xr"]
    profile = [(-400.0, rise), (-width - xr, rise), (-width, 0.0), (-width, floor)]
    profile += [(0.0, floor), (0.0, 0.0), (xr, rise), (400.0, rise)]
    tops = (-mc.REFERENCE_VOID[0], floor, -legacy_cases.LOWER)
    layers = tuple(G.Layer(f"layer {i + 1}", z, M.STANDARD_ROCK) for i, z in enumerate(tops))
    ground = G.Ground("ground", profile, -200.0, -400.0, 400.0, SOIL, layers=layers)
    assert [piece.owner for piece in ground.pieces()] == ["ground"] * 2 + ["layer 1"] * 2 + [
        "layer 2",
        "layer 3",
    ]
    scene = G.Scene([ground], altitude_m=mc.ALTITUDE)
    scene.check()
    origin = np.array([-width - mc.REFERENCE_VOID[0], 0.0, (floor - legacy_cases.LOWER) / 2])  # 17 m deep
    rng = np.random.default_rng(6)
    elevation, azimuth = np.radians(rng.uniform(10.0, 80.0, 300)), np.radians(rng.uniform(-20.0, 20.0, 300))
    u = np.stack(
        [np.cos(elevation) * np.cos(azimuth), np.cos(elevation) * np.sin(azimuth), np.sin(elevation)], axis=1
    )
    table = scene.scan(origin, u)
    bands = zip(("ground", *(layer.name for layer in layers)), (np.inf, *tops), (*tops, -200.0), strict=True)
    for name, z_hi, z_lo in bands:
        expected = [path_under_profile(profile, origin, v, z_lo, z_hi) for v in u]
        np.testing.assert_allclose(table.length_in_m(name), expected, rtol=0, atol=ATOL_PATH)
    assert table.length_in_m("layer 1").max() > 20  # both banks are in it
    assert np.all(table.length_in_m("layer 2") >= 5 - ATOL_PATH)  # 5 m of it lie above the lines of sight
    assert np.all(table.length_in_m("layer 3") == 0)  # the lines of sight start above it


# ------------------------------------------------ analytic chords of sphere, cylinder and box


@pytest.mark.parametrize("layer_tops", [(), PIER_LAYERS])
def test_chords_of_sphere_cylinder_and_box(layer_tops):
    """The pier of the spike: chords of its cap, piles and void against analytic ones to 1e-9 m.

    With layers the cap stays whole and every pile is in three parts; a part ends on the face of
    its layer. Largest differences met: 1e-14 m for the cap, 3e-12 m for the piles, 4e-13 m for
    the void.
    """
    scene, cap, piles = pier_scene(layer_tops)
    origin = SPIKE.detector(SPIKE.PIER_NICHE)
    u = slope_grid(201, 1.2).reshape(-1, 3)
    table = scene.scan(origin, u)
    box = legacy.box_chord(u, origin, np.array(cap.lo_m), np.array(cap.hi_m))
    cylinders = sum(cylinder_chord(u, origin, pile) for pile in piles)
    sphere = legacy.sphere_chord(u, origin, np.array(PIER_VOID[:3]), PIER_VOID[3])
    for name, chord in (("cap", box), ("piles", cylinders), ("void", sphere)):
        assert np.count_nonzero(chord) > 100, name
        np.testing.assert_allclose(table.length_in_m(name), chord, rtol=0, atol=ATOL_PATH, err_msg=name)
    assert scene.hosts == {"void": "layer 1" if layer_tops else "soil"}
    assert len(scene.mulder.media) == (34 if layer_tops else 14)  # the spike scene had 14 media
    # everything that is not air adds up to the path behind a wall with this niche
    solid = table.length_m[:, [i for i, b in enumerate(table.bodies) if b.kind not in ("air", "cut")]]
    floor = legacy.Niche(hn=NICHE.hn - SPIKE.EPS)  # the detector of the spike is EPS above the floor
    prelim = legacy.soil_length(legacy.Model(niche=floor), u, SPIKE.PIER_NICHE - SPIKE.EPS)
    np.testing.assert_allclose(solid.sum(axis=1), prelim, rtol=0, atol=ATOL_PATH)
    scene.check()


def test_a_structure_in_the_air_and_an_anomaly_inside_a_structure():
    """The pile cap of legacy_cases in the ground, a block of its size standing on it, a bubble in the cap.

    The block is a daughter of the world, the bubble a daughter of the cap, placed relative to
    its centre.
    """
    cap = G.Box(CAP_LO, CAP_HI)
    block = G.Box((CAP_LO[0], CAP_LO[1], 0.0), (CAP_HI[0], CAP_HI[1], CAP_HI[2] - CAP_LO[2]))
    centre, radius = (CAP_LO + CAP_HI) / 2 + np.array([1.0, 0.5, 0.0]), mc.REFERENCE_VOID[3]
    extra = [
        G.Structure("cap", cap, M.CONCRETE),
        G.Structure("block", block, M.CONCRETE),
        G.Anomaly("bubble", G.Sphere(centre, radius), drho_rel=-1.0),
    ]
    depth = legacy_cases.LOWER
    scene = wall_scene(depth, extra=extra)
    assert scene.hosts == {"bubble": "cap"}
    origin = mc.detector(depth)
    axis = (centre - origin) / np.linalg.norm(centre - origin)
    # 1152 directions over the cap and the block, 23 m away, and as many within the bubble and around it
    u = np.concatenate([mc.cone(axis, 0.3)[0], mc.cone(axis, 0.04)[0]])
    table = scene.scan(origin, u)
    bubble = legacy.sphere_chord(u, origin, centre, radius)
    whole = legacy.box_chord(u, origin, CAP_LO, CAP_HI)
    above = legacy.box_chord(u, origin, np.array(block.lo_m), np.array(block.hi_m))
    assert min(np.count_nonzero(c) for c in (bubble, whole, above)) > 20
    np.testing.assert_allclose(table.length_in_m("bubble"), bubble, rtol=0, atol=ATOL_PATH)
    np.testing.assert_allclose(table.length_in_m("cap"), whole - bubble, rtol=0, atol=ATOL_PATH)
    np.testing.assert_allclose(table.length_in_m("block"), above, rtol=0, atol=ATOL_PATH)
    assert list(scene.locate([centre, block.centre_m, origin, (5.0, 0.0, -5.0), (-0.5, 0.0, 0.0)])) == [
        "bubble",
        "block",
        "niche",
        "ground",
        "air",
    ]
    scene.check()


# ------------------------------------------------ shapes of a scan


def test_a_scan_keeps_the_shape_of_its_lines_of_sight():
    """Directions on a grid give a table on that grid; points are broadcast or given one by one."""
    depth = legacy_cases.UPPER
    scene = wall_scene(depth)
    grid = slope_grid(7)
    table = scene.scan(mc.detector(depth), grid)
    assert table.shape == (7, 7) and table.length_m.shape == (7, 7, 3)
    assert table.position_m.shape == table.direction.shape == (7, 7, 3)
    flat = scene.scan(np.broadcast_to(mc.detector(depth), (49, 3)), grid.reshape(-1, 3))
    np.testing.assert_array_equal(flat.length_m, table.length_m.reshape(49, 3))
    one = scene.scan(mc.detector(depth), [0.0, 0.0, 1.0])
    assert one.shape == () and one.length_in_m("ground") == pytest.approx(depth - NICHE.hn, abs=ATOL_PATH)
    # two detectors, one direction each
    points = np.array([mc.detector(depth), [5.0, 1.0, -mc.REFERENCE_VOID[0]]])
    both = scene.scan(points, [[0.0, 0.0, 1.0]] * 2)
    np.testing.assert_allclose(
        both.length_in_m("ground"), [depth - NICHE.hn, mc.REFERENCE_VOID[0]], atol=ATOL_PATH
    )
    with pytest.raises(ValueError, match="read-only"):
        table.length_m[0, 0, 0] = 1.0
    none = scene.scan(mc.detector(depth), np.zeros((0, 3)))  # no direction: an empty table
    assert none.shape == (0,) and none.opacity_g_cm2().shape == (0,)


# ------------------------------------------------ materials and densities given to Mulder


def test_media_carry_the_materials_and_bulk_densities_of_the_bodies():
    """A soil gives its mixture and rho_d (1 + w); concrete its own density; an anomaly scales its host."""
    scene, _, _ = pier_scene()
    by_name = {body.name: body for body in scene.table_bodies}
    assert by_name["soil"].material == register(SOIL.material)
    assert by_name["soil"].density_g_cm3 == pytest.approx(1.60 * 1.15, rel=1e-12)
    assert by_name["layer 2"].density_g_cm3 == pytest.approx(1.70 * 1.10, rel=1e-12)
    assert by_name["cap"].material == by_name["piles"].material == register(M.CONCRETE)
    assert by_name["piles"].density_g_cm3 == 2.300
    assert (by_name["void"].host, by_name["void"].material) == ("layer 1", by_name["layer 1"].material)
    assert by_name["void"].density_g_cm3 == 0.0  # drho_rel = -1
    assert by_name["niche"].material == by_name["air"].material == register(M.AIR)
    # kg/m3 for Mulder; a void is 1e-9 kg/m3 there, not zero (`test_a_void_is_empty_for_mulder_too`)
    known = {(body.material, round(max(1e3 * body.density_g_cm3, 1e-9), 9)) for body in scene.table_bodies}
    for medium in scene.mulder.media:
        assert (medium.material, round(medium.density, 9)) in known
    denser = G.Scene(
        [G.Ground("ground", wall_profile("flat"), -400.0, -1500.0, 1500.0, M.CONCRETE, density_g_cm3=2.4)],
        altitude_m=0.0,
    )
    assert denser.table_bodies[0].density_g_cm3 == 2.4  # as the preliminary calculation took concrete


def test_a_void_is_empty_for_mulder_too():
    """Flux above 1 GeV along the axis of the reference void: +8.7 % against the same scene without it.

    Mulder reads a density of exactly zero as "not set" and takes the density of the material,
    2.27 g/cm3 for the loam of the spike: the flux through a void would then be 2.1 % below that
    through soil of 1.8 g/cm3. The medium of a void therefore gets 1e-9 kg/m3, while its
    `PathTable` density stays zero. The frozen scene has air of 1.2 kg/m3 in the void and 20 m
    of sea-level air above the brow: 1.8e-4 less flux.
    """
    void = G.Anomaly("void", G.Sphere(mc.REFERENCE_VOID[:3], mc.REFERENCE_VOID[3]), drho_rel=-1.0)
    where = {"latitude_deg": mc.SITE["latitude"], "longitude_deg": mc.SITE["longitude"]}
    loam = {"material": spike_loam(), "density_g_cm3": 1e-3 * mc.SOIL}
    origin = mc.detector(mc.VOID_DEPTH)
    axis = np.array(mc.REFERENCE_VOID[:3]) - origin
    u = (axis / np.linalg.norm(axis))[None, :]

    def flux(scene):
        return mc.integral_many_above(mulder.Fluxmeter(geometry=scene.mulder), origin, u, 1.0)[0]

    scene = wall_scene(mc.VOID_DEPTH, extra=[void], **loam, **where)
    (medium,) = [m for m in scene.mulder.media if m.density < 1.0]
    assert medium.density == 1e-9 and scene.table_bodies[2].density_g_cm3 == 0.0
    through_void, through_soil = flux(scene), flux(wall_scene(mc.VOID_DEPTH, **loam, **where))
    assert through_void / through_soil - 1 == pytest.approx(0.0871, abs=0.0005)
    frozen = mc.build_scene(mc.VOID_DEPTH, sphere=mc.REFERENCE_VOID)  # 20 m of default air above the brow
    frozen.media[mc.medium_index(frozen, "Anomaly")].density = mc.AIR
    reference = mc.integral_many_above(mulder.Fluxmeter(geometry=frozen), origin, u, 1.0)[0]
    assert through_void == pytest.approx(reference, rel=3e-4)  # met: 1.8e-4
    medium.density = 0.0
    assert flux(scene) / through_soil - 1 == pytest.approx(-0.021, abs=0.002)  # what the floor is for


# ------------------------------------------------ the air of a scene


def test_scene_air_has_the_density_of_the_atmosphere_at_the_site():
    """1.0997 kg/m3 at 1150 m in the US standard atmosphere, against 1.205 of Mulder's air; a thin box."""
    scene = wall_scene(mc.VOID_DEPTH)
    assert scene.air_density_g_cm3 == pytest.approx(1.0997e-3, rel=1e-4)
    assert scene.air_density_g_cm3 == pytest.approx(
        1e-3 * FROZEN["scene_air/1150m/top1m/site"]["air_density_kg_m3"], rel=1e-12
    )
    air = [m for m in scene.mulder.media if m.material == register(M.AIR)]
    assert len(air) == 2 and all(
        m.density == 1e3 * scene.air_density_g_cm3 for m in air
    )  # the box and the niche
    assert scene.hi_m[2] == 1.0 and scene.lo_m == (-1.0, -1501.0, -401.0)  # 1 m around the bodies
    assert wall_scene(mc.VOID_DEPTH, "rising", margin_m=3.0).hi_m[2] == 1.5 + 3.0
    assert G.air_density_g_cm3(0.0) == pytest.approx(1.2311e-3, rel=1e-4)
    winter = wall_scene(mc.VOID_DEPTH, atmosphere="midlatitude-winter")
    assert winter.air_density_g_cm3 > scene.air_density_g_cm3
    for bad in ({"atmosphere": "martian"}, {"altitude_m": 1e5}, {"altitude_m": float("nan")}):
        with pytest.raises(ValueError):
            wall_scene(mc.VOID_DEPTH, **bad)


def test_open_sky_through_the_scene_air_is_that_of_the_atmosphere():
    """Above 1 GeV at 45 degrees towards the canyon: within 1e-5 of Mulder without a local geometry.

    The scene has 1 m of air above the brow at the density of the site and gives -2.5e-7 against
    the atmosphere alone (frozen `scene_air/1150m/atmosphere`). The frozen scene with the same
    box and Mulder's own air is at -2.6e-6: the air of the library has the composition of
    Mulder's, but their CSDA ranges differ by 3e-5 to 7e-5, which is +2.4e-6 here. The tolerance
    is four times that, and thirty times below what 20 m of sea-level air does.
    """
    scene = wall_scene(
        mc.VOID_DEPTH,
        material=spike_loam(),
        latitude_deg=mc.SITE["latitude"],
        longitude_deg=mc.SITE["longitude"],
    )
    elevation = np.radians(mc.AIR_ELEVATION)
    u = np.array([[-np.cos(elevation), 0.0, np.sin(elevation)]])
    meter = mulder.Fluxmeter(geometry=scene.mulder)
    flux = mc.integral_many_above(meter, mc.detector(mc.VOID_DEPTH), u, 1.0, k_max=1e5, n=801)[0]
    atmosphere = FROZEN["scene_air/1150m/atmosphere"]["open_sky/K>1GeV"]
    assert flux == pytest.approx(atmosphere, rel=1e-5)
    assert abs(flux / atmosphere - 1) < 1e-6  # met: -2.5e-7
    assert flux == pytest.approx(FROZEN["scene_air/1150m/top1m/site"]["open_sky/K>1GeV"], rel=1e-5)
    k = np.array([0.3, 1.0, 3.0, 10.0])  # kinetic energy, GeV
    own = np.asarray(mulder.Physics().compile("Air").range(k)) * mc.default_density("Air") / 10  # g/cm2
    np.testing.assert_allclose(M.csda_range(M.AIR, k), own, rtol=1e-4)
    assert abs(flux / FROZEN["scene_air/1150m/top20m/default"]["open_sky/K>1GeV"] - 1) > 2e-4


def test_what_the_air_of_a_local_geometry_costs():
    """The frozen `scene_air` cases: sea-level air and a thick box against the air of a scene.

    Relative to 1 m of air at the density of the site, above 1 GeV at 1150 m:
    20 m of Mulder's default air: void counts -0.012 %, open sky -0.030 %;
    1000 m of it: -0.91 % and -1.55 %; 1000 m at the site density: -0.31 % and -0.55 %;
    20 m at the site density: -0.0002 % and -0.0008 %. The effect of the void moves by less than
    0.03 percentage points in every case.
    """
    base = FROZEN["scene_air/1150m/top1m/site"]
    expected = {  # case: (void counts, open sky), relative to the base
        "top20m/default": (-1.24e-4, -2.97e-4),
        "top1000m/default": (-9.09e-3, -1.551e-2),
        "top1000m/site": (-3.11e-3, -5.48e-3),
        "top20m/site": (-2.0e-6, -7.7e-6),
    }
    for name, (void, sky) in expected.items():
        case = FROZEN[f"scene_air/1150m/{name}"]
        assert case["void/N0"] / base["void/N0"] - 1 == pytest.approx(void, abs=5e-6 + 5e-3 * abs(void))
        assert case["open_sky/K>1GeV"] / base["open_sky/K>1GeV"] - 1 == pytest.approx(
            sky, abs=5e-7 + 5e-3 * abs(sky)
        )
        assert abs(case["void/effect_percent"] - base["void/effect_percent"]) < 0.03
    # 20 m of default air is the scene of the frozen reference void
    assert (
        FROZEN["scene_air/1150m/top20m/default"]["void/N0"]
        == FROZEN["reference_void/1150m/cone1/K>1GeV"]["N0"]
    )
    atmosphere = FROZEN["scene_air/1150m/atmosphere"]["open_sky/K>1GeV"]
    assert base["open_sky/K>1GeV"] / atmosphere - 1 == pytest.approx(-2.6e-6, abs=1e-6)


# ------------------------------------------------ speed


@pytest.mark.validation
def test_scan_speed_on_the_pier():
    """At least 5e4 directions/s on the pier of the spike, with and without layers.

    Linux VM (4 cores, arm64): 4.0e5 directions/s with 14 media, 2.7e5 with 34; macOS on the
    same MacBook Air (Python 3.14): 5.9e4 and 5.4e4, 18 % and 8 % above the limit. The best of
    three runs after a warm-up, 40 401 directions each; the checks of `Scene.scan` are included.

    Marked `validation`: with 8 % to spare on macOS a timing does not belong in every run.
    Run it with `python -m pytest tests/unit/test_geometry_scene.py -m validation`.
    """
    origin = SPIKE.detector(SPIKE.PIER_NICHE)
    u = slope_grid(201, 1.2).reshape(-1, 3)
    for layer_tops in ((), PIER_LAYERS):
        scene, _, _ = pier_scene(layer_tops)
        scene.scan(origin, u[:1000])
        best = float("inf")
        for _ in range(3):
            start = time.perf_counter()
            scene.scan(origin, u)
            best = min(best, time.perf_counter() - start)
        assert len(u) / best > 5e4, f"{len(u) / best:.2e} directions/s with layers at {layer_tops}"


# ------------------------------------------------ the example of the note


def test_example_of_the_note():
    """The example of `docs/notes/geometry.md`: the reference void behind a wall, hypotheses on one scan."""
    soil = M.loam(1.60, 15.0)  # bulk density 1.84 g/cm3
    bank = G.Ground(
        "bank",
        profile_xz_m=[(0, -200), (0, 0), (400, 0)],
        z_bottom_m=-200,
        y_min_m=-400,
        y_max_m=400,
        material=soil,
    )
    niche = G.Cut("niche", G.Box((0, -0.15, -10.001), (0.5, 0.15, -9.7)))
    void = G.Anomaly("void", G.Sphere((3, 0, -3), 0.5), drho_rel=-1.0)
    scene = G.Scene([bank, niche, void], altitude_m=1150.0)

    detector = (0.25, 0.0, -10.0)
    towards_void = np.array([2.75, 0.0, 7.0]) / np.hypot(2.75, 7.0)
    table = scene.scan(detector, towards_void)

    assert table.names == ("bank", "niche", "void", "air")
    assert table.length_in_m("void") == pytest.approx(1.0, abs=ATOL_PATH)  # its diameter
    assert table.length_in_m("bank") == pytest.approx(9.4217, abs=1e-4)
    assert table.opacity_g_cm2() == pytest.approx(1733.6, abs=0.05)  # with the void
    assert table.opacity_g_cm2(drho_rel={"void": 0.0}) == pytest.approx(1917.6, abs=0.05)  # without it
    assert table.opacity_g_cm2({"bank": 1.70}, {"void": -0.2}) == pytest.approx(1737.7, abs=0.05)
