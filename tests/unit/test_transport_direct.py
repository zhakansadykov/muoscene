"""`transport.direct` and `Site.differential_flux`: the flux of Mulder through a scene, in ENU.

Roadmap 0.1-08. Sizes are those the repository already has: the bank and the niche of
`tests/regression/legacy_cases.py`, the loam and the site of `mulder_cases.py`. The direct flux
is slow, so the file asks Mulder for fewer than a hundred lines of sight.
"""

import functools
import json
import math
import pathlib

import numpy as np
import pytest

from muoscene import geometry as G
from muoscene.flux import LEGACY_RESPONSE, EnergyResponse, MulderTable
from muoscene.materials.composition import Material, Mixture
from muoscene.site import Site
from muoscene.transport import DirectCheck, compare_with_direct, direct_flux, energy_nodes, seed_of
from tests.regression import legacy
from tests.regression import mulder_cases as mc

FROZEN = json.loads(
    (pathlib.Path(mc.__file__).parent / "data" / "mulder_reference.json").read_text(encoding="utf-8")
)["cases"]
NICHE = legacy.Niche()
DEPTH = mc.VOID_DEPTH  # the niche of the reference void, 10 m below the brow
HARD = EnergyResponse.threshold(1.0)


def spike_loam():
    """The loam of the spike: the `RegLoam` of the frozen cases as a muoscene material."""
    quartz = Material.from_formula("spike quartz", "SiO2", 2.65, i_ev="mulder")
    calcite = Material.from_formula("spike calcite", "CaCO3", 2.71, i_ev="mulder")
    alumina = Material.from_formula("spike alumina", "Al2O3", 3.95, i_ev="mulder")
    water = Material.from_formula("spike water", "H2O", 1.02, i_ev="mulder")  # Mulder's Water
    return Mixture.of("spike loam", {quartz: 0.60, calcite: 0.15, alumina: 0.12, water: 0.13})


LOAM = spike_loam()


@functools.cache
def site(azimuth=0.0):
    return Site(mc.SITE["latitude"], mc.SITE["longitude"], mc.ALTITUDE, azimuth)


@functools.cache
def bank(azimuth=0.0):
    """The flat bank of legacy_cases at 1150 m with the niche of the reference void; the wall is x = 0."""
    ground = G.Ground(
        "ground",
        [(0.0, -mc.GROUND[2]), (0.0, 0.0), (mc.GROUND[0], 0.0)],
        z_bottom_m=-mc.GROUND[2],
        y_min_m=-mc.GROUND[1] / 2,
        y_max_m=mc.GROUND[1] / 2,
        material=LOAM,
        density_g_cm3=mc.SOIL / 1000.0,
    )
    niche = G.Cut(
        "niche", G.Box((0.0, -NICHE.w / 2, -DEPTH - mc.EPS), (NICHE.xn, NICHE.w / 2, -DEPTH + NICHE.hn))
    )
    return site(azimuth).scene([ground, niche])


DETECTOR = (NICHE.xd, 0.0, -DEPTH)  # in the axes of the scene


def unit(azimuth_deg, elevation_deg):
    """Unit vector in ENU at an azimuth from north, clockwise, and an elevation."""
    a, e = math.radians(azimuth_deg), math.radians(elevation_deg)
    return np.array([math.cos(e) * math.sin(a), math.cos(e) * math.cos(a), math.sin(e)])


# ------------------------------------------------------------------ energies and seeds


def test_energy_nodes_start_at_the_threshold():
    """240 nodes a decade over the first half decade from the threshold, 30 beyond: 226 from 1 GeV."""
    k = energy_nodes(HARD)
    assert k[0] == 1.0 and k[-1] == 1.0e4 and k.size == 226
    assert k[120] == pytest.approx(10**0.5, rel=1e-14)
    assert np.diff(np.log(k[:121])) == pytest.approx(math.log(10) / 240)
    assert np.diff(np.log(k[120:])) == pytest.approx(math.log(10) / 30)
    assert energy_nodes(HARD, first_nodes_per_decade=30).size == 121  # the nodes before 22 November 2026
    soft = EnergyResponse((0.2, 0.5, 2.0), (0.3, 0.8, 1.0))
    k = energy_nodes(soft)
    assert k[0] == 0.2 and {0.5, 2.0} <= set(k.tolist()) and np.all(np.diff(k) > 0)
    assert energy_nodes(HARD, nodes_per_decade=240).size == 961
    with pytest.raises(TypeError, match="EnergyResponse"):
        energy_nodes(LEGACY_RESPONSE)
    with pytest.raises(ValueError, match="k_max_gev"):
        energy_nodes(HARD, k_max_gev=0.5)
    with pytest.raises(ValueError, match="nodes_per_decade"):
        energy_nodes(HARD, nodes_per_decade=2)
    with pytest.raises(ValueError, match="first_nodes_per_decade"):
        energy_nodes(HARD, first_nodes_per_decade=2)


def test_the_seed_of_a_generator_is_what_it_was_made_with():
    assert seed_of(np.random.default_rng(5)) == 5
    assert seed_of(np.random.default_rng(2**70 + 1)) == 2**70 + 1
    assert seed_of(np.random.default_rng()) >= 0  # the entropy numpy drew: the sample can be drawn again
    with pytest.raises(TypeError, match="Generator"):
        seed_of(np.random.RandomState(5))
    with pytest.raises(ValueError, match="spawned"):
        seed_of(np.random.default_rng(np.random.SeedSequence(5).spawn(1)[0]))
    with pytest.raises(ValueError, match="does not know its seed"):
        seed_of(np.random.default_rng([1, 2]))


# ------------------------------------------------------------------ Mulder through the scene


def test_the_open_sky_through_a_scene_is_that_of_mulder():
    """The frozen `scene_air/1150m/atmosphere`: towards the canyon at 45 degrees, above 1 GeV.

    Tolerance 2e-4: 121 energy nodes here against 801 there and the air of the scene; met 5e-5.
    """
    frozen = FROZEN[f"scene_air/{mc.ALTITUDE:.0f}m/atmosphere"]["open_sky/K>1GeV"]
    flux = direct_flux(site(), bank(), DETECTOR, unit(270.0, mc.AIR_ELEVATION), HARD)
    assert flux.shape == ()
    assert flux == pytest.approx(frozen, rel=2e-4)


def test_points_and_directions_are_in_enu():
    """A site whose scene is turned by 90 degrees: the massif, +x of the scene, lies to the south.

    The flux along a direction in ENU is that of the scene which is not turned along the same
    direction in its own axes, within 1e-5: without a field the azimuth of the frame reaches
    Mulder only through the shape of the Earth, 1.8e-6. The massif is where the site says: a third
    of the open sky at 45 degrees, under 14 m of soil.
    """
    turned, plain = site(90.0), site(0.0)
    k = energy_nodes(HARD)[::10]
    position = turned.enu_from_scene(DETECTOR)
    south, north, east = unit(180.0, 45.0), unit(0.0, 45.0), unit(90.0, 45.0)
    flux = turned.differential_flux(bank(90.0), position, np.array([south, north, east]), k)
    assert flux.shape == (3, k.size)
    own = plain.differential_flux(
        bank(0.0), DETECTOR, turned.scene_from_enu(np.array([south, north, east])), k
    )
    assert flux == pytest.approx(own, rel=1e-5)
    into, towards = np.trapezoid(flux[0] * k, np.log(k)), np.trapezoid(flux[1] * k, np.log(k))
    assert 0.30 < into / towards < 0.37  # south is under the massif, north is the canyon: 0.335
    assert np.all(flux[2] < flux[1])  # along the wall: through the roof of the niche


def test_the_integral_starts_at_the_threshold():
    """Trapezoids from the first energy of the response: the default nodes against 240 a decade within 2e-4.

    Met 7e-5 above 1 GeV and 3e-5 above 0.2 GeV under the massif with 30 nodes a decade everywhere;
    2.3e-5 and 4.8e-5 with the default since 22 November 2026, 240 a decade over the first half
    decade. A response that rises from
    0.2 to 1 GeV gives a flux between the two thresholds.
    """
    into = unit(90.0, 60.0)
    values = {}
    for k_min in (0.2, 1.0):
        response = EnergyResponse.threshold(k_min)
        coarse = direct_flux(site(), bank(), DETECTOR, into, response)
        fine = direct_flux(site(), bank(), DETECTOR, into, response, nodes_per_decade=240)
        assert coarse == pytest.approx(fine, rel=2e-4)
        values[k_min] = float(coarse)
    rising = direct_flux(site(), bank(), DETECTOR, into, EnergyResponse((0.2, 1.0), (0.0, 1.0)))
    assert values[1.0] < rising < values[0.2]


def test_the_check_at_the_edge_of_the_cut_off():
    """Open sky above 1 GeV just above the elevation where 1 GeV meets the cut-off at 1150 m.

    A film 0.5 m above a flat ground, lines towards the south at 15.5, 16 and 17 degrees. The
    table against the direct flux: with 30 nodes a decade from the threshold the direct flux was
    0.62, 0.14 and 0.05 % above the table and refused a right run at 15.5 degrees (found with
    `examples/open_sky`, 22 November 2026); with 240 nodes a decade over the first half decade it
    is within 1e-3, met 0.02 %. With 2000 nodes a decade the table is 0.010 % low at 15.5 degrees.
    """
    plain = site()
    flat = G.Ground("ground", [(-100.0, 0.0), (100.0, 0.0)], -20.0, -100.0, 100.0, LOAM, mc.SOIL / 1000.0)
    scene = plain.scene([flat])
    table = MulderTable(mc.ALTITUDE + 0.5, LOAM, cache=False)
    for elevation, before in ((15.5, 6.0e-3), (16.0, 1.3e-3), (17.0, 4.0e-4)):
        u = unit(180.0, elevation)
        fast = float(table.integral(HARD, u[2], 0.0))
        direct = float(direct_flux(plain, scene, (0.0, 0.0, 0.5), u, HARD))
        assert fast == pytest.approx(direct, rel=1e-3)
        k = energy_nodes(HARD, first_nodes_per_decade=30)
        coarse = np.trapezoid(plain.differential_flux(scene, (0.0, 0.0, 0.5), u, k) * k, np.log(k))
        assert coarse / fast - 1 > before  # what the check had before


def test_what_mulder_would_answer_with_a_crash_is_refused():
    k = np.array([1.0, 10.0])
    up = unit(0.0, 90.0)
    with pytest.raises(ValueError, match="outside the box"):
        site().differential_flux(bank(), (0.0, 0.0, 50.0), up, k)
    with pytest.raises(ValueError, match="unit vectors"):
        site().differential_flux(bank(), DETECTOR, 2.0 * up, k)
    with pytest.raises(ValueError, match="towards the sky"):
        site().differential_flux(bank(), DETECTOR, -up, k)
    with pytest.raises(ValueError, match="another site"):
        site(90.0).differential_flux(bank(), DETECTOR, up, k)
    for energy in (1.0, [], [1.0, -1.0], [[1.0, 2.0]], [math.nan]):
        with pytest.raises(ValueError, match="kinetic_energy_gev"):
            site().differential_flux(bank(), DETECTOR, up, energy)


# ------------------------------------------------------------------ the comparison


def test_a_comparison_stores_the_deviation_and_refuses_beyond_its_limit():
    lines = np.array([unit(90.0, 60.0), unit(90.0, 80.0), unit(270.0, 45.0)])
    direct = direct_flux(site(), bank(), DETECTOR, lines, HARD)
    check = compare_with_direct(site(), bank(), DETECTOR, lines, HARD, direct * [1.001, 0.998, 1.0])
    assert isinstance(check, DirectCheck) and check.n_directions == 3
    assert check.deviation == pytest.approx([1e-3, -2e-3, 0.0], abs=1e-12)
    assert check.max == pytest.approx(2e-3) and check.worst == 1
    assert check.mean == pytest.approx(-1e-3 / 3) and check.rms == pytest.approx(math.sqrt(5e-6 / 3))
    assert check.limit == 5e-3 and not check.direct_m2_s_sr.flags.writeable
    with pytest.raises(RuntimeError, match="0.60 %"):
        compare_with_direct(site(), bank(), DETECTOR, lines, HARD, direct * [1.0, 1.006, 1.0])
    loose = compare_with_direct(site(), bank(), DETECTOR, lines, HARD, direct * 1.01, limit=None)
    assert loose.max == pytest.approx(0.01) and loose.limit is None
    with pytest.raises(ValueError, match="one flux"):
        compare_with_direct(site(), bank(), DETECTOR, lines, HARD, direct[:2])
    with pytest.raises(ValueError, match="limit"):
        compare_with_direct(site(), bank(), DETECTOR, lines, HARD, direct, limit=0.0)
    with pytest.raises(RuntimeError, match="nothing to compare"):
        compare_with_direct(site(), bank(), DETECTOR, lines, HARD, [1.0, math.nan, 1.0], limit=None)
