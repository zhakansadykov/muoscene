"""`transport.table`: the equivalent opacity of a path and the domain of a flux table.

Roadmap 0.1-08. The numbers that depend on Mulder are the frozen cases `threshold_slab/…` and
`scene_air/1150m/atmosphere` of `tests/regression/mulder_cases.py`: the uniform slab and the
open sky the roadmap asks the fast path to match. Materials are the loam and the concrete of
those cases; depths and zenith angles are theirs.
"""

import functools
import json
import math
import pathlib

import numpy as np
import pytest

from muoscene import materials as M
from muoscene.flux import LEGACY_RESPONSE, EnergyResponse, LegacyTable, MulderTable
from muoscene.geometry import air_density_g_cm3
from muoscene.materials.composition import Material, Mixture
from muoscene.transport import (
    BELOW_CUTOFF,
    DEFINED,
    LOW_ELEVATION,
    TOO_OPAQUE,
    equivalent_opacity,
    flux_domain,
)
from tests.regression import mulder_cases as mc

FROZEN = json.loads(
    (pathlib.Path(mc.__file__).parent / "data" / "mulder_reference.json").read_text(encoding="utf-8")
)["cases"]


def spike_loam():
    """The loam of the spike: the `RegLoam` of the frozen cases as a muoscene material."""
    quartz = Material.from_formula("spike quartz", "SiO2", 2.65, i_ev="mulder")
    calcite = Material.from_formula("spike calcite", "CaCO3", 2.71, i_ev="mulder")
    alumina = Material.from_formula("spike alumina", "Al2O3", 3.95, i_ev="mulder")
    water = Material.from_formula("spike water", "H2O", 1.02, i_ev="mulder")  # Mulder's Water
    return Mixture.of("spike loam", {quartz: 0.60, calcite: 0.15, alumina: 0.12, water: 0.13})


LOAM = spike_loam()
CONCRETE = Material.from_elements("spike concrete", mc.CONCRETE, mc.CONCRETE_DENSITY / 1000.0, i_ev="mulder")
NAMES = {M.register(LOAM): LOAM, M.register(CONCRETE): CONCRETE}
LOAM_NAME, CONCRETE_NAME = NAMES
AIR = air_density_g_cm3(mc.ALTITUDE)  # 1.0997e-3 g/cm³


@functools.cache
def table(altitude_m, material=LOAM):
    """Flux table of a detector at an altitude; not cached on disk."""
    return MulderTable(altitude_m, material, cache=False)


# ------------------------------------------------------------------ weights


@pytest.mark.parametrize("k_min", [0.2, 1.0])
def test_one_material_is_exact(k_min):
    """A path through concrete alone: X_eq = R_ref(K_top) − R_ref(K_min) with R_c(K_top) = R_c(K_min) + X.

    Tolerance 3e-5: the range of Mulder and its inverse round-trip to 6e-6 (`docs/notes/materials.md`).
    """
    x = np.array([30.0, 345.0, 2000.0, 2.0e4])
    got = equivalent_opacity(table(mc.ALTITUDE), {CONCRETE_NAME: x}, NAMES, k_min)
    top = M.energy_to_cross(CONCRETE, np.asarray(M.csda_range(CONCRETE, k_min)) + x)
    exact = np.asarray(M.csda_range(LOAM, top)) - M.csda_range(LOAM, k_min)
    assert got.matter_g_cm2 == pytest.approx(exact, rel=3e-5)
    assert np.array_equal(got.lookup_g_cm2, got.matter_g_cm2)  # no air asked for
    assert 0.99 < got.matter_g_cm2[2] / x[2] < 1.0  # concrete weighs 0.995 against this loam at 7 GeV
    assert not got.matter_g_cm2.flags.writeable


def test_without_a_threshold_it_is_the_weight_of_the_conventions():
    """At the first energy of the table, 1 MeV, X_eq = R_ref(R_m⁻¹(X)) of `materials.equivalent_opacity`."""
    x = np.array([345.0, 2000.0, 2.0e4])
    sky = table(mc.ALTITUDE)
    got = equivalent_opacity(sky, {CONCRETE_NAME: x}, NAMES, sky.range_energy_gev[0])
    assert got.matter_g_cm2 == pytest.approx(M.equivalent_opacity(CONCRETE, LOAM, x), rel=1e-4)


def test_the_reference_material_weighs_one():
    x = np.array([0.0, 12.5, 3680.0])
    got = equivalent_opacity(table(mc.ALTITUDE), {LOAM_NAME: x}, NAMES, 1.0)
    assert np.array_equal(got.matter_g_cm2, x)
    assert got.iterations == 0


def test_the_choice_of_the_reference_material_does_not_matter():
    """Loam with concrete in it, through a table of either material: the same flux.

    Met 2e-7 above 1 GeV; a weight taken at one energy would differ between the two tables.
    """
    loam, concrete = np.array([3000.0, 3600.0, 5400.0]), np.array([345.0, 0.0, 1840.0])
    sine = np.array([1.0, 0.8, 0.75])
    response = EnergyResponse.threshold(1.0)
    flux = []
    for material in (LOAM, CONCRETE):
        sky = table(mc.ALTITUDE - 20.0, material)
        x = equivalent_opacity(sky, {LOAM_NAME: loam, CONCRETE_NAME: concrete}, NAMES, 1.0)
        flux.append(sky.integral(response, sine, x.matter_g_cm2))
    assert flux[0] == pytest.approx(flux[1], rel=1e-4)


def test_displaced_air():
    """The air a body displaces leaves the opacity of the table: ρ_air L, weighted like a material."""
    sky = table(mc.ALTITUDE)
    x, length = np.array([3680.0, 3680.0, 0.01]), np.array([20.0, 0.0, 30.0])
    got = equivalent_opacity(sky, {LOAM_NAME: x}, NAMES, 1.0, solid_length_m=length, air_density_g_cm3=AIR)
    weight = got.displaced_air_g_cm2[0] / (AIR * 20.0 * 100.0)
    assert weight == pytest.approx(1.155, abs=5e-3)  # a gas has no density effect
    assert got.lookup_g_cm2[0] == pytest.approx(3680.0 - got.displaced_air_g_cm2[0], rel=1e-14)
    assert got.lookup_g_cm2[1] == 3680.0
    assert got.lookup_g_cm2[2] == 0.0  # 3.3 g/cm² of air against 0.01 of matter: not below zero


def test_what_the_equivalent_opacity_refuses():
    sky = table(mc.ALTITUDE)
    with pytest.raises(TypeError, match="FluxTable"):
        equivalent_opacity(LegacyTable(), {LOAM_NAME: 1.0}, NAMES, 1.0)
    with pytest.raises(ValueError, match="no material"):
        equivalent_opacity(sky, {"ms-unknown": 1.0}, NAMES, 1.0)
    with pytest.raises(ValueError, match="non-negative"):
        equivalent_opacity(sky, {CONCRETE_NAME: -1.0}, NAMES, 1.0)
    with pytest.raises(ValueError, match="range table"):
        equivalent_opacity(sky, {CONCRETE_NAME: 1.0}, NAMES, 1.0e-4)
    with pytest.raises(ValueError, match="both"):
        equivalent_opacity(sky, {CONCRETE_NAME: 1.0}, NAMES, 1.0, solid_length_m=1.0)


def test_the_alternatives_that_were_not_taken():
    """What the other ways of weighting cost, on 3000 g/cm² of loam under 1840 g/cm² of concrete.

    The judge is the chain of ranges in the true order — loam at the detector, concrete above —
    which is continuous transport itself; flux above the threshold at a sine of elevation of
    0.75, relative to it. The iteration taken: +0.023 % above 1 GeV, +0.034 % above 0.2 GeV.
    The same chain in the wrong order: +0.048 and +0.075 %. A weight at one energy: +0.17,
    +0.06, −0.08 % at 1, 7, 100 GeV. No weights at all: −0.17 %. The numbers hold the decision
    of 21 November 2026 (`docs/notes/forward.md`).
    """
    sky = table(mc.ALTITUDE - 20.0)
    loam, concrete, sine = 3000.0, 1840.0, 0.75

    def chain(k_min, order):
        k = k_min
        for material, x in order:
            k = M.energy_to_cross(material, M.csda_range(material, k) + x)
        return M.csda_range(LOAM, k) - M.csda_range(LOAM, k_min)

    for k_min, ours, wrong, fixed, plain in (
        (1.0, 2.3e-4, 4.8e-4, (1.75e-3, 5.9e-4, -8.0e-4), -1.72e-3),
        (0.2, 3.4e-4, 7.5e-4, (1.80e-3, 5.8e-4, -8.8e-4), -1.86e-3),
    ):
        response = EnergyResponse.threshold(k_min)
        true = sky.integral(response, sine, chain(k_min, [(LOAM, loam), (CONCRETE, concrete)]))

        def off(x, true=true, response=response):
            return sky.integral(response, sine, x) / true - 1

        x = {LOAM_NAME: np.array([loam]), CONCRETE_NAME: np.array([concrete])}
        assert off(equivalent_opacity(sky, x, NAMES, k_min).matter_g_cm2[0]) == pytest.approx(ours, abs=2e-5)
        assert off(chain(k_min, [(CONCRETE, concrete), (LOAM, loam)])) == pytest.approx(wrong, abs=2e-5)
        for energy, expected in zip((1.0, 7.0, 100.0), fixed, strict=True):
            weighted = loam + M.weight(CONCRETE, LOAM, energy) * concrete
            assert off(weighted) == pytest.approx(expected, abs=3e-5)
        assert off(loam + concrete) == pytest.approx(plain, abs=3e-5)


def test_a_table_at_the_altitude_of_the_surface_is_wrong_for_the_open_sky():
    """Seen from a niche 20 m below the brow the open sky above 1 GeV would be 0.31 to 0.39 % high.

    Sines of elevation of 1, 0.71 and 0.5. Under soil such a table is 0.03–0.05 % high
    (`docs/notes/flux.md`); one table at the altitude of the detector serves both sides.
    """
    response = EnergyResponse.threshold(1.0)
    sine = np.array([1.0, math.sqrt(0.5), 0.5])
    high = table(mc.ALTITUDE).integral(response, sine, 0.0) / table(mc.ALTITUDE - 20.0).integral(
        response, sine, 0.0
    )
    assert high - 1 == pytest.approx([3.10e-3, 3.45e-3, 3.90e-3], abs=3e-5)


# ------------------------------------------------------------------ against Mulder


@pytest.mark.parametrize("depth", mc.THRESHOLD_DEPTHS)
def test_a_uniform_slab_matches_mulder(depth):
    """The frozen `threshold_slab` cases: Mulder with the real soil above a detector at 1150 m − depth.

    A table at the altitude of the detector, entered with the soil less the air it displaces.
    Tolerance 3e-4: what is left is the decay inside the soil; met −0.003 to +0.021 % over
    0.5, 2 and 20 m, zenith angles of 0 and 45 degrees and both thresholds. Without the air term
    the same table is −0.074 to −0.078 % under 20 m (`docs/notes/flux.md`).
    """
    sky = table(mc.ALTITUDE - depth)
    for zenith in mc.THIN_ZENITH:
        sine = math.cos(math.radians(zenith))
        length = depth / sine
        x = {LOAM_NAME: np.array([mc.SOIL / 1000.0 * length * 100.0])}
        got = equivalent_opacity(sky, x, NAMES, 1.0, solid_length_m=length, air_density_g_cm3=AIR)
        frozen = FROZEN[f"threshold_slab/{mc.ALTITUDE:.0f}m/{depth:g}m/zenith{zenith:g}"]
        for k_min in mc.THRESHOLDS:
            response = EnergyResponse.threshold(k_min)
            assert flux_domain(sky, response, sine, got.lookup_g_cm2)[0] == DEFINED
            with_air = sky.integral(response, sine, got.lookup_g_cm2)[0] / frozen[f"K>{k_min:g}GeV"] - 1
            without = sky.integral(response, sine, got.matter_g_cm2)[0] / frozen[f"K>{k_min:g}GeV"] - 1
            assert abs(with_air) < 3e-4
            assert without < with_air  # the table alone is the lower one
            if depth == 20.0:
                assert -9e-4 < without < -6e-4 and abs(with_air) < abs(without) / 2


def test_the_open_sky_matches_mulder():
    """The frozen `scene_air/1150m/atmosphere`: Mulder without a geometry, above 1 GeV at 45 degrees."""
    sky = table(mc.ALTITUDE - mc.VOID_DEPTH)
    response, sine = EnergyResponse.threshold(1.0), math.sin(math.radians(mc.AIR_ELEVATION))
    assert flux_domain(sky, response, sine, 0.0) == DEFINED
    frozen = FROZEN[f"scene_air/{mc.ALTITUDE:.0f}m/atmosphere"]["open_sky/K>1GeV"]
    assert sky.integral(response, sine, 0.0) == pytest.approx(frozen, rel=1e-4)


# ------------------------------------------------------------------ domain


def test_where_a_table_has_a_flux_at_1150_m():
    """Sea-level reference at 1150 m: 0.2 GeV holds only under soil, 1 GeV from 15.7 degrees up."""
    sky = table(mc.ALTITUDE)
    soft, hard = EnergyResponse.threshold(0.2), EnergyResponse.threshold(1.0)
    sine = np.sin(np.radians([90.0, 45.0, 16.0, 15.2, 13.0]))
    assert flux_domain(sky, soft, sine, 0.0).tolist() == [BELOW_CUTOFF] * 4 + [LOW_ELEVATION]
    assert flux_domain(sky, hard, sine, 0.0).tolist() == [DEFINED] * 3 + [BELOW_CUTOFF, LOW_ELEVATION]
    # no cut-off is left under 138 g/cm² by the vertical; under 100 g/cm² 0.2 GeV is still above what is left
    assert flux_domain(sky, soft, 1.0, [50.0, 100.0, 140.0, 3.0e4, 3.1e4]).tolist() == [
        BELOW_CUTOFF,
        DEFINED,
        DEFINED,
        DEFINED,
        TOO_OPAQUE,
    ]
    assert flux_domain(sky, soft, sine, 600.0).tolist() == [DEFINED] * 4 + [LOW_ELEVATION]
    assert flux_domain(sky, hard, 1.0, 0.0).dtype == np.int8


def test_the_domain_is_where_the_table_answers():
    """`table.integral` succeeds on the lines that are DEFINED and on no others: 400 random lines."""
    sky = table(mc.ALTITUDE)
    rng = np.random.default_rng(8)
    sine = rng.uniform(0.2, 1.0, 400)
    opacity = np.where(rng.uniform(size=400) < 0.3, 0.0, np.exp(rng.uniform(0.0, math.log(4.0e4), 400)))
    for response in (EnergyResponse.threshold(0.2), EnergyResponse.threshold(1.0)):
        code = flux_domain(sky, response, sine, opacity)
        assert {DEFINED, BELOW_CUTOFF, LOW_ELEVATION, TOO_OPAQUE} == set(code.tolist())
        assert np.all(sky.integral(response, sine[code == DEFINED], opacity[code == DEFINED]) > 0)
        for i in np.flatnonzero(code != DEFINED):
            with pytest.raises(ValueError):
                sky.integral(response, sine[i], opacity[i])


def test_the_domain_of_the_preliminary_calculation():
    legacy = LegacyTable()
    code = flux_domain(legacy, LEGACY_RESPONSE, [1.0, 0.25, 0.2, 1.0], [0.0, 3.0e4, 0.0, 3.1e4])
    assert code.tolist() == [DEFINED, DEFINED, LOW_ELEVATION, TOO_OPAQUE]


def test_what_the_domain_refuses():
    sky = table(mc.ALTITUDE)
    with pytest.raises(TypeError, match="EnergyResponse"):
        flux_domain(sky, LEGACY_RESPONSE, 1.0, 0.0)
    with pytest.raises(TypeError, match="LEGACY_RESPONSE"):
        flux_domain(LegacyTable(), EnergyResponse.threshold(1.0), 1.0, 0.0)
    with pytest.raises(TypeError, match="FluxTable"):
        flux_domain(object(), EnergyResponse.threshold(1.0), 1.0, 0.0)
    with pytest.raises(ValueError, match="non-negative"):
        flux_domain(sky, EnergyResponse.threshold(1.0), 1.0, -1.0)
    with pytest.raises(ValueError, match="highest energy"):
        flux_domain(sky, EnergyResponse((1.0, 2.0e4), (0.5, 1.0)), 1.0, 0.0)
