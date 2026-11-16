"""Compositions, mean excitation energies and soils (no physics tables are compiled here)."""

import math

import mulder.materials as mm
import pytest

from muoscene import materials as M
from muoscene.materials.backend import mulder_name
from muoscene.materials.composition import (
    Material,
    Mixture,
    Soil,
    formula_mass_fractions,
    oxide_mass_fractions,
    parse_formula,
)

# ------------------------------------------------ formulas


@pytest.mark.parametrize(
    ("formula", "atoms"),
    [
        ("H2O", {"H": 2, "O": 1}),
        ("SiO2", {"Si": 1, "O": 2}),
        ("CaCO3", {"Ca": 1, "C": 1, "O": 3}),
        ("Al2Si2O5(OH)4", {"Al": 2, "Si": 2, "O": 9, "H": 4}),
        ("Ca(Mg(CO3))2", {"Ca": 1, "Mg": 2, "C": 2, "O": 6}),
        ("Na0.5Ca0.5Al1.5Si2.5O8", {"Na": 0.5, "Ca": 0.5, "Al": 1.5, "Si": 2.5, "O": 8}),
    ],
)
def test_parse_formula(formula, atoms):
    assert parse_formula(formula) == pytest.approx(atoms)


@pytest.mark.parametrize("bad", ["", "h2o", "H2O)", "(H2O", "2H", "Xx2O", "H2O!"])
def test_bad_formulas_are_refused(bad):
    with pytest.raises(ValueError):
        parse_formula(bad)


def test_water_fractions_are_those_of_mulder():
    ours = formula_mass_fractions("H2O")
    theirs = dict(mm.Material("Water").composition)
    assert ours == pytest.approx(theirs, rel=1e-14)


def test_oxides():
    fractions = oxide_mass_fractions({"SiO2": 50.0, "CaO": 25.0, "CO2": 25.0})
    assert sum(fractions.values()) == pytest.approx(1.0)
    assert fractions["Si"] == pytest.approx(0.5 * 28.0855 / (28.0855 + 2 * 15.999))
    with pytest.raises(ValueError, match="loss on ignition"):
        oxide_mass_fractions({"SiO2": 90.0, "LOI": 10.0})
    with pytest.raises(ValueError):
        oxide_mass_fractions({"SiO2": 90.0, "CaO": -1.0})


# ------------------------------------------------ Z/A and I against PDG and ICRU 37


@pytest.mark.parametrize(
    ("material", "pdg"),
    [
        (M.WATER, 0.55509),
        (M.CONCRETE, 0.50274),
        (M.QUARTZ, 0.49930),
        (M.CALCITE, 0.49955),
        (M.AIR, 0.49919),
        (M.STANDARD_ROCK, 0.50000),
    ],
)
def test_z_over_a_as_in_pdg(material, pdg):
    assert material.z_over_a == pytest.approx(pdg, abs=3e-5)


@pytest.mark.parametrize(
    ("formula", "tabulated", "tolerance"),
    [
        ("H2O", 75.0, 0.5),  # ICRU 37, the value the rule was built on
        ("CaCO3", 136.4, 0.3),  # PDG
        ("SiO2", 139.2, 2.5),  # PDG; the rule gives 141.0
        ("Al2O3", 145.2, 3.0),  # PDG; the rule gives 142.6
    ],
)
def test_berger_seltzer_rule(formula, tabulated, tolerance):
    assert M.berger_seltzer_mean_excitation(formula_mass_fractions(formula)) == pytest.approx(
        tabulated, abs=tolerance
    )


def test_mulder_default_i_of_water_is_too_low():
    """Plain Bragg additivity of elemental values, what PUMAS does when I is not given: 69 eV."""
    water = Material.from_formula("water", "H2O", 1.0, i_ev="mulder")
    assert water.i_ev is None
    assert water.mean_excitation_ev == pytest.approx(69.0, abs=0.2)


def test_skeletons_from_their_sources():
    for skeleton in (M.SILICATE_SKELETON, M.LOESS_SKELETON):
        assert sum(skeleton.elements.values()) == pytest.approx(1.0)
        assert 0.495 < skeleton.z_over_a < 0.497
        assert 140 < skeleton.i_ev < 155  # iron and titanium raise it above quartz
        assert skeleton.source
    assert M.LOESS_SKELETON.elements["Ca"] > 2 * M.SILICATE_SKELETON.elements["Ca"]  # carbonate


def test_every_library_material_has_a_source_and_an_explicit_i():
    for material in M.library.MATERIALS.values():
        assert material.source and material.i_ev is not None


def test_bad_materials_are_refused():
    with pytest.raises(ValueError):
        Material.from_elements("x", {"Si": 1.0}, 0.0)
    with pytest.raises(ValueError):
        Material.from_elements("x", {}, 1.0)
    with pytest.raises(ValueError):
        Material.from_elements("x", {"Si": 1.0}, 1.0, i_ev="guess")
    with pytest.raises(ValueError, match="unknown element"):
        Material.from_elements("x", {"Qq": 1.0}, 1.0)
    with pytest.raises(ValueError, match="listed twice"):
        Material("x", (("Si", 0.5), ("Si", 0.2), ("O", 0.5)), 2.0, 100.0)


# ------------------------------------------------ mixtures and soils


def test_mixture_density_has_no_pores():
    mix = Mixture.of("half", {M.QUARTZ: 0.5, M.WATER: 0.5})
    assert mix.density_g_cm3 == pytest.approx(1 / (0.5 / 2.65 + 0.5 / 1.0))
    assert sum(mix.elements.values()) == pytest.approx(1.0)
    assert Mixture.of("same", {M.QUARTZ: 1.0, M.WATER: 1.0}).parts == mix.parts
    with pytest.raises(TypeError):
        Mixture.of("nested", {mix: 1.0})


def test_soil_quantities():
    s = M.loam(1.60, 15.0)
    assert s.bulk_density_g_cm3 == pytest.approx(1.84)
    assert s.porosity == pytest.approx(1 - 1.60 / 2.70)
    assert s.volumetric_water_content == pytest.approx(0.24)
    assert s.saturation == pytest.approx(0.24 / (1 - 1.60 / 2.70))
    assert s.water_mass_fraction == pytest.approx(0.15 / 1.15)
    parts = dict(s.material.parts)
    assert parts[M.WATER] == pytest.approx(0.15 / 1.15)
    assert parts[M.SILICATE_SKELETON] == pytest.approx(1 / 1.15)
    assert M.loam(1.60, 0.0).material is M.SILICATE_SKELETON


SATURATED_W = {dry: 100 * (1 - dry / 2.70) / dry for dry in (1.2, 1.45, 1.6, 1.85, 2.1)}  # %, water 1 g/cm3


@pytest.mark.parametrize(
    ("dry", "w"),
    [(dry, w) for dry, limit in SATURATED_W.items() for w in (0.0, 3.0, 12.5, 20.0, limit) if w <= limit],
)
def test_moisture_round_trip(dry, w):
    """Dry density and water content -> bulk density and water mass fraction -> back."""
    s = M.loam(dry, w)
    back = Soil.from_bulk(s.skeleton, s.bulk_density_g_cm3, s.water_mass_fraction, s.water)
    assert back.dry_density_g_cm3 == pytest.approx(dry, rel=1e-13)
    assert back.water_content_percent == pytest.approx(w, rel=1e-13, abs=1e-12)
    assert mulder_name(back.material) == mulder_name(s.material)  # no new Mulder material


def test_impossible_soils_are_refused():
    with pytest.raises(ValueError, match="exceeds saturation"):
        M.loam(1.50, 30.0)  # saturation at 29.6 %
    M.loam(1.50, 29.6)
    with pytest.raises(ValueError, match="particle density"):
        M.sand(2.65, 0.0)
    with pytest.raises(ValueError):
        M.loam(1.6, -1.0)


# ------------------------------------------------ names in Mulder


def test_names_follow_the_definition():
    a = Material.from_formula("Quartz sand!", "SiO2", 2.65, i_ev=139.2)
    b = Material.from_formula("Quartz sand!", "SiO2", 2.65, i_ev=139.2)
    assert mulder_name(a) == mulder_name(b)
    assert mulder_name(a).startswith("ms-quartz-sand-")
    for other in (
        Material.from_formula("Quartz sand!", "SiO2", 2.66, i_ev=139.2),
        Material.from_formula("Quartz sand!", "SiO2", 2.65, i_ev=139.3),
        Material.from_formula("Quartz sand!", "SiO2", 2.65, i_ev="mulder"),
    ):
        assert mulder_name(other) != mulder_name(a)
    assert mulder_name(M.loam(1.6, 12.0).material) != mulder_name(M.loam(1.6, 18.0).material)


def test_names_are_bare_toml_keys():
    """Mulder writes names into its disk cache as bare TOML keys; a dot would defeat the cache."""
    for material in [*M.library.MATERIALS.values(), M.loess(1.5, 10.0).material]:
        name = mulder_name(material)
        assert all(c.isascii() and (c.isalnum() or c == "-") for c in name), name


def test_dry_density_and_saturation_formula():
    s = M.sand(1.70, 10.0)
    assert s.saturation == pytest.approx(0.1 * 1.70 / (1 - 1.70 / 2.65))
    assert math.isclose(s.porosity, 1 - 1.70 / 2.65)
