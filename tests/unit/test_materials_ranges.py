"""Materials in Mulder: ranges against PDG, the thresholds of the spike, weights, equivalent opacity.

These tests compile Mulder physics tables: about a minute on a machine that has none cached.
"""

import json
import pathlib

import mulder.materials as mm
import numpy as np
import pytest

from muoscene import materials as M
from muoscene.materials.backend import mulder_density_kg_m3, register
from muoscene.materials.composition import Material, Mixture
from tests.regression import mulder_cases

K_PDG = np.array([0.3, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0, 1000.0])  # kinetic energy, GeV
PDG_RANGE = {  # CSDA range, g/cm^2: PDG muon tables (D. E. Groom et al., ADNDT 78 (2001) 183)
    # pdg.lbl.gov/2024/AtomicNuclearProperties/MUE/muE_standard_rock.txt
    "standard rock": [1.523e2, 5.534e2, 1.599e3, 4.920e3, 1.362e4, 4.084e4, 1.041e5, 2.453e5],
    # pdg.lbl.gov/2022/AtomicNuclearProperties/MUE/muE_shielding_concrete.txt
    "concrete": [1.506e2, 5.460e2, 1.577e3, 4.855e3, 1.346e4, 4.038e4, 1.030e5, 2.434e5],
    # pdg.web.cern.ch/2020/AtomicNuclearProperties/MUE/muE_water_liquid.pdf
    "water": [1.291e2, 4.732e2, 1.378e3, 4.279e3, 1.195e4, 3.642e4, 9.614e4, 2.433e5],
}
LIBRARY = {"standard rock": M.STANDARD_ROCK, "concrete": M.CONCRETE, "water": M.WATER}


@pytest.mark.parametrize("name", sorted(PDG_RANGE))
def test_ranges_agree_with_pdg(name):
    """Within 0.6 % from 0.3 GeV to 1 TeV (largest: water, -0.52 % at 1 GeV)."""
    ratio = M.csda_range(LIBRARY[name], K_PDG) / np.array(PDG_RANGE[name])
    np.testing.assert_allclose(ratio, 1.0, atol=0.006)


def test_mulder_default_i_shortens_the_range_of_water():
    plain = Material.from_formula("water", "H2O", 1.0, i_ev="mulder")
    ratio = M.csda_range(plain, [0.3, 1.0, 10.0]) / M.csda_range(M.WATER, [0.3, 1.0, 10.0])
    assert np.all(ratio < 1)
    assert ratio[1] == pytest.approx(0.995, abs=0.002)


def test_registration_is_idempotent_and_keeps_densities():
    soil = M.loam(1.60, 15.0).material
    name = register(soil)
    assert register(soil) == name
    assert name in mm.Composite.all()
    assert mulder_density_kg_m3(soil) == pytest.approx(1e3 * soil.density_g_cm3, rel=1e-12)
    assert mulder_density_kg_m3(M.QUARTZ) == pytest.approx(2650.0)


# ------------------------------------------------ the spike, through the muoscene interface


def spike_materials():
    """The loam and concrete of the spike (roadmap 0.1-02), with Mulder's default I as there."""
    quartz = Material.from_formula("spike quartz", "SiO2", 2.65, i_ev="mulder")
    calcite = Material.from_formula("spike calcite", "CaCO3", 2.71, i_ev="mulder")
    alumina = Material.from_formula("spike alumina", "Al2O3", 3.95, i_ev="mulder")
    water = Material.from_formula("spike water", "H2O", 1.02, i_ev="mulder")  # Mulder's Water
    loam = Mixture.of("spike loam", {quartz: 0.60, calcite: 0.15, alumina: 0.12, water: 0.13})
    concrete = Material.from_elements("spike concrete", mulder_cases.CONCRETE, 2.30, i_ev="mulder")
    return loam, concrete


def test_spike_thresholds_are_reproduced():
    frozen = json.loads(
        (pathlib.Path(mulder_cases.__file__).parent / "data" / "mulder_reference.json").read_text(
            encoding="utf-8"
        )
    )["cases"]
    loam, concrete = spike_materials()
    x = np.array(mulder_cases.OPACITIES)
    for material, case in ((loam, "threshold/RegLoam"), (concrete, "threshold/RegConcrete")):
        np.testing.assert_allclose(M.energy_to_cross(material, x), frozen[case], rtol=mulder_cases.RTOL)


def test_concrete_weighs_0995_against_loam_at_7_gev():
    loam, concrete = spike_materials()
    assert M.weight(concrete, loam, 7.0) == pytest.approx(0.995, abs=0.0005)
    # library: concrete against a loam with 15 % of water is closer still
    assert M.weight(M.CONCRETE, M.loam(1.60, 15.0).material, 7.0) == pytest.approx(0.9995, abs=0.0005)


# ------------------------------------------------ range functions


def test_range_and_inverse_range_round_trip():
    k = np.logspace(-1, 3, 9)
    for material in (M.STANDARD_ROCK, M.loam(1.6, 15.0).material):
        # Mulder interpolates its range tables: the round trip holds to 6e-6
        np.testing.assert_allclose(M.energy_to_cross(material, M.csda_range(material, k)), k, rtol=1e-5)
    assert isinstance(M.csda_range(M.WATER, 10.0), float)
    with pytest.raises(ValueError):
        M.energy_to_cross(M.WATER, -1.0)


def test_weight_and_equivalent_opacity():
    x = np.array([500.0, 3600.0, 3e4])
    np.testing.assert_allclose(M.equivalent_opacity(M.CONCRETE, M.CONCRETE, x), x, rtol=1e-6)
    assert M.weight(M.WATER, M.WATER, 5.0) == pytest.approx(1.0)
    # one material: the equivalent opacity needs the same energy in the reference
    rock_eq = M.equivalent_opacity(M.CONCRETE, M.STANDARD_ROCK, x)
    np.testing.assert_allclose(
        M.energy_to_cross(M.STANDARD_ROCK, rock_eq), M.energy_to_cross(M.CONCRETE, x), rtol=1e-6
    )
    # and the weight at that energy is the ratio of the opacities
    k = M.energy_to_cross(M.CONCRETE, x)
    np.testing.assert_allclose(M.weight(M.CONCRETE, M.STANDARD_ROCK, k), rock_eq / x, rtol=1e-6)


def test_water_raises_the_energy_needed_at_fixed_opacity():
    """More electrons per gram: +1.6 % per 10 % of water content at 1000 g/cm^2 (docs/notes/materials.md)."""
    dry = M.energy_to_cross(M.SILICATE_SKELETON, 1000.0)
    wet = M.energy_to_cross(M.loam(1.50, 10.0).material, 1000.0)
    assert wet / dry - 1 == pytest.approx(0.0167, abs=0.001)
