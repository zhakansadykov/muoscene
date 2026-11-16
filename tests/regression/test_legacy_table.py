"""`LegacyTable` against the preliminary calculation of `legacy.py` and its frozen numbers.

The roadmap asks for parity to 1e-6. The table repeats the arithmetic of `legacy.py` expression
by expression, so the two agree to the last bit; the tolerance below leaves room for a platform
whose mathematical library rounds differently.
"""

import json
import pathlib

import numpy as np
import pytest

from muoscene.flux import LEGACY_RESPONSE, EnergyResponse, LegacyTable
from muoscene.flux.parametrizations import COS_GRID, E_GRID, guan_intensity
from tests.regression import legacy, legacy_cases

FROZEN = json.loads(
    (pathlib.Path(__file__).resolve().parent / "data" / "legacy_reference.json").read_text(encoding="utf-8")
)["cases"]
RTOL = 1e-12
FROZEN_RTOL = 1e-10  # as in test_legacy_reference: last-bit differences between platforms


def preliminary(cos_zenith, opacity, model=None):
    """Integral intensity of the preliminary calculation under an opacity, m^-2 s^-1 sr^-1."""
    model = model or legacy.Model()
    return legacy.integral_intensity(legacy.threshold(model, opacity), cos_zenith) * 1e4


def test_parity_with_the_preliminary_calculation():
    table = LegacyTable()
    rng = np.random.default_rng(4)
    c = np.r_[rng.uniform(0.25, 1.0, 5000), 0.25, 1.0, COS_GRID[7], 0.5 * (COS_GRID[7] + COS_GRID[8])]
    x = np.r_[np.exp(rng.uniform(np.log(0.1), np.log(3.0e4), 5000)), 0.0, 3.0e4, 1800.0, 0.0]
    ours = table.integral(LEGACY_RESPONSE, c, x)
    np.testing.assert_allclose(ours, preliminary(c, x), rtol=RTOL, atol=0)
    np.testing.assert_array_equal(ours, preliminary(c, x))  # the same arithmetic: equal to the last bit
    assert isinstance(table.integral(LEGACY_RESPONSE, 1.0, 1800.0), float)
    assert table.integral(LEGACY_RESPONSE, np.ones((3, 1)), np.array([0.0, 500.0])).shape == (3, 2)


def test_parity_for_other_parameters():
    model = legacy.Model(a_loss=2.5e-3, b_loss=3.5e-6, offset_gev=0.2)
    table = LegacyTable(a_loss_gev_cm2_g=2.5e-3, b_loss_cm2_g=3.5e-6, offset_gev=0.2)
    c, x = np.array([1.0, 0.8, 0.4]), np.array([0.0, 3600.0, 2.0e4])
    np.testing.assert_allclose(table.integral(LEGACY_RESPONSE, c, x), preliminary(c, x, model), rtol=RTOL)
    np.testing.assert_allclose(table.threshold_gev(x), legacy.threshold(model, x), rtol=RTOL)


def test_frozen_numbers():
    table = LegacyTable()
    opacities = np.array([0.0, 100.0, 500.0, 1000.0, 1800.0, 3600.0, 5400.0, 1e4, 3e4])  # of legacy_cases
    np.testing.assert_allclose(table.threshold_gev(opacities), FROZEN["threshold"], rtol=FROZEN_RTOL)
    # the frozen intensities are given against energy: the opacities that need those energies
    energies = np.array([0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 30.0, 100.0, 1000.0])
    cosines = np.array([0.25, 0.5, 0.7071, 0.9, 1.0])
    frozen = 1e4 * np.array(FROZEN["integral_intensity"]).reshape(cosines.size, energies.size)
    m = legacy.Model()
    opacity = np.log1p((energies - m.offset_gev) * m.b_loss / m.a_loss) / m.b_loss
    inside = opacity <= 3.0e4  # 100 GeV needs 4.3e4 g/cm2
    assert inside.sum() == 7
    ours = table.integral(LEGACY_RESPONSE, cosines[:, None], opacity[inside])
    np.testing.assert_allclose(
        ours, frozen[:, inside], rtol=1e-9
    )  # the inversion of the threshold costs 1e-13


def test_vertical_intensity_and_the_spectrum():
    """60.5 m^-2 s^-1 sr^-1 above 1 GeV of total energy; the open sky of this table starts at 0.1 GeV."""
    m = legacy.Model()
    one_gev = np.log1p((1.0 - m.offset_gev) * m.b_loss / m.a_loss) / m.b_loss
    table = LegacyTable()
    assert table.threshold_gev(one_gev) == pytest.approx(1.0)
    assert table.integral(LEGACY_RESPONSE, 1.0, one_gev) == pytest.approx(60.5, abs=0.05)
    assert table.threshold_gev(0.0) == 0.1
    np.testing.assert_array_equal(guan_intensity(E_GRID, 0.6), legacy.dnde(E_GRID, 0.6))
    assert legacy_cases.RHO * 100 * 10.0 == 1800.0  # 10 m of the soil of the preliminary calculation, g/cm2


def test_only_its_own_response():
    """The preliminary calculation has no detector energy response: nothing else is taken, nothing assumed."""
    table = LegacyTable()
    for response in (EnergyResponse.threshold(1.0), None, 1.0, "legacy"):
        with pytest.raises(ValueError, match="no detector energy response"):
            table.integral(response, 1.0, 1800.0)
    with pytest.raises(TypeError):
        table.integral(1.0, 1800.0)
    assert repr(LEGACY_RESPONSE) == "LEGACY_RESPONSE"


@pytest.mark.parametrize(
    ("c", "x"),
    [(0.2, 100.0), (1.01, 100.0), (float("nan"), 100.0), (0.5, -1.0), (0.5, 3.1e4), (0.5, float("nan"))],
)
def test_outside_the_table_is_refused(c, x):
    """The preliminary calculation clips cos(zenith) at 0.25 and never meets more than 3e4 g/cm2."""
    with pytest.raises(ValueError, match="outside the table"):
        LegacyTable().integral(LEGACY_RESPONSE, c, x)


@pytest.mark.parametrize(
    "parameters",
    [
        {"a_loss_gev_cm2_g": 0.0},
        {"b_loss_cm2_g": -1e-6},
        {"offset_gev": 0.05},  # below 0.1 GeV, where the table of the spectrum starts
        {"opacity_max_g_cm2": 1.0e6},  # the threshold would pass 10 TeV, where it ends
        {"opacity_max_g_cm2": 0.0},
    ],
)
def test_parameters_are_checked(parameters):
    with pytest.raises(ValueError):
        LegacyTable(**parameters)
