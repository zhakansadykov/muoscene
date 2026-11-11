import pytest

from muoscene import units


def test_opacity_of_one_metre_of_water():
    assert units.density_length_to_opacity(1.0, 1.0) == pytest.approx(units.MWE)


def test_opacity_of_soil():
    # 20 m of soil at 1.8 g/cm^3
    assert units.density_length_to_opacity(1.8, 20.0) == pytest.approx(3600.0)


def test_constants():
    assert units.MUON_MASS == pytest.approx(0.10566, abs=1e-5)
    assert units.DAY == 24 * 3600
