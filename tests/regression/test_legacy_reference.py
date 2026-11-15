"""The preliminary calculation, frozen on a synthetic canyon (`data/legacy_reference.json`).

The port in `legacy.py` was checked to be bit-identical to the original scripts. These tests keep
it so: any change of its numbers fails here. Regenerate the file only for a deliberate change,
with `python -m tests.regression.make_legacy_reference`, and record the reason in CHANGELOG.md.
"""

import json
import pathlib

import numpy as np
import pytest

from tests.regression import legacy, legacy_cases

FROZEN = json.loads(
    (pathlib.Path(__file__).resolve().parent / "data" / "legacy_reference.json").read_text(encoding="utf-8")
)
RTOL = 1e-10  # same arithmetic; leaves room only for last-bit differences between platforms


def assert_same(actual, expected, where=""):
    if isinstance(expected, dict):
        assert isinstance(actual, dict) and actual.keys() == expected.keys(), where
        for key in expected:
            assert_same(actual[key], expected[key], f"{where}.{key}")
    elif expected is None or isinstance(expected, bool):
        assert actual == expected, where
    else:
        np.testing.assert_allclose(actual, expected, rtol=RTOL, atol=0, err_msg=where)


def test_frozen_file_lists_every_case():
    assert sorted(FROZEN["cases"]) == sorted(legacy_cases.CASES)


def test_parameters_unchanged():
    assert json.loads(json.dumps(legacy_cases.parameters())) == FROZEN["parameters"]


@pytest.mark.parametrize("name", sorted(legacy_cases.CASES))
def test_case(name):
    assert_same(legacy_cases.compute(name), FROZEN["cases"][name], name)


# ------------------------------------------------ what the frozen numbers mean


def test_reference_void_is_the_number_of_the_application():
    """1 m void 3 m from the wall and 3 m under a flat brow, niche 10 m, 240 cm^2, 60 days: 10.3 sigma."""
    result = FROZEN["cases"]["reference_void"]
    assert result["inside"]
    assert result["sig"] == pytest.approx(10.31, abs=0.01)
    assert result["eff"] == pytest.approx(6.17, abs=0.01)


def test_vertical_intensity_above_1_gev():
    """Guan et al. (2015) at sea level: 60.5 m^-2 s^-1 sr^-1 above 1 GeV of total energy."""
    value = legacy.integral_intensity(np.array([1.0]), np.array([1.0]))[0] * 1e4
    assert value == pytest.approx(60.5, abs=0.05)


def test_threshold_is_monotonic_and_starts_below_the_muon_mass():
    e = np.array(FROZEN["cases"]["threshold"])
    assert np.all(np.diff(e) > 0)
    assert e[0] == pytest.approx(0.1)  # below 0.1057 GeV: a known simplification, docs/legacy.md


def test_paths_are_symmetric_along_the_wall():
    u = legacy_cases.direction_set()
    mirrored = u * np.array([1.0, -1.0, 1.0])
    m = legacy.Model()
    for g in legacy_cases.BANKS.values():
        for depth in (legacy_cases.UPPER, legacy_cases.LOWER):
            np.testing.assert_array_equal(
                legacy.soil_length(m, u, depth, **g), legacy.soil_length(m, mirrored, depth, **g)
            )


def test_rising_brow_adds_ground_towards_the_massif():
    u = legacy_cases.direction_set()
    flat = np.array(FROZEN["cases"]["soil_length/flat/22m"])
    rising = np.array(FROZEN["cases"]["soil_length/rising/22m"])
    into = u[:, 0] > 0
    assert np.all(rising[into] >= flat[into] - 1e-12)
    assert np.all(rising[~into] == flat[~into])


def test_canyon_side_directions_cross_only_the_lip_above_the_niche():
    """Towards the canyon the ground ends at the wall face: only the lip of the niche roof remains."""
    m = legacy.Model()
    u = legacy_cases.direction_set()
    out = u[:, 0] < 0
    to_face = m.niche.xd / -u[out, 0]  # distance to the plane of the wall face
    for bank in legacy_cases.BANKS:
        for depth in (legacy_cases.UPPER, legacy_cases.LOWER):
            path = np.array(FROZEN["cases"][f"soil_length/{bank}/{depth:.0f}m"])[out]
            assert np.all(path <= to_face + 1e-12)
            assert path[u[out, 0] < -0.5].max() <= 0.5  # metres, against 12-22 m towards the massif


def test_open_sky_track_density_does_not_depend_on_the_bank():
    cases = FROZEN["cases"]
    assert cases["track_density/flat/0m"] == cases["track_density/rising/0m"]
    assert 3e4 < cases["track_density/flat/22m"] < 4e4  # (3-4)e4 per cm^2 in 60 days, as in the application


def test_tilting_towards_the_massif_reduces_the_open_sky_reference():
    cases = FROZEN["cases"]
    for bank in legacy_cases.BANKS:
        assert cases[f"reference_fraction/{bank}/22m/beta20"] < cases[f"reference_fraction/{bank}/22m/beta0"]
