"""Mulder and Calzone on the synthetic scenes of the spike, frozen in `data/mulder_reference.json`.

A release of Mulder or Calzone that changes any of these numbers by more than 1e-5 fails here:
energy loss, cross-sections, the reference flux, the atmosphere, the geomagnetic field, the
navigation. If the change is understood and accepted, regenerate the file with
`python -m tests.regression.make_mulder_reference` and record why in CHANGELOG.md.

The second half of the file checks what the frozen numbers mean: the statements of
`docs/notes/mulder.md` and the parity with the preliminary calculation of `legacy.py`.
"""

import json
import pathlib

import numpy as np
import pytest

from tests.regression import legacy, legacy_cases, mulder_cases

DATA = pathlib.Path(__file__).resolve().parent / "data"
FROZEN = json.loads((DATA / "mulder_reference.json").read_text(encoding="utf-8"))
CASES = FROZEN["cases"]
RTOL = mulder_cases.RTOL  # above the spread between platforms, below any change of physics
ATOL_PATH = 1e-9  # m, for paths in the ground (zero towards the canyon)


def assert_close(actual, expected, atol, where=""):
    if isinstance(expected, dict):
        assert isinstance(actual, dict) and actual.keys() == expected.keys(), where
        for key in expected:
            assert_close(actual[key], expected[key], atol, f"{where}.{key}")
    elif isinstance(expected, str):
        assert actual == expected, where
    else:
        np.testing.assert_allclose(actual, expected, rtol=RTOL, atol=atol, err_msg=where)


def test_frozen_file_lists_every_case():
    assert sorted(CASES) == sorted(mulder_cases.CASES)


def test_parameters_unchanged():
    assert json.loads(json.dumps(mulder_cases.parameters())) == FROZEN["parameters"]


@pytest.mark.parametrize("name", sorted(mulder_cases.CASES))
def test_case(name):
    atol = ATOL_PATH if name.startswith("wall/") else 0.0
    assert_close(mulder_cases.compute(name), CASES[name], atol, name)


# ------------------------------------------------ what the frozen numbers mean


def gain(case, other, key):
    return CASES[case][key] / CASES[other][key] - 1


def test_reference_is_guan_with_total_energy():
    ref = CASES["reference/differential"]
    k = np.array(ref["kinetic_energy_gev"])
    for elevation in (90, 60, 30):
        guan = legacy.dnde(k + mulder_cases.MU, np.full(k.shape, np.sin(np.radians(elevation)))) * 1e4
        np.testing.assert_allclose(ref[f"elev{elevation}"], guan, rtol=2e-6)


def test_vertical_flux_at_sea_level():
    sky = CASES["open_sky/0m/elev90"]
    assert sky["E>1GeV"] == pytest.approx(60.50, abs=0.01)
    assert sky["all"] == pytest.approx(88.42, abs=0.01)


def test_altitude_gain_of_the_open_sky_at_1150_m():
    case, sea = "open_sky/1150m/elev90", "open_sky/0m/elev90"
    assert gain(case, sea, "E>1GeV") == pytest.approx(0.198, abs=0.001)
    assert gain(case, sea, "E>5GeV") == pytest.approx(0.096, abs=0.001)
    assert gain(case, sea, "E>10GeV") == pytest.approx(0.063, abs=0.001)
    under = CASES["slab/1150m/20m/zenith0"] / CASES["slab/0m/20m/zenith0"] - 1
    assert under == pytest.approx(0.076, abs=0.001)


@pytest.mark.parametrize(("altitude", "expected"), [(1150, (0.28, 0.38, 0.53)), (3340, (0.69, 0.97, 1.40))])
def test_altitude_cutoff(altitude, expected):
    for elevation, k_cut in zip((90, 45, 30), expected, strict=True):
        assert CASES[f"cutoff/{altitude}m/elev{elevation}"] == pytest.approx(k_cut, abs=0.006)
    assert CASES["season/us-standard"]["open_sky/cutoff_gev"] == pytest.approx(
        CASES["cutoff/1150m/elev90"], rel=1e-4
    )


def test_preliminary_threshold_is_too_high():
    loam = np.array(CASES["threshold/RegLoam"]) + mulder_cases.MU
    excess = legacy.threshold(legacy.Model(), np.array(mulder_cases.OPACITIES)) / loam - 1
    np.testing.assert_allclose(excess[:6], [0.16, 0.16, 0.133, 0.089, 0.063, 0.025], atol=0.002)
    assert excess[-1] < 0  # beyond 2e4 g/cm^2 the exponential of a + bE overtakes the real loss


def test_soils_and_concrete_differ_by_a_few_percent():
    loam = np.array(CASES["threshold/RegLoam"])
    for name, worst in (("Rock", 0.025), ("RegConcrete", 0.008)):
        assert np.all(np.abs(np.array(CASES[f"threshold/{name}"]) / loam - 1) < worst)


def test_equivalent_opacity_weight_of_concrete():
    i = mulder_cases.RANGE_ENERGIES.index(7.0)
    weight = CASES["range/RegLoam"][i] / CASES["range/RegConcrete"][i]
    assert weight == pytest.approx(0.995, abs=0.0005)


def test_slab_against_the_preliminary_calculation():
    m = legacy.Model()
    for altitude, (lo, hi) in (("0m", (1.025, 1.125)), ("1150m", (1.08, 1.26))):
        for depth in mulder_cases.DEPTHS:
            for zenith in mulder_cases.ZENITH:
                c = np.cos(np.radians(zenith))
                x = mulder_cases.SOIL / 10 * depth / c
                prelim = legacy.integral_intensity(np.array([legacy.threshold(m, x)]), np.array([c]))[0] * 1e4
                assert lo < CASES[f"slab/{altitude}/{depth:.0f}m/zenith{zenith:g}"] / prelim < hi


def test_wall_paths_equal_the_preliminary_calculation():
    """Equal to round-off, except along the wall and towards the canyon under the rising brow.

    There the preliminary calculation takes the brow as flat; Calzone follows the mesh, which is
    rise * xd / xr higher above the detector.
    """
    m, u = legacy.Model(), legacy_cases.direction_set()
    inward = u[:, 0] > 0
    for bank, g in legacy_cases.BANKS.items():
        for depth in (legacy_cases.UPPER, legacy_cases.LOWER):
            ours = np.array(CASES[f"wall/{bank}/{depth:.0f}m"])
            prelim = legacy.soil_length(m, u, depth, **g)
            if g["rise"] == 0:
                np.testing.assert_allclose(ours, prelim, rtol=0, atol=ATOL_PATH)
                continue
            np.testing.assert_allclose(ours[inward], prelim[inward], rtol=0, atol=ATOL_PATH)
            extra = ours[~inward] - prelim[~inward]
            assert np.all(extra >= -ATOL_PATH)
            assert np.all(extra <= g["rise"] * m.niche.xd / g["xr"] / u[~inward, 2] + ATOL_PATH)


def test_reference_void():
    """10.3 sigma of the application becomes 10.7 sigma at sea level and 12.0 sigma at 1150 m."""
    prelim = json.loads((DATA / "legacy_reference.json").read_text(encoding="utf-8"))["cases"][
        "reference_void"
    ]
    sea, site = CASES["reference_void/0m/cone1"], CASES["reference_void/1150m/cone1"]
    assert sea["sig"] == pytest.approx(10.73, abs=0.01)
    assert site["sig"] == pytest.approx(12.05, abs=0.01)
    assert prelim["sig"] < sea["sig"] < site["sig"]
    assert sea["effect_percent"] == pytest.approx(prelim["eff"], abs=0.15)


def test_season_changes_the_flux_by_less_than_one_percent():
    summer, winter = CASES["season/midlatitude-summer"], CASES["season/midlatitude-winter"]
    assert winter["open_sky/cutoff_gev"] > summer["open_sky/cutoff_gev"]  # denser winter air
    for key in ("open_sky/all", "open_sky/E>1GeV", "under_20m/all", "under_20m/E>1GeV"):
        assert abs(winter[key] / summer[key] - 1) < 0.01


def test_geomagnetic_field_east_against_west():
    g = CASES["geomagnet/1150m/elev45"]
    assert g["east/all"] / g["west/all"] - 1 == pytest.approx(0.0076, abs=0.0005)
    assert g["east/E>1GeV"] / g["west/E>1GeV"] - 1 == pytest.approx(0.0027, abs=0.0005)


def test_physics_models_are_the_defaults_of_the_spike():
    assert CASES["physics/defaults"] == {
        "bremsstrahlung": "SSR19",
        "pair_production": "SSR19",
        "photonuclear": "DRSS01",
        "cutoff": 0.05,
        "elastic_ratio": 0.05,
    }
