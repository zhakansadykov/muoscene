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


def test_detector_threshold_lowers_the_significance():
    """12.05 sigma at 1150 m counts every muon; above 0.2 GeV it is 11.52, above 1 GeV 9.37.

    The preliminary calculation has no detector threshold; the conventions default to 1 GeV.
    """
    for altitude, expected in (("0m", (10.30, 8.48)), ("1150m", (11.52, 9.37))):
        plain = CASES[f"reference_void/{altitude}/cone1"]
        low, high = (CASES[f"reference_void/{altitude}/cone1/K>{k:g}GeV"] for k in mulder_cases.THRESHOLDS)
        assert (low["sig"], high["sig"]) == pytest.approx(expected, abs=0.01)
        assert high["N0"] < low["N0"] < plain["N0"]
        assert high["effect_percent"] < low["effect_percent"] < plain["effect_percent"]


def test_cutoff_survives_under_thin_soil():
    """Half a metre of soil loses less energy than the air down to sea level: the cut-off stays."""
    name = "thin_slab/1150m/0.5m/zenith{}"
    assert CASES[name.format(0)]["cutoff_gev"] == pytest.approx(0.124, abs=0.002)
    assert CASES[name.format(45)]["cutoff_gev"] == pytest.approx(0.162, abs=0.002)
    for zenith in (0, 45):
        case = CASES[name.format(zenith)]
        assert 0 < case["cutoff_gev"] < CASES["cutoff/1150m/elev90" if zenith == 0 else "cutoff/1150m/elev45"]
        assert case["E>1GeV"] < case["all"]


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


def test_below_sea_level_the_flux_is_lower_and_has_no_cutoff():
    """At -430 m the reference is carried down through 430 m more air: 56.69 against 60.50 above 1 GeV."""
    low, sea = CASES["below_sea_level/-430m/elev90"], CASES["open_sky/0m/elev90"]
    assert low["cutoff_gev"] == 0
    assert low["E>1GeV"] == pytest.approx(56.69, abs=0.01)
    for key in ("E>1GeV", "E>5GeV", "E>10GeV"):
        assert low[key] < sea[key]
    assert low["E>10GeV"] / sea["E>10GeV"] > low["E>1GeV"] / sea["E>1GeV"]  # the soft muons lose most
    assert low["K>1GeV"] < low["E>1GeV"] < low["K>0.2GeV"]  # E > 1 GeV is K > 0.894 GeV


def test_detector_threshold_under_a_slab():
    """Under 20 m of soil a threshold of 1 GeV takes 15 % of the muons above 0.2 GeV, under 0.5 m a third."""
    for depth, loss in ((0.5, 0.333), (2.0, 0.288), (20.0, 0.149)):
        case = CASES[f"threshold_slab/1150m/{depth:g}m/zenith0"]
        assert 1 - case["K>1GeV"] / case["K>0.2GeV"] == pytest.approx(loss, abs=0.001)
    for zenith in (0, 45):
        thick, mid, thin = (CASES[f"threshold_slab/1150m/{d}m/zenith{zenith}"] for d in ("20", "2", "0.5"))
        for key in ("K>0.2GeV", "K>1GeV"):
            assert thick[key] < mid[key] < thin[key]
        # the lower threshold is above the cut-off that half a metre of soil leaves
        assert CASES[f"thin_slab/1150m/0.5m/zenith{zenith}"]["cutoff_gev"] < mulder_cases.THRESHOLDS[0]
        assert thin["K>1GeV"] < CASES[f"thin_slab/1150m/0.5m/zenith{zenith}"]["E>1GeV"]
