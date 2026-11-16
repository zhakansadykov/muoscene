"""`MulderTable`: the criteria of roadmap 0.1-05 against Mulder and its frozen numbers.

The reference material is the loam of the spike (roadmap 0.1-02), with Mulder's default mean
excitation energies as there, so that the frozen cases of `tests/regression/mulder_cases.py`
apply. Tables are built without a cache or into `tmp_path`: nothing is written elsewhere.
A table over the default grid takes 2 s (Linux VM, 4 cores); most tests build a band of two
elevation nodes around the direction they need, 0.03 s.

A threshold "E > 1 GeV" is one in total energy, as in the spike note: K > 1 GeV - m.
"""

import functools
import json
import pathlib

import mulder
import numpy as np
import pytest

from muoscene import materials as M
from muoscene.flux import EnergyResponse, FluxTable, MulderTable, atmosphere_models, default_cache_dir
from muoscene.flux.atmosphere import fluxmeter, open_sky_cutoff, states
from muoscene.flux.tables import CUTOFF_MARGIN_GEV, DENSE_SLAB_KG_M3, K_MAX_GEV, check_sample, slab_flux
from muoscene.materials.backend import register
from muoscene.materials.composition import Material, Mixture
from tests.regression import mulder_cases as mc

FROZEN = json.loads(
    (pathlib.Path(mc.__file__).parent / "data" / "mulder_reference.json").read_text(encoding="utf-8")
)["cases"]
SITE = {"latitude_deg": mc.SITE["latitude"], "longitude_deg": mc.SITE["longitude"]}
SIN45 = float(np.sin(np.radians(45.0)))
STEP = 0.005  # spacing of the default grid in the sine of elevation
ABOVE_1_GEV_TOTAL = EnergyResponse.threshold(1.0 - mc.MU)
FROM_1_MEV = EnergyResponse.threshold(1.0e-3)  # every muon the table holds where there is no cut-off
# The frozen integrals are trapezoids on 801 nodes up to 1e5 GeV, which sit 1.6e-5 above their
# limit; the tables integrate their own interpolation, 0.6e-5 below it. Open-sky numbers of a table
# are 2.2e-5 to 3.5e-5 below the frozen ones; the tolerance leaves a factor of two.
RTOL = 7.0e-5


def spike_loam():
    """The loam of the spike: the `RegLoam` of the frozen cases as a muoscene material."""
    quartz = Material.from_formula("spike quartz", "SiO2", 2.65, i_ev="mulder")
    calcite = Material.from_formula("spike calcite", "CaCO3", 2.71, i_ev="mulder")
    alumina = Material.from_formula("spike alumina", "Al2O3", 3.95, i_ev="mulder")
    water = Material.from_formula("spike water", "H2O", 1.02, i_ev="mulder")  # Mulder's Water
    return Mixture.of("spike loam", {quartz: 0.60, calcite: 0.15, alumina: 0.12, water: 0.13})


LOAM = spike_loam()


def nodes(sine):
    """Two nodes of the default grid around a sine of elevation."""
    low = np.floor(min(sine, 1.0 - STEP) / STEP + 1e-9) * STEP
    return (round(low, 6), round(low + STEP, 6))


@functools.cache
def band(altitude_m, sine, **options):
    """Table at an altitude over the two default nodes around `sine`, not cached on disk."""
    return MulderTable(altitude_m, LOAM, cache=False, sin_elevation=nodes(sine), **(SITE | options))


@pytest.fixture(scope="module")
def site():
    """The table of a detector at 1150 m over the default grid."""
    return MulderTable(mc.ALTITUDE, LOAM, cache=False, **SITE)


def sum_from_the_first_node(table, sine, k_min=1.0, k_max=1.0e4, n=61):
    """The integral as the spike took it: every node, the first one included, gets a whole step."""
    k = np.logspace(np.log10(k_min), np.log10(k_max), n)
    return float(np.sum(table.differential(k, sine, 0.0) * k * np.gradient(np.log(k))))


# ------------------------------------------------ open sky at sea level


def test_vertical_flux_at_sea_level():
    """60.50 m^-2 s^-1 sr^-1 above 1 GeV of total energy; 88.42 over all energies."""
    sea = band(0.0, 1.0)
    frozen = FROZEN["open_sky/0m/elev90"]
    above = sea.integral(ABOVE_1_GEV_TOTAL, 1.0, 0.0)
    assert above == pytest.approx(60.50, abs=0.01)
    assert above == pytest.approx(frozen["E>1GeV"], rel=RTOL)
    for e in (5.0, 10.0):
        value = sea.integral(EnergyResponse.threshold(e - mc.MU), 1.0, 0.0)
        assert value == pytest.approx(frozen[f"E>{e:g}GeV"], rel=RTOL)
    # the table starts at 1 MeV: the first MeV, 0.042 at 41.7 per GeV, is all that is missing of 88.42
    first_mev = sea.differential(1.0e-3, 1.0, 0.0) * 1.0e-3
    assert sea.integral(FROM_1_MEV, 1.0, 0.0) + first_mev == pytest.approx(frozen["all"], rel=RTOL)
    assert sea.cutoff(1.0, 0.0) == 0.0
    oblique = band(0.0, SIN45).integral(ABOVE_1_GEV_TOTAL, SIN45, 0.0)
    assert oblique == pytest.approx(FROZEN["open_sky/0m/elev45"]["E>1GeV"], rel=RTOL)


def test_a_sum_from_the_first_node_is_two_to_three_percent_high(site):
    """Why the integral is a trapezoid sum from the threshold: 61 nodes from 1 GeV, as in the spike.

    Vertical: +2.85 % at sea level, +3.29 % at 1150 m; at 1150 m it is +2.45 % at 45 degrees,
    +1.88 % at 30 degrees and +2.02 % at 17.5 degrees, where the cut-off is close to 1 GeV.
    """
    above_1_gev = EnergyResponse.threshold(1.0)
    sea = band(0.0, 1.0)
    assert sum_from_the_first_node(sea, 1.0) / sea.integral(above_1_gev, 1.0, 0.0) - 1 == pytest.approx(
        0.0285, abs=0.0005
    )
    for sine, excess in ((1.0, 0.0329), (0.71, 0.0245), (0.5, 0.0188), (0.3, 0.0202)):
        ratio = sum_from_the_first_node(site, sine) / site.integral(above_1_gev, sine, 0.0)
        assert ratio - 1 == pytest.approx(excess, abs=0.0005)


# ------------------------------------------------ altitude cut-off


def test_cutoff_numbers_of_the_spike(site):
    """0.28, 0.38, 0.53 GeV at 1150 m and 0.69, 0.97, 1.40 GeV at 3340 m for 90, 45 and 30 degrees."""
    for elevation, expected in ((90, 0.28), (45, 0.38), (30, 0.53)):
        sine = float(np.sin(np.radians(elevation)))
        assert site.cutoff(sine, 0.0) == pytest.approx(expected, abs=0.006)
        # at a node the cut-off is the bisection itself; between nodes it is interpolated to 3e-7
        assert site.cutoff(sine, 0.0) == pytest.approx(FROZEN[f"cutoff/1150m/elev{elevation}"], rel=2e-6)
    for elevation, expected in ((90, 0.69), (45, 0.97), (30, 1.40)):
        sine = float(np.sin(np.radians(elevation)))
        high = band(mc.HIGH, sine)
        assert high.cutoff(sine, 0.0) == pytest.approx(expected, abs=0.006)
        assert high.cutoff(sine, 0.0) == pytest.approx(FROZEN[f"cutoff/3340m/elev{elevation}"], rel=2e-6)


def test_the_table_holds_only_energies_above_the_cutoff(site):
    assert site.energy_min(1.0, 0.0) == pytest.approx(site.cutoff(1.0, 0.0) + CUTOFF_MARGIN_GEV)
    with pytest.raises(ValueError, match="cut-off of 0.2834 GeV"):
        site.differential(0.2, 1.0, 0.0)
    with pytest.raises(ValueError, match="detector threshold 0.2 GeV is below 0.2934 GeV"):
        site.integral(EnergyResponse.threshold(0.2), 1.0, 0.0)
    assert site.differential(0.3, 1.0, 0.0) > 0
    # towards the cut-off the flux rises as (K - K_cut)^-0.4: from 67 to 170 per GeV between 100 and 10 MeV
    near, far = (site.differential(site.cutoff(1.0, 0.0) + d, 1.0, 0.0) for d in (0.0101, 0.1))
    assert np.log(near / far) / np.log(0.0101 / 0.1) == pytest.approx(-0.40, abs=0.02)
    assert (far, near) == pytest.approx((67.2, 169.3), rel=0.01)
    # a threshold of 1 GeV holds from 15.7 degrees up: the cut-off is 1.0 GeV at 15.2 degrees
    assert np.all(site.integral(EnergyResponse.threshold(1.0), site.sin_elevation[4:], 0.0) > 0)
    assert site.sin_elevation[4] == pytest.approx(0.27)
    assert site.cutoff(0.25, 0.0) == pytest.approx(1.058, abs=0.001)
    with pytest.raises(ValueError, match="detector threshold 1 GeV"):
        site.integral(EnergyResponse.threshold(1.0), 0.26, 0.0)  # 15.1 degrees
    assert site.cutoff(0.26, 0.0) == pytest.approx(1.016, abs=0.001)
    assert site.cutoff(0.27, 0.0) == pytest.approx(0.977, abs=0.001)
    # under soil the cut-off shrinks and ends at 138 g/cm2 by the vertical, 569 g/cm2 at 14.5 degrees
    assert site.cutoff(1.0, 137.0) > 0
    assert site.cutoff(1.0, 139.0) == 0
    assert site.cutoff(0.25, 568.0) > 0
    assert site.cutoff(0.25, 571.0) == 0
    assert site.integral(EnergyResponse.threshold(0.2), 1.0, 150.0) > 0
    assert site.energy_min(1.0, 1000.0) == 1.0e-3


def test_cutoff_survives_under_thin_soil():
    """Frozen `thin_slab`: 0.5 m of soil at 1150 m leave a cut-off of 0.12 GeV, 0.16 GeV at 45 degrees."""
    for zenith in mc.THIN_ZENITH:
        sine = float(np.cos(np.radians(zenith)))
        opacity = mc.SOIL / 10 * mc.THIN_DEPTH / sine  # g/cm2
        frozen = FROZEN[f"thin_slab/1150m/{mc.THIN_DEPTH:g}m/zenith{zenith:g}"]
        table = band(mc.ALTITUDE - mc.THIN_DEPTH, sine)
        # the table is at the altitude of the detector, 0.5 m below the surface of the frozen case:
        # half a metre of air less under it is 1.1e-4 GeV of cut-off by the vertical, 9e-4 of its value
        assert table.cutoff(sine, opacity) == pytest.approx(frozen["cutoff_gev"], rel=1.5e-3)
        at_surface = band(mc.ALTITUDE, sine).cutoff(sine, opacity)  # range table against bisection
        assert at_surface == pytest.approx(frozen["cutoff_gev"], rel=2e-5)
        # -8e-5: the half metre of air, and the trapezoids of the frozen number (see RTOL)
        assert table.integral(ABOVE_1_GEV_TOTAL, sine, opacity) == pytest.approx(frozen["E>1GeV"], rel=1.5e-4)
        with pytest.raises(ValueError, match="cut-off"):
            table.differential(0.1, sine, opacity)
        with pytest.raises(ValueError, match="cut-off"):
            table.integral(EnergyResponse.threshold(0.1), sine, opacity)
        assert table.integral(EnergyResponse.threshold(0.2), sine, opacity) > 0  # 0.2 GeV holds under it


# ------------------------------------------------ under soil


def test_detector_thresholds_under_soil():
    """Frozen `threshold_slab`: above 0.2 and 1 GeV under 0.5, 2 and 20 m of real soil at 1150 m.

    The table of a detector is made at its altitude, so `depth` of soil above it comes with as
    much air that is not there: -0.007 to -0.013 % under 0.5 m, -0.018 % under 2 m and
    -0.074 to -0.077 % under 20 m.
    """
    for depth, low, high in ((0.5, -1.5e-4, -0.5e-4), (2.0, -2.2e-4, -1.5e-4), (20.0, -8.5e-4, -7.0e-4)):
        for zenith in mc.THIN_ZENITH:
            sine = float(np.cos(np.radians(zenith)))
            opacity = mc.SOIL / 10 * depth / sine
            frozen = FROZEN[f"threshold_slab/1150m/{depth:g}m/zenith{zenith:g}"]
            table = band(mc.ALTITUDE - depth, sine)
            for k_min in mc.THRESHOLDS:
                ratio = (
                    table.integral(EnergyResponse.threshold(k_min), sine, opacity) / frozen[f"K>{k_min:g}GeV"]
                )
                assert low < ratio - 1 < high


@pytest.mark.parametrize("surface", [0.0, mc.ALTITUDE])
def test_slab_numbers_of_the_spike(surface):
    """The 24 frozen cases `slab`: every muon under 10, 20 and 30 m of soil, zenith angles 0 to 60 degrees.

    A table at the altitude of the detector is 0.05 to 0.10 % below them: the air it puts in the
    place of the soil. A table at the altitude of the surface is 0.03 to 0.05 % above them: the
    frozen numbers are for soil of 1.8 g/cm3, where a few muons decay on the way, and a table
    takes no decay in the opacity.
    """
    for depth in mc.DEPTHS:
        for zenith in mc.ZENITH:
            sine = float(np.cos(np.radians(zenith)))
            opacity = mc.SOIL / 10 * depth / sine
            frozen = FROZEN[f"slab/{surface:.0f}m/{depth:.0f}m/zenith{zenith:g}"]
            at_detector = band(surface - depth, sine).integral(FROM_1_MEV, sine, opacity) / frozen - 1
            at_surface = band(surface, sine).integral(FROM_1_MEV, sine, opacity) / frozen - 1
            assert -1.1e-3 < at_detector < -5.0e-4  # measured: -0.053 to -0.104 %
            assert 2.5e-4 < at_surface < 5.5e-4  # measured: +0.030 to +0.052 %


def test_more_than_3e4_g_cm2_is_refused(site):
    assert site.integral(EnergyResponse.threshold(1.0), 0.5, 3.0e4) > 0
    with pytest.raises(ValueError, match="continuous transport is not valid beyond"):
        site.integral(EnergyResponse.threshold(1.0), 0.5, 3.01e4)
    with pytest.raises(ValueError, match="opacity"):
        site.differential(10.0, 0.5, 1.0e5)
    with pytest.raises(ValueError, match="sin elevation 0.2 is outside"):
        site.integral(EnergyResponse.threshold(1.0), 0.2, 100.0)


# ------------------------------------------------ against Mulder itself


def test_table_is_within_02_percent_of_mulder_at_random_points(site):
    """Differential flux under a dense slab at 20 000 random energies, elevations and opacities.

    Up to 1 TeV at the detector, 16 150 points: rms 0.0035 %, largest difference 0.045 %. Above
    1 TeV, 3850 points: rms 0.044 %, largest 0.79 %, because there the flux of Mulder under a
    thick slab is not smooth in energy (next test).
    """
    energy, sine, opacity = check_sample(site, 20000, seed=1, k_max_gev=K_MAX_GEV)
    ratio = site.differential(energy, sine, opacity) / slab_flux(
        site.meta["inputs"], LOAM, energy, sine, opacity
    )
    low = energy <= 1.0e3
    assert low.sum() == 16150
    assert np.abs(ratio[low] - 1).max() < 2.0e-3  # the criterion of the roadmap
    assert np.abs(ratio[low] - 1).max() == pytest.approx(4.5e-4, abs=1.0e-4)
    assert np.std(ratio[low]) == pytest.approx(3.5e-5, abs=1.0e-5)
    assert np.abs(ratio[~low] - 1).max() == pytest.approx(7.9e-3, abs=1.0e-3)
    assert np.std(ratio[~low]) == pytest.approx(4.4e-4, abs=1.0e-4)
    # the same comparison is made when a table is built, and kept with it
    check = site.meta["check"]
    assert check["points"] == 64 and check["max"] < 2.0e-3 and check["rms"] < 1.0e-4


def test_mulder_ripples_above_1_tev_under_a_thick_slab(site):
    """Why the comparison of a new table stops at 1 TeV: the flux of Mulder wavers around 2 TeV.

    Vertical, under 2.5e4 g/cm2, 1.2 to 3.5 TeV at the detector. Largest deviation of ln flux
    from a polynomial of degree 4 in ln K: 0.34 % for Mulder under the slab, 0.09 % for the
    table, 0.0003 % for Mulder in the open sky at the energies above the slab. The table is up
    to 0.52 % from Mulder there. Muons above 1 TeV at the detector are 0.09 % of the flux above
    1 GeV under that slab, so the ripple is 3e-6 of an integral.
    """
    k = np.logspace(np.log10(1200.0), np.log10(3500.0), 47)
    one, opacity = np.ones(k.size), np.full(k.size, 2.5e4)
    direct = slab_flux(site.meta["inputs"], LOAM, k, one, opacity)
    ours = site.differential(k, one, opacity)
    top = M.energy_to_cross(LOAM, M.csda_range(LOAM, k) + 2.5e4)
    sky = np.asarray(mulder.Fluxmeter().flux(energy=top, **states(one * mc.ALTITUDE, one, **SITE)))

    def waviness(energy, flux):
        fit = np.polynomial.Polynomial.fit(np.log(energy), np.log(flux), 4)
        return float(np.abs(np.log(flux) - fit(np.log(energy))).max())

    assert waviness(k, direct) == pytest.approx(3.4e-3, abs=3.0e-4)
    assert waviness(k, ours) == pytest.approx(8.8e-4, abs=1.0e-4)
    assert waviness(top, sky) < 1.0e-5
    assert np.abs(ours / direct - 1).max() == pytest.approx(5.2e-3, abs=5.0e-4)
    above_1_tev = site.integral(EnergyResponse.threshold(1.0e3), 1.0, 2.5e4)
    share = above_1_tev / site.integral(EnergyResponse.threshold(1.0), 1.0, 2.5e4)
    assert share == pytest.approx(8.9e-4, abs=1.0e-4)


def test_integrals_are_those_of_mulder_under_a_slab(site):
    """Above 0.2 and 1 GeV in 12 random directions each: within 6e-5 of trapezoids on 1601 nodes."""
    rng = np.random.default_rng(8)
    sky = fluxmeter(atmosphere="us-standard")
    layer = mulder.Layer(mc.ALTITUDE, material=register(LOAM), density=DENSE_SLAB_KG_M3)
    slab = fluxmeter(layer, atmosphere="us-standard")
    worst = 0.0
    for k_min in mc.THRESHOLDS:
        done = 0
        while done < 12:
            sine = float(rng.uniform(0.25, 1.0))
            opacity = 0.0 if rng.random() < 0.15 else float(np.exp(rng.uniform(0.0, np.log(3.0e4))))
            if site.energy_min(sine, opacity) > k_min:
                with pytest.raises(ValueError):
                    site.integral(EnergyResponse.threshold(k_min), sine, opacity)
                continue
            depth = opacity * 10 / DENSE_SLAB_KG_M3 * sine
            state = {k: float(v) for k, v in states(mc.ALTITUDE - depth, sine, **SITE).items()}
            direct = mc.integral(sky if opacity == 0 else slab, k_min, n=1601, **state)
            ours = site.integral(EnergyResponse.threshold(k_min), sine, opacity)
            worst = max(worst, abs(ours / direct - 1))
            done += 1
    assert worst == pytest.approx(5.6e-5, abs=2.0e-5)


def test_muons_above_the_table_do_not_matter():
    """A table ends at 1e4 GeV: 3e-6 of the muons above 1 GeV under 3e4 g/cm2, 2e-8 in the open sky."""
    layer = mulder.Layer(mc.ALTITUDE, material=register(LOAM), density=DENSE_SLAB_KG_M3)
    for meter, opacity, share in (
        (fluxmeter(atmosphere="us-standard"), 0.0, 1.7e-8),
        (fluxmeter(layer, atmosphere="us-standard"), 3.0e4, 3.2e-6),
    ):
        depth = opacity * 10 / DENSE_SLAB_KG_M3
        state = {k: float(v) for k, v in states(mc.ALTITUDE - depth, 1.0).items()}
        beyond = mc.integral(meter, K_MAX_GEV, k_max=1e7, **state)
        assert beyond / mc.integral(meter, 1.0, k_max=1e7, n=3201, **state) == pytest.approx(share, rel=0.05)


def test_a_grid_too_coarse_is_refused_by_the_check():
    """Two elevation nodes 30 degrees apart: the comparison with Mulder under a slab fails the table."""
    with pytest.raises(RuntimeError, match="differs from Mulder under a slab"):
        MulderTable(mc.ALTITUDE, LOAM, cache=False, sin_elevation=(0.5, 1.0))
    unchecked = MulderTable(mc.ALTITUDE, LOAM, cache=False, sin_elevation=(0.5, 1.0), check_points=0)
    assert unchecked.meta["check"] == {"points": 0}


# ------------------------------------------------ below sea level


def test_below_sea_level_the_flux_is_lower_and_has_no_cutoff():
    """A site at -430 m: 56.69 against 60.50 m^-2 s^-1 sr^-1 above 1 GeV of total energy."""
    low, sea = band(mc.BELOW, 1.0), band(0.0, 1.0)
    frozen = FROZEN["below_sea_level/-430m/elev90"]
    assert low.cutoff(1.0, 0.0) == 0.0 == frozen["cutoff_gev"]
    assert not np.any(low.cutoff_gev)
    assert low.energy_min(1.0, 0.0) == 1.0e-3
    assert low.differential(1.0e-3, 1.0, 0.0) > 0
    for response in (ABOVE_1_GEV_TOTAL, EnergyResponse.threshold(0.2), EnergyResponse.threshold(1.0)):
        assert low.integral(response, 1.0, 0.0) < sea.integral(response, 1.0, 0.0)
    assert low.integral(ABOVE_1_GEV_TOTAL, 1.0, 0.0) == pytest.approx(56.69, abs=0.01)
    assert low.integral(ABOVE_1_GEV_TOTAL, 1.0, 0.0) == pytest.approx(frozen["E>1GeV"], rel=RTOL)
    for k_min in mc.THRESHOLDS:
        value = low.integral(EnergyResponse.threshold(k_min), 1.0, 0.0)
        assert value == pytest.approx(frozen[f"K>{k_min:g}GeV"], rel=RTOL)


# ------------------------------------------------ atmosphere, date, reference


def test_atmosphere_model():
    """Winter air at 1150 m: the frozen case `season`; the cut-off rises from 0.283 to 0.296 GeV."""
    assert "us-standard" in atmosphere_models() and "midlatitude-winter" in atmosphere_models()
    frozen = FROZEN["season/midlatitude-winter"]
    winter = band(mc.ALTITUDE + 0.01, 1.0, atmosphere="midlatitude-winter")
    assert winter.integral(ABOVE_1_GEV_TOTAL, 1.0, 0.0) == pytest.approx(frozen["open_sky/E>1GeV"], rel=RTOL)
    assert winter.cutoff(1.0, 0.0) == pytest.approx(frozen["open_sky/cutoff_gev"], rel=2e-5)
    assert winter.cutoff(1.0, 0.0) > band(mc.ALTITUDE + 0.01, 1.0).cutoff(1.0, 0.0)
    assert winter.meta["inputs"]["atmosphere"] == "midlatitude-winter"
    with pytest.raises(ValueError, match="unknown atmosphere 'martian'; Mulder has"):
        MulderTable(mc.ALTITUDE, LOAM, cache=False, atmosphere="martian")


def test_geomagnetic_field_of_a_date():
    """East against west at 45 degrees and 1150 m, IGRF: +0.27 % above 1 GeV, the frozen case `geomagnet`."""
    frozen = FROZEN["geomagnet/1150m/elev45"]
    value = {}
    for side, azimuth in (("east", 90.0), ("west", 270.0)):
        table = band(mc.ALTITUDE + 0.01, SIN45, date=mc.DATE, azimuth_deg=azimuth)
        value[side] = table.integral(ABOVE_1_GEV_TOTAL, SIN45, 0.0)
        assert value[side] == pytest.approx(frozen[f"{side}/E>1GeV"], rel=RTOL)
        assert table.meta["inputs"]["date"] == mc.DATE
    assert value["east"] / value["west"] - 1 == pytest.approx(0.0027, abs=0.0001)
    without = band(mc.ALTITUDE + 0.01, SIN45).integral(ABOVE_1_GEV_TOTAL, SIN45, 0.0)
    assert value["west"] < without < value["east"]
    for options in (
        {"date": mc.DATE},  # no azimuth
        {"date": mc.DATE, "azimuth_deg": 90.0, "latitude_deg": None},
        {"azimuth_deg": 90.0},  # an azimuth without a field
        {"date": "15.10.2028", "azimuth_deg": 90.0},
    ):
        with pytest.raises(ValueError):
            MulderTable(mc.ALTITUDE, LOAM, cache=False, **(SITE | options))


def test_latitude_matters_little_without_a_field():
    """The shape of the Earth alone: 1.6e-6 between the equator and 43 degrees north."""
    here = band(mc.ALTITUDE, SIN45).integral(ABOVE_1_GEV_TOTAL, SIN45, 0.0)
    equator = MulderTable(mc.ALTITUDE, LOAM, cache=False, sin_elevation=nodes(SIN45))
    assert here / equator.integral(ABOVE_1_GEV_TOTAL, SIN45, 0.0) - 1 == pytest.approx(1.6e-6, abs=1e-6)


def test_reference_follows_the_conventions():
    """Version 0.1.0 has the sea-level reference only; a reference at the site altitude is roadmap 0.1-05b."""
    auto = MulderTable(mc.ALTITUDE, LOAM, cache=False, sin_elevation=nodes(1.0))
    named = MulderTable(mc.ALTITUDE, LOAM, cache=False, sin_elevation=nodes(1.0), reference="sea-level")
    assert auto.key == named.key
    assert auto.meta["inputs"]["reference"] == "GCCLY15"
    assert auto.meta["inputs"]["mulder"] == mulder.config.VERSION
    for reference in ("site-altitude", auto, np.ones((3, 3))):  # by name, or as a table of some kind
        with pytest.raises(NotImplementedError, match="0.1-05b"):
            MulderTable(mc.ALTITUDE, LOAM, cache=False, reference=reference)
    with pytest.raises(ValueError, match="reference must be"):
        MulderTable(mc.ALTITUDE, LOAM, cache=False, reference="GCCLY15 at the site")


def test_cutoff_search_is_the_bisection_of_the_frozen_cases():
    cutoff = open_sky_cutoff(
        mulder.Fluxmeter(), mc.ALTITUDE, [1.0, 0.5], latitude_deg=43.35, longitude_deg=79.08
    )
    # the same search to the last bit on one machine; the frozen numbers come from another
    np.testing.assert_allclose(
        cutoff, [FROZEN["cutoff/1150m/elev90"], FROZEN["cutoff/1150m/elev30"]], rtol=mc.RTOL
    )
    np.testing.assert_array_equal(open_sky_cutoff(mulder.Fluxmeter(), -10.0, [1.0, 0.5]), [0.0, 0.0])


@pytest.mark.parametrize(
    "arguments",
    [
        {"material": "loam"},
        {"material": M.loam(1.6, 15.0)},  # a soil is not a material: its material is
        {"altitude_m": float("nan")},
        {"altitude_m": 2.0e4},
        {"sin_elevation": (1.0, 0.5)},
        {"sin_elevation": (0.5,)},
        {"sin_elevation": (0.0, 0.5)},
        {"nodes_per_decade": 5},
        {"check_points": -1},
        {"latitude_deg": 100.0},
    ],
)
def test_arguments_are_checked(arguments):
    arguments = {"altitude_m": mc.ALTITUDE, "material": LOAM} | arguments
    with pytest.raises((TypeError, ValueError)):
        MulderTable(arguments.pop("altitude_m"), arguments.pop("material"), cache=False, **arguments)


def test_example_of_the_note():
    """The example of `docs/notes/flux.md`: 20 m below a brow at 1150 m, thresholds of 0.2 and 1 GeV."""
    soil = M.loam(1.60, 15.0)  # dry density in g/cm3, water content in %
    table = MulderTable(1130.0, soil.material, cache=False, sin_elevation=nodes(1.0))
    opacity = 20.0 * soil.bulk_density_g_cm3 * 100  # g/cm2 by the vertical
    assert opacity == pytest.approx(3680.0)
    above_1_gev = table.integral(EnergyResponse.threshold(1.0), 1.0, opacity)
    above_02_gev = table.integral(EnergyResponse.threshold(0.2), 1.0, opacity)
    assert above_1_gev == pytest.approx(10.84, abs=0.01)
    assert above_02_gev == pytest.approx(12.73, abs=0.01)
    with pytest.raises(ValueError, match="cut-off of 0.2794 GeV"):
        table.integral(EnergyResponse.threshold(0.2), 1.0, 0.0)


def test_a_soil_of_the_library_as_reference_material():
    """Loam of 1.6 g/cm3 with 15 % of water: the table passes its comparison with Mulder under that soil."""
    soil = M.loam(1.60, 15.0)
    table = MulderTable(mc.ALTITUDE, soil.material, cache=False, sin_elevation=nodes(1.0), **SITE)
    assert table.meta["check"]["max"] < 2.0e-3
    assert table.meta["material_name"] == soil.material.name
    # in the open sky the material does not matter (5e-8: the energy grid ends where it needs to);
    # under soil it does, through its range
    assert table.integral(ABOVE_1_GEV_TOTAL, 1.0, 0.0) == pytest.approx(
        band(mc.ALTITUDE, 1.0).integral(ABOVE_1_GEV_TOTAL, 1.0, 0.0), rel=1e-6
    )
    under = table.integral(EnergyResponse.threshold(1.0), 1.0, 3600.0)
    assert under != pytest.approx(
        band(mc.ALTITUDE, 1.0).integral(EnergyResponse.threshold(1.0), 1.0, 3600.0), rel=1e-3
    )


# ------------------------------------------------ disk cache


def test_cache_by_hash_of_the_inputs(tmp_path):
    options = {"sin_elevation": nodes(1.0), "cache_dir": tmp_path} | SITE
    first = MulderTable(mc.ALTITUDE, LOAM, **options)
    assert not first.from_cache
    assert first.cache_path == tmp_path / f"{first.key[:24]}.npz"
    assert [p.name for p in tmp_path.iterdir()] == [first.cache_path.name]
    second = MulderTable(mc.ALTITUDE, LOAM, **options)
    assert second.from_cache and second.key == first.key
    np.testing.assert_array_equal(second.flux_m2_s_sr_gev, first.flux_m2_s_sr_gev)
    np.testing.assert_array_equal(second.cutoff_gev, first.cutoff_gev)
    assert second.meta == first.meta and second.meta["check"]["points"] == 64
    assert second.integral(ABOVE_1_GEV_TOTAL, 1.0, 90.0) == first.integral(ABOVE_1_GEV_TOTAL, 1.0, 90.0)
    assert isinstance(FluxTable.load(first.cache_path), FluxTable)  # a cached table is a plain file

    # every input is part of the key
    changed = [
        MulderTable(mc.ALTITUDE + 1.0, LOAM, **options),
        MulderTable(mc.ALTITUDE, M.STANDARD_ROCK, **options),
        MulderTable(mc.ALTITUDE, LOAM, **(options | {"atmosphere": "midlatitude-summer"})),
        MulderTable(mc.ALTITUDE, LOAM, **(options | {"date": mc.DATE, "azimuth_deg": 0.0})),
        MulderTable(mc.ALTITUDE, LOAM, **(options | {"latitude_deg": 44.0})),
        MulderTable(mc.ALTITUDE, LOAM, **(options | {"sin_elevation": nodes(0.5)})),
        MulderTable(mc.ALTITUDE, LOAM, **(options | {"nodes_per_decade": 50})),
        MulderTable(mc.ALTITUDE, LOAM, **(options | {"check_points": 16})),
    ]
    keys = {first.key} | {table.key for table in changed}
    assert len(keys) == 9 and not any(table.from_cache for table in changed)
    assert len(list(tmp_path.iterdir())) == 9


def test_a_damaged_cache_file_is_rebuilt(tmp_path):
    options = {"sin_elevation": nodes(1.0), "cache_dir": tmp_path}
    first = MulderTable(mc.ALTITUDE, LOAM, **options)
    first.cache_path.write_bytes(b"not a table")
    again = MulderTable(mc.ALTITUDE, LOAM, **options)
    assert not again.from_cache
    assert MulderTable(mc.ALTITUDE, LOAM, **options).from_cache
    # a file made from other inputs under the same name is not trusted either
    other = MulderTable(mc.ALTITUDE + 5.0, LOAM, **options)
    first.cache_path.write_bytes(other.cache_path.read_bytes())
    assert not MulderTable(mc.ALTITUDE, LOAM, **options).from_cache
    # nor one that reads but does not make a table
    with np.load(first.cache_path, allow_pickle=False) as data:
        stored = {name: data[name] for name in data.files}
    np.savez(first.cache_path, **(stored | {"cutoff_gev": -stored["cutoff_gev"]}))
    assert not MulderTable(mc.ALTITUDE, LOAM, **options).from_cache
    assert MulderTable(mc.ALTITUDE, LOAM, **options).from_cache


def test_without_a_cache_nothing_is_written(tmp_path, monkeypatch):
    monkeypatch.setenv("MUOSCENE_CACHE", str(tmp_path / "cache"))
    table = MulderTable(mc.ALTITUDE, LOAM, cache=False, sin_elevation=nodes(1.0))
    assert table.cache_path is None and not table.from_cache
    assert list(tmp_path.iterdir()) == []
    # the default directory is the user cache, never the working directory
    cached = MulderTable(mc.ALTITUDE, LOAM, sin_elevation=nodes(1.0))
    assert cached.cache_path.parent == tmp_path / "cache" / "flux"
    assert cached.cache_path.is_file()


def test_a_table_on_disk_is_under_a_megabyte(site, tmp_path):
    site.save(tmp_path / "site.npz")
    assert 0.7e6 < (tmp_path / "site.npz").stat().st_size < 0.9e6  # 151 elevations by 607 energies
    again = FluxTable.load(tmp_path / "site.npz")
    assert again.meta["check"] == site.meta["check"]
    assert again.integral(ABOVE_1_GEV_TOTAL, 0.8, 500.0) == site.integral(ABOVE_1_GEV_TOTAL, 0.8, 500.0)


def test_default_cache_directory(tmp_path, monkeypatch):
    monkeypatch.setenv("MUOSCENE_CACHE", str(tmp_path / "explicit"))
    assert default_cache_dir() == tmp_path / "explicit"
    monkeypatch.delenv("MUOSCENE_CACHE")
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "xdg"))
    assert default_cache_dir() == tmp_path / "xdg" / "muoscene"
    monkeypatch.delenv("XDG_CACHE_HOME")
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    monkeypatch.setenv("USERPROFILE", str(tmp_path / "home"))
    assert default_cache_dir() == tmp_path / "home" / ".cache" / "muoscene"
