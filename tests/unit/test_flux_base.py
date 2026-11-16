"""`FluxTable` and `EnergyResponse` on an analytic table: no Mulder here.

The table of these tests has a constant stopping power A, so that the range is K / A, the energy
above an opacity X is K + A X and the Jacobian is one; its open-sky flux is
exp(B s) (K + 1)^-G in the sine of elevation s. Everything the table returns has a closed form.
With a cut-off C / s the flux is exp(B s) d^-1/2 (1 + d)^-G in d = K - C / s, which has the
integrable singularity of a sea-level reference seen from altitude.
"""

import numpy as np
import pytest

from muoscene.flux import OPACITY_MAX_G_CM2, EnergyResponse, FluxTable

A = 2.0e-3  # stopping power, GeV per g/cm^2
B = 1.5  # slope of ln flux in the sine of elevation
G = 2.7  # spectral index
C = 0.28  # cut-off times sine of elevation, GeV
K_MIN, TOP = 1.0e-3, 2.0e4  # GeV
MARGIN = 0.01  # GeV, where the grid of a table with a cut-off starts above it
SINES = np.linspace(0.25, 1.0, 16)
# Interpolation is linear in logarithms on 100 nodes per decade: the flux of these tables is
# reproduced to 4.5e-5 and its integrals to 2.4e-5, so the tolerance below leaves a factor of two.
RTOL = 1.0e-4


def table(cutoff=False, **changes):
    first = MARGIN if cutoff else K_MIN
    d = np.logspace(np.log10(first), np.log10(TOP), int(np.log10(TOP / first) * 100) + 1)
    k = np.logspace(np.log10(K_MIN), np.log10(1.1 * TOP), 1500)
    shape = d**-0.5 * (1 + d) ** -G if cutoff else (d + 1) ** -G
    arguments = {
        "altitude_m": 1000.0 if cutoff else 0.0,
        "sin_elevation": SINES,
        "cutoff_gev": C / SINES if cutoff else np.zeros(SINES.size),
        "energy_above_cutoff_gev": d,
        "flux_m2_s_sr_gev": np.exp(B * SINES)[:, None] * shape,
        "range_energy_gev": k,
        "range_g_cm2": k / A,
        "stopping_power_gev_cm2_g": np.full(k.size, A),
        "meta": {"made by": "tests/unit/test_flux_base.py"},
    }
    return FluxTable(**(arguments | changes))


def above(k, s, x, top):
    """Closed form of the integral flux between the kinetic energy k and `top`, without a cut-off."""
    return np.exp(B * s) * ((k + A * x + 1) ** (1 - G) - (top + 1) ** (1 - G)) / (G - 1)


@pytest.fixture(scope="module")
def plain():
    return table()


@pytest.fixture(scope="module")
def cut():
    return table(cutoff=True)


# ------------------------------------------------ differential flux


def test_flux_under_an_opacity_is_the_open_sky_flux_at_the_energy_above_it(plain):
    rng = np.random.default_rng(1)
    s = rng.uniform(0.25, 1.0, 2000)
    x = np.where(rng.random(2000) < 0.2, 0.0, np.exp(rng.uniform(0.0, np.log(OPACITY_MAX_G_CM2), 2000)))
    k = np.exp(rng.uniform(np.log(K_MIN), np.log(1.0e4), 2000))
    expected = np.exp(B * s) * (k + A * x + 1) ** -G
    np.testing.assert_allclose(plain.differential(k, s, x), expected, rtol=RTOL)


def test_jacobian_is_the_ratio_of_stopping_powers(plain):
    """With S = A (1 + K) the range is ln(1 + K) / A: the flux gains (1 + K_top) / (1 + K)."""
    k = plain.range_energy_gev
    # such a muon crosses 4950 g/cm2 at most: the table is made for 2000 g/cm2
    curved = table(
        range_g_cm2=np.log1p(k) / A, stopping_power_gev_cm2_g=A * (1 + k), opacity_max_g_cm2=2000.0
    )
    energy, s, x = np.array([0.5, 3.0, 40.0]), 0.7, np.array([200.0, 1500.0, 900.0])
    top = np.expm1(np.log1p(energy) + A * x)
    expected = np.exp(B * s) * (top + 1) ** -G * (1 + top) / (1 + energy)
    np.testing.assert_allclose(curved.differential(energy, s, x), expected, rtol=RTOL)


def test_scalars_give_floats_and_arrays_broadcast(plain):
    assert isinstance(plain.differential(1.0, 0.5, 100.0), float)
    assert isinstance(plain.integral(EnergyResponse.threshold(1.0), 0.5, 100.0), float)
    assert isinstance(plain.cutoff(0.5, 0.0), float)
    grid = plain.differential(np.array([1.0, 2.0, 5.0])[:, None], np.array([0.3, 0.9]), 0.0)
    assert grid.shape == (3, 2)
    assert plain.integral(EnergyResponse.threshold(1.0), np.array([0.3, 0.9]), np.zeros((4, 1))).shape == (
        4,
        2,
    )


# ------------------------------------------------ integral over an energy response


def test_integral_starts_at_the_threshold_itself(plain):
    """Thresholds between the nodes of the grid: the closed form, not the sum from the next node."""
    s, x = np.array([1.0, 0.61, 0.25]), np.array([0.0, 730.0, 3.0e4])
    top = plain.energy_above_cutoff_gev[-1]
    for k_min in (K_MIN, 0.2, 1.0, 1.0123, 37.0):
        value = plain.integral(EnergyResponse.threshold(k_min), s, x)
        np.testing.assert_allclose(value, above(k_min, s, x, top), rtol=RTOL)
    # one step of the grid is 2.3 % in energy and 2.0 % in this integral: no threshold is rounded to a node
    step = plain.integral(EnergyResponse.threshold(1.0), 1.0, 0.0) / plain.integral(
        EnergyResponse.threshold(10**0.01), 1.0, 0.0
    )
    assert step - 1 == pytest.approx(above(1.0, 1.0, 0.0, top) / above(10**0.01, 1.0, 0.0, top) - 1, rel=1e-3)


def test_integral_over_a_response_is_the_integral_of_efficiency_times_flux(plain):
    response = EnergyResponse((0.5, 1.0, 3.0), (0.0, 0.6, 0.95))
    for s, x in ((1.0, 0.0), (0.4, 2500.0)):
        k = np.logspace(np.log10(0.5), np.log10(plain.k_max_gev), 20001)
        brute = np.trapezoid(response(k) * plain.differential(k, s, x) * k, np.log(k))
        beyond = 0.95 * (above(plain.k_max_gev, s, x, plain.energy_above_cutoff_gev[-1]))
        assert plain.integral(response, s, x) == pytest.approx(brute + beyond, rel=RTOL)


def test_constant_efficiency_scales_the_integral(plain):
    whole = plain.integral(EnergyResponse.threshold(1.0), 0.8, 500.0)
    assert plain.integral(EnergyResponse.threshold(1.0, efficiency=0.9), 0.8, 500.0) == pytest.approx(
        0.9 * whole
    )
    flat = EnergyResponse((1.0, 5.0), (0.9, 0.9))
    assert plain.integral(flat, 0.8, 500.0) == pytest.approx(0.9 * whole)


def test_the_response_is_never_assumed(plain):
    with pytest.raises(TypeError, match="EnergyResponse"):
        plain.integral(1.0, 0.8, 500.0)
    with pytest.raises(TypeError):
        plain.integral(None, 0.8, 500.0)
    with pytest.raises(TypeError):
        plain.integral(0.8, 500.0)  # the response has no default


def test_energy_response():
    step = EnergyResponse.threshold(1.0)
    assert step.k_min_gev == 1.0
    np.testing.assert_array_equal(step([0.5, 1.0, 7.0]), [0.0, 1.0, 1.0])
    ramp = EnergyResponse((0.2, 1.0), (0.0, 0.8))
    np.testing.assert_allclose(ramp([0.1, 0.2, 0.6, 1.0, 50.0]), [0.0, 0.0, 0.4, 0.8, 0.8])
    assert isinstance(ramp(0.6), float)
    assert EnergyResponse([0.2, 1.0], np.array([0.0, 0.8])) == ramp  # sequences become tuples
    for energies, efficiency in (
        ((), ()),
        ((1.0, 0.5), (0.5, 1.0)),  # not increasing
        ((0.0, 1.0), (0.5, 1.0)),  # zero energy
        ((1.0,), (1.2,)),  # efficiency above one
        ((1.0,), (-0.1,)),
        ((1.0, 2.0), (1.0,)),
        ((float("nan"),), (1.0,)),
    ):
        with pytest.raises(ValueError):
            EnergyResponse(energies, efficiency)


# ------------------------------------------------ cut-off


def test_cutoff_in_the_open_sky_and_under_an_opacity(cut):
    s = np.array([1.0, 0.5, 0.3275, 0.25])  # 0.3275 lies between nodes
    np.testing.assert_allclose(cut.cutoff(s, 0.0), C / s, rtol=1e-12)
    np.testing.assert_allclose(cut.energy_min(s, 0.0), C / s + MARGIN, rtol=1e-12)
    # an opacity X takes A X off the cut-off; where less than the first energy of the table is left, none
    np.testing.assert_allclose(
        cut.cutoff(1.0, [50.0, 100.0, 139.4, 139.6, 500.0]), [0.18, 0.08, 1.2e-3, 0, 0]
    )
    np.testing.assert_allclose(cut.energy_min(1.0, [50.0, 100.0, 500.0]), [0.19, 0.09, K_MIN], rtol=1e-9)
    assert table().cutoff(0.5, 0.0) == 0.0
    assert table().energy_min(0.5, 0.0) == K_MIN


def test_flux_above_the_cutoff(cut):
    rng = np.random.default_rng(2)
    s = rng.uniform(0.25, 1.0, 2000)
    x = np.where(rng.random(2000) < 0.3, 0.0, np.exp(rng.uniform(0.0, np.log(OPACITY_MAX_G_CM2), 2000)))
    low = cut.energy_min(s, x)
    k = low * np.exp(rng.uniform(0.0, 1.0, 2000) * np.log(1.0e4 / low))
    d = k + A * x - C / s
    np.testing.assert_allclose(cut.differential(k, s, x), np.exp(B * s) * d**-0.5 * (1 + d) ** -G, rtol=RTOL)
    # the integral from a threshold just above the cut-off, against a fine sum in ln(K - K_cut)
    s, x, k_min = 0.5, 0.0, C / 0.5 + 0.011
    d = np.logspace(np.log10(0.011), np.log10(cut.energy_above_cutoff_gev[-1]), 200001)
    fine = np.trapezoid(np.exp(B * s) * d**-0.5 * (1 + d) ** -G * d, np.log(d))
    assert cut.integral(EnergyResponse.threshold(k_min), s, x) == pytest.approx(fine, rel=RTOL)


def test_below_the_cutoff_is_refused(cut):
    """Not zero and not an extrapolation: the table holds no energy there."""
    with pytest.raises(ValueError, match="cut-off of 0.28 GeV"):
        cut.differential(0.2, 1.0, 0.0)
    with pytest.raises(ValueError, match="cut-off"):
        cut.differential(0.285, 1.0, 0.0)  # above the cut-off, inside the margin
    assert cut.differential(0.2901, 1.0, 0.0) > 0
    with pytest.raises(ValueError, match="detector threshold 0.2 GeV is below 0.29 GeV"):
        cut.integral(EnergyResponse.threshold(0.2), 1.0, 0.0)
    with pytest.raises(ValueError, match="detector threshold"):
        cut.integral(EnergyResponse.threshold(1.0), [1.0, 0.5, 0.25], 0.0)  # 1.12 GeV at sin 0.25
    assert np.all(cut.integral(EnergyResponse.threshold(1.0), [1.0, 0.5, 0.3], 0.0) > 0)
    # under 100 g/cm2 the cut-off is 0.08 GeV: 0.2 GeV holds there and not in the open sky
    assert cut.integral(EnergyResponse.threshold(0.2), 1.0, 100.0) > 0
    with pytest.raises(ValueError):
        cut.integral(EnergyResponse((0.2, 1.0), (0.0, 1.0)), 1.0, 0.0)  # a response that only starts below


# ------------------------------------------------ outside the grid


@pytest.mark.parametrize(
    ("k", "s", "x", "message"),
    [
        (1.0, 0.5, 3.0001e4, "opacity"),
        (1.0, 0.5, -1.0, "opacity"),
        (1.0, 0.5, float("nan"), "opacity"),
        (1.0, 0.2, 0.0, "sin elevation"),
        (1.0, 1.01, 0.0, "sin elevation"),
        (1.0, float("nan"), 0.0, "sin elevation"),
        (5.0e-4, 0.5, 0.0, "below"),
        (0.0, 0.5, 100.0, "below"),
        (-1.0, 0.5, 100.0, "below"),
        (1.999e4, 0.5, 0.0, "above"),  # in the grid, but 3e4 g/cm2 cost 60 GeV: served up to 19 940 GeV
        (float("nan"), 0.5, 0.0, "finite"),
        (float("inf"), 0.5, 0.0, "finite"),
    ],
)
def test_outside_the_grid_is_refused(plain, k, s, x, message):
    with pytest.raises(ValueError, match=message):
        plain.differential(k, s, x)
    if message in ("opacity", "sin elevation"):
        with pytest.raises(ValueError, match=message):
            plain.integral(EnergyResponse.threshold(1.0), s, x)
        with pytest.raises(ValueError, match=message):
            plain.cutoff(s, x)


def test_one_bad_value_refuses_the_whole_request(plain):
    with pytest.raises(ValueError, match="opacity 35000 g/cm2"):
        plain.integral(EnergyResponse.threshold(1.0), 0.5, np.array([10.0, 3.5e4, 20.0]))
    with pytest.raises(ValueError, match="detector threshold 2e\\+04 GeV is above"):
        plain.integral(EnergyResponse.threshold(2.0e4), 0.5, 0.0)
    assert plain.differential(plain.k_max_gev, 0.5, OPACITY_MAX_G_CM2) > 0  # the corner of the domain


def test_edges_of_the_grid_are_inside(plain):
    assert plain.differential(K_MIN, 0.25, 0.0) > 0
    assert plain.differential(1.0, 1.0, OPACITY_MAX_G_CM2) > 0
    assert plain.k_min_gev == K_MIN
    assert 1.0e4 < plain.k_max_gev < TOP  # what is left of the top of the grid under 3e4 g/cm2


# ------------------------------------------------ construction and files


@pytest.mark.parametrize(
    "changes",
    [
        {"sin_elevation": SINES[::-1]},
        {"sin_elevation": np.linspace(0.0, 1.0, 16)},
        {"cutoff_gev": np.r_[np.zeros(15), 0.3]},  # a cut-off at one elevation only
        {"cutoff_gev": np.zeros(3)},
        {"flux_m2_s_sr_gev": np.zeros((16, 731))},
        {"flux_m2_s_sr_gev": np.ones((16, 7))},
        {"range_g_cm2": np.ones(1500)},
        {"stopping_power_gev_cm2_g": np.zeros(1500)},
        {"range_energy_gev": np.logspace(-3, 3, 1500)},  # ends below the flux grid
        {"range_energy_gev": np.logspace(-2, 4.5, 1500)},  # starts above the flux grid
        {"opacity_max_g_cm2": 0.0},
        {"opacity_max_g_cm2": 2.0e7},  # more than the top of the grid can cross
        {"altitude_m": float("nan")},
    ],
)
def test_inconsistent_tables_are_refused(changes):
    with pytest.raises(ValueError):
        table(**changes)


def test_arrays_cannot_be_changed(plain):
    with pytest.raises(ValueError):
        plain.flux_m2_s_sr_gev[0, 0] = 1.0
    with pytest.raises(ValueError):
        plain.cutoff_gev[0] = 1.0


def test_save_and_load(tmp_path, cut):
    path = tmp_path / "table.npz"
    cut.save(path)
    assert [p.name for p in tmp_path.iterdir()] == ["table.npz"]  # no temporary file left
    again = FluxTable.load(path)
    assert again.meta == cut.meta == {"made by": "tests/unit/test_flux_base.py"}
    assert again.altitude_m == cut.altitude_m and again.opacity_max_g_cm2 == cut.opacity_max_g_cm2
    for name in ("sin_elevation", "cutoff_gev", "energy_above_cutoff_gev", "flux_m2_s_sr_gev", "range_g_cm2"):
        np.testing.assert_array_equal(getattr(again, name), getattr(cut, name))
    assert again.integral(EnergyResponse.threshold(1.0), 0.7, 300.0) == cut.integral(
        EnergyResponse.threshold(1.0), 0.7, 300.0
    )
    with np.load(path, allow_pickle=False) as data:  # plain arrays: readable without pickle
        assert int(data["format"]) == 1
        stored = {name: data[name] for name in data.files}
    np.savez(tmp_path / "newer.npz", **(stored | {"format": np.array(99)}))
    with pytest.raises(ValueError, match="file layout 99"):
        FluxTable.load(tmp_path / "newer.npz")
