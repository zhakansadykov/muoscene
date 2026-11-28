"""`muoscene.viz.curves`: the flux of a table against the zenith angle, and where it has no number.

Roadmap 0.1-12a. The tables are those of `test_flux_base.py`, analytic, so that no Mulder is
needed: without a cut-off the integral above K has a closed form; with a cut-off C / sin(elevation)
a threshold K has a number while C / s plus the margin of the grid (10 MeV) is at most K, that is
up to the zenith angle arccos(C / (K − 0.01)), and the table ends at sin(elevation) = 0.25, at
75.52°. The numbers of a table made by Mulder are checked in `test_example_open_sky.py`.
"""

import math

import matplotlib.text
import numpy as np
import pytest

from muoscene.flux import EnergyResponse
from muoscene.transport import BELOW_CUTOFF, DEFINED, LOW_ELEVATION
from muoscene.viz import curves, style
from tests.unit.test_flux_base import MARGIN, RTOL, TOP, C, above, table

ONE, LOW = EnergyResponse.threshold(1.0), EnergyResponse.threshold(0.2)
TABLE_END_DEG = math.degrees(math.acos(0.25))


@pytest.fixture(scope="module")
def plain():
    return table()


@pytest.fixture(scope="module")
def cut():
    return table(cutoff=True)


def test_the_curve_is_the_integral_of_the_table(plain):
    """Without a cut-off: the closed form up to the end of the table, then no number."""
    curve = curves.flux_curve(plain, ONE)
    s = np.cos(np.radians(curve.zenith_deg))
    inside = curve.zenith_deg <= TABLE_END_DEG
    np.testing.assert_array_equal(curve.status == DEFINED, inside)
    assert set(curve.status[~inside].tolist()) == {LOW_ELEVATION}
    np.testing.assert_allclose(curve.flux_m2_s_sr.compressed(), above(1.0, s[inside], 0.0, TOP), rtol=RTOL)
    assert curve.edge_deg == pytest.approx(TABLE_END_DEG, abs=1e-5)
    under = curves.flux_curve(plain, ONE, [0.0, 30.0], opacity_g_cm2=500.0)
    np.testing.assert_allclose(
        under.flux_m2_s_sr, above(1.0, np.cos(np.radians([0.0, 30.0])), 500.0, TOP), rtol=RTOL
    )


def test_where_the_cut_off_ends_the_curve(cut):
    """With a cut-off: 1 GeV up to arccos(C / (1 − 0.01)); 0.2 GeV nowhere, as at 1150 m."""
    edge = math.degrees(math.acos(C / (1.0 - MARGIN)))
    curve = curves.flux_curve(cut, ONE)
    assert curve.edge_deg == pytest.approx(edge, abs=1e-5)
    assert np.all(curve.status[curve.zenith_deg <= edge] == DEFINED)
    beyond = curve.status[(curve.zenith_deg > edge) & (curve.zenith_deg <= TABLE_END_DEG)]
    assert beyond.size and np.all(beyond == BELOW_CUTOFF)
    assert np.all(curve.status[curve.zenith_deg > TABLE_END_DEG] == LOW_ELEVATION)
    np.testing.assert_allclose(
        curve.flux_m2_s_sr.compressed(),
        cut.integral(ONE, np.cos(np.radians(curve.zenith_deg[curve.status == 0])), 0.0),
    )
    nothing = curves.flux_curve(cut, LOW)
    assert nothing.edge_deg is None and nothing.last_defined_deg is None
    assert np.all(np.ma.getmaskarray(nothing.flux_m2_s_sr))


def test_the_figure_draws_the_curves_and_says_where_they_stop(cut):
    others = [("a formula", [0.0, 30.0, 60.0], [70.0, 50.0, 20.0])]
    figure = curves.flux_figure(cut, [ONE, LOW], others=others)
    (ax,) = figure.axes
    lines = {line.get_gid(): line for line in ax.get_lines()}
    curve = curves.flux_curve(cut, ONE)
    np.testing.assert_array_equal(lines["muoscene-curve-0"].get_ydata(), curve.flux_m2_s_sr.compressed())
    np.testing.assert_array_equal(lines["muoscene-curve-1"].get_ydata(), [])
    assert ax.get_yscale() == "log"
    text = " ".join(t.get_text() for t in figure.findobj(matplotlib.text.Text)).replace("\n", " ")
    assert f"K > 1 GeV: a number up to {curve.edge_deg:.2f}°" in text
    assert "K > 0.2 GeV: no number at any angle (below the cut-off" in text
    assert "a formula" in text and "table at 1000 m" in text and "muoscene" in text
    shaded = {p.get_gid() for p in ax.patches}
    assert {f"muoscene-reason-{BELOW_CUTOFF}", f"muoscene-reason-{LOW_ELEVATION}"} <= shaded


def test_wrong_input_is_refused(cut, tmp_path):
    with pytest.raises(ValueError, match="up to, not including, 90"):
        curves.flux_curve(cut, ONE, [0.0, 90.0])
    with pytest.raises(ValueError, match="non-negative"):
        curves.flux_curve(cut, ONE, opacity_g_cm2=-1.0)
    with pytest.raises(TypeError, match="EnergyResponse"):
        curves.flux_figure(cut, [1.0])
    with pytest.raises(ValueError, match="same length"):
        curves.flux_figure(cut, ONE, others=[("x", [0.0, 1.0], [1.0])])
    path = style.save(curves.flux_figure(cut, ONE, annotate=False), tmp_path / "flux.pdf")
    assert path.stat().st_size > 5_000
