"""Discrete transport on the reference void: how much of its signal stays in the geometric cone.

The frozen numbers of `tests/regression` are continuous transport; of the discrete mode they
hold only the tables it draws from. This test runs the Monte Carlo itself, with scattering and
straggling, on the scene of `mulder_cases.reference_void` at 1150 m: the excess of counts in the
cone of the void, discrete over continuous. With muons above 0.1 GeV at the detector the share
is 0.74 +- 0.02 (`docs/numbers.md`); 20 000 states give it to +- 0.05, so the bounds below are
three standard deviations wide. Without scattering the share would be 1.

Run with `pytest -m validation` (about 40 s).
"""

import mulder
import numpy as np
import pytest

from tests.regression import mulder_cases as mc

pytestmark = pytest.mark.validation

STATES = 20000
K_MIN, K_MAX = 0.1, 1e3  # kinetic energy at the detector, GeV
SEED = 11


def excess(mode):
    """Relative excess of counts in the cone of the void and its standard error."""
    geometry = mc.build_scene(mc.VOID_DEPTH, altitude=mc.ALTITUDE, sphere=mc.REFERENCE_VOID)
    index = mc.medium_index(geometry, "Anomaly")
    d, centre, radius = mc.detector(mc.VOID_DEPTH), np.array(mc.REFERENCE_VOID[:3]), mc.REFERENCE_VOID[3]
    axis = (centre - d) / np.linalg.norm(centre - d)
    alpha = np.arcsin(radius / np.linalg.norm(centre - d))
    e1 = np.cross(axis, [1.0, 0.0, 0.0])
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(axis, e1)
    rng = np.random.default_rng(SEED)
    energy = np.exp(rng.uniform(np.log(K_MIN), np.log(K_MAX), STATES))
    cb = rng.uniform(np.cos(alpha), 1.0, STATES)
    sb, ph = np.sqrt(1 - cb**2), rng.uniform(0, 2 * np.pi, STATES)
    u = cb[:, None] * axis + (sb * np.cos(ph))[:, None] * e1 + (sb * np.sin(ph))[:, None] * e2
    position = np.broadcast_to(d, u.shape)
    values = []
    for density in (mc.SOIL, mc.AIR):
        geometry.media[index].density = density
        meter = mulder.Fluxmeter(geometry=geometry, mode=mode, seed=SEED)
        extra = {"events": 1} if mode == "discrete" else {}
        flux = np.asarray(meter.flux(position=position, direction=u, energy=energy, **extra))
        values.append((flux[:, 0] if flux.ndim == 2 else flux) * energy * u[:, 2])
    base, void = values
    ratio = void.sum() / base.sum() - 1
    error = np.sqrt(STATES) * np.std((void - base) - ratio * base) / base.sum()
    return ratio, error


def test_scattering_moves_a_quarter_of_the_signal_out_of_the_cone():
    continuous, _ = excess("continuous")
    discrete, error = excess("discrete")
    assert continuous == pytest.approx(0.064, abs=0.002)
    share = discrete / continuous
    assert error / continuous < 0.07
    assert 0.59 < share < 0.89
