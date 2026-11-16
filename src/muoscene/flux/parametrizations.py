"""flux.parametrizations — `LegacyTable`: Guan et al. (2015) at sea level with a + bE, for regression tests.

The flux model of the preliminary calculation behind the 2027 grant application, with the same
arithmetic as `tests/regression/legacy.py` (`docs/legacy.md`): the spectrum of Guan et al. (2015)
at sea level in total energy, tabulated as the integral intensity above E on 301 values of
cos(zenith) from 0.25 (nearest node) and 2001 values of ln E from 0.1 GeV to 10 TeV (linear
interpolation of the logarithm), and the total energy needed to cross an opacity X,
E(X) = a/b (exp(bX) - 1) + offset. It has no altitude, no cut-off and no detector energy
response: it counts every muon above E(X) at the surface. It is a reference for parity tests,
not a model offered to users.

Units: opacity in g/cm², total energy in GeV, integral flux in m⁻² s⁻¹ sr⁻¹ (the preliminary
calculation works in cm⁻²).
"""

from __future__ import annotations

import functools
from dataclasses import dataclass

import numpy as np

from .base import EDGE, OPACITY_MAX_G_CM2

P1, P2, P3, P4, P5 = 0.102573, -0.068287, 0.958633, 0.0407253, 0.817285
COS_GRID = np.linspace(0.25, 1.0, 301)
E_GRID = np.logspace(-1, 4, 2001)
LN_E = np.log(E_GRID)


def cos_star(c):
    """Effective cosine of the zenith angle of Guan et al. (2015): it accounts for the curved atmosphere."""
    return np.sqrt((c**2 + P1**2 + P2 * c**P3 + P4 * c**P5) / (1 + P1**2 + P2 + P4))


def guan_intensity(total_energy_gev, cos_zenith):
    """Differential intensity of Guan et al. (2015) at sea level, cm⁻² s⁻¹ sr⁻¹ GeV⁻¹; total energy in GeV."""
    E = total_energy_gev
    cs = cos_star(cos_zenith)
    Es = E * (1 + 3.64 / (E * cs**1.29))
    return 0.14 * Es**-2.7 * (1 / (1 + 1.1 * E * cs / 115) + 0.054 / (1 + 1.1 * E * cs / 850))


@functools.cache
def _ln_integral() -> np.ndarray:
    """Logarithm of the integral intensity above each energy of `E_GRID`, per node of `COS_GRID`."""
    table = np.empty((COS_GRID.size, E_GRID.size))
    for i, c in enumerate(COS_GRID):
        f = guan_intensity(E_GRID, c) * E_GRID  # integral over ln E
        cum = np.concatenate([np.cumsum(((f[1:] + f[:-1]) / 2 * np.diff(LN_E))[::-1])[::-1], [0.0]])
        table[i] = np.log(np.maximum(cum, 1e-30))
    table.flags.writeable = False
    return table


class _LegacyResponse:
    """The energy response of the preliminary calculation: none, every muon it transports is counted."""

    __slots__ = ()

    def __repr__(self) -> str:
        return "LEGACY_RESPONSE"


LEGACY_RESPONSE = _LegacyResponse()
"""What `LegacyTable.integral` takes as its response: the preliminary calculation has no other."""


@dataclass(frozen=True)
class LegacyTable:
    """Integral flux of the preliminary calculation against elevation and opacity.

    a_loss_gev_cm2_g   a of the energy loss a + bE, GeV per g/cm²
    b_loss_cm2_g       b of the energy loss, per g/cm²
    offset_gev         total energy added to the threshold, GeV; 0.1 GeV, below the muon mass
    opacity_max_g_cm2  largest opacity served, g/cm²

    The defaults are the values of the preliminary calculation (`legacy.Model`).
    """

    a_loss_gev_cm2_g: float = 2.2e-3
    b_loss_cm2_g: float = 4.0e-6
    offset_gev: float = 0.1
    opacity_max_g_cm2: float = OPACITY_MAX_G_CM2

    def __post_init__(self) -> None:
        if not (self.a_loss_gev_cm2_g > 0 and self.b_loss_cm2_g > 0 and self.opacity_max_g_cm2 > 0):
            raise ValueError("a_loss_gev_cm2_g, b_loss_cm2_g and opacity_max_g_cm2 must be positive")
        if not self.offset_gev >= E_GRID[0]:
            raise ValueError(
                f"offset_gev must be at least {E_GRID[0]:g} GeV, where the table of the spectrum starts"
            )
        if self.threshold_gev(self.opacity_max_g_cm2) > E_GRID[-1]:
            raise ValueError(
                "the threshold under opacity_max_g_cm2 is beyond the table of the spectrum,"
                f" {E_GRID[-1]:g} GeV"
            )

    def threshold_gev(self, opacity_g_cm2):
        """Total energy (GeV) a muon needs at the surface to be counted under the opacity (g/cm²)."""
        x = np.asarray(opacity_g_cm2, dtype=float)
        value = self.a_loss_gev_cm2_g / self.b_loss_cm2_g * np.expm1(self.b_loss_cm2_g * x) + self.offset_gev
        return float(value) if x.ndim == 0 else value

    def integral(self, response, sin_elevation, opacity_g_cm2):
        """Integral flux under the opacity, m⁻² s⁻¹ sr⁻¹: the intensity above `threshold_gev`.

        `response` must be `LEGACY_RESPONSE`: the preliminary calculation has no detector energy
        response, and a threshold given here could not be honoured. The sine of elevation (the
        cosine of the zenith angle) must lie in [0.25, 1], the opacity in [0, opacity_max_g_cm2].
        """
        if response is not LEGACY_RESPONSE:
            raise ValueError(
                "the preliminary calculation has no detector energy response: pass LEGACY_RESPONSE,"
                " or use a MulderTable for a response in kinetic energy"
            )
        c, x = np.broadcast_arrays(
            np.asarray(sin_elevation, dtype=float), np.asarray(opacity_g_cm2, dtype=float)
        )
        bad = ~((c >= COS_GRID[0] * (1 - EDGE)) & (c <= COS_GRID[-1] * (1 + EDGE)))
        if np.any(bad):
            raise ValueError(
                f"sin elevation {c[bad].flat[0]:g} is outside the table, {COS_GRID[0]:g} to {COS_GRID[-1]:g}"
            )
        bad = ~((x >= 0) & (x <= self.opacity_max_g_cm2 * (1 + EDGE)))
        if np.any(bad):
            raise ValueError(
                f"opacity {x[bad].flat[0]:g} g/cm2 is outside the table,"
                f" 0 to {self.opacity_max_g_cm2:g} g/cm2"
            )
        e_th = self.a_loss_gev_cm2_g / self.b_loss_cm2_g * np.expm1(self.b_loss_cm2_g * x) + self.offset_gev
        table = _ln_integral()
        # as in the preliminary calculation: nearest node in cos(zenith), linear in ln E
        ci = np.clip(
            np.rint((c - COS_GRID[0]) / (COS_GRID[1] - COS_GRID[0])).astype(int), 0, COS_GRID.size - 1
        )
        u = np.clip(np.log(e_th), LN_E[0], LN_E[-1])
        j = np.clip(np.searchsorted(LN_E, u) - 1, 0, LN_E.size - 2)
        t = (u - LN_E[j]) / (LN_E[j + 1] - LN_E[j])
        value = np.exp(table[ci, j] * (1 - t) + table[ci, j + 1] * t) * 1e4
        scalar = np.ndim(sin_elevation) == 0 and np.ndim(opacity_g_cm2) == 0
        return float(value) if scalar else value
