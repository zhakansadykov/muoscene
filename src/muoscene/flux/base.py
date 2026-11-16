"""flux.base — `FluxTable`: grid, interpolation, integral over a detector energy response.

In continuous transport the energy a muon loses is a function of the opacity it crosses, so the
flux under an opacity X of one material follows from the open-sky flux at the same place:

    phi(K, X) = phi_sky(K_top) S(K_top) / S(K),    R(K_top) = R(K) + X,

with R the CSDA range and S the stopping power of the material. The muons above K under X are
the muons above K_top in the open sky, so an integral above a threshold needs no Jacobian.

A `FluxTable` therefore holds the open-sky differential flux at the detector altitude on a grid
of kinetic energy and sine of elevation, and the range and stopping power of its reference
material; the opacity enters through them. `MulderTable` checks this against Mulder under a
slab every time it builds a table (`docs/notes/flux.md`).

With a sea-level reference above sea level the open-sky flux exists only above the altitude
cut-off K_cut(elevation), where it has an integrable singularity. The energy nodes of each
elevation are therefore K_cut + d on one grid of d, the energy above the cut-off, which starts
a margin above zero: the table holds no energy at or below the cut-off, and a request there is
refused. Under an opacity the cut-off left at the detector is lower, and under enough of it
there is none.

Units: kinetic energy at the detector in GeV, opacity in g/cm², differential flux in
m⁻² s⁻¹ sr⁻¹ GeV⁻¹, integral flux in m⁻² s⁻¹ sr⁻¹. Outside the grid every method raises
ValueError: nothing is extrapolated and nothing comes back as zero or nan.
"""

from __future__ import annotations

import json
import os
import pathlib
from dataclasses import dataclass

import numpy as np

OPACITY_MAX_G_CM2 = 3.0e4
"""Largest opacity served, g/cm²: the limit of continuous transport in the conventions."""
FORMAT = 1
"""Version of the file layout written by `FluxTable.save`."""
REFINE = 4
"""Trapezoids per interval of the energy grid in the integral over energy."""
RESPONSE_STEP = 0.02
"""Largest step in ln K of the trapezoids between two nodes of an energy response."""
EDGE = 1e-9
"""Relative slack at the edges of the grids, for values that are on an edge up to round-off."""


@dataclass(frozen=True)
class EnergyResponse:
    """Efficiency of a detector against the kinetic energy of the muon at the detector.

    `kinetic_energy_gev` are increasing nodes in GeV and `efficiency` the values there, between
    0 and 1. The efficiency is zero below the first node, linear between nodes and constant
    above the last one. A sharp threshold is one node: `EnergyResponse.threshold(1.0)`.
    """

    kinetic_energy_gev: tuple[float, ...]
    efficiency: tuple[float, ...]

    def __post_init__(self) -> None:
        k = np.atleast_1d(np.asarray(self.kinetic_energy_gev, dtype=float))
        e = np.atleast_1d(np.asarray(self.efficiency, dtype=float))
        if k.ndim != 1 or k.shape != e.shape or k.size == 0:
            raise ValueError("an energy response needs as many efficiencies as energies, at least one")
        if not (np.all(np.isfinite(k)) and np.all(k > 0) and np.all(np.diff(k) > 0)):
            raise ValueError("the energies of a response must be positive and increasing")
        if not np.all((e >= 0) & (e <= 1)):
            raise ValueError("an efficiency lies between 0 and 1")
        object.__setattr__(self, "kinetic_energy_gev", tuple(k.tolist()))
        object.__setattr__(self, "efficiency", tuple(e.tolist()))

    @classmethod
    def threshold(cls, k_min_gev: float, efficiency: float = 1.0) -> EnergyResponse:
        """A sharp threshold: `efficiency` at and above the kinetic energy `k_min_gev`, zero below."""
        return cls((k_min_gev,), (efficiency,))

    @property
    def k_min_gev(self) -> float:
        """Lowest kinetic energy the detector responds to, GeV."""
        return self.kinetic_energy_gev[0]

    def __call__(self, kinetic_energy_gev):
        """Efficiency at the kinetic energy (GeV)."""
        k = np.asarray(kinetic_energy_gev, dtype=float)
        value = np.interp(k, self.kinetic_energy_gev, self.efficiency, left=0.0)
        return float(value) if k.ndim == 0 else value


def _frozen(values, name: str) -> np.ndarray:
    array = np.array(values, dtype=float)
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name}: values must be finite")
    array.flags.writeable = False
    return array


def _increasing(array: np.ndarray, name: str) -> None:
    if array.ndim != 1 or array.size < 2 or not np.all(np.diff(array) > 0):
        raise ValueError(f"{name}: at least two increasing values are needed")


class FluxTable:
    """Differential muon flux at a detector against kinetic energy, elevation and opacity.

    altitude_m                 altitude of the detector above sea level, m
    sin_elevation              nodes of the sine of the elevation of the line of sight
    cutoff_gev                 open-sky cut-off K_cut at each of them, GeV; zeros when there is none
    energy_above_cutoff_gev    nodes of d = K - K_cut, GeV
    flux_m2_s_sr_gev           open-sky flux at K_cut + d, m⁻² s⁻¹ sr⁻¹ GeV⁻¹; shape (elevations, energies)
    range_energy_gev, range_g_cm2, stopping_power_gev_cm2_g
                               CSDA range (g/cm²) and stopping power (GeV cm²/g) of the reference
                               material against kinetic energy (GeV)
    opacity_max_g_cm2          largest opacity served, g/cm²
    meta                       how the table was made: plain values that go to JSON

    A table has a cut-off at every elevation or at none. All interpolation is linear in
    logarithms along energy and linear in the sine of elevation.
    """

    def __init__(
        self,
        *,
        altitude_m: float,
        sin_elevation,
        cutoff_gev,
        energy_above_cutoff_gev,
        flux_m2_s_sr_gev,
        range_energy_gev,
        range_g_cm2,
        stopping_power_gev_cm2_g,
        opacity_max_g_cm2: float = OPACITY_MAX_G_CM2,
        meta: dict | None = None,
    ) -> None:
        self.altitude_m = float(altitude_m)
        self.sin_elevation = _frozen(sin_elevation, "sin_elevation")
        self.cutoff_gev = _frozen(cutoff_gev, "cutoff_gev")
        self.energy_above_cutoff_gev = _frozen(energy_above_cutoff_gev, "energy_above_cutoff_gev")
        self.flux_m2_s_sr_gev = _frozen(flux_m2_s_sr_gev, "flux_m2_s_sr_gev")
        self.range_energy_gev = _frozen(range_energy_gev, "range_energy_gev")
        self.range_g_cm2 = _frozen(range_g_cm2, "range_g_cm2")
        self.stopping_power_gev_cm2_g = _frozen(stopping_power_gev_cm2_g, "stopping_power_gev_cm2_g")
        self.opacity_max_g_cm2 = float(opacity_max_g_cm2)
        self.meta = json.loads(json.dumps(meta or {}))

        s, c, d = self.sin_elevation, self.cutoff_gev, self.energy_above_cutoff_gev
        if not np.isfinite(self.altitude_m):
            raise ValueError("altitude_m must be finite")
        _increasing(s, "sin_elevation")
        _increasing(d, "energy_above_cutoff_gev")
        _increasing(self.range_energy_gev, "range_energy_gev")
        _increasing(self.range_g_cm2, "range_g_cm2")
        if s[0] <= 0 or s[-1] > 1:
            raise ValueError("sin_elevation: nodes must lie in (0, 1]")
        if c.shape != s.shape or np.any(c < 0) or (np.any(c > 0) and not np.all(c > 0)):
            raise ValueError("cutoff_gev: one value per elevation, all zero or all positive")
        if self.flux_m2_s_sr_gev.shape != (s.size, d.size) or np.any(self.flux_m2_s_sr_gev <= 0):
            raise ValueError("flux_m2_s_sr_gev: positive values of shape (elevations, energies) are needed")
        if self.range_g_cm2.shape != self.range_energy_gev.shape or self.range_g_cm2[0] <= 0:
            raise ValueError("range_g_cm2: one positive value per energy")
        stopping = self.stopping_power_gev_cm2_g
        if stopping.shape != self.range_energy_gev.shape or np.any(stopping <= 0):
            raise ValueError("stopping_power_gev_cm2_g: one positive value per energy")
        if not self.opacity_max_g_cm2 > 0:
            raise ValueError("opacity_max_g_cm2 must be positive")
        if self.range_energy_gev[0] <= 0 or d[0] < self.range_energy_gev[0]:
            raise ValueError("the energy grid must start at or above the first energy of the range table")
        top = c.min() + d[-1]  # highest open-sky energy every elevation has
        if top > self.range_energy_gev[-1] * (1 + EDGE):
            raise ValueError("the range table must reach the highest energy of the flux grid")

        self._ln_e = np.log(self.range_energy_gev)
        self._ln_r = np.log(self.range_g_cm2)
        self._ln_s = np.log(stopping)
        self._ln_d = np.log(d)
        self._ln_flux = np.log(self.flux_m2_s_sr_gev)
        self._has_cutoff = bool(c[0] > 0)
        self._cutoff_times_sine = c * s
        self._cumulative = None
        self.k_min_gev = float(self.range_energy_gev[0])
        """Lowest kinetic energy at the detector that the table serves anywhere, GeV."""
        residual = float(self._range(min(top, self.range_energy_gev[-1]))) - self.opacity_max_g_cm2
        if residual <= self.range_g_cm2[0]:
            raise ValueError("the flux grid ends below the energy needed to cross opacity_max_g_cm2")
        self.k_max_gev = float(self._energy(residual))
        """Highest kinetic energy at the detector served under every opacity, GeV."""

    def __repr__(self) -> str:
        cut = (
            f"cut-off {self.cutoff_gev.min():.3g}-{self.cutoff_gev.max():.3g} GeV"
            if self._has_cutoff
            else "no cut-off"
        )
        return (
            f"{type(self).__name__}(altitude {self.altitude_m:g} m, sin elevation "
            f"{self.sin_elevation[0]:g}-{self.sin_elevation[-1]:g}, "
            f"K {self.k_min_gev:g}-{self.k_max_gev:.4g} GeV, {cut})"
        )

    # ------------------------------------------------------------- reference material

    def _range(self, kinetic_energy_gev):
        return np.exp(np.interp(np.log(kinetic_energy_gev), self._ln_e, self._ln_r))

    def _energy(self, range_g_cm2):
        return np.exp(np.interp(np.log(range_g_cm2), self._ln_r, self._ln_e))

    def _stopping(self, kinetic_energy_gev):
        return np.exp(np.interp(np.log(kinetic_energy_gev), self._ln_e, self._ln_s))

    def _top(self, k, x):
        """Kinetic energy above the opacity `x` of a muon that has `k` below it."""
        inside = np.clip(k, self.range_energy_gev[0], self.range_energy_gev[-1])
        return np.where(x > 0, self._energy(self._range(inside) + x), k)

    # ------------------------------------------------------------- checks

    def _directions(self, sin_elevation, opacity_g_cm2, *more):
        """Broadcast the arguments and refuse directions outside the table."""
        arrays = np.broadcast_arrays(
            *(np.asarray(a, dtype=float) for a in (sin_elevation, opacity_g_cm2, *more))
        )
        s, x = arrays[0], arrays[1]
        lo, hi = self.sin_elevation[0], self.sin_elevation[-1]
        bad = ~((s >= lo * (1 - EDGE)) & (s <= hi * (1 + EDGE)))
        if np.any(bad):
            raise ValueError(
                f"sin elevation {s[bad].flat[0]:g} is outside the table, {lo:g} to {hi:g}"
                f" ({np.degrees(np.arcsin(lo)):.1f} to {np.degrees(np.arcsin(hi)):.1f} degrees)"
            )
        bad = ~((x >= 0) & (x <= self.opacity_max_g_cm2 * (1 + EDGE)))
        if np.any(bad):
            raise ValueError(
                f"opacity {x[bad].flat[0]:g} g/cm2 is outside the table, 0 to {self.opacity_max_g_cm2:g}"
                " g/cm2: continuous transport is not valid beyond it (conventions, domain of validity)"
            )
        return [np.clip(s, lo, hi), np.clip(x, 0.0, self.opacity_max_g_cm2), *arrays[2:]]

    def _energies(self, k, s, x, what: str) -> None:
        """Refuse kinetic energies `k` that the table does not hold in the directions (s, x)."""
        if not np.all(np.isfinite(k)):
            raise ValueError(f"{what}: the kinetic energy must be finite")
        low = self._energy_min(s, x)
        bad = k < low * (1 - EDGE)
        if np.any(bad):
            i = np.flatnonzero(bad)[0]
            where = f"at sin elevation {s.flat[i]:.4g} under {x.flat[i]:.4g} g/cm2"
            message = (
                f"{what} {k.flat[i]:.4g} GeV is below {low.flat[i]:.4g} GeV,"
                f" the lowest energy of the table {where}"
            )
            if self._has_cutoff and low.flat[i] > self.k_min_gev:
                message += (
                    f": the reference flux has a cut-off of {self._cutoff(s, x).flat[i]:.4g} GeV there."
                    " Raise the detector threshold or use a reference at the site altitude"
                )
            raise ValueError(message)
        bad = k > self.k_max_gev * (1 + EDGE)
        if np.any(bad):
            raise ValueError(
                f"{what} {k[bad].flat[0]:.4g} GeV is above {self.k_max_gev:.4g} GeV,"
                " the highest energy of the table"
            )

    # ------------------------------------------------------------- cut-off

    def _cutoff_open(self, s):
        if not self._has_cutoff:
            return np.zeros(np.shape(s))
        return np.interp(s, self.sin_elevation, self._cutoff_times_sine) / s

    def _left(self, energy_open, x):
        """Energy under the opacity `x` of a muon with `energy_open` above it; 0 if it is stopped."""
        inside = np.clip(energy_open, self.range_energy_gev[0], self.range_energy_gev[-1])
        residual = self._range(inside) - x
        reaches = (energy_open >= self.range_energy_gev[0]) & (residual > self.range_g_cm2[0])
        under = np.where(reaches, self._energy(np.where(reaches, residual, self.range_g_cm2[0])), 0.0)
        return np.where(x > 0, under, energy_open)

    def _cutoff(self, s, x):
        return self._left(self._cutoff_open(s), x)

    def _energy_min(self, s, x):
        first = self._cutoff_open(s) + self.energy_above_cutoff_gev[0]
        return np.maximum(self._left(first, x), self.k_min_gev)

    def cutoff(self, sin_elevation, opacity_g_cm2):
        """Cut-off of the reference flux left at the detector, GeV; 0 where there is none.

        In the open sky it is the energy lost in the air between the site and the altitude of the
        reference; under an opacity it is what is left of that energy after the opacity, and a
        cut-off below the first energy of the range table (1 MeV for Mulder) counts as none.
        """
        s, x = self._directions(sin_elevation, opacity_g_cm2)
        return _like(self._cutoff(s, x), sin_elevation, opacity_g_cm2)

    def energy_min(self, sin_elevation, opacity_g_cm2):
        """Lowest kinetic energy (GeV) the table holds for the direction: requests below it are refused."""
        s, x = self._directions(sin_elevation, opacity_g_cm2)
        return _like(self._energy_min(s, x), sin_elevation, opacity_g_cm2)

    # ------------------------------------------------------------- differential flux

    def _cell(self, s):
        nodes = self.sin_elevation
        i = np.clip(np.searchsorted(nodes, s, side="right") - 1, 0, nodes.size - 2)
        return i, (s - nodes[i]) / (nodes[i + 1] - nodes[i])

    def _node(self, grid, u):
        k = np.clip(np.searchsorted(grid, u, side="right") - 1, 0, grid.size - 2)
        return k, np.clip((u - grid[k]) / (grid[k + 1] - grid[k]), 0.0, 1.0)

    def differential(self, kinetic_energy_gev, sin_elevation, opacity_g_cm2):
        """Differential flux at the detector, m⁻² s⁻¹ sr⁻¹ GeV⁻¹.

        `kinetic_energy_gev` is the kinetic energy at the detector, `opacity_g_cm2` the opacity of
        the reference material along the line of sight. The arguments broadcast together.
        """
        s, x, k = self._directions(sin_elevation, opacity_g_cm2, kinetic_energy_gev)
        self._energies(k, s, x, "kinetic energy")
        top = self._top(k, x)
        d = np.clip(
            top - self._cutoff_open(s), self.energy_above_cutoff_gev[0], self.energy_above_cutoff_gev[-1]
        )
        i, a = self._cell(s)
        j, t = self._node(self._ln_d, np.log(d))
        table = self._ln_flux
        low = table[i, j] * (1 - t) + table[i, j + 1] * t
        high = table[i + 1, j] * (1 - t) + table[i + 1, j + 1] * t
        jacobian = np.where(x > 0, self._stopping(top) / self._stopping(np.maximum(k, self.k_min_gev)), 1.0)
        return _like(
            np.exp(low * (1 - a) + high * a) * jacobian, kinetic_energy_gev, sin_elevation, opacity_g_cm2
        )

    # ------------------------------------------------------------- integral flux

    def _integrals(self):
        """Open-sky integral flux above each node of the refined energy grid, per elevation."""
        if self._cumulative is None:
            n = self._ln_d.size
            position = np.arange((n - 1) * REFINE + 1) / REFINE
            j = np.minimum(position.astype(int), n - 2)
            t = position - j
            fine = self._ln_d[j] * (1 - t) + self._ln_d[j + 1] * t
            ln_flux = self._ln_flux[:, j] * (1 - t) + self._ln_flux[:, j + 1] * t
            integrand = np.exp(ln_flux + fine)  # flux times d: the integral runs in ln d
            steps = 0.5 * (integrand[:, 1:] + integrand[:, :-1]) * np.diff(fine)
            above = np.concatenate(
                [np.cumsum(steps[:, ::-1], axis=1)[:, ::-1], np.zeros((steps.shape[0], 1))], axis=1
            )
            self._cumulative = fine, ln_flux, integrand, above
        return self._cumulative

    def _above(self, k, s, x):
        """Integral flux above the kinetic energy `k` at the detector, trapezoids from `k` itself."""
        fine, ln_flux, integrand, above = self._integrals()
        d = np.clip(
            self._top(k, x) - self._cutoff_open(s),
            self.energy_above_cutoff_gev[0],
            self.energy_above_cutoff_gev[-1],
        )
        u = np.log(d)
        i, a = self._cell(s)
        j, t = self._node(fine, u)
        out = 0.0
        for row, weight in ((i, 1 - a), (i + 1, a)):
            first = np.exp(ln_flux[row, j] * (1 - t) + ln_flux[row, j + 1] * t + u)
            total = 0.5 * (first + integrand[row, j + 1]) * (fine[j + 1] - u) + above[row, j + 1]
            out = out + weight * np.log(np.maximum(total, 1e-300))
        return np.exp(out)

    def integral(self, response: EnergyResponse, sin_elevation, opacity_g_cm2):
        """Flux integrated over the energy response of a detector, m⁻² s⁻¹ sr⁻¹.

        The integral of efficiency times differential flux over the kinetic energy at the
        detector, by trapezoids in the logarithm of energy that start at the first node of the
        response, up to the top of the table. The response is always given: there is no default
        threshold. A response that starts below `energy_min` of a direction is refused.
        """
        if not isinstance(response, EnergyResponse):
            raise TypeError("the energy response of the detector must be an EnergyResponse")
        s, x = self._directions(sin_elevation, opacity_g_cm2)
        nodes = np.array(response.kinetic_energy_gev)
        efficiency = np.array(response.efficiency)
        self._energies(np.full(s.shape, nodes[0]), s, x, "the detector threshold")
        self._energies(np.full(s.shape, nodes[-1]), s, x, "the last energy of the response")
        out = efficiency[0] * self._above(nodes[0], s, x)
        # by parts: the integral of the flux above K against the slope of the efficiency
        for m in range(nodes.size - 1):
            slope = (efficiency[m + 1] - efficiency[m]) / (nodes[m + 1] - nodes[m])
            if slope == 0:
                continue
            width = np.log(nodes[m + 1] / nodes[m])
            k = nodes[m] * np.exp(np.linspace(0.0, width, max(2, int(np.ceil(width / RESPONSE_STEP)) + 1)))
            values = [self._above(kk, s, x) * kk for kk in k]
            for n in range(k.size - 1):
                out = out + slope * 0.5 * (values[n] + values[n + 1]) * np.log(k[n + 1] / k[n])
        return _like(out, sin_elevation, opacity_g_cm2)

    # ------------------------------------------------------------- files

    def _arrays(self) -> dict:
        return {
            "altitude_m": np.array(self.altitude_m),
            "sin_elevation": self.sin_elevation,
            "cutoff_gev": self.cutoff_gev,
            "energy_above_cutoff_gev": self.energy_above_cutoff_gev,
            "flux_m2_s_sr_gev": self.flux_m2_s_sr_gev,
            "range_energy_gev": self.range_energy_gev,
            "range_g_cm2": self.range_g_cm2,
            "stopping_power_gev_cm2_g": self.stopping_power_gev_cm2_g,
            "opacity_max_g_cm2": np.array(self.opacity_max_g_cm2),
        }

    def save(self, path) -> None:
        """Write the table to `path` (NumPy .npz, no pickle); the file appears complete or not at all."""
        path = pathlib.Path(path)
        temporary = path.with_name(f"{path.name}.{os.getpid()}.tmp")
        with open(temporary, "wb") as handle:
            np.savez(handle, format=np.array(FORMAT), meta=np.array(json.dumps(self.meta)), **self._arrays())
        os.replace(temporary, path)

    @staticmethod
    def read(path) -> dict:
        """Arguments of `FluxTable` stored in a file written by `save`."""
        with np.load(path, allow_pickle=False) as data:
            if int(data["format"]) != FORMAT:
                raise ValueError(f"{path}: file layout {int(data['format'])}, this version reads {FORMAT}")
            out = {name: data[name] for name in data.files if name not in ("format", "meta")}
            out["meta"] = json.loads(str(data["meta"]))
        out["altitude_m"] = float(out["altitude_m"])
        out["opacity_max_g_cm2"] = float(out["opacity_max_g_cm2"])
        return out

    @staticmethod
    def load(path) -> FluxTable:
        """The table stored in a file written by `save`: a plain `FluxTable`, whatever class wrote it."""
        return FluxTable(**FluxTable.read(path))


def _like(value, *arguments):
    """A float when every argument is a scalar, otherwise an array."""
    value = np.asarray(value, dtype=float)
    return float(value) if all(np.ndim(a) == 0 for a in arguments) else value
