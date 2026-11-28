"""viz.curves — the flux of a table against the zenith angle above given detector thresholds.

`flux_curve` gives the integral flux of a `FluxTable` above an energy response along a range of
zenith angles, at one opacity (zero for the open sky), and the reason of every angle where the
table has no number; `flux_figure` draws one curve per response on a logarithmic axis. Decided in
roadmap 0.1-12a (`docs/notes/figures.md`):

- Where the table has no number the curve stops, and the figure says why and from which angle:
  with the sea-level reference carried up to 1150 m a threshold of 1 GeV has a flux up to a
  zenith angle of 74.3° (15.7° of elevation, the cut-off of the reference), one of 0.2 GeV none
  in the open sky; beyond 75.5° (14.5° of elevation) the table itself ends. Nothing is drawn
  there in place of the missing number.
- Other curves — a formula another program computes with, a measured spectrum — are given as
  arrays and drawn dashed on the same axes (`others`); muoscene computes none of them.

matplotlib is imported inside the functions.

Layer 5.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from importlib import metadata

import numpy as np

from .. import __version__
from ..flux import EnergyResponse, FluxTable
from ..transport import DEFINED, flux_domain
from . import style

ZENITH_DEG = tuple(np.arange(0.0, 90.0, 0.25).tolist())
"""Zenith angles of a curve by default, degrees: every quarter of a degree below the horizon."""


@dataclass(frozen=True)
class FluxCurve:
    """Integral flux of a table above an energy response against the zenith angle.

    zenith_deg       zenith angles, degrees
    flux_m2_s_sr     integral flux above the response, m⁻² s⁻¹ sr⁻¹; masked where the table has
                     no number
    status           0 where it has one, otherwise the reason: the status codes of an
                     `ExpectedMap` (2 below the cut-off of the reference, 3 elevation outside the
                     table, 4 opacity beyond it)
    response         the energy response
    opacity_g_cm2    the opacity of the reference material the flux is taken under, g/cm²
    edge_deg         the zenith angle, to 1e-6 degrees, where the number ends after the last angle
                     of `zenith_deg` that has one; None when no angle has one or every angle does
    """

    zenith_deg: np.ndarray
    flux_m2_s_sr: np.ma.MaskedArray
    status: np.ndarray
    response: EnergyResponse
    opacity_g_cm2: float

    edge_deg: float | None = None

    @property
    def last_defined_deg(self) -> float | None:
        """Largest zenith angle of `zenith_deg` with a number; None when there is none."""
        good = self.status == DEFINED
        return float(self.zenith_deg[good].max()) if np.any(good) else None


def flux_curve(table: FluxTable, response: EnergyResponse, zenith_deg=ZENITH_DEG, opacity_g_cm2: float = 0.0):
    """The integral flux of `table` above `response` at the zenith angles `zenith_deg` (0 to 90, degrees).

    The table is entered at sin(elevation) = cos(zenith) and the opacity `opacity_g_cm2` of its
    reference material. Angles where it has no number are masked, with their reason in `status`.
    """
    if not isinstance(table, FluxTable):
        raise TypeError("flux_curve takes a FluxTable")
    if not isinstance(response, EnergyResponse):
        raise TypeError("the energy response must be an EnergyResponse")
    zenith = np.asarray(zenith_deg, dtype=float)
    if zenith.ndim != 1 or zenith.size == 0 or not np.all(np.isfinite(zenith)):
        raise ValueError("zenith_deg: a non-empty sequence of finite angles")
    if np.any(zenith < 0) or np.any(zenith >= 90):
        raise ValueError("zenith_deg: angles from 0 up to, not including, 90 degrees")
    if not (math.isfinite(opacity_g_cm2) and opacity_g_cm2 >= 0):
        raise ValueError("opacity_g_cm2 must be finite and non-negative")
    sine = np.cos(np.radians(zenith))
    status = flux_domain(table, response, sine, np.full(sine.shape, float(opacity_g_cm2)))
    good = status == DEFINED
    flux = np.zeros(sine.shape)
    if np.any(good):
        flux[good] = table.integral(response, sine[good], float(opacity_g_cm2))
    zenith = zenith.copy()
    zenith.flags.writeable = False
    edge = _edge(table, response, zenith, status, float(opacity_g_cm2))
    flux = np.ma.MaskedArray(flux, mask=~good)
    return FluxCurve(zenith, flux, status, response, float(opacity_g_cm2), edge)


def _edge(table, response, zenith, status, opacity) -> float | None:
    """Bisection between the last angle with a number and the next one without."""
    good = np.flatnonzero(status == DEFINED)
    if good.size == 0 or good[-1] == zenith.size - 1:
        return None
    low, high = float(zenith[good[-1]]), float(zenith[good[-1] + 1])

    def defined(z: float) -> bool:
        code = flux_domain(table, response, np.array([math.cos(math.radians(z))]), np.array([opacity]))
        return bool(code[0] == DEFINED)

    while high - low > 1e-6:
        middle = 0.5 * (low + high)
        low, high = (middle, high) if defined(middle) else (low, middle)
    return low


def _response_label(response: EnergyResponse) -> str:
    nodes = response.kinetic_energy_gev
    if len(nodes) == 1:
        return f"K > {nodes[0]:g} GeV"
    return f"response {nodes[0]:g}–{nodes[-1]:g} GeV"


@style.styled
def flux_figure(
    table: FluxTable,
    responses,
    *,
    opacity_g_cm2: float = 0.0,
    others=(),
    zenith_deg=ZENITH_DEG,
    width="half",
    annotate: bool = True,
):
    """Integral flux of a table against the zenith angle, one curve per energy response.

    responses       an `EnergyResponse` or a sequence of them: the thresholds of the detectors
    opacity_g_cm2   opacity of the reference material of the table, g/cm²; 0 is the open sky
    others          curves computed elsewhere, drawn dashed: (label, zenith_deg, flux_m2_s_sr)
    annotate        False leaves out the lines above and under the plot, for a caption to carry
    The figure states the altitude of the table, its reference and its material, and in its
    footer the versions of muoscene and Mulder. Returns a `matplotlib.figure.Figure`.
    """
    import matplotlib
    from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter

    if isinstance(responses, EnergyResponse):
        responses = (responses,)
    responses = tuple(responses)
    if not responses or not all(isinstance(r, EnergyResponse) for r in responses):
        raise TypeError("responses: an EnergyResponse or a sequence of them")
    extra = []
    for other in others:
        label, z, f = other
        z, f = np.asarray(z, dtype=float), np.asarray(f, dtype=float)
        if z.shape != f.shape or z.ndim != 1 or not (np.all(np.isfinite(z)) and np.all(np.isfinite(f))):
            raise ValueError(f"{label}: zenith angles and fluxes of the same length, finite")
        extra.append((str(label), z, f))
    curves = [flux_curve(table, r, zenith_deg, opacity_g_cm2) for r in responses]
    figure = style.new_figure(width, aspect=1.25 if annotate else 0.75)
    styles = ("-", "-.", ":", (0, (5, 1, 1, 1)))
    greys = ("0.0", "0.35", "0.55", "0.2")
    with matplotlib.rc_context(style.rc()):
        ax = figure.add_subplot()
        notes = []
        for n, curve in enumerate(curves):
            good = curve.status == DEFINED
            line = ax.plot(
                curve.zenith_deg[good],
                curve.flux_m2_s_sr.compressed(),
                linestyle=styles[n % len(styles)],
                color=greys[n % len(greys)],
                label=_response_label(curve.response),
            )[0]
            line.set_gid(f"muoscene-curve-{n}")
            notes.append(_where_it_stops(curve))
        for label, z, f in extra:
            ax.plot(z, f, linestyle="--", color="tab:red", linewidth=0.9, label=label)
        shades = _shade(ax, curves[0].zenith_deg, curves)
        ax.set_xlim(0.0, 90.0)
        ax.set_xticks(np.arange(0.0, 91.0, 15.0))
        ax.set_xlabel("zenith angle, degrees")
        ax.set_ylabel("flux above the threshold, m⁻² s⁻¹ sr⁻¹" if not annotate else "flux, m⁻² s⁻¹ sr⁻¹")
        if any(np.any(c.status == DEFINED) for c in curves) or extra:
            ax.set_yscale("log")
            ax.yaxis.set_major_locator(LogLocator(subs=(1.0, 2.0, 5.0)))
            ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
            ax.yaxis.set_minor_formatter(NullFormatter())
        handles, _ = ax.get_legend_handles_labels()
        ax.legend(handles=handles + shades, loc="lower left", fontsize=style.FOOTER_PT)
        ax.grid(True, which="major", linewidth=0.3, color="0.85")
    where = "open sky" if opacity_g_cm2 == 0 else f"under {opacity_g_cm2:g} g/cm² of its reference material"
    inputs = table.meta.get("inputs", {})
    material = table.meta.get("material_name")
    reference = inputs.get("reference")
    header = [
        f"Flux above the threshold · table at {table.altitude_m:g} m · {where}",
        " · ".join(
            x
            for x in (
                f"reference {reference}" if reference else "",
                f"{inputs['atmosphere']} atmosphere" if inputs.get("atmosphere") else "",
                f"material {material}" if material and opacity_g_cm2 > 0 else "",
            )
            if x
        ),
        *notes,
    ]
    footer = [f"muoscene {__version__} · Mulder {_version('mulder')} · continuous transport"]
    from .maps import _texts

    if annotate:
        _texts(figure, header, footer)
    return figure


def _where_it_stops(curve: FluxCurve) -> str:
    label = _response_label(curve.response)
    last = curve.last_defined_deg
    if last is None:
        reasons = sorted(set(curve.status.tolist()) - {DEFINED})
        why = ", ".join(style.LOOKS[r].short for r in reasons)
        return f"{label}: no number at any angle ({why})"
    beyond = curve.status[curve.zenith_deg > last]
    if beyond.size == 0:
        return f"{label}: a number at every angle drawn"
    reason = style.LOOKS[int(beyond[0])].short
    edge = curve.edge_deg if curve.edge_deg is not None else last
    return f"{label}: a number up to {edge:.2f}°, then {reason}"


def _shade(ax, zenith, curves) -> list:
    """Shade the zenith angles of every reason for no number, once, for the first curve that has it."""
    from matplotlib.patches import Patch

    handles = []
    seen = set()
    for curve in curves:
        status = curve.status
        for code in sorted(set(status.tolist()) - {DEFINED}):
            if code in seen:
                continue
            hit = zenith[status == code]
            look = style.LOOKS[int(code)]
            ax.axvspan(
                hit.min(),
                90.0 if hit.max() >= zenith.max() else hit.max(),
                facecolor=str(look.grey),
                edgecolor=style.HATCH_COLOUR,
                hatch=look.hatch or None,
                linewidth=0.0,
                alpha=0.55,
                zorder=0,
            ).set_gid(f"muoscene-reason-{int(code)}")
            seen.add(code)
            handles.append(
                Patch(
                    facecolor=str(look.grey),
                    edgecolor=style.HATCH_COLOUR,
                    hatch=look.hatch or None,
                    linewidth=0.0,
                    alpha=0.55,
                    label=f"{_response_label(curve.response)}: {look.short}",
                )
            )
    return handles


def _version(name: str) -> str:
    try:
        return metadata.version(name)
    except metadata.PackageNotFoundError:
        return "?"
