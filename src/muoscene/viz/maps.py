"""viz.maps — angular maps of a detector in the slopes of its frame: an `ExpectedMap` or an array on its grid.

What is drawn, decided in roadmap 0.1-12a (`docs/notes/figures.md`):

- The axes are the slopes tx, ty of the detector (`docs/conventions.md`, "Detector pose"), bins
  as they are, one unit of slope the same length on both axes. tx grows to the right and ty
  upwards: the plane of the slopes is seen from the side of the sky, against the normal of the
  detector — the mirror image of what the detector sees. For a detector with its three angles at
  zero east is to the right and north is up, as on a plan, and the binned tracks of a scan, whose
  axes are those of the film, go onto the same axes without a flip.
- Where the detector looks, without reading the code: dotted lines of equal zenith angle every
  15 degrees with their values, the horizon as a thick line, the zenith (a circle) and the normal
  of the detector at (0, 0) (a cross), and the letters N, E, S, W where the lines of sight towards
  those azimuths leave the acceptance circle — or reach the horizon inside it. In the plane of
  the slopes, a gnomonic projection, every vertical plane is a straight line through the zenith:
  for a roll of zero the zenith is at (0, −tan β) and the horizon is the line ty = cot β.
- A map has a number in a bin or a reason why not. Every reason has its grey and its hatch
  (`muoscene.viz.style.LOOKS`) and a line in the legend with the number of its bins; a bin
  without a number is never drawn as zero or left empty, and a nan where a number is expected is
  refused.
- Scales. `scale="auto"` takes a logarithmic scale when every value is positive and the largest
  is at least `LOG_SPAN` times the smallest, a linear one from zero otherwise; "linear" and "log"
  force it, and `vmin`, `vmax` fix the limits, so that two maps — data and model — share one
  scale. A ratio is drawn on RdBu_r, symmetric in the logarithm about one, at least ±1 % wide.

What each quantity of an `ExpectedMap` is (`map_figure`):

- "counts": `counts`, the expected number of tracks in a bin; a bin with status above 0 has no
  number for its own reason (not accepted, below the cut-off, outside the table).
- "opacity": `opacity_gcm2`, the mean equivalent opacity of the matter over the accepted lines of
  sight, g/cm². It is geometry: it exists in every accepted bin, whatever the flux table holds.
- "transmission": `transmission`, the counts of a bin over the counts the same lines of sight of
  the same detector would give with no matter — with the same energy response, area, efficiency
  and exposure, so that everything but the matter cancels. It has no number where the counts
  have none, nor where the open sky has none (`open_sky_status`): above 0.2 GeV at 1150 m with
  the sea-level reference, nowhere in the open sky.
- `ratio_figure`: the counts under one hypothesis over the counts under another, bin by bin.

The figures are `matplotlib.figure.Figure` objects not attached to any window; nothing is shown
on a screen, and files are written by `muoscene.viz.style.save`. `array_figure` takes any array on
a `DirectionGrid` with a mask or status codes: the binned tracks of roadmap 0.1-15 and the ratio
of data to model of 0.1-16 are drawn by it. matplotlib is imported inside the functions.

Layer 5.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, field

import numpy as np

from ..detectors import DirectionGrid, Pose
from ..model import grid_of
from ..transport import DEFINED, NOT_ACCEPTED
from . import style

QUANTITIES = {
    "counts": ("counts", "expected tracks per bin"),
    "opacity": ("opacity_gcm2", "equivalent opacity, g/cm²"),
    "transmission": ("transmission", "transmission: tracks / open sky"),
}
"""Quantities of `map_figure`: the variable of the map and the label of its colour bar."""
LOG_SPAN = 30.0
"""With scale "auto", a logarithmic scale when the largest value is at least this times the smallest."""
RATIO_SPAN_MIN = 0.01
"""The scale of a ratio spans at least 1 ± this."""
ZENITH_STEP_DEG = 15.0
"""Lines of equal zenith angle are drawn every this many degrees."""
COMPASS = (("N", 0.0), ("E", 90.0), ("S", 180.0), ("W", 270.0))
NARROW_MM = 120.0
"""A figure narrower than this gets the short labels of the reasons and a legend of one column."""
MARGIN = 0.14
"""Room around the grid for the letters of the compass, as a share of its half-width."""
AXES_NOTE = "slopes tx, ty in the frame of the detector; the plane of the slopes seen from the sky"
"""The footer line that says what the axes of a map are."""
VALUES_GID = "muoscene-values"
"""gid of the mesh of values; overlays of bins without a number are "muoscene-reason-<code>"."""


@dataclass(frozen=True)
class Drawn:
    """What `draw_map` put on an Axes.

    values      the array drawn, masked where a bin has no number; shape (bins in ty, bins in tx)
    reasons     0 where a bin has a number, otherwise the code of its reason (`style.LOOKS`)
    bins        number of bins for every reason drawn
    scale       "linear", "log" or "ratio"
    limits      lower and upper end of the colour scale
    zenith      slopes of the zenith, None when there is no pose or it lies in the plane
    letters     slopes where each letter of the compass is put
    """

    values: np.ma.MaskedArray
    reasons: np.ndarray
    bins: dict
    scale: str
    limits: tuple[float, float]
    zenith: tuple[float, float] | None = None
    letters: dict = field(default_factory=dict)


# ------------------------------------------------------------------ public


@style.styled
def draw_map(
    ax,
    values,
    grid: DirectionGrid,
    *,
    mask=None,
    status=None,
    pose: Pose | None = None,
    tan_max: float | None = None,
    kind: str = "sequential",
    label: str = "",
    scale: str = "auto",
    vmin: float | None = None,
    vmax: float | None = None,
    masked_label: str = "masked",
    legend: bool = True,
) -> Drawn:
    """Draw an array over the bins of `grid` on a matplotlib Axes, with a colour bar and a legend.

    values      array of shape `grid.shape` (bins in ty, bins in tx), in the unit of `label`
    mask        bins without a number, True where masked; drawn with the look "masked" and the
                legend line `masked_label`
    status      status codes of an `ExpectedMap` (0 to 4): every code above 0 is a bin without a
                number with its own look. A bin is drawn when its status is 0 and it is not masked
    pose        the pose of the detector, for the zenith, the horizon and the compass; None for none
    tan_max     radius of the acceptance circle in the slopes, drawn when given
    kind        "sequential" or "ratio" (a scale symmetric in the logarithm about one)
    scale       "auto", "linear" or "log" for kind "sequential"
    vmin, vmax  limits of the colour scale; for a ratio only vmax, the scale is 1/vmax to vmax

    A drawn bin must hold a finite number (and a positive one on a logarithmic scale or in a
    ratio): anything else is refused with ValueError, never drawn as zero.
    """
    if not isinstance(grid, DirectionGrid):
        raise TypeError("grid must be a DirectionGrid")
    reasons = np.zeros(grid.shape, dtype=np.int16)
    if status is not None:
        status = _shaped(status, grid, "status")
        if not np.all(np.isin(status, (0, 1, 2, 3, 4))):
            raise ValueError("status holds the codes 0 to 4 of an ExpectedMap")
        reasons = status.astype(np.int16)
    if mask is not None:
        mask = _shaped(mask, grid, "mask").astype(bool)
        reasons = np.where((reasons == 0) & mask, style.MASKED, reasons).astype(np.int16)
    return _draw(
        ax,
        values,
        grid,
        reasons,
        pose=pose,
        tan_max=tan_max,
        kind=kind,
        label=label,
        scale=scale,
        vmin=vmin,
        vmax=vmax,
        looks={style.MASKED: style.Look(masked_label, *_look_of(style.MASKED))},
        legend=legend,
    )


@style.styled
def array_figure(
    values,
    grid: DirectionGrid,
    *,
    mask=None,
    status=None,
    pose: Pose | None = None,
    tan_max: float | None = None,
    kind: str = "sequential",
    label: str = "",
    title: str = "",
    notes=(),
    footer=(),
    masked_label: str = "masked",
    scale: str = "auto",
    vmin: float | None = None,
    vmax: float | None = None,
    width="half",
):
    """A figure of one array on a `DirectionGrid`: binned tracks, a ratio of data to model, any map.

    title, notes, footer   the first line, lines under it, small lines at the bottom
    width                  a name of `muoscene.viz.style.WIDTHS_MM` or mm
    The other arguments are those of `draw_map`. Returns a `matplotlib.figure.Figure`.
    """
    figure = style.new_figure(width, aspect=1.2)
    ax = figure.add_subplot()
    draw_map(
        ax,
        values,
        grid,
        mask=mask,
        status=status,
        pose=pose,
        tan_max=tan_max,
        kind=kind,
        label=label,
        scale=scale,
        vmin=vmin,
        vmax=vmax,
        masked_label=masked_label,
    )
    _texts(figure, [title, *notes], [*footer, AXES_NOTE])
    return figure


@style.styled
def map_figure(
    expected,
    detector: str,
    quantity: str = "counts",
    hypothesis: str | None = None,
    *,
    scale: str = "auto",
    vmin: float | None = None,
    vmax: float | None = None,
    width="half",
    annotate: bool = True,
):
    """One quantity of an `ExpectedMap` for one detector under one hypothesis.

    quantity     "counts", "opacity" or "transmission" (see the docstring of the module)
    hypothesis   a hypothesis of the map; by default the first one
    annotate     False leaves out the lines above and under the map, for a figure whose caption
                 carries them (a paper); the legend of the bins without a number stays
    The figure states the detector, the hypothesis, the energy response, the exposure, the area
    and efficiency, where the normal points, the altitude, and in its footer the versions of
    muoscene and of the backends and the check of the fast path against Mulder.
    """
    info = _Station(expected, detector)
    hypothesis = _hypothesis(expected, hypothesis)
    values, reasons, label = _quantity(expected, detector, hypothesis, quantity)
    figure = style.new_figure(width, aspect=1.2)
    ax = figure.add_subplot()
    drawn = _draw(
        ax,
        values,
        info.grid,
        reasons,
        pose=info.pose,
        tan_max=info.tan_max,
        label=label,
        scale=scale,
        vmin=vmin,
        vmax=vmax,
    )
    title = f"{detector} · {QUANTITIES[quantity][1].split(':')[0]} · hypothesis {hypothesis}"
    if annotate:
        _texts(figure, [title, *info.lines(), *_summary(quantity, drawn)], info.footer())
    return figure


@style.styled
def ratio_figure(
    expected,
    detector: str,
    hypothesis: str,
    reference: str | None = None,
    *,
    vmax: float | None = None,
    width="half",
    annotate: bool = True,
):
    """The counts of one detector under `hypothesis` over those under `reference`, bin by bin.

    reference   by default the first hypothesis of the map. A bin has no number when either map
                has none there (the reason of `hypothesis` first), or when nothing is expected
                under `reference`.
    annotate    as in `map_figure`
    """
    info = _Station(expected, detector)
    reference = _hypothesis(expected, reference)
    hypothesis = _hypothesis(expected, hypothesis)
    if hypothesis == reference:
        raise ValueError("a ratio needs two different hypotheses")
    values, reasons = _ratio(expected, detector, hypothesis, reference)
    figure = style.new_figure(width, aspect=1.2)
    ax = figure.add_subplot()
    _draw(
        ax,
        values,
        info.grid,
        reasons,
        pose=info.pose,
        tan_max=info.tan_max,
        kind="ratio",
        label=f"tracks {hypothesis} / {reference}",
        vmax=vmax,
    )
    title = f"{detector} · ratio of expected tracks · {hypothesis} / {reference}"
    if annotate:
        _texts(figure, [title, *info.lines()], info.footer())
    return figure


@style.styled
def sheet_figure(expected, detector: str, panels=None, *, width="report", annotate: bool = True):
    """Several maps of one detector on one figure: what a run gives for it.

    panels   a sequence of ("counts" | "opacity" | "transmission", hypothesis) and
             ("ratio", hypothesis, reference); a bare name takes the first hypothesis. By default
             counts, opacity and transmission under the first hypothesis, then the ratio of every
             other hypothesis to the first.
    annotate as in `map_figure`
    """
    info = _Station(expected, detector)
    names = [str(h) for h in expected["hypothesis"].values]
    if panels is None:
        panels = [(q, names[0]) for q in QUANTITIES] + [("ratio", h, names[0]) for h in names[1:]]
    specs = [_panel(p, names) for p in panels]
    if not specs:
        raise ValueError("a sheet needs at least one panel")
    columns = min(3, len(specs))
    rows = math.ceil(len(specs) / columns)
    w = style.width_mm(width)
    figure = style.new_figure(w, height_mm=rows * (w / columns * 0.82 + 20.0) + 30.0)
    axes = figure.subplots(rows, columns, squeeze=False)
    for ax in axes.flat[len(specs) :]:
        ax.set_axis_off()
    summary, codes = [], set()
    for ax, spec in zip(axes.flat, specs, strict=False):
        if spec[0] == "ratio":
            _, hypothesis, reference = spec
            values, reasons = _ratio(expected, detector, hypothesis, reference)
            label, kind = f"tracks {hypothesis} / {reference}", "ratio"
        else:
            quantity, hypothesis = spec
            values, reasons, label = _quantity(expected, detector, hypothesis, quantity)
            kind = "sequential"
        drawn = _draw(
            ax,
            values,
            info.grid,
            reasons,
            pose=info.pose,
            tan_max=info.tan_max,
            kind=kind,
            label=label,
            short=True,
        )
        ax.set_title(_panel_title(spec), fontsize=style.FONT_PT)
        n = list(axes.flat).index(ax)
        if n % columns:
            ax.set_ylabel("")
        if n + columns < len(specs):
            ax.set_xlabel("")
            legend = ax.get_legend()
            if legend is not None:  # no axis label under it: the legend moves up
                legend.set_bbox_to_anchor((0.0, -0.13), transform=ax.transAxes)
        if spec[0] == "counts" and not summary:
            summary = _summary("counts", drawn, hypothesis)
        codes |= set(drawn.bins)
    if codes:
        handles = [_patch(style.LOOKS[c], style.LOOKS[c].label) for c in sorted(codes)]
        spare = list(axes.flat)[len(specs) :]
        if spare:  # an empty place of the grid holds what the hatches mean
            spare[0].legend(
                handles=handles, loc="center", title="bins without a number", title_fontsize=style.TICK_PT
            )
        else:
            figure.legend(
                handles=handles,
                loc="outside upper center",
                ncols=2,
                title="bins without a number",
                title_fontsize=style.TICK_PT,
            )
    if annotate:
        _texts(figure, [f"{detector}", *info.lines(), *summary], info.footer())
    return figure


def _patch(look: style.Look, label: str):
    """A legend entry with the grey and the hatch of a reason."""
    from matplotlib.patches import Patch

    return Patch(
        facecolor=str(look.grey),
        edgecolor=style.HATCH_COLOUR,
        linewidth=0.0,
        hatch=look.hatch or None,
        label=label,
    )


# ------------------------------------------------------------------ the maps of an ExpectedMap


class _Station:
    """What a figure says about one detector of a map, read from its attributes."""

    def __init__(self, expected, detector: str) -> None:
        names = [str(d) for d in expected["detector"].values]
        if detector not in names:
            raise ValueError(f"no detector {detector!r} in the map: {', '.join(names)}")
        self.name = detector
        self.attrs = expected.attrs
        stations = json.loads(expected.attrs["stations"])
        self.station = stations[detector]
        pose = self.station["pose"]
        self.pose = Pose(tuple(pose["position_m"]), pose["azimuth_deg"], pose["tilt_deg"], pose["roll_deg"])
        self.tan_max = float(self.station["detector"]["tan_max"])
        self.grid = grid_of(expected)
        self.exposure_days = float(expected["exposure_days"].sel(detector=detector).values)
        self.altitude_m = json.loads(expected.attrs["site"])["altitude_m"] + self.pose.position_m[2]
        checks = json.loads(expected.attrs.get("fast_vs_direct", "{}"))
        self.check = checks.get(detector, {})

    def response(self) -> str:
        response = self.station["detector"]["energy_response"]
        if response == "legacy":
            return "response of the preliminary calculation"
        nodes = response["kinetic_energy_gev"]
        if len(nodes) == 1:
            return f"K > {nodes[0]:g} GeV at the detector"
        return f"response from {nodes[0]:g} to {nodes[-1]:g} GeV"

    def lines(self) -> list[str]:
        d = self.station["detector"]
        area = d.get("area_cm2")
        size = f"{area:g} cm²" if area is not None else d["class"]
        pose = self.pose
        return [
            f"{self.response()} · {self.exposure_days:g} days · {size} · efficiency {d['efficiency']:g}",
            f"normal: azimuth {pose.azimuth_deg:.0f}°, tilt {pose.tilt_deg:.0f}°"
            f" · altitude {self.altitude_m:.1f} m",
        ]

    def footer(self) -> list[str]:
        a = self.attrs
        out = [
            f"muoscene {a.get('muoscene_version', '?')} · Mulder {a.get('mulder_version', '?')}"
            f" · Calzone {a.get('calzone_version', '?')} · map of {str(a.get('created', '?'))[:10]}"
            f" · {a.get('transport', '')}",
            AXES_NOTE,
        ]
        check = self.check
        if check.get("compared"):
            out.append(
                f"fast path against Mulder: largest {check['max'] * 100:.3f} % on {check['n_directions']}"
                " lines of sight"
            )
        elif check:
            out.append(f"fast path not compared with Mulder: {check.get('why', '')}")
        return out


def _hypothesis(expected, hypothesis: str | None) -> str:
    names = [str(h) for h in expected["hypothesis"].values]
    if hypothesis is None:
        return names[0]
    if hypothesis not in names:
        raise ValueError(f"no hypothesis {hypothesis!r} in the map: {', '.join(names)}")
    return hypothesis


def _quantity(expected, detector: str, hypothesis: str, quantity: str):
    """Values, reasons and label of one quantity of a map."""
    if quantity not in QUANTITIES:
        raise ValueError(f"quantity is one of {', '.join(QUANTITIES)}, not {quantity!r}")
    variable, label = QUANTITIES[quantity]
    chosen = expected.sel(detector=detector, hypothesis=hypothesis)
    values = np.asarray(chosen[variable].values, dtype=float)
    status = np.asarray(chosen["status"].values).astype(np.int16)
    if quantity == "counts":
        reasons = status
    elif quantity == "opacity":
        accepted = np.asarray(expected["accepted_solid_angle_sr"].sel(detector=detector).values) > 0
        reasons = np.where(accepted, 0, NOT_ACCEPTED).astype(np.int16)
    else:
        sky = np.asarray(expected["open_sky_status"].sel(detector=detector).values).astype(np.int16)
        reasons = np.where((status == DEFINED) & (sky != DEFINED), style.OPEN_SKY + sky, status)
        reasons = np.where((reasons == 0) & ~np.isfinite(values), style.ZERO_BELOW, reasons)
    return values, reasons.astype(np.int16), label


def _ratio(expected, detector: str, hypothesis: str, reference: str):
    top = expected.sel(detector=detector, hypothesis=hypothesis)
    bottom = expected.sel(detector=detector, hypothesis=reference)
    s_top = np.asarray(top["status"].values).astype(np.int16)
    s_bottom = np.asarray(bottom["status"].values).astype(np.int16)
    reasons = np.where(s_top != DEFINED, s_top, s_bottom)
    below = np.asarray(bottom["counts"].values, dtype=float)
    above = np.asarray(top["counts"].values, dtype=float)
    reasons = np.where((reasons == 0) & ~(below > 0), style.ZERO_BELOW, reasons).astype(np.int16)
    with np.errstate(divide="ignore", invalid="ignore"):
        values = np.where(reasons == 0, above / np.where(below > 0, below, 1.0), np.nan)
    return values, reasons


def _panel(panel, names):
    if isinstance(panel, str):
        panel = (panel,)
    panel = tuple(panel)
    if not panel:
        raise ValueError("an empty panel")
    if panel[0] == "ratio":
        if len(panel) not in (2, 3):
            raise ValueError('a ratio panel is ("ratio", hypothesis) or ("ratio", hypothesis, reference)')
        reference = panel[2] if len(panel) == 3 else names[0]
        for h in (panel[1], reference):
            if h not in names:
                raise ValueError(f"no hypothesis {h!r} in the map: {', '.join(names)}")
        if panel[1] == reference:
            raise ValueError("a ratio needs two different hypotheses")
        return ("ratio", panel[1], reference)
    if panel[0] not in QUANTITIES or len(panel) > 2:
        raise ValueError(f"a panel is (quantity, hypothesis) with a quantity of {', '.join(QUANTITIES)}")
    hypothesis = panel[1] if len(panel) == 2 else names[0]
    if hypothesis not in names:
        raise ValueError(f"no hypothesis {hypothesis!r} in the map: {', '.join(names)}")
    return (panel[0], hypothesis)


def _panel_title(spec) -> str:
    if spec[0] == "ratio":
        return f"ratio {spec[1]} / {spec[2]}"
    return f"{spec[0]} · {spec[1]}"


def _summary(quantity: str, drawn: Drawn, hypothesis: str | None = None) -> list[str]:
    if quantity != "counts":
        return []
    n = int(np.count_nonzero(drawn.reasons == 0))
    total = float(drawn.values.sum()) if n else 0.0
    under = f" under {hypothesis}" if hypothesis else ""
    if n == 0:
        return [f"no bin has a number of tracks{under}"]
    return [f"{style.number(total)} tracks{under} in the {n} bins with a number"]


# ------------------------------------------------------------------ drawing


def _shaped(array, grid: DirectionGrid, name: str) -> np.ndarray:
    out = np.asarray(array)
    if out.shape != grid.shape:
        raise ValueError(f"{name} has the shape {out.shape}; the grid has {grid.shape} (bins in ty, in tx)")
    return out


def _look_of(code: int):
    look = style.LOOKS[code]
    return look.grey, look.hatch


def _draw(
    ax,
    values,
    grid,
    reasons,
    *,
    pose=None,
    tan_max=None,
    kind="sequential",
    label="",
    scale="auto",
    vmin=None,
    vmax=None,
    looks=None,
    legend=True,
    short=False,
) -> Drawn:
    import matplotlib
    from matplotlib import colors
    from matplotlib.collections import PolyCollection
    from matplotlib.patches import Circle, Patch

    looks = {**style.LOOKS, **(looks or {})}
    values = np.asarray(_shaped(values, grid, "values"), dtype=float)
    reasons = np.asarray(reasons, dtype=np.int16)
    unknown = set(np.unique(reasons).tolist()) - {0} - set(looks)
    if unknown:
        raise ValueError(f"no look for the reasons {sorted(unknown)}")
    if kind not in ("sequential", "ratio"):
        raise ValueError('kind is "sequential" or "ratio"')
    if scale not in ("auto", "linear", "log"):
        raise ValueError('scale is "auto", "linear" or "log"')
    if tan_max is not None and not (
        isinstance(tan_max, int | float) and math.isfinite(tan_max) and tan_max > 0
    ):
        raise ValueError("tan_max must be a positive number")
    drawn = reasons == 0
    if not np.all(np.isfinite(values[drawn])):
        raise ValueError(
            "a bin to be drawn holds nan or inf: give the bins without a number in mask or status;"
            " nothing is drawn in their place"
        )
    shown = values[drawn]
    with matplotlib.rc_context(style.rc()):
        if kind == "ratio":
            if np.any(shown <= 0):
                raise ValueError(
                    "a ratio is drawn on a logarithmic scale: every drawn value must be positive"
                )
            if vmax is not None:
                if not (math.isfinite(vmax) and vmax > 1):
                    raise ValueError("vmax of a ratio must be above 1")
                span = math.log(vmax)
            else:
                span = max(
                    float(np.max(np.abs(np.log(shown)))) if shown.size else 0.0, math.log1p(RATIO_SPAN_MIN)
                )
            norm = colors.LogNorm(math.exp(-span), math.exp(span))
            cmap, used = style.DIVERGING, "ratio"
        else:
            used = _scale(shown, scale)
            low, high = _limits(shown, used, vmin, vmax)
            norm = colors.LogNorm(low, high) if used == "log" else colors.Normalize(low, high)
            cmap = style.SEQUENTIAL
        masked = np.ma.MaskedArray(np.where(drawn, values, np.nan), mask=~drawn)
        mesh = ax.pcolormesh(
            grid.tx_edges, grid.ty_edges, masked, cmap=cmap, norm=norm, shading="flat", antialiased=False
        )
        mesh.set_gid(VALUES_GID)
        bins, handles = {}, []
        narrow = ax.figure.get_figwidth() * style.MM_PER_INCH < NARROW_MM
        for code in sorted(set(np.unique(reasons).tolist()) - {0}):
            j, i = np.nonzero(reasons == code)
            x0, x1 = grid.tx_edges[i], grid.tx_edges[i + 1]
            y0, y1 = grid.ty_edges[j], grid.ty_edges[j + 1]
            verts = np.stack([np.stack(c, axis=-1) for c in ((x0, y0), (x1, y0), (x1, y1), (x0, y1))], axis=1)
            look = looks[code]
            overlay = PolyCollection(
                verts,
                facecolors=str(look.grey),
                edgecolors=style.HATCH_COLOUR,
                linewidths=0.0,
                hatch=look.hatch or None,
                antialiased=False,
            )
            overlay.set_gid(f"muoscene-reason-{code}")
            ax.add_collection(overlay)
            bins[code] = int(j.size)
            handles.append(
                Patch(
                    facecolor=str(look.grey),
                    edgecolor=style.HATCH_COLOUR,
                    linewidth=0.0,
                    hatch=look.hatch or None,
                    label=f"{j.size}"
                    if short
                    else f"{(look.short or look.label) if narrow else look.label} ({j.size} bins)",
                )
            )
        if kind == "ratio" and shown.size >= 4:
            _ratio_contours(ax, grid, masked, norm)
        half_x = 0.5 * (grid.tx_edges[-1] - grid.tx_edges[0])
        half_y = 0.5 * (grid.ty_edges[-1] - grid.ty_edges[0])
        margin = MARGIN * max(half_x, half_y)
        ax.set_xlim(grid.tx_edges[0] - margin, grid.tx_edges[-1] + margin)
        ax.set_ylim(grid.ty_edges[0] - margin, grid.ty_edges[-1] + margin)
        ax.set_aspect("equal")
        ax.set_xlabel("slope tx")
        ax.set_ylabel("slope ty")
        circle = None
        if tan_max is not None:
            circle = Circle((0.0, 0.0), tan_max, fill=False, linewidth=0.6, edgecolor=style.LINE_COLOUR)
            ax.add_patch(circle)
        zenith, letters = (None, {})
        if pose is not None:
            zenith, letters = _sky(ax, pose, tan_max, circle)
        if not shown.size:
            ax.text(0.5, 0.5, "no bin has a number", transform=ax.transAxes, ha="center", va="center")
        bar = ax.figure.colorbar(mesh, ax=ax, shrink=0.82, pad=0.02, fraction=0.06)
        bar.set_label(label)
        if kind == "ratio":
            ticks = _ratio_ticks(span)
            bar.set_ticks(ticks)
            bar.set_ticklabels(_tick_labels(ticks))
            bar.minorticks_off()
        if legend and handles and short:  # the bins of each reason; what the hatches mean is in the sheet
            ax.legend(
                handles=handles,
                loc="upper left",
                bbox_to_anchor=(0.0, -0.30),
                ncols=len(handles) if len(handles) <= 3 else 2,
                handlelength=1.2,
                handletextpad=0.3,
                columnspacing=0.8,
                borderaxespad=0.0,
                title="no number, bins:",
                title_fontsize=style.FOOTER_PT,
                fontsize=style.FOOTER_PT,
                alignment="left",
            )
        elif legend and handles:
            ax.figure.legend(
                handles=handles,
                loc="outside upper center",
                ncols=1 if narrow else 2,
                handlelength=1.6,
                title="bins without a number",
                title_fontsize=style.TICK_PT,
            )
    limits = (float(norm.vmin), float(norm.vmax))
    return Drawn(masked, reasons, bins, used, limits, zenith, letters)


def _scale(shown: np.ndarray, scale: str) -> str:
    if scale != "auto":
        return scale
    if shown.size and np.all(shown > 0) and shown.max() >= LOG_SPAN * shown.min():
        return "log"
    return "linear"


def _limits(shown, used, vmin, vmax):
    for name, value in (("vmin", vmin), ("vmax", vmax)):
        if value is not None and not math.isfinite(value):
            raise ValueError(f"{name} must be finite")
    if used == "log":
        if shown.size and np.any(shown <= 0):
            raise ValueError("a logarithmic scale needs every drawn value positive")
        low = vmin if vmin is not None else (float(shown.min()) if shown.size else 1.0)
        high = vmax if vmax is not None else (float(shown.max()) if shown.size else 10.0)
        if not 0 < low < high:
            raise ValueError("a logarithmic scale needs 0 < vmin < vmax")
        return low, high
    low = vmin if vmin is not None else (min(0.0, float(shown.min())) if shown.size else 0.0)
    high = vmax if vmax is not None else (float(shown.max()) if shown.size else 1.0)
    if high <= low:
        high = low + (abs(low) if low else 1.0)
    return low, high


def _ratio_ticks(span: float) -> list[float]:
    from matplotlib.ticker import MaxNLocator

    top = math.exp(span)
    if top <= 2.0:
        ticks = MaxNLocator(nbins=4, steps=[1, 2, 2.5, 5, 10]).tick_values(1.0 / top, top)
        ticks = [float(t) for t in ticks if 1.0 / top - 1e-12 <= t <= top + 1e-12]
        return sorted(set(ticks) | {1.0})
    nice = [k for k in (1.5, 2.0, 3.0, 5.0, 10.0, 20.0, 50.0, 100.0, 1e3, 1e4) if k <= top][-2:]
    return sorted({1.0, *nice, *(1.0 / k for k in nice)})


def _tick_labels(ticks) -> list[str]:
    """Labels with as many decimals as the closest two ticks need to differ."""
    ticks = sorted(ticks)
    step = min((b - a for a, b in zip(ticks, ticks[1:], strict=False)), default=1.0)
    decimals = max(0, -math.floor(math.log10(step))) if step > 0 else 3
    if all(abs(t - round(t, decimals)) < 1e-9 * max(1.0, abs(t)) for t in ticks):
        return [f"{t:.{decimals}f}" for t in ticks]
    return [f"{t:.3g}" for t in ticks]


def _ratio_contours(ax, grid, masked, norm) -> None:
    """Contours of a ratio: solid above one, dashed below, so that the sign reads in black and white."""
    levels = [t for t in _ratio_ticks(math.log(norm.vmax)) if t != 1.0]
    low, high = float(masked.min()), float(masked.max())
    levels = [t for t in levels if low < t < high]
    if not levels:
        return
    x, y = np.meshgrid(grid.tx_centers, grid.ty_centers)
    ax.contour(
        x,
        y,
        masked,
        levels=levels,
        colors=style.LINE_COLOUR,
        linewidths=0.5,
        linestyles=["dashed" if t < 1 else "solid" for t in levels],
    )


def _sky(ax, pose: Pose, tan_max, circle):
    """Lines of equal zenith angle, the horizon, the zenith, the normal and the compass."""
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    tx, ty = np.meshgrid(np.linspace(x0, x1, 241), np.linspace(y0, y1, 241))
    elevation = np.degrees(np.arcsin(np.clip(pose.direction(tx, ty)[..., 2], -1.0, 1.0)))
    if tan_max is not None:  # lines and their labels inside the acceptance circle only
        elevation = np.ma.MaskedArray(elevation, mask=np.hypot(tx, ty) > tan_max * 1.001)
    levels = [90.0 - z for z in np.arange(ZENITH_STEP_DEG, 90.0 - 1e-9, ZENITH_STEP_DEG)][::-1]
    levels = [e for e in levels if elevation.min() < e < elevation.max()]
    clip = circle if circle is not None else ax.patch
    if levels:
        lines = ax.contour(
            tx, ty, elevation, levels=levels, colors=style.LINE_COLOUR, linewidths=0.45, linestyles="dotted"
        )
        lines.set_clip_path(clip)
        ax.clabel(lines, fmt=lambda e: f"{90.0 - e:.0f}°", fontsize=style.FOOTER_PT, inline=True)
    if elevation.min() < 0.0 < elevation.max():
        horizon = ax.contour(tx, ty, elevation, levels=[0.0], colors=style.LINE_COLOUR, linewidths=1.3)
        horizon.set_clip_path(clip)
        ax.clabel(horizon, fmt=lambda e: "horizon", fontsize=style.FOOTER_PT, inline=True)
    zenith = None
    try:
        sx, sy = pose.slopes(np.array([0.0, 0.0, 1.0]))
        zenith = (float(sx), float(sy))
    except ValueError:  # a vertical detector: the zenith lies in its plane
        pass
    ax.plot([0.0], [0.0], marker="+", markersize=6, color=style.LINE_COLOUR, linestyle="none")
    if zenith is not None and math.hypot(*zenith) < 1e-9:
        ax.annotate("zenith = normal", (0.0, 0.0), xytext=(3, 3), textcoords="offset points", fontsize=6)
    else:
        ax.annotate("normal", (0.0, 0.0), xytext=(3, 3), textcoords="offset points", fontsize=6)
        if zenith is not None and x0 <= zenith[0] <= x1 and y0 <= zenith[1] <= y1:
            ax.plot(
                [zenith[0]],
                [zenith[1]],
                marker="o",
                markersize=4,
                markerfacecolor="none",
                color=style.LINE_COLOUR,
                linestyle="none",
            )
            ax.annotate("zenith", zenith, xytext=(3, 3), textcoords="offset points", fontsize=6)
    letters = {}
    radius = tan_max if tan_max is not None else min(x1 - x0, y1 - y0) / 2.0
    for letter, azimuth in COMPASS:
        point = compass_point(pose, azimuth, radius)
        if point is None:
            continue
        letters[letter] = point
        r = math.hypot(*point)
        where = (point[0] * 1.08, point[1] * 1.08) if r > 0.97 * radius else point
        ax.text(
            *where,
            letter,
            ha="center",
            va="center",
            fontsize=style.FONT_PT,
            fontweight="bold",
            color=style.LINE_COLOUR,
        )
    return zenith, letters


def compass_point(pose: Pose, azimuth_deg: float, radius: float, steps: int = 900):
    """Where the lines of sight towards `azimuth_deg` leave the circle of `radius` in the slopes.

    The lines of sight of one azimuth, from the zenith down to the horizon, lie on a straight line
    of the plane of the slopes (a gnomonic projection). The point returned is the one nearest the
    horizon that is in front of the detector and within the circle, in steps of 0.1 degree of
    zenith angle from 0.1; None when no such line is.
    """
    zenith = np.radians(np.linspace(0.0, 90.0, steps, endpoint=False)[1:])  # the zenith has no azimuth
    a = math.radians(azimuth_deg)
    u = np.stack([np.sin(zenith) * math.sin(a), np.sin(zenith) * math.cos(a), np.cos(zenith)], axis=-1)
    w = u @ pose.axes.T
    front = w[:, 2] > 1e-9
    if not np.any(front):
        return None
    sx = np.where(front, w[:, 0] / np.where(front, w[:, 2], 1.0), np.inf)
    sy = np.where(front, w[:, 1] / np.where(front, w[:, 2], 1.0), np.inf)
    inside = front & (np.hypot(sx, sy) <= radius)
    if not np.any(inside):
        return None
    k = np.flatnonzero(inside)[-1]
    return (float(sx[k]), float(sy[k]))


def _texts(figure, header, footer) -> None:
    """The lines above the maps and the small lines under them, broken at " · " to the width of the figure."""
    width = figure.get_figwidth() * style.MM_PER_INCH * 0.96
    header = [h for h in header if h]
    if header:
        lines = [part for h in header for part in _wrap(h, width, style.FONT_PT)]
        figure.suptitle("\n".join(lines), x=0.02, ha="left", fontsize=style.FONT_PT)
    if footer:
        lines = [part for f in footer for part in _wrap(f, width, style.FOOTER_PT)]
        figure.supxlabel("\n".join(lines), x=0.02, ha="left", fontsize=style.FOOTER_PT, color="0.25")


CHARACTER_EM = 0.56
"""Mean width of a character of DejaVu Sans, in em, a little generous: for breaking lines."""
PT_MM = 25.4 / 72.0


def _wrap(line: str, width_mm: float, size_pt: float) -> list[str]:
    """Break a line at " · ", and a part still too long at spaces, to fit `width_mm` at `size_pt`."""
    import textwrap

    most = max(10, int(width_mm / (CHARACTER_EM * size_pt * PT_MM)))
    if len(line) <= most:
        return [line]
    out, current = [], ""
    for part in line.split(" · "):
        joined = f"{current} · {part}" if current else part
        if len(joined) <= most:
            current = joined
            continue
        if current:
            out.append(current)
        pieces = textwrap.wrap(part, most) or [part]
        out += pieces[:-1]
        current = pieces[-1]
    out.append(current)
    return out
