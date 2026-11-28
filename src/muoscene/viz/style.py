"""viz.style — one style for every figure: sizes, fonts, colours, the look of bins without a number.

Decisions of roadmap 0.1-12a (`docs/notes/figures.md`):

- Sizes are physical, so that the type keeps its size when a journal or a report scales the
  figure to its width: `WIDTHS_MM` gives a half width (two figures side by side), a full width
  and the width of a report. JINST is set in one column on A4 and its manual gives no text width;
  the full width of 160 mm is A4 less margins of 25 mm, to be checked against `\\the\\textwidth`
  of the paper when it is written.
- Type is DejaVu Sans, which comes with matplotlib: the same letters and the same layout on
  Windows, macOS and Linux. 8 pt for text, 7 pt for ticks, 6 pt for the footer.
- Colours. Sequential quantities use viridis: its lightness rises monotonically, so it reads in
  black and white and with the common forms of colour blindness. Ratios use RdBu_r on a scale
  symmetric in the logarithm about one; in black and white the sign is carried by contours,
  solid above one and dashed below.
- Bins without a number are never a colour of the scale, never zero and never empty: each reason
  has its grey and its hatch (`LOOKS`), and a figure lists the reasons it shows, with the number
  of bins, in its legend.
- Files: PNG at 200 dpi (the JINST manual asks for 150 to 600), PDF with the fonts embedded as
  TrueType, SVG with the text as paths.

matplotlib is imported inside the functions: importing `muoscene.viz` does not import it.

Layer 5.
"""

from __future__ import annotations

import pathlib
from dataclasses import dataclass

MM_PER_INCH = 25.4
WIDTHS_MM = {"half": 80.0, "full": 160.0, "report": 180.0}
"""Figure widths, mm: one of two figures side by side, the text width of a paper, a report."""
FONT_PT = 8.0
"""Size of text, pt."""
TICK_PT = 7.0
"""Size of tick labels and legends, pt."""
FOOTER_PT = 6.0
"""Size of the footer with the versions, pt."""
DPI = 200
"""Resolution of PNG files."""
FORMATS = (".png", ".pdf", ".svg")
"""File types `save` writes, by the suffix of the path."""
SEQUENTIAL = "viridis"
DIVERGING = "RdBu_r"
HATCH_COLOUR = "0.30"
LINE_COLOUR = "0.15"


@dataclass(frozen=True)
class Look:
    """How the bins without a number for one reason are drawn.

    label   what the legend says, followed by the number of bins
    grey    face colour, a grey level between 0 (black) and 1 (white)
    hatch   matplotlib hatch pattern; empty for none
    short   the label in the narrow legends of a sheet of several maps
    """

    label: str
    grey: float
    hatch: str
    short: str = ""


# Reasons a bin of a figure has no number. Codes 1 to 4 are the status codes of an ExpectedMap
# (docs/conventions.md, "Expected map"); 12 to 14 the same reasons for the open sky, the
# denominator of the transmission; 20 a bin masked by the caller; 21 a ratio with nothing below.
NOT_ACCEPTED, BELOW_CUTOFF, LOW_ELEVATION, TOO_OPAQUE = 1, 2, 3, 4
OPEN_SKY = 10
MASKED = 20
ZERO_BELOW = 21
LOOKS = {
    NOT_ACCEPTED: Look("not accepted", 0.93, "", "not accepted"),
    BELOW_CUTOFF: Look("below the cut-off of the reference flux", 0.80, "////", "below the cut-off"),
    LOW_ELEVATION: Look("elevation outside the flux table", 0.80, "....", "elevation off the table"),
    TOO_OPAQUE: Look("opacity beyond the flux table", 0.62, "xxxx", "too opaque for the table"),
    OPEN_SKY + BELOW_CUTOFF: Look(
        "open sky below the cut-off: no transmission", 0.80, "\\\\\\\\", "open sky below the cut-off"
    ),
    OPEN_SKY + LOW_ELEVATION: Look(
        "open sky: elevation outside the flux table", 0.80, "oo", "open sky: elevation off the table"
    ),
    OPEN_SKY + TOO_OPAQUE: Look(
        "open sky: opacity beyond the flux table", 0.62, "++", "open sky: too opaque"
    ),
    MASKED: Look("masked", 0.80, "||||", "masked"),
    ZERO_BELOW: Look("nothing expected in the denominator", 0.70, "--", "nothing below"),
}
"""The look of every reason for a bin without a number."""


def styled(function):
    """Run a function that draws with the settings of `rc`: every artist it makes gets them."""
    import functools

    @functools.wraps(function)
    def wrapper(*arguments, **keywords):
        import matplotlib

        with matplotlib.rc_context(rc()):
            return function(*arguments, **keywords)

    return wrapper


def width_mm(width) -> float:
    """A width given by name (`WIDTHS_MM`) or in mm."""
    if isinstance(width, str):
        if width not in WIDTHS_MM:
            raise ValueError(f"width is one of {sorted(WIDTHS_MM)} or a number of mm, not {width!r}")
        return WIDTHS_MM[width]
    if isinstance(width, bool) or not isinstance(width, int | float) or not 20.0 <= width <= 400.0:
        raise ValueError("width in mm must be a number between 20 and 400")
    return float(width)


def rc() -> dict:
    """matplotlib settings of every figure; used with `matplotlib.rc_context` when drawing and saving."""
    return {
        "font.family": "DejaVu Sans",
        "mathtext.fontset": "dejavusans",
        "font.size": FONT_PT,
        "axes.titlesize": FONT_PT,
        "axes.labelsize": FONT_PT,
        "xtick.labelsize": TICK_PT,
        "ytick.labelsize": TICK_PT,
        "legend.fontsize": TICK_PT,
        "legend.frameon": False,
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "lines.linewidth": 1.0,
        "hatch.linewidth": 0.5,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "path",
        "savefig.dpi": DPI,
        "figure.dpi": 100,
        "image.cmap": SEQUENTIAL,
    }


def new_figure(width="half", height_mm: float | None = None, *, aspect: float = 1.0):
    """An empty `matplotlib.figure.Figure` with constrained layout, not attached to any window.

    width        a name of `WIDTHS_MM` or mm
    height_mm    height, mm; by default `aspect` times the width
    """
    import matplotlib
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    w = width_mm(width)
    h = w * aspect if height_mm is None else float(height_mm)
    if not 20.0 <= h <= 600.0:
        raise ValueError("height in mm must be between 20 and 600")
    with matplotlib.rc_context(rc()):
        figure = Figure(figsize=(w / MM_PER_INCH, h / MM_PER_INCH), layout="constrained")
    FigureCanvasAgg(figure)  # drawing and text sizes without a display, on every system
    return figure


def save(figure, path, *, dpi: int = DPI) -> pathlib.Path:
    """Write a figure to `path`, PNG, PDF or SVG by its suffix; the folder is created. Returns the path."""
    import matplotlib
    from matplotlib.figure import Figure

    if not isinstance(figure, Figure):
        raise TypeError("save takes a matplotlib Figure")
    path = pathlib.Path(path)
    if path.suffix.lower() not in FORMATS:
        raise ValueError(f"a figure is written as {', '.join(FORMATS)}, not {path.suffix or path.name!r}")
    if isinstance(dpi, bool) or not isinstance(dpi, int) or not 72 <= dpi <= 1200:
        raise ValueError("dpi must be an integer between 72 and 1200")
    path.parent.mkdir(parents=True, exist_ok=True)
    with matplotlib.rc_context(rc()):
        figure.savefig(path, dpi=dpi, metadata=_metadata(path.suffix.lower()))
    return path


def _metadata(suffix: str) -> dict:
    """No creation date in the file, so that the same figure gives the same bytes on one system."""
    if suffix == ".pdf":
        return {"CreationDate": None, "Creator": "muoscene"}
    if suffix == ".svg":
        return {"Date": None, "Creator": "muoscene"}
    return {"Software": "muoscene"}


def superscript(n: int) -> str:
    """An integer in superscript digits: 6 → '⁶', −3 → '⁻³'."""
    return str(n).translate(str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹"))


def number(x: float, digits: int = 4) -> str:
    """A number for a label: plain below 10⁴, otherwise 1.771·10⁶."""
    if x == 0 or 1e-3 <= abs(x) < 1e4:
        return f"{x:.{digits}g}"
    mantissa, exponent = f"{x:.{digits - 1}e}".split("e")
    return f"{mantissa}·10{superscript(int(exponent))}"
