"""`muoscene.viz.maps`: what a map figure draws is the map, and every bin without a number says why.

Roadmap 0.1-12a. The maps here are made by hand with the variables and attributes of an
`ExpectedMap` (`docs/conventions.md`, "Expected map"), so that no Mulder is needed: a film tilted
by 60 degrees to the south, whose bins are defined, not accepted, below the cut-off, below the
table and behind too much matter. What a figure draws is read back from its artists: the mesh
of values, the overlay of every reason and the texts. Pixels are never compared: fonts and
antialiasing differ between systems. Figures are written to a temporary folder.
"""

import json
import subprocess
import sys

import matplotlib.text
import numpy as np
import pytest
import xarray as xr

from muoscene.detectors import DirectionGrid, Pose
from muoscene.viz import maps, style

GRID = DirectionGrid.square()
TILTED = Pose((0.0, 0.0, 0.5), 180.0, 60.0)


def fake_map(pose=TILTED, hypotheses=("none", "void"), k_min=1.0, check=True):
    """An ExpectedMap of one detector made by hand, with every status code present for the tilted pose."""
    tx, ty = np.meshgrid(GRID.tx_centers, GRID.ty_centers)
    u = pose.direction(tx, ty)
    elevation = np.degrees(np.arcsin(u[..., 2]))
    status = np.zeros(GRID.shape, np.int8)
    status[elevation < 15.7] = 2
    status[elevation < 14.5] = 3
    status[(np.abs(tx - 0.3) < 0.12) & (np.abs(ty + 0.45) < 0.12) & (status == 0)] = 4
    status[(np.hypot(tx, ty) > 1.0) | (elevation <= 0)] = 1
    counts = np.where(status == 0, 1.0e4 * u[..., 2] ** 4 + 1.0, np.nan)
    counts = np.where(status == 1, 0.0, counts)
    void = np.where(status == 0, counts * (1.0 + 0.1 * np.exp(-((tx - 0.2) ** 2 + ty**2) / 0.02)), counts)
    opacity = np.where(status == 1, np.nan, np.where(ty < -0.3, 0.0, 300.0 * (ty + 0.3)))
    sky_status = np.where(status == 4, 0, status).astype(np.int8)
    sky_status[(status == 0) & (tx < -0.6)] = 2  # under matter there is a number, in the open sky none
    transmission = np.where((status == 0) & (sky_status == 0), 0.5 + 0.1 * ty, np.nan)
    n = len(hypotheses)
    per_map = ("detector", "hypothesis", "ty", "tx")
    data = {
        "solid_angle_sr": (("ty", "tx"), GRID.solid_angle_sr),
        "accepted_solid_angle_sr": (
            ("detector", "ty", "tx"),
            np.where(status == 1, 0.0, GRID.solid_angle_sr)[None],
        ),
        "opacity_gcm2": (per_map, np.stack([opacity] * n)[None]),
        "counts": (per_map, np.stack([counts, void, counts][:n])[None]),
        "transmission": (per_map, np.stack([transmission] * n)[None]),
        "status": (per_map, np.stack([status] * n)[None]),
        "open_sky_status": (("detector", "ty", "tx"), sky_status[None]),
        "reference": (("detector", "ty", "tx"), (status == 0)[None]),
        "zenith_deg": (("detector", "ty", "tx"), (90.0 - elevation)[None]),
        "azimuth_deg": (("detector", "ty", "tx"), np.zeros(GRID.shape)[None]),
        "exposure_days": (("detector",), [45.0]),
    }
    station = {
        "detector": {
            "class": "EmulsionStack",
            "tan_max": 1.0,
            "efficiency": 0.9,
            "area_cm2": 120.0,
            "energy_response": {"kinetic_energy_gev": [k_min], "efficiency": [1.0]},
        },
        "pose": {
            "position_m": list(pose.position_m),
            "azimuth_deg": pose.azimuth_deg,
            "tilt_deg": pose.tilt_deg,
            "roll_deg": pose.roll_deg,
        },
        "exposure_days": 45.0,
        "table": {"class": "MulderTable"},
    }
    checks = {"film": {"compared": True, "max": 6.5e-4, "n_directions": 80}} if check else {}
    return xr.Dataset(
        data,
        coords={
            "detector": ["film"],
            "hypothesis": list(hypotheses),
            "ty": GRID.ty_centers,
            "tx": GRID.tx_centers,
            "ty_edges": ("ty_edge", GRID.ty_edges),
            "tx_edges": ("tx_edge", GRID.tx_edges),
        },
        attrs={
            "muoscene_version": "0.0.1.dev0",
            "mulder_version": "0.3.7",
            "calzone_version": "1.2.2",
            "created": "2026-11-28T04:00:00+00:00",
            "transport": "continuous, straight lines of sight",
            "site": json.dumps({"latitude_deg": 43.35, "longitude_deg": 79.08, "altitude_m": 1150.0}),
            "stations": json.dumps({"film": station}),
            "fast_vs_direct": json.dumps(checks),
        },
    )


@pytest.fixture(scope="module")
def expected():
    return fake_map()


def mesh_of(figure):
    (mesh,) = [c for ax in figure.axes for c in ax.collections if c.get_gid() == maps.VALUES_GID]
    return mesh


def overlays(ax):
    """Number of bins drawn for every reason on an Axes."""
    out = {}
    for c in ax.collections:
        gid = c.get_gid() or ""
        if gid.startswith("muoscene-reason-"):
            out[int(gid.rsplit("-", 1)[1])] = len(c.get_paths())
    return out


def texts(figure) -> str:
    return "\n".join(t.get_text() for t in figure.findobj(matplotlib.text.Text))


def legend_lines(figure) -> list[str]:
    return [t.get_text() for legend in figure.findobj(matplotlib.legend.Legend) for t in legend.get_texts()]


def bincount(status) -> dict:
    values, counts = np.unique(status, return_counts=True)
    return {int(v): int(c) for v, c in zip(values, counts, strict=True) if v != 0}


@pytest.mark.parametrize("quantity", ["counts", "opacity", "transmission"])
def test_what_is_drawn_is_the_map(expected, quantity):
    """The mesh holds the variable of the map in every bin with a number and nothing elsewhere."""
    figure = maps.map_figure(expected, "film", quantity)
    drawn = mesh_of(figure).get_array()
    variable = maps.QUANTITIES[quantity][0]
    values = expected[variable].sel(detector="film", hypothesis="none").values
    status = expected["status"].sel(detector="film", hypothesis="none").values
    sky = expected["open_sky_status"].sel(detector="film").values
    if quantity == "counts":
        reasons = status
    elif quantity == "opacity":
        reasons = np.where(status == 1, 1, 0)
    else:
        reasons = np.where((status == 0) & (sky != 0), style.OPEN_SKY + sky, status)
    assert drawn.shape == GRID.shape
    np.testing.assert_array_equal(np.ma.getmaskarray(drawn), reasons != 0)
    np.testing.assert_array_equal(drawn.compressed(), values[reasons == 0])
    ax = figure.axes[0]
    assert overlays(ax) == bincount(reasons)


def test_every_reason_has_its_look_and_a_line_in_the_legend(expected):
    """Five reasons in the transmission of the film: each with its overlay, its label and its bins."""
    figure = maps.map_figure(expected, "film", "transmission", width="full")
    status = expected["status"].sel(detector="film", hypothesis="none").values
    sky = expected["open_sky_status"].sel(detector="film").values
    reasons = np.where((status == 0) & (sky != 0), style.OPEN_SKY + sky, status)
    found = overlays(figure.axes[0])
    assert set(found) == {1, 2, 3, 4, 12}
    lines = legend_lines(figure)
    for code, n in found.items():
        assert f"{style.LOOKS[code].label} ({n} bins)" in lines
        assert n == int(np.count_nonzero(reasons == code))
    looks = {(style.LOOKS[c].grey, style.LOOKS[c].hatch) for c in found}
    assert len(looks) == len(found)  # no two reasons look alike


def test_the_figure_states_what_was_computed(expected):
    """Detector, hypothesis, threshold, exposure, area, efficiency, pose, versions, check, axes."""
    text = texts(maps.map_figure(expected, "film", "counts"))
    for part in (
        "film",
        "hypothesis none",
        "K > 1 GeV at the detector",
        "45 days",
        "120 cm²",
        "efficiency 0.9",
        "azimuth 180°",
        "tilt 60°",
        "altitude 1150.5 m",
        "muoscene 0.0.1.dev0",
        "Mulder 0.3.7",
        "largest 0.065 % on 80",
        "seen from the sky",
    ):
        assert part in text.replace("\n", " "), part
    text = text.replace("\n", " ")
    counts = expected["counts"].sel(detector="film", hypothesis="none").values
    status = expected["status"].sel(detector="film", hypothesis="none").values
    total = style.number(float(counts[status == 0].sum()))
    assert f"{total} tracks in the {int(np.count_nonzero(status == 0))} bins with a number" in text
    bare = texts(maps.map_figure(expected, "film", "counts", annotate=False))
    assert "45 days" not in bare and "not accepted" in bare


def test_scales(expected):
    """Counts spanning more than 30 are logarithmic, opacity from zero linear; limits can be fixed."""
    counts = maps.map_figure(expected, "film", "counts")
    assert mesh_of(counts).norm.__class__.__name__ == "LogNorm"
    opacity = maps.map_figure(expected, "film", "opacity")
    norm = mesh_of(opacity).norm
    assert norm.__class__.__name__ == "Normalize" and norm.vmin == 0.0
    fixed = maps.map_figure(expected, "film", "counts", scale="log", vmin=1.0, vmax=1.0e5)
    assert (mesh_of(fixed).norm.vmin, mesh_of(fixed).norm.vmax) == (1.0, 1.0e5)


def test_a_ratio_is_symmetric_about_one(expected):
    figure = maps.ratio_figure(expected, "film", "void")
    mesh = mesh_of(figure)
    assert mesh.norm.vmin * mesh.norm.vmax == pytest.approx(1.0, rel=1e-12)
    top = expected["counts"].sel(detector="film", hypothesis="void").values
    bottom = expected["counts"].sel(detector="film", hypothesis="none").values
    status = expected["status"].sel(detector="film", hypothesis="none").values
    with np.errstate(invalid="ignore"):  # not accepted: zero over zero
        ratio = top / bottom
    np.testing.assert_allclose(mesh.get_array().compressed(), ratio[status == 0], rtol=1e-15)
    assert mesh.norm.vmax == pytest.approx(np.nanmax(ratio), rel=1e-12)
    (bar,) = [ax for ax in figure.axes if ax.get_label() == "<colorbar>"]
    assert 1.0 in [float(t.get_text()) for t in bar.get_yticklabels()]
    small = maps.ratio_figure(fake_map(hypotheses=("none", "void", "same")), "film", "same")
    assert mesh_of(small).norm.vmax == pytest.approx(1.0 + maps.RATIO_SPAN_MIN)
    with pytest.raises(ValueError, match="two different"):
        maps.ratio_figure(expected, "film", "none")


@pytest.mark.parametrize(
    "pose", [Pose((0, 0, 0)), Pose((0, 0, 0), 300.0, 20.0), TILTED, Pose((0, 0, 0), 30, 45, 25)]
)
def test_the_sky_is_marked_where_it_is(pose):
    """The zenith and the letters of the compass, checked against the directions of their slopes.

    The zenith of a pose with a roll of zero is at (0, −tan β). Every letter sits on a slope whose
    direction on the side of the normal (`Pose.direction`) has the azimuth of the letter, within
    0.3 degrees (900 steps of 0.1 degree in zenith angle), and is above the horizon.
    """
    figure = maps.map_figure(fake_map(pose), "film", "opacity")
    ax = figure.axes[0]
    tx, ty = pose.slopes(np.array([0.0, 0.0, 1.0]))
    if pose.roll_deg == 0:
        assert (tx, ty) == pytest.approx((0.0, -np.tan(np.radians(pose.tilt_deg))), abs=1e-12)
    letters = {t.get_text(): t.get_position() for t in ax.texts if t.get_text() in "NESW" and t.get_text()}
    assert letters, "no letter of the compass"
    for letter, azimuth in maps.COMPASS:
        point = maps.compass_point(pose, azimuth, 1.0)
        if point is None:
            assert letter not in letters
            continue
        u = pose.direction(*point)
        assert u[2] > 0
        found = np.degrees(np.arctan2(u[0], u[1])) % 360.0
        assert min(abs(found - azimuth), 360.0 - abs(found - azimuth)) < 0.3, (letter, found)
        assert letter in letters
    if pose.tilt_deg == 0:
        assert maps.compass_point(pose, 0.0, 1.0) == pytest.approx((0.0, 1.0), abs=2e-3)  # north up
        assert maps.compass_point(pose, 90.0, 1.0) == pytest.approx((1.0, 0.0), abs=2e-3)  # east right


def test_an_array_on_the_grid(expected):
    """Any array with a mask or status codes, as the binned tracks of roadmap 0.1-15 will be."""
    rng = np.random.default_rng(3)
    values = rng.poisson(50.0, GRID.shape).astype(float)
    mask = np.zeros(GRID.shape, bool)
    mask[:3] = True
    status = expected["status"].sel(detector="film", hypothesis="none").values
    figure = maps.array_figure(
        values,
        GRID,
        mask=mask,
        status=status,
        pose=TILTED,
        tan_max=1.0,
        label="tracks per bin",
        title="binned tracks",
        masked_label="not scanned",
    )
    reasons = np.where((status == 0) & mask, style.MASKED, status)
    drawn = mesh_of(figure).get_array()
    np.testing.assert_array_equal(np.ma.getmaskarray(drawn), reasons != 0)
    np.testing.assert_array_equal(drawn.compressed(), values[reasons == 0])
    assert f"not scanned ({int(np.count_nonzero(reasons == style.MASKED))} bins)" in legend_lines(figure)
    bare = maps.array_figure(values, GRID)
    assert not np.any(np.ma.getmaskarray(mesh_of(bare).get_array()))


def test_a_missing_number_is_never_drawn():
    values = np.ones(GRID.shape)
    values[5, 5] = np.nan
    with pytest.raises(ValueError, match="nan or inf"):
        maps.array_figure(values, GRID)
    mask = np.zeros(GRID.shape, bool)
    mask[5, 5] = True
    maps.array_figure(values, GRID, mask=mask)
    with pytest.raises(ValueError, match="shape"):
        maps.array_figure(np.ones((3, 3)), GRID)
    with pytest.raises(ValueError, match="codes 0 to 4"):
        maps.array_figure(np.ones(GRID.shape), GRID, status=np.full(GRID.shape, 7))
    with pytest.raises(ValueError, match="positive"):
        maps.array_figure(-np.ones(GRID.shape), GRID, kind="ratio")


def test_wrong_names_are_refused(expected):
    with pytest.raises(ValueError, match="no detector"):
        maps.map_figure(expected, "other")
    with pytest.raises(ValueError, match="no hypothesis"):
        maps.map_figure(expected, "film", "counts", "other")
    with pytest.raises(ValueError, match="quantity"):
        maps.map_figure(expected, "film", "flux")
    with pytest.raises(ValueError, match="panel"):
        maps.sheet_figure(expected, "film", [("counts", "none", "void")])


def test_a_sheet_of_maps(expected):
    """Counts, opacity, transmission and the ratio: four panels, and the legend in the fifth place."""
    figure = maps.sheet_figure(expected, "film")
    panels = [ax for ax in figure.axes if ax.get_label() != "<colorbar>" and ax.axison]
    assert [ax.get_title() for ax in panels] == [
        "counts · none",
        "opacity · none",
        "transmission · none",
        "ratio void / none",
    ]
    counts = expected["counts"].sel(detector="film", hypothesis="none").values
    status = expected["status"].sel(detector="film", hypothesis="none").values
    (mesh,) = [c for c in panels[0].collections if c.get_gid() == maps.VALUES_GID]
    np.testing.assert_array_equal(mesh.get_array().compressed(), counts[status == 0])
    lines = legend_lines(figure)
    for code in (1, 2, 3, 4, 12):
        assert style.LOOKS[code].label in lines
    assert str(int(np.count_nonzero(status == 2))) in lines


def test_files_are_written(expected, tmp_path):
    """PNG, PDF and SVG, folders made on the way; other types refused. Nothing goes to a screen."""
    figure = maps.map_figure(expected, "film", "counts")
    assert figure.canvas.__class__.__name__ == "FigureCanvasAgg"
    for suffix in (".png", ".pdf", ".svg"):
        path = style.save(figure, tmp_path / "a" / f"map{suffix}")
        assert path.stat().st_size > 10_000
    assert (tmp_path / "a" / "map.png").read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"
    with pytest.raises(ValueError, match="written as"):
        style.save(figure, tmp_path / "map.jpg")
    with pytest.raises(ValueError, match="width"):
        style.new_figure("column")


def test_importing_the_figures_does_not_import_matplotlib():
    code = (
        "import sys\n"
        "import muoscene.viz, muoscene.viz.maps, muoscene.viz.curves, muoscene.viz.style\n"
        "assert 'matplotlib' not in sys.modules, sorted(m for m in sys.modules if 'matplotlib' in m)[:3]\n"
    )
    run = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, timeout=120)
    assert run.returncode == 0, run.stderr


def test_numbers_and_lines_of_text():
    assert style.number(1.771e6) == "1.771·10⁶"
    assert style.number(950.0) == "950"
    assert maps._wrap("a" * 30 + " · " + "b" * 30, 60.0, 8.0) == ["a" * 30, "b" * 30]
    assert maps._wrap("a" * 30, 30.0, 8.0) == ["a" * 18, "a" * 12]  # 30 mm hold 18 characters of 8 pt
    assert maps._tick_labels([0.995, 1.0, 1.005]) == ["0.995", "1.000", "1.005"]
