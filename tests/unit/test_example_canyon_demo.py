"""`examples/canyon_demo`: it runs, its detectors see what its scene says, a wrong line is refused.

Roadmap 0.1-10. The example is loaded from its file and run as it stands, with the flux tables
cached in a temporary folder and the check against Mulder on 80 lines of sight per detector and
threshold instead of 256: two maps of six detectors, 2.4e5 lines scanned and 960 given to Mulder
for the direct flux. The run as a user starts it, `python run.py`, is marked `validation`.

The counts depend on Mulder; tolerances are those of the platform spread of the frozen cases
(3e-6) widened to 1e-3. The number of bins next to the cut-off of the reference may differ by
one between platforms: those counts are compared within one bin.
"""

import importlib.util
import json
import os
import pathlib
import re
import subprocess
import sys

import numpy as np
import pytest

from muoscene.model import expected_counts
from muoscene.transport import DEFINED, NOT_ACCEPTED

PATH = pathlib.Path(__file__).resolve().parents[2] / "examples" / "canyon_demo" / "run.py"


def load():
    spec = importlib.util.spec_from_file_location("canyon_demo_run", PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


EXAMPLE = load()
NAMES = ["L12n", "L12w", "L22n", "L22w", "R12", "R22"]


@pytest.fixture(scope="module")
def result(tmp_path_factory):
    with pytest.MonkeyPatch.context() as patch:
        patch.setenv("MUOSCENE_CACHE", str(tmp_path_factory.mktemp("cache")))
        return EXAMPLE.run(check_directions=80)


def statuses(expected, name):
    status = expected["status"].sel(detector=name, hypothesis="none").values
    return np.bincount(status.ravel(), minlength=5).tolist()


def test_where_the_detectors_are(result):
    """Niches and pins where the table of the example puts them; one flux table per altitude."""
    setup = result.setup
    assert [row.name for row in setup.rows] == NAMES
    assert setup.scene_points["L12n"] == (-0.25, -1.0, -12.0)  # in a niche, 0.25 m behind the left wall
    assert setup.scene_points["L12w"] == (0.15, 1.0, -12.0)  # on the wall, 0.15 m in front of it
    assert setup.scene_points["R22"] == (29.85, 0.0, -17.0)  # 22 m below the right edge at +5 m
    for name, point in setup.scene_points.items():
        assert setup.enu[name] == pytest.approx(tuple(setup.site.enu_from_scene(point)), abs=1e-12)
    assert sorted(result.tables) == [1128.0, 1133.0, 1138.0, 1143.0]
    names = {body.name for body in setup.scene.bodies}
    assert {"niche L12n", "niche L22n", "niche R12"} <= names and "niche L12w" not in names
    lines = EXAMPLE.echo(setup)
    assert any(line.split()[:2] == ["R22", "wall"] and line.rstrip().endswith("15 deg") for line in lines)


def test_above_1_gev_every_bin_has_a_number(result):
    """1324 accepted bins for every detector; the lowest lines, 15 degrees from R22, are under soil."""
    expected = result.maps[1.0]
    for name in NAMES:
        defined, refused, *rest = statuses(expected, name)
        assert refused == 276 and defined + sum(rest) == 1324
        assert defined >= 1323, name


def test_above_0_2_gev_the_open_sky_has_no_number(result):
    """Below the cut-off of the sea-level reference: more of it for a film on the wall than in a niche."""
    expected = result.maps[0.2]
    defined = {name: statuses(expected, name)[DEFINED] for name in NAMES}
    met = {"L12n": 1292, "L12w": 942, "L22n": 990, "L22w": 622, "R12": 1276, "R22": 400}
    for name, value in met.items():
        assert abs(defined[name] - value) <= 1, (name, defined[name])
    status = expected["status"].sel(detector="L22w").values
    counts = expected["counts"].sel(detector="L22w").values
    assert np.all(np.isnan(counts[status > NOT_ACCEPTED]))  # a missing number, never a zero


def test_the_anomalies_are_seen_where_the_scene_puts_them(result):
    """The void behind the left wall from the left detectors, the decompaction from R12 alone.

    Above 1 GeV: the void adds 1947, 1614, 360 and 345 tracks to L12n, L12w, L22n and L22w, up to
    +12.6, +12.1, +10.6 and +10.8 % in a bin; R22 across the canyon sees 16 tracks of it, up to
    3.0 %; the decompaction adds 1189 tracks to R12, up to +3.6 %.
    """
    expected = result.maps[1.0]
    met = {
        ("L12n", "void"): (1946.8, 0.1258),
        ("L12w", "void"): (1614.4, 0.1207),
        ("L22n", "void"): (360.1, 0.1056),
        ("L22w", "void"): (345.4, 0.1080),
        ("R22", "void"): (15.6, 0.0304),
        ("R12", "decompaction"): (1188.7, 0.0364),
    }
    for name in NAMES:
        for hypothesis in ("void", "decompaction"):
            more, _, ratio = EXAMPLE._effect(expected, name, hypothesis, "none")
            if (name, hypothesis) in met:
                tracks, largest = met[(name, hypothesis)]
                assert more == pytest.approx(tracks, rel=1e-2), (name, hypothesis)
                assert float(np.ma.max(np.abs(ratio))) == pytest.approx(largest, abs=1e-3), (name, hypothesis)
            else:
                assert abs(more) < 1.0 and float(np.ma.max(np.abs(ratio))) < 1e-3, (name, hypothesis)
    tracks = {
        "L12n": 1.7707e6,
        "L12w": 2.1469e6,
        "L22n": 2.2928e6,
        "L22w": 2.7511e6,
        "R12": 1.7326e6,
        "R22": 4.9432e6,
    }
    for name, value in tracks.items():
        assert float(expected_counts(expected, name, "none").sum()) == pytest.approx(value, rel=1e-3), name


def test_a_niche_against_pins(result):
    """On the wall a film collects 21 and 20 % more tracks above 1 GeV: part of its cone is open sky.

    382 and 662 of its bins see the open sky, against 0 and 52 in the niche. The void changes a bin
    by about as much, and adds 17 and 4 % fewer tracks.
    """
    expected = result.maps[1.0]
    pairs = [(n.name, w.name) for n, w in EXAMPLE._pairs(result.setup.rows)]
    assert pairs == [("L12n", "L12w"), ("L22n", "L22w")]
    for niche, wall in pairs:
        n0, n1 = (float(expected_counts(expected, name, "none").sum()) for name in (niche, wall))
        assert 1.18 < n1 / n0 < 1.22
        sky = [int(expected["reference"].sel(detector=name).values.sum()) for name in (niche, wall)]
        assert sky == {"L12n": [0, 382], "L22n": [52, 662]}[niche]
        e0, e1 = (EXAMPLE._effect(expected, name, "void", "none") for name in (niche, wall))
        assert abs(float(np.ma.max(e0[2])) - float(np.ma.max(e1[2]))) < 0.01
        assert e1[0] < e0[0]


def test_every_run_checks_itself_against_mulder(result):
    """80 lines per detector and threshold; met at most 0.065 % above 1 GeV and 0.14 % above 0.2 GeV."""
    for threshold, limit in ((1.0, 1e-3), (0.2, 2e-3)):
        checks = json.loads(result.maps[threshold].attrs["fast_vs_direct"])
        for name in NAMES:
            assert checks[name]["compared"] and checks[name]["n_directions"] == 80
            assert checks[name]["max"] < limit, (threshold, name, checks[name]["max"])
        assert result.maps[threshold].attrs["check_seed"] == EXAMPLE.SEED


def test_the_report(result):
    lines = EXAMPLE.report(result)
    text = "\n".join(lines)
    assert text.count("Above 1 GeV") == 1 and text.count("Above 0.2 GeV") == 1
    assert "Bins without a number: none." in text and "below the cut-off of the reference flux" in text
    assert "L12n / L12w" in text and "L22n / L22w" in text
    assert any(line.strip().startswith("R12") and "not seen" in line for line in lines)


@pytest.mark.parametrize(
    ("line", "message"),
    [
        (("A", "middle", 0.0, 12.0, "niche", 300.0, 0.0, 60), "detector A: bank is 'left' or 'right'"),
        (("A", "left", 0.0, 12.0, "pins", 300.0, 0.0, 60), "detector A: mount is 'niche' or 'wall'"),
        (("A", "left", 0.0, 45.0, "wall", 300.0, 0.0, 60), "detector A: .*not on the left wall"),
        (("A", "left", 299.9, 12.0, "niche", 300.0, 0.0, 60), "detector A: .*outside the model"),
        (("A", "left", 0.0, 12.0, "niche", 300.0, 0.0, 0), "detector A: exposure_days must be positive"),
        (("A", "left", 0.0, 12.0, "niche", 300.0, 0.0), "line 1 of DETECTORS: 8 values are needed"),
    ],
)
def test_a_wrong_line_is_refused_without_a_traceback(line, message, capsys):
    assert EXAMPLE.main([line]) == 2
    error = capsys.readouterr().err
    assert error.startswith("canyon_demo: ") and "Traceback" not in error
    assert re.search(message, error), error


def test_a_film_tilted_too_far_is_named(monkeypatch, tmp_path):
    """Tilted by 50 degrees into the massif the lines of a film reach the far end of the model.

    The refusal of the scan names the detector, the second of two, and says what to change; no
    line is given to Mulder (`check_directions=0`).
    """
    monkeypatch.setenv("MUOSCENE_CACHE", str(tmp_path))
    rows = [
        ("ok", "left", 0.0, 22.0, "wall", 300.0, 0.0, 60),
        ("B", "left", 5.0, 22.0, "wall", 300.0, 50.0, 60),
    ]
    with pytest.raises(
        EXAMPLE.InputError, match=r"^detector B: .*cross a limit of the model.*lower tilt_deg"
    ):
        EXAMPLE.run(rows, check_directions=0)


@pytest.mark.validation
def test_the_example_as_a_user_runs_it(tmp_path):
    """`python run.py`: 256 lines per detector and threshold, 3072 given to Mulder; 4 s on the Linux VM."""
    environment = dict(os.environ, MUOSCENE_CACHE=str(tmp_path))
    run = subprocess.run(
        [sys.executable, str(PATH)],
        env=environment,
        capture_output=True,
        text=True,
        timeout=600,
        cwd=tmp_path,
    )
    assert run.returncode == 0, run.stderr
    assert "(256 lines)" in run.stdout and "Done in" in run.stdout


def test_the_figures_of_the_example(result, tmp_path):
    """The sheet of L12n above each threshold: counts, opacity, transmission and the two ratios;
    what the counts panel draws is `expected_counts`. Two sheets of the twelve a run writes; the
    folder of the example gets nothing."""
    from muoscene.viz import maps

    before = sorted(p.name for p in PATH.parent.iterdir())
    written = EXAMPLE.figures(result, tmp_path, names=["L12n"])
    assert sorted(p.name for p in written) == ["canyon_demo_L12n_0.2GeV.png", "canyon_demo_L12n_1GeV.png"]
    assert all(p.stat().st_size > 10_000 for p in written)
    assert sorted(p.name for p in PATH.parent.iterdir()) == before
    figure = maps.sheet_figure(result.maps[1.0], "L12n")
    titles = [ax.get_title() for ax in figure.axes if ax.get_title()]
    assert titles == [
        "counts · none",
        "opacity · none",
        "transmission · none",
        "ratio void / none",
        "ratio decompaction / none",
    ]
    (mesh,) = [c for c in figure.axes[0].collections if c.get_gid() == maps.VALUES_GID]
    np.testing.assert_array_equal(mesh.get_array(), expected_counts(result.maps[1.0], "L12n", "none"))
