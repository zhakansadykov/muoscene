"""`examples/open_sky`: the first check of an installation holds, as a user runs it and inside.

Roadmap 0.1-10. Under the open sky the expected map is the flux table itself, so the counts of
every bin with a number equal flux × projected area × solid angle × efficiency × time to
round-off. The example is run once as a user starts it, `python run.py`, in a fresh interpreter
(512 lines of sight given to Mulder: no line has a flux above 0.2 GeV), and once loaded from its
file with the check on 64 lines per film (128 more). Flux tables are cached in a temporary folder.
"""

import importlib.util
import os
import pathlib
import subprocess
import sys

import numpy as np
import pytest

from muoscene.flux import EnergyResponse
from muoscene.transport import DEFINED

PATH = pathlib.Path(__file__).resolve().parents[2] / "examples" / "open_sky" / "run.py"


def load():
    spec = importlib.util.spec_from_file_location("open_sky_run", PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


EXAMPLE = load()


@pytest.fixture(scope="module")
def result(tmp_path_factory):
    with pytest.MonkeyPatch.context() as patch:
        patch.setenv("MUOSCENE_CACHE", str(tmp_path_factory.mktemp("cache")))
        return EXAMPLE.run(check_directions=64)


def test_as_a_user_runs_it(tmp_path):
    """A fresh interpreter, from another folder: every check holds and the exit code is 0."""
    environment = dict(os.environ, MUOSCENE_CACHE=str(tmp_path))
    run = subprocess.run(
        [sys.executable, str(PATH)],
        env=environment,
        capture_output=True,
        text=True,
        timeout=600,
        cwd=tmp_path,
    )
    assert run.returncode == 0, run.stdout + run.stderr
    assert "checks passed: 5 of 5" in run.stdout and "FAILED" not in run.stdout


def test_every_check_holds(result):
    found = EXAMPLE.checks(*result)
    assert len(found) == 5 and all(holds for holds, _ in found), found


def test_what_the_sea_level_reference_gives_at_1150_m(result):
    """1 GeV: every bin of the flat film; for the film tilted by 60 degrees to the south 820 bins.

    Of the tilted film 476 bins are not accepted (below the horizon or outside the circle), 270
    look below 14.5 degrees, where the table ends, and 34 below 15.7 degrees, where 1 GeV is
    below the cut-off. 0.2 GeV: no bin of either film has a number. Bins next to the cut-off are
    compared within one.
    """
    _, _, _, tables, maps = result
    flat = np.bincount(maps[1.0]["status"].sel(detector="flat").values.ravel(), minlength=5).tolist()
    tilted = np.bincount(maps[1.0]["status"].sel(detector="tilted").values.ravel(), minlength=5).tolist()
    assert flat == [1324, 276, 0, 0, 0]
    assert abs(tilted[0] - 820) <= 1 and tilted[1] == 476 and abs(tilted[2] - 34) <= 1
    assert abs(tilted[3] - 270) <= 1 and sum(tilted) == 1600
    assert int(np.count_nonzero(maps[0.2]["status"].values == DEFINED)) == 0
    (table,) = tables.values()
    assert float(table.cutoff(1.0, 0.0)) == pytest.approx(0.284, abs=1e-3)
    vertical = float(table.integral(EnergyResponse.threshold(1.0), 1.0, 0.0))
    assert vertical == pytest.approx(69.22, rel=1e-3)  # muons m-2 s-1 sr-1 above 1 GeV at 1150.5 m
    tracks = float(np.nansum(maps[1.0]["counts"].sel(detector="flat").values))
    assert tracks == pytest.approx(4.77e6, rel=1e-3)


# ------------------------------------------------------------------ figures (roadmap 0.1-12a)


def drawn_mesh(figure):
    from muoscene.viz import maps

    (mesh,) = [c for ax in figure.axes for c in ax.collections if c.get_gid() == maps.VALUES_GID]
    return mesh.get_array()


def test_the_figures_of_the_example(result, tmp_path):
    """Five files — the counts of the two films above the two thresholds and the flux of the table —
    written where they are asked for; the folder of the example gets nothing."""
    before = sorted(p.name for p in PATH.parent.iterdir())
    written = EXAMPLE.figures(result, tmp_path / "figures")
    assert sorted(p.name for p in written) == [
        "open_sky_flat_0.2GeV.png",
        "open_sky_flat_1GeV.png",
        "open_sky_flux.png",
        "open_sky_tilted_0.2GeV.png",
        "open_sky_tilted_1GeV.png",
    ]
    assert all(p.parent == tmp_path / "figures" and p.stat().st_size > 10_000 for p in written)
    assert sorted(p.name for p in PATH.parent.iterdir()) == before


def test_a_map_figure_draws_the_map(result):
    """The tilted film above 1 GeV: the counts in the 820 bins with a number, the other 780 masked
    and named in the legend by their reason."""
    import matplotlib.legend

    from muoscene.viz import maps

    _, _, _, _, expected = result
    figure = maps.map_figure(expected[1.0], "tilted", "counts")
    drawn = drawn_mesh(figure)
    status = expected[1.0]["status"].sel(detector="tilted", hypothesis="open sky").values
    counts = expected[1.0]["counts"].sel(detector="tilted", hypothesis="open sky").values
    np.testing.assert_array_equal(np.ma.getmaskarray(drawn), status != DEFINED)
    np.testing.assert_array_equal(drawn.compressed(), counts[status == DEFINED])
    lines = [t.get_text() for legend in figure.findobj(matplotlib.legend.Legend) for t in legend.get_texts()]
    for code in (1, 2, 3):
        n = int(np.count_nonzero(status == code))
        assert any(line.endswith(f"({n} bins)") for line in lines), (code, lines)


def test_the_flux_curve_against_mulder(result):
    """The curve of the table of the example against the open sky of Mulder frozen in
    `mulder_reference.json`, and where a threshold of 1 GeV ends.

    Above total energies of 1, 5 and 10 GeV at elevations of 90 and 45 degrees, met within 6e-5.
    The table is at 1150.5 m and the frozen cases at 1150 m: half a metre moves the open sky by
    at most 1e-4 (0.02 % a metre, `muoscene.model.expected.TABLE_ALTITUDE_M`), and the table
    agrees with the frozen numbers within 7e-5 (`docs/notes/flux.md`): tolerance 2e-4. The curve
    of 1 GeV ends at 74.528 degrees of zenith angle, 15.47 of elevation, where the cut-off of the
    reference plus the 10 MeV of the grid reaches 1 GeV; the frozen cut-offs agree with the table
    within 2e-6, which moves that angle by less than 0.005 degree. 0.2 GeV has no number at all.
    """
    import json

    from muoscene.viz import curves
    from tests.regression import mulder_cases as mc

    (table,) = result[3].values()
    frozen = json.loads(
        (pathlib.Path(mc.__file__).parent / "data" / "mulder_reference.json").read_text(encoding="utf-8")
    )["cases"]
    for elevation in (90, 45):
        for total in (1.0, 5.0, 10.0):
            curve = curves.flux_curve(table, EnergyResponse.threshold(mc.k_above(total)), [90.0 - elevation])
            expected = frozen[f"open_sky/1150m/elev{elevation}"][f"E>{total:g}GeV"]
            assert float(curve.flux_m2_s_sr[0]) == pytest.approx(expected, rel=2e-4), (elevation, total)
    one = curves.flux_curve(table, EnergyResponse.threshold(1.0))
    assert one.edge_deg == pytest.approx(74.528, abs=0.005)
    assert float(one.flux_m2_s_sr[0]) == pytest.approx(69.22, rel=1e-3)
    assert curves.flux_curve(table, EnergyResponse.threshold(0.2)).edge_deg is None
