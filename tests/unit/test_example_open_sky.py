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
