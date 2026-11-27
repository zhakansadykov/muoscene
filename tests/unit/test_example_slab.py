"""`examples/slab`: the films under a uniform slab hold the known answers, and those of Mulder frozen.

Roadmap 0.1-10. The opacity of every bin is density × depth / cos θ to round-off; the flux along
the vertical and along 45 degrees is that of Mulder. Run with the loam of the frozen cases at
their density, the example reproduces the frozen `threshold_slab` cases of
`tests/regression/mulder_cases.py`: under 0.5, 2 and 20 m at 1150 m, above 0.2 and 1 GeV,
tolerance 3e-4 as in `test_transport_table.py`. Budget: two runs; 384 and 0 random lines of the
maps and 2 × 192 along the two lines given to Mulder; 6 × 2e4 lines scanned per run.

The film "2 m" stands at the centre of the slab, and 4 of the 16 lines of its pencil, 114 of the
lines of its map and 106 of those of the film "20 m" leave the ground through the diagonal of
its top: on macOS the flux of Mulder along such a line never returned and the scan gave wrong
paths along some (roadmap 0.1-10b). Before anything else a child process gives Mulder the 16
lines of the pencil with a timeout, and every test of this file needs it to have answered: a
regression fails here within minutes instead of hanging the run.
"""

import functools
import importlib.util
import json
import pathlib
import subprocess
import sys

import numpy as np
import pytest

from muoscene.materials.composition import Material, Mixture
from tests.regression import mulder_cases as mc

PATH = pathlib.Path(__file__).resolve().parents[2] / "examples" / "slab" / "run.py"
FROZEN = json.loads(
    (pathlib.Path(mc.__file__).parent / "data" / "mulder_reference.json").read_text(encoding="utf-8")
)["cases"]


def load():
    spec = importlib.util.spec_from_file_location("slab_run", PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


EXAMPLE = load()


PENCIL_CHILD = """
import importlib.util, json, sys
from muoscene import geometry as G
from muoscene import materials as M
from muoscene.detectors import Pose
from muoscene.flux import EnergyResponse
from muoscene.site import Site
from muoscene.transport import energy_nodes
spec = importlib.util.spec_from_file_location("slab_run", sys.argv[1])
ex = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ex)
site = Site(**ex.SITE)
half, depth = ex.GROUND["half_size_m"], ex.GROUND["depth_m"]
soil = M.loam(ex.SOIL["dry_density_g_cm3"], ex.SOIL["water_percent"])
scene = site.scene([G.Ground("ground", [(-half, 0.0), (half, 0.0)], -depth, -half, half, soil)])
lines = ex.film(1.0).lines_of_sight(Pose((0.0, 0.0, -2.0)), ex.PENCIL).direction
flux = site.differential_flux(scene, (0.0, 0.0, -2.0), lines, energy_nodes(EnergyResponse.threshold(1.0)))
print(json.dumps(flux.tolist()))
"""
PENCIL_TIMEOUT_S = 300
"""The child takes 15 s on the Linux VM with an empty cache of Mulder (the tables of the loam are
compiled) and 1 s with it; macOS runs 5 to 7 times slower (`docs/notes/geometry.md`)."""
MIRRORS = ((3, 0), (6, 5), (9, 10), (12, 15))
"""Lines of the pencil through the diagonal and lines of the same zenith off it: (tx, ty) and (-tx, ty)."""


@pytest.fixture(scope="module")
def pencil():
    """Flux of Mulder along the 16 lines of the pencil of the film "2 m", from a child process."""
    try:
        run = subprocess.run(
            [sys.executable, "-c", PENCIL_CHILD, str(PATH)],
            capture_output=True,
            text=True,
            timeout=PENCIL_TIMEOUT_S,
        )
    except subprocess.TimeoutExpired:
        pytest.fail(
            f"the flux of Mulder along the pencil of the film at the centre of the slab did not return in"
            f" {PENCIL_TIMEOUT_S} s: a line of sight through an edge of a mesh hangs it (roadmap 0.1-10b)"
        )
    assert run.returncode == 0, run.stderr[-2000:]
    return np.array(json.loads(run.stdout))


def test_the_lines_through_the_diagonal_of_the_top_answer(pencil):
    """Every line has a flux, and the four through the diagonal have that of their mirror: met 6e-15."""
    assert pencil.shape[0] == 16 and np.all(np.isfinite(pencil)) and np.all(pencil[:, 0] > 0)
    for on, off in MIRRORS:
        np.testing.assert_allclose(pencil[on], pencil[off], rtol=1e-9)


def spike_loam():
    """The loam of the spike: the `RegLoam` of the frozen cases as a muoscene material."""
    quartz = Material.from_formula("spike quartz", "SiO2", 2.65, i_ev="mulder")
    calcite = Material.from_formula("spike calcite", "CaCO3", 2.71, i_ev="mulder")
    alumina = Material.from_formula("spike alumina", "Al2O3", 3.95, i_ev="mulder")
    water = Material.from_formula("spike water", "H2O", 1.02, i_ev="mulder")  # Mulder's Water
    return Mixture.of("spike loam", {quartz: 0.60, calcite: 0.15, alumina: 0.12, water: 0.13})


@pytest.fixture(scope="module")
def cache(tmp_path_factory):
    return str(tmp_path_factory.mktemp("cache"))


@functools.cache
def runs(cache, frozen):
    with pytest.MonkeyPatch.context() as patch:
        patch.setenv("MUOSCENE_CACHE", cache)
        if frozen:
            return EXAMPLE.run(check_directions=0, soil=spike_loam(), density_g_cm3=mc.SOIL / 1000.0)
        return EXAMPLE.run(check_directions=64)


def test_every_check_holds(cache, pencil):
    found = EXAMPLE.checks(runs(cache, False))
    assert len(found) == 4 and all(holds for holds, _ in found), found


def test_the_flux_along_two_lines_is_that_of_mulder(cache, pencil):
    """With the loam of the example: within 3e-4 of Mulder through the same scene, met 0.024 %."""
    for (threshold, depth, zenith), (_, deviation) in EXAMPLE.along(runs(cache, False)).items():
        assert abs(deviation) < 3e-4, (threshold, depth, zenith, deviation)


@pytest.mark.parametrize("depth", mc.THRESHOLD_DEPTHS)
def test_the_frozen_slab(cache, depth, pencil):
    """The loam of the frozen cases at 1.8 g/cm³: the frozen `threshold_slab` cases within 3e-4.

    Met −0.003 to +0.019 %; the largest is under 20 m, the decay of muons inside the soil, which a
    table does not have (`docs/notes/forward.md`).
    """
    flux = EXAMPLE.along(runs(cache, True))
    for zenith in mc.THIN_ZENITH:
        frozen = FROZEN[f"threshold_slab/{mc.ALTITUDE:.0f}m/{depth:g}m/zenith{zenith:g}"]
        for k_min in mc.THRESHOLDS:
            got, _ = flux[(k_min, depth, zenith)]
            assert got == pytest.approx(frozen[f"K>{k_min:g}GeV"], rel=3e-4), (depth, zenith, k_min)
