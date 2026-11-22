"""`examples/slab`: the films under a uniform slab hold the known answers, and those of Mulder frozen.

Roadmap 0.1-10. The opacity of every bin is density × depth / cos θ to round-off; the flux along
the vertical and along 45 degrees is that of Mulder. Run with the loam of the frozen cases at
their density, the example reproduces the frozen `threshold_slab` cases of
`tests/regression/mulder_cases.py`: under 0.5, 2 and 20 m at 1150 m, above 0.2 and 1 GeV,
tolerance 3e-4 as in `test_transport_table.py`. Budget: two runs; 384 and 0 random lines of the
maps and 2 × 192 along the two lines given to Mulder; 6 × 2e4 lines scanned per run.
"""

import functools
import importlib.util
import json
import pathlib

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


def test_every_check_holds(cache):
    found = EXAMPLE.checks(runs(cache, False))
    assert len(found) == 4 and all(holds for holds, _ in found), found


def test_the_flux_along_two_lines_is_that_of_mulder(cache):
    """With the loam of the example: within 3e-4 of Mulder through the same scene, met 0.024 %."""
    for (threshold, depth, zenith), (_, deviation) in EXAMPLE.along(runs(cache, False)).items():
        assert abs(deviation) < 3e-4, (threshold, depth, zenith, deviation)


@pytest.mark.parametrize("depth", mc.THRESHOLD_DEPTHS)
def test_the_frozen_slab(cache, depth):
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
