"""`expected_map`: the criteria of roadmap 0.1-08 on the invented pier of the spike.

The demonstration scene is the pier of `docs/notes/spike/mulder_spike.py` — a pile cap, 3 × 3
piles and a void in front of them, seen from a niche 20 m below the brow — in the loam and the
concrete of the frozen cases; `canyon_demo` comes with roadmap 0.1-10. Checked here:

- the fast path is within 0.5 % of the direct flux of Mulder, in the run itself and against the
  frozen cases `pier/1150m/…`, which Mulder computed through a geometry built without muoscene;
- a hypothesis without a new scan equals a rebuilt scene;
- bins the flux table has no flux for are refused or marked, never filled;
- the layout of the `ExpectedMap` of the conventions.

The open sky and the uniform slab against Mulder are in `test_transport_table.py` and
`test_transport_direct.py`. Budget of the file: 4e5 lines of sight scanned, 900 given to Mulder
for the direct flux.
"""

import functools
import importlib.util
import json
import math
import pathlib

import numpy as np
import pytest
import xarray as xr

from muoscene import __version__
from muoscene import geometry as G
from muoscene import materials as M
from muoscene.detectors import DirectionGrid, EmulsionStack, LinesOfSight, Pose
from muoscene.flux import LEGACY_RESPONSE, EnergyResponse, LegacyTable, MulderTable
from muoscene.materials.composition import Material, Mixture
from muoscene.model import Hypothesis, HypothesisSet, Station, expected_counts, expected_map, grid_of
from muoscene.model.expected import VARIABLES
from muoscene.site import Site
from muoscene.transport import (
    BELOW_CUTOFF,
    DEFINED,
    LOW_ELEVATION,
    NOT_ACCEPTED,
    equivalent_opacity,
)
from tests.regression import mulder_cases as mc

ROOT = pathlib.Path(__file__).resolve().parents[2]
FROZEN = json.loads(
    (pathlib.Path(mc.__file__).parent / "data" / "mulder_reference.json").read_text(encoding="utf-8")
)["cases"]
SITE = Site(mc.SITE["latitude"], mc.SITE["longitude"], mc.ALTITUDE)
SEA = Site(mc.SITE["latitude"], mc.SITE["longitude"], 0.0)
GRID = DirectionGrid.square()
DAYS = 60.0
SECONDS = DAYS * 86400.0


def _spike():
    """The spike script as a module: the invented pier is defined there and nowhere else."""
    path = ROOT / "docs" / "notes" / "spike" / "mulder_spike.py"
    spec = importlib.util.spec_from_file_location("mulder_spike", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


SPIKE = _spike()
PIER_X0, PIER_SIZE, PIER_CAP, PILE_RADIUS, PILE_LENGTH, PILES, PILE_PITCH = SPIKE.PIER
DETECTOR = tuple(SPIKE.detector(SPIKE.PIER_NICHE).tolist())  # 1 mm above the floor of the niche


def spike_loam():
    """The loam of the spike: the `RegLoam` of the frozen cases as a muoscene material."""
    quartz = Material.from_formula("spike quartz", "SiO2", 2.65, i_ev="mulder")
    calcite = Material.from_formula("spike calcite", "CaCO3", 2.71, i_ev="mulder")
    alumina = Material.from_formula("spike alumina", "Al2O3", 3.95, i_ev="mulder")
    water = Material.from_formula("spike water", "H2O", 1.02, i_ev="mulder")  # Mulder's Water
    return Mixture.of("spike loam", {quartz: 0.60, calcite: 0.15, alumina: 0.12, water: 0.13})


LOAM = spike_loam()
CONCRETE = Material.from_elements("spike concrete", mc.CONCRETE, mc.CONCRETE_DENSITY / 1000.0, i_ev="mulder")
SOIL = mc.SOIL / 1000.0  # g/cm³


@functools.cache
def pier(drho_rel=-1.0, soil_g_cm3=SOIL, site=SITE):
    """The invented pier in a flat bank of loam behind a wall at x = 0, with its niche and its void."""
    ground = G.Ground(
        "soil",
        [(0.0, -200.0), (0.0, 0.0), (400.0, 0.0)],
        -200.0,
        -400.0,
        400.0,
        LOAM,
        density_g_cm3=soil_g_cm3,
    )
    cap = G.Box((PIER_X0, -PIER_SIZE / 2, -PIER_CAP), (PIER_X0 + PIER_SIZE, PIER_SIZE / 2, 0.0))
    centre = (PIER_X0 + PIER_SIZE / 2, 0.0)
    piles = G.pile_grid(centre, -PIER_CAP, PILE_LENGTH, PILE_RADIUS, PILES, PILES, PILE_PITCH)
    niche = G.Box((0.0, -0.15, -SPIKE.PIER_NICHE), (0.5, 0.15, -SPIKE.PIER_NICHE + 0.3))
    void = G.Sphere(SPIKE.PIER_VOID[:3], SPIKE.PIER_VOID[3])
    return site.scene(
        [
            ground,
            G.Structure("cap", cap, CONCRETE),
            G.Structure("piles", piles, CONCRETE),
            G.Cut("niche", niche),
            G.Anomaly("void", void, drho_rel=drho_rel),
        ]
    )


@functools.cache
def table(altitude_m=mc.ALTITUDE - SPIKE.PIER_NICHE):
    """Flux table at the altitude of the niche of the pier, loam of the spike; not cached on disk."""
    return MulderTable(altitude_m, LOAM, cache=False)


def stack(k_min_gev=1.0, **changes):
    """Six films of 120 mm × 100 mm with a sharp threshold."""
    arguments = {
        "area_cm2": 120.0,
        "n_films": 6,
        "tan_max": 1.0,
        "response": EnergyResponse.threshold(k_min_gev),
        "efficiency": 0.9,
    }
    return EmulsionStack(**(arguments | changes))


def station(k_min_gev=1.0, name=None, pose=None, **changes):
    name = name or f"K>{k_min_gev:g}GeV"
    return Station(name, stack(k_min_gev, **changes), pose or Pose(DETECTOR), table(), DAYS)


BOTH = HypothesisSet([Hypothesis.without("void"), Hypothesis.void("void")])


@functools.cache
def demo():
    """The map of the demonstration: two thresholds, with and without the void, compared with Mulder."""
    stations = [station(1.0), station(0.2)]
    return expected_map(SITE, pier(), stations, GRID, BOTH, rng=np.random.default_rng(8), undefined="mask")


def quiet(*arguments, **options):
    """A map without the comparison with Mulder, for what does not need it."""
    options.setdefault("rng", np.random.default_rng(0))
    options.setdefault("undefined", "mask")
    return expected_map(*arguments, check_directions=0, **options)


# ------------------------------------------------------------------ against Mulder


def test_the_demo_is_within_half_a_percent_of_mulder():
    """The criterion of the roadmap, in the run itself: 256 random lines of sight per detector, seed 8.

    Met above 1 GeV: mean +0.006 %, rms 0.014 %, largest 0.032 %. Above 0.2 GeV, on the lines
    that have a flux: mean +0.020 %, rms 0.022 %, largest 0.049 % (0.061 % before the denser
    energy nodes of the direct flux at its threshold, 22 November 2026). The limits asserted are three
    times what was met.
    """
    expected = demo()
    assert expected.attrs["check_seed"] == 8 and expected.attrs["check_directions"] == 256
    checks = json.loads(expected.attrs["fast_vs_direct"])
    hard, soft = checks["K>1GeV"], checks["K>0.2GeV"]
    for check in (hard, soft):
        assert check["compared"] and check["n_directions"] == 256 and check["limit"] == 5e-3
        assert check["hypothesis"] == "the scene as built"
        assert check["max"] < 5e-3  # the roadmap
        assert abs(check["mean"]) <= check["rms"] <= check["max"]
    assert hard["max"] < 1e-3 and hard["rms"] < 4e-4
    assert soft["max"] < 2e-3 and soft["rms"] < 6e-4
    assert expected.attrs["fast_vs_direct_max"] == max(hard["max"], soft["max"])


@pytest.mark.parametrize(("case", "drho_rel"), [("soil", 0.0), ("void", mc.AIR / mc.SOIL - 1.0)])
def test_the_fast_path_against_the_frozen_pier(case, drho_rel):
    """The frozen `pier/1150m/…`: Mulder through a pier built without muoscene, 77 lines of sight.

    Tolerance 1e-3; met within 0.033 % above 1 GeV and 0.050 % above 0.2 GeV on the side of the
    ground. The void of the frozen case holds air of 1.2 kg/m³, which a hypothesis can say.
    """
    frozen = FROZEN[f"pier/{mc.ALTITUDE:.0f}m/{case}"]
    u = mc.pier_directions()
    pose = Pose(DETECTOR)
    paths = SITE.scan(pier(), pose.position_m, u)
    hypothesis = Hypothesis(case, drho_rel={"void": drho_rel})
    materials = {M.register(LOAM): LOAM, M.register(CONCRETE): CONCRETE}
    solid = paths.length_m[:, [i for i, b in enumerate(paths.bodies) if b.kind not in ("cut", "air")]]
    ground = u[:, 0] > 0
    for k_min, mask, key in ((1.0, np.ones(len(u), bool), "K>1GeV"), (0.2, ground, "ground/K>0.2GeV")):
        opacity = equivalent_opacity(
            table(),
            hypothesis.opacity_by_material_g_cm2(paths),
            materials,
            k_min,
            solid_length_m=solid.sum(axis=1),
            air_density_g_cm3=pier().air_density_g_cm3,
        )
        lines = LinesOfSight(pose, u[mask], np.ones(mask.sum()))
        flux = stack(k_min).flux_m2_s_sr(table(), lines, opacity.lookup_g_cm2[mask])
        assert flux == pytest.approx(np.array(frozen[key]), rel=1e-3)


def test_the_void_of_the_frozen_pier_is_seen():
    """One of the 77 lines crosses the void, +12.6 % of flux above 1 GeV; in the map its bins gain 10.9 %."""
    soil, void = (np.array(FROZEN[f"pier/{mc.ALTITUDE:.0f}m/{case}"]["K>1GeV"]) for case in ("soil", "void"))
    effect = void / soil - 1
    assert np.count_nonzero(effect > 1e-6) == 1 and effect.max() == pytest.approx(0.126, abs=1e-3)
    counts = demo()["counts"].sel(detector="K>1GeV")
    gain = (counts.sel(hypothesis="void") / counts.sel(hypothesis="none") - 1).values
    more = float((counts.sel(hypothesis="void") - counts.sel(hypothesis="none")).sum())
    assert more == pytest.approx(572.0, abs=3.0)  # tracks the void adds to the 2.398e6 of the map
    assert np.nanmax(gain) == pytest.approx(0.109, abs=2e-3) and np.nanmin(gain) == 0.0
    j, i = np.unravel_index(np.nanargmax(gain), gain.shape)
    assert GRID.tx_centers[i] == pytest.approx(0.425) and abs(GRID.ty_centers[j]) == pytest.approx(0.025)


# ------------------------------------------------------------------ hypotheses


def test_a_hypothesis_without_a_new_scan_equals_a_rebuilt_scene():
    """Counts under a hypothesis on the scan of one scene against the scene built with those densities.

    No void, a decompaction of 20 %, and a drier bank of 1.7 g/cm³ with the void: paths are the
    same numbers, so the counts agree to 1e-12 (met 3e-15).
    """
    cases = {
        "none": (Hypothesis.without("void"), pier(0.0)),
        "decompaction": (Hypothesis.decompaction("void", -0.2), pier(-0.2)),
        "drier": (Hypothesis("drier", {"soil": 1.7}), pier(-1.0, 1.7)),
    }
    one = quiet(SITE, pier(), station(), GRID, [h for h, _ in cases.values()])
    for name, (_, scene) in cases.items():
        rebuilt = quiet(SITE, scene, station(), GRID, Hypothesis.as_built())
        ours, theirs = one["counts"].sel(hypothesis=name).values, rebuilt["counts"].values[0]
        assert ours == pytest.approx(theirs, rel=1e-12)
        assert one["opacity_gcm2"].sel(hypothesis=name).values == pytest.approx(
            rebuilt["opacity_gcm2"].values[0], rel=1e-12, nan_ok=True
        )
    drier, none = (expected_counts(one, "K>1GeV", name).sum() for name in ("drier", "none"))
    # 5.6 % less soil and the void: 1.7 % more tracks over the map, most of which come from the canyon
    assert drier / none == pytest.approx(1.017, abs=1e-3)
    with pytest.raises(ValueError, match="hypothesis 'cavity'"):
        quiet(SITE, pier(), station(), GRID, Hypothesis("cavity", drho_rel={"hollow": -1.0}))


# ------------------------------------------------------------------ bins without a flux


def test_bins_without_a_flux_are_marked_not_filled():
    """Sea-level reference at 1150 m, a horizontal film in the niche, 1324 bins in its acceptance.

    Above 1 GeV every bin has counts. Above 0.2 GeV 356 bins have a line of sight below the
    cut-off — the canyon and the thin soil above the niche — and the open sky has no flux in any.
    """
    expected = demo()
    status = expected["status"]
    assert status.dtype == np.int8 and status.attrs["flag_meanings"].split()[BELOW_CUTOFF].startswith("below")
    hard, soft = (status.sel(detector=name, hypothesis="none").values for name in ("K>1GeV", "K>0.2GeV"))
    assert np.bincount(hard.ravel(), minlength=5).tolist() == [1324, 276, 0, 0, 0]
    assert np.bincount(soft.ravel(), minlength=5).tolist() == [968, 276, 356, 0, 0]
    counts = expected["counts"].sel(detector="K>0.2GeV", hypothesis="none").values
    assert np.all(np.isnan(counts[soft == BELOW_CUTOFF]))  # the mark of a missing number, never a zero
    assert np.all(counts[soft == DEFINED] > 0) and np.all(counts[soft == NOT_ACCEPTED] == 0.0)
    assert np.all(
        np.isfinite(expected["opacity_gcm2"].sel(detector="K>0.2GeV").values[:, soft != NOT_ACCEPTED])
    )
    # the open sky: defined above 1 GeV from 45 degrees up, nowhere above 0.2 GeV
    sky = expected["open_sky_status"]
    assert np.all(sky.sel(detector="K>1GeV").values[hard == DEFINED] == DEFINED)
    assert np.all(sky.sel(detector="K>0.2GeV").values[hard == DEFINED] == BELOW_CUTOFF)
    assert np.all(np.isnan(expected["transmission"].sel(detector="K>0.2GeV").values))
    # a reference bin is geometry: the same 52 bins whether or not the table has a flux there
    reference = expected["reference"]
    assert reference.dtype == bool and int(reference.sel(detector="K>1GeV").sum()) == 52
    assert np.array_equal(reference.sel(detector="K>1GeV").values, reference.sel(detector="K>0.2GeV").values)
    assert np.all(soft[reference.sel(detector="K>0.2GeV").values] == BELOW_CUTOFF)
    transmission = expected["transmission"].sel(detector="K>1GeV", hypothesis="none").values
    assert np.all(transmission[reference.sel(detector="K>1GeV").values] == 1.0)
    inside = hard == DEFINED
    assert np.all((transmission[inside] > 0.1) & (transmission[inside] <= 1.0))
    assert transmission[inside].min() == pytest.approx(0.154, abs=2e-3)  # 28 m of soil at 45 degrees
    # the masked counts: a sum cannot take a missing number for a zero
    masked = expected_counts(expected, "K>0.2GeV", "none")
    assert masked.mask.sum() == 356 and masked.sum() == pytest.approx(np.nansum(counts))
    assert not np.ma.is_masked(expected_counts(expected, "K>1GeV", "none"))


def test_a_run_refuses_what_it_has_no_flux_for():
    """`undefined="refuse"`: the conventions as they stand. The message says how many lines and why."""
    refused = quiet(SITE, pier(), station(1.0), GRID, BOTH, undefined="refuse")
    assert int((refused["status"] > NOT_ACCEPTED).sum()) == 0
    with pytest.raises(ValueError, match=r"K>0.2GeV.*5162 of 20388.*below cutoff.*45.0 to 67.9 degrees"):
        quiet(SITE, pier(), station(0.2), GRID, BOTH, undefined="refuse")
    for undefined in (None, "nan", True):
        with pytest.raises(ValueError, match="always given"):
            quiet(SITE, pier(), station(), GRID, BOTH, undefined=undefined)
    with pytest.raises(TypeError, match="undefined"):
        expected_map(SITE, pier(), station(), GRID, BOTH, rng=np.random.default_rng(0))


def test_a_tilted_detector_and_the_edge_of_the_table():
    """A film tilted by 30 degrees towards the canyon (azimuth 270) looks down to 15 degrees above 1 GeV.

    Six bins hold lines between 15.1 and 15.7 degrees of elevation, below the cut-off of 1 GeV;
    tilted by 60 degrees part of the cone is below the horizon and below the table.
    """
    low = quiet(SITE, pier(), station(pose=Pose(DETECTOR, 270.0, 30.0)), GRID, Hypothesis.as_built())
    assert np.bincount(low["status"].values.ravel(), minlength=5).tolist() == [1318, 276, 6, 0, 0]
    assert float(low["zenith_deg"].max()) < 90.0
    lower = quiet(SITE, pier(), station(pose=Pose(DETECTOR, 270.0, 60.0)), GRID, Hypothesis.as_built())
    status, zenith = lower["status"].values[0, 0], lower["zenith_deg"].values[0]
    assert {DEFINED, NOT_ACCEPTED, BELOW_CUTOFF, LOW_ELEVATION} == set(status.ravel().tolist())
    assert zenith.max() > 90.0  # bin centres on the side of the normal that look below the horizon
    accepted = lower["accepted_solid_angle_sr"].values[0]
    below = zenith > 92.0
    assert np.all(accepted[below] == 0.0) and np.all(status[below] == NOT_ACCEPTED)
    assert np.all(lower["counts"].values[0, 0][below] == 0.0)
    lines = stack().lines_of_sight(Pose(DETECTOR, 270.0, 60.0), GRID)
    assert accepted.sum() == pytest.approx(lines.solid_angle_total_sr, rel=1e-13)
    assert accepted.sum() < 0.9 * low["accepted_solid_angle_sr"].values.sum()
    # the azimuth of the bin centres: the axis of the film points to the west
    j, i = 20, 20
    assert lower["azimuth_deg"].values[0, j, i] == pytest.approx(270.0, abs=3.0)


# ------------------------------------------------------------------ layout


def test_the_layout_of_the_map():
    expected = demo()
    assert isinstance(expected, xr.Dataset)
    assert set(VARIABLES) == set(expected.data_vars)
    conventions = {
        "solid_angle_sr": ("ty", "tx"),
        "opacity_gcm2": ("detector", "hypothesis", "ty", "tx"),
        "counts": ("detector", "hypothesis", "ty", "tx"),
        "transmission": ("detector", "hypothesis", "ty", "tx"),
        "reference": ("detector", "ty", "tx"),
        "zenith_deg": ("detector", "ty", "tx"),
        "azimuth_deg": ("detector", "ty", "tx"),
    }
    for name, dims in conventions.items():
        assert expected[name].dims == dims, name
    assert expected["status"].dims == expected["counts"].dims
    assert (
        expected["accepted_solid_angle_sr"].dims
        == expected["open_sky_status"].dims
        == ("detector", "ty", "tx")
    )
    assert expected["detector"].values.tolist() == ["K>1GeV", "K>0.2GeV"]
    assert expected["hypothesis"].values.tolist() == ["none", "void"]
    assert np.array_equal(expected["tx"].values, GRID.tx_centers) and expected.sizes["ty"] == 40
    assert grid_of(expected) == GRID
    # the solid angle of a bin is that of its rectangle; what a detector accepts of it is apart
    assert np.array_equal(expected["solid_angle_sr"].values, GRID.solid_angle_sr)
    accepted = expected["accepted_solid_angle_sr"].sel(detector="K>1GeV").values
    assert accepted.sum() == pytest.approx(2 * math.pi * (1 - math.cos(math.atan(1.0))), rel=1e-13)
    assert (
        np.all(accepted <= GRID.solid_angle_sr * (1 + 1e-12))
        and np.count_nonzero(accepted < 0.999 * GRID.solid_angle_sr) == 424
    )
    # the centre of a bin in ENU: a horizontal film looks up, +tx is east and +ty north
    zenith, azimuth = expected["zenith_deg"].values[0], expected["azimuth_deg"].values[0]
    assert zenith[20, 20] == pytest.approx(math.degrees(math.atan(math.hypot(0.025, 0.025))))
    assert azimuth[20, 39] == pytest.approx(math.degrees(math.atan2(0.975, 0.025)))  # east, a little north
    assert azimuth[39, 20] == pytest.approx(math.degrees(math.atan2(0.025, 0.975)))
    assert expected["exposure_days"].values.tolist() == [DAYS, DAYS]
    attrs = expected.attrs
    assert attrs["muoscene_version"] == __version__ and attrs["undefined"] == "mask"
    assert (
        attrs["mulder_version"] == mc.versions()["mulder"]
        and attrs["calzone_version"] == mc.versions()["calzone"]
    )
    stations = json.loads(attrs["stations"])
    assert stations["K>0.2GeV"]["detector"]["energy_response"] == {
        "kinetic_energy_gev": [0.2],
        "efficiency": [1.0],
    }
    assert stations["K>1GeV"]["exposure_days"] == DAYS and stations["K>1GeV"]["pose"]["position_m"] == list(
        DETECTOR
    )
    assert stations["K>1GeV"]["table"]["altitude_m"] == mc.ALTITUDE - SPIKE.PIER_NICHE
    assert json.loads(attrs["hypotheses"]) == {
        "none": {"densities_g_cm3": {}, "drho_rel": {"void": 0.0}},
        "void": {"densities_g_cm3": {}, "drho_rel": {"void": -1.0}},
    }
    assert json.loads(attrs["site"])["altitude_m"] == mc.ALTITUDE and attrs["created"].endswith("+00:00")


def test_counts_are_the_sum_of_the_detector_over_its_lines():
    """The map against the pieces it is made of, for the scene as built, above 1 GeV; 2.40e6 tracks."""
    films, pose = stack(), Pose(DETECTOR)
    lines = films.lines_of_sight(pose, GRID)
    paths = SITE.scan(pier(), pose.position_m, lines.direction)
    solid = paths.length_m[:, [i for i, b in enumerate(paths.bodies) if b.kind not in ("cut", "air")]]
    opacity = equivalent_opacity(
        table(),
        paths.opacity_by_material_g_cm2(),
        {M.register(LOAM): LOAM, M.register(CONCRETE): CONCRETE},
        1.0,
        solid_length_m=solid.sum(axis=1),
        air_density_g_cm3=pier().air_density_g_cm3,
    )
    by_hand = films.counts(table(), lines, opacity.lookup_g_cm2, DAYS)
    expected = quiet(SITE, pier(), station(), GRID, Hypothesis.as_built())
    assert expected["counts"].values[0, 0] == pytest.approx(by_hand, rel=1e-13)
    assert by_hand.sum() == pytest.approx(2.398e6, rel=1e-3)
    mean = np.bincount(lines.bin, lines.solid_angle_sr * opacity.matter_g_cm2, 1600) / np.maximum(
        np.bincount(lines.bin, lines.solid_angle_sr, 1600), 1e-300
    )
    reported = expected["opacity_gcm2"].values[0, 0].ravel()
    seen = np.isfinite(reported)
    assert reported[seen] == pytest.approx(mean[seen], rel=1e-13)  # the matter, without the air term
    assert json.loads(expected.attrs["fast_vs_direct"])["K>1GeV"] == {
        "compared": False,
        "why": "check_directions is 0",
    }
    assert math.isnan(expected.attrs["fast_vs_direct_max"])


def test_the_open_sky_above_a_site_at_sea_level():
    """A film 0.5 m above the brow at sea level: every bin is a reference bin, also above 0.2 GeV.

    No cut-off at sea level. The counts are flux × area × cos × solid angle × efficiency × time
    of the table itself: 5.29e6 tracks above 0.2 GeV in the cone of 45 degrees in 60 days.
    """
    scene = pier(site=SEA)
    pose = Pose((-0.5, 0.0, 0.5))
    sky = MulderTable(0.5, LOAM, cache=False)
    place = Station("above the brow", stack(0.2), pose, sky, DAYS)
    expected = expected_map(
        SEA, scene, place, GRID, BOTH, rng=np.random.default_rng(3), undefined="refuse", check_directions=64
    )
    status = expected["status"].values[0, 0]
    assert (
        np.all(expected["reference"].values[0] == (status == DEFINED)) and (status == DEFINED).sum() == 1324
    )
    assert np.all(expected["transmission"].values[0, 0][status == DEFINED] == 1.0)
    assert np.all(expected["opacity_gcm2"].values[0, 0][status == DEFINED] == 0.0)
    lines = stack(0.2).lines_of_sight(pose, GRID)
    flux = sky.integral(EnergyResponse.threshold(0.2), lines.direction[:, 2], 0.0)
    total = float(np.sum(flux * 120.0e-4 * lines.cos_normal * lines.solid_angle_sr) * 0.9 * SECONDS)
    assert float(expected["counts"].sum()) == pytest.approx(2 * total, rel=1e-13)  # two hypotheses, one sky
    assert total == pytest.approx(5.288e6, rel=1e-3)
    check = json.loads(expected.attrs["fast_vs_direct"])["above the brow"]
    assert (
        check["n_directions"] == 64 and check["max"] < 2e-4
    )  # met 5e-5: the energy nodes of the direct flux


def test_the_preliminary_calculation_goes_through_the_same_map():
    """A `LegacyTable` with `LEGACY_RESPONSE`: a plain sum of density times length; Mulder is not asked."""
    scene = pier(site=SEA)
    films = stack(response=LEGACY_RESPONSE)
    place = Station("legacy", films, Pose(DETECTOR), LegacyTable(), DAYS)
    expected = expected_map(SEA, scene, place, GRID, BOTH, rng=np.random.default_rng(1), undefined="refuse")
    lines = films.lines_of_sight(Pose(DETECTOR), GRID)
    paths = SEA.scan(scene, DETECTOR, lines.direction)
    by_hand = films.counts(LegacyTable(), lines, paths.opacity_g_cm2(drho_rel={"void": 0.0}), DAYS)
    assert expected["counts"].sel(hypothesis="none").values[0] == pytest.approx(by_hand, rel=1e-13)
    check = json.loads(expected.attrs["fast_vs_direct"])["legacy"]
    assert check["compared"] is False and "preliminary" in check["why"]
    assert json.loads(expected.attrs["stations"])["legacy"]["detector"]["energy_response"] == "legacy"
    with pytest.raises(ValueError, match="LEGACY_RESPONSE goes with a LegacyTable"):
        Station("mixed", films, Pose(DETECTOR), table(), DAYS)


# ------------------------------------------------------------------ what is refused


def test_what_a_station_and_a_map_refuse():
    pose = Pose(DETECTOR)
    with pytest.raises(ValueError, match="name"):
        Station("", stack(), pose, table(), DAYS)
    with pytest.raises(TypeError, match="Detector"):
        Station("a", "films", pose, table(), DAYS)
    with pytest.raises(TypeError, match="Pose"):
        Station("a", stack(), DETECTOR, table(), DAYS)
    with pytest.raises(TypeError, match="FluxTable"):
        Station("a", stack(), pose, None, DAYS)
    for days in (0.0, -3.0, math.inf, True):
        with pytest.raises(ValueError, match="exposure_days"):
            Station("a", stack(), pose, table(), days)
    # the table of a detector is at its altitude, for the atmosphere of its scene, without a date
    with pytest.raises(ValueError, match="at 1130 m and the detector at 1140"):
        quiet(SITE, pier(), station(pose=Pose((0.25, 0.0, -9.999))), GRID, BOTH)
    winter = pier.__wrapped__(site=SITE).bodies
    with pytest.raises(ValueError, match="atmosphere"):
        quiet(SITE, SITE.scene(winter, atmosphere="midlatitude-winter"), station(), GRID, BOTH)
    with pytest.raises(ValueError, match="different names"):
        quiet(SITE, pier(), [station(), station()], GRID, BOTH)
    with pytest.raises(TypeError, match="at least one Station"):
        quiet(SITE, pier(), [], GRID, BOTH)
    with pytest.raises(TypeError, match="DirectionGrid"):
        quiet(SITE, pier(), station(), None, BOTH)
    turned = Site(SITE.latitude_deg, SITE.longitude_deg, SITE.altitude_m, 90.0)
    with pytest.raises(ValueError, match="another site"):
        quiet(turned, pier(), station(), GRID, BOTH)
    with pytest.raises(ValueError, match="does not know its seed"):
        quiet(SITE, pier(), station(), GRID, BOTH, rng=np.random.default_rng([1, 2]))
    with pytest.raises(ValueError, match="check_directions"):
        expected_map(
            SITE,
            pier(),
            station(),
            GRID,
            BOTH,
            rng=np.random.default_rng(0),
            undefined="mask",
            check_directions=-1,
        )
    # a limit the fast path cannot meet refuses the run
    with pytest.raises(RuntimeError, match="above the limit"):
        expected_map(
            SITE,
            pier(),
            station(),
            GRID,
            BOTH,
            rng=np.random.default_rng(0),
            undefined="mask",
            check_directions=16,
            check_limit=1e-5,
        )


def test_example_of_the_note():
    """The example of `docs/notes/forward.md`."""
    site = Site(latitude_deg=43.35, longitude_deg=79.08, altitude_m=1150.0, scene_azimuth_deg=30.0)
    soil = M.loam(1.60, 15.0)  # bulk density 1.84 g/cm³
    bank = G.Ground("bank", [(0, -200), (0, 0), (400, 0)], -200, -400, 400, soil)  # the wall is x = 0
    niche = G.Cut("niche", G.Box((0, -0.15, -10.001), (0.5, 0.15, -9.7)))
    void = G.Anomaly("void", G.Sphere((3, 0, -3), 0.5), drho_rel=-1.0)
    scene = site.scene([bank, niche, void])

    films = EmulsionStack(
        area_cm2=120.0, n_films=6, tan_max=1.0, response=EnergyResponse.threshold(1.0), efficiency=0.9
    )
    position = tuple(site.enu_from_scene((0.25, 0.0, -10.0)))
    pose = Pose(position, azimuth_deg=120.0, tilt_deg=20.0)  # the massif lies at the azimuth 120
    place = Station("niche", films, pose, MulderTable(1140.0, soil.material), exposure_days=60)
    hypotheses = [Hypothesis.without("void"), Hypothesis.void("void")]
    expected = expected_map(
        site,
        scene,
        place,
        DirectionGrid.square(),
        hypotheses,
        rng=np.random.default_rng(1),
        undefined="refuse",
    )

    none, with_void = (expected_counts(expected, "niche", name) for name in ("none", "void"))
    assert none.sum() == pytest.approx(1.968e6, rel=1e-3)
    assert with_void.sum() - none.sum() == pytest.approx(928.0, rel=1e-2)
    assert np.sum(with_void > none) == 12
    assert (with_void / none - 1).max() == pytest.approx(0.076, abs=1e-3)
    assert expected.attrs["fast_vs_direct_max"] == pytest.approx(1.5e-4, abs=5e-5)
