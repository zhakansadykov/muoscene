"""Track density of `EmulsionStack` against the preliminary calculation (roadmap 0.1-07).

The frozen cases `track_density/…` and `track_density_tilted/…` of `data/legacy_reference.json`
are Monte Carlo sums of the preliminary calculation: 200 000 and 400 000 directions drawn with
the seed 5. muoscene sums the same integrand over a grid of lines of sight, in the conditions
of the preliminary calculation: `LegacyTable` with `LEGACY_RESPONSE`, straight paths, its wall
and its niche. The path in the ground comes from a scene of `muoscene.geometry` through
`Site.scan`.

Two questions are kept apart.

1. Is it the same integrand? On the directions of the preliminary calculation themselves, with
   its solid angle per direction, the density is its own to the last bits: 4e-16 with the path
   of `legacy.soil_length`, and the same through the scene under a flat brow. Under the rising
   brow the scene differs by −3.0e-6 and −3.5e-6: the preliminary calculation takes the brow as
   flat towards the canyon and along the wall, up to 0.33 m of path on these directions.
2. How far is the frozen number from the integral? By the Monte Carlo error of the frozen
   number. Its standard deviation, from 40 other seeds (100 to 139) with the same number of
   directions, is `SIGMA` below: 0.28 to 0.39 % and 0.21 to 0.22 %
   (`test_monte_carlo_error_of_the_frozen_densities`, marked `validation`). The grid sum is
   converged to 4e-6 (`test_the_grid_sum_is_converged`, marked `validation`) and lies within
   0.03 % of the frozen horizontal cases and within 0.10 % of the tilted ones, 0.1 and 0.5
   standard deviations. The tolerance is three of them.
"""

import json
import math
import pathlib

import numpy as np
import pytest

from muoscene import geometry as G
from muoscene import materials as M
from muoscene.detectors import EmulsionStack, LinesOfSight, Pose
from muoscene.flux import LEGACY_RESPONSE, LegacyTable
from muoscene.site import Site
from tests.regression import legacy, legacy_cases
from tests.regression import mulder_cases as mc

FROZEN = json.loads(
    (pathlib.Path(__file__).resolve().parent / "data" / "legacy_reference.json").read_text(encoding="utf-8")
)["cases"]
MODEL = legacy.Model()
NICHE = MODEL.niche
DAYS = 60.0
WHERE = (mc.SITE["latitude"], mc.SITE["longitude"], 0.0)  # the preliminary calculation is at sea level
HORIZON = math.sqrt(15.0)  # cos(zenith) > 0.25, the cut of `legacy.track_density`, as a tangent
# Relative standard deviation of the frozen Monte Carlo numbers: 40 seeds, 100 to 139.
SIGMA = {
    "track_density/flat/0m": 0.0028,
    "track_density/flat/12m": 0.0034,
    "track_density/flat/22m": 0.0039,
    "track_density/rising/0m": 0.0028,
    "track_density/rising/12m": 0.0035,
    "track_density/rising/22m": 0.0039,
    "track_density_tilted/flat/22m/beta0": 0.0021,
    "track_density_tilted/flat/22m/beta20": 0.0022,
    "track_density_tilted/rising/22m/beta0": 0.0021,
    "track_density_tilted/rising/22m/beta20": 0.0022,
}
STACK = EmulsionStack(
    area_cm2=240.0, n_films=1, tan_max=MODEL.tmax, response=LEGACY_RESPONSE, efficiency=MODEL.eff
)


class ClippedLegacyTable:
    """`LegacyTable` as the tilted films of the preliminary calculation use it.

    `legacy.track_density_tilted` counts directions down to cos(zenith) = 0.05 and takes the
    intensity below 0.25 at 0.25; `LegacyTable` refuses them. The opacity limit is raised
    because grazing lines cross up to 8.5e4 g/cm² there; the preliminary calculation has none.
    """

    def __init__(self):
        self.table = LegacyTable(opacity_max_g_cm2=1.0e5)

    def integral(self, response, sin_elevation, opacity_g_cm2):
        return self.table.integral(response, np.maximum(sin_elevation, 0.25), opacity_g_cm2)


def canyon(site, bank, depth):
    """The bank of legacy_cases behind a wall at x = 0 of the scene, soil of 1.8 g/cm³, with the niche."""
    rise, xr = legacy_cases.BANKS[bank]["rise"], legacy_cases.BANKS[bank]["xr"]
    brow = [(xr, rise)] if rise > 0 else []
    ground = G.Ground(
        "ground",
        [(0.0, -mc.GROUND[2]), (0.0, 0.0), *brow, (mc.GROUND[0], rise)],
        z_bottom_m=-mc.GROUND[2],
        y_min_m=-mc.GROUND[1] / 2,
        y_max_m=mc.GROUND[1] / 2,
        material=M.STANDARD_ROCK,  # a path does not depend on the material
        density_g_cm3=MODEL.rho,
    )
    niche = G.Cut(
        "niche", G.Box((0.0, -NICHE.w / 2, -depth - mc.EPS), (NICHE.xn, NICHE.w / 2, -depth + NICHE.hn))
    )
    return site.scene([ground, niche])


def density(site, bank, depth, tilt, table, **lines):
    """Track density of a film on the floor of the niche, tilted towards the massif (+x of the scene)."""
    position = site.enu_from_scene([NICHE.xd, 0.0, -depth])
    pose = Pose(tuple(position), azimuth_deg=site.scene_azimuth_deg + 90.0, tilt_deg=tilt)
    sight = STACK.track_lines(pose, **lines)
    opacity = site.scan(canyon(site, bank, depth), position, sight.direction).opacity_g_cm2()
    return STACK.track_density(table, sight, opacity, DAYS)


def test_sigma_is_given_for_every_frozen_density():
    assert sorted(SIGMA) == sorted(name for name in FROZEN if name.startswith("track_density"))


@pytest.mark.parametrize("bank", sorted(legacy_cases.BANKS))
@pytest.mark.parametrize("depth", [legacy_cases.UPPER, legacy_cases.LOWER])
def test_horizontal_films(bank, depth):
    """`track_density/…`: horizontal film in the niche, cos(zenith) > 0.25; met within 0.03 %."""
    name = f"track_density/{bank}/{depth:.0f}m"
    ours = density(
        Site(*WHERE), bank, depth, 0.0, LegacyTable(), sin_elevation_min=0.0, tan_max=HORIZON, n_bins=120
    )
    assert str(ours).count("LEGACY_RESPONSE") == 1 and ours.n_directions == 56592
    assert ours.tracks_per_cm2 == pytest.approx(FROZEN[name], rel=3 * SIGMA[name])
    assert abs(ours.tracks_per_cm2 / FROZEN[name] - 1) < 5e-4  # met: +2.6e-4 to +3.2e-4


def test_open_sky():
    """`track_density/…/0m`: no ground; exact value from the table of the preliminary calculation.

    Its intensity is constant between the nodes of cos(zenith), so the integral of I cos θ over
    the sky is a sum over the cells of the nodes: 69 444.8 tracks/cm² in 60 days. The grid sum
    is 2.1e-5 below it, the steps of the table sampled at points; the frozen Monte Carlo number
    is 69 434.6, 1.5e-4 below, half of its standard deviation.
    """
    from muoscene.flux.parametrizations import COS_GRID

    table = LegacyTable()
    step = COS_GRID[1] - COS_GRID[0]
    edges = np.clip(np.append(COS_GRID - step / 2, COS_GRID[-1] + step / 2), 0.25, 1.0)
    cells = math.pi * np.diff(edges**2)  # the integral of cos θ dΩ between two cosines
    exact = float(np.sum(table.integral(LEGACY_RESPONSE, COS_GRID, 0.0) * cells)) * 1e-4 * DAYS * 86400.0
    assert exact == pytest.approx(69444.8, abs=0.1)
    lines = STACK.track_lines(Pose((0.0, 0.0, 0.0)), sin_elevation_min=0.0, tan_max=HORIZON)
    ours = STACK.track_density(table, lines, 0.0, DAYS).tracks_per_cm2
    assert ours == pytest.approx(exact, rel=5e-5)  # met: −2.1e-5
    for bank in legacy_cases.BANKS:
        name = f"track_density/{bank}/0m"
        assert ours == pytest.approx(FROZEN[name], rel=3 * SIGMA[name])
        assert FROZEN[name] / exact - 1 == pytest.approx(-1.5e-4, abs=2e-5)


@pytest.mark.parametrize("bank", sorted(legacy_cases.BANKS))
@pytest.mark.parametrize("tilt", [0.0, 20.0])
def test_tilted_films(bank, tilt):
    """`track_density_tilted/…`: both sides of the film, tan 10 from its normal, cos(zenith) > 0.05.

    Met: −0.10 % at 0 degrees and −0.06 % at 20 degrees, with 0.21 % of Monte Carlo error.
    At 20 degrees 6.7 % of the lines look at the sky from behind the film.
    """
    name = f"track_density_tilted/{bank}/{legacy_cases.LOWER:.0f}m/beta{tilt:g}"
    ours = density(Site(*WHERE), bank, legacy_cases.LOWER, tilt, ClippedLegacyTable(), sin_elevation_min=0.05)
    assert ours.both_sides and ours.tan_max == 10.0 and ours.sin_elevation_min == 0.05
    assert ours.tracks_per_cm2 == pytest.approx(FROZEN[name], rel=3 * SIGMA[name])
    assert ours.tracks_per_cm2 / FROZEN[name] - 1 == pytest.approx(
        -1.0e-3 if tilt == 0 else -0.6e-3, abs=1.5e-4
    )


@pytest.mark.parametrize("bank", sorted(legacy_cases.BANKS))
@pytest.mark.parametrize("tilt", [0.0, 20.0])
def test_the_same_integrand_on_the_directions_of_the_preliminary_calculation(bank, tilt):
    """The directions and the solid angle of `legacy.track_density_tilted` given to the stack as they are."""
    depth, g = legacy_cases.LOWER, legacy_cases.BANKS[bank]
    n = 400000
    rng = np.random.default_rng(5)
    c = rng.uniform(0.0, 1.0, n)
    s = np.sqrt(1 - c**2)
    p = rng.uniform(0, 2 * np.pi, n)
    u = np.stack([s * np.cos(p), s * np.sin(p), c], 1)
    ca = u @ legacy.normal(tilt)
    u = u[(c > 0.05) & (np.sqrt(np.maximum(1 - ca**2, 0)) <= 10 * np.abs(ca))]
    site = Site(*WHERE)
    pose = Pose((NICHE.xd, 0.0, -depth), azimuth_deg=90.0, tilt_deg=tilt)
    lines = LinesOfSight(pose, u, np.full(len(u), 2 * np.pi / n), both_sides=True, sin_elevation_min=0.05)
    frozen = FROZEN[f"track_density_tilted/{bank}/{depth:.0f}m/beta{tilt:g}"]
    table = ClippedLegacyTable()
    by_legacy_path = MODEL.rho * 100 * legacy.soil_length(MODEL, u, depth, g["rise"], g["xr"])
    own = STACK.track_density(table, lines, by_legacy_path, DAYS).tracks_per_cm2
    assert own == pytest.approx(frozen, rel=1e-12)  # met: 4e-16
    by_scene = site.scan(canyon(site, bank, depth), pose.position_m, u).opacity_g_cm2()
    through_scene = STACK.track_density(table, lines, by_scene, DAYS).tracks_per_cm2
    if bank == "flat":
        np.testing.assert_allclose(by_scene, by_legacy_path, rtol=0, atol=1e-9 * 100 * MODEL.rho)  # 1e-9 m
        assert through_scene == pytest.approx(frozen, rel=1e-12)
    else:
        # the preliminary calculation takes the rising brow as flat towards the canyon and along the wall
        assert np.abs(by_scene - by_legacy_path).max() / (100 * MODEL.rho) == pytest.approx(0.334, abs=0.001)
        assert through_scene / frozen - 1 == pytest.approx(-3.0e-6 if tilt == 0 else -3.5e-6, abs=3e-7)


@pytest.mark.parametrize("azimuth", [37.0, 180.0])
def test_a_turned_scene_gives_the_same_density(azimuth):
    """The wall at another azimuth, the film turned with it: the same lines in the scene, the same number."""
    table = ClippedLegacyTable()
    options = {"sin_elevation_min": 0.05, "n_bins": 120}
    straight = density(Site(*WHERE), "rising", legacy_cases.LOWER, 20.0, table, **options)
    turned = density(Site(*WHERE, azimuth), "rising", legacy_cases.LOWER, 20.0, table, **options)
    assert turned.n_directions == straight.n_directions
    assert turned.tracks_per_cm2 == pytest.approx(straight.tracks_per_cm2, rel=1e-9)  # met: 2e-16


@pytest.mark.validation
def test_the_grid_sum_is_converged():
    """120, 240 and 480 bins a side: the tilted density moves by 2.6e-4 and then by 4e-6.

    1.2e6 lines of sight, 20 s on macOS: marked `validation`.
    """
    table = ClippedLegacyTable()
    values = [
        density(
            Site(*WHERE), "flat", legacy_cases.LOWER, 20.0, table, sin_elevation_min=0.05, n_bins=n
        ).tracks_per_cm2
        for n in (120, 240, 480)
    ]
    assert abs(values[0] / values[2] - 1) < 5e-4 and abs(values[1] / values[2] - 1) < 2e-5


def test_the_back_side_of_a_tilted_film():
    """What the acceptance of the conventions would leave out of the density: 0.48 % at 20°, 9.4 % at 45°.

    `legacy.track_density_tilted` counts both sides of the film; with the front side alone the
    density at 20 degrees is 0.48 % lower and the frozen number would be missed by two of its
    standard deviations.
    """
    table = ClippedLegacyTable()
    for tilt, share in ((20.0, -0.0048), (45.0, -0.094)):
        options = {"sin_elevation_min": 0.05, "n_bins": 120}
        both = density(Site(*WHERE), "flat", legacy_cases.LOWER, tilt, table, **options)
        front = density(Site(*WHERE), "flat", legacy_cases.LOWER, tilt, table, both_sides=False, **options)
        assert front.tracks_per_cm2 / both.tracks_per_cm2 - 1 == pytest.approx(share, rel=0.03)


@pytest.mark.parametrize(
    ("depth", "low_share", "canyon_share"),
    [(legacy_cases.UPPER, 0.0084, 0.775), (legacy_cases.LOWER, 0.0086, 0.875)],
)
def test_where_the_tracks_of_a_film_in_a_niche_come_from(depth, low_share, canyon_share):
    """Horizontal film, conditions of the preliminary calculation, tan 10 and cos(zenith) above 0.05.

    Elevations below 16 degrees, which a threshold of 1 GeV cannot have at 1150 m with the
    sea-level reference, hold 0.84 and 0.86 % of the density; the side of the canyon holds
    78 and 87 % of it.
    """
    site = Site(*WHERE)
    position = site.enu_from_scene([NICHE.xd, 0.0, -depth])
    lines = STACK.track_lines(Pose(tuple(position)), sin_elevation_min=0.05)
    opacity = site.scan(canyon(site, "flat", depth), position, lines.direction).opacity_g_cm2()
    table = ClippedLegacyTable()
    whole = STACK.track_density(table, lines, opacity, DAYS).tracks_per_cm2

    def share(mask):
        return STACK.track_density(table, lines.select(mask), opacity[mask], DAYS).tracks_per_cm2 / whole

    assert share(lines.direction[:, 2] <= math.sin(math.radians(16.0))) == pytest.approx(low_share, abs=1e-4)
    assert share(lines.direction[:, 0] < 0) == pytest.approx(canyon_share, abs=1e-3)


@pytest.mark.validation
def test_monte_carlo_error_of_the_frozen_densities():
    """`SIGMA`: the spread of the estimators of the preliminary calculation over the seeds 100 to 139."""
    seeds = range(100, 140)
    low = legacy.Model(area_cm2=legacy_cases.AREA_LOWER)
    for bank, g in legacy_cases.BANKS.items():
        for depth in (0.0, legacy_cases.UPPER, legacy_cases.LOWER):
            values = [
                legacy.track_density(MODEL, depth, g["rise"], g["xr"], legacy_cases.T60, seed=s)
                for s in seeds
            ]
            name = f"track_density/{bank}/{depth:.0f}m"
            assert np.std(values, ddof=1) / np.mean(values) == pytest.approx(SIGMA[name], abs=1e-4)
            assert abs(FROZEN[name] / np.mean(values) - 1) < SIGMA[name]
        for tilt in (0.0, 20.0):
            values = [
                legacy.track_density_tilted(
                    low, legacy_cases.LOWER, g["rise"], g["xr"], tilt, legacy_cases.T60, seed=s
                )
                for s in seeds
            ]
            name = f"track_density_tilted/{bank}/{legacy_cases.LOWER:.0f}m/beta{tilt:g}"
            assert np.std(values, ddof=1) / np.mean(values) == pytest.approx(SIGMA[name], abs=1e-4)
            assert abs(FROZEN[name] / np.mean(values) - 1) < SIGMA[name]
