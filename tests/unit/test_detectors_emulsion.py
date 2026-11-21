"""`EmulsionStack`: detector interface, lines of sight, and the track density above a threshold.

Roadmap 0.1-07: every output states its energy threshold; open sky below the cut-off of the
sea-level reference is refused, not returned as zero. The numbers that depend on Mulder are
compared with the frozen cases `track_density/1150m/…` of `tests/regression/mulder_cases.py`,
on their own lines of sight. Parity with the preliminary calculation is in
`tests/regression/test_legacy_track_density.py`.

Sizes are those the repository already has: the bank and the niche of
`tests/regression/legacy_cases.py`, the loam and the site of `mulder_cases.py`.
"""

import functools
import json
import math
import pathlib

import numpy as np
import pytest

from muoscene import geometry as G
from muoscene.detectors import (
    Detector,
    DirectionGrid,
    EmulsionStack,
    LinesOfSight,
    Pose,
    TrackDensity,
    lines_of_sight,
)
from muoscene.flux import LEGACY_RESPONSE, EnergyResponse, LegacyTable, MulderTable
from muoscene.materials.composition import Material, Mixture
from muoscene.site import Site
from tests.regression import legacy, legacy_cases
from tests.regression import mulder_cases as mc

FROZEN = json.loads(
    (pathlib.Path(mc.__file__).parent / "data" / "mulder_reference.json").read_text(encoding="utf-8")
)["cases"]
NICHE = legacy.Niche()
SITE = Site(mc.SITE["latitude"], mc.SITE["longitude"], mc.ALTITUDE)
ORIGIN = (0.0, 0.0, 0.0)
DAYS = 60.0
SECONDS = DAYS * 86400.0
ELEVATION = math.radians(mc.TRACK_ELEVATION)  # 16 degrees: 1 GeV is above the cut-off from 15.7
SIN_MIN, COT = math.sin(ELEVATION), 1 / math.tan(ELEVATION)


def spike_loam():
    """The loam of the spike: the `RegLoam` of the frozen cases as a muoscene material."""
    quartz = Material.from_formula("spike quartz", "SiO2", 2.65, i_ev="mulder")
    calcite = Material.from_formula("spike calcite", "CaCO3", 2.71, i_ev="mulder")
    alumina = Material.from_formula("spike alumina", "Al2O3", 3.95, i_ev="mulder")
    water = Material.from_formula("spike water", "H2O", 1.02, i_ev="mulder")  # Mulder's Water
    return Mixture.of("spike loam", {quartz: 0.60, calcite: 0.15, alumina: 0.12, water: 0.13})


LOAM = spike_loam()


def stack(k_min_gev=1.0, **changes):
    """A stack of six films of 120 mm × 100 mm with a sharp threshold."""
    arguments = {
        "area_cm2": 120.0,
        "n_films": 6,
        "tan_max": 1.0,
        "response": EnergyResponse.threshold(k_min_gev),
        "efficiency": 0.9,
    }
    return EmulsionStack(**(arguments | changes))


@functools.cache
def table(altitude_m):
    """Flux table of a detector at an altitude, loam of the spike; not cached on disk."""
    return MulderTable(altitude_m, LOAM, cache=False)


@functools.cache
def bank(depth):
    """The flat bank of legacy_cases at 1150 m, loam of 1.8 g/cm³ as in the frozen cases, with the niche."""
    ground = G.Ground(
        "ground",
        [(0.0, -mc.GROUND[2]), (0.0, 0.0), (mc.GROUND[0], 0.0)],
        z_bottom_m=-mc.GROUND[2],
        y_min_m=-mc.GROUND[1] / 2,
        y_max_m=mc.GROUND[1] / 2,
        material=LOAM,
        density_g_cm3=mc.SOIL / 1000.0,
    )
    niche = G.Cut(
        "niche", G.Box((0.0, -NICHE.w / 2, -depth - mc.EPS), (NICHE.xn, NICHE.w / 2, -depth + NICHE.hn))
    )
    return SITE.scene([ground, niche])


class Flat:
    """A sky of one flux in every direction, m⁻² s⁻¹ sr⁻¹, that remembers the response it was asked with."""

    def __init__(self, flux=100.0):
        self.flux, self.asked = flux, []

    def integral(self, response, sin_elevation, opacity_g_cm2):
        self.asked.append(response)
        return np.full(np.broadcast_shapes(np.shape(sin_elevation), np.shape(opacity_g_cm2)), self.flux)


# ------------------------------------------------------------------ the detector


def test_a_stack_has_no_defaults():
    """Every argument is given by name; the energy response above all."""
    with pytest.raises(TypeError, match="response"):
        EmulsionStack(area_cm2=120.0, n_films=6, tan_max=1.0, efficiency=0.9)
    with pytest.raises(TypeError):
        EmulsionStack(120.0, 6, 1.0, EnergyResponse.threshold(1.0), 0.9)
    for response in (None, 1.0, "1 GeV", (1.0,)):
        with pytest.raises(TypeError, match="no default threshold"):
            stack(response=response)
    assert stack(response=LEGACY_RESPONSE).response is LEGACY_RESPONSE
    ramp = EnergyResponse((0.2, 1.0), (0.0, 1.0))
    assert stack(response=ramp).response == ramp
    s = stack()
    assert isinstance(s, Detector) and s.response.k_min_gev == 1.0 and s.n_films == 6
    with pytest.raises(AttributeError):
        s.area_cm2 = 1.0


@pytest.mark.parametrize(
    "changes",
    [
        {"area_cm2": 0.0},
        {"area_cm2": math.nan},
        {"n_films": 0},
        {"n_films": 6.0},
        {"n_films": True},
        {"tan_max": 0.0},
        {"tan_max": math.inf},
        {"efficiency": 0.0},
        {"efficiency": 1.1},
        {"efficiency": "0.9"},
    ],
)
def test_a_stack_checks_its_numbers(changes):
    with pytest.raises(ValueError):
        stack(**changes)


def test_projected_area_is_the_area_times_the_cosine():
    """Conventions: the projected sensitive area is A (u · z_d). 120 cm² at 60 degrees: 60 cm²."""
    s = stack()
    assert s.sensitive_area_cm2(0.3, -0.2) == 120.0
    assert s.projected_area_cm2(0.0, 0.0) == 120.0
    assert s.projected_area_cm2(math.sqrt(3.0), 0.0) == pytest.approx(60.0, rel=1e-15)
    pose = Pose(ORIGIN, azimuth_deg=200.0, tilt_deg=35.0, roll_deg=10.0)
    tx, ty = np.random.default_rng(1).uniform(-2, 2, (2, 50))
    u = pose.direction(tx, ty)
    np.testing.assert_allclose(s.projected_area_cm2(tx, ty), 120.0 * (u @ pose.normal), rtol=1e-14)
    assert s.sensitive_area_cm2(tx, ty).shape == (50,)
    # the acceptance is a circle in the slopes, apart from the area
    assert s.accepts(0.6, 0.8) and not s.accepts(0.6, 0.81) and s.accepts(0.0, -1.0)
    np.testing.assert_array_equal(s.accepts(np.array([0.0, 2.0]), 0.0), [True, False])


def test_a_detector_whose_area_depends_on_the_slopes():
    """The interface holds a telescope of two planes a × b at the distance D (roadmap 0.4).

    A track crosses both planes if it enters the first within (a − D|tx|) × (b − D|ty|). With
    D = 0 the geometric factor within the cone is that of a plane, a b π sin² θ_max; with
    D > 0 it is compared with a numerical integral of A_eff / (1 + tx² + ty²)².
    """
    from scipy import integrate

    class Telescope(Detector):
        def __init__(self, a_cm, b_cm, base_cm, tan_max):
            self.a_cm, self.b_cm, self.base_cm = a_cm, b_cm, base_cm
            self.tan_max, self.response, self.efficiency = tan_max, EnergyResponse.threshold(1.0), 1.0
            self._check_detector()

        def sensitive_area_cm2(self, tx, ty):
            along = np.maximum(self.a_cm - self.base_cm * np.abs(tx), 0.0)
            across = np.maximum(self.b_cm - self.base_cm * np.abs(ty), 0.0)
            return along * across

    def factor(telescope, grid, n):
        lines = telescope.lines_of_sight(Pose(ORIGIN), grid, supersample=n)
        return float(np.sum(telescope.projected_area_cm2(lines.tx, lines.ty) * lines.solid_angle_sr))

    grid = DirectionGrid.square(1.0, 0.05)
    plane = Telescope(60.0, 40.0, 0.0, 1.0)
    assert factor(plane, grid, 4) == pytest.approx(60 * 40 * math.pi / 2, rel=5e-6)  # met: 2.6e-6
    telescope = Telescope(60.0, 40.0, 30.0, 2.5)  # sees up to tan 2 along tx and 4/3 along ty
    assert telescope.sensitive_area_cm2(0.0, 0.0) == 2400.0 and telescope.sensitive_area_cm2(2.0, 0.0) == 0.0
    expected, _ = integrate.dblquad(
        lambda ty, tx: (60 - 30 * tx) * (40 - 30 * ty) / (1 + tx * tx + ty * ty) ** 2,
        0.0,
        2.0,
        0.0,
        4.0 / 3.0,
    )
    wide = DirectionGrid.square(2.5, 0.05)
    assert factor(telescope, wide, 4) == pytest.approx(4 * expected, rel=1e-4)  # met: 2.9e-5


# ------------------------------------------------------------------ lines of sight


def test_accepted_lines_of_a_horizontal_stack():
    s = stack()
    grid = DirectionGrid.square()
    lines = s.lines_of_sight(Pose((1.0, 2.0, -3.0)), grid)
    assert isinstance(lines, LinesOfSight) and len(lines) == 20388 and lines.supersample == 4
    assert lines.solid_angle_total_sr == pytest.approx(2 * math.pi * (1 - math.sqrt(0.5)), rel=1e-13)
    assert np.all(lines.front) and not lines.both_sides and lines.tan_max == 1.0 and lines.grid is grid
    assert lines.direction[:, 2].min() >= math.sqrt(0.5) * (1 - 1e-12)  # the cone of 45 degrees
    np.testing.assert_allclose(lines.cos_normal, lines.direction[:, 2], rtol=1e-15)
    np.testing.assert_allclose(lines.direction[:, 0] / lines.direction[:, 2], lines.tx, atol=1e-15)
    per_bin = np.bincount(lines.bin, lines.solid_angle_sr, minlength=1600).reshape(40, 40)
    np.testing.assert_allclose(per_bin, grid.accepted_solid_angle_sr(1.0), rtol=0, atol=5e-14)
    for array in (lines.direction, lines.solid_angle_sr, lines.tx, lines.cos_normal, lines.bin):
        with pytest.raises(ValueError):
            array[0] = 0  # read-only
    assert "20388 lines" in repr(lines)


def test_lines_of_a_tilted_stack_and_its_back_side():
    """Tilt 30 degrees, cone 45: every accepted line is above the horizon. Tilt 60: some are not.

    At a tilt of 60 degrees the lines more than 30 degrees from the normal on its low side have
    their sky side behind the plane of the film. The acceptance leaves them out; `both_sides`
    keeps them, with the direction that looks at the sky.
    """
    s = stack()
    grid = DirectionGrid.square()
    cone = 2 * math.pi * (1 - math.sqrt(0.5))
    gentle = s.lines_of_sight(Pose(ORIGIN, azimuth_deg=90.0, tilt_deg=30.0), grid)
    assert gentle.solid_angle_total_sr == pytest.approx(cone, rel=1e-13) and np.all(
        gentle.direction[:, 2] > 0
    )
    steep = Pose(ORIGIN, azimuth_deg=90.0, tilt_deg=60.0)
    front = s.lines_of_sight(steep, grid)
    both = lines_of_sight(steep, grid, tan_max=1.0, both_sides=True)
    assert both.solid_angle_total_sr == pytest.approx(cone, rel=1e-13)  # every line has a sky side
    assert np.all(front.front) and np.all(front.direction @ steep.normal > 0)
    back = both.select(~both.front)
    assert len(both) == len(front) + len(back) and len(back) > 0
    assert np.all(back.direction[:, 2] > 0) and np.all(back.direction @ steep.normal < 0)
    assert np.all(back.direction[:, 0] < 0)  # the film leans east: its back side looks west
    assert np.all(back.cos_normal > 0) and back.both_sides
    # the part of the cone that dips below the horizon: 0.222 sr of its 1.840 sr, by whole sub-samples
    assert back.solid_angle_total_sr == pytest.approx(0.222, abs=0.001)
    assert front.solid_angle_total_sr == pytest.approx(cone - back.solid_angle_total_sr, rel=1e-13)
    # an elevation limit leaves out whole sub-samples
    high = lines_of_sight(steep, grid, tan_max=1.0, sin_elevation_min=0.5)
    assert 0 < len(high) < len(front) and high.direction[:, 2].min() > 0.5 and high.sin_elevation_min == 0.5


def test_lines_given_directly_and_what_is_refused():
    pose = Pose(ORIGIN, azimuth_deg=90.0, tilt_deg=60.0)
    up = np.array([[0.0, 0.0, 1.0], [0.6, 0.0, 0.8]])
    lines = LinesOfSight(pose, up, [0.1, 0.2])
    assert len(lines) == 2 and lines.solid_angle_total_sr == pytest.approx(0.3) and lines.bin is None
    assert lines.tx == pytest.approx([0.0, 0.0], abs=1e-15)
    assert lines.select(np.array([False, True])).solid_angle_total_sr == 0.2
    west = np.array([[-0.8, 0.0, 0.6]])  # 67 degrees from the normal on its low side... behind the film
    with pytest.raises(ValueError, match="both_sides"):
        LinesOfSight(pose, west, [0.1])
    assert not LinesOfSight(pose, west, [0.1], both_sides=True).front[0]
    for bad, weights, options, message in (
        (up * 2.0, [0.1, 0.2], {}, "unit vectors"),
        (-up, [0.1, 0.2], {}, "towards the sky"),
        (up, [0.1, 0.2], {"sin_elevation_min": 0.9}, "towards the sky"),
        (up, [0.1, -0.2], {}, "positive"),
        (up, [0.1], {}, "shape"),
        (up, [0.1, 0.2], {"tan_max": 0.5}, "outside tan_max"),
        (up, [0.1, 0.2], {"bin": [0, 1]}, "needs a grid"),
        (up, [0.1, 0.2], {"sin_elevation_min": 1.0}, "sin_elevation_min"),
    ):
        with pytest.raises(ValueError, match=message):
            LinesOfSight(pose, bad, weights, **options)
    with pytest.raises(TypeError):
        LinesOfSight("pose", up, [0.1, 0.2])
    with pytest.raises(TypeError):
        lines_of_sight(pose, "grid", tan_max=1.0)
    with pytest.raises(ValueError, match="mask"):
        lines.select([1, 0])


# ------------------------------------------------------------------ track density


def test_track_density_of_a_uniform_sky_by_hand():
    """Flux F in every direction, horizontal film, cone θ_max: F T π sin² θ_max per unit area.

    With F = 100 m⁻² s⁻¹ sr⁻¹, 60 days and the cone of the default track lines (tan 10):
    100 × 1e-4 × 5 184 000 s × π × 100/101 = 161 247.7 tracks/cm².
    """
    sky = Flat(100.0)
    s = stack()
    lines = s.track_lines(Pose(ORIGIN), sin_elevation_min=0.0)
    density = s.track_density(sky, lines, 0.0, DAYS)
    assert isinstance(density, TrackDensity)
    by_hand = 100.0 * 1e-4 * SECONDS * math.pi * 100 / 101
    assert by_hand == pytest.approx(161247.7, abs=0.05)
    assert density.tracks_per_cm2 == pytest.approx(by_hand, rel=1e-4)  # met: −2.9e-5 with 2.3e5 lines
    # neither the area nor the efficiency of the stack enters: the film holds every track
    assert stack(area_cm2=600.0, efficiency=0.5).track_density(sky, lines, 0.0, DAYS) == density
    assert s.track_density(sky, lines, 0.0, 30.0).tracks_per_cm2 == pytest.approx(by_hand / 2, rel=1e-4)
    assert s.track_density(sky, lines, np.zeros(len(lines)), DAYS) == density


@pytest.mark.parametrize("tilt", [0.0, 20.0, 45.0, 90.0])
def test_a_film_has_two_sides(tilt):
    """Under a uniform sky a film collects the same density at any tilt, if both its sides count.

    The integral of |u · z_d| over the upper hemisphere is π for every normal, π 100/101 within
    tan 10 of it. From the front side alone it is less by the sky behind the film.
    """
    sky = Flat(100.0)
    s = stack()
    pose = Pose(ORIGIN, azimuth_deg=130.0, tilt_deg=tilt)
    both = s.track_density(sky, s.track_lines(pose, sin_elevation_min=0.0), 0.0, DAYS)
    by_hand = 100.0 * 1e-4 * SECONDS * math.pi * 100 / 101
    assert both.tracks_per_cm2 == pytest.approx(
        by_hand, rel=1e-4
    )  # met: 2.9e-5; the horizon cuts whole sub-samples
    assert both.both_sides
    front = s.track_density(sky, s.track_lines(pose, sin_elevation_min=0.0, both_sides=False), 0.0, DAYS)
    assert not front.both_sides
    if tilt == 0.0:
        assert front.tracks_per_cm2 == both.tracks_per_cm2
    else:
        assert front.tracks_per_cm2 < both.tracks_per_cm2
        # for a vertical film exactly half of the sky is behind it
        assert (tilt != 90.0) or front.tracks_per_cm2 == pytest.approx(by_hand / 2, rel=1e-3)


def test_every_density_states_its_threshold():
    """The table is asked with the response of the stack, and the result carries it."""
    sky = Flat()
    for response in (
        EnergyResponse.threshold(0.2),
        EnergyResponse.threshold(1.0),
        EnergyResponse((0.2, 1.0), (0.0, 1.0)),
    ):
        s = stack(response=response)
        lines = s.track_lines(Pose(ORIGIN, 90.0, 20.0), sin_elevation_min=SIN_MIN, n_bins=24)
        density = s.track_density(sky, lines, 0.0, DAYS)
        assert sky.asked[-1] is response and density.response is response
        assert density.exposure_days == DAYS and density.sin_elevation_min == SIN_MIN
        assert density.tan_max == 10.0 and density.n_directions == len(lines)
        assert density.solid_angle_sr == pytest.approx(lines.solid_angle_total_sr)
    text = str(stack(1.0).track_density(sky, lines, 0.0, DAYS))
    assert "K > 1 GeV" in text and "60 days" in text and "above 16.0 degrees" in text and "both sides" in text
    assert "K > 0.2 GeV" in str(stack(0.2).track_density(sky, lines, 0.0, DAYS))
    legacy_text = str(stack(response=LEGACY_RESPONSE).track_density(LegacyTable(), lines, 0.0, DAYS))
    assert "LEGACY_RESPONSE" in legacy_text


def test_what_a_track_density_refuses():
    s = stack()
    lines = s.track_lines(Pose(ORIGIN), sin_elevation_min=SIN_MIN, n_bins=12)
    sky = Flat()
    for days in (0.0, -1.0, math.nan, "60"):
        with pytest.raises(ValueError, match="exposure_days"):
            s.track_density(sky, lines, 0.0, days)
    with pytest.raises(ValueError, match="opacity_g_cm2"):
        s.track_density(sky, lines, np.zeros(3), DAYS)
    with pytest.raises(TypeError, match="LinesOfSight"):
        s.track_density(sky, lines.direction, 0.0, DAYS)
    with pytest.raises(ValueError, match="sin_elevation_min"):
        s.track_lines(Pose(ORIGIN), sin_elevation_min=1.0)
    with pytest.raises(TypeError):
        s.track_lines(Pose(ORIGIN))  # the elevation limit is always given
    with pytest.raises(
        ValueError, match="no detector energy response"
    ):  # a threshold is not for a LegacyTable
        s.track_density(LegacyTable(), lines, 0.0, DAYS)


# ------------------------------------------------------------------ with Mulder, at 1150 m


def test_open_sky_below_the_cutoff_is_refused():
    """At 1150 m the sea-level reference has a cut-off: 0.28 GeV by the vertical, 0.98 GeV at 15.7 degrees.

    A threshold of 0.2 GeV is below it in the whole open sky, one of 1 GeV below 15.7 degrees.
    The density is then refused with the cut-off in the message; it is never zero.
    """
    sky = table(mc.ALTITUDE)
    pose = Pose(ORIGIN)
    low = stack(0.2)
    with pytest.raises(ValueError) as error:
        low.track_density(
            sky, low.track_lines(pose, sin_elevation_min=SIN_MIN, tan_max=COT, n_bins=24), 0.0, DAYS
        )
    message = str(error.value)
    assert (
        "no flux for the energy response" in message and "the detector threshold 0.2 GeV is below" in message
    )
    assert "cut-off" in message and "0.1-05b" in message and "Nothing is put in its place" in message
    high = stack(1.0)
    with pytest.raises(ValueError, match="the detector threshold 1 GeV is below"):  # down to 15 degrees
        lines = high.track_lines(pose, sin_elevation_min=math.sin(math.radians(15.0)), n_bins=48)
        high.track_density(sky, lines, 0.0, DAYS)
    with pytest.raises(ValueError, match="outside the table"):  # the table itself ends at 14.5 degrees
        high.track_density(sky, high.track_lines(pose, sin_elevation_min=0.0, n_bins=24), 0.0, DAYS)
    # from 16 degrees a threshold of 1 GeV holds
    lines = high.track_lines(pose, sin_elevation_min=SIN_MIN, tan_max=COT)
    density = high.track_density(sky, lines, 0.0, DAYS)
    assert density.tracks_per_cm2 == pytest.approx(60039.0, rel=1e-4)
    assert lines.solid_angle_total_sr == pytest.approx(2 * math.pi * (1 - SIN_MIN), rel=1e-13)
    # at sea level there is no cut-off, and the lower threshold holds too
    sea = table(0.0)
    assert high.track_density(sea, lines, 0.0, DAYS).tracks_per_cm2 == pytest.approx(49843.0, rel=1e-4)
    assert low.track_density(sea, lines, 0.0, DAYS).tracks_per_cm2 == pytest.approx(64431.0, rel=1e-4)


@pytest.mark.parametrize("depth", [legacy_cases.UPPER, legacy_cases.LOWER])
def test_track_density_in_a_niche_against_mulder(depth):
    """The frozen cases `track_density/1150m/…`: Mulder through the scene, on the same lines of sight.

    Horizontal film on the floor of the niche, 60 days. Tolerance 0.15 %: a table at the altitude
    of the detector has air where the soil above it is, −0.05 to −0.10 % (`docs/notes/flux.md`),
    and the scene of the frozen cases has 20 m of air at sea-level density above the brow,
    −0.03 % in the open sky (`docs/notes/geometry.md`). Met: −0.04 and −0.06 % on the side of
    the ground, +0.03 % on the side of the canyon.
    """
    frozen = FROZEN[f"track_density/{mc.ALTITUDE:.0f}m/{depth:.0f}m"]
    u, d_omega = mc.sky_grid()
    pose = Pose((NICHE.xd, 0.0, -depth))
    opacity = SITE.scan(bank(depth), pose.position_m, u).opacity_g_cm2()
    ground = u[:, 0] > 0
    sky = table(mc.ALTITUDE - depth)

    def density(k_min, mask):
        lines = LinesOfSight(pose, u[mask], np.full(mask.sum(), d_omega), sin_elevation_min=SIN_MIN)
        return stack(k_min).track_density(sky, lines, opacity[mask], DAYS).tracks_per_cm2

    assert density(0.2, ground) == pytest.approx(frozen["ground/K>0.2GeV"], rel=1.5e-3)
    assert density(1.0, ground) == pytest.approx(frozen["ground/K>1GeV"], rel=1.5e-3)
    assert density(1.0, ~ground) == pytest.approx(frozen["canyon/K>1GeV"], rel=1.5e-3)
    assert density(1.0, ground) / frozen["ground/K>1GeV"] - 1 < -2e-4  # the table is the lower one
    # the side of the canyon above 0.2 GeV does not exist with a sea-level reference
    with pytest.raises(ValueError, match="cut-off"):
        density(0.2, ~ground)


@pytest.mark.parametrize(
    ("depth", "total", "ground_1", "ground_02", "canyon_1"),
    [
        (legacy_cases.UPPER, 34758.0, 8514.0, 10081.0, 26244.0),
        (legacy_cases.LOWER, 30346.0, 4242.0, 4794.0, 26104.0),
    ],
)
def test_track_density_in_a_niche_at_1150_m(depth, total, ground_1, ground_02, canyon_1):
    """The numbers of `docs/numbers.md`: horizontal film in the niche, 60 days, elevations above 16 degrees.

    2.2e5 lines of sight (240 bins a side in equal angles, 2 × 2 sub-samples); with 480 bins a
    side the numbers move by less than 3e-5. The side of the canyon holds 76 % and 86 % of the
    tracks above 1 GeV: the density of a film in a niche is mostly open sky. Above 0.2 GeV only
    the side of the ground exists until roadmap 0.1-05b.
    """
    pose = Pose((NICHE.xd, 0.0, -depth))
    sky = table(mc.ALTITUDE - depth)
    lines = stack().track_lines(pose, sin_elevation_min=SIN_MIN, tan_max=COT)
    assert len(lines) == 224820
    opacity = SITE.scan(bank(depth), pose.position_m, lines.direction).opacity_g_cm2()
    ground = lines.direction[:, 0] > 0
    assert opacity[ground].min() > 2000.0  # g/cm²: no cut-off is left under 138 to 569 g/cm²
    whole = stack(1.0).track_density(sky, lines, opacity, DAYS)
    assert whole.tracks_per_cm2 == pytest.approx(total, rel=1e-4)
    assert "K > 1 GeV" in str(whole) and "above 16.0 degrees" in str(whole)
    towards_ground, towards_canyon = lines.select(ground), lines.select(~ground)
    in_ground = stack(1.0).track_density(sky, towards_ground, opacity[ground], DAYS).tracks_per_cm2
    in_canyon = stack(1.0).track_density(sky, towards_canyon, opacity[~ground], DAYS).tracks_per_cm2
    assert in_ground == pytest.approx(ground_1, rel=1e-4) and in_canyon == pytest.approx(canyon_1, rel=1e-4)
    assert in_ground + in_canyon == pytest.approx(whole.tracks_per_cm2, rel=1e-12)
    soft = stack(0.2).track_density(sky, towards_ground, opacity[ground], DAYS)
    assert soft.tracks_per_cm2 == pytest.approx(ground_02, rel=1e-4) and soft.response.k_min_gev == 0.2
    with pytest.raises(ValueError, match="cut-off"):
        stack(0.2).track_density(sky, lines, opacity, DAYS)


def test_example_of_the_note():
    """The example of `docs/notes/detectors.md`: a niche 10 m below a brow at 1150 m, a turned scene."""
    from muoscene import materials as M

    site = Site(latitude_deg=43.35, longitude_deg=79.08, altitude_m=1150.0, scene_azimuth_deg=30.0)
    soil = M.loam(1.60, 15.0)
    ground = G.Ground("bank", [(0, -200), (0, 0), (400, 0)], -200, -400, 400, soil)
    niche = G.Cut("niche", G.Box((0, -0.15, -10.001), (0.5, 0.15, -9.7)))
    scene = site.scene([ground, niche])
    films = stack()
    pose = Pose(tuple(site.enu_from_scene((0.25, 0.0, -10.0))))
    assert pose.position_m == pytest.approx((0.2165, -0.125, -10.0), abs=5e-5)
    lines = films.track_lines(pose, sin_elevation_min=SIN_MIN, tan_max=COT)
    paths = site.scan(scene, pose.position_m, lines.direction)
    sky = MulderTable(1140.0, soil.material, cache=False)
    density = films.track_density(sky, lines, paths.opacity_g_cm2(), exposure_days=60)
    assert density.tracks_per_cm2 == pytest.approx(36212.0, rel=1e-4)
    assert str(density).startswith(
        "3.621e+04 tracks/cm² in 60 days, K > 1 GeV, elevations above 16.0 degrees"
    )
    massif = site.scene_from_enu(lines.direction)[:, 0] > 0
    under_soil = lines.select(massif)
    soft = stack(0.2)
    with pytest.raises(ValueError, match="cut-off of 0.2427 GeV"):
        soft.track_density(sky, lines, paths.opacity_g_cm2(), exposure_days=60)
    assert soft.track_density(
        sky, under_soil, paths.opacity_g_cm2()[massif], 60
    ).tracks_per_cm2 == pytest.approx(11941.0, rel=1e-4)
    assert films.track_density(
        sky, under_soil, paths.opacity_g_cm2()[massif], 60
    ).tracks_per_cm2 == pytest.approx(9947.0, rel=1e-4)
    assert films.projected_area_cm2(0.5, 0.0) == pytest.approx(107.3, abs=0.05)
