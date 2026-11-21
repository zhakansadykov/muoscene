"""`Detector.counts`: expected tracks per bin from a flux table and the opacity along the lines of sight.

Roadmap 0.1-08, the border drawn in 0.1-07: the detector calls the table with its own energy
response; the forward model brings the opacity. No Mulder here: the sky is a table of one flux.
"""

import math

import numpy as np
import pytest

from muoscene.detectors import DirectionGrid, EmulsionStack, LinesOfSight, Pose, lines_of_sight
from muoscene.flux import EnergyResponse

ORIGIN = (0.0, 0.0, 0.0)
DAYS = 60.0
SECONDS = DAYS * 86400.0
GRID = DirectionGrid.square()


def stack(**changes):
    """Six films of 120 mm × 100 mm with a threshold of 1 GeV."""
    arguments = {
        "area_cm2": 120.0,
        "n_films": 6,
        "tan_max": 1.0,
        "response": EnergyResponse.threshold(1.0),
        "efficiency": 0.9,
    }
    return EmulsionStack(**(arguments | changes))


class Flat:
    """A sky of one flux in every direction, m⁻² s⁻¹ sr⁻¹, that remembers the response it was asked with."""

    def __init__(self, flux=100.0):
        self.flux, self.asked = flux, []

    def integral(self, response, sin_elevation, opacity_g_cm2):
        self.asked.append(response)
        return np.full(np.broadcast_shapes(np.shape(sin_elevation), np.shape(opacity_g_cm2)), self.flux)


def test_counts_of_a_uniform_sky_by_hand():
    """F A ε T ∫ cos θ dΩ over the cone of 45 degrees is F A ε T π sin² θ_max: 2.20e6 tracks.

    The grid sum is 2.6e-6 low with 4 × 4 sub-samples (`test_sub_samples_integrate_a_smooth_function`).
    """
    films, sky = stack(), Flat(100.0)
    lines = films.lines_of_sight(Pose(ORIGIN), GRID)
    counts = films.counts(sky, lines, 0.0, DAYS)
    assert counts.shape == GRID.shape
    expected = 100.0 * 120.0e-4 * 0.9 * SECONDS * math.pi * 0.5
    assert counts.sum() == pytest.approx(expected, rel=1e-5)
    assert sky.asked == [films.response]  # the detector asks with its own response
    # one bin on the axis side: flux × area × cos θ × solid angle of the bin, to the size of the bin
    j, i = 20, 23
    cos = 1 / math.sqrt(1 + GRID.tx_centers[i] ** 2 + GRID.ty_centers[j] ** 2)
    one = 100.0 * 120.0e-4 * cos * GRID.solid_angle_sr[j, i] * 0.9 * SECONDS
    assert counts[j, i] == pytest.approx(one, rel=1e-3)
    # the corners of the grid are outside the circle: nothing is accepted, nothing is expected
    assert counts[0, 0] == 0.0 and counts[-1, -1] == 0.0
    assert np.count_nonzero(counts) == 1324


def test_counts_are_the_track_density_times_area_and_efficiency():
    """On the accepted lines of a film the two sums differ by the area and the efficiency alone."""
    films, sky = stack(), Flat(37.0)
    for pose in (Pose(ORIGIN), Pose(ORIGIN, 250.0, 20.0, 15.0)):
        lines = films.lines_of_sight(pose, GRID)
        density = films.track_density(sky, lines, 0.0, DAYS).tracks_per_cm2
        assert films.counts(sky, lines, 0.0, DAYS).sum() == pytest.approx(density * 120.0 * 0.9, rel=1e-13)


def test_counts_follow_the_opacity_line_by_line():
    """The table is asked once, with the opacity of every line; a bin is the sum over its lines."""

    class Slope:
        def integral(self, response, sin_elevation, opacity_g_cm2):
            return 1.0 + np.asarray(opacity_g_cm2, dtype=float)

    films = stack()
    lines = films.lines_of_sight(Pose(ORIGIN), GRID, supersample=2)
    opacity = np.arange(len(lines), dtype=float)
    counts = films.counts(Slope(), lines, opacity, DAYS)
    rate = (1.0 + opacity) * 120.0e-4 * lines.cos_normal * lines.solid_angle_sr * 0.9 * SECONDS
    by_hand = np.zeros(GRID.shape).ravel()
    np.add.at(by_hand, lines.bin, rate)
    assert counts.ravel() == pytest.approx(by_hand, rel=1e-13)


def test_a_tilted_detector_counts_only_what_is_in_front_of_it_and_above_the_horizon():
    """At a tilt of 60 degrees part of the cone looks below the horizon: those bins get nothing."""
    films, sky = stack(), Flat()
    pose = Pose(ORIGIN, 90.0, 60.0)
    lines = films.lines_of_sight(pose, GRID)
    counts = films.counts(sky, lines, 0.0, DAYS)
    below = GRID.directions(pose)[..., 2] < -0.05
    assert below.any() and np.all(counts[below] == 0.0)
    assert (
        0.5 < counts.sum() / films.counts(sky, films.lines_of_sight(Pose(ORIGIN), GRID), 0.0, DAYS).sum() < 1
    )


def test_what_counts_refuse():
    films, sky, pose = stack(), Flat(), Pose(ORIGIN, 90.0, 60.0)
    lines = films.lines_of_sight(Pose(ORIGIN), GRID)
    for days in (0.0, -1.0, math.nan, True):
        with pytest.raises(ValueError, match="exposure_days"):
            films.counts(sky, lines, 0.0, days)
    with pytest.raises(TypeError, match="LinesOfSight"):
        films.counts(sky, lines.direction, 0.0, DAYS)
    # the back side of a film belongs to its track density, not to its map
    both = lines_of_sight(pose, GRID, tan_max=1.0, both_sides=True)
    with pytest.raises(ValueError, match="acceptance"):
        films.counts(sky, both, 0.0, DAYS)
    wide = lines_of_sight(Pose(ORIGIN), DirectionGrid.square(2.0, 0.1), tan_max=2.0)
    with pytest.raises(ValueError, match="acceptance"):
        films.counts(sky, wide, 0.0, DAYS)
    bare = LinesOfSight(Pose(ORIGIN), lines.direction, lines.solid_angle_sr, tan_max=1.0)
    with pytest.raises(ValueError, match="grid"):
        films.counts(sky, bare, 0.0, DAYS)
    with pytest.raises(ValueError, match="one per line"):
        films.counts(sky, lines, np.zeros(3), DAYS)
    empty = lines.select(np.zeros(len(lines), dtype=bool))
    assert np.array_equal(films.counts(sky, empty, 0.0, DAYS), np.zeros(GRID.shape))
