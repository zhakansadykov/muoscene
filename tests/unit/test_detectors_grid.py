"""`DirectionGrid`: exact solid angles of rectangular bins, also where the acceptance circle cuts them.

The criterion of roadmap 0.1-07: the solid angles of the bins add up to the analytic solid
angle of the acceptance, 2π (1 − cos θ_max). The numbers quoted in `docs/notes/detectors.md`
for the alternative that was not taken, the share of sub-samples inside the circle, are
computed here as well.
"""

import math

import numpy as np
import pytest
from scipy import integrate

from muoscene.detectors import (
    DEFAULT_BIN_WIDTH,
    DEFAULT_SUPERSAMPLE,
    DEFAULT_TAN_MAX,
    DirectionGrid,
    Pose,
    clipped_solid_angle_sr,
    rectangle_solid_angle_sr,
)

SUM_RTOL = 1e-13  # sums of up to 2.4e4 exact terms; met: 4e-15


def cone_sr(tan_max):
    """Solid angle of the cone with the half-angle arctan(tan_max): 2π (1 − cos θ_max)."""
    return 2 * math.pi * (1 - 1 / math.sqrt(1 + tan_max**2))


def element(ty, tx):
    """Solid angle per unit area of slope space."""
    return (1 + tx**2 + ty**2) ** -1.5


# ------------------------------------------------------------------ rectangles


def test_the_formula_of_the_conventions():
    """Ω = arctan(a b / sqrt(1 + a² + b²)) for [0, a] × [0, b]; a quadrant of slope space is π/2."""
    assert rectangle_solid_angle_sr(0.0, 1.0, 0.0, 1.0) == pytest.approx(math.atan(1 / math.sqrt(3)))
    assert rectangle_solid_angle_sr(0.0, 1.0, 0.0, 1.0) == pytest.approx(math.pi / 6, rel=1e-15)
    assert rectangle_solid_angle_sr(0.0, 1e9, 0.0, 1e9) == pytest.approx(math.pi / 2, rel=1e-8)
    assert rectangle_solid_angle_sr(-1e9, 1e9, -1e9, 1e9) == pytest.approx(2 * math.pi, rel=1e-8)
    value, _ = integrate.dblquad(element, 0.0, 0.7, 0.0, 1.9, epsabs=1e-13, epsrel=1e-13)
    assert rectangle_solid_angle_sr(0.0, 0.7, 0.0, 1.9) == pytest.approx(value, rel=1e-12)


@pytest.mark.parametrize(
    "box",
    [
        (-0.7, 0.0, -1.9, 0.0),
        (-0.7, 0.0, 0.0, 1.9),
        (-0.3, 0.9, -1.1, 0.2),
        (-2.0, -1.5, 0.4, 0.45),
        (1.0, 1.05, -3.0, -2.95),
    ],
)
def test_negative_slopes(box):
    """Four corner terms hold for slopes of either sign: against a numerical integral, and by symmetry."""
    tx0, tx1, ty0, ty1 = box
    value, _ = integrate.dblquad(element, tx0, tx1, ty0, ty1, epsabs=1e-14, epsrel=1e-13)
    ours = rectangle_solid_angle_sr(tx0, tx1, ty0, ty1)
    assert ours > 0
    assert ours == pytest.approx(value, rel=1e-10)
    assert rectangle_solid_angle_sr(-tx1, -tx0, ty0, ty1) == pytest.approx(ours, rel=1e-12)
    assert rectangle_solid_angle_sr(-tx1, -tx0, -ty1, -ty0) == pytest.approx(ours, rel=1e-12)
    assert rectangle_solid_angle_sr(ty0, ty1, tx0, tx1) == pytest.approx(ours, rel=1e-12)


def test_rectangles_add_up():
    whole = rectangle_solid_angle_sr(-0.3, 0.9, -1.1, 0.2)
    parts = (
        rectangle_solid_angle_sr(-0.3, 0.0, -1.1, 0.0)
        + rectangle_solid_angle_sr(0.0, 0.9, -1.1, 0.0)
        + rectangle_solid_angle_sr(-0.3, 0.0, 0.0, 0.2)
        + rectangle_solid_angle_sr(0.0, 0.9, 0.0, 0.2)
    )
    assert parts == pytest.approx(whole, rel=1e-15)
    grid = DirectionGrid([-0.3, 0.0, 0.4, 0.9], [-1.1, 0.0, 0.2])
    assert grid.solid_angle_sr.shape == (2, 3)
    assert grid.solid_angle_sr.sum() == pytest.approx(whole, rel=1e-15)
    assert grid.solid_angle_sr[0, 0] == pytest.approx(rectangle_solid_angle_sr(-0.3, 0.0, -1.1, 0.0))


# ------------------------------------------------------------------ the circle


def test_the_default_grid_adds_up_to_the_cone():
    """tan_max = 1, bins of 0.05: 40 × 40 bins, 148 of them cut by the circle, which touches four edges."""
    assert (DEFAULT_TAN_MAX, DEFAULT_BIN_WIDTH, DEFAULT_SUPERSAMPLE) == (1.0, 0.05, 4)
    grid = DirectionGrid.square()
    assert grid.shape == (40, 40)
    assert (grid.tx_edges[0], grid.tx_edges[-1], grid.tx_edges[20]) == (-1.0, 1.0, 0.0)
    accepted = grid.accepted_solid_angle_sr(DEFAULT_TAN_MAX)
    assert accepted.sum() == pytest.approx(cone_sr(1.0), rel=SUM_RTOL)
    assert cone_sr(1.0) == pytest.approx(1.840302, abs=1e-6)
    cut = (accepted > 0) & (accepted < grid.solid_angle_sr * (1 - 1e-12))
    assert cut.sum() == 148
    assert np.count_nonzero(accepted == 0) == 276 and np.all(accepted <= grid.solid_angle_sr * (1 + 1e-12))
    inside = np.isclose(accepted, grid.solid_angle_sr, rtol=1e-12, atol=0)
    assert inside.sum() == 1600 - 148 - 276 == 1176


@pytest.mark.parametrize(
    ("tan_max", "bin_width"),
    [
        (1.0, 0.1),
        (1.0, 0.07),  # not a multiple of the bin: the grid ends at 1.05
        (0.6, 0.05),
        (math.sqrt(15.0), 0.05),  # elevations above 14.5 degrees, the lowest of a flux table
        (math.hypot(0.3, 0.4), 0.1),  # through the corners of bins
        (0.02, 0.05),  # a circle within the four central bins
        (10.0, 0.5),
    ],
)
def test_bins_add_up_to_the_cone(tan_max, bin_width):
    grid = DirectionGrid.square(max(tan_max, bin_width), bin_width)
    assert grid.accepted_solid_angle_sr(tan_max).sum() == pytest.approx(cone_sr(tan_max), rel=SUM_RTOL)
    for n in (1, 3, 4):
        sub = grid.subsample(n, tan_max)
        assert sub.solid_angle_sr.sum() == pytest.approx(cone_sr(tan_max), rel=SUM_RTOL)


def test_a_circle_inside_one_bin_and_uneven_edges():
    grid = DirectionGrid([-2.0, -0.3, 0.5, 0.6, 3.0], [-1.0, -0.2, 0.4, 2.5])
    for tan_max in (0.1, 0.19, 0.7, 1.5, 4.0):
        accepted = grid.accepted_solid_angle_sr(tan_max)
        covered = cone_sr(tan_max) if tan_max < 0.2 else None
        if covered is not None:  # the circle lies in the bin around the axis
            assert accepted[1, 1] == pytest.approx(covered, rel=1e-14) and accepted.sum() == accepted[1, 1]
        sub = grid.subsample(5, tan_max)
        per_bin = np.bincount(sub.bin, sub.solid_angle_sr, minlength=12).reshape(3, 4)
        np.testing.assert_allclose(per_bin, accepted, rtol=0, atol=5e-14)
    np.testing.assert_allclose(grid.accepted_solid_angle_sr(100.0), grid.solid_angle_sr, rtol=1e-15)


def test_cut_rectangles_against_a_numerical_integral():
    """The contour integral against quadrature over the part of a rectangle in the circle: within 1e-11 sr.

    The integral over ty has a closed form, ty / ((1 + tx²) sqrt(1 + tx² + ty²)); the one over tx
    is done by scipy between the points where the circle meets the sides.
    """

    def inner(tx, ty):
        return ty / ((1 + tx * tx) * math.sqrt(1 + tx * tx + ty * ty))

    rng = np.random.default_rng(21)
    worst, cut = 0.0, 0
    for _ in range(200):
        tx0, tx1 = np.sort(rng.uniform(-1.5, 1.5, 2))
        ty0, ty1 = np.sort(rng.uniform(-1.5, 1.5, 2))
        t = rng.uniform(0.3, 1.5)

        def strip(tx, t=t, ty0=ty0, ty1=ty1):
            reach = math.sqrt(max(t * t - tx * tx, 0.0))
            low, high = max(ty0, -reach), min(ty1, reach)
            return inner(tx, high) - inner(tx, low) if high > low else 0.0

        lo, hi = max(tx0, -t), min(tx1, t)
        expected = 0.0
        if hi > lo:
            kinks = [s * math.sqrt(t * t - y * y) for y in (ty0, ty1) if abs(y) < t for s in (-1, 1)]
            kinks = sorted(k for k in kinks if lo < k < hi)
            expected, _ = integrate.quad(
                strip, lo, hi, points=kinks or None, epsabs=1e-13, epsrel=1e-12, limit=200
            )
        ours = float(clipped_solid_angle_sr(tx0, tx1, ty0, ty1, t))
        cut += 0 < ours < float(rectangle_solid_angle_sr(tx0, tx1, ty0, ty1)) * (1 - 1e-9)
        worst = max(worst, abs(ours - expected))
    assert cut > 100  # most of the rectangles are cut by their circle
    assert worst < 1e-11  # the tolerance asked of the quadrature; met: 4e-14 sr


def test_the_share_of_sub_samples_inside_the_circle_is_not_exact():
    """The alternative that was not taken, on the default grid.

    Sub-samples whose centre is inside the circle, each with the whole solid angle of its
    rectangle: the sum is off by +3.4e-3, −3.5e-4, +4.3e-5 and +2.0e-4 of the cone for 1, 2, 4
    and 8 sub-samples per side, and by 7.9e-4 with 4 per side for bins of 0.07.
    """

    def by_share(grid, n, tan_max):
        sub = grid.subsample(n)
        return sub.solid_angle_sr[np.hypot(sub.tx, sub.ty) <= tan_max].sum() / cone_sr(tan_max) - 1

    grid = DirectionGrid.square()
    errors = [by_share(grid, n, 1.0) for n in (1, 2, 4, 8)]
    assert errors == pytest.approx([3.36e-3, -3.47e-4, 4.26e-5, 2.01e-4], rel=0.01)
    assert by_share(DirectionGrid.square(1.0, 0.07), 4, 1.0) == pytest.approx(7.88e-4, rel=0.01)
    assert abs(grid.subsample(1, 1.0).solid_angle_sr.sum() / cone_sr(1.0) - 1) < SUM_RTOL


# ------------------------------------------------------------------ sub-samples


def test_sub_samples_lie_in_their_bin_and_in_the_circle():
    grid = DirectionGrid.square()
    n_ty, n_tx = grid.shape
    for n in (1, 2, DEFAULT_SUPERSAMPLE):
        sub = grid.subsample(n, DEFAULT_TAN_MAX)
        j, i = np.divmod(sub.bin, n_tx)
        assert np.all((sub.tx >= grid.tx_edges[i]) & (sub.tx <= grid.tx_edges[i + 1]))
        assert np.all((sub.ty >= grid.ty_edges[j]) & (sub.ty <= grid.ty_edges[j + 1]))
        assert np.all(np.hypot(sub.tx, sub.ty) <= DEFAULT_TAN_MAX) and np.all(sub.solid_angle_sr > 0)
        per_bin = np.bincount(sub.bin, sub.solid_angle_sr, minlength=n_ty * n_tx).reshape(grid.shape)
        # 5e-14 sr: the contour integral of a cut rectangle is a difference of terms of order
        # one, good to 1e-15 sr, and a bin holds up to 16 of them; met: 1.2e-14 sr, in the
        # bins where the circle touches the grid. The smallest cut bin holds 8.6e-6 sr
        np.testing.assert_allclose(per_bin, grid.accepted_solid_angle_sr(1.0), rtol=0, atol=5e-14)
    # with 4 per side: 16 in a bin inside the circle, fewer in a bin on its edge, none outside
    count = np.bincount(sub.bin, minlength=n_ty * n_tx).reshape(grid.shape)
    assert count[20, 20] == 16 and count[0, 0] == 0 and 0 < count[0, 13] < 16
    assert len(sub) == 20388


def test_sub_samples_without_a_circle():
    grid = DirectionGrid([-0.2, 0.0, 0.5], [0.0, 0.1])
    sub = grid.subsample(2)
    assert len(sub) == 8
    np.testing.assert_allclose(sub.tx, [-0.15, -0.05, 0.125, 0.375] * 2, atol=1e-16)
    np.testing.assert_allclose(sub.ty, [0.025] * 4 + [0.075] * 4, atol=1e-16)
    np.testing.assert_array_equal(sub.bin, [0, 0, 1, 1, 0, 0, 1, 1])
    per_bin = np.bincount(sub.bin, sub.solid_angle_sr)
    np.testing.assert_allclose(per_bin, grid.solid_angle_sr[0], rtol=1e-14)
    with pytest.raises(ValueError):
        sub.tx[0] = 1.0  # read-only


def test_sub_samples_integrate_a_smooth_function():
    """The integral of cos θ over the cone of 45 degrees is π sin² θ_max = π/2.

    On the default grid the error is −5.3e-5 with one direction per bin and −2.6e-6 with 4 × 4;
    the solid angle itself is exact, what is left is the function taken at one point of each
    rectangle.
    """
    grid = DirectionGrid.square()
    errors = []
    for n in (1, 4):
        sub = grid.subsample(n, 1.0)
        cosine = 1 / np.sqrt(1 + sub.tx**2 + sub.ty**2)
        errors.append(float((cosine * sub.solid_angle_sr).sum() / (math.pi / 2) - 1))
    assert errors == pytest.approx([-5.27e-5, -2.59e-6], rel=0.01)


# ------------------------------------------------------------------ the grid itself


def test_edges_centres_and_directions():
    grid = DirectionGrid([-0.2, 0.0, 0.2, 0.6], [-0.1, 0.1, 0.5])
    assert grid.shape == (2, 3)
    np.testing.assert_allclose(grid.tx_centers, [-0.1, 0.1, 0.4])
    np.testing.assert_allclose(grid.ty_centers, [0.0, 0.3])
    for array in (grid.tx_edges, grid.ty_edges, grid.solid_angle_sr):
        with pytest.raises(ValueError):
            array[0] = 9.0  # read-only
    u = grid.directions(Pose((0.0, 0.0, 0.0)))
    assert u.shape == (2, 3, 3)
    # horizontal detector with no angles: the slopes are east and north over up
    np.testing.assert_allclose(
        u[1, 2], np.array([0.4, 0.3, 1.0]) / math.sqrt(1 + 0.4**2 + 0.3**2), atol=1e-15
    )
    tilted = Pose((0.0, 0.0, 0.0), azimuth_deg=90.0, tilt_deg=30.0)
    np.testing.assert_allclose(
        DirectionGrid([-0.1, 0.1], [-0.1, 0.1]).directions(tilted)[0, 0], tilted.normal
    )
    assert grid == DirectionGrid([-0.2, 0.0, 0.2, 0.6], [-0.1, 0.1, 0.5]) and hash(grid) == hash(
        DirectionGrid([-0.2, 0.0, 0.2, 0.6], [-0.1, 0.1, 0.5])
    )
    assert grid != DirectionGrid.square() and "3 bins in tx" in repr(grid)


def test_equal_angle_grid():
    grid = DirectionGrid.equal_angle(10.0, 240)
    assert grid.shape == (240, 240)
    assert (grid.tx_edges[0], grid.tx_edges[120], grid.tx_edges[-1]) == (-10.0, 0.0, 10.0)
    steps = np.diff(np.degrees(np.arctan(grid.tx_edges)))
    np.testing.assert_allclose(steps, 2 * math.degrees(math.atan(10.0)) / 240, rtol=1e-12)
    assert grid.accepted_solid_angle_sr(10.0).sum() == pytest.approx(cone_sr(10.0), rel=SUM_RTOL)
    assert DirectionGrid.equal_angle(1.0, 3).shape == (3, 3)


@pytest.mark.parametrize(
    "edges",
    [
        ([0.0], [0.0, 1.0]),
        ([0.0, 1.0, 1.0], [0.0, 1.0]),
        ([1.0, 0.0], [0.0, 1.0]),
        ([0.0, math.nan], [0.0, 1.0]),
        ([[0.0, 1.0]], [0.0, 1.0]),
    ],
)
def test_edges_are_checked(edges):
    with pytest.raises(ValueError, match="increasing"):
        DirectionGrid(*edges)


def test_what_a_grid_refuses():
    grid = DirectionGrid.square()
    for n in (0, -1, 2.0, True):
        with pytest.raises(ValueError, match="sub-samples"):
            grid.subsample(n)
    for tan_max in (0.0, -1.0, math.inf, math.nan, "1"):
        with pytest.raises(ValueError, match="tan_max"):
            grid.accepted_solid_angle_sr(tan_max)
        with pytest.raises(ValueError, match="tan_max"):
            grid.subsample(2, tan_max)
    with pytest.raises(ValueError, match="bin_width"):
        DirectionGrid.square(1.0, 0.0)
    with pytest.raises(ValueError, match="bin_width"):
        DirectionGrid.square(1.0, 2.0)
    with pytest.raises(ValueError, match="n_bins"):
        DirectionGrid.equal_angle(1.0, 1)
