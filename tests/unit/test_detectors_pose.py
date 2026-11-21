"""`Pose`: the axes and the slopes of a detector by the conventions ("Detector pose").

Round trips — angles to axes to angles, direction to slopes to direction — hold to round-off.
They do not catch a sign that is wrong both ways, so the axes and the slopes are also compared
with numbers worked out by hand and with the formulas of the conventions written out here.
"""

import math

import numpy as np
import pytest

from muoscene.detectors import Pose

ORIGIN = (0.0, 0.0, 0.0)


def unit(azimuth_deg, elevation_deg):
    """Unit vector in ENU: azimuth from north, clockwise (conventions, "Directions")."""
    a, e = math.radians(azimuth_deg), math.radians(elevation_deg)
    return np.array([math.cos(e) * math.sin(a), math.cos(e) * math.cos(a), math.sin(e)])


def test_zero_angles_are_enu():
    pose = Pose((1.0, 2.0, -3.0))
    np.testing.assert_array_equal(pose.axes, np.eye(3))
    assert pose.position_m == (1.0, 2.0, -3.0)
    assert (pose.azimuth_deg, pose.tilt_deg, pose.roll_deg) == (0.0, 0.0, 0.0)
    # a direction 45 degrees up towards the east has the slopes (1, 0), towards the north (0, 1)
    assert pose.slopes(unit(90.0, 45.0)) == pytest.approx((1.0, 0.0), abs=1e-15)
    assert pose.slopes(unit(0.0, 45.0)) == pytest.approx((0.0, 1.0), abs=1e-15)


def test_axes_worked_out_by_hand():
    """A film tilted by 30 degrees towards the east (azimuth 90).

    Normal: (sin 30, 0, cos 30), leaning east. x_d = (cos 90, −sin 90, 0): south, the azimuth
    90 + 90. y_d = z_d × x_d = (cos 30, 0, −sin 30): east and down, along the slope of the film.
    The vertical has tx = 0 and ty = −sin 30 / cos 30 = −tan 30; the normal has the slopes (0, 0);
    a line 30 degrees from the vertical towards the west lies 60 degrees from the normal, in
    the plane of z_d and y_d, on the side of −y_d: ty = −tan 60.
    """
    pose = Pose(ORIGIN, azimuth_deg=90.0, tilt_deg=30.0)
    half, root = 0.5, math.sqrt(3) / 2
    np.testing.assert_allclose(pose.axes, [[0, -1, 0], [root, 0, -half], [half, 0, root]], atol=1e-15)
    np.testing.assert_allclose(pose.normal, unit(90.0, 60.0), atol=1e-15)
    assert pose.slopes([0.0, 0.0, 1.0]) == pytest.approx((0.0, -math.tan(math.radians(30))), abs=1e-15)
    assert pose.slopes(pose.normal) == pytest.approx((0.0, 0.0), abs=1e-15)
    assert pose.slopes(unit(270.0, 60.0)) == pytest.approx((0.0, -math.tan(math.radians(60))), abs=1e-14)
    assert pose.slopes(unit(180.0, 45.0))[0] > 0  # towards the south: along x_d


def test_roll_is_counter_clockwise_seen_from_the_tip_of_the_normal():
    """A horizontal film rolled by 90 degrees: x_d points north, y_d west."""
    pose = Pose(ORIGIN, roll_deg=90.0)
    np.testing.assert_array_equal(pose.axes, [[0, 1, 0], [-1, 0, 0], [0, 0, 1]])
    # at zero tilt the azimuth turns the same axes clockwise: x_d at the azimuth 90 + 90, south
    np.testing.assert_array_equal(Pose(ORIGIN, azimuth_deg=90.0).axes, [[0, -1, 0], [1, 0, 0], [0, 0, 1]])
    np.testing.assert_allclose(Pose(ORIGIN, 40.0, 0.0, 40.0).axes, np.eye(3), atol=2e-16)


@pytest.mark.parametrize("seed", range(5))
def test_axes_follow_the_formulas_of_the_conventions(seed):
    rng = np.random.default_rng(seed)
    alpha, beta, gamma = rng.uniform(0, 360), rng.uniform(0, 90), rng.uniform(-180, 180)
    a, b, g = np.radians([alpha, beta, gamma])
    z = np.array([np.sin(b) * np.sin(a), np.sin(b) * np.cos(a), np.cos(b)])
    x0 = np.array([np.cos(a), -np.sin(a), 0.0])
    y0 = np.cross(z, x0)
    pose = Pose(ORIGIN, alpha, beta, gamma)
    np.testing.assert_allclose(pose.axes[2], z, atol=1e-15)
    np.testing.assert_allclose(pose.axes[0], np.cos(g) * x0 + np.sin(g) * y0, atol=1e-15)
    np.testing.assert_allclose(pose.axes[1], -np.sin(g) * x0 + np.cos(g) * y0, atol=1e-15)
    np.testing.assert_allclose(pose.axes @ pose.axes.T, np.eye(3), atol=1e-15)
    np.testing.assert_allclose(np.cross(pose.axes[0], pose.axes[1]), pose.axes[2], atol=1e-15)
    u = unit(rng.uniform(0, 360), rng.uniform(1, 89))
    tx, ty = pose.slopes(u)
    assert tx == pytest.approx((u @ pose.axes[0]) / (u @ z), rel=1e-13)
    assert ty == pytest.approx((u @ pose.axes[1]) / (u @ z), rel=1e-13)


def test_angles_to_axes_to_angles():
    """Round trip to round-off: 6e-14 degrees over 2000 poses with a tilt of at least 0.01 degrees.

    Near a tilt of zero the azimuth is ill-conditioned by 1 / sin(tilt): 1e-16 / 1.7e-4.
    """
    rng = np.random.default_rng(11)
    worst = 0.0
    for _ in range(2000):
        angles = (rng.uniform(0, 360), rng.uniform(0.01, 90), rng.uniform(-180, 180))
        pose = Pose((0.5, -1.0, 2.0), *angles)
        back = Pose.from_axes(pose.position_m, pose.axes)
        assert back.position_m == pose.position_m
        difference = np.abs(np.array([back.azimuth_deg, back.tilt_deg, back.roll_deg]) - angles)
        worst = max(worst, np.minimum(difference, 360 - difference).max())
        np.testing.assert_allclose(back.axes, pose.axes, atol=1e-15)
    assert worst < 1e-11
    for angles in ((0.0, 90.0, 0.0), (270.0, 90.0, 180.0), (45.0, 45.0, -90.0)):
        back = Pose.from_axes(ORIGIN, Pose(ORIGIN, *angles).axes)
        assert (back.azimuth_deg, back.tilt_deg, back.roll_deg) == pytest.approx(angles, abs=1e-13)


def test_a_horizontal_detector_keeps_one_angle():
    """At zero tilt the azimuth and the roll are one rotation: (α, 0, γ) comes back as (0, 0, γ − α)."""
    pose = Pose(ORIGIN, azimuth_deg=70.0, tilt_deg=0.0, roll_deg=25.0)
    back = Pose.from_axes(ORIGIN, pose.axes)
    assert (back.azimuth_deg, back.tilt_deg) == (0.0, 0.0)
    assert back.roll_deg == pytest.approx(-45.0, abs=1e-13)
    np.testing.assert_allclose(back.axes, pose.axes, atol=1e-15)


def test_direction_to_slopes_to_direction():
    """Round trip to round-off for 4000 directions of 4 poses: 3e-16 in the direction, 2e-15 in the slopes."""
    rng = np.random.default_rng(12)
    for angles in ((0.0, 0.0, 0.0), (120.0, 20.0, 0.0), (300.0, 75.0, 40.0), (10.0, 90.0, -130.0)):
        pose = Pose(ORIGIN, *angles)
        tx, ty = rng.uniform(-3, 3, (2, 1000))
        u = pose.direction(tx, ty)
        assert u.shape == (1000, 3)
        np.testing.assert_allclose(np.linalg.norm(u, axis=1), 1.0, rtol=4e-16)
        assert np.all(u @ pose.normal > 0)
        np.testing.assert_allclose(u @ pose.normal, 1 / np.sqrt(1 + tx**2 + ty**2), rtol=1e-15)
        back = pose.slopes(u)
        np.testing.assert_allclose(back[0], tx, atol=2e-14)
        np.testing.assert_allclose(back[1], ty, atol=2e-14)
        np.testing.assert_allclose(pose.direction(*back), u, atol=1e-15)
        # a pair of slopes is a line: the opposite direction and any length give the same pair
        np.testing.assert_allclose(pose.slopes(-2.5 * u), back, atol=2e-14)
    assert Pose(ORIGIN).direction(0.0, 0.0).shape == (3,)
    assert Pose(ORIGIN).direction(np.zeros((2, 3)), 0.5).shape == (2, 3, 3)


def test_the_sky_side_of_a_line():
    """Film tilted by 60 degrees towards the east; lines in the vertical east-west plane (tx = 0).

    ty = −tan 60 is the vertical. ty = +tan 40 lies 40 degrees from the normal further east:
    100 degrees from the vertical, below the horizon, so its sky side is the opposite direction,
    80 degrees from the vertical towards the west, behind the plane of the film.
    """
    pose = Pose(ORIGIN, azimuth_deg=90.0, tilt_deg=60.0)
    ty = np.array([-math.tan(math.radians(60)), 0.0, math.tan(math.radians(40))])
    u, front = pose.sky_direction(np.zeros(3), ty)
    np.testing.assert_array_equal(front, [True, True, False])
    np.testing.assert_allclose(u[0], [0, 0, 1], atol=1e-15)
    np.testing.assert_allclose(u[1], pose.normal, atol=1e-15)
    np.testing.assert_allclose(u[2], unit(270.0, 10.0), atol=1e-15)
    assert np.all(u[:, 2] > 0) and u[2] @ pose.normal < 0
    np.testing.assert_allclose(pose.direction(0.0, ty[2]), -u[2], atol=1e-15)
    # for a horizontal film every line looks at the sky from the front
    assert np.all(Pose(ORIGIN).sky_direction(np.linspace(-5, 5, 11), 0.3)[1])


def test_angles_are_kept_in_their_ranges():
    pose = Pose(ORIGIN, azimuth_deg=-90.0, tilt_deg=10, roll_deg=270.0)
    assert (pose.azimuth_deg, pose.tilt_deg, pose.roll_deg) == (270.0, 10.0, -90.0)
    assert Pose(ORIGIN, 720.0, 0.0, -180.0).roll_deg == 180.0
    assert Pose(ORIGIN, 720.0, 0.0, 540.0).azimuth_deg == 0.0
    assert Pose(ORIGIN, 1.0, 2.0, 3.0) == Pose([0, 0, 0], 361.0, 2.0, 363.0)


@pytest.mark.parametrize(
    "arguments",
    [
        ((0.0, 0.0),),
        ((0.0, 0.0, math.nan),),
        (ORIGIN, 0.0, -1.0),
        (ORIGIN, 0.0, 90.5),
        (ORIGIN, math.inf),
        (ORIGIN, 0.0, 0.0, "1"),
        (ORIGIN, True),
    ],
)
def test_a_pose_checks_its_numbers(arguments):
    with pytest.raises(ValueError):
        Pose(*arguments)


def test_what_a_pose_refuses():
    pose = Pose(ORIGIN, azimuth_deg=90.0, tilt_deg=30.0)
    with pytest.raises(AttributeError):
        pose.tilt_deg = 0.0
    with pytest.raises(ValueError):
        pose.axes[0, 0] = 2.0  # read-only
    with pytest.raises(ValueError, match="plane of the detector"):
        pose.slopes(pose.axes[0])
    with pytest.raises(ValueError, match="three finite"):
        pose.slopes([1.0, 0.0])
    with pytest.raises(ValueError, match="finite"):
        pose.direction(math.nan, 0.0)
    with pytest.raises(ValueError, match="orthonormal"):
        Pose.from_axes(ORIGIN, 2.0 * np.eye(3))
    with pytest.raises(ValueError, match="right-handed"):
        Pose.from_axes(ORIGIN, np.diag([1.0, -1.0, 1.0]))
    with pytest.raises(ValueError, match="below the horizon"):
        Pose.from_axes(ORIGIN, np.diag([1.0, -1.0, -1.0]))
    with pytest.raises(ValueError, match="three rows"):
        Pose.from_axes(ORIGIN, np.eye(2))
