"""`geometry.canyon`: the trough of two banks, its layers, piers, niches and points on the walls.

Roadmap 0.1-10. The sizes are those of `examples/canyon_demo`, invented for it (decided on
22 November 2026): a floor 40 m below the level left brow, 30 m between vertical walls, a right
brow 5 m higher that rises by 1.5 m over 8 m, layer tops at −8 and −18 m, a pier on each bank.
Paths are checked against lengths worked out by hand; nothing here asks Mulder for a flux.
"""

import dataclasses
import math

import numpy as np
import pytest

from muoscene import geometry as G
from muoscene import materials as M
from muoscene.detectors import DirectionGrid, EmulsionStack, Pose
from muoscene.flux import EnergyResponse
from muoscene.geometry.canyon import FLOOR_GAP_M, GROUND
from muoscene.site import Site

ROOT2 = math.sqrt(2.0)
SITE = Site(43.35, 79.08, 1150.0, 30.0)


def soil(dry_density_g_cm3):
    return M.loam(dry_density_g_cm3, 15.0)


LEFT_PIER = G.Pier("left pier", "left", 0.0, 8.0, (6.0, 6.0, 1.5), (3, 3), 0.8, 8.0, 2.2, M.CONCRETE)
RIGHT_PIER = G.Pier("right pier", "right", 0.0, 9.0, (6.0, 6.0, 2.0), (3, 3), 1.0, 10.0, 2.0, M.CONCRETE)
CANYON = G.Canyon(
    floor_z_m=-40.0,
    width_m=30.0,
    left=G.Bank(0.0, 0.0, 0.0),
    right=G.Bank(5.0, 1.5, 8.0),
    soil=soil(1.60),
    layers=(G.Layer("layer 1", -8.0, soil(1.70)), G.Layer("layer 2", -18.0, soil(1.80))),
    piers=(LEFT_PIER, RIGHT_PIER),
    behind_m=300.0,
    half_length_m=300.0,
    z_bottom_m=-150.0,
    niche_m=(0.5, 0.3, 0.3),
)


def scene(*extra):
    return G.Scene(CANYON.bodies(extra), altitude_m=1150.0)


def lengths(table):
    """Path length in every body by name, m, for one line of sight."""
    return {name: float(table.length_in_m(name).ravel()[0]) for name in table.names}


def test_the_profile_is_a_trough():
    assert CANYON.profile_xz_m == (
        (-300.0, 0.0),
        (0.0, 0.0),
        (0.0, -40.0),
        (30.0, -40.0),
        (30.0, 5.0),
        (38.0, 6.5),
        (330.0, 6.5),
    )
    ground = CANYON.ground
    assert ground.name == GROUND and [layer.name for layer in ground.layers] == ["layer 1", "layer 2"]
    assert (ground.y_min_m, ground.y_max_m, ground.z_bottom_m) == (-300.0, 300.0, -150.0)
    # above the floor the top of the ground and layer 1 are a piece under each bank; layer 2 goes on
    # under the floor of the canyon from one bank to the other, one piece
    owners = [piece.owner for piece in ground.pieces()]
    assert owners.count(GROUND) == owners.count("layer 1") == 2 and owners.count("layer 2") == 1


def test_the_piers_are_level_with_the_surface():
    """Caps: x from −14 to −8 m on the left, 39 to 45 m on the right, tops at the surface; 3 × 3 piles."""
    cap, piles, right_cap, right_piles = CANYON.bodies()[1:5]
    assert cap.name == "left pier cap" and cap.shapes[0].aabb() == ((-14.0, -3.0, -1.5), (-8.0, 3.0, 0.0))
    assert right_cap.shapes[0].aabb() == ((39.0, -3.0, 4.5), (45.0, 3.0, 6.5))
    left_axes = sorted({round(p.centre_xy_m[0], 9) for p in piles.shapes})
    right_axes = sorted({round(p.centre_xy_m[1], 9) for p in right_piles.shapes})
    assert len(piles.shapes) == len(right_piles.shapes) == 9
    assert left_axes == [-13.2, -11.0, -8.8] and right_axes == [-2.0, 0.0, 2.0]
    assert {(p.z_bottom_m, p.z_top_m, p.radius_m) for p in piles.shapes} == {(-9.5, -1.5, 0.4)}
    assert {(p.z_bottom_m, p.z_top_m, p.radius_m) for p in right_piles.shapes} == {(-5.5, 4.5, 0.5)}


def test_paths_across_both_walls_and_the_layers():
    """Lines at 45 degrees into each massif from a detector on the wall, 50 m along, away from the piers.

    Left, from 22 m below the brow: 0.15 m of canyon air, then 3.85, 10 and 8 m of rise in layer 2,
    layer 1 and the top of the ground. Right, from 12 m below its edge (z = −7): the top of the
    ground up to the level part of the brow at +6.5 m, 13.35 m of rise. Within 1e-9 m.
    """
    here = scene()
    left = here.scan(CANYON.mount("left", 50.0, 22.0, 0.15), np.array([-1.0, 0.0, 1.0]) / ROOT2)
    got = lengths(left)
    assert got["layer 2"] == pytest.approx(3.85 * ROOT2, abs=1e-9)
    assert got["layer 1"] == pytest.approx(10.0 * ROOT2, abs=1e-9)
    assert got[GROUND] == pytest.approx(8.0 * ROOT2, abs=1e-9)
    right = here.scan(CANYON.mount("right", 50.0, 12.0, 0.15), np.array([1.0, 0.0, 1.0]) / ROOT2)
    got = lengths(right)
    assert got[GROUND] == pytest.approx(13.35 * ROOT2, abs=1e-9)
    assert got["layer 1"] == got["layer 2"] == 0.0


@pytest.mark.parametrize(
    ("x", "expected"),
    [
        (-11.0, {"layer 2": 12.0, "layer 1": 8.5, "left pier piles": 8.0, "left pier cap": 1.5, GROUND: 0.0}),
        (
            42.0,
            {"layer 2": 12.0, "layer 1": 10.0, GROUND: 2.5, "right pier piles": 10.0, "right pier cap": 2.0},
        ),
    ],
)
def test_a_vertical_line_through_the_axis_of_a_pile(x, expected):
    """From 30 m below the origin straight up through the middle pile and the cap of each pier.

    The piles of the left pier cross the top of layer 1 at −8 m and are split there; their
    parts add up to the 8 m of a pile.
    """
    got = lengths(scene().scan((x, 0.0, -30.0), (0.0, 0.0, 1.0)))
    for name, length in expected.items():
        assert got[name] == pytest.approx(length, abs=1e-9), name


def test_a_niche_and_its_detector():
    """A niche of 0.5 × 0.3 × 0.3 m: the detector 0.25 m behind the wall and 1 mm above the floor."""
    cut, point = CANYON.niche("niche L12n", "left", -1.0, 12.0)
    assert point == (-0.25, -1.0, -12.0)
    assert cut.shape.lo_m == pytest.approx((-0.5, -1.15, -12.0 - FLOOR_GAP_M))
    assert cut.shape.hi_m == pytest.approx((0.0, -0.85, -12.0 - FLOOR_GAP_M + 0.3))
    cut, point = CANYON.niche("niche R12", "right", 50.0, 12.0)
    assert point == (30.25, 50.0, -7.0) and cut.shape.lo_m[0] == 30.0 and cut.shape.hi_m[0] == 30.5
    got = lengths(scene(cut).scan(point, (0.0, 0.0, 1.0)))
    assert got["niche R12"] == pytest.approx(0.3 - FLOOR_GAP_M, abs=1e-9)
    assert got[GROUND] == pytest.approx(5.0 + 1.5 * 0.25 / 8.0 - (-7.0 + 0.3 - FLOOR_GAP_M), abs=1e-9)


def test_points_on_the_walls():
    assert CANYON.mount("left", 1.0, 12.0, 0.15) == (0.15, 1.0, -12.0)
    assert CANYON.mount("right", 0.0, 22.0, 0.15) == (29.85, 0.0, -17.0)
    assert CANYON.point("left", 0.0, 6.0, 6.0) == (-6.0, 0.0, -6.0)  # the void of the example
    assert CANYON.point("right", 0.0, 5.0, 6.0) == (36.0, 0.0, 0.0)  # its decompaction zone
    assert CANYON.into("left") == -1.0 and CANYON.wall_x_m("right") == 30.0


def test_the_scene_of_the_example_builds_and_has_no_overlaps():
    cuts = [
        CANYON.niche(f"niche {i}", bank, along, depth)[0]
        for i, (bank, along, depth) in enumerate(
            [("left", -1.0, 12.0), ("left", -1.0, 22.0), ("right", 0.0, 12.0)]
        )
    ]
    void = G.Anomaly("void", G.Sphere(CANYON.point("left", 0.0, 6.0, 6.0), 1.0), -1.0)
    loose = G.Anomaly("decompaction", G.Sphere(CANYON.point("right", 0.0, 5.0, 6.0), 1.5), -0.2)
    built = SITE.scene(CANYON.bodies([*cuts, void, loose]))
    built.check()
    assert built.hosts == {"void": GROUND, "decompaction": GROUND}


def test_how_far_a_film_on_the_wall_can_tilt_into_the_massif():
    """A film 22 m down the left wall, 0.15 m in front of it, tilted into the massif (azimuth 300).

    At 40 degrees every line of the default grid, down to 5.2 degrees of elevation, leaves the
    ground through its surface; at 44
    degrees lines below 4 degrees of elevation reach the far end of the model and the scan
    refuses them. The flux table starts at 14.5 degrees: what the scan passes at 40 degrees,
    the table does not all have (`examples/canyon_demo/README.md`).
    """
    films = EmulsionStack(
        area_cm2=120.0, n_films=6, tan_max=1.0, response=EnergyResponse.threshold(1.0), efficiency=0.9
    )
    built = SITE.scene(CANYON.bodies())
    position = tuple(SITE.enu_from_scene(CANYON.mount("left", 0.0, 22.0, 0.15)))
    lines = films.lines_of_sight(Pose(position, 300.0, 40.0), DirectionGrid.square())
    SITE.scan(built, position, lines.direction)
    assert math.degrees(math.asin(lines.direction[:, 2].min())) == pytest.approx(5.16, abs=0.01)
    assert np.unique(lines.bin[lines.direction[:, 2] < 0.25]).size == 154  # bins below 14.5 degrees
    lines = films.lines_of_sight(Pose(position, 300.0, 44.0), DirectionGrid.square())
    with pytest.raises(
        ValueError, match=r"cross a limit of the model .* at elevations of [0-9.]+ to 3\.\d+ degrees"
    ):
        SITE.scan(built, position, lines.direction)


def test_what_a_canyon_refuses():
    with pytest.raises(ValueError, match="level brow"):
        dataclasses.replace(CANYON, piers=(dataclasses.replace(RIGHT_PIER, distance_m=3.0),))
    with pytest.raises(ValueError, match="do not fit under the cap"):
        dataclasses.replace(LEFT_PIER, pile_pitch_m=3.0)
    with pytest.raises(ValueError, match="two positive integers"):
        dataclasses.replace(LEFT_PIER, piles=(3,))
    with pytest.raises(ValueError, match="bank"):
        dataclasses.replace(LEFT_PIER, bank="middle")
    with pytest.raises(ValueError, match="below the bottom"):
        dataclasses.replace(CANYON, piers=(dataclasses.replace(LEFT_PIER, pile_length_m=200.0),))
    with pytest.raises(ValueError, match="rise_over_m = 0"):
        G.Bank(0.0, 0.0, 5.0)
    with pytest.raises(ValueError, match="positive rise_over_m"):
        G.Bank(0.0, 1.0, 0.0)
    with pytest.raises(ValueError, match="above the floor"):
        dataclasses.replace(CANYON, floor_z_m=1.0)
    with pytest.raises(ValueError, match="below the floor"):
        dataclasses.replace(CANYON, z_bottom_m=-40.0)
    with pytest.raises(ValueError, match="end inside the model"):
        dataclasses.replace(CANYON, right=G.Bank(5.0, 1.5, 300.0))
    with pytest.raises(ValueError, match="not on the left wall"):
        CANYON.niche("n", "left", 0.0, 45.0)
    with pytest.raises(ValueError, match="does not fit in the wall"):
        CANYON.niche("n", "left", 0.0, 0.2)
    with pytest.raises(ValueError, match="outside the model"):
        CANYON.niche("n", "left", 299.9, 12.0)
    with pytest.raises(ValueError, match="beyond the other wall"):
        CANYON.mount("left", 0.0, 12.0, 30.0)
    with pytest.raises(ValueError, match="positive"):
        CANYON.mount("left", 0.0, 12.0, 0.0)
    with pytest.raises(ValueError, match="not on the right wall"):
        CANYON.mount("right", 0.0, -1.0, 0.15)
    with pytest.raises(ValueError, match="bank"):
        CANYON.mount("middle", 0.0, 12.0, 0.15)
