"""geometry.canyon — `Canyon`: the bodies of a synthetic canyon with two banks, piers and detectors.

A generator of scenes for the examples (`examples/canyon_demo`), for the tests and for the tasks
that need one canyon to talk about. It knows nothing a `Ground`, a `Structure` and a `Cut` do not
already know: it computes their sizes from a handful of numbers and gives them back as bodies for
`Site.scene`. No size has a default; the numbers of the demonstration live in the example.

Axes of the scene (decided on 22 November 2026, roadmap 0.1-10): x across the canyon, from the
left wall to the right one; y along the canyon; z up. The left wall is the plane x = 0 with its
massif at x < 0, the right wall the plane x = `width_m` with its massif beyond. "Left" and
"right" are as seen looking along +y. The two banks are one `Ground` whose profile is a trough:
the brow of each bank, a vertical wall down to a flat floor, the wall of the other bank. Each
brow has an edge, the top of its wall, and may rise or fall linearly over a distance from the
edge before it levels off.

A place on a wall is given in the terms of the canyon, not in ENU: the bank, the distance along
the canyon, the depth below the edge of the brow of that bank and, for points in the ground,
the distance behind the wall (`Canyon.point`). `Site` turns it into ENU. Moving a detector is
one change of these numbers, and its niche moves with it.

A detector stands on a wall in one of two ways (decided on 22 November 2026):

- `Canyon.niche` — a niche is cut into the wall where the detector is: a box of air of the
  given depth into the wall, width and height, open to the canyon. The detector is in its
  middle, half the depth of the niche behind the wall, `FLOOR_GAP_M` above its floor;
- `Canyon.mount` — no niche: the detector hangs in the air of the canyon in front of the wall,
  its centre `standoff_m` from it, as on pins drilled into the rock. The pins, the cassette and
  any cover are not modelled.

A pier is a pile cap whose top is level with the surface of its brow and a grid of vertical
piles under it. The surface under the cap must be level: a cap on the slope of a brow would
cross the surface, which version 0.1 refuses (`docs/notes/geometry.md`).

What the scan needs from the sizes: every line of sight of a detector must leave the ground
through its surface before the limits of the model (`behind_m`, `half_length_m`, `z_bottom_m`),
or `Scene.scan` refuses it. With 300 m behind each wall and along the canyon all the detectors
of the demonstration pass; a film 22 m down the left wall tilted into the massif by 40 degrees
still passes, down to 5.2 degrees of elevation, and one tilted by 44 degrees is refused: its
lines below 4 degrees reach the far end of the model (`tests/unit/test_geometry_canyon.py`).

Units: metres; densities in g/cm³.

Layer 1.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from ..materials.composition import Material, Mixture, Soil
from .primitives import Body, Box, Cut, Ground, Layer, Structure, pile_grid

BANKS = ("left", "right")
"""Banks of a canyon, as seen looking along +y of the scene."""
GROUND = "ground"
"""Name of the ground of a canyon in a `PathTable`; its layers keep their own names."""
FLOOR_GAP_M = 1.0e-3
"""A detector in a niche is this much above its floor: a scan refuses a point on a boundary
(`checks.BOUNDARY_M`), and the frozen cases of the spike keep the same millimetre."""
LEVEL_M = 1.0e-9
"""The surface under a pile cap is level when it varies by less than this."""


def _number(value, what: str) -> float:
    if isinstance(value, bool) or not isinstance(value, int | float) or not math.isfinite(value):
        raise ValueError(f"{what} must be a finite number")
    return float(value)


def _positive(value, what: str) -> float:
    value = _number(value, what)
    if not value > 0:
        raise ValueError(f"{what} must be positive")
    return value


def _bank(bank, what: str) -> str:
    if bank not in BANKS:
        raise ValueError(f"{what}: the bank is {' or '.join(map(repr, BANKS))}, not {bank!r}")
    return bank


@dataclass(frozen=True)
class Bank:
    """The brow of one bank: the edge on top of its wall and how the surface goes on from it.

    edge_z_m      elevation of the edge, the top of the wall, m
    rise_m        the surface rises by this much (falls, if negative) going away from the canyon ...
    rise_over_m   ... over this distance from the edge, m, and is level beyond it; 0 with a
                  rise of 0, a level brow
    """

    edge_z_m: float
    rise_m: float
    rise_over_m: float

    def __post_init__(self) -> None:
        object.__setattr__(self, "edge_z_m", _number(self.edge_z_m, "edge_z_m"))
        object.__setattr__(self, "rise_m", _number(self.rise_m, "rise_m"))
        over = _number(self.rise_over_m, "rise_over_m")
        if self.rise_m == 0.0 and over != 0.0:
            raise ValueError("a level brow (rise_m = 0) has rise_over_m = 0")
        if self.rise_m != 0.0 and not over > 0:
            raise ValueError("a brow that rises or falls needs a positive rise_over_m")
        object.__setattr__(self, "rise_over_m", over)

    @property
    def level_z_m(self) -> float:
        """Elevation of the brow beyond its rise, m."""
        return self.edge_z_m + self.rise_m


@dataclass(frozen=True)
class Pier:
    """A pier on a brow: a pile cap level with the surface and a grid of vertical piles under it.

    name             the cap is the body `<name> cap`, the piles `<name> piles`
    bank             "left" or "right"
    along_m          axis of the pier along the canyon (y), m
    distance_m       from the wall to the near face of the cap, m
    cap_m            sizes of the cap across the canyon, along it, and its thickness, m
    piles            number of piles across and along the canyon
    pile_diameter_m  diameter of a pile, m
    pile_length_m    length of a pile under the cap, m
    pile_pitch_m     distance between the axes of neighbouring piles, both ways, m
    material         `Material`, `Mixture` or `Soil` of the cap and the piles
    density_g_cm3    bulk density; by default that of the material
    """

    name: str
    bank: str
    along_m: float
    distance_m: float
    cap_m: tuple[float, float, float]
    piles: tuple[int, int]
    pile_diameter_m: float
    pile_length_m: float
    pile_pitch_m: float
    material: Material | Mixture | Soil
    density_g_cm3: float | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.name, str) or not self.name.strip():
            raise ValueError("a pier needs a name")
        _bank(self.bank, self.name)
        object.__setattr__(self, "along_m", _number(self.along_m, f"{self.name}: along_m"))
        distance = _number(self.distance_m, f"{self.name}: distance_m")
        if distance < 0:
            raise ValueError(f"{self.name}: distance_m, from the wall to the cap, must not be negative")
        object.__setattr__(self, "distance_m", distance)
        cap = tuple(self.cap_m)
        if len(cap) != 3:
            raise ValueError(f"{self.name}: cap_m holds the size across, the size along and the thickness")
        object.__setattr__(self, "cap_m", tuple(_positive(v, f"{self.name}: cap_m") for v in cap))
        piles = tuple(self.piles)
        if len(piles) != 2 or not all(
            isinstance(n, int) and not isinstance(n, bool) and n >= 1 for n in piles
        ):
            raise ValueError(f"{self.name}: piles is two positive integers, across and along the canyon")
        object.__setattr__(self, "piles", piles)
        for field in ("pile_diameter_m", "pile_length_m", "pile_pitch_m"):
            object.__setattr__(self, field, _positive(getattr(self, field), f"{self.name}: {field}"))
        across, along, _ = self.cap_m
        for n, size, direction in ((piles[0], across, "across"), (piles[1], along, "along")):
            if (n - 1) * self.pile_pitch_m + self.pile_diameter_m > size + 1e-9:
                raise ValueError(f"{self.name}: the piles do not fit under the cap {direction} the canyon")
        if not isinstance(self.material, Material | Mixture | Soil):
            raise TypeError(
                f"{self.name}: the material is a Material, a Mixture or a Soil of muoscene.materials"
            )


@dataclass(frozen=True)
class Canyon:
    """A canyon between two banks, in the axes of a scene.

    floor_z_m      elevation of the floor of the canyon, m
    width_m        distance between the walls, m: the left wall is x = 0, the right x = width_m
    left, right    the brows of the two banks (`Bank`)
    soil           `Material`, `Mixture` or `Soil` of the ground from the surface down to the
                   first layer
    layers         `Layer`s from the top down, horizontal, through both banks
    piers          `Pier`s
    behind_m       extent of the model behind each wall, m
    half_length_m  extent of the model along the canyon on either side of y = 0, m
    z_bottom_m     flat bottom of the model, m
    niche_m        a niche: its depth into the wall, its width and its height, m

    The ground is the body `GROUND`. `bodies` gives the scene, `point`, `niche` and `mount`
    give places on and behind the walls.
    """

    floor_z_m: float
    width_m: float
    left: Bank
    right: Bank
    soil: Material | Mixture | Soil
    layers: tuple[Layer, ...]
    piers: tuple[Pier, ...]
    behind_m: float
    half_length_m: float
    z_bottom_m: float
    niche_m: tuple[float, float, float]

    def __post_init__(self) -> None:
        object.__setattr__(self, "floor_z_m", _number(self.floor_z_m, "floor_z_m"))
        for field in ("width_m", "behind_m", "half_length_m"):
            object.__setattr__(self, field, _positive(getattr(self, field), field))
        object.__setattr__(self, "z_bottom_m", _number(self.z_bottom_m, "z_bottom_m"))
        if not (isinstance(self.left, Bank) and isinstance(self.right, Bank)):
            raise TypeError("left and right are the Bank of each side")
        for name in BANKS:
            bank = getattr(self, name)
            if not bank.edge_z_m > self.floor_z_m:
                raise ValueError(f"the edge of the {name} brow must lie above the floor of the canyon")
            if not bank.rise_over_m < self.behind_m:
                raise ValueError(f"the rise of the {name} brow must end inside the model (behind_m)")
        if not self.z_bottom_m < self.floor_z_m:
            raise ValueError("z_bottom_m, the bottom of the model, must lie below the floor of the canyon")
        object.__setattr__(self, "layers", tuple(self.layers))
        object.__setattr__(self, "piers", tuple(self.piers))
        if not all(isinstance(pier, Pier) for pier in self.piers):
            raise TypeError("piers are Pier objects")
        niche = tuple(self.niche_m)
        if len(niche) != 3:
            raise ValueError("niche_m holds the depth of a niche into the wall, its width and its height")
        object.__setattr__(self, "niche_m", tuple(_positive(v, "niche_m") for v in niche))
        object.__setattr__(self, "_ground", self._make_ground())
        object.__setattr__(self, "_piers", tuple(body for pier in self.piers for body in self._pier(pier)))

    # ------------------------------------------------------------------ the ground

    @property
    def profile_xz_m(self) -> tuple[tuple[float, float], ...]:
        """The profile of the trough in the x-z plane, m, from the far end of the left bank."""
        left, right, w, behind = self.left, self.right, self.width_m, self.behind_m
        points = [(-behind, left.level_z_m)]
        if left.rise_m:
            points.append((-left.rise_over_m, left.level_z_m))
        points += [(0.0, left.edge_z_m), (0.0, self.floor_z_m), (w, self.floor_z_m), (w, right.edge_z_m)]
        if right.rise_m:
            points.append((w + right.rise_over_m, right.level_z_m))
        points.append((w + behind, right.level_z_m))
        return tuple(points)

    def _make_ground(self) -> Ground:
        return Ground(
            GROUND,
            self.profile_xz_m,
            self.z_bottom_m,
            -self.half_length_m,
            self.half_length_m,
            self.soil,
            layers=self.layers,
        )

    @property
    def ground(self) -> Ground:
        """The two banks as one `Ground`."""
        return self._ground

    # ------------------------------------------------------------------ places

    def wall_x_m(self, bank: str) -> float:
        """x of the wall of a bank, m."""
        return 0.0 if _bank(bank, "wall_x_m") == "left" else self.width_m

    def into(self, bank: str) -> float:
        """Sign of x going from the wall of a bank into its massif: −1 on the left, +1 on the right."""
        return -1.0 if _bank(bank, "into") == "left" else 1.0

    def point(
        self, bank: str, along_m: float, depth_m: float, behind_m: float = 0.0
    ) -> tuple[float, float, float]:
        """A point in the axes of the scene, m, given in the terms of the canyon.

        bank       "left" or "right"
        along_m    y, along the canyon, m
        depth_m    below the edge of the brow of that bank, m (negative: above it)
        behind_m   from the wall into the massif, m; negative in front of the wall, in the canyon
        """
        _bank(bank, "point")
        along = _number(along_m, "along_m")
        depth = _number(depth_m, "depth_m")
        behind = _number(behind_m, "behind_m")
        edge = getattr(self, bank).edge_z_m
        return (self.wall_x_m(bank) + self.into(bank) * behind, along, edge - depth)

    def _on_wall(self, bank: str, along_m: float, depth_m: float, half_width_m: float, what: str) -> None:
        if not abs(_number(along_m, f"{what}: along_m")) + half_width_m < self.half_length_m:
            raise ValueError(
                f"{what}: along_m = {along_m:g} m is outside the model, which reaches"
                f" {self.half_length_m:g} m either way along the canyon"
            )
        depth = _number(depth_m, f"{what}: depth_m")
        edge = getattr(self, bank).edge_z_m
        if not 0 < depth < edge - self.floor_z_m:
            raise ValueError(
                f"{what}: depth_m = {depth:g} m is not on the {bank} wall, which goes from the edge of its"
                f" brow down {edge - self.floor_z_m:g} m to the floor of the canyon"
            )

    def niche(
        self, name: str, bank: str, along_m: float, depth_m: float
    ) -> tuple[Cut, tuple[float, float, float]]:
        """A niche in a wall and the point of the detector in it.

        name      name of the niche in a `PathTable`
        bank      "left" or "right"
        along_m   where along the canyon, m
        depth_m   the detector below the edge of the brow of that bank, m

        Returns the `Cut` and the point, in the axes of the scene: half the depth of the niche
        behind the wall, in the middle of its width, `FLOOR_GAP_M` above its floor.
        """
        _bank(bank, name)
        deep, wide, high = self.niche_m
        self._on_wall(bank, along_m, depth_m, wide / 2, name)
        x, y, z = self.point(bank, along_m, depth_m, deep / 2)
        floor = z - FLOOR_GAP_M
        if not floor > self.floor_z_m or not floor + high < getattr(self, bank).edge_z_m:
            raise ValueError(f"{name}: the niche does not fit in the wall between the floor and the brow")
        wall = self.wall_x_m(bank)
        x_lo, x_hi = sorted((wall, wall + self.into(bank) * deep))
        return Cut(name, Box((x_lo, y - wide / 2, floor), (x_hi, y + wide / 2, floor + high))), (x, y, z)

    def mount(
        self, bank: str, along_m: float, depth_m: float, standoff_m: float
    ) -> tuple[float, float, float]:
        """The point of a detector on a wall without a niche, in the axes of the scene, m.

        bank         "left" or "right"
        along_m      where along the canyon, m
        depth_m      below the edge of the brow of that bank, m
        standoff_m   from the wall to the centre of the detector, in the air of the canyon, m
        """
        _bank(bank, "mount")
        standoff = _positive(standoff_m, "standoff_m")
        if not standoff < self.width_m:
            raise ValueError(f"standoff_m = {standoff:g} m puts the detector beyond the other wall")
        self._on_wall(bank, along_m, depth_m, 0.0, "mount")
        return self.point(bank, along_m, depth_m, -standoff)

    # ------------------------------------------------------------------ piers and bodies

    def _pier(self, pier: Pier) -> tuple[Structure, Structure]:
        across, along, thick = pier.cap_m
        near = self.wall_x_m(pier.bank) + self.into(pier.bank) * pier.distance_m
        x_lo, x_hi = sorted((near, near + self.into(pier.bank) * across))
        if not abs(pier.along_m) + along / 2 < self.half_length_m or not (
            -self.behind_m < x_lo and x_hi < self.width_m + self.behind_m
        ):
            raise ValueError(f"{pier.name}: the pile cap is outside the model")
        low, high = self._ground.surface_range_m(x_lo, x_hi)
        if math.isnan(low) or high - low > LEVEL_M:
            raise ValueError(
                f"{pier.name}: the surface under the pile cap goes from {low:g} to {high:g} m; a cap level"
                " with the surface needs a level brow there: move it further from the wall"
            )
        top = high
        bottom = top - thick - pier.pile_length_m
        if not bottom > self.z_bottom_m:
            raise ValueError(f"{pier.name}: the piles reach below the bottom of the model")
        y0, y1 = pier.along_m - along / 2, pier.along_m + along / 2
        cap = Box((x_lo, y0, top - thick), (x_hi, y1, top))
        piles = pile_grid(
            ((x_lo + x_hi) / 2, pier.along_m),
            top - thick,
            pier.pile_length_m,
            pier.pile_diameter_m / 2,
            pier.piles[0],
            pier.piles[1],
            pier.pile_pitch_m,
        )
        return (
            Structure(f"{pier.name} cap", cap, pier.material, pier.density_g_cm3),
            Structure(f"{pier.name} piles", piles, pier.material, pier.density_g_cm3),
        )

    def bodies(self, extra=()) -> tuple[Body, ...]:
        """The bodies of the scene: the ground, the cap and the piles of every pier, then `extra`.

        extra   further bodies: the niches of the detectors (`niche`), anomalies, other structures
        """
        return (self._ground, *self._piers, *tuple(extra))
