"""detectors.base — `Pose`, `LinesOfSight` and the interface of a detector.

`Pose` is the reference point of a detector in ENU and its three angles, with the axes and the
slopes of the conventions ("Detector pose"):

    z_d = (sin β sin α, sin β cos α, cos β)      normal; α azimuth, β tilt
    x_d = (cos α, −sin α, 0)                      horizontal, at the azimuth α + 90°
    y_d = z_d × x_d
    tx = (u · x_d) / (u · z_d),   ty = (u · y_d) / (u · z_d)

and a roll γ that turns x_d and y_d about z_d, counter-clockwise seen from its tip.

`LinesOfSight` are directions in ENU, from the detector towards the sky, each with the solid
angle it stands for: what a scene is scanned along and what a flux is summed over.
`lines_of_sight` makes them from a `DirectionGrid` and a pose. A pair of slopes is a line, and
its sky side is the direction with u_z > 0. For a tilted detector that direction can lie
behind its plane (u · z_d < 0): the muon then crosses the detector from the back. Such lines
are outside the acceptance of the conventions and are left out, unless `both_sides` is asked
for, as the track density of an emulsion film does.

`Detector` is what the forward model needs from any detector: the circle of accepted slopes,
the energy response, the efficiency, and the sensitive area a track of given slopes can cross,
A for a film and max(a − D|tx|, 0) · max(b − D|ty|, 0) for a telescope of two planes a × b at
the distance D. The projected area is that area times cos θ = 1 / sqrt(1 + tx² + ty²).

Where the task ends (decided on 21 November 2026): a detector takes a flux table and the
opacity along its lines of sight and calls the table with its own energy response, so that
every number it gives belongs to a stated threshold. Where the opacity comes from — the scan,
the density hypotheses, the equivalent opacity — is the forward model, roadmap 0.1-08.

Units: metres and degrees for a pose, cm² for areas, g/cm² for opacity, m⁻² s⁻¹ sr⁻¹ for flux.
"""

from __future__ import annotations

import abc
import math
from dataclasses import dataclass

import numpy as np

from ..flux import LEGACY_RESPONSE, EnergyResponse
from ..site.frames import cos_sin_deg
from .grid import DEFAULT_SUPERSAMPLE, DirectionGrid

ORTHONORMAL = 1e-9
"""Axes given to `Pose.from_axes` must be orthonormal and right-handed within this."""
LEVEL = 1e-12
"""A normal whose horizontal part is shorter than this is vertical: its azimuth is taken as 0."""
IN_PLANE = 1e-12
"""A direction with |u · z_d| below this lies in the plane of the detector and has no slopes."""


def _number(value, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, int | float) or not math.isfinite(value):
        raise ValueError(f"{name} must be a finite number")
    return float(value)


def _like(value, *arguments):
    value = np.asarray(value, dtype=float)
    return float(value) if all(np.ndim(a) == 0 for a in arguments) else value


@dataclass(frozen=True)
class Pose:
    """Where a detector is and how it is turned.

    position_m    reference point (centre of the sensitive volume) in ENU, m
    azimuth_deg   α, azimuth of the horizontal projection of the normal, from true north,
                  clockwise; kept in [0, 360)
    tilt_deg      β, angle between the normal and the vertical, 0 to 90; 0 is a horizontal film
    roll_deg      γ, rotation about the normal, counter-clockwise seen from its tip; kept in
                  (−180, 180]

    With the three angles at zero the detector frame is ENU. At a tilt of zero the azimuth and
    the roll turn about the same axis, clockwise and counter-clockwise: (α, 0, γ) and
    (0, 0, γ − α) are one orientation, and `from_axes` gives the second.
    """

    position_m: tuple[float, float, float]
    azimuth_deg: float = 0.0
    tilt_deg: float = 0.0
    roll_deg: float = 0.0

    def __post_init__(self) -> None:
        position = np.asarray(self.position_m, dtype=float)
        if position.shape != (3,) or not np.all(np.isfinite(position)):
            raise ValueError("position_m must be three finite coordinates in ENU, m")
        object.__setattr__(self, "position_m", tuple(position.tolist()))
        azimuth = _number(self.azimuth_deg, "azimuth_deg") % 360.0
        tilt = _number(self.tilt_deg, "tilt_deg")
        roll = _number(self.roll_deg, "roll_deg")
        if not 0.0 <= tilt <= 90.0:
            raise ValueError("tilt_deg must lie between 0 (horizontal) and 90 (vertical)")
        roll = 0.0 - ((180.0 - roll) % 360.0 - 180.0)  # into (−180, 180]
        object.__setattr__(self, "azimuth_deg", azimuth)
        object.__setattr__(self, "tilt_deg", tilt)
        object.__setattr__(self, "roll_deg", roll)
        (ca, sa), (cb, sb), (cg, sg) = cos_sin_deg(azimuth), cos_sin_deg(tilt), cos_sin_deg(roll)
        z = np.array([sb * sa, sb * ca, cb])
        x0 = np.array([ca, -sa, 0.0])
        y0 = np.array([cb * sa, cb * ca, -sb])  # z × x0
        axes = np.array([cg * x0 + sg * y0, -sg * x0 + cg * y0, z])
        axes.flags.writeable = False
        object.__setattr__(self, "_axes", axes)

    @property
    def axes(self) -> np.ndarray:
        """The axes of the detector in ENU: rows x_d, y_d, z_d; shape (3, 3), read-only."""
        return self._axes

    @property
    def normal(self) -> np.ndarray:
        """The normal z_d of the detector in ENU."""
        return self._axes[2]

    @classmethod
    def from_axes(cls, position_m, axes) -> Pose:
        """The pose whose detector axes in ENU are the rows x_d, y_d, z_d of `axes`.

        The axes must be orthonormal and right-handed, and the normal must not point below the
        horizon. For a vertical normal the azimuth is 0 and the roll carries the rotation.
        """
        a = np.asarray(axes, dtype=float)
        if a.shape != (3, 3) or not np.all(np.isfinite(a)):
            raise ValueError("axes: three rows x_d, y_d, z_d of three finite components")
        if np.abs(a @ a.T - np.eye(3)).max() > ORTHONORMAL or np.abs(np.cross(a[0], a[1]) - a[2]).max() > (
            ORTHONORMAL
        ):
            raise ValueError("axes must be orthonormal and right-handed: z_d = x_d × y_d")
        x, z = a[0], a[2]
        if z[2] < -ORTHONORMAL:
            raise ValueError("the normal of a detector does not point below the horizon: tilt is at most 90")
        horizontal = math.hypot(z[0], z[1])
        tilt = math.atan2(horizontal, max(z[2], 0.0))
        azimuth = math.atan2(z[0], z[1]) if horizontal > LEVEL else 0.0
        ca, sa, cb, sb = math.cos(azimuth), math.sin(azimuth), math.cos(tilt), math.sin(tilt)
        x0 = np.array([ca, -sa, 0.0])
        y0 = np.array([cb * sa, cb * ca, -sb])
        roll = math.atan2(float(x @ y0), float(x @ x0))
        return cls(position_m, math.degrees(azimuth), math.degrees(tilt), math.degrees(roll))

    def slopes(self, direction):
        """Slopes (tx, ty) in the detector frame of directions given in ENU; shape (..., 3) in.

        A direction and its opposite have the same slopes. The directions need not be unit
        vectors; one in the plane of the detector has no slopes and is refused.
        """
        u = np.asarray(direction, dtype=float)
        if u.shape[-1:] != (3,) or not np.all(np.isfinite(u)):
            raise ValueError("direction needs three finite components along its last axis")
        w = u @ self._axes.T
        length = np.linalg.norm(u, axis=-1)
        if np.any(np.abs(w[..., 2]) <= IN_PLANE * length):
            raise ValueError("a direction in the plane of the detector has no slopes")
        return w[..., 0] / w[..., 2], w[..., 1] / w[..., 2]

    def direction(self, tx, ty) -> np.ndarray:
        """Unit vector in ENU of the slopes (tx, ty), on the side of the normal; shape (..., 3).

        u = (tx x_d + ty y_d + z_d) / sqrt(1 + tx² + ty²), so u · z_d > 0. For a tilted detector
        it can point below the horizon; `sky_direction` gives the side that looks at the sky.
        """
        tx, ty = np.broadcast_arrays(np.asarray(tx, dtype=float), np.asarray(ty, dtype=float))
        if not (np.all(np.isfinite(tx)) and np.all(np.isfinite(ty))):
            raise ValueError("slopes must be finite")
        x, y, z = self._axes
        u = tx[..., None] * x + ty[..., None] * y + z
        return u / np.sqrt(1.0 + tx * tx + ty * ty)[..., None]

    def sky_direction(self, tx, ty) -> tuple[np.ndarray, np.ndarray]:
        """The sky side of the lines with slopes (tx, ty): unit vectors in ENU with u_z ≥ 0, and `front`.

        `front` is true where the sky side is on the side of the normal (u · z_d > 0) and false
        where it is behind the plane of the detector. A line in the horizontal plane keeps the
        side of the normal and u_z = 0: it carries no muons.
        """
        u = self.direction(tx, ty)
        front = u[..., 2] >= 0.0
        return np.where(front[..., None], u, -u), front


class LinesOfSight:
    """Lines of sight of a detector: directions towards the sky and the solid angle of each.

    pose                the pose of the detector
    direction           unit vectors in ENU from the detector towards the sky (u_z > 0); shape (m, 3)
    solid_angle_sr      solid angle each direction stands for, sr; shape (m,)
    bin                 flat index of the bin of each line in `grid` (j · n_tx + i), or None
    grid, supersample   the `DirectionGrid` the lines were made from and its sub-sampling, or None
    tan_max             the lines lie within sqrt(tx² + ty²) ≤ tan_max of the detector frame, or None
    sin_elevation_min   every line is above this sine of elevation
    both_sides          whether lines whose sky side is behind the detector plane are included

    Derived, shape (m,): `tx`, `ty` (slopes in the detector frame), `cos_normal` (|u · z_d|,
    the cosine of the angle to the normal) and `front` (u · z_d > 0). The arrays are read-only.
    """

    def __init__(
        self,
        pose: Pose,
        direction,
        solid_angle_sr,
        *,
        bin=None,
        grid: DirectionGrid | None = None,
        supersample: int | None = None,
        tan_max: float | None = None,
        sin_elevation_min: float = 0.0,
        both_sides: bool = False,
    ) -> None:
        if not isinstance(pose, Pose):
            raise TypeError("pose must be a Pose")
        u = np.array(direction, dtype=float)
        weight = np.array(solid_angle_sr, dtype=float)
        if u.ndim != 2 or u.shape[1] != 3 or weight.shape != u.shape[:1]:
            raise ValueError("direction must have the shape (m, 3) and solid_angle_sr the shape (m,)")
        if not (np.all(np.isfinite(u)) and np.all(np.isfinite(weight)) and np.all(weight > 0)):
            raise ValueError("directions must be finite and solid angles finite and positive")
        if np.any(np.abs(np.linalg.norm(u, axis=1) - 1.0) > 1e-9):
            raise ValueError("direction must hold unit vectors")
        floor = _number(sin_elevation_min, "sin_elevation_min")
        if not 0.0 <= floor < 1.0:
            raise ValueError("sin_elevation_min must lie in [0, 1)")
        if np.any(u[:, 2] <= 0.0) or np.any(u[:, 2] < floor * (1.0 - 1e-12)):
            raise ValueError(
                "a line of sight points from the detector towards the sky, above sin_elevation_min"
            )
        tx, ty = pose.slopes(u)
        along = u @ pose.normal
        front = along > 0.0
        if not both_sides and not np.all(front):
            raise ValueError(
                "lines whose sky side is behind the plane of the detector need both_sides=True:"
                " they are outside the acceptance of a detector"
            )
        if tan_max is not None and np.any(np.hypot(tx, ty) > tan_max * (1.0 + 1e-9)):
            raise ValueError("a line of sight lies outside tan_max")
        if bin is not None:
            bin = np.array(bin, dtype=np.int64)
            if grid is None or bin.shape != weight.shape:
                raise ValueError("bin needs a grid and one index per line of sight")
            if bin.size and (bin.min() < 0 or bin.max() >= grid.shape[0] * grid.shape[1]):
                raise ValueError("bin: an index outside the grid")
            bin.flags.writeable = False
        self.pose = pose
        self.grid = grid
        self.supersample = supersample
        self.tan_max = None if tan_max is None else float(tan_max)
        self.sin_elevation_min = floor
        self.both_sides = bool(both_sides)
        self.bin = bin
        self.direction, self.solid_angle_sr = u, weight
        self.tx, self.ty = np.array(tx), np.array(ty)
        self.cos_normal, self.front = np.abs(along), front
        for array in (self.direction, self.solid_angle_sr, self.tx, self.ty, self.cos_normal, self.front):
            array.flags.writeable = False

    def __len__(self) -> int:
        return self.solid_angle_sr.size

    def __repr__(self) -> str:
        return (
            f"LinesOfSight({len(self)} lines, {self.solid_angle_total_sr:.4g} sr,"
            f" sin elevation above {self.sin_elevation_min:g},"
            f" {'both sides' if self.both_sides else 'front side'})"
        )

    @property
    def solid_angle_total_sr(self) -> float:
        """Solid angle of all the lines together, sr."""
        return float(self.solid_angle_sr.sum())

    def select(self, mask) -> LinesOfSight:
        """The lines of sight where the boolean `mask` (shape (m,)) is true, with what they were made from.

        For sums over a part of the field of view: the lines under soil, one side of a wall.
        """
        mask = np.asarray(mask)
        if mask.dtype != bool or mask.shape != self.solid_angle_sr.shape:
            raise ValueError("mask must hold one boolean per line of sight")
        return LinesOfSight(
            self.pose,
            self.direction[mask],
            self.solid_angle_sr[mask],
            bin=None if self.bin is None else self.bin[mask],
            grid=self.grid,
            supersample=self.supersample,
            tan_max=self.tan_max,
            sin_elevation_min=self.sin_elevation_min,
            both_sides=self.both_sides,
        )


def lines_of_sight(
    pose: Pose,
    grid: DirectionGrid,
    *,
    tan_max: float,
    supersample: int = DEFAULT_SUPERSAMPLE,
    both_sides: bool = False,
    sin_elevation_min: float = 0.0,
) -> LinesOfSight:
    """Lines of sight of a detector at `pose` through the sub-samples of `grid`.

    tan_max             lines within sqrt(tx² + ty²) ≤ tan_max; solid angles on the edge are exact
    supersample         sub-samples per side of a bin
    both_sides          keep the lines whose sky side is behind the plane of the detector
    sin_elevation_min   leave out lines at or below this sine of elevation; 0 is the horizon

    The horizon and the back side are cut by whole sub-samples, not exactly: for a horizontal
    detector give the elevation limit as `tan_max`, which is exact.
    """
    floor = _number(sin_elevation_min, "sin_elevation_min")
    if not 0.0 <= floor < 1.0:
        raise ValueError("sin_elevation_min must lie in [0, 1)")
    if not isinstance(grid, DirectionGrid):
        raise TypeError("grid must be a DirectionGrid")
    sub = grid.subsample(supersample, tan_max)
    u, front = pose.sky_direction(sub.tx, sub.ty)
    keep = (u[:, 2] > floor) & (front | bool(both_sides))
    return LinesOfSight(
        pose,
        u[keep],
        sub.solid_angle_sr[keep],
        bin=sub.bin[keep],
        grid=grid,
        supersample=supersample,
        tan_max=tan_max,
        sin_elevation_min=floor,
        both_sides=both_sides,
    )


class Detector(abc.ABC):
    """What the forward model needs from a detector.

    A detector has `tan_max` (it accepts slopes with sqrt(tx² + ty²) ≤ tan_max on the side of
    its normal), `response` (efficiency against the kinetic energy at the detector, always
    given: there is no default threshold) and `efficiency` (a constant factor on the counts),
    and says which area a track of given slopes can cross (`sensitive_area_cm2`).
    """

    tan_max: float
    response: EnergyResponse
    efficiency: float

    def _check_detector(self) -> None:
        """Checks shared by the detectors; called by their constructors."""
        tan_max = _number(self.tan_max, "tan_max")
        efficiency = _number(self.efficiency, "efficiency")
        if not tan_max > 0:
            raise ValueError("tan_max must be positive")
        if not 0 < efficiency <= 1:
            raise ValueError("efficiency must lie in (0, 1]")
        if not (isinstance(self.response, EnergyResponse) or self.response is LEGACY_RESPONSE):
            raise TypeError(
                "the energy response of a detector is always given: an EnergyResponse, for instance"
                " EnergyResponse.threshold(1.0); there is no default threshold"
            )

    @abc.abstractmethod
    def sensitive_area_cm2(self, tx, ty):
        """Area (cm²) in the plane of the detector that a track with the slopes (tx, ty) can cross."""

    def projected_area_cm2(self, tx, ty):
        """Sensitive area seen along the slopes (tx, ty), cm²: the area times cos θ = 1 / sqrt(1 + tx² + ty²).

        It does not include the acceptance: outside `tan_max` it is the area a track would
        meet, not zero (`accepts`).
        """
        tx_, ty_ = np.asarray(tx, dtype=float), np.asarray(ty, dtype=float)
        return _like(self.sensitive_area_cm2(tx_, ty_) / np.sqrt(1.0 + tx_ * tx_ + ty_ * ty_), tx, ty)

    def accepts(self, tx, ty):
        """Whether the slopes (tx, ty) lie in the acceptance, sqrt(tx² + ty²) ≤ tan_max."""
        value = np.hypot(np.asarray(tx, dtype=float), np.asarray(ty, dtype=float)) <= self.tan_max
        return bool(value) if value.ndim == 0 else value

    def lines_of_sight(
        self, pose: Pose, grid: DirectionGrid, supersample: int = DEFAULT_SUPERSAMPLE
    ) -> LinesOfSight:
        """The accepted lines of sight of the detector at `pose` through the sub-samples of `grid`.

        Within `tan_max`, on the side of the normal, above the horizon (conventions,
        "Acceptance").
        """
        return lines_of_sight(pose, grid, tan_max=self.tan_max, supersample=supersample)

    def flux_m2_s_sr(self, table, lines: LinesOfSight, opacity_g_cm2):
        """Flux over the energy response of the detector along every line of sight, m⁻² s⁻¹ sr⁻¹.

        table           a flux table of the detector altitude: anything with
                        `integral(response, sin_elevation, opacity_g_cm2)`, a `FluxTable`
        lines           the lines of sight
        opacity_g_cm2   opacity of the reference material of the table along each of them,
                        g/cm²; one value per line, or one for all (0 for the open sky)

        The table is called with the response of the detector, so the result belongs to its
        threshold. Where the table has no flux for that response — below the cut-off of a
        sea-level reference in the open sky, below its lowest elevation, beyond its largest
        opacity — the ValueError of the table is raised again with the lines of sight in it:
        nothing comes back as zero.
        """
        if not isinstance(lines, LinesOfSight):
            raise TypeError("lines must be LinesOfSight")
        opacity = np.asarray(opacity_g_cm2, dtype=float)
        if opacity.ndim not in (0, 1) or (opacity.ndim == 1 and opacity.shape != (len(lines),)):
            raise ValueError("opacity_g_cm2: one value, or one per line of sight")
        sine = lines.direction[:, 2]
        try:
            return np.asarray(table.integral(self.response, sine, opacity), dtype=float)
        except ValueError as error:
            low = math.degrees(math.asin(float(sine.min()))) if len(lines) else float("nan")
            raise ValueError(
                f"no flux for the energy response {self.response!r} on these {len(lines)} lines of sight,"
                f" the lowest at an elevation of {low:.1f} degrees: {error}."
                " Nothing is put in its place: for a sum over the rest leave those lines out"
                " (sin_elevation_min, or LinesOfSight.select for the lines under soil). The reference"
                " at the site altitude is roadmap 0.1-05b"
            ) from error
