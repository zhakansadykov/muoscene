"""detectors.grid — `DirectionGrid`: rectangular bins in the slopes (tx, ty) with exact solid angles.

The direction variables of every detector are the slopes in its own frame (conventions,
"Detector pose"), and bins are rectangles in them. In slope space the solid angle element is
dΩ = dtx dty / (1 + tx² + ty²)^(3/2), so the rectangle [0, a] × [0, b] holds

    Ω = arctan(a b / sqrt(1 + a² + b²)),

and any rectangle follows from four such corner terms, for slopes of either sign.

The acceptance of a detector is a circle, sqrt(tx² + ty²) ≤ tan_max, and the bins on its edge
lie in it only partly. Their accepted solid angle is exact as well (decided on 21 November
2026): the solid angle of a region around the axis is the contour integral of (1 − cos θ) dφ
along its boundary, which has a closed form on a straight side and on an arc of the circle.
The alternative, the share of sub-samples that fall inside the circle, puts the sum over the
default grid off by 3.4e-3, 3.5e-4, 4.3e-5 and 2.0e-4 for 1, 2, 4 and 8 sub-samples per side,
not even monotonically, and a single edge bin by up to 1/(2n) (`docs/notes/detectors.md`).

Sub-samples. A bin is divided into n × n rectangles of equal sides in the slopes; each carries
its own exact solid angle, cut by the circle when there is one, and one direction: its centre,
or, when the centre is outside the circle, a point of the part that is inside. The weights of
a bin therefore add up to its accepted solid angle, also on the edge.

A grid needs no scene and no detector: tracks are binned in it and synthetic tracks are drawn
from it by its edges alone.

Units: slopes are dimensionless (tangents), solid angles in steradians.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

DEFAULT_TAN_MAX = 1.0
"""Half-width of the default grid in the slopes: the cone of the preliminary calculation, 45 degrees."""
DEFAULT_BIN_WIDTH = 0.05
"""Bin width of the default grid in the slopes: 2.9 degrees on the axis, 40 × 40 bins (21 November 2026).

Of these bins 1176 lie in the circle of the default `tan_max`, 148 are cut by it and 276 lie outside."""
DEFAULT_SUPERSAMPLE = 4
"""Sub-samples per side of a bin by default: 16 directions per bin (21 November 2026)."""
GAP = 1e-12
"""Two points of a boundary closer than this, in the slopes, are one point."""


def _frozen(values, name: str) -> np.ndarray:
    array = np.array(values, dtype=float)
    if array.ndim != 1 or array.size < 2 or not np.all(np.isfinite(array)) or not np.all(np.diff(array) > 0):
        raise ValueError(f"{name}: at least two finite, increasing values are needed")
    array.flags.writeable = False
    return array


def _tan_max(tan_max) -> float:
    if isinstance(tan_max, bool) or not isinstance(tan_max, int | float) or not (0 < tan_max < math.inf):
        raise ValueError("tan_max must be a positive finite number")
    return float(tan_max)


def _corner(a, b):
    """Solid angle of the rectangle between the axis and the corner (a, b), signed like a · b."""
    return np.arctan(a * b / np.sqrt(1.0 + a * a + b * b))


def rectangle_solid_angle_sr(tx0, tx1, ty0, ty1):
    """Solid angle (sr) of the rectangle [tx0, tx1] × [ty0, ty1] in slope space; arrays broadcast."""
    tx0, tx1, ty0, ty1 = (np.asarray(v, dtype=float) for v in (tx0, tx1, ty0, ty1))
    return _corner(tx1, ty1) - _corner(tx0, ty1) - _corner(tx1, ty0) + _corner(tx0, ty0)


def _side(px: float, py: float, qx: float, qy: float) -> float:
    """Integral of (1 − cos θ) dφ along the straight line from P to Q in slope space.

    On a line at the distance d from the axis, with ψ the polar angle counted from the foot of
    the perpendicular, cos θ = cos ψ / sqrt(cos² ψ + d²), and the integral is
    ψ − arcsin(sin ψ / sqrt(1 + d²)). A line through the axis sweeps no angle.
    """
    cross = px * qy - py * qx
    length = math.hypot(qx - px, qy - py)
    if length == 0.0 or cross == 0.0:
        return 0.0
    d = abs(cross) / length
    ex, ey = (qx - px) / length, (qy - py) / length
    k = math.sqrt(1.0 + d * d)

    def primitive(along: float) -> float:
        return math.atan2(along, d) - math.asin(along / (math.hypot(along, d) * k))

    value = primitive(qx * ex + qy * ey) - primitive(px * ex + py * ey)
    return value if cross > 0 else -value


def _rectangle_in_circle(x0: float, x1: float, y0: float, y1: float, t: float) -> float:
    """Solid angle (sr) of the part of the rectangle [x0, x1] × [y0, y1] within the circle of radius t.

    The boundary of the part is walked counter-clockwise: the pieces of the sides that lie in
    the circle, joined by arcs of the circle where the walk leaves it and comes back.
    """
    corners = ((x0, y0), (x1, y0), (x1, y1), (x0, y1))
    pieces = []
    for i in range(4):
        (px, py), (qx, qy) = corners[i], corners[(i + 1) % 4]
        dx, dy = qx - px, qy - py
        a = dx * dx + dy * dy
        b = 2.0 * (px * dx + py * dy)
        c = px * px + py * py - t * t
        discriminant = b * b - 4.0 * a * c
        if discriminant <= 0.0:  # the line of the side misses the circle or touches it
            continue
        root = math.sqrt(discriminant)
        s0, s1 = max((-b - root) / (2.0 * a), 0.0), min((-b + root) / (2.0 * a), 1.0)
        if (s1 - s0) * math.sqrt(a) > GAP:
            pieces.append((px + s0 * dx, py + s0 * dy, px + s1 * dx, py + s1 * dy))
    cap = 1.0 - 1.0 / math.sqrt(1.0 + t * t)
    if not pieces:  # no side enters the circle: it lies in the rectangle, or apart from it
        return 2.0 * math.pi * cap if (x0 < 0.0 < x1 and y0 < 0.0 < y1) else 0.0
    total = 0.0
    for i, (px, py, qx, qy) in enumerate(pieces):
        total += _side(px, py, qx, qy)
        nx, ny = pieces[(i + 1) % len(pieces)][:2]
        if math.hypot(nx - qx, ny - qy) > GAP * max(1.0, t):
            total += cap * ((math.atan2(ny, nx) - math.atan2(qy, qx)) % (2.0 * math.pi))
    return total


def clipped_solid_angle_sr(tx0, tx1, ty0, ty1, tan_max: float) -> np.ndarray:
    """Solid angle (sr) of rectangles in slope space within the circle sqrt(tx² + ty²) ≤ tan_max.

    Exact: the rectangle formula where a rectangle lies in the circle, zero where it lies
    outside, and the contour integral for those the circle cuts. Arrays broadcast.
    """
    t = _tan_max(tan_max)
    x0, x1, y0, y1 = np.broadcast_arrays(*(np.asarray(v, dtype=float) for v in (tx0, tx1, ty0, ty1)))
    out = np.array(rectangle_solid_angle_sr(x0, x1, y0, y1), dtype=float)
    far = np.maximum(x0 * x0, x1 * x1) + np.maximum(y0 * y0, y1 * y1)
    near = np.clip(0.0, x0, x1) ** 2 + np.clip(0.0, y0, y1) ** 2
    out[near >= t * t] = 0.0
    cut = (far > t * t) & (near < t * t)
    flat = out.reshape(-1)
    for index in np.flatnonzero(cut.reshape(-1)):
        box = (float(v.reshape(-1)[index]) for v in (x0, x1, y0, y1))
        flat[index] = max(_rectangle_in_circle(*box, t), 0.0)
    return out


@dataclass(frozen=True)
class Subsamples:
    """Directions inside the bins of a grid, each with the solid angle it stands for.

    tx, ty           slopes of the directions; shape (m,)
    solid_angle_sr   solid angle of each, sr; shape (m,). Those of one bin add up to its
                     accepted solid angle
    bin              flat index of the bin of each, row-major over (ty, tx): j · n_tx + i
    """

    tx: np.ndarray
    ty: np.ndarray
    solid_angle_sr: np.ndarray
    bin: np.ndarray

    def __post_init__(self) -> None:
        for name in ("tx", "ty", "solid_angle_sr", "bin"):
            array = np.array(getattr(self, name))
            array.flags.writeable = False
            object.__setattr__(self, name, array)

    def __len__(self) -> int:
        return self.tx.size


class DirectionGrid:
    """Rectangular bins in the slopes (tx, ty) of a detector frame.

    tx_edges, ty_edges   increasing bin edges along tx and ty; they need not be equally spaced

    Arrays over the bins have the shape (number of bins in ty, number of bins in tx), as in an
    `ExpectedMap`. The arrays are read-only.
    """

    def __init__(self, tx_edges, ty_edges) -> None:
        self.tx_edges = _frozen(tx_edges, "tx_edges")
        self.ty_edges = _frozen(ty_edges, "ty_edges")
        self._solid_angle = None

    @classmethod
    def square(cls, tan_max: float = DEFAULT_TAN_MAX, bin_width: float = DEFAULT_BIN_WIDTH) -> DirectionGrid:
        """Square bins of `bin_width` in the slopes, symmetric about the axis, that cover ±`tan_max`.

        The defaults give 40 × 40 bins of 0.05 up to a slope of 1. When `tan_max` is not a
        multiple of `bin_width` the grid ends at the next multiple.
        """
        t = _tan_max(tan_max)
        if not (isinstance(bin_width, int | float) and 0 < bin_width <= t):
            raise ValueError("bin_width must be positive and at most tan_max")
        half = math.ceil(t / bin_width - 1e-9)
        edges = bin_width * np.arange(-half, half + 1)
        return cls(edges, edges)

    @classmethod
    def equal_angle(cls, tan_max: float, n_bins: int) -> DirectionGrid:
        """`n_bins` bins a side whose edges are equally spaced in angle, up to the slope `tan_max`.

        For wide fields of view: equal steps in the slope would spend most bins near the
        horizon of the detector, where the solid angle of a bin is small.
        """
        t = _tan_max(tan_max)
        if isinstance(n_bins, bool) or not isinstance(n_bins, int) or n_bins < 2:
            raise ValueError("n_bins must be an integer of at least 2")
        half = math.atan(t)
        edges = np.tan(np.linspace(-half, half, n_bins + 1))
        edges[0], edges[-1] = -t, t
        if n_bins % 2 == 0:
            edges[n_bins // 2] = 0.0
        return cls(edges, edges)

    def __repr__(self) -> str:
        n_ty, n_tx = self.shape
        return (
            f"DirectionGrid({n_tx} bins in tx from {self.tx_edges[0]:g} to {self.tx_edges[-1]:g},"
            f" {n_ty} bins in ty from {self.ty_edges[0]:g} to {self.ty_edges[-1]:g})"
        )

    def __eq__(self, other) -> bool:
        if not isinstance(other, DirectionGrid):
            return NotImplemented
        return np.array_equal(self.tx_edges, other.tx_edges) and np.array_equal(self.ty_edges, other.ty_edges)

    def __hash__(self) -> int:
        return hash((self.tx_edges.tobytes(), self.ty_edges.tobytes()))

    @property
    def shape(self) -> tuple[int, int]:
        """Number of bins along ty and along tx."""
        return self.ty_edges.size - 1, self.tx_edges.size - 1

    @property
    def tx_centers(self) -> np.ndarray:
        """Middle of every bin along tx, in the slope."""
        return 0.5 * (self.tx_edges[1:] + self.tx_edges[:-1])

    @property
    def ty_centers(self) -> np.ndarray:
        """Middle of every bin along ty, in the slope."""
        return 0.5 * (self.ty_edges[1:] + self.ty_edges[:-1])

    @property
    def solid_angle_sr(self) -> np.ndarray:
        """Solid angle of every bin, sr; shape (bins in ty, bins in tx). Exact, by four corner terms."""
        if self._solid_angle is None:
            x, y = self.tx_edges[None, :], self.ty_edges[:, None]
            value = rectangle_solid_angle_sr(x[:, :-1], x[:, 1:], y[:-1], y[1:])
            value.flags.writeable = False
            self._solid_angle = value
        return self._solid_angle

    def accepted_solid_angle_sr(self, tan_max: float) -> np.ndarray:
        """Solid angle of the part of every bin within sqrt(tx² + ty²) ≤ `tan_max`, sr.

        Exact also for the bins the circle cuts. Over a grid that covers the circle the sum is
        2π (1 − cos θ_max) with tan θ_max = `tan_max`.
        """
        x, y = self.tx_edges[None, :], self.ty_edges[:, None]
        return clipped_solid_angle_sr(x[:, :-1], x[:, 1:], y[:-1], y[1:], tan_max)

    def directions(self, pose) -> np.ndarray:
        """Unit vectors in ENU of the bin centres for a detector `pose`; shape (bins in ty, bins in tx, 3).

        The direction on the side of the normal of the detector (`Pose.direction`): for a
        tilted detector it can point below the horizon.
        """
        tx, ty = np.meshgrid(self.tx_centers, self.ty_centers)
        return pose.direction(tx, ty)

    def subsample(self, n: int = DEFAULT_SUPERSAMPLE, tan_max: float | None = None) -> Subsamples:
        """n × n directions in every bin with their solid angles.

        n         sub-samples per side of a bin
        tan_max   with it, only what lies within sqrt(tx² + ty²) ≤ tan_max: the solid angle of
                  a sub-sample is cut by the circle exactly, its direction is moved inside the
                  circle when its centre is outside, and sub-samples outside are left out

        The sub-samples of a bin are rectangles of equal sides in the slopes. Order: row-major
        over (ty, tx) of the fine grid.
        """
        if isinstance(n, bool) or not isinstance(n, int) or n < 1:
            raise ValueError("the number of sub-samples per side must be a positive integer")
        x_edges, y_edges = _refined(self.tx_edges, n), _refined(self.ty_edges, n)
        x0, x1 = x_edges[None, :-1], x_edges[None, 1:]
        y0, y1 = y_edges[:-1, None], y_edges[1:, None]
        tx, ty = np.broadcast_arrays(0.5 * (x0 + x1), 0.5 * (y0 + y1))
        tx, ty = tx.copy(), ty.copy()
        if tan_max is None:
            weight = np.array(rectangle_solid_angle_sr(x0, x1, y0, y1))
        else:
            t = _tan_max(tan_max)
            weight = clipped_solid_angle_sr(x0, x1, y0, y1, t)
            outside = (tx * tx + ty * ty > t * t) & (weight > 0)
            if np.any(outside):
                # from the point of the rectangle nearest to the axis halfway to where the line
                # towards its centre leaves the circle: inside the rectangle and inside the circle
                px = np.broadcast_to(np.clip(0.0, x0, x1), tx.shape)[outside]
                py = np.broadcast_to(np.clip(0.0, y0, y1), ty.shape)[outside]
                dx, dy = tx[outside] - px, ty[outside] - py
                a, b, c = dx * dx + dy * dy, 2.0 * (px * dx + py * dy), px * px + py * py - t * t
                s = (-b + np.sqrt(np.maximum(b * b - 4.0 * a * c, 0.0))) / (2.0 * a)
                tx[outside], ty[outside] = px + 0.5 * s * dx, py + 0.5 * s * dy
        n_ty, n_tx = self.shape
        row = np.repeat(np.arange(n_ty), n)[:, None]
        column = np.repeat(np.arange(n_tx), n)[None, :]
        index = np.broadcast_to(row * n_tx + column, tx.shape)
        keep = weight > 0
        return Subsamples(tx[keep], ty[keep], weight[keep], index[keep])


def _refined(edges: np.ndarray, n: int) -> np.ndarray:
    """Every interval of `edges` divided into n equal parts; the edges themselves are kept exactly."""
    inner = edges[:-1, None] + np.diff(edges)[:, None] * (np.arange(n) / n)
    return np.append(inner.reshape(-1), edges[-1])
