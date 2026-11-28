"""data.binning — tracks into counts per bin of a `DirectionGrid`, with the cuts of the model.

`bin_tracks` puts the tracks of a `TrackSet` into the bins of the grid of an expected map and
counts every track it leaves out, by the reason (roadmap 0.1-15, decided on 28 November 2026,
`docs/notes/tracks.md`):

- Edges. A bin holds its lower edges and not its upper ones; the last bin along each axis holds
  both: a track on an edge between two bins belongs to the bin above it, as in
  `numpy.histogram2d`. Edges are compared as the grid holds them, without a tolerance.
- Circle. A track with sqrt(tx² + ty²) ≤ `tan_max` is inside the acceptance, the circle included,
  as `Detector.accepts`; one beyond it is counted in `outside_circle`.
- Grid. A track inside the circle but beyond the edges of the grid is counted in `outside_grid`.
- Side. With the pose of the detector, a track counts only when its line looks at the sky in
  front of the plane of the detector and above the horizon: the direction of its slopes on the
  side of the normal (`Pose.direction`) has u_z > 0. A track whose line reaches the sky only
  from behind the plane, or along the horizon, is counted in `behind`. A film records such a
  muon, but the acceptance of a detector (`docs/conventions.md`, "Acceptance") and so the model
  leave it out. Without a pose the side is not looked at: tracks in the frame of a scan have none.
- Bins of the model. The bins the model accepts are those with an accepted line of sight: with a
  pose, the lines of sight that `Detector.lines_of_sight` gives the forward model
  (`muoscene.detectors.lines_of_sight` with the same `tan_max` and `supersample`), so that
  `accepted` is false exactly where the status of the map is 1, "not accepted"; without a pose,
  the bins with a part inside the circle. A track that passes the cuts above and still falls in a
  bin the model does not accept — on the circle at the corner of a bin that touches it, or next
  to the horizon where the model cuts by whole sub-samples — is counted in `not_accepted`.

Every track is in exactly one of `counts`, `outside_circle`, `outside_grid`, `behind` and
`not_accepted`: nothing is lost without a number. The model cuts the horizon and the back of the
plane by whole sub-samples and the tracks exactly, so in the bins those lines cross the two
differ by a part of a sub-sample (`docs/notes/tracks.md`).

Layer 4: may import only `muoscene` modules of layer <= 4 (see `muoscene._layers`).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from ..detectors import DEFAULT_SUPERSAMPLE, DirectionGrid, Pose, TrackSet, lines_of_sight

NOT_ACCEPTED = 1
"""Status of a bin the model does not accept, as in an `ExpectedMap`; 0 for one it does."""


def _frozen(values, dtype) -> np.ndarray:
    array = np.array(values, dtype=dtype)
    array.flags.writeable = False
    return array


@dataclass(frozen=True, eq=False)
class BinnedTracks:
    """Counts of tracks per bin of a grid and the number of tracks left out, by the reason.

    counts          tracks per bin, shape (bins in ty, bins in tx), integers; zero in the bins
                    the model does not accept
    accepted        whether the model accepts the bin: it has an accepted line of sight
    grid            the `DirectionGrid`
    tan_max         radius of the acceptance circle in the slopes
    pose            the pose of the detector, or None when the side was not looked at
    supersample     sub-samples per side of a bin of the lines of sight, with a pose; None without
    frame           frame of the tracks binned, "detector" or "scan"
    detector        name of the detector
    n_tracks        tracks in the set binned
    outside_circle, outside_grid, behind, not_accepted
                    tracks left out (see the module docstring)

    Arrays are read-only. `status` gives the codes for `muoscene.viz.maps.array_figure`.
    """

    counts: np.ndarray
    accepted: np.ndarray
    grid: DirectionGrid
    tan_max: float
    pose: Pose | None
    supersample: int | None
    frame: str
    detector: str
    n_tracks: int
    outside_circle: int
    outside_grid: int
    behind: int
    not_accepted: int

    def __post_init__(self) -> None:
        counts = _frozen(self.counts, np.int64)
        accepted = _frozen(self.accepted, bool)
        if counts.shape != self.grid.shape or accepted.shape != self.grid.shape:
            raise ValueError(f"counts and accepted must have the shape of the grid, {self.grid.shape}")
        if np.any(counts < 0) or np.any(counts[~accepted] != 0):
            raise ValueError("counts are not negative and are zero in the bins the model does not accept")
        left_out = (self.outside_circle, self.outside_grid, self.behind, self.not_accepted)
        if any(n < 0 for n in left_out) or int(counts.sum()) + sum(left_out) != self.n_tracks:
            raise ValueError("every track is binned or counted as left out, once")
        object.__setattr__(self, "counts", counts)
        object.__setattr__(self, "accepted", accepted)

    @property
    def binned(self) -> int:
        """Number of tracks in the bins."""
        return int(self.counts.sum())

    @property
    def status(self) -> np.ndarray:
        """0 where the model accepts the bin, 1 ("not accepted") where not; int8, the shape of the grid."""
        return np.where(self.accepted, 0, NOT_ACCEPTED).astype(np.int8)

    def summary(self) -> dict:
        """The numbers of the binning as plain values: tracks, binned, left out by reason, bins."""
        return {
            "detector": self.detector,
            "frame": self.frame,
            "tracks": self.n_tracks,
            "binned": self.binned,
            "outside_circle": self.outside_circle,
            "outside_grid": self.outside_grid,
            "behind": self.behind,
            "not_accepted": self.not_accepted,
            "bins": int(self.accepted.size),
            "bins_accepted": int(self.accepted.sum()),
            "tan_max": self.tan_max,
        }


def _index(edges: np.ndarray, values: np.ndarray) -> np.ndarray:
    """Bin of each value: lower edge included, upper excluded, the last edge in the last bin; −1 outside."""
    index = np.searchsorted(edges, values, side="right") - 1
    index[values == edges[-1]] = edges.size - 2
    index[(values < edges[0]) | (values > edges[-1])] = -1
    return index


def bin_tracks(
    tracks: TrackSet,
    grid: DirectionGrid,
    *,
    tan_max: float,
    pose: Pose | None = None,
    supersample: int = DEFAULT_SUPERSAMPLE,
) -> BinnedTracks:
    """Count the tracks of `tracks` in the bins of `grid` with the cuts of the model.

    tracks        a `TrackSet`; with a pose it must be in the frame of the detector
    grid          the grid of the expected map (`muoscene.model.grid_of`)
    tan_max       radius of the acceptance circle of the detector in the slopes (`Detector.tan_max`)
    pose          the pose of the detector, to leave out tracks whose line looks at the sky from
                  behind its plane or along the horizon and to find the bins the model accepts;
                  None looks only at the circle
    supersample   sub-samples per side of a bin with which the map was computed (its attribute
                  "supersample"); used only with a pose

    Returns a `BinnedTracks`: the counts, the bins the model accepts and the tracks left out.
    """
    if not isinstance(tracks, TrackSet):
        raise TypeError("tracks must be a TrackSet")
    if not isinstance(grid, DirectionGrid):
        raise TypeError("grid must be a DirectionGrid")
    if isinstance(tan_max, bool) or not isinstance(tan_max, int | float) or not 0 < tan_max < math.inf:
        raise ValueError("tan_max must be a positive finite number")
    tan_max = float(tan_max)
    if pose is not None:
        if not isinstance(pose, Pose):
            raise TypeError("pose must be a Pose or None")
        if tracks.frame != "detector":
            raise ValueError(
                "tracks in the frame of a scan have no side: bin them without a pose, or read them"
                " again with the scan-frame transform of the exposure (muoscene.io.tracks.ScanFrame)"
            )
        lines = lines_of_sight(pose, grid, tan_max=tan_max, supersample=supersample)
        n_ty, n_tx = grid.shape
        accepted = np.bincount(lines.bin, minlength=n_ty * n_tx).reshape(grid.shape) > 0
    else:
        supersample = None
        accepted = grid.accepted_solid_angle_sr(tan_max) > 0

    tx, ty = tracks.tx, tracks.ty
    inside = np.hypot(tx, ty) <= tan_max
    i, j = _index(grid.tx_edges, tx), _index(grid.ty_edges, ty)
    on_grid = inside & (i >= 0) & (j >= 0)
    if pose is not None:
        front = pose.direction(tx, ty)[..., 2] > 0.0 if len(tracks) else np.zeros(0, bool)
    else:
        front = np.ones(len(tracks), bool)
    candidate = on_grid & front
    ok = candidate.copy()
    ok[candidate] = accepted[j[candidate], i[candidate]]
    n_ty, n_tx = grid.shape
    counts = np.bincount(j[ok] * n_tx + i[ok], minlength=n_ty * n_tx).reshape(grid.shape)
    return BinnedTracks(
        counts=counts,
        accepted=accepted,
        grid=grid,
        tan_max=tan_max,
        pose=pose,
        supersample=supersample,
        frame=tracks.frame,
        detector=tracks.detector,
        n_tracks=len(tracks),
        outside_circle=int(np.count_nonzero(~inside)),
        outside_grid=int(np.count_nonzero(inside & ~on_grid)),
        behind=int(np.count_nonzero(on_grid & ~front)),
        not_accepted=int(np.count_nonzero(candidate & ~ok)),
    )
