"""detectors.emulsion — `EmulsionStack`: films, cone, energy response, track density.

A stack of nuclear emulsion films is a detector without electronics: it records every charged
track that crosses it from the moment it is poured to the moment it is developed. Two numbers
follow from the flux: the tracks in the accepted cone, which make the angular map (the forward
model, roadmap 0.1-08), and the density of all tracks on a film, which the scanning tolerates
only up to a limit. This module gives the detector and the second number.

Track density. A film has no front: a muon that comes from the sky crosses it whichever side it
enters, and for a tilted film part of the sky lies behind its plane. The density is therefore
summed over both sides (decided on 21 November 2026; the back side is 0.5 % of the density of
the preliminary calculation at a tilt of 20 degrees and 9.4 % at 45 degrees), over a cone much
wider than the one accepted for the map, with |u · z_d| as the projection. It counts the muons
the energy response lets through and nothing else: no efficiency of reconstruction, no soft
component (electrons and positrons, roadmap 0.1-05b), no tracks collected before and after the
exposure (version 0.3). With a sea-level reference above sea level the open sky exists only
above the altitude cut-off: at 1150 m a threshold of 1 GeV holds from 15.7 degrees of elevation
up, and one of 0.2 GeV only under soil. A density asked for outside that is refused by the flux
table; it is not returned as zero and not filled in.

Units: cm² for areas, days for the exposure, tracks per cm² for the density.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from ..flux import EnergyResponse
from ..units import DAY
from .base import Detector, LinesOfSight, Pose, _number, lines_of_sight
from .grid import DirectionGrid

M2_PER_CM2 = 1.0e-4
TRACK_TAN_MAX = 10.0
"""Cone of the track density by default: 84.3 degrees from the normal, as in the preliminary calculation."""
TRACK_BINS = 240
"""Bins a side of the equal-angle grid of the track density: 0.7 degrees each."""
TRACK_SUPERSAMPLE = 2
"""Sub-samples per side of those bins: 2.3e5 lines of sight for a horizontal film."""


@dataclass(frozen=True)
class TrackDensity:
    """Density of muon tracks on a film and what it holds for.

    tracks_per_cm2      muons crossing a cm² of the film during the exposure
    response            the energy response it was summed over: its threshold
    exposure_days       length of the exposure, days
    solid_angle_sr      solid angle of the lines of sight it was summed over, sr
    sin_elevation_min   only lines above this sine of elevation were counted
    tan_max             only lines within this tangent of the angle to the normal, or None
    both_sides          whether the sky behind the plane of the film was counted
    n_directions        number of lines of sight in the sum
    """

    tracks_per_cm2: float
    response: EnergyResponse
    exposure_days: float
    solid_angle_sr: float
    sin_elevation_min: float
    tan_max: float | None
    both_sides: bool
    n_directions: int

    def __str__(self) -> str:
        threshold = getattr(self.response, "k_min_gev", None)
        above = repr(self.response) if threshold is None else f"K > {threshold:g} GeV"
        cone = "" if self.tan_max is None else f", within tan {self.tan_max:g} of the normal"
        return (
            f"{self.tracks_per_cm2:.4g} tracks/cm² in {self.exposure_days:g} days, {above},"
            f" elevations above {math.degrees(math.asin(self.sin_elevation_min)):.1f} degrees{cone},"
            f" {'both sides' if self.both_sides else 'front side'} of the film"
        )


@dataclass(frozen=True, kw_only=True)
class EmulsionStack(Detector):
    """A stack of emulsion films at one point.

    area_cm2     sensitive area a muon along the normal meets, cm²: that of one film (a film of
                 120 mm × 100 mm has 120 cm²) times the number of stacks laid side by side
    n_films      films of a stack, one above the other along the normal; they do not add to the
                 area. Kept for the track selection and the scanning cost of later versions
    tan_max      accepted cone: tracks with sqrt(tx² + ty²) ≤ tan_max
    response     efficiency against the kinetic energy of the muon at the detector; always
                 given, there is no default threshold. `LEGACY_RESPONSE` for a `LegacyTable`
    efficiency   constant share of the tracks in the cone that are reconstructed, in (0, 1]

    Every argument is given by name; none has a default.
    """

    area_cm2: float
    n_films: int
    tan_max: float
    response: EnergyResponse
    efficiency: float

    def __post_init__(self) -> None:
        area = _number(self.area_cm2, "area_cm2")
        if not area > 0:
            raise ValueError("area_cm2 must be positive")
        if isinstance(self.n_films, bool) or not isinstance(self.n_films, int) or self.n_films < 1:
            raise ValueError("n_films must be a positive integer")
        self._check_detector()
        object.__setattr__(self, "area_cm2", area)
        object.__setattr__(self, "tan_max", float(self.tan_max))
        object.__setattr__(self, "efficiency", float(self.efficiency))

    def sensitive_area_cm2(self, tx, ty):
        """Area a track with the slopes (tx, ty) can cross, cm²: that of the films, whatever the slopes."""
        shape = np.broadcast_shapes(np.shape(tx), np.shape(ty))
        return self.area_cm2 if shape == () else np.full(shape, self.area_cm2)

    def track_lines(
        self,
        pose: Pose,
        *,
        sin_elevation_min: float,
        tan_max: float = TRACK_TAN_MAX,
        n_bins: int = TRACK_BINS,
        supersample: int = TRACK_SUPERSAMPLE,
        both_sides: bool = True,
    ) -> LinesOfSight:
        """Lines of sight for the track density of a film at `pose`: a wide cone, both sides.

        sin_elevation_min   lines at or below this sine of elevation are left out. It is always
                            given, because the density holds only for what was summed: 0 for
                            everything above the horizon, if the flux table reaches that low
        tan_max             cone of the sum, tangent of the angle to the normal
        n_bins, supersample bins a side of an equal-angle grid and sub-samples per side of a bin
        both_sides          count the sky behind the plane of the film as well
        """
        grid = DirectionGrid.equal_angle(tan_max, n_bins)
        return lines_of_sight(
            pose,
            grid,
            tan_max=tan_max,
            supersample=supersample,
            both_sides=both_sides,
            sin_elevation_min=sin_elevation_min,
        )

    def track_density(self, table, lines: LinesOfSight, opacity_g_cm2, exposure_days: float) -> TrackDensity:
        """Muon tracks per cm² of a film over an exposure, above the energy response of the stack.

        table           flux table at the altitude of the detector (`FluxTable.integral`)
        lines           lines of sight of the sum (`track_lines`)
        opacity_g_cm2   opacity of the reference material of the table along each line, g/cm²;
                        one value per line, or one for all (0 for the open sky)
        exposure_days   length of the exposure, days

        The sum over the lines of flux × |u · z_d| × solid angle × time. The efficiency of the
        stack is not applied: it is the share of tracks reconstructed, and the film holds them
        all. Lines for which the table has no flux for this response are refused with a
        ValueError (`Detector.flux_m2_s_sr`).
        """
        days = _number(exposure_days, "exposure_days")
        if not days > 0:
            raise ValueError("exposure_days must be positive")
        flux = self.flux_m2_s_sr(table, lines, opacity_g_cm2)
        per_second = float(np.sum(flux * lines.cos_normal * lines.solid_angle_sr)) * M2_PER_CM2
        return TrackDensity(
            tracks_per_cm2=per_second * days * DAY,
            response=self.response,
            exposure_days=days,
            solid_angle_sr=lines.solid_angle_total_sr,
            sin_elevation_min=lines.sin_elevation_min,
            tan_max=lines.tan_max,
            both_sides=lines.both_sides,
            n_directions=len(lines),
        )
