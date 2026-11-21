"""site.frames — `Site`: origin of the world frame, local ENU and the azimuth of the axes of a scene.

The world frame of the conventions is local ENU: x east, y north, z up, metres, with its origin
at a latitude, a longitude and an altitude. A scene of `muoscene.geometry` has no rotations: a
canyon wall that does not run north-south is described along its own axes. `Site` holds the one
angle between the two frames, the azimuth of the +y axis of the scene from true north,
clockwise, and is the only place where it is applied: `Site.scan` takes points and directions
in ENU, turns them into the axes of the scene, scans, and gives the path table back in ENU. No
caller rotates a vector itself (decided on 21 November 2026, `docs/notes/detectors.md`).

The local frame of Mulder has the same angle under the name `azimuth` (`Site.mulder_frame`):
it is how Mulder knows where north is for the geomagnetic field. A scene made by `Site.scene`
carries it.

Survey coordinate systems (pyproj), grid north and magnetic declination come with version 0.2.

Units: degrees for the angles of a site, metres for points; directions are unit vectors.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import mulder
import numpy as np

from ..geometry.pathtable import PathTable
from ..geometry.scene import ALTITUDE_MAX_M, Scene


def cos_sin_deg(angle_deg: float) -> tuple[float, float]:
    """Cosine and sine of an angle in degrees, exact at the multiples of 90 degrees."""
    quarter, rest = divmod(float(angle_deg) % 360.0, 90.0)
    if rest == 0.0:
        return ((1.0, 0.0), (0.0, 1.0), (-1.0, 0.0), (0.0, -1.0))[int(quarter)]
    angle = math.radians(angle_deg)
    return math.cos(angle), math.sin(angle)


@dataclass(frozen=True)
class Site:
    """Origin of the world frame and the orientation of the axes of a scene in it.

    latitude_deg, longitude_deg   of the origin, degrees (WGS84)
    altitude_m                    of the origin above sea level, m; z = 0 of the world frame
    scene_azimuth_deg             azimuth of the +y axis of the scene from true north, clockwise,
                                  degrees, kept in [0, 360). At 0 the axes of the scene are those
                                  of ENU; at 90 its +y points east and its +x south

    The two frames share their origin and their vertical. In ENU the axes of the scene are

        x_s = (cos ψ, −sin ψ, 0),   y_s = (sin ψ, cos ψ, 0),   z_s = (0, 0, 1)

    with ψ the azimuth: the same formulas as for the axes of a detector with no tilt
    (conventions, "Detector pose").
    """

    latitude_deg: float
    longitude_deg: float
    altitude_m: float
    scene_azimuth_deg: float = 0.0

    def __post_init__(self) -> None:
        for name in ("latitude_deg", "longitude_deg", "altitude_m", "scene_azimuth_deg"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int | float) or not math.isfinite(value):
                raise ValueError(f"{name} must be a finite number")
            object.__setattr__(self, name, float(value))
        if abs(self.latitude_deg) > 90 or abs(self.longitude_deg) > 360:
            raise ValueError("latitude_deg or longitude_deg is out of range")
        if abs(self.altitude_m) > ALTITUDE_MAX_M:
            raise ValueError(f"altitude_m must lie within {ALTITUDE_MAX_M:g} m of sea level")
        object.__setattr__(self, "scene_azimuth_deg", self.scene_azimuth_deg % 360.0)

    # ------------------------------------------------------------------ frames

    @property
    def scene_axes(self) -> np.ndarray:
        """The axes of the scene in ENU: rows x_s, y_s, z_s; shape (3, 3)."""
        c, s = cos_sin_deg(self.scene_azimuth_deg)
        return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])

    def scene_from_enu(self, vectors) -> np.ndarray:
        """Points (m) or directions given in ENU, in the axes of the scene; shape (..., 3)."""
        return _vectors(vectors) @ self.scene_axes.T

    def enu_from_scene(self, vectors) -> np.ndarray:
        """Points (m) or directions given in the axes of the scene, in ENU; shape (..., 3)."""
        return _vectors(vectors) @ self.scene_axes

    def mulder_frame(self):
        """The frame of the scene for Mulder: a `mulder.LocalFrame` with its y axis at the azimuth."""
        where = {"latitude": self.latitude_deg, "longitude": self.longitude_deg, "altitude": self.altitude_m}
        if self.scene_azimuth_deg:
            where["azimuth"] = self.scene_azimuth_deg
        return mulder.LocalFrame(**where)

    # ------------------------------------------------------------------ scenes

    def scene(self, bodies, **options) -> Scene:
        """A `Scene` of `bodies`, described along the axes of the scene, at this site.

        The altitude, latitude, longitude and azimuth of the scene are those of the site;
        `options` are the other arguments of `Scene` (`atmosphere`, `margin_m`, `air`).
        """
        for name in ("altitude_m", "latitude_deg", "longitude_deg", "azimuth_deg"):
            if name in options:
                raise ValueError(f"{name} of a scene made by a site is that of the site")
        return Scene(
            bodies,
            altitude_m=self.altitude_m,
            latitude_deg=self.latitude_deg,
            longitude_deg=self.longitude_deg,
            azimuth_deg=self.scene_azimuth_deg,
            **options,
        )

    def scan(self, scene: Scene, position_m, direction) -> PathTable:
        """Path length in every body of a scene along lines of sight given in ENU.

        scene        a scene of this site (`Site.scene`): same altitude, origin and azimuth
        position_m   points in ENU, m; shape (..., 3) or (3,)
        direction    unit vectors in ENU from the detector towards the sky; shape (..., 3)

        The points and directions are turned into the axes of the scene and scanned by
        `Scene.scan`, which refuses what it refuses. The `PathTable` returned holds the points
        and directions in ENU, as given.
        """
        self._check(scene)
        position, direction = _vectors(position_m, "position_m"), _vectors(direction, "direction")
        table = scene.scan(self.scene_from_enu(position), self.scene_from_enu(direction))
        shape = (*table.shape, 3)
        return PathTable(
            table.bodies, np.broadcast_to(position, shape), np.broadcast_to(direction, shape), table.length_m
        )

    def _check(self, scene: Scene) -> None:
        if not isinstance(scene, Scene):
            raise TypeError("a site scans a muoscene.geometry.Scene")
        ours = (self.altitude_m, self.latitude_deg, self.longitude_deg, self.scene_azimuth_deg)
        theirs = (scene.altitude_m, scene.latitude_deg, scene.longitude_deg, scene.azimuth_deg)
        if ours != theirs:
            raise ValueError(
                "the scene was built for another site: altitude, latitude, longitude and azimuth are"
                f" {theirs} in the scene and {ours} in the site. Build it with Site.scene"
            )


def _vectors(values, name: str = "vectors") -> np.ndarray:
    array = np.asarray(values, dtype=float)
    if array.shape[-1:] != (3,):
        raise ValueError(f"{name} needs three components along its last axis")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must be finite")
    return array
