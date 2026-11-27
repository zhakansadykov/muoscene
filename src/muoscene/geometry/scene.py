"""geometry.scene — `Scene`: built geometry, media, densities, `scan` through Mulder.

A `Scene` takes the bodies of `primitives`, has `builder` turn them into a Calzone geometry,
wraps it as a local geometry of Mulder at the altitude of the site and gives every medium its
material and bulk density. `Scene.scan` then returns the path length in every body for every
line of sight, a `PathTable`.

Air. Mulder gives the air of a local geometry one density, by default that of sea level
(1.205 kg/m³), and continues with its atmosphere model only above the box. The air of a scene —
its box and its cuts — therefore takes the density of the atmosphere model at the altitude of the
site, and the box ends `margin_m` above the highest body: 1000 m of sea-level air above the brow
would lower the counts of the reference void at 1150 m by 0.9 %, and 1000 m of air of the site
density still by 0.3 % (`docs/numbers.md`).

What a scan refuses, because Mulder 0.3.7 answers it with a crash, an endless loop or a wrong
number: a point outside the box of the scene, a direction that is not finite or not a unit
vector, a point on a boundary between bodies, a line of sight cut by the limits of a ground, and
a line that meets an edge of a mesh at a grazing angle. A line that meets an edge or a vertex of
a mesh otherwise is moved off it by a few nanometres before Mulder sees it (`edges`, roadmap
0.1-10b); `traced_lines` gives the lines as Mulder is given them.

Units at this boundary: metres on both sides; g/cm³ become kg/m³ of Mulder.
"""

from __future__ import annotations

import math

import mulder
import numpy as np

from ..materials.backend import register
from ..materials.composition import Material, Mixture
from ..materials.library import AIR
from . import checks, edges
from .builder import assemble, prism_mesh
from .pathtable import BodyInfo, PathTable
from .primitives import Anomaly, Cut, Ground, Prism, Structure

KG_M3_PER_G_CM3 = 1.0e3
UNIT_TOLERANCE = 1e-9
"""A direction whose length differs from one by more than this is refused."""
VOID_KG_M3 = 1.0e-9
"""Density given to Mulder for a medium of zero density. Mulder reads a zero as "not set" and
takes the density of the material: the flux through a void would be that through dense soil."""
CHUNK = 65536
"""Lines of sight scanned at a time: bounds the memory of the lengths per medium."""
ALTITUDE_MAX_M = 9000.0


def air_density_g_cm3(altitude_m: float, atmosphere: str = "us-standard") -> float:
    """Density of the air (g/cm³) at an altitude above sea level (m) in an atmosphere model of Mulder."""
    if atmosphere not in mulder.Atmosphere.models:
        raise ValueError(
            f"unknown atmosphere {atmosphere!r}; Mulder has {', '.join(mulder.Atmosphere.models)}"
        )
    if not (math.isfinite(altitude_m) and abs(altitude_m) <= ALTITUDE_MAX_M):
        raise ValueError(f"altitude_m must lie within {ALTITUDE_MAX_M:g} m of sea level")
    return float(mulder.Atmosphere(atmosphere).density(float(altitude_m))) / KG_M3_PER_G_CM3


class Scene:
    """Bodies in the world frame, built for scanning.

    bodies          `Ground`, `Structure`, `Cut` and `Anomaly` objects; the order is that of the
                    columns of a `PathTable` (the layers of a ground follow it, the air is last)
    altitude_m      altitude of the origin of the world frame above sea level, m
    atmosphere      density profile of the air, a name of Mulder; give the flux tables the same
    latitude_deg, longitude_deg
                    of the origin, degrees; both or none (then those of Mulder, 0 and 0). They
                    do not enter a scan; Mulder uses them for the geomagnetic field
    azimuth_deg     azimuth of the +y axis of the scene from true north, clockwise, degrees;
                    0 when the axes of the scene are east, north, up. It does not enter a scan
                    either, and nothing is rotated here: it only tells Mulder where north is,
                    for the geomagnetic field. Points and directions in ENU go through
                    `muoscene.site.Site.scan`, which turns them into the axes of the scene
    margin_m        the box of the scene is this much larger than its bodies on every side, m;
                    a detector must lie inside the box
    air             material of the air of the scene and of its cuts

    Attributes: `bodies`, `altitude_m`, `latitude_deg`, `longitude_deg`, `azimuth_deg`,
    `atmosphere`, `air_density_g_cm3`, `lo_m` and `hi_m`
    (corners of the box, m), `hosts` (the host of every anomaly), `calzone` (the
    `calzone.Geometry`; its materials are placeholders) and `mulder` (the `mulder.LocalGeometry`
    with the materials and densities of the scene, for direct transport).
    """

    def __init__(
        self,
        bodies,
        *,
        altitude_m: float,
        atmosphere: str = "us-standard",
        latitude_deg: float | None = None,
        longitude_deg: float | None = None,
        azimuth_deg: float = 0.0,
        margin_m: float = 1.0,
        air: Material | Mixture = AIR,
    ) -> None:
        self.bodies = tuple(bodies)
        self.altitude_m = float(altitude_m)
        self.atmosphere = atmosphere
        self.air_density_g_cm3 = air_density_g_cm3(self.altitude_m, atmosphere)
        if (latitude_deg is None) != (longitude_deg is None):
            raise ValueError("give both latitude_deg and longitude_deg, or neither")
        if latitude_deg is not None and not (abs(latitude_deg) <= 90 and abs(longitude_deg) <= 360):
            raise ValueError("latitude_deg or longitude_deg is out of range")
        if not math.isfinite(azimuth_deg):
            raise ValueError("azimuth_deg must be finite")
        self.latitude_deg = None if latitude_deg is None else float(latitude_deg)
        self.longitude_deg = None if longitude_deg is None else float(longitude_deg)
        self.azimuth_deg = float(azimuth_deg) % 360.0
        if not isinstance(air, Material | Mixture):
            raise TypeError("the air of a scene is a Material or a Mixture")
        assembly = assemble(self.bodies, margin_m)
        self.lo_m, self.hi_m = assembly.lo_m, assembly.hi_m
        self.hosts = dict(assembly.hosts)
        self.calzone = assembly.geometry
        self._volumes = assembly.volumes
        self._grounds = tuple(body for body in self.bodies if isinstance(body, Ground))
        self._edges = edges.mesh_edges(
            prism_mesh(volume.shape) for volume in self._volumes if isinstance(volume.shape, Prism)
        )

        described = _describe(self.bodies, self.hosts, register(air), self.air_density_g_cm3)
        self.table_bodies = tuple(described)
        """`BodyInfo` of every column of a `PathTable` of this scene."""
        where = (
            {}
            if latitude_deg is None
            else {"latitude": float(latitude_deg), "longitude": float(longitude_deg)}
        )
        if self.azimuth_deg:
            # the local frame of Mulder has its y axis at this azimuth from geographic north
            where["azimuth"] = self.azimuth_deg
        frame = mulder.LocalFrame(altitude=self.altitude_m, **where)
        self.mulder = mulder.LocalGeometry(self.calzone, frame=frame)
        by_volume = {volume.name: volume.body for volume in self._volumes}
        info = {body.name: body for body in described}
        names = [body.name for body in described]
        columns = []
        for medium in self.mulder.media:
            body = info[by_volume[medium.description.split(".")[-1]]]
            medium.material = body.material
            medium.density = max(KG_M3_PER_G_CM3 * body.density_g_cm3, VOID_KG_M3)
            columns.append(names.index(body.name))
        self._medium_body = np.array(columns)
        # media sorted by body, so that the lengths of one body are summed in one pass
        self._order = np.argsort(self._medium_body, kind="stable")
        sorted_bodies = self._medium_body[self._order]
        self._starts = np.flatnonzero(np.diff(sorted_bodies, prepend=-1))
        self._present = sorted_bodies[self._starts]

    # ------------------------------------------------------------------ scanning

    def scan(self, position_m, direction) -> PathTable:
        """Path length in every body along lines of sight.

        position_m   points in the world frame, m; shape (..., 3) or (3,)
        direction    unit vectors from the detector towards the sky (u_z > 0); shape (..., 3)
        The two are broadcast against each other. Returns a `PathTable` of that shape.

        ValueError for a point outside the box of the scene or closer than 1 µm to a boundary
        between bodies, for a direction that is not a unit vector pointing upwards, for a line of
        sight that crosses a limit of the model inside a ground (`Ground`), and for one that meets
        an edge of a mesh at less than 0.57 degrees to a face there (`edges`). The table holds
        the points as given; Mulder scans them as `traced_lines` gives them.
        """
        position, direction, o, u, _ = self._traced(position_m, direction)
        shape = position.shape[:-1]
        lengths = np.zeros((len(o), len(self.table_bodies)))
        for start in range(0, len(o), CHUNK):
            part = slice(start, start + CHUNK)
            per_medium = np.asarray(self.mulder.scan(position=o[part], direction=u[part]), dtype=float)
            per_medium = per_medium.reshape(-1, len(self._medium_body))
            lengths[part, self._present] = np.add.reduceat(per_medium[:, self._order], self._starts, axis=1)
        return PathTable(self.table_bodies, position, direction, lengths.reshape(*shape, lengths.shape[1]))

    def traced_lines(self, position_m, direction) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """The lines of sight as Mulder is given them, after the checks of `scan`.

        position_m, direction   as for `scan`
        Returns the points (m) and the unit directions, broadcast to one shape (..., 3), and a
        boolean array (...) of the lines whose point was moved off an edge of a mesh, by 3 or
        9 nm across the line (`edges.move_off_edges`); the other points are those given, bit for
        bit. `scan` and `muoscene.site.Site.differential_flux` both trace these lines. Raises
        what `scan` raises.
        """
        position, direction, o, u, moved = self._traced(position_m, direction)
        return o.reshape(position.shape), u.reshape(position.shape), moved.reshape(position.shape[:-1])

    def _traced(self, position_m, direction):
        """The checked lines broadcast as given, then flat as Mulder takes them: points moved off edges."""
        position, direction = self._lines(position_m, direction)
        o, u = position.reshape(-1, 3), direction.reshape(-1, 3)
        self._refuse_cut_lines(o, u)
        moved_o, moved = edges.move_off_edges(self._edges, o, u)
        return position, direction, moved_o, u, moved

    def locate(self, position_m) -> np.ndarray:
        """Name of the body at every point (m, shape (..., 3)); an array of strings of shape (...)."""
        position = self._points(position_m)
        index = np.asarray(self.mulder.locate(position=position.reshape(-1, 3))).reshape(-1)
        names = np.array([body.name for body in self.table_bodies])
        return names[self._medium_body[index]].reshape(position.shape[:-1])

    def distance_to_boundary_m(self, position_m) -> tuple[float, str]:
        """Distance (m) of one point to the nearest boundary between bodies, and the body it bounds.

        Inside a cut only the surface of the cut counts: the ground it is cut from is not there.
        """
        distance, index = self._boundaries(np.asarray(position_m, dtype=float).reshape(1, 3))
        return float(distance[0]), self._volumes[index[0]].body

    def _boundaries(self, points: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Distance of every point (shape (n, 3)) to the nearest boundary and the index of its volume."""
        signed = np.array([volume.shape.signed_distance_m(points) for volume in self._volumes[1:]])
        nearest = np.abs(signed).argmin(axis=0)
        distance = np.abs(signed).min(axis=0)
        for i, volume in enumerate(self._volumes[1:]):
            if volume.kind == "cut":
                inside = signed[i] < 0
                nearest, distance = np.where(inside, i, nearest), np.where(inside, -signed[i], distance)
        return distance, nearest + 1

    def check(self, resolution: int | None = None) -> None:
        """Run the overlap check of Calzone on the built geometry; ValueError if volumes overlap.

        The description is checked analytically when the scene is built; this is the independent
        check of Geant4, which samples `resolution` points on the surface of every volume.
        """
        checks.check_overlaps(self.calzone, resolution)

    # ------------------------------------------------------------------ helpers

    def _points(self, position_m) -> np.ndarray:
        position = np.asarray(position_m, dtype=float)
        if position.shape[-1:] != (3,):
            raise ValueError("position_m needs three coordinates along its last axis")
        if not np.all(np.isfinite(position)):
            raise ValueError("position_m must be finite")
        if not (np.all(position > self.lo_m) and np.all(position < self.hi_m)):
            raise ValueError(
                f"a point lies outside the box of the scene, {self.lo_m} to {self.hi_m} m:"
                " move it or build the scene with a larger margin_m"
            )
        return position

    def _lines(self, position_m, direction) -> tuple[np.ndarray, np.ndarray]:
        position = self._points(position_m)
        direction = np.asarray(direction, dtype=float)
        if direction.shape[-1:] != (3,):
            raise ValueError("direction needs three components along its last axis")
        if not np.all(np.isfinite(direction)):
            raise ValueError("direction must be finite")
        norm = np.linalg.norm(direction, axis=-1, keepdims=True)
        if not np.all(np.abs(norm - 1.0) <= UNIT_TOLERANCE):
            raise ValueError("direction must hold unit vectors")
        direction = direction / norm
        if not np.all(direction[..., 2] > 0):
            raise ValueError("a direction points from the detector towards the sky: u_z must be positive")
        shape = np.broadcast_shapes(position.shape, direction.shape)
        points = np.unique(position.reshape(-1, 3), axis=0)
        distance, index = self._boundaries(points)
        if distance.size and distance.min() < checks.BOUNDARY_M:
            worst = int(distance.argmin())
            raise ValueError(
                f"the point {tuple(points[worst].tolist())} m lies on a boundary of"
                f" {self._volumes[index[worst]].body!r} (within {checks.BOUNDARY_M:g} m):"
                " the medium there is ambiguous, move it off"
            )
        return (
            np.ascontiguousarray(np.broadcast_to(position, shape)),
            np.ascontiguousarray(np.broadcast_to(direction, shape)),
        )

    def _refuse_cut_lines(self, o: np.ndarray, u: np.ndarray) -> None:
        for ground in self._grounds:
            cut = checks.cut_by_limits(ground, o, u)
            if np.any(cut):
                lowest = math.degrees(math.asin(float(u[cut, 2].min())))
                highest = math.degrees(math.asin(float(u[cut, 2].max())))
                raise ValueError(
                    f"{int(cut.sum())} of {len(o)} lines of sight cross a limit of the model inside ground"
                    f" {ground.name!r} (its bottom, its faces across y or the ends of its profile), at"
                    f" elevations of {lowest:.2f} to {highest:.2f} degrees: the path beyond is not modelled."
                    " Extend the ground or leave these directions out"
                )


def _describe(bodies, hosts: dict[str, str], air: str, air_density_g_cm3: float) -> list[BodyInfo]:
    """`BodyInfo` of every column: bodies in the given order, layers after their ground, the air last."""
    out: list[BodyInfo] = []
    for body in bodies:
        if isinstance(body, Ground):
            out.append(BodyInfo(body.name, "ground", register(body.composition), body.bulk_density_g_cm3))
            for layer in body.layers:
                out.append(
                    BodyInfo(layer.name, "layer", register(layer.composition), layer.bulk_density_g_cm3)
                )
        elif isinstance(body, Structure):
            out.append(BodyInfo(body.name, "structure", register(body.composition), body.bulk_density_g_cm3))
        elif isinstance(body, Cut):
            out.append(BodyInfo(body.name, "cut", air, air_density_g_cm3))
    solid = {body.name: body for body in out}
    for body in bodies:
        if isinstance(body, Anomaly):
            host = solid[hosts[body.name]]
            density = host.density_g_cm3 * (1.0 + body.drho_rel)
            out.append(BodyInfo(body.name, "anomaly", host.material, density, host.name, body.drho_rel))
    out.append(BodyInfo(checks.AIR, "air", air, air_density_g_cm3))
    return out
