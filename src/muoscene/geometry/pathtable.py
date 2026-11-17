"""geometry.pathtable — `PathTable`: path length per body and direction; opacity under a density hypothesis.

A scan of a scene gives, for every line of sight, the length of the path in every body. The
opacity along a line of sight is the sum over the bodies of density times length, so a change of
densities — another water content of a layer, a void in place of soil — is a new sum over the
same table and never a new scan (conventions, "Bodies, densities and anomalies").

A hypothesis is two optional mappings: `densities_g_cm3` gives bulk densities to grounds, layers
and structures by name, and `drho_rel` gives relative density changes to anomalies. The density
of an anomaly is that of its host under the same hypothesis times (1 + drho_rel).

Units: lengths in metres, densities in g/cm³, opacity in g/cm² (X = ρ L × 100).
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import dataclass

import numpy as np

from ..units import density_length_to_opacity

SOLID = ("ground", "layer", "structure", "anomaly")
"""Kinds of bodies that count in an opacity; the air of the scene and of the cuts does not."""
KINDS = (*SOLID, "cut", "air")


@dataclass(frozen=True)
class BodyInfo:
    """What a `PathTable` knows about a body.

    name            name given in the scene; "air" is the air of the scene box
    kind            "ground", "layer", "structure", "anomaly", "cut" or "air"
    material        name of the material in Mulder (`materials.mulder_name`); for an anomaly
                    that of its host
    density_g_cm3   bulk density in the scene; for an anomaly that of its host times (1 + drho_rel)
    host            for an anomaly, the name of the body it lies in
    drho_rel        for an anomaly, its relative density change in the scene
    """

    name: str
    kind: str
    material: str
    density_g_cm3: float
    host: str | None = None
    drho_rel: float | None = None


class PathTable:
    """Path lengths of lines of sight through the bodies of a scene.

    bodies       `BodyInfo` of every body, in the order of the last axis of `length_m`
    position_m   start of every line of sight in the world frame, m; shape (..., 3)
    direction    unit vectors from the detector towards the sky; shape (..., 3)
    length_m     path length in every body, m; shape (..., number of bodies)

    The leading shape (...) is that of the scan. The arrays are read-only.
    """

    def __init__(self, bodies, position_m, direction, length_m) -> None:
        self.bodies = tuple(bodies)
        self.position_m = _frozen(position_m)
        self.direction = _frozen(direction)
        self.length_m = _frozen(length_m)
        names = [body.name for body in self.bodies]
        if len(set(names)) != len(names) or not names:
            raise ValueError("a path table needs bodies with different names")
        if any(body.kind not in KINDS for body in self.bodies):
            raise ValueError(f"the kind of a body is one of {', '.join(KINDS)}")
        for body in self.bodies:
            if (body.kind == "anomaly") != (body.host is not None and body.drho_rel is not None):
                raise ValueError(f"{body.name}: an anomaly, and only an anomaly, has a host and a drho_rel")
            if body.kind == "anomaly" and (
                body.host not in names or self.body(body.host).kind not in SOLID[:3]
            ):
                raise ValueError(f"{body.name}: the host of an anomaly is a ground, a layer or a structure")
        shape = self.length_m.shape[:-1]
        if self.length_m.shape[-1:] != (len(names),):
            raise ValueError("length_m needs one column for every body")
        if self.position_m.shape != (*shape, 3) or self.direction.shape != (*shape, 3):
            raise ValueError("position_m and direction need the shape of length_m with three components")
        if not np.all(self.length_m >= 0):
            raise ValueError("path lengths must be finite and non-negative")

    @property
    def names(self) -> tuple[str, ...]:
        """Names of the bodies, in the order of the last axis of `length_m`."""
        return tuple(body.name for body in self.bodies)

    @property
    def shape(self) -> tuple[int, ...]:
        """Shape of the scan: that of the lines of sight."""
        return self.length_m.shape[:-1]

    def body(self, name: str) -> BodyInfo:
        """`BodyInfo` of a body; ValueError for an unknown name."""
        for body in self.bodies:
            if body.name == name:
                return body
        raise ValueError(f"no body named {name!r}; the table has {', '.join(self.names)}")

    def length_in_m(self, name: str) -> np.ndarray:
        """Path length in one body for every line of sight, m."""
        self.body(name)
        return self.length_m[..., self.names.index(name)]

    def densities_g_cm3(
        self,
        densities_g_cm3: Mapping[str, float] | None = None,
        drho_rel: Mapping[str, float] | None = None,
    ) -> np.ndarray:
        """Bulk density of every body under a hypothesis, g/cm³; one value per body.

        densities_g_cm3   bulk densities of grounds, layers and structures by name, g/cm³
        drho_rel          relative density changes of anomalies by name, at least -1
        Bodies that are not named keep the values of the scene.
        """
        densities = dict(densities_g_cm3 or {})
        changes = dict(drho_rel or {})
        for name, value in densities.items():
            if self.body(name).kind not in SOLID[:3]:
                raise ValueError(
                    f"{name}: a density is given to a ground, a layer or a structure;"
                    " an anomaly takes drho_rel, the air is that of the scene"
                )
            if not (math.isfinite(value) and value >= 0):
                raise ValueError(f"{name}: a density must be finite and non-negative")
        for name, value in changes.items():
            if self.body(name).kind != "anomaly":
                raise ValueError(f"{name}: drho_rel is given to an anomaly")
            if not (math.isfinite(value) and value >= -1):
                raise ValueError(f"{name}: drho_rel must be at least -1 (a void)")
        plain = {
            b.name: float(densities.get(b.name, b.density_g_cm3)) for b in self.bodies if b.kind != "anomaly"
        }
        out = []
        for body in self.bodies:
            if body.kind == "anomaly":
                out.append(plain[body.host] * (1.0 + float(changes.get(body.name, body.drho_rel))))
            else:
                out.append(plain[body.name])
        return np.array(out)

    def opacity_by_body_g_cm2(
        self,
        densities_g_cm3: Mapping[str, float] | None = None,
        drho_rel: Mapping[str, float] | None = None,
    ) -> np.ndarray:
        """Opacity of every body along every line of sight, g/cm²: density times length.

        Shape (..., number of bodies), the air and the cuts included, each at its own density.
        """
        return density_length_to_opacity(self.densities_g_cm3(densities_g_cm3, drho_rel), self.length_m)

    def opacity_g_cm2(
        self,
        densities_g_cm3: Mapping[str, float] | None = None,
        drho_rel: Mapping[str, float] | None = None,
        *,
        air: bool = False,
    ) -> np.ndarray:
        """Opacity along every line of sight under a hypothesis, g/cm²; shape of the scan.

        The sum of density times length over grounds, layers, structures and anomalies. The air
        of the scene and of the cuts is part of the atmosphere and is left out unless `air` is
        true.
        """
        columns = [i for i, body in enumerate(self.bodies) if air or body.kind in SOLID]
        return self.opacity_by_body_g_cm2(densities_g_cm3, drho_rel)[..., columns].sum(axis=-1)

    def opacity_by_material_g_cm2(
        self,
        densities_g_cm3: Mapping[str, float] | None = None,
        drho_rel: Mapping[str, float] | None = None,
    ) -> dict[str, np.ndarray]:
        """Opacity in every solid material along every line of sight, g/cm².

        Keys are the names of the materials in Mulder; an anomaly counts with the material of
        its host. This is what the equivalent opacity of the conventions is summed over.
        """
        opacity = self.opacity_by_body_g_cm2(densities_g_cm3, drho_rel)
        out: dict[str, np.ndarray] = {}
        for i, body in enumerate(self.bodies):
            if body.kind in SOLID:
                out[body.material] = out.get(body.material, 0.0) + opacity[..., i]
        return out


def _frozen(values) -> np.ndarray:
    array = np.array(values, dtype=float)
    if not np.all(np.isfinite(array)):
        raise ValueError("a path table holds finite numbers")
    array.flags.writeable = False
    return array
