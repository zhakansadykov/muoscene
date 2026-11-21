"""model.hypotheses — `Hypothesis` and `HypothesisSet`: what a map is computed for, without a new scan.

A hypothesis is a statement about densities: other bulk densities of grounds, layers and
structures, and other relative density changes of anomalies — no anomaly (0), a void (−1), a
decompaction zone (−0.2). Shapes and materials are those of the scene, so the opacity under any
hypothesis is a new sum over the one `PathTable` of a scan
(`PathTable.opacity_by_material_g_cm2`) and never a new scan. A hypothesis that changes a shape
or a composition — another water content changes the material — is another scene.

A body that a hypothesis does not name keeps what the scene gave it: `Hypothesis.as_built` names
none and is the scene as it was built, the one Mulder holds and the fast path is checked on.

Variants of systematics, the `variant` dimension of the conventions, come with version 0.3.

Units: densities in g/cm³; `drho_rel` is dimensionless, at least −1.

Layer 3.
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Iterator, Mapping
from dataclasses import dataclass, field
from types import MappingProxyType

from ..geometry import PathTable


def _frozen(values, what: str, lowest: float) -> Mapping[str, float]:
    if not isinstance(values, Mapping):
        raise TypeError(f"{what}: a mapping from the name of a body to a number")
    out = {}
    for name, value in values.items():
        if not isinstance(name, str) or not name:
            raise ValueError(f"{what}: bodies are given by name")
        if isinstance(value, bool) or not isinstance(value, int | float) or not math.isfinite(value):
            raise ValueError(f"{what}[{name!r}] must be a finite number")
        if value < lowest:
            raise ValueError(f"{what}[{name!r}] must be at least {lowest:g}")
        out[name] = float(value)
    return MappingProxyType(dict(sorted(out.items())))


@dataclass(frozen=True, eq=False)
class Hypothesis:
    """Densities a map is computed for.

    name              name of the hypothesis: its label along the `hypothesis` dimension of a map
    densities_g_cm3   bulk densities of grounds, layers and structures by name, g/cm³
    drho_rel          relative density changes of anomalies by name: 0 no anomaly, −1 a void

    Bodies that are not named keep the values of the scene; the density of an anomaly is that of
    its host under the same hypothesis times (1 + drho_rel). The mappings are read-only.
    """

    name: str
    densities_g_cm3: Mapping[str, float] = field(default_factory=dict)
    drho_rel: Mapping[str, float] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not isinstance(self.name, str) or not self.name.strip():
            raise ValueError("a hypothesis needs a name")
        object.__setattr__(self, "densities_g_cm3", _frozen(self.densities_g_cm3, "densities_g_cm3", 0.0))
        object.__setattr__(self, "drho_rel", _frozen(self.drho_rel, "drho_rel", -1.0))

    def _key(self):
        return self.name, tuple(self.densities_g_cm3.items()), tuple(self.drho_rel.items())

    def __eq__(self, other) -> bool:
        if not isinstance(other, Hypothesis):
            return NotImplemented
        return self._key() == other._key()

    def __hash__(self) -> int:
        return hash(self._key())

    @classmethod
    def as_built(cls, name: str = "as built") -> Hypothesis:
        """The scene as it was built: no density is changed."""
        return cls(name)

    @classmethod
    def without(cls, anomalies: Iterable[str], name: str = "none") -> Hypothesis:
        """No anomaly: the named anomalies have the density of their hosts (drho_rel = 0)."""
        return cls(name, drho_rel={anomaly: 0.0 for anomaly in _names(anomalies)})

    @classmethod
    def void(cls, anomalies: Iterable[str], name: str = "void") -> Hypothesis:
        """The named anomalies are voids (drho_rel = −1)."""
        return cls(name, drho_rel={anomaly: -1.0 for anomaly in _names(anomalies)})

    @classmethod
    def decompaction(
        cls, anomalies: Iterable[str], drho_rel: float, name: str = "decompaction"
    ) -> Hypothesis:
        """The named anomalies are zones of the relative density change `drho_rel`, e.g. −0.2."""
        return cls(name, drho_rel={anomaly: drho_rel for anomaly in _names(anomalies)})

    def opacity_by_material_g_cm2(self, paths: PathTable) -> dict:
        """Opacity in every solid material along the lines of a `PathTable` under this hypothesis, g/cm²."""
        try:
            return paths.opacity_by_material_g_cm2(dict(self.densities_g_cm3), dict(self.drho_rel))
        except ValueError as error:
            raise ValueError(f"hypothesis {self.name!r}: {error}") from error

    def opacity_g_cm2(self, paths: PathTable):
        """Plain opacity, the sum of density times length, under this hypothesis, g/cm²."""
        try:
            return paths.opacity_g_cm2(dict(self.densities_g_cm3), dict(self.drho_rel))
        except ValueError as error:
            raise ValueError(f"hypothesis {self.name!r}: {error}") from error

    def to_dict(self) -> dict:
        """Plain values for a run record."""
        return {"densities_g_cm3": dict(self.densities_g_cm3), "drho_rel": dict(self.drho_rel)}


def _names(anomalies) -> tuple[str, ...]:
    if isinstance(anomalies, str):
        return (anomalies,)
    names = tuple(anomalies)
    if not names:
        raise ValueError("name at least one anomaly")
    return names


class HypothesisSet:
    """The hypotheses of one map, in a fixed order and under different names.

    hypotheses   one `Hypothesis` or several

    `names` are the labels of the `hypothesis` dimension. A set is read by index or by name,
    iterated and measured; it cannot be changed.
    """

    def __init__(self, hypotheses) -> None:
        if isinstance(hypotheses, Hypothesis):
            hypotheses = (hypotheses,)
        self._items = tuple(hypotheses)
        if not self._items or not all(isinstance(h, Hypothesis) for h in self._items):
            raise TypeError("a hypothesis set holds at least one Hypothesis")
        names = [h.name for h in self._items]
        if len(set(names)) != len(names):
            raise ValueError("the hypotheses of a set have different names")

    @property
    def names(self) -> tuple[str, ...]:
        """Names of the hypotheses, in order."""
        return tuple(h.name for h in self._items)

    def __iter__(self) -> Iterator[Hypothesis]:
        return iter(self._items)

    def __len__(self) -> int:
        return len(self._items)

    def __getitem__(self, key) -> Hypothesis:
        if isinstance(key, str):
            for hypothesis in self._items:
                if hypothesis.name == key:
                    return hypothesis
            raise KeyError(key)
        return self._items[key]

    def __eq__(self, other) -> bool:
        if not isinstance(other, HypothesisSet):
            return NotImplemented
        return self._items == other._items

    def __hash__(self) -> int:
        return hash(self._items)

    def __repr__(self) -> str:
        return f"HypothesisSet({', '.join(self.names)})"

    def check(self, paths: PathTable) -> None:
        """ValueError if a hypothesis names a body the path table does not have or cannot change."""
        for hypothesis in self._items:
            try:
                paths.densities_g_cm3(dict(hypothesis.densities_g_cm3), dict(hypothesis.drho_rel))
            except ValueError as error:
                raise ValueError(f"hypothesis {hypothesis.name!r}: {error}") from error
