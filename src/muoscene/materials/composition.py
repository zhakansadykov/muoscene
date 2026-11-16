"""materials.composition — what a material is made of: elements, mixtures and moist soils.

Three kinds of material:

- `Material` — a base material: element mass fractions, the density of the solid and the mean
  excitation energy I. It becomes a Mulder base material.
- `Mixture` — base materials mixed by mass, each keeping its own density and I. It becomes a
  Mulder composite. Mulder computes the density effect of each part at the density of that part.
- `Soil` — a mineral skeleton and pore water, given as in an engineering-geological survey: dry
  density and gravimetric water content. Its material is a `Mixture` of skeleton and water; its
  bulk density goes to the geometry.

Why a soil is a mixture and not one material at its bulk density: the density effect of
ionisation loss builds up over tens of nanometres (c/ω_p is 6 nm in quartz), far below the size
of sand and silt grains, so the grains, not the bulk, set it. The same quartz defined at
1.6 g/cm³ instead of 2.65 g/cm³ has a CSDA range shorter by 1.1–1.5 % per g/cm² at 1–100 GeV
(`docs/notes/materials.md`).

Mean excitation energy: a measured or tabulated value where one exists; otherwise the rule of
Berger and Seltzer (1982), adopted in ICRU Report 37 — Bragg additivity with I = 19.2, 81, 82,
106, 112 and 180 eV for H, C, N, O, F and Cl bound in condensed compounds and 1.13 times the
elemental value for any other element. Mulder's own default (Bragg additivity of the elemental
values) puts I of water at 69 eV instead of 75–80 eV and shortens its range by 0.5 % at 1 GeV.

Element data (Z, A, elemental I) come from Mulder, so that the fractions computed here are
those Mulder uses. Densities in g/cm³, energies in eV for I.
"""

from __future__ import annotations

import math
import re
from collections.abc import Mapping
from dataclasses import dataclass
from functools import cache

import mulder.materials as mm

BERGER_SELTZER_EV = {"H": 19.2, "C": 81.0, "N": 82.0, "O": 106.0, "F": 112.0, "Cl": 180.0}
"""I of light elements bound in condensed compounds (Berger and Seltzer 1982, Table 5.1), eV."""
BERGER_SELTZER_FACTOR = 1.13
"""Factor on the elemental I of every other element in a compound (Berger and Seltzer 1982)."""


@dataclass(frozen=True)
class Element:
    symbol: str
    Z: int
    A: float  # g/mol
    i_ev: float  # mean excitation energy of the element


@cache
def element(symbol: str) -> Element:
    """Element data as Mulder has them; ValueError for an unknown symbol."""
    try:
        e = mm.Element(symbol)
    except Exception as error:  # noqa: BLE001 - Mulder raises a generic error for an unknown symbol
        raise ValueError(f"unknown element {symbol!r}") from error
    return Element(symbol, int(e.Z), float(e.A), float(e.I) * 1e9)


# ------------------------------------------------------------------ formulas

_TOKEN = re.compile(r"([A-Z][a-z]?)|(\()|(\))|(\d+(?:\.\d+)?)")


def parse_formula(formula: str) -> dict[str, float]:
    """Atoms per formula unit: "Al2Si2O5(OH)4" -> {"Al": 2, "Si": 2, "O": 9, "H": 4}."""
    tokens, pos = [], 0
    text = formula.replace(" ", "")
    while pos < len(text):
        m = _TOKEN.match(text, pos)
        if m is None:
            raise ValueError(f"cannot read formula {formula!r} at {text[pos:]!r}")
        tokens.append(m)
        pos = m.end()

    def group(i: int, depth: int) -> tuple[dict[str, float], int]:
        atoms: dict[str, float] = {}
        while i < len(tokens):
            symbol, opening, closing, _ = tokens[i].groups()
            if closing:
                if depth == 0:
                    raise ValueError(f"unbalanced ')' in {formula!r}")
                return atoms, i + 1
            if opening:
                inner, i = group(i + 1, depth + 1)
            elif symbol:
                element(symbol)
                inner, i = {symbol: 1.0}, i + 1
            else:
                raise ValueError(f"a number without an element in {formula!r}")
            count = 1.0
            if i < len(tokens) and tokens[i].group(4):
                count, i = float(tokens[i].group(4)), i + 1
            for s, n in inner.items():
                atoms[s] = atoms.get(s, 0.0) + n * count
        if depth:
            raise ValueError(f"unbalanced '(' in {formula!r}")
        return atoms, i

    if not tokens:
        raise ValueError("empty formula")
    atoms, _ = group(0, 0)
    return atoms


def formula_mass_fractions(formula: str) -> dict[str, float]:
    """Element mass fractions of a chemical formula."""
    atoms = parse_formula(formula)
    mass = {s: n * element(s).A for s, n in atoms.items()}
    total = sum(mass.values())
    return {s: m / total for s, m in mass.items()}


def oxide_mass_fractions(oxides: Mapping[str, float]) -> dict[str, float]:
    """Element mass fractions of a rock or soil analysis given as oxides (any mass unit, e.g. wt%).

    The analysis is normalised to its own total. Volatiles (CO2, H2O) count only if they are
    given as such; a loss on ignition of unknown nature is not accepted.
    """
    if any(k.strip().upper() == "LOI" for k in oxides):
        raise ValueError("give the loss on ignition as CO2 and H2O, or leave it out")
    out: dict[str, float] = {}
    total = sum(oxides.values())
    if total <= 0 or any(v < 0 for v in oxides.values()):
        raise ValueError("oxide contents must be non-negative with a positive total")
    for formula, amount in oxides.items():
        for s, f in formula_mass_fractions(formula).items():
            out[s] = out.get(s, 0.0) + f * amount / total
    return out


# ------------------------------------------------------------------ I and Z/A


def z_over_a(mass_fractions: Mapping[str, float]) -> float:
    """<Z/A> = sum of w_j Z_j / A_j, mol/g."""
    return sum(w * element(s).Z / element(s).A for s, w in mass_fractions.items())


def bragg_mean_excitation(mass_fractions: Mapping[str, float], i_ev: Mapping[str, float]) -> float:
    """Bragg additivity: ln I = sum w_j (Z/A)_j ln I_j / <Z/A>."""
    num = den = 0.0
    for s, w in mass_fractions.items():
        za = w * element(s).Z / element(s).A
        num += za * math.log(i_ev[s])
        den += za
    return math.exp(num / den)


def berger_seltzer_mean_excitation(mass_fractions: Mapping[str, float]) -> float:
    """Mean excitation energy (eV) of a condensed compound or mixture without a measured value."""
    values = {s: BERGER_SELTZER_EV.get(s, BERGER_SELTZER_FACTOR * element(s).i_ev) for s in mass_fractions}
    return bragg_mean_excitation(mass_fractions, values)


# ------------------------------------------------------------------ materials


def _normalised(fractions: Mapping[str, float], what: str) -> tuple[tuple[str, float], ...]:
    if not fractions:
        raise ValueError(f"{what}: empty composition")
    if any(w < 0 for w in fractions.values()):
        raise ValueError(f"{what}: negative mass fraction")
    total = sum(fractions.values())
    if total <= 0:
        raise ValueError(f"{what}: mass fractions sum to zero")
    return tuple(sorted((k, w / total) for k, w in fractions.items() if w > 0))


@dataclass(frozen=True)
class Material:
    """A base material.

    `i_ev` is the mean excitation energy; None leaves it to Mulder (Bragg additivity of the
    elemental values, too low for condensed compounds) and exists to reproduce old numbers.
    """

    name: str
    mass_fractions: tuple[tuple[str, float], ...]
    density_g_cm3: float
    i_ev: float | None
    source: str = ""

    def __post_init__(self) -> None:
        symbols = [s for s, _ in self.mass_fractions]
        if len(set(symbols)) != len(symbols):
            raise ValueError(f"{self.name}: an element is listed twice")
        object.__setattr__(self, "mass_fractions", _normalised(dict(self.mass_fractions), self.name))
        for s, _ in self.mass_fractions:
            element(s)
        if not self.density_g_cm3 > 0:
            raise ValueError(f"{self.name}: density must be positive")
        if self.i_ev is not None and not self.i_ev > 0:
            raise ValueError(f"{self.name}: I must be positive")

    @classmethod
    def from_elements(
        cls,
        name: str,
        mass_fractions: Mapping[str, float],
        density_g_cm3: float,
        i_ev: float | str = "berger-seltzer",
        source: str = "",
    ) -> Material:
        """From element mass fractions; `i_ev` a value in eV, "berger-seltzer" or "mulder"."""
        fractions = dict(_normalised(mass_fractions, name))
        return cls(name, tuple(fractions.items()), density_g_cm3, _resolve_i(i_ev, fractions), source)

    @classmethod
    def from_formula(
        cls,
        name: str,
        formula: str,
        density_g_cm3: float,
        i_ev: float | str = "berger-seltzer",
        source: str = "",
    ) -> Material:
        """From a chemical formula, e.g. "CaCO3"."""
        return cls.from_elements(name, formula_mass_fractions(formula), density_g_cm3, i_ev, source)

    @classmethod
    def from_oxides(
        cls,
        name: str,
        oxides: Mapping[str, float],
        density_g_cm3: float,
        i_ev: float | str = "berger-seltzer",
        source: str = "",
    ) -> Material:
        """From a chemical analysis in oxides, e.g. {"SiO2": 66.6, "Al2O3": 15.4, ...} in wt%."""
        return cls.from_elements(name, oxide_mass_fractions(oxides), density_g_cm3, i_ev, source)

    @property
    def elements(self) -> dict[str, float]:
        return dict(self.mass_fractions)

    @property
    def z_over_a(self) -> float:
        return z_over_a(self.elements)

    @property
    def mean_excitation_ev(self) -> float:
        """I as Mulder will use it."""
        if self.i_ev is not None:
            return self.i_ev
        return bragg_mean_excitation(self.elements, {s: element(s).i_ev for s in self.elements})


def _resolve_i(i_ev: float | str, fractions: Mapping[str, float]) -> float | None:
    if i_ev == "berger-seltzer":
        return berger_seltzer_mean_excitation(fractions)
    if i_ev == "mulder":
        return None
    if isinstance(i_ev, str):
        raise ValueError(f'i_ev must be a number, "berger-seltzer" or "mulder", not {i_ev!r}')
    return float(i_ev)


@dataclass(frozen=True)
class Mixture:
    """Base materials mixed by mass; the density is that of the solid mixture, without pores."""

    name: str
    parts: tuple[tuple[Material, float], ...]
    source: str = ""

    def __post_init__(self) -> None:
        if not self.parts:
            raise ValueError(f"{self.name}: empty mixture")
        if not all(isinstance(m, Material) for m, _ in self.parts):
            raise TypeError(f"{self.name}: a mixture is made of base materials")
        if any(w < 0 for _, w in self.parts):
            raise ValueError(f"{self.name}: negative mass fraction")
        total = sum(w for _, w in self.parts)
        if total <= 0:
            raise ValueError(f"{self.name}: mass fractions sum to zero")
        merged: dict[Material, float] = {}
        for m, w in self.parts:
            if w > 0:
                merged[m] = merged.get(m, 0.0) + w / total
        object.__setattr__(self, "parts", tuple(merged.items()))

    @classmethod
    def of(cls, name: str, parts: Mapping[Material, float], source: str = "") -> Mixture:
        return cls(name, tuple(parts.items()), source)

    @property
    def density_g_cm3(self) -> float:
        return 1.0 / sum(w / m.density_g_cm3 for m, w in self.parts)

    @property
    def elements(self) -> dict[str, float]:
        out: dict[str, float] = {}
        for m, w in self.parts:
            for s, f in m.mass_fractions:
                out[s] = out.get(s, 0.0) + w * f
        return out

    @property
    def z_over_a(self) -> float:
        return z_over_a(self.elements)

    @property
    def mean_excitation_ev(self) -> float:
        """Bragg combination of the I of the parts (for reports; Mulder treats the parts separately)."""
        num = den = 0.0
        for m, w in self.parts:
            za = w * m.z_over_a
            num += za * math.log(m.mean_excitation_ev)
            den += za
        return math.exp(num / den)


AnyMaterial = Material | Mixture


@dataclass(frozen=True)
class Soil:
    """A mineral skeleton with pore water, as an engineering-geological survey describes it.

    dry_density_g_cm3   ρ_d, mass of solids per bulk volume
    water_content_percent  w, mass of water per mass of solids (gravimetric, ГОСТ 5180), %
    The density of the skeleton material is the particle density ρ_s; pore air is neglected
    (its mass is below 1e-3 of the soil).
    """

    skeleton: Material | Mixture
    dry_density_g_cm3: float
    water_content_percent: float
    water: Material

    def __post_init__(self) -> None:
        if not self.dry_density_g_cm3 > 0:
            raise ValueError("dry density must be positive")
        if not self.water_content_percent >= 0:
            raise ValueError("water content must be non-negative")
        if self.dry_density_g_cm3 >= self.skeleton.density_g_cm3:
            raise ValueError(
                f"dry density {self.dry_density_g_cm3:g} g/cm3 is not below the particle density"
                f" {self.skeleton.density_g_cm3:g} g/cm3 of {self.skeleton.name}"
            )
        if self.saturation > 1 + 1e-9:
            limit = 100 * self.porosity * self.water.density_g_cm3 / self.dry_density_g_cm3
            raise ValueError(
                f"water content {self.water_content_percent:g} % exceeds saturation ({limit:.1f} %)"
                f" at dry density {self.dry_density_g_cm3:g} g/cm3 and particle density"
                f" {self.skeleton.density_g_cm3:g} g/cm3"
            )

    @classmethod
    def from_bulk(
        cls,
        skeleton: Material | Mixture,
        bulk_density_g_cm3: float,
        water_mass_fraction: float,
        water: Material,
    ) -> Soil:
        """From bulk density and the mass fraction of water in the moist soil (the inverse view)."""
        if not 0 <= water_mass_fraction < 1:
            raise ValueError("water mass fraction must be in [0, 1)")
        w = water_mass_fraction / (1 - water_mass_fraction)
        return cls(skeleton, bulk_density_g_cm3 / (1 + w), 100 * w, water)

    @property
    def w(self) -> float:
        """Water content as a fraction."""
        return self.water_content_percent / 100

    @property
    def bulk_density_g_cm3(self) -> float:
        return self.dry_density_g_cm3 * (1 + self.w)

    @property
    def porosity(self) -> float:
        return 1 - self.dry_density_g_cm3 / self.skeleton.density_g_cm3

    @property
    def volumetric_water_content(self) -> float:
        return self.w * self.dry_density_g_cm3 / self.water.density_g_cm3

    @property
    def saturation(self) -> float:
        return self.volumetric_water_content / self.porosity

    @property
    def water_mass_fraction(self) -> float:
        return self.w / (1 + self.w)

    @property
    def material(self) -> Material | Mixture:
        """Material of the moist soil: skeleton and water mixed by mass; the skeleton if dry."""
        if self.w == 0:
            return self.skeleton
        parts = self.skeleton.parts if isinstance(self.skeleton, Mixture) else ((self.skeleton, 1.0),)
        scaled = [(m, f * (1 - self.water_mass_fraction)) for m, f in parts]
        name = f"{self.skeleton.name}, w = {self.water_content_percent:g} %"
        return Mixture(name, (*scaled, (self.water, self.water_mass_fraction)))
