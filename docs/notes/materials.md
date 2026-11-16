# Materials: choices and their size

Roadmap task 0.1-04, 16 November 2026. Mulder 0.3.7. Every number below is a ratio of CSDA ranges, or of the kinetic energies needed to cross an opacity X, computed by Mulder for materials defined through `muoscene.materials`. Opacities 500–3·10⁴ g/cm² cover the domain of version 0.1. The tests in `tests/unit/test_materials_*.py` hold the statements that matter.

## Summary

| Choice | Size | What muoscene does |
|---|---|---|
| Mean excitation energy I | Mulder's default (Bragg additivity of elemental values) gives water 69 eV instead of 75–80 eV: range −1.5 % at 0.1 GeV, −0.5 % at 1 GeV, −0.1 % at 10 GeV; for the silicate skeleton the threshold rises by 0.6 % at 500 g/cm² and 0.05 % at 3·10⁴ | I set explicitly: tabulated (PDG) where it exists, otherwise the Berger–Seltzer rule |
| Density effect in a porous soil | a soil defined as one material at its bulk density has a range 0.5–0.75 % shorter at 1–100 GeV and a threshold 0.5–0.8 % higher than the same soil as grains plus water; quartz at 1.6 instead of 2.65 g/cm³: −1.1 to −1.5 % | a soil is a mixture of skeleton (at particle density) and water; bulk density goes only to the geometry |
| Water content, at fixed opacity | threshold +0.9 % at w = 5 %, +1.6 % at 10 %, +3.0 % at 20 %, +4.0 % at 28 % (X = 500 g/cm²; a sixth less at 3·10⁴) — on top of the opacity itself, which grows as 1 + w | water is part of the composition; w and dry density are inputs |
| Chemistry of the skeleton | against the upper continental crust: loess +0.02 to +0.23 %, quartz +0.7 to +1.2 %, standard rock +0.8 to +1.0 %, calcite +0.9 to +1.1 %, concrete +2.0 to +2.3 % | a skeleton from a chemical analysis if the site has one; the library otherwise |
| Particle density of the skeleton | 2.75 against 2.65 g/cm³: −0.09 to −0.13 % | 2.70 by default; the survey value when there is one |
| Carbonate CO₂ left out of a volatile-free loess analysis | +0.26 % at 500 g/cm² to 0 at 3·10⁴ | analysis used as published |
| Mulder against the PDG range tables | standard rock within 0.34 %, concrete 0.24 %, water 0.52 %, 0.1 GeV – 1 TeV (with the PDG values of I); with Mulder's default I concrete −1.0 % and water −1.5 % at 0.1 GeV | checked by a test at 0.6 % |
| Concrete against loam | 0.995 at 7 GeV against the loam of the spike; 0.9991–0.9996 at 1–100 GeV against loam with 15 % water; 1.017–1.023 against dry silicate skeleton | weight w(K) = R_ref(K) / R_m(K) |

Net for a soil at a site: the water content is the material property that matters most, through bulk density (+10 % of opacity per 10 % of w) and then through composition (+1.6 % of threshold per 10 % of w). Composition of the skeleton is a 0.1–1 % effect; "standard rock" is not a good stand-in for a silicate soil at the 1 % level. All of these choices are at or below the 0.2–0.5 % by which Mulder differs from the PDG range tables, and two orders of magnitude below what a detector threshold and scattering do to a significance (`docs/numbers.md`).

## Mean excitation energy

PUMAS, when not given I, combines the elemental values by Bragg additivity. For condensed compounds that is about 10 % low (water: 69 eV against 75 eV of ICRU 37 and 79.7 eV of PDG; concrete: 121.7 eV against 135.2 eV). muoscene therefore always passes I:

- tabulated values from PDG for water (79.7 eV), standard rock (136.4 eV), shielding concrete (135.2 eV), dry air (85.7 eV), SiO₂ (139.2 eV), CaCO₃ (136.4 eV);
- otherwise the rule of M. J. Berger and S. M. Seltzer, *Stopping powers and ranges of electrons and positrons*, NBSIR 82-2550 (1982; second edition 82-2550-A, 1983), Section 5, adopted in ICRU Report 37 (1984): Bragg additivity with I = 19.2, 81, 82, 106, 112 and 180 eV for H, C, N, O, F and Cl bound in condensed compounds and 1.13 times the elemental value for other elements.

The rule against tabulated values: water 75.3 eV (ICRU 37: 75.0), CaCO₃ 136.4 (136.4), SiO₂ 141.0 (139.2), Al₂O₃ 142.6 (145.2), shielding concrete 136.3 (135.2). A 4 % change of I (water, 75 against 78 eV) moves the range by 0.15 % at 1 GeV and 0.05 % at 10 GeV.

## Density effect and porous soils

The density effect of ionisation loss depends on the electron density of the medium. Mulder computes it for each base material at that material's density. The screening length c/ω_p is 6 nm in quartz and 9 nm in water, and the field of the muon adjusts to the medium over a few times that, tens of nanometres (an order-of-magnitude estimate from the formation length 2cω/ω_p² at the frequencies of atomic excitation). Sand and silt grains are larger by two to five orders of magnitude, so the grains, not the bulk, set the effect. Clay platelets with their bound water are thinner than that scale; for them the homogeneous limit is closer, and it differs from the mixture by 0.15–0.20 % of stopping power at 1–100 GeV (quartz with 15 % of water as one material at the pore-free density against the composite). Mulder adds the stopping powers of the parts of a composite by mass, each at its own density: checked to round-off. A soil is therefore a `Mixture` of its mineral skeleton at the particle density ρ_s and water at 1 g/cm³ (a Mulder composite); the bulk density ρ = ρ_d (1 + w) goes only into the geometry. Pore air is neglected: its mass is below 10⁻³ of the soil.

| Soil as | Range at 1 / 10 / 100 GeV | Threshold at 500 / 3600 / 3·10⁴ g/cm² |
|---|---|---|
| one material at bulk density 1.84 g/cm³ (loam, ρ_d = 1.60, w = 15 %) against grains and water | −0.50 / −0.73 / −0.72 % | +0.48 / +0.76 / +0.82 % |

The concrete of the library keeps the PDG definition (one material at 2.30 g/cm³), as do the range tables it is compared with.

## Soil inputs

As in an engineering-geological survey (ГОСТ 5180): dry density ρ_d, gravimetric water content w (mass of water per mass of solids), and if measured the particle density ρ_s. Derived: bulk density ρ_d (1 + w), porosity n = 1 − ρ_d / ρ_s, volumetric water content θ = w ρ_d / ρ_w, saturation θ / n. A water content above saturation is refused. The round trip (ρ_d, w) → (ρ, water mass fraction) → (ρ_d, w) is exact to 10⁻¹³.

## Library

| Entry | Composition | Density, g/cm³ | I, eV | Source |
|---|---|---|---|---|
| water | H₂O | 1.000 | 79.7 | PDG |
| standard rock | Z = 11, A = 22 | 2.650 | 136.4 | PDG; Groom et al. (2001) |
| concrete | PDG shielding concrete (10 elements) | 2.300 | 135.2 | PDG |
| air | PDG dry air | 1.205·10⁻³ | 85.7 | PDG |
| quartz | SiO₂ | 2.65 | 139.2 | PDG (I) |
| calcite | CaCO₃ | 2.71 | 136.4 | PDG (I) |
| silicate skeleton (loam) | upper continental crust, 10 oxides | 2.70 | 148.2 (rule) | R. L. Rudnick, S. Gao, Treatise on Geochemistry, vol. 3 (2003), Table 3, which prints SiO₂ 66.6 and Al₂O₃ 15.4 |
| loess skeleton | Chinese Loess Plateau average, 9 oxides, volatile-free | 2.70 | 149.7 (rule) | Z. L. Ding et al., Geochim. Cosmochim. Acta 65 (2001) 901, as compiled in Y. Zhang et al., Atmosphere 15 (2024) 261, Table 3: an average of seven samples, not checked against the original paper |
| sand skeleton | quartz | 2.65 | 139.2 | — |

Particle densities: 2.65 g/cm³ for quartz sand and 2.70 g/cm³ for silty and clayey soils and loess, within the usual 2.60–2.80 g/cm³ of soil solids; a site survey replaces them. Iron is total iron as FeO. A survey or a site analysis replaces any of the skeletons through `Material.from_oxides`.

## Mulder in practice

- Names. Mulder keeps one registry of materials per process and refuses a second definition under a taken name, so muoscene names each material by a hash of its definition (`ms-<label>-<hash>`); definitions with the same label that are equal to 12 significant digits share a name, and the same definition under another label is a second material that compiles again. Mulder writes the names into its disk cache as bare TOML keys: a name with a dot never hits the cache again, and each process recompiles (10 s for the silicate skeleton). Names are therefore letters, digits and hyphens.
- Cost. A new soil compiles once per machine: 11 s for the 11-element silicate skeleton with water, 2.7 s for quartz sand (Linux VM, arm64); afterwards it loads from the cache in milliseconds. A dozen moisture hypotheses cost two minutes once.
- `inverse_range` interpolates the range table: range and inverse range round-trip to 6·10⁻⁶.
