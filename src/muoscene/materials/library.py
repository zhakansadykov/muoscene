"""materials.library — named materials with their sources.

Every entry gives its composition, density and mean excitation energy I with the source. Soils
are a skeleton (the solid mineral part, at its particle density) and pore water; their dry
density and water content are site data and are given by the user (`loam`, `loess`, `sand`,
`soil`). The names of soils choose the chemistry of the skeleton; grain size does not enter the
muon physics.

PDG — Particle Data Group, Atomic and Nuclear Properties of Materials (pdg.lbl.gov), with the
muon tables of D. E. Groom, N. V. Mokhov, S. I. Striganov, Atomic Data and Nuclear Data Tables
78 (2001) 183. I of skeletons given by a chemical analysis follows Berger and Seltzer (1982),
see `composition`.
"""

from __future__ import annotations

from .composition import Material, Mixture, Soil

WATER = Material.from_formula(
    "water",
    "H2O",
    density_g_cm3=1.000,
    i_ev=79.7,
    source="PDG, water (liquid): I = 79.7 eV, Z/A = 0.55509",
)

STANDARD_ROCK = Material.from_elements(
    "standard rock",
    {"Rk": 1.0},  # Z = 11, A = 22, the pseudo-element of Mulder
    density_g_cm3=2.650,
    i_ev=136.4,
    source="PDG, standard rock: Z = 11, A = 22, 2.650 g/cm3, I = 136.4 eV (Groom et al. 2001)",
)

CONCRETE = Material.from_elements(
    "concrete",
    {
        "H": 0.010000,
        "C": 0.001000,
        "O": 0.529107,
        "Na": 0.016000,
        "Mg": 0.002000,
        "Al": 0.033872,
        "Si": 0.337021,
        "K": 0.013000,
        "Ca": 0.044000,
        "Fe": 0.014000,
    },
    density_g_cm3=2.300,
    i_ev=135.2,
    source="PDG, shielding concrete: 2.300 g/cm3, I = 135.2 eV, Z/A = 0.50274",
)

AIR = Material.from_elements(
    "air",
    {"C": 0.000124, "N": 0.755267, "O": 0.231781, "Ar": 0.012827},
    density_g_cm3=1.20479e-3,
    i_ev=85.7,
    source="PDG, air (dry, 1 atm): 1.205e-3 g/cm3 at 20 °C, I = 85.7 eV; scene air gets the site density",
)

QUARTZ = Material.from_formula(
    "quartz",
    "SiO2",
    density_g_cm3=2.65,
    i_ev=139.2,
    source="SiO2; I = 139.2 eV (PDG, silicon dioxide); 2.65 g/cm3, density of quartz",
)

CALCITE = Material.from_formula(
    "calcite",
    "CaCO3",
    density_g_cm3=2.71,
    i_ev=136.4,
    source="CaCO3; I = 136.4 eV (PDG, calcium carbonate); 2.71 g/cm3, density of calcite",
)

SILICATE_SKELETON = Material.from_oxides(
    "silicate soil skeleton",
    {
        "SiO2": 66.62,
        "TiO2": 0.64,
        "Al2O3": 15.40,
        "FeO": 5.04,  # total iron as FeO
        "MnO": 0.10,
        "MgO": 2.48,
        "CaO": 3.59,
        "Na2O": 3.27,
        "K2O": 2.80,
        "P2O5": 0.15,
    },
    density_g_cm3=2.70,
    source="upper continental crust, R. L. Rudnick, S. Gao, Treatise on Geochemistry (2003; 2nd ed. 2014), "
    "Table 3; particle density 2.70 g/cm3, typical of silty and clayey soils (2.60-2.80); "
    "I by Berger-Seltzer",
)

LOESS_SKELETON = Material.from_oxides(
    "loess skeleton",
    {
        "SiO2": 62.89,
        "TiO2": 0.58,
        "Al2O3": 12.59,
        "FeO": 4.33,  # total iron as FeO
        "MnO": 0.09,
        "MgO": 3.84,
        "CaO": 10.91,
        "Na2O": 1.67,
        "K2O": 2.32,
    },
    density_g_cm3=2.70,
    source="Chinese Loess Plateau average, Z. L. Ding et al., Geochim. Cosmochim. Acta 65 (2001) 901, "
    "as compiled in Y. Zhang et al., Atmosphere 15 (2024) 261, Table 3; volatile-free analysis, the CO2 of "
    "carbonate is not included; particle density 2.70 g/cm3; I by Berger-Seltzer",
)

SAND_SKELETON = QUARTZ

SKELETONS = {"loam": SILICATE_SKELETON, "loess": LOESS_SKELETON, "sand": SAND_SKELETON}

MATERIALS = {
    "water": WATER,
    "standard_rock": STANDARD_ROCK,
    "concrete": CONCRETE,
    "air": AIR,
    "quartz": QUARTZ,
    "calcite": CALCITE,
    "silicate_skeleton": SILICATE_SKELETON,
    "loess_skeleton": LOESS_SKELETON,
}


def soil(skeleton: Material | Mixture, dry_density_g_cm3: float, water_content_percent: float) -> Soil:
    """A soil of the given skeleton, dry density (g/cm³) and gravimetric water content (%)."""
    return Soil(skeleton, dry_density_g_cm3, water_content_percent, WATER)


def loam(dry_density_g_cm3: float, water_content_percent: float) -> Soil:
    """Loam (суглинок) and other silicate soils: skeleton of the upper continental crust."""
    return soil(SILICATE_SKELETON, dry_density_g_cm3, water_content_percent)


def loess(dry_density_g_cm3: float, water_content_percent: float) -> Soil:
    """Loess: carbonate-bearing silicate skeleton."""
    return soil(LOESS_SKELETON, dry_density_g_cm3, water_content_percent)


def sand(dry_density_g_cm3: float, water_content_percent: float) -> Soil:
    """Quartz sand."""
    return soil(SAND_SKELETON, dry_density_g_cm3, water_content_percent)
