"""Units and physical constants used across the package.

Internal units: metres, g/cm^2 for opacity, GeV for energy, seconds for time.
Conversions to days and metres water equivalent are explicit helpers here.

Layer 0.
"""

# lengths
M = 1.0
CM = 1e-2
# opacity: g/cm^2; 1 m of water equivalent = 100 g/cm^2
MWE = 100.0
# energy: GeV
GEV = 1.0
MEV = 1e-3
MUON_MASS = 0.1056583755  # GeV (PDG 2024)
# time: s
S = 1.0
DAY = 86400.0


def density_length_to_opacity(rho_g_cm3: float, length_m: float) -> float:
    """Opacity X = rho * L in g/cm^2 for density in g/cm^3 and length in metres."""
    return rho_g_cm3 * length_m * 100.0
