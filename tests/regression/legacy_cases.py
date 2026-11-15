"""Cases frozen in `data/legacy_reference.json`: the preliminary calculation on a synthetic canyon.

The geometry is synthetic and differs from any real site: a flat bank and a bank whose brow
rises by 1.5 m over 8 m, niches 12 m and 22 m below the brow, a pile cap 6 m x 6 m x 2 m whose
near face is 9 m from the wall. One case, `reference_void`, is the textbook example of the
preliminary calculation: a 1 m void 3 m from the wall and 3 m below a flat brow, seen from a
niche 10 m below it with 240 cm^2 of horizontal films for 60 days.

Each case is a name, a function of `legacy` and its arguments; `compute` returns plain Python
numbers and lists, so the results go to JSON unchanged.
"""

import numpy as np

from tests.regression import legacy

DAY = 86400.0
T60 = 60 * DAY
RHO = legacy.Model().rho

BANKS = {
    "flat": {"rise": 0.0, "xr": 1.0},
    "rising": {"rise": 1.5, "xr": 8.0},
}
UPPER, LOWER = 12.0, 22.0
AREA_UPPER, AREA_LOWER = 240.0, 600.0
VOID, DECOMPACTION = -RHO, -0.2 * RHO
CONCRETE = 2.4 - RHO


def pile_cap(rise):
    """Pile cap at the surface: 6 m along the wall normal, 6 m wide, 2 m thick."""
    return np.array([9.0, -3.0, rise - 2.0]), np.array([15.0, 3.0, rise])


def direction_set():
    """169 directions on a slope grid up to tan 1.2, plus vertical, grazing and canyon-side ones."""
    t = np.linspace(-1.2, 1.2, 13)
    tx, ty = np.meshgrid(t, t)
    v = np.stack([tx.ravel(), ty.ravel(), np.ones(tx.size)], axis=1)
    extra = np.array(
        [
            [0.0, 0.0, 1.0],
            [3.0, 0.0, 1.0],
            [-3.0, 0.0, 1.0],
            [0.0, 3.0, 1.0],
            [1.0, 1.0, 0.05],
            [-0.2, 0.0, 1.0],
        ]
    )
    v = np.concatenate([v, extra])
    return v / np.linalg.norm(v, axis=1)[:, None]


def plain(value):
    """Numbers, arrays, tuples and dicts of them as JSON-ready Python values."""
    if isinstance(value, dict):
        return {k: plain(v) for k, v in value.items()}
    if isinstance(value, tuple | list):
        return [plain(v) for v in value]
    if isinstance(value, np.ndarray):
        return [plain(v) for v in value.tolist()]
    if isinstance(value, np.bool_ | bool):
        return bool(value)
    if value is None or isinstance(value, str):
        return value
    return float(value)


def _cases():
    m = legacy.Model()
    up = legacy.Model(area_cm2=AREA_UPPER)
    low = legacy.Model(area_cm2=AREA_LOWER)
    u = direction_set()
    out = {}

    energies = np.array([0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 30.0, 100.0, 1000.0])
    cosines = np.array([0.25, 0.5, 0.7071, 0.9, 1.0])
    e, c = np.meshgrid(energies, cosines)
    out["integral_intensity"] = lambda: legacy.integral_intensity(e.ravel(), c.ravel())
    opacities = np.array([0.0, 100.0, 500.0, 1000.0, 1800.0, 3600.0, 5400.0, 1e4, 3e4])
    out["threshold"] = lambda: legacy.threshold(m, opacities)

    for bank, g in BANKS.items():
        for depth in (UPPER, LOWER):
            out[f"soil_length/{bank}/{depth:.0f}m"] = lambda g=g, depth=depth: legacy.soil_length(
                m, u, depth, **g
            )
    d = legacy.detector(m, UPPER)
    out["sphere_chord"] = lambda: legacy.sphere_chord(u, d, np.array([3.0, 0.0, -3.0]), 1.0)
    out["box_chord"] = lambda: legacy.box_chord(u, legacy.detector(m, LOWER), *pile_cap(0.0))

    out["reference_void"] = lambda: legacy.evaluate_sphere(up, 10.0, np.array([3.0, 0.0, -3.0]), 0.5, VOID)

    targets_upper = [
        ((3.0, 0.0, -3.0), 0.5, VOID),
        ((5.0, 0.0, -5.0), 0.5, VOID),
        ((3.0, 0.0, -3.0), 1.0, DECOMPACTION),
    ]
    targets_lower = [
        ((8.0, 0.0, -7.0), 1.0, VOID),
        ((8.0, 0.0, -12.0), 1.0, VOID),
        ((8.0, 0.0, -12.0), 2.0, DECOMPACTION),
    ]
    outside = ((25.0, 0.0, -3.0), 0.5, VOID)
    for bank, g in BANKS.items():
        for model, depth, targets in ((up, UPPER, [*targets_upper, outside]), (low, LOWER, targets_lower)):
            for centre, r, drho in targets:
                name = f"sphere/{bank}/{depth:.0f}m/{centre[0]:g},{centre[2]:g}/r{r:g}/drho{drho:+.2f}"
                out[name] = lambda model=model, depth=depth, centre=centre, r=r, drho=drho, g=g: (
                    legacy.evaluate_sphere(model, depth, np.array(centre), r, drho, **g)
                )
        lo, hi = pile_cap(g["rise"])
        out[f"box/{bank}/{LOWER:.0f}m"] = lambda lo=lo, hi=hi, g=g: legacy.evaluate_box(
            low, LOWER, lo, hi, CONCRETE, **g
        )
        for depth in (0.0, UPPER, LOWER):
            out[f"track_density/{bank}/{depth:.0f}m"] = lambda depth=depth, g=g: legacy.track_density(
                m, depth, g["rise"], g["xr"], T60
            )

        for beta in (0.0, 20.0):
            for centre, r, drho in targets_lower:
                where = f"{centre[0]:g},{centre[2]:g}/r{r:g}/drho{drho:+.2f}"
                name = f"sphere_tilted/{bank}/{LOWER:.0f}m/beta{beta:g}/{where}"
                out[name] = lambda centre=centre, r=r, drho=drho, g=g, beta=beta: (
                    legacy.evaluate_sphere_tilted(
                        low, LOWER, np.array(centre), r, drho, g["rise"], g["xr"], beta
                    )
                )
            out[f"box_tilted/{bank}/{LOWER:.0f}m/beta{beta:g}"] = lambda lo=lo, hi=hi, g=g, beta=beta: (
                legacy.evaluate_box_tilted(low, LOWER, lo, hi, CONCRETE, g["rise"], g["xr"], beta)
            )
            out[f"track_density_tilted/{bank}/{LOWER:.0f}m/beta{beta:g}"] = lambda g=g, beta=beta: (
                legacy.track_density_tilted(low, LOWER, g["rise"], g["xr"], beta, T60)
            )
            out[f"reference_fraction/{bank}/{LOWER:.0f}m/beta{beta:g}"] = lambda g=g, beta=beta: (
                legacy.reference_fraction(low, LOWER, g["rise"], g["xr"], beta)
            )
    return out


CASES = _cases()


def compute(name):
    """Result of one case as plain Python numbers."""
    return plain(CASES[name]())


def parameters():
    """Everything a reader needs to reproduce the cases by hand."""
    m = legacy.Model()
    return {
        "model": {k: getattr(m, k) for k in ("rho", "eff", "tmax", "neigh", "a_loss", "b_loss", "offset_gev")}
        | {"niche": vars(m.niche)},
        "banks": BANKS,
        "niche_depths_m": [UPPER, LOWER],
        "areas_cm2": {"upper": AREA_UPPER, "lower": AREA_LOWER},
        "exposure_days": 60,
        "pile_cap": "x from 9 to 15 m, y from -3 to 3 m, 2 m thick at the surface; density 2.4 g/cm^3",
        "direction_set": "slopes from -1.2 to 1.2 in 13 steps on both axes, plus six extra directions",
    }
