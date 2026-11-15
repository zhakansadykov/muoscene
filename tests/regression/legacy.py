"""Physics of the preliminary sensitivity calculation, ported as pure functions.

The preliminary calculation (October 2026) gave the sensitivity numbers of the 2027 grant
application: emulsion detectors in niches cut into a vertical canyon wall, spherical voids and
decompaction zones, a box for a pile cap. It was three scripts with module-level parameters,
mutated between calls. Here the same arithmetic takes its parameters explicitly (`Model`), so
that the numbers can be frozen and reproduced. Expressions and their order are kept as in the
scripts: the port is checked to give bit-identical results against them (`docs/legacy.md`).

Geometry, metres: vertical wall at x = 0, ground at x > 0, canyon at x < 0. The ground surface
is z = 0 at the wall, rising linearly by `rise` over `xr` and flat beyond. A niche of depth `xn`,
width `w` and height `hn` is cut into the wall at depth `h` below the brow; the detector sits on
its floor at x = `xd`. Directions point from the detector to the sky.

Physics: Guan et al. (2015) sea-level spectrum in total energy, tabulated as the integral
intensity above E on a grid of cos(zenith) from 0.25 (clipped below) and ln E; energy threshold
E(X) = a/b (exp(bX) - 1) + 0.1 GeV for opacity X; straight lines; no scattering. Significance:
the best of five cones around the target, with a background term from eight neighbouring bins.

Units: metres, g/cm^3, g/cm^2, GeV (total energy), seconds, cm^2 for detector area; intensities
in cm^-2 s^-1 sr^-1.
"""

from dataclasses import dataclass, field

import numpy as np

# ---------------------------------------------------------------- spectrum and its table
P1, P2, P3, P4, P5 = 0.102573, -0.068287, 0.958633, 0.0407253, 0.817285


def cos_star(c):
    return np.sqrt((c**2 + P1**2 + P2 * c**P3 + P4 * c**P5) / (1 + P1**2 + P2 + P4))


def dnde(E, c):
    """Guan et al. (2015) differential intensity, cm^-2 s^-1 sr^-1 GeV^-1; total energy E in GeV."""
    cs = cos_star(c)
    Es = E * (1 + 3.64 / (E * cs**1.29))
    return 0.14 * Es**-2.7 * (1 / (1 + 1.1 * E * cs / 115) + 0.054 / (1 + 1.1 * E * cs / 850))


COS_GRID = np.linspace(0.25, 1.0, 301)
E_GRID = np.logspace(-1, 4, 2001)
LN_E = np.log(E_GRID)


def _integral_table():
    table = np.empty((COS_GRID.size, E_GRID.size))
    for i, c in enumerate(COS_GRID):
        f = dnde(E_GRID, c) * E_GRID  # integral over ln E
        cum = np.concatenate([np.cumsum(((f[1:] + f[:-1]) / 2 * np.diff(LN_E))[::-1])[::-1], [0.0]])
        table[i] = np.log(np.maximum(cum, 1e-30))
    return table


LN_INTEGRAL = _integral_table()


def integral_intensity(e_th, c):
    """Integral intensity above total energy `e_th`, cm^-2 s^-1 sr^-1.

    Nearest node in cos(zenith), clipped to [0.25, 1]; linear interpolation of the logarithm
    in ln E, clipped to [0.1, 1e4] GeV.
    """
    ci = np.clip(np.rint((c - COS_GRID[0]) / (COS_GRID[1] - COS_GRID[0])).astype(int), 0, COS_GRID.size - 1)
    x = np.clip(np.log(e_th), LN_E[0], LN_E[-1])
    j = np.clip(np.searchsorted(LN_E, x) - 1, 0, LN_E.size - 2)
    t = (x - LN_E[j]) / (LN_E[j + 1] - LN_E[j])
    return np.exp(LN_INTEGRAL[ci, j] * (1 - t) + LN_INTEGRAL[ci, j + 1] * t)


# ---------------------------------------------------------------- parameters
@dataclass(frozen=True)
class Niche:
    xn: float = 0.5  # depth into the wall, m
    w: float = 0.3  # width, m
    hn: float = 0.3  # height, m
    xd: float = 0.25  # detector distance from the wall, m


@dataclass(frozen=True)
class Model:
    """Parameters of the preliminary calculation; defaults are the values it used."""

    rho: float = 1.8  # ground density, g/cm^3
    area_cm2: float = 240.0  # emulsion area at a point
    t_exp_s: float = 60 * 86400.0  # exposure
    eff: float = 0.9  # fraction of reconstructed tracks
    tmax: float = 1.0  # acceptance, tan of the angle to the film normal
    neigh: int = 8  # neighbouring bins in the background estimate
    a_loss: float = 2.2e-3  # GeV per g/cm^2
    b_loss: float = 4.0e-6  # per g/cm^2
    offset_gev: float = 0.1
    niche: Niche = field(default_factory=Niche)


def threshold(m: Model, X):
    """Total energy needed to cross opacity X (g/cm^2), GeV."""
    return m.a_loss / m.b_loss * np.expm1(m.b_loss * X) + m.offset_gev


# ---------------------------------------------------------------- geometry
def soil_length(m: Model, u, h, rise=0.0, xr=1.0):
    """Path length in the ground (m) from the detector on the niche floor; u: (N, 3) unit vectors."""
    xd, xn, w, hn = m.niche.xd, m.niche.xn, m.niche.w, m.niche.hn
    ax, ay, c = u[:, 0], u[:, 1], u[:, 2]
    big = 1e9
    t_roof = hn / c
    t_side = np.where(np.abs(ay) > 1e-12, (w / 2) / np.maximum(np.abs(ay), 1e-12), big)
    t_back = np.where(ax > 1e-12, (xn - xd) / np.maximum(ax, 1e-12), big)
    t_front = np.where(ax < -1e-12, xd / np.maximum(-ax, 1e-12), big)
    t_exit = np.minimum.reduce([t_roof, t_side, t_back, t_front])
    via_front = t_front <= t_exit
    if rise > 0:
        den = c - rise * ax / xr
        t1 = np.where(den > 1e-9, (h + rise * xd / xr) / np.maximum(den, 1e-9), big)
        x1 = xd + t1 * ax
        t2 = (h + rise) / c
        t_surf = np.where((t1 < big) & (x1 <= xr), t1, t2)
        t_surf = np.where(ax <= 0, h / c, t_surf)
    else:
        t_surf = h / c
    t_end = np.where(ax < 0, np.minimum(t_surf, t_front), t_surf)
    return np.where(via_front, 0.0, np.maximum(t_end - t_exit, 0.0))


def sphere_chord(u, D, C, r):
    """Chord (m) of rays from D along u through the sphere of centre C and radius r."""
    v = C - D
    tca = u @ v
    d2 = (v @ v) - tca**2
    ok = (tca > 0) & (d2 < r * r)
    return np.where(ok, 2 * np.sqrt(np.maximum(r * r - d2, 0)), 0.0)


def box_chord(u, D, lo, hi):
    """Chord (m) of rays from D along u through the box [lo, hi]."""
    with np.errstate(divide="ignore", invalid="ignore"):
        inv = 1.0 / u
        t1 = (lo - D) * inv
        t2 = (hi - D) * inv
    tmin = np.nanmax(np.minimum(t1, t2), axis=1)
    tmax = np.nanmin(np.maximum(t1, t2), axis=1)
    tmin = np.maximum(tmin, 0)
    return np.clip(tmax - tmin, 0, None)


def cone_dirs(n, alpha, N, rng):
    """N random directions uniform in solid angle within a cone of half-angle alpha around n."""
    cb = rng.uniform(np.cos(alpha), 1.0, N)
    sb = np.sqrt(1 - cb**2)
    ps = rng.uniform(0, 2 * np.pi, N)
    n = n / np.linalg.norm(n)
    t = np.array([1.0, 0, 0]) if abs(n[0]) < 0.9 else np.array([0, 1.0, 0])
    e1 = np.cross(n, t)
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(n, e1)
    return cb[:, None] * n + (sb * np.cos(ps))[:, None] * e1 + (sb * np.sin(ps))[:, None] * e2


def detector(m: Model, h):
    return np.array([m.niche.xd, 0.0, -h])


# ---------------------------------------------------------------- horizontal films
def counts(m: Model, u, h, rise, xr, extra_X=None, d_omega=1.0):
    """Expected tracks over directions u, each carrying solid angle d_omega; horizontal films."""
    c = u[:, 2]
    acc = (c > 0) & (np.hypot(u[:, 0], u[:, 1]) <= m.tmax * c)
    X = m.rho * 100 * soil_length(m, u, h, rise, xr)
    if extra_X is not None:
        X = X + extra_X
    intensity = integral_intensity(threshold(m, np.maximum(X, 0)), c)
    return np.sum(np.where(acc, intensity * c, 0.0)) * d_omega * m.area_cm2 * m.t_exp_s * m.eff


CONE_FACTORS = (0.6, 0.8, 1.0, 1.3, 1.7)


def evaluate_sphere(m: Model, h, C, r, drho, rise=0.0, xr=1.0, N=40000, seed=1):
    """Significance of a sphere (density change drho, g/cm^3) seen by horizontal films."""
    rng = np.random.default_rng(seed)
    D = detector(m, h)
    n = C - D
    dist = np.linalg.norm(n)
    n /= dist
    zen = np.degrees(np.arccos(n[2]))
    if np.hypot(n[0], n[1]) > m.tmax * n[2]:
        return dict(zen=zen, dist=dist, inside=False)
    ang_r = np.arcsin(min(r / dist, 1))
    best = None
    for k in CONE_FACTORS:
        alpha = max(np.radians(1.0), k * ang_r)
        u = cone_dirs(n, alpha, N, rng)
        d_omega = 2 * np.pi * (1 - np.cos(alpha)) / N
        n0 = counts(m, u, h, rise, xr, None, d_omega)
        n1 = counts(m, u, h, rise, xr, drho * 100 * sphere_chord(u, D, C, r), d_omega)
        s = abs(n1 - n0) / np.sqrt(n0) / np.sqrt(1 + 1 / m.neigh) if n0 > 0 else 0
        if best is None or s > best["sig"]:
            best = dict(
                zen=zen,
                dist=dist,
                inside=True,
                alpha=np.degrees(alpha),
                N0=n0,
                eff=(n1 / n0 - 1) * 100,
                sig=s,
            )
    return best


def evaluate_box(m: Model, h, lo, hi, drho, rise, xr, N=60000, seed=3):
    """Significance of a box (pile cap) seen by horizontal films; no background term."""
    rng = np.random.default_rng(seed)
    D = detector(m, h)
    C = (lo + hi) / 2
    n = C - D
    n /= np.linalg.norm(n)
    corners = np.array([[a, b, c] for a in (lo[0], hi[0]) for b in (lo[1], hi[1]) for c in (lo[2], hi[2])])
    v = corners - D
    v /= np.linalg.norm(v, axis=1)[:, None]
    alpha = np.max(np.arccos(np.clip(v @ n, -1, 1)))
    u = cone_dirs(n, alpha, N, rng)
    d_omega = 2 * np.pi * (1 - np.cos(alpha)) / N
    n0 = counts(m, u, h, rise, xr, None, d_omega)
    n1 = counts(m, u, h, rise, xr, drho * 100 * box_chord(u, D, lo, hi), d_omega)
    zen = np.degrees(np.arccos(n[2]))
    return dict(
        zen=zen, alpha=np.degrees(alpha), N0=n0, eff=(n1 / n0 - 1) * 100, sig=abs(n1 - n0) / np.sqrt(n0)
    )


def track_density(m: Model, h, rise, xr, T, tmax=10.0, N=200000, seed=5):
    """Tracks per cm^2 of a horizontal film over T seconds, directions with cos(zenith) > 0.25."""
    rng = np.random.default_rng(seed)
    c = rng.uniform(0.0, 1.0, N)
    s = np.sqrt(1 - c**2)
    p = rng.uniform(0, 2 * np.pi, N)
    u = np.stack([s * np.cos(p), s * np.sin(p), c], 1)
    u = u[(c > 0.25) & (s <= tmax * c)]
    d_omega = 2 * np.pi / N
    X = m.rho * 100 * soil_length(m, u, h, rise, xr) if h > 0 else np.zeros(len(u))
    intensity = integral_intensity(threshold(m, X), u[:, 2])
    return np.sum(intensity * u[:, 2]) * d_omega * T


# ---------------------------------------------------------------- films tilted towards the massif
def normal(beta):
    """Film normal tilted by beta degrees from the vertical towards +x (into the massif)."""
    b = np.radians(beta)
    return np.array([np.sin(b), 0.0, np.cos(b)])


def counts_tilted(m: Model, u, h, rise, xr, beta, extra_X=None, d_omega=1.0):
    n = normal(beta)
    c = u[:, 2]
    ca = u @ n
    acc = (c > 0) & (ca > 0) & (np.sqrt(np.maximum(1 - ca**2, 0)) <= m.tmax * ca)
    X = m.rho * 100 * soil_length(m, u, h, rise, xr)
    if extra_X is not None:
        X = X + extra_X
    intensity = integral_intensity(threshold(m, np.maximum(X, 0)), np.maximum(c, 0.25))
    return np.sum(np.where(acc, intensity * ca, 0.0)) * d_omega * m.area_cm2 * m.t_exp_s * m.eff


def inside_tilted(m: Model, n_dir, beta):
    ca = n_dir @ normal(beta)
    return ca > 0 and np.sqrt(max(1 - ca**2, 0)) <= m.tmax * ca and n_dir[2] > 0


def evaluate_sphere_tilted(m: Model, h, C, r, drho, rise, xr, beta, N=40000, seed=1):
    """Significance of a sphere for tilted films; None when the target is outside the cone."""
    rng = np.random.default_rng(seed)
    D = detector(m, h)
    n = C - D
    dist = np.linalg.norm(n)
    n /= dist
    if not inside_tilted(m, n, beta):
        return None
    ang_r = np.arcsin(min(r / dist, 1))
    best = None
    for k in CONE_FACTORS:
        alpha = max(np.radians(1.0), k * ang_r)
        u = cone_dirs(n, alpha, N, rng)
        d_omega = 2 * np.pi * (1 - np.cos(alpha)) / N
        n0 = counts_tilted(m, u, h, rise, xr, beta, None, d_omega)
        n1 = counts_tilted(m, u, h, rise, xr, beta, drho * 100 * sphere_chord(u, D, C, r), d_omega)
        s = abs(n1 - n0) / np.sqrt(n0) / np.sqrt(1 + 1 / m.neigh) if n0 > 0 else 0
        if best is None or s > best:
            best = s
    return best


def evaluate_box_tilted(m: Model, h, lo, hi, drho, rise, xr, beta, N=60000, seed=3):
    """Significance and expected tracks of a box for tilted films."""
    rng = np.random.default_rng(seed)
    D = detector(m, h)
    C = (lo + hi) / 2
    n = C - D
    n /= np.linalg.norm(n)
    corners = np.array([[a, b, c] for a in (lo[0], hi[0]) for b in (lo[1], hi[1]) for c in (lo[2], hi[2])])
    v = corners - D
    v /= np.linalg.norm(v, axis=1)[:, None]
    alpha = np.max(np.arccos(np.clip(v @ n, -1, 1)))
    u = cone_dirs(n, alpha, N, rng)
    d_omega = 2 * np.pi * (1 - np.cos(alpha)) / N
    n0 = counts_tilted(m, u, h, rise, xr, beta, None, d_omega)
    n1 = counts_tilted(m, u, h, rise, xr, beta, drho * 100 * box_chord(u, D, lo, hi), d_omega)
    return (abs(n1 - n0) / np.sqrt(n0) if n0 > 0 else 0.0), n0


def track_density_tilted(m: Model, h, rise, xr, beta, T, N=400000, seed=5):
    """Tracks per cm^2 of a tilted film over T seconds, both sides, cos(zenith) > 0.05."""
    rng = np.random.default_rng(seed)
    c = rng.uniform(0.0, 1.0, N)
    s = np.sqrt(1 - c**2)
    p = rng.uniform(0, 2 * np.pi, N)
    u = np.stack([s * np.cos(p), s * np.sin(p), c], 1)
    ca = u @ normal(beta)
    keep = (c > 0.05) & (np.sqrt(np.maximum(1 - ca**2, 0)) <= 10 * np.abs(ca))
    u, ca = u[keep], ca[keep]
    X = m.rho * 100 * soil_length(m, u, h, rise, xr) if h > 0 else np.zeros(len(u))
    intensity = integral_intensity(threshold(m, X), np.maximum(u[:, 2], 0.25))
    return np.sum(intensity * np.abs(ca)) * (2 * np.pi / N) * T


def reference_fraction(m: Model, h, rise, xr, beta, N=200000, seed=7):
    """Fraction of accepted tracks arriving from the canyon side (x < 0): the open-sky reference."""
    rng = np.random.default_rng(seed)
    u = cone_dirs(normal(beta), np.radians(45), N, rng)
    u = u[u[:, 2] > 0]
    ca = u @ normal(beta)
    X = m.rho * 100 * soil_length(m, u, h, rise, xr)
    w = integral_intensity(threshold(m, X), np.maximum(u[:, 2], 0.25)) * ca
    return w[u[:, 0] < 0].sum() / w.sum()
