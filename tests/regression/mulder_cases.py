"""Cases frozen in `data/mulder_reference.json`: Mulder and Calzone on the synthetic scenes of the spike.

The numbers of `docs/notes/mulder.md` (roadmap 0.1-02) in a form that a new release of Mulder or
Calzone can be checked against: the GCCLY15 reference, open-sky flux at sea level and at 1150 m,
energy thresholds and ranges, a uniform soil slab, the altitude cut-off, atmosphere models, the
geomagnetic field, the tables behind the transport, the path in the wall with a niche of
`legacy_cases`, and the reference void.

Energies: K is the kinetic energy at the detector (Mulder), E = K + m the total energy (the
preliminary calculation and the spike note). Units: m and kg/m^3 (Mulder), g/cm^2 for opacity;
fluxes in m^-2 s^-1 sr^-1 (GeV^-1 for differential ones). Every case uses continuous transport
and is deterministic. The physics of discrete transport enters through the frozen tables of
cross-sections and elastic scattering; discrete transport itself is checked by the validation
script of roadmap 0.1-13.

"All energies" above sea level starts at the altitude cut-off of the sea-level reference, where
the differential flux has an integrable singularity, so these integrals run in ln(K - K_cut).
A plain grid in ln K across the cut-off is wrong by up to 1.3 % even with 100 nodes per decade.

Materials are registered with Mulder under their own names (prefix `Reg`), so that a library of
muoscene materials can define `Loam` differently in the same process.
"""

import functools
import pathlib
import tempfile

import calzone
import mulder
import mulder.materials as mm
import numpy as np

from tests.regression import legacy, legacy_cases

MU = 0.1056583755  # GeV
CM = 100.0  # Calzone lengths are centimetres
SITE = {"latitude": 43.35, "longitude": 79.08}  # as in the spike; matters for the geomagnetic case only
ALTITUDE, HIGH = 1150.0, 3340.0  # m
SOIL = 1800.0  # kg/m^3
AIR = 1.2  # kg/m^3, density of the void
EPS = 1e-3  # m, the niche floor is this much below the detector
DATE = "2028-10-15"
ZENITH = (0.0, 30.0, 45.0, 60.0)
DEPTHS = (10.0, 20.0, 30.0)
OPACITIES = (500.0, 1000.0, 1800.0, 3600.0, 5400.0, 1e4, 2e4, 3e4)  # g/cm^2
RANGE_ENERGIES = (1.0, 7.0, 10.0, 100.0, 1000.0)  # K, GeV
RANGE_MATERIALS = ("Rock", "RegLoam", "RegConcrete", "Water")
TABLE_ENERGIES = (0.01, 0.1, 1.0, 10.0, 100.0, 1000.0)  # K, GeV
K_GRID = np.logspace(-3, 4, 71)  # kinetic energy at the detector for many directions, GeV
REFERENCE_VOID = (3.0, 0.0, -3.0, 0.5)  # x, y, z, radius, m; the `reference_void` of legacy_cases
VOID_DEPTH = 10.0
# Relative tolerance of the frozen numbers. Spread between platforms: 3e-9 between Linux arm64 and
# x86-64; on macOS arm64 up to 3.1e-6 in the elastic-scattering tables and 4e-9 elsewhere. The
# mildest change of physics tried, the pair-production model of Mulder, moves thresholds by 8e-5.
RTOL = 1e-5
GROUND = np.array([1500.0, 3000.0, 400.0])  # m: into the massif, along the wall, down from the brow
WORLD = np.array([4000.0, 4000.0, 420.0])  # m; the air of the scene ends 20 m above the brow
WORLD_Z = -190.0  # m, centre of the world box

BASE = {"RegQuartz": ("SiO2", 2650.0), "RegCalcite": ("CaCO3", 2710.0), "RegAlumina": ("Al2O3", 3950.0)}
LOAM = {"RegQuartz": 0.60, "RegCalcite": 0.15, "RegAlumina": 0.12, "Water": 0.13}  # mass fractions
CONCRETE = {"O": 0.529, "Si": 0.337, "Ca": 0.044, "Al": 0.034, "Na": 0.016, "Fe": 0.014}
CONCRETE |= {"K": 0.013, "H": 0.010, "Mg": 0.002, "C": 0.001}
CONCRETE_DENSITY = 2300.0


@functools.cache
def define_materials():
    """Register the materials of the cases with Mulder (once per process)."""
    for name, (formula, density) in BASE.items():
        mm.Material.define(name, composition=formula, density=density)
    mm.Composite.define("RegLoam", composition=LOAM)
    mm.Material.define("RegConcrete", composition=CONCRETE, density=CONCRETE_DENSITY)
    return True


def default_density(name):
    define_materials()
    try:
        return mm.Material(name).density
    except ValueError:  # not a base material
        return mm.Composite(name).density


# ------------------------------------------------------------------ integrals over energy


def integral(meter, k_min, k_max=1e5, n=801, **state):
    """Integral flux above the kinetic energy `k_min` (above the cut-off), trapezoids in ln K."""
    k = np.logspace(np.log10(k_min), np.log10(k_max), n)
    return float(np.trapezoid(np.asarray(meter.flux(energy=k, **state)) * k, np.log(k)))


def cutoff_of(meter, steps=60, **state):
    """Lowest K (GeV) with a non-zero flux, by bisection in ln K; 0 if the flux is positive at 1e-7 GeV."""
    if float(meter.flux(energy=1e-7, **state)) > 0:
        return 0.0
    lo, hi = np.log(1e-7), np.log(10.0)
    assert float(meter.flux(energy=np.exp(hi), **state)) > 0
    for _ in range(steps):
        mid = 0.5 * (lo + hi)
        if float(meter.flux(energy=np.exp(mid), **state)) > 0:
            hi = mid
        else:
            lo = mid
    return float(np.exp(hi))


def integral_all(meter, k_max=1e5, n=5601, **state):
    """Integral flux over all energies: from 1e-7 GeV, or from the cut-off in ln(K - K_cut)."""
    k_cut = cutoff_of(meter, **state)
    d = np.logspace(-9 if k_cut else -7, np.log10(k_max), n)
    flux = np.asarray(meter.flux(energy=k_cut + d, **state))
    return float(np.trapezoid(flux * d, np.log(d)))


def integral_many(meter, position, directions, k=K_GRID):
    """Integral flux for many directions at once on the grid `K_GRID`, rectangles in ln K."""
    n, m = len(directions), len(k)
    flux = meter.flux(
        position=np.broadcast_to(position, (n * m, 3)),
        direction=np.repeat(directions, m, axis=0),
        energy=np.tile(k, n),
    )
    return (np.asarray(flux).reshape(n, m) * k * np.gradient(np.log(k))).sum(axis=1)


def k_above(e_total):
    """Kinetic energy for a threshold given as total energy."""
    return e_total - MU


# ------------------------------------------------------------------ flux and physics


def reference_differential():
    """GCCLY15 at sea level, differential in K."""
    reference = mulder.Reference()
    k = np.array([0.5, 1.0, 10.0, 100.0, 1000.0])
    out = {"kinetic_energy_gev": k}
    for e in (90.0, 60.0, 30.0):
        out[f"elev{e:g}"] = np.asarray(reference.flux(energy=k, elevation=e, altitude=0.0))
    return out


def open_sky(altitude, elevation):
    """Open-sky integral flux over all energies and above total energies of 1, 5 and 10 GeV."""
    meter = mulder.Fluxmeter()
    state = dict(SITE, altitude=altitude, azimuth=0.0, elevation=elevation)
    out = {"all": integral_all(meter, **state)}
    return out | {f"E>{e:g}GeV": integral(meter, k_above(e), **state) for e in (1.0, 5.0, 10.0)}


def cutoff(altitude, elevation):
    """Altitude cut-off K_cut (GeV) of the sea-level reference, open sky, US standard atmosphere."""
    return cutoff_of(mulder.Fluxmeter(), **dict(SITE, altitude=altitude, azimuth=0.0, elevation=elevation))


def threshold(name):
    """Kinetic energy (GeV) that a muon needs to cross OPACITIES, from the CSDA range."""
    define_materials()
    compiled = mulder.Physics().compile(name)
    x = np.array(OPACITIES) * 10 / default_density(name)  # m at the default density
    return np.asarray(compiled.inverse_range(x))


def csda_range(name):
    """CSDA range (g/cm^2) at RANGE_ENERGIES."""
    define_materials()
    compiled = mulder.Physics().compile(name)
    return np.asarray(compiled.range(np.array(RANGE_ENERGIES))) * default_density(name) / 10


def physics_defaults():
    """Default physics models of Mulder (PUMAS)."""
    p = mulder.Physics()
    names = ("bremsstrahlung", "pair_production", "photonuclear", "cutoff", "elastic_ratio")
    return {a: getattr(p, a) for a in names}


def physics_tables(name):
    """Tables behind the transport at TABLE_ENERGIES, in the units of Mulder.

    Stopping power (GeV/m), proper time and transport path (m) at the default density, macroscopic
    cross-section of the discrete processes, elastic scattering (angle, path).
    (`CompiledMaterial.magnetic_gyration` crashes Mulder 0.3.7 and is left out.)
    """
    define_materials()
    compiled = mulder.Physics().compile(name)
    k = np.array(TABLE_ENERGIES)
    out = {
        f: np.asarray(getattr(compiled, f)(k)) for f in ("stopping_power", "proper_time", "transport_path")
    }
    out["cross_section"] = np.asarray(compiled.cross_section(k))
    elastic = np.asarray(compiled.elastic_scattering(k))
    return out | {f"elastic_{f}": elastic[f] for f in elastic.dtype.names}


def slab(altitude, depth, zenith):
    """Integral flux over all energies under `depth` of soil, surface at `altitude`."""
    define_materials()
    meter = mulder.Fluxmeter(mulder.Layer(altitude, material="RegLoam", density=SOIL))
    return integral_all(meter, **dict(SITE, altitude=altitude - depth, azimuth=0.0, elevation=90.0 - zenith))


def season(model):
    """Vertical integral flux at 1150 m for an atmosphere model: open sky and under 20 m of soil.

    Only the air between the site and sea level changes with the model; the reference does not.
    """
    define_materials()
    meter = mulder.Fluxmeter(mulder.Layer(ALTITUDE, material="RegLoam", density=SOIL), atmosphere=model)
    out = {}
    for where, altitude in (("open_sky", ALTITUDE + 0.01), ("under_20m", ALTITUDE - 20.0)):
        state = dict(SITE, altitude=altitude, azimuth=0.0, elevation=90.0)
        out[f"{where}/all"] = integral_all(meter, **state)
        out[f"{where}/E>1GeV"] = integral(meter, k_above(1.0), **state)
        out[f"{where}/cutoff_gev"] = cutoff_of(meter, **state)
    return out


def geomagnet():
    """Open-sky integral flux at 1150 m and 45 degrees elevation towards east and west, IGRF on DATE."""
    define_materials()
    meter = mulder.Fluxmeter(mulder.Layer(ALTITUDE, material="RegLoam", density=SOIL), date=DATE)
    out = {}
    for side, azimuth in (("east", 90.0), ("west", 270.0)):
        state = dict(SITE, altitude=ALTITUDE + 0.01, azimuth=azimuth, elevation=45.0)
        out[f"{side}/all"] = integral_all(meter, **state)
        out[f"{side}/E>1GeV"] = integral(meter, k_above(1.0), **state)
    return out


# ------------------------------------------------------------------ scenes


def write_prism(path, polygon, y0, y1):
    """Closed OBJ mesh (metres, with face normals): a convex x-z polygon extruded along y."""
    n = len(polygon)
    v = np.array([(x, y0, z) for x, z in polygon] + [(x, y1, z) for x, z in polygon], float)
    faces = []
    for i in range(n):
        j = (i + 1) % n
        faces += [(i, j, n + j), (i, n + j, n + i)]
    for i in range(1, n - 1):
        faces += [(0, i + 1, i), (n, n + i, n + i + 1)]
    centre = v.mean(axis=0)
    lines = ["v {:.6f} {:.6f} {:.6f}".format(*p) for p in v]
    rows = []
    for k, (a, b, c) in enumerate(faces, start=1):
        normal = np.cross(v[b] - v[a], v[c] - v[a])
        normal /= np.linalg.norm(normal)
        if normal @ (v[a] - centre) < 0:
            b, c, normal = c, b, -normal
        lines.append("vn {:.6f} {:.6f} {:.6f}".format(*normal))
        rows.append(f"f {a + 1}//{k} {b + 1}//{k} {c + 1}//{k}")
    pathlib.Path(path).write_text("\n".join(lines + rows) + "\n", encoding="utf-8")


def detector(depth):
    """Detector of the preliminary calculation: on the niche floor, `depth` below the brow."""
    return np.array([legacy.Niche().xd, 0.0, -depth])


def build_scene(depth, bank="flat", altitude=ALTITUDE, sphere=None):
    """Ground behind a vertical wall at x = 0 with a niche; brow flat (boxes) or rising (OBJ mesh).

    The ground reaches 1500 m into the massif and 1500 m along the wall each way, so that the
    grazing directions of `legacy_cases.direction_set` stay inside it. The air of the scene ends
    20 m above the brow: Mulder gives it sea-level density, and a thicker layer of it lowers the
    flux at altitude (spike, E12). The niche is that of the preliminary calculation with its floor
    EPS below the detector. sphere = (x, y, z, radius).
    """
    define_materials()
    niche = legacy.Niche()
    rise, xr = legacy_cases.BANKS[bank]["rise"], legacy_cases.BANKS[bank]["xr"]
    size = GROUND
    centre = np.array([size[0] / 2, 0.0, -size[2] / 2])
    with tempfile.TemporaryDirectory() as tmp:
        if rise > 0:
            mesh = pathlib.Path(tmp) / "brow.obj"
            polygon = [(0, -size[2]), (size[0], -size[2]), (size[0], rise), (xr, rise), (0, 0)]
            write_prism(mesh, polygon, -size[1] / 2, size[1] / 2)
            ground = {"mesh": {"path": str(mesh), "units": "m"}}
        else:
            ground = {"box": list(size * CM), "position": list(centre * CM)}
        ground |= {"material": "G4_SILICON_DIOXIDE", "subtract": ["Niche"]}
        if sphere is not None:
            assert rise == 0
            x, y, z, r = sphere
            where = (np.array([x, y, z]) - centre) * CM  # daughters are placed relative to the mother centre
            ground["Anomaly"] = {"sphere": {"radius": r * CM}, "position": list(where), "material": "G4_AIR"}
        z_niche = -depth - EPS + (niche.hn + EPS) / 2
        hole = {
            "box": [niche.xn * CM, niche.w * CM, (niche.hn + EPS) * CM],
            "position": [niche.xn / 2 * CM, 0.0, z_niche * CM],
            "material": "G4_AIR",
        }
        world = {"box": list(WORLD * CM), "position": [0.0, 0.0, WORLD_Z * CM], "material": "G4_AIR"}
        world |= {"Ground": ground, "Niche": hole}
        frame = mulder.LocalFrame(altitude=altitude, **SITE)
        geometry = mulder.LocalGeometry(calzone.Geometry({"Environment": world}), frame=frame)
    for medium in geometry.media:
        if medium.description.split(".")[-1] in ("Environment", "Niche"):
            medium.material = "Air"
        else:
            medium.material, medium.density = "RegLoam", SOIL
    return geometry


def medium_index(geometry, name):
    return next(i for i, m in enumerate(geometry.media) if m.description.endswith(name))


def wall(bank, depth):
    """Path in the ground (m) from the detector along `legacy_cases.direction_set()`."""
    geometry = build_scene(depth, bank)
    u = legacy_cases.direction_set()
    position = np.broadcast_to(detector(depth), u.shape)
    grammage = np.asarray(geometry.scan(position=position, direction=u, output="grammage"))
    return grammage[:, medium_index(geometry, "Ground")] / SOIL


def cone(axis, alpha, n_r=24, n_phi=48):
    """Directions within `alpha` of `axis`, uniform in cosine and azimuth, and the solid angle of each."""
    e1 = np.cross(axis, [1.0, 0.0, 0.0])
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(axis, e1)
    edges = np.linspace(np.cos(alpha), 1.0, n_r + 1)
    mid, phi = 0.5 * (edges[1:] + edges[:-1]), (np.arange(n_phi) + 0.5) * 2 * np.pi / n_phi
    cb, ph = np.meshgrid(mid, phi, indexing="ij")
    sb = np.sqrt(1 - cb**2)
    u = cb[..., None] * axis + (sb * np.cos(ph))[..., None] * e1 + (sb * np.sin(ph))[..., None] * e2
    return u.reshape(-1, 3), (edges[1] - edges[0]) * 2 * np.pi / n_phi


def reference_void(altitude, factor=1.0):
    """`reference_void` of legacy_cases with Mulder: counts, effect and cone significance.

    The cone is `factor` times the angular radius of the void; 1.0 is the best of
    `legacy.CONE_FACTORS` for the preliminary calculation and for Mulder.
    """
    m = legacy.Model()
    geometry = build_scene(VOID_DEPTH, altitude=altitude, sphere=REFERENCE_VOID)
    index = medium_index(geometry, "Anomaly")
    meter = mulder.Fluxmeter(geometry=geometry)
    d, centre, radius = detector(VOID_DEPTH), np.array(REFERENCE_VOID[:3]), REFERENCE_VOID[3]
    axis = (centre - d) / np.linalg.norm(centre - d)
    u, d_omega = cone(axis, max(np.radians(1.0), factor * np.arcsin(radius / np.linalg.norm(centre - d))))
    scale = 1e-4 * m.area_cm2 * m.t_exp_s * m.eff
    counts = []
    for density in (SOIL, AIR):
        geometry.media[index].density = density
        counts.append(float(np.sum(integral_many(meter, d, u) * u[:, 2]) * d_omega * scale))
    n0, n1 = counts
    sig = abs(n1 - n0) / np.sqrt(n0) / np.sqrt(1 + 1 / m.neigh)
    return {"N0": n0, "N1": n1, "effect_percent": 100 * (n1 / n0 - 1), "sig": sig}


# ------------------------------------------------------------------ registry


def _cases():
    out = {"reference/differential": reference_differential}
    for altitude in (0.0, ALTITUDE):
        for elevation in (90.0, 45.0):
            out[f"open_sky/{altitude:.0f}m/elev{elevation:g}"] = functools.partial(
                open_sky, altitude, elevation
            )
    for altitude in (ALTITUDE, HIGH):
        for elevation in (90.0, 45.0, 30.0):
            out[f"cutoff/{altitude:.0f}m/elev{elevation:g}"] = functools.partial(cutoff, altitude, elevation)
    for name in RANGE_MATERIALS:
        out[f"threshold/{name}"] = functools.partial(threshold, name)
        out[f"range/{name}"] = functools.partial(csda_range, name)
    out["physics/defaults"] = physics_defaults
    for name in ("RegLoam", "RegConcrete", "Air"):
        out[f"physics/{name}"] = functools.partial(physics_tables, name)
    for altitude in (0.0, ALTITUDE):
        for depth in DEPTHS:
            for zenith in ZENITH:
                name = f"slab/{altitude:.0f}m/{depth:.0f}m/zenith{zenith:g}"
                out[name] = functools.partial(slab, altitude, depth, zenith)
    for model in ("us-standard", "midlatitude-summer", "midlatitude-winter"):
        out[f"season/{model}"] = functools.partial(season, model)
    out["geomagnet/1150m/elev45"] = geomagnet
    for bank in legacy_cases.BANKS:
        for depth in (legacy_cases.UPPER, legacy_cases.LOWER):
            out[f"wall/{bank}/{depth:.0f}m"] = functools.partial(wall, bank, depth)
    for altitude in (0.0, ALTITUDE):
        out[f"reference_void/{altitude:.0f}m/cone1"] = functools.partial(reference_void, altitude)
    return out


CASES = _cases()


def compute(name):
    """Result of one case as plain Python numbers."""
    return legacy_cases.plain(CASES[name]())


def versions():
    return {"mulder": mulder.config.VERSION, "calzone": calzone.VERSION}


def parameters():
    """Everything a reader needs to reproduce the cases by hand."""
    return {
        "site": SITE | {"altitudes_m": [0.0, ALTITUDE, HIGH], "date_geomagnet": DATE},
        "muon_mass_gev": MU,
        "materials": {
            "base": {k: {"formula": f, "density_kg_m3": d} for k, (f, d) in BASE.items()},
            "RegLoam": {"mass_fractions": LOAM, "bulk_density_kg_m3": SOIL},
            "RegConcrete": {"mass_fractions": CONCRETE, "density_kg_m3": CONCRETE_DENSITY},
            "Rock, Water, Air": "Mulder defaults",
        },
        "opacities_g_cm2": list(OPACITIES),
        "range_kinetic_energies_gev": list(RANGE_ENERGIES),
        "table_kinetic_energies_gev": list(TABLE_ENERGIES),
        "slab": {"depths_m": list(DEPTHS), "zenith_deg": list(ZENITH), "soil": "RegLoam"},
        "integration": {
            "above a threshold": "trapezoids in ln K, 801 nodes up to 1e5 GeV",
            "all energies": "trapezoids in ln(K - K_cut), 5601 nodes from 1e-9 GeV; 1e-7 GeV if no cut-off",
            "many directions": "rectangles in ln K, 71 nodes from 1e-3 to 1e4 GeV",
        },
        "wall": "niche and directions of legacy_cases; ground 1500 x 3000 x 400 m; floor EPS below it",
        "reference_void": {
            "centre_m": list(REFERENCE_VOID[:3]),
            "radius_m": REFERENCE_VOID[3],
            "niche_depth_m": VOID_DEPTH,
            "void_density_kg_m3": AIR,
            "cone": "1.0 x the angular radius, 24 x 48 directions uniform in cosine and azimuth",
        },
        "eps_m": EPS,
    }
