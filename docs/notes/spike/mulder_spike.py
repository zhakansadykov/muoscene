"""Mulder and Calzone spike (roadmap task 0.1-02): the numbers behind docs/notes/mulder.md.

Compares Mulder 0.3.7 with the preliminary calculation used for the 2027 grant application
(Guan et al. 2015 spectrum at sea level, energy loss a + bE, straight paths, cone significance).
The scene is synthetic: a flat bank behind a vertical wall, a niche, a spherical void, a pier.

Run:
    pip install mulder calzone
    python docs/notes/spike/mulder_spike.py                 all experiments, about 6 minutes
    python docs/notes/spike/mulder_spike.py speed           only the speed test
    python docs/notes/spike/mulder_spike.py --events 40000 scattering
    python docs/notes/spike/mulder_spike.py figure          needs matplotlib
"""

import argparse
import pathlib
import platform
import tempfile
import time

import calzone
import mulder
import mulder.materials as mm
import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
MU = 0.1056583755  # GeV
CM = 100.0  # Calzone lengths are centimetres
SITE = {"latitude": 43.35, "longitude": 79.08}
ALTITUDE = 1150.0  # m
SOIL = 1800.0  # kg/m^3
EPS = 1e-3  # m, the detector sits this much above the niche floor

# preliminary calculation
P1, P2, P3, P4, P5 = 0.102573, -0.068287, 0.958633, 0.0407253, 0.817285
A_LEGACY, B_LEGACY, OFFSET_LEGACY = 2.2e-3, 4.0e-6, 0.1
AREA, EXPOSURE, EFFICIENCY, NEIGHBOURS = 240.0, 60 * 86400.0, 0.9, 8


def define_materials() -> None:
    mm.Material.define("Quartz", composition="SiO2", density=2650.0)
    mm.Material.define("Calcite", composition="CaCO3", density=2710.0)
    mm.Material.define("Alumina", composition="Al2O3", density=3950.0)
    mm.Composite.define("Loam", composition={"Quartz": 0.60, "Calcite": 0.15, "Alumina": 0.12, "Water": 0.13})
    concrete = {"O": 0.529, "Si": 0.337, "Ca": 0.044, "Al": 0.034, "Na": 0.016, "Fe": 0.014}
    concrete |= {"K": 0.013, "H": 0.010, "Mg": 0.002, "C": 0.001}
    mm.Material.define("Concrete", composition=concrete, density=2300.0)


def guan(energy, cos_theta):
    """Guan et al. (2015), m^-2 s^-1 sr^-1 GeV^-1; total energy in GeV."""
    c = np.sqrt((cos_theta**2 + P1**2 + P2 * cos_theta**P3 + P4 * cos_theta**P5) / (1 + P1**2 + P2 + P4))
    e = energy * (1 + 3.64 / (energy * c**1.29))
    return 1e4 * 0.14 * e**-2.7 * (1 / (1 + 1.1 * energy * c / 115) + 0.054 / (1 + 1.1 * energy * c / 850))


def guan_integral(e_min, cos_theta, e_max=1e5, n=2001):
    e = np.logspace(np.log10(e_min), np.log10(e_max), n)
    return np.trapezoid(guan(e, cos_theta) * e, np.log(e))


def legacy_threshold(opacity):
    """Total energy needed to cross `opacity` (g/cm^2) in the preliminary calculation."""
    return A_LEGACY / B_LEGACY * np.expm1(B_LEGACY * opacity) + OFFSET_LEGACY


def legacy_path(u, depth, rise=0.0, x_rise=1.0, x_det=0.25, niche=(0.5, 0.3, 0.3)):
    """Path in the ground (m) of the preliminary calculation: wall x = 0, ground x > 0, brow rising
    linearly by `rise` over `x_rise`, niche of (depth, width, height), detector on its floor."""
    xn, w, hn = niche
    ax, ay, c = u[:, 0], u[:, 1], u[:, 2]
    big = 1e9
    t_roof = hn / c
    t_side = np.where(np.abs(ay) > 1e-12, (w / 2) / np.maximum(np.abs(ay), 1e-12), big)
    t_back = np.where(ax > 1e-12, (xn - x_det) / np.maximum(ax, 1e-12), big)
    t_front = np.where(ax < -1e-12, x_det / np.maximum(-ax, 1e-12), big)
    t_exit = np.minimum.reduce([t_roof, t_side, t_back, t_front])
    if rise > 0:
        den = c - rise * ax / x_rise
        t1 = np.where(den > 1e-9, (depth + rise * x_det / x_rise) / np.maximum(den, 1e-9), big)
        t_surf = np.where((t1 < big) & (x_det + t1 * ax <= x_rise), t1, (depth + rise) / c)
        t_surf = np.where(ax <= 0, depth / c, t_surf)
    else:
        t_surf = depth / c
    t_end = np.where(ax < 0, np.minimum(t_surf, t_front), t_surf)
    return np.where(t_front <= t_exit, 0.0, np.maximum(t_end - t_exit, 0.0))


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


def build_scene(depth, altitude=ALTITUDE, sphere=None, pier=None, brow_mesh=None):
    """Local geometry: ground x in [0, 400] m behind a vertical wall, niche at `depth` below the brow.

    sphere = (x, y, z, radius); pier = (x0, cap size, cap height, pile radius, pile length, n, pitch);
    brow_mesh = path of an OBJ mesh replacing the ground box (then no sphere or pier).
    """
    centre = np.array([200.0, 0.0, -100.0])

    def rel(x, y, z):
        return [float(q) * CM for q in np.array([x, y, z]) - centre]

    if brow_mesh is None:
        ground = {"box": [400 * CM, 800 * CM, 200 * CM], "position": list(centre * CM)}
    else:
        ground = {"mesh": {"path": str(brow_mesh), "units": "m"}}
    ground |= {"material": "G4_SILICON_DIOXIDE", "subtract": ["Niche"]}
    if sphere is not None:
        x, y, z, r = sphere
        ground["Anomaly"] = {"sphere": {"radius": r * CM}, "position": rel(x, y, z), "material": "G4_AIR"}
    if pier is not None:
        x0, size, cap, radius, length, n, pitch = pier
        ground["Cap"] = {"box": [size * CM, size * CM, cap * CM], "material": "G4_CONCRETE"}
        ground["Cap"]["position"] = rel(x0 + size / 2, 0, -cap / 2)
        for i in range(n):
            for j in range(n):
                px = x0 + size / 2 + (i - (n - 1) / 2) * pitch
                py = (j - (n - 1) / 2) * pitch
                pile = {"cylinder": {"length": length * CM, "radius": radius * CM}, "material": "G4_CONCRETE"}
                pile["position"] = rel(px, py, -cap - length / 2)
                ground[f"Pile{i * n + j}"] = pile
    niche = {"box": [0.5 * CM, 0.3 * CM, 0.3 * CM], "position": [0.25 * CM, 0, (-depth + 0.15) * CM]}
    niche["material"] = "G4_AIR"
    tree = {"Environment": {"box": [1000 * CM, 1000 * CM, 600 * CM], "material": "G4_AIR"}}
    tree["Environment"] |= {"Ground": ground, "Niche": niche}
    frame = mulder.LocalFrame(altitude=altitude, **SITE)
    geometry = mulder.LocalGeometry(calzone.Geometry(tree), frame=frame)
    for medium in geometry.media:
        name = medium.description.split(".")[-1]
        if name in ("Environment", "Niche"):
            medium.material = "Air"
        elif name in ("Ground", "Anomaly"):
            medium.material, medium.density = "Loam", SOIL
        else:
            medium.material, medium.density = "Concrete", 2400.0
    return geometry


def medium_index(geometry, name):
    return [i for i, m in enumerate(geometry.media) if m.description.endswith(name)][0]


def detector(depth):
    return np.array([0.25, 0.0, -depth + EPS])


def slope_grid(n, t_max=1.0):
    tx, ty = np.meshgrid(np.linspace(-t_max, t_max, n), np.linspace(-t_max, t_max, n))
    norm = np.sqrt(1 + tx**2 + ty**2)
    return np.stack([tx / norm, ty / norm, 1 / norm], axis=-1).reshape(-1, 3)


K_GRID = np.logspace(-3, 4, 71)  # kinetic energy at the detector, GeV
K_WEIGHT = K_GRID * np.gradient(np.log(K_GRID))


def integral_flux(meter, k_min=1e-3, k_max=1e5, n=801, **state):
    """Integral over kinetic energy at the detector of the Mulder differential flux, m^-2 s^-1 sr^-1."""
    k = np.logspace(np.log10(k_min), np.log10(k_max), n)
    return np.trapezoid(np.asarray(meter.flux(energy=k, **state)) * k, np.log(k))


def integral_flux_many(meter, position, directions, k=K_GRID):
    """The same for many directions at once, on a coarser energy grid."""
    n, m = len(directions), len(k)
    flux = meter.flux(
        position=np.broadcast_to(position, (n * m, 3)),
        direction=np.repeat(directions, m, axis=0),
        energy=np.tile(k, n),
    )
    return (np.asarray(flux).reshape(n, m) * k * np.gradient(np.log(k))).sum(axis=1)


def exp_reference(args):
    print("\n## E2 reference flux (GCCLY15 against the Guan formula, total energy = kinetic + mass)")
    reference = mulder.Reference()
    worst = 0.0
    for elevation in (90.0, 60.0, 30.0):
        for k in (0.5, 1.0, 10.0, 100.0):
            ours = guan(k + MU, np.sin(np.radians(elevation)))
            theirs = float(reference.flux(energy=k, elevation=elevation, altitude=0.0))
            worst = max(worst, abs(theirs / ours - 1))
    print(f"largest relative difference of the differential flux: {worst:.1e}")
    meter = mulder.Fluxmeter()
    state = dict(SITE, azimuth=0.0, elevation=90.0)
    print(
        f"vertical, sea level: I(E > 1 GeV) = {integral_flux(meter, k_min=1 - MU, altitude=0.0, **state):.2f}"
        f" (formula {guan_integral(1.0, 1.0):.2f}),"
        f" all energies {integral_flux(meter, altitude=0.0, **state):.2f}"
        f" (formula {guan_integral(MU * 1.0001, 1.0):.2f}) m^-2 s^-1 sr^-1"
    )
    print("altitude  elevation  all energies  E>1 GeV  E>5 GeV  E>10 GeV")
    for altitude in (0.0, ALTITUDE):
        for elevation in (90.0, 45.0):
            s = dict(SITE, azimuth=0.0, elevation=elevation, altitude=altitude)
            row = [integral_flux(meter, k_min=e - MU if e else 1e-3, **s) for e in (0, 1.0, 5.0, 10.0)]
            print(f"{altitude:8.0f}  {elevation:9.0f}  " + "  ".join(f"{v:8.3f}" for v in row))


def exp_range(args):
    print("\n## E3 energy threshold: total energy to cross an opacity X, GeV")
    names = ("Rock", "Loam", "Concrete", "Water")
    compiled = mulder.Physics().compile(*names)
    density = {"Rock": 2650.0, "Loam": mm.Composite("Loam").density, "Concrete": 2300.0}
    density["Water"] = mm.Material("Water").density
    print("X, g/cm2   a+bE   " + "  ".join(f"{n:>8s}" for n in names) + "   a+bE over Loam")
    for x in (500.0, 1000.0, 1800.0, 3600.0, 5400.0, 10000.0):
        row = [float(c.inverse_range(x * 10 / density[c.name])) + MU for c in compiled]
        legacy = legacy_threshold(x)
        print(
            f"{x:8.0f} {legacy:7.3f}  "
            + "  ".join(f"{v:8.3f}" for v in row)
            + f"   {100 * (legacy / row[1] - 1):+5.1f} %"
        )


def exp_slab(args):
    print("\n## E4 uniform soil slab at sea level: integral flux, m^-2 s^-1 sr^-1")
    layer = mulder.Layer(0.0, material="Loam", density=SOIL)
    continuous = mulder.Fluxmeter(layer)
    print("depth  zenith  X, g/cm2   Mulder   a+bE and Guan   ratio")
    for depth in (10.0, 20.0, 30.0):
        for zenith in (0.0, 30.0, 45.0, 60.0):
            c = np.cos(np.radians(zenith))
            x = SOIL / 10 * depth / c
            s = dict(SITE, altitude=-depth, azimuth=0.0, elevation=90.0 - zenith)
            ours = integral_flux(continuous, **s)
            legacy = guan_integral(legacy_threshold(x), c)
            print(f"{depth:5.0f}  {zenith:6.0f}  {x:8.0f}  {ours:7.3f}  {legacy:14.3f}  {ours / legacy:6.3f}")
    k = np.logspace(-1, 3, 25)
    w = k * np.gradient(np.log(k))
    events = max(200, args.events // 20)
    print(f"continuous against mixed and discrete, {events} events per energy:")
    for depth, zenith in ((20.0, 0.0), (20.0, 45.0)):
        s = dict(SITE, altitude=-depth, azimuth=0.0, elevation=90.0 - zenith)
        base = np.sum(np.asarray(continuous.flux(energy=k, **s)) * w)
        line = f"depth {depth:.0f} m, zenith {zenith:.0f}: continuous {base:.3f}"
        for mode in ("mixed", "discrete"):
            r = np.asarray(mulder.Fluxmeter(layer, mode=mode, seed=1).flux(energy=k, events=events, **s))
            value, error = np.sum(r[:, 0] * w), np.sqrt(np.sum((r[:, 1] * w) ** 2))
            line += f"; {mode} {100 * (value / base - 1):+.2f} ± {100 * error / base:.2f} %"
        print(line)


def exp_wall(args):
    print("\n## E5 wall with a niche: path in the ground against the preliminary calculation")
    u = slope_grid(201)
    mesh = pathlib.Path(tempfile.mkdtemp()) / "brow.obj"
    rise, x_rise = 2.0, 9.0
    write_prism(mesh, [(0, -200), (400, -200), (400, rise), (x_rise, rise), (0, 0)], -400, 400)
    for label, brow, kwargs in (
        ("flat brow, boxes", None, {}),
        ("rising brow, mesh", mesh, {"rise": rise, "x_rise": x_rise}),
    ):
        for depth in (10.0, 25.0):
            geometry = build_scene(depth, brow_mesh=brow)
            start = time.perf_counter()
            position = np.broadcast_to(detector(depth), u.shape)
            grammage = np.asarray(geometry.scan(position=position, direction=u, output="grammage"))
            rate = len(u) / (time.perf_counter() - start)
            ours = grammage[:, medium_index(geometry, "Ground")] / SOIL
            legacy = legacy_path(u, depth - EPS, niche=(0.5, 0.3, 0.3 - EPS), **kwargs)
            diff = np.abs(ours - legacy)
            inward = u[:, 0] > 0
            print(
                f"{label}, niche {depth:.0f} m: {rate:.1e} directions/s; max |dL| = {diff.max():.1e} m;"
                f" towards the massif {diff[inward].max():.1e} m;"
                f" differing by > 1 um: {np.sum(diff > 1e-6)} of {len(u)}"
            )


def cone(axis, alpha, n_r=24, n_phi=48):
    t = np.array([1.0, 0.0, 0.0])
    e1 = np.cross(axis, t)
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(axis, e1)
    edges = np.linspace(np.cos(alpha), 1.0, n_r + 1)
    cb, ph = np.meshgrid(
        0.5 * (edges[1:] + edges[:-1]), (np.arange(n_phi) + 0.5) * 2 * np.pi / n_phi, indexing="ij"
    )
    sb = np.sqrt(1 - cb**2)
    u = cb[..., None] * axis + (sb * np.cos(ph))[..., None] * e1 + (sb * np.sin(ph))[..., None] * e2
    return u.reshape(-1, 3), (edges[1] - edges[0]) * 2 * np.pi / n_phi


VOID = (3.0, 0.0, -3.0, 0.5)  # a 1 m void 3 m from the wall and 3 m below the brow
DEPTH = 10.0


def exp_significance(args):
    print("\n## E6 a 1 m void seen from a niche 10 m below the brow: counts and significance")
    print(
        f"horizontal films {AREA:.0f} cm2, {EXPOSURE / 86400:.0f} days, efficiency {EFFICIENCY},"
        " cone estimator"
    )
    centre, radius = np.array(VOID[:3]), VOID[3]
    d = detector(DEPTH)
    axis = (centre - d) / np.linalg.norm(centre - d)
    angular = np.arcsin(radius / np.linalg.norm(centre - d))
    scenes = {}
    for label, altitude in (("Mulder, sea level", 0.0), (f"Mulder, {ALTITUDE:.0f} m", ALTITUDE)):
        geometry = build_scene(DEPTH, altitude=altitude, sphere=VOID)
        scenes[label] = (geometry, medium_index(geometry, "Anomaly"), mulder.Fluxmeter(geometry=geometry))
    scale = 1e-4 * AREA * EXPOSURE * EFFICIENCY

    def sigma(n0, n1):
        return abs(n1 - n0) / np.sqrt(n0) / np.sqrt(1 + 1 / NEIGHBOURS)

    def chord(u):
        v = centre - d
        tca = u @ v
        d2 = v @ v - tca**2
        return np.where((tca > 0) & (d2 < radius**2), 2 * np.sqrt(np.maximum(radius**2 - d2, 0)), 0.0)

    best = {}
    print("cone   preliminary: N0, effect, sigma   " + "   ".join(f"{k}: N0, effect, sigma" for k in scenes))
    for k in (0.6, 0.8, 1.0, 1.3, 1.7):
        u, weight = cone(axis, max(np.radians(1.0), k * angular))
        path = legacy_path(u, DEPTH - EPS, niche=(0.5, 0.3, 0.3 - EPS))
        counts = []
        for void in (False, True):
            x = SOIL / 10 * (path - chord(u) * void)
            flux = np.array(
                [guan_integral(legacy_threshold(xi), ci, n=801) for xi, ci in zip(x, u[:, 2], strict=True)]
            )
            counts.append(np.sum(flux * u[:, 2]) * weight * scale)
        line = (
            f"{k:4.1f}   {counts[0]:8.0f} {100 * (counts[1] / counts[0] - 1):+6.2f} % {sigma(*counts):5.2f}"
        )
        best["preliminary"] = max(best.get("preliminary", 0), sigma(*counts))
        for label, (geometry, index, meter) in scenes.items():
            counts = []
            for density in (SOIL, 1.2):
                geometry.media[index].density = density
                counts.append(np.sum(integral_flux_many(meter, d, u) * u[:, 2]) * weight * scale)
            line += f"   {counts[0]:8.0f} {100 * (counts[1] / counts[0] - 1):+6.2f} % {sigma(*counts):5.2f}"
            best[label] = max(best.get(label, 0), sigma(*counts))
        print(line)
    print("best significance: " + "; ".join(f"{k} {v:.2f}" for k, v in best.items()))


def exp_scattering(args):
    print(
        "\n## E6 scattering: excess of counts inside a cone around the 1 m void, continuous against discrete"
    )
    centre, radius = np.array(VOID[:3]), VOID[3]
    d = detector(DEPTH)
    axis = (centre - d) / np.linalg.norm(centre - d)
    angular = np.arcsin(radius / np.linalg.norm(centre - d))
    geometry = build_scene(DEPTH, sphere=VOID)
    index = medium_index(geometry, "Anomaly")
    n = args.events
    t = np.array([1.0, 0.0, 0.0])
    e1 = np.cross(axis, t)
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(axis, e1)
    for k in (1.0, 1.7):
        rng = np.random.default_rng(11)
        energy = np.exp(rng.uniform(np.log(1e-1), np.log(1e3), n))
        cb = rng.uniform(np.cos(k * angular), 1.0, n)
        sb, ph = np.sqrt(1 - cb**2), rng.uniform(0, 2 * np.pi, n)
        u = cb[:, None] * axis + (sb * np.cos(ph))[:, None] * e1 + (sb * np.sin(ph))[:, None] * e2
        position = np.broadcast_to(d, u.shape)
        line = f"cone {k:.1f} x the angular radius ({np.degrees(k * angular):.1f} deg), {n} states:"
        for mode in ("continuous", "discrete"):
            values = []
            for density in (SOIL, 1.2):
                geometry.media[index].density = density
                meter = mulder.Fluxmeter(geometry=geometry, mode=mode, seed=5)
                start = time.perf_counter()
                extra = {"events": 1} if mode == "discrete" else {}
                flux = np.asarray(meter.flux(position=position, direction=u, energy=energy, **extra))
                rate = n / (time.perf_counter() - start)
                values.append((flux[:, 0] if flux.ndim == 2 else flux) * energy * u[:, 2])
            base, void = values
            ratio = void.sum() / base.sum() - 1
            error = np.sqrt(n) * np.std((void - base) - ratio * base) / base.sum()
            line += f" {mode} {100 * ratio:+.2f} ± {100 * error:.2f} % ({rate:.0f} states/s);"
        print(line)


def exp_season(args):
    print(f"\n## E9 season at {ALTITUDE:.0f} m: open sky against 20 m of soil, vertical, m^-2 s^-1 sr^-1")
    layer = mulder.Layer(ALTITUDE, material="Loam", density=SOIL)
    rows = {}
    for model in ("us-standard", "midlatitude-summer", "midlatitude-winter"):
        meter = mulder.Fluxmeter(layer, atmosphere=model)
        s = dict(SITE, azimuth=0.0, elevation=90.0)
        rows[model] = (
            integral_flux(meter, altitude=ALTITUDE + 0.01, **s),
            integral_flux(meter, altitude=ALTITUDE - 20, **s),
        )
        print(
            f"{model:20s} open sky {rows[model][0]:8.3f}   under 20 m {rows[model][1]:7.4f}"
            f"   ratio {rows[model][1] / rows[model][0]:.5f}"
        )
    summer, winter = rows["midlatitude-summer"], rows["midlatitude-winter"]
    print(
        f"winter against summer: open sky {100 * (winter[0] / summer[0] - 1):+.2f} %, under 20 m"
        f" {100 * (winter[1] / summer[1] - 1):+.2f} %,"
        f" ratio {100 * (winter[1] / winter[0] / (summer[1] / summer[0]) - 1):+.2f} %"
    )
    meter = mulder.Fluxmeter(layer, date="2028-10-15")
    east, west = (
        integral_flux(meter, altitude=ALTITUDE + 0.01, azimuth=a, elevation=45.0, **SITE)
        for a in (90.0, 270.0)
    )
    print(
        "geomagnetic field (IGRF, 2028-10-15), open sky at 45 deg elevation:"
        f" east {east:.2f}, west {west:.2f}"
        f" ({100 * (east / west - 1):+.1f} %)"
    )


def exp_speed(args):
    print("\n## E7 speed")
    print(
        f"{platform.platform()}, Python {platform.python_version()},"
        f" mulder {mulder.config.VERSION}, calzone {calzone.VERSION}"
    )
    pier = (10.0, 7.0, 2.0, 0.5, 10.0, 3, 2.5)
    start = time.perf_counter()
    geometry = build_scene(25.0, sphere=(8.0, 0.0, -6.0, 1.0), pier=pier)
    print(
        f"scene with a pier: {len(geometry.media)} media,"
        f" built in {1e3 * (time.perf_counter() - start):.0f} ms"
    )
    d = detector(25.0)
    u = slope_grid(401, 1.2)
    position = np.broadcast_to(d, u.shape)
    geometry.scan(position=position[:100], direction=u[:100], output="grammage")  # warm-up
    best = float("inf")
    for _ in range(3):
        start = time.perf_counter()
        geometry.scan(position=position, direction=u, output="grammage")
        best = min(best, time.perf_counter() - start)
    print(f"scan (opacity per medium): {len(u) / best:.2e} directions/s (best of 3)")
    meter = mulder.Fluxmeter(geometry=geometry)
    integral_flux_many(meter, d, u[:10])
    start = time.perf_counter()
    integral_flux_many(meter, d, u[::16])
    elapsed = time.perf_counter() - start
    n = len(u[::16])
    print(
        f"continuous flux: {n * len(K_GRID) / elapsed:.2e} states/s"
        f" = {n / elapsed:.2e} directions/s with {len(K_GRID)} energies"
    )
    meter = mulder.Fluxmeter(geometry=geometry, mode="discrete", seed=1)
    n = 2000
    start = time.perf_counter()
    meter.flux(
        position=np.broadcast_to(d, (n, 3)),
        direction=u[:: len(u) // n][:n],
        energy=np.full(n, 10.0),
        events=1,
    )
    print(f"discrete flux: {n / (time.perf_counter() - start):.0f} states/s")


def exp_cutoff(args):
    print("\n## E13 altitude cut-off: with the sea-level reference the flux at altitude is zero below K_cut")
    meter = mulder.Fluxmeter()
    k = np.logspace(-2, 0.5, 2001)
    for altitude in (ALTITUDE, 3340.0):
        for elevation in (90.0, 45.0, 30.0):
            state = dict(SITE, altitude=altitude, azimuth=0.0, elevation=elevation)
            flux = np.asarray(meter.flux(energy=k, **state))
            print(
                f"altitude {altitude:.0f} m, elevation {elevation:.0f}:"
                f" K_cut = {k[np.argmax(flux > 0)]:.2f} GeV"
            )


def exp_hybrid(args):
    print(
        "\n## E12 fast path (opacity from scan, flux from a table built by Mulder) against direct Mulder flux"
    )
    k = np.logspace(0.0, 4.0, 61)  # kinetic energy above 1 GeV, above the altitude cut-off
    depth = 25.0
    geometry = build_scene(depth, sphere=(8.0, 0.0, -6.0, 1.0), pier=(10.0, 7.0, 2.0, 0.5, 10.0, 3, 2.5))
    geometry.media[medium_index(geometry, "Anomaly")].density = 1.2
    d = detector(depth)
    altitude = ALTITUDE + d[2]

    names = ("Loam", "Concrete")
    compiled = dict(zip(names, mulder.Physics().compile(*names), strict=True))
    default = {"Loam": mm.Composite("Loam").density, "Concrete": mm.Material("Concrete").density}
    k_ref = 7.0
    grammage = {name: float(compiled[name].range(k_ref)) * default[name] for name in names}
    weight = {name: grammage["Loam"] / grammage[name] for name in names}
    print(
        f"equivalent-opacity weight of concrete relative to loam at {k_ref:.0f} GeV: {weight['Concrete']:.4f}"
    )

    dense = 1.0e5  # kg/m^3: 2e4 g/cm^2 is 2 m thick, so the table stays at the detector altitude
    meter = mulder.Fluxmeter(mulder.Layer(altitude, material="Loam", density=dense))
    sines = np.linspace(0.30, 1.0, 36)
    opacity = np.concatenate([[0.0], np.logspace(0.5, np.log10(3.0e4), 70)])  # g/cm^2
    s, x = np.meshgrid(sines, opacity, indexing="ij")
    n, m = s.size, len(k)
    start = time.perf_counter()
    flux = meter.flux(
        latitude=np.full(n * m, SITE["latitude"]),
        longitude=np.full(n * m, SITE["longitude"]),
        altitude=np.repeat((altitude - x * 10 / dense * s).ravel() - 1e-6, m),
        azimuth=np.zeros(n * m),
        elevation=np.repeat(np.degrees(np.arcsin(s)).ravel(), m),
        energy=np.tile(k, n),
    )
    table = np.log((np.asarray(flux).reshape(n, m) * k * np.gradient(np.log(k))).sum(axis=1).reshape(s.shape))
    print(
        f"table of {len(sines)} elevations x {len(opacity)} opacities"
        f" built in {time.perf_counter() - start:.1f} s"
    )
    log_x = np.log(opacity[1:])

    def lookup(sin_e, x_eq):
        """Bilinear in (sin elevation, ln opacity); linear in opacity below the first non-zero node."""
        i = np.clip(np.searchsorted(sines, sin_e) - 1, 0, len(sines) - 2)
        a = (sin_e - sines[i]) / (sines[i + 1] - sines[i])
        lx = np.log(np.maximum(x_eq, opacity[1]))
        j = np.clip(np.searchsorted(log_x, lx) - 1, 0, len(log_x) - 2)
        b = (lx - log_x[j]) / (log_x[j + 1] - log_x[j])
        inner = (1 - a) * ((1 - b) * table[i, j + 1] + b * table[i, j + 2])
        inner += a * ((1 - b) * table[i + 1, j + 1] + b * table[i + 1, j + 2])
        c = x_eq / opacity[1]
        edge = (1 - c) * ((1 - a) * table[i, 0] + a * table[i + 1, 0]) + c * (
            (1 - a) * table[i, 1] + a * table[i + 1, 1]
        )
        return np.exp(np.where(x_eq < opacity[1], edge, inner))

    rng = np.random.default_rng(3)
    tx, ty = rng.uniform(-1.2, 1.2, (2, 3000))
    norm = np.sqrt(1 + tx**2 + ty**2)
    u = np.stack([tx / norm, ty / norm, 1 / norm], axis=1)
    u = u[u[:, 2] > 0.32]
    thickness = np.asarray(geometry.scan(position=np.broadcast_to(d, u.shape), direction=u))  # m per medium
    material = [medium.material for medium in geometry.media]
    density = np.array([medium.density or mm.Material(medium.material).density for medium in geometry.media])
    solid = np.array([name != "Air" for name in material])
    w = np.array([weight.get(name, 1.0) for name in material])
    x_eq = (thickness * density * w)[:, solid].sum(axis=1) / 10
    start = time.perf_counter()
    fast = lookup(u[:, 2], x_eq)
    t_fast = time.perf_counter() - start
    start = time.perf_counter()
    exact = integral_flux_many(mulder.Fluxmeter(geometry=geometry), d, u, k)
    t_exact = time.perf_counter() - start
    concrete = thickness[:, [name == "Concrete" for name in material]].sum(axis=1) > 0.5
    sky = x_eq < 1.0
    for label, selection in (
        ("all", np.full(len(u), True)),
        ("open sky", sky),
        ("soil only", ~sky & ~concrete),
        ("through concrete", concrete),
    ):
        r = fast[selection] / exact[selection] - 1
        print(
            f"{label:17s} {selection.sum():5d} directions: mean {100 * r.mean():+.3f} %,"
            f" rms {100 * r.std():.3f} %,"
            f" max {100 * np.abs(r).max():.3f} %"
        )
    print(
        f"fast path {len(u) / t_fast:.1e} directions/s, direct {len(u) / t_exact:.0f} directions/s,"
        f" speed-up {t_exact / t_fast:.0f}"
    )


def exp_figure(args):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    geometry = build_scene(25.0, sphere=(8.0, 0.0, -6.0, 1.0), pier=(10.0, 7.0, 2.0, 0.5, 10.0, 3, 2.5))
    n, t_max = 401, 1.2
    u = slope_grid(n, t_max)
    grammage = np.asarray(
        geometry.scan(position=np.broadcast_to(detector(25.0), u.shape), direction=u, output="grammage")
    )
    concrete = [i for i, m in enumerate(geometry.media) if m.material == "Concrete"]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8), dpi=110)
    panels = (
        ("Opacity from a niche 25 m below the brow, g/cm²", grammage.sum(axis=1), "viridis"),
        ("Concrete only: pile cap and nine piles, g/cm²", grammage[:, concrete].sum(axis=1), "magma"),
    )
    for ax, (title, values, cmap) in zip(axes, panels, strict=True):
        image = ax.imshow(
            values.reshape(n, n) / 10, origin="lower", extent=[-t_max, t_max, -t_max, t_max], cmap=cmap
        )
        fig.colorbar(image, ax=ax)
        ax.set(title=title, xlabel="tx (into the massif →)", ylabel="ty (along the wall)")
    fig.tight_layout()
    fig.savefig(HERE / "pier_opacity.png")
    print("written", HERE / "pier_opacity.png")

    layers = {0.0: mulder.Layer(0.0, material="Loam", density=SOIL)}
    layers[ALTITUDE] = mulder.Layer(ALTITUDE, material="Loam", density=SOIL)
    depths = np.array([2.0, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 60.0])
    fig, ax = plt.subplots(figsize=(6.4, 4.2), dpi=110)
    for altitude, layer in layers.items():
        meter = mulder.Fluxmeter(layer)
        for zenith, style in ((0.0, "-"), (45.0, "--")):
            c = np.cos(np.radians(zenith))
            state = dict(SITE, azimuth=0.0, elevation=90.0 - zenith)
            ours = [integral_flux(meter, altitude=altitude - d, **state) for d in depths]
            legacy = [guan_integral(legacy_threshold(SOIL / 10 * d / c), c) for d in depths]
            label = f"{altitude:.0f} m a.s.l., zenith {zenith:.0f}°"
            ax.plot(SOIL / 10 * depths / c, np.array(ours) / np.array(legacy), style, label=label)
    ax.axhline(1.0, color="0.6", lw=0.8)
    ax.set(xscale="log", xlabel="opacity along the line of sight, g/cm²")
    ax.set(ylabel="Mulder / preliminary calculation", title="Integral flux under a uniform soil slab")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(HERE / "slab_ratio.png")
    print("written", HERE / "slab_ratio.png")


EXPERIMENTS = {
    "reference": exp_reference,
    "range": exp_range,
    "slab": exp_slab,
    "wall": exp_wall,
    "significance": exp_significance,
    "scattering": exp_scattering,
    "season": exp_season,
    "speed": exp_speed,
    "hybrid": exp_hybrid,
    "cutoff": exp_cutoff,
    "figure": exp_figure,
}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("experiments", nargs="*", help="any of: " + ", ".join(EXPERIMENTS))
    parser.add_argument("--events", type=int, default=20000, help="Monte Carlo states for the discrete mode")
    args = parser.parse_args()
    unknown = set(args.experiments) - set(EXPERIMENTS)
    if unknown:
        parser.error(f"unknown experiments: {sorted(unknown)}")
    define_materials()
    for name in args.experiments or [n for n in EXPERIMENTS if n != "figure"]:
        EXPERIMENTS[name](args)


if __name__ == "__main__":
    main()
