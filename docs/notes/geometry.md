# Geometry

Roadmap task 0.1-06, 17 November 2026, Mulder 0.3.7, Calzone 1.2.2; the section on lines of sight through an edge of a mesh: roadmap task 0.1-10b, 27 November 2026. Module `muoscene.geometry`; conventions, revision 4, sections "Frames" and "Bodies, densities and anomalies". Every number below comes from the named test of `tests/unit/test_geometry_scene.py`, `test_geometry_pathtable.py` or `test_geometry_checks.py`, from a frozen case of `tests/regression/data/mulder_reference.json`, or from the command printed next to it; the headline ones are in `docs/numbers.md`. Sizes in tests are those the repository already had: the wall, banks, niches and pile cap of `tests/regression/legacy_cases.py`, the reference void of `mulder_cases.py`, the invented pier of `docs/notes/spike/mulder_spike.py`.

## What a scene is

A list of bodies in the world frame — local ENU, metres — and the altitude of its origin.

| Body | Shape | What it is |
|---|---|---|
| `Ground` | a profile of the surface in the x–z plane, extruded along y between two limits, with a flat bottom | everything under the profile; `Layer`s divide it by elevation, each from its top down to the next |
| `Structure` | a `Box`, a vertical `Cylinder`, or several (`pile_grid`) | one material, one column of the path table |
| `Cut` | a `Box` | a niche or an excavation: air, with precedence over everything else |
| `Anomaly` | a `Sphere` | scales the density of the body it lies in by (1 + `drho_rel`); −1 is a void |

A body takes a `Material`, a `Mixture` or a `Soil` of `muoscene.materials`. A soil gives its mixture of skeleton and water to Mulder and its bulk density ρ_d (1 + w) to the geometry, and takes no other density; a material or a mixture defaults to its own density. The geometry knows a material only by its name in Mulder (`materials.register`).

The profile is a list of points (x, z) with x not decreasing; two points with one x are a vertical wall. Overhangs cannot be given — a niche is a cut. The bottom of a ground, its two faces across y and the vertical faces under the first and the last point of its profile are limits of the model, not surfaces of the object, so a canyon wall is a segment of the profile: the bank of the frozen cases is `[(0, -400), (0, 0), (1500, 0)]` over a bottom at −400 m.

Version 0.1 has no rotations: boxes follow the axes, cylinders are vertical, the ground runs along y. A scene whose wall does not run north–south is described in its own axes, and the frame is turned in roadmap 0.1-07 (decided on 17 November 2026). Since that task `Scene` takes `azimuth_deg`, the azimuth of its +y axis from true north: it changes no scan and only tells Mulder where north is, for the geomagnetic field; points and directions in ENU go through `muoscene.site.Site.scan` (`docs/notes/detectors.md`).

## How it is built

`Scene` has `builder.assemble` make a Calzone geometry, wraps it as a `mulder.LocalGeometry` at the altitude of the site and sets the material and the density of every medium.

- **Ground and layers.** Every connected piece of ground between two layer tops is one closed mesh; the pieces stand side by side in the world and touch along the layer tops. Two banks of a trench give two pieces of each layer above its floor. The point where a layer top cuts a slope is the same number in both pieces (`test_pieces_of_a_ground_share_their_vertices`). Layers are not nested, because Calzone subtracts only volumes of one mother, and a nested layer could not subtract the niche.
- **Cuts** are boxes of air in the world, subtracted by every piece of ground they reach.
- **Structures.** A shape lies in one ground or in the air. In the ground it is split at the tops of the layers it crosses, and each part is a daughter of the piece of ground around it; in the air it is a daughter of the world.
- **Anomalies** are daughters of the one volume they lie in, their host: a part of a structure, or a piece of a ground or of a layer.
- **Units and frames.** Metres become Calzone centimetres. A daughter is placed relative to the origin of its mother: the centre of a box or a cylinder, the origin of the scene for a mesh and for the world, whose own position only moves its box.
- **Meshes.** Calzone reads a mesh only from a file: prisms are written as OBJ files with one normal per face into a temporary directory of the system, removed as soon as the geometry is built. Calzone names are letters and digits (`V0`, `V1`, …); the names of bodies live in muoscene.
- **Air.** The box of the scene and the cuts are the air of the library at the density of the atmosphere model at the altitude of the site; the box is `margin_m` (1 m) larger than the bodies on every side. See below.

Where bodies meet, the conventions give precedence to cuts, then structures, then layers, then the ground. What version 0.1 does with each case:

| Case | What happens |
|---|---|
| a cut in a ground or in layers, also across a layer top or out of the wall | subtracted |
| a structure across layer tops | split |
| structures, cuts or anomalies that touch | accepted; the overlap check of Geant4 (`Scene.check`) passes |
| structures that overlap, cuts that overlap, grounds that overlap in plan | refused |
| a cut through a structure | refused: shape the structure around the cut |
| a structure across the surface or a limit of a ground | refused: give the part in the ground and the part in the air as two structures |
| an anomaly across a surface, a layer top, a structure, a cut or another anomaly | refused, not split |

All refusals are `ValueError`s raised when the scene is built, with the names of the bodies (`tests/unit/test_geometry_checks.py`).

## Scan and path table

`Scene.scan(position_m, direction)` gives a `PathTable`: for every line of sight the length of the path, in metres, in every body. Directions are unit vectors from the detector towards the sky; points and directions are broadcast, and the table keeps their shape. The columns are the bodies in the order given, the layers of a ground after it and the air of the scene box last; the parts of a split structure and the pieces of a layer are summed into one column.

Opacity is a sum over that table (conventions): `opacity_g_cm2(densities_g_cm3, drho_rel)` adds density times length over grounds, layers, structures and anomalies, with other bulk densities for named bodies and other relative changes for named anomalies. An anomaly has the density of its host under the same hypothesis times (1 + `drho_rel`). The air of the box and of the cuts is part of the atmosphere and is left out unless asked for. `opacity_by_material_g_cm2` gives the same per material of Mulder, which is what the equivalent opacity of the conventions is summed over (roadmap 0.1-08). An unknown name, a density given to an anomaly or to the air, a `drho_rel` below −1 are refused.

| Check | Result | Test |
|---|---|---|
| path in the ground behind the wall with a niche against the frozen `wall/…` cases | within 1.1·10⁻¹³ m; the criterion is 10⁻⁹ m | `test_wall_paths_equal_the_frozen_ones` |
| the same against the preliminary calculation, flat brow | within 10⁻⁹ m | `test_wall_paths_equal_the_preliminary_calculation` |
| one to three layers against the whole ground, also with tops through the niches | sum within 1.1·10⁻¹³ m | `test_layers_add_up_to_the_ground` |
| a layer between 3 and 8 m of depth | 5 m / u_z | `test_a_layer_holds_its_thickness_over_the_sine_of_elevation` |
| two banks of a trench, three layer tops, one at its floor: every layer against the analytic path | within 2·10⁻¹⁴ m | `test_two_banks_and_layers_against_the_analytic_path` |
| chords of the cap (box), the nine piles (cylinders, in three parts each) and the void (sphere) of the invented pier | within 1·10⁻¹⁴, 3·10⁻¹² and 4·10⁻¹³ m | `test_chords_of_sphere_cylinder_and_box` |
| a structure in the air, an anomaly inside a structure | within 10⁻⁹ m | `test_a_structure_in_the_air_and_an_anomaly_inside_a_structure` |
| opacity of the table against the grammage Mulder sums along the same lines | within 10⁻¹² | `test_opacity_is_density_times_length_as_mulder_sums_it` |
| a hypothesis on the table against a scene rebuilt with those densities | lengths identical, opacity within 10⁻¹² | `test_a_hypothesis_without_a_new_scan_equals_a_rebuilt_scene` |

## What a scan refuses

| Refused | Why |
|---|---|
| a point outside the box of the scene | Mulder 0.3.7 ends the interpreter with a segmentation fault |
| a direction with a nan | Mulder never returns |
| a direction that is not a unit vector (10⁻⁹) | Mulder divides the lengths by its length: half the path for a vector of length two |
| a direction with u_z ≤ 0 | conventions: directions point to the sky |
| a point within 1 µm of a boundary between bodies | the medium there is ambiguous; the frozen cases keep the detector 1 mm above the niche floor. Inside a cut only the surface of the cut counts |
| a line of sight that crosses a limit of a ground inside it | beyond its bottom, its faces across y and the ends of its profile the model has nothing: the path would be short without a sign of it. The message gives the number of such lines and their elevations |
| a line of sight within 1 nm of an edge of a mesh at less than 0.57° to a face there, or along the edge | Mulder 0.3.7 may never return or give a wrong path along a line through an edge (next section), and a move of a few nanometres changes the path of such a line by more than a micrometre. Any other line through an edge is moved off it, not refused |

The first three are reproduced by a geometry of two boxes (the third prints 30 and 5 m where the first prints 60 and 10 m):

```bash
python -c "
import calzone, mulder
box = {'box': [1000.0] * 3, 'material': 'G4_WATER'}
world = {'box': [10000.0] * 3, 'material': 'G4_AIR', 'A': box, 'B': box | {'position': [1000.0, 0.0, 0.0]}}
geometry = mulder.LocalGeometry(calzone.Geometry({'World': world}))
print(geometry.scan(position=[[0.0, 0.0, -20.0]], direction=[[0.0, 0.0, 1.0]]))
print(geometry.scan(position=[[0.0, 0.0, -20.0]], direction=[[0.0, 0.0, 2.0]]))
# geometry.scan(position=[[0.0, 0.0, 60.0]], direction=[[0.0, 0.0, 1.0]])            # segmentation fault
# geometry.scan(position=[[0.0, 0.0, -20.0]], direction=[[float('nan'), 0.0, 1.0]])  # never returns
"
```

## Lines of sight through an edge of a mesh

Roadmap 0.1-10b, decided on 27 November 2026; module `muoscene.geometry.edges`. The pieces of a ground and of its layers are closed meshes of triangles (`builder.prism_mesh`): every side face of a prism is two triangles that share a diagonal, the caps are triangulated. Where a line of sight meets an edge or a vertex of these triangles exactly, Calzone 1.2.2 and Mulder 0.3.7 fail without an error:

| Line of sight | What Mulder gives | Where |
|---|---|---|
| leaves the ground through the diagonal of a face | the flux never returns: the navigator of Geant4 turns in `G4Navigator::CheckOverlapsIterative`, sampling points of the tessellated solid. The scan gives no path in the ground along such a line: 0 m instead of 2 m along the lines 3 and 9 of the pencil of `examples/slab`, and along 57 of the 114 lines of the map of its film under 2 m through the diagonal of the top, up to 2.83 m lost | macOS 27.3, arm64, Python 3.14.5 (27 November 2026); not seen on Linux |
| touches a concave edge, as at the foot of a canyon wall, passing from the ground under the floor into the bank | the path beyond the edge is taken for air: 29.48 m of solid instead of 51.84 m in `test_lines_through_edges_get_their_analytic_path` | Linux and macOS |
| leaves through the diagonal of a face, meshes traversed with `bvh` | no path in the ground at all | Linux and macOS; the builder uses `voxels` |

Such lines come from symmetric scenes, not by chance: a detector at the centre of a square ground looks through the diagonal of its top along every direction with tx = −ty. The trap is narrow. On macOS the flux along the lines 3 and 9 of the pencil of `examples/slab` never returned with their point moved 1e-14 m across the diagonal and was right from 1e-13 m; on Linux the line through the foot of the step below has the wrong path 1e-13 m off the edge and the right one from 1e-12 m. Of all the lines tried on macOS — the pencil and other lines through the diagonal of the slab, the edges of an invented canyon of one layer, the edges and a vertex of a box structure and the edges of a niche — every one that failed passed within 1e-11 m of an edge of a mesh; the box and the niche did not fail.

What is done. `Scene.scan` and `muoscene.site.Site.differential_flux` give Mulder every line through `Scene.traced_lines`: a line within 1 nm (`edges.TOLERANCE_M`) of an edge or a vertex of a mesh has its point moved by 3 nm across the line and the edge; where that lands near another edge — a vertex, or the crossing of the diagonals with which two pieces split a face they share at the top of a layer — the next of 16 fixed moves of 3 and 9 nm across the line is taken. The direction is kept, the path table keeps the point as given, and the moves are the same on every system. A line at less than 0.57° to a face at the edge it meets is refused (previous section). The move changes the path by at most its size over |u · n| of a face the line crosses: 0.9 µm at the edge at worst, 3 nm over the sine of the elevation at a flat top. The edges are those of every mesh as written for Calzone: the top of a ground and of each layer, the walls and slopes of the profile, the faces where pieces meet at the top of a layer, and the bottom, the faces across y and the faces under the ends of the profile, which are limits of the model; the parts of meshes inside a cut are kept. Boxes, cylinders and spheres are solids of Geant4, not meshes; cylinders and spheres were not tried.

| Check | Result | Test |
|---|---|---|
| the pencil of the film under 2 m of `examples/slab` | lines 3, 6, 9 and 12, those with tx = −ty, moved by 3 nm across; the others bit for bit as given | `test_the_lines_of_the_pencil_through_the_diagonal_are_moved` |
| the map of that film | 114 of 20 388 lines moved; none left within 1 nm of an edge, by a distance found independently; a second pass moves nothing | `test_after_the_move_every_line_is_off_every_edge` |
| lines through the edges of an invented canyon: a brow, the top of a layer on a wall, the crossing of two diagonals, a slope, both feet of the walls, an edge and a vertex at the limits of the model | the analytic path within 10⁻⁷ m, met 9·10⁻⁹ m; without the move the two feet get 29.48 and 22.36 m instead of 51.84 and 44.72 m | `test_lines_through_edges_get_their_analytic_path` |
| a line 0.006° from the wall through its brow; lines through the diagonal of a top at 0.3° and 1° of elevation | refused, refused, moved | `test_a_line_that_grazes_a_face_at_an_edge_is_refused` |
| the lines of `examples/canyon_demo` | none within 0.3 mm of an edge, none moved | `test_the_lines_of_the_canyon_demo_pass_no_closer_than_0_3_mm` |
| the points `Site.differential_flux` gives Mulder | those of `traced_lines`, in the axes of a scene turned by 30° | `test_the_flux_of_mulder_goes_along_the_traced_lines` |
| the flux of Mulder along the 16 lines of the pencil, in a child process with a timeout of 300 s | every line answers; the four through the diagonal have the flux of their mirror within 6·10⁻¹⁵ on Linux | `test_the_lines_through_the_diagonal_of_the_top_answer` |
| the maps of `examples/slab` against those without the move | 56 of 7944 bins change, by at most 3.4·10⁻¹¹ (the film under 20 m, whose moves are not horizontal); `examples/canyon_demo` bit for bit | the command below |
| the cost | 1.2 and 2.2 % of the time of a scan on macOS, with 18 and 54 edges; 5.34·10⁴ directions/s with 34 media against 5.4·10⁴ before; 8 to 10 % on the Linux VM | the command of "Run on another machine" |

Calzone and Mulder alone, without muoscene, on a step — a polygon extruded along y, 24 triangles — and on a box of 12 triangles whose top is split along x = y. On Linux and macOS the step gives 10.770330 m of rock instead of 32.310989 m on the line through its concave edge and on that line moved 1e-13 m off it, and the right path from 1e-12 m; on Linux the box gives the right path with `voxels` and none with `bvh`, and its flux 0.914936 m⁻² s⁻¹ sr⁻¹ GeV⁻¹; on macOS it gives none with either, and its flux along the diagonal never returns (27 November 2026). On macOS the same box 40 or 4 m wide answered: of those tried, the one of 400 m is the smallest that hangs.

```python
import pathlib, tempfile
import calzone, mulder, numpy as np


def geometry(vertices, faces, world_m, algorithm="voxels"):
    """Calzone and Mulder alone: one closed mesh of rock (OBJ with a normal per face) in a box of air."""
    v = np.array(vertices, dtype=float)
    lines = [f"v {x!r} {y!r} {z!r}" for x, y, z in vertices]
    for a, b, c in faces:
        n = np.cross(v[b] - v[a], v[c] - v[a])
        lines.append("vn {!r} {!r} {!r}".format(*(n / np.linalg.norm(n)).tolist()))
    lines += [f"f {a + 1}//{k} {b + 1}//{k} {c + 1}//{k}" for k, (a, b, c) in enumerate(faces, start=1)]
    with tempfile.TemporaryDirectory() as folder:
        path = pathlib.Path(folder) / "mesh.obj"
        path.write_text("\n".join(lines) + "\n")
        mesh = {"mesh": {"path": str(path), "units": "m", "algorithm": algorithm}, "material": "G4_WATER"}
        world = {"box": [100.0 * world_m] * 3, "material": "G4_AIR", "Mesh": mesh}  # Calzone takes cm
        g = mulder.LocalGeometry(calzone.Geometry({"World": world}))
    g.media[0].material, g.media[1].material = "Air", "Rock"
    return g


def prism(polygon, y0, y1):
    """A counter-clockwise polygon (x, z) extruded along y: each side two triangles, the caps a fan."""
    n = len(polygon)
    vertices = [(x, y0, z) for x, z in polygon] + [(x, y1, z) for x, z in polygon]
    faces = [f for i in range(n) for f in ((i, n + i, n + (i + 1) % n), (i, n + (i + 1) % n, (i + 1) % n))]
    faces += [f for k in range(1, n - 1) for f in ((0, k, k + 1), (n, n + k + 1, n + k))]
    return vertices, faces


# a step with its concave edge at x = z = 0; the line goes from under the step through that edge into its top
step = geometry(*prism([(-20, -20), (20, -20), (20, 0), (0, 0), (0, 10), (-20, 10)], -20.0, 20.0), 60.0)
o, u = np.array([10.0, 0.0, -4.0]), np.array([-10.0, 0.0, 4.0]) / np.hypot(10.0, 4.0)
for offset in (0.0, 1e-13, 1e-12, 3e-9):
    rock = float(
        np.asarray(step.scan(position=[o + offset * np.array([u[2], 0.0, -u[0]])], direction=[u])).ravel()[1]
    )
    print(f"step, line {offset:g} m off the edge: {rock:.6f} m of rock, analytic {32.310989:.6f}")
# a box 400 m wide, 60 m deep, its top split along x = y; the line from 2 m under its centre leaves through it
box = [
    (-200, -200, -60),
    (200, -200, -60),
    (200, 200, -60),
    (-200, 200, -60),
    (-200, -200, 0),
    (200, -200, 0),
    (200, 200, 0),
    (-200, 200, 0),
]
faces = [
    (0, 2, 1),
    (0, 3, 2),
    (4, 5, 6),
    (4, 6, 7),
    (0, 1, 5),
    (0, 5, 4),
    (1, 2, 6),
    (1, 6, 5),
    (2, 3, 7),
    (2, 7, 6),
    (3, 0, 4),
    (3, 4, 7),
]
u = [-0.000749999578125356, -0.000749999578125356, 0.9999994375004746]
for algorithm in ("voxels", "bvh"):
    g = geometry(box, faces, 402.0, algorithm)
    print(
        f"box, {algorithm}: {float(np.asarray(g.scan(position=[[0.0, 0.0, -2.0]], direction=[u])).ravel()[1]):.9f} m of rock, analytic {2.0 / u[2]:.9f}"
    )
meter = mulder.Fluxmeter(geometry=geometry(box, faces, 402.0))
meter.mode = "continuous"
print(
    "box, flux at 10 GeV:",
    float(np.asarray(meter.flux(position=[0.0, 0.0, -2.0], direction=u, energy=10.0)).ravel()[0]),
)  # macOS: never returns
```

The maps of `examples/slab` with and without the move, in one process, on Linux, where the paths without it are right for these films:

```bash
python -c "
import importlib.util, numpy as np
from muoscene.geometry import edges
spec = importlib.util.spec_from_file_location('slab', 'examples/slab/run.py')
slab = importlib.util.module_from_spec(spec)
spec.loader.exec_module(slab)
def counts():
    maps = slab.run(check_directions=0)['maps']
    return np.concatenate([m['counts'].values.ravel() for m in maps.values()])
moved = counts()
edges.move_off_edges = lambda e, o, u: (np.array(o, dtype=float), np.zeros(len(o), dtype=bool))  # as before 0.1-10b
before = counts()
ok = np.isfinite(before) & (before > 0)
change = np.abs(moved[ok] / before[ok] - 1.0)
print(f'{np.count_nonzero(change)} of {ok.sum()} bins change, by at most {change.max():.1e}')
"
```

Not done: a spatial index of the edges. The check takes every edge for every line of sight, which is cheap for the prisms of version 0.1 (18 to 120 edges) and not for the meshes of terrain of version 0.2.

## Air of the scene

Mulder gives the air of a local geometry one density — by default that of its `Air`, 1.205 kg/m³, sea level — and follows its atmosphere model only above the box. At 1150 m the US standard atmosphere has 1.0997 kg/m³. Both the density and the height of the box matter. Frozen cases `scene_air/1150m/…`: the reference void and the open sky at 45° from its niche towards the canyon, above 1 GeV, relative to a box that ends 1 m above the brow with air at the density of the site (`test_what_the_air_of_a_local_geometry_costs`):

| Air above the brow | Counts of the reference void | Open sky |
|---|---|---|
| 1000 m at Mulder's default density | −0.91 % | −1.55 % |
| 1000 m at the density of the site | −0.31 % | −0.55 % |
| 20 m at the default density (the scene of the frozen reference void) | −0.012 % | −0.030 % |
| 20 m at the density of the site | −0.0002 % | −0.0008 % |

A uniform kilometre of air at the density of the site is still too heavy, because the real air thins upwards. A scene therefore takes the density of its atmosphere model (`atmosphere`, the same names as for the flux tables) at the altitude of the site for its box and its cuts, and ends the box `margin_m` above the highest body. Against Mulder without a local geometry the open sky is then within 2.6·10⁻⁶ with Mulder's own air (`scene_air/1150m/atmosphere`) and within 2.5·10⁻⁷ with the air of the library, which has the same composition and CSDA ranges 3·10⁻⁵ to 7·10⁻⁵ longer (`test_open_sky_through_the_scene_air_is_that_of_the_atmosphere`). The effect of the void moves by less than 0.03 percentage points in every row. A detector must lie inside the box: one that stands higher than 1 m above every body needs a larger `margin_m`.

A void is not a zero for Mulder: it reads a density of exactly zero as "not set" and takes the density of the material, so that the flux above 1 GeV along the axis of the reference void comes out 2.1 % below that through the soil instead of 8.7 % above it. The medium of a body of zero density therefore gets 10⁻⁹ kg/m³; in the path table it stays zero (`test_a_void_is_empty_for_mulder_too`).

## Calzone in practice

- `subtract` works between volumes of one mother only.
- Meshes come from files only, and two traversal algorithms exist. The builder asks for `voxels`. With `bvh` lines of sight that cross the face shared by two touching meshes lose path: in a ground with one layer top 8085 of 400 000 lines differ by up to 5.9 m and the sum is 0.57 % low, with three layer tops 9044 lines, up to 11.4 m and 0.65 %; a single mesh agrees to round-off.

```bash
python -c "
import numpy as np
from muoscene import geometry as G, materials as M
from muoscene.geometry import builder
def paths(algorithm, tops):
    builder.MESH_ALGORITHM = algorithm
    layers = [G.Layer(f'layer {i}', z, M.STANDARD_ROCK) for i, z in enumerate(tops)]
    profile = [(0, -400), (0, 0), (1500, 0)]
    scene = G.Scene([G.Ground('ground', profile, -400, -1500, 1500, M.STANDARD_ROCK, layers=layers)], altitude_m=0.0)
    t = np.random.default_rng(1).uniform(-1.2, 1.2, (400000, 2))
    u = np.column_stack([t, np.ones(len(t))])
    return scene.scan((0.25, 0.0, -20.0), u / np.linalg.norm(u, axis=1)[:, None]).length_m[:, :-1].sum(axis=1)
for tops in ((), (-3.0,), (-3.0, -6.0, -9.0)):
    bvh, voxels = paths('bvh', tops), paths('voxels', tops)
    print(tops, bvh.sum() / voxels.sum() - 1, np.count_nonzero(np.abs(bvh - voxels) > 1e-6), np.abs(bvh - voxels).max())
"
```

- The overlap check of Calzone (Geant4 samples points on every surface) passes volumes that touch and reports those that overlap; `Scene.check` turns its report into a `ValueError`.
- The materials of the Calzone volumes are placeholders, `G4_AIR` and `G4_WATER`: the materials and densities of a scene are on the media of `Scene.mulder`.

## Speed

`Scene.scan`, checks included, on the invented pier of the spike seen from its niche, 40 401 directions, best of three (`test_scan_speed_on_the_pier`, which asserts the 5·10⁴ directions/s of the roadmap and is marked `validation`: `python -m pytest tests/unit/test_geometry_scene.py -m validation`): 4.0·10⁵ directions/s with 14 media and 2.7·10⁵ with 34 media (two layer tops, every pile in three parts) on a Linux VM with 4 cores, arm64. On macOS itself, on the same MacBook Air (Python 3.14), the command below gives 5.9·10⁴ and 5.4·10⁴ directions/s: the roadmap is met by 18 % and 8 %, and Linux on the same machine is 5 to 7 times faster. The bare `scan` of Mulder gave 6.5·10⁴ there in the spike, a first call included (`docs/notes/mulder.md`, E7): the numbers are close, which points at Mulder and Calzone on macOS rather than at the checks of `Scene.scan`; the two were not timed apart. A scene of this size is built in 2 ms. Since roadmap 0.1-10b the scan also moves lines of sight off the edges of the meshes (above): on 27 November 2026 the command below gave 5.85·10⁴ and 5.34·10⁴ directions/s on macOS, of which the check of the 18 and 54 edges takes 1.2 and 2.2 % — the roadmap is met by 17 % and 7 % — and 3.6·10⁵ and 2.3·10⁵ on the Linux VM, where the check takes 8 to 10 % (it varies from run to run there).

## Example

The reference void behind a wall, in loam of 1.60 g/cm³ dry density and 15 % of water (`test_example_of_the_note`):

```python
import numpy as np

from muoscene import geometry as G
from muoscene import materials as M

soil = M.loam(1.60, 15.0)  # bulk density 1.84 g/cm³
bank = G.Ground(
    "bank",
    profile_xz_m=[(0, -200), (0, 0), (400, 0)],  # a vertical wall at x = 0, a flat brow
    z_bottom_m=-200,
    y_min_m=-400,
    y_max_m=400,
    material=soil,
)
niche = G.Cut("niche", G.Box((0, -0.15, -10.001), (0.5, 0.15, -9.7)))
void = G.Anomaly("void", G.Sphere((3, 0, -3), 0.5), drho_rel=-1.0)
scene = G.Scene([bank, niche, void], altitude_m=1150.0)

detector = (0.25, 0.0, -10.0)  # 1 mm above the floor of the niche
towards_void = np.array([2.75, 0.0, 7.0]) / np.hypot(2.75, 7.0)
table = scene.scan(detector, towards_void)

table.names  # ("bank", "niche", "void", "air")
table.length_in_m("void")  # 1.0 m, its diameter; 9.4217 m in the bank
table.opacity_g_cm2()  # 1733.6 g/cm² with the void
table.opacity_g_cm2(drho_rel={"void": 0.0})  # 1917.6 g/cm² without it
table.opacity_g_cm2({"bank": 1.70}, {"void": -0.2})  # 1737.7 g/cm²: a drier bank, a decompaction zone
```

## Not in version 0.1

Rotations of bodies (the frame is turned in 0.1-07); closed meshes, terrain models, survey import and layer horizons that are not planes (0.2), and with them a spatial index of the edges of meshes for `edges`; cuts through structures; structures across the surface; anomalies of other shapes and anomalies split at boundaries; the voxel overlay.

## Run on another machine

```bash
# python -m: these tests import tests.regression, which a bare pytest on one file does not find
python -m pytest tests/unit/test_geometry_scene.py tests/unit/test_geometry_pathtable.py tests/unit/test_geometry_checks.py
python -c "
import sys, time
import numpy as np
sys.path.insert(0, 'tests/unit')
import test_geometry_scene as T
from muoscene.geometry import edges
u = T.slope_grid(401, 1.2).reshape(-1, 3)
origin = T.SPIKE.detector(T.SPIKE.PIER_NICHE)
o = np.broadcast_to(np.asarray(origin, dtype=float), u.shape)
for tops in ((), T.PIER_LAYERS):
    scene = T.pier_scene(tops)[0]
    scene.scan(origin, u[:1000])
    best = check = float('inf')
    for _ in range(3):
        start = time.perf_counter()
        scene.scan(origin, u)
        best = min(best, time.perf_counter() - start)
        start = time.perf_counter()
        edges.move_off_edges(scene._edges, o, u)  # the part of a scan that 0.1-10b added
        check = min(check, time.perf_counter() - start)
    print(f'{len(scene.mulder.media)} media: {len(u) / best:.1e} directions/s; the edges {100 * check / best:.1f} % of it')
"
```

With these 160 801 directions: 4.0·10⁵ and 2.6·10⁵ directions/s on the Linux VM, 5.9·10⁴ and 5.4·10⁴ on macOS (MacBook Air, arm64, Python 3.14; 17 November 2026); with the check of the edges, 3.6·10⁵ and 2.3·10⁵ (the edges 8 to 10 %) and 5.85·10⁴ and 5.34·10⁴ (1.2 and 2.2 %; 27 November 2026).
