# Geometry

Roadmap task 0.1-06, 17 November 2026, Mulder 0.3.7, Calzone 1.2.2. Module `muoscene.geometry`; conventions, revision 4, sections "Frames" and "Bodies, densities and anomalies". Every number below comes from the named test of `tests/unit/test_geometry_scene.py`, `test_geometry_pathtable.py` or `test_geometry_checks.py`, from a frozen case of `tests/regression/data/mulder_reference.json`, or from the command printed next to it; the headline ones are in `docs/numbers.md`. Sizes in tests are those the repository already had: the wall, banks, niches and pile cap of `tests/regression/legacy_cases.py`, the reference void of `mulder_cases.py`, the invented pier of `docs/notes/spike/mulder_spike.py`.

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

Version 0.1 has no rotations: boxes follow the axes, cylinders are vertical, the ground runs along y. A scene whose wall does not run north–south is described in its own axes, and the frame is turned in roadmap 0.1-07 (decided on 17 November 2026).

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

`Scene.scan`, checks included, on the invented pier of the spike seen from its niche, 40 401 directions, best of three (`test_scan_speed_on_the_pier`, which asserts the 5·10⁴ directions/s of the roadmap and is marked `validation`: `python -m pytest tests/unit/test_geometry_scene.py -m validation`): 4.0·10⁵ directions/s with 14 media and 2.7·10⁵ with 34 media (two layer tops, every pile in three parts) on a Linux VM with 4 cores, arm64. On macOS itself, on the same MacBook Air (Python 3.14), the command below gives 5.9·10⁴ and 5.4·10⁴ directions/s: the roadmap is met by 18 % and 8 %, and Linux on the same machine is 5 to 7 times faster. The bare `scan` of Mulder gave 6.5·10⁴ there in the spike, a first call included (`docs/notes/mulder.md`, E7): the numbers are close, which points at Mulder and Calzone on macOS rather than at the checks of `Scene.scan`; the two were not timed apart. A scene of this size is built in 2 ms.

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

Rotations of bodies (the frame is turned in 0.1-07); closed meshes, terrain models, survey import and layer horizons that are not planes (0.2); cuts through structures; structures across the surface; anomalies of other shapes and anomalies split at boundaries; the voxel overlay.

## Run on another machine

```bash
# python -m: these tests import tests.regression, which a bare pytest on one file does not find
python -m pytest tests/unit/test_geometry_scene.py tests/unit/test_geometry_pathtable.py tests/unit/test_geometry_checks.py
python -c "
import sys, time
sys.path.insert(0, 'tests/unit')
import test_geometry_scene as T
u = T.slope_grid(401, 1.2).reshape(-1, 3)
origin = T.SPIKE.detector(T.SPIKE.PIER_NICHE)
for tops in ((), T.PIER_LAYERS):
    scene = T.pier_scene(tops)[0]
    scene.scan(origin, u[:1000])
    best = float('inf')
    for _ in range(3):
        start = time.perf_counter()
        scene.scan(origin, u)
        best = min(best, time.perf_counter() - start)
    print(f'{len(scene.mulder.media)} media: {len(u) / best:.1e} directions/s')
"
```

With these 160 801 directions: 4.0·10⁵ and 2.6·10⁵ directions/s on the Linux VM, 5.9·10⁴ and 5.4·10⁴ on macOS (MacBook Air, arm64, Python 3.14; 17 November 2026).
