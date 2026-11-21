# Site and detectors

Roadmap task 0.1-07, 21 November 2026, Mulder 0.3.7, Calzone 1.2.2. Modules `muoscene.site` (layer 1) and `muoscene.detectors` (layer 2); conventions, revision 4, sections "Frames", "Detector pose", "Direction grid" and "Energy and flux". Every number below comes from the named test of `tests/unit/test_site_frames.py`, `test_detectors_pose.py`, `test_detectors_grid.py`, `test_detectors_emulsion.py` or `tests/regression/test_legacy_track_density.py`, or from a frozen case of `tests/regression/data/`; the headline ones are in `docs/numbers.md`. Sizes in tests and in the example are those the repository already had: the bank and the niche of `tests/regression/legacy_cases.py`, the loam and the site of `mulder_cases.py`, the example of `docs/notes/geometry.md`.

## What is there

| Object | Module | What it is |
|---|---|---|
| `Site` | `site.frames` | the origin of the world frame (latitude, longitude, altitude) and the azimuth of the axes of a scene; `Site.scan` |
| `Pose` | `detectors.base` | reference point of a detector in ENU and its azimuth, tilt and roll; axes, slopes, directions |
| `DirectionGrid` | `detectors.grid` | rectangular bins in the slopes with exact solid angles, also within an acceptance circle; sub-samples |
| `LinesOfSight` | `detectors.base` | directions towards the sky, each with its solid angle: what a scene is scanned along and a flux is summed over |
| `Detector` | `detectors.base` | the interface: accepted circle, energy response, efficiency, sensitive area against the slopes |
| `EmulsionStack`, `TrackDensity` | `detectors.emulsion` | films at one point; the density of muon tracks above the energy response, with what it holds for |

## Decisions of 21 November 2026

| Question | Decided | Alternatives and their price |
|---|---|---|
| Where the turn between ENU and the axes of a scene lives | in `Site`: it holds the azimuth and `Site.scan` is the only place that applies it. `Scene` gets `azimuth_deg`, which it only passes to the local frame of Mulder | `Scene.scan` taking ENU: every test of the geometry rewritten. Without `azimuth_deg` in `Scene`, Mulder would take the y axis of a turned scene for north: 0.25 % between east and west above 1 GeV with the geomagnetic field of a date |
| Bins on the edge of the acceptance circle | exact solid angle of the part of every bin, and of every sub-sample, within the circle | the share of sub-samples inside: the sum over the default grid off by 3.4·10⁻³ to 4.3·10⁻⁵, not monotonic, one edge bin by up to 1/(2n). Bins kept or dropped whole: a jagged acceptance |
| Default grid | slopes up to 1, bins of 0.05 (40 × 40), 4 × 4 sub-samples | none was defined before; the cone is that of the preliminary calculation |
| Where a detector ends and the forward model begins | a detector takes a flux table and the opacity along its lines of sight and calls the table with its own energy response; where the opacity comes from is roadmap 0.1-08 | a ready flux per line: the threshold of a result could not be vouched for. The scene inside the detector: 0.1-08 would rewrite it |
| Track density and the two sides of a film | summed over both sides; the acceptance of the conventions (u · z_d > 0) stays for the counts of the map | front side only: −0.48 % at a tilt of 20°, −9.4 % at 45° in the conditions of the preliminary calculation, which counts both (`test_the_back_side_of_a_tilted_film`) |
| Parity with the preliminary calculation | through a scene and `Site.scan`, with `LegacyTable` and `LEGACY_RESPONSE`; `legacy.soil_length` only to tell the terms of the difference apart | `legacy.soil_length` alone would not test the chain a user runs |
| Area and films of a stack | `area_cm2` is the sensitive area (a film of 120 mm × 100 mm has 120 cm²; stacks side by side add up); `n_films` counts the films along the normal and does not multiply the area | — |

Two of them are not in the conventions yet and are proposed for their revision 5: the azimuth of the axes of a scene ("Frames" calls the world frame the local frame of Mulder, which holds only for a scene that is not turned), and the two sides of a film in a track density ("Acceptance").

## Site and the axes of a scene

A scene of `muoscene.geometry` has no rotations (`docs/notes/geometry.md`): a canyon wall is the plane x = 0 of the scene wherever it runs. `Site.scene_azimuth_deg` ψ is the azimuth of the +y axis of the scene from true north, clockwise; the two frames share origin and vertical:

    x_s = (cos ψ, −sin ψ, 0),   y_s = (sin ψ, cos ψ, 0),   z_s = (0, 0, 1)        in ENU

At ψ = 0 the scene is ENU; at 90° its +y points east and its +x south. `Site.scene(bodies)` builds the `Scene` with the altitude, origin and azimuth of the site; `Site.scan(scene, position_m, direction)` takes points and directions in ENU, turns them, scans, and returns a `PathTable` that holds them in ENU as given. A scene of another site is refused. Poses, grids and detectors never see the axes of a scene.

| Check | Result | Test |
|---|---|---|
| paths worked out by hand in a scene turned by 30°: into the massif (azimuth 120°) through the roof and through the back wall of the niche, towards the canyon (300°) through the lip, along the wall, and at 45° to it | within 2·10⁻¹⁵ m of 13.5100, 23.7113, 0.1536, 16.7584, 16.7584 and 16.6706 m; with the azimuth counted the other way the paths differ by more than 0.5 m | `test_paths_worked_out_by_hand_in_a_turned_scene` |
| the turn against `LocalFrame(azimuth=ψ)` of Mulder | vectors agree to 9·10⁻¹⁶ | `test_the_azimuth_is_that_of_the_local_frame_of_mulder` |
| `Site.scan` against `Scene.scan` along the axes of the scene | paths within 10⁻¹¹ m; identical at ψ = 0 | `test_scan_is_the_scan_of_the_scene_along_its_own_axes` |
| a scene turned by 90° tells Mulder where north is: open-sky flux above 1 GeV at 45° along ±y of the scene, with the field of 15 October 2028 | equal to the geographic east and west within 10⁻⁵; east is 0.25 % above west | `test_a_scene_of_a_site_tells_mulder_where_north_is` |
| `Scene(azimuth_deg=…)` changes no path | identical | `test_scene_azimuth_changes_no_scan` |

Mulder knows north through the `azimuth` of its `LocalFrame`, with the same definition as ψ. It matters for the geomagnetic field only: a scan does not depend on it. `MulderTable(azimuth_deg=…)` is a geographic azimuth already, because a table is built without a local geometry.

## Pose

Axes and slopes are those of the conventions; the roll turns x_d and y_d about the normal, counter-clockwise seen from its tip. The azimuth is kept in [0°, 360°), the tilt lies in [0°, 90°], the roll is kept in (−180°, 180°]. At a tilt of zero the azimuth and the roll turn about the same axis in opposite senses: (α, 0, γ) and (0, 0, γ − α) are one orientation, and `Pose.from_axes` returns the second.

| Check | Result | Test |
|---|---|---|
| a film tilted by 30° towards the east: axes, and the slopes of the vertical, (0, −tan 30°) | as worked out by hand, to 10⁻¹⁵ | `test_axes_worked_out_by_hand` |
| angles → axes → angles, 2000 poses with a tilt of at least 0.01° | within 6·10⁻¹⁴ degrees | `test_angles_to_axes_to_angles` |
| direction → slopes → direction, 4000 directions | within 3·10⁻¹⁶; slopes within 2·10⁻¹⁵ | `test_direction_to_slopes_to_direction` |

`Pose.direction(tx, ty)` is the direction on the side of the normal. `Pose.sky_direction(tx, ty)` is the side of the same line that looks at the sky (u_z ≥ 0), with a flag that tells whether it is in front of the plane of the detector or behind it.

## Direction grid

Bins are rectangles in (tx, ty); the solid angle of a rectangle is four corner terms arctan(a b / sqrt(1 + a² + b²)), for slopes of either sign (`test_negative_slopes`: against a numerical integral to 10⁻¹⁰, and by symmetry).

The acceptance is the circle sqrt(tx² + ty²) ≤ `tan_max`. The solid angle of a region around the axis is ∮ (1 − cos θ) dφ along its boundary. On a straight side at the distance d from the axis, with ψ the polar angle from the foot of the perpendicular, the integral is ψ − arcsin(sin ψ / sqrt(1 + d²)); on an arc of the circle it is (1 − cos θ_max) Δφ. The part of a rectangle within the circle is therefore exact (`clipped_solid_angle_sr`): within 4·10⁻¹⁴ sr of quadrature for 200 random rectangles (`test_cut_rectangles_against_a_numerical_integral`).

| Grid | Sum of the accepted solid angles against 2π (1 − cos θ_max) | Test |
|---|---|---|
| default: `tan_max` 1, bins of 0.05; 1176 bins inside the circle, 148 cut, 276 outside; the circle touches the grid at four points | within 10⁻¹³; met 2·10⁻¹⁵ | `test_the_default_grid_adds_up_to_the_cone` |
| bins of 0.1 and 0.07, `tan_max` 0.6, 0.5 through corners, sqrt(15), 10, and a circle of 0.02 inside the four central bins | within 10⁻¹³ | `test_bins_add_up_to_the_cone` |
| the alternative, the share of sub-samples inside the circle, default grid | +3.4·10⁻³, −3.5·10⁻⁴, +4.3·10⁻⁵, +2.0·10⁻⁴ for 1, 2, 4, 8 sub-samples per side; 7.9·10⁻⁴ for bins of 0.07 and 4 per side | `test_the_share_of_sub_samples_inside_the_circle_is_not_exact` |

`subsample(n, tan_max)` divides every bin into n × n rectangles of equal sides in the slopes. Each has its exact solid angle within the circle and one direction: its centre, or, when the centre is outside the circle, a point halfway from its corner nearest to the axis to where the line towards its centre leaves the circle. The weights of a bin add up to its accepted solid angle within 5·10⁻¹⁴ sr (`test_sub_samples_lie_in_their_bin_and_in_the_circle`). What is left is the integrand taken at one point of each rectangle: the integral of cos θ over the cone of 45° is 5.3·10⁻⁵ low with one direction per bin of the default grid and 2.6·10⁻⁶ low with 4 × 4 (`test_sub_samples_integrate_a_smooth_function`).

`DirectionGrid.square(tan_max, bin_width)` is the grid of the maps; `DirectionGrid.equal_angle(tan_max, n_bins)` has edges equally spaced in angle, for wide fields of view. A grid gives `tx_edges`, `ty_edges`, `shape` (bins in ty, bins in tx), centres, `solid_angle_sr`, `accepted_solid_angle_sr(tan_max)` and `directions(pose)`, and needs neither a scene nor a detector.

## Detector and lines of sight

A `Detector` has `tan_max`, `response` (an `EnergyResponse`, always given) and `efficiency`, and implements `sensitive_area_cm2(tx, ty)`, the area in its plane a track of given slopes can cross: A for a film, max(a − D|tx|, 0) · max(b − D|ty|, 0) for two planes a × b at the distance D. `projected_area_cm2` is that area times cos θ = 1 / sqrt(1 + tx² + ty²). A telescope written against the interface in `test_a_detector_whose_area_depends_on_the_slopes` has the geometric factor a b π sin² θ_max at D = 0 within 2.6·10⁻⁶ and that of a numerical integral within 2.9·10⁻⁵ at D = 30 cm.

`Detector.lines_of_sight(pose, grid, supersample)` gives the accepted lines: within `tan_max`, in front of the plane, above the horizon. `lines_of_sight(pose, grid, tan_max=…, both_sides=…, sin_elevation_min=…)` gives any other set, and `LinesOfSight(pose, direction, solid_angle_sr)` takes directions made elsewhere. The horizon, an elevation limit and the back side are cut by whole sub-samples; for a horizontal detector an elevation limit given as `tan_max` is exact.

`Detector.flux_m2_s_sr(table, lines, opacity_g_cm2)` calls `table.integral` with the response of the detector. This is the border with the forward model (roadmap 0.1-08): the detector does not know where the opacity comes from — a scan, a density hypothesis, the equivalent opacity of several materials — and the forward model does not choose the threshold. Counts per bin are the same sum with the projected area, the efficiency and `LinesOfSight.bin`, and belong to 0.1-08 with the `ExpectedMap`.

## Track density

`EmulsionStack.track_density(table, lines, opacity_g_cm2, exposure_days)` is the sum over the lines of sight of flux × |u · z_d| × solid angle × time, per cm² of film. It returns a `TrackDensity`: the number with the energy response, the exposure, the cone, the lowest elevation, the sides of the film and the number of lines it was summed over. `EmulsionStack.track_lines(pose, sin_elevation_min=…)` makes the lines: by default tan 10 from the normal (84.3°, as in the preliminary calculation), both sides, 240 bins a side in equal angles with 2 × 2 sub-samples — 2.2·10⁵ lines for a horizontal film; with 480 bins a side the numbers below move by less than 3·10⁻⁵. The elevation limit is always given.

What it counts and what it does not:

- muons the energy response lets through; no efficiency of reconstruction and no area — the film holds every track;
- both sides of a film. Under a uniform sky the density is then the same at any tilt, F T π · 100/101 within tan 10 of the normal (`test_a_film_has_two_sides`, within 2.9·10⁻⁵);
- no soft component (electrons and positrons) and no muons below the threshold: with it, a budget against the scanning limit is a lower bound. Roadmap 0.1-05b;
- no tracks collected before and after the exposure: version 0.3.

**Thresholds and the cut-off.** With the sea-level reference at 1150 m the open sky exists above 0.28 GeV by the vertical and 0.98 GeV at 15.7° (`docs/notes/flux.md`). A film in a niche of a wall looks at open sky with most of its field of view, so

- above 1 GeV the density is defined for elevations from 15.7° up; the numbers below are from 16°. In the conditions of the preliminary calculation 0.84 and 0.86 % of the density of a horizontal film in the niches at 12 and 22 m come from below 16° (`test_where_the_tracks_of_a_film_in_a_niche_come_from`);
- above 0.2 GeV it is defined only for lines of sight under soil. The side of the canyon, 76 to 86 % of the tracks above 1 GeV, has no flux above 0.2 GeV until 0.1-05b: `track_density` raises the `ValueError` of the table, with the cut-off and the direction in it, and nothing is returned as zero (`test_open_sky_below_the_cutoff_is_refused`). `LinesOfSight.select` takes the lines under soil.

Horizontal film on the floor of the niche of `legacy_cases` in a flat bank of loam of 1.8 g/cm³ at 1150 m, 60 days, elevations above 16°, tracks per cm² (`test_track_density_in_a_niche_at_1150_m`):

| Niche below the brow | Above 1 GeV, all | Side of the canyon, above 1 GeV | Side of the ground, above 1 GeV | Side of the ground, above 0.2 GeV | All, above 0.2 GeV |
|---|---|---|---|---|---|
| 12 m | 34 758 | 26 244 | 8 514 | 10 081 | refused |
| 22 m | 30 346 | 26 104 | 4 242 | 4 794 | refused |

Open sky, horizontal film, 60 days, above 16°: 60 039 tracks/cm² above 1 GeV at 1150 m, 49 843 at sea level; above 0.2 GeV 64 431 at sea level and refused at 1150 m (`test_open_sky_below_the_cutoff_is_refused`).

Against Mulder through the scene, on the 1152 lines of sight of the frozen cases `track_density/1150m/12m` and `…/22m` (`test_track_density_in_a_niche_against_mulder`): −0.04 and −0.06 % on the side of the ground, above both thresholds, +0.03 % on the side of the canyon. The first is the table at the altitude of the detector, which has air where the soil is (`docs/notes/flux.md`); the second is the 20 m of sea-level air in the scene of the frozen cases (`docs/notes/geometry.md`). The frozen numbers themselves are sums over a coarse grid, 24 × 48 directions, and differ from the table above by up to 0.9 % on the side of the canyon: they are a reference for the same lines of sight, not for the integral.

Found on the way: the open-sky flux of Mulder through a local geometry is not the same number in every process. If a fluxmeter without a local geometry ran before — any of the open-sky, slab or cut-off cases — the side of the canyon of `track_density/1150m/22m` comes out 5.7·10⁻¹⁴ lower, and 1.3·10⁻¹⁴ higher at 12 m. The frozen file holds the numbers of the order of `mulder_cases.CASES`, which `make_mulder_reference --check` reproduces bit for bit; the tolerance of the frozen cases is 10⁻⁵.

```bash
python -c "
from tests.regression import mulder_cases as mc
print(repr(mc.compute('track_density/1150m/22m')['canyon/K>1GeV']))  # 26350.248433470006
"
python -c "
from tests.regression import mulder_cases as mc
mc.compute('open_sky/0m/elev90')
print(repr(mc.compute('track_density/1150m/22m')['canyon/K>1GeV']))  # 26350.248433468514
"
```

## Parity with the preliminary calculation

The frozen `track_density/…` and `track_density_tilted/…` of `legacy_reference.json` are Monte Carlo sums: 200 000 directions (150 000 kept) and 400 000, seed 5. muoscene sums the same integrand over a grid, with `LegacyTable`, `LEGACY_RESPONSE`, straight paths and the wall and the niche of the preliminary calculation as a scene. Tilted films there count directions down to cos θ = 0.05 with the intensity taken at 0.25, which `LegacyTable` refuses: the test wraps the table to repeat it.

| Frozen case | Frozen | muoscene | Difference | Monte Carlo σ of the frozen number |
|---|---|---|---|---|
| `track_density/flat/0m` (open sky) | 69 435 | 69 443 | +0.013 % | 0.28 % |
| `track_density/flat/12m` | 38 143 | 38 155 | +0.031 % | 0.34 % |
| `track_density/flat/22m` | 33 694 | 33 705 | +0.032 % | 0.39 % |
| `track_density/rising/12m` | 37 523 | 37 533 | +0.027 % | 0.35 % |
| `track_density/rising/22m` | 33 420 | 33 430 | +0.031 % | 0.39 % |
| `track_density_tilted/flat/22m/beta0` | 33 963 | 33 929 | −0.100 % | 0.21 % |
| `track_density_tilted/flat/22m/beta20` | 27 135 | 27 119 | −0.059 % | 0.22 % |
| `track_density_tilted/rising/22m/beta0` | 33 687 | 33 652 | −0.102 % | 0.21 % |
| `track_density_tilted/rising/22m/beta20` | 26 816 | 26 799 | −0.062 % | 0.22 % |

σ is the relative standard deviation of the estimator of the preliminary calculation over the seeds 100 to 139 with the same number of directions (`test_monte_carlo_error_of_the_frozen_densities`, marked `validation`); the tests allow three of them (`test_horizontal_films`, `test_tilted_films`, `test_open_sky`). What the difference is made of:

- **the integrand is the same.** On the directions of the preliminary calculation themselves, with its solid angle per direction, the density equals the frozen one to 4·10⁻¹⁶ with the path of `legacy.soil_length`, and through the scene as well under a flat brow, where the paths agree to 2·10⁻¹³ m (`test_the_same_integrand_on_the_directions_of_the_preliminary_calculation`);
- **the geometry differs under the rising brow** by −3.0·10⁻⁶ and −3.5·10⁻⁶ of the density: the preliminary calculation takes the brow as flat towards the canyon and along the wall, up to 0.33 m of path on these directions (0.15 m on those of `docs/legacy.md`);
- **the grid sum is converged** to 4·10⁻⁶ between 240 and 480 bins a side (`test_the_grid_sum_is_converged`, marked `validation`), and for the open sky lies 2.1·10⁻⁵ below the exact integral of the table of the preliminary calculation, 69 444.8 (`test_open_sky`);
- **the rest is the Monte Carlo error of the frozen number**: 0.05 to 0.5 σ.

A scene turned by 37° or 180°, with the film turned with it, gives the same density to 2·10⁻¹⁶ (`test_a_turned_scene_gives_the_same_density`).

## Example

Six films of 120 mm × 100 mm on the floor of a niche 10 m below a brow at 1150 m; the wall runs along the azimuth 30°, the massif lies at 120° (`test_example_of_the_note`):

```python
import math

from muoscene import geometry as G
from muoscene import materials as M
from muoscene.detectors import EmulsionStack, Pose
from muoscene.flux import EnergyResponse, MulderTable
from muoscene.site import Site

site = Site(latitude_deg=43.35, longitude_deg=79.08, altitude_m=1150.0, scene_azimuth_deg=30.0)
soil = M.loam(1.60, 15.0)  # bulk density 1.84 g/cm³
bank = G.Ground("bank", [(0, -200), (0, 0), (400, 0)], -200, -400, 400, soil)  # the wall is x = 0
niche = G.Cut("niche", G.Box((0, -0.15, -10.001), (0.5, 0.15, -9.7)))
scene = site.scene([bank, niche])

films = EmulsionStack(
    area_cm2=120.0, n_films=6, tan_max=1.0, response=EnergyResponse.threshold(1.0), efficiency=0.9
)
pose = Pose(tuple(site.enu_from_scene((0.25, 0.0, -10.0))))  # (0.2165, -0.125, -10.0) in ENU, horizontal
films.projected_area_cm2(0.5, 0.0)  # 107.3 cm² at a slope of 0.5

elevation = math.radians(16.0)  # 1 GeV is above the cut-off from 15.7 degrees
lines = films.track_lines(pose, sin_elevation_min=math.sin(elevation), tan_max=1 / math.tan(elevation))
paths = site.scan(scene, pose.position_m, lines.direction)  # points and directions in ENU
table = MulderTable(1140.0, soil.material)  # at the altitude of the detector
density = films.track_density(table, lines, paths.opacity_g_cm2(), exposure_days=60)
print(density)  # 3.621e+04 tracks/cm² in 60 days, K > 1 GeV, elevations above 16.0 degrees, ...

soft = EmulsionStack(
    area_cm2=120.0, n_films=6, tan_max=1.0, response=EnergyResponse.threshold(0.2), efficiency=0.9
)
soft.track_density(table, lines, paths.opacity_g_cm2(), 60)  # ValueError: ... a cut-off of 0.2427 GeV there
massif = site.scene_from_enu(lines.direction)[:, 0] > 0  # the lines of sight under soil
soft.track_density(table, lines.select(massif), paths.opacity_g_cm2()[massif], 60)  # 11 941 tracks/cm²
```

The opacity here is that of one material. With several materials the forward model gives the equivalent opacity (roadmap 0.1-08).

## Not in version 0.1

Survey coordinate systems, grid north and magnetic declination (0.2); tables of angular efficiency, latent-image fading, the track-density budget with transport and storage (0.3); `TrackSet`, the straw telescope and the two-plane telescope (0.4). Counts per bin, hypotheses and the `ExpectedMap` are roadmap 0.1-08; the open sky below the cut-off and the soft component are 0.1-05b.

## Run on another machine

```bash
# python -m: these tests import tests.regression, which a bare pytest on one file does not find
python -m pytest tests/unit/test_site_frames.py tests/unit/test_detectors_pose.py tests/unit/test_detectors_grid.py tests/unit/test_detectors_emulsion.py tests/regression/test_legacy_track_density.py
python -m pytest tests/regression/test_legacy_track_density.py -m validation    # Monte Carlo error of the frozen densities, convergence of the grid: 20 s
```

On the Linux VM with 4 cores the five files take 27 s. Not run on macOS yet, where a scan is 5 to 7 times slower (`docs/notes/geometry.md`).
