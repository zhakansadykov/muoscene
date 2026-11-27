# Forward model

Roadmap task 0.1-08, 21 November 2026, Mulder 0.3.7, Calzone 1.2.2. Modules `muoscene.transport` (layer 2) and `muoscene.model` (layer 3); conventions, revision 4, sections "Energy and flux" and "Expected map". Every number below comes from the named test of `tests/unit/test_transport_table.py`, `test_transport_direct.py`, `test_detectors_counts.py` or `test_model_expected.py`, or from a frozen case of `tests/regression/data/mulder_reference.json`; the headline ones are in `docs/numbers.md`. The demonstration scene is the invented pier of `docs/notes/spike/mulder_spike.py` in the loam and the concrete of the frozen cases; `examples/canyon_demo` comes with roadmap 0.1-10.

## What is there

| Object | Module | What it is |
|---|---|---|
| `expected_map` | `model.expected` | detectors × hypotheses × bins → an `ExpectedMap`, an `xarray.Dataset` |
| `Station` | `model.expected` | a detector with its pose, the flux table of its altitude and its exposure |
| `Hypothesis`, `HypothesisSet` | `model.hypotheses` | densities a map is computed for, without a new scan |
| `equivalent_opacity`, `flux_domain` | `transport.table` | the opacity a table is entered with; where a table has a flux |
| `direct_flux`, `compare_with_direct` | `transport.direct` | the flux of Mulder through the scene; the check of the fast path |
| `Detector.counts` | `detectors.base` | tracks per bin: flux × projected area × efficiency × solid angle × time |
| `Site.differential_flux` | `site.frames` | Mulder through a scene along lines given in ENU |

The path of a run: the accepted lines of sight of a detector through the grid → `Site.scan` → `PathTable` → opacity by material under every hypothesis → equivalent opacity → the flux table, called by the detector with its own energy response → counts per bin. The map of the pier — two detectors, two hypotheses, 40 × 40 bins with 4 × 4 sub-samples, the comparison with Mulder included — takes 0.8 s on a Linux VM with 4 cores.

## Decisions of 21 November 2026

| Question | Decided | Alternatives and their price |
|---|---|---|
| Energy at which the weights of the equivalent opacity are taken | none is fixed: w_m = ΔR_ref / ΔR_m between the first energy of the detector response and the energy needed above the path, found by iteration | on 3000 g/cm² of loam under 1840 g/cm² of concrete, flux relative to the chain of ranges in the true order: this one +0.023 %; the chain in the wrong order +0.048 %; a weight at 1, 7, 100 GeV +0.17, +0.06, −0.08 %; no weights −0.17 % (`test_the_alternatives_that_were_not_taken`) |
| Altitude of the table of a detector under soil | that of the detector; the air the solid bodies displace, ρ_air times the path in them, is taken out of the opacity the table is entered with | the table alone: −0.074 to −0.078 % under 20 m; a table at the altitude of the surface: +0.03 to +0.05 % under soil and +0.31 to +0.39 % in the open sky from a niche 20 m below the brow (`test_a_table_at_the_altitude_of_the_surface_is_wrong_for_the_open_sky`) |
| Bins the table has no flux for | `undefined` is always given: `"refuse"` refuses the run, `"mask"` computes the rest; a bin with one line of sight without a flux is undefined as a whole, with the reason in `status` and nan in `counts` | refusal alone: no map above 0.2 GeV for a film in a niche at 1150 m, where 356 of its 1324 bins look at the canyon or through thin soil. Zeros in their place: a missing number taken for a count |
| How the direct flux reaches Mulder | `Site.differential_flux`: turns ENU into the axes of the scene and scans the lines first, so that `Scene.scan` refuses what crashes Mulder; the integral over the response is in `transport.direct` | the whole of it in `transport.direct`: the refusals of the scan repeated there |
| Sample of the comparison | 256 accepted lines of sight per detector that have a flux, drawn without replacement; the scene as built; refused above 0.5 % | every hypothesis compared: a Mulder geometry per hypothesis. The seed is that of the generator, written to the map |
| `solid_angle_sr` of a bin on the edge of the circle | the whole rectangle; what a detector accepts of it is `accepted_solid_angle_sr` | the accepted part under the name of the conventions: it depends on the detector, and the variable has no `detector` dimension |

The first three change what the conventions say and are proposed for their next revision: the weight with a threshold, the displaced air ("the atmosphere is part of the table and is not counted"), and the mask (the domain of validity refuses a run whose response reaches below the cut-off in any accepted direction, while numbers are asked for 0.2 GeV, which holds only under soil). The conventions were not changed by this task; all three entered their revision 5 on 27 November 2026.

## Equivalent opacity

A `PathTable` gives the opacity X_m in every material along a line of sight, not their order. The table of a detector is entered with

    X = Σ w_m X_m − w_air ρ_air L_solid,    w_m = (R_ref(K_top) − R_ref(K_min)) / (R_m(K_top) − R_m(K_min)),

with R the CSDA range, K_min the first energy of the response, R_ref(K_top) = R_ref(K_min) + Σ w_m X_m and L_solid the length of the path in solid bodies. For K_min → 0 the weight is R_ref(K) / R_m(K) of the conventions.

| Check | Result | Test |
|---|---|---|
| a path through concrete alone against R_ref(K_top) − R_ref(K_min) | within 3·10⁻⁵, the round trip of the ranges of Mulder | `test_one_material_is_exact` |
| the same at K_min = 1 MeV against `materials.equivalent_opacity` | within 10⁻⁴ | `test_without_a_threshold_it_is_the_weight_of_the_conventions` |
| loam with concrete through a table of either material | the same flux to 2·10⁻⁷ | `test_the_choice_of_the_reference_material_does_not_matter` |
| weight of the displaced air, 20 m of soil by the vertical, above 1 GeV | 1.155: a gas has no density effect | `test_displaced_air` |
| uniform slab against Mulder with the real soil, 0.5, 2 and 20 m, zenith 0° and 45°, above 0.2 and 1 GeV | −0.003 to +0.021 %; without the air term −0.007 to −0.078 % | `test_a_uniform_slab_matches_mulder`; `threshold_slab/…` |
| open sky against Mulder without a geometry, above 1 GeV at 45° | table −0.004 %, direct flux through the scene −0.005 % | `test_the_open_sky_matches_mulder`, `test_the_open_sky_through_a_scene_is_that_of_mulder`; `scene_air/1150m/atmosphere` |

What is left under soil, +0.01 to +0.02 %, is the decay of muons inside it, which a table does not have. A void of a scene is empty, also of air, for Mulder and for the table alike. With a `LegacyTable` the opacity is the plain sum of density times length, as in the preliminary calculation.

## Bins without a flux

`flux_domain` gives one code per line of sight; a bin takes the reason that holds most of its undefined solid angle.

| `status` | Meaning | `counts` |
|---|---|---|
| 0 | defined | a number |
| 1 | not accepted: outside the circle, behind the plane of the detector, below the horizon | 0 |
| 2 | the response starts below the lowest energy of the table there: the cut-off of the reference flux | nan |
| 3 | the elevation is outside the table (below 14.5°) | nan |
| 4 | the opacity is beyond the table (3·10⁴ g/cm²) | nan |

A horizontal film in the niche of the pier at 1150 m, 1324 bins in its acceptance (`test_bins_without_a_flux_are_marked_not_filled`): above 1 GeV all are defined; above 0.2 GeV 968 are and 356 are below the cut-off, and the open sky has no flux in any bin, so that `transmission` is undefined everywhere while the counts under soil exist. `reference` is geometry and marks the same 52 bins for both thresholds. Tilted by 30° towards the canyon the film has 6 bins below the cut-off of 1 GeV, between 15.1° and 15.7° of elevation (`test_a_tilted_detector_and_the_edge_of_the_table`). `expected_counts` gives the counts of one detector and hypothesis as a masked array.

## The map

Variables of the conventions: `solid_angle_sr`, `opacity_gcm2`, `counts`, `transmission`, `reference`, `zenith_deg`, `azimuth_deg`. Added: `accepted_solid_angle_sr`, `status`, `open_sky_status` (detector, ty, tx) and `exposure_days` (detector); the bin edges are the coordinates `tx_edges` and `ty_edges`, from which `grid_of` makes the `DirectionGrid`. `zenith_deg` and `azimuth_deg` are those of the bin centre on the side of the normal: above 90° the centre looks below the horizon, and the counts of the bin come from its lines of sight above it. `opacity_gcm2` is the matter, Σ w_m X_m, without the air term.

Attributes: versions of muoscene, Mulder and Calzone; creation time; `undefined`; the site; `stations` (detector, energy response, pose, exposure, table) and `hypotheses` as JSON; `check_seed`, `check_directions`, `fast_vs_direct` (per detector: mean, rms and largest deviation, where the largest is, the limit) and `fast_vs_direct_max`. The configuration hash and the file come with roadmap 0.1-11.

## Against Mulder

| Comparison | Result | Test or frozen case |
|---|---|---|
| the pier, in the run: 256 random lines per detector, seed 8, above 1 GeV | mean +0.006 %, rms 0.014 %, largest 0.032 % | `test_the_demo_is_within_half_a_percent_of_mulder` |
| the same above 0.2 GeV, on the lines that have a flux | mean +0.020 %, rms 0.022 %, largest 0.049 % | the same |
| the pier built without muoscene, 77 lines on a grid of slopes | within 0.033 % above 1 GeV; within 0.050 % above 0.2 GeV on the 34 lines into the massif | `test_the_fast_path_against_the_frozen_pier`; `pier/1150m/soil`, `…/void` |
| the void of the pier | +12.6 % on the line through it; +10.9 % in its two bins; +572 tracks of 2.398·10⁶ in the map | `test_the_void_of_the_frozen_pier_is_seen` |
| a hypothesis on one scan against a scene rebuilt with those densities | counts within 10⁻¹² | `test_a_hypothesis_without_a_new_scan_equals_a_rebuilt_scene` |
| energy nodes of the direct flux, 30 against 240 per decade | 7·10⁻⁵ above 1 GeV, 3·10⁻⁵ above 0.2 GeV | `test_the_integral_starts_at_the_threshold` |
| the open sky above 1 GeV at 15.5, 16 and 17 degrees, where 1 GeV is just above the cut-off at 1150 m: the table against the direct flux | with 30 nodes a decade from the threshold 0.62, 0.14 and 0.05 % below it, and a right run refused; with 240 a decade over the first half decade, the default since 22 November 2026 (roadmap 0.1-10), within 0.021 % | `test_the_check_at_the_edge_of_the_cut_off` |
| a scene turned by 90° against the same scene in its own axes | within 1.8·10⁻⁶ | `test_points_and_directions_are_in_enu` |

The criterion of the roadmap is 0.5 %; the spike had rms 0.17 % and a largest deviation of 0.39 % with its coarse table (`docs/notes/mulder.md`, E12).

## Example

Six films in a niche 10 m below a brow at 1150 m, tilted by 20° into the massif, and a void of 1 m behind the wall (`test_example_of_the_note`):

```python
import numpy as np

from muoscene import geometry as G
from muoscene import materials as M
from muoscene.detectors import DirectionGrid, EmulsionStack, Pose
from muoscene.flux import EnergyResponse, MulderTable
from muoscene.model import Hypothesis, Station, expected_counts, expected_map
from muoscene.site import Site

site = Site(latitude_deg=43.35, longitude_deg=79.08, altitude_m=1150.0, scene_azimuth_deg=30.0)
soil = M.loam(1.60, 15.0)  # bulk density 1.84 g/cm³
bank = G.Ground("bank", [(0, -200), (0, 0), (400, 0)], -200, -400, 400, soil)  # the wall is x = 0
niche = G.Cut("niche", G.Box((0, -0.15, -10.001), (0.5, 0.15, -9.7)))
void = G.Anomaly("void", G.Sphere((3, 0, -3), 0.5), drho_rel=-1.0)
scene = site.scene([bank, niche, void])

films = EmulsionStack(
    area_cm2=120.0, n_films=6, tan_max=1.0, response=EnergyResponse.threshold(1.0), efficiency=0.9
)
position = tuple(site.enu_from_scene((0.25, 0.0, -10.0)))
pose = Pose(position, azimuth_deg=120.0, tilt_deg=20.0)  # the massif lies at the azimuth 120
station = Station("niche", films, pose, MulderTable(1140.0, soil.material), exposure_days=60)
hypotheses = [Hypothesis.without("void"), Hypothesis.void("void")]
expected = expected_map(
    site, scene, station, DirectionGrid.square(), hypotheses, rng=np.random.default_rng(1), undefined="refuse"
)

none, with_void = (expected_counts(expected, "niche", name) for name in ("none", "void"))
none.sum()  # 1.968e6 tracks in 60 days above 1 GeV
with_void.sum() - none.sum()  # 928 more with the void, in 12 bins
(with_void / none - 1).max()  # +7.6 % in the bin that looks at its centre
expected.attrs["fast_vs_direct_max"]  # 1.5e-4: the fast path against Mulder on 256 lines of sight
```

Whether 928 tracks are seen is the business of roadmap 0.1-09.

## Not here

Significance, χ² and exposure to kσ (roadmap 0.1-09); configuration, the run folder and NetCDF (0.1-11); figures (0.1-12); scattering in maps, the geomagnetic field in maps — a table with a date holds one azimuth and is refused — and systematic variants (0.3); the open sky below the cut-off (0.1-05b).

## Run on another machine

```bash
# python -m: these tests import tests.regression, which a bare pytest on one file does not find
python -m pytest tests/unit/test_transport_table.py tests/unit/test_transport_direct.py tests/unit/test_detectors_counts.py tests/unit/test_model_hypotheses.py tests/unit/test_model_expected.py
```

On the Linux VM with 4 cores the five files take 21 s and give Mulder about 900 lines of sight for the direct flux. Not run on macOS or Windows yet.
