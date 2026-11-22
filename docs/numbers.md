# Headline numbers and their conditions

The numbers that other documents quote, each with the conditions it holds under and the place it comes from. Other documents refer here instead of repeating them. Frozen cases are those of `tests/regression/data/`; commands are those of `docs/notes/spike/mulder_spike.py`; tests are those of `tests/unit/test_flux_tables.py`, in the geometry section of `tests/unit/test_geometry_scene.py`, and in the sections on the site and the detectors, on the forward model and on the examples of the files named there. Last update: 22 November 2026, Mulder 0.3.7, Calzone 1.2.2.

## The reference void

A 1 m void 3 m from the wall and 3 m below a flat brow, seen from a niche 10 m below the brow; 240 cm² of horizontal films, 60 days, efficiency 0.9; counts in the cone of the angular radius of the void (3.8°); significance |N₁ − N₀| / √N₀ / √(1 + 1/8). Statistics of one point only: no systematics, no fading, no angular efficiency, no soft component.

Continuous transport, straight paths:

| Model | Altitude | Detector threshold | N₀ | Effect | Significance | Frozen case |
|---|---|---|---|---|---|---|
| preliminary calculation: Guan at sea level, a + bE | 0 m | none (E > 0.1 GeV) | 31 405 | +6.17 % | 10.31σ | `legacy_reference.json`, `reference_void` |
| Mulder | 0 m | none (K > 1 MeV) | 35 089 | +6.08 % | 10.73σ | `reference_void/0m/cone1` |
| Mulder | 0 m | K > 0.2 GeV | 33 604 | +5.96 % | 10.30σ | `…/K>0.2GeV` |
| Mulder | 0 m | K > 1 GeV | 27 439 | +5.43 % | 8.48σ | `…/K>1GeV` |
| Mulder | 1150 m | none (K > 1 MeV) | 39 162 | +6.46 % | 12.05σ | `reference_void/1150m/cone1` |
| Mulder | 1150 m | K > 0.2 GeV | 37 411 | +6.32 % | 11.52σ | `…/K>0.2GeV` |
| Mulder | 1150 m | K > 1 GeV | 30 221 | +5.71 % | 9.37σ | `…/K>1GeV` |

The spike script gives 10.33σ and 12.02σ where this table has 10.31σ and 12.05σ: it puts the directions of the preliminary calculation on a regular grid instead of drawing them, and it has 300 m instead of 20 m of scene air, at sea-level density, above the brow. The grant application quotes the first row as 10.3σ.

Scattering — the share of the excess that stays in the geometric cone, discrete over continuous transport, at 1150 m; means over the seeds 11, 12 and 13 with 50 000 states each (`--events 50000 --seed N [--threshold T] scattering`):

| Muons at the detector | Continuous | Discrete | Share in the cone | Share in a cone 1.7 times wider |
|---|---|---|---|---|
| above 0.1 GeV | +6.39 % | +4.72 ± 0.13 % | 0.74 ± 0.02 | 0.93 ± 0.06 |
| above 0.2 GeV | +6.30 % | +4.62 ± 0.13 % | 0.73 ± 0.02 | 0.99 ± 0.06 |
| above 1 GeV | +5.70 % | +5.23 ± 0.13 % | 0.92 ± 0.02 | 0.99 ± 0.06 |

Significance with scattering, the continuous one times the share:

| Detector threshold | Continuous, 1150 m | With scattering |
|---|---|---|
| none | 12.05σ | 8.9 ± 0.2σ |
| 0.2 GeV | 11.52σ | 8.4 ± 0.2σ |
| 1 GeV | 9.37σ | 8.6 ± 0.2σ |

The error is the Monte Carlo error of the share. The share without a threshold is the one measured above 0.1 GeV; muons below 0.2 GeV are 4.5 % of the counts. A threshold and scattering are not independent losses: the threshold removes the muons that scatter most.

## Open sky and altitude

| Quantity | Value | Conditions | Frozen case |
|---|---|---|---|
| vertical intensity at sea level | 60.50 m⁻² s⁻¹ sr⁻¹ | total energy above 1 GeV, GCCLY15 | `open_sky/0m/elev90` |
| the same over all energies | 88.42 m⁻² s⁻¹ sr⁻¹ | from K = 0 | `open_sky/0m/elev90` |
| gain at 1150 m, vertical | +19.8 %, +9.6 %, +6.3 % | above 1, 5, 10 GeV of total energy; decay along the air path included | `open_sky/1150m/elev90` |
| the same over all energies | +20.0 % | above the cut-off of 0.28 GeV | `open_sky/1150m/elev90` |
| gain under 20 m of soil | +7.6 % | all energies at the detector | `slab/…/20m/zenith0` |
| cut-off of the sea-level reference at 1150 m | 0.28, 0.38, 0.53 GeV | elevation 90°, 45°, 30°; 0.62 GeV at 25°, 1.0 GeV at 15° | `cutoff/1150m/…` |
| the same at 3340 m | 0.69, 0.97, 1.40 GeV | elevation 90°, 45°, 30° | `cutoff/3340m/…` |
| cut-off left under 0.5 m of soil at 1150 m | 0.12, 0.16 GeV | zenith 0°, 45° | `thin_slab/1150m/0.5m/…` |

## Flux tables

`MulderTable` of a detector at 1150 m with the loam of the spike as reference material, unless stated (`docs/notes/flux.md`). Continuous transport; thresholds are kinetic energies at the detector.

| Quantity | Value | Conditions | Test or frozen case |
|---|---|---|---|
| vertical intensity 430 m below sea level | 56.69 m⁻² s⁻¹ sr⁻¹ | total energy above 1 GeV, against 60.50 at sea level; no cut-off | `below_sea_level/-430m/elev90` |
| vertical flux under soil of 1.8 g/cm³ | 96.70 and 64.50 under 0.5 m, 74.33 and 52.92 under 2 m, 13.02 and 11.09 m⁻² s⁻¹ sr⁻¹ under 20 m | above 0.2 and 1 GeV; Mulder with the real soil | `threshold_slab/1150m/…/zenith0` |
| table against Mulder under a dense slab, differential flux | rms 0.0035 %, largest 0.045 % | 16 150 random energies, elevations and opacities, up to 1 TeV and from 10 MeV above the cut-off | `test_table_is_within_02_percent_of_mulder_at_random_points` |
| the same above 1 TeV | rms 0.044 %, largest 0.79 % | 3850 points; under 2.5·10⁴ g/cm² the flux of Mulder wavers by 0.34 % around 2 TeV, the table by 0.09 % | the same and `test_mulder_ripples_above_1_tev_under_a_thick_slab` |
| table against Mulder under a dense slab, integral flux | within 6·10⁻⁵ | above 0.2 and 1 GeV, 12 random directions each | `test_integrals_are_those_of_mulder_under_a_slab` |
| table at the altitude of the detector against real soil above it | −0.05 to −0.10 % | every muon under 10–30 m at zenith angles of 0–60°: the table has air in place of that soil; a table at the altitude of the surface is +0.03 to +0.05 %, it takes no decay in the soil | `test_slab_numbers_of_the_spike` |
| the same above 0.2 and 1 GeV | −0.007 to −0.013 % under 0.5 m, −0.018 % under 2 m, −0.074 to −0.077 % under 20 m | zenith 0° and 45° | `test_detector_thresholds_under_soil` |
| a sum that gives its first node a whole step | +2.85 % at sea level, +3.29 % at 1150 m | vertical, 61 nodes from 1 GeV to 10 TeV, against trapezoids from the threshold | `test_a_sum_from_the_first_node_is_two_to_three_percent_high` |
| where a threshold of 1 GeV holds in the open sky | from 15.7° of elevation (sine 0.27) | the cut-off is 0.98 GeV there and 1.02 GeV at 15.1°; a table starts 10 MeV above its cut-off | `test_the_table_holds_only_energies_above_the_cutoff` |
| opacity under which no cut-off is left | 138 g/cm² by the vertical, 569 g/cm² at 14.5° | loam of the spike | the same |
| muons above the top of a table, 10 TeV at the detector | 3·10⁻⁶ under 3·10⁴ g/cm², 2·10⁻⁸ in the open sky | of the muons above 1 GeV, vertical | `test_muons_above_the_table_do_not_matter` |

## Geometry

`muoscene.geometry` (`docs/notes/geometry.md`). Paths are lengths in metres along straight lines of sight; the frozen cases `scene_air/1150m/…` are continuous transport above 1 GeV at the detector.

| Quantity | Value | Conditions | Test or frozen case |
|---|---|---|---|
| path in the ground behind the wall with a niche, a scene of muoscene against the frozen one | within 1.1·10⁻¹³ m | flat and rising brow, niches 12 and 22 m below it, the 175 directions of `legacy_cases` | `test_wall_paths_equal_the_frozen_ones` |
| layers against the whole ground | sum within 1.1·10⁻¹³ m | one to three layer tops, also through the niches | `test_layers_add_up_to_the_ground` |
| chords of a box, of cylinders and of a sphere against analytic ones | within 1·10⁻¹⁴, 3·10⁻¹² and 4·10⁻¹³ m | cap, nine piles and void of the invented pier, 40 401 directions; with two layer tops every pile is in three parts | `test_chords_of_sphere_cylinder_and_box` |
| scan speed | 4.0·10⁵ and 2.7·10⁵ directions/s on Linux, 5.9·10⁴ and 5.4·10⁴ on macOS | the invented pier, 14 and 34 media, checks of `Scene.scan` included; a Linux VM with 4 cores and macOS itself (Python 3.14) on the same MacBook Air, arm64 | `test_scan_speed_on_the_pier` (asserts 5·10⁴; marked `validation`); on macOS the command of `docs/notes/geometry.md`, 17 November 2026 |
| 1000 m of scene air above the brow at Mulder's default density, 1.205 kg/m³ | counts of the reference void −0.91 %, open sky at 45° −1.55 % | 1150 m; against 1 m of air at the density of the site, 1.0997 kg/m³ | `scene_air/1150m/top1000m/default` |
| the same at the density of the site | −0.31 % and −0.55 % | a uniform kilometre is heavier than the real air above the site | `scene_air/1150m/top1000m/site` |
| 20 m of air at the default density, the scene of the frozen reference void | −0.012 % and −0.030 % | | `scene_air/1150m/top20m/default` |
| the air of a scene: 1 m above the highest body, at the density of the site | open sky within 2.6·10⁻⁶ of Mulder without a local geometry | 2.5·10⁻⁷ with the air of the library | `scene_air/1150m/atmosphere`, `test_open_sky_through_the_scene_air_is_that_of_the_atmosphere` |
| a void given to Mulder as a density of zero | flux 2.1 % below that through soil, instead of 8.7 % above | along the axis of the reference void; Mulder takes the density of the material, so a void gets 10⁻⁹ kg/m³ | `test_a_void_is_empty_for_mulder_too` |

## Site and detectors

`muoscene.site` and `muoscene.detectors` (`docs/notes/detectors.md`). Tests are those of `tests/unit/test_site_frames.py`, `test_detectors_pose.py`, `test_detectors_grid.py`, `test_detectors_emulsion.py` and `tests/regression/test_legacy_track_density.py`. Track densities are muon tracks per cm² of film in 60 days, continuous transport, straight paths, both sides of the film; they hold no soft component and no muons below the threshold.

| Quantity | Value | Conditions | Test or frozen case |
|---|---|---|---|
| sum of the bin solid angles against 2π (1 − cos θ_max) | within 2·10⁻¹⁵ | default grid: slopes up to 1, 40 × 40 bins of 0.05, 148 of them cut by the acceptance circle; exact part of every bin within the circle | `test_the_default_grid_adds_up_to_the_cone` |
| the same with edge bins taken by the share of sub-samples inside the circle | +3.4·10⁻³, −3.5·10⁻⁴, +4.3·10⁻⁵, +2.0·10⁻⁴ | 1, 2, 4, 8 sub-samples per side; the alternative that was not taken | `test_the_share_of_sub_samples_inside_the_circle_is_not_exact` |
| pose round trips | 6·10⁻¹⁴ degrees; 3·10⁻¹⁶ | angles → axes → angles for tilts from 0.01°; direction → slopes → direction | `test_angles_to_axes_to_angles`, `test_direction_to_slopes_to_direction` |
| east against west with the frame of a scene turned by 90° | +0.25 % | open sky above 1 GeV at 45°, 1150 m, field of 15 October 2028; a scene that kept its y axis for north would give north and south | `test_a_scene_of_a_site_tells_mulder_where_north_is` |
| track density of the preliminary calculation, grid sum against the frozen Monte Carlo number | +0.01 to +0.03 % horizontal films, −0.06 to −0.10 % tilted | sea level, every muon above a + bE (`LegacyTable`), niches 12 and 22 m, flat and rising brow, tilts of 0° and 20°; Monte Carlo σ of the frozen numbers 0.28–0.39 % and 0.21–0.22 % (seeds 100–139) | `test_horizontal_films`, `test_tilted_films`, `test_open_sky`; `legacy_reference.json`, `track_density/…`, `track_density_tilted/…` |
| the same on the directions of the preliminary calculation themselves | 4·10⁻¹⁶; −3.0·10⁻⁶ and −3.5·10⁻⁶ under the rising brow through a scene | 400 000 directions, seed 5; the preliminary calculation takes the rising brow as flat towards the canyon and along the wall, up to 0.33 m of path | `test_the_same_integrand_on_the_directions_of_the_preliminary_calculation` |
| the back side of a tilted film | 0.48 % of the track density at a tilt of 20° | conditions of the preliminary calculation, niche 22 m; lines of sight within tan 10 of the normal, cos θ above 0.05 | `test_the_back_side_of_a_tilted_film` |
| track density in the open sky above 1 GeV | 60 039 at 1150 m, 49 843 at sea level | horizontal film, elevations above 16°; kinetic energy at the detector | `test_open_sky_below_the_cutoff_is_refused` |
| the same above 0.2 GeV | 64 431 at sea level; refused at 1150 m | the sea-level reference has a cut-off of 0.28 GeV by the vertical at 1150 m | the same |
| track density in a niche above 1 GeV | 34 758 at 12 m, 30 346 at 22 m below the brow | horizontal film on the floor of the niche of `legacy_cases`, flat bank of loam of 1.8 g/cm³ at 1150 m, elevations above 16°; table at the altitude of the detector; 2.2·10⁵ lines of sight, converged to 3·10⁻⁵ | `test_track_density_in_a_niche_at_1150_m` |
| of it, from the side of the canyon | 26 244 (76 %) and 26 104 (86 %) | open sky but for the lip of the niche | the same |
| track density in a niche above 0.2 GeV, side of the ground only | 10 081 at 12 m, 4794 at 22 m; 8514 and 4242 above 1 GeV | the side of the canyon has no flux above 0.2 GeV with a sea-level reference: the whole is refused until roadmap 0.1-05b | the same |
| table and scene against Mulder through the scene, track density | −0.04 and −0.06 % on the side of the ground, +0.03 % on the side of the canyon | the 1152 lines of sight of the frozen cases, above 0.2 and 1 GeV | `test_track_density_in_a_niche_against_mulder`; `track_density/1150m/12m`, `…/22m` |

## Forward model

`muoscene.transport` and `muoscene.model` (`docs/notes/forward.md`). Tests are those of `tests/unit/test_transport_table.py`, `test_transport_direct.py` and `test_model_expected.py`. Continuous transport, straight lines of sight; thresholds are kinetic energies at the detector; the scene is the invented pier of the spike at 1150 m, seen from its niche 20 m below the brow, in loam of 1.8 g/cm³ with concrete of 2.3 g/cm³.

| Quantity | Value | Conditions | Test or frozen case |
|---|---|---|---|
| fast path against the direct flux of Mulder, in the run | mean +0.006 %, rms 0.014 %, largest 0.032 % | above 1 GeV; 256 random accepted lines of sight, seed 8; table at the altitude of the detector, equivalent opacity with the displaced air taken out | `test_the_demo_is_within_half_a_percent_of_mulder` |
| the same above 0.2 GeV | mean +0.020 %, rms 0.022 %, largest 0.049 % | on the lines of sight that have a flux: under soil that leaves no cut-off | the same |
| fast path against Mulder through a pier built without muoscene | within 0.033 % above 1 GeV, 0.050 % above 0.2 GeV | 77 lines on a grid of slopes; above 0.2 GeV the 34 into the massif | `test_the_fast_path_against_the_frozen_pier`; `pier/1150m/soil`, `…/void` |
| uniform slab, table against Mulder with the real soil | −0.003 to +0.021 %; −0.007 to −0.078 % without the displaced air | 0.5, 2 and 20 m of soil, zenith 0° and 45°, above 0.2 and 1 GeV | `test_a_uniform_slab_matches_mulder`; `threshold_slab/…` |
| a table at the altitude of the surface, open sky seen 20 m below it | +0.31 to +0.39 % | above 1 GeV, sines of elevation 1 to 0.5 | `test_a_table_at_the_altitude_of_the_surface_is_wrong_for_the_open_sky` |
| weights of the equivalent opacity, the iteration taken | +0.023 % in flux above 1 GeV, +0.034 % above 0.2 GeV | 3000 g/cm² of loam under 1840 g/cm² of concrete, against the chain of ranges in the true order | `test_the_alternatives_that_were_not_taken` |
| the same with a weight at one energy | +0.17, +0.06, −0.08 % | weight taken at 1, 7, 100 GeV; −0.17 % with no weights; +0.048 % for the chain in the wrong order | the same |
| a hypothesis without a new scan against a rebuilt scene | counts within 10⁻¹² | no void, a decompaction of 20 %, a bank of 1.7 g/cm³ | `test_a_hypothesis_without_a_new_scan_equals_a_rebuilt_scene` |
| expected tracks in the map of the pier | 2.398·10⁶ above 1 GeV; the void adds 572, 10.9 % in its two bins | 120 cm² of horizontal film, 60 days, efficiency 0.9, cone of 45°, 1324 bins of 0.05 in the slopes | `test_counts_are_the_sum_of_the_detector_over_its_lines`, `test_the_void_of_the_frozen_pier_is_seen` |
| bins with a flux above 0.2 GeV | 968 of 1324; the open sky in none | sea-level reference at 1150 m: the canyon and the thin soil above the niche are below the cut-off; all 1324 above 1 GeV | `test_bins_without_a_flux_are_marked_not_filled` |
| a film tilted by 30° towards the canyon, above 1 GeV | 6 of 1324 bins without a flux | lines of sight between 15.1° and 15.7° of elevation | `test_a_tilted_detector_and_the_edge_of_the_table` |
| the check of the fast path where 1 GeV is just above the cut-off | the table within 0.021 % of the direct flux; it was 0.62 % below the direct flux of 30 energy nodes a decade, which refused right runs | open sky at 1150.5 m, 15.5° of elevation; 240 nodes a decade over the first half decade from the threshold since 22 November 2026 (roadmap 0.1-10) | `tests/unit/test_transport_direct.py`, `test_the_check_at_the_edge_of_the_cut_off` |

## Examples

`examples/canyon_demo`, `examples/open_sky` and `examples/slab` (roadmap 0.1-10), all at 1150 m with the sea-level reference, continuous transport; films of 120 cm², efficiency 0.9, a cone of 45° in 1324 bins of 0.05 in the slopes; loam with 15 % of water. The sizes of the canyon are in `examples/canyon_demo/README.md`. Tests are those of `tests/unit/test_example_canyon_demo.py`, `test_example_open_sky.py` and `test_example_slab.py`; the run as a user starts it is `python run.py` in each folder.

| Quantity | Value | Conditions | Test |
|---|---|---|---|
| tracks above 1 GeV, synthetic canyon | 1.771·10⁶ and 2.293·10⁶ in the niches 12 and 22 m below the left brow; 2.147·10⁶ and 2.751·10⁶ on the wall beside them | 60 days; the upper films tilted by 20° into the massif, the lower horizontal; a niche 0.5 × 0.3 × 0.3 m, a film on the wall 0.15 m in front of it | `test_the_anomalies_are_seen_where_the_scene_puts_them`, `test_a_niche_against_pins` |
| the void of 2 m behind the left wall, above 1 GeV | +1947 and +360 tracks in the niches, up to +12.6 and +10.6 % in a bin; +1614 and +345 on the wall, up to +12.1 and +10.8 % | the same films | `test_the_anomalies_are_seen_where_the_scene_puts_them` |
| open-sky bins, niche against wall | 0 and 52 against 382 and 662 | the same films | `test_a_niche_against_pins` |
| bins with a number above 0.2 GeV | 1292 and 990 in the niches, 942 and 622 on the wall, of 1324 | the same films; the others are below the cut-off of the reference | `test_above_0_2_gev_the_open_sky_has_no_number` |
| the decompaction of 3 m behind the right wall, −20 %, above 1 GeV | +1189 tracks, up to +3.6 % in a bin | a film tilted by 20° in a niche 12 m below the right edge | `test_the_anomalies_are_seen_where_the_scene_puts_them` |
| fast path against Mulder in the demo | at most 0.066 % above 1 GeV, 0.14 % above 0.2 GeV | 256 lines per film and threshold, seed 10 (`python run.py`); 0.065 and 0.14 % on 80 lines in the test | `test_every_run_checks_itself_against_mulder` |
| open sky: counts against the flux table | within 7·10⁻¹⁶; 69.22 m⁻² s⁻¹ sr⁻¹ above 1 GeV at the vertical; no bin above 0.2 GeV | films 0.5 m above the ground; cut-off 0.284 GeV at the vertical | `test_every_check_holds`, `test_what_the_sea_level_reference_gives_at_1150_m` |
| slab: the flux along the vertical and 45° against the frozen `threshold_slab` | −0.003 to +0.019 % | the loam of the frozen cases at 1.8 g/cm³, 0.5, 2 and 20 m, above 0.2 and 1 GeV | `test_the_frozen_slab` |
| slab: opacity against ρ d / cos θ | within 8·10⁻¹⁶ in every bin | loam of 1.84 g/cm³ | `test_every_check_holds` |
| a film on the left wall tilted into the massif | at 40° 154 bins below 14.5° of elevation, where the table ends; from 44° refused, its lines below 4° leave the model | 22 m below the brow, 0.15 m in front of the wall, a model 300 m behind the wall and along it | `tests/unit/test_geometry_canyon.py`, `test_how_far_a_film_on_the_wall_can_tilt_into_the_massif` |

## Accuracy of the base

| Comparison | Size |
|---|---|
| Mulder CSDA ranges against the PDG tables, 0.3 GeV – 1 TeV | within 0.34 % for standard rock, 0.24 % for concrete, 0.52 % for water; Mulder is the shorter |
| other models of pair production, bremsstrahlung, photonuclear loss in Mulder | thresholds change by 8·10⁻⁵, 1.6·10⁻⁴, 1.3·10⁻³ |
| fast path against direct Mulder flux, invented pier, above 1 GeV | spike, with its coarse table: rms 0.17 %, maximum 0.39 %; muoscene: rms 0.014 %, maximum 0.032 % ("Forward model") |
| a + bE of the preliminary calculation against loam | threshold 16 % too high at 500–1000 g/cm², 9 % at 3600, 6 % at 5400 |
| platforms: Linux arm64, Linux x86-64, macOS arm64 | 3·10⁻⁹, and 3.1·10⁻⁶ in the elastic-scattering tables on macOS |
