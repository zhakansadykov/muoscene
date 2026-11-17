# Headline numbers and their conditions

The numbers that other documents quote, each with the conditions it holds under and the place it comes from. Other documents refer here instead of repeating them. Frozen cases are those of `tests/regression/data/`; commands are those of `docs/notes/spike/mulder_spike.py`; tests are those of `tests/unit/test_flux_tables.py` and, in the geometry section, of `tests/unit/test_geometry_scene.py`. Last update: 17 November 2026, Mulder 0.3.7, Calzone 1.2.2.

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

## Accuracy of the base

| Comparison | Size |
|---|---|
| Mulder CSDA ranges against the PDG tables, 0.3 GeV – 1 TeV | within 0.34 % for standard rock, 0.24 % for concrete, 0.52 % for water; Mulder is the shorter |
| other models of pair production, bremsstrahlung, photonuclear loss in Mulder | thresholds change by 8·10⁻⁵, 1.6·10⁻⁴, 1.3·10⁻³ |
| fast path against direct Mulder flux, invented pier, above 1 GeV | rms 0.17 %, maximum 0.39 % |
| a + bE of the preliminary calculation against loam | threshold 16 % too high at 500–1000 g/cm², 9 % at 3600, 6 % at 5400 |
| platforms: Linux arm64, Linux x86-64, macOS arm64 | 3·10⁻⁹, and 3.1·10⁻⁶ in the elastic-scattering tables on macOS |
