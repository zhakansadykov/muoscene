# Headline numbers and their conditions

The numbers that other documents quote, each with the conditions it holds under and the place it comes from. Other documents refer here instead of repeating them. Frozen cases are those of `tests/regression/data/`; commands are those of `docs/notes/spike/mulder_spike.py`. Last update: 16 November 2026, Mulder 0.3.7, Calzone 1.2.2.

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

## Accuracy of the base

| Comparison | Size |
|---|---|
| Mulder CSDA ranges against the PDG tables, 0.3 GeV – 1 TeV | within 0.34 % for standard rock, 0.24 % for concrete, 0.52 % for water; Mulder is the shorter |
| other models of pair production, bremsstrahlung, photonuclear loss in Mulder | thresholds change by 8·10⁻⁵, 1.6·10⁻⁴, 1.3·10⁻³ |
| fast path against direct Mulder flux, invented pier, above 1 GeV | rms 0.17 %, maximum 0.39 % |
| a + bE of the preliminary calculation against loam | threshold 16 % too high at 500–1000 g/cm², 9 % at 3600, 6 % at 5400 |
| platforms: Linux arm64, Linux x86-64, macOS arm64 | 3·10⁻⁹, and 3.1·10⁻⁶ in the elastic-scattering tables on macOS |
