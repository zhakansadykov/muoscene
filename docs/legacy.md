# The preliminary calculation

The sensitivity numbers of the 2027 grant application came from three short Python scripts written in October 2026: `sens.py` (spectrum, threshold, path in the wall, cone significance), `scheme.py` (box anomaly, track density) and `tilt.py` (tilted films). They are not part of the repository: their parameters are those of a real site. Their physics lives on in two places.

- `tests/regression/legacy.py` — the same arithmetic as pure functions with explicit parameters, frozen on a synthetic canyon (`data/legacy_reference.json`). It is the "legacy mode" of the conventions: a reference for parity tests, not a model offered to users.
- `tests/regression/mulder_cases.py` — the same questions answered by Mulder and Calzone (`data/mulder_reference.json`), so that every simplification below has a measured size.

## Physics of the preliminary calculation

| Part | Model |
|---|---|
| Spectrum | Guan et al. (2015) at sea level, total energy E; integral intensity above E tabulated on 301 values of cos θ from 0.25 (nearest node, clipped below) and 2001 values of ln E from 0.1 GeV (linear interpolation of ln I) |
| Energy threshold | E(X) = a/b (exp(bX) − 1) + 0.1 GeV, a = 2.2·10⁻³ GeV cm²/g, b = 4·10⁻⁶ cm²/g, for opacity X in g/cm² |
| Geometry | vertical wall at x = 0, ground x > 0; brow flat or rising linearly by `rise` over `xr`; niche 0.5 m deep, 0.3 m wide, 0.3 m high, detector on its floor 0.25 m from the face; one soil density (1.8 g/cm³); anomalies as density offsets on sphere or box chords |
| Transport | straight lines, no scattering |
| Detector | horizontal or tilted films of given area, constant efficiency 0.9, angular cut tan θ ≤ 1 |
| Significance | counts in the best of five cones (0.6–1.7 times the angular radius of the target), Monte Carlo directions with fixed seeds, \|N₁ − N₀\| / √N₀ / √(1 + 1/8): the background is estimated from eight neighbouring cones |

## Port and its verification

| Original | Port |
|---|---|
| `sens._tab`, `sens.I_int` | `LN_INTEGRAL`, `integral_intensity` |
| `sens.Eth_of_X` | `threshold` |
| `sens.soil_length`, `sens.chord`, `scheme.box_chord` | `soil_length`, `sphere_chord`, `box_chord` |
| `sens.cone_dirs`, `sens.counts`, `sens.evaluate` | `cone_dirs`, `counts`, `evaluate_sphere` |
| `scheme.eval_box`, `scheme.track_density` | `evaluate_box`, `track_density` |
| `tilt.eval_sphere`, `tilt.eval_box`, `tilt.track_density`, `tilt.ref_flux` | `evaluate_sphere_tilted`, `evaluate_box_tilted`, `track_density_tilted`, `reference_fraction` |
| module-level parameters (`A_PT` and others), changed between calls | fields of `Model` and `Niche` |

On 15 November 2026 every function of the port was compared with the original on the same inputs — random directions and opacities, five anomalies, two detector areas, tilts of 0°, 20° and 45°, flat and rising brows — with `numpy.array_equal`: no differences. The originals had these SHA-256 sums (also recorded in `legacy_reference.json`):

```
sens.py    2f93533b628a8ac34733cce62121cb762d36a723576de9c5b2dc427e9001f18d
scheme.py  e3febef1bee2530d85fb97f87760cce50369f42e66df817b00ab17847ffc4776
tilt.py    ea5d0583b1d90452bec69ab3b90dbe3b669bd7615f7e481d29eb49a475e1abd5
```

## Frozen numbers

```bash
python -m tests.regression.make_legacy_reference --check   # 55 cases, about 2 s
python -m tests.regression.make_mulder_reference --check   # 76 cases, about 25 s (a minute the first time)
pytest tests/regression
```

`legacy_reference.json` is bit-identical on regeneration (fixed seeds). The synthetic canyon of `legacy_cases.py` has a flat bank and a bank whose brow rises by 1.5 m over 8 m, niches 12 m and 22 m below the brow, films of 240 and 600 cm² for 60 days, tilts of 0° and 20°, spheres (void, decompaction, concrete) and a 6 × 6 × 2 m pile cap. One case, `reference_void`, is the textbook example of the application: a 1 m void 3 m from the wall and 3 m under a flat brow, seen from a niche 10 m below it — 10.31σ.

`mulder_reference.json` holds continuous-transport results: the reference flux, open sky at 0 and 1150 m, altitude cut-offs, thresholds and ranges of four materials, the tables behind discrete transport, a soil slab, three atmosphere models, the geomagnetic field, the paths in the wall of the synthetic canyon, a slab thinner than the air below the site, the reference void — without a detector threshold and above 0.2 and 1 GeV — and, for the flux tables, open sky at a site below sea level and the flux above the two thresholds under a soil slab, and, for the geometry, the reference void and the open sky against the height and the density of the air of the local geometry. The tests compare a fresh computation with it at a relative tolerance of 10⁻⁵. The same numbers computed on Linux arm64 and x86-64 agree to 3·10⁻⁹; on macOS arm64 (macOS 27.3, Python 3.14) the elastic-scattering tables differ by up to 3.1·10⁻⁶, all other numbers by at most 4·10⁻⁹. Switching the pair-production, bremsstrahlung or photonuclear model of Mulder changes the thresholds by 8·10⁻⁵, 1.6·10⁻⁴ and 1.3·10⁻³ and fails the tests. Regenerate only for an understood change of Mulder or Calzone and record it in `CHANGELOG.md`.

## Simplifications and their size

Sizes from `mulder_reference.json` and `docs/notes/mulder.md`; K is the kinetic energy, E = K + m.

| Simplification | Size |
|---|---|
| a + bE threshold instead of the CSDA range of a real soil | threshold too high by 16 % at 500–1000 g/cm², 13 % at 1800, 9 % at 3600, 6 % at 5400, 2.5 % at 10⁴; too low beyond 2·10⁴ (−1 to −3 %). Flux under 10–30 m of soil underestimated by 3–12 % at sea level |
| sea-level spectrum at a site at 1150 m | open sky +19.8 % above 1 GeV, +20.0 % over all energies; under 20 m of soil +7.6 %, so the through-mass to open-sky ratio shifts by −10 %. The reference void, every muon counted: 10.73σ at sea level, 12.05σ at 1150 m |
| straight paths, no scattering | a quarter of the signal of the reference void leaves its geometric cone (discrete mode, share 0.74 ± 0.02): 12.05σ becomes 8.9 ± 0.2σ |
| one material: concrete as soil of density 2.4 | composition changes the threshold by 0.3–0.7 % (concrete against loam); equivalent-opacity weight 0.995 |
| brow taken as flat towards the canyon and along the wall | path short by up to rise · x_d / x_r / cos θ under a rising brow: 0.15 m in the synthetic canyon |
| threshold at zero opacity 0.1 GeV, below the muon mass | open-sky intensity +0.27 % |
| nearest cos θ node of the table | integral intensity within 0.95 % (0.4 % for cos θ above 0.7) |
| no altitude cut-off, no atmosphere model, no geomagnetic field | season (winter against summer) and east against west below 1 % (E9 of the spike note) |
| constant efficiency, no energy response, no angular acceptance | a detector threshold of 1 GeV alone lowers the reference void from 12.05σ to 9.37σ at 1150 m (0.2 GeV: 11.52σ); the response of the emulsion is not known yet; angular acceptance comes in roadmap 0.1-07 and 0.3 |
| Gaussian significance with a background from neighbouring cones, no systematics | replaced by likelihood-ratio significance in 0.1-09 and systematics in 0.3 |

Net effect on the reference void: 10.3σ in the application, 10.7σ with Mulder at sea level, 12.05σ at the altitude of the site — all three counting every muon — and 8.4–8.9σ with scattering for detector thresholds between none and 1 GeV (`docs/numbers.md`). Statistics of one point only.
