# Mulder and Calzone: spike results

Roadmap task 0.1-02. Mulder 0.3.7 and Calzone 1.2.2 (V. Niess, LGPL-3.0-or-later), 12 November 2026. All numbers are reproduced by `docs/notes/spike/mulder_spike.py`; the scene is synthetic. "Preliminary calculation" is the one used for the 2027 grant application: Guan et al. (2015) spectrum at sea level, energy loss a + bE, straight paths, cone significance.

Status: complete. Decision of 13 November 2026 at the end; corrections of 15 and 16 November above the practical notes. The numbers that other documents quote are collected, with their conditions, in `docs/numbers.md`. The numbers are frozen as regression tests in `tests/regression/data/mulder_reference.json` (`docs/legacy.md`).

## Summary

| # | Question | Result |
|---|---|---|
| E1 | Installation | `pip install mulder calzone`: wheels for Linux and macOS, x86-64 and arm64, Python ≥ 3.8 (abi3); 128 MB with numpy. Geometry needs no Geant4 data download. Physics tables compile in 3–10 s on first use and are cached in `~/.cache/mulder`. |
| E2 | Reference flux | `Reference("GCCLY15")` is the Guan et al. (2015) formula with total energy E = K + m: differential flux agrees to 1e-6. Units m⁻² s⁻¹ sr⁻¹ GeV⁻¹, differential in kinetic energy at the detector; integral flux is obtained by integrating over energy. Vertical, sea level: 60.50 m⁻² s⁻¹ sr⁻¹ above 1 GeV, 88.42 in total. |
| E2 | Altitude | The reference is defined at sea level; `Fluxmeter` transports through the atmosphere, energy loss and decay included. At 1150 m the open-sky flux is +19.8 % above 1 GeV, +9.6 % above 5 GeV, +6.3 % above 10 GeV (total energy), +20.0 % over all energies from the cut-off. Under 20 m of soil it is +7.6 %, so the through-mass to open-sky ratio shifts by about 10 % against a sea-level model. Below the altitude cut-off the flux is not defined (E13). |
| E3 | Energy threshold | a + bE of the preliminary calculation is too high: against loam by 16 % at 500–1000 g/cm², 9 % at 3600, 6 % at 5400, 2.5 % at 10⁴ g/cm². Standard rock needs 1.4–2.5 % less energy than loam, concrete 0.3–0.7 % less (500–3·10⁴ g/cm²). |
| E4 | Uniform slab | Integral flux under 10–30 m of soil is 3–12 % above the preliminary calculation at sea level and 9–25 % above it at 1150 m for zenith angles 0–60° (`spike/slab_ratio.png`). Discrete against continuous mode: −0.5 to −0.9 ± 0.2 %. |
| E5 | Wall with a niche | Calzone geometry: ground box or closed OBJ mesh, niche by `subtract`. Path in the ground agrees with the preliminary calculation to 1e-12 m. With a rising brow 603 of 40 401 directions differ by up to 0.08 m: the preliminary calculation takes the brow as flat for directions towards the canyon and along the wall; Mulder is right. |
| E6 | Hypotheses | A sphere placed as a daughter volume is a separate medium; its density is changed without rebuilding the geometry. |
| E6 | Significance | 1 m void 3 m from the wall and 3 m below the brow, niche at 10 m, 240 cm², 60 days: preliminary 10.33σ (the application quotes 10.3σ); Mulder continuous at sea level 10.73σ; at 1150 m 12.02σ (12.05σ with 20 m instead of 300 m of scene air above the brow, `tests/regression`; see E12). All of these count every muon down to 1 MeV at the detector, as the preliminary calculation does; above 0.2 GeV at the detector the last number is 11.52σ, above 1 GeV it is 9.37σ. |
| E6 | Scattering | Discrete mode, muons above 0.1 GeV at the detector, 50 000 states: the excess inside the geometric cone of the void (3.8°) is +4.87 ± 0.22 % against +6.39 % in continuous mode with seed 11, and +4.72 ± 0.13 % over the seeds 11, 12 and 13 — a share of 0.74 ± 0.02; in a cone 1.7 times wider the share is 0.93 ± 0.06. Scattering moves a quarter of the signal out of the geometric cone: 12.05σ at 1150 m becomes 8.9 ± 0.2σ. Above 1 GeV at the detector the share is 0.92 ± 0.02 and 9.37σ becomes 8.6 ± 0.2σ: a threshold removes the muons that scatter most (`mulder_spike.py --events 50000 --seed N [--threshold 1] scattering`). |
| E7 | Speed (Linux VM, 4 cores, arm64) | `scan` (opacity per medium): 6.0e5 directions/s with 14 media. Continuous flux: 4.5e4–1.7e5 states/s, i.e. 640–2400 directions/s with 71 energies. Discrete flux: 1.2e3 states/s above 0.1 GeV at the detector, 3.6e3 above 1 GeV, 2.8e3 at 10 GeV. |
| E7 | Speed (MacBook Air, macOS 27.3, arm64, Python 3.14) | `scan`: 6.5e4 directions/s (single run, first call included). Continuous flux: 2.7e4 states/s, i.e. 375 directions/s with 71 energies. Discrete flux: 2.4e4 states/s at 10 GeV. Measured on 12 November with the earlier pier scene, also of 14 media. Opacity is fast enough on both machines; the energy-integrated flux, 375–2400 directions/s, is 4 to 27 times below the 1e4 directions/s needed for planning. |
| E8 | DEM terrain with a local geometry | `Fluxmeter` takes either layers or a local geometry, not both. Above a local geometry the atmosphere is still traversed: open sky at 45° seen from a niche 10 m below the brow is 50.59 against 50.61 m⁻² s⁻¹ sr⁻¹ for layers only (`mulder_spike.py local`). Calzone converts a DEM to a closed mesh, so the whole scene can be local. |
| E9 | Season | Winter against summer mid-latitude atmosphere at 1150 m, all energies: open sky −0.36 %, under 20 m +0.41 %, their ratio +0.78 %; above 1 GeV open sky +0.78 %, under 20 m +0.38 %. The cut-off rises from 0.279 to 0.296 GeV. Only the air column between the site and sea level changes; the reference flux does not, so this is not a full temperature-effect model. |
| E9 | Geomagnetic field | IGRF by date: east against west at 45° elevation +0.76 % in the open-sky integral flux over all energies, +0.27 % above 1 GeV. This is the bending of muons in the air between the site and the sea-level reference, not the east–west effect of the primaries. |
| E10 | Camera | A camera inside the niche sees only the niche walls. The opacity map from `scan` is the view of the detector (`spike/pier_opacity.png`: the pile cap and nine piles of an invented pier from a niche 20 m below the brow). |
| E11 | Dependency | Mulder: 11 releases since March 2023, eight in December 2025 – June 2026, version 0.x, developed by V. Niess. Calzone: 15 releases since December 2024, last in September 2026. LGPL-3: usable as an imported dependency of BSD-3 code. |
| E12 | Fast path | Opacity from `scan` with equivalent-opacity weights (concrete 0.995 relative to loam) and a table of the integral flux above 1 GeV built by Mulder under a dense slab at the detector altitude (36 elevations × 71 opacities, 2 s), against direct Mulder flux for 3000 random directions through the pier scene: rms 0.17 %, maximum 0.39 %; through concrete rms 0.04 %. Open sky is offset by +0.36 % because the air of the scene box has sea-level density. Fast path 8e6 directions/s, direct 2900 directions/s. Both sides are sums over a grid that starts at 1 GeV and gives its first node a whole step, which overestimates the integral above 1 GeV by 2–3 %, equally: the comparison holds, the absolute values of the table do not. |
| E13 | Altitude cut-off | With the sea-level reference the flux at altitude is zero below K_cut, the energy lost in the air down to sea level: 0.28, 0.38, 0.53 GeV at 1150 m and 0.69, 0.97, 1.40 GeV at 3340 m for elevations of 90°, 45°, 30°. Under soil thinner than about 130 g/cm² the cut-off survives at the detector: 0.12 GeV under 0.5 m at 1150 m. The differential flux has an integrable singularity at the cut-off: a coarse energy grid across it is wrong by up to a factor of two. Integrate from above the cut-off, or use a reference at the site altitude. |

## Corrections (15 November 2026)

Found while freezing these numbers as regression tests (roadmap 0.1-03, `tests/regression/mulder_cases.py`). Above sea level the differential flux of the sea-level reference has an integrable singularity at the altitude cut-off (E13). The spike integrated "all energies" on a grid in ln K across it, 100 nodes per decade, which is wrong by up to 1.3 %. The integrals now run in ln(K − K_cut). Changed: E2 total at sea level 88.38 → 88.42 (integral from K = 0 instead of 1 MeV); E8 49.97 / 50.08 → 50.59 / 50.61 (E8 is now in the script); E9 season −1.45 % / +1.9 % → −0.36 % / +0.78 %; E9 geomagnetic field +1.0 % → +0.76 %. Integrals above 1 GeV and everything under soil are unaffected. E3 and E4 ranges are restated for the full set of frozen cases.

## Corrections (16 November 2026)

Independent review of the spike, of the frozen numbers and of the documents.

- Detector threshold. The significances of E6 count every muon down to 1 MeV at the detector and did not say so. The same example above 0.2 and 1 GeV is added and frozen (`reference_void/.../K>0.2GeV` and `K>1GeV`): 11.52σ and 9.37σ at 1150 m.
- Scattering. The estimate is now a mean over three seeds. The earlier "12.0σ to about 9.2σ" came from one seed and carried no error: it becomes 8.9 ± 0.2σ without a threshold and 8.6 ± 0.2σ above 1 GeV.
- Pier. The pier of E7, E10 and E12 is replaced by an invented one (`PIER` in the script). The figure and the numbers of E7 on Linux and of E12 are recomputed: rms 0.19 % → 0.17 %, maximum 0.40 % → 0.39 %.
- Decay. Mulder applies muon decay along the air path between the site and sea level: the flux at 1150 m is the reference at the shifted energy times the ratio of stopping powers times exp(Δτ/τ), 1.21 at 1 GeV and 1.04 at 5 GeV. The altitude gains of E2 contain it.
- Cut-off at 45° and 1150 m: 0.38 GeV (0.3848), not 0.39; under half a metre of soil a cut-off of 0.12 GeV survives at the detector (frozen as `thin_slab`).
- E7: the energy-integrated flux is 4 to 27 times below the speed needed, not 30–40 times.
- E12: the sums overestimate the integral above 1 GeV by 2–3 % on both sides of the comparison.
- `spike/slab_ratio.png` is not touched by the integration error corrected on 15 November: under 2 m of soil and more the two integrators agree to 0.02 %.

## Practical notes

- Calzone lengths are centimetres, Mulder uses metres and kg/m³, muoscene conventions use metres and g/cm³.
- Positions of daughter volumes are relative to the centre of the mother volume.
- OBJ meshes need face normals (`vn`); `units` is set per mesh.
- Calzone media carry Geant4 material names; they are replaced through `LocalGeometry.media[i].material` and `.density` by Mulder materials defined with `mulder.materials.Material.define` and `Composite.define`.
- The default water of Mulder has 1020 kg/m³.
- Mulder energies are kinetic.
- A detector exactly on a boundary is ambiguous; the spike puts it 1 mm above the niche floor.

## What the preliminary calculation gets wrong, in numbers

| Simplification | Effect on the 1 m void example |
|---|---|
| a + bE threshold | flux under soil underestimated by 3–12 %; significance +4 % when corrected |
| sea-level spectrum | at 1150 m significance +12 % more; the reference-to-through ratio shifts by about 10 % |
| no scattering | a quarter of the signal leaves the geometric cone (share 0.74 ± 0.02 without a detector threshold) |
| net | 10.3σ becomes 8.4–8.9σ for detector thresholds between none and 1 GeV (`docs/numbers.md`) |

## Run on another machine

```bash
pip install mulder calzone
python docs/notes/spike/mulder_spike.py speed
```

## Decision (13 November 2026)

**Hybrid on top of Mulder and Calzone.** muoscene implements no ray tracer, no range tables and no flux model of its own.

| Part | From Mulder and Calzone | Written in muoscene |
|---|---|---|
| geometry | solids, meshes, subtraction, navigation, `scan` | description of engineering scenes, builder to a Calzone tree, checks, `PathTable` |
| materials | compositions, muon ranges | library of soils and concrete, moisture, equivalent-opacity weights |
| flux | reference flux, atmosphere, geomagnetic field, continuous transport | `FluxTable` built by Mulder for each detector; legacy table for regression |
| transport | direct continuous and discrete flux on the scene | fast path: equivalent opacity and table lookup; scattering kernel |
| the rest | — | detectors, forward model, planning and statistics, synthetic data, comparison with data, inversion, visualisation |

Reasons:

1. Paths agree with the preliminary calculation to 1e-12 m (E5): an own tracer would duplicate Calzone.
2. Direct flux with the energy integral runs at 375–2900 directions/s where planning needs 1e4 or more; the fast path gives 8e6 directions/s within 0.2 % rms of it (E12).
3. The physics of the preliminary calculation is replaced by a maintained code: energy loss for real compositions (E3), altitude (E2), scattering (E6).
4. The cost is two binary dependencies with a 0.x interface from one developer.

How the cost is contained: backend imports are confined to `geometry`, `materials`, `flux`, `transport` and `site` and checked by a test; versions are capped; the numbers of this spike are frozen as regression tests, so a release that changes the physics is noticed; every run compares the fast path with the direct one on a sample of directions.

Consequences:

- The energy variable is the kinetic energy at the detector, and every detector has an energy response.
- With the sea-level reference, comparisons use muons above the altitude cut-off. The soft part, which matters for the open-sky reference and for the track-density budget, needs a reference flux at the site altitude (tabulated from a model or measured); moved from 0.3 to 0.1 on 15 November (roadmap 0.1-05b, after the release of 0.1.0).
- Version 0.1 uses continuous transport; scattering is quantified with the discrete mode and applied in 0.3.
- The whole scene is a local geometry; terrain from a DEM enters as a Calzone mesh.
- The air of the scene box takes the density of the atmosphere at the site altitude.
- The paper claims the scene description, detector responses, planning and comparison with data; Mulder, PUMAS and Calzone are cited as the base.
