# Mulder and Calzone: spike results

Roadmap task 0.1-02. Mulder 0.3.7 and Calzone 1.2.2 (V. Niess, LGPL-3.0-or-later), 12 November 2026. All numbers are reproduced by `docs/notes/spike/mulder_spike.py`; the scene is synthetic. "Preliminary calculation" is the one used for the 2027 grant application: Guan et al. (2015) spectrum at sea level, energy loss a + bE, straight paths, cone significance.

Status: measurements are complete; the decision is at the end. The decision — build on Mulder, hybrid, or own core — and the rewrite of tasks 03–07 follow in a separate step.

## Summary

| # | Question | Result |
|---|---|---|
| E1 | Installation | `pip install mulder calzone`: wheels for Linux and macOS, x86-64 and arm64, Python ≥ 3.8 (abi3); 128 MB with numpy. Geometry needs no Geant4 data download. Physics tables compile in 3–10 s on first use and are cached in `~/.cache/mulder`. |
| E2 | Reference flux | `Reference("GCCLY15")` is the Guan et al. (2015) formula with total energy E = K + m: differential flux agrees to 1e-6. Units m⁻² s⁻¹ sr⁻¹ GeV⁻¹, differential in kinetic energy at the detector; integral flux is obtained by integrating over energy. Vertical, sea level: 60.50 m⁻² s⁻¹ sr⁻¹ above 1 GeV, 88.38 in total. |
| E2 | Altitude | The reference is defined at sea level; `Fluxmeter` transports through the atmosphere. At 1150 m the open-sky flux is +19.8 % above 1 GeV, +9.6 % above 5 GeV, +6.3 % above 10 GeV (total energy). Under 20 m of soil it is +7.6 %, so the through-mass to open-sky ratio shifts by about 10 % against a sea-level model. Below the altitude cut-off the flux is not defined (E13). |
| E3 | Energy threshold | a + bE of the preliminary calculation is too high: against loam by 16 % at 500–1000 g/cm², 9 % at 3600, 6 % at 5400, 2.5 % at 10⁴ g/cm². Loam, concrete and standard rock differ from each other by 0.5–2 %. |
| E4 | Uniform slab | Integral flux under 10–30 m of soil is 3–12 % above the preliminary calculation at sea level and 11–25 % above it at 1150 m (`spike/slab_ratio.png`). Discrete against continuous mode: −0.5 to −0.9 ± 0.2 %. |
| E5 | Wall with a niche | Calzone geometry: ground box or closed OBJ mesh, niche by `subtract`. Path in the ground agrees with the preliminary calculation to 1e-12 m. With a rising brow 603 of 40 401 directions differ by up to 0.08 m: the preliminary calculation takes the brow as flat for directions towards the canyon; Mulder is right. |
| E6 | Hypotheses | A sphere placed as a daughter volume is a separate medium; its density is changed without rebuilding the geometry. |
| E6 | Significance | 1 m void 3 m from the wall and 3 m below the brow, niche at 10 m, 240 cm², 60 days: preliminary 10.33σ (the application quotes 10.3σ); Mulder continuous at sea level 10.73σ; at 1150 m 12.02σ. |
| E6 | Scattering | Discrete mode: the excess inside the geometric cone of the void (3.8°) is +4.87 ± 0.22 % against +6.39 % in continuous mode (−24 %); in a cone 1.7 times wider it is +2.02 ± 0.21 % against +2.21 %. Scattering moves a quarter of the signal out of the geometric cone; the significance at 1150 m drops from 12.0σ to about 9.2σ. |
| E7 | Speed (Linux VM, 4 cores, arm64) | `scan` (opacity per medium): 5.6e5 directions/s with 14 media. Continuous flux: 4.6e4–1.7e5 states/s, i.e. 650–2400 directions/s with 71 energies. Discrete flux: 1.3e3–2.5e3 states/s. |
| E7 | Speed (MacBook Air, macOS 27.3, arm64, Python 3.14) | `scan`: 6.5e4 directions/s (single run, first call included). Continuous flux: 2.7e4 states/s, i.e. 375 directions/s with 71 energies. Discrete flux: 2.4e4 states/s at 10 GeV. Opacity is fast enough on both machines; the energy-integrated flux is 30–40 times below the 1e4 directions/s needed for planning. |
| E8 | DEM terrain with a local geometry | `Fluxmeter` takes either layers or a local geometry, not both. Above a local geometry the atmosphere is still traversed: open sky seen from the niche is 49.97 against 50.08 m⁻² s⁻¹ sr⁻¹ for layers only. Calzone converts a DEM to a closed mesh, so the whole scene can be local. |
| E9 | Season | Winter against summer mid-latitude atmosphere at 1150 m: open sky above the cut-off −1.45 %, under 20 m +0.41 %, their ratio +1.9 %. Only the air column between the site and sea level changes; the reference flux does not, so this is not a full temperature-effect model. |
| E9 | Geomagnetic field | IGRF by date: east against west at 45° elevation +1.0 % in the open-sky integral flux. |
| E10 | Camera | A camera inside the niche sees only the niche walls. The opacity map from `scan` is the view of the detector (`spike/pier_opacity.png`: pile cap and nine piles from a niche 25 m below the brow). |
| E11 | Dependency | Mulder: 11 releases since March 2023, eight in December 2025 – June 2026, version 0.x, developed by V. Niess. Calzone: 15 releases since December 2024, last in September 2026. LGPL-3: usable as an imported dependency of BSD-3 code. |
| E12 | Fast path | Opacity from `scan` with equivalent-opacity weights (concrete 0.995 relative to loam) and a table of the integral flux above 1 GeV built by Mulder under a dense slab at the detector altitude (36 elevations × 71 opacities, 2 s), against direct Mulder flux for 3000 random directions through the pier scene: rms 0.19 %, maximum 0.40 %; through concrete rms 0.04 %. Open sky is offset by +0.37 % because the air of the scene box has sea-level density. Fast path 8e6 directions/s, direct 2900 directions/s. |
| E13 | Altitude cut-off | With the sea-level reference the flux at altitude is zero below K_cut, the energy lost in the air down to sea level: 0.28, 0.39, 0.53 GeV at 1150 m and 0.69, 0.97, 1.40 GeV at 3340 m for elevations of 90°, 45°, 30°. The differential flux has an integrable singularity at the cut-off: a coarse energy grid across it is wrong by up to a factor of two. Integrate from above the cut-off, or use a reference at the site altitude. |

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
| no scattering | signal inside the geometric cone overestimated by 24 % |
| net | 10.3σ becomes about 9.2σ |

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
- With the sea-level reference, comparisons use muons above the altitude cut-off. The soft part, which matters for the open-sky reference and for the track-density budget, needs a reference flux at the site altitude (tabulated from MCEq or measured); planned for 0.3.
- Version 0.1 uses continuous transport; scattering is quantified with the discrete mode and applied in 0.3.
- The whole scene is a local geometry; terrain from a DEM enters as a Calzone mesh.
- The air of the scene box takes the density of the atmosphere at the site altitude.
- The paper claims the scene description, detector responses, planning and comparison with data; Mulder, PUMAS and Calzone are cited as the base.
