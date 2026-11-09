# Roadmap

| Version | Target | Content |
|---|---|---|
| 0.1 | Dec 2026 | parity with the preliminary calculation on a synthetic canyon; first 3D scene |
| 0.2 | May 2027 | real sites: coordinate systems, meshes, survey import, layers, structures, voxels |
| 0.3 | Aug 2027 | physics and planning: range tables, altitude, season, scattering, emulsion fading and track-density budget, systematics, profile likelihood, sensitivity volume, scheme optimisation |
| 0.4 | Nov 2027 | data: synthetic tracks, track readers, comparison and fits; straw telescope |
| 1.0 | Dec 2027 – Mar 2028 | validation suite, documentation, Zenodo DOI, PyPI, software paper |
| 2.0 | 2029 | inversion on the voxel grid; real exposure data |

One task is one working session and one or a few commits. A task is done when it has code, tests, docstrings, a `CHANGELOG.md` line and its box ticked here.

## 0.1 — parity and first scene

Goal: `muoscene plan examples/canyon_demo/run.yaml` reproduces the preliminary sensitivity calculation in legacy mode on a synthetic canyon and renders the scene. Conventions — `docs/conventions.md`.

- [ ] **0.1-01 Conventions.** Frames, directions, pose, slopes, units, energy threshold, anomalies, expected-map layout, randomness, legacy mode. *Done when:* `docs/conventions.md` is approved.
- [ ] **0.1-02 Prior-art spike: Mulder and Calzone (next).** Mulder 0.3.7 already provides layered DEM terrain, local geometry through Calzone (meshes, boxes, cylinders, spheres, subtraction), ENU frames from latitude and longitude, the GCCLY15 (Guan et al. 2015) reference flux with altitude, seasonal atmospheres, IGRF, ranges for arbitrary compositions and CSDA / Monte Carlo transport. Experiments: E1 install on 3.12 and 3.14 (Geant4 data needed for geometry?); E2 reference flux against our Guan formula; E3 `inverse_range` against the legacy a + bE threshold; E4 uniform slab, continuous versus discrete; E5 legacy wall-with-niche geometry in Calzone, grammage parity with the legacy path to 1e-6; E6 density hypotheses without rebuilding; E7 speed of `scan` and continuous `flux` on the Mac (threshold 1e4 directions/s); E8 DEM terrain together with a local geometry; E9 season and altitude effect on the reference-to-through-mass ratio; E10 camera view from a niche; E11 dependency risk (release cadence, API stability, LGPL-3 with BSD-3). Decision: A — build on Mulder and Calzone; B — hybrid (grammage and ranges from Mulder, own fast flux); C — own core, Mulder for validation only. *Done when:* `docs/notes/mulder.md` with numbers and the decision; tasks 03–07 and `docs/conventions.md` rewritten accordingly.
- [ ] **0.1-03 Legacy reference.** Port the physics of the preliminary scripts (spectrum, integral table, threshold, path in the wall with a niche, sphere and box chords, cone counts and significance, tilted films, track density) to `tests/regression/legacy.py` as pure functions with explicit parameters; generate `tests/regression/data/legacy_reference.json` for a synthetic parameter set (two banks with flat and rising brows, niches at two depths, tilts 0° and 20°, sphere and box anomalies). *Done when:* regeneration is bit-identical with fixed seeds; no parameters of a real site.
- [ ] **0.1-04 Flux.** `FluxModel` protocol (integral and differential intensity); Guan et al. (2015) and Gaisser; integral table in (cos θ, ln E) with interpolation of ln I over cos θ ∈ [0, 1]; legacy option clipping cos θ at 0.25. *Done when:* vertical Guan intensity is 60.5 m⁻² s⁻¹ sr⁻¹ above 1 GeV and 88.4 above m_μ (±0.5 %); legacy parity to 1e-6.
- [ ] **0.1-05 Materials and threshold.** `Material`; library (standard rock, water, concrete, soil, air); `RangeModel` protocol; `LinearAB` (legacy) and `RangeTable` (PDG muon tables for standard rock, water and concrete, log-log interpolation, Z/A scaling for other materials); data files with citations. *Done when:* legacy parity; tables monotonic and round-trip.
- [ ] **0.1-06 Geometry I.** `RayBatch`; per-ray intervals and priority resolution; box, cylinder, sphere, extruded polygon (canyon cross-section), pile grid, horizontal layer inside a parent body; `Scene.trace` → `PathTable` (length per body, opacity). *Done when:* analytic chords; overlap tests; legacy wall-with-niche path parity to 1e-9 m; randomised check against point sampling along rays.
- [ ] **0.1-07 Site and detectors I.** `Site` (ENU, origin metadata; CRS conversion in 0.2); `Pose` per conventions; `DirectionGrid` with exact bin solid angles and sub-sampling; `EmulsionStack` (area, films, cone, constant efficiency, projected area, track density over the whole sky). *Done when:* grid solid angles sum to the analytic cone; pose round trips; legacy track-density parity within Monte Carlo tolerance.
- [ ] **0.1-08 Forward model.** Grid → rays → trace → hypotheses applied to stored paths → threshold → flux → counts; `ExpectedMap`. *Done when:* open sky matches the analytic integral; uniform slab transmission; hypotheses without retracing equal retracing a modified scene.
- [ ] **0.1-09 Significance and exposure.** Legacy cone estimator; binned χ² and Asimov likelihood-ratio significance; days to kσ; track density versus days against the scanning limit; Markdown table in the layout of the preliminary calculation. *Done when:* legacy parity within 2 %; unit tests of the formulas.
- [ ] **0.1-10 Examples.** `examples/canyon_demo` (generator: canyon depth and widths, two banks with different brows, two niche levels, a pier on each bank with a pile cap and 3 × 3 piles, two horizontal layers) with `site.yaml`, `detectors.yaml`, `run.yaml`; `examples/open_sky`, `examples/slab`. *Done when:* the parameters are documented and differ from any real site.
- [ ] **0.1-11 Configuration, CLI, provenance.** pydantic schema for every 0.1 object (`extra="forbid"`, units in field names); loader with `include` and relative paths; `from_config` builders; `muoscene check` and `muoscene plan`; run folder with resolved configuration, input hashes, version, git commit, NetCDF and Markdown report. *Done when:* examples run from the command line; broken configurations give readable errors.
- [ ] **0.1-12 Visualisation I.** Angular maps (counts, opacity, transmission, significance per hypothesis) in one style; PyVista scene (canyon, niches, pile caps, piles, anomalies, detectors with view cones coloured by significance) exported to HTML. *Done when:* the HTML opens in a browser and the scene matches the configuration.
- [ ] **0.1-13 Release 0.1.0.** README quickstart, `docs/legacy.md`, CHANGELOG, roadmap; CI green; tag `v0.1.0` (no Zenodo, no PyPI yet).

Order: 02 → 01 → 03 → 04 → 05 → 06 → 07 → 08 → 09 → 10 → 11 → 12 → 13. Conventions (01) are finalised after the spike (02): energy variable and density units may change. Tasks 03–07 are independent of each other once 01 is fixed.

## 0.2 — real sites

Coordinate systems with pyproj; true north, meridian convergence and magnetic declination; closed-mesh bodies with a ray-casting backend; survey import (LAS/LAZ, CSV points, LandXML, DXF, GeoTIFF) and point cloud → closed mesh; height-field horizons as layer boundaries; structures from meshes; voxel overlay and ray-by-voxel path matrix; template of a private site workspace outside the repository.

## 0.3 — physics and planning

Range tables by default; altitude correction; seasonal correction with optional MCEq; geomagnetic cut-off as a systematic of the open-sky reference; multiple-scattering kernel; angular efficiency tables; latent-image fading from temperature and humidity logs; track-density budget including transport and storage; systematics framework and profile-likelihood significance; sensitivity volume per voxel; optimisation of niches and tilts under film and scanning budgets; required accuracy of orientation and wall survey; report on how the improved physics changes the preliminary numbers.

## 0.4 — data

Synthetic tracks in the scanning output format; readers for scanning text output and FEDRA ROOT; binning with the same grid and cuts; data/model ratio, pulls, χ², region fits with nuisance parameters; straw telescope (level-1 response) and generic two-plane telescope; first comparison with an archival emulsion exposure.

## 1.0 — release and paper

Validation against Mulder/PUMAS and MUTE; documentation (guide, theory, configuration reference generated from the schema, API); one script per paper figure; Zenodo DOI; PyPI; software paper.

## 2.0 — inversion

Joint regularised inversion on the voxel grid with borehole priors; comparison with real exposure data.
