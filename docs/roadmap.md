# Roadmap

| Version | Target | Content |
|---|---|---|
| 0.1 | Dec 2026 | scenes of primitives through Calzone and Mulder, flux tables at the site altitude, parity with the preliminary calculation on a synthetic canyon, first 3D scene |
| 0.2 | May 2027 | real sites: coordinate systems, meshes, survey import, layer horizons, voxels |
| 0.3 | Aug 2027 | physics and planning: reference at the site altitude, scattering, emulsion fading and track-density budget, systematics, profile likelihood, sensitivity volume, scheme optimisation |
| 0.4 | Nov 2027 | data: synthetic tracks, track readers, comparison and fits; straw telescope |
| 1.0 | Dec 2027 – Mar 2028 | validation suite, documentation, Zenodo DOI, PyPI, software paper |
| 2.0 | 2029 | inversion on the voxel grid; real exposure data |

One task is one working session and one or a few commits. A task is done when it has code, tests, docstrings, a `CHANGELOG.md` line and its box ticked here.

Architecture (decided 13 November 2026, `docs/notes/mulder.md`): geometry navigation and muon physics come from Calzone and Mulder; muoscene adds the scene description, a fast flux path through tables built by Mulder, detector responses, planning, and comparison with data.

## 0.1 — scenes, flux tables, parity

Goal: `muoscene plan examples/canyon_demo/run.yaml` computes expected maps and significance for a synthetic canyon with a pier at the site altitude, reproduces the preliminary sensitivity calculation in legacy mode, and renders the scene. Conventions — `docs/conventions.md`.

- [x] **0.1-02 Prior-art spike: Mulder and Calzone.** Thirteen experiments and the decision: hybrid on top of Mulder and Calzone. `docs/notes/mulder.md`, `docs/notes/spike/mulder_spike.py`.
- [ ] **0.1-01 Conventions.** Revised after the spike: backends and units, kinetic energy, detector energy response, altitude cut-off, flux tables, equivalent opacity, bodies and anomalies. *Done when:* `docs/conventions.md` is approved.
- [ ] **0.1-03 Reference numbers.** Port the physics of the preliminary scripts (spectrum, threshold, path in the wall with a niche, sphere and box chords, cone counts and significance, tilted films, track density) to `tests/regression/legacy.py` as pure functions with explicit parameters; freeze `tests/regression/data/legacy_reference.json` for a synthetic parameter set (two banks with flat and rising brows, niches at two depths, tilts 0° and 20°, sphere and box anomalies) and `mulder_reference.json` with the numbers of the spike (reference flux, altitude, thresholds, slab, cut-off). *Done when:* regeneration is bit-identical with fixed seeds; no parameters of a real site; a new Mulder version that changes the physics fails the test.
- [ ] **0.1-04 Materials.** Material descriptions (formula or mass fractions, mixtures, moisture to composition and bulk density) registered with Mulder; library with sources: loam, loess, sand, concrete, standard rock, water, air; equivalent-opacity weights from Mulder ranges. *Done when:* the threshold table of the spike is reproduced; concrete weighs 0.995 relative to loam at 7 GeV; moisture round trip.
- [ ] **0.1-05 Flux tables.** `FluxTable` (differential flux on kinetic energy, sine of elevation and opacity; interpolation; integral over a detector energy response); `MulderTable` (detector altitude, atmosphere model, optional date for the geomagnetic field, dense slab of the reference material, disk cache by input hash; altitude cut-off recorded and enforced); `LegacyTable` (Guan at sea level with a + bE). *Done when:* 60.50 m⁻² s⁻¹ sr⁻¹ above 1 GeV at sea level; slab and cut-off numbers of the spike; legacy parity to 1e-6; interpolation within 0.2 % of direct Mulder at random points.
- [ ] **0.1-06 Geometry.** Scene description: ground as an extruded cross-section, horizontal layers, box, cylinder, pile grid, niche, spherical anomaly; builder to a Calzone tree (centimetres, positions relative to mothers, subtraction of cuts, splitting of structures at layer boundaries, air at the site density); `Scene.scan` → `PathTable`; density hypotheses on the `PathTable`; checks (anomaly within one host, `Geometry.check` for overlaps). *Done when:* wall-with-niche path parity with the reference to 1e-9 m; analytic chords of sphere, cylinder and box; layers add up to the ground; at least 5e4 directions/s for the demo pier.
- [ ] **0.1-07 Site and detectors I.** `Site` (origin by latitude, longitude, altitude; ENU; CRS conversion in 0.2); `Pose` per conventions; `DirectionGrid` with exact bin solid angles and sub-sampling; `EmulsionStack` (area, films, cone, energy threshold, constant efficiency, projected area, track density above the threshold). *Done when:* grid solid angles sum to the analytic cone; pose round trips; legacy track-density parity within Monte Carlo tolerance.
- [ ] **0.1-08 Forward model.** Grid → scan → `PathTable` → equivalent opacity per hypothesis → table lookup → counts; `ExpectedMap`; the fast path is compared with direct Mulder flux on a random sample of directions in every run and the deviation is stored. *Done when:* deviation within 0.5 % on the demo; open sky and uniform slab match Mulder; a hypothesis without a new scan equals a rebuilt scene.
- [ ] **0.1-09 Significance and exposure.** Legacy cone estimator; binned χ² and Asimov likelihood-ratio significance; days to kσ; track density versus days against the scanning limit; Markdown table in the layout of the preliminary calculation. *Done when:* 10.3σ of the reference void reproduced in legacy mode within 2 %; unit tests of the formulas.
- [ ] **0.1-10 Examples.** `examples/canyon_demo` (generator: canyon depth and widths, two banks with different brows, two niche levels, a pier on each bank with a pile cap and 3 × 3 piles, two horizontal layers) with `site.yaml`, `detectors.yaml`, `run.yaml`; `examples/open_sky`, `examples/slab`. *Done when:* the parameters are documented and differ from any real site.
- [ ] **0.1-11 Configuration, pipeline, CLI, provenance.** pydantic schema for every 0.1 object (`extra="forbid"`, units in field names); loader with `include` and relative paths; `from_config` builders; orchestration in `muoscene.pipeline`; `muoscene check` and `muoscene plan`; run folder with resolved configuration, input hashes, package and backend versions, git commit, NetCDF and Markdown report. *Done when:* examples run from the command line; broken configurations give readable errors.
- [ ] **0.1-12 Visualisation I.** Angular maps (counts, opacity, transmission, significance per hypothesis) in one style; PyVista scene (canyon, niches, pile caps, piles, anomalies, detectors with view cones coloured by significance) exported to HTML. *Done when:* the HTML opens in a browser and the scene matches the configuration.
- [ ] **0.1-13 Release 0.1.0.** README quickstart, `docs/legacy.md`, a validation script for the scattering factor of the demo void in discrete mode, CHANGELOG, roadmap; CI green; tag `v0.1.0` (no Zenodo, no PyPI yet).

Order: 01 → 03 → 04 → 05 → 06 → 07 → 08 → 09 → 10 → 11 → 12 → 13. Tasks 03–07 are independent of each other once 01 is fixed. Version 0.1 uses continuous transport only; scattering is quantified but not applied to maps.

## 0.2 — real sites

Coordinate systems with pyproj; true north, meridian convergence, magnetic declination from IGRF; mesh bodies (OBJ and STL with face normals) and DEM to a closed mesh through Calzone; survey import (LAS/LAZ, CSV points, LandXML, DXF, GeoTIFF) and point cloud to a closed mesh; layer horizons as nested height-field solids; structures crossing non-planar horizons; voxel overlay and ray-by-voxel path matrix; template of a private site workspace outside the repository.

## 0.3 — physics and planning

Tabulated reference flux at the site altitude (MCEq or measured) for muons below the altitude cut-off, needed for the open-sky reference and the track-density budget; seasonal atmospheres; geomagnetic field as a systematic of the open-sky reference; scattering — discrete-mode maps for final numbers and a kernel calibrated on them for planning loops; angular efficiency tables; latent-image fading from temperature and humidity logs; track-density budget including transport and storage; systematics framework and profile-likelihood significance; sensitivity volume per voxel; optimisation of niches and tilts under film and scanning budgets; required accuracy of orientation and wall survey; report on how the full physics changes the preliminary numbers.

## 0.4 — data

`TrackSet` in `detectors.tracks`; synthetic tracks in the scanning output format; readers for scanning text output and FEDRA ROOT; binning with the same grid and cuts; data/model ratio, pulls, χ², region fits with nuisance parameters; straw telescope (level-1 response) and generic two-plane telescope; first comparison with an archival emulsion exposure.

## 1.0 — release and paper

Validation: analytic cases, the preliminary calculation, continuous against discrete transport, MUTE for intensities at depth, an archival exposure; documentation (guide, theory, configuration reference generated from the schema, API); one script per paper figure; Zenodo DOI; PyPI; software paper.

## 2.0 — inversion

Joint regularised inversion on the voxel grid with borehole priors; comparison with real exposure data.
