# Roadmap

| Version | Target | Content |
|---|---|---|
| 0.1 | early Jan 2027 | scenes of primitives through Calzone and Mulder, open-sky reference at the site altitude, flux tables, parity with the preliminary calculation on a synthetic canyon, first 3D scene |
| 0.2 | May 2027 | real sites: coordinate systems, meshes, survey import, layer horizons, voxels |
| 0.3 | Aug 2027 | physics and planning: scattering, mixed-mode tables for thick overburden, emulsion fading and track-density budget, systematics, profile likelihood, sensitivity volume, scheme optimisation |
| 0.4 | Nov 2027 | data: synthetic tracks, track readers, comparison and fits; straw telescope |
| 1.0 | Dec 2027 – Mar 2028 | validation suite, documentation, Zenodo DOI, PyPI, software paper |
| 2.0 | 2029 | inversion on the voxel grid; real exposure data |

One task is one working session and one or a few commits. A task is done when it has code, tests, docstrings, a `CHANGELOG.md` line and its box ticked here.

Architecture (decided 13 November 2026, `docs/notes/mulder.md`): geometry navigation and muon physics come from Calzone and Mulder; muoscene adds the scene description, a fast flux path through tables built by Mulder, detector responses, planning, and comparison with data.

## 0.1 — scenes, flux tables, parity

Goal: `muoscene plan examples/canyon_demo/run.yaml` computes expected maps and significance for a synthetic canyon with a pier at the site altitude, reproduces the preliminary sensitivity calculation in legacy mode, and renders the scene. Conventions — `docs/conventions.md`.

- [x] **0.1-02 Prior-art spike: Mulder and Calzone.** Thirteen experiments and the decision: hybrid on top of Mulder and Calzone. `docs/notes/mulder.md`, `docs/notes/spike/mulder_spike.py`.
- [x] **0.1-01 Conventions.** Approved on 15 November 2026: kinetic energy at the detector, default threshold 1 GeV, reference at the site altitude in 0.1, runs outside the domain of validity refused.
- [x] **0.1-03 Reference numbers.** The preliminary scripts ported to `tests/regression/legacy.py` (bit-identical to the originals) and frozen on a synthetic canyon in `legacy_reference.json` (55 cases, 10.31σ for the reference void); the spike frozen in `mulder_reference.json` (57 cases at a relative tolerance of 10⁻⁵, above the spread between Linux and macOS, including the tables behind discrete transport and the wall paths of the synthetic canyon); `docs/legacy.md` with the sizes of the simplifications. Found on the way: integrals over all energies across the altitude cut-off were wrong by up to 1.3 % in the spike (season and geomagnetic numbers corrected in `docs/notes/mulder.md`); Mulder gives the air of a local geometry sea-level density: 1000 m of it above the brow lowers the counts of the reference void at 1150 m by 1 % (0.1-06: air at the site density, or a thin air layer).
- [x] **0.1-04 Materials.** `muoscene.materials`: materials from formulas, element fractions or oxide analyses with an explicit mean excitation energy (tabulated, or the Berger–Seltzer rule); mixtures; soils from dry density and gravimetric water content as a mixture of skeleton and water, so that Mulder takes the density effect at the particle density; names in Mulder derived from the definition; CSDA ranges, energy to cross, weights w(K) and equivalent opacity. Library with sources: water, standard rock, concrete, air, quartz, calcite, silicate and loess skeletons. Mulder ranges agree with the PDG tables within 0.52 %; the spike thresholds are reproduced exactly; concrete weighs 0.9952 against the spike loam at 7 GeV; moisture round trip to 1e-13. Sizes of the choices: `docs/notes/materials.md`.
- [ ] **0.1-05 Flux tables.** `FluxTable` (differential flux on kinetic energy, sine of elevation and opacity; interpolation; integral over a detector energy response); `MulderTable` (detector altitude, atmosphere model, optional date for the geomagnetic field, dense slab of the reference material, disk cache by input hash, reference chosen per conventions; altitude cut-off recorded and enforced for a sea-level reference); `LegacyTable` (Guan at sea level with a + bE). *Done when:* 60.50 m⁻² s⁻¹ sr⁻¹ above 1 GeV at sea level; slab and cut-off numbers of the spike; legacy parity to 1e-6; interpolation within 0.2 % of direct Mulder at random points; a site below sea level (−430 m) gives a flux below the sea-level one and no cut-off; a direction with more than 3e4 g/cm² is refused; whether Mulder applies muon decay along the air path is established and recorded.
- [ ] **0.1-05b Reference at the site altitude.** Open-sky differential flux in K and elevation for a given altitude, latitude and date, computed with MCEq (optional `mceq` extra; tables cached and shipped for a grid of altitudes so that the core runs without MCEq), passed to Mulder as a tabulated reference; normalisation against GCCLY15 at sea level; altitude dependence checked against published muon measurements at mountain altitudes (sources chosen in the task and listed in `docs/notes/reference.md`); a measured open-sky spectrum can replace the table. *Done when:* at sea level the table agrees with GCCLY15 within its stated uncertainty above 1 GeV; at 1150 m and 3340 m there is no cut-off and the flux above 1 GeV agrees with the sea-level reference carried up within the size of the decay and production effects, which are quantified; the choice of reference follows the conventions automatically.
- [ ] **0.1-06 Geometry.** Scene description: ground as an extruded cross-section, horizontal layers, box, cylinder, pile grid, niche, spherical anomaly; builder to a Calzone tree (centimetres, positions relative to mothers, subtraction of cuts, splitting of structures at layer boundaries, air at the site density); `Scene.scan` → `PathTable`; density hypotheses on the `PathTable`; checks (anomaly within one host, `Geometry.check` for overlaps). *Done when:* wall-with-niche path parity with the reference to 1e-9 m; analytic chords of sphere, cylinder and box; layers add up to the ground; at least 5e4 directions/s for the demo pier.
- [ ] **0.1-07 Site and detectors I.** `Site` (origin by latitude, longitude, altitude; ENU; CRS conversion in 0.2); `Pose` per conventions; `DirectionGrid` with exact bin solid angles and sub-sampling; `EmulsionStack` (area, films, cone, energy threshold, constant efficiency, projected area, track density above the threshold). *Done when:* grid solid angles sum to the analytic cone; pose round trips; legacy track-density parity within Monte Carlo tolerance.
- [ ] **0.1-08 Forward model.** Grid → scan → `PathTable` → equivalent opacity per hypothesis → table lookup → counts; `ExpectedMap`; the fast path is compared with direct Mulder flux on a random sample of directions in every run and the deviation is stored. *Done when:* deviation within 0.5 % on the demo; open sky and uniform slab match Mulder; a hypothesis without a new scan equals a rebuilt scene.
- [ ] **0.1-09 Significance and exposure.** Legacy cone estimator; binned χ² and Asimov likelihood-ratio significance; days to kσ; track density versus days against the scanning limit; Markdown table in the layout of the preliminary calculation. *Done when:* 10.3σ of the reference void reproduced in legacy mode within 2 %; unit tests of the formulas.
- [ ] **0.1-10 Examples.** `examples/canyon_demo` (generator: canyon depth and widths, two banks with different brows, two niche levels, a pier on each bank with a pile cap and 3 × 3 piles, two horizontal layers) with `site.yaml`, `detectors.yaml`, `run.yaml`; `examples/open_sky`, `examples/slab`. *Done when:* the parameters are documented and differ from any real site.
- [ ] **0.1-11 Configuration, pipeline, CLI, provenance.** pydantic schema for every 0.1 object (`extra="forbid"`, units in field names); loader with `include` and relative paths; `from_config` builders; orchestration in `muoscene.pipeline`; `muoscene check` and `muoscene plan`; run folder with resolved configuration, input hashes, package and backend versions, git commit, NetCDF and Markdown report. *Done when:* examples run from the command line; broken configurations give readable errors.
- [ ] **0.1-12 Visualisation I.** Angular maps (counts, opacity, transmission, significance per hypothesis) in one style; PyVista scene (canyon, niches, pile caps, piles, anomalies, detectors with view cones coloured by significance) exported to HTML. *Done when:* the HTML opens in a browser and the scene matches the configuration.
- [ ] **0.1-13 Release 0.1.0.** README quickstart, `docs/legacy.md`, a validation script for the scattering factor of the demo void in discrete mode, CHANGELOG, roadmap; CI green; tag `v0.1.0` (no Zenodo, no PyPI yet).

Order: 01 → 03 → 04 → 05 → 05b → 06 → 07 → 08 → 09 → 10 → 11 → 12 → 13. Tasks 03–07 are independent of each other; 05b needs 05. Version 0.1 uses continuous transport only; scattering is quantified but not applied to maps.

## 0.2 — real sites

Coordinate systems with pyproj; true north, meridian convergence, magnetic declination from IGRF; mesh bodies (OBJ and STL with face normals) and DEM to a closed mesh through Calzone; survey import (LAS/LAZ, CSV points, LandXML, DXF, GeoTIFF) and point cloud to a closed mesh; layer horizons as nested height-field solids; structures crossing non-planar horizons; voxel overlay and ray-by-voxel path matrix; template of a private site workspace outside the repository.

## 0.3 — physics and planning

Seasonal atmospheres and seasonal site-altitude references; geomagnetic field as a systematic of the open-sky reference; mixed-mode tables for opacity beyond 3e4 g/cm²; scattering — discrete-mode maps for final numbers and a kernel calibrated on them for planning loops; angular efficiency tables; latent-image fading from temperature and humidity logs; track-density budget including transport and storage; systematics framework and profile-likelihood significance; sensitivity volume per voxel; optimisation of niches and tilts under film and scanning budgets; required accuracy of orientation and wall survey; report on how the full physics changes the preliminary numbers.

## 0.4 — data

`TrackSet` in `detectors.tracks`; synthetic tracks in the scanning output format; readers for scanning text output and FEDRA ROOT; binning with the same grid and cuts; data/model ratio, pulls, χ², region fits with nuisance parameters; straw telescope (level-1 response) and generic two-plane telescope; first comparison with an archival emulsion exposure.

## 1.0 — release and paper

Validation: analytic cases, the preliminary calculation, continuous against discrete transport, MUTE for intensities at depth, an archival exposure; documentation (guide, theory, configuration reference generated from the schema, API); one script per paper figure; Zenodo DOI; PyPI; software paper.

## 2.0 — inversion

Joint regularised inversion on the voxel grid with borehole priors; comparison with real exposure data.
