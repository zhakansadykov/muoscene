# Conventions

Status: approved on 15 November 2026, after the Mulder and Calzone spike (`docs/notes/mulder.md`); revision 5 of 27 November 2026. Revision 2 (16 November): equivalent-opacity weights depend on energy, a soil is a mixture of skeleton and water. Revision 3 (16 November, after the review): decay along the air path, numbers quoted for two detector thresholds, the site-altitude reference from 0.1.1. Revision 4 (16 November, roadmap 0.1-05): a flux table holds the open-sky flux on a grid of energy and elevation and takes the opacity through the range of its reference material. Revision 5 (27 November, after roadmap tasks 0.1-07 to 0.1-10b): the axes of a scene, the tilt of a detector, the two sides of a film in a track density, the solid angle of the part of a bin in the acceptance, lines of sight through an edge of a mesh, the weights of the equivalent opacity between the threshold and the energy above the path, the air the solid bodies displace, the comparison with Mulder in every run, bins without a flux marked on request, the variables of the expected map. Every module follows these rules; a change here is a breaking change, is made on purpose in a commit of its own and is recorded in `CHANGELOG.md`.

## Backends and units

Scenes are built with Calzone and scanned by Mulder; reference flux, atmosphere and energy loss come from Mulder. Only `geometry`, `materials`, `flux`, `transport` and `site` talk to them (`calzone` only `geometry`), and units are converted at that boundary.

| Quantity | muoscene | Mulder | Calzone |
|---|---|---|---|
| length | m | m | cm |
| density | g/cm³ | kg/m³ | — |
| opacity X | g/cm² | kg/m² | — |
| energy | GeV, kinetic | GeV, kinetic | — |
| flux | m⁻² s⁻¹ sr⁻¹ (per GeV when differential) | the same | — |
| detector area | cm² | — | — |
| time | s; exposure by dates or in days | — | — |
| angles | rad internally, degrees in configuration | degrees | degrees |

X = ρ L × 100 for ρ in g/cm³ and L in m; 1 m w.e. = 100 g/cm². Configuration field names carry their unit: `length_m`, `area_cm2`, `tilt_deg`, `exposure_days`, `density_g_cm3`, `k_min_gev`.

## Frames

**World frame — local ENU.** x points east, y north, z up; metres. The origin is the site origin given by latitude, longitude and altitude; z = 0 is that altitude above sea level. For a scene that is not turned it is the Mulder local frame of the scene. Survey coordinate systems (UTM, a contractor's local grid, WGS84) are converted only in `muoscene.site`, at input and output.

**Axes of a scene.** A scene is described along its own axes, which share the origin and the vertical of the world frame and are turned about the vertical. The site holds the angle: `scene_azimuth_deg` ψ is the azimuth of the +y axis of the scene from true north, clockwise, so that in ENU x_s = (cos ψ, −sin ψ, 0) and y_s = (sin ψ, cos ψ, 0). The Mulder local frame of a scene is the frame of these axes, `LocalFrame(azimuth=ψ)`. Points and directions pass between ENU and the axes of a scene only in `muoscene.site` (`Site.scan`, `Site.differential_flux`); poses, grids, detectors and output are in ENU.

**Directions.** A direction `u` is a unit vector pointing from the detector towards the sky, i.e. opposite to the muon momentum (the observation direction of Mulder). Zenith angle θ = arccos(u_z); elevation is 90° − θ. Azimuth φ is measured from true north, clockwise (towards east):

    u = (sin θ sin φ, sin θ cos φ, cos θ)

Only directions with u_z > 0 carry muons; upward-going particles are ignored.

**North.** Azimuths in configuration and output are relative to true north. Grid north (UTM meridian convergence) and magnetic north (declination for the date, IGRF) are converted in `muoscene.site`; a value given relative to them must say so explicitly.

## Detector pose

A detector has a reference point (centre of the sensitive volume) in ENU and three angles:

- `azimuth_deg` α — azimuth of the horizontal projection of the detector normal;
- `tilt_deg` β — angle between the detector normal and the vertical; 0 is a horizontal film or plane;
- `roll_deg` γ — rotation about the normal.

The tilt lies between 0° and 90°. At a tilt of zero the azimuth and the roll turn about the same axis in opposite senses, and the pose is written with an azimuth of zero.

With γ = 0 the detector axes in ENU are

    z_d = (sin β sin α, sin β cos α, cos β)      normal
    x_d = (cos α, −sin α, 0)                      horizontal, pointing to azimuth α + 90°
    y_d = z_d × x_d

and roll rotates x_d and y_d about z_d by γ (counter-clockwise seen from the tip of z_d). For α = β = γ = 0 the detector frame coincides with ENU.

**Slopes.** The direction variables of every detector are the slopes in its own frame:

    tx = (u · x_d) / (u · z_d),   ty = (u · y_d) / (u · z_d)

They are the slopes an emulsion scan reports (dx/dz, dy/dz) and the ratios of hit offsets to the base of a two-plane telescope. A slope pair describes a line, so u and −u give the same (tx, ty); the sky side is fixed by u_z > 0.

**Scan frame.** Film axes used by the scanning laboratory are related to (x_d, y_d) by a fixed rotation or flip recorded in the exposure log; `muoscene.io.tracks` applies it when reading tracks.

**Acceptance.** A detector accepts directions with u · z_d > 0 and sqrt(tx² + ty²) ≤ `tan_max`. The projected sensitive area is A (u · z_d). A detector point is never placed exactly on a boundary between bodies. A line of sight that passes within 1 nm of an edge or a vertex of a mesh is traced 3 nm aside (9 nm, or in another direction, where that lands near a second edge); one that runs at less than 0.57° to a face at such an edge, or along the edge, is refused (`docs/notes/geometry.md`).

A track density is not a count in the acceptance: a film records a muon whichever side it enters. It is summed over the lines of sight on both sides of the film, with |u · z_d| as the projection, and states the cone and the lowest elevation it was summed over.

## Direction grid

Bins are rectangles in (tx, ty). The solid angle of the rectangle [0, a] × [0, b] in slope space is

    Ω = arctan(a b / sqrt(1 + a² + b²))

and a general bin follows from four corner terms. Expected counts in a bin are integrated over sub-samples (n × n directions per bin, `supersample` in the configuration). The part of a bin within the acceptance circle has an exact solid angle as well, the contour integral of (1 − cos θ) dφ along its boundary; a sub-sample carries the exact solid angle of its rectangle within the circle and a direction inside it.

## Bodies, densities and anomalies

A scene is a list of bodies; each has a shape, a material (composition) and a bulk density. A soil is given by its mineral skeleton, dry density and gravimetric water content; its material is a mixture of the skeleton at the particle density and water, and its bulk density ρ_d (1 + w) goes to the geometry. Every material carries an explicit mean excitation energy: tabulated, or by the Berger–Seltzer rule (`docs/notes/materials.md`). Where bodies overlap, precedence is: cuts (niches, excavations), then structures, then layers, then the ground; the builder realises it by subtracting cuts and by splitting structures at layer boundaries. The air of the scene has the density of the atmosphere at the site altitude.

An anomaly is a body that scales the density of its host by (1 + `drho_rel`): −1 is a void, −0.2 a 20 % decompaction. It lies inside exactly one host body; the builder refuses or splits an anomaly that crosses a boundary. Voxel overlays follow the same rule.

`PathTable` stores the path length in every body for every direction. Opacity under a hypothesis is the sum of length times density, so density hypotheses and density systematics never need a new scan.

## Energy and flux

**Energy** is the kinetic energy K of the muon at the detector, in GeV.

**Detector energy response.** Each detector has an efficiency ε(K); the default is a threshold `k_min_gev = 1.0`. The response of an emulsion stack depends on the track selection of the scanning laboratory and is not known yet; until it is, numbers are quoted for thresholds of 0.2 and 1 GeV (`docs/numbers.md`); with a sea-level reference the lower one holds only under soil, where no cut-off is left. Counts in a bin are the exposure times the efficiency times the projected area times the integral over the bin and over K of ε(K) φ(K, direction). Momentum, where needed (scattering), is p = sqrt(K² + 2 K m_μ).

**Reference flux and the altitude cut-off.** The reference flux is the open-sky flux at the site, from which Mulder transports muons to the detector through the scene. Two kinds are supported:

- *Sea-level reference* — GCCLY15 (Guan et al. 2015), the Mulder default. Mulder carries it through the air between sea level and the site. Above sea level this defines the flux only above K_cut(elevation), the energy lost in that air: 0.28, 0.38, 0.53 GeV at 1150 m and 0.69, 0.97, 1.40 GeV at 3340 m for elevations of 90°, 45°, 30°. Mulder applies muon decay along that air path: the flux at altitude is the sea-level flux at the shifted energy times the ratio of stopping powers times exp(Δτ/τ), 1.21 at 1 GeV and 1.04 at 5 GeV at 1150 m. Production of muons in that air layer is not modelled. Below sea level the reference is carried down through the extra air and has no cut-off.
- *Site-altitude reference* — a table of the open-sky differential flux in K and elevation computed for the altitude of the site (a model of the atmospheric flux such as PARMA, cross-checked with MCEq, or a measured open-sky spectrum), passed to Mulder as a tabulated reference. It has no cut-off and includes production and decay in the air. It is the default when the site is more than 100 m above sea level, from version 0.1.1 on; version 0.1.0 has the sea-level reference only.

Where the energy response reaches below the cut-off of a sea-level reference there is no number: the run is refused or, on request, the bins are marked as undefined with the reason; nothing is put in their place.

**Flux table.** For each detector a `FluxTable` gives the differential flux against K, sine of elevation and opacity of the reference material. It holds the open-sky flux at the detector altitude on a grid of K and sine of elevation, computed by Mulder with the configured atmosphere and, optionally, the geomagnetic field of a date (then for one azimuth), together with the CSDA range R and the stopping power S of the reference material. Under an opacity X the flux is that of the open sky at the energy K_top with R(K_top) = R(K) + X, times S(K_top) / S(K): this is what continuous transport does, and every new table is compared with Mulder under a dense slab of the reference material at the detector altitude and refused if they differ (`docs/notes/flux.md`). Tables are cached on disk, outside the repository, by a hash of their inputs. With a sea-level reference a table holds only energies from 10 MeV above the cut-off, also under an opacity that leaves some of it, and a request below is refused. An integral over the energy response is a trapezoid sum that starts at the threshold itself; the response is always given, a table has no default threshold.

**Equivalent opacity.** Path lengths in solid bodies are combined as X_eq = Σ w_m ρ L, with weights from the CSDA muon ranges R_ref in the reference material and R_m in material m; at one kinetic energy their ratio R_ref(K) / R_m(K) is 0.995 for concrete relative to the loam of the spike at 7 GeV. For a detector whose response starts at K_min the weight is taken between K_min and the energy K_top a muon needs above the path to arrive with K_min: w_m = (R_ref(K_top) − R_ref(K_min)) / (R_m(K_top) − R_m(K_min)), with R_ref(K_top) = R_ref(K_min) + X_eq, found by iteration. It is exact for a path through one material, needs no order of the materials along the path and does not depend on which material is the reference. The atmosphere is part of the table: a table holds the open sky at the altitude of the detector, so the air that the solid bodies of a path displace, ρ_air times the length of the path in them, weighted like a material, is taken out of the opacity the table is entered with.

**Fast and direct paths.** Expected maps use the fast path: scan, equivalent opacity, table lookup. Every run also computes the direct Mulder continuous flux, for the scene as built, on a random sample of 256 accepted lines of sight per detector and stores the deviation with the seed of the sample; a run is refused above 0.5 %; on the invented pier it is 0.014 % rms above 1 GeV. Mulder discrete flux, which includes scattering, is used for final numbers and to calibrate the scattering kernel.

## Domain of validity

Every run checks where it stands and writes it to the run report. Outside the domain it is refused, not approximated.

| Condition | Version 0.1 | Beyond |
|---|---|---|
| site below sea level | sea-level reference carried down | — |
| site above sea level | 0.1.0: sea-level reference, only above K_cut; from 0.1.1: site-altitude reference, no cut-off | — |
| sea-level reference above sea level | above K_cut only; a run says what to do with accepted directions below it: refuse, or compute the rest and mark the bins | — |
| opacity up to 3e4 g/cm² (about 100 m of rock) along accepted directions | continuous transport | — |
| larger opacity | refused | tables in mixed mode, version 0.3 |
| targets of a few metres and less | continuous transport, scattering not applied; the report gives the scattering factor measured for the configuration when available | discrete mode and a scattering kernel, version 0.3 |

## Expected map

`ExpectedMap` is an `xarray.Dataset` with dimensions `detector`, `hypothesis`, `ty`, `tx` (later `variant` for systematics):

| Variable | Dimensions | Meaning |
|---|---|---|
| `solid_angle_sr` | ty, tx | solid angle of the bin: the whole rectangle in the slopes |
| `accepted_solid_angle_sr` | detector, ty, tx | solid angle of the accepted lines of sight of the bin: within the circle, in front of the detector, above the horizon |
| `opacity_gcm2` | detector, hypothesis, ty, tx | mean equivalent opacity over sub-samples |
| `counts` | detector, hypothesis, ty, tx | expected number of tracks |
| `transmission` | detector, hypothesis, ty, tx | flux through matter over open-sky flux |
| `status` | detector, hypothesis, ty, tx | 0 defined; 1 not accepted (counts are zero); 2 below the cut-off of the reference flux; 3 elevation outside the flux table; 4 opacity beyond the flux table (counts are nan for 2 to 4) |
| `open_sky_status` | detector, ty, tx | the same for the counts without matter, the denominator of `transmission` |
| `reference` | detector, ty, tx | bin sees open sky: no accepted line of sight of the bin has an equivalent opacity above `reference_opacity_g_cm2`, zero by default, under any hypothesis |
| `zenith_deg`, `azimuth_deg` | detector, ty, tx | bin centre in ENU on the side of the detector normal; a zenith angle above 90° means that the centre looks below the horizon; the azimuth of a vertical centre is 0 |
| `exposure_days` | detector | length of the exposure |

Attributes carry the package and backend versions, the configuration hash, creation time, exposure, and the deviation of the fast path from the direct one; also the site, the stations with their energy responses, the hypotheses, and the seed and the size of the sample of the comparison.

## Tracks

`TrackSet` (`muoscene.detectors.tracks`) holds measured and synthetic tracks in one layout: slopes `tx`, `ty` and position in the detector frame, film index, quality, and the source (data or synthetic). It lives with the detectors, not with file readers, because synthetic data and comparison with data use it without touching files; readers in `muoscene.io.tracks` convert scanning output into it and apply the scan-frame transform. Model and data are binned in the same `DirectionGrid` with the same cuts.

## Randomness

Every stochastic function takes a `numpy.random.Generator`; Mulder discrete runs take an explicit seed. Seeds are part of the configuration and are written to the run folder.

## Legacy mode

The preliminary calculation used for the 2027 grant application is reproduced by `LegacyTable` (Guan et al. 2015 at sea level, energy loss a + bE with a 0.1 GeV offset) and by the cone significance estimator with the best of five cone sizes. They are used only by `tests/regression`. `docs/legacy.md` lists its simplifications and their size.
