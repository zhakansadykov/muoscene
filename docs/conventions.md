# Conventions

Status: approved on 15 November 2026, after the Mulder and Calzone spike (`docs/notes/mulder.md`); revision 4 of 16 November 2026. Revision 2 (16 November): equivalent-opacity weights depend on energy, a soil is a mixture of skeleton and water. Revision 3 (16 November, after the review): decay along the air path, numbers quoted for two detector thresholds, the site-altitude reference from 0.1.1. Revision 4 (16 November, roadmap 0.1-05): a flux table holds the open-sky flux on a grid of energy and elevation and takes the opacity through the range of its reference material. Every module follows these rules; a change here is a breaking change, is made on purpose in a commit of its own and is recorded in `CHANGELOG.md`.

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

**World frame — local ENU.** x points east, y north, z up; metres. The origin is the site origin given by latitude, longitude and altitude; z = 0 is that altitude above sea level. It is the Mulder local frame of the scene. Survey coordinate systems (UTM, a contractor's local grid, WGS84) are converted only in `muoscene.site`, at input and output.

**Directions.** A direction `u` is a unit vector pointing from the detector towards the sky, i.e. opposite to the muon momentum (the observation direction of Mulder). Zenith angle θ = arccos(u_z); elevation is 90° − θ. Azimuth φ is measured from true north, clockwise (towards east):

    u = (sin θ sin φ, sin θ cos φ, cos θ)

Only directions with u_z > 0 carry muons; upward-going particles are ignored.

**North.** Azimuths in configuration and output are relative to true north. Grid north (UTM meridian convergence) and magnetic north (declination for the date, IGRF) are converted in `muoscene.site`; a value given relative to them must say so explicitly.

## Detector pose

A detector has a reference point (centre of the sensitive volume) in ENU and three angles:

- `azimuth_deg` α — azimuth of the horizontal projection of the detector normal;
- `tilt_deg` β — angle between the detector normal and the vertical; 0 is a horizontal film or plane;
- `roll_deg` γ — rotation about the normal.

With γ = 0 the detector axes in ENU are

    z_d = (sin β sin α, sin β cos α, cos β)      normal
    x_d = (cos α, −sin α, 0)                      horizontal, pointing to azimuth α + 90°
    y_d = z_d × x_d

and roll rotates x_d and y_d about z_d by γ (counter-clockwise seen from the tip of z_d). For α = β = γ = 0 the detector frame coincides with ENU.

**Slopes.** The direction variables of every detector are the slopes in its own frame:

    tx = (u · x_d) / (u · z_d),   ty = (u · y_d) / (u · z_d)

They are the slopes an emulsion scan reports (dx/dz, dy/dz) and the ratios of hit offsets to the base of a two-plane telescope. A slope pair describes a line, so u and −u give the same (tx, ty); the sky side is fixed by u_z > 0.

**Scan frame.** Film axes used by the scanning laboratory are related to (x_d, y_d) by a fixed rotation or flip recorded in the exposure log; `muoscene.io.tracks` applies it when reading tracks.

**Acceptance.** A detector accepts directions with u · z_d > 0 and sqrt(tx² + ty²) ≤ `tan_max`. The projected sensitive area is A (u · z_d). A detector point is never placed exactly on a boundary between bodies.

## Direction grid

Bins are rectangles in (tx, ty). The solid angle of the rectangle [0, a] × [0, b] in slope space is

    Ω = arctan(a b / sqrt(1 + a² + b²))

and a general bin follows from four corner terms. Expected counts in a bin are integrated over sub-samples (n × n directions per bin, `supersample` in the configuration).

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

A run whose energy response reaches below the cut-off of a sea-level reference is refused.

**Flux table.** For each detector a `FluxTable` gives the differential flux against K, sine of elevation and opacity of the reference material. It holds the open-sky flux at the detector altitude on a grid of K and sine of elevation, computed by Mulder with the configured atmosphere and, optionally, the geomagnetic field of a date (then for one azimuth), together with the CSDA range R and the stopping power S of the reference material. Under an opacity X the flux is that of the open sky at the energy K_top with R(K_top) = R(K) + X, times S(K_top) / S(K): this is what continuous transport does, and every new table is compared with Mulder under a dense slab of the reference material at the detector altitude and refused if they differ (`docs/notes/flux.md`). Tables are cached on disk, outside the repository, by a hash of their inputs. With a sea-level reference a table holds only energies from 10 MeV above the cut-off, also under an opacity that leaves some of it, and a request below is refused. An integral over the energy response is a trapezoid sum that starts at the threshold itself; the response is always given, a table has no default threshold.

**Equivalent opacity.** Path lengths in solid bodies are combined as X_eq = Σ w_m ρ L, where w_m(K) = R_ref(K) / R_m(K) is the ratio of CSDA muon ranges in the reference material and in material m at kinetic energy K (0.995 for concrete relative to the loam of the spike at 7 GeV). For a path through one material X_eq = R_ref(R_m⁻¹(X)) exactly; the energy at which the weights of a path through several materials are taken is set by the forward model (roadmap 0.1-08). The atmosphere is part of the table and is not counted.

**Fast and direct paths.** Expected maps use the fast path: scan, equivalent opacity, table lookup. Every run also computes the direct Mulder continuous flux for a random sample of directions and stores the deviation; in the spike it is 0.2 % rms. Mulder discrete flux, which includes scattering, is used for final numbers and to calibrate the scattering kernel.

## Domain of validity

Every run checks where it stands and writes it to the run report. Outside the domain it is refused, not approximated.

| Condition | Version 0.1 | Beyond |
|---|---|---|
| site below sea level | sea-level reference carried down | — |
| site above sea level | 0.1.0: sea-level reference, only above K_cut; from 0.1.1: site-altitude reference, no cut-off | — |
| sea-level reference above sea level | allowed only above K_cut for all accepted directions | — |
| opacity up to 3e4 g/cm² (about 100 m of rock) along accepted directions | continuous transport | — |
| larger opacity | refused | tables in mixed mode, version 0.3 |
| targets of a few metres and less | continuous transport, scattering not applied; the report gives the scattering factor measured for the configuration when available | discrete mode and a scattering kernel, version 0.3 |

## Expected map

`ExpectedMap` is an `xarray.Dataset` with dimensions `detector`, `hypothesis`, `ty`, `tx` (later `variant` for systematics):

| Variable | Dimensions | Meaning |
|---|---|---|
| `solid_angle_sr` | ty, tx | solid angle of the bin |
| `opacity_gcm2` | detector, hypothesis, ty, tx | mean equivalent opacity over sub-samples |
| `counts` | detector, hypothesis, ty, tx | expected number of tracks |
| `transmission` | detector, hypothesis, ty, tx | flux through matter over open-sky flux |
| `reference` | detector, ty, tx | bin sees open sky (opacity below the reference threshold) |
| `zenith_deg`, `azimuth_deg` | detector, ty, tx | bin centre in ENU |

Attributes carry the package and backend versions, the configuration hash, creation time, exposure, and the deviation of the fast path from the direct one.

## Tracks

`TrackSet` (`muoscene.detectors.tracks`) holds measured and synthetic tracks in one layout: slopes `tx`, `ty` and position in the detector frame, film index, quality, and the source (data or synthetic). It lives with the detectors, not with file readers, because synthetic data and comparison with data use it without touching files; readers in `muoscene.io.tracks` convert scanning output into it and apply the scan-frame transform. Model and data are binned in the same `DirectionGrid` with the same cuts.

## Randomness

Every stochastic function takes a `numpy.random.Generator`; Mulder discrete runs take an explicit seed. Seeds are part of the configuration and are written to the run folder.

## Legacy mode

The preliminary calculation used for the 2027 grant application is reproduced by `LegacyTable` (Guan et al. 2015 at sea level, energy loss a + bE with a 0.1 GeV offset) and by the cone significance estimator with the best of five cone sizes. They are used only by `tests/regression`. `docs/legacy.md` lists its simplifications and their size.
