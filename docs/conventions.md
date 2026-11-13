# Conventions

Status: revised after the Mulder and Calzone spike (`docs/notes/mulder.md`, 13 November 2026). Every module follows these rules; a change here is a breaking change and is recorded in `CHANGELOG.md`.

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

A scene is a list of bodies; each has a shape, a material (composition) and a bulk density. Moisture changes bulk density and composition in `muoscene.materials`. Where bodies overlap, precedence is: cuts (niches, excavations), then structures, then layers, then the ground; the builder realises it by subtracting cuts and by splitting structures at layer boundaries. The air of the scene has the density of the atmosphere at the site altitude.

An anomaly is a body that scales the density of its host by (1 + `drho_rel`): −1 is a void, −0.2 a 20 % decompaction. It lies inside exactly one host body; the builder refuses or splits an anomaly that crosses a boundary. Voxel overlays follow the same rule.

`PathTable` stores the path length in every body for every direction. Opacity under a hypothesis is the sum of length times density, so density hypotheses and density systematics never need a new scan.

## Energy and flux

**Energy** is the kinetic energy K of the muon at the detector, in GeV.

**Detector energy response.** Each detector has an efficiency ε(K); the default is a threshold `k_min_gev`. Counts in a bin are the exposure times the efficiency times the projected area times the integral over the bin and over K of ε(K) φ(K, direction).

**Reference flux and the altitude cut-off.** The default reference is GCCLY15 (Guan et al. 2015) at sea level; Mulder carries it through the atmosphere to the detector. For a detector above sea level this defines the flux only above K_cut(elevation), the energy lost in the air down to sea level: 0.28, 0.39 and 0.53 GeV at 1150 m for elevations of 90°, 45° and 30°. A run whose energy response reaches below the cut-off needs a tabulated reference at the site altitude; without one it is refused.

**Flux table.** For each detector a `FluxTable` holds the differential flux on a grid of K, sine of elevation and opacity, computed by Mulder under a dense slab of the reference material at the detector altitude, with the configured atmosphere and, optionally, the geomagnetic field of a date. Tables are cached on disk by a hash of their inputs.

**Equivalent opacity.** Path lengths in solid bodies are combined as X_eq = Σ w_m ρ L, where w_m is the ratio of the muon range in the reference material to that in material m at a reference energy (0.995 for concrete relative to loam). The atmosphere is part of the table and is not counted.

**Fast and direct paths.** Expected maps use the fast path: scan, equivalent opacity, table lookup. Every run also computes the direct Mulder continuous flux for a random sample of directions and stores the deviation; in the spike it is 0.2 % rms. Mulder discrete flux, which includes scattering, is used for final numbers and to calibrate the scattering kernel.

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
