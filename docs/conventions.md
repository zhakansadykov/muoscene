# Conventions

Status: draft for version 0.1, to be finalised after the Mulder and Calzone spike (roadmap task 0.1-02): the energy variable (total or kinetic) and the density unit may change. Every module follows these rules; a change here is a breaking change and is recorded in `CHANGELOG.md`.

## Frames

**World frame — local ENU.** x points east, y north, z up; metres. The origin is the site origin given by latitude, longitude and altitude; z = 0 is that altitude above sea level. All computations run in this frame. Survey coordinate systems (UTM, a contractor's local grid, WGS84) are converted only in `muoscene.site`, at input and output.

**Directions.** A direction `u` is a unit vector pointing from the detector towards the sky, i.e. opposite to the muon momentum. Zenith angle θ = arccos(u_z). Azimuth φ is measured from true north, clockwise (towards east):

    u = (sin θ sin φ, sin θ cos φ, cos θ)

Only directions with u_z > 0 carry muons; upward-going particles are ignored.

**North.** Azimuths in configuration and output are relative to true north. Grid north (UTM meridian convergence) and magnetic north (declination for the date) are converted in `muoscene.site`; a value given relative to them must say so explicitly.

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

**Acceptance.** A detector accepts directions with u · z_d > 0 and sqrt(tx² + ty²) ≤ `tan_max`. The projected sensitive area is A (u · z_d).

## Direction grid

Bins are rectangles in (tx, ty). The solid angle of the rectangle [0, a] × [0, b] in slope space is

    Ω = arctan(a b / sqrt(1 + a² + b²))

and a general bin follows from four corner terms. Expected counts in a bin are integrated over sub-samples (n × n directions per bin, `supersample` in the configuration).

## Units

| Quantity | Internal unit | Notes |
|---|---|---|
| length | m | |
| area of detectors | cm² | configuration fields end in `_cm2` |
| density | g/cm³ | |
| opacity X | g/cm² | X = ρ L × 100 for ρ in g/cm³ and L in m; 1 m w.e. = 100 g/cm² |
| energy | GeV | total energy E ≥ m_μ |
| flux | cm⁻² s⁻¹ sr⁻¹ | reports may show m⁻² s⁻¹ sr⁻¹ |
| time | s | exposure in configuration is given by dates or in days |
| angles | rad internally | configuration fields end in `_deg` |

Configuration field names carry their unit: `length_m`, `area_cm2`, `tilt_deg`, `exposure_days`.

## Energy threshold

E_min(X) is the total energy with which a muon traverses opacity X along the line of sight and still reaches the sensitive volume: the CSDA range of the material sequence equals X, plus the detector's own packaging (`e_floor_gev` or an equivalent opacity). The flux through matter is I(> E_min(X), cos θ).

## Materials, density and anomalies

A body has a material (composition, Z/A, mean excitation energy) and a bulk density. Moisture changes bulk density and composition in `muoscene.materials`.

Anomalies are bodies that scale the density of whatever lies inside them by (1 + `drho_rel`): −1 is a void, −0.2 a 20 % decompaction. Voxel overlays work the same way. Density changes are applied to stored path lengths without retracing rays.

Overlapping bodies are resolved by priority: cuts (niches, excavations) 100, structures 80, layers 40, ground 20; outside every body is air. Air opacity is ignored by default.

## Expected map

`ExpectedMap` is an `xarray.Dataset` with dimensions `detector`, `hypothesis`, `ty`, `tx` (later `variant` for systematics):

| Variable | Dimensions | Meaning |
|---|---|---|
| `solid_angle_sr` | ty, tx | solid angle of the bin |
| `opacity_gcm2` | detector, hypothesis, ty, tx | mean opacity over sub-samples |
| `counts` | detector, hypothesis, ty, tx | expected number of tracks |
| `transmission` | detector, hypothesis, ty, tx | flux through matter over open-sky flux |
| `reference` | detector, ty, tx | bin sees open sky (opacity below the reference threshold) |
| `zenith_deg`, `azimuth_deg` | detector, ty, tx | bin centre in ENU |

Attributes carry the package version, configuration hash, creation time and exposure.

## Tracks

`TrackSet` (`muoscene.detectors.tracks`) holds measured and synthetic tracks in one layout: slopes `tx`, `ty` and position in the detector frame, film index, quality, and the source (data or synthetic). It lives with the detectors, not with file readers, because synthetic data and comparison with data use it without touching files; readers in `muoscene.io.tracks` convert scanning output into it and apply the scan-frame transform. Model and data are binned in the same `DirectionGrid` with the same cuts.

## Randomness

Every stochastic function takes a `numpy.random.Generator`. Seeds are part of the configuration and are written to the run folder.

## Legacy mode

The preliminary calculation used for the 2027 grant application (spectrum clipped at cos θ = 0.25, energy loss a + bE with a 0.1 GeV offset, cone-based significance with the best of five cone sizes) is reproduced by explicit `legacy` options, used only by `tests/regression`. Defaults do not inherit its simplifications; `docs/legacy.md` lists them.
