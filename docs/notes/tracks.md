# Tracks

Roadmap task 0.1-15, 28 November 2026. Modules `muoscene.detectors.tracks` (layer 2), `muoscene.data.binning` (layer 4) and `muoscene.io.tracks` (layer 5); conventions, revision 6, sections "Detector pose" (slopes, scan frame, acceptance), "Direction grid" and "Tracks"; the file formats are in `docs/formats.md`. Every number below comes from the named test of `tests/unit/test_detectors_tracks.py`, `test_io_tracks.py` or `test_data_binning.py`, or from the command next to it. No file, line or number of a real exposure is in the repository.

## What is there

| Object | Module | What it is |
|---|---|---|
| `TrackSet` | `detectors.tracks` | the tracks of one detector: slopes and, where the scan gives them, positions, in the frame of the detector or of the scan, with the area scanned, the length of the exposure, the cuts applied and where they come from; immutable; `.npz` files with a format version |
| `read_fedra_text` | `io.tracks` | the text output of the FEDRA macro `evaluate_VT_nseg.C` into a `TrackSet` |
| `ScanFrame` | `io.tracks` | the scan-frame transform of an exposure: which detector axis each axis of the scan runs along, with its sign, and where the reference point of the detector is on the stage |
| `bin_tracks`, `BinnedTracks` | `data.binning` | counts of tracks per bin of a `DirectionGrid` with the cuts of the model, the bins the model accepts, and every track left out with its reason |

What a user does with them, from a file to a figure (`muoscene.viz.maps.array_figure`, roadmap 0.1-12a):

```python
from muoscene.data import bin_tracks
from muoscene.detectors import DirectionGrid
from muoscene.io.tracks import read_fedra_text
from muoscene.viz import maps, style

tracks = read_fedra_text("tracks.txt", detector="film A")
binned = bin_tracks(tracks, DirectionGrid.square(), tan_max=1.0)
figure = maps.array_figure(
    binned.counts, binned.grid, status=binned.status, tan_max=1.0, label="tracks per bin"
)
style.save(figure, "tracks.png")
```

Without the pose of the detector and its scan-frame transform the figure is in the axes of the scan, with no zenith and no compass; with them (`scan_frame=` to the reader, `pose=` to `bin_tracks` and to `array_figure`) it is on the slopes of the detector, as the maps of the model.

## Decisions of 28 November 2026

| Question | Decided | Alternatives and their price |
|---|---|---|
| What a `TrackSet` holds | slopes `tx`, `ty`; positions `x_m`, `y_m` or none; `frame`; `source` (data or synthetic); `detector`; `area_cm2`, `exposure_days` or none; `cuts` and `meta`, plain values that go to JSON. No number of segments, first film or quality: the text output of the first exposures has none; they come with the reader that brings them, and with them version 2 of the file format (decided by the author) | three fields that are always empty: a reader, a figure and a comparison would each have to know that they mean nothing, and the format would claim columns no file has |
| Frames | "detector" — the axes of the detector, positions from its reference point — or "scan" — the axes and origin of the stage, as in the file. A set does not change its frame: the transform is applied once, by the reader, and recorded in `meta`; only a set in the frame of the detector is compared with a model, and `bin_tracks` refuses a pose for one in the frame of the scan (decided by the author) | a method that turns a set: two places where a transform can be applied, and a set whose history has to be read from its `meta` to know whether it already was |
| The scan-frame transform | `ScanFrame(origin_um=, x=, y=, z=)`: the detector axis each axis of the scan runs along, with its sign, z along or against the normal, and the stage coordinates of the reference point; all 16 combinations, the 8 left-handed ones included, because the transform is what the exposure log says (decided by the author) | a turn by 0, 90, 180 or 270° and a flip: the same 16 matrices, but it hides that a film scanned from its other side turns positions and slopes differently — x → x, y → −y for a point, tx → −tx, ty → ty for a slope (`test_scan_frames_by_hand`) |
| What Theta is | the tangent of the angle from the normal, by default: in FEDRA `EdbSegP::Theta()` is sqrt(TX² + TY²), and the group wrote so; `theta="radians"` reads it as an angle and refuses π/2 or more. Phi is atan2(TY, TX), as `EdbSegP::Phi()` | a guess from the values: they cannot tell a tangent below 1.57 from an angle below π/2 |
| Reading without a loop over lines | the labels of a track line are replaced by numbers, `numpy.loadtxt` reads eight columns and four of them must be the labels in their order; a loop over the lines runs only when that fails, to name the line. 10⁶ tracks, 65 MB, in 1.2 s on the Linux VM (the command below) | a regular expression or a split per line in Python: the same result, and a loop over every line of every file |
| A line that cannot be read | refused with the name of the file and the number of the line, never skipped: a header that promises a number of tracks is checked against it, and after the tracks only `(int) <n>` may follow | skipping bad lines with a warning: a count that silently loses tracks is a wrong count |
| A track on the edge of a bin | the bin above the edge; the last edge belongs to the last bin, as in `numpy.histogram2d`; edges compared as the grid holds them, without a tolerance. So 0.05 opens the bin [0.05, 0.10) (`test_edges_belong_to_the_bin_above`) | the nearest bin within a tolerance: one track could go to two bins on two systems |
| A track on the acceptance circle | inside, as `Detector.accepts`: sqrt(tx² + ty²) ≤ `tan_max` (`test_the_circle_is_inside_and_what_is_left_out_is_counted`) | — |
| The side of the plane | with a pose, a track counts when its line looks at the sky in front of the detector and above the horizon: `Pose.direction(tx, ty)` has u_z > 0; otherwise it is counted in `behind`. A film records such a muon, but the acceptance of the conventions leaves it out, and so does the model. Without a pose only the circle applies: a set in the frame of a scan has no side | counting a track when its sub-sample is one the model uses: data and model would cover the same region of slopes, but tracks behind the plane would be counted, and the counts would depend on `supersample` |
| Which bins are accepted | those where `muoscene.detectors.lines_of_sight`, the function the forward model calls, gives a line with the same `tan_max` and `supersample`; without a pose, those with a part in the circle. So `accepted` is false exactly where the status of the map is 1, on both films of `examples/open_sky` (`test_the_bins_left_out_are_those_the_model_does_not_accept`) | a second rule written for the tracks: it would agree with the model only until one of them changes |
| Nothing lost | every track is in exactly one of the counts, `outside_circle`, `outside_grid`, `behind` and `not_accepted`; `BinnedTracks` refuses numbers that do not add up, and counts are zero in the bins the model does not accept | — |

The model cuts the horizon and the back of the plane by whole sub-samples (`lines_of_sight`), the tracks are cut exactly. In the bins those lines cross the two differ by up to a row of sub-samples, a quarter of a bin at the default `supersample` of 4. On the film of `examples/open_sky` tilted by 60° to the south the horizon, ty = cot 60° = 0.5774, crosses 34 accepted bins: the model counts their sub-samples up to ty = 0.575, the tracks are counted up to the horizon. Whether such bins enter a comparison is decided by roadmap 0.1-16.

## Tests

| What | Number | Test |
|---|---|---|
| tracks drawn bin by bin from a map made by hand, one Poisson draw per bin, uniform in the part of each bin inside the circle, then shuffled | the counts come back exactly in all 1600 bins | `test_tracks_drawn_in_their_bins_come_back_to_them` |
| tracks drawn over the cone of 45° independently of the grid — an isotropic sky, halved in a block of 8 × 8 bins so that a swap of tx and ty or of a sign shows — against the exact accepted solid angle of each bin; Pearson's χ² over the bins expecting at least 10 tracks, seed 20261128 | 194 697 tracks; χ² = 1315.5 for K = 1296 bins, 0.38 standard deviations above K, where the test fails beyond 4; tx and ty swapped: 192, tx mirrored: 183 standard deviations | `test_tracks_drawn_over_the_cone_agree_with_the_map_within_poisson_errors` |
| how often that test would fail by chance for another seed | 6·10⁻⁵ with a normal tail, 8·10⁻⁵ with the skew of a χ² of 1300 degrees of freedom (Wilson–Hilferty) | — |
| the bins `bin_tracks` does not accept and the bins of status 1 of the map, for both films of `examples/open_sky` above 1 GeV | the same; 276 and 476 bins; no track above the horizon falls in a bin the model leaves out | `test_the_bins_left_out_are_those_the_model_does_not_accept` |
| every one of the 16 scan frames against the product of its matrix with a point (dx, dy, 0) and a direction (tx, ty, 1) | within 1e-15; 8 are right-handed | `test_sixteen_scan_frames` |
| a file with a decimal comma, with CR LF, with tabs and double spaces | the same numbers, bit for bit, as with a point, LF and single spaces | `test_comma_line_ends_and_blanks_do_not_change_the_numbers` |

The time to read a file, on the Linux VM, from a folder outside the repository:

```bash
python - <<'EOF'
import time

import numpy as np

from muoscene.io.tracks import read_fedra_text

rng = np.random.default_rng(1)
n = 10**6
x, y = rng.uniform(1e4, 1e5, (2, n))
theta, phi = rng.uniform(0, 1.5, n), rng.uniform(-3.14159, 3.14159, n)
lines = [f"VT param: X = {a:g} Y = {b:g} Theta = {c:g} Phi = {d:g}" for a, b, c, d in zip(x, y, theta, phi)]
with open("invented.txt", "w") as file:
    file.write("\n".join([f"N VTs = {n}", *lines, "(int) 1"]) + "\n")
start = time.perf_counter()
tracks = read_fedra_text("invented.txt", detector="d")
print(len(tracks), round(time.perf_counter() - start, 2), "s")  # 1000000 1.2 s
EOF
```

## Not here

Synthetic tracks drawn from an expected map and their files in the format of a scan (`muoscene.synth`, roadmap 0.4; task EZ-05); the reader of FEDRA ROOT and of the output of the processing program of the group, if the group asks for one (0.4); the ratio of data to model, its errors and its normalisation (0.1-16, `muoscene.data.compare`); the efficiency against the angle and the cuts of the scan, which the group has not given; the schema of a `ScanFrame` in a configuration (0.1-11).
