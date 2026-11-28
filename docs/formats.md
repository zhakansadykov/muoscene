# File formats

Roadmap task 0.1-15, 28 November 2026. The files muoscene reads and writes, part by part and column by column, with what of each its authors have confirmed. Every number on this page is invented. The decisions behind the readers are in `docs/notes/tracks.md`.

## Text output of the FEDRA macro `evaluate_VT_nseg.C`

Read by `muoscene.io.tracks.read_fedra_text`. The first exposures of the scanning group came in this form on 27 November 2026: one text file per detector, written by a macro of the emulsion software FEDRA run in ROOT. What it looks like, with invented numbers:

```text
Load FEDRA libs

Processing evaluate_VT_nseg.C++...
N VTs = 3
VT param: X = 12000 Y = 15000 Theta = 0.12 Phi = 0.5
VT param: X = 12500.5 Y = 9000 Theta = 0.4 Phi = -2.1
VT param: X = 30000 Y = 21000 Theta = 1.3 Phi = 1.2
(int) 1
```

### Parts

| Part | Lines | What the reader does |
|---|---|---|
| header | every line before the first line that begins with `VT param:`; in the files at hand four lines: `Load FEDRA libs`, an empty line, `Processing evaluate_VT_nseg.C++...` and `N VTs = <n>` | keeps them in `TrackSet.meta["header"]`; the number of `N VTs = <n>` must be the number of track lines. A file without a header is read |
| track lines | `VT param: X = <x> Y = <y> Theta = <θ> Phi = <φ>`, one line per track, nothing between them | reads them; every line from the first to the last track line must have this form |
| after the tracks | empty lines and `(int) <n>`, which ROOT prints when the macro returns | keeps `(int) <n>` in `meta["trailer"]`; any other line is refused |

Numbers have a decimal point or a decimal comma — both kinds came — but one file has one kind. Runs of spaces and tabs separate words; lines end in LF or CR LF.

### Columns

| Column | Meaning | Unit | Confirmed |
|---|---|---|---|
| `X`, `Y` | where the track crosses the plane of the scan, in the axes of the scanning stage | µm, from the origin of the stage | the unit, by the group; the origin and the film whose plane it is, not |
| `Theta` | the tangent of the angle between the track and the normal of the film: sqrt(TX² + TY²) of the slopes TX = dx/dz, TY = dy/dz, as `EdbSegP::Theta()` in FEDRA | — | by the group (the letter of 27 November 2026); asked again, because a range in the name of one file reads more naturally as a tangent than as an angle and nothing else in the file says which. The reader takes `theta="radians"` for the other meaning |
| `Phi` | the angle of the track in the plane of the film from the x axis of the scan, counted towards its y axis: atan2(TY, TX), as `EdbSegP::Phi()` in FEDRA | rad | "from the x axis", by the group; the sense of the count is that of FEDRA, not confirmed |

So the slopes in the axes of the scan are tx = tan θ cos φ, ty = tan θ sin φ, with tan θ the column `Theta` itself. A line is a line: a track has no side, and the file does not say whether the muon went up or down.

### Not in the file

The scan-frame transform — how the axes of the stage lie on the detector and on which side of the film the scan looks — is not in the file and is given to the reader (`ScanFrame`), from the exposure log; without it a set stays in the frame of the scan and is not compared with a model. The area scanned, the dates of the exposure and the cuts already applied are given to the reader as arguments when they are known. The number of segments of a track, its first film and a quality of its fit are not in this format at all. Whether a correction for the efficiency is in the tracks, what the cuts are that the names of the files hint at, and the origin of the coordinates are asked of the group.

### What is refused

A ValueError that names the file and the line: an empty file; a file whose header promises tracks that are not there, or a number of them that the file does not have; a line among the track lines that is not a track line, an empty one included; labels out of their order; a value that is not a finite number; decimal points and decimal commas in one file; a carriage return without a line feed; a line after the tracks other than `(int) <n>`; a negative `Theta`, or with `theta="radians"` a `Theta` of π/2 or more.

## `TrackSet` files, `.npz`, version 1

Written by `TrackSet.save`, read by `TrackSet.load`. A NumPy `.npz` archive without pickled objects:

| Array | Shape | Meaning |
|---|---|---|
| `tx`, `ty` | (n,), float64 | slopes, in the frame of the set |
| `x_m`, `y_m` | (n,), float64 | positions, m; only when the set has them |
| `header` | (), text | JSON: `format` (`"muoscene.TrackSet"`), `version` (1), `n_tracks`, `frame`, `source`, `detector`, `positions` (true or false), `area_cm2`, `exposure_days`, `cuts`, `meta` |

Another format, another version, another set of arrays or a number of tracks that differs from `n_tracks` is refused. Version 2 will add the columns of a reader that brings them (segments, first film, quality); version 1 stays readable.

## Example

What `tests/unit/test_io_tracks.py::test_the_example_of_formats_md` runs, in an empty folder:

```python
from muoscene.data import bin_tracks
from muoscene.detectors import DirectionGrid
from muoscene.io.tracks import ScanFrame, read_fedra_text

lines = [
    "Load FEDRA libs",
    "",
    "Processing evaluate_VT_nseg.C++...",
    "N VTs = 3",
    "VT param: X = 12000 Y = 15000 Theta = 0.12 Phi = 0.5",
    "VT param: X = 12500.5 Y = 9000 Theta = 0.4 Phi = -2.1",
    "VT param: X = 30000 Y = 21000 Theta = 1.3 Phi = 1.2",
    "(int) 1",
]
with open("invented.txt", "w") as file:
    file.write("\n".join(lines) + "\n")

tracks = read_fedra_text("invented.txt", detector="film A")  # in the frame of the scan
binned = bin_tracks(tracks, DirectionGrid.square(), tan_max=1.0)
print(binned.summary())  # 3 tracks: 2 in the bins, 1 outside the circle (tan θ = 1.3)

# the film scanned from its other side, its reference point at X = 20 mm, Y = 15 mm of the stage
frame = ScanFrame(origin_um=(20000.0, 15000.0), y="-y", z="-z")
in_detector = read_fedra_text("invented.txt", detector="film A", scan_frame=frame)
```
