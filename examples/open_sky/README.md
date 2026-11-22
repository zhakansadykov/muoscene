# open_sky

Two films under the open sky at 1150 m, on a flat ground: one horizontal, one tilted by 60° to the south. Nothing stands between them and the sky, so the answer is known: the expected map is the flux table itself. This is also the first check of an installation.

## Run

With muoscene installed (until the installation page of roadmap 0.1-14: `pip install -e .` in the repository), from this folder:

```bash
python run.py
```

The first run builds a flux table and compiles the physics tables of Mulder, up to a minute; later runs take about a second. It ends with `checks passed: 5 of 5` and the exit code 0; a check that fails is marked `FAILED` and the exit code is 1.

## What it checks

1. In every bin with a number the counts are the flux of the table times the projected area, the solid angle, the efficiency and the exposure, to round-off (met 7·10⁻¹⁶).
2. Every bin the films accept sees the open sky.
3. The fast path agrees with Mulder, computed along 256 random lines of sight per film (met within 0.025 %).
4. Above 0.2 GeV no bin has a number: the sea-level reference flux, carried up to 1150 m, has no muons below its cut-off, 0.284 GeV at the vertical. A reference at the site altitude is roadmap 0.1-05b.
5. The tilted film, above 1 GeV: nothing is counted below the horizon (476 bins are not accepted); 270 bins look below 14.5° of elevation, where the flux table ends, and 34 below the elevation where 1 GeV meets the cut-off; those 304 bins have no number. 820 bins have one.

It also prints the cut-off and the vertical flux above 1 GeV of the table: 69.22 muons m⁻² s⁻¹ sr⁻¹ at 1150.5 m. The numbers are checked by `tests/unit/test_example_open_sky.py`.

## Change it

`DETECTORS` at the top of `run.py`: name, position in ENU (x east, y north, m), height above the ground, azimuth and tilt of the normal, exposure. A film must stay below 1 m above the ground, the top of the scene.
