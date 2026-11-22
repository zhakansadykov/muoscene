# slab

Three horizontal films under 0.5, 2 and 20 m of flat ground at 1150 m: loam with 15 % of water, dry density 1.60 g/cm³, 1.84 g/cm³ with the water. The answer is known twice over: the opacity along every line of sight is the density times the depth over the cosine of the zenith angle, and the flux under the slab is the one Mulder gives.

## Run

With muoscene installed (until the installation page of roadmap 0.1-14: `pip install -e .` in the repository), from this folder:

```bash
python run.py
```

The first run builds three flux tables and compiles the physics tables of Mulder, up to a minute; later runs take a few seconds. It ends with `checks passed: 4 of 4` and the exit code 0; a check that fails is marked `FAILED` and the exit code is 1.

## What it checks

1. The opacity of every bin is density × depth / cos θ, averaged over the bin (met 8·10⁻¹⁶).
2. Under 0.5 m of soil and more every bin has a number above 0.2 GeV as well: the soil removes the cut-off of the reference flux.
3. The fast path agrees with Mulder along 256 random lines of sight per film (met within 0.026 %).
4. The flux along the vertical and along 45°, from one narrow bin around each, agrees with Mulder along the same lines (met within 0.024 %). These are the directions of the numbers frozen for Mulder in the tests of the program: with the loam of those numbers the example reproduces them within 0.019 % (`tests/unit/test_example_slab.py`).

Above 1 GeV the vertical flux is 64.43, 52.72 and 10.85 muons m⁻² s⁻¹ sr⁻¹ under 0.5, 2 and 20 m; above 0.2 GeV 96.59, 74.07 and 12.74.

## Change it

`SOIL`, `GROUND` and `DETECTORS` at the top of `run.py`: the soil, the size of the ground, and the name, position (ENU, m), depth and exposure of every film. A film stays inside the ground; at 45° from it a line runs as far sideways as the film is deep, so the ground must reach further than that.
