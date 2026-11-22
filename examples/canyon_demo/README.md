# canyon_demo

A synthetic canyon at 1150 m with a pier on each bank, a void and a zone of decompacted soil near the piers, and six emulsion detectors on the walls: four in niches cut into the rock and two on pins in front of it. Every size is invented and differs from any real site.

## The scene

A section across the canyon, looking along it (towards the azimuth of 30°, the `+y` axis of the scene); z is up from the edge of the left brow, x across from the left wall. Not to scale.

```
  z, m        left bank                          canyon                right bank
  +6.5                                                         ________[cap]_______  level brow
  +5.0                                                    edge /         |||         (rises 1.5 m over 8 m)
     0  ______________[cap]________________ edge          |   (decompaction) |||
                       |||    (void)       |  air          |                 |||
  -7.0                 |||                 |           R12 n|
  -8.0  - - - - - - - -|||- - - layer 1 - -|               |- - - - - - - - - - - - -
 -12.0                 |||           L12n n|w L12w         |
 -17.0                                     |           R22 w|
 -18.0  - - - - - - - - - - - - layer 2 - -|               |- - - - - - - - - - - - -
 -22.0                               L22n n|w L22w         |
 -40.0                                     |_______________|  floor, 30 m wide
  x, m   -300   -14   -8   -6              0               30   36  39   45     330
```

- Ground: loam with 15 % of water, dry density 1.60 g/cm³ (1.84 g/cm³ with the water) down to z = −8 m, 1.70 g/cm³ down to −18 m, 1.80 g/cm³ below; the model reaches 300 m behind each wall, 300 m along the canyon either way and down to −150 m.
- Left pier: a pile cap 6 × 6 × 1.5 m level with the surface, its near face 8 m from the wall, on 3 × 3 piles 0.8 m across and 8 m long at 2.2 m. Right pier: a cap 6 × 6 × 2 m, 9 m from the wall, on 3 × 3 piles 1.0 m across and 10 m long at 2.0 m. Concrete of 2.30 g/cm³.
- Void: a sphere 2 m across, 6 m behind the left wall and 6 m below the brow. Decompaction: a sphere 3 m across with 20 % less density, 6 m behind the right wall and 5 m below its edge.
- Detectors: one film of 120 cm² each (six films in a stack), accepted within 45° of the normal, efficiency 0.9.

| Detector | Bank | Along, m | Below the edge, m | Mount | Normal: azimuth, tilt | Days |
|---|---|---|---|---|---|---|
| L12n | left | −1 | 12 | niche | 300°, 20° into the massif | 60 |
| L12w | left | +1 | 12 | wall | 300°, 20° | 60 |
| L22n | left | −1 | 22 | niche | horizontal | 60 |
| L22w | left | +1 | 22 | wall | horizontal | 60 |
| R12 | right | 0 | 12 | niche | 120°, 20° into the massif | 60 |
| R22 | right | 0 | 22 | wall | 300°, 30° across the canyon | 90 |

A niche is 0.5 m deep into the wall, 0.3 m wide and 0.3 m high; the film lies in its middle. A film on the wall hangs 0.15 m in front of it; the pins, the cassette and any cover are not modelled.

## Run

With muoscene installed (until the installation page of roadmap 0.1-14: `pip install -e .` in the repository), from this folder:

```bash
python run.py
```

The first run builds four flux tables and compiles the physics tables of Mulder, a minute or two; later runs take a few seconds. Nothing is written except the cache of the tables, in the cache folder of the user. The script refers to nothing outside its own file: copy the folder anywhere and run it there.

## Change it

Everything is at the top of `run.py`; nothing under the line "nothing to edit below" needs changing.

- **Move a detector:** change its line in `DETECTORS` — the bank, the distance along the canyon, the depth below the edge of the brow. Its niche moves with it.
- **Pins instead of a niche:** change `"niche"` to `"wall"` in its line.
- **Turn it:** the azimuth of the normal is from true north, clockwise: 300 points into the left massif, and from the right wall across the canyon; 120 into the right massif. A tilt up to 30° is computed in full. Tilted 40° into the massif, 154 bins look below 14.5° of elevation, where the flux table ends, and have no number; from 44° the program refuses the detector, because its lowest lines run along the wall out of the model.
- **Other anomalies or densities:** `ANOMALIES` and `HYPOTHESES`; sizes of the canyon, the brows, the layers and the piers in `CANYON`, `BROWS`, `LAYERS` and `PIERS`.

A wrong line is refused with one line that names the detector, for example `canyon_demo: detector L12w: mount is 'niche' or 'wall', not 'pins'`, and the exit code 2.

## What the numbers mean

First the detectors as the program read them: their place in ENU (x east, y north, z up from the site, m), their altitude, the direction of their normal and the lowest elevation they look at. Then, for a threshold of 1 GeV and of 0.2 GeV — the kinetic energy a muon needs at the detector to leave a track:

- **bins with a number:** the map has 40 × 40 bins in the slopes of the film, 1324 of them within its cone. A bin has no number where the flux table cannot give one: below the cut-off of the reference flux (the open sky and thin soil above 0.2 GeV — the reference at the site altitude is roadmap 0.1-05b), below 14.5° of elevation, or behind more than 3·10⁴ g/cm². Such a bin is left out, never counted as zero.
- **tracks in them:** the expected number of muon tracks in the bins with a number, for the exposure of the detector, with no anomaly.
- **fast path against Mulder:** every run computes the flux along 256 random lines of sight per detector with Mulder itself and compares; above 0.5 % the run is refused.
- **anomaly against none:** how many tracks the void or the decompaction adds, in how many bins the counts change by more than 1 %, and the largest change with the direction of its bin. Whether that is enough to be seen in data is significance, roadmap 0.1-09.
- **niche against wall:** the two films of a pair side by side: their tracks, what the anomaly adds, their bins with a number and their open-sky bins.

Above 1 GeV (`python run.py`; the counts are checked by `tests/unit/test_example_canyon_demo.py`):

| Detector | Tracks | Void | Decompaction | Largest change in a bin |
|---|---|---|---|---|
| L12n | 1.771·10⁶ | +1947 | — | +12.6 % |
| L12w | 2.147·10⁶ | +1614 | — | +12.1 % |
| L22n | 2.293·10⁶ | +360 | — | +10.6 % |
| L22w | 2.751·10⁶ | +345 | — | +10.8 % |
| R12 | 1.733·10⁶ | — | +1189 | +3.6 % |
| R22 | 4.943·10⁶ (90 days) | +16 | — | +3.0 % |

On the wall a film collects about 20 % more tracks than in a niche — part of its cone is open sky (382 and 662 bins against 0 and 52) — while the void changes its bins by about as much. Above 0.2 GeV the open sky has no number at 1150 m: 942 and 622 bins of the films on the wall have one, against 1292 and 990 in the niches. The fast path agrees with Mulder within 0.066 % above 1 GeV and 0.14 % above 0.2 GeV, on 256 lines per detector.

## Not modelled

Scattering of muons (version 0.3), the measured energy response of the emulsion (a sharp threshold here), the soft component of the open sky (0.1-05b), fading of the latent image and the track density a scan tolerates (0.1-09 and 0.3), the pins and cassettes of a film on the wall.
