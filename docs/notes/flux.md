# Flux tables

Roadmap task 0.1-05, 16 November 2026, Mulder 0.3.7. Module `muoscene.flux`; conventions, revision 4. Every number below comes from the named test of `tests/unit/test_flux_tables.py` or from a frozen case of `tests/regression/data/mulder_reference.json`; the headline ones are in `docs/numbers.md`. Reference material of the tests: the loam of the spike; site at 1150 m.

## What a table is

A `FluxTable` answers two questions for a detector: the differential flux at a kinetic energy K, a sine of elevation and an opacity X of the reference material (`differential`), and that flux integrated over the energy response of the detector (`integral`). It holds

- the open-sky differential flux at the altitude of the detector on a grid of energy and sine of elevation: 100 energies per decade, sines from 0.25 (14.5°) to 1 in steps of 0.005;
- the CSDA range R(K) and the stopping power S(K) of the reference material, 200 energies per decade from 1 MeV.

The opacity is not a grid dimension. In continuous transport the energy a muon loses depends only on the opacity it crosses, so

    φ(K, X) = φ_sky(K_top) · S(K_top) / S(K),    R(K_top) = R(K) + X,

and the muons above K under X are the muons above K_top in the open sky: an integral above a threshold needs no Jacobian. A grid in opacity filled by Mulder under a slab, as revision 3 of the conventions had it, would repeat this information, and it cannot be interpolated near the altitude cut-off, where the flux is singular.

`MulderTable` builds the table with Mulder: atmosphere model, optionally the geomagnetic field of a date, the GCCLY15 reference at sea level. One table takes 2 s on a Linux VM with 4 cores and is 0.8 MB on disk. A reference at the site altitude (roadmap 0.1-05b) will enter through the argument `reference`, which until then accepts only the sea-level one.

## Check against Mulder

Every new table is compared with Mulder under a dense slab of its reference material (10⁵ kg/m³, top at the altitude of the table) at 64 random energies up to 1 TeV, elevations and opacities, and is refused if the largest difference exceeds 0.2 %. The result is kept in `table.meta["check"]`. A grid with two elevation nodes 30° apart fails it (`test_a_grid_too_coarse_is_refused_by_the_check`).

| Comparison | Result | Test |
|---|---|---|
| differential flux, 16 150 random points up to 1 TeV | rms 0.0035 %, largest 0.045 % | `test_table_is_within_02_percent_of_mulder_at_random_points` |
| the same, 3850 points from 1 to 10 TeV | rms 0.044 %, largest 0.79 % | the same |
| integral above 0.2 and 1 GeV, 12 random directions each | within 6·10⁻⁵ | `test_integrals_are_those_of_mulder_under_a_slab` |
| open sky and thin soil against the frozen numbers | within 7·10⁻⁵; 1.5·10⁻⁴ under half a metre of soil | `test_vertical_flux_at_sea_level`, `test_cutoff_survives_under_thin_soil` |

Above 1 TeV the comparison is limited by Mulder: under 2.5·10⁴ g/cm² its flux between 1.2 and 3.5 TeV deviates from a smooth curve by up to 0.34 %, while the open-sky flux it starts from is smooth to 0.0003 % and the table to 0.09 % (`test_mulder_ripples_above_1_tev_under_a_thick_slab`). Muons above 1 TeV are 0.09 % of the flux above 1 GeV under that slab, so the effect on an integral is 3·10⁻⁶. The check of a new table therefore stops at 1 TeV.

## Altitude cut-off

With the sea-level reference above sea level the open-sky flux exists only above K_cut(elevation), the energy lost in the air down to sea level, and rises towards it as (K − K_cut)^(−0.4) between 100 and 10 MeV above it. The table finds K_cut for every elevation node by bisection, stores it and takes the energy nodes of each elevation as K_cut + d on one grid of d that starts 10 MeV above zero. Under an opacity the cut-off left at the detector follows from the range table.

- At 1150 m the table gives 0.28, 0.38, 0.53 GeV at 90°, 45°, 30° and 1.06 GeV at 14.5°; at 3340 m 0.69, 0.97, 1.40 GeV (`test_cutoff_numbers_of_the_spike`).
- Soil shortens it: 0.12 GeV is left under 0.5 m by the vertical, 0.16 GeV at 45° (frozen `thin_slab`), and nothing under more than 138 g/cm² by the vertical or 569 g/cm² at 14.5° (`test_the_table_holds_only_energies_above_the_cutoff`).
- A request below `energy_min(sin_elevation, opacity)`, the cut-off plus 10 MeV, is refused with the cut-off in the message; it is neither extrapolated nor returned as zero. At 1150 m a detector threshold of 1 GeV holds in the open sky from 15.7° of elevation; one of 0.2 GeV holds only under soil.
- 430 m below sea level there is no cut-off and the flux is lower: 56.69 against 60.50 m⁻² s⁻¹ sr⁻¹ above 1 GeV of total energy (`test_below_sea_level_the_flux_is_lower_and_has_no_cutoff`).

## Integral over the energy response

`integral(response, sin_elevation, opacity_g_cm2)` always takes the response: an `EnergyResponse` is a list of energies and efficiencies, zero below the first, linear in between, constant above the last, and `EnergyResponse.threshold(1.0)` is a sharp threshold. The integral is a trapezoid sum in ln K that starts at the first energy of the response, four trapezoids to an interval of the grid. A sum that starts at a node and gives it a whole step, as the spike took it, is 2.85 % high at sea level and 3.29 % at 1150 m for 61 nodes from 1 GeV (`test_a_sum_from_the_first_node_is_two_to_three_percent_high`). A table ends at 10 TeV at the detector; the muons above are 3·10⁻⁶ of those above 1 GeV under 3·10⁴ g/cm² (`test_muons_above_the_table_do_not_matter`).

## What a table leaves out

| Left out | Size | Test |
|---|---|---|
| the table of a detector is at its altitude: the soil above it comes with as much air that is not there | −0.05 to −0.10 % under 10–30 m of soil, every muon counted; −0.07 to −0.08 % above 0.2 and 1 GeV under 20 m | `test_slab_numbers_of_the_spike`, `test_detector_thresholds_under_soil` |
| decay inside the opacity | +0.03 to +0.05 % under 10–30 m of soil of 1.8 g/cm³ | `test_slab_numbers_of_the_spike` |
| latitude and longitude without a geomagnetic field | 1.6·10⁻⁶ between the equator and 43° north | `test_latitude_matters_little_without_a_field` |
| azimuth: with a date a table is made for one azimuth | east against west at 45°: +0.27 % above 1 GeV | `test_geomagnetic_field_of_a_date` |

Not in any table of version 0.1.0: production of muons in the air between sea level and the site, and electrons and positrons (roadmap 0.1-05b); scattering and energy-loss fluctuations (version 0.3). Opacities beyond 3·10⁴ g/cm² are refused.

## Cache

Tables are cached as `.npz` files, without pickle, in `$MUOSCENE_CACHE/flux`, else `$XDG_CACHE_HOME/muoscene/flux`, else `~/.cache/muoscene/flux`; `cache_dir` overrides the directory and `cache=False` turns the cache off. The file name is a hash of the altitude, position, azimuth, atmosphere, date, reference, definition of the material, grids, version of Mulder and version of the algorithm. A file that cannot be read, was made from other inputs or does not make a table is rebuilt (`test_cache_by_hash_of_the_inputs`, `test_a_damaged_cache_file_is_rebuilt`).

## Example

A detector 20 m below a brow at 1150 m, in loam of 1.60 g/cm³ dry density and 15 % of water (`test_example_of_the_note`):

```python
from muoscene import materials as M
from muoscene.flux import EnergyResponse, MulderTable

soil = M.loam(1.60, 15.0)
table = MulderTable(1130.0, soil.material)
opacity = 20.0 * soil.bulk_density_g_cm3 * 100  # 3680 g/cm² by the vertical
table.integral(EnergyResponse.threshold(1.0), 1.0, opacity)  # 10.84 m⁻² s⁻¹ sr⁻¹ above 1 GeV
table.integral(EnergyResponse.threshold(0.2), 1.0, opacity)  # 12.73 above 0.2 GeV: holds under soil
table.integral(EnergyResponse.threshold(0.2), 1.0, 0.0)  # ValueError: the cut-off is 0.28 GeV
```

## Run on another machine

```bash
pytest tests/unit/test_flux_tables.py --durations=3    # the setup of the first test is one table: 2 s
python -c "
import time, numpy as np
from muoscene import materials as M
from muoscene.flux import EnergyResponse, MulderTable
table = MulderTable(1150.0, M.STANDARD_ROCK, cache=False)
s = np.random.default_rng(0).uniform(0.3, 1.0, 10**6)
x = np.random.default_rng(1).uniform(0.0, 3.0e4, 10**6)
t = time.perf_counter()
table.integral(EnergyResponse.threshold(1.0), s, x)
print(f'{1e6 / (time.perf_counter() - t):.1e} directions/s')
"
```

On the Linux VM: 2.5·10⁶ directions/s. Not measured on macOS yet.
