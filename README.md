# muoscene

[![tests](https://github.com/zhakansadykov/muoscene/actions/workflows/tests.yml/badge.svg)](https://github.com/zhakansadykov/muoscene/actions/workflows/tests.yml)

**Muography scene** — plan and check absorption muography of engineering-scale objects.

A site model — terrain and a near-vertical canyon wall, soil layers, piles and pile caps, niches cut into the wall — with detectors placed in it (nuclear emulsion stacks; straw telescopes) gives expected angular maps, significance versus exposure days with the emulsion track-density ceiling and fading, and the volume of ground actually seen from all detectors together. Measured tracks are compared with the model by the same code.

Geometry navigation and muon physics are not reimplemented: scenes are built with [Calzone](https://github.com/niess/calzone) and scanned by [Mulder](https://github.com/niess/mulder) (V. Niess; PUMAS inside), which also provides the reference flux, the atmosphere and energy loss. muoscene adds the description of engineering scenes, a fast flux path through tables built by Mulder, detector responses, planning statistics and comparison with data. The reasons and the numbers are in `docs/notes/mulder.md`.

> Status: skeleton (0.0.1.dev0). Nothing is computed yet. Roadmap — `docs/roadmap.md`.

## Install (development)

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
pre-commit install
pytest
```

Optional extras: `mesh`, `survey`, `viz3d`, `data`, `mceq`, `docs`.

## Layout

| Package | Responsibility |
|---|---|
| `units`, `config` | units and constants; configuration schema, YAML loading, registry |
| `site`, `materials`, `geometry`, `flux` | site frame and north; materials and equivalent-opacity weights; scene description, Calzone builder and scan; flux tables built by Mulder |
| `transport`, `detectors` | fast path from paths to flux, direct Mulder flux for checks, scattering; detector responses |
| `model` | forward model: ExpectedMap |
| `plan`, `synth`, `data`, `invert` | significance and planning; synthetic tracks; comparison with data; inversion (v2) |
| `io`, `viz`, `report`, `provenance`, `pipeline`, `cli` | files in and out; 3D scene and figures; reports; run folders; orchestration of the stages; commands |

Three import rules are checked by `tests/unit/test_layers.py`: a module imports only modules of its own or a lower layer (`src/muoscene/_layers.py`); packages installed with extras are imported inside functions, so every module imports with the core dependencies alone; the computing core (layers 0–4) never imports graphics; Mulder and Calzone are imported only by the modules that wrap them.

## Site data stay outside

Geometry and survey of real objects, their configurations and runs live outside this repository. Refused at commit and checked again in CI (`scripts/check_private.py`): files with the marker line `# PRIVATE-SITE`, names `*.private.*`, anything under a `private/` directory, files over 500 kB, and survey or CAD formats (`.dwg`, `.dxf`, `.las`, `.laz`, `.e57`, `.ifc`, `.rvt`, `.landxml`) outside `examples/` and `tests/`. Put the marker line at the top of every configuration of a real object. Examples use a synthetic canyon (`examples/canyon_demo`).

## Citation

See `CITATION.cff`. `.zenodo.json` is generated from it by `scripts/make_zenodo_json.py`; a test and a pre-commit hook keep the two in step.

## License

BSD 3-Clause.
