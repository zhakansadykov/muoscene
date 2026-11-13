# Changelog

## Unreleased

- Repository skeleton: package layout by layers, configuration of tests, linting and private-data guard.
- Draft conventions (`docs/conventions.md`) and detailed roadmap for 0.1 (`docs/roadmap.md`).
- Design review of the skeleton: `TrackSet` moved from `io` to `detectors.tracks` (layer 2), because `synth` and `data` (layer 4) could not import it from layer 5; `pipeline` module for orchestration shared by the CLI and Python users.
- Import rules tightened: relative imports are resolved by the layer test; optional dependencies must be imported lazily; the computing core must not import graphics; module docstrings must state their layer.
- Private-data guard: scans whole files for the marker, refuses survey and CAD formats outside `examples/` and `tests/`, runs in CI on all tracked files, covered by tests.
- Metadata: `.zenodo.json` checked against `CITATION.cff` by a test, a pre-commit hook and CI; licence, authors and repository URL checked across `pyproject.toml` and `CITATION.cff`.
- Tooling: ruff pinned to one version for pre-commit, CI and development; local hooks run in their own environment; CI runs the format check, keeps all matrix jobs (`fail-fast: false`) and cancels superseded runs.
- Mulder and Calzone spike (roadmap 0.1-02): results in `docs/notes/mulder.md`, reproducible by `docs/notes/spike/mulder_spike.py`.
- Architecture decision: hybrid on top of Mulder and Calzone. They become core dependencies (`mulder>=0.3.7,<0.4`, `calzone>=1.2,<2`); numba is dropped; the `mc` extra is removed. Stub modules renamed to the new design: `geometry.builder`, `geometry.pathtable`, `transport.table`, `transport.direct`, `flux.atmosphere`. A test confines backend imports to the wrapping modules.
- Conventions revised: kinetic energy at the detector, detector energy response, altitude cut-off, flux tables, equivalent opacity, bodies and anomalies. Roadmap for 0.1 rewritten.
