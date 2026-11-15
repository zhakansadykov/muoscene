"""Freeze Mulder and Calzone on the synthetic scenes of the spike into data/mulder_reference.json.

Run from the repository root (about a minute; the first run also compiles physics tables):
    python -m tests.regression.make_mulder_reference           write the file
    python -m tests.regression.make_mulder_reference --check   compare a fresh run with the frozen file:
        the ten largest relative changes; exit 1 if one exceeds mulder_cases.RTOL

Regenerate only on purpose — a new Mulder or Calzone release whose changes are understood — and
record the reason and the largest relative changes in CHANGELOG.md.
"""

import json
import pathlib
import platform
import sys

import numpy as np

from tests.regression import mulder_cases

PATH = pathlib.Path(__file__).resolve().parent / "data" / "mulder_reference.json"


def build():
    return {
        "about": "Numbers of the Mulder and Calzone spike (docs/notes/mulder.md) on synthetic scenes, "
        "continuous transport. K is the kinetic energy at the detector, E = K + m. Units: m, kg/m^3, "
        "g/cm^2 for opacity, GeV, m^-2 s^-1 sr^-1 (GeV^-1 for differential fluxes); physics tables in "
        "the units of Mulder.",
        "generated_with": {"python": platform.python_version(), "numpy": np.__version__}
        | mulder_cases.versions(),
        "parameters": mulder_cases.parameters(),
        "cases": {name: mulder_cases.compute(name) for name in mulder_cases.CASES},
    }


def largest_change(new, old):
    """Largest relative difference between two nested results; inf if their structure differs."""
    if isinstance(old, dict):
        if not isinstance(new, dict) or new.keys() != old.keys():
            return float("inf")
        return max((largest_change(new[k], old[k]) for k in old), default=0.0)
    if isinstance(old, str | bool) or old is None:
        return 0.0 if new == old else float("inf")
    a, b = np.asarray(new, float), np.asarray(old, float)
    if a.shape != b.shape:
        return float("inf")
    scale = np.maximum(np.abs(b), 1e-300)
    return float(np.max(np.where(a == b, 0.0, np.abs(a - b) / scale), initial=0.0))


def main(argv):
    data = build()
    if "--check" in argv:
        frozen = json.loads(PATH.read_text(encoding="utf-8"))
        if frozen["cases"] == data["cases"]:
            print("bit-identical")
            return 0
        changes = {k: largest_change(data["cases"].get(k), v) for k, v in frozen["cases"].items()}
        for name, change in sorted(changes.items(), key=lambda kv: -kv[1])[:10]:
            print(f"{change:.2e}  {name}")
        worst = max(changes.values())
        print(f"largest relative change {worst:.1e}, tolerance {mulder_cases.RTOL:.0e}")
        return 0 if worst <= mulder_cases.RTOL else 1
    PATH.write_text(json.dumps(data, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    print(f"written {PATH.name}: {len(data['cases'])} cases, {PATH.stat().st_size // 1024} kB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
