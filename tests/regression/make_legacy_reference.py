"""Freeze the preliminary calculation on the synthetic canyon into data/legacy_reference.json.

Run from the repository root:
    python -m tests.regression.make_legacy_reference           write the file
    python -m tests.regression.make_legacy_reference --check   exit 1 unless a fresh run is bit-identical
"""

import json
import pathlib
import platform
import sys

import numpy as np

from tests.regression import legacy_cases

PATH = pathlib.Path(__file__).resolve().parent / "data" / "legacy_reference.json"
ORIGINAL_SHA256 = {
    "sens.py": "2f93533b628a8ac34733cce62121cb762d36a723576de9c5b2dc427e9001f18d",
    "scheme.py": "e3febef1bee2530d85fb97f87760cce50369f42e66df817b00ab17847ffc4776",
    "tilt.py": "ea5d0583b1d90452bec69ab3b90dbe3b669bd7615f7e481d29eb49a475e1abd5",
}


def build():
    return {
        "about": "Preliminary sensitivity calculation of the 2027 grant application, ported to "
        "tests/regression/legacy.py, on a synthetic canyon. Units: m, g/cm^3, g/cm^2, GeV (total "
        "energy), s, cm^2; intensities in cm^-2 s^-1 sr^-1; densities of tracks per cm^2.",
        "ported_from_sha256": ORIGINAL_SHA256,
        "generated_with": {"python": platform.python_version(), "numpy": np.__version__},
        "parameters": legacy_cases.parameters(),
        "cases": {name: legacy_cases.compute(name) for name in legacy_cases.CASES},
    }


def render(data):
    return json.dumps(data, indent=1, sort_keys=True) + "\n"


def main(argv):
    data = build()
    if "--check" in argv:
        frozen = json.loads(PATH.read_text(encoding="utf-8"))
        same = frozen["cases"] == data["cases"]
        print("bit-identical" if same else "differs from the frozen file")
        return 0 if same else 1
    PATH.write_text(render(data), encoding="utf-8")
    print(f"written {PATH.name}: {len(data['cases'])} cases, {PATH.stat().st_size // 1024} kB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
