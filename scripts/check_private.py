"""Refuse to commit site data. Run by pre-commit on staged files and by CI on all tracked files.

A path is refused when

- a line of the file consists of the marker `# PRIVATE-SITE` (put it at the top of every
  configuration of a real object);
- its name matches `*.private.*` or it lies under a directory named `private`;
- it is larger than 500 kB;
- it has a survey or CAD suffix (.dwg, .dxf, .las, .laz, .e57, .ifc, .rvt, .landxml) outside
  `examples/` and `tests/`, where only synthetic geometry lives.

Site data of real objects live outside the repository (e.g. ~/Data/muoscene-<site>).

Usage: python scripts/check_private.py PATH [PATH ...]
"""

import pathlib
import sys

MARKER = b"# PRIVATE-SITE"
MAX_BYTES = 500 * 1024
SITE_SUFFIXES = frozenset({".dwg", ".dxf", ".las", ".laz", ".e57", ".ifc", ".rvt", ".landxml"})
SYNTHETIC_DIRS = frozenset({"examples", "tests"})


def problems(path: pathlib.Path) -> list[str]:
    """Reasons why `path` must not be committed; empty if it may."""
    try:
        parts = path.resolve().relative_to(pathlib.Path.cwd().resolve()).parts
    except ValueError:
        parts = path.parts
    out = []
    if "private" in parts:
        out.append("inside a private/ directory")
    if ".private." in path.name:
        out.append("name matches *.private.*")
    if path.suffix.lower() in SITE_SUFFIXES and not SYNTHETIC_DIRS.intersection(parts):
        out.append(f"survey or CAD format {path.suffix} outside examples/ and tests/")
    if path.is_file():
        if path.stat().st_size > MAX_BYTES:
            out.append(f"larger than {MAX_BYTES // 1024} kB")
        else:
            try:
                lines = path.read_bytes().splitlines()
            except OSError:
                lines = []
            if any(line.strip() == MARKER for line in lines):
                out.append("contains the PRIVATE-SITE marker")
    return out


def main(argv: list[str]) -> int:
    refused = 0
    for arg in argv:
        for message in problems(pathlib.Path(arg)):
            print(f"refused: {arg}: {message}")
            refused += 1
    return 1 if refused else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
