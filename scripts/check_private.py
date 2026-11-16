"""Refuse to commit site data. Run by pre-commit on staged files and by CI on all tracked files.

A path is refused when

- a line of the file consists of the marker `# PRIVATE-SITE` (put it at the top of every
  configuration of a real object);
- its name matches `*.private.*` or it lies under a directory named `private`;
- it is larger than 500 kB;
- it lies outside `examples/` and `tests/`, where only synthetic geometry lives, and its type is
  not on the short list of what the repository is made of: source code, Markdown, YAML, TOML,
  JSON, plain text, and figures under `docs/`. Drawings, scans, point clouds, terrain models and
  spreadsheets are refused whatever their format; survey and CAD formats (.dwg, .dxf, .las,
  .laz, .e57, .ifc, .rvt, .landxml) are named in the message.

A new kind of file is added to `ALLOWED_SUFFIXES` on purpose, never by habit.

Site data of real objects live outside the repository (e.g. ~/Data/muoscene-<site>).

Usage: python scripts/check_private.py PATH [PATH ...]
"""

import pathlib
import sys

MARKER = b"# PRIVATE-SITE"
MAX_BYTES = 500 * 1024
SITE_SUFFIXES = frozenset({".dwg", ".dxf", ".las", ".laz", ".e57", ".ifc", ".rvt", ".landxml"})
SYNTHETIC_DIRS = frozenset({"examples", "tests"})
ALLOWED_SUFFIXES = frozenset({".py", ".md", ".toml", ".yaml", ".yml", ".json", ".cff", ".txt", ".typed"})
ALLOWED_NAMES = frozenset({"LICENSE", ".gitignore"})
FIGURE_SUFFIXES = frozenset({".png", ".svg"})
FIGURE_DIR = "docs"


def allowed_type(path: pathlib.Path, parts: tuple[str, ...]) -> bool:
    """Whether a file outside the synthetic directories is of a kind the repository is made of."""
    suffix = path.suffix.lower()
    if path.name in ALLOWED_NAMES or suffix in ALLOWED_SUFFIXES:
        return True
    return suffix in FIGURE_SUFFIXES and parts[:1] == (FIGURE_DIR,)


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
    synthetic = parts[:1] != () and parts[0] in SYNTHETIC_DIRS
    if path.suffix.lower() in SITE_SUFFIXES and not synthetic:
        out.append(f"survey or CAD format {path.suffix} outside examples/ and tests/")
    elif not synthetic and not allowed_type(path, parts):
        out.append(f"file type {path.suffix or path.name!r} is not allowed outside examples/ and tests/")
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
