"""Refuse to commit site data. Run by pre-commit on staged files.

Rejects: files with a line consisting of the marker `# PRIVATE-SITE` in the first 4 kB,
names matching *.private.*, anything under a `private/` directory, and files larger than
500 kB. Site data of real objects live outside the repository (e.g. ~/Data/muoscene-<site>).
"""

import pathlib
import sys

MARKER = b"# PRIVATE-SITE"
MAX_BYTES = 500 * 1024


def problems(path: pathlib.Path) -> list[str]:
    out = []
    parts = path.parts
    if "private" in parts:
        out.append("inside a private/ directory")
    if ".private." in path.name:
        out.append("name matches *.private.*")
    if path.is_file():
        if path.stat().st_size > MAX_BYTES:
            out.append(f"larger than {MAX_BYTES // 1024} kB")
        try:
            head = path.read_bytes()[:4096].splitlines()
            if any(line.strip() == MARKER for line in head):
                out.append("contains the PRIVATE-SITE marker")
        except OSError:
            pass
    return out


def main(argv: list[str]) -> int:
    if pathlib.Path(__file__).name in [pathlib.Path(a).name for a in argv]:
        argv = [a for a in argv if pathlib.Path(a).name != pathlib.Path(__file__).name]
    bad = 0
    for a in argv:
        for msg in problems(pathlib.Path(a)):
            print(f"refused: {a}: {msg}")
            bad += 1
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
