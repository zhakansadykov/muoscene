"""Build .zenodo.json from CITATION.cff.

Title, authors, licence and keywords live only in CITATION.cff; never edit .zenodo.json by hand.

Usage:
    python scripts/make_zenodo_json.py           write .zenodo.json
    python scripts/make_zenodo_json.py --check   exit 1 if .zenodo.json is out of date
"""

import json
import pathlib
import sys

import yaml

ROOT = pathlib.Path(__file__).resolve().parents[1]


def build(cff: dict) -> dict:
    """Zenodo metadata record for a parsed CITATION.cff."""
    creators = []
    for author in cff["authors"]:
        creator = {"name": f"{author['family-names']}, {author['given-names']}"}
        if "orcid" in author:
            creator["orcid"] = author["orcid"].rsplit("/", 1)[-1]
        if "affiliation" in author:
            creator["affiliation"] = author["affiliation"]
        creators.append(creator)
    record = {
        "title": cff["title"],
        "upload_type": "software",
        "creators": creators,
        "license": cff["license"],
        "keywords": list(cff.get("keywords", [])),
    }
    if "abstract" in cff:
        record["description"] = cff["abstract"]
    return record


def render(record: dict) -> str:
    return json.dumps(record, ensure_ascii=False, indent=2) + "\n"


def main(argv: list[str]) -> int:
    cff = yaml.safe_load((ROOT / "CITATION.cff").read_text(encoding="utf-8"))
    text = render(build(cff))
    target = ROOT / ".zenodo.json"
    if "--check" in argv:
        if not target.exists() or target.read_text(encoding="utf-8") != text:
            print(".zenodo.json is out of date: run python scripts/make_zenodo_json.py")
            return 1
        return 0
    target.write_text(text, encoding="utf-8")
    print(f"written .zenodo.json with {len(cff['authors'])} creator(s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
