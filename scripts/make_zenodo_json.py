"""Build .zenodo.json from CITATION.cff. Authors live only in CITATION.cff; never edit .zenodo.json by hand.

Usage: python scripts/make_zenodo_json.py
"""

import json
import pathlib

import yaml

root = pathlib.Path(__file__).resolve().parents[1]
cff = yaml.safe_load((root / "CITATION.cff").read_text(encoding="utf-8"))
creators = []
for a in cff["authors"]:
    c = {"name": f"{a['family-names']}, {a['given-names']}"}
    if "orcid" in a:
        c["orcid"] = a["orcid"].rsplit("/", 1)[-1]
    if "affiliation" in a:
        c["affiliation"] = a["affiliation"]
    creators.append(c)
z = {
    "title": cff["title"],
    "upload_type": "software",
    "creators": creators,
    "license": cff["license"],
    "keywords": cff.get("keywords", []),
}
(root / ".zenodo.json").write_text(json.dumps(z, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
print("written .zenodo.json with", len(creators), "creator(s)")
