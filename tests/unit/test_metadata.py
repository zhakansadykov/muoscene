"""Package metadata stay consistent across pyproject.toml, CITATION.cff and .zenodo.json."""

import importlib.util
import pathlib
import re
import tomllib

import yaml

import muoscene

ROOT = pathlib.Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("make_zenodo_json", ROOT / "scripts" / "make_zenodo_json.py")
make_zenodo_json = importlib.util.module_from_spec(spec)
spec.loader.exec_module(make_zenodo_json)

PYPROJECT = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
CFF = yaml.safe_load((ROOT / "CITATION.cff").read_text(encoding="utf-8"))


def test_zenodo_json_follows_citation():
    expected = make_zenodo_json.render(make_zenodo_json.build(CFF))
    assert (ROOT / ".zenodo.json").read_text(encoding="utf-8") == expected


def test_licence_is_the_same_everywhere():
    assert PYPROJECT["project"]["license"] == CFF["license"]
    assert "BSD 3-Clause" in (ROOT / "LICENSE").read_text(encoding="utf-8")


def test_authors_are_the_same_everywhere():
    in_pyproject = {a["name"] for a in PYPROJECT["project"]["authors"]}
    in_citation = {f"{a['given-names']} {a['family-names']}" for a in CFF["authors"]}
    assert in_pyproject == in_citation


def test_repository_url_is_the_same_everywhere():
    assert PYPROJECT["project"]["urls"]["Repository"] == CFF["repository-code"]


def test_version_is_valid():
    assert re.fullmatch(r"\d+\.\d+\.\d+((a|b|rc)\d+)?(\.dev\d+)?", muoscene.__version__)
