"""The guard that keeps site data out of the repository."""

import importlib.util
import pathlib

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("check_private", ROOT / "scripts" / "check_private.py")
check_private = importlib.util.module_from_spec(spec)
spec.loader.exec_module(check_private)

MARKER = check_private.MARKER.decode()


@pytest.fixture
def repo(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    return tmp_path


def write(path: str, text: str = "x: 1\n") -> pathlib.Path:
    p = pathlib.Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(text, encoding="utf-8")
    return p


def test_clean_file_passes(repo):
    assert check_private.problems(write("examples/canyon_demo/site.yaml")) == []


@pytest.mark.parametrize("line", [MARKER, f"  {MARKER}  ", f"{MARKER}\r"])
def test_marker_line_is_refused(repo, line):
    p = write("site.yaml", f"{line}\nsite: x\n")
    assert any("marker" in m for m in check_private.problems(p))


def test_marker_far_from_the_top_is_refused(repo):
    p = write("site.yaml", "# comment\n" * 2000 + MARKER + "\n")
    assert any("marker" in m for m in check_private.problems(p))


def test_marker_inside_a_sentence_passes(repo):
    p = write("README.md", f"Files with the marker line `{MARKER}` are refused.\n")
    assert check_private.problems(p) == []


def test_private_name_is_refused(repo):
    assert check_private.problems(write("site.private.yaml"))


def test_private_directory_is_refused(repo):
    assert check_private.problems(write("data/private/site.yaml"))


def test_large_file_is_refused(repo):
    p = write("big.bin", "0" * (check_private.MAX_BYTES + 1))
    assert any("larger" in m for m in check_private.problems(p))


@pytest.mark.parametrize("name", ["survey.dwg", "wall.LAS", "pier.ifc", "docs/terrain.dxf"])
def test_survey_and_cad_formats_are_refused(repo, name):
    assert any("survey or CAD" in m for m in check_private.problems(write(name)))


@pytest.mark.parametrize("name", ["examples/canyon_demo/wall.dxf", "tests/data/points.las"])
def test_synthetic_geometry_passes(repo, name):
    assert check_private.problems(write(name)) == []


def test_missing_path_does_not_crash(repo):
    assert check_private.problems(pathlib.Path("deleted.yaml")) == []


def test_main_exit_codes(repo, capsys):
    clean = write("a.yaml")
    private = write("b.yaml", f"{MARKER}\n")
    assert check_private.main([str(clean)]) == 0
    assert check_private.main([str(clean), str(private)]) == 1
    assert "refused: b.yaml" in capsys.readouterr().out
