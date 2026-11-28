"""`muoscene.io.tracks`: the text output of the FEDRA macro `evaluate_VT_nseg.C` and the scan-frame transform.

Roadmap 0.1-15. Every file is written by the test into a temporary folder from invented numbers,
in the layout of `docs/formats.md`: no line of a real exposure is here. A line that cannot be
read is refused with its number, never skipped.
"""

import hashlib
import itertools
import math
import pathlib
import re

import numpy as np
import pytest

from muoscene.io.tracks import ScanFrame, read_fedra_text

HEADER = ["Load FEDRA libs", "", "Processing evaluate_VT_nseg.C++..."]
FORMATS_MD = pathlib.Path(__file__).resolve().parents[2] / "docs" / "formats.md"
TRACKS_MD = pathlib.Path(__file__).resolve().parents[2] / "docs" / "notes" / "tracks.md"


def invented(n=50, seed=11, theta_max=1.2):
    """Columns X, Y (µm), Theta, Phi of invented tracks, as text with six significant digits, like ROOT."""
    rng = np.random.default_rng(seed)
    columns = (
        rng.uniform(1.0e3, 9.0e4, n),
        rng.uniform(1.0e3, 9.0e4, n),
        rng.uniform(0.0, theta_max, n),
        rng.uniform(-math.pi, math.pi, n),
    )
    return [[f"{value:g}" for value in row] for row in zip(*columns, strict=True)]


def lines_of(rows, header=True, declared=True, trailer=True):
    out = [*HEADER, f"N VTs = {len(rows)}"] if header else []
    if header and not declared:
        out.pop()
    out += [f"VT param: X = {x} Y = {y} Theta = {t} Phi = {p}" for x, y, t, p in rows]
    return out + (["(int) 1"] if trailer else [])


def write(path, lines, ending="\n"):
    path.write_bytes((ending.join(lines) + ending).encode("ascii"))
    return path


def numbers(rows):
    return np.array([[float(v.replace(",", ".")) for v in row] for row in rows]).T


def test_a_file_as_the_macro_writes_it(tmp_path):
    rows = invented()
    path = write(tmp_path / "tracks.txt", lines_of(rows))
    tracks = read_fedra_text(path, detector="film A")
    x, y, theta, phi = numbers(rows)
    assert len(tracks) == len(rows) and tracks.frame == "scan" and tracks.source == "data"
    assert np.array_equal(tracks.x_m, x * 1e-6) and np.array_equal(tracks.y_m, y * 1e-6)
    assert np.array_equal(tracks.tx, theta * np.cos(phi)) and np.array_equal(tracks.ty, theta * np.sin(phi))
    assert np.allclose(np.hypot(tracks.tx, tracks.ty), theta, rtol=1e-15, atol=0)
    meta = tracks.meta
    assert meta["header"] == [*HEADER, f"N VTs = {len(rows)}"] and meta["declared_tracks"] == len(rows)
    assert meta["trailer"] == ["(int) 1"] and meta["first_track_line"] == 5
    assert meta["decimal"] == "." and meta["line_ends"] == "LF" and meta["theta"] == "tangent"
    assert meta["scan_frame"] is None and meta["exposure"] is None
    raw = path.read_bytes()
    assert meta["file"] == {
        "name": "tracks.txt",
        "bytes": len(raw),
        "sha256": hashlib.sha256(raw).hexdigest(),
    }


def test_comma_line_ends_and_blanks_do_not_change_the_numbers(tmp_path):
    rows = invented()
    reference = read_fedra_text(write(tmp_path / "point.txt", lines_of(rows)), detector="d")
    comma = [[v.replace(".", ",") for v in row] for row in rows]
    variants = {
        "comma": write(tmp_path / "comma.txt", lines_of(comma)),
        "crlf": write(tmp_path / "crlf.txt", lines_of(rows), ending="\r\n"),
        "blanks": write(
            tmp_path / "blanks.txt",
            [line.replace(" = ", "\t=  ").replace("Y", " Y") + "  " for line in lines_of(rows)],
        ),
    }
    for name, path in variants.items():
        tracks = read_fedra_text(path, detector="d")
        for field in ("tx", "ty", "x_m", "y_m"):
            assert np.array_equal(getattr(tracks, field), getattr(reference, field)), (name, field)
    assert read_fedra_text(variants["comma"], detector="d").meta["decimal"] == ","
    assert read_fedra_text(variants["crlf"], detector="d").meta["line_ends"] == "CR LF"


def test_header_and_trailer_are_optional(tmp_path):
    rows = invented(5)
    for header, declared, trailer in itertools.product((True, False), (True, False), (True, False)):
        path = write(tmp_path / "t.txt", lines_of(rows, header, declared, trailer))
        tracks = read_fedra_text(path, detector="d")
        assert len(tracks) == 5
        assert tracks.meta["declared_tracks"] == (5 if header and declared else None)
    bare = read_fedra_text(
        write(tmp_path / "bare.txt", lines_of(rows, header=False, trailer=False)), detector="d"
    )
    assert bare.meta["header"] == [] and bare.meta["first_track_line"] == 1 and bare.meta["trailer"] == []


def test_no_tracks(tmp_path):
    empty = read_fedra_text(write(tmp_path / "none.txt", [*HEADER, "N VTs = 0", "(int) 1"]), detector="d")
    assert len(empty) == 0 and empty.meta["declared_tracks"] == 0 and empty.meta["decimal"] is None
    (tmp_path / "zero.txt").write_bytes(b"")
    (tmp_path / "blank.txt").write_bytes(b" \n\n\t\n")
    for name in ("zero.txt", "blank.txt"):
        with pytest.raises(ValueError, match="empty"):
            read_fedra_text(tmp_path / name, detector="d")
    with pytest.raises(ValueError, match="no track lines"):
        read_fedra_text(write(tmp_path / "header.txt", HEADER), detector="d")
    with pytest.raises(ValueError, match="N VTs = 3"):
        read_fedra_text(write(tmp_path / "promise.txt", [*HEADER, "N VTs = 3"]), detector="d")


def test_the_count_of_the_header_is_checked(tmp_path):
    rows = invented(7)
    lines = lines_of(rows)
    lines[3] = "N VTs = 8"
    with pytest.raises(ValueError, match="N VTs = 8.*7 track lines"):
        read_fedra_text(write(tmp_path / "t.txt", lines), detector="d")


GOOD = "VT param: X = 1200.5 Y = 3400 Theta = 0.25 Phi = 1.1"
BROKEN = {
    "a letter in a number": "VT param: X = 12o0.5 Y = 3400 Theta = 0.25 Phi = 1.1",
    "a value missing": "VT param: X = 1200.5 Y = 3400 Theta = 0.25 Phi =",
    "a label missing": "VT param: X = 1200.5 Y = 3400 Theta = 0.25",
    "labels out of order": "VT param: X = 1200.5 Phi = 1.1 Theta = 0.25 Y = 3400",
    "another line": "Warning in <TFile::Init>: no keys recovered",
    "an empty line": "",
    "not a finite number": "VT param: X = 1200.5 Y = nan Theta = 0.25 Phi = 1.1",
    "a decimal comma among points": "VT param: X = 1200,5 Y = 3400 Theta = 0.25 Phi = 1.1",
    "a byte outside ASCII": "VT param: X = 1200.5 Y = 3400 Theta = 0.25 Phi = 1.1µ",
}


@pytest.mark.parametrize("kind", sorted(BROKEN))
def test_a_line_that_cannot_be_read_is_refused_with_its_number(tmp_path, kind):
    lines = lines_of(invented(40))
    number = 4 + 23  # the 23rd track line
    lines[number - 1] = BROKEN[kind]
    path = tmp_path / "t.txt"
    path.write_bytes(("\n".join(lines) + "\n").encode("utf-8"))
    with pytest.raises(ValueError, match=rf"t\.txt.*line {number}\b"):
        read_fedra_text(path, detector="d")


def test_a_good_line_is_good(tmp_path):
    tracks = read_fedra_text(write(tmp_path / "t.txt", [GOOD]), detector="d")
    assert len(tracks) == 1 and tracks.x_m[0] == 1200.5e-6 and tracks.y_m[0] == 3400e-6


def test_a_line_after_the_tracks_is_refused_with_its_number(tmp_path):
    lines = lines_of(invented(10))  # 4 lines of header, 10 track lines, (int) 1
    lines += ["", "end of the macro"]
    with pytest.raises(ValueError, match=r"line 17\b.*after the last track line"):
        read_fedra_text(write(tmp_path / "t.txt", lines), detector="d")
    lone = GOOD.replace(" Y", "\rY")  # a carriage return without a line feed
    with pytest.raises(ValueError, match=r"line 5\b.*carriage return"):
        read_fedra_text(write(tmp_path / "u.txt", [*HEADER, GOOD, lone, GOOD]), detector="d")


def test_what_theta_means(tmp_path):
    rows = invented(30, theta_max=1.5)
    path = write(tmp_path / "t.txt", lines_of(rows))
    _, _, theta, phi = numbers(rows)
    radians = read_fedra_text(path, detector="d", theta="radians")
    assert np.array_equal(radians.tx, np.tan(theta) * np.cos(phi)) and radians.meta["theta"] == "radians"
    rows[9][2] = "1.5708"  # just above π/2
    rows[20][2] = "1.6"
    with pytest.raises(ValueError, match=r"line 14\b.*Theta = 1.5708.*2 track lines"):
        read_fedra_text(write(tmp_path / "u.txt", lines_of(rows)), detector="d", theta="radians")
    assert len(read_fedra_text(tmp_path / "u.txt", detector="d")) == 30  # as tangents they are fine
    rows[3][2] = "-0.1"
    with pytest.raises(ValueError, match=r"line 8\b.*not negative"):
        read_fedra_text(write(tmp_path / "v.txt", lines_of(rows)), detector="d")
    with pytest.raises(ValueError, match="theta"):
        read_fedra_text(path, detector="d", theta="degrees")


def explicit(frame: ScanFrame, x_um, y_um, tx, ty):
    """The transform by its matrix: points (dx, dy, 0) and directions (tx, ty, 1) of the scan."""
    a = frame.matrix
    point = a @ np.array([x_um - frame.origin_um[0], y_um - frame.origin_um[1], np.zeros_like(tx)])
    direction = a @ np.array([tx, ty, np.ones_like(tx)])
    return point[0] * 1e-6, point[1] * 1e-6, direction[0] / direction[2], direction[1] / direction[2]


FRAMES = [
    ScanFrame(origin_um=(5.0e4, 4.0e4), x=x, y=y, z=z)
    for x, y, z in itertools.product(("+x", "-x", "+y", "-y"), ("+x", "-x", "+y", "-y"), ("+z", "-z"))
    if x[1] != y[1]
]


def test_sixteen_scan_frames():
    assert len(FRAMES) == 16 and sum(frame.right_handed for frame in FRAMES) == 8
    rng = np.random.default_rng(3)
    x, y = rng.uniform(0, 1e5, 20), rng.uniform(0, 1e5, 20)
    tx, ty = rng.uniform(-1, 1, 20), rng.uniform(-1, 1, 20)
    for frame in FRAMES:
        a = frame.matrix
        assert sorted(np.abs(a).sum(axis=0)) == [1, 1, 1] and abs(abs(np.linalg.det(a)) - 1) < 1e-15
        assert frame.right_handed == (np.linalg.det(a) > 0)
        got = frame.apply(x, y, tx, ty)
        for value, expected in zip(got, explicit(frame, x, y, tx, ty), strict=True):
            assert np.allclose(value, expected, rtol=0, atol=1e-15)


def test_scan_frames_by_hand():
    same = ScanFrame(origin_um=(0.0, 0.0))
    assert same.apply(1e4, 2e4, 0.3, -0.1) == pytest.approx((0.01, 0.02, 0.3, -0.1), abs=1e-15)
    turned = ScanFrame(origin_um=(0.0, 0.0), x="+y", y="-x")  # turned by 90° about the normal
    assert turned.apply(1e4, 0.0, 0.3, 0.0) == pytest.approx((0.0, 0.01, 0.0, 0.3), abs=1e-15)
    assert turned.right_handed
    over = ScanFrame(origin_um=(1e4, 1e4), y="-y", z="-z")  # the film scanned from its other side
    assert over.apply(2e4, 3e4, 0.3, 0.2) == pytest.approx((0.01, -0.02, -0.3, 0.2), abs=1e-15)
    assert over.right_handed
    mirror = ScanFrame(origin_um=(0.0, 0.0), y="-y")
    assert not mirror.right_handed and mirror.apply(0.0, 1e4, 0.0, 0.3) == pytest.approx((0, -0.01, 0, -0.3))
    for wrong in (
        {"x": "+x", "y": "-x"},
        {"z": "z"},
        {"x": "x"},
        {"origin_um": (1.0,)},
        {"origin_um": (0, np.nan)},
    ):
        with pytest.raises(ValueError):
            ScanFrame(**({"origin_um": (0.0, 0.0)} | wrong))


def test_a_file_read_into_the_detector_frame(tmp_path):
    rows = invented(20)
    path = write(tmp_path / "t.txt", lines_of(rows))
    frame = ScanFrame(origin_um=(5.0e4, 4.0e4), x="-y", y="+x", z="-z")
    tracks = read_fedra_text(path, detector="d", scan_frame=frame, area_cm2=90.0, cuts={"films": "unknown"})
    x, y, theta, phi = numbers(rows)
    expected = explicit(frame, x, y, theta * np.cos(phi), theta * np.sin(phi))
    for value, wanted in zip((tracks.x_m, tracks.y_m, tracks.tx, tracks.ty), expected, strict=True):
        assert np.allclose(value, wanted, rtol=0, atol=1e-15)
    assert tracks.frame == "detector" and tracks.area_cm2 == 90.0 and tracks.cuts == {"films": "unknown"}
    assert tracks.meta["scan_frame"] == {
        "x": "-y",
        "y": "+x",
        "z": "-z",
        "origin_um": [5.0e4, 4.0e4],
        "right_handed": frame.right_handed,
    }


def test_the_exposure(tmp_path):
    path = write(tmp_path / "t.txt", lines_of(invented(3)))
    tracks = read_fedra_text(path, detector="d", exposure_start="2022-03-01", exposure_end="2022-03-31")
    assert tracks.exposure_days == 30.0 and tracks.meta["exposure"] == {
        "start": "2022-03-01",
        "end": "2022-03-31",
    }
    timed = read_fedra_text(
        path, detector="d", exposure_start="2022-03-01T12:00", exposure_end="2022-03-02T00:00"
    )
    assert timed.exposure_days == 0.5
    assert read_fedra_text(path, detector="d", exposure_days=30).exposure_days == 30.0
    for wrong in (
        {"exposure_start": "2022-03-01"},
        {"exposure_start": "2022-03-31", "exposure_end": "2022-03-01"},
        {"exposure_start": "2022-03-01", "exposure_end": "2022-03-31", "exposure_days": 31},
        {"exposure_start": "1 March 2022", "exposure_end": "2022-03-31"},
        {"exposure_days": -1},
        {"source": "measured"},
    ):
        with pytest.raises(ValueError):
            read_fedra_text(path, detector="d", **wrong)


def test_the_example_of_formats_md(tmp_path, monkeypatch):
    """The Python block of `docs/formats.md` runs as written and gives the numbers it states."""
    blocks = re.findall(r"```python\n(.*?)```", FORMATS_MD.read_text(encoding="utf-8"), flags=re.S)
    assert len(blocks) == 1
    monkeypatch.chdir(tmp_path)
    space = {}
    exec(compile(blocks[0], str(FORMATS_MD), "exec"), space)
    assert space["binned"].binned == 2 and space["binned"].outside_circle == 1
    assert space["tracks"].frame == "scan" and space["in_detector"].frame == "detector"


def test_the_example_of_notes_tracks_md(tmp_path, monkeypatch):
    """The Python block of `docs/notes/tracks.md` runs as written on a file of this format; it draws."""
    blocks = re.findall(r"```python\n(.*?)```", TRACKS_MD.read_text(encoding="utf-8"), flags=re.S)
    assert len(blocks) == 1
    monkeypatch.chdir(tmp_path)
    write(tmp_path / "tracks.txt", lines_of(invented(200)))
    space = {}
    exec(compile(blocks[0], str(TRACKS_MD), "exec"), space)
    assert space["binned"].n_tracks == 200 and (tmp_path / "tracks.png").stat().st_size > 0
