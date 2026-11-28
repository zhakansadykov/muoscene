"""`muoscene.detectors.tracks.TrackSet`: its input is checked, nothing changes it, a file gives it back.

Roadmap 0.1-15. Sets are made by hand with invented numbers; files are written to a temporary
folder.
"""

import json

import numpy as np
import pytest

from muoscene.detectors import TrackSet
from muoscene.detectors.tracks import FORMAT, FORMAT_VERSION


def three(**changes):
    values = {
        "tx": [0.1, -0.2, 0.0],
        "ty": [0.0, 0.3, -0.4],
        "x_m": [0.01, 0.02, 0.03],
        "y_m": [-0.01, 0.0, 0.01],
        "frame": "detector",
        "source": "synthetic",
        "detector": "film A",
        "area_cm2": 120.0,
        "exposure_days": 30.0,
        "cuts": {"min_films": 4},
        "meta": {"seed": 7, "file": {"name": "invented.txt"}, "list": [1, 2.5, None, "x"]},
    }
    values |= changes
    tx, ty = values.pop("tx"), values.pop("ty")
    return TrackSet(tx, ty, **values)


def same(a: TrackSet, b: TrackSet) -> bool:
    arrays = all(
        (getattr(a, k) is None and getattr(b, k) is None) or np.array_equal(getattr(a, k), getattr(b, k))
        for k in ("tx", "ty", "x_m", "y_m")
    )
    fields = ("frame", "source", "detector", "area_cm2", "exposure_days", "cuts", "meta")
    return arrays and all(getattr(a, k) == getattr(b, k) for k in fields)


def test_a_set_holds_what_it_was_given():
    tracks = three()
    assert len(tracks) == 3
    assert tracks.tx.tolist() == [0.1, -0.2, 0.0] and tracks.ty.tolist() == [0.0, 0.3, -0.4]
    assert tracks.frame == "detector" and tracks.source == "synthetic" and tracks.detector == "film A"
    assert tracks.area_cm2 == 120.0 and tracks.exposure_days == 30.0
    assert tracks.cuts == {"min_films": 4} and tracks.meta["list"] == [1, 2.5, None, "x"]
    assert "3 tracks of 'film A'" in repr(tracks)


def test_nothing_changes_a_set():
    given = np.array([0.1, 0.2])
    tracks = TrackSet(given, [0.0, 0.0], frame="scan", source="data", detector="d")
    given[0] = 5.0  # the set holds a copy
    assert tracks.tx[0] == 0.1
    with pytest.raises(ValueError):
        tracks.tx[0] = 1.0
    with pytest.raises(AttributeError):
        tracks.frame = "detector"
    with pytest.raises(AttributeError):
        del tracks.tx
    meta = three().meta
    meta["seed"] = 8
    cuts = three().cuts
    cuts["min_films"] = 1
    assert three().meta["seed"] == 7 and three().cuts == {"min_films": 4}


def test_an_empty_set_and_one_without_positions():
    empty = TrackSet([], [], frame="detector", source="data", detector="d")
    assert len(empty) == 0 and empty.x_m is None and empty.cuts == {} and empty.meta == {}
    bare = three(x_m=None, y_m=None, area_cm2=None, exposure_days=None, cuts=None, meta=None)
    assert bare.x_m is None and bare.y_m is None and bare.area_cm2 is None and bare.exposure_days is None


@pytest.mark.parametrize(
    "changes",
    [
        {"tx": [[0.1, 0.2, 0.3]]},  # not one-dimensional
        {"tx": [0.1, 0.2]},  # lengths differ
        {"ty": [0.0, np.nan, 0.0]},
        {"tx": [0.0, np.inf, 0.0]},
        {"x_m": None},  # one of the two positions
        {"y_m": [0.0, 0.0]},
        {"x_m": [0.0, np.nan, 0.0]},
        {"frame": "film"},
        {"source": "measured"},
        {"detector": ""},
        {"detector": "  "},
        {"detector": 3},
        {"area_cm2": 0.0},
        {"area_cm2": -1.0},
        {"area_cm2": np.inf},
        {"area_cm2": True},
        {"exposure_days": 0},
        {"exposure_days": "30"},
        {"cuts": {"min_films": np.nan}},
        {"cuts": {1: "a"}},
        {"meta": {"array": np.zeros(2)}},
        {"meta": ["not", "a", "dictionary"]},
    ],
)
def test_wrong_input_is_refused(changes):
    with pytest.raises(ValueError):
        three(**changes)


def test_a_file_gives_the_set_back(tmp_path):
    for tracks in (
        three(),
        three(x_m=None, y_m=None, area_cm2=None, cuts=None),
        TrackSet([], [], frame="scan", source="data", detector="d"),
    ):
        path = tracks.save(tmp_path / "tracks.npz")
        back = TrackSet.load(path)
        assert same(tracks, back), (tracks, back)
        assert back.tx.flags.writeable is False


def write(path, header, **arrays):
    with path.open("wb") as handle:
        np.savez(handle, header=np.array(json.dumps(header)), **arrays)


def test_other_files_are_refused(tmp_path):
    good = tmp_path / "good.npz"
    three().save(good)
    with np.load(good) as data:
        header, arrays = (
            json.loads(str(data["header"][()])),
            {k: data[k] for k in data.files if k != "header"},
        )
    cases = {
        "version": header | {"version": FORMAT_VERSION + 1},
        "format": header | {"format": FORMAT + "x"},
        "count": header | {"n_tracks": 4},
    }
    for name, changed in cases.items():
        path = tmp_path / f"{name}.npz"
        write(path, changed, **arrays)
        with pytest.raises(ValueError):
            TrackSet.load(path)
    extra = tmp_path / "extra.npz"
    write(extra, header, **arrays, quality=np.zeros(3))
    with pytest.raises(ValueError, match="arrays"):
        TrackSet.load(extra)
    nothing = tmp_path / "nothing.npz"
    with nothing.open("wb") as handle:
        np.savez(handle, tx=np.zeros(1))
    with pytest.raises(ValueError, match="no header"):
        TrackSet.load(nothing)
    text = tmp_path / "text.npz"
    text.write_text("not a zip file")
    with pytest.raises(ValueError):
        TrackSet.load(text)
    with pytest.raises(ValueError, match=".npz"):
        three().save(tmp_path / "tracks.dat")
