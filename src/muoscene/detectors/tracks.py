"""detectors.tracks — `TrackSet`, the common container of measured and synthetic tracks.

A track is a line through a film: its slopes (tx, ty) = (dx/dz, dy/dz) and, where the scan gives
it, the point where it crosses the plane of the film. `TrackSet` holds the tracks of one detector
with what a comparison with a model needs besides (roadmap 0.1-15): the frame they are given in,
the area scanned, the length of the exposure, the cuts already applied and where they come from.
It is defined here, in layer 2, because `synth` and `data` (layer 4) produce and consume it; the
readers in `muoscene.io.tracks` only fill it.

Decided on 28 November 2026 (`docs/notes/tracks.md`):

- Two frames. "detector": the axes x_d, y_d, z_d of the detector (`docs/conventions.md`, "Detector
  pose"), positions from its reference point; only such a set is compared with a model. "scan":
  the axes and the origin of the scan, as in the file, while the scan-frame transform of the
  exposure is not known. A set does not change its frame: the transform is applied by the reader
  (`muoscene.io.tracks.ScanFrame`), and a set in the scan frame is read again once it is known.
- Fields are what the scanning outputs at hand give. The text output of the FEDRA macro of the
  first exposures has no number of segments, first film or quality of a track; such columns come
  with the reader that brings them, and with them version 2 of the file format.
- A set is immutable: its arrays are read-only, and `cuts` and `meta` are given out as fresh copies.
- Files: `save` writes `.npz` with a header in JSON that names the format and its version; `load`
  refuses another version and anything else in the file. No pickled objects are written or read.

Units: slopes are dimensionless (tangents), positions in metres, areas in cm², times in days.

Layer 2: may import only `muoscene` modules of layer <= 2 (see `muoscene._layers`).
"""

from __future__ import annotations

import json
import math
import pathlib

import numpy as np

FORMAT = "muoscene.TrackSet"
"""Name of the file format written by `TrackSet.save`."""
FORMAT_VERSION = 1
"""Version of that format; `TrackSet.load` reads only this one."""
FRAMES = ("detector", "scan")
"""Frames a set can be given in: that of the detector, or that of the scan, as in the file."""
SOURCES = ("data", "synthetic")
"""Where the tracks of a set come from."""
FIELDS = (
    "tx",
    "ty",
    "x_m",
    "y_m",
    "frame",
    "source",
    "detector",
    "area_cm2",
    "exposure_days",
    "_cuts",
    "_meta",
)


def _slopes(values, name: str) -> np.ndarray:
    array = np.array(values, dtype=float)
    if array.ndim != 1:
        raise ValueError(f"{name} must be one-dimensional: one value per track")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must hold finite numbers")
    array.flags.writeable = False
    return array


def _positive(value, name: str) -> float | None:
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, int | float) or not (0 < value < math.inf):
        raise ValueError(f"{name} must be a positive finite number or None")
    return float(value)


def _plain(value, name: str) -> str:
    """The JSON text of a dictionary of plain values: str keys, numbers, strings, booleans, None, lists."""
    if value is None:
        value = {}
    if not isinstance(value, dict) or not all(isinstance(key, str) for key in value):
        raise ValueError(f"{name} must be a dictionary with string keys")
    try:
        return json.dumps(value, allow_nan=False, sort_keys=True, ensure_ascii=False)
    except (TypeError, ValueError) as error:
        raise ValueError(
            f"{name} must hold plain values that go to JSON, finite numbers only: {error}"
        ) from None


class TrackSet:
    """Tracks of one detector: their slopes and, where the scan gives them, their positions.

    tx, ty          slopes dx/dz, dy/dz of the tracks in `frame`; shape (n,), finite
    frame           "detector" or "scan" (see the module docstring)
    source          "data" or "synthetic"
    detector        name of the detector, a non-empty string
    x_m, y_m        where each track crosses the plane of the film, m, in `frame`; shape (n,);
                    both None when the scan does not give it
    area_cm2        area scanned, cm²; None when not known
    exposure_days   length of the exposure, days; None when not known
    cuts            cuts already applied to the tracks, name → plain value; {} when none is known
    meta            where the set comes from: file, its hash, the reader and its arguments, the
                    scan-frame transform; a dictionary of plain values that go to JSON

    A slope pair is a line, so a track has no side: a film records a muon whichever side it
    enters (`docs/conventions.md`, "Detector pose"). An empty set (n = 0) is allowed.
    """

    __slots__ = FIELDS

    def __init__(
        self,
        tx,
        ty,
        *,
        frame: str,
        source: str,
        detector: str,
        x_m=None,
        y_m=None,
        area_cm2: float | None = None,
        exposure_days: float | None = None,
        cuts: dict | None = None,
        meta: dict | None = None,
    ) -> None:
        tx, ty = _slopes(tx, "tx"), _slopes(ty, "ty")
        if tx.shape != ty.shape:
            raise ValueError(f"tx and ty must have the same length: {tx.size} and {ty.size}")
        if (x_m is None) != (y_m is None):
            raise ValueError("x_m and y_m are given together or not at all")
        if x_m is not None:
            x_m, y_m = _slopes(x_m, "x_m"), _slopes(y_m, "y_m")
            if x_m.shape != tx.shape or y_m.shape != tx.shape:
                raise ValueError(f"x_m and y_m must have one value per track: {tx.size}")
        if frame not in FRAMES:
            raise ValueError(f"frame must be one of {FRAMES}, not {frame!r}")
        if source not in SOURCES:
            raise ValueError(f"source must be one of {SOURCES}, not {source!r}")
        if not isinstance(detector, str) or not detector.strip():
            raise ValueError("detector must be the name of the detector, a non-empty string")
        values = {
            "tx": tx,
            "ty": ty,
            "x_m": x_m,
            "y_m": y_m,
            "frame": frame,
            "source": source,
            "detector": detector,
            "area_cm2": _positive(area_cm2, "area_cm2"),
            "exposure_days": _positive(exposure_days, "exposure_days"),
            "_cuts": _plain(cuts, "cuts"),
            "_meta": _plain(meta, "meta"),
        }
        for name, value in values.items():
            object.__setattr__(self, name, value)

    def __setattr__(self, name, value) -> None:
        raise AttributeError("a TrackSet is immutable: make a new one")

    def __delattr__(self, name) -> None:
        raise AttributeError("a TrackSet is immutable: make a new one")

    def __len__(self) -> int:
        return self.tx.size

    def __repr__(self) -> str:
        positions = "with positions" if self.x_m is not None else "without positions"
        return (
            f"TrackSet({len(self)} tracks of {self.detector!r}, {self.source}, {self.frame} frame,"
            f" {positions}, area {self.area_cm2} cm², exposure {self.exposure_days} days)"
        )

    @property
    def cuts(self) -> dict:
        """Cuts already applied to the tracks: a fresh copy, name → plain value."""
        return json.loads(self._cuts)

    @property
    def meta(self) -> dict:
        """Where the set comes from: a fresh copy of a dictionary of plain values."""
        return json.loads(self._meta)

    def _header(self) -> dict:
        return {
            "format": FORMAT,
            "version": FORMAT_VERSION,
            "n_tracks": len(self),
            "frame": self.frame,
            "source": self.source,
            "detector": self.detector,
            "positions": self.x_m is not None,
            "area_cm2": self.area_cm2,
            "exposure_days": self.exposure_days,
            "cuts": self.cuts,
            "meta": self.meta,
        }

    def save(self, path) -> pathlib.Path:
        """Write the set to `path`, a file ending in `.npz`; returns the path.

        The file holds the arrays `tx`, `ty` (and `x_m`, `y_m` when there are positions) and
        `header`, the JSON text of the format, its version and the other fields.
        """
        path = pathlib.Path(path)
        if path.suffix != ".npz":
            raise ValueError(f"a TrackSet is written to a file ending in .npz, not {path.name!r}")
        arrays = {"tx": self.tx, "ty": self.ty}
        if self.x_m is not None:
            arrays |= {"x_m": self.x_m, "y_m": self.y_m}
        header = np.array(json.dumps(self._header(), allow_nan=False, ensure_ascii=False))
        with path.open("wb") as handle:
            np.savez(handle, header=header, **arrays)
        return path

    @classmethod
    def load(cls, path) -> TrackSet:
        """Read a set written by `save`; another format, another version or other arrays are refused."""
        path = pathlib.Path(path)
        try:
            with np.load(path, allow_pickle=False) as data:
                arrays = {name: data[name] for name in data.files}
        except (OSError, ValueError) as error:
            raise ValueError(f"{path.name}: not a file written by TrackSet.save: {error}") from None
        if "header" not in arrays or arrays["header"].shape != () or arrays["header"].dtype.kind != "U":
            raise ValueError(f"{path.name}: not a file written by TrackSet.save: no header")
        try:
            header = json.loads(str(arrays.pop("header")[()]))
        except json.JSONDecodeError as error:
            raise ValueError(f"{path.name}: the header is not JSON: {error}") from None
        if not isinstance(header, dict) or header.get("format") != FORMAT:
            raise ValueError(f"{path.name}: not a file of the format {FORMAT}")
        if header.get("version") != FORMAT_VERSION:
            raise ValueError(
                f"{path.name}: format version {header.get('version')!r};"
                f" this muoscene reads version {FORMAT_VERSION}"
            )
        expected = {"tx", "ty"} | ({"x_m", "y_m"} if header.get("positions") else set())
        if set(arrays) != expected:
            raise ValueError(f"{path.name}: arrays {sorted(arrays)}, expected {sorted(expected)}")
        tracks = cls(
            arrays["tx"],
            arrays["ty"],
            x_m=arrays.get("x_m"),
            y_m=arrays.get("y_m"),
            frame=header.get("frame"),
            source=header.get("source"),
            detector=header.get("detector"),
            area_cm2=header.get("area_cm2"),
            exposure_days=header.get("exposure_days"),
            cuts=header.get("cuts"),
            meta=header.get("meta"),
        )
        if len(tracks) != header.get("n_tracks"):
            raise ValueError(f"{path.name}: {len(tracks)} tracks, the header says {header.get('n_tracks')!r}")
        return tracks
