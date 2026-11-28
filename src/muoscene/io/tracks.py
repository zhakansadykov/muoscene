"""io.tracks — readers of scanning output into `muoscene.detectors.tracks.TrackSet`.

The first reader (roadmap 0.1-15, decided on 27 November 2026) is for the form in which the first
exposures of the scanning group came: the text output of the FEDRA macro `evaluate_VT_nseg.C`
(`read_fedra_text`). Its format, column by column, and what of it the group has confirmed are in
`docs/formats.md`; the decisions are in `docs/notes/tracks.md`. In short:

- Lines before the first track line are the header and are kept in `meta`; a line
  `N VTs = <n>` among them is checked against the number of tracks read. After the last track
  line only empty lines and the line `(int) <n>` that ROOT prints when a macro returns may follow.
- A track line is `VT param: X = <x> Y = <y> Theta = <θ> Phi = <φ>`: the point where the track
  crosses the plane of the scan, in micrometres, and its direction. Theta is the tangent of the
  angle from the normal of the film by default — in FEDRA `EdbSegP::Theta()` is sqrt(TX² + TY²)
  — or that angle in radians (`theta="radians"`); Phi is the angle from the x axis of the scan
  towards its y axis, in radians, as FEDRA's `EdbSegP::Phi()` = atan2(TY, TX). So
  tx = tan θ cos φ and ty = tan θ sin φ in the axes of the scan.
- Numbers have a decimal point or a decimal comma; a file with both is refused. Runs of blanks
  are one separator; line ends LF or CR LF.
- A line that cannot be read is refused with its number, never skipped; so is an empty file,
  a header that promises tracks that are not there, and a Theta outside its meaning (negative,
  or at least π/2 as radians).
- The file is read without a loop over its lines in Python: the labels are replaced by numbers
  and `numpy.loadtxt` reads eight columns, four of which must be the labels in their order. Only
  when that fails is the file walked line by line, to name the line.

The scan-frame transform (`ScanFrame`) takes the axes of the scan to those of the detector
(`docs/conventions.md`, "Detector pose", "Scan frame"). Without it a set stays in the frame of
the scan (`TrackSet.frame` "scan") and cannot be compared with a model. With it, positions are
turned and moved to the reference point of the detector, and slopes are turned and, where the z
axis of the scan runs against the normal of the detector, change sign: a point (X, Y) on the film
goes where the axes go, a direction (tx, ty, 1) goes to (a, b, ±1) and its slopes are (±a, ±b).

Units: positions in the file in micrometres, in a `TrackSet` in metres; angles in radians.

Layer 5: may import only `muoscene` modules of layer <= 5 (see `muoscene._layers`).
"""

from __future__ import annotations

import datetime
import hashlib
import io
import math
import pathlib
import re
from dataclasses import dataclass

import numpy as np

from .. import __version__
from ..detectors.tracks import SOURCES, TrackSet

MICROMETRE_M = 1.0e-6
"""One micrometre in metres: positions of FEDRA are in micrometres."""
THETA_MEANINGS = ("tangent", "radians")
"""What the column Theta can mean: the tangent of the angle from the normal, or the angle in radians."""
TRACK_LINE = re.compile(
    rb"VT param: X = (\S+) Y = (\S+) Theta = (\S+) Phi = (\S+)",
)
"""A track line of `evaluate_VT_nseg.C`, after runs of blanks are made one space."""
LABELS = ((b"VT param: X = ", b"-1 "), (b" Y = ", b" -2 "), (b" Theta = ", b" -3 "), (b" Phi = ", b" -4 "))
"""The labels of a track line and the numbers that stand for them while it is read."""
DECLARED = re.compile(rb"^\s*N VTs\s*=\s*(\d+)\s*$")
"""The line of the header that gives the number of tracks."""
RETURNED = re.compile(rb"^\(int\) -?\d+$")
"""The line ROOT prints when the macro returns, after the tracks."""
AXES = ("+x", "-x", "+y", "-y")
NORMALS = ("+z", "-z")


@dataclass(frozen=True, kw_only=True)
class ScanFrame:
    """How the axes of a scan lie in the frame of the detector: the scan-frame transform of an exposure.

    origin_um   the scan coordinates (X, Y) of the reference point of the detector, µm
    x, y        the detector axis the x and the y axis of the scan run along: "+x", "-x", "+y" or
                "-y"; the two must lie on different axes
    z           "+z" when the z axis of the scan runs along the normal z_d of the detector, "-z"
                when against it (the film was scanned from the other side)

    Every one of the 16 combinations is allowed, also those that make the axes of the scan a
    left-handed triple in the detector frame (`right_handed` False): the transform is what the
    exposure log says, not what it ought to be. It is never guessed: without it a set stays in
    the frame of the scan.
    """

    origin_um: tuple[float, float]
    x: str = "+x"
    y: str = "+y"
    z: str = "+z"

    def __post_init__(self) -> None:
        origin = np.asarray(self.origin_um, dtype=float)
        if origin.shape != (2,) or not np.all(np.isfinite(origin)):
            raise ValueError("origin_um must be two finite scan coordinates X, Y of the reference point, µm")
        object.__setattr__(self, "origin_um", tuple(origin.tolist()))
        if self.x not in AXES or self.y not in AXES:
            raise ValueError(f"x and y must be one of {AXES}")
        if self.x[1] == self.y[1]:
            raise ValueError(f"the x and the y axis of the scan cannot both run along {self.x[1]}")
        if self.z not in NORMALS:
            raise ValueError(f"z must be one of {NORMALS}")

    @property
    def matrix(self) -> np.ndarray:
        """3 × 3 matrix whose columns are the axes of the scan in the detector frame; read-only."""
        unit = {"x": (1.0, 0.0, 0.0), "y": (0.0, 1.0, 0.0), "z": (0.0, 0.0, 1.0)}
        columns = [
            (1.0 if axis[0] == "+" else -1.0) * np.array(unit[axis[1]]) for axis in (self.x, self.y, self.z)
        ]
        out = np.array(columns).T
        out.flags.writeable = False
        return out

    @property
    def right_handed(self) -> bool:
        """Whether the axes of the scan are a right-handed triple in the detector frame."""
        return bool(np.linalg.det(self.matrix) > 0)

    def to_dict(self) -> dict:
        """The transform as plain values, for `TrackSet.meta`."""
        return {
            "x": self.x,
            "y": self.y,
            "z": self.z,
            "origin_um": list(self.origin_um),
            "right_handed": self.right_handed,
        }

    def apply(self, x_um, y_um, tx, ty) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Positions (µm, from the scan origin) and slopes in the scan → x_m, y_m, tx, ty of the detector.

        Positions go where the axes go and are moved to the reference point; slopes are turned and
        divided by the z component of the direction (tx, ty, 1) in the detector frame, ±1.
        """
        a = self.matrix
        dx = np.asarray(x_um, dtype=float) - self.origin_um[0]
        dy = np.asarray(y_um, dtype=float) - self.origin_um[1]
        x_m = (a[0, 0] * dx + a[0, 1] * dy) * MICROMETRE_M
        y_m = (a[1, 0] * dx + a[1, 1] * dy) * MICROMETRE_M
        tx, ty = np.asarray(tx, dtype=float), np.asarray(ty, dtype=float)
        sign = a[2, 2]  # the z component of the direction (tx, ty, 1) of the scan in the detector frame
        tx_d = (a[0, 0] * tx + a[0, 1] * ty) / sign
        ty_d = (a[1, 0] * tx + a[1, 1] * ty) / sign
        return x_m, y_m, tx_d, ty_d


def _days(start, end, exposure_days):
    """Length of the exposure from its dates, checked against `exposure_days`; and the dates as text."""
    if exposure_days is not None and (
        isinstance(exposure_days, bool)
        or not isinstance(exposure_days, int | float)
        or not 0 < exposure_days < math.inf
    ):
        raise ValueError("exposure_days must be a positive finite number or None")
    if start is None and end is None:
        return exposure_days, None
    if start is None or end is None:
        raise ValueError("exposure_start and exposure_end are given together or not at all")

    def moment(value, name):
        if isinstance(value, datetime.datetime | datetime.date):
            value = value.isoformat()
        if not isinstance(value, str):
            raise ValueError(f"{name} must be a date or a time in ISO 8601, such as 2022-03-01")
        try:
            return value, datetime.datetime.fromisoformat(value)
        except ValueError:
            raise ValueError(f"{name} must be a date or a time in ISO 8601, not {value!r}") from None

    (start_text, t0), (end_text, t1) = moment(start, "exposure_start"), moment(end, "exposure_end")
    if (t0.tzinfo is None) != (t1.tzinfo is None):
        raise ValueError("exposure_start and exposure_end must both have a time zone or both have none")
    days = (t1 - t0).total_seconds() / 86400.0
    if not days > 0:
        raise ValueError("exposure_end must be after exposure_start")
    if exposure_days is not None and not math.isclose(days, exposure_days, rel_tol=1e-9):
        raise ValueError(f"exposure_days {exposure_days} does not match the dates, {days} days")
    return days, {"start": start_text, "end": end_text}


def _line_of(text: bytes, offset: int) -> int:
    """Number, from 1, of the line of `text` that holds the byte `offset`."""
    return text.count(b"\n", 0, offset) + 1


def _show(line: bytes) -> str:
    return repr(line[:100].decode("utf-8", errors="replace"))


def _bad_line(body: bytes, first: int, name: str) -> ValueError:
    """The error naming the first line of `body` that is not a track line; `first`: number of line one."""
    for number, line in enumerate(body.split(b"\n"), start=first):
        match = TRACK_LINE.fullmatch(line)
        if match is None:
            return ValueError(f"{name}, line {number}: not a track line of evaluate_VT_nseg.C: {_show(line)}")
        for value in match.groups():
            try:
                float(value.replace(b",", b"."))
            except ValueError:
                return ValueError(
                    f"{name}, line {number}: {value.decode(errors='replace')!r} is not a number"
                )
    return ValueError(f"{name}: the track lines could not be read")


def read_fedra_text(
    path,
    *,
    detector: str,
    theta: str = "tangent",
    scan_frame: ScanFrame | None = None,
    source: str = "data",
    area_cm2: float | None = None,
    exposure_days: float | None = None,
    exposure_start=None,
    exposure_end=None,
    cuts: dict | None = None,
) -> TrackSet:
    """Read the text output of the FEDRA macro `evaluate_VT_nseg.C` into a `TrackSet`.

    path            the file
    detector        name of the detector the tracks belong to
    theta           what the column Theta is: "tangent" (default) of the angle from the normal of
                    the film, or that angle in "radians"
    scan_frame      the scan-frame transform of the exposure (`ScanFrame`); None leaves the set in
                    the frame of the scan, with positions from the origin of the scan
    source          "data" or "synthetic"
    area_cm2        area scanned, cm², when known
    exposure_days   length of the exposure, days, when known
    exposure_start, exposure_end
                    dates (or times) of the exposure in ISO 8601, when known: they give the length,
                    which must then agree with `exposure_days`
    cuts            cuts the group applied before writing the file, as far as they are known

    `meta` records the name, the size and the SHA-256 of the file, the header and the lines after
    the tracks, the decimal mark, the meaning given to Theta and Phi, the transform and the version
    of muoscene. Every refusal is a ValueError that names the file and, where there is one, the line.
    """
    path = pathlib.Path(path)
    name = path.name
    if theta not in THETA_MEANINGS:
        raise ValueError(f"theta must be one of {THETA_MEANINGS}, not {theta!r}")
    if scan_frame is not None and not isinstance(scan_frame, ScanFrame):
        raise TypeError("scan_frame must be a ScanFrame or None")
    if source not in SOURCES:
        raise ValueError(f"source must be one of {SOURCES}, not {source!r}")
    days, dates = _days(exposure_start, exposure_end, exposure_days)
    raw = path.read_bytes()
    if not raw.strip():
        raise ValueError(f"{name}: the file is empty")
    text = raw.replace(b"\r\n", b"\n")
    lone = text.find(b"\r")
    if lone >= 0:
        raise ValueError(f"{name}, line {_line_of(text, lone)}: a carriage return without a line feed")

    first = re.search(rb"(?m)^VT param:", text)
    if first is None:
        header_text, body, tail = text, b"", b""
    else:
        start = first.start()
        last = text.rfind(b"\nVT param:", start)
        last = start if last < 0 else last + 1
        end = text.find(b"\n", last)
        end = len(text) if end < 0 else end
        header_text, body, tail = text[:start], text[start:end], text[end + 1 :]
    header = header_text.split(b"\n")
    if header and header[-1] == b"":
        header.pop()
    first_line = len(header) + 1
    declared = None
    for line in header:
        found = DECLARED.match(line)
        if found:
            declared = int(found.group(1))
    trailer = tail.split(b"\n")
    after = first_line + body.count(b"\n") + 1
    for number, line in enumerate(trailer, start=after):
        if line.strip() and not RETURNED.match(line.strip()):
            raise ValueError(f"{name}, line {number}: after the last track line: {_show(line)}")

    if body:
        if b"\t" in body or b"  " in body or b"\n " in body or b" \n" in body or body.endswith(b" "):
            body = re.sub(rb"[ \t]+", b" ", body)
            body = re.sub(rb" ?\n ?", b"\n", body).strip(b" ")
        blank = body.find(b"\n\n")
        if blank >= 0:
            raise ValueError(
                f"{name}, line {first_line + _line_of(body, blank)}: an empty line among the track lines"
            )
        comma, point = body.find(b","), body.find(b".")
        if comma >= 0 and point >= 0:
            raise ValueError(
                f"{name}: decimal points (line {first_line - 1 + _line_of(body, point)}) and decimal commas"
                f" (line {first_line - 1 + _line_of(body, comma)}) in one file"
            )
        decimal = "," if comma >= 0 else "."
        numbers = body.replace(b",", b".") if comma >= 0 else body
        for label, mark in LABELS:
            numbers = numbers.replace(label, mark)
        try:
            values = np.loadtxt(io.StringIO(numbers.decode("ascii")), dtype=float, comments=None, ndmin=2)
        except (ValueError, UnicodeDecodeError):
            raise _bad_line(body, first_line, name) from None
        if values.shape[1] != 8 or values.shape[0] != body.count(b"\n") + 1:
            raise _bad_line(body, first_line, name)
        wrong = np.flatnonzero(np.any(values[:, 0::2] != [-1.0, -2.0, -3.0, -4.0], axis=1))
        if wrong.size:
            line = body.split(b"\n")[wrong[0]]
            raise ValueError(f"{name}, line {first_line + wrong[0]}: labels out of order: {_show(line)}")
        bad = np.flatnonzero(~np.all(np.isfinite(values[:, 1::2]), axis=1))
        if bad.size:
            raise ValueError(f"{name}, line {first_line + bad[0]}: a value is not a finite number")
        x_um, y_um, angle, phi = values[:, 1], values[:, 3], values[:, 5], values[:, 7]
    else:
        decimal = None
        x_um = y_um = angle = phi = np.zeros(0)

    if theta == "tangent":
        bad = np.flatnonzero(angle < 0.0)
        meaning = "a tangent is not negative"
    else:
        bad = np.flatnonzero((angle < 0.0) | (angle >= 0.5 * math.pi))
        meaning = "an angle from the normal in radians lies in [0, π/2)"
    if bad.size:
        raise ValueError(
            f"{name}, line {first_line + bad[0]}: Theta = {angle[bad[0]]:g}, and {meaning};"
            f" {bad.size} track lines like it"
        )
    n = x_um.size
    if n == 0 and declared is None:
        raise ValueError(f"{name}: no track lines, and no line 'N VTs = 0' in the header")
    if declared is not None and declared != n:
        raise ValueError(f"{name}: the header says N VTs = {declared}, the file has {n} track lines")

    tangent = angle if theta == "tangent" else np.tan(angle)
    tx, ty = tangent * np.cos(phi), tangent * np.sin(phi)
    if scan_frame is None:
        frame, x_m, y_m = "scan", x_um * MICROMETRE_M, y_um * MICROMETRE_M
    else:
        frame = "detector"
        x_m, y_m, tx, ty = scan_frame.apply(x_um, y_um, tx, ty)

    meta = {
        "reader": "muoscene.io.tracks.read_fedra_text",
        "muoscene_version": __version__,
        "file": {"name": name, "bytes": len(raw), "sha256": hashlib.sha256(raw).hexdigest()},
        "format": "text output of the FEDRA macro evaluate_VT_nseg.C",
        "header": [line.decode("utf-8", errors="replace") for line in header],
        "declared_tracks": declared,
        "trailer": [line.strip().decode("utf-8", errors="replace") for line in trailer if line.strip()],
        "first_track_line": first_line if n else None,
        "decimal": decimal,
        "line_ends": "CR LF" if b"\r\n" in raw else "LF",
        "theta": theta,
        "phi": "radians from the x axis of the scan towards its y axis",
        "positions": "micrometres from the origin of the scan in the file",
        "scan_frame": None if scan_frame is None else scan_frame.to_dict(),
        "exposure": dates,
    }
    return TrackSet(
        tx,
        ty,
        x_m=x_m,
        y_m=y_m,
        frame=frame,
        source=source,
        detector=detector,
        area_cm2=area_cm2,
        exposure_days=days,
        cuts=cuts,
        meta=meta,
    )
