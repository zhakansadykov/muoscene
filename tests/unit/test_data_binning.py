"""`muoscene.data.binning`: tracks into the bins of a grid with the cuts of the model, nothing lost.

Roadmap 0.1-15. The rules — an edge belongs to the bin above it, the circle is inside, a track
whose line sees the sky only from behind the plane is left out and counted — are checked on
tracks placed by hand. Tracks drawn from known maps and binned back must give the map: exactly,
when each is drawn inside its own bin, and within Poisson errors, when they are drawn over the
whole cone independently of the grid. The bins `bin_tracks` takes as accepted are those of
status 1 of a map of `examples/open_sky`, computed through Mulder; every other map here is made
by hand.
"""

import importlib.util
import math
import pathlib

import numpy as np
import pytest

from muoscene.data import BinnedTracks, bin_tracks
from muoscene.detectors import DirectionGrid, Pose, TrackSet
from muoscene.model import grid_of

GRID = DirectionGrid.square()
TILTED = Pose((0.0, 0.0, 0.0), 180.0, 60.0)  # looks south; the horizon is the line ty = cot 60°
SEED = 20261128
OPEN_SKY = pathlib.Path(__file__).resolve().parents[2] / "examples" / "open_sky" / "run.py"


def tracks(tx, ty, frame="detector"):
    return TrackSet(tx, ty, frame=frame, source="synthetic", detector="d")


def where(binned):
    """(ty bin, tx bin) of every bin with a track, with the count."""
    return {(int(j), int(i)): int(binned.counts[j, i]) for j, i in np.argwhere(binned.counts)}


def test_edges_belong_to_the_bin_above():
    edges = GRID.tx_edges
    binned = bin_tracks(
        tracks([edges[21], edges[0], edges[-1], 0.01, -0.01], [edges[19], edges[0], edges[-1], 0.0, -0.0]),
        GRID,
        tan_max=2.0,
    )
    # edges[21] = 0.05 opens bin 21, edges[19] = −0.05 bin 19; the first edge is in bin 0, the last
    # in bin 39; (0.01, 0) and (−0.01, −0) sit on ty = 0, the lower edge of bin 20
    assert where(binned) == {(19, 21): 1, (0, 0): 1, (39, 39): 1, (20, 20): 1, (20, 19): 1}
    assert binned.binned == 5 and binned.outside_grid == binned.outside_circle == 0


def test_the_circle_is_inside_and_what_is_left_out_is_counted():
    small = DirectionGrid.square(0.5)
    binned = bin_tracks(
        tracks([1.0, 1.0 + 1e-12, 0.7, 0.0, 0.6], [0.0, 0.0, 0.0, 0.0, 0.8]), small, tan_max=1.0
    )
    # (1, 0) is on the circle but beyond the grid; so is (0.7, 0); (0.6, 0.8) has sqrt(tx² + ty²) = 1
    assert binned.outside_circle == 1 and binned.outside_grid == 3 and binned.binned == 1
    on = bin_tracks(tracks([1.0, 0.0], [0.0, -1.0]), GRID, tan_max=1.0)
    assert on.binned == 2 and where(on) == {(20, 39): 1, (0, 20): 1}
    summary = binned.summary()
    assert summary["tracks"] == 5 == sum(summary[k] for k in ("binned", "outside_circle", "outside_grid"))


def test_behind_the_plane_and_on_the_horizon():
    horizon = 1.0 / math.tan(math.radians(60.0))
    ty = [0.0, 0.5, horizon - 1e-9, horizon + 1e-9, 0.9, -0.9]
    binned = bin_tracks(tracks(np.zeros(6), ty), GRID, tan_max=1.0, pose=TILTED)
    assert binned.binned == 4 and binned.behind == 2 and binned.not_accepted == 0
    assert TILTED.direction(0.0, horizon + 1e-9)[2] < 0 < TILTED.direction(0.0, horizon - 1e-9)[2]
    level = Pose((0.0, 0.0, 0.0))
    rng = np.random.default_rng(SEED)
    many = tracks(rng.uniform(-1.2, 1.2, 5000), rng.uniform(-1.2, 1.2, 5000))
    with_pose, without = bin_tracks(many, GRID, tan_max=1.0, pose=level), bin_tracks(many, GRID, tan_max=1.0)
    assert np.array_equal(with_pose.counts, without.counts) and np.array_equal(
        with_pose.accepted, without.accepted
    )
    assert with_pose.behind == 0 and with_pose.supersample == 4 and without.supersample is None


def test_what_comes_back():
    binned = bin_tracks(tracks([0.1, 0.2], [0.0, 0.0]), GRID, tan_max=1.0)
    assert (
        isinstance(binned, BinnedTracks)
        and binned.counts.dtype == np.int64
        and binned.counts.shape == (40, 40)
    )
    assert int(binned.accepted.sum()) == 1324  # 1176 bins in the circle and 148 cut by it (detectors.grid)
    assert np.array_equal(binned.status, np.where(binned.accepted, 0, 1)) and binned.status.dtype == np.int8
    with pytest.raises(ValueError):
        binned.counts[20, 22] = 0
    empty = bin_tracks(tracks([], []), GRID, tan_max=1.0, pose=TILTED)
    assert empty.n_tracks == empty.binned == 0


def test_wrong_input_is_refused():
    scan = tracks([0.1], [0.0], frame="scan")
    assert bin_tracks(scan, GRID, tan_max=1.0).frame == "scan"
    with pytest.raises(ValueError, match="frame of a scan"):
        bin_tracks(scan, GRID, tan_max=1.0, pose=TILTED)
    for wrong in (0.0, -1.0, math.inf, True, "1"):
        with pytest.raises(ValueError):
            bin_tracks(scan, GRID, tan_max=wrong)
    with pytest.raises(TypeError):
        bin_tracks([0.1], GRID, tan_max=1.0)
    with pytest.raises(TypeError):
        bin_tracks(scan, GRID.tx_edges, tan_max=1.0)
    with pytest.raises(ValueError, match="once"):  # a count that does not add up
        BinnedTracks(
            np.zeros(GRID.shape), np.ones(GRID.shape, bool), GRID, 1.0, None, None, "scan", "d", 1, 0, 0, 0, 0
        )


def in_bins(rng, counts, grid, tan_max):
    """Tracks drawn bin by bin: n tracks in a bin, uniform in its part inside the circle."""
    tx, ty = [], []
    edges_x, edges_y = grid.tx_edges, grid.ty_edges
    for (j, i), n in np.ndenumerate(counts):
        got_x, got_y = np.zeros(0), np.zeros(0)
        while got_x.size < n:
            x = rng.uniform(edges_x[i], edges_x[i + 1], 4 * n + 8)
            y = rng.uniform(edges_y[j], edges_y[j + 1], 4 * n + 8)
            keep = np.hypot(x, y) <= tan_max
            got_x, got_y = np.append(got_x, x[keep]), np.append(got_y, y[keep])
        tx.append(got_x[:n])
        ty.append(got_y[:n])
    return np.concatenate(tx), np.concatenate(ty)


def test_tracks_drawn_in_their_bins_come_back_to_them():
    """A map made by hand, one Poisson draw per bin, tracks inside the bins: the counts come back exactly."""
    rng = np.random.default_rng(SEED)
    tx, ty = np.meshgrid(GRID.tx_centers, GRID.ty_centers)
    shape = 40.0 * (1.0 + tx + 0.5 * ty**2)  # not symmetric in tx, ty or their signs
    expected = np.where(GRID.accepted_solid_angle_sr(1.0) > 0, shape, 0.0)
    drawn = rng.poisson(expected)
    x, y = in_bins(rng, drawn, GRID, 1.0)
    order = rng.permutation(x.size)
    binned = bin_tracks(tracks(x[order], y[order]), GRID, tan_max=1.0)
    assert np.array_equal(binned.counts, drawn) and binned.binned == drawn.sum() > 5e4


def test_tracks_drawn_over_the_cone_agree_with_the_map_within_poisson_errors():
    """Tracks drawn over the whole cone, independently of the grid, and binned back give the map.

    The map: an isotropic sky over the cone of 45°, halved in a block of 8 × 8 bins at tx in
    [0.2, 0.6), ty in [−0.5, −0.1), so that a swap of tx and ty or of a sign shows. The expected
    count of a bin is the number drawn times its exact accepted solid angle over that of the cone.
    The statistic is Pearson's χ² over the K bins expecting at least 10 tracks; for Poisson counts
    its mean is K and its variance Σ(2 + 1/λ). The test fails when it is more than 4 standard
    deviations from K: for a random seed about 6·10⁻⁵ with a normal tail, 8·10⁻⁵ with the skew
    of a χ² of 1300 degrees of freedom (Wilson–Hilferty). The seed is fixed, so the outcome is too:
    K = 1296, χ² = 1315.5, 0.38 standard deviations above K; the two controls, tx and ty swapped
    and tx mirrored, are 192 and 183 standard deviations above it.
    """
    rng = np.random.default_rng(SEED)
    theta_max = math.atan(1.0)
    n = rng.poisson(2.0e5)
    cos_theta = rng.uniform(math.cos(theta_max), 1.0, n)
    phi = rng.uniform(0.0, 2.0 * math.pi, n)
    slope = np.sqrt(1.0 - cos_theta**2) / cos_theta
    tx, ty = slope * np.cos(phi), slope * np.sin(phi)
    e_x, e_y = GRID.tx_edges, GRID.ty_edges
    shaded = (tx >= e_x[24]) & (tx < e_x[32]) & (ty >= e_y[10]) & (ty < e_y[18])
    keep = ~shaded | (rng.uniform(size=n) < 0.5)
    binned = bin_tracks(tracks(tx[keep], ty[keep]), GRID, tan_max=1.0)

    omega = GRID.accepted_solid_angle_sr(1.0)
    assert omega.sum() == pytest.approx(2.0 * math.pi * (1.0 - math.cos(theta_max)), rel=1e-12)
    weight = np.ones(GRID.shape)
    weight[10:18, 24:32] = 0.5
    expected = 2.0e5 * omega / omega.sum() * weight
    assert binned.outside_circle == 0 and binned.outside_grid == 0 and binned.binned == keep.sum()

    def chi2(counts):
        big = expected >= 10.0
        k = int(big.sum())
        value = float((((counts - expected) ** 2) / np.where(big, expected, 1.0))[big].sum())
        sigma = math.sqrt(float((2.0 + 1.0 / expected[big]).sum()))
        return value, k, sigma

    value, k, sigma = chi2(binned.counts)
    assert k > 1250 and abs(value - k) <= 4.0 * sigma, (value, k, sigma)
    small = (expected > 0) & (expected < 10.0)
    assert abs(binned.counts[small].sum() - expected[small].sum()) <= 4.0 * math.sqrt(expected[small].sum())
    for control in (binned.counts.T, binned.counts[:, ::-1]):
        wrong, _, _ = chi2(control)
        assert wrong - k > 20.0 * sigma, wrong


def load_open_sky():
    spec = importlib.util.spec_from_file_location("open_sky_binning", OPEN_SKY)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def open_sky(tmp_path_factory):
    with pytest.MonkeyPatch.context() as patch:
        patch.setenv("MUOSCENE_CACHE", str(tmp_path_factory.mktemp("cache")))
        return load_open_sky().run(check_directions=0)


def test_the_bins_left_out_are_those_the_model_does_not_accept(open_sky):
    """For both films of `examples/open_sky`, the bins `bin_tracks` does not accept are those of status 1."""
    _, _, poses, _, maps = open_sky
    expected = maps[1.0]
    grid = grid_of(expected)
    rng = np.random.default_rng(SEED)
    many = tracks(rng.uniform(-1.0, 1.0, 200_000), rng.uniform(-1.0, 1.0, 200_000))
    for name, pose in poses.items():
        status = expected["status"].sel(detector=name).values[0]
        binned = bin_tracks(
            many, grid, tan_max=1.0, pose=pose, supersample=int(expected.attrs["supersample"])
        )
        assert np.array_equal(~binned.accepted, status == 1), name
        assert np.all(binned.counts[status == 1] == 0)
        inside = np.hypot(many.tx, many.ty) <= 1.0
        sky = pose.direction(many.tx, many.ty)[:, 2] > 0
        assert binned.behind == np.count_nonzero(inside & ~sky)
        assert binned.binned + binned.not_accepted == np.count_nonzero(inside & sky)
        if name == "flat":
            assert binned.behind == 0 and binned.not_accepted == 0
        else:
            # 476 bins not accepted (docs/notes/figures.md). The horizon, ty = cot 60° = 0.577, leaves
            # two of the four rows of sub-samples of the bins it crosses above it, so no track above
            # it falls in a bin the model leaves out. In the 34 bins it crosses the model counts
            # whole sub-samples, up to ty = 0.575, and the tracks are cut at the horizon itself
            assert int((status == 1).sum()) == 476 and binned.not_accepted == 0
            i = np.searchsorted(grid.tx_edges, many.tx, side="right") - 1
            j = np.searchsorted(grid.ty_edges, many.ty, side="right") - 1
            crossed = np.zeros(grid.shape, bool)
            crossed[j[inside & ~sky], i[inside & ~sky]] = True
            assert int((crossed & binned.accepted).sum()) == 34
