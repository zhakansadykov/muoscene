"""transport.direct — direct path: Mulder continuous flux through the built scene; check of the fast path.

The fast path — scan, equivalent opacity, table lookup — is compared with Mulder itself in
every run: on a random sample of the lines of sight of a detector the flux of the table is set
against the continuous flux Mulder carries through the scene as it was built, integrated over
the energy response of the detector, and the deviation is stored with the result.

How the direct flux reaches Mulder (decided on 21 November 2026): through `Site.differential_flux`,
which takes points and directions in ENU, turns them into the axes of the scene and has
`Scene.scan` refuse what Mulder answers with a crash. Nothing here rotates a vector.

The integral over the response is a trapezoid sum in ln K that starts at the first energy of the
response itself, 30 nodes per decade up to 10 TeV with the nodes of the response among them:
against 240 nodes per decade it is 7e-5 off above 1 GeV and 3e-5 above 0.2 GeV under 12 m of
soil. A sum that gives its first node a whole step is 2–3 % high (`docs/notes/flux.md`).

In the first half decade from the threshold the nodes are 240 to a decade (22 November 2026,
roadmap 0.1-10). With a sea-level reference the flux of the open sky falls steeply from a
threshold just above its cut-off, and 30 nodes there put the direct flux above 1 GeV at 15.5
degrees of elevation at 1150 m 0.62 % above the table, which is right: the check refused a run
it should have passed. With the denser start the difference is 0.02 %, against 0.01 % with 2000
nodes a decade (`tests/unit/test_transport_direct.py`); a line takes 226 nodes instead of 121.

Mulder discrete flux, with scattering, comes with version 0.3.

Units: kinetic energy in GeV, flux in m⁻² s⁻¹ sr⁻¹; points in m and directions as unit vectors in ENU.

Layer 2.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from ..flux import EnergyResponse
from ..geometry import Scene
from ..site import Site

K_MAX_GEV = 1.0e4
"""Top of the integral over energy: the highest energy at the detector of a `MulderTable`."""
NODES_PER_DECADE = 30
"""Energy nodes per decade of the integral over the response."""
FIRST_DECADES = 0.5
"""From the first energy of the response up by this many decades the nodes are denser:"""
FIRST_NODES_PER_DECADE = 240
"""Nodes per decade there: near the cut-off of a sea-level reference the flux falls steeply from
the threshold."""
LIMIT = 5.0e-3
"""Largest deviation of the fast path from the direct one a run accepts: the 0.5 % of the roadmap."""


def energy_nodes(
    response: EnergyResponse,
    k_max_gev: float = K_MAX_GEV,
    nodes_per_decade: int = NODES_PER_DECADE,
    first_nodes_per_decade: int = FIRST_NODES_PER_DECADE,
) -> np.ndarray:
    """Kinetic energies (GeV) of the integral over a response: from its first node to `k_max_gev`.

    Equal steps in ln K: at least `first_nodes_per_decade` to a decade over the first
    `FIRST_DECADES` from the first node, at least `nodes_per_decade` beyond, and the other nodes
    of the response among them.
    """
    if not isinstance(response, EnergyResponse):
        raise TypeError("the direct flux is integrated over an EnergyResponse")
    for name, value in (
        ("nodes_per_decade", nodes_per_decade),
        ("first_nodes_per_decade", first_nodes_per_decade),
    ):
        if isinstance(value, bool) or not isinstance(value, int) or value < 4:
            raise ValueError(f"{name} must be an integer of at least 4")
    first = response.k_min_gev
    if not (math.isfinite(k_max_gev) and k_max_gev > first):
        raise ValueError("k_max_gev must be above the first energy of the response")
    top = min(first * 10.0**FIRST_DECADES, k_max_gev)
    dense = max(nodes_per_decade, first_nodes_per_decade)
    count = int(math.ceil(math.log10(top / first) * dense)) + 1
    k = np.logspace(math.log10(first), math.log10(top), max(count, 2))
    if top < k_max_gev:
        count = int(math.ceil(math.log10(k_max_gev / top) * nodes_per_decade)) + 1
        k = np.concatenate([k, np.logspace(math.log10(top), math.log10(k_max_gev), max(count, 2))[1:]])
    k[0], k[-1] = first, k_max_gev
    more = [e for e in response.kinetic_energy_gev[1:] if e < k_max_gev]
    return np.unique(np.concatenate([k, more]))


def direct_flux(
    site: Site,
    scene: Scene,
    position_m,
    direction,
    response: EnergyResponse,
    *,
    k_max_gev: float = K_MAX_GEV,
    nodes_per_decade: int = NODES_PER_DECADE,
) -> np.ndarray:
    """Flux of Mulder through a scene over the energy response of a detector, m⁻² s⁻¹ sr⁻¹.

    site, scene     the site and a scene of it, with the densities it was built with
    position_m      points in ENU, m; shape (..., 3) or (3,)
    direction       unit vectors in ENU from the detector towards the sky; shape (..., 3)
    response        efficiency against the kinetic energy at the detector

    Continuous transport along straight lines of sight; the integral of efficiency times
    differential flux by trapezoids in ln K over `energy_nodes`. Returns one value per line of
    sight. Mulder gives a flux of zero below the cut-off of the reference: a response that
    starts there is integrated across it, coarsely, and the caller keeps such lines out.
    """
    k = energy_nodes(response, k_max_gev, nodes_per_decade)
    flux = site.differential_flux(scene, position_m, direction, k)
    return np.trapezoid(flux * (np.asarray(response(k)) * k), np.log(k), axis=-1)


def seed_of(rng: np.random.Generator) -> int:
    """The seed a generator was made from, to be written with a result.

    That of `numpy.random.default_rng(seed)`, or the entropy numpy drew for a generator made
    without one. A generator that does not know its seed is refused: a sample that cannot be
    drawn again is not stored.
    """
    if not isinstance(rng, np.random.Generator):
        raise TypeError("randomness comes from a numpy.random.Generator")
    sequence = getattr(rng.bit_generator, "seed_seq", None)
    entropy = getattr(sequence, "entropy", None)
    if isinstance(entropy, bool) or not isinstance(entropy, int):
        raise ValueError("the generator does not know its seed: make it with numpy.random.default_rng(seed)")
    if getattr(sequence, "spawn_key", ()):
        raise ValueError(
            "a spawned generator has no single seed: make it with numpy.random.default_rng(seed)"
        )
    return entropy


@dataclass(frozen=True)
class DirectCheck:
    """The fast path against the direct flux of Mulder on a sample of lines of sight.

    fast_m2_s_sr     flux of the fast path along each line, m⁻² s⁻¹ sr⁻¹
    direct_m2_s_sr   flux of Mulder through the scene along the same lines
    mean, rms, max   of the relative deviation fast / direct − 1: its mean, its root mean
                     square and its largest absolute value
    worst            index of the line with the largest deviation
    limit            deviation above which the run was to be refused, or None
    """

    fast_m2_s_sr: np.ndarray
    direct_m2_s_sr: np.ndarray
    mean: float
    rms: float
    max: float
    worst: int
    limit: float | None

    @property
    def n_directions(self) -> int:
        """Number of lines of sight in the sample."""
        return self.fast_m2_s_sr.size

    @property
    def deviation(self) -> np.ndarray:
        """Relative deviation of the fast path along each line: fast / direct − 1."""
        return self.fast_m2_s_sr / self.direct_m2_s_sr - 1.0


def compare_with_direct(
    site: Site,
    scene: Scene,
    position_m,
    direction,
    response: EnergyResponse,
    fast_flux_m2_s_sr,
    *,
    limit: float | None = LIMIT,
    k_max_gev: float = K_MAX_GEV,
    nodes_per_decade: int = NODES_PER_DECADE,
) -> DirectCheck:
    """Compare the flux of the fast path along lines of sight with Mulder through the scene.

    direction            unit vectors in ENU; shape (n, 3), at least one
    fast_flux_m2_s_sr    flux of the fast path along each, over the same response; shape (n,)
    limit                largest relative deviation accepted; None only records

    RuntimeError if the largest deviation exceeds `limit`: the table and the equivalent
    opacity do not stand for this scene, and a map made with them would be wrong without a sign
    of it.
    """
    fast = np.array(fast_flux_m2_s_sr, dtype=float)
    direction = np.asarray(direction, dtype=float)
    if fast.ndim != 1 or fast.size == 0 or direction.shape != (fast.size, 3):
        raise ValueError("one flux of the fast path for each of at least one line of sight")
    if limit is not None and not (math.isfinite(limit) and limit > 0):
        raise ValueError("limit must be positive, or None")
    direct = direct_flux(
        site, scene, position_m, direction, response, k_max_gev=k_max_gev, nodes_per_decade=nodes_per_decade
    )
    if not (np.all(np.isfinite(fast)) and np.all(direct > 0)):
        raise RuntimeError("no flux along a line of sight of the sample: nothing to compare with")
    ratio = fast / direct - 1.0
    worst = int(np.argmax(np.abs(ratio)))
    fast.flags.writeable = False
    direct.flags.writeable = False
    check = DirectCheck(
        fast,
        direct,
        float(ratio.mean()),
        float(np.sqrt(np.mean(ratio**2))),
        float(np.abs(ratio).max()),
        worst,
        limit,
    )
    if limit is not None and check.max > limit:
        raise RuntimeError(
            f"the fast path differs from the direct flux of Mulder by {100 * check.max:.2f} % on a line"
            f" at an elevation of {math.degrees(math.asin(float(direction[worst, 2]))):.1f} degrees"
            f" (rms {100 * check.rms:.2f} % over {fast.size} lines), above the limit of {100 * limit:g} %"
        )
    return check
