"""Muon flux at the detector as a function of what lies along the line of sight.

`FluxTable` holds the differential flux on a grid of kinetic energy, sine of elevation and
opacity, and integrates it over a detector energy response. `MulderTable` builds it with Mulder
for the detector altitude, an atmosphere model and optionally the geomagnetic field of a date;
the reference is GCCLY15 (Guan et al. 2015) at sea level, or a tabulated reference at the site
altitude, which is needed for muons below the altitude cut-off. `LegacyTable` reproduces the
preliminary calculation (sea level, a + bE) for regression tests. Tables are cached on disk.

Layer 1: may import only `muoscene` modules of layer <= 1 (see `muoscene._layers`).
"""
