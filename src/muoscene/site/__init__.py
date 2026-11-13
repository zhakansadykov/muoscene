"""Site and time.

Local ENU frame (east, north, up; metres) anchored at the site origin given by latitude,
longitude and altitude; it is the Mulder local frame of the scene. Conversions from the survey
CRS (UTM, a contractor's local grid) and WGS84 via pyproj happen only here, at the I/O boundary.
True north, grid convergence and magnetic declination for a given date (IGRF through Mulder).
Exposure windows and season helpers.

Layer 1: may import only `muoscene` modules of layer <= 1 (see `muoscene._layers`).
"""
