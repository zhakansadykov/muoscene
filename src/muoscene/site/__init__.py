"""Site and time.

Local ENU frame (east, north, up; metres) anchored at the site origin given by latitude,
longitude and altitude. Conversions from the survey CRS (UTM, a contractor's local grid) and
WGS84 via pyproj happen only here, at the I/O boundary. True north, grid convergence and
magnetic declination for a given date. Exposure windows and season helpers.

Layer 1: may import only `muoscene` modules of layer <= 1 (see `muoscene._layers`).
"""
