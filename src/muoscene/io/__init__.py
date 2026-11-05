"""Input and output.

Survey import (LAS/LAZ, CSV points, LandXML, DXF, GeoTIFF); point cloud -> closed mesh;
track readers (MISIS text, FEDRA ROOT via uproot) -> TrackSet; temperature/humidity loggers;
export to NetCDF, VTK, CSV. The core receives objects, not files.

Layer 5: may import only `muoscene` modules of layer <= 5 (see `muoscene._layers`).
"""
