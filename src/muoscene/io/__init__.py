"""Input and output.

Survey import (LAS/LAZ, CSV points, LandXML, DXF, GeoTIFF); point cloud -> closed mesh;
track readers that fill `muoscene.detectors.tracks.TrackSet` — the text output of the FEDRA macro
`evaluate_VT_nseg.C` (`read_fedra_text`, `ScanFrame`; `docs/formats.md`), later FEDRA ROOT via uproot;
temperature/humidity loggers; export to NetCDF, VTK, CSV. The core receives objects, not files.

Layer 5: may import only `muoscene` modules of layer <= 5 (see `muoscene._layers`).
"""
