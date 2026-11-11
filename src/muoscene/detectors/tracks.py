"""detectors.tracks — `TrackSet`, the common container of measured and synthetic tracks.

Slopes (tx, ty) and positions in the detector frame, film index, quality and the source (data
or synthetic). It is defined here, in layer 2, because `synth` and `data` (layer 4) produce and
consume it; the readers in `muoscene.io.tracks` only fill it. Not implemented yet.
"""
