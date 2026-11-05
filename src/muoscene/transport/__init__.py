"""Muon transport.

Fast path: energy threshold from opacity via range tables. Multiple-scattering blur kernel
(Highland). Optional adapter to PUMAS / Mulder used for validation, not as the default engine.

Layer 2: may import only `muoscene` modules of layer <= 2 (see `muoscene._layers`).
"""
