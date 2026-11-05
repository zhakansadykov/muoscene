"""Detectors.

Detector interface: pose in ENU; direction grid in the detector frame (tan theta_x, tan theta_y,
solid angles) and its transform to ENU; acceptance, efficiency eps(direction, time), blur kernel,
background. Emulsion stack: cone, angular efficiency, latent-image fading from temperature and
humidity logs, track-density budget (exposure plus transport and storage), number of films and
scanning cost. Straw telescope (level-1 response). Generic two-plane telescope.

Layer 2: may import only `muoscene` modules of layer <= 2 (see `muoscene._layers`).
"""
