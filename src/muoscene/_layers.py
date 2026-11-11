"""Import layers. A module may import only muoscene modules of the same or a lower layer.

Checked by tests/unit/test_layers.py.
"""

LAYERS = {
    "units": 0,
    "config": 0,
    "site": 1,
    "materials": 1,
    "geometry": 1,
    "flux": 1,
    "transport": 2,
    "detectors": 2,
    "model": 3,
    "plan": 4,
    "synth": 4,
    "data": 4,
    "invert": 4,
    "io": 5,
    "viz": 5,
    "report": 5,
    "provenance": 5,
    "pipeline": 5,
    "cli": 5,
}
