"""muoscene — muography scene: plan and check absorption muography of engineering-scale objects.

A site model (terrain and canyon wall, soil layers, structures, niches) with detectors placed in it
(nuclear emulsion, straw telescopes) gives expected angular maps, significance versus exposure
days and the volume actually seen, and compares them with scanned data using the same code.

Importing the package gives the process a `HOME` when it has none, as on Windows: Mulder 0.3.7
finds its cache of physics tables through `$HOME` alone and otherwise fails every calculation
with "bad cache (could not resolve $HOME)". The value is the home directory of the user, so the
cache of Mulder is `.cache/mulder` under it on every system, next to that of muoscene.
"""

import os as _os
import pathlib as _pathlib

if "HOME" not in _os.environ:
    try:
        _os.environ["HOME"] = str(_pathlib.Path.home())
    except RuntimeError:  # no home directory at all: Mulder says so itself when it needs its cache
        pass

__version__ = "0.0.1.dev0"
