import importlib

from muoscene._layers import LAYERS


def test_all_layer_modules_import():
    for name in LAYERS:
        importlib.import_module(f"muoscene.{name}")
