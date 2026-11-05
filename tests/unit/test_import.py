import importlib

from muoscene._layers import LAYERS


def test_all_layer_modules_import():
    for name in LAYERS:
        importlib.import_module(f"muoscene.{name}")


def test_cli_help(capsys):
    from muoscene.cli import main

    assert main([]) == 0
