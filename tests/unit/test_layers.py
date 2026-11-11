"""Import rules: layers, optional dependencies and graphics."""

import ast
import importlib
import pathlib
import re

import pytest

from muoscene._layers import LAYERS

SRC = pathlib.Path(__file__).resolve().parents[2] / "src" / "muoscene"

# Installed only with extras: importing any muoscene module must work without them.
OPTIONAL = {
    "trimesh",
    "open3d",
    "laspy",
    "rasterio",
    "ezdxf",
    "pyvista",
    "uproot",
    "awkward",
    "mulder",
    "calzone",
    "MCEq",
}
# The computing core (layers 0-4) never draws.
GRAPHICS = {"matplotlib", "pyvista", "vtk"}
CORE_MAX_LAYER = 4


def top_module(path: pathlib.Path) -> str | None:
    """Top-level muoscene module a source file belongs to; None for package-level files."""
    name = path.relative_to(SRC).parts[0].removesuffix(".py")
    return name if name in LAYERS else None


def package_of(path: pathlib.Path) -> tuple[str, ...]:
    """Package of a source file relative to muoscene: ("geometry",) for geometry/scene.py."""
    return path.relative_to(SRC).parent.parts


def imported_modules(tree: ast.AST, package: tuple[str, ...]) -> set[str]:
    """Top-level muoscene modules imported anywhere in a syntax tree, relative imports included."""
    out: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                parts = alias.name.split(".")
                if parts[0] == "muoscene" and len(parts) > 1:
                    out.add(parts[1])
            continue
        if not isinstance(node, ast.ImportFrom):
            continue
        module = tuple(node.module.split(".")) if node.module else ()
        if node.level == 0:
            if not module or module[0] != "muoscene":
                continue
            target = module[1:]
        else:
            up = node.level - 1
            if up > len(package):
                continue
            target = package[: len(package) - up] + module
        if target:
            out.add(target[0])
        else:
            out.update(alias.name for alias in node.names)
    return {name for name in out if name in LAYERS}


def third_party_imports(tree: ast.AST) -> tuple[set[str], set[str]]:
    """Top-level names of absolute imports: (at module level, anywhere)."""
    at_module_level: set[str] = set()
    anywhere: set[str] = set()

    def visit(node: ast.AST, in_function: bool) -> None:
        for child in ast.iter_child_nodes(node):
            names: list[str] = []
            if isinstance(child, ast.Import):
                names = [alias.name.split(".")[0] for alias in child.names]
            elif isinstance(child, ast.ImportFrom) and child.level == 0 and child.module:
                names = [child.module.split(".")[0]]
            anywhere.update(names)
            if not in_function:
                at_module_level.update(names)
            inner = in_function or isinstance(child, ast.FunctionDef | ast.AsyncFunctionDef | ast.Lambda)
            visit(child, inner)

    visit(tree, False)
    return at_module_level, anywhere


def sources():
    for path in sorted(SRC.rglob("*.py")):
        yield path, ast.parse(path.read_text(encoding="utf-8"))


def test_layers_respected():
    bad = []
    for path, tree in sources():
        own = top_module(path)
        if own is None:
            continue
        for dep in imported_modules(tree, package_of(path)):
            if LAYERS[dep] > LAYERS[own]:
                bad.append(
                    f"{path.relative_to(SRC)}: {own} (layer {LAYERS[own]}) "
                    f"imports {dep} (layer {LAYERS[dep]})"
                )
    assert not bad, "\n".join(bad)


def test_optional_dependencies_are_imported_lazily():
    bad = []
    for path, tree in sources():
        at_module_level, _ = third_party_imports(tree)
        for name in sorted(at_module_level & OPTIONAL):
            bad.append(f"{path.relative_to(SRC)}: {name} imported at module level")
    assert not bad, "\n".join(bad)


def test_core_does_not_draw():
    bad = []
    for path, tree in sources():
        own = top_module(path)
        if own is None or LAYERS[own] > CORE_MAX_LAYER:
            continue
        _, anywhere = third_party_imports(tree)
        for name in sorted(anywhere & GRAPHICS):
            bad.append(f"{path.relative_to(SRC)}: {own} (layer {LAYERS[own]}) imports {name}")
    assert not bad, "\n".join(bad)


def test_every_package_has_a_layer():
    names = {
        p.name.removesuffix(".py")
        for p in SRC.iterdir()
        if not p.name.startswith("_") and p.name != "py.typed"
    }
    assert names == set(LAYERS), names ^ set(LAYERS)


def test_docstrings_state_the_layer():
    bad = []
    for name, layer in LAYERS.items():
        doc = importlib.import_module(f"muoscene.{name}").__doc__ or ""
        found = re.search(r"Layer (\d+)", doc)
        if found is None or int(found.group(1)) != layer:
            bad.append(name)
    assert not bad, bad


@pytest.mark.parametrize(
    ("source", "package", "expected"),
    [
        ("import muoscene.viz", ("geometry",), {"viz"}),
        ("import muoscene.viz.maps as m", (), {"viz"}),
        ("from muoscene import viz, flux", ("model",), {"viz", "flux"}),
        ("from muoscene.plan.significance import z", ("model",), {"plan"}),
        ("from . import scene", ("geometry",), {"geometry"}),
        ("from .scene import Scene", ("geometry",), {"geometry"}),
        ("from .. import viz", ("geometry",), {"viz"}),
        ("from ..viz import maps", ("geometry",), {"viz"}),
        ("from ..viz.maps import f", ("geometry", "sub"), {"geometry"}),
        ("from ...viz.maps import f", ("geometry", "sub"), {"viz"}),
        ("from . import viz", (), {"viz"}),
        ("from .viz.maps import f", (), {"viz"}),
        ("def f():\n    from muoscene.viz import maps\n    return maps", ("model",), {"viz"}),
        ("import numpy as np\nfrom pathlib import Path", ("model",), set()),
        ("from muoscene import __version__", (), set()),
    ],
)
def test_import_scanner(source, package, expected):
    assert imported_modules(ast.parse(source), package) == expected


def test_third_party_scanner():
    tree = ast.parse(
        "import numpy as np\n"
        "from scipy import optimize\n"
        "def f():\n"
        "    import pyvista\n"
        "    return pyvista\n"
        "class A:\n"
        "    import trimesh\n"
    )
    at_module_level, anywhere = third_party_imports(tree)
    assert at_module_level == {"numpy", "scipy", "trimesh"}
    assert anywhere == {"numpy", "scipy", "trimesh", "pyvista"}
