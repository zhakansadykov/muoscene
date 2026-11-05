"""A module may import only muoscene modules of its own or a lower layer."""

import ast
import pathlib

from muoscene._layers import LAYERS

SRC = pathlib.Path(__file__).resolve().parents[2] / "src" / "muoscene"


def top_module(path: pathlib.Path) -> str | None:
    rel = path.relative_to(SRC)
    name = rel.parts[0].removesuffix(".py")
    return name if name in LAYERS else None


def imported_modules(tree: ast.AST) -> set[str]:
    out = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for a in node.names:
                parts = a.name.split(".")
                if parts[0] == "muoscene" and len(parts) > 1:
                    out.add(parts[1])
        elif isinstance(node, ast.ImportFrom) and node.module and node.level == 0:
            parts = node.module.split(".")
            if parts[0] == "muoscene":
                if len(parts) > 1:
                    out.add(parts[1])
                else:
                    out.update(a.name for a in node.names)
    return {m for m in out if m in LAYERS}


def test_layers_respected():
    bad = []
    for path in SRC.rglob("*.py"):
        own = top_module(path)
        if own is None:
            continue
        for dep in imported_modules(ast.parse(path.read_text(encoding="utf-8"))):
            if LAYERS[dep] > LAYERS[own]:
                bad.append(
                    f"{path.relative_to(SRC)}: {own} (layer {LAYERS[own]}) "
                    f"imports {dep} (layer {LAYERS[dep]})"
                )
    assert not bad, "\n".join(bad)


def test_every_package_has_a_layer():
    names = {
        p.name.removesuffix(".py")
        for p in SRC.iterdir()
        if not p.name.startswith("_") and p.name != "py.typed"
    }
    assert names <= set(LAYERS), names - set(LAYERS)
