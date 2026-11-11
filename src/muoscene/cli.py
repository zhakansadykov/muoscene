"""Command-line interface: `muoscene check|import-survey|scene|plan|synth|compare|invert`.

A thin wrapper: arguments are parsed here, the work is done by `muoscene.pipeline`.

Layer 5.
"""

import argparse
import sys

from muoscene import __version__

COMMANDS = {
    "check": "validate configuration and scene geometry",
    "import-survey": "convert a survey into a closed mesh in the site ENU frame",
    "scene": "render the 3D scene only (HTML and VTK)",
    "plan": "expected maps, significance, exposure days, sensitivity volume, report",
    "synth": "synthetic tracks from expected maps",
    "compare": "measured tracks against the model",
    "invert": "density inversion on the voxel grid (version 2)",
}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog="muoscene", description="Muography scene: planning and data comparison."
    )
    parser.add_argument("--version", action="version", version=f"muoscene {__version__}")
    sub = parser.add_subparsers(dest="command", metavar="command")
    for name, help_text in COMMANDS.items():
        p = sub.add_parser(name, help=help_text)
        p.add_argument("config", nargs="?", help="run configuration (YAML)")
    args = parser.parse_args(argv)
    if args.command is None:
        parser.print_help()
        return 0
    print(f"muoscene {args.command}: not implemented in {__version__}", file=sys.stderr)
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
