#! /usr/bin/env python3
"""Dispatcher for the `mycotools gff` subcommand.

Routes `mycotools gff <TOOL> ...` to a gff manipulation/rendering module."""
from mycotools.lib.subcmd import Dispatcher

# subcommand name/alias -> submodule within this package (mycotools.gff.<module>)
SUBCOMMANDS = {
    "add": "add",
    "svg": "svg",
}

DESCRIPTION = """Manipulate and render gff3 files

Tools (all following arguments are forwarded to the tool):
  add    add an entry to an existing gff
  svg    render a gff locus to SVG

Examples:
  mycotools gff add -h
  mycotools gff svg -h"""

_dispatcher = Dispatcher(
    "mycotools gff", "mycotools.gff", SUBCOMMANDS, DESCRIPTION,
    metavar="TOOL", arg_help="gff tool (see below)",
)
main = _dispatcher.main
cli = _dispatcher.cli


if __name__ == "__main__":
    cli()
