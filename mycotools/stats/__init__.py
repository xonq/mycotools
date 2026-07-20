#! /usr/bin/env python3
"""Dispatcher for the `mycotools stats` subcommand.

Routes `mycotools stats <KIND> ...` to a statistics module."""
from mycotools.lib.subcmd import Dispatcher

# subcommand name/alias -> submodule within this package (mycotools.stats.<module>)
SUBCOMMANDS = {
    "annotation": "annotation",
    "assembly": "assembly",
}

DESCRIPTION = """Summary statistics for annotations and assemblies

Kinds (all following arguments are forwarded to the kind):
  annotation   gene annotation statistics (gff3/gtf or MTDB)
  assembly     genome assembly statistics

Examples:
  mycotools stats annotation -h
  mycotools stats assembly -h"""

_dispatcher = Dispatcher(
    "mycotools stats", "mycotools.stats", SUBCOMMANDS, DESCRIPTION,
    metavar="KIND", arg_help="statistic to compute (see below)",
)
main = _dispatcher.main
cli = _dispatcher.cli


if __name__ == "__main__":
    cli()
