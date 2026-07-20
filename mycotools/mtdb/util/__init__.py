#! /usr/bin/env python3
"""Dispatcher for the `mtdb util` subcommand.

Routes `mtdb util <TOOL> ...` to a database-access helper. These tools read the
primary MTDB (they do not curate it), so they live under `mtdb` alongside the
database lifecycle commands."""
from mycotools.lib.subcmd import Dispatcher

# subcommand name/alias -> submodule within this package (mycotools.mtdb.util.<module>)
SUBCOMMANDS = {
    "files": "files",
    "name": "name",
}

DESCRIPTION = """Database-access helpers

Tools (all following arguments are forwarded to the tool):
  files   symlink/copy selected files from the database
  name    convert MTDB ome codes to taxonomic names

Examples:
  mtdb util files -h
  mtdb util name -h"""

_dispatcher = Dispatcher(
    "mtdb util", "mycotools.mtdb.util", SUBCOMMANDS, DESCRIPTION,
    metavar="TOOL", arg_help="utility (see below)",
)
main = _dispatcher.main
cli = _dispatcher.cli


if __name__ == "__main__":
    cli()
