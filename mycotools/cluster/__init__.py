#! /usr/bin/env python3
"""Dispatcher for the `mycotools cluster` subcommand.

Routes `mycotools cluster <METHOD> ...` to a per-method clustering module."""
from mycotools.lib.subcmd import Dispatcher

# subcommand name/alias -> submodule within this package (mycotools.cluster.<module>)
SUBCOMMANDS = {
    "homology": "homology",
    "hgs": "homology",
    "fasta": "fasta",
    "fa": "fasta",
}

DESCRIPTION = """Group sequences into clusters / homology groups

Methods (all following arguments are forwarded to the method):
  homology  (hgs)   circumscribe database sequences into homology groups
  fasta     (fa)    iterative sequence-similarity clustering of a fasta

Examples:
  mycotools cluster homology -h
  mycotools cluster fasta -h"""

_dispatcher = Dispatcher(
    "mycotools cluster", "mycotools.cluster", SUBCOMMANDS, DESCRIPTION,
    metavar="METHOD", arg_help="clustering method (see below)",
)
main = _dispatcher.main
cli = _dispatcher.cli


if __name__ == "__main__":
    cli()
