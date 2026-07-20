#! /usr/bin/env python3
"""Dispatcher for the `mycotools phylo` subcommand.

Routes `mycotools phylo <TOOL> ...` to a phylogenetics module."""
from mycotools.lib.subcmd import Dispatcher

# subcommand name/alias -> submodule within this package (mycotools.phylo.<module>)
SUBCOMMANDS = {
    "crap": "crap",
    "tree": "tree",
    "synteny": "synteny",
    "tools": "tools",
}

DESCRIPTION = """Build phylogenies and phylogenetic pipelines

Tools (all following arguments are forwarded to the tool):
  crap       Cluster Reconstruction and Phylogeny (CRAP) pipeline
  tree       build a phylogeny from a fasta (align -> trim -> infer)
  synteny    build a microsynteny tree
  tools      manipulate an existing phylogeny (root/prune/rename/strip support)

Examples:
  mycotools phylo crap -h
  mycotools phylo tools -h"""

_dispatcher = Dispatcher(
    "mycotools phylo",
    "mycotools.phylo",
    SUBCOMMANDS,
    DESCRIPTION,
    metavar="TOOL",
    arg_help="phylogenetics tool (see below)",
)
main = _dispatcher.main
cli = _dispatcher.cli


if __name__ == "__main__":
    cli()
