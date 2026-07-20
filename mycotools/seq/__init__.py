#! /usr/bin/env python3
"""Dispatcher for the `mycotools seq` subcommand.

Routes `mycotools seq <TOOL> ...` to a sequence/coordinate transform module."""
from mycotools.lib.subcmd import Dispatcher

# subcommand name/alias -> submodule within this package (mycotools.seq.<module>)
SUBCOMMANDS = {
    "translate": "translate",
    "coords": "coords",
    "gff": "gff",
    "convert": "convert",
    "mass": "mass",
}

DESCRIPTION = """Sequence and coordinate transforms

Tools (all following arguments are forwarded to the tool):
  translate   translate a nucleotide fasta to protein
  coords      extract subsequence(s) from a fasta by coordinate
  gff         extract sequences from a gff3 (+ optional assembly)
  convert     convert between sequence file formats
  mass        compute protein masses from an amino-acid fasta

Examples:
  mycotools seq translate -h
  mycotools seq coords -h"""

_dispatcher = Dispatcher(
    "mycotools seq", "mycotools.seq", SUBCOMMANDS, DESCRIPTION,
    metavar="TOOL", arg_help="sequence tool (see below)",
)
main = _dispatcher.main
cli = _dispatcher.cli


if __name__ == "__main__":
    cli()
