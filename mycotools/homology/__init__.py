#! /usr/bin/env python3
"""Dispatcher for the `mycotools homology` subcommand.

Routes `mycotools homology <METHOD> ...` to a homology-search module: `db`
searches a query against the database (BLAST/mmseqs/hmmer), `fasta` runs
hmmsearch/nhmmer on a fasta and returns a fasta of hits."""
from mycotools.lib.subcmd import Dispatcher

# subcommand name/alias -> submodule within this package (mycotools.homology.<module>)
SUBCOMMANDS = {
    "db": "db",
    "fasta": "fasta",
    "fa": "fasta",
}

DESCRIPTION = """Search query sequence(s) against the database or a fasta

Methods (all following arguments are forwarded to the method):
  db              search a query against the database (BLAST/mmseqs/hmmer)
  fasta   (fa)    hmmsearch/nhmmer a fasta and return a fasta of hits

Examples:
  mycotools homology db -h
  mycotools homology fasta -h"""

_dispatcher = Dispatcher(
    "mycotools homology", "mycotools.homology", SUBCOMMANDS, DESCRIPTION,
    metavar="METHOD", arg_help="homology-search method (see below)",
)
main = _dispatcher.main
cli = _dispatcher.cli


if __name__ == "__main__":
    cli()
