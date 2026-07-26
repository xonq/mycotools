#! /usr/bin/env python3
"""Dispatcher for the `mycotools download` subcommand.

Routes `mycotools download <SOURCE> ...` to a per-source retrieval module. Note
that `mtdb update` already assimilates JGI/NCBI genomes into the primary MTDB;
these commands are for fetching genomes/annotations directly as files."""
from mycotools.lib.subcmd import Dispatcher

# subcommand name/alias -> submodule within this package (mycotools.download.<module>)
SUBCOMMANDS = {
    "jgi": "jgi",
    "ncbi": "ncbi",
}

DESCRIPTION = """Download genomes/annotations from external repositories

Sources (all following arguments are forwarded to the source's downloader):
  jgi     download from JGI MycoCosm
  ncbi    download from NCBI GenBank/RefSeq

Examples:
  mycotools download jgi -h
  mycotools download ncbi -h"""

_dispatcher = Dispatcher(
    "mycotools download",
    "mycotools.download",
    SUBCOMMANDS,
    DESCRIPTION,
    metavar="SOURCE",
    arg_help="download source (see below)",
)
main = _dispatcher.main
cli = _dispatcher.cli


if __name__ == "__main__":
    cli()
