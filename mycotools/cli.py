#! /usr/bin/env python3
"""Top-level dispatcher for the `mycotools` command.

`mycotools` is the downstream-analysis entrypoint; the database lifecycle
(update/extract/predb/manage/accession/files) lives under the separate `mtdb`
command. Routes `mycotools <GROUP> ...` to a group dispatcher, or, for a leaf
tool like `rename`, directly to that tool's module."""
from mycotools.lib.subcmd import Dispatcher

# group name/alias -> submodule within the mycotools package. Groups (download,
# homology, cluster, phylo, stats, seq, gff) are subpackage dispatchers;
# `rename` is a single leaf module.
SUBCOMMANDS = {
    "download": "download",
    "d": "download",
    "homology": "homology",
    "h": "homology",
    "cluster": "cluster",
    "c": "cluster",
    "phylo": "phylo",
    "p": "phylo",
    "rename": "rename",
    "r": "rename",
    "stats": "stats",
    "seq": "seq",
    "gff": "gff",
}

DESCRIPTION = """Mycotools downstream-analysis toolkit

Groups (all following arguments are forwarded to the group/tool):
  download   (d)      download genomes/annotations from JGI or NCBI
  homology   (h)      search query sequence(s) against the database or a fasta
  cluster    (c)      cluster sequences / circumscribe homology groups
  phylo      (p)      phylogenies and phylogenetic pipelines (crap/tree/synteny)
  rename     (r)      substitute MTDB ome codes with taxonomic names in a file
  stats               annotation / assembly statistics
  seq                 sequence & coordinate transforms
  gff                 manipulate / render gff3 files

Database operations (build/manage/query the MTDB) live under the `mtdb` command;
run `mtdb -h`.

Examples:
  mycotools download jgi -h
  mycotools phylo crap -h
  mycotools homology -h"""

_dispatcher = Dispatcher(
    "mycotools", "mycotools", SUBCOMMANDS, DESCRIPTION,
    metavar="GROUP", arg_help="analysis group or tool (see below)",
)
main = _dispatcher.main
cli = _dispatcher.cli


if __name__ == "__main__":
    cli()
