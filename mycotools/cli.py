#! /usr/bin/env python3
"""Top-level dispatcher for the `mycotools` command.

`mycotools` is the downstream-analysis entrypoint; the database lifecycle
(update/extract/predb/manage/accession/util) lives under the separate `mtdb`
command. Routes `mycotools <GROUP> ...` to a group dispatcher (or, for `search`,
directly to the search module)."""
from mycotools.lib.subcmd import Dispatcher

# group name/alias -> submodule within the mycotools package. Groups (download,
# cluster, phylo, stats, seq, gff) are subpackage dispatchers; `search` is a
# single leaf module.
SUBCOMMANDS = {
    "download": "download",
    "dl": "download",
    "search": "search",
    "s": "search",
    "cluster": "cluster",
    "clus": "cluster",
    "phylo": "phylo",
    "p": "phylo",
    "stats": "stats",
    "seq": "seq",
    "gff": "gff",
}

DESCRIPTION = """Mycotools downstream-analysis toolkit

Groups (all following arguments are forwarded to the group/tool):
  download   (dl)     download genomes/annotations from JGI or NCBI
  search     (s)      search query sequence(s) against the database
  cluster    (clus)   cluster sequences / circumscribe homology groups
  phylo      (p)      phylogenies and phylogenetic pipelines (crap/tree/synteny)
  stats               annotation / assembly statistics
  seq                 sequence & coordinate transforms
  gff                 manipulate / render gff3 files

Database operations (build/manage/query the MTDB) live under the `mtdb` command;
run `mtdb -h`.

Examples:
  mycotools download jgi -h
  mycotools phylo crap -h
  mycotools search -h"""

_dispatcher = Dispatcher(
    "mycotools", "mycotools", SUBCOMMANDS, DESCRIPTION,
    metavar="GROUP", arg_help="analysis group or tool (see below)",
)
main = _dispatcher.main
cli = _dispatcher.cli


if __name__ == "__main__":
    cli()
