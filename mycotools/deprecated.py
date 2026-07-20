#! /usr/bin/env python3
"""Backwards-compatibility shims for the pre-v2 flat command names.

Every tool that moved under the nested `mycotools`/`mtdb` dispatchers keeps its
old flat console command (e.g. `jgiDwnld`, `acc2fa`) working for one transition
period. Each old command resolves to a shim here that prints a deprecation
notice to stderr and then hands off to the tool's relocated `cli()` unchanged;
`sys.argv` is untouched, so argument parsing behaves exactly as before.

These direct entrypoints are TEMPORARY and will be removed in a subsequent
release -- use the nested invocation printed in the warning instead."""
import sys
import importlib


def _make_shim(old, new, target):
    """Build a console-script shim: warn, then delegate to `target`'s cli().

    old    : the deprecated flat command name
    new    : the nested invocation the user should switch to
    target : importable module exposing `cli()` (the relocated tool)
    """

    def shim():
        sys.stderr.write(
            f"WARNING: `{old}` is a deprecated entrypoint and will be removed in "
            f"a subsequent release. Use `{new}` instead.\n"
        )
        importlib.import_module(target).cli()

    shim.__name__ = old
    shim.__qualname__ = old
    shim.__doc__ = f"Deprecated alias for `{new}`."
    return shim


# --- relocated under the `mycotools` analysis entrypoint ---------------------
jgiDwnld = _make_shim("jgiDwnld", "mycotools download jgi", "mycotools.download.jgi")
ncbiDwnld = _make_shim("ncbiDwnld", "mycotools download ncbi", "mycotools.download.ncbi")
db2search = _make_shim("db2search", "mycotools homology db", "mycotools.homology.db")
fa2hmmer2fa = _make_shim(
    "fa2hmmer2fa", "mycotools homology fasta", "mycotools.homology.fasta"
)
db2hgs = _make_shim("db2hgs", "mycotools cluster db", "mycotools.cluster.db")
fa2clus = _make_shim("fa2clus", "mycotools cluster fasta", "mycotools.cluster.fasta")
crap = _make_shim("crap", "mycotools phylo crap", "mycotools.phylo.crap")
fa2tree = _make_shim("fa2tree", "mycotools phylo tree", "mycotools.phylo.tree")
db2microsyntree = _make_shim(
    "db2microsyntree", "mycotools phylo synteny", "mycotools.phylo.synteny"
)
annotationStats = _make_shim(
    "annotationStats", "mycotools stats annotation", "mycotools.stats.annotation"
)
assemblyStats = _make_shim(
    "assemblyStats", "mycotools stats assembly", "mycotools.stats.assembly"
)
fna2faa = _make_shim("fna2faa", "mycotools seq translate", "mycotools.seq.translate")
coords2fa = _make_shim("coords2fa", "mycotools seq coords", "mycotools.seq.coords")
gff2seq = _make_shim("gff2seq", "mycotools seq gff", "mycotools.seq.gff")
bioreform = _make_shim("bioreform", "mycotools seq convert", "mycotools.seq.convert")
fa2mass = _make_shim("fa2mass", "mycotools seq mass", "mycotools.seq.mass")
add2gff = _make_shim("add2gff", "mycotools gff add", "mycotools.gff.add")
gff2svg = _make_shim("gff2svg", "mycotools gff svg", "mycotools.gff.svg")

# --- relocated under the `mtdb` database entrypoint --------------------------
db2files = _make_shim("db2files", "mtdb files", "mycotools.mtdb.files")

# --- relocated under the `mycotools` analysis entrypoint ---------------------
ome2name = _make_shim("ome2name", "mycotools rename", "mycotools.rename")

# --- restored acc2* aliases (retired in v2, now `mtdb accession <FORMAT>`) ----
acc2fa = _make_shim("acc2fa", "mtdb accession fa", "mycotools.mtdb.acc2.fa")
acc2gff = _make_shim("acc2gff", "mtdb accession gff", "mycotools.mtdb.acc2.gff")
acc2gbk = _make_shim("acc2gbk", "mtdb accession gbk", "mycotools.mtdb.acc2.gbk")
acc2locus = _make_shim("acc2locus", "mtdb accession locus", "mycotools.mtdb.acc2.locus")
