#! /usr/bin/env python3

# NEED to add a log option
# list of DBs, ability to change names, quickly change between
# list update dates
# report storage information
# report taxonomy data
# NEED to add option to export NCBI/JGI credentials
# NEED to pay attention to old ome versions

import re
import sys
import logging
import argparse
import importlib
from pathlib import Path
from mycotools.lib.kontools import format_path, setup_logging
from mycotools.lib.dbtools import (
    primary_db,
    mtdb_disconnect,
    mtdb_initialize,
    mtdb,
    parse_user_config,
)

logger = logging.getLogger(__name__)

# subcommand name/alias -> submodule within this package (mycotools.mtdb.<module>).
# `accession`/`a` targets the acc2 subpackage, which further dispatches by format.
SUBCOMMANDS = {
    "extract": "extract",
    "e": "extract",
    "update": "update",
    "u": "update",
    "predb2mtdb": "predb",
    "p": "predb",
    "manage": "manage",
    "m": "manage",
    "accession": "acc2",
    "a": "acc2",
    "files": "files",
    "f": "files",
}

DESCRIPTION = """MycotoolsDB (MTDB) utility

Run without arguments to print the primary MTDB path.

Subcommands (all following arguments are forwarded to the subcommand):
  extract     (e)   extract a sub-.mtdb file
  update      (u)   update / initialize the primary MTDB
  predb2mtdb  (p)   add local genomes to the primary MTDB
  manage      (m)   MTDB management utility
  accession   (a)   retrieve data for accession(s) by format (fa/gff/gbk/locus)
  files       (f)   symlink/copy selected files from the database

Ome lookup:
  mtdb <OME>[.gff3|.fna|.faa]   print an ome's row, or a specific file path"""


def get_version():
    """Return the `-v`/`--version` output string."""
    from importlib.metadata import version

    return f'Mycotools version {version("mycotools")}'


def build_parser():
    """Build the base `mtdb` argument parser.

    Subcommands and ome-lookup share one positional (`target`) followed by a
    REMAINDER: this lets subcommand arguments pass through verbatim to their
    subcommand module (native argparse subparsers cannot, as they intercept
    forwarded flags and reject arbitrary ome positionals)."""
    parser = argparse.ArgumentParser(
        prog="mtdb",
        description=DESCRIPTION,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("-v", "--version", action="version", version=get_version())
    ops = parser.add_mutually_exclusive_group()
    ops.add_argument(
        "-i", "--interface", metavar="DBPATH", help="link/initialize the MTDB at path"
    )
    ops.add_argument("-u", "--unlink", action="store_true", help="unlink from the MTDB")
    ops.add_argument(
        "-l",
        "--list",
        action="store_true",
        help="list historically linked primary MTDBs",
    )
    parser.add_argument(
        "target",
        nargs="?",
        metavar="SUBCOMMAND|OME",
        help="subcommand (see below) or ome code to look up",
    )
    parser.add_argument("rest", nargs=argparse.REMAINDER, help=argparse.SUPPRESS)
    return parser


def delegate(module_name, args):
    """Run a subcommand's CLI in-process and return its exit code.

    The subcommand modules live in this package (`mycotools.mtdb.<module>`) and
    parse `sys.argv` themselves, so their argv is swapped in for the call and the
    `SystemExit` they raise (argparse errors, explicit exits) is translated back
    to an exit code. Imports are deferred so a bare `mtdb` invocation stays light
    - `update` in particular pulls in the JGI/NCBI download stack."""
    module = importlib.import_module(f"mycotools.mtdb.{module_name}")
    saved_argv = sys.argv
    sys.argv = [f"mtdb {module_name}"] + list(args)
    try:
        module.cli()
        return 0
    except SystemExit as exc:
        if exc.code is None:
            return 0
        return exc.code if isinstance(exc.code, int) else 1
    finally:
        sys.argv = saved_argv


def link_mtdb(path):
    """Link (interface) to the MTDB rooted at `path`."""
    mtdb_loc = format_path(path)
    if not Path(mtdb_loc).is_dir() or not Path(mtdb_loc + "mtdb").is_dir():
        raise FileNotFoundError("invalid MTDB path")
    mtdb_initialize(mtdb_loc)


def list_links(config):
    """Print the historically linked primary MTDBs recorded in the user config."""
    for dbtype in config.get("log", {}):
        print(dbtype, flush=True)
        for mtdb_loc, login_time in config["log"][dbtype].items():
            print(f"\t{mtdb_loc} {login_time}", flush=True)
        print()


def lookup_omes(omes):
    """Print the database row, or a specific file path, for ome code(s)."""
    db = mtdb(primary_db()).set_index()
    for ome_prep in omes:
        if ome_prep in db:
            print(ome_prep + "\t" + "\t".join(str(v) for v in db[ome_prep].values()))
            return
        ome = re.sub(r"\.\w+[\w\d]$", "", ome_prep)
        ext_srch = re.search(r"^\d+\.?\d*\.(.*$)", ome_prep[6:])
        extension = ext_srch[1] if ext_srch is not None else None
        if ome in db:
            try:
                print(db[ome][extension] if extension else {"ome": ome, **db[ome]})
            except KeyError:
                raise KeyError("Invalid extension " + extension)
        else:
            for ref_ome, row in db.items():
                if ref_ome.startswith(ome + "."):
                    print(row[extension] if extension else {"ome": ome, **row})
                    break
            else:
                raise KeyError("Invalid ome " + ome)


def print_primary():
    """Print the primary MTDB path; return an exit code."""
    path = primary_db()
    if path:
        print(path, flush=True)
        return 0
    logger.error("Link a MycotoolsDB via `mtdb -i <MTDB_DIR>`")
    return 1


def main(argv=sys.argv):
    setup_logging()
    config = parse_user_config()
    args = build_parser().parse_args(argv[1:])

    # 1. forward recognized subcommands to their standalone tool
    if args.target in SUBCOMMANDS:
        sys.exit(delegate(SUBCOMMANDS[args.target], args.rest))
    # 2. terminal base operations
    if args.unlink:
        mtdb_disconnect()
        sys.exit(0)
    if args.list:
        list_links(config)
        sys.exit(0)
    # 3. ome-code lookup (skipped when linking, which reports the path afterward)
    if args.target and not args.interface:
        lookup_omes(args.target.replace('"', "").replace("'", "").split())
        sys.exit(0)
    # 4. optionally (re)link, then always report the primary MTDB path
    if args.interface:
        link_mtdb(args.interface)
    sys.exit(print_primary())


def cli():
    main(sys.argv)


if __name__ == "__main__":
    cli()
