#! /usr/bin/env python3

import os
import sys
import logging
import argparse
from mycotools.lib.dbtools import loginCheck, primaryDB, mtdb, encrypt_pw
from mycotools.lib.kontools import format_path, read_json, collect_files, setup_logging
from pathlib import Path

logger = logging.getLogger(__name__)


def rm_outdated(omes, yes=False):
    """Remove outdated genomes after compiling them"""

    biofiles, to_del = [], []
    # compile the files
    biofiles.extend(
        [f"{os.environ['MYCOGFF3']}/{x}" for x in [p.name for p in Path(os.environ["MYCOGFF3"]).iterdir()]]
    )
    biofiles.extend(
        [f"{os.environ['MYCOFAA']}/{x}" for x in [p.name for p in Path(os.environ["MYCOFAA"]).iterdir()]]
    )
    biofiles.extend(
        [f"{os.environ['MYCOFNA']}/{x}" for x in [p.name for p in Path(os.environ["MYCOFNA"]).iterdir()]]
    )

    # remove each biofile
    for i in biofiles:
        ome = None
        ome_prep = Path(i).name
        if ome_prep.endswith(".gff3"):
            ome = ome_prep[:-5]
        elif ome_prep.endswith(".faa"):
            ome = ome_prep[:-4]
        elif ome_prep.endswith(".fna"):
            ome = ome_prep[:-4]
        else:  # safer to preserve independent placements
            continue
        if ome not in omes:
            to_del.append(i)

    if to_del:
        if yes:
            data = "y"
        else:
            data = input(f"\n{len(to_del)} omes to be deleted.\n" + "Continue [y/N]? ")
        if data.lower() in {"yes", "y"}:
            for i in to_del:
                Path(i).unlink()
        else:
            raise KeyError("cache removal stopped")


def restrictions(
    db, restr_list, mtdb_config=format_path("~/.mycotools/config.json"), yes=False
):
    mtdb_config = read_json(mtdb_config)
    restr_path = mtdb_config[mtdb_config["active"]]["MYCODB"] + "../log/failed.tsv"

    try:
        with open(restr_path, "r") as raw:
            restricted = [x.rstrip().split() for x in raw]
    except FileNotFoundError:
        restricted = []

    accs = set(x[0] for x in restricted)
    for r, s, reason in restr_list:
        if s.lower() in {"ncbi", "jgi"} and r not in accs:
            restricted.append([r, s.lower(), str(reason)])
            logger.info("%s %s", r, s)

    in_db = [x[0] for x in restricted if x[0] in db]
    while in_db:
        if not yes:
            check = input("Some restrictions are in the MTDB. Delete them? [y/N]: ")
            if check.lower() in {"yes", "y"}:
                break
            else:
                sys.exit(1)

    with open(restr_path, "w") as out:
        out.write("\n".join(["\t".join(x) for x in restricted]))


def cli():
    parser = argparse.ArgumentParser(
        description="Primary MycotoolsDB management utility"
    )
    parser.add_argument(
        "-c", "--clear_cache", action="store_true", help="Clear MycotoolsDB legacy data"
    )
    parser.add_argument(
        "-p",
        "--password",
        action="store_true",
        help="Encrypt NCBI/JGI passwords to expedite access",
    )
    parser.add_argument(
        "-r",
        "--restrict",
        help="Restrict assembly accessions file, formatted: "
        + "<ACCESSION>\t<SOURCE>\t[REASON]",
    )
    parser.add_argument("-y", "--yes", help="Answer yes", action="store_true")
    args = parser.parse_args()
    setup_logging(verbose=getattr(args, "verbose", False))

    db = mtdb(primaryDB()).set_index("assembly_acc")

    if args.password:
        ncbi_email, ncbi_api, jgi_email, jgi_pwd = loginCheck()
        encrypt_pw(ncbi_email, ncbi_api, jgi_email, jgi_pwd)
    if args.restrict:
        restrict_path = format_path(args.restrict)
        with open(restrict_path, "r") as raw:
            restricted = [x.rstrip().split("\t") for x in raw]
        for v in restricted:
            if len(v) < 3:
                v = v + [None]
        restrictions(db, restricted, yes=args.yes)
    if args.clear_cache:
        rm_outdated(mtdb(primaryDB())["ome"], args.yes)

    sys.exit(0)


if __name__ == "__main__":
    cli()
