#! /usr/bin/env python3

import os
import sys
import logging
import argparse
from mycotools.lib.dbtools import (
    login_check,
    primary_db,
    mtdb,
    encrypt_pw,
    get_login,
    store_login,
)
from mycotools.lib import mtdb_sql
from mycotools.lib.kontools import format_path, read_json, setup_logging
from pathlib import Path

logger = logging.getLogger(__name__)


def ome_list(db_path):
    """Every ome in the database, without materializing the rest of it."""
    if mtdb_sql.is_sqlite(db_path):
        return mtdb_sql.omes(db_path)
    return list(mtdb(db_path)["ome"])


def rm_outdated(omes, yes=False):
    """Remove outdated genomes after compiling them"""

    biofiles, to_del = [], []
    # compile the files
    biofiles.extend(
        [
            f"{os.environ['MYCOGFF3']}/{x}"
            for x in [p.name for p in Path(os.environ["MYCOGFF3"]).iterdir()]
        ]
    )
    biofiles.extend(
        [
            f"{os.environ['MYCOFAA']}/{x}"
            for x in [p.name for p in Path(os.environ["MYCOFAA"]).iterdir()]
        ]
    )
    biofiles.extend(
        [
            f"{os.environ['MYCOFNA']}/{x}"
            for x in [p.name for p in Path(os.environ["MYCOFNA"]).iterdir()]
        ]
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


def migrate(yes=False):
    """Convert a flat-file primary MTDB into the SQLite backend.

    The `.mtdb` it was built from is left in place -- `primary_db()` prefers
    `mtdb.db` once it exists, so the flat file becomes an inert snapshot that
    can be deleted, kept for provenance, or handed to an older Mycotools."""
    db_path = primary_db()
    if not db_path:
        logger.error("Link a MycotoolsDB via `mtdb -i <MTDB_DIR>`")
        return 1
    if mtdb_sql.is_sqlite(db_path):
        logger.info("Primary MTDB is already SQLite: %s", db_path)
        return 0

    target = format_path("$MYCODB/" + mtdb_sql.PRIMARY_DB_NAME)
    db = mtdb(db_path)
    n = len(db["ome"])
    if not yes:
        check = input(f"Convert {n} genomes in {db_path} to {target}? [y/N]: ")
        if check.lower() not in {"yes", "y"}:
            return 1

    db.to_sql(target)
    migrated = mtdb_sql.count(target)
    if migrated != n:
        logger.error("migrated %d of %d genomes; %s left in place", migrated, n, db_path)
        return 1
    logger.info("Migrated %d genomes -> %s", migrated, target)
    logger.info("%s is now a snapshot and is no longer read", db_path)
    return 0


def cli():
    parser = argparse.ArgumentParser(
        description="Primary MycotoolsDB management utility"
    )
    parser.add_argument(
        "-m",
        "--migrate",
        action="store_true",
        help="Convert a flat-file primary MTDB to the SQLite backend",
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
        "-s",
        "--store",
        action="store_true",
        help="Store NCBI/JGI credentials WITHOUT a password (unencrypted, chmod 600)",
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

    if args.migrate:
        sys.exit(migrate(args.yes))

    if args.password and args.store:
        logger.error("--password and --store are mutually exclusive")
        sys.exit(1)
    if args.password:
        ncbi_api, jgi_email, jgi_pwd = login_check()
        encrypt_pw(ncbi_api, jgi_email, jgi_pwd)
    if args.store:
        ncbi_api, jgi_email, jgi_pwd = get_login(ncbi=True, jgi=True)
        store_login(ncbi_api, jgi_email, jgi_pwd)
    if args.restrict:
        restrict_path = format_path(args.restrict)
        with open(restrict_path, "r") as raw:
            restricted = [x.rstrip().split("\t") for x in raw]
        for v in restricted:
            if len(v) < 3:
                v = v + [None]
        # loaded here rather than up front so the credential and migration
        # operations do not pay for reading the whole database
        db = mtdb(primary_db()).set_index("assembly_acc")
        restrictions(db, restricted, yes=args.yes)
    if args.clear_cache:
        rm_outdated(set(ome_list(primary_db())), args.yes)

    sys.exit(0)


if __name__ == "__main__":
    cli()
