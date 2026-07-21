#! /usr/bin/env python3

import logging
import os
import sys
import argparse
from datetime import datetime
from shutil import copy as cp
from mycotools.lib import mtdb_sql
from mycotools.lib.dbtools import primary_db, mtdb
from mycotools.lib.kontools import format_path, prep_output, setup_logging
from pathlib import Path

logger = logging.getLogger(__name__)


def soft_main(filetypes, db, output_path, print_link=False, verbose=False):
    """Symlink or print files from each file_type"""

    db = db.set_index("ome")
    # symlink files
    if not print_link:
        # make the directories for each requested file type
        for ftype in filetypes:
            if not Path(output_path + ftype).is_dir():
                Path(output_path + ftype).mkdir()
        # grab the files for each genome code
        for ome, row in db.items():
            for ftype in filetypes:
                if Path(row[ftype]).is_file():
                    sym_path = f"{output_path}{ftype}/{ome}.{ftype}"
                    try:
                        Path(sym_path).symlink_to(row[ftype])
                    except FileExistsError:
                        if Path(sym_path).is_symlink():
                            Path(sym_path).unlink()
                            Path(sym_path).symlink_to(row[ftype])
                        else:
                            logger.debug("" + ome + " " + ftype + " exists")
                else:
                    logger.debug("" + ome + " " + ftype)
    # simply print the link for each file
    else:
        for ome, row in db.items():
            for ftype in filetypes:
                print(row[ftype], flush=True)


def hard_main(filetypes, db, output_path):
    """Hard copy files from filetypes to their filetype output directory"""

    db = db.set_index("ome")
    # create the directories to output each file type
    for ftype in filetypes:
        if not Path(output_path + ftype).is_dir():
            Path(output_path + ftype).mkdir()

    # copy each file by genome
    for ome, row in db.items():
        for ftype in filetypes:
            try:
                cp(row[ftype], output_path + ftype + "/" + Path(row[ftype]).name)
            except FileNotFoundError:
                logger.error("" + ome + " " + ftype)


def mtdb_main(db, output_path, og_mtdb_path):
    """Create a MycotoolsDB directory with the files wanted for copy"""

    # generate the base directory for output
    if not output_path:
        output_path = str(Path.cwd()) + "/"
    if not Path(output_path).is_dir():
        Path(output_path).mkdir()
    mtdb_dir = output_path + "mycotoolsdb/"
    if not Path(mtdb_dir).is_dir():
        Path(mtdb_dir).mkdir()

    # generate the MTDB hierarchy subdirectories
    sub_dirs = [
        f"{mtdb_dir}log/",
        f"{mtdb_dir}config/",
        f"{mtdb_dir}mtdb/",
        f"{mtdb_dir}data/",
    ]
    for dir_ in sub_dirs:
        if not Path(dir_).is_dir():
            Path(dir_).mkdir()

    # copy the og_mtdb configuration
    cp(og_mtdb_path + "config/mtdb.json", f"{mtdb_dir}config/mtdb.json")

    # output the database: the generated hierarchy is meant to be linked with
    # `mtdb -i`, so its primary uses the SQLite backend, with a dated `.mtdb`
    # snapshot alongside it for portability
    cdate = datetime.now().strftime("%Y%m%d")
    db.to_sql(f"{mtdb_dir}mtdb/{mtdb_sql.PRIMARY_DB_NAME}")
    db.df2db(f"{mtdb_dir}log/{cdate}.mtdb", headers=True)

    # output the files
    hard_main(["gff3", "faa", "fna"], db, f"{mtdb_dir}data/")


def cli():

    parser = argparse.ArgumentParser(
        description="Symlinks/copies selected files from database"
    )
    parser.add_argument(
        "-d", "--mtdb", default=primary_db(), help="DEFAULT: primary_db"
    )
    parser.add_argument("-a", "--assembly", action="store_true", help="Grab assemblies")
    parser.add_argument("-p", "--proteome", action="store_true", help="Grab proteomes")
    parser.add_argument("-g", "--gff", action="store_true", help="Grab gff`s")
    parser.add_argument("--print", action="store_true", help="Print paths, no copy")
    parser.add_argument("--hard", action="store_true", help="Hard copy files")
    parser.add_argument(
        "-n", "--new_mtdb", action="store_true", help="Create MTDB directory hierarchy"
    )
    parser.add_argument("-o", "--output", default=str(Path.cwd()))
    args = parser.parse_args()
    setup_logging(verbose=getattr(args, "verbose", False))

    if not args.assembly and not args.proteome and not args.gff and not args.new_mtdb:
        logger.error("--assembly/--proteome/--gff/--new_mtdb required")
        sys.exit(4)
    if args.new_mtdb:
        args.hard = False
        args.print = False

    db_path = format_path(args.mtdb)
    args.output = format_path(args.output, force_dir=True)
    output_path = prep_output(args.output, cd=False)

    filetypes = []
    if args.proteome:
        filetypes.append("faa")
    if args.gff:
        filetypes.append("gff3")
    if args.assembly:
        filetypes.append("fna")

    db = mtdb(db_path)
    if args.new_mtdb:
        mtdb_main(db, output_path, format_path(os.environ["MYCODB"] + "/../"))
    elif args.print or not args.hard:
        soft_main(filetypes, db, output_path, print_link=args.print)
    else:
        hard_main(filetypes, db, output_path)

    sys.exit(0)


if __name__ == "__main__":
    cli()
