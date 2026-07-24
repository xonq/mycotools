#! /usr/bin/env python3

# NEED multiple lineages from command line
# NEED stdin acceptance for most of these arguments

import os
import sys
import logging
import argparse
from mycotools.lib.kontools import (
    file2list,
    format_path,
    setup_logging,
    mk_output,
)
from mycotools.lib import mtdb_sql
from mycotools.lib.dbtools import mtdb, primary_db, db_stem
from mycotools.mtdb.files import mtdb_main as gen_full_mtdb
from pathlib import Path

logger = logging.getLogger(__name__)


def load_db(db_path, omes_set=(), aa_set=(), lineage_list=()):
    """Load only the rows an extraction can possibly need.

    A SQLite primary database can answer "which genomes" before anything is
    materialized, so an ome list, an assembly-accession list, or a set of
    lineages narrows the read to an index seek. Anything else -- and any
    `.mtdb` flat file -- falls back to reading the whole database."""
    if not mtdb_sql.is_sqlite(db_path):
        return mtdb(db_path)
    if omes_set:
        return mtdb(mtdb_sql.select_omes(db_path, omes_set))
    if aa_set:
        return mtdb(mtdb_sql.select_column(db_path, "assembly_acc", aa_set))
    if lineage_list:
        # resolve lineages against the normalized taxonomy table; a
        # species/strain lineage is not answerable there, so read everything
        ranks = {mtdb_sql.infer_rank(db_path, lin) for lin in lineage_list}
        if not ranks.intersection({None, "species", "strain"}):
            genera = mtdb_sql.genera_for_lineages(db_path, lineage_list)
            return mtdb(mtdb_sql.select_column(db_path, "genus", genera))
    return mtdb(db_path)


def main(
    db,
    rank=None,
    x_number=0,
    lineage_list=[],
    omes_set=set(),
    by_rank=False,
    source=None,
    nonpublished=False,
    inverse=False,
    aa_set=set(),
    seed=None,
):
    """Python entry point for extract_mtdb"""

    db = db.set_index("ome")
    if x_number > 0:
        db = db.extract_unique(x_number, rank=rank, seed=seed)

    # extract each taxonomic entry based on the classification specified
    if lineage_list:
        new_db = db.extract_tax(lineage_list)
    # if an ome list is specified then open it, store each entry in a list and pull each ome
    elif omes_set:
        new_db = db.extract_ome(omes_set)
    elif aa_set:
        new_db = db.extract_ome(aa_set, "assembly_acc")
    # if none of these are specified then create a `new_db` variable to work for later
    else:
        new_db = db

    # if there is a source specified, extract it or the opposite
    if source:
        new_db = new_db.extract_source(source)

    # if you want publisheds, then just pull those out
    if not nonpublished:
        new_db = new_db.extract_pub()

    if inverse:
        new_omes = set(new_db.keys())
        inv_db = mtdb().set_index()
        for ome, row in db.items():
            if ome not in new_omes:
                inv_db[ome] = row
        new_db = inv_db

    if by_rank:
        lineages = set(v["taxonomy"][rank] for k, v in new_db.items())
        dbs = {}
        for lineage in lineages:
            if lineage:
                dbs[lineage.lower()] = new_db.extract_tax([lineage]).reset_index()
            else:
                dbs["unclassified"] = new_db.extract_tax([""]).reset_index()
        return dbs
    else:
        return new_db.reset_index()


def cli():
    ranks = [
        "kingdom",
        "phylum",
        "subphylum",
        "class",
        "order",
        "family",
        "genus",
        "species",
        "strain",
    ]
    parser = argparse.ArgumentParser(
        description="Extracts a MycotoolsDB from arguments. E.g.\t`mtdb extract "
        + "-l Atheliaceae`"
    )

    ex_opt = parser.add_argument_group("Extraction parameters")
    ex_opt.add_argument("-l", "--lineage", default="")
    ex_opt.add_argument("-s", "--source")
    ex_opt.add_argument(
        "-n", "--nonpublished", action="store_true", help="Include restricted"
    )
    ex_opt.add_argument("-r", "--rank", help=f"[-a|-b] Taxonomic rank: {ranks}")
    ex_opt.add_argument(
        "-a",
        "--allowed_rank",
        default=0,
        type=int,
        help="[-r] Number of randomly sampled --rank allowed",
    )
    ex_opt.add_argument(
        "-b",
        "--by_rank",
        action="store_true",
        help="[-r] Output MTDBs for each lineage in --rank",
    )
    ex_opt.add_argument(
        "-i",
        "--inverse",
        action="store_true",
        help="Inverse [source|lineage(s)|nonpublished]",
    )
    ex_opt.add_argument(
        "--seed",
        type=int,
        help="[-a] Seed the random sample for a reproducible selection",
    )
    ex_opt.add_argument("-ol", "--ome", help="File w/list of omes")
    ex_opt.add_argument(
        "-al", "--assembly_list", help="File w/list of assembly accessions"
    )
    ex_opt.add_argument("-ll", "--lineages", help="File w/list of lineages")

    out_opt = parser.add_argument_group("Output parameters")
    out_opt.add_argument(
        "-m", "--new_mtdb", action="store_true", help="Create MTDB directory hierarchy"
    )
    out_opt.add_argument("-p", "--paths", help="Output with paths", action="store_true")
    out_opt.add_argument("--headers", action="store_true")
    out_opt.add_argument("-d", "--mtdb", help="- for stdin", default=primary_db())
    out_opt.add_argument("-o", "--output")

    args = parser.parse_args()
    setup_logging(verbose=getattr(args, "verbose", False))
    db_path = format_path(args.mtdb)

    if args.lineage or args.lineages:
        logger.warning(
            "extracting taxonomy is subject to " + "errors in NCBI's hierarchy"
        )

    # these arguments require one another
    if args.by_rank and not args.rank:
        logger.error("--by_rank requires --rank")
        sys.exit(10)
    elif args.allowed_rank and not args.rank:
        logger.error("--allowed_rank requres --rank")
        sys.exit(11)
    elif args.rank and not args.allowed_rank and not args.by_rank:
        logger.error("--rank requires --allowed_rank or --by_rank")
    elif args.rank:
        args.rank = args.rank.lower()
        if args.rank not in set(ranks):
            logger.error(f"--rank not in {ranks}")
            sys.exit(12)

    args.lineage = args.lineage.replace('"', "").replace("'", "")

    output = ""
    if args.output:
        # `-o` names the directory to write into; the file inside it is named
        # for the source database and the filters applied. The extension is
        # always `.mtdb` -- an extract is an interchange file regardless of
        # which backend it was read from.
        out_dir = format_path(args.output, force_dir=True)
        Path(out_dir).mkdir(parents=True, exist_ok=True)
        tag = ""
        if args.lineage:
            tag += "_" + args.lineage
        if args.lineages:
            tag += "_taxonomy"
        if args.source:
            tag += "_" + args.source.lower()
        if not args.nonpublished:
            tag += "_pub"
        output = f"{out_dir}{db_stem(db_path)}{tag}.mtdb"

    if args.ome:
        omes = set(file2list(format_path(args.ome)))
    else:
        omes = set()

    if args.assembly_list:
        aa_set = set(file2list(format_path(args.assembly_list)))
    else:
        aa_set = set()

    if args.lineages:
        lineage_list = file2list(format_path(args.lineages))
    elif args.lineage:
        lineage_list = [args.lineage]
    else:
        lineage_list = []

    if args.mtdb == "-":
        db = mtdb.from_string(sys.stdin.read())
    elif args.inverse:
        # the inverse needs every row to subtract from, so no narrowing
        db = mtdb(db_path)
    else:
        db = load_db(db_path, omes, aa_set, lineage_list)

    new_db = main(
        db,
        lineage_list=lineage_list,
        omes_set=omes,
        source=args.source,
        rank=args.rank,
        x_number=args.allowed_rank,
        by_rank=args.by_rank,
        nonpublished=args.nonpublished,
        inverse=args.inverse,
        aa_set=aa_set,
        seed=args.seed,
    )
    if args.new_mtdb:
        gen_full_mtdb(
            new_db, format_path(output), format_path(os.environ["MYCODB"] + "/../")
        )
    elif args.output or args.by_rank:
        if isinstance(new_db, mtdb):
            new_db.df2db(output, headers=bool(args.headers), paths=args.paths)
        else:
            out_dir = mk_output(output or str(Path.cwd()), "extract_mtdb")
            prefix = db_stem(db_path)
            for lineage, db in new_db.items():
                out_f = f"{out_dir}{prefix}.{lineage}.mtdb"
                db.df2db(out_f, headers=bool(args.headers))
    else:
        new_db.df2db(headers=bool(args.headers), paths=args.paths)

    sys.exit(0)


if __name__ == "__main__":
    cli()
