#! /usr/bin/env python3

import logging
import os
import re
import sys
import copy
import datetime
import argparse
import pandas as pd
import numpy as np
from mycotools.lib.kontools import intro, outro, setup_logging
from mycotools.lib.dbtools import db2df, df2db, read_log, log_editor
from mycotools.download.jgi import main as jgi_dwnld
from pathlib import Path

logger = logging.getLogger(__name__)


def compile_log(log_path):

    log = {}
    if not Path(log_path).is_file():
        with open(log_path, "w") as out:
            out.write("#assembly_acc\tfna\tgff3\tfaa")
    else:
        log = read_log(log_path)

    return log


def jgi_redundancy_check(db, jgi_df, duplicates={}, ome_col="portal", jgi2ncbi={}):
    '''everything in jgi_df (pd.DataFrame()) should be have a valid biosample
    at this point.
    db will have the index set at "assembly_acc"'''

    db["version"].fillna(0.0)
    jgi_df["version"] = jgi_df.loc[:, "version"].replace(np.nan, "0.0")
    preOmes, biosamples = set(copy.deepcopy(db.index)), set(db["biosample"])
    updates, old_omes, todel_db, todel_jgi = {}, {}, [], []

    for i, row in jgi_df.iterrows():
        #        version = float(row['version'].replace('v',''))
        version = float(row["version"])
        ome, biosample = (
            copy.deepcopy(row[ome_col]),
            row["biosample"],
        )  # implied to exist
        if ome in duplicates:
            todel_jgi.append(i)
        elif ome in preOmes:
            db.at[ome, "published"] = copy.deepcopy(row["publication(s)"])
            if pd.isnull(db["version"][ome]):
                db.at[ome, "version"] = 0
            # update most current pub status
            if not db["version"][ome] or pd.isnull(db["version"][ome]):
                db.at[ome, "version"] = 0
                db_version = 0
            else:
                try:
                    db_version = float(db["version"][ome])
                except (ValueError, AttributeError):
                    db_version = float(db["version"][ome].replace("v", ""))
                    db.at[ome, "version"] = db_version
            if version > db_version and db["source"][ome] == "jgi":
                db_organism = (
                    db.loc[ome, "genus"]
                    + "_"
                    + db.loc[ome, "species"]
                    + "_"
                    + db.loc[ome, "strain"]
                )
                organism = jgi_df["name"].replace(" ", "_")
                updates[ome] = [db.loc[ome, "ome"], db_organism, organism]
                old_omes[db.loc[ome, "ome"]] = copy.deepcopy(db.loc[ome])
                # retain row information
                #                jgi_df.at[i, 'ome'] = db.at[ome, 'ome'] # keep the same ome code.
                # I think this ought to be updated to a new code
                todel_db.append(ome)
            else:
                todel_jgi.append(i)
        elif biosample in biosamples and biosample.rstrip():
            # genome_code isn't here, but the biosample is
            check_db = db[db["biosample"] == biosample]
            checkOme = list(check_db["ome"])[0]
            if len(check_db) > 1:  # more than one entry to the biosample, keep as it is
                todel_jgi.append(i)
            else:  # this is an NCBI accession that can be updated to JGI
                index = list(check_db.index)[0]
                db_organism = (
                    db.loc[index, "genus"]
                    + "_"
                    + db.loc[index, "species"]
                    + "_"
                    + db.loc[index, "strain"]
                )
                organism = jgi_df["name"].replace(" ", "_")
                updates[ome] = [checkOme, db_organism, organism]
                todel_db.append(index)
        elif ome.lower() in jgi2ncbi:  # NCBI can be updated to JGI
            low_ome = ome.lower()
            ncbi_acc = jgi2ncbi[low_ome]
            if ncbi_acc in set(db["assembly_acc"]):
                check_db = db[db["assembly_acc"] == ncbi_acc]
                checkOme = list(check_db["ome"])[0]
                index = list(check_db.index)[0]
                db_organism = (
                    db.loc[index, "genus"]
                    + "_"
                    + db.loc[index, "species"]
                    + "_"
                    + db.loc[index, "strain"]
                )
                organism = jgi_df["name"].replace(" ", "_")
                updates[ome] = [checkOme, db_organism, organism]
                todel_db.append(index)
            else:
                organism = jgi_df["name"].replace(" ", "_")
                updates[ome] = [None, None, organism]
        else:
            organism = jgi_df["name"].replace(" ", "_")
            updates[ome] = [None, None, organism]

    todel_db.sort(reverse=True)
    todel_jgi.sort(reverse=True)
    for i in todel_db:
        db = db.drop(i)
    for i in todel_jgi:
        jgi_df = jgi_df.drop(i)

    db = db.reset_index()

    return jgi_df, db, updates, old_omes


def logged_file(log, ome, typ, output):
    """Return the path of an ome's `typ` file recorded in a previous run's log,
    if that file is still on disk; otherwise None. JGI files arrive gzipped and
    curation decompresses them in place, so either form is accepted."""
    basename = log.get(ome, {}).get(typ, "na")
    if basename in {"na", "error", "pending", ""}:
        return None
    path = f"{output}/{typ}/{basename}"
    for candidate in (path, re.sub(r"\.gz$", "", path)):
        if Path(candidate).is_file():
            return candidate
    return None


def runjgi_dwnld(
    jgi_df,
    user,
    pwd,
    ome_col,
    output,
    log,
    log_path,
    dwnlds,
    failed,
    rerun,
    masked,
    spacer="\t\t",
    restore_wait=None,
):
    """Download the MycoCosm portals in `jgi_df` via the JGI Data Portal API,
    skipping omes a previous run already completed and recording each outcome in
    the resume log. Returns (jgi_df, log, failed, deferred).

    Tape-archived genomes are skipped on the first pass and revisited once every
    portal has been visited - nearly every MycoCosm genome needs a restore, so
    waiting on each in turn would stall the update. `restore_wait` caps how many
    minutes any one genome is waited on there; None waits for as long as JGI
    takes.

    Omes without an assembly or gff3 are dropped from `jgi_df`. Genuine failures
    (portal absent from MycoCosm, no such file type, corrupt download) are logged
    as `error` and appended to `failed` as [ome, version], as the legacy per-ome
    loop did. Omes whose files are merely awaiting a JGI tape restore are instead
    logged as `pending` and returned in `deferred` - they are retried on the next
    run rather than blacklisted."""

    # resume: an ome needs no download when every requested file is logged and
    # present - or, without `rerun`, was previously logged as unobtainable
    todwnld_i, preexisting = [], {}
    for i, row in jgi_df.iterrows():
        ome = row[ome_col]
        paths = {typ: logged_file(log, ome, typ, output) for typ in dwnlds}
        settled = [
            typ
            for typ, path in paths.items()
            if path or (not rerun and log.get(ome, {}).get(typ) == "error")
        ]
        if len(settled) == len(dwnlds):
            preexisting[i] = {t: p for t, p in paths.items() if p}
        else:
            todwnld_i.append(i)

    if preexisting:
        logger.debug(f"{spacer}{len(preexisting)} preexisting download(s)")
    for i, paths in preexisting.items():
        for typ, path in paths.items():
            jgi_df.at[i, typ + "_path"] = path

    api_deferred = set()
    if todwnld_i:
        post_df, api_failed = jgi_dwnld(
            jgi_df.loc[todwnld_i].copy(),
            output,
            user,
            pwd,
            assembly="fna" in dwnlds,
            proteome="faa" in dwnlds,
            gff3="gff3" in dwnlds,
            masked=masked,
            spacer=spacer,
            ome_col=ome_col,
            deferred=api_deferred,
            restore_wait=restore_wait,
            defer_tape=True,
        )
    else:
        post_df, api_failed = None, set()

    todel, deferred = [], []
    if post_df is not None:
        for i, row in post_df.iterrows():
            ome = row[ome_col]
            if ome not in log:
                log[ome] = {"fna": "na", "gff3": "na", "faa": "na"}
            # a file JGI has yet to stage to disk is pending, not failed
            unobtained = "pending" if ome in api_deferred else "error"
            for typ in dwnlds:
                path = row.get(typ + "_path")
                if isinstance(path, str) and path:
                    jgi_df.at[i, typ + "_path"] = path
                    log[ome][typ] = Path(path).name
                    logger.debug(f"{spacer}{ome} {typ}: {Path(path).name}")
                else:
                    log[ome][typ] = unobtained
                    logger.debug(f"{spacer}{ome} {typ}: {unobtained.upper()}")
            # JGI metadata may fill in curation fields the MycoCosm table lacks
            for col in ("genus", "species", "strain"):
                if col in post_df.columns:
                    jgi_df.at[i, col] = row[col]
            log_editor(
                log_path,
                ome,
                ome
                + "\t"
                + log[ome]["fna"]
                + "\t"
                + log[ome]["gff3"]
                + "\t"
                + log[ome]["faa"],
            )
            if ome in api_deferred:
                deferred.append([ome, jgi_df["version"][i]])
                todel.append(i)
            elif ome in api_failed:
                failed.append([ome, jgi_df["version"][i]])
                todel.append(i)

    for i in todel:
        jgi_df = jgi_df.drop(i)

    if deferred:
        logger.info(
            f"{spacer}{len(deferred)} genome(s) awaiting JGI tape restore; "
            "they will be retried on the next run"
        )

    return jgi_df, log, failed, deferred


def main(
    jgi_df,
    ref_db,
    output,
    user,
    pwd,
    date=datetime.datetime.now().strftime("%Y%m%d"),
    assembly=True,
    proteome=False,
    gff3=True,
    update=True,
    repeatmasked=True,
    nonpublished=False,
    rerun=False,
    failed_dict={},
    duplicates={},
    spacer="\t\t",
    jgi2ncbi={},
    restore_wait=None,
):

    if not nonpublished:
        jgi_df = jgi_df[jgi_df["is published"] == "Y"]
    if "assembly_acc" in jgi_df.columns:
        ome_col = "assembly_acc"
    elif "portal" in jgi_df.columns:
        ome_col = "portal"
    else:
        logger.debug(spacer + "invalid MycoCosm tsv headers")
        sys.exit(3)

    toDel = []

    if not rerun:
        jgi_df = jgi_df.set_index(ome_col)
        jgi_omes = set(jgi_df.index)
        for failure in failed_dict:
            if failure in jgi_omes:
                if "v" in str(jgi_df["version"][failure]):
                    version = float(jgi_df["version"][failure].replace("v", ""))
                    if failed_dict[failure]["version"] != "":
                        old_vers = float(
                            failed_dict[failure]["version"].replace("v", "")
                        )
                        if not version > old_vers:
                            toDel.append(failure)
                else:
                    toDel.append(failure)
        for failed in toDel:
            jgi_df = jgi_df.drop(failed)
        jgi_df = jgi_df.reset_index()

    logger.debug(spacer + "Redundancy check")
    if isinstance(ref_db, pd.DataFrame):
        ref_db["index"] = ref_db["assembly_acc"].copy()
        ref_db = ref_db.set_index("index")
        jgi_df, new_ref_db, updates, old_rows = jgi_redundancy_check(
            ref_db, jgi_df, ome_col=ome_col, jgi2ncbi=jgi2ncbi
        )
        output_str = "ome\tdb_organism\tjgi_organism\tassembly_acc\n"
        if isinstance(updates, dict):
            for ome, update in updates.items():
                output_str += "\t".join([str(x) for x in update]) + "\t" + ome + "\n"
            with open(output + "/jgiUpdates.tsv", "w") as out:
                out.write(output_str)
        else:
            updates.to_csv(f"{output}/jgiUpdates.tsv", sep="\t")
        #        update_check = {i[-1]: i[0:3] for i in updates if i[0]}
        logger.debug(spacer + "" + str(len(jgi_df)) + " genomes to assimilate")
    else:
        new_ref_db, updates = None, {}

    failed = []
    log_path = output + "/jgi2db.log"
    log = compile_log(log_path)
    if not rerun:
        prev_omes = set(jgi_df[ome_col])
        for ome in log:
            if ome not in prev_omes:
                continue
            if log[ome]["fna"] == "error" or log[ome]["gff3"] == "error":
                drop_index = list(jgi_df[jgi_df[ome_col] == ome].index)[0]
                failed.append([ome, jgi_df["version"][drop_index]])
                jgi_df = jgi_df.drop(drop_index)

    logger.debug(spacer + "Downloading JGI data")
    dwnlds = []
    if assembly:
        dwnlds.append("fna")
    if gff3:
        dwnlds.append("gff3")
    if proteome:
        dwnlds.append("faa")

    for typ in dwnlds:
        if not Path(output + "/" + typ).is_dir():
            Path(output + "/" + typ).mkdir()

    jgi_df, log, failed, deferred = runjgi_dwnld(
        jgi_df,
        user,
        pwd,
        ome_col,
        output,
        log,
        log_path,
        dwnlds,
        failed,
        rerun,
        repeatmasked,
        spacer,
        restore_wait=restore_wait,
    )

    jgi_df = jgi_df.rename(
        columns={
            "published(s)": "published",
            ome_col: "assembly_acc",
            "fna_path": "assemblyPath",
            "gff3_path": "gffPath",
        }
    )
    jgi_df["source"] = "jgi"

    # add back entries whose attempted update did not complete, whether it
    # failed outright or is still awaiting a JGI tape restore
    for ome_d in failed + deferred:
        ome = ome_d[0]
        if ome in old_rows:  # if there is an old row to add back
            new_ref_db = new_ref_db.append(old_rows[ome])

    old_assembly_accs = set(new_ref_db["assembly_acc"])
    jgi_df = jgi_df.set_index("assembly_acc")
    new_ref_db = new_ref_db.set_index("assembly_acc")
    for ome in updates:  # apply the old ome code to update the version
        if ome in old_assembly_accs:
            jgi_df.at[ome, "ome"] = new_ref_db.loc[ome, "ome"]

    jgi_df = jgi_df.astype(str)
    jgi_premtdb_df = jgi_df.reset_index()
    for i, row in jgi_premtdb_df.iterrows():
        if not pd.isnull(row["publication(s)"]) and row["publication(s)"]:
            jgi_premtdb_df.at[i, "published"] = row["publication(s)"]
        elif not pd.isnull(row["is published"]) and row["is published"]:
            jgi_premtdb_df.at[i, "published"] = 1
        elif not pd.isnull(row["is public"]) and row["is public"]:
            jgi_premtdb_df.at[i, "published"] = 1

    return jgi_premtdb_df, new_ref_db.reset_index(), failed, deferred


def cli():

    mycoCosmURL = (
        "https://mycocosm.jgi.doe.gov/ext-api/mycocosm/catalog/"
        + "download-group?flt=&seq=all&pub=all&grp=fungi&srt=released&ord=desc"
    )
    parser = argparse.ArgumentParser(
        description=" Downloads or imports MycoCosm table."
        + " Downloads files and outputs a mycotools `.db` of novel JGI downloads relative to"
        + " a reference `.db`. Can also begin a new `.db`"
    )
    parser.add_argument(
        "-l",
        "--login",
        help=r'Login file: "<JGI username>\t<JGI Password>\n<NCBI API key>"',
    )
    parser.add_argument("-d", "--database", help="Existing myctools `.db` to reference")
    parser.add_argument(
        "-m",
        "--mycocosm",
        default=mycoCosmURL,
        help="MycoCosm master table path / URL. This may change from time to time "
        + "DEFAULT: "
        + mycoCosmURL,
    )
    parser.add_argument(
        "-u",
        "--update",
        help="Download potential JGI updates for existing omes "
        "- NOT FUNCTIONAL FOR DATABASE SAFETY",
        default=False,
        action="store_true",
    )
    parser.add_argument(
        "-o", "--output", help="Output directory. DEFAULT: Current date"
    )
    parser.add_argument(
        "-a",
        "--assembly",
        default=True,
        action="store_false",
        help="Do not download assemblies.",
    )
    parser.add_argument(
        "-p",
        "--proteome",
        default=True,
        action="store_false",
        help="Do not download proteomes.",
    )
    parser.add_argument(
        "-g", "--gff", default=True, action="store_false", help="Do not download gff3s."
    )
    parser.add_argument(
        "-r",
        "--repeatmasked",
        default=True,
        action="store_false",
        help="Download nonmasked assemblies.",
    )

    args = parser.parse_args()
    setup_logging(verbose=getattr(args, "verbose", False))

    args_dict = {
        "Preexisting db": args.database,
        "MycoCosm Table": args.mycocosm,
        "Assemblies": args.assembly,
        "RepeatMasked": args.repeatmasked,
        "Proteomes": args.proteome,
        ".gff3's": args.gff,
    }

    start_time = intro("jgi2db", args_dict)
    if args.output:
        output = os.path.abspath(args.output)
    else:
        output = start_time.strftime("%Y%m%d") + "_jgi2db"
    if not Path(output).is_dir():
        Path(output).mkdir()

    if args.login:
        with open(args.login, "r") as raw:
            prep = raw.read()
        data = [x.split("\t") for x in prep.split("\n")]
        apikey = None
        if len(data) > 1 and data[1]:
            # NCBI API key is the last field of the second line; a legacy
            # leading NCBI email column (now unused) is tolerated
            api_field = data[1][-1]
            if api_field != "":
                apikey = api_field

    ref_db = db2df(format_path(args.database))
    jgi_df = main(args.mycocosm, refdb, output)

    df2db(jgi_df, output + "/new.db")
    logger.debug(
        "Success! "
        + str(len(jgi_df))
        + " added to database\n \
            Run updateDB to confirm and finish update."
    )

    outro(start_time)


if __name__ == "__main__":
    cli()
