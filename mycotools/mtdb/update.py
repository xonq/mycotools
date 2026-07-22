#! /usr/bin/env python3

# NEED reinit implementation
# NEED to update introduction
# NEED a verbose option
# NEED revert version option (ome-by-ome/list of omes)
# NEED to reference a manually curated duplicate check
# NEED a prohibit option to import prohibited JGI/NCBI IDs and option to update
# NEED to finish --save implementation
# NEED an only one of a strain feature in config
# NEED to pull failed JGI from NCBI
# will require logging whatever NCBI omes directly overlap MycoCosm
# NEED to remove overlap when rerunning failed genomes

import logging
import os
import re
import sys
import json
import shutil
import zipfile
import requests
import argparse
import subprocess
import numpy as np
import pandas as pd
from tqdm import tqdm
from Bio import Entrez
from datetime import datetime
from collections import defaultdict
from mycotools.lib import mtdb_sql
from mycotools.lib.dbtools import (
    db2df,
    df2db,
    gather_taxonomy,
    assimilate_tax,
    primary_db,
    login_check,
    log_editor,
    mtdb,
    mtdb_initialize,
)
from mycotools.lib.kontools import (
    intro,
    outro,
    format_path,
    prep_output,
    collect_files,
    read_json,
    write_json,
    split_input,
    find_execs,
    setup_logging,
    atomic_write,
)
from mycotools.lib.biotools import fa2dict, gff2list, dict2fa, list2gff
from mycotools.download.ncbi import (
    esearch_ncbi,
    esummary_ncbi,
    run_datasets,
    compile_organism_names,
    main as ncbiDwnld,
)
from mycotools.download.jgi import main as jgiDwnld
from mycotools.utils.ncbi2db import main as ncbi2db
from mycotools.utils.jgi2db import main as jgi2db
from mycotools.mtdb.predb import main as predb2mtdb
from mycotools.mtdb.predb import predb_headers, read_predb, gen_omes
from pathlib import Path

logger = logging.getLogger(__name__)


def _read_ledger(file_path, comment="#"):
    """Return the meaningful lines of a ledger/sidecar file.

    Yields the shared read half of the ledger parsers below: returns a list of
    right-stripped lines, skipping comment lines (those beginning with
    `comment`) and lines that are empty once stripped. Returns an empty list if
    the file does not exist, so callers can treat a missing ledger as empty.
    """
    if not Path(file_path).is_file():
        return []
    with open(file_path, "r") as raw:
        return [
            line.rstrip()
            for line in raw
            if not line.startswith(comment) and line.strip()
        ]


def validate_t_and_c(config, discrepancy=False):
    """Validate that user understands and accepts the conditions of
    use-restricted data, and that responsibility over confirming the
    desigination of use-restriction is the user's responsibility"""

    # if there is a configuration json, just query that
    if config and not discrepancy:
        try:
            if config["nonpublished"].lower() in {"yes", "y"}:
                nonpublished = "yes"
            else:
                nonpublished = False  # weird situation for prokaryotes, how to
                # handle?
        except AttributeError:
            nonpublished = False

    # if there isnt a configuration, alert the user to use-restriction policies
    else:
        logger.debug(
            "Please review JGI use-restricted data policy here: "
            + "https://jgi.doe.gov/user-programs/pmo-overview/policies/"
            + "\nPlease review GenBank use-restricted data policy here: "
            + "https://ncbi.nlm.nih.gov/genbank/"
            + "\nPlease review how Mycotools handles use-restricted data here:"
            + " https://github.com/xonq/mycotools/blob/master/MTDB.md"
        )
        check = ""
        if check.lower() not in {"y", "yes"}:
            check = input(
                "\nWARNING: This is a permanent configuration "
                + "option for the primary MTDB. Do you agree to honor "
                + "the terms and conditions "
                + "of use-restricted JGI and GenBank data; acknowledge that "
                + "Mycotools may not comprehensively determine use-restricted "
                + "designations; and acknowledge that you will validate "
                + "any use-restricted assignments in your local MTDB "
                + "prior to publication?"
                + "\n\nPlease type [y]es/[N]o if you acknowledge these terms: "
            )
        if check.lower() not in {"y", "yes"}:
            logger.info("Rerun without --nonpublished")
            sys.exit(1)

        nonpublished = "yes"

    return nonpublished


def gen_config(
    branch="fungi",
    forbidden="",
    repo=None,
    rogue=False,
    nonpublished=False,
    jgi=False,
    rank2lineages={},
):

    config = {
        "forbidden": forbidden,
        "repository": repo,
        "branch": branch,
        "nonpublished": nonpublished,
        "rogue": rogue,
        "jgi": jgi,
        "lineage_constraints": rank2lineages,
    }

    return config


def add_vars(init_dir, dbtype):
    """Initialize the environmental variables for MTDB"""
    mtdb_initialize(init_dir, dbtype, init=True)


def init_db(
    init_dir,
    branch,
    envs,
    dbtype,
    date=None,
    rogue=False,
    nonpublished=False,
    jgi=True,
    repo=None,
    rank2lineages={},
):
    """Initialize database in `init_dir`"""

    new_dirs = [
        init_dir + "data/",
        init_dir + "config/",
        init_dir + "data/fna/",
        init_dir + "data/faa/",
        init_dir + "mtdb",
        init_dir + "data/gff3/",
        init_dir + "log/",
        init_dir + "data/db/",
    ]
    output = prep_output(init_dir, cd=True)
    if not output.endswith("/"):
        output += "/"
    for new_dir in new_dirs:
        if not Path(new_dir).is_dir():
            Path(new_dir).mkdir()

    config = gen_config(
        branch=branch,
        rogue=rogue,
        forbidden="$MYCODB/log/forbidden.tsv",
        nonpublished=nonpublished,
        jgi=jgi,
        repo=repo,
        rank2lineages=rank2lineages,
    )
    write_json(config, init_dir + "config/mtdb.json", indent=1)

    if not rogue:
        # this is a relic, and needs to be adjusted to a central reference if
        # that is ever created
        if not Path(init_dir + "mtdb").is_dir():
            # NEED TO CHANGE FROM SSH TO LINK ONCE OPEN (config['repository'])
            git_exit = subprocess.call(
                [
                    "git",
                    "clone",
                    "git@gitlab.com:xonq/mtdb",
                    init_dir + "mtdb",
                    #                '-b', branch
                ]
            )
            if git_exit != 0:
                logger.error("git clone failed.")
                sys.exit(2)
        else:
            logger.info("mycotoolsdb directory already exists")
        # NEED TO ADD GITIGNORE TO GIT
        if not primary_db():
            logger.error("no YYYYmmdd.mtdb in " + format_path(envs["MYCODB"]))
            sys.exit(3)
    else:
        new_db_path = output + "mtdb/" + mtdb_sql.PRIMARY_DB_NAME
        if not Path(new_db_path).is_file():
            mtdb().to_sql(new_db_path)

    return output, config


def parse_dups(file_path):
    """Retrieve a file containing replicated genomes and ignore these.
    This is important for dereplication of discrepant genus naming between NCBI
    and MycoCosm due to not adhering to conserved genus naming standards and
    updating relic genus names. It appears that MycoCosm will name genera by
    their anamorph occassionally, and not the consensus name - though I assume
    this is also to some extent present in NCBI. Ultimately, a manually curated
    file is necessary for this and should be held in a central repository."""
    duplicates = {}
    for line in _read_ledger(file_path):
        data = [x.rstrip() for x in line.split("\t") if x]
        if data:
            duplicates[data[0]] = [data[1], data[2], data[3]]
    return duplicates


# def add_dups(
#   dup_code, dup_entry, file_path
#  ):
#    edit = dup_code + '\t' + '\t'.join(dup_entry)
#   log_editor(file_path, dup_code, edit)


def acq_forbid_omes(file_path):
    """Parse a file with forbidden ome accessions - ome codes that have been
    used before and are no longer valid"""
    return set(_read_ledger(file_path))


def write_forbid_omes(omes, file_path):
    """Add to a file of forbidden ome accessions so that these are not ever
    used again, even if the codename is removed from the database"""
    # union with any pre-existing relics; atomic_write guards against losing
    # the old data if the write is cancelled midway
    new_relics = set(_read_ledger(file_path)).union(set(omes))
    with atomic_write(file_path) as out:
        out.write("\n".join([str(x) for x in sorted(new_relics)]))


def parse_failed(file_path=None, rerun=False):
    """Parse a file that stores the failed accessions and metadata of the
    attempted acquisition. Return a dictionary that contains the failed
    accession and its metadata."""
    if not Path(file_path).is_file() or rerun:
        with open(file_path, "w") as out:
            out.write("#code\tsource\tversion\tattempt_date")
        return {}
    prev_failed = {}
    for line in _read_ledger(file_path):
        data = [x.rstrip() for x in line.split("\t")]
        while len(data) < 4:
            data.append("")
        prev_failed[data[0]] = {
            "source": data[1],
            "version": data[2],
            "attempt_date": data[3],
        }
    return prev_failed


def parse_jgi2ncbi(file_path):
    """Parse previously collected NCBI to JGI data to limit querying"""
    if not Path(file_path).is_file():
        with open(file_path, "w") as out:
            out.write("#ncbi_acc\tmycocosm_portal")
        return {}
    jgi2ncbi = {}
    for line in _read_ledger(file_path):
        d = line.split("\t")
        jgi2ncbi[d[1].lower()] = d[0]
    return jgi2ncbi


def parse_true_ncbi(file_path):
    """Parse accessions considered to be unique to NCBI"""
    if not Path(file_path).is_file():
        with open(file_path, "w") as out:
            out.write("#ncbi_acc")
        return set()
    return set(_read_ledger(file_path))


def add_true_ncbi(true_ncbi, file_path=None):
    """Add to a ledger of accessions considered to be unique to NCBI"""
    with atomic_write(file_path) as out:
        out.write("#ncbi_acc\n" + "\n".join([str(x) for x in true_ncbi]))


def add_jgi2ncbi(jgi2ncbi, file_path=None):
    """Add to a ledger that seeks to associated NCBI accessions with JGI. This
    is prone to failure given that the field JGI uses to supply their genome
    accession is either absent from some NCBI entries, or is in a different
    field"""
    with atomic_write(file_path) as out:
        out.write("#ncbi_acc\tmycocosm_portal\n")
        for jgi, ncbi in jgi2ncbi.items():
            out.write(ncbi + "\t" + jgi + "\n")


def add_failed(code, source, version, date, file_path):
    """Add genomes to the failed acquisition file, including metadata on where
    it was downloaded, the version attempted to download, and the date of the
    run"""

    version = version.replace("-", "").replace(" 00:00:00", "")
    edit = code + "\t" + source + "\t" + version + "\t" + date
    log_editor(file_path, code, edit)


def dwnld_mycocosm(
    out_file,
    mycocosm_url="https://mycocosm.jgi.doe.gov/ext-api/mycocosm/catalog/"
    + "download-group?flt=&seq=all&pub=all&grp=fungi&srt="
    + "released&ord=desc",
):
    """Download the MycoCosm genome data spreadsheet, format to UTF-8 and
    return a Pandas dataframe of the data"""

    check_curl = find_execs(["curl"], verbose=False)

    if not Path(out_file).is_file():
        for attempt in range(3):
            if check_curl:
                curl_cmd = subprocess.call(
                    ["curl", mycocosm_url, "-o", out_file + ".tmp"],
                    stdout=subprocess.PIPE,
                )
                if not curl_cmd:
                    shutil.move(out_file + ".tmp", out_file)
                    break
            if curl_cmd:
                logger.error("failed to retrieve MycoCosm table")
            else:
                resp = requests.get(url)
                with open(out_file + ".tmp", "wb") as f:
                    f.write(resp.content)
                shutil.move(out_file + ".tmp", out_file)

    try:
        jgi_df = pd.read_csv(out_file, encoding="cp1252")
    except UnicodeDecodeError:
        jgi_df = pd.read_csv(out_file, encoding="latin1")
    except UnicodeDecodeError:
        jgi_df = pd.read_csv(out_file, encoding="utf-8")
    jgi_df.columns = [x.replace('"', "").replace('"', "") for x in jgi_df.columns]

    return jgi_df


def dwnld_ncbi_metadata(
    ncbi_file,
    ncbi_url="https://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/",
    group="eukaryotes",
):
    """Download the NCBI genome data spreadsheet for prokaryotes or
    eukaryotes, and return a Pandas dataframe"""

    ncbi_url = ncbi_url + group + ".txt"
    if not Path(ncbi_file).is_file():
        getTbl = subprocess.call(["curl", ncbi_url, "-o", ncbi_file + ".tmp"])
        shutil.move(ncbi_file + ".tmp", ncbi_file)
    ncbi_df = pd.read_csv(ncbi_file, sep="\t")

    return ncbi_df


def dwnld_data_reports(
    accs, out_dir, api=None, chunk=100, max_attempts=3, exit_code=11, label="GenBank"
):
    """Acquire NCBI assembly data reports; return acc2org, acc2meta, failed.

    `datasets` is called on `chunk` accessions at a time rather than on the
    whole list at once. A from-scratch initialization queries >16,000
    accessions, and one request that large is both slow enough to be dropped
    mid-transfer and all-or-nothing when it is: a single failure discards every
    record. Each chunk downloads into its own directory and is left there, so an
    interrupted run resumes at the first chunk it had not finished."""
    acc2org, acc2meta, failed = {}, {}, []
    empty_chunks = []

    # a run started before chunking landed leaves a single whole-list download
    # here; parse it rather than re-acquiring everything
    legacy_dir = out_dir + "ncbi_dataset/"
    if Path(legacy_dir + "data/assembly_data_report.jsonl").is_file():
        return compile_organism_names(legacy_dir)

    acc_chunks = [accs[i : i + chunk] for i in range(0, len(accs), chunk)]
    # `datasets` draws its own per-call bar, which is meaningless here because
    # it restarts every chunk; it is silenced below in favor of one bar over the
    # whole acquisition
    for chunk_i, acc_chunk in enumerate(
        tqdm(
            acc_chunks,
            total=len(acc_chunks),
            desc=f"{label} reports",
            unit=" chunk",
            disable=len(acc_chunks) < 2,
        )
    ):
        chunk_dir = f"{out_dir}chunk_{chunk_i}/"
        acc_file = f"{chunk_dir}assembly_accs.txt"
        unzip_dir = f"{chunk_dir}ncbi_dataset/"
        zip_path = f"{chunk_dir}ncbi_dataset.zip"
        if not Path(chunk_dir).is_dir():
            Path(chunk_dir).mkdir(parents=True)

        # only reuse a completed chunk that covers exactly these accessions;
        # otherwise the chunk boundaries have shifted since it was written and
        # its contents no longer correspond to this index
        expected = "\n".join(acc_chunk)
        cached = (
            Path(unzip_dir).is_dir()
            and Path(acc_file).is_file()
            and Path(acc_file).read_text() == expected
        )
        if not cached:
            shutil.rmtree(unzip_dir, ignore_errors=True)
            with open(acc_file, "w") as out:
                out.write(expected)
            # every message in this loop is debug: a retry is routine, it would
            # overdraw the bar above, and a chunk that never succeeds is
            # accounted for in the single summary once the bar has finished
            attempts = 0
            while attempts < max_attempts:
                if attempts:
                    logger.debug(f"Reattempting {label} chunk {chunk_i + 1}")
                    if Path(zip_path).is_file():
                        Path(zip_path).unlink()
                attempts += 1
                run_datasets(
                    None,
                    acc_file,
                    chunk_dir,
                    True,
                    api=api,
                    verbose=False,
                    # datasets' own stderr would overdraw the chunk bar above;
                    # it is kept at debug, and the failures below are what the
                    # user is told about
                    mute_stderr=True,
                )
                try:
                    with zipfile.ZipFile(zip_path, "r") as zip_ref:
                        zip_ref.extractall(chunk_dir)
                    Path(zip_path).unlink()
                    break
                except zipfile.BadZipFile:
                    # a truncated archive is corruption rather than absence, and
                    # exhausting the attempts on it is fatal -- so this one is
                    # said out loud
                    if attempts == max_attempts:
                        logger.error(f"{label} chunk {chunk_i + 1} download corrupted")
                        sys.exit(exit_code)
                    logger.debug(f"datasets download corrupted - {attempts}")
                except FileNotFoundError:
                    logger.debug(f"datasets failed - {attempts}")

        # datasets exited without writing a report on every attempt: the chunk
        # holds no retrievable genomes, which is not fatal to the remainder
        if not Path(unzip_dir + "data/assembly_data_report.jsonl").is_file():
            empty_chunks.append(chunk_i + 1)
            continue

        c_acc2org, c_acc2meta, c_failed = compile_organism_names(unzip_dir)
        acc2org.update(c_acc2org)
        acc2meta.update(c_acc2meta)
        failed.extend(c_failed)
        logger.debug(
            f"\t\t{label} chunk {chunk_i + 1}/{len(acc_chunks)}: "
            + f"{len(c_acc2meta)} genomes"
        )

    # one line once the bar is done, rather than a burst of them through it;
    # empty chunks are routine for the RefSeq pass, where most accessions are
    # speculative, so this reports the scale and leaves the detail to -v
    if empty_chunks:
        logger.warning(
            f"{len(empty_chunks)}/{len(acc_chunks)} {label} chunk(s) returned "
            + f"no genomes ({len(empty_chunks) * chunk} accessions at most); "
            + "rerun with -v for the datasets output"
        )
        logger.debug(f"empty {label} chunks: {empty_chunks}")

    return acc2org, acc2meta, failed


def prep_taxa_cols(
    df,
    taxonomy_dir,
    col="#Organism/Name",
    api=None,
    acc2org={},
    max_attempts=3,
    chunk=100,
):

    skip_prep = list(acc2org.keys())
    gca_prep = [x.upper().replace("GCF", "GCA") for x in skip_prep]
    gcf_prep = [x.upper().replace("GCA", "GCF") for x in skip_prep]
    skip = set(gca_prep + gcf_prep)
    if not Path(taxonomy_dir).is_dir():
        Path(taxonomy_dir).mkdir()

    gb_accs = [x for x in list(df["assembly_acc"]) if x not in skip]
    acc2org_n, acc2meta, org_failed = dwnld_data_reports(
        gb_accs,
        taxonomy_dir,
        api=api,
        chunk=chunk,
        max_attempts=max_attempts,
        exit_code=11,
        label="GenBank",
    )
    logger.info(
        "%s %s",
        f"\t\t{len(acc2meta) + len(org_failed)}",
        "genomes queried from GenBank",
    )
    if len(df["assembly_acc"]):
        logger.debug(f'\t\t{len(org_failed)/len(df["assembly_acc"])*100}% failed')

    # check for RefSeq for failed entries
    refseq_dir = taxonomy_dir + "refseq/"
    if not Path(refseq_dir).is_dir():
        Path(refseq_dir).mkdir()
    missing_accs = sorted(
        set(df["assembly_acc"]).difference(
            set(acc2org_n.keys()).union(set(acc2org.keys()))
        )
    )
    reattempt_acc = []
    for acc in missing_accs:
        if acc.upper().startswith("GCA"):
            reattempt_acc.append(acc.upper().replace("GCA_", "GCF_"))
        elif acc.upper().startswith("GCF"):
            reattempt_acc.append(acc.upper().replace("GCF_", "GCA_"))
    logger.debug(f"Checking RefSeq for {len(reattempt_acc)} entries")
    acc2org_rs, acc2meta_rs, org_failed_2 = dwnld_data_reports(
        reattempt_acc,
        refseq_dir,
        api=api,
        chunk=chunk,
        max_attempts=max_attempts,
        exit_code=10,
        label="RefSeq",
    )
    logger.debug(f"{len(acc2meta_rs)} genome(s) queried from RefSeq")

    acc2org, acc2meta = {**acc2org, **acc2org_n, **acc2org_rs}, {
        **acc2meta,
        **acc2meta_rs,
    }

    df["strain"] = ""
    todel = set()
    for i, row in df.iterrows():
        acc = row["assembly_acc"]
        if acc in acc2org:
            df.at[i, "genus"] = acc2org[acc]["genus"]
            df.at[i, "species"] = acc2org[acc]["species"]
            df.at[i, "strain"] = acc2org[acc]["strain"]
        else:
            todel.add(i)

    df = df[~df.index.isin(todel)]

    return df, acc2meta, acc2org


def prep_jgi_cols(jgi_df, name_col="name"):

    jgi_df["publication(s)"] = jgi_df["publication(s)"].astype(str)
    for i, row in jgi_df.iterrows():
        taxa = re.sub(r"[^ a-zA-Z0-9\.]", "", row[name_col]).split()
        jgi_df.at[i, "genus"] = taxa[0].replace(".", "")
        jgi_df.at[i, "publication(s)"] = jgi_df.at[i, "publication(s)"].replace(
            '"""', ""
        )
        if len(taxa) > 1:
            jgi_df.at[i, "species"] = taxa[1]
            if len(taxa) > 2:
                jgi_df.at[i, "strain"] = "".join(taxa[2:])
                vers_search = re.search(r"v(\d+\.\d+$)", jgi_df["strain"][i])
                if vers_search is not None:
                    version = vers_search[1]
                    jgi_df.at[i, "version"] = float(version)
                    jgi_df.at[i, "strain"] = jgi_df["strain"][i][: -len(vers_search[0])]
                else:
                    jgi_df.at[i, "version"] = 0
                jgi_df.at[i, "strain"] = jgi_df.at[i, "strain"].replace(".", "")
        else:
            jgi_df.at[i, "species"] = "sp."

    jgi_df = jgi_df.sort_values(by="version", ascending=False)
    jgi_df = jgi_df.drop_duplicates("portal")

    return jgi_df


def clean_ncbi_df(
    ncbi_df, update_path, kingdom="Fungi", api=None, max_attempts=3, chunk=100
):
    ncbi_df = ncbi_df.astype(str).replace(np.nan, "")

    acc2org_path = update_path + "../gca2org.tsv"
    acc2org = {}
    if Path(acc2org_path).is_file():
        with open(acc2org_path, "r") as raw:
            for line in raw:
                d = line.split("\t")
                acc2org[d[0]] = {
                    "genus": d[1],
                    "species": d[2],
                    "strain": d[3].rstrip(),
                }

    if kingdom.lower() == "fungi":
        # extract group of interest (case sensitive to first letter)
        ncbi_df = ncbi_df[ncbi_df["Group"] == "Fungi"]
    elif kingdom.lower() in {"plants", "viridiplantae"}:
        ncbi_df = ncbi_df[ncbi_df["Group"] == "Plants"]
    elif kingdom.lower() in {"animals", "metazoa"}:
        ncbi_df = ncbi_df[ncbi_df["Group"] == "Animals"]
    elif kingdom.lower() in {"bacteria", "prokaryotes"}:
        ncbi_df = ncbi_df[ncbi_df["Group"] != "Archaea"]
    elif kingdom.lower() in {"archaea"}:
        ncbi_df = ncbi_df[ncbi_df["Group"] == "Archaea"]

    ncbi_df = ncbi_df[ncbi_df["assembly_acc"].str.startswith(("GCA", "GCF"))]
    ncbi_df, acc2meta, acc2org = prep_taxa_cols(
        ncbi_df, update_path + "taxonomy/", api=api, acc2org=acc2org, chunk=chunk
    )

    with atomic_write(acc2org_path) as out:
        for acc, org in acc2org.items():
            org_meta = f'{org["genus"]}\t{org["species"]}\t{org["strain"]}'
            out.write(f"{acc}\t{org_meta}\n")

    # remove entries without sufficient metadata
    ncbi_df = ncbi_df.dropna(subset=["genus"])

    # remove assembly- & annotation-lacking entries
    ncbi_df = ncbi_df[~ncbi_df["Genes"].isin({"-", "", "0"})]
    ncbi_df = ncbi_df[~ncbi_df["Proteins"].isin({"-", "", "0"})]
    ncbi_df = ncbi_df[~ncbi_df["assembly_acc"].isin({"-", "", "0"})]

    # sort by version, keep most recent version of duplicate assemblies
    ncbi_df["version"] = pd.to_datetime(ncbi_df["Modify Date"])
    ncbi_df = ncbi_df.sort_values(by="version", ascending=False)
    ncbi_df = ncbi_df.drop_duplicates("assembly_acc")

    # check for explicit new versions of assemblies
    ncbi_df["assembly_base"] = [x[: x.find(".")] for x in list(ncbi_df["assembly_acc"])]
    ncbi_df["assembly_vers"] = [x[x.find(".") :] for x in list(ncbi_df["assembly_acc"])]
    ncbi_df = ncbi_df.sort_values(by="assembly_vers", ascending=False)
    ncbi_df = ncbi_df.drop_duplicates("assembly_base")

    return ncbi_df, acc2meta


def exec_rm_overlap(ncbi_df, todel_i):
    """Execute the removal of JGI redundancy from NCBI genomes and create a
    DataFrame that represents the overlapping genomes"""
    todel = ncbi_df.index.intersection(todel_i)
    ncbi_jgi_overlap = pd.DataFrame(columns=ncbi_df.columns)
    ncbi_jgi_overlap = ncbi_df.loc[todel]
    ncbi_df = ncbi_df.drop(todel)
    return ncbi_jgi_overlap, ncbi_df


def rm_ncbi_overlap(ncbi_df, mycocosm_df, jgi2ncbi, fails=set(), acc2meta={}, api=3):
    """Acquire MycoCosm assembly accessions from NCBI.
    pd.DataFrame(ncbi_df) = post clean ncbi_df;
    set(mycocosm_omes) = set of lower-cased mycocosm genome codes;
    api = number of iterations before sleeping (3/10);
    First, obtain the genome UID by esearching via Entrez;
    Next, use the genome UID to obtain the genome summary via Entrez;
    Finally, relate the accession to the MycoCosm dataframe"""
    version_comp = re.compile(r"\s[Vv]\d+\.\d+")
    mycocosm_omes = set([x.lower() for x in list(mycocosm_df["portal"])])

    todel = []
    ncbi2jgi, jgi2biosample = {v: k for k, v in jgi2ncbi.items()}, {}
    ass_count = 0

    # could vectorize these
    jgi_names = {
        f"{v['genus']}_{v['species']}_{v['strain']}" for k, v in mycocosm_df.iterrows()
    }
    jgi_gen_sp = {f"{v['genus']}_{v['species']}" for k, v in mycocosm_df.iterrows()}
    for i, row in tqdm(ncbi_df.iterrows(), total=len(ncbi_df)):
        if row["assembly_acc"] in ncbi2jgi:
            jgi2biosample[ncbi2jgi[row["assembly_acc"]]] = row["BioSample Accession"]
            todel.append(i)
        #        elif row['assembly_acc'] in fails:
        #           pass
        #            todel.append(i)
        elif row["assembly_acc"] not in fails:
            if row["assembly_acc"] in acc2meta:
                ass_name = acc2meta[row["assembly_acc"]]["accession"]
                submitter = acc2meta[row["assembly_acc"]]["submitter"]
            else:
                ass_uid = esearch_ncbi(row["assembly_acc"], "assembly", "assembly")
                if not ass_uid:
                    fails.add(row["assembly_acc"])
                    continue
                ncbi_df.at[i, "uid"] = ass_uid
                # this is a huge bottleneck, taking approximately 1s per query
                summary = esummary_ncbi(max(ass_uid), "assembly")
                ass_name_prep = summary["DocumentSummarySet"]["DocumentSummary"][0][
                    "AssemblyName"
                ]
                ass_name = version_comp.sub("", ass_name_prep)
                submitter = summary["DocumentSummarySet"]["DocumentSummary"][0][
                    "SubmitterOrganization"
                ]
            if ass_name.lower() in mycocosm_omes:
                todel.append(i)
                jgi2ncbi[ass_name.lower()] = row["assembly_acc"]
                jgi2biosample[ass_name.lower()] = row["BioSample Accession"]
            elif any(
                x in submitter.lower()
                for x in ["joint genome institute", "jgi", "joint_genome_institute"]
            ):  # very crude, but mycocosm does not give the option to be
                # systematic
                todel.append(i)
                jgi2ncbi[ass_name.lower() + f"${ass_count}"] = row["assembly_acc"]
                jgi2biosample[ass_name.lower() + f"${ass_count}"] = row[
                    "BioSample Accession"
                ]
                ass_count += 1
            elif f"{row['genus']}_{row['species']}_{row['strain']}" in jgi_names:
                todel.append(i)
                jgi2ncbi[f"{row['genus']}_{row['species']}_{row['strain']}"] = row[
                    "assembly_acc"
                ]
                jgi2biosample[
                    f"{row['genus']}_{row['species']}_{row['strain']}_{ass_count}"
                ] = row["BioSample Accession"]
                ass_count += 1
            elif (
                not row["strain"].rstrip()
                and f'{row["genus"]}_{row["species"]}' in jgi_gen_sp
            ):
                # unfortunately we cannot validate that this is either a JGI or
                # NCBI unique specimen, but to err on the side of caution we
                # should remove it.
                todel.append(i)
                gen_sp = f'{row["genus"]}_{row["species"]}'
                jgi2ncbi[gen_sp] = row["assembly_acc"]
                jgi2biosample[gen_sp] = row["BioSample Accession"]
            else:
                fails.add(row["assembly_acc"])
    #    for i in reversed(todel):
    #       ncbi_jgi_overlap = pd.concat([ncbi_jgi_overlap, ncbi_df.loc[i]])
    #        ncbi_df = ncbi_df.drop(i)

    ncbi_df, ncbi_jgi_overlap = exec_rm_overlap(ncbi_df, todel)

    return ncbi_df, jgi2ncbi, jgi2biosample, fails, ncbi_jgi_overlap, todel


def write_primary(db, date, update_path=None):
    """Install `db` as the primary MTDB.

    The primary is the SQLite database at `$MYCODB/mtdb.db`, written atomically
    so an interrupted update can never leave a partial database where
    `primary_db()` would pick it up. The outgoing primary is archived under
    `log/<date>/` first, alongside a `.mtdb` snapshot of the new one -- every
    historical primary stays readable with nothing but a text editor."""
    new_path = format_path("$MYCODB/" + mtdb_sql.PRIMARY_DB_NAME)
    prior = primary_db(verbose=False)

    # copy, rather than move, so a failed write leaves the old primary in place
    if update_path and prior and Path(prior).is_file():
        archive = update_path + Path(prior).name
        if format_path(prior) != format_path(archive):
            shutil.copy(prior, archive)

    db.to_sql(new_path)

    if update_path:
        db.df2db(update_path + date + ".mtdb", headers=True)
    # a dated flat primary predates the SQLite backend; it has been archived, so
    # drop it rather than leave a stale database beside the real one
    if prior and Path(prior).is_file() and format_path(prior) != format_path(new_path):
        Path(prior).unlink()
    return new_path


def mk_wrk_dirs(update_path):
    """Make the download directories in the update path"""
    wrk_dirs = ["faa/", "fna/", "gff3/"]
    for wrk_dir in wrk_dirs:
        if not Path(update_path + wrk_dir).is_dir():
            Path(update_path + wrk_dir).mkdir()


def prepare_ref_db(ref_db, date):
    """Prepare MTDBs of the data in the reference MTDB that are from the two
    download sources"""
    ref_db["acquisition_date"] = [date for x in ref_db["acquisition_date"]]
    ref_db = ref_db.set_index()

    jgi = mtdb(
        {k: v for k, v in ref_db.items() if v["source"].lower() == "jgi"}, index="ome"
    )
    ncbi = mtdb(
        {k: v for k, v in ref_db.items() if v["source"].lower() == "ncbi"}, index="ome"
    )
    if set(ref_db.keys()).difference(set(jgi.keys()).union(set(ncbi.keys()))):
        logger.warning(
            '\tWARNING: reference entries that are not labeled "jgi/ncbi" are excluded'
        )

    return jgi.mtdb2pd(), ncbi.mtdb2pd()


def internal_redundancy_check(db):
    """Check the inputted database for overlapping assembly accessions and
    dereplicate, including those with different versions of the same
    accession"""
    db = mtdb.pd2mtdb(db).set_index()
    ncbi_db = {k: v for k, v in db.items() if v["source"] == "ncbi"}
    jgi_db = {k: v for k, v in db.items() if v["source"] == "jgi"}

    ncbi_accs = defaultdict(list)
    for ome, row in ncbi_db.items():
        ass_acc = row["assembly_acc"][: row["assembly_acc"].find(".")]
        ncbi_accs[ass_acc].append(ome)
    red_ncbi = {k: v for k, v in ncbi_accs.items() if len(v) > 1}
    for ass_acc, omes in red_ncbi.items():
        omes = sorted(omes, reverse=True)  # large ome number to small
        accs = defaultdict(list)
        for ome in omes:
            try:
                acc_ver = int(
                    ncbi_db[ome]["assembly_acc"][
                        ncbi_db[ome]["assembly_acc"].find(".") + 1 :
                    ]
                )
            except ValueError:
                acc_ver = 0
            accs[acc_ver].append(ome)
        max_ver = max(accs.keys())
        for ver, omes in accs.items():
            if ver != max_ver:
                for ome in omes:
                    del db[ome]
            else:
                for ome in omes[1:]:
                    del db[ome]

    jgi_accs = defaultdict(list)
    for ome, row in jgi_db.items():
        ass_acc = row["assembly_acc"]
        jgi_accs[ass_acc].append(ome)
    red_jgi = {k: v for k, v in jgi_accs.items() if len(v) > 1}
    for ass_acc, omes in red_jgi.items():
        omes = sorted(omes, reverse=True)
        accs = defaultdict(list)
        for ome in omes:
            acc_ver = jgi_db[ome]["version"].replace("v", "").replace("V", "")
            try:
                ver = int(float(acc_ver))
            except ValueError:
                ver = 0
            accs[ver].append(ome)
        max_ver = max(accs.keys())
        for ver, omes in accs.items():
            if ver != max_ver:
                for ome in omes:
                    del db[ome]
            else:
                for ome in omes[1:]:
                    del db[ome]

    return db2df(db.reset_index())


def read_prev_tax(tax_path):
    """Open a genus to taxonomy JSON path"""
    tax_dicts = {}
    for line in _read_ledger(tax_path):
        data = line.split("\t")
        tax_dicts[data[0]] = json.loads(data[1])
    return tax_dicts


def ref_update(
    ref_db,
    update_path,
    date,
    rerun,
    jgi_email,
    jgi_pwd,
    config,
    ncbi_api,
    cpus=1,
    check_MD5=True,
    jgi=True,
    group="eukaryotes",
    kingdom="Fungi",
    remove=True,
    taxonomy=True,
    chunk=100,
    tape_wait=None,
):
    """Initialize/Update the primary MTDB based on a reference database
    acquired external from any existing primary MTDB"""
    # NEED to mark none for new databases' refdb
    # initialize update
    logger.info("Initializing run")
    mk_wrk_dirs(update_path)

    jgi_df, ncbi_df = prepare_ref_db(ref_db, date)

    # run JGI
    if jgi and len(jgi_df) > 0:
        logger.info("Assimilating MycoCosm")
        jgi_predb_path = update_path + date + ".jgi.predb2.mtdb"

        if not Path(jgi_predb_path).is_file():
            logger.info("Downloading MycoCosm data")
            jgi_deferred = set()
            post_jgi_df, jgi_dwnld_failed = jgiDwnld(
                jgi_df,
                update_path,
                jgi_email,
                jgi_pwd,
                deferred=jgi_deferred,
                restore_wait=tape_wait,
                defer_tape=True,
            )
            jgi_predb = post_jgi_df.rename(
                columns={
                    "published(s)": "published",
                    "fna_path": "assemblyPath",
                    "gff3_path": "gffPath",
                }
            )

            logger.info("Curating MycoCosm data")
            jgi_premtdb = jgi_predb.fillna("").to_dict(orient="list")
            # jgiDwnld reports bare portal ids; the failed ledger records
            # [accession, version] pairs. Genomes only awaiting a JGI tape
            # restore are retryable, so they are not recorded as failures
            versions = dict(zip(jgi_df["assembly_acc"], jgi_df["version"]))
            jgi_failed = [
                [acc, versions.get(acc, "")]
                for acc in jgi_dwnld_failed
                if acc not in jgi_deferred
            ]
            if jgi_deferred:
                logger.info(
                    f"\t{len(jgi_deferred)} genome(s) awaiting JGI tape restore; "
                    "they will be retried on the next run"
                )
            # a downloaded assembly path is required to curate; if no JGI genome
            # was successfully retrieved (e.g. all portals failed) skip curation
            # rather than raising a KeyError and aborting the whole run
            if "assemblyPath" in jgi_premtdb:
                jgi_mtdb, jgi_failed1 = predb2mtdb(
                    jgi_premtdb,
                    mtdb(),
                    update_path,
                    #                                            forbidden = forbid_omes,
                    cpus=cpus,
                    remove=remove,
                    spacer="\t\t",
                )
                jgi_failed.extend(jgi_failed1)
            else:
                logger.warning(
                    "No JGI assemblies downloaded; skipping MycoCosm curation"
                )
                jgi_mtdb = mtdb()
            jgi_mtdb.df2db(jgi_predb_path)
            for failure in jgi_failed:
                add_failed(
                    failure[0],
                    "jgi",
                    str(failure[1]),
                    date,
                    format_path("$MYCODB/../log/failed.tsv"),
                )

        else:
            jgi_mtdb = mtdb(jgi_predb_path)

    else:
        jgi_mtdb = mtdb()
    new_db = jgi_mtdb.mtdb2pd()

    logger.info("Assimilating NCBI")
    if not Path(update_path + date + ".ncbi.predb").is_file():
        logger.info("Downloading NCBI data")
        ncbi_predb, ncbi_failed1 = ncbiDwnld(
            assembly=True,
            proteome=False,
            gff3=True,
            ncbi_df=ncbi_df,
            remove=True,
            output_path=update_path,
            column="assembly_acc",
            ncbi_column="genome",
            check_MD5=check_MD5,
            verbose=True,
            chunk=chunk,
        )

        for failure in ncbi_failed1:
            add_failed(
                failure[0],
                "ncbi",
                str(failure[1]),
                date,
                format_path("$MYCODB/../log/failed.tsv"),
            )

        #        for dup in new_dups:
        #           add_dups(dup, new_dups[dup], format_path('$MYCODB/../log/duplicates.tsv'))
        ncbi_predb.to_csv(update_path + date + ".ncbi.predb", sep="\t", index=None)
    else:
        #        refdbncbi = mtdb(update_path + date + '.ncbi.ref.mtdb')
        ncbi_predb = pd.read_csv(update_path + date + ".ncbi.predb", sep="\t")

    logger.info("Curating NCBI data")
    if not Path(update_path + date + ".ncbi.predb2.mtdb").is_file():
        for key in ncbi_predb.columns:
            ncbi_predb[key] = ncbi_predb[key].fillna("")
        ncbi_predb["version"] = ncbi_predb["version"].astype(str)
        ncbi_premtdb = ncbi_predb.to_dict(orient="list")
        ncbi_mtdb, ncbi_failed2 = predb2mtdb(
            ncbi_premtdb,
            mtdb(),
            update_path,
            #                                       forbidden = forbid_omes,
            cpus=cpus,
            remove=remove,
            spacer="\t\t",
        )
        for failure in ncbi_failed2:
            add_failed(
                failure[0],
                "ncbi",
                str(failure[1]),
                date,
                format_path("$MYCODB/../log/failed.tsv"),
            )
        ncbi_mtdb.df2db(update_path + date + ".ncbi.predb2.mtdb")
    else:
        ncbi_mtdb = mtdb(f"{update_path}{date}.ncbi.predb2.mtdb")

    try:
        ncbi_db = db2df(update_path + date + ".ncbi.predb2.mtdb")
    except pd.errors.EmptyDataError:
        ncbi_db = pd.DataFrame({x: [] for x in refdbncbi.keys()})
    if len(ncbi_db) > 0:
        #        df2db(ncbi_db, ncbi_db_path)
        new_db = pd.concat([new_db, ncbi_db])

    logger.info("Assimilating NCBI taxonomy data")
    new_mtdb = mtdb.pd2mtdb(new_db)

    if kingdom.lower() == "fungi":
        rank = "kingdom"
    else:
        rank = "superkingdom"

    if jgi_mtdb and ncbi_mtdb:
        update_mtdb = mtdb(
            {**jgi_mtdb.set_index(), **ncbi_mtdb.set_index()}, index="ome"
        )
    elif ncbi_mtdb:
        update_mtdb = ncbi_mtdb
    elif jgi_mtdb:
        update_mtdb = jgi_mtdb
    else:
        logger.info("No updates")
        sys.exit(0)

    if taxonomy:  # already completed
        tax_path = f"{update_path}../taxonomy.tsv"
        tax_dicts = read_prev_tax(tax_path)
        tax_dicts = gather_taxonomy(
            new_mtdb,
            api_key=ncbi_api,
            king=kingdom,
            rank=rank,
            output_path=tax_path,
            tax_dicts=tax_dicts,
        )
        new_mtdb, genus_dicts = assimilate_tax(new_mtdb, tax_dicts)

        for ome, row in update_mtdb.items():
            if row["genus"] in genus_dicts:
                row["taxonomy"] = genus_dicts[row["genus"]]

    return new_mtdb, update_mtdb


def extract_constraint_lineages(
    df, ncbi_api, kingdom, lineage_constraints, tax_dicts, tax_path
):
    """Extract genera from NCBI and JGI Pandas dataframes
    that hit a dictionary of lineages of interest"""

    # begin extracting lineages of interest and store tax_dicts for later
    if "taxonomy" not in df.columns:
        df["taxonomy"] = [{} for x in range(len(df))]
    if kingdom.lower() in {"fungi", "plants", "animals", "metazoa", "viridiplantae"}:
        query_rank = "kingdom"
    else:
        query_rank = "superkingdom"

    # skip querying ncbi if we are only gathering the genus
    if not set(lineage_constraints.keys()).difference({"genus"}):
        passing_tax = set(x[0].upper() + x[1:] for x in lineage_constraints["genus"])
        df = df[df["genus"].isin(passing_tax)]
        return tax_dicts, df

    tax_dicts = gather_taxonomy(
        df,
        api_key=ncbi_api,
        king=kingdom,
        rank=query_rank,
        tax_dicts=tax_dicts,
        output_path=tax_path,
    )

    lineage_constraints = {k: set(v) for k, v in lineage_constraints.items()}

    # extract genera that pass
    passing_tax = set()
    if "genus" in lineage_constraints:
        # add constraint genera immediately
        passing_tax = set(x[0].upper() + x[1:] for x in lineage_constraints["genus"])
        lineage_constraints = {
            k: v for k, v in lineage_constraints.items() if k != "genus"
        }

    for genus, tax in tax_dicts.items():
        for rank, lineages in lineage_constraints.items():
            if rank in tax:
                if tax[rank].lower() in lineages:
                    passing_tax.add(genus)

    df = df[df["genus"].isin(passing_tax)]
    return tax_dicts, df


def taxonomy_update(
    orig_db,
    update_path,
    date,
    config,
    ncbi_api,
    rank="kingdom",
    group="fungi",
):
    """Reset the taxonomy for the entire database and overwrite the previous
    tax path data to accomodate new taxonomy"""
    taxless_db = orig_db.reset_index()
    taxless_db["taxonomy"] = [{} for x in taxless_db["taxonomy"]]
    tax_path = f"{update_path}../taxonomy.tsv"
    gca_path = f"{update_path}../gca2org.tsv"
    if Path(tax_path).is_file():
        Path(tax_path).rename(update_path + "old_taxonomy.tsv")
    if Path(gca_path).is_file():
        Path(gca_path).rename(update_path + "old_gca2org.tsv")
    tax_dicts = gather_taxonomy(
        taxless_db, api_key=ncbi_api, king=group, rank=rank, output_path=tax_path
    )
    tax_db, genus_dicts = assimilate_tax(taxless_db, tax_dicts)
    if not isinstance(tax_db, mtdb):
        return tax_db, mtdb.pd2mtdb(tax_db)
    else:
        return tax_db.mtdb2pd(), tax_db


def rogue_update(
    db,
    update_path,
    date,
    rerun,
    jgi_email,
    jgi_pwd,
    config,
    ncbi_api,
    cpus=1,
    check_MD5=True,
    jgi=True,
    group="eukaryotes",
    kingdom="Fungi",
    remove=True,
    lineage_constraints={},
    chunk=100,
    tape_wait=None,
):
    """Initialize/update a standalone primary MTDB"""
    # NEED to mark none for new databases' refdb
    # initialize update
    logger.info("Initializing run")
    mk_wrk_dirs(update_path)
    prev_failed = parse_failed(
        rerun=rerun, file_path=format_path("$MYCODB/../log/failed.tsv")
    )
    duplicates = parse_dups(format_path("$MYCODB/../log/duplicates.tsv"))
    forbid_omes = acq_forbid_omes(file_path=format_path("$MYCODB/../log/relics.txt"))

    db = internal_redundancy_check(db)

    # prepare ncbi_df
    if ncbi_api:
        api = 10
    else:
        api = 3
    ncbi_db_path = update_path + date + ".ncbi.mtdb"
    pre_ncbi_df0 = dwnld_ncbi_metadata(update_path + date + ".ncbi.tsv", group=group)
    pre_ncbi_df1 = pre_ncbi_df0.rename(columns={"Assembly Accession": "assembly_acc"})
    logger.info("Acquiring NCBI metadata")
    ncbi_df, acc2meta = clean_ncbi_df(
        pre_ncbi_df1, update_path, kingdom=kingdom, api=ncbi_api, chunk=chunk
    )

    # begin extracting lineages of interest and store tax_dicts for later
    tax_path = f"{update_path}../taxonomy.tsv"
    tax_dicts = read_prev_tax(tax_path)
    #    tax_dicts = {v['genus']: v['taxonomy'] for k, v in db.iterrows() \
    #                if any(y for x, y in v['taxonomy'].items() \
    #                      if x not in {'genus', 'species', 'strain'})}
    if lineage_constraints:
        lineage_path = update_path + date + ".ncbi.posttax.df"
        if not Path(lineage_path).is_file():
            logger.info("Extracting lineages from NCBI")
            # NEED to transition to datasets
            tax_dicts, ncbi_df = extract_constraint_lineages(
                ncbi_df, ncbi_api, kingdom, lineage_constraints, tax_dicts, tax_path
            )
            ncbi_df.to_csv(lineage_path, sep="\t", index=None)
        else:
            ncbi_df = pd.read_csv(lineage_path, sep="\t")
            ncbi_df["version"] = pd.to_datetime(ncbi_df["version"])

    old_len = len(db["ome"])
    new_len = len(db["ome"])
    if old_len - new_len:
        logger.debug("" + str(old_len - new_len) + " redundant entries removed")

    # run JGI
    if jgi:
        logger.info("Assimilating MycoCosm (1 download/minute)")
        jgi_db_path = update_path + date + ".jgi.mtdb"
        mycocosm_path = update_path + date + ".mycocosm.csv"

        # acquire the mycocosm master table
        jgi_df = dwnld_mycocosm(mycocosm_path)
        jgi_df["biosample"] = ""
        jgi_df = prep_jgi_cols(jgi_df, "name")

        # extract JGI lineages of interest and store tax_dicts for later
        if lineage_constraints:
            lineage_path = update_path + date + ".jgi.posttax.df"
            if not Path(lineage_path).is_file():
                logger.info("Extracting lineages from MycoCosm")
                tax_dicts, jgi_df = extract_constraint_lineages(
                    jgi_df, ncbi_api, kingdom, lineage_constraints, tax_dicts, tax_path
                )
                jgi_df.to_csv(lineage_path, sep="\t", index=None)
            else:
                jgi_df = pd.read_csv(lineage_path, sep="\t")

        logger.info("Searching NCBI for MycoCosm overlap")
        jgi_ncbi_overlap_file = f"{update_path}/redundant_ncbi.tsv"
        jgi2ncbi = parse_jgi2ncbi(update_path + "../jgi2ncbi.tsv")
        ncbi_df = ncbi_df.set_index("assembly_acc", drop=False)
        if Path(jgi_ncbi_overlap_file).is_file():
            with open(jgi_ncbi_overlap_file, "r") as raw:
                todel_i = [x.rstrip() for x in raw]
            ncbi_df, ncbi_jgi_overlap = exec_rm_overlap(ncbi_df, todel_i)
        else:
            true_ncbi = parse_true_ncbi(update_path + "../supported_ncbi.tsv")
            ncbi_df, jgi2ncbi, jgi2biosample, true_ncbi, ncbi_jgi_overlap, todel_i = (
                rm_ncbi_overlap(ncbi_df, jgi_df, jgi2ncbi, true_ncbi, acc2meta, api=api)
            )

            logger.debug("" + str(len(jgi2ncbi)) + " overlapping genomes")
            with atomic_write(jgi_ncbi_overlap_file) as out:
                out.write("\n".join([x for x in todel_i]))
            add_true_ncbi(true_ncbi, update_path + "../supported_ncbi.tsv")
            add_jgi2ncbi(jgi2ncbi, update_path + "../jgi2ncbi.tsv")
            for i, row in jgi_df.iterrows():
                if (
                    row["portal"].lower() in jgi2ncbi
                    and row["portal"].lower() in jgi2biosample
                ):
                    jgi_df.at[i, "biosample"] = jgi2biosample[row["portal"].lower()]

        logger.info("Downloading MycoCosm data")
        jgi_predb_path = update_path + date + ".jgi.predb2.mtdb"
        jgi_predb, db, jgi_failed, jgi_deferred = jgi2db(
            jgi_df,
            db,
            update_path,
            jgi_email,
            jgi_pwd,
            date=date,
            nonpublished=config["nonpublished"],
            rerun=rerun,
            failed_dict=prev_failed,
            jgi2ncbi=jgi2ncbi,
            repeatmasked=True,
            restore_wait=tape_wait,
        )  # download JGI files and ready predb

        # get ncbi hits that hit jgi runs which did not yield data - genomes
        # awaiting a JGI tape restore included, so NCBI covers them until the
        # next run retrieves the MycoCosm copy
        failed_ncbi2jgi = {
            jgi2ncbi[f[0]]: f[0] for f in jgi_failed + jgi_deferred if f[0] in jgi2ncbi
        }
        ncbi_jgi_overlap = ncbi_jgi_overlap[
            ncbi_jgi_overlap["assembly_acc"].isin(failed_ncbi2jgi)
        ]
        ncbi_df = pd.concat([ncbi_df, ncbi_jgi_overlap])

        refdbjgi = mtdb.pd2mtdb(db)
        if not Path(jgi_predb_path).is_file():
            logger.info("Curating MycoCosm data")
            jgi_premtdb = jgi_predb.fillna("").to_dict(orient="list")
            if "assemblyPath" in jgi_premtdb:
                jgi_mtdb, jgi_failed1 = predb2mtdb(
                    jgi_premtdb,
                    refdbjgi,
                    update_path,
                    forbidden=forbid_omes,
                    cpus=cpus,
                    remove=remove,
                    spacer="\t\t",
                )
                jgi_failed.extend(jgi_failed1)
            else:
                jgi_mtdb = mtdb()
            jgi_mtdb.df2db(jgi_predb_path)

            # only genuine failures are blacklisted; genomes pending a JGI tape
            # restore are absent from jgi_failed so the next run retries them
            for failure in jgi_failed:
                add_failed(
                    failure[0],
                    "jgi",
                    str(failure[1]),
                    date,
                    format_path("$MYCODB/../log/failed.tsv"),
                )
        else:
            jgi_mtdb = mtdb(jgi_predb_path)

        try:
            jgi_db = db2df(jgi_predb_path)
        except pd.errors.EmptyDataError:  # empty JGI predb
            jgi_db = pd.DataFrame({x: [] for x in refdbjgi.keys()})

        new_db_path = update_path + date + ".checkpoint.jgi.mtdb"
        if not Path(new_db_path).is_file():
            if len(jgi_db) > 0:
                df2db(jgi_db, jgi_db_path)
                if not db is None:
                    new_db = pd.concat([jgi_db, db])
                else:
                    new_db = jgi_db
            else:
                new_db = db
            df2db(new_db, new_db_path)
        else:
            new_db = db2df(new_db_path)
    else:
        jgi_mtdb = None
        new_db = db

    logger.info("Assimilating NCBI (10 download/second w/API key, 3 w/o)")
    new_db["version"] = new_db["version"].astype(str)
    if not Path(update_path + date + ".ncbi.predb").is_file():
        #    if not os.path.isfile(update_path + date + '.ncbi.predb'):
        logger.info("Downloading NCBI data")
        ncbi_predb, new_db, ncbi_failed1 = ncbi2db(
            update_path,
            ncbi_df,
            ref_db=new_db,
            date=date,
            failed_dict=prev_failed,
            rerun=rerun,
            duplicates=duplicates,
            check_MD5=check_MD5,
            chunk=chunk,
        )

        for failure in ncbi_failed1:
            add_failed(
                failure[0],
                "ncbi",
                str(failure[1]),
                date,
                format_path("$MYCODB/../log/failed.tsv"),
            )
        #        for dup in new_dups:
        #           add_dups(dup, new_dups[dup], format_path('$MYCODB/../log/duplicates.tsv'))
        refdbncbi = mtdb.pd2mtdb(new_db)
        refdbncbi.df2db(update_path + date + ".ncbi.ref.mtdb")
        ncbi_predb.to_csv(update_path + date + ".ncbi.predb", sep="\t", index=None)
    else:
        refdbncbi = mtdb(update_path + date + ".ncbi.ref.mtdb")
        ncbi_predb = pd.read_csv(update_path + date + ".ncbi.predb", sep="\t")

    logger.info("Curating NCBI data")
    if not Path(update_path + date + ".ncbi.predb2.mtdb").is_file():
        for key in ncbi_predb.columns:
            ncbi_predb[key] = ncbi_predb[key].fillna("")
        ncbi_predb["version"] = ncbi_predb["version"].astype(str)
        ncbi_premtdb = ncbi_predb.to_dict(orient="list")
        for col in predb_headers:
            if col not in ncbi_premtdb:
                ncbi_premtdb[col] = ["" for x in ncbi_predb[key]]
        ncbi_premtdb["restriction"] = ["0" for x in ncbi_premtdb["assemblyPath"]]
        ncbi_mtdb, ncbi_failed2 = predb2mtdb(
            ncbi_premtdb,
            refdbncbi,
            update_path,
            forbidden=forbid_omes,
            cpus=cpus,
            remove=remove,
            spacer="\t\t",
        )
        for failure in ncbi_failed2:
            add_failed(
                failure[0],
                "ncbi",
                str(failure[1]),
                date,
                format_path("$MYCODB/../log/failed.tsv"),
            )
        ncbi_mtdb.df2db(update_path + date + ".ncbi.predb2.mtdb")

    try:
        ncbi_db = db2df(update_path + date + ".ncbi.predb2.mtdb")
        ncbi_mtdb = mtdb(update_path + date + ".ncbi.predb2.mtdb")
    except pd.errors.EmptyDataError:
        ncbi_db = pd.DataFrame({x: [] for x in refdbncbi.keys()})
        ncbi_mtdb = mtdb()
    if len(ncbi_db) > 0:
        df2db(ncbi_db, ncbi_db_path)
        new_db = pd.concat([new_db, ncbi_db])

    new_mtdb = mtdb.pd2mtdb(new_db)

    logger.info("Assimilating NCBI taxonomy data")
    if kingdom.lower() == "fungi":
        rank = "kingdom"
    else:
        rank = "superkingdom"

    tax_path = f"{update_path}../taxonomy.tsv"
    tax_dicts = gather_taxonomy(
        new_mtdb,
        api_key=ncbi_api,
        king=kingdom,
        rank=rank,
        tax_dicts=tax_dicts,
        output_path=tax_path,
    )
    new_mtdb, genus_dicts = assimilate_tax(new_mtdb, tax_dicts)

    if jgi_mtdb and ncbi_mtdb:
        update_mtdb = mtdb(
            {**jgi_mtdb.set_index(), **ncbi_mtdb.set_index()}, index="ome"
        )
    elif ncbi_mtdb:
        update_mtdb = ncbi_mtdb
    elif jgi_mtdb:
        update_mtdb = jgi_mtdb
    else:
        logger.info("No updates")
        sys.exit(0)

    return new_mtdb, update_mtdb


def rm_raw_data(out_dir):
    """Remove raw data after completion"""
    for i in ["faa", "gff3", "gff", "xml", "fna"]:
        if Path(out_dir + i).is_dir():
            shutil.rmtree(out_dir + i)


def gen_algn_db(update_path, omes):
    """Generate an alignment database for the complete primary MTDB"""
    date = Path(os.path.abspath(update_path)).name
    fas = collect_files(os.environ["MYCOFAA"] + "/", ".faa")
    fas = [x for x in fas if Path(x).name[:-6] in omes]
    mkdb_base = "cat " + " ".join(fas)
    mkdb_blast = (
        mkdb_base
        + " | makeblastdb -in -"
        + " -out "
        + os.environ["MYCOGFF3"]
        + "../db/"
        + date
        + ".db -parse_seqids -dbtype prot -title "
        + date
        + ".db"
    )
    #  mkdb_mmseqs = mkdb_base + ' | mmseqs createdb stdin ' + \
    #       format_path('$MYCOFAA/' + date + '.mmseqs.db') + '; ' + \
    #        'mmseqs createdb ' + format_path('$MYCOFAA/' + date + \
    #     '.mmseqs.db') + ' tmp'
    with open(update_path + date + "_makeblastdb.sh", "w") as out:
        out.write(mkdb_blast)
    # with open(update_path + date + '_mmseqsdb.sh', 'w') as out:
    #   out.write(mkdb_mmseqs)

    logger.debug(
        "OPTIONAL: To generate blastdb | mmseqsdb, run the following"
        + "\nbash "
        + update_path
        + date
        + "_makeblastdb.sh"
    )
    # bash ' + update_path \
    #   + date + '_mmseqsdb.sh')


def check_add_mtdb(orig_mtdb, add_mtdb, update_path, overwrite=True):
    """Check the original MTDB for overlapping omes and curate if necessary"""
    orig_mtdb = orig_mtdb.set_index()
    add_mtdb = add_mtdb.set_index()
    orig_aa2ome = {v["assembly_acc"]: k for k, v in orig_mtdb.items()}

    #    failed_aas = []
    overwrite_omes = []
    for assembly_acc, row in add_mtdb.items():
        if assembly_acc in orig_aa2ome:
            #           failed_aas.append(assembly_acc)
            if overwrite:
                overwrite_omes.append(orig_aa2ome[assembly_acc])
            else:
                overwrite_omes.append(row["ome"])

    #    if failed_aas:
    if overwrite:
        for ome in overwrite_omes:
            del orig_mtdb[ome]
    else:
        for ome in overwrite_omes:
            del add_mtdb[ome]
    #        eprint('\nERROR: assembly accessions ("assembly_acc") must be ' \
    #              'unique between databases: ', flush = True)
    #      eprint(', '.join(failed_aas), flush = True)
    #     sys.exit(123)

    orig_omes = set(orig_mtdb.keys())
    new_omes = set(add_mtdb.keys())
    inter_omes = orig_omes.intersection(new_omes)

    # if there are overlapping omes between the addition MTDB and existing
    if inter_omes:
        # prepare to create new omes for the overlapping names
        need_ome_mtdb = mtdb({k: add_mtdb[k] for k in sorted(inter_omes)}, index="ome")
        # these genomes have unique codenames and can be submitted as-is
        fine_ome_mtdb = mtdb(
            {k: v for k, v in add_mtdb.items() if k not in inter_omes}, index="ome"
        )
        need_ome_mtdb = need_ome_mtdb.reset_index()
        # delete the previous ome codes of those that need new ones
        need_ome_mtdb["ome"] = ["" for x in need_ome_mtdb["assembly_acc"]]
        # the reference database should now contain the original and fine
        # codenames
        ref_ome_mtdb = mtdb({**fine_ome_mtdb, **orig_mtdb}, index="ome")
        # generate new codenames

        part_ome_mtdb, failed = gen_omes(need_ome_mtdb, ref_ome_mtdb.reset_index())
        overlap_mtdb = add_mtdb.set_index("assembly_acc")

        # prepare a database of new genome codes and previously fine ones
        new_ome_mtdb = mtdb(
            {**part_ome_mtdb.set_index("ome"), **fine_ome_mtdb}, index="ome"
        )
        # choose the assembly accession column to reference for old names
        new_ome_mtdb = new_ome_mtdb.set_index("assembly_acc")
        old_ome2new_ome = {
            v["ome"]: new_ome_mtdb[k]["ome"]
            for k, v in overlap_mtdb.items()
            if v["ome"] in inter_omes
        }
        new_ome2old_ome = {v: k for k, v in old_ome2new_ome.items()}
        new_ome_mtdb = new_ome_mtdb.set_index("ome")
        for k, v in old_ome2new_ome.items():
            logger.debug(f"{k} converted to {v}")

        # create directories for new files
        fna_dir, gff_dir, faa_dir = (
            f"{update_path}fna/",
            f"{update_path}gff3/",
            f"{update_path}faa/",
        )
        for path_ in [fna_dir, gff_dir, faa_dir]:
            if not Path(path_).is_dir():
                Path(path_).mkdir()

        # convert the file header names to the new omes
        for ome, old_ome in new_ome2old_ome.items():
            row = new_ome_mtdb[ome]
            if row["assembly_acc"] == old_ome:
                new_ome_mtdb[ome]["assembly_acc"] = ome

            old_ome = new_ome2old_ome[ome]
            fna = fa2dict(row["fna"])
            new_fna = {k.replace(old_ome + "_", ome + "_"): v for k, v in fna.items()}
            with open(f"{fna_dir}{ome}.fna", "w") as out:
                out.write(dict2fa(new_fna))
            gff = gff2list(row["gff3"])
            for entry in gff:
                entry["seqid"] = entry["seqid"].replace(f"{old_ome}_", f"{ome}_")
                entry["attributes"] = entry["attributes"].replace(
                    f"{old_ome}_", f"{ome}_"
                )
            with open(f"{gff_dir}{ome}.gff3", "w") as out:
                out.write(list2gff(gff))
            faa = fa2dict(row["faa"])
            new_faa = {k.replace(old_ome + "_", ome + "_"): v for k, v in faa.items()}
            with open(f"{faa_dir}{ome}.faa", "w") as out:
                out.write(dict2fa(new_faa))
            row["fna"] = f"{fna_dir}{ome}.fna"
            row["faa"] = f"{faa_dir}{ome}.faa"
            row["gff3"] = f"{gff_dir}{ome}.gff3"

        # assembly accessions must be unique
        new_ome_mtdb = new_ome_mtdb.reset_index()

        return new_ome_mtdb
    else:
        return add_mtdb.reset_index()


def db2primary(addDB, refDB, save=False, combined=False):
    """Finalize an update by converting the updated MTDB into the primary
    MTDB"""
    if save:
        move_ns = shutil.copy
    else:
        move_ns = shutil.move

    addDB = addDB.reset_index()
    refDB = refDB.reset_index()

    refOmes = set(refDB["ome"])
    addOmes = set(addDB["ome"])
    base_ome2update_ome = {re.search(r"^[^\d]+\d+", x)[0]: x for x in refDB["ome"] if x}
    updates = {}
    refDB = refDB.set_index()
    if refOmes.intersection(addOmes) and not combined:
        logger.info(refOmes.intersection(addOmes))
        raise KeyError(
            "ERROR: ome codes exist in database. Rerun `mtdb predb` or remove manually"
        )
    for i, ome in enumerate(addDB["ome"]):
        base_ome = re.search(r"^[^\d]+\d+", ome)[0]
        if base_ome in base_ome2update_ome:
            update_ome = base_ome2update_ome[base_ome]
            updates[update_ome] = ome
            del refDB[update_ome]
        if Path(addDB["gff3"][i]).is_file():
            move_ns(addDB["gff3"][i], format_path("$MYCOGFF3/" + ome + ".gff3"))
        elif not Path(format_path("$MYCOGFF3/" + ome + ".gff3")).is_file():
            raise FileNotFoundError(f"{ome} missing gff3 for unknown reason")
        if Path(addDB["fna"][i]).is_file():
            move_ns(addDB["fna"][i], format_path("$MYCOFNA/" + ome + ".fna"))
        elif not Path(format_path("$MYCOFNA/" + ome + ".fna")).is_file():
            raise FileNotFoundError(f"{ome} missing fna for unknown reason")
        if Path(addDB["faa"][i]).is_file():
            move_ns(addDB["faa"][i], format_path("$MYCOFAA/" + ome + ".faa"))
        elif not Path(format_path("$MYCOFAA/" + ome + ".faa")).is_file():
            raise FileNotFoundError(f"{ome} missing faa for unknown reason")
        addDB["gff3"][i] = os.environ["MYCOGFF3"] + ome + ".gff3"
        addDB["fna"][i] = os.environ["MYCOFNA"] + ome + ".fna"
        addDB["faa"][i] = os.environ["MYCOFAA"] + ome + ".faa"
    addDB = addDB.set_index()
    for ome, row in addDB.items():
        refDB[ome] = row

    return refDB.reset_index(), updates


def control_flow(
    init,
    update,
    reference,
    add,
    taxonomy,
    predb,
    save,
    nonpublished,
    ncbi_only,
    lineage,
    rank,
    kingdom,
    failed,
    forbidden,
    resume,
    no_md5,
    cpu,
    ncbi_api=None,
    overwrite=True,
    chunk=100,
    tape_wait=None,
):

    abbr2king = {
        "a": "animals",
        "r": "archaea",
        "f": "fungi",
        "b": "bacteria",
        "p": "plants",
    }

    kingdom = kingdom.lower()
    if kingdom not in abbr2king:
        if kingdom not in set(abbr2king.values()):
            logger.error("invalid --kingdom")
            sys.exit(431)
    else:
        kingdom = abbr2king[kingdom]

    if not init and not update and not reference and not add and not taxonomy:
        logger.error("--update/--init/--reference/--add must be specified")
        sys.exit(15)
    elif reference and not init:
        logger.error("--reference requires a --init directory")
        sys.exit(14)
    elif lineage and not rank:
        logger.error("--lineage requires --rank")
        sys.exit(16)
    elif lineage and not init:
        logger.error("--lineage requires --init")
        sys.exit(17)
    elif predb and not init:
        logger.error("--predb requires --init")
        sys.exit(18)
    elif predb and lineage:
        logger.error("--predb and --lineage are incompatible")
        sys.exit(20)
    elif reference:
        if add:
            logger.error("--add and --reference are incompatible")
            sys.exit(13)
        elif predb:
            logger.error("--reference and --predb are incompatible")
            sys.exit(19)
        else:
            ref_db = mtdb(format_path(reference), add_paths=False)

    if predb:
        predb_path = format_path(predb)

    #    if rogue:
    rogue_bool = True
    if ncbi_only:
        jgi = False
    else:
        jgi = True

    # acquire the lineages inputted
    rank2lineages = {}
    permitted_ranks = {"phylum", "subphylum", "class", "order", "family", "genus"}
    if lineage:
        lineage_constraints = split_input(lineage)
        rank_constraints = split_input(rank)
        if len(lineage_constraints) != len(rank_constraints):
            logger.error("--lineage must be same length as --rank")
            sys.exit(18)
        for rank_c in rank_constraints:
            if rank_c.lower() not in permitted_ranks:
                logger.error(f"accepted ranks: {permitted_ranks}")
                sys.exit(22)
        rank2lineages = defaultdict(set)
        for i, v in enumerate(lineage_constraints):
            rank2lineages[rank_constraints[i]].add(v.lower())
        rank2lineages = {
            k.lower(): sorted(v)
            for k, v in sorted(rank2lineages.items(), key=lambda x: x[0])
        }

    # parse and check configuration nonpublished arguments
    config = {}
    if "MYCODB" in os.environ:
        config_path = format_path("$MYCODB/../config/mtdb.json")
        if Path(config_path).is_file():
            config = read_json(format_path(config_path))
            # for LEGACY installs:
            if "lineage_constraints" not in config:
                config["lineage_constraints"] = {}
                write_json(config, config_path)
        elif not init:
            logger.error("corrupted MycotoolsDB - no configuration found")
            sys.exit(21)
        if not init:  # is MYCODB initialized?
            #            rogue_bool = config['rogue']
            #                nonpublished = config['nonpublished']
            if bool(nonpublished) and not bool(config["nonpublished"]):
                config["nonpublished"] = validate_t_and_c(config, discrepancy=True)
                write_json(config, config_path)
            if bool(config["jgi"]) and bool(ncbi_only):  # and not overwrite:
                logger.error("--ncbi_only specified after initialization")
                sys.exit(173)
        elif init:
            if format_path(init) != format_path(os.environ["MYCODB"] + "../../"):
                logger.error("MTDB linked. Unlink via `mtdb -u`")
                sys.exit(175)

    # nonfungi is nonpublished by default because it is all GenBank
    if kingdom != "fungi":
        nonpublished = True
    # archaic placeholder for reference / rogue DB setup
    elif nonpublished and rogue_bool:
        nonpublished = validate_t_and_c(config)
    else:
        nonpublished = False

    #    branch = 'stable'
    db_path = primary_db()
    if not resume or add:
        date = datetime.now().strftime("%Y%m%d")
    else:
        date = str(resume)

    if not ncbi_api:
        ncbi_api, jgi_email, jgi_pwd = login_check()
    if ncbi_api:
        Entrez.api_key = ncbi_api

    if init:
        dbtype = kingdom
        init_dir = format_path(init)
        if Path(init_dir).is_dir():
            init_dir += "mycotoolsdb/"
        if not init_dir.endswith("/"):
            init_dir += "/"
        envs = {
            "MYCOFNA": init_dir + "data/fna",
            "MYCOFAA": init_dir + "data/faa",
            "MYCOGFF3": init_dir + "data/gff3",
            "MYCODB": init_dir + "mtdb/",
        }
        os.environ["MYCODB"] = init_dir + "mtdb/"
        output, config = init_db(
            init_dir,
            dbtype,
            envs,
            dbtype,
            date=date,
            rogue=rogue_bool,
            nonpublished=nonpublished,
            jgi=jgi,
            repo=format_path(reference),
            rank2lineages=rank2lineages,
        )
        for env in envs:
            os.environ[env] = envs[env]
        orig_db = db2df(mtdb())  # initialize a new database
        update_path = output + "log/" + date + "/"
        if not Path(update_path).is_dir():
            Path(update_path).mkdir()
        mtdb_initialize(
            init_dir, init=True
        )  # init_dir + 'config/mtdb.json', init = True)
    else:
        try:
            output = format_path("$MYCODB/..")
        except KeyError:
            logger.error("MTDB not linked. Link via `mtdb -i <DB_PATH>`")
            sys.exit(50)
        update_path = output + "log/" + date + "/"
        if not Path(update_path).is_dir():
            Path(update_path).mkdir()
        if not True:  # config['rogue']: # NEED TO MAKE THIS wget a particular URL
            old_db = db2df(db_path)
            shutil.move(db_path, update_path + Path(db_path).name)
            git_pull = subprocess.call(
                [
                    "git",
                    "pull",
                    "-C",
                    output + "mtdb",
                    config["repository"],
                    "-B",
                    branch,
                ],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            new_db = db2df(primary_db())
            orig_db = pd.concat([old_db, new_db.loc[~new_db["ome"].isin(old_db.index)]])
        else:
            orig_db = db2df(primary_db())

    orig_db = orig_db.dropna(subset=["ome"])

    if config["branch"] in {"prokaryote", "bacteria"}:
        jgi = False
        group = "prokaryotes"
        king = "bacteria"  # NEED to make DB tools pull from this
        rank = "superkingdom"
    elif config["branch"] in {"plants"}:
        jgi = False
        group = "eukaryotes"
        king = "viridiplantae"
        rank = "kingdom"
    #    elif config['branch'] in {'protists'}:
    #       jgi = False
    #      group = 'eukaryotes'
    #     king = 'protists'
    #    rank = 'kingdom'
    elif config["branch"] in {"animals"}:
        jgi = False
        group = "eukaryotes"
        king = "metazoa"
        rank = "kingdom"
    elif config["branch"] in {"archaea"}:
        jgi = False
        group = "prokaryotes"
        king = "Archaea"
        rank = "superkingdom"
    else:
        jgi = not ncbi_only
        group = "eukaryotes"
        king = "fungi"
        rank = "kingdom"

    if add or predb:  # add predb2mtdb 2 master database
        if predb:
            add_predb = read_predb(predb_path)
            addDB, init_failed = predb2mtdb(
                add_predb, orig_db, update_path, cpus=cpu, remove=False, spacer="\t\t"
            )
            if init_failed:
                if not failed:
                    logger.error("some genomes failed curation")
                    sys.exit(23)
                else:
                    logger.warning("some genomes failed curation")

        else:
            addDB = mtdb(format_path(add))
        # we need full Paths for an addDB
        gff_fail, fna_fail, faa_fail = False, False, False
        if not all(Path(format_path(x)).is_file() for x in addDB.reset_index()["gff3"]):
            logger.error("some GFF paths do not exist")
            gff_fail = [
                x
                for x in addDB.reset_index()["gff3"]
                if not Path(format_path(x)).is_file()
            ]
            logger.debug(",".join(gff_fail))
        if not all(Path(format_path(x)).is_file() for x in addDB.reset_index()["fna"]):
            logger.error("some FNA paths do not exist")
            fna_fail = [
                x
                for x in addDB.reset_index()["fna"]
                if not Path(format_path(x)).is_file()
            ]
            logger.debug(",".join(fna_fail))
        if not all(Path(format_path(x)).is_file() for x in addDB.reset_index()["faa"]):
            logger.error("some FAA paths do not exist")
            faa_fail = [
                x
                for x in addDB.reset_index()["faa"]
                if not Path(format_path(x)).is_file()
            ]
            logger.debug(",".join(faa_fail))
        if gff_fail or fna_fail or faa_fail:
            sys.exit(124)

        addDB["aquisition_date"] = [date for x in addDB["ome"]]
        # make date the acquisition time
        orig_mtdb = mtdb(primary_db())
        update_path = format_path("$MYCODB/../" + "log/" + date + "/")
        if not Path(update_path).is_dir():
            Path(update_path).mkdir()
        shutil.copy(primary_db(), update_path)

        tax_path = f"{update_path}../taxonomy.tsv"
        tax_dicts = read_prev_tax(tax_path)
        tax_dicts = gather_taxonomy(
            addDB,
            api_key=ncbi_api,
            king=king,
            rank=rank,
            tax_dicts=tax_dicts,
            output_path=tax_path,
        )
        addDB, genus_dicts = assimilate_tax(addDB, tax_dicts)
        addDB = check_add_mtdb(orig_mtdb, addDB, update_path, overwrite)

        write_forbid_omes(set(addDB["ome"]), format_path("$MYCODB/../log/relics.txt"))

        new_mtdb, update_omes = db2primary(addDB, orig_mtdb, save=True)
        return write_primary(new_mtdb, date, update_path)

    if taxonomy:
        new_db, update_mtdb = taxonomy_update(
            orig_db,
            update_path,
            date,
            config,
            ncbi_api,
            rank=rank,
            group=king,
        )
        write_primary(update_mtdb, date, update_path)
        sys.exit(0)
    elif reference:
        if any(not x for x in ref_db["published"]) and not nonpublished:
            logger.warning(
                "nonpublished data detected in reference and will be ignored"
            )

        new_mtdb, update_mtdb = ref_update(
            ref_db,
            update_path,
            date,
            failed,
            jgi_email,
            jgi_pwd,
            config,
            ncbi_api,
            cpus=cpu,
            check_MD5=not bool(no_md5),
            jgi=jgi,
            group=group,
            kingdom=king,
            remove=not save,
            taxonomy=True,
            chunk=chunk,
            tape_wait=tape_wait,
        )
    else:
        new_mtdb, update_mtdb = rogue_update(
            orig_db,
            update_path,
            date,
            failed,
            jgi_email,
            jgi_pwd,
            config,
            ncbi_api,
            cpus=cpu,
            check_MD5=not bool(no_md5),
            jgi=jgi,
            group=group,
            kingdom=king,
            remove=not save,
            lineage_constraints=config["lineage_constraints"],
            chunk=chunk,
            tape_wait=tape_wait,
        )

    if not update_mtdb:
        logger.info("No new data acquired")

    if not save:  # add the predb2mtdb and remove files
        #        df2db(db, format_path('$MYCODB/' + date + '.mtdb'))
        # output new database and new list of omes

        logger.info("Moving data into database")
        write_forbid_omes(
            set(new_mtdb["ome"]), format_path("$MYCODB/../log/relics.txt")
        )

        full_mtdb, update_omes = db2primary(
            update_mtdb, new_mtdb, save=False, combined=True
        )
        write_primary(full_mtdb, date, update_path)
        rm_raw_data(update_path)
        logger.info("MTDB update complete")
    #        gen_algn_db(
    #           update_path, set(full_mtdb['ome'])
    #          )
    else:
        # NEED to: insert note aboutrunning updatedb on predb
        new_mtdb.df2db(format_path(update_path + date + ".mtdb"))
        logger.info(
            f"Update ready for `mtdb u -a` at "
            + f'{format_path(update_path + date + ".mtdb")}'
        )
        # output new database and new list of omes

    return primary_db()


def main():

    abbr2king = {
        "a": "[a]nimals",  #'r': 'a[r]chaea',
        "f": "[f]ungi",
        "b": "[b]acteria",
        "p": "[p]lants",
    }
    parser = argparse.ArgumentParser(
        description="Initializes or updates "
        + "MycotoolsDB (MTDB) derived from a kingdom of interest. Animals "
        + "= metazoa, plants = viridiplantae"
    )

    init_args = parser.add_argument_group("MTDB Initializiation")
    init_args.add_argument("-i", "--init", help="Initialize MTDB in dir")
    init_args.add_argument(
        "-r",
        "--reference",
        help="[-i]: Initialize primary MTDB using a reference .mtdb",
    )
    init_args.add_argument(
        "-p",
        "--predb",
        help="[-i]: Initialize primary MTDB using a reference predb .tsv",
    )
    init_args.add_argument(
        "-k",
        "--kingdom",
        default="fungi",
        help="[-i]: Kingdom - " + str(sorted(abbr2king.values())) + "; DEFAULT: fungi",
    )

    upd_args = parser.add_argument_group("MTDB Updating")
    upd_args.add_argument("-u", "--update", action="store_true")
    upd_args.add_argument(
        "-a", "--add", help=".mtdb with full paths to add to database"
    )
    upd_args.add_argument(
        "-t",
        "--taxonomy",
        action="store_true",
        help="Remove old taxonomy metadata, update taxonomy, and exit",
    )
    upd_args.add_argument(
        "--keep",
        action="store_true",
        help="[-a] Keep original MTDB data when adding overlapping accessions",
    )
    upd_args.add_argument(
        "--save",
        action="store_true",
        help="[-u] Do not integrate/delete new data; -a to complete",
    )

    #    init_args.add_argument('--reinit', action = 'store_true', help = 'Redownload all web data')
    #    parser.add_argument('--rogue', action = 'store_true',
    #       help = 'De novo MTDB') # currently required

    conf_args = parser.add_argument_group("Configuration")
    conf_args.add_argument(
        "--nonpublished",
        action="store_true",
        help="[FUNGI]: Include MycoCosm restricted-use",
    )
    conf_args.add_argument(
        "--ncbi_only", help="[FUNGI, -i]: Forego MycoCosm", action="store_true"
    )
    conf_args.add_argument(
        "-l", "--lineage", help="[-i, -r]: Lineage(s) to initialize MTDB with"
    )
    conf_args.add_argument(
        "-rk",
        "--rank",
        help="[-i, -l]: Rank(s) that positionally " + "correspond to -l",
    )
    conf_args.add_argument("--failed", action="store_true", help="Rerun/ignore failed")
    conf_args.add_argument("--forbidden", action="store_true", help="Rerun forbidden")

    #    conf_args.add_argument('--deviate', action = 'store_true', help = 'Deviate' \
    #       + ' from existing config without prompting')

    run_args = parser.add_argument_group("Runtime")
    run_args.add_argument("--resume", type=int, help="Resume previous date (YYYYmmdd)")
    run_args.add_argument(
        "--no_md5",
        action="store_true",
        help="Skip NCBI MD5" + " (expedite large reruns)",
    )
    run_args.add_argument(
        "--chunk",
        type=int,
        default=100,
        help="Accessions to download per datasets call; DEFAULT: 100",
    )
    run_args.add_argument(
        "--tape_wait",
        type=int,
        default=None,
        help="[FUNGI]: Maximum minutes to wait for a MycoCosm genome's tape "
        + "restore before deferring it to a later run; DEFAULT: wait indefinitely",
    )
    run_args.add_argument("-c", "--cpu", type=int, default=1)
    args = parser.parse_args()
    setup_logging(verbose=getattr(args, "verbose", False))

    args_dict = {
        "Primary MTDB": primary_db(verbose=False),
        "Update": args.update,
        "Initialize": args.init,
        "Add": format_path(args.add),  #'Rogue': rogue_bool,
        "Include Restricted": bool(args.nonpublished),
        "Resume": args.resume,
        "Retry failed": args.failed,
        "Retry forbidden": args.forbidden,
        "Save raw data": args.save,
        "Chunk": args.chunk,
        "Max tape wait": (
            "indefinite" if args.tape_wait is None else f"{args.tape_wait} minute(s)"
        ),
    }

    find_execs(["datasets"], exit={"datasets"})
    start_time = intro("Update MycotoolsDB", args_dict)

    control_flow(
        args.init,
        args.update,
        args.reference,
        args.add,
        args.taxonomy,
        args.predb,
        args.save,
        args.nonpublished,
        args.ncbi_only,
        args.lineage,
        args.rank,
        args.kingdom,
        args.failed,
        args.forbidden,
        args.resume,
        args.no_md5,
        args.cpu,
        overwrite=not args.keep,
        chunk=args.chunk,
        tape_wait=args.tape_wait,
    )

    outro(start_time)


def cli():
    main()


if __name__ == "__main__":
    cli()
