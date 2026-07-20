#! /usr/bin/env python3

# NEED a db check to ensure the log is relevant to the input
# NEED to convert to datasets
# NEED to consider refseq genomes with annotations when genbank doesn't have them

import logging
import os
import re
import sys
import time
import shutil
import urllib
import urllib.request
import requests
import argparse
import subprocess
import pandas as pd
from tqdm import tqdm
from Bio import Entrez
from datetime import datetime
from mycotools.lib.kontools import (
    intro,
    outro,
    format_path,
    prep_output,
    findExecs,
    setup_logging,
)
from mycotools.lib.dbtools import log_editor, loginCheck, mtdb
from pathlib import Path

logger = logging.getLogger(__name__)


def ncbidb2df(data, stdin=False):
    import pandas as pd

    if isinstance(data, mtdb):
        db_df = pd.DataFrame(data.reset_index())
    elif not stdin:
        data = format_path(data)
        db_df = pd.read_csv(data, sep="\t")
        if "ome" not in set(db_df.columns) and "assembly_acc" not in set(db_df.columns):
            db_df = pd.read_csv(data, sep="\t", header=None)
    else:
        db_df = pd.read_csv(StringIO(data), sep="\t")
        if "ome" not in set(db_df.columns) and "assembly_acc" not in set(db_df.columns):
            db_df = pd.read_csv(StringIO(data), sep="\t", header=None)

    db_df = db_df.fillna("")

    return db_df


def prepare_folders(output_path, gff, prot, assem, transcript):

    file_types = []
    if assem:
        if not Path(output_path + "fna").exists():
            Path(output_path + "fna").mkdir()
        file_types.append("fna")
    if gff:
        if not Path(output_path + "gff3").exists():
            Path(output_path + "gff3").mkdir()
        file_types.append("gff3")
    if prot:
        if not Path(output_path + "faa").exists():
            Path(output_path + "faa").mkdir()
        file_types.append("faa")
    if transcript:
        if not Path(output_path + "transcript").exists():
            Path(output_path + "transcript").mkdir()
        file_types.append("transcript")

    return file_types


def compile_log(output_path, remove=False):

    acc2log = {}
    if not Path(output_path + "ncbiDwnld.fallback.log").is_file():
        with open(output_path + "ncbiDwnld.fallback.log", "w") as out:
            out.write(
                "#ome\tassembly_acc\tassembly\tproteome\tgff3\ttranscript\t"
                + "fna_md5\tfaa_md5\tgff3_md5\ttrans_md5\tgenome_id(s)\tgenus\tspecies\tstrain"
            )

        # too risky, too many things can go wrong and then users would be in a
        # loop, but necessary for huge downloads
    else:
        with open(output_path + "ncbiDwnld.fallback.log", "r") as raw:
            for line in raw:
                if not line.startswith("#"):
                    data = [x.rstrip() for x in line.split("\t")]
                    while len(data) < 13:
                        data.append("")
                    acc2log[data[0]] = {
                        "assembly_acc": str(data[1]),
                        "fna": str(data[2]),
                        "faa": str(data[3]),
                        "gff3": str(data[4]),
                        "transcript": str(data[5]),
                        "fna_md5": str(data[6]),
                        "faa_md5": str(data[7]),
                        "gff3_md5": str(data[8]),
                        "transcript_md5": str(data[9]),
                        "genome_id": data[10],
                        "genus": str(data[11]),
                        "species": str(data[12]),
                        "strain": str(data[13]),
                    }

    return acc2log


def wait_for_ncbi(count, api=False):
    if count >= 2:
        if not api:
            time.sleep(1)
            count = 0
        elif count >= 7:
            time.sleep(1)
            count = 0
    return count


def esearch_ncbi(accession, column, database="assembly"):
    search_term, esc_count = accession + "[" + column + "]", 0
    while esc_count < 3:
        try:
            handle = Entrez.esearch(db=database, term=search_term)
            genome_ids = Entrez.read(handle)["IdList"]
            break
        except (RuntimeError, urllib.error.HTTPError) as e:
            time.sleep(1)
            esc_count += 1
    else:
        logger.error("%s %s %s", "\tERROR:", accession, "failed to search NCBI")
        return None
    return genome_ids


def esummary_ncbi(ID, database):

    esc_count = 0
    while esc_count < 10:
        esc_count += 1
        try:
            handle = Entrez.esummary(db=database, id=ID, report="full")
            record = Entrez.read(handle, validate=False)
        except urllib.error.HTTPError:
            time.sleep(0.1)
            continue
        if database == "assembly":
            try:  # is it populated with an FTP?
                ftp_path = str(
                    record["DocumentSummarySet"]["DocumentSummary"][0][
                        "FtpPath_GenBank"
                    ]
                )
            except IndexError:  # wait a sec and retry
                time.sleep(1)
                continue
        break
    else:  # too many failed attempts
        if esc_count >= 10:
            raise urllib.error.HTTPError("\tERROR: FTP request failed")

    return record


# collects paths to download proteomes and assemblies
def collect_ftps(
    ncbi_df,
    acc2log,
    api_key=0,
    column="assembly_acc",
    ncbi_column="Assembly Accession",
    database="assembly",
    output_path="",
    verbose=True,
    remove=False,
    spacer="\t\t",
):

    count, failed = 0, []

    # for each row in the assembly, grab the accession number, form the search term for Entrez, use Entrez,
    if ncbi_column in {"assembly", "genome", "uid"}:
        out_df = ncbi_df[ncbi_df.index.isin(set(acc2log.keys()))]
        ncbi_df = ncbi_df[~ncbi_df.index.isin(set(acc2log.keys()))]
    for accession, row in tqdm(ncbi_df.iterrows(), total=len(ncbi_df)):
        if accession in acc2log:  # add all rows that have indices associated with this
            # query type
            out_df = pd.concat([out_df, row.to_frame().T])
            icount = 1
            test = str(accession) + "_" + str(icount)
            while test in acc2log:
                count += 1
                #                row['assembly_acc'] = acc2log[test]['genome_id']
                if "ome" in row.keys():
                    row["ome"] = None  # haven't assigned a mycotools ID yet
                out_df = pd.concat([out_df, row.to_frame().T])
                sys.exit()
                test = str(accession) + "_" + str(icount)
            continue

        elif pd.isnull(row[column]) or not row[column]:  # ignore blank entries
            acc2log[str(accession)] = {
                "assembly_acc": accession,
                "fna": "",
                "faa": "",
                "gff3": "",
                "transcript": "",
                "fna_md5": "",
                "faa_md5": "",
                "gff3_md5": "",
                "transcript_md5": "",
                "genome_id": "",
                "genus": "",
                "species": "",
                "strain": "",
            }
            failed.append([accession, datetime.strftime(row["version"], "%Y%m%d")])
            continue

        if ncbi_column not in {"uid"}:  # we already have the uid, no worries
            genome_id = esearch_ncbi(accession, ncbi_column, database="assembly")
        else:
            genome_id = [accession]

        if not genome_id:  # No IDs retrieved
            if "ome" in row.keys():
                accession = row["ome"]
            logger.info(spacer + "" + accession + " failed to find genome ID")
            try:
                failed.append([accession, datetime.strftime(row["version"], "%Y%m%d")])
            except TypeError:  # if the row can't be formatted as a date
                failed.append([accession, row["version"]])
            continue

        if ncbi_column in {
            "Assembly Accession",
            "assembly",
            "genome",
            "uid",
        }:  # be confident it is the most
            # recent assembly UID
            genome_id = [str(max([int(i) for i in genome_id]))]

        icount = 0
        for ID in genome_id:
            if icount:
                new_acc = str(accession) + "$" + str(icount)
            else:
                new_acc = accession
            # obtain the path from a summary of the ftp directory and create the standard paths for proteomes and assemblies
            acc2log[str(new_acc)] = {
                "assembly_acc": accession,
                "fna": "",
                "faa": "",
                "gff3": "",
                "transcript": "",
                "fna_md5": "",
                "faa_md5": "",
                "gff3_md5": "",
                "transcript_md5": "",
                "genome_id": ID,
                "genus": "",
                "species": "",
                "strain": "",
            }
            record = esummary_ncbi(ID, database)

            strain = ""  # populate a fallback strain
            record_info = record["DocumentSummarySet"]["DocumentSummary"][0]
            org = record_info["SpeciesName"].split()
            genus = org[0]
            if len(org) > 2:
                species = org[1]
                strain1 = "".join(org[2:])
            elif len(org) == 2:
                species = org[1]
                strain1 = ""
            else:
                species = "sp."
                strain1 = ""

            try:
                for attr in record_info["Biosource"]["InfraspeciesList"]:
                    if attr["Sub_type"].lower() == "strain":
                        strain = attr["Sub_value"]
                        if strain.lower() in {"missing", "none"}:
                            strain = ""
                        break
            except KeyError:
                if strain1:
                    strain = strain1

            strain = re.sub(r"[^a-zA-Z0-9]", "", strain)

            ftp_path = str(record_info["FtpPath_GenBank"])

            if not ftp_path:
                logger.info(spacer + "" + new_acc + " failed to return any FTP path")
                try:
                    failed.append(
                        [accession, datetime.strftime(row["version"], "%Y%m%d")]
                    )
                except TypeError:
                    failed.append([accession, str(row["version"])])
                continue

            ass_md5, gff_md5, trans_md5, prot_md5, md5s = "", "", "", "", {}
            basename = Path(ftp_path).name

            dwnld = 0
            for attempt in range(3):
                try:
                    r = requests.head(
                        ftp_path.replace("ftp://", "https://"), allow_redirects=True
                    )
                    break
                except:
                    time.sleep(1)

            if r.status_code != 200:
                dwnld = -1
            else:
                md5_path = ftp_path.replace("ftp://", "https://") + "/md5checksums.txt"

                dwnld = subprocess.call(
                    [
                        "curl",
                        md5_path,
                        "-o",
                        output_path + ".tmpmd5",
                        "--connect-timeout",
                        "5",
                    ],
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                )

                count += 1
            if dwnld == 0:
                with open(output_path + ".tmpmd5", "r") as raw:
                    for line in raw:
                        data = line.rstrip().split("  ")
                        #                        data = line.rstrip().split()
                        if data and len(data) == 2:
                            try:
                                md5s[ftp_path + "/" + Path(data[1]).name] = data[0]
                            except IndexError:  # 404 error or something else
                                md5s = {}
                                break
            else:
                md5s = {}

            tranname = Path(ftp_path.replace("/GCA", "/GCF")).name
            #          tranname = os.path.basename(ftp_path)
            assembly = ftp_path + "/" + basename + "_genomic.fna.gz"
            if assembly in md5s:
                ass_md5 = md5s[assembly]
            else:
                assembly = ""
            proteome = ftp_path + "/" + basename + "_protein.faa.gz"
            if proteome in md5s:
                prot_md5 = md5s[proteome]
            else:
                proteome = ""
            gff3 = ftp_path + "/" + basename + "_genomic.gff3.gz"
            test_gff3 = re.sub(r"\.gff3\.gz$", ".gff.gz", gff3)
            if gff3 in md5s:
                gff_md5 = md5s[gff3]
            elif test_gff3 in md5s:
                gff3 = test_gff3
                gff_md5 = md5s[gff3]
            else:
                gff3 = ""

            transcript = (
                ftp_path.replace("/GCA", "/GCF") + "/" + tranname + "_rna.fna.gz"
            )
            #           transcript = f'{ftp_path}/{tranname}_rna.fna.gz'
            if transcript in md5s:
                trans_md5 = md5s[transcript]

            if (not assembly or not gff3) and remove:
                try:
                    failed.append(
                        [accession, datetime.strftime(row["version"], "%Y%m%d")]
                    )
                except TypeError:
                    failed.append(
                        [accession, datetime.strftime(datetime.now(), "%Y%m%d")]
                    )

            log_editor(
                output_path + "ncbiDwnld.fallback.log",
                str(new_acc),
                str(accession)
                + "\t"
                + accession
                + "\t"
                + assembly
                + "\t"
                + proteome
                + "\t"
                + gff3
                + "\t"
                + transcript
                + "\t"
                + ass_md5
                + "\t"
                + prot_md5
                + "\t"
                + gff_md5
                + "\t"
                + trans_md5
                + "\t"
                + ID
                + f"\t{genus}\t{species}\t{strain}",
            )
            acc2log[str(new_acc)] = {
                "assembly_acc": accession,
                "fna": assembly,
                "fna_md5": ass_md5,
                "faa": proteome,
                "faa_md5": prot_md5,
                "gff3": gff3,
                "gff3_md5": gff_md5,
                "transcript": transcript,
                "transcript_md5": trans_md5,
                "genome_id": ID,
                "genus": genus,
                "species": species,
                "strain": strain,
            }
            row["dwnld_id"] = ID
            if icount:
                if "ome" in row:
                    row["ome"] = None
            out_df = pd.concat([out_df, row.to_frame().T])
            icount += 1

            # if no API key is used, we can only generate 3 queries per second, otherwise we can use 10
            count = wait_for_ncbi(count, api_key)

    return acc2log, failed, out_df


# download the file depending on the type inputted
def download_files(
    acc_prots, acc, file_types, output_dir, count, remove=False, spacer="\t\t"
):

    dwnlds = {}
    for file_type in file_types:
        ftp_link = acc_prots[file_type]
        dwnlds[file_type] = -1
        if file_type == "fna":
            file_path = output_dir + "fna/" + Path(acc_prots[file_type]).name
        elif file_type == "gff3":
            file_path = output_dir + "gff3/" + Path(acc_prots[file_type]).name
        elif file_type == "faa":
            file_path = output_dir + "faa/" + Path(acc_prots[file_type]).name
        elif file_type == "transcript":
            file_path = output_dir + "transcript/" + Path(acc_prots[file_type]).name

        if Path(file_path).is_file():
            count += 1
            md5_cmd = subprocess.run(
                ["md5sum", file_path], stdout=subprocess.PIPE
            )  # check the md5
            md5_res = md5_cmd.stdout.decode("utf-8")
            md5_find = re.search(r"\w+", md5_res)
            md5 = md5_find[0]

            if md5 == acc_prots[file_type + "_md5"]:
                logger.info(f"{spacer}\t{file_type}: {Path(file_path).name}")
                dwnlds[file_type] = 0
                continue
        elif Path(file_path[:-3]).is_file():
            logger.info(f"{spacer}\t{file_type}: {Path(file_path).name}")
            dwnlds[file_type] = 0
            continue

        if ftp_link == "":
            dwnlds[file_type] = 15
            continue
        #       esc_count = 0
        for esc_count in range(3):
            count += 1
            dwnld = subprocess.call(
                ["curl", ftp_link, "-o", file_path + ".tmp", "--connect-timeout", "5"],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            if not dwnld:
                Path(file_path + ".tmp").rename(file_path)
                break
            else:
                time.sleep(1)
                count = 0

        if dwnld:
            logger.error(f"{spacer}\t\tERROR: {file_type} failed")
            dwnlds[file_type] = 69
            acc_prots[file_type] = ""
            log_editor(
                output_dir + "ncbiDwnld.fallback.log",
                str(acc),
                str(acc)
                + "\t"
                + str(acc_prots["assembly_acc"])
                + "\t"
                + str(acc_prots["fna"])
                + "\t"
                + str(acc_prots["faa"])
                + "\t"
                + str(acc_prots["gff3"])
                + "\t"
                + str(acc_prots["transcript"])
                + "\t"
                + str(acc_prots["fna_md5"])
                + "\t"
                + str(acc_prots["faa_md5"])
                + "\t"
                + str(acc_prots["gff3_md5"])
                + "\t"
                + str(acc_prots["transcript_md5"])
                + "\t"
                + f"{acc_prots['genome_id']}\t{acc_prots['genus']}\t"
                + f"{acc_prots['species']}\t{acc_prots['strain']}",
            )
            if remove and file_type in {"fna", "gff3"}:
                break
            continue

        if not Path(file_path).is_file():
            dwnlds[file_type] = 1
            logger.error(f"{spacer}\t\tERROR: {file_type} missing")
            if remove and file_type in {"fna", "gff3"}:
                break
        else:
            dwnlds[file_type] = 0
            if Path(file_path).stat().st_size < 150:
                logger.error(f"{spacer}\t{file_type}: ERROR, file too small")
                dwnlds[file_type] = 420
                if remove and file_type in {"fna", "gff3"}:
                    break
        logger.info(f"{spacer}\t{file_type}: {Path(file_path).name}")

    return dwnlds, count


def dwnld_mngr(ncbi_df, data, acc, file_types, output_path, count, remove, api, spacer):
    fail = []
    exits, count = download_files(
        data, acc, file_types, output_path, count, remove=remove, spacer=spacer
    )
    count = wait_for_ncbi(count, api)
    t_acc = acc
    try:
        if exits["fna"] != 0:
            if "$" in acc:
                t_acc = acc[: acc.find("$")]
            fail = [t_acc, ncbi_df["version"][t_acc]]
            return fail, count
    except KeyError:
        pass
    try:
        if exits["gff3"] != 0:
            if "$" in acc:  # is this correct?
                t_acc = acc[: acc.find("$")]
            fail = [t_acc, ncbi_df["version"][t_acc]]
            return fail, count
    except KeyError:
        pass
    return fail, count


def dwnld_mngr_no_MD5(
    ncbi_df, data, acc, file_types, output_path, count, remove, api, spacer
):
    run, fail = False, []
    for file_type in file_types:
        file_path = output_path + file_type + "/" + Path(data[file_type]).name
        if not Path(file_path).is_file():
            run = True
            break

    if run:
        exits, count = download_files(
            data, acc, file_types, output_path, count, remove=remove, spacer=spacer
        )
        count = wait_for_ncbi(count, api)
    else:
        exits = {x: 0 for x in file_types}
        # we aren't checking md5s, so assume the exit is 0
        return fail, count
    try:
        if exits["fna"] != 0:
            if "$" in acc:
                t_acc = acc[: acc.find("$")]
            else:
                t_acc = acc
            fail = [t_acc, ncbi_df["version"][t_acc]]
            return fail, count
    except KeyError:
        pass
    try:
        if exits["gff3"] != 0:
            if "$" in acc:  # dont know if this should be here
                t_acc = acc[: acc.find("$")]
            else:
                t_acc = acc

            fail = [t_acc, ncbi_df["version"][t_acc]]
            return fail, count
    except KeyError:
        pass
    return fail, count


def main(
    api=None,
    assembly=True,
    proteome=False,
    gff3=True,
    transcript=False,
    ncbi_df=False,
    remove=False,
    output_path=str(Path.cwd()),
    verbose=False,
    column="assembly_acc",
    ncbi_column="Assembly",
    check_MD5=True,
    spacer="\t\t",
):

    # initialize run directory and information
    output_path = format_path(output_path)
    file_types = prepare_folders(output_path, gff3, proteome, assembly, transcript)
    acc2log = compile_log(output_path, remove)

    # check if ncbi_df is a dataframe, and import if not
    if not isinstance(ncbi_df, pd.DataFrame) and Path(ncbi_df).is_file():
        ncbi_df = ncbidb2df(ncbi_df)
    if len(ncbi_df.index) == 0:
        ncbi_df = pd.DataFrame({i: [v] for i, v in enumerate(list(ncbi_df.keys()))})

    # make the modify date from a standard NCBI table the version if it does
    # not otherwise exist, else there isn't a version to reference
    if "Modify Date" in ncbi_df.keys() and not "version" in ncbi_df.keys():
        ncbi_df["version"] = pd.to_datetime(ncbi_df["Modify Date"])
    elif "version" not in ncbi_df.keys():
        ncbi_df["version"] = ""

    ncbi_df = ncbi_df.set_index(pd.Index(list(ncbi_df[column])))
    # preserve the original column, but index ncbi_df on it as well
    logger.debug("" + spacer + "Assembling NCBI ftp directories")
    acc2log, failed, ncbi_df = collect_ftps(
        ncbi_df,
        acc2log,
        remove=remove,
        ncbi_column=ncbi_column,
        column=column,
        api_key=api,
        output_path=output_path,
        verbose=verbose,
        spacer=spacer,
    )

    if remove:
        acc2log = {
            o: acc2log[o]
            for o in acc2log
            if all(acc2log[o][p] for p in ["fna", "gff3"])
        }
    new_df = pd.DataFrame()

    logger.debug(f"{spacer}Downloading {len(acc2log)} NCBI files")
    count = 0
    if "strain" not in ncbi_df.columns:
        ncbi_df["strain"] = ""
    if check_MD5:
        for acc, data in acc2log.items():
            logger.info(spacer + "" + str(acc))
            if data:
                fail, count = dwnld_mngr(
                    ncbi_df,
                    data,
                    acc,
                    file_types,
                    output_path,
                    count,
                    remove,
                    api,
                    spacer,
                )
                if fail:
                    failed.append(fail)
                else:
                    ncbi_df.at[acc, "assemblyPath"] = (
                        output_path + "fna/" + Path(acc2log[acc]["fna"]).name
                    )
                    ncbi_df.at[acc, "faa"] = (
                        output_path + "faa/" + Path(acc2log[acc]["faa"]).name
                    )
                    ncbi_df.at[acc, "gffPath"] = (
                        output_path + "gff3/" + Path(acc2log[acc]["gff3"]).name
                    )
                    ncbi_df.at[acc, "genus"] = acc2log[acc]["genus"]
                    ncbi_df.at[acc, "species"] = acc2log[acc]["species"]
                    if not ncbi_df.loc[acc, "strain"] or pd.isnull(
                        ncbi_df.loc[acc, "strain"]
                    ):
                        ncbi_df.at[acc, "strain"] = acc2log[acc]["strain"]
                    new_df = pd.concat([new_df, ncbi_df.loc[acc].to_frame().T])
            else:
                check = ncbi_df[ncbi_df[column] == acc[: acc.find("$")]]
                # check for entries in the inputted table that match the accession
                # provided without version modification
                db_vers = datetime.strftime(row["version"], "%Y%m%d")
                failed.append([acc, db_vers])
    else:  # there is no checking md5, this is for efficient, so make it
        # efficient by avoiding conditional expressions
        for acc, data in acc2log.items():
            logger.info(spacer + "" + str(acc))
            fail, count = dwnld_mngr_no_MD5(
                ncbi_df, data, acc, file_types, output_path, count, remove, api, spacer
            )
            if fail:
                failed.append(fail)
            else:
                for file_type in file_types:
                    ncbi_df.at[acc, file_type] = (
                        output_path + file_type + "/" + Path(data[file_type]).name
                    )
                new_df = pd.concat([new_df, ncbi_df.loc[acc].to_frame().T])

    if "fna" in new_df.keys():
        new_df = new_df.rename(columns={"fna": "assemblyPath"})
    if "gff3" in new_df.keys():
        new_df = new_df.rename(columns={"gff3": "gffPath"})
    new_df = new_df.reset_index()
    return new_df, failed


def get_SRA(assembly_acc, fastqdump="fastq-dump", pe=True):

    handle = Entrez.esearch(db="SRA", term=assembly_acc)
    ids = Entrez.read(handle)["IdList"]
    for id in ids:
        handle = Entrez.esummary(db="SRA", id=id, report="full")
        records = Entrez.read(handle, validate=False)
        for record in records:
            srr = re.search(r'Run acc="(S\w+\d+)"', record["Runs"])[1]
            logger.debug("" + srr)
            cmd, count = 1, 0
            if pe:
                while cmd and count < 3:
                    count += 1
                    cmd = subprocess.call(
                        ["prefetch", srr, "--max-size", "10t"], stdout=subprocess.PIPE
                    )
                    if cmd:
                        continue
                    cmd = subprocess.call(["vdb-validate", srr], stdout=subprocess.PIPE)
                    if cmd:
                        continue
                    cmd = subprocess.call(
                        [fastqdump, "--split-3", "--gzip", srr], stdout=subprocess.PIPE
                    )
                    if Path(srr + "_1.fastq").is_file():
                        #                  if os.path.isfile(srr + '_1.fastq'):
                        #                        cmd = subprocess.call(['gzip', f'{srr}_1.fastq'])
                        #                       cmd = subprocess.call(['gzip', f'{srr}_2.fastq'])
                        shutil.move(
                            srr + "_1.fastq.gz", assembly_acc + "_" + srr + "_1.fq.gz"
                        )
                        shutil.move(
                            srr + "_2.fastq.gz", assembly_acc + "_" + srr + "_2.fq.gz"
                        )
                    else:
                        #                        cmd = subprocess.call(['gzip', f'{srr}.fastq'])
                        logger.warning("file failed or not paired-end")
            else:
                while cmd and count < 3:
                    count += 1
                    cmd = subprocess.call(["prefetch", srr], stdout=subprocess.PIPE)
                    if cmd:
                        continue
                    cmd = subprocess.call(["vdb-validate", srr], stdout=subprocess.PIPE)
                    if cmd:
                        continue
                    cmd = subprocess.call(
                        [fastqdump, srr, "--gzip"], stdout=subprocess.PIPE
                    )
                    if cmd:
                        continue
                    #                    cmd = subprocess.call(['gzip', f'{srr}.fastq'])
                    if Path(srr + ".fastq.gz").is_file():
                        shutil.move(
                            srr + ".fastq.gz", assembly_acc + "_" + srr + ".fq.gz"
                        )
                    else:
                        logger.error("file failed")


def goSRA(df, output=str(Path.cwd()) + "/", pe=True):

    print()
    sra_dir = output + "sra/"
    if not Path(sra_dir).is_dir():
        Path(sra_dir).mkdir()
    os.chdir(sra_dir)
    fastqdump = findExecs("fastq-dump", exit=set("fastq-dump"))
    count = 0

    if "sra" in df.keys():
        row_key = "sra"
    else:
        row_key = "assembly_acc"

    for i, row in df.iterrows():
        logger.debug("" + row[row_key])
        get_SRA(row[row_key], fastqdump[0])
        count += 1
        if count >= 10:
            time.sleep(1)
            count = 0


def cli():
    parser = argparse.ArgumentParser(
        description="GenBank downloading utility. Downloads "
        + "accession by accession, files without MD5s are excluded"
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Space delimited accession; tab delimited file with -c",
    )
    parser.add_argument("-a", "--assembly", action="store_true")
    parser.add_argument("-p", "--proteome", action="store_true")
    parser.add_argument("-g", "--gff3", action="store_true")
    parser.add_argument("-t", "--transcript", action="store_true")
    parser.add_argument("-s", "--sra", action="store_true", help="Download SRAs only")
    parser.add_argument(
        "-pe",
        "--paired",
        action="store_true",
        help="Download paired-end SRAs. (REQUIRES -s)",
    )
    parser.add_argument(
        "-c", "--column", help='Accession column num/name; DEFAULT ["assembly_acc" | 0]'
    )
    parser.add_argument(
        "-n",
        "--ncbi_column",
        help="NCBI database associated with column. "
        + '{"assembly", "biosample", "bioproject", "genome" ...}; '
        + "DEFAULT: attempt to decipher",
    )
    parser.add_argument("-o", "--output", help="Output directory")
    parser.add_argument("-e", "--email", help="NCBI email")
    parser.add_argument("--api", help="NCBI API key for high query rate")
    args = parser.parse_args()
    setup_logging(verbose=getattr(args, "verbose", False))

    if args.email:
        ncbi_email = args.email
        Entrez.email = ncbi_email
        if args.api:
            ncbi_api = args.api
            Entrez.api_key = ncbi_api
        else:
            ncbi_api = None
    else:
        ncbi_email, ncbi_api, jgi_email, jgi_pwd = loginCheck(jgi=False)
        Entrez.email = ncbi_email
        if ncbi_api:
            Entrez.api_key = ncbi_api

    if not args.output:
        output = str(Path.cwd()) + "/"
    else:
        output = format_path(args.output)

    args_dict = {
        "NCBI Table": args.input,
        "email": ncbi_email,
        "Assemblies": args.assembly,
        "Proteomes": args.proteome,
        ".gff3's": args.gff3,
        "Transcripts": args.transcript,
        "SRA": args.sra,
    }

    start_time = intro("Download NCBI files", args_dict)
    #    if not args.assembly and not args.proteome and not args.gff3 and not args.sra and not args.transcript:
    #       eprint('\nERROR: You must choose at least one download option\nExit code 37', flush = True)
    #      sys.exit( 37 )

    if args.sra:
        if Path(format_path(args.input)).is_file():
            goSRA(
                pd.read_csv(format_path(args.input), sep="\t"), output, pe=args.paired
            )
        else:
            goSRA(pd.DataFrame({"sra": [args.input.rstrip()]}), output, pe=args.paired)
    else:
        if Path(format_path(args.input)).is_file():
            ncbi_df = pd.read_csv(args.input, sep="\t", header=None)
            if not args.column:
                if "assembly_acc" in ncbi_df.keys():
                    column = "assembly_acc"
                    ncbi_column = "Assembly Accession"
                elif "Assembly Accession" in ncbi_df.keys():
                    column = "assembly_acc"
                    ncbi_column = "Assembly Accession"
                else:
                    column = 0
            else:
                try:
                    column = ncbi_df.columns[int(0)]
                    ncbi_column = column
                except ValueError:  # not an integer
                    pass
            if not args.ncbi_column:
                if args.column is not None:
                    if column.lower() in {"assembly"}:
                        ncbi_column = "assembly"
                    elif column.lower() in {
                        "genome",
                        "assembly accession",
                        "assembly_acc",
                    }:
                        ncbi_column = "genome"
                    elif column.lower() in {"biosample", "biosample accession"}:
                        ncbi_column = "biosample"
                    else:
                        ncbi_column = column.lower()
                else:
                    ncbi_column = "genome"
            else:
                ncbi_column = args.ncbi_column.lower()
        else:
            ncbi_df = pd.DataFrame(
                {"assembly_acc": args.input.replace('"', "").replace("'", "").split()}
            )
            column = "assembly_acc"
            ncbi_column = "assembly"

        ncbi_df = ncbi_df.drop_duplicates(column)

        new_df, failed = main(
            assembly=args.assembly,
            column=column,
            ncbi_column=ncbi_column,
            proteome=args.proteome,
            gff3=args.gff3,
            transcript=args.transcript,
            ncbi_df=ncbi_df,
            output_path=output,
            verbose=True,
            spacer="",
        )
        new_df = new_df.rename(columns={"index": "#assembly_accession"})
        new_df["source"] = "ncbi"
        new_df["useRestriction (yes/no)"] = "no"
        #        if 'index' in new_df.columns:
        #           del new_df['index']
        if 0 in new_df.columns:
            del new_df[0]

        new_df.to_csv(args.input + ".predb", sep="\t", index=None)

    outro(start_time)


if __name__ == "__main__":
    cli()
