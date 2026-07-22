#! /usr/bin/env python3

# NEED a db check to ensure the log is relevant to the input
# NEED to consider refseq genomes with annotations when genbank doesn't have them

import os
import re
import sys
import json
import math
import time
import shutil
import urllib
import logging
import zipfile
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
    mk_output,
    find_execs,
    read_json,
    split_input,
    setup_logging,
)
from mycotools.lib.dbtools import clean_api_key, log_editor, login_check, mtdb
from pathlib import Path

logger = logging.getLogger(__name__)

pd.options.mode.chained_assignment = None


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


def compile_log(output_path):

    acc2log = {}
    if not Path(output_path).is_file():
        with open(output_path, "w") as out:
            out.write("#acc\tassembly_acc\n")
    else:
        with open(output_path, "r") as raw:
            for line in raw:
                if not line.startswith("#"):
                    data = line.rstrip().split("\t")
                    acc2log[data[0]] = data[1]

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
    search_term, esc_count = f"{accession}[{column}]", 0
    while esc_count < 3:
        try:
            handle = Entrez.esearch(db=database, term=search_term)
            genome_ids = Entrez.read(handle)["IdList"]
            break
        except (RuntimeError, urllib.error.HTTPError) as e:
            time.sleep(1)
            esc_count += 1
    else:
        logger.error(f"{accession} failed to search NCBI")
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
def collect_assembly_accs(
    ncbi_df,
    acc2log,
    api_key=0,
    column="assembly_acc",
    ncbi_column="Assembly Accession",
    database="assembly",
    output_path="",
    verbose=True,
    spacer="\t\t",
):

    count, failed = 0, []

    # for each row in the assembly, grab the accession number, form the search term for Entrez, use Entrez,
    if ncbi_column in {"assembly", "genome", "uid"}:
        out_df = ncbi_df[ncbi_df.index.isin(set(acc2log.keys()))]
        ncbi_df = ncbi_df[~ncbi_df.index.isin(set(acc2log.keys()))]
        if acc2log:
            out_df["assembly_acc"] = pd.Series(acc2log)
    for accession, row in tqdm(ncbi_df.iterrows(), total=len(ncbi_df)):
        if accession in acc2log:  # add all rows that have indices associated with this
            # query type
            out_df = pd.concat([out_df, row.to_frame().T])
            icount = 1
            test = str(accession) + "_" + str(icount)
            while test in acc2log:
                count += 1
                if "ome" in row.keys():
                    row["ome"] = None  # haven't assigned a mycotools ID yet
                out_df = pd.concat([out_df, row.to_frame().T])
                #                sys.exit()
                test = str(accession) + "_" + str(icount)
            continue
        elif accession.startswith(("GCA_", "GCF_")):
            row["assembly_acc"] = accession
            out_df = pd.concat([out_df, row.to_frame().T])
            acc2log[accession] = accession
            log_editor(
                output_path + "ncbiDwnld.log",
                accession,
                str(accession) + "\t" + accession,
            )

            continue

        elif pd.isnull(row[column]) or not row[column]:  # ignore blank entries
            acc2log[str(accession)] = ""
            failed.append([accession, datetime.strftime(row["version"], "%Y%m%d")])
            continue

        if ncbi_column not in {"uid"}:  # we already have the uid, no worries
            genome_id = esearch_ncbi(accession, ncbi_column, database="assembly")
        else:
            genome_id = [accession]

        if not genome_id:  # No IDs retrieved
            if "ome" in row.keys():
                accession = row["ome"]
            logger.error(spacer + "\t" + accession + " failed to find genome ID")
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

            record = esummary_ncbi(ID, database)
            record_info = record["DocumentSummarySet"]["DocumentSummary"][0]
            assemblyID = record_info["AssemblyAccession"]

            log_editor(
                output_path + "ncbiDwnld.log",
                str(new_acc),
                str(accession) + "\t" + assemblyID,
            )
            acc2log[str(new_acc)] = assemblyID
            row["assembly_acc"] = assemblyID
            out_df = pd.concat([out_df, row.to_frame().T])
            icount += 1

            count = wait_for_ncbi(count, api_key)

    return acc2log, failed, out_df


def run_datasets(
    include,
    accs_file,
    output_path,
    annotated,
    api=None,
    verbose=False,
    mute_stderr=False,
):
    """Run NCBI datasets to download genomes or metadata

    `mute_stderr` keeps datasets' own progress bar and error text off the
    terminal entirely; it is captured and logged at debug instead. It exists for
    callers that draw their own progress bar and would otherwise be overdrawn,
    so it belongs to the caller rather than the CLI and is not exposed as a
    flag. The caller stays responsible for reporting the failure itself."""
    dataset_scaf = [
        "datasets",
        "download",
        "genome",
        "accession",
        "--inputfile",
        accs_file,
    ]  # , '--no-progressbar']
    if include:
        dataset_scaf.extend(["--include", include])
    else:
        dataset_scaf.append("--dehydrated")
    api = clean_api_key(api)
    if api:
        dataset_scaf += ["--api-key", api]
    if annotated:
        dataset_scaf.append("--annotated")

    cwd = str(Path.cwd())
    os.chdir(output_path)
    if verbose:
        dataset_call = subprocess.call(dataset_scaf)
    else:
        # capture stderr rather than discard it: datasets reports why it failed
        # there, and a silenced call that fails leaves no other explanation
        proc = subprocess.run(
            dataset_scaf, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True
        )
        dataset_call = proc.returncode
        # datasets' progress bar shares this stream, so collapse each line to
        # its last carriage-return frame rather than reprinting every redraw
        # (note: splitlines() would break on \r, hence the explicit split)
        err = "\n".join(
            x.split("\r")[-1].rstrip()
            for x in proc.stderr.split("\n")
            if x.split("\r")[-1].strip()
        )
        if err:
            # the bar lands here too, so it is only worth surfacing on failure
            # -- and not even then once a caller has muted it
            loud = dataset_call and not mute_stderr
            (logger.error if loud else logger.debug)(err)
    os.chdir(cwd)

    return dataset_call


def compile_organism_names(unzip_path, spacer="\t"):
    acc2org, acc2meta = {}, {}
    failed = []
    with open(unzip_path + "data/assembly_data_report.jsonl", "r") as raw:
        for line in raw:
            data = json.loads(line.rstrip())
            if "accession" in data:
                acc = data["accession"]
                try:
                    org0 = data["organism"]["organismName"]
                except KeyError:
                    failed.append(acc)
                org1 = re.sub(r"[^ a-zA-Z0-9]", "", org0)
                org = org1.split()
                genus = org[0]
                if len(org) > 2:
                    species = org[1]
                    strain = "".join(org[2:])
                elif len(org) == 2:
                    species = org[1]
                    strain = ""
                else:
                    species = "sp."
                    strain = ""
                try:
                    if "infraspecificNames" in data["organism"]:
                        if "strain" in data["organism"]["infraspecificNames"]:
                            strain = data["organism"]["infraspecificNames"]["strain"]
                        elif "isolate" in data["organism"]["infraspecificNames"]:
                            strain = data["organism"]["infraspecificNames"]["isolate"]
                    elif not strain:
                        for attr in data["assemblyInfo"]["biosample"]["attributes"]:
                            if attr["name"].lower() == "strain":
                                strain = attr["value"]
                                if strain.lower() in {"missing", "none"}:
                                    strain = ""
                                break
                except KeyError:
                    failed.append(acc)
                    pass
                strain = re.sub(r"[^a-zA-Z0-9]", "", strain)

                try:
                    acc2meta[acc] = {
                        "accession": data["assemblyInfo"]["assemblyName"],
                        "submitter": data["annotationInfo"]["provider"],
                    }
                except KeyError:
                    failed.append(acc)
                    pass
                #                if not strain:
                #                   eprint(f'{spacer}WARNING: {acc} no strain metadata', flush = True)
                acc2org[acc] = {"genus": genus, "species": species, "strain": strain}
    return acc2org, acc2meta, failed


def resolve_paired_accs(
    accs, acc_file, output_path, api=None, spacer="\t\t", summary_chunk=500
):
    """Map accessions onto their counterpart in the other NCBI repository.

    GenBank and RefSeq version their assemblies independently, so the
    counterpart of GCA_017499595.2 is GCF_017499595.1 -- swapping the prefix
    and keeping the version names GCF_017499595.2, which does not exist. NCBI
    reports the pair it actually holds, and reports none for an assembly that
    was never mirrored, so an accession missing from the result has nothing to
    reattempt rather than a counterpart that failed.

    These are metadata records rather than genomes, so the chunk is sized well
    above the genome chunk for the same reason the data reports are. A chunk
    that cannot be resolved is skipped rather than fatal: the accessions in it
    simply keep their existing failure, which is the outcome without a
    counterpart anyway."""

    acc2pair = {}
    accs = [str(x) for x in accs]
    for i in range(0, len(accs), summary_chunk):
        batch = accs[i : i + summary_chunk]
        with open(acc_file, "w") as out:
            out.write("\n".join(batch))
        cmd = [
            "datasets",
            "summary",
            "genome",
            "accession",
            "--inputfile",
            acc_file,
            "--as-json-lines",
        ]
        api_key = clean_api_key(api)
        if api_key:
            cmd += ["--api-key", api_key]

        cwd = str(Path.cwd())
        os.chdir(output_path)
        try:
            proc = subprocess.run(
                cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
            )
        # an absent datasets is the caller's problem to report, not a reason to
        # end a run that has already downloaded everything it could
        except FileNotFoundError:
            logger.debug(f"{spacer}\tdatasets is unavailable to resolve pairs")
            os.chdir(cwd)
            return acc2pair
        os.chdir(cwd)
        if proc.returncode:
            logger.debug(
                f"{spacer}\tcould not resolve paired accessions: "
                + f"{proc.stderr.rstrip()}"
            )
            continue

        for line in proc.stdout.split("\n"):
            if not line.strip():
                continue
            try:
                report = json.loads(line)
            except json.JSONDecodeError:
                continue
            if report.get("accession") and report.get("paired_accession"):
                acc2pair[report["accession"]] = report["paired_accession"]

    return acc2pair


def parse_datasets(datasets_path, unzip_base, req_files, spacer="\t"):
    """Unzip, identify complete downloads, parse file outputs and metadata,
    report missing data to check alternative repository"""
    try:
        with zipfile.ZipFile(datasets_path, "r") as zip_ref:
            zip_ref.extractall(unzip_base)
    # a dropped transfer leaves a truncated archive, and a call that died before
    # it opened the file leaves none at all; both are the same failure here
    except (zipfile.BadZipFile, FileNotFoundError):
        return False, False, False
    Path(datasets_path).unlink()
    unzip_path = unzip_base + "ncbi_dataset/"

    type2ncbi = {
        "fna": "GENOMIC_NUCLEOTIDE_FASTA",
        "gff3": "GFF3",
        "faa": "PROTEIN_FASTA",
        "rna": "RNA_NUCLEOTIDE_FASTA",
    }
    ncbi2type = {v: k for k, v in type2ncbi.items()}

    failed = []
    data_dict = read_json(unzip_path + "data/dataset_catalog.json")
    acc2data = {}
    basename = unzip_path + "data/"
    for data in data_dict["assemblies"]:
        if "accession" in data:
            files_p = {x["fileType"]: x["filePath"] for x in data["files"]}

            files = {ncbi2type[k]: basename + v for k, v in files_p.items()}
            # if there are missing files
            if req_files.difference(set(files.keys())):
                failed.append(data["accession"])
                for t, f_ in files.items():
                    if Path(f_).is_file():
                        Path(f_).unlink()
            else:
                acc2data[data["accession"]] = files

    acc2org, acc2meta, org_failed = compile_organism_names(unzip_path, spacer)

    return acc2data, acc2org, failed


def download_datasets(
    accs,
    acc_file,
    include,
    req_files,
    annotated,
    output_path,
    api=None,
    verbose=False,
    spacer="\t\t",
    max_attempts=3,
    min_accs=10,
    max_dead=None,
):
    """Write a chunk of accessions to acc_file, download them via NCBI datasets,
    and parse the output. Returns acc2files, acc2org, failed (each False only if
    nothing in `accs` could be retrieved).

    NCBI drops these transfers mid-stream, and the larger the archive the likelier
    it is to be dropped -- so a batch that fails every attempt is halved and its
    halves are downloaded separately rather than retried whole. Retrying whole
    restarts a multi-GB transfer from zero at the size that was already failing
    and loses the entire batch when it fails again; halving both shrinks the
    transfer and keeps whatever the other half retrieved.

    Splitting is abandoned once `max_dead` batches have come back with nothing
    and nothing at all has been retrieved, because a repository serving nothing
    at any size is down rather than overloaded and halving it only multiplies the
    calls made against it. The default tolerance is the depth of the split tree:
    halving descends the first half before trying its sibling, so the first batch
    that can succeed may be that many failures away."""
    zip_path = output_path + "ncbi_dataset.zip"
    if max_dead is None:
        max_dead = 2 + math.ceil(math.log2(max(len(accs) / max(min_accs, 1), 1)))

    def attempt_batch(batch):
        """Download one batch, retrying a dropped transfer at the same size.
        Returns the parsed result, or None if every attempt failed."""
        for attempt in range(1, max_attempts + 1):
            if attempt > 1:
                # datasets leaves a truncated archive behind when the stream is
                # reset, and a call that dies before reopening it would hand
                # that same partial file back to the parser
                Path(zip_path).unlink(missing_ok=True)
                # a reset usually means NCBI is loaded; retrying instantly adds
                # to the load that caused it
                time.sleep(3 * 2 ** (attempt - 2))
                logger.debug(f"{spacer}\tAttempt {attempt} ({len(batch)} accessions)")

            with open(acc_file, "w") as out:
                out.write("\n".join([str(x) for x in batch]))

            code = run_datasets(
                include,
                acc_file,
                output_path,
                api=api,
                verbose=verbose,
                annotated=annotated,
                # a failed attempt is routine now that it is retried and split,
                # so datasets' error text is debug detail; the failures that
                # survive the splitting are what the caller is told about
                mute_stderr=True,
            )
            if code:
                continue

            parsed = parse_datasets(zip_path, output_path, req_files, spacer)
            if parsed[0] is not False:
                return parsed
        Path(zip_path).unlink(missing_ok=True)
        return None

    logger.debug(f"{spacer}Downloading data")
    acc2files, acc2org, failed = {}, {}, []
    queue, retrieved, dead_streak = [list(accs)], False, 0
    while queue:
        batch = queue.pop(0)
        parsed = attempt_batch(batch)
        if parsed is not None:
            retrieved = True
            dead_streak = 0
            acc2files.update(parsed[0])
            acc2org.update(parsed[1])
            failed.extend(parsed[2])
            continue

        # once anything has come back the repository is serving, and every later
        # failure is that batch's problem rather than grounds to stop splitting
        if not retrieved:
            dead_streak += 1
            if dead_streak >= max_dead:
                logger.debug(
                    f"{spacer}\t{dead_streak} batches returned nothing; abandoning "
                    + f"{sum(len(b) for b in queue) + len(batch)} accessions"
                )
                break
        # below this size the transfer is no longer what is failing, so the
        # accessions are abandoned to the caller rather than split again. It
        # derives what it never received from acc2files, so they need no
        # accounting here
        if len(batch) <= min_accs:
            continue

        mid = len(batch) // 2
        logger.debug(
            f"{spacer}\t{len(batch)} accessions failed {max_attempts} attempts; "
            + "splitting"
        )
        queue[:0] = [batch[:mid], batch[mid:]]

    if not retrieved:
        return False, False, False

    return acc2files, acc2org, failed


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
    chunk=25,
):

    # initialize run directory and information
    output_path = format_path(output_path)
    #    file_types = prepare_folders(output_path, gff3, proteome,
    #                                assembly, transcript)

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

    # Uppercase before the index is taken from it, not after. The index is what
    # downloads and failures are matched on downstream, and NCBI reports its
    # accessions uppercase, so normalizing the column afterwards leaves a
    # lowercase index that matches neither and drops the row out of the run
    if "assembly_acc" in ncbi_df.keys():
        ncbi_df["assembly_acc"] = [str(x).upper() for x in ncbi_df["assembly_acc"]]

    # preserve the original column, but index ncbi_df on it as well
    ncbi_df = ncbi_df.set_index(pd.Index(list(ncbi_df[column])))

    ## CHANGE TO ACCOMODATE BIOSAMPLE/OTHER NCBICOLUMNS
    if ncbi_column.lower() != "assembly":
        logger.debug(f"{spacer}Assembling NCBI ftp directories")
        acc2log = compile_log(output_path + "ncbiDwnld.log")
        acc2log, failed, ncbi_df = collect_assembly_accs(
            ncbi_df,
            acc2log,
            ncbi_column=ncbi_column,
            column=column,
            api_key=api,
            output_path=output_path,
            verbose=verbose,
            spacer=spacer,
        )
        column = "assembly_acc"

        # collect_assembly_accs supplies the column here, so it needs the same
        # normalization before this index is taken from it
        ncbi_df["assembly_acc"] = [str(x).upper() for x in ncbi_df["assembly_acc"]]
        ncbi_df = ncbi_df.set_index(pd.Index(list(ncbi_df[column])))
    new_df = pd.DataFrame()

    ## GUARANTEE ASSEMBLY ACCESSIONS ARE LABELED THIS COLUMN NAME
    acc_file = output_path + "assembly_accs.txt"

    include = ""
    req_files = set()
    if assembly:
        include += "genome,"
        req_files.add("fna")
    if proteome:
        include += "protein,"
        req_files.add("faa")
    if gff3:
        include += "gff3,"
        req_files.add("gff3")
    if transcript:
        include += "rna,"
        req_files.add("rna")
    include = include[:-1]

    # Download via datasets
    if proteome or gff3:
        annotated = True
    else:
        annotated = False

    # Chunk the accessions so datasets is called on `chunk` accessions at a time
    all_accs = [str(x) for x in list(ncbi_df["assembly_acc"])]
    acc_chunks = [all_accs[i : i + chunk] for i in range(0, len(all_accs), chunk)]

    # Run downloads chunk-by-chunk, accumulating results
    acc2files, acc2org, failed = {}, {}, []
    dead_chunks, consecutive_dead = [], 0
    for chunk_i, acc_chunk in enumerate(acc_chunks):
        if len(acc_chunks) > 1:
            logger.debug(
                f"{spacer}Chunk {chunk_i + 1}/{len(acc_chunks)} "
                + f"({len(acc_chunk)} accessions)"
            )
        c_acc2files, c_acc2org, c_failed = download_datasets(
            acc_chunk,
            acc_file,
            include,
            req_files,
            annotated,
            output_path,
            api=api,
            verbose=verbose,
            spacer=spacer,
        )
        if c_acc2files is False:
            # one chunk NCBI will not serve is not worth discarding the chunks
            # that succeeded; its accessions fall out of the set difference below
            # and are reported as failures. Several in a row is the repository or
            # the connection being gone, which the remaining chunks cannot fix
            dead_chunks.append(chunk_i + 1)
            consecutive_dead += 1
            if consecutive_dead >= 3:
                logger.error(
                    f"{spacer}ncbiDwnld failed {consecutive_dead} consecutive chunks"
                )
                sys.exit(10)
            logger.warning(
                f"{spacer}chunk {chunk_i + 1}/{len(acc_chunks)} failed; continuing"
            )
            continue
        consecutive_dead = 0
        acc2files.update(c_acc2files)
        acc2org.update(c_acc2org)
        failed.extend(c_failed)

    if dead_chunks:
        logger.warning(
            f"{spacer}{len(dead_chunks)}/{len(acc_chunks)} chunk(s) failed to "
            + "download; their accessions are reported as failures"
        )
        logger.debug(f"{spacer}failed chunks: {dead_chunks}")

    failed.extend(
        sorted(set(ncbi_df["assembly_acc"]).difference(set(acc2files.keys())))
    )

    # Attempt the counterpart repository for whatever this one did not serve
    acc2pair = {}
    if failed:
        logger.debug(f"{spacer}Attempting alternative repository for failed downloads")
        acc_file_re = output_path + "assembly_accs.reattempt.txt"
        acc2pair = resolve_paired_accs(
            failed, acc_file_re, output_path, api=api, spacer=spacer
        )
        pair2acc = {v: k for k, v in acc2pair.items()}
        reattempt_acc = sorted(pair2acc)
        unpaired = len(failed) - len(reattempt_acc)
        if unpaired:
            logger.debug(
                f"{spacer}\t{unpaired} failed accession(s) are not mirrored in "
                + "the alternative repository"
            )

        # Chunk the reattempt accessions as well
        reattempt_chunks = [
            reattempt_acc[i : i + chunk]
            for i in range(0, len(reattempt_acc), chunk)
        ]
        recovered = set()
        for acc_chunk in reattempt_chunks:
            c_acc2files, c_acc2org, c_failed = download_datasets(
                acc_chunk,
                acc_file_re,
                include,
                req_files,
                annotated,
                output_path,
                api=api,
                verbose=verbose,
                spacer=spacer,
            )
            if c_acc2files is False:
                continue
            acc2files = {**acc2files, **c_acc2files}
            acc2org = {**acc2org, **c_acc2org}
            recovered.update(pair2acc[x] for x in c_acc2files if x in pair2acc)

        # only what was actually retrieved leaves the failure list. A counterpart
        # NCBI never served is absent from the archive rather than named in
        # c_failed, so subtracting what came back is the only accounting that
        # sees it -- rebuilding the list from c_failed instead dropped every
        # unmirrored accession out of the run without a trace
        failed = sorted(set(failed).difference(recovered))

    # Parse download output, add to df
    # Report failed
    failed_set = set(failed)
    rep_failed = []

    def report_failure(acc, row):
        try:
            rep_failed.append([acc, datetime.strftime(row["version"], "%Y%m%d")])
        except TypeError:
            rep_failed.append([acc, row["version"]])

    for acc, row in ncbi_df.iterrows():
        # the counterpart NCBI actually holds, which versions independently of
        # acc; None when the assembly is not mirrored at all
        check_acc = acc2pair.get(acc)
        if acc in failed_set:
            report_failure(acc, row)
        elif acc in acc2files:
            try:
                for file_t, file_p in acc2files[acc].items():
                    ncbi_df.at[acc, file_t] = file_p
                for tax, name in acc2org[acc].items():
                    ncbi_df.at[acc, tax] = name
                new_df = pd.concat([new_df, ncbi_df.loc[acc].to_frame().T])
            except AttributeError:  # multiple entries
                logger.warning(f"{spacer}{acc} is redundant")
                for acc1, row1 in ncbi_df.loc[acc].iterrows():
                    for file_t, file_p in acc2files[acc].items():
                        row1[file_t] = file_p
                    for tax, name in acc2org[acc].items():
                        row1[tax] = name
                    new_df = pd.concat([new_df, row1.to_frame().T])

        elif check_acc in acc2files:
            try:
                for file_t, file_p in acc2files[check_acc].items():
                    ncbi_df.at[acc, file_t] = file_p
                for tax, name in acc2org[check_acc].items():
                    ncbi_df.at[acc, tax] = name
                new_df = pd.concat([new_df, ncbi_df.loc[acc].to_frame().T])
            except AttributeError:  # multiple entries
                logger.debug(f"{spacer}\t{check_acc} is redundant")
                for acc1, row1 in ncbi_df.loc[acc].iterrows():
                    for file_t, file_p in acc2files[check_acc].items():
                        row1[file_t] = file_p
                    for tax, name in acc2org[check_acc].items():
                        row1[tax] = name
                    new_df = pd.concat([new_df, row1.to_frame().T])

        else:
            # nothing retrieved it and nothing named it a failure. Without this
            # the accession leaves the run in neither new_df nor rep_failed,
            # which is how an initialization silently lost 214 of 525 genomes
            logger.debug(f"{spacer}\t{acc} was neither retrieved nor reported")
            report_failure(acc, row)

    if "fna" in new_df.keys():
        new_df = new_df.rename(columns={"fna": "assemblyPath"})
    if "gff3" in new_df.keys():
        new_df = new_df.rename(columns={"gff3": "gffPath"})
    new_df = new_df.reset_index()
    return new_df, rep_failed


def get_sra(assembly_acc, fastqdump="fastq-dump", pe=True):

    handle = Entrez.esearch(db="SRA", term=assembly_acc)
    ids = Entrez.read(handle)["IdList"]
    for id in ids:
        handle = Entrez.esummary(db="SRA", id=id, report="full")
        records = Entrez.read(handle, validate=False)
        for record in records:
            srr = re.search(r'Run acc="(S\w+\d+)"', record["Runs"])[1]
            logger.info("\t\t" + srr)
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
                    if Path(f"{srr}_1.fastq.gz").is_file():
                        #                  if os.path.isfile(srr + '_1.fastq'):
                        #                        cmd = subprocess.call(['gzip', f'{srr}_1.fastq'])
                        #                       cmd = subprocess.call(['gzip', f'{srr}_2.fastq'])
                        shutil.move(
                            f"{srr}_1.fastq.gz", f"{assembly_acc}_{srr}_1.fq.gz"
                        )
                        shutil.move(
                            f"{srr}_2.fastq.gz", f"{assembly_acc}_{srr}_2.fq.gz"
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
                    if Path(f"{srr}.fastq.gz").is_file():
                        shutil.move(f"{srr}.fastq.gz", f"{assembly_acc}_{srr}.fq.gz")
                    else:
                        logger.error("file failed")


def go_sra(df, output=str(Path.cwd()) + "/", pe=True, column="sra"):

    print()
    sra_dir = output + "sra/"
    if not Path(sra_dir).is_dir():
        Path(sra_dir).mkdir()
    os.chdir(sra_dir)
    fastqdump = find_execs("fastq-dump", exit={"fastq-dump"})
    count = 0

    for i, row in df.iterrows():
        logger.info("\t" + row[column])
        get_sra(row[column], fastqdump[0])
        count += 1
        if count >= 10:
            time.sleep(1)
            count = 0


def cli():
    parser = argparse.ArgumentParser(
        description="GenBank/RefSeq downloading utility. Downloads "
        + "accession by accession"
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
    parser.add_argument("--api", help="NCBI API key for high query rate")
    parser.add_argument(
        "--chunk",
        type=int,
        default=25,
        help="Accessions to download per datasets call; larger chunks are "
        + "likelier to be reset mid-transfer by NCBI; DEFAULT: 25",
    )
    args = parser.parse_args()
    setup_logging(verbose=getattr(args, "verbose", False))

    if args.api:
        ncbi_api = args.api
    else:
        ncbi_api, jgi_email, jgi_pwd = login_check(jgi=False)
    if ncbi_api:
        Entrez.api_key = ncbi_api

    if not args.output:
        output = mk_output(None, "download_ncbi")
    else:
        output = format_path(args.output)

    find_execs("datasets", exit={"datasets"})

    args_dict = {
        "NCBI Table": args.input,
        "Assemblies": args.assembly,
        "Proteomes": args.proteome,
        ".gff3's": args.gff3,
        "Transcripts": args.transcript,
        "SRA": args.sra,
        "Chunk": args.chunk,
    }

    start_time = intro("Download NCBI files", args_dict)
    if args.sra:
        if Path(format_path(args.input)).is_file():
            if not args.column:
                go_sra(
                    pd.read_csv(format_path(args.input), sep="\t", names=["sra"]),
                    output,
                    pe=args.paired,
                    column="sra",
                )
            else:
                go_sra(
                    pd.read_csv(format_path(args.input), sep="\t"),
                    output,
                    pe=args.paired,
                    column=args.column,
                )
        else:
            go_sra(
                pd.DataFrame({"sra": split_input(args.input)}),
                output,
                pe=args.paired,
                column="sra",
            )
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
            ncbi_df = pd.DataFrame({"assembly_acc": split_input(args.input)})
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
            chunk=args.chunk,
        )
        new_df = new_df.rename(columns={"index": "#assembly_accession"})
        new_df["source"] = "ncbi"
        new_df["useRestriction (yes/no)"] = "no"
        if 0 in new_df.columns:
            del new_df[0]

        new_df.to_csv(output + "ncbiDwnld.predb", sep="\t", index=None)
        if failed:
            logger.error(",".join([str(x[0]) for x in failed]))

    outro(start_time)


if __name__ == "__main__":
    cli()
