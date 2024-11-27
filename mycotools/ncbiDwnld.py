#! /usr/bin/env python3

# NEED a db check to ensure the log is relevant to the input
# NEED to convert to datasets
# NEED to consider refseq genomes with annotations when genbank doesn't have them

import os
import re
import sys
import gzip
import json
import time
import shutil
import urllib
import zipfile
import argparse
import subprocess
import numpy as np
import pandas as pd
from contextlib import closing
from tqdm import tqdm
from Bio import Entrez
from datetime import datetime
from mycotools.lib.kontools import intro, outro, format_path, prep_output, \
                                   mkOutput, eprint, vprint, findExecs, \
                                   read_json, split_input
from mycotools.lib.dbtools import log_editor, loginCheck, mtdb, read_tax

pd.options.mode.chained_assignment = None

def ncbidb2df(data, stdin = False):
    import pandas as pd, pandas
    columns = mtdb.columns
    if isinstance(data, mtdb):
        db_df = pd.DataFrame(data.reset_index())
    elif not stdin:
        data = format_path( data )
        db_df = pd.read_csv( data, sep='\t' )
        if 'ome' not in set( db_df.columns ) and 'assembly_acc' not in set( db_df.columns ):
            db_df = pd.read_csv( data, sep = '\t', header = None )
    else:
        db_df = pd.read_csv( StringIO(data), sep='\t' )
        if 'ome' not in set( db_df.columns ) and 'assembly_acc' not in set( db_df.columns ):
            db_df = pd.read_csv( StringIO(data), sep = '\t', header = None )

    db_df = db_df.fillna('')

    return db_df


def prepare_folders(output_path, gff, prot, assem, transcript):

    file_types = []
    if assem:
        if not os.path.exists(output_path + 'fna'):
            os.mkdir(output_path + 'fna')
        file_types.append('fna')
    if gff:
        if not os.path.exists(output_path + 'gff3'):
            os.mkdir(output_path + 'gff3')
        file_types.append('gff3')
    if prot:
        if not os.path.exists(output_path + 'faa'):
            os.mkdir(output_path + 'faa')
        file_types.append('faa') 
    if transcript:
        if not os.path.exists(output_path + 'transcript'):
            os.mkdir(output_path + 'transcript')
        file_types.append('transcript')

    return file_types  


def compile_log(output_path):

    acc2log = {}
    if not os.path.isfile(output_path):
        with open(output_path, 'w' ) as out:
            out.write('#acc\tassembly_acc\n')
    else:
        with open(output_path, 'r') as raw:
            for line in raw:
                if not line.startswith('#'):
                    data = line.rstrip().split('\t')
                    acc2log[data[0]] = data[1]

    return acc2log

def wait_for_ncbi(count, api = False):
    if count >= 2:
        if not api:
            time.sleep(1)
            count = 0
        elif count >= 7:
            time.sleep(1)
            count = 0
    return count

def esearch_ncbi(accession, column, database = 'assembly'):
    search_term, esc_count = f'{accession}[{column}]', 0
    while esc_count < 3:
        try:
            handle = Entrez.esearch(db=database, term=search_term)
            genome_ids = Entrez.read(handle)['IdList']
            break
        except (RuntimeError, urllib.error.HTTPError) as e:
            time.sleep(1)
            esc_count += 1
    else:
        print('\tERROR:', accession, 'failed to search NCBI')
        return None
    return genome_ids

def esummary_ncbi(ID, database):

    esc_count = 0
    while esc_count < 10:
        esc_count += 1
        try:
            handle = Entrez.esummary(db=database, id=ID, report="full")
            record = Entrez.read(handle, validate = False)
        except urllib.error.HTTPError:
            time.sleep(0.1)
            continue
        if database == 'assembly':
            try: # is it populated with an FTP?
                ftp_path = str(record['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_GenBank'])
            except IndexError: # wait a sec and retry
                time.sleep(1)
                continue
        break
    else: # too many failed attempts
        if esc_count >= 10:
            raise urllib.error.HTTPError('\tERROR: FTP request failed')

    return record


# collects paths to download proteomes and assemblies
def collect_assembly_accs(
    ncbi_df, acc2log, api_key=0, column = 'assembly_acc',
    ncbi_column='Assembly Accession', database="assembly", output_path = '',
    verbose=True, spacer = '\t\t'
    ):

    count, failed = 0, []

# for each row in the assembly, grab the accession number, form the search term for Entrez, use Entrez,
    if ncbi_column in {'assembly', 'genome', 'uid'}:
        out_df = ncbi_df[ncbi_df.index.isin(set(acc2log.keys()))]
        ncbi_df = ncbi_df[~ncbi_df.index.isin(set(acc2log.keys()))]
        if acc2log:
            out_df['assembly_acc'] = pd.Series(acc2log)
    for accession, row in tqdm(ncbi_df.iterrows(), total = len(ncbi_df)):
        if accession in acc2log: # add all rows that have indices associated with this
            # query type
            out_df = pd.concat([out_df, row.to_frame().T])
            icount = 1
            test = str(accession) + '_' + str(icount)
            while test in acc2log:
                count += 1
                if 'ome' in row.keys():
                    row['ome'] = None # haven't assigned a mycotools ID yet
                out_df = pd.concat([out_df, row.to_frame().T])
#                sys.exit()
                test = str(accession) + '_' + str(icount)
            continue
        elif accession.startswith(('GCA_', 'GCF_')):
            row['assembly_acc'] = accession
            out_df = pd.concat([out_df, row.to_frame().T])
            acc2log[accession] = accession
            log_editor(output_path + 'ncbiDwnld.log', accession, 
                       str(accession) + '\t' + accession)

            continue

        elif pd.isnull(row[column]) or not row[column]: # ignore blank entries
            acc2log[str(accession)] = ''
            failed.append([accession, datetime.strftime(row['version'], '%Y%m%d')])
            continue

        if ncbi_column not in {'uid'}: # we already have the uid, no worries
            genome_id = esearch_ncbi(accession, ncbi_column, database = 'assembly')
        else:
            genome_id = [accession]

        if not genome_id: # No IDs retrieved
            if 'ome' in row.keys():
                accession = row['ome']
            eprint(spacer + '\t' + accession + ' failed to find genome ID', flush = True)
            try:
                failed.append([accession, datetime.strftime(row['version'], '%Y%m%d')])
            except TypeError: # if the row can't be formatted as a date
                failed.append([accession, row['version']])
            continue

        if ncbi_column in {'Assembly Accession', 'assembly', 
                           'genome', 'uid'}: # be confident it is the most
        # recent assembly UID
            genome_id = [str(max([int(i) for i in genome_id]))]

        icount = 0
        for ID in genome_id:
            if icount:
                new_acc = str(accession) + '$' + str(icount)
            else:
                new_acc = accession

            record = esummary_ncbi(ID, database)
            record_info = record['DocumentSummarySet']['DocumentSummary'][0]
            assemblyID = record_info['AssemblyAccession']
           
            esc_count = 0
            log_editor(output_path + 'ncbiDwnld.log', str(new_acc), 
                       str(accession) + '\t' + assemblyID)
            acc2log[str(new_acc)] = assemblyID
            row['assembly_acc'] = assemblyID
            out_df = pd.concat([out_df, row.to_frame().T])
            icount += 1
    
            count = wait_for_ncbi(count, api_key)
    
    return acc2log, failed, out_df

def run_datasets(include, accs_file, output_path, 
                 annotated, api = None, verbose = False):
    """Run NCBI datasets to download genomes or metadata"""
    dataset_scaf = ['datasets', 'download', 'genome', 'accession',
                    '--inputfile', accs_file] #, '--no-progressbar']
    if include:
        dataset_scaf.extend(['--include', include])
    else:
        dataset_scaf.append('--dehydrated')
    if api:
        dataset_scaf += ['--api-key', api]
    if annotated:
        dataset_scaf.append('--annotated')

    cwd = os.getcwd()
    os.chdir(output_path)
    if verbose:
        v = None
    else:
        v = subprocess.DEVNULL
    dataset_call = subprocess.call(dataset_scaf, stdout = v, stderr = v)
    os.chdir(cwd)

    return dataset_call

def compile_organism_names(unzip_path, spacer = '\t'):
    acc2org, acc2meta = {}, {}
    failed = []
    with open(unzip_path + 'data/assembly_data_report.jsonl', 'r') as raw:
        for line in raw:
            data = json.loads(line.rstrip())
            if 'accession' in data:
                acc = data['accession']
                try:
                    org0 = data['organism']['organismName']
                except KeyError:
                    failed.append(acc)
                org1 = re.sub(r'[^ a-zA-Z0-9]', '', org0) 
                org = org1.split()
                genus = org[0] 
                if len(org) > 2:
                    species = org[1]
                    strain = ''.join(org[2:])
                elif len(org) == 2:
                    species = org[1]
                    strain = ''
                else:
                    species = 'sp.'
                    strain = ''
                try:
                    if 'infraspecificNames' in data['organism']:
                        if 'strain' in data['organism']['infraspecificNames']:
                            strain = data['organism']['infraspecificNames']['strain']
                        elif 'isolate' in data['organism']['infraspecificNames']:
                            strain = data['organism']['infraspecificNames']['isolate']
                    elif not strain:
                        for attr in data['assemblyInfo']['biosample']['attributes']:
                            if attr['name'].lower() == 'strain':
                                strain = attr['value']
                                if strain.lower() in {'missing', 'none'}:
                                    strain = ''
                                break
                except KeyError:
                    failed.append(acc)
                    pass
                strain = re.sub(r'[^a-zA-Z0-9]', '', strain) 

                try: 
                    acc2meta[acc] = {'accession': data['assemblyInfo']['assemblyName'],
                                     'submitter': data['annotationInfo']['provider']}
                except KeyError:
                    failed.append(acc)
                    pass
#                if not strain:
 #                   eprint(f'{spacer}WARNING: {acc} no strain metadata', flush = True)
                acc2org[acc] = {'genus': genus, 'species': species, 'strain': strain}
    return acc2org, acc2meta, failed

def parse_datasets(datasets_path, unzip_base, req_files, spacer = '\t'):
    """Unzip, identify complete downloads, parse file outputs and metadata, 
    report missing data to check alternative repository"""
    try:
        with zipfile.ZipFile(datasets_path, 'r') as zip_ref:
            zip_ref.extractall(unzip_base) 
    except zipfile.BadZipFile:
        return False, False, False
    os.remove(datasets_path)
    unzip_path = unzip_base + 'ncbi_dataset/'

    type2ncbi = {'fna': 'GENOMIC_NUCLEOTIDE_FASTA',
                 'gff3': 'GFF3',
                 'faa': 'PROTEIN_FASTA',
                 'rna': 'RNA_NUCLEOTIDE_FASTA'}
    ncbi2type = {v: k for k, v in type2ncbi.items()}

    failed = []
    data_dict = read_json(unzip_path + 'data/dataset_catalog.json') 
    acc2data = {}
    basename = unzip_path + 'data/'
    for data in data_dict['assemblies']:
        if 'accession' in data:
            files_p = {x['fileType']: x['filePath'] \
                     for x in data['files']}
             
            files = {ncbi2type[k]: basename + v for k, v in files_p.items()}
            # if there are missing files
            if req_files.difference(set(files.keys())):
                failed.append(data['accession'])
                for t, f_ in files.items():
                    if os.path.isfile(f_):
                        os.remove(f_)
            else:
                acc2data[data['accession']] = files

    acc2org, acc2meta, org_failed = compile_organism_names(unzip_path, spacer)

    return acc2data, acc2org, failed

def main( 
    api = None, 
    assembly = True, proteome = False, gff3 = True, transcript = False,
    ncbi_df = False, remove = False, output_path = os.getcwd(), verbose = False,
    column = 'assembly_acc', ncbi_column = 'Assembly', check_MD5 = True,
    spacer = '\t\t'
    ):

    # initialize run directory and information
    output_path = format_path(output_path)
#    file_types = prepare_folders(output_path, gff3, proteome, 
 #                                assembly, transcript)

    # check if ncbi_df is a dataframe, and import if not
    if not isinstance(ncbi_df, pd.DataFrame) and os.path.isfile(ncbi_df):
        ncbi_df = ncbidb2df(ncbi_df)
    if len(ncbi_df.index) == 0:
        ncbi_df = pd.DataFrame(
            {i: [v] for i, v in enumerate(list(ncbi_df.keys()))}
            )

    # make the modify date from a standard NCBI table the version if it does
    # not otherwise exist, else there isn't a version to reference
    if 'Modify Date' in ncbi_df.keys() and not 'version' in ncbi_df.keys():
        ncbi_df['version'] = pd.to_datetime(ncbi_df['Modify Date'])
    elif 'version' not in ncbi_df.keys():
        ncbi_df['version'] = ''

    # preserve the original column, but index ncbi_df on it as well
    ncbi_df = ncbi_df.set_index(pd.Index(list(ncbi_df[column])))

    ## CHANGE TO ACCOMODATE BIOSAMPLE/OTHER NCBICOLUMNS
    if ncbi_column.lower() != 'assembly':
        vprint(f'{spacer}Assembling NCBI ftp directories', v = verbose, flush = True)
        acc2log = compile_log(output_path + 'ncbiDwnld.log')
        acc2log, failed, ncbi_df = collect_assembly_accs( 
            ncbi_df, acc2log,
            ncbi_column = ncbi_column, column = column, api_key=api,
            output_path = output_path, verbose = verbose,
            spacer = spacer
            )
        column = 'assembly_acc'

        ncbi_df = ncbi_df.set_index(pd.Index(list(ncbi_df[column])))
    new_df = pd.DataFrame()

    ## GUARANTEE ASSEMBLY ACCESSIONS ARE LABELED THIS COLUMN NAME
    acc_file = output_path + 'assembly_accs.txt'
    ncbi_df['assembly_acc'] = list([x.upper() for x in ncbi_df['assembly_acc']])
    with open(acc_file, 'w') as out:
        out.write('\n'.join([str(x) for x in list(ncbi_df['assembly_acc'])]))

    include = ''
    req_files = set()
    if assembly:
        include += 'genome,'
        req_files.add('fna')
    if proteome:
        include += 'protein,'
        req_files.add('faa')
    if gff3:
        include += 'gff3,'
        req_files.add('gff3')
    if transcript:
        include += 'rna,'
        req_files.add('rna')
    include = include[:-1]

    # Download via datasets
    if proteome or gff3:
        annotated = True
    else:
        annotated = False

    # Run downloads
    count = 0
    while count < 3:
        if not count:
            vprint(f'{spacer}Downloading data', v = verbose, flush = True)
            count += 1
        else:
            count += 1
            vprint(f'{spacer}\tAttempt {count}', v = verbose, flush = True)

        run_datasets(include, acc_file, output_path, api = api,
                     verbose = verbose, annotated = annotated)

            # Parse download output, add to df
        acc2files, acc2org, failed = parse_datasets(output_path + 'ncbi_dataset.zip', 
                       output_path, req_files, spacer)
        if acc2files == False and acc2org == False and failed == False:
            continue
        else:
            break

    if acc2files == False and acc2org == False and failed == False:
        eprint(f'{spacer}ERROR: ncbiDwnld failed {count} attempts', 
            flush = True)
        # maybe add a fallback to the old methodology here
        eprint(f'{spacer}Consider --fallback',
            flush = True)
        sys.exit(10)

    failed.extend(sorted(set(ncbi_df['assembly_acc']).difference(set(acc2files.keys()))))

    # Attempt RefSeq accessions 
    if failed:
        vprint(f'{spacer}Attempting alternative repository for failed downloads', 
               v = verbose, flush = True)
        reattempt_acc = []
        for acc in failed:
            if acc.upper().startswith('GCA'):
                reattempt_acc.append(acc.upper().replace('GCA_', 'GCF_'))
            elif acc.upper().startswith('GCF'):
                reattempt_acc.append(acc.upper().replace('GCF_', 'GCA_'))    
        acc_file_re = output_path + 'assembly_accs.reattempt.txt'
        with open(acc_file_re, 'w') as out:
            out.write('\n'.join(reattempt_acc))

        run_datasets(include, acc_file_re, output_path, verbose = verbose,
                     annotated = annotated)
        acc2files_r, acc2org_r, failed_r = parse_datasets(output_path + 'ncbi_dataset.zip',
                                              output_path, req_files)
        acc2files = {**acc2files, **acc2files_r}
        acc2org = {**acc2org, **acc2org_r}

        failed = []
        for acc in failed_r:
            if acc.upper().startswith('GCA'):
                failed.append(acc.upper().replace('GCA_', 'GCF_'))
            elif acc.upper().startswith('GCF'):
                failed.append(acc.upper().replace('GCF_', 'GCA_'))

    # Parse download output, add to df
    # Report failed
    failed_set = set(failed)
    rep_failed = []
    for acc, row in ncbi_df.iterrows():
        if acc.startswith('GCA_'):
            check_acc = acc.replace('GCA_', 'GCF_')
        elif acc.startswith('GCF_'):
            check_acc = acc.replace('GCF_', 'GCA_')
        if acc in failed_set:
            try:
                rep_failed.append([acc, datetime.strftime(row['version'], '%Y%m%d')])
            except TypeError:
                rep_failed.append([acc, row['version']])
        elif acc in acc2files:
            try:
                for file_t, file_p in acc2files[acc].items():
                    ncbi_df.at[acc, file_t] = file_p
                for tax, name in acc2org[acc].items():
                    ncbi_df.at[acc, tax] = name
                new_df = pd.concat([new_df, ncbi_df.loc[acc].to_frame().T])
            except AttributeError: # multiple entries
                eprint(f'{spacer}WARNING: {acc} is redundant', flush = True)
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
            except AttributeError: # multiple entries
                vprint(f'{spacer}\tWARNING: {check_acc} is redundant', 
                       v = verbose, e = True, flush = True)
                for acc1, row1 in ncbi_df.loc[acc].iterrows():
                    for file_t, file_p in acc2files[check_acc].items():
                        row1[file_t] = file_p
                    for tax, name in acc2org[check_acc].items():
                        row1[tax] = name
                    new_df = pd.concat([new_df, row1.to_frame().T])
                    

    if 'fna' in new_df.keys():
        new_df = new_df.rename(columns = {'fna': 'assemblyPath'})
    if 'gff3' in new_df.keys():
        new_df = new_df.rename(columns = {'gff3': 'gffPath'})
    new_df = new_df.reset_index()
    return new_df, rep_failed


def get_SRA(assembly_acc, fastqdump = 'fastq-dump', pe = True):

    handle = Entrez.esearch(db='SRA', term=assembly_acc)
    ids = Entrez.read(handle)['IdList']
    for id in ids:
        handle = Entrez.esummary(db='SRA', id = id, report='full')
        records = Entrez.read(handle, validate = False)
        for record in records:
            srr = re.search(r'Run acc="(S\w+\d+)"', record['Runs'])[1]
            print('\t\t' + srr, flush = True)
            cmd, count = 1, 0
            if pe:
                while cmd and count < 3:
                    count += 1
                    cmd = subprocess.call(['prefetch', srr, '--max-size', '10t'], 
                                           stdout = subprocess.PIPE)
                    if cmd:
                        continue
                    cmd = subprocess.call(['vdb-validate', srr], 
                                          stdout = subprocess.PIPE)
                    if cmd:
                        continue
                    cmd = subprocess.call([
                        fastqdump, '--split-3', '--gzip', srr], 
                        stdout = subprocess.PIPE)
                    if os.path.isfile(f'{srr}_1.fastq.gz'):
  #                  if os.path.isfile(srr + '_1.fastq'):
#                        cmd = subprocess.call(['gzip', f'{srr}_1.fastq'])
 #                       cmd = subprocess.call(['gzip', f'{srr}_2.fastq'])
                        shutil.move(f'{srr}_1.fastq.gz', 
                                    f'{assembly_acc}_{srr}_1.fq.gz')
                        shutil.move(f'{srr}_2.fastq.gz', 
                                    f'{assembly_acc}_{srr}_2.fq.gz')
                    else:
#                        cmd = subprocess.call(['gzip', f'{srr}.fastq'])
                        print('\t\t\tWARNING: file failed or not paired-end', flush = True)
            else:
                while cmd and count < 3:
                    count += 1
                    cmd = subprocess.call(['prefetch', srr], stdout = subprocess.PIPE)
                    if cmd:
                        continue
                    cmd = subprocess.call(['vdb-validate', srr], 
                                          stdout = subprocess.PIPE)
                    if cmd:
                        continue
                    cmd = subprocess.call([
                        fastqdump, srr, '--gzip'], 
                        stdout = subprocess.PIPE)
                    if cmd:
                        continue
#                    cmd = subprocess.call(['gzip', f'{srr}.fastq'])
                    if os.path.isfile(f'{srr}.fastq.gz'):
                        shutil.move(f'{srr}.fastq.gz', f'{assembly_acc}_{srr}.fq.gz')
                    else:
                        print('\t\t\tERROR: file failed', flush = True)


def goSRA(df, output = os.getcwd() + '/', pe = True, column = 'sra'):

    print()
    sra_dir = output + 'sra/'
    if not os.path.isdir(sra_dir):
        os.mkdir(sra_dir)
    os.chdir(sra_dir)
    fastqdump = findExecs('fastq-dump', exit = {'fastq-dump'})
    count = 0

    for i, row in df.iterrows():
        print('\t' + row[column], flush = True)
        get_SRA(row[column], fastqdump[0])
        count +=1
        if count >= 10:
            time.sleep(1)
            count = 0


def cli():
    parser = argparse.ArgumentParser(
        description = 'GenBank/RefSeq downloading utility. Downloads ' \
                    + 'accession by accession')
    parser.add_argument('-i', '--input', required = True, \
    help = 'Space delimited accession; tab delimited file with -c')
    parser.add_argument('-a', '--assembly', action = 'store_true')
    parser.add_argument('-p', '--proteome', action = 'store_true')
    parser.add_argument('-g', '--gff3', action = 'store_true')
    parser.add_argument('-t', '--transcript', action = 'store_true')
    parser.add_argument('-s', '--sra', action = 'store_true', \
        help = 'Download SRAs only' )
    parser.add_argument('-pe', '--paired', action = 'store_true', \
        help = 'Download paired-end SRAs. (REQUIRES -s)')
    parser.add_argument('-c', '--column', \
        help = 'Accession column num/name; DEFAULT ["assembly_acc" | 0]')
    parser.add_argument('-n', '--ncbi_column', \
        help = 'NCBI database associated with column. ' \
            + '{"assembly", "biosample", "bioproject", "genome" ...}; ' \
            + 'DEFAULT: attempt to decipher')
    parser.add_argument('-o', '--output', help = 'Output directory' )
    parser.add_argument('-e', '--email', help = 'NCBI email')
    parser.add_argument('--api', help = 'NCBI API key for high query rate')
    parser.add_argument('--fallback', action = 'store_true', 
        help = 'Fallback mode if datasets fails')
    args = parser.parse_args()

    if args.email:
        ncbi_email = args.email
        Entrez.email = ncbi_email
        if args.api:
            ncbi_api = args.api
            Entrez.api_key = ncbi_api
        else:
            ncbi_api = None
    else:
        ncbi_email, ncbi_api, jgi_email, jgi_pwd = loginCheck(jgi = False) 
        Entrez.email = ncbi_email
        if ncbi_api:
            Entrez.api_key = ncbi_api

    if not args.output:
        output = mkOutput(None, 'ncbiDwnld')
    else:
        output = format_path(args.output)

    findExecs('datasets', exit = {'datasets'})

    args_dict = {
        'NCBI Table': args.input,
        'email': ncbi_email,
        'Assemblies': args.assembly,
        'Proteomes': args.proteome,
        ".gff3's": args.gff3,
        'Transcripts': args.transcript,
        'SRA': args.sra}

    start_time = intro('Download NCBI files',args_dict)
    if args.sra:
        if os.path.isfile(format_path(args.input)):
            if not args.column:
                goSRA(pd.read_csv(format_path(args.input), sep = '\t',
                                  names = ['sra']), output, 
                      pe = args.paired, column = 'sra')
            else:
                goSRA(pd.read_csv(format_path(args.input), sep = '\t'), 
                      output, 
                      pe = args.paired, column = args.column)
        else:
            goSRA(pd.DataFrame({'sra': split_input(args.input)}), output, 
                  pe = args.paired, column = 'sra')
    else:
        if os.path.isfile(format_path(args.input)):
            ncbi_df = pd.read_csv(args.input, sep = '\t', header = None)
            if not args.column:
                if 'assembly_acc' in ncbi_df.keys():
                    column = 'assembly_acc'
                    ncbi_column = 'Assembly Accession'
                elif 'Assembly Accession' in ncbi_df.keys():
                    column = 'assembly_acc'
                    ncbi_column = 'Assembly Accession'
                else:
                    column = 0
            else:
                try:
                    column = ncbi_df.columns[int(0)]
                    ncbi_column = column
                except ValueError: # not an integer
                    pass
            if not args.ncbi_column:
                if args.column is not None:
                    if column.lower() in {'assembly'}:
                        ncbi_column = 'assembly'
                    elif column.lower() in \
                        {'genome', 'assembly accession', 'assembly_acc'}:
                        ncbi_column = 'genome'
                    elif column.lower() in {'biosample', 'biosample accession'}:
                        ncbi_column = 'biosample'
                    else:
                        ncbi_column = column.lower()
                else:
                    ncbi_column = 'genome'
            else:
                ncbi_column = args.ncbi_column.lower()
        else:
            ncbi_df = pd.DataFrame({'assembly_acc': split_input(args.input)})
            column = 'assembly_acc'
            ncbi_column = 'assembly'

        ncbi_df = ncbi_df.drop_duplicates(column)
        if args.fallback:
            from mycotools.ncbi_dwnld_fallback import main as main_fallback
            new_df, failed = main_fallback( 
                assembly = args.assembly, column = column, 
                ncbi_column = ncbi_column, proteome = args.proteome, 
                gff3 = args.gff3, transcript = args.transcript, ncbi_df = ncbi_df,
                output_path = output, verbose = True, spacer = ''
                )
        else:
            new_df, failed = main( 
                assembly = args.assembly, column = column, 
                ncbi_column = ncbi_column, proteome = args.proteome, 
                gff3 = args.gff3, transcript = args.transcript, ncbi_df = ncbi_df,
                output_path = output, verbose = True, spacer = ''
                )
        new_df = new_df.rename(columns = {'index': '#assembly_accession'})
        new_df['source'] = 'ncbi'
        new_df['useRestriction (yes/no)'] = 'no'
        if 0 in new_df.columns:
            del new_df[0]

        new_df.to_csv(output + 'ncbiDwnld.predb', sep = '\t', index = None)
        if failed:
            eprint('ERROR: ' + ','.join([str(x[0]) for x in failed]), flush = True)
           

    outro(start_time)


if __name__ == '__main__':
    cli()
