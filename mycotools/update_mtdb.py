#! /usr/bin/env python3

#NEED reinit implementation
#NEED to update introduction
#NEED a verbose option
#NEED revert version option (ome-by-ome/list of omes)
#NEED to reference a manually curated duplicate check
    # NEED a prohibit option to import prohibited JGI/NCBI IDs and option to update
#NEED to finish --save implementation
#NEED an only one of a strain feature in config
#NEED to pull failed JGI from NCBI
    # will require logging whatever NCBI omes directly overlap MycoCosm
#NEED to remove overlap when rerunning failed genomes

import os
import re
import sys
import time
import json
import base64
import shutil
import getpass
import hashlib
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
from mycotools.lib.dbtools import db2df, df2db, gather_taxonomy, assimilate_tax, \
    primaryDB, loginCheck, log_editor, mtdb, \
    mtdb_initialize
from mycotools.lib.kontools import intro, outro, format_path, eprint, \
                                   prep_output, collect_files, read_json, \
                                   write_json, split_input, findExecs
from mycotools.lib.biotools import fa2dict, gff2list, dict2fa, list2gff
from mycotools.ncbiDwnld import esearch_ncbi, esummary_ncbi, run_datasets, \
                                   compile_organism_names, main as ncbiDwnld
from mycotools.jgiDwnld import main as jgiDwnld
from mycotools.utils.ncbi2db import main as ncbi2db
from mycotools.utils.jgi2db import main as jgi2db
from mycotools.predb2mtdb import main as predb2mtdb
from mycotools.predb2mtdb import predb_headers, read_predb, gen_omes
from mycotools.assemblyStats import main as assStats
from mycotools.annotationStats import main as annStats


def validate_t_and_c(config, discrepancy = False):
    """Validate that user understands and accepts the conditions of
    use-restricted data, and that responsibility over confirming the
    desigination of use-restriction is the user's responsibility"""

    # if there is a configuration json, just query that
    if config and not discrepancy:
        try:
            if config['nonpublished'].lower() in {'yes', 'y'}:
                nonpublished = 'yes'
            else:
                nonpublished = False # weird situation for prokaryotes, how to
                # handle?
        except AttributeError:
            nonpublished = False

    # if there isnt a configuration, alert the user to use-restriction policies
    else:
        print('\nPlease review JGI use-restricted data policy here: ' \
            + 'https://jgi.doe.gov/user-programs/pmo-overview/policies/' \
            + '\nPlease review GenBank use-restricted data policy here: ' \
            + 'https://ncbi.nlm.nih.gov/genbank/'
            + '\nPlease review how Mycotools handles use-restricted data here:' \
            + ' https://github.com/xonq/mycotools/blob/master/MTDB.md', flush = True)
        check = ''
        if check.lower() not in {'y', 'yes'}:
            check = input('\nWARNING: This is a permanent configuration ' \
                  + 'option for the primary MTDB. Do you agree to honor ' \
                  + 'the terms and conditions ' \
                  + 'of use-restricted JGI and GenBank data; acknowledge that ' \
                  + 'Mycotools may not comprehensively determine use-restricted ' \
                  + 'designations; and acknowledge that you will validate ' \
                  + 'any use-restricted assignments in your local MTDB ' \
                  + 'prior to publication?' \
                  + '\n\nPlease type [y]es/[N]o if you acknowledge these terms: ')
        if check.lower() not in {'y', 'yes'}:
            print('\nRerun without --nonpublished', flush = True)
            sys.exit(1)

        nonpublished = 'yes'

    return nonpublished

def gen_config(
    branch = 'fungi', forbidden = '',
    repo = None,
    rogue = False, nonpublished = False, jgi = False,
    rank2lineages = {}
    ):

    config = {
        'forbidden': forbidden,
        'repository': repo,
        'branch': branch,
        'nonpublished': nonpublished,
        'rogue': rogue,
        'jgi': jgi,
        'lineage_constraints': rank2lineages
    }

    return config


def add_vars(init_dir, dbtype):
    """Initialize the environmental variables for MTDB"""
    mtdb_initialize(init_dir, dbtype, init = True)


def initDB(init_dir, branch, envs, dbtype, date = None, 
           rogue = False, nonpublished = False, jgi = True,
           repo = None, rank2lineages = {}):
    '''Initialize database in `init_dir`'''

    new_dirs = [
        init_dir + 'data/', init_dir + 'config/', 
        init_dir + 'data/fna/', init_dir + 'data/faa/',
        init_dir + 'mtdb',
        init_dir + 'data/gff3/', 
        init_dir + 'log/',
        init_dir + 'data/db/'
    ]
    output = prep_output(init_dir, cd = True)
    if not output.endswith('/'):
        output += '/'
    for new_dir in new_dirs:
        if not os.path.isdir(new_dir):
            os.mkdir(new_dir)

    config = gen_config(branch = branch, rogue = rogue, 
                        forbidden = '$MYCODB/log/forbidden.tsv',
                        nonpublished = nonpublished, jgi = jgi,
                        repo = repo, rank2lineages = rank2lineages)
    write_json(config, init_dir + 'config/mtdb.json', indent = 1)

    if not rogue:
        # this is a relic, and needs to be adjusted to a central reference if
        # that is ever created
        if not os.path.isdir( init_dir + 'mtdb' ):
# NEED TO CHANGE FROM SSH TO LINK ONCE OPEN (config['repository'])
            git_exit = subprocess.call( [
                'git', 'clone', 'git@gitlab.com:xonq/mtdb', init_dir + 'mtdb'
#                '-b', branch
                ] )
            if git_exit != 0:
                eprint('\nERROR: git clone failed.', flush = True)
                sys.exit(2)
        else:
            print('\nmycotoolsdb directory already exists', flush = True)
# NEED TO ADD GITIGNORE TO GIT
        if not primaryDB():
            eprint('\nERROR: no YYYYmmdd.mtdb in ' + format_path(envs['MYCODB']), flush = True)
            sys.exit(3)
    else:
        new_db_path = output + 'mtdb/' + date + '.mtdb'
        if not os.path.isfile(new_db_path):
            with open(output + 'mtdb/' + date + '.mtdb', 'w') as out:
                out.write(''.join(['\t' for x in mtdb.columns]))


    return output, config


def parse_forbidden(forbidden_path):
    """Read the log file containing information on forbidden genomes"""
    # NEED to be rewritten to reference a central repository forbidden file
    
    log_path = format_path(forbidden_path)
    if os.path.isfile(log_path):
        log_dict = readLog(log_path)
    else:
        log_dict = {}

    return log_dict


def add_forbidden(
    tag, source, file_path = None,
    flag = 'failed download'
    ):
    edit = tag + '\t' + source + '\t' + flag
    log_editor(file_path, tag, edit)


def parse_dups(file_path):
    """Retrieve a file containing replicated genomes and ignore these.
    This is important for dereplication of discrepant genus naming between NCBI
    and MycoCosm due to not adhering to conserved genus naming standards and
    updating relic genus names. It appears that MycoCosm will name genera by
    their anamorph occassionally, and not the consensus name - though I assume
    this is also to some extent present in NCBI. Ultimately, a manually curated
    file is necessary for this and should be held in a central repository."""
    duplicates = {}
    if os.path.isfile(file_path):
        with open(file_path, 'r') as raw:
            for line in raw:
                if not line.startswith('#'):
                    data = [x.rstrip() for x in line.split('\t') if x]
                    if data:
                        duplicates[data[0]] = [data[1], data[2], data[3]]

    return duplicates


#def add_dups(
 #   dup_code, dup_entry, file_path
  #  ):
#    edit = dup_code + '\t' + '\t'.join(dup_entry)
 #   log_editor(file_path, dup_code, edit)

def acq_forbid_omes(file_path):
    """Parse a file with forbidden ome accessions - ome codes that have been
    used before and are no longer valid"""
    if not os.path.isfile(file_path):
        return set()
    with open(file_path, 'r') as raw:
        relics = set([x.rstrip() for x in raw])
    return relics

def write_forbid_omes(omes, file_path):
    """Add to a file of forbidden ome accessions so that these are not ever
    used again, even if the codename is removed from the database"""
    if os.path.isfile(file_path):
        with open(file_path, 'r') as raw:
            old_relics = set([x.rstrip() for x in raw])
        new_relics = old_relics.union(set(omes))
    else:
        new_relics = set(omes)

    with open(file_path + '.tmp', 'w') as out: # be cautious because if it 
    # cancels then we lose the old data
        out.write('\n'.join([str(x) for x in sorted(new_relics)]))
    shutil.move(file_path + '.tmp', file_path)

def parse_failed(file_path = None, rerun = False):
    """Parse a file that stores the failed accessions and metadata of the
    attempted acquisition. Return a dictionary that contains the failed
    accession and its metadata."""
    prev_failed = {}
    if not os.path.isfile(file_path) or rerun:
        with open(file_path, 'w') as out:
            out.write('#code\tsource\tversion\tattempt_date')
    else:
        with open(file_path, 'r') as raw:
            for line in raw:
                if not line.startswith('#'):
                    data = [x.rstrip() for x in line.split('\t')]
                    while len(data) < 4:
                        data.append('')
                    prev_failed[data[0]] = {
                        'source': data[1], 'version': data[2],
                        'attempt_date': data[3]
                        }

    return prev_failed

def parse_jgi2ncbi(file_path):
    """Parse previously collected NCBI to JGI data to limit querying"""
    jgi2ncbi = {}
    if not os.path.isfile(file_path):
        with open(file_path, 'w') as out:
            out.write('#ncbi_acc\tmycocosm_portal')
    else:
        with open(file_path, 'r') as raw:
            for line in raw:
                if not line.startswith('#'):
                    d = line.rstrip().split('\t')
                    ncbi, jgi = d[0], d[1].lower()
                    jgi2ncbi[jgi] = ncbi

    return jgi2ncbi

def parse_true_ncbi(file_path):
    """Parse accessions considered to be unique to NCBI"""
    true_ncbi = set()
    if not os.path.isfile(file_path):
        with open(file_path, 'w') as out:
            out.write('#ncbi_acc')
    else:
        with open(file_path, 'r') as raw:
            true_ncbi = set([x.rstrip() for x in raw if not x.startswith('#')])

    return true_ncbi

def add_true_ncbi(true_ncbi, file_path = None):
    """Add to a ledger of accessions considered to be unique to NCBI"""
    with open(file_path, 'w') as out:
        out.write('#ncbi_acc\n' + '\n'.join([str(x) for x in list(true_ncbi)]))

def add_jgi2ncbi(jgi2ncbi, file_path = None):
    """Add to a ledger that seeks to associated NCBI accessions with JGI. This
    is prone to failure given that the field JGI uses to supply their genome
    accession is either absent from some NCBI entries, or is in a different
    field"""
    with open(file_path, 'w') as out:
        out.write('#ncbi_acc\tmycocosm_portal\n')
        for jgi, ncbi in jgi2ncbi.items():
            out.write(ncbi + '\t' + jgi + '\n')

def add_failed(
    code, source, version, date,
    file_path    
    ):
    """Add genomes to the failed acquisition file, including metadata on where
    it was downloaded, the version attempted to download, and the date of the
    run"""

    version = version.replace('-','').replace(' 00:00:00','') 
    edit = code + '\t' + source + '\t' + version + '\t' + date
    log_editor(file_path, code, edit)


def dwnld_mycocosm(
    out_file,
    mycocosm_url = 'https://mycocosm.jgi.doe.gov/ext-api/mycocosm/catalog/' \
                 + 'download-group?flt=&seq=all&pub=all&grp=fungi&srt=' \
                 + 'released&ord=desc'
    ):
    """Download the MycoCosm genome data spreadsheet, format to UTF-8 and
    return a Pandas dataframe of the data"""

    check_curl = findExecs(['curl'], verbose = False)

    if not os.path.isfile(out_file):
        for attempt in range(3):
            if check_curl:
                curl_cmd = subprocess.call(
                    ['curl', mycocosm_url, '-o', out_file + '.tmp'],
                    stdout = subprocess.PIPE
                    )
                if not curl_cmd:
                    shutil.move(out_file + '.tmp', out_file)
                    break
            if curl_cmd:
                eprint('\nERROR: failed to retrieve MycoCosm table', flush = True)
            else:
                resp = requests.get(url)
                with open(out_file + '.tmp', 'wb') as f:
                    f.write(resp.content)
                shutil.move(out_file + '.tmp', out_file)
                
    try:
        jgi_df = pd.read_csv(out_file, encoding = 'cp1252')
    except UnicodeDecodeError:
        jgi_df = pd.read_csv(out_file, encoding = 'latin1')
    except UnicodeDecodeError:
        jgi_df = pd.read_csv(out_file, encoding = 'utf-8')
    jgi_df.columns = [x.replace('"','').replace('"','') for x in jgi_df.columns]
    
    return jgi_df


def dwnld_ncbi_metadata( 
    ncbi_file,
    ncbi_url = 'https://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/',
    group = 'eukaryotes'
    ):
    """Download the NCBI genome data spreadsheet for prokaryotes or
    eukaryotes, and return a Pandas dataframe"""

    ncbi_url = ncbi_url + group + '.txt'
    if not os.path.isfile(ncbi_file):
        getTbl = subprocess.call( 
            ['curl', ncbi_url, '-o', ncbi_file + '.tmp']
            )
        shutil.move(ncbi_file + '.tmp', ncbi_file)
    ncbi_df = pd.read_csv(ncbi_file, sep = '\t')
 
    return ncbi_df

def prep_taxa_cols(df, taxonomy_dir, col = '#Organism/Name', api = None,
                   acc2org = {}, max_attempts = 3):

    skip_prep = list(acc2org.keys())
    gca_prep = [x.upper().replace('GCF', 'GCA') for x in skip_prep]
    gcf_prep = [x.upper().replace('GCA', 'GCF') for x in skip_prep]
    skip = set(gca_prep + gcf_prep)
    if not os.path.isdir(taxonomy_dir):
        os.mkdir(taxonomy_dir)
    aa_file = taxonomy_dir + 'assembly_accs.genbank.txt'
    with open(aa_file, 'w') as out:
        out.write('\n'.join([x for x in list(df['assembly_acc']) if x not in skip]))

    attempts, datasets_cmd = 0, 0
    if not os.path.isdir(taxonomy_dir + 'ncbi_dataset'):
        datasets_path = taxonomy_dir + 'ncbi_dataset.zip'
        while attempts < max_attempts:
            if attempts:
                eprint('\t\t\tReattempting', flush = True)
                if os.path.isfile(datasets_path):
                    os.remove(datasets_path)
            attempts += 1
            datasets_cmd = run_datasets(None, aa_file, taxonomy_dir, True, 
                                         api = api, verbose = True)
            try:
                with zipfile.ZipFile(datasets_path, 'r') as zip_ref:
                    zip_ref.extractall(taxonomy_dir)
                os.remove(datasets_path)
                break
            except zipfile.BadZipFile:
                eprint(f'\t\tERROR: datasets download corrupted - {attempts}',
                        flush = True)
                if attempts == max_attempts:
                    sys.exit(11)
            except FileNotFoundError:
                eprint(f'\t\tERROR: datasets failed - {attempts}',
                        flush = True)

    if datasets_cmd:
        eprint(f'\t\tWARNING: datasets failed, assuming no genomes found',
               flush = True)
        acc2org_n, acc2meta = {}, {}
    else:
        acc2org_n, acc2meta, org_failed = compile_organism_names(taxonomy_dir + 'ncbi_dataset/')
        print(f'\t\t{len(acc2meta) + len(org_failed)}', 
              'genomes queried from GenBank', flush = True)
        print(f'\t\t{len(org_failed)/len(df["assembly_acc"])*100}% failed', flush = True)
    

    # check for RefSeq for failed entries
    refseq_dir = taxonomy_dir + 'refseq/'
    if not os.path.isdir(refseq_dir):
        os.mkdir(refseq_dir)
    missing_accs = \
        sorted(set(df['assembly_acc']).difference(
                       set(acc2org_n.keys()).union(set(acc2org.keys()))))
    reattempt_acc = []
    for acc in missing_accs:
        if acc.upper().startswith('GCA'):
            reattempt_acc.append(acc.upper().replace('GCA_', 'GCF_'))
        elif acc.upper().startswith('GCF'):
            reattempt_acc.append(acc.upper().replace('GCF_', 'GCA_'))
    acc_file_re = refseq_dir + 'assembly_accs.refseq.txt'
    with open(acc_file_re, 'w') as out:
        out.write('\n'.join(reattempt_acc))

    # attempt to download the ncbi_datasets zip file until allowed attempts are
    # exhausted
    rs_datasets_cmd = 0
    if not os.path.isdir(refseq_dir + 'ncbi_dataset'):
        print(f'\t\tChecking RefSeq for {len(reattempt_acc)} entries', flush = True)
        rs_datasets_path = refseq_dir + 'ncbi_dataset.zip'
        attempts = 0
        while attempts < max_attempts:
            if attempts:
                eprint('\t\t\t\tReattempting', flush = True)
                if os.path.isfile(rs_datasets_path):
                    os.remove(rs_datasets_path)
            attempts += 1
            rs_datasets_cmd = run_datasets(None, acc_file_re, refseq_dir, True, 
                                   api = api, verbose = True)
            try:
                with zipfile.ZipFile(rs_datasets_path, 'r') as zip_ref:
                    zip_ref.extractall(refseq_dir)
                os.remove(rs_datasets_path)
                break
            except zipfile.BadZipFile:
                eprint(f'\t\t\tERROR: datasets download corrupted - {attempts}',
                        flush = True)
                if attempts == max_attempts:
                    sys.exit(10)
            except FileNotFoundError:
                eprint(f'\t\tERROR: datasets failed - {attempts}',
                        flush = True)

    if rs_datasets_cmd:
        eprint(f'\t\tWARNING: datasets failed, assuming no genomes found',
               flush = True)
        acc2org_rs, acc2meta_rs = {}, {}
    else:
        acc2org_rs, acc2meta_rs, org_failed_2 = compile_organism_names(refseq_dir + 'ncbi_dataset/')
        print(f'\t\t{len(acc2meta_rs)} genome(s) queried from RefSeq', flush = True)
   
    acc2org, acc2meta = {**acc2org, **acc2org_n, **acc2org_rs}, {**acc2meta, **acc2meta_rs}

    df['strain'] = ''
    todel = set()
    for i, row in df.iterrows():
        acc = row['assembly_acc']
        if acc in acc2org:
            df.at[i, 'genus'] = acc2org[acc]['genus']
            df.at[i, 'species'] = acc2org[acc]['species']
            df.at[i, 'strain'] = acc2org[acc]['strain']
        else:
            todel.add(i)

    df = df[~df.index.isin(todel)]

    return df, acc2meta, acc2org


def prep_jgi_cols(jgi_df, name_col = 'name'):

    jgi_df['publication(s)'] = jgi_df['publication(s)'].astype(str)
    for i, row in jgi_df.iterrows():
        taxa = re.sub(r'[^ a-zA-Z0-9\.]', '', row[name_col]).split()
        jgi_df.at[i, 'genus'] = taxa[0].replace('.','')
        jgi_df.at[i, 'publication(s)'] = jgi_df.at[i, 'publication(s)'].replace('"""','')
        if len(taxa) > 1:
            jgi_df.at[i, 'species'] = taxa[1]
            if len(taxa) > 2:
                jgi_df.at[i, 'strain'] = ''.join(taxa[2:])
                vers_search = re.search(r'v(\d+\.\d+$)', jgi_df['strain'][i])
                if vers_search is not None:
                    version = vers_search[1]
                    jgi_df.at[i, 'version'] = float(version)
                    jgi_df.at[i, 'strain'] = jgi_df['strain'][i][:-len(vers_search[0])]
                else:
                    jgi_df.at[i, 'version'] = 0
                jgi_df.at[i, 'strain'] = jgi_df.at[i, 'strain'].replace('.','')
        else:
            jgi_df.at[i, 'species'] = 'sp.'

    jgi_df = jgi_df.sort_values(by = 'version', ascending = False)
    jgi_df = jgi_df.drop_duplicates('portal')

    return jgi_df


def clean_ncbi_df(ncbi_df, update_path, kingdom = 'Fungi', 
                  api = None, max_attempts = 3):
    ncbi_df = ncbi_df.astype(str).replace(np.nan, '')

    acc2org_path = update_path + '../gca2org.tsv'
    acc2org = {}
    if os.path.isfile(acc2org_path):
        with open(acc2org_path, 'r') as raw:
            for line in raw:
                d = line.split('\t')
                acc2org[d[0]] = {'genus': d[1], 'species': d[2],
                                 'strain': d[3].rstrip()}

    if kingdom.lower() == 'fungi':
        # extract group of interest (case sensitive to first letter)
        ncbi_df = ncbi_df[ncbi_df['Group'] == 'Fungi']
    elif kingdom.lower() in {'plants', 'viridiplantae'}:
        ncbi_df = ncbi_df[ncbi_df['Group'] == 'Plants']
    elif kingdom.lower() in {'animals', 'metazoa'}:
        ncbi_df = ncbi_df[ncbi_df['Group'] == 'Animals']
    elif kingdom.lower() in {'bacteria', 'prokaryotes'}:
        ncbi_df = ncbi_df[ncbi_df['Group'] != 'Archaea']
    elif kingdom.lower() in {'archaea'}:
        ncbi_df = ncbi_df[ncbi_df['Group'] == 'Archaea']

    ncbi_df = ncbi_df[ncbi_df['assembly_acc'].str.startswith(('GCA', 'GCF'))]
    ncbi_df, acc2meta, acc2org = prep_taxa_cols(ncbi_df, 
                                       update_path + 'taxonomy/', 
                                       api = api, acc2org = acc2org)

    with open(acc2org_path + '.tmp', 'w') as out:
        for acc, org in acc2org.items():
            org_meta = f'{org["genus"]}\t{org["species"]}\t{org["strain"]}'
            out.write(f'{acc}\t{org_meta}\n')
    os.rename(acc2org_path + '.tmp', acc2org_path)

    # remove entries without sufficient metadata
    ncbi_df = ncbi_df.dropna(subset = ['genus'])

    # remove assembly- & annotation-lacking entries
    ncbi_df = ncbi_df[~ncbi_df['Genes'].isin({'-', '', '0'})]
    ncbi_df = ncbi_df[~ncbi_df['Proteins'].isin({'-', '' ,'0'})]
    ncbi_df = ncbi_df[~ncbi_df['assembly_acc'].isin({'-', '', '0'})]

    # sort by version, keep most recent version of duplicate assemblies
    ncbi_df['version'] = pd.to_datetime(ncbi_df['Modify Date'])
    ncbi_df = ncbi_df.sort_values(by = 'version', ascending = False)
    ncbi_df = ncbi_df.drop_duplicates('assembly_acc')

    # check for explicit new versions of assemblies
    ncbi_df['assembly_base'] = [
        x[:x.find('.')] for x in list(ncbi_df['assembly_acc'])
        ]
    ncbi_df['assembly_vers'] = [
        x[x.find('.'):] for x in list(ncbi_df['assembly_acc'])
        ]
    ncbi_df = ncbi_df.sort_values(by = 'assembly_vers', ascending = False)
    ncbi_df = ncbi_df.drop_duplicates('assembly_base')

    return ncbi_df, acc2meta


def exec_rm_overlap(ncbi_df, todel_i):
    """Execute the removal of JGI redundancy from NCBI genomes and create a
    DataFrame that represents the overlapping genomes"""
    todel = ncbi_df.index.intersection(todel_i)
    ncbi_jgi_overlap = pd.DataFrame(columns = ncbi_df.columns)
    ncbi_jgi_overlap = ncbi_df.loc[todel]
    ncbi_df = ncbi_df.drop(todel)
    return ncbi_jgi_overlap, ncbi_df


def rm_ncbi_overlap(ncbi_df, mycocosm_df, jgi2ncbi, 
                    fails = set(), acc2meta = {}, api = 3):
    """Acquire MycoCosm assembly accessions from NCBI.
    pd.DataFrame(ncbi_df) = post clean ncbi_df;
    set(mycocosm_omes) = set of lower-cased mycocosm genome codes;
    api = number of iterations before sleeping (3/10);
    First, obtain the genome UID by esearching via Entrez;
    Next, use the genome UID to obtain the genome summary via Entrez;
    Finally, relate the accession to the MycoCosm dataframe"""
    version_comp = re.compile(r'\s[Vv]\d+\.\d+')
    mycocosm_omes = set(
        [x.lower() for x in list(mycocosm_df['portal'])]
        )

    todel = []
    ncbi2jgi, jgi2biosample = {v: k for k, v in jgi2ncbi.items()}, {}
    ass_count = 0

    # could vectorize these
    jgi_names = {f"{v['genus']}_{v['species']}_{v['strain']}" \
                 for k, v in mycocosm_df.iterrows()}
    jgi_gen_sp = {f"{v['genus']}_{v['species']}" \
                 for k, v in mycocosm_df.iterrows()}
    for i, row in tqdm(ncbi_df.iterrows(), total = len(ncbi_df)):
        if row['assembly_acc'] in ncbi2jgi:
            jgi2biosample[ncbi2jgi[row['assembly_acc']]] = \
                row['BioSample Accession']
            todel.append(i)
#        elif row['assembly_acc'] in fails:
 #           pass
#            todel.append(i)
        elif row['assembly_acc'] not in fails:
            if row['assembly_acc'] in acc2meta:
                ass_name = acc2meta[row['assembly_acc']]['accession']
                submitter = acc2meta[row['assembly_acc']]['submitter']
            else:
                ass_uid = esearch_ncbi(row['assembly_acc'],
                                       'assembly', 'assembly')
                if not ass_uid:
                    fails.add(row['assembly_acc'])
                    continue
                ncbi_df.at[i, 'uid'] = ass_uid
                # this is a huge bottleneck, taking approximately 1s per query
                summary = esummary_ncbi(max(ass_uid), 'assembly')
                ass_name_prep = \
                    summary['DocumentSummarySet']['DocumentSummary'][0]['AssemblyName']
                ass_name = version_comp.sub('', ass_name_prep)
                submitter = \
                    summary['DocumentSummarySet']['DocumentSummary'][0]['SubmitterOrganization']
            if ass_name.lower() in mycocosm_omes:
                todel.append(i)
                jgi2ncbi[ass_name.lower()] = row['assembly_acc']
                jgi2biosample[ass_name.lower()] = row['BioSample Accession']
            elif any(x in submitter.lower() for x in [
                'joint genome institute', 'jgi', 'joint_genome_institute'
                ]): # very crude, but mycocosm does not give the option to be
                # systematic
                todel.append(i)
                jgi2ncbi[ass_name.lower() + f'${ass_count}'] = row['assembly_acc']
                jgi2biosample[ass_name.lower() + f'${ass_count}'] = \
                    row['BioSample Accession']
                ass_count += 1
            elif f"{row['genus']}_{row['species']}_{row['strain']}" \
                 in jgi_names:
                todel.append(i)
                jgi2ncbi[f"{row['genus']}_{row['species']}_{row['strain']}"] = \
                    row['assembly_acc']
                jgi2biosample[f"{row['genus']}_{row['species']}_{row['strain']}_{ass_count}"] = \
                    row['BioSample Accession']
                ass_count += 1
            elif not row['strain'].rstrip() \
                 and f'{row["genus"]}_{row["species"]}' in jgi_gen_sp:
                # unfortunately we cannot validate that this is either a JGI or
                # NCBI unique specimen, but to err on the side of caution we
                # should remove it. 
                todel.append(i)
                gen_sp = f'{row["genus"]}_{row["species"]}'
                jgi2ncbi[gen_sp] = row['assembly_acc']
                jgi2biosample[gen_sp] = row['BioSample Accession']
            else:
                fails.add(row['assembly_acc'])
#    for i in reversed(todel):
 #       ncbi_jgi_overlap = pd.concat([ncbi_jgi_overlap, ncbi_df.loc[i]])
#        ncbi_df = ncbi_df.drop(i)

    ncbi_df, ncbi_jgi_overlap = exec_rm_overlap(ncbi_df, todel)

    return ncbi_df, jgi2ncbi, jgi2biosample, fails, ncbi_jgi_overlap, todel


def mk_wrk_dirs(update_path):
    """Make the download directories in the update path"""
    wrk_dirs = ['faa/', 'fna/', 'gff3/']
    for wrk_dir in wrk_dirs:
        if not os.path.isdir(update_path + wrk_dir):
            os.mkdir(update_path + wrk_dir)


def prepare_ref_db(ref_db, date):
    """Prepare MTDBs of the data in the reference MTDB that are from the two
    download sources"""
    ref_db['acquisition_date'] = [date for x in ref_db['acquisition_date']]
    ref_db = ref_db.set_index()
    
    jgi = mtdb({k: v for k,v in ref_db.items() \
                if v['source'].lower() == 'jgi'}, 
               index = 'ome')
    ncbi = mtdb({k: v for k,v in ref_db.items() \
                if v['source'].lower() == 'ncbi'},
               index = 'ome')
    if set(ref_db.keys()).difference(set(jgi.keys()).union(set(ncbi.keys()))):
        eprint('\tWARNING: reference entries that are not labeled "jgi/ncbi" are excluded')

    return jgi.mtdb2pd(), ncbi.mtdb2pd()


def internal_redundancy_check(db):
    """Check the inputted database for overlapping assembly accessions and
    dereplicate, including those with different versions of the same
    accession"""
    db = mtdb.pd2mtdb(db).set_index()
    ncbi_db = {k: v for k,v in db.items() if v['source'] == 'ncbi'}
    jgi_db = {k: v for k,v in db.items() if v['source'] == 'jgi'}

    ncbi_accs = defaultdict(list)
    for ome, row in ncbi_db.items():
        ass_acc = row['assembly_acc'][:row['assembly_acc'].find('.')]
        ncbi_accs[ass_acc].append(ome)
    red_ncbi = {k: v for k, v in ncbi_accs.items() if len(v) > 1}
    for ass_acc, omes in red_ncbi.items():
        omes = sorted(omes, reverse = True) # large ome number to small
        accs = defaultdict(list)
        for ome in omes:
            try:
                acc_ver = int(
                    ncbi_db[ome]['assembly_acc'][ncbi_db[ome]['assembly_acc'].find('.')+1:]
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
        ass_acc = row['assembly_acc']
        jgi_accs[ass_acc].append(ome)
    red_jgi = {k: v for k,v in jgi_accs.items() if len(v) > 1}
    for ass_acc, omes in red_jgi.items():
        omes = sorted(omes, reverse = True)
        accs = defaultdict(list)
        for ome in omes:
            acc_ver = jgi_db[ome]['version'].replace('v','').replace('V','')
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
    if os.path.isfile(tax_path):
        with open(tax_path, 'r') as raw:
            for line in raw:
                data = line.rstrip().split('\t')
                tax_dicts[data[0]] = json.loads(data[1])
    return tax_dicts


def ref_update(
    ref_db, update_path, date, rerun, jgi_email, jgi_pwd,
    config, ncbi_email, ncbi_api, cpus = 1, check_MD5 = True,
    jgi = True, group = 'eukaryotes', kingdom = 'Fungi',
    remove = True, taxonomy = True, ncbi_fallback = False
    ):
    """Initialize/Update the primary MTDB based on a reference database
    acquired external from any existing primary MTDB"""
# NEED to mark none for new databases' refdb
    # initialize update
    print('\nInitializing run', flush = True)
    mk_wrk_dirs(update_path)

    jgi_df, ncbi_df = prepare_ref_db(ref_db, date)


    # run JGI
    if jgi and len(jgi_df) > 0:
        print('\nAssimilating MycoCosm', flush = True)
        jgi_db_path = update_path + date + '.jgi.mtdb'
        jgi_predb_path = update_path + date + '.jgi.predb2.mtdb'

        if not os.path.isfile(jgi_predb_path):
            print('\tDownloading MycoCosm data', flush = True)
            post_jgi_df, jgi_failed = jgiDwnld(jgi_df, update_path, jgi_email, jgi_pwd)
            jgi_predb = post_jgi_df.rename(columns = {'published(s)': 'published',
                                                      'fna_path': 'assemblyPath', 
                                                      'gff3_path': 'gffPath'})


            print('\tCurating MycoCosm data', flush = True)
            jgi_premtdb = jgi_predb.fillna('').to_dict(orient='list')
            jgi_mtdb, jgi_failed1 = predb2mtdb(jgi_premtdb, mtdb(), update_path,
#                                            forbidden = forbid_omes, 
                                            cpus = cpus, 
                                            remove = remove, spacer = '\t\t')
            jgi_failed = list(jgi_failed)
            jgi_failed.extend(jgi_failed1)
            jgi_mtdb.df2db(jgi_predb_path)
            for failure in jgi_failed:
                add_failed(failure[0], 'jgi', str(failure[1]), date,
                          format_path('$MYCODB/../log/failed.tsv'))

        else:
            jgi_mtdb = mtdb(jgi_predb_path)

    else:
        jgi_mtdb = mtdb()
    new_db = jgi_mtdb.mtdb2pd()

    print('\nAssimilating NCBI', flush = True)
    if not os.path.isfile(update_path + date + '.ncbi.predb'):
        print('\tDownloading NCBI data', flush = True)
        if ncbi_fallback:
            from mycotools.ncbi_dwnld_fallback \
                import main as ncbi_dwnld_fallback
            ncbi_predb, ncbi_failed1 = ncbi_dwnld_fallback(
                assembly = True, proteome = False, gff3 = True,
                ncbi_df = ncbi_df, remove = True, output_path = update_path,
                column = 'assembly_acc', ncbi_column = 'genome', 
                check_MD5 = check_MD5, verbose = True
                )

        else:
            ncbi_predb, ncbi_failed1 = ncbiDwnld(
                assembly = True, proteome = False, gff3 = True,
                ncbi_df = ncbi_df, remove = True, output_path = update_path,
                column = 'assembly_acc', ncbi_column = 'genome', 
                check_MD5 = check_MD5, verbose = True
                )

        for failure in ncbi_failed1:
            add_failed(failure[0], 'ncbi', str(failure[1]), date,
                      format_path('$MYCODB/../log/failed.tsv'))

#        for dup in new_dups:
 #           add_dups(dup, new_dups[dup], format_path('$MYCODB/../log/duplicates.tsv'))
        ncbi_predb.to_csv(update_path + date + '.ncbi.predb',
                          sep = '\t', index = None)
    else:
#        refdbncbi = mtdb(update_path + date + '.ncbi.ref.mtdb')
        ncbi_predb = pd.read_csv(update_path + date + '.ncbi.predb',
                                 sep = '\t')

    print('\tCurating NCBI data', flush = True)
    if not os.path.isfile(update_path + date + '.ncbi.predb2.mtdb'):
        for key in ncbi_predb.columns:
            ncbi_predb[key] = ncbi_predb[key].fillna('')
        ncbi_predb['version'] = ncbi_predb['version'].astype(str)
        ncbi_premtdb = ncbi_predb.to_dict(orient='list')
        ncbi_mtdb, ncbi_failed2 = predb2mtdb(ncbi_premtdb, mtdb(), update_path,
    #                                       forbidden = forbid_omes, 
                                           cpus = cpus,
                                           remove = remove, spacer = '\t\t')
        for failure in ncbi_failed2:
            add_failed(failure[0], 'ncbi', str(failure[1]), date,
                      format_path('$MYCODB/../log/failed.tsv'))
        ncbi_mtdb.df2db(update_path + date + '.ncbi.predb2.mtdb')
    else:
        ncbi_mtdb = mtdb(f'{update_path}{date}.ncbi.predb2.mtdb')

    try:
        ncbi_db = db2df(update_path + date + '.ncbi.predb2.mtdb')
    except pd.errors.EmptyDataError:
        ncbi_db = pd.DataFrame({x: [] for x in refdbncbi.keys()})
    if len(ncbi_db) > 0:
#        df2db(ncbi_db, ncbi_db_path)
        new_db = pd.concat([new_db, ncbi_db])

    print('\nAssimilating NCBI taxonomy data', flush = True)
    new_mtdb = mtdb.pd2mtdb(new_db)

    if kingdom.lower() == 'fungi':
        rank = 'kingdom'
    else:
        rank = 'superkingdom'

    if jgi_mtdb and ncbi_mtdb:
        update_mtdb = mtdb({**jgi_mtdb.set_index(), **ncbi_mtdb.set_index()}, 
                            index = 'ome')
    elif ncbi_mtdb:
        update_mtdb = ncbi_mtdb
    elif jgi_mtdb:
        update_mtdb = jgi_mtdb
    else:
        eprint('\nNo updates', flush = True)
        sys.exit(0)

    if taxonomy: #already completed
        tax_path = f'{update_path}../taxonomy.tsv'
        tax_dicts = read_prev_tax(tax_path)
        tax_dicts = gather_taxonomy(new_mtdb, api_key = ncbi_api, 
                                king=kingdom, rank = rank, 
                                output_path = tax_path, tax_dicts = tax_dicts)
        new_mtdb, genus_dicts = assimilate_tax(new_mtdb, tax_dicts) 
        dupFiles = {'fna': {}, 'faa': {}, 'gff3': {}}
    
        for ome, row in update_mtdb.items():
            if row['genus'] in genus_dicts:
                row['taxonomy'] = genus_dicts[row['genus']]
    
    return new_mtdb, update_mtdb 


def extract_constraint_lineages(df, ncbi_api, kingdom,
                                lineage_constraints, tax_dicts,
                                tax_path):
    """Extract genera from NCBI and JGI Pandas dataframes 
    that hit a dictionary of lineages of interest"""

    # begin extracting lineages of interest and store tax_dicts for later
    if 'taxonomy' not in df.columns:
        df['taxonomy'] = [{} for x in range(len(df))]
    if kingdom.lower() in {'fungi', 'plants', 'animals', 'metazoa',
                           'viridiplantae'}:
        query_rank = 'kingdom'
    else:
        query_rank = 'superkingdom'

    # skip querying ncbi if we are only gathering the genus
    if not set(lineage_constraints.keys()).difference({'genus'}):
        passing_tax = set(x[0].upper() + x[1:] \
                          for x in lineage_constraints['genus'])
        df = df[df['genus'].isin(passing_tax)]
        return tax_dicts, df

    tax_dicts = gather_taxonomy(df, api_key = ncbi_api,
                                king = kingdom, rank = query_rank,
                                tax_dicts = tax_dicts, output_path = tax_path)

    lineage_constraints = {k: set(v) for k, v in lineage_constraints.items()}

    # extract genera that pass
    passing_tax = set()
    if 'genus' in lineage_constraints:
        # add constraint genera immediately
        passing_tax = set(x[0].upper() + x[1:] \
                          for x in lineage_constraints['genus'])
        lineage_constraints = {k: v for k, v in lineage_constraints.items() \
                               if k != 'genus'}

    for genus, tax in tax_dicts.items():
        for rank, lineages in lineage_constraints.items():
            if rank in tax:
                if tax[rank].lower() in lineages:
                    passing_tax.add(genus)

    df = df[df['genus'].isin(passing_tax)]
    return tax_dicts, df


def taxonomy_update(orig_db, update_path, date, 
                    config, ncbi_email, ncbi_api,
                    rank = 'kingdom', group = 'fungi'):
    """Reset the taxonomy for the entire database and overwrite the previous
    tax path data to accomodate new taxonomy"""
    taxless_db = orig_db.reset_index()
    taxless_db['taxonomy'] = [{} for x in taxless_db['taxonomy']]
    tax_path = f'{update_path}../taxonomy.tsv'
    gca_path = f'{update_path}../gca2org.tsv'
    if os.path.isfile(tax_path):
        os.rename(tax_path, update_path + 'old_taxonomy.tsv')
    if os.path.isfile(gca_path):
        os.rename(gca_path, update_path + 'old_gca2org.tsv')
    tax_dicts = gather_taxonomy(taxless_db, api_key = ncbi_api,
                                    king=group, rank = rank, 
                                    output_path = tax_path)
    tax_db, genus_dicts = assimilate_tax(taxless_db, tax_dicts)
    if not isinstance(tax_db, mtdb):
        return tax_db, mtdb.pd2mtdb(tax_db)
    else:
        return tax_db.mtdb2pd(), tax_db



def rogue_update(
    db, update_path, date, rerun, jgi_email, jgi_pwd,
    config, ncbi_email, ncbi_api, cpus = 1, check_MD5 = True,
    jgi = True, group = 'eukaryotes', kingdom = 'Fungi',
    remove = True, lineage_constraints = {}, ncbi_fallback = False
    ):
    """Initialize/update a standalone primary MTDB"""
# NEED to mark none for new databases' refdb
    # initialize update
    print('\nInitializing run', flush = True)
    mk_wrk_dirs(update_path)
    prev_failed = parse_failed(rerun = rerun, file_path = format_path('$MYCODB/../log/failed.tsv'))
    duplicates = parse_dups(format_path('$MYCODB/../log/duplicates.tsv'))
    forbid_omes = acq_forbid_omes(
        file_path = format_path('$MYCODB/../log/relics.txt')
        )

    db = internal_redundancy_check(db)

    # prepare ncbi_df
    if ncbi_api:
        api = 10
    else:
        api = 3
    ncbi_db_path = update_path + date + '.ncbi.mtdb'
    pre_ncbi_df0 = dwnld_ncbi_metadata(update_path + date + '.ncbi.tsv',
                                   group = group) 
    pre_ncbi_df1 = pre_ncbi_df0.rename(columns={'Assembly Accession': 'assembly_acc'})
    print('\tAcquiring NCBI metadata', flush = True)
    ncbi_df, acc2meta = clean_ncbi_df(pre_ncbi_df1, update_path, 
                            kingdom = kingdom, api = ncbi_api)

    # begin extracting lineages of interest and store tax_dicts for later
    tax_path = f'{update_path}../taxonomy.tsv'
    tax_dicts = read_prev_tax(tax_path)
#    tax_dicts = {v['genus']: v['taxonomy'] for k, v in db.iterrows() \
 #                if any(y for x, y in v['taxonomy'].items() \
  #                      if x not in {'genus', 'species', 'strain'})}
    if lineage_constraints:
        lineage_path = update_path + date + '.ncbi.posttax.df'
        if not os.path.isfile(lineage_path):
            print('\nExtracting lineages from NCBI', flush = True)
            # NEED to transition to datasets
            tax_dicts, ncbi_df = extract_constraint_lineages(ncbi_df, 
                                                      ncbi_api, kingdom,
                                                      lineage_constraints,
                                                      tax_dicts, tax_path)
            ncbi_df.to_csv(lineage_path, sep = '\t', index = None)
        else:
            ncbi_df = pd.read_csv(lineage_path, sep = '\t')

    old_len = len(db['ome'])
    new_len = len(db['ome'])
    if old_len - new_len:
        print('\t' + str(old_len-new_len) + ' redundant entries removed',
              flush = True)

    # run JGI
    if jgi:
        print('\nAssimilating MycoCosm (1 download/minute)', flush = True)
        jgi_db_path = update_path + date + '.jgi.mtdb'
        mycocosm_path = update_path + date + '.mycocosm.csv'

        # acquire the mycocosm master table
        jgi_df = dwnld_mycocosm(mycocosm_path)
        jgi_df['biosample'] = ''
        jgi_df = prep_jgi_cols(jgi_df, 'name')

        # extract JGI lineages of interest and store tax_dicts for later
        if lineage_constraints:
            lineage_path = update_path + date + '.jgi.posttax.df'
            if not os.path.isfile(lineage_path):
                print('\tExtracting lineages from MycoCosm', flush = True)
                tax_dicts, jgi_df = extract_constraint_lineages(jgi_df, 
                                                        ncbi_api, kingdom,
                                                        lineage_constraints,
                                                        tax_dicts, tax_path)
                jgi_df.to_csv(lineage_path, sep = '\t', index = None)
            else:
                jgi_df = pd.read_csv(lineage_path, sep = '\t')

        print('\tSearching NCBI for MycoCosm overlap', flush = True)
        jgi_ncbi_overlap_file = f'{update_path}/redundant_ncbi.tsv'
        jgi2ncbi = parse_jgi2ncbi(update_path + '../jgi2ncbi.tsv')
        ncbi_df = ncbi_df.set_index('assembly_acc', drop = False)
        if os.path.isfile(jgi_ncbi_overlap_file):
            with open(jgi_ncbi_overlap_file, 'r') as raw:
                todel_i = [x.rstrip() for x in raw]
            ncbi_df, ncbi_jgi_overlap = exec_rm_overlap(ncbi_df, todel_i)
        else:
            true_ncbi = parse_true_ncbi(update_path + '../supported_ncbi.tsv')
            ncbi_df, jgi2ncbi, jgi2biosample, \
            true_ncbi, ncbi_jgi_overlap, todel_i \
                = rm_ncbi_overlap(ncbi_df, jgi_df, jgi2ncbi, 
                                  true_ncbi, acc2meta, api = api)

            print('\t\t' + str(len(jgi2ncbi)) + ' overlapping genomes',
                 flush = True)
            with open(jgi_ncbi_overlap_file + '.tmp', 'w') as out:
                out.write('\n'.join([x for x in todel_i]))
            os.rename(jgi_ncbi_overlap_file + '.tmp', jgi_ncbi_overlap_file)
            add_true_ncbi(true_ncbi, update_path + '../supported_ncbi.tsv')
            add_jgi2ncbi(jgi2ncbi, update_path + '../jgi2ncbi.tsv')
            for i, row in jgi_df.iterrows():
                if row['portal'].lower() in jgi2ncbi \
                    and row['portal'].lower() in jgi2biosample:
                    jgi_df.at[i, 'biosample'] = jgi2biosample[
                        row['portal'].lower()
                        ]

        print('\tDownloading MycoCosm data', flush = True)
        jgi_predb_path = update_path + date + '.jgi.predb2.mtdb'
        jgi_predb, db, jgi_failed = jgi2db( 
            jgi_df, db, update_path, jgi_email, jgi_pwd,
            date = date, nonpublished = config['nonpublished'],
            rerun = rerun, failed_dict = prev_failed,
            jgi2ncbi = jgi2ncbi, repeatmasked = True
            ) # download JGI files and ready predb

        # get ncbi hits that hit failed jgi runs
        failed_ncbi2jgi = {jgi2ncbi[f[0]]: f[0] \
                           for f in jgi_failed \
                           if f[0] in jgi2ncbi}
        ncbi_jgi_overlap = \
            ncbi_jgi_overlap[ncbi_jgi_overlap['assembly_acc'].isin(failed_ncbi2jgi)]
        ncbi_df = pd.concat([ncbi_df, ncbi_jgi_overlap])

        refdbjgi = mtdb.pd2mtdb(db)
        if not os.path.isfile(jgi_predb_path):
            print('\tCurating MycoCosm data', flush = True)
            jgi_premtdb = jgi_predb.fillna('').to_dict(orient='list')
            if 'assemblyPath' in jgi_premtdb:
                jgi_mtdb, jgi_failed1 = predb2mtdb(jgi_premtdb, refdbjgi, update_path,
                                                forbidden = forbid_omes, cpus = cpus, 
                                                remove = remove, spacer = '\t\t')
                jgi_failed.extend(jgi_failed1)
            else:
                jgi_mtdb = mtdb()
            jgi_mtdb.df2db(jgi_predb_path)

            for failure in jgi_failed:
                add_failed(failure[0], 'jgi', str(failure[1]), date,
                          format_path('$MYCODB/../log/failed.tsv'))
        else:
            jgi_mtdb =  mtdb(jgi_predb_path)

        try:
            jgi_db = db2df(jgi_predb_path)
        except pd.errors.EmptyDataError: # empty JGI predb
            jgi_db = pd.DataFrame({x: [] for x in refdbjgi.keys()})
    
        new_db_path = update_path + date + '.checkpoint.jgi.mtdb'
        if not os.path.isfile(new_db_path): 
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
        new_dups = duplicates

    print('\nAssimilating NCBI (10 download/second w/API key, 3 w/o)', flush = True)
    new_db['version'] = new_db['version'].astype(str)
    if not os.path.isfile(update_path + date + '.ncbi.predb'):
#    if not os.path.isfile(update_path + date + '.ncbi.predb'):
        print('\tDownloading NCBI data', flush = True)
        ncbi_predb, new_db, ncbi_failed1 = ncbi2db( 
            update_path, ncbi_df, ref_db = new_db, 
            date = date, failed_dict = prev_failed, 
            rerun = rerun, duplicates = duplicates,
            check_MD5 = check_MD5, fallback = ncbi_fallback
            )

        for failure in ncbi_failed1:
            add_failed(failure[0], 'ncbi', str(failure[1]), date,
                      format_path('$MYCODB/../log/failed.tsv'))
#        for dup in new_dups:
 #           add_dups(dup, new_dups[dup], format_path('$MYCODB/../log/duplicates.tsv'))
        refdbncbi = mtdb.pd2mtdb(new_db)
        refdbncbi.df2db(update_path + date + '.ncbi.ref.mtdb')
        ncbi_predb.to_csv(update_path + date + '.ncbi.predb',
                          sep = '\t', index = None)
    else:
        refdbncbi = mtdb(update_path + date + '.ncbi.ref.mtdb')
        ncbi_predb = pd.read_csv(update_path + date + '.ncbi.predb',
                                 sep = '\t')

    print('\tCurating NCBI data', flush = True)
    if not os.path.isfile(update_path + date + '.ncbi.predb2.mtdb'):
        for key in ncbi_predb.columns:
            ncbi_predb[key] = ncbi_predb[key].fillna('')
        ncbi_predb['version'] = ncbi_predb['version'].astype(str)
        ncbi_premtdb = ncbi_predb.to_dict(orient='list')
        for col in predb_headers:
            if col not in ncbi_premtdb:
                ncbi_premtdb[col] = ['' for x in ncbi_predb[key]]
        ncbi_premtdb['restriction'] = ['0' for x in ncbi_premtdb['assemblyPath']]
        ncbi_mtdb, ncbi_failed2 = predb2mtdb(ncbi_premtdb, refdbncbi, update_path,
                                           forbidden = forbid_omes, cpus = cpus,
                                           remove = remove, spacer = '\t\t')
        for failure in ncbi_failed2:
            add_failed(failure[0], 'ncbi', str(failure[1]), date,
                      format_path('$MYCODB/../log/failed.tsv'))
        ncbi_mtdb.df2db(update_path + date + '.ncbi.predb2.mtdb')

    try:
        ncbi_db = db2df(update_path + date + '.ncbi.predb2.mtdb')
        ncbi_mtdb = mtdb(update_path + date + '.ncbi.predb2.mtdb')
    except pd.errors.EmptyDataError:
        ncbi_db = pd.DataFrame({x: [] for x in refdbncbi.keys()})
        ncbi_mtdb = mtdb()
    if len(ncbi_db) > 0:
        df2db(ncbi_db, ncbi_db_path)
        new_db = pd.concat([new_db, ncbi_db])

    new_mtdb = mtdb.pd2mtdb(new_db)

    print('\nAssimilating NCBI taxonomy data', flush = True)
    if kingdom.lower() == 'fungi':
        rank = 'kingdom'
    else:
        rank = 'superkingdom'

    tax_path = f'{update_path}../taxonomy.tsv'
    tax_dicts = gather_taxonomy(new_mtdb, api_key = ncbi_api, 
                                king=kingdom, rank = rank, 
                                tax_dicts = tax_dicts,
                                output_path = tax_path)
    new_mtdb, genus_dicts = assimilate_tax(new_mtdb, tax_dicts) 
    dupFiles = {'fna': {}, 'faa': {}, 'gff3': {}}

    if jgi_mtdb and ncbi_mtdb:
        update_mtdb = mtdb({**jgi_mtdb.set_index(), **ncbi_mtdb.set_index()}, 
                            index = 'ome')
    elif ncbi_mtdb:
        update_mtdb = ncbi_mtdb
    elif jgi_mtdb:
        update_mtdb = jgi_mtdb
    else:
        eprint('\nNo updates', flush = True)
        sys.exit(0)

    return new_mtdb, update_mtdb 


def rm_raw_data(out_dir):
    """Remove raw data after completion"""
    for i in ['faa', 'gff3', 'gff', 'xml', 'fna']:
        if os.path.isdir(out_dir + i):
            shutil.rmtree(out_dir + i)


def gen_algn_db(update_path, omes):
    """Generate an alignment database for the complete primary MTDB"""
    date = os.path.basename(os.path.abspath(update_path))
    fas = collect_files(os.environ['MYCOFAA'] + '/', '.faa')
    fas = [x for x in fas if os.path.basename(x)[:-6] in omes]
    mkdb_base = 'cat ' + ' '.join(fas) 
    mkdb_blast = mkdb_base + ' | makeblastdb -in -' + \
        ' -out ' + os.environ['MYCOGFF3'] + '../db/' + date + \
        '.db -parse_seqids -dbtype prot -title ' + date + '.db'
  #  mkdb_mmseqs = mkdb_base + ' | mmseqs createdb stdin ' + \
 #       format_path('$MYCOFAA/' + date + '.mmseqs.db') + '; ' + \
#        'mmseqs createdb ' + format_path('$MYCOFAA/' + date + \
   #     '.mmseqs.db') + ' tmp'
    with open(update_path + date + '_makeblastdb.sh', 'w') as out:
        out.write(mkdb_blast)
    #with open(update_path + date + '_mmseqsdb.sh', 'w') as out:
     #   out.write(mkdb_mmseqs)

    print('\nOPTIONAL: To generate blastdb | mmseqsdb, run the following' + \
        '\nbash ' + update_path + date + '_makeblastdb.sh')
     #bash ' + update_path \
     #   + date + '_mmseqsdb.sh')
  

def check_add_mtdb(orig_mtdb, add_mtdb, update_path, overwrite = True):
    """Check the original MTDB for overlapping omes and curate if necessary"""
    orig_mtdb = orig_mtdb.set_index()
    add_mtdb = add_mtdb.set_index()
    orig_aa2ome = {v['assembly_acc']: k for k, v in orig_mtdb.items()}

#    failed_aas = []
    overwrite_omes = []
    for assembly_acc, row in add_mtdb.items():
        if assembly_acc in orig_aa2ome:
 #           failed_aas.append(assembly_acc)
            if overwrite:
                overwrite_omes.append(orig_aa2ome[assembly_acc])
            else:
                overwrite_omes.append(row['ome'])

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
        need_ome_mtdb = mtdb({k: add_mtdb[k] for k in sorted(inter_omes)},
                             index = 'ome')
        # these genomes have unique codenames and can be submitted as-is
        fine_ome_mtdb = mtdb({k: v for k, v in add_mtdb.items() \
                         if k not in inter_omes}, index = 'ome')
        need_ome_mtdb = need_ome_mtdb.reset_index()
        # delete the previous ome codes of those that need new ones
        need_ome_mtdb['ome'] = ['' for x in need_ome_mtdb['assembly_acc']]
        # the reference database should now contain the original and fine
        # codenames
        ref_ome_mtdb = mtdb({**fine_ome_mtdb, **orig_mtdb}, index = 'ome')
        # generate new codenames

        part_ome_mtdb, failed = gen_omes(need_ome_mtdb, ref_ome_mtdb.reset_index())
        overlap_mtdb = add_mtdb.set_index('assembly_acc')

        # prepare a database of new genome codes and previously fine ones
        new_ome_mtdb = mtdb({**part_ome_mtdb.set_index('ome'), **fine_ome_mtdb}, 
                            index = 'ome')
        # choose the assembly accession column to reference for old names
        new_ome_mtdb = new_ome_mtdb.set_index('assembly_acc')
        old_ome2new_ome = {v['ome']: new_ome_mtdb[k]['ome'] \
                           for k, v in overlap_mtdb.items() \
                           if v['ome'] in inter_omes}
        new_ome2old_ome = {v: k for k, v in old_ome2new_ome.items()}
        new_ome_mtdb = new_ome_mtdb.set_index('ome')
        for k, v in old_ome2new_ome.items():
            print(f'\t{k} converted to {v}',  flush = True)

        # create directories for new files
        fna_dir, gff_dir, faa_dir = f'{update_path}fna/', \
                                    f'{update_path}gff3/', \
                                    f'{update_path}faa/'
        for path_ in [fna_dir, gff_dir, faa_dir]:
            if not os.path.isdir(path_):
                os.mkdir(path_)

        # convert the file header names to the new omes
        for ome, old_ome in new_ome2old_ome.items():
            row = new_ome_mtdb[ome]
            if row['assembly_acc'] == old_ome:
                new_ome_mtdb[ome]['assembly_acc'] = ome

            old_ome = new_ome2old_ome[ome]
            fna = fa2dict(row['fna'])
            new_fna = {k.replace(old_ome + '_', ome + '_'): v \
                       for k, v in fna.items()}
            with open(f'{fna_dir}{ome}.fna', 'w') as out:
                out.write(dict2fa(new_fna))
            gff = gff2list(row['gff3'])
            for entry in gff:
                entry['seqid'] = entry['seqid'].replace(f'{old_ome}_',
                                                        f'{ome}_')
                entry['attributes'] = entry['attributes'].replace(f'{old_ome}_',
                                                        f'{ome}_')
            with open(f'{gff_dir}{ome}.gff3', 'w') as out:
                out.write(list2gff(gff))
            faa = fa2dict(row['faa'])
            new_faa = {k.replace(old_ome + '_', ome + '_'): v \
                       for k, v in faa.items()}
            with open(f'{faa_dir}{ome}.faa', 'w') as out:
                out.write(dict2fa(new_faa))
            row['fna'] = f'{fna_dir}{ome}.fna'
            row['faa'] = f'{faa_dir}{ome}.faa'
            row['gff3'] = f'{gff_dir}{ome}.gff3'

        # assembly accessions must be unique
        new_ome_mtdb = new_ome_mtdb.reset_index()

        return new_ome_mtdb
    else:
        return add_mtdb.reset_index()

 

def db2primary(addDB, refDB, save = False, combined = False):
    """Finalize an update by converting the updated MTDB into the primary
    MTDB"""
    if save:
        move_ns = shutil.copy
    else:
        move_ns = shutil.move

    addDB = addDB.reset_index()
    refDB = refDB.reset_index()

    refOmes = set(refDB['ome'])
    addOmes = set(addDB['ome'])
    base_ome2update_ome = {re.search(r'^[^\d]+\d+', x)[0]: x \
                           for x in refDB['ome'] if x}
    updates = {}
    refDB = refDB.set_index()
    if refOmes.intersection(addOmes) and not combined:
        eprint(refOmes.intersection(addOmes), flush = True)
        raise KeyError('ERROR: ome codes exist in database. Rerun predb2mtdb or remove manually')
    for i, ome in enumerate(addDB['ome']):
        base_ome = re.search(r'^[^\d]+\d+', ome)[0]
        if base_ome in base_ome2update_ome:
            update_ome = base_ome2update_ome[base_ome]
            updates[update_ome] = ome
            del refDB[update_ome]
        if os.path.isfile(addDB['gff3'][i]):
            move_ns(addDB['gff3'][i], format_path('$MYCOGFF3/' + ome + '.gff3'))
        elif not os.path.isfile(format_path('$MYCOGFF3/' + ome + '.gff3')):
            raise FileNotFoundError(f'{ome} missing gff3 for unknown reason')
        if os.path.isfile(addDB['fna'][i]):
            move_ns(addDB['fna'][i], format_path('$MYCOFNA/' + ome + '.fna'))
        elif not os.path.isfile(format_path('$MYCOFNA/' + ome + '.fna')):
            raise FileNotFoundError(f'{ome} missing fna for unknown reason')
        if os.path.isfile(addDB['faa'][i]):
            move_ns(addDB['faa'][i], format_path('$MYCOFAA/' + ome + '.faa'))
        elif not os.path.isfile(format_path('$MYCOFAA/' + ome + '.faa')):
            raise FileNotFoundError(f'{ome} missing faa for unknown reason')
        addDB['gff3'][i] = os.environ['MYCOGFF3'] + ome + '.gff3'
        addDB['fna'][i] = os.environ['MYCOFNA'] + ome + '.fna'
        addDB['faa'][i] = os.environ['MYCOFAA'] + ome + '.faa'
    addDB = addDB.set_index()
    for ome, row in addDB.items():
        refDB[ome] = row

    return refDB.reset_index(), updates


def control_flow(init, update, reference, add, taxonomy,
                 predb, save, nonpublished,
                 ncbi_only, lineage, rank, kingdom, failed,
                 forbidden, resume, no_md5, cpu, ncbi_email = False,
                 ncbi_api = None, overwrite = True, fallback = False):

    abbr2king = {'a': 'animals', 'r': 'archaea', 'f': 'fungi', 
                 'b': 'bacteria', 'p': 'plants'} 

    kingdom = kingdom.lower()
    if kingdom not in abbr2king:
        if kingdom not in set(abbr2king.values()):
            eprint('\nERROR: invalid --kingdom', flush = True)
            sys.exit(431)
    else:
        kingdom = abbr2king[kingdom]

    if not init \
        and not update \
        and not reference \
        and not add \
        and not taxonomy:
        eprint('\nERROR: --update/--init/--reference/--add must be specified', flush = True)
        sys.exit(15)
    elif reference and not init:
        eprint('\nERROR: --reference requires a --init directory', flush = True)
        sys.exit(14)
    elif lineage and not rank:
        eprint('\nERROR: --lineage requires --rank')
        sys.exit(16)
    elif lineage and not init:
        eprint('\nERROR: --lineage requires --init')
        sys.exit(17)
    elif predb and not init:
        eprint('\nERROR: --predb requires --init')
        sys.exit(18)
    elif predb and lineage:
        eprint('\nERROR: --predb and --lineage are incompatible')
        sys.exit(20)
    elif reference:
        if add:
            eprint('\nERROR: --add and --reference are incompatible')
            sys.exit(13)
        elif predb:
            eprint('\nERROR: --reference and --predb are incompatible')
            sys.exit(19)
        else:
            ref_db = mtdb(format_path(reference), add_paths = False)

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
    permitted_ranks = {'phylum', 'subphylum', 'class', 
                       'order', 'family', 'genus'}
    if lineage:
        lineage_constraints = split_input(lineage)
        rank_constraints = split_input(rank)
        if len(lineage_constraints) != len(rank_constraints):
            eprint('\nERROR: --lineage must be same length as --rank')
            sys.exit(18)
        for rank_c in rank_constraints:
            if rank_c.lower() not in permitted_ranks:
                eprint(f'\nERROR: accepted ranks: {permitted_ranks}')
                sys.exit(22)
        rank2lineages = defaultdict(set)
        for i, v in enumerate(lineage_constraints):
            rank2lineages[rank_constraints[i]].add(v.lower())
        rank2lineages = {k.lower(): sorted(v) for k, v \
                         in sorted(rank2lineages.items(), 
                                   key = lambda x: x[0])}

    # parse and check configuration nonpublished arguments
    config = {}
    if 'MYCODB' in os.environ:
        config_path = format_path('$MYCODB/../config/mtdb.json')
        if os.path.isfile(config_path):
            config = read_json(format_path(config_path))
            # for LEGACY installs:
            if 'lineage_constraints' not in config:
                config['lineage_constraints'] = {}
                write_json(config, config_path)
        elif not init:
            eprint('\nERROR: corrupted MycotoolsDB - no configuration found')
            sys.exit(21)
        if not init: # is MYCODB initialized?
    #            rogue_bool = config['rogue']
#                nonpublished = config['nonpublished']
                if bool(nonpublished) and not bool(config['nonpublished']):
                    config['nonpublished'] = validate_t_and_c(config, discrepancy = True)
                    write_json(config, config_path)
                if bool(config['jgi']) and bool(ncbi_only): #and not overwrite:
                    eprint('\nERROR: --ncbi_only specified after initialization',
                           flush = True)
                    sys.exit(173)
        elif init:
            if format_path(init) != \
               format_path(os.environ['MYCODB'] + '../../'):
                eprint('\nERROR: MTDB linked. Unlink via `mtdb -u`')
                sys.exit(175)

    # nonfungi is nonpublished by default because it is all GenBank
    if kingdom != 'fungi':
        nonpublished = True
    # archaic placeholder for reference / rogue DB setup
    elif nonpublished and rogue_bool: 
        nonpublished = validate_t_and_c(config)
    else:
        nonpublished = False

#    branch = 'stable'
    db_path = primaryDB()
    if not resume or add:
        date = datetime.now().strftime('%Y%m%d')
    else:
        date = str(resume)

    if not ncbi_email:
        ncbi_email, ncbi_api, jgi_email, jgi_pwd = loginCheck()
    Entrez.email = ncbi_email
    if ncbi_api:
        Entrez.api_key = ncbi_api


    if init:
        dbtype = kingdom
        init_dir = format_path(init)
        if os.path.isdir(init_dir):
            init_dir += 'mycotoolsdb/'
        if not init_dir.endswith( '/' ):
             init_dir += '/'
        envs = { 
            'MYCOFNA': init_dir + 'data/fna', 
            'MYCOFAA': init_dir + 'data/faa', 
            'MYCOGFF3': init_dir + 'data/gff3', 
            'MYCODB': init_dir + 'mtdb/'
            }
        os.environ['MYCODB'] = init_dir + 'mtdb/'
        output, config = initDB( 
            init_dir, dbtype, envs, dbtype, date = date, 
            rogue = rogue_bool, nonpublished = nonpublished,
            jgi = jgi, repo = format_path(reference),
            rank2lineages = rank2lineages
            )
        for env in envs:
            os.environ[env] = envs[env]
        orig_db = db2df(mtdb()) # initialize a new database
        update_path = output + 'log/' + date + '/'
        if not os.path.isdir(update_path):
            os.mkdir(update_path)
        mtdb_initialize(init_dir, init = True) #init_dir + 'config/mtdb.json', init = True)
    else:
        try:
            output = format_path('$MYCODB/..')
        except KeyError:
            eprint('\nERROR: MTDB not linked. Link via `mtdb -i <DB_PATH>`',
                    flush = True)
            sys.exit(50)
        update_path = output + 'log/' + date + '/'
        if not os.path.isdir(update_path):
            os.mkdir(update_path)
        if not True: #config['rogue']: # NEED TO MAKE THIS wget a particular URL
            old_db = db2df(db_path)
            shutil.move(db_path, update_path + os.path.basename(db_path))
            git_pull = subprocess.call([
                'git', 'pull', '-C', output + 'mtdb', config['repository'], '-B', branch
                 ], stdout = subprocess.PIPE, stderr = subprocess.PIPE 
                )
            new_db = db2df(primaryDB())
            orig_db = pd.concat([old_db, new_db.loc[~new_db['ome'].isin(old_db.index)]])
        else:
            orig_db = db2df(primaryDB())

    orig_db = orig_db.dropna(subset = ['ome'])

    if config['branch'] in {'prokaryote', 'bacteria'}:
        jgi = False
        group = 'prokaryotes'
        king = 'bacteria' # NEED to make DB tools pull from this
        rank = 'superkingdom'
    elif config['branch'] in {'plants'}:
        jgi = False
        group = 'eukaryotes'
        king = 'viridiplantae'
        rank = 'kingdom'
#    elif config['branch'] in {'protists'}:
 #       jgi = False
  #      group = 'eukaryotes'
   #     king = 'protists'
    #    rank = 'kingdom'
    elif config['branch'] in {'animals'}:
        jgi = False
        group = 'eukaryotes'
        king = 'metazoa'
        rank = 'kingdom'
    elif config['branch'] in {'archaea'}:
        jgi = False
        group = 'prokaryotes'
        king = 'Archaea'
        rank = 'superkingdom'
    else:
        jgi = not ncbi_only
        group = 'eukaryotes'
        king = 'fungi'
        rank = 'kingdom'


    if add or predb: # add predb2mtdb 2 master database
        if predb:
            add_predb = read_predb(predb_path)
            addDB, init_failed = predb2mtdb(add_predb, orig_db, update_path,
                                           cpus = cpu, remove = False, 
                                           spacer = '\t\t')
            if init_failed:
                if not failed:
                    eprint('\nERROR: some genomes failed curation', flush = True)
                    sys.exit(23)
                else:
                    eprint('\nWARNING: some genomes failed curation',
                           flush = True)
 
        else:
            addDB = mtdb(format_path(add))
        # we need full Paths for an addDB
        gff_fail, fna_fail, faa_fail = False, False, False
        if not all(os.path.isfile(format_path(x)) \
                   for x in addDB.reset_index()['gff3']):
            eprint('\nERROR: some GFF paths do not exist', flush = True)
            gff_fail = [x for x in addDB.reset_index()['gff3'] \
                        if not os.path.isfile(format_path(x))]
            print(','.join(gff_fail), flush = True)
        if not all(os.path.isfile(format_path(x)) \
                   for x in addDB.reset_index()['fna']):
            eprint('\nERROR: some FNA paths do not exist', flush = True)
            fna_fail = [x for x in addDB.reset_index()['fna'] \
                        if not os.path.isfile(format_path(x))]
            print(','.join(fna_fail), flush = True)
        if not all(os.path.isfile(format_path(x)) \
                   for x in addDB.reset_index()['faa']):
            eprint('\nERROR: some FAA paths do not exist', flush = True)
            faa_fail = [x for x in addDB.reset_index()['faa'] \
                        if not os.path.isfile(format_path(x))]
            print(','.join(faa_fail), flush = True)
        if gff_fail or fna_fail or faa_fail:
            sys.exit(124)

        addDB['aquisition_date'] = [date for x in addDB['ome']] 
        # make date the acquisition time
        orig_mtdb = mtdb(primaryDB())
        update_path = format_path('$MYCODB/../' + 'log/' + date + '/')
        if not os.path.isdir(update_path):
            os.mkdir(update_path)
        shutil.copy(primaryDB(), update_path)

        tax_path = f'{update_path}../taxonomy.tsv'
        tax_dicts = read_prev_tax(tax_path)
        tax_dicts = gather_taxonomy(addDB, api_key = ncbi_api, 
                                    king=king, rank = rank, 
                                    tax_dicts = tax_dicts, 
                                    output_path = tax_path)
        addDB, genus_dicts = assimilate_tax(addDB, tax_dicts) 
        addDB = check_add_mtdb(orig_mtdb, addDB, update_path, overwrite)

        write_forbid_omes(set(addDB['ome']), format_path('$MYCODB/../log/relics.txt'))

        new_mtdb, update_omes = db2primary(addDB, orig_mtdb, save = True)
        new_db_path = format_path('$MYCODB/' + date + '.mtdb')

        new_mtdb.df2db(new_db_path)

        if new_db_path != db_path:
            if db_path:
                os.remove(db_path)
        return new_db_path


    if taxonomy:
        new_db, update_mtdb = taxonomy_update(orig_db, update_path, date, 
                                              config, ncbi_email, ncbi_api,
                                              rank = rank, group = king)
        new_path = format_path('$MYCODB/' + date + '.mtdb')
        update_mtdb.df2db(new_path)
        sys.exit(0)
    elif reference:
        if any(not x for x in ref_db['published']) and not nonpublished:
            eprint('\nWARNING: nonpublished data detected in reference and will be ignored', 
                   flush = True)
        
        new_mtdb, update_mtdb = ref_update(
            ref_db, update_path, date, failed, jgi_email, jgi_pwd,
            config, ncbi_email, ncbi_api, cpus = cpu, check_MD5 = not bool(no_md5),
            jgi = jgi, group = group, kingdom = king,
            remove = not save, taxonomy = True, ncbi_fallback = fallback
            )
    else:
        new_mtdb, update_mtdb = rogue_update(
            orig_db, update_path, date, failed, jgi_email, jgi_pwd,
            config, ncbi_email, ncbi_api, cpus = cpu, 
            check_MD5 = not bool(no_md5), jgi = jgi, group = group,
            kingdom = king, remove = not save, 
            lineage_constraints = config['lineage_constraints'],
            ncbi_fallback = fallback
            )


    if not update_mtdb:
        eprint('\nNo new data acquired', flush = True)

    if not save: # add the predb2mtdb and remove files
#        df2db(db, format_path('$MYCODB/' + date + '.mtdb'))
        # output new database and new list of omes

        eprint('\nMoving data into database', flush = True)
        write_forbid_omes(set(new_mtdb['ome']), format_path('$MYCODB/../log/relics.txt'))

        new_path = format_path('$MYCODB/' + date + '.mtdb')
        if format_path(db_path) == new_path:
            shutil.copy(db_path, db_path + '.tmp')
        full_mtdb, update_omes = db2primary(update_mtdb, new_mtdb, save = False, 
                                            combined = True)
        full_mtdb.df2db(new_path + '.tmp')
        try:
            shutil.move(primaryDB(), update_path + os.path.basename(primaryDB()))
            # move master database to log if it exists
        except FileNotFoundError:
            pass
        shutil.move(new_path + '.tmp', new_path)
        rm_raw_data(update_path)
        eprint('\nMTDB update complete', flush = True)
#        gen_algn_db(
 #           update_path, set(full_mtdb['ome'])
  #          )
    else:
        # NEED to: insert note aboutrunning updatedb on predb
        new_mtdb.df2db(format_path(update_path + date + '.mtdb'))
        eprint(f'\nUpdate ready for `mtdb u -a` at ' \
            +  f'{format_path(update_path + date + ".mtdb")}')
        # output new database and new list of omes

    return primaryDB()


def main():

    abbr2king = {'a': '[a]nimals', #'r': 'a[r]chaea', 
                 'f': '[f]ungi', 
                 'b': '[b]acteria', 'p': '[p]lants'} 
    parser = argparse.ArgumentParser(description = "Initializes or updates " \
        + "MycotoolsDB (MTDB) derived from a kingdom of interest. Animals " \
        + "= metazoa, plants = viridiplantae")

    init_args = parser.add_argument_group('MTDB Initializiation')
    init_args.add_argument('-i', '--init', 
                           help = 'Initialize MTDB in dir')
    init_args.add_argument('-r', '--reference', 
        help = '[-i]: Initialize primary MTDB using a reference .mtdb')
    init_args.add_argument('-p', '--predb', 
        help = '[-i]: Initialize primary MTDB using a reference predb .tsv')
    init_args.add_argument('-k', '--kingdom', default = 'fungi',
        help = '[-i]: Kingdom - ' + str(sorted(abbr2king.values())) \
             + '; DEFAULT: fungi')


    upd_args = parser.add_argument_group("MTDB Updating")
    upd_args.add_argument('-u', '--update', action = 'store_true')
    upd_args.add_argument('-a', '--add', help = '.mtdb with full paths to add to database')
    upd_args.add_argument('-t', '--taxonomy', action = 'store_true',
        help = 'Remove old taxonomy metadata, update taxonomy, and exit')
    upd_args.add_argument('--keep', action = 'store_true',
        help = '[-a] Keep original MTDB data when adding overlapping accessions')
    upd_args.add_argument('--save', action = 'store_true', 
        help = '[-u] Do not integrate/delete new data; -a to complete')

#    init_args.add_argument('--reinit', action = 'store_true', help = 'Redownload all web data')
#    parser.add_argument('--rogue', action = 'store_true', 
 #       help = 'De novo MTDB') # currently required

    conf_args = parser.add_argument_group('Configuration')
    conf_args.add_argument('--nonpublished', action = 'store_true', 
        help = '[FUNGI]: Include MycoCosm restricted-use')
    conf_args.add_argument('--ncbi_only', help = '[FUNGI, -i]: Forego MycoCosm', 
        action = 'store_true')
    conf_args.add_argument('-l', '--lineage', 
                           help = '[-i, -r]: Lineage(s) to initialize MTDB with')
    conf_args.add_argument('-rk', '--rank',
                           help = '[-i, -l]: Rank(s) that positionally ' \
                                + 'correspond to -l')
    conf_args.add_argument('--failed', action = 'store_true', 
        help = 'Rerun/ignore failed')
    conf_args.add_argument('--forbidden', action = 'store_true', 
        help = 'Rerun forbidden')

#    conf_args.add_argument('--deviate', action = 'store_true', help = 'Deviate' \
 #       + ' from existing config without prompting')

    run_args = parser.add_argument_group('Runtime')
    run_args.add_argument('--resume', type = int, help = 'Resume previous date (YYYYmmdd)')
    run_args.add_argument('--no_md5', action = 'store_true', help = 'Skip NCBI MD5'
        + ' (expedite large reruns)')
    run_args.add_argument('--fallback', action = 'store_false', default = True, 
        help = '[ALPHA] use NCBI datasets utility for downloading NCBI data')
    run_args.add_argument('-c', '--cpu', type = int, default = 1)
    args = parser.parse_args()

    args_dict = { 
        'Primary MTDB': primaryDB(verbose = False), 'Update': args.update, 'Initialize': args.init,
        'Add': format_path(args.add), #'Rogue': rogue_bool, 
        'Include Restricted': bool(args.nonpublished), 'Resume': args.resume,
        'Retry failed': args.failed, 'Retry forbidden': args.forbidden,
        'Save raw data': args.save
        }

    findExecs(['datasets'], exit = {'datasets'})
    start_time = intro('Update MycotoolsDB', args_dict)

    control_flow(args.init, args.update, args.reference, args.add, args.taxonomy,
                 args.predb, args.save, args.nonpublished,
                 args.ncbi_only, args.lineage, args.rank, args.kingdom, args.failed,
                 args.forbidden, args.resume, args.no_md5, args.cpu, 
                 ncbi_email = None, overwrite = not args.keep, 
                 fallback = args.fallback)

    outro(start_time)


def cli():
    main()


if __name__ == '__main__':
    cli()
