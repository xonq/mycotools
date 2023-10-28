#! /usr/bin/env python3

#NEED to add taxonomy for --add options
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
import base64
import shutil
import getpass
import hashlib
import argparse
import subprocess
import numpy as np
import pandas as pd
from Bio import Entrez
from collections import defaultdict
from mycotools.lib.dbtools import db2df, df2db, gather_taxonomy, assimilate_tax, \
    primaryDB, loginCheck, log_editor, mtdb, mtdb_connect, \
    mtdb_initialize
from mycotools.lib.kontools import intro, outro, format_path, eprint, prep_output, collect_files, read_json, write_json
from mycotools.lib.biotools import fa2dict, gff2list
from mycotools.ncbiDwnld import esearch_ncbi, esummary_ncbi, main as ncbiDwnld
from mycotools.jgiDwnld import main as jgiDwnld
from mycotools.utils.ncbi2db import main as ncbi2db
from mycotools.utils.jgi2db import main as jgi2db
from mycotools.predb2mtdb import main as predb2mtdb
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
    rogue = False, nonpublished = False, jgi = False
    ):

    config = {
        'forbidden': forbidden,
        'repository': repo,
        'branch': branch,
        'nonpublished': nonpublished,
        'rogue': rogue,
        'jgi': jgi
    }

    return config


def add_vars(init_dir, dbtype):
    """Initialize the environmental variables for MTDB"""
    mtdb_initialize(init_dir, dbtype, init = True)


def initDB(init_dir, branch, envs, dbtype, date = None, 
           rogue = False, nonpublished = False, jgi = True,
           repo = None):
    '''Initialize database in `init_dir`'''

    new_dirs = [
        init_dir + 'data/', init_dir + 'config/', 
        init_dir + 'data/fna/', init_dir + 'data/faa/',
        init_dir + 'mtdb',
        init_dir + 'data/gff3/', 
        init_dir + 'log/',
        init_dir + 'data/db/'
    ]
    output = prep_output( init_dir, cd = True )
    if not output.endswith('/'):
        output += '/'
    for new_dir in new_dirs:
        if not os.path.isdir( new_dir ):
            os.mkdir( new_dir )

    config = gen_config(branch = branch, rogue = rogue, 
                        forbidden = '$MYCODB/log/forbidden.tsv',
                        nonpublished = nonpublished, jgi = jgi,
                        repo = repo)
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

def parse_ncbi2jgi(file_path):
    """Parse previously collected NCBI to JGI data to limit querying"""
    ncbi2jgi = {}
    if not os.path.isfile(file_path):
        with open(file_path, 'w') as out:
            out.write('#ncbi_acc\tmycocosm_portal')
    else:
        with open(file_path, 'r') as raw:
            for line in raw:
                if not line.startswith('#'):
                    d = line.rstrip().split('\t')
                    ncbi, jgi = d[0], d[1].lower()
                    ncbi2jgi[ncbi] = jgi

    return ncbi2jgi

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

def add_ncbi2jgi(jgi2ncbi, file_path = None):
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
    mycocosm_url = 'https://mycocosm.jgi.doe.gov/ext-api/mycocosm/catalog/' + \
        'download-group?flt=&seq=all&pub=all&grp=fungi&srt=released&ord=desc'
    ):
    """Download the MycoCosm genome data spreadsheet, format to UTF-8 and
    return a Pandas dataframe of the data"""

    if not os.path.isfile(out_file):
        subprocess.call(
            ['curl', mycocosm_url, '-o', out_file + '.tmp'],
            stdout = subprocess.PIPE
            )
        shutil.move(out_file + '.tmp', out_file)
    try:
        jgi_df = pd.read_csv( out_file, encoding = 'cp1252' )
    except UnicodeDecodeError:
        jgi_df = pd.read_csv( out_file, encoding = 'utf-8' )
    jgi_df.columns = [x.replace('"','').replace('"','') for x in jgi_df.columns]
    return jgi_df


def dwnld_ncbi_table( 
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

def prep_ncbi_names(ncbi_df):

    ncbi_df['strain'] = ''
    for i, row in ncbi_df.iterrows():
        organism = row['#Organism/Name'].split(' ')
        organism = [ 
            x.replace('[','').replace(']','').replace('(','').replace(')','') for x in organism 
            ]
        ncbi_df.at[i, 'genus'] = organism[0]
        if len(organism) > 1:
            ncbi_df.at[i, 'species'] = organism[1]
            if len(organism) > 2:
                ncbi_df.at[i, 'strain'] = ''.join(organism[2:])
                ncbi_df.at[i, 'strain'] = re.sub(r'[^a-zA-Z0-9]', '', ncbi_df['strain'][i])
        else:
            ncbi_df.at[i, 'species'] = 'sp.'

    return ncbi_df



def clean_ncbi_df(ncbi_df, kingdom = 'Fungi'):
    ncbi_df = ncbi_df.astype(str).replace(np.nan, '')
    ncbi_df['strain'] = ''

    if kingdom.lower() == 'fungi':
        # extract group of interest (case sensitive to first letter)
        ncbi_df = ncbi_df[ncbi_df['Group'] == 'Fungi']

    # remove entries without sufficient metadata
    ncbi_df = ncbi_df.dropna(subset = ['genus'])

    # remove assembly- & annotation-lacking entries
    ncbi_df = ncbi_df[~ncbi_df['Genes'].isin({'-', '', '0'})]
    ncbi_df = ncbi_df[~ncbi_df['Proteins'].isin({'-', '' ,'0'})]
    ncbi_df = ncbi_df[~ncbi_df['Assembly Accession'].isin({'-', '', '0'})]

    # sort by version, keep most recent version of duplicate assemblies
    ncbi_df['version'] = pd.to_datetime(ncbi_df['Modify Date'])
    ncbi_df = ncbi_df.sort_values(by = 'version', ascending = False)
    ncbi_df = ncbi_df.drop_duplicates('Assembly Accession')

    # check for explicit new versions of assemblies
    ncbi_df['assembly_base'] = [
        x[:x.find('.')] for x in list(ncbi_df['Assembly Accession'])
        ]
    ncbi_df['assembly_vers'] = [
        x[x.find('.'):] for x in list(ncbi_df['Assembly Accession'])
        ]
    ncbi_df = ncbi_df.sort_values(by = 'assembly_vers', ascending = False)
    ncbi_df = ncbi_df.drop_duplicates('assembly_base')

    return ncbi_df

def rm_ncbi_overlap(ncbi_df, mycocosm_df, ncbi2jgi, fails = set(), api = 3):
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
    jgi2ncbi, jgi2biosample = {v: k for k, v in ncbi2jgi.items()}, {}
    ass_count = 0

    for i, row in ncbi_df.iterrows():
        if row['assembly_acc'] in ncbi2jgi:
            jgi2biosample[ncbi2jgi[row['assembly_acc']]] = \
                row['BioSample Accession']
            todel.append(i)
#        elif row['assembly_acc'] in fails:
 #           pass
#            todel.append(i)
        elif row['assembly_acc'] not in fails:
            ass_uid = esearch_ncbi(row['assembly_acc'],
                                   'assembly', 'assembly')
            ncbi_df.at[i, 'uid'] = ass_uid
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
            else:
                fails.add(row['assembly_acc'])

    ncbi_jgi_overlap = pd.DataFrame(columns = ncbi_df.columns)
    ncbi_jgi_overlap = ncbi_df.loc[todel]
    ncbi_df = ncbi_df.drop(todel)
#    for i in reversed(todel):
 #       ncbi_jgi_overlap = pd.concat([ncbi_jgi_overlap, ncbi_df.loc[i]])
#        ncbi_df = ncbi_df.drop(i)
    
    return ncbi_df, jgi2ncbi, jgi2biosample, fails, ncbi_jgi_overlap

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


def ref_update(
    ref_db, update_path, date, rerun, jgi_email, jgi_pwd,
    config, ncbi_email, ncbi_api, cpus = 1, check_MD5 = True,
    jgi = True, group = 'eukaryotes', kingdom = 'Fungi',
    remove = True
    ):
    """Initialize/Update the primary MTDB based on a reference database
    acquired external from any existing primary MTDB"""
# NEED to mark none for new databases' refdb
    # initialize update
    print('\nInitializing run', flush = True)
    mk_wrk_dirs(update_path)
#    prev_failed = parse_failed(rerun = rerun, file_path = format_path('$MYCODB/../log/failed.tsv'))
#    duplicates = parse_dups(format_path('$MYCODB/../log/duplicates.tsv'))

#    forbid_omes = acq_forbid_omes(
 #       file_path = format_path('$MYCODB/../log/relics.txt')
  #      )

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

    
#        new_db_path = update_path + date + '.checkpoint.jgi.mtdb'
 #       if not os.path.isfile(new_db_path): 
  #          if len(jgi_db) > 0:
   #             df2db( jgi_db, jgi_db_path )
#                if not db is None:
 #                   new_db = pd.concat([jgi_db, db])
  #              else:
   #                 new_db = jgi_db
#            else:
 #               new_db = db
  #          df2db(new_db, new_db_path)
   #     else:
    #        new_db = db2df(new_db_path)
    else:
        jgi_mtdb = mtdb()
#        new_dups = duplicates
    new_db = jgi_mtdb.mtdb2pd()

    print('\nAssimilating NCBI', flush = True)
    if not os.path.isfile(update_path + date + '.ncbi.predb'):
        print('\tDownloading NCBI data', flush = True)
        ncbi_predb, ncbi_failed1 = ncbiDwnld(
            assembly = True, proteome = False, gff3 = True,
            ncbi_df = ncbi_df, remove = True, output_path = update_path,
            column = 'assembly_acc', ncbi_column = 'genome', check_MD5 = check_MD5
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

    try:
        ncbi_db = db2df(update_path + date + '.ncbi.predb2.mtdb')
    except pd.errors.EmptyDataError:
        ncbi_db = pd.DataFrame({x: [] for x in refdbncbi.keys()})
    if len(ncbi_db) > 0:
#        df2db(ncbi_db, ncbi_db_path)
        new_db = pd.concat([new_db, ncbi_db])

    print('\nAssimilating NCBI taxonomy data', flush = True)
    if kingdom.lower() == 'fungi':
        rank = 'kingdom'
    else:
        rank = 'superkingdom'

    tax_dicts = gather_taxonomy(new_db, api_key = ncbi_api, 
                                king=kingdom, rank = rank)
    new_db, genus_dicts = assimilate_tax(new_db, tax_dicts) 
    dupFiles = { 'fna': {}, 'faa': {}, 'gff3': {} }

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

    for ome, row in update_mtdb.items():
        if row['genus'] in genus_dicts:
            row['taxonomy'] = genus_dicts[row['genus']]

    return new_db, update_mtdb 


def rogue_update(
    db, update_path, date, rerun, jgi_email, jgi_pwd,
    config, ncbi_email, ncbi_api, cpus = 1, check_MD5 = True,
    jgi = True, group = 'eukaryotes', kingdom = 'Fungi',
    remove = True
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
    pre_ncbi_df0 = dwnld_ncbi_table(update_path + date + '.ncbi.tsv',
                                   group = group) 
    pre_ncbi_df1 = prep_ncbi_names(pre_ncbi_df0)
    ncbi_df = clean_ncbi_df(pre_ncbi_df1, kingdom = kingdom)
    ncbi_df = ncbi_df.rename(columns={'Assembly Accession': 'assembly_acc'})

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
        jgi_df = dwnld_mycocosm(mycocosm_path)
        jgi_df['biosample'] = ''
        # acquire the mycocosm master table

        print('\tSearching NCBI for MycoCosm overlap', flush = True)
        ncbi2jgi = parse_ncbi2jgi(update_path + '../ncbi2jgi.tsv')
        true_ncbi = parse_true_ncbi(update_path + '../true_ncbi.tsv')
        ncbi_df, jgi2ncbi, jgi2biosample, true_ncbi, ncbi_jgi_overlap = \
            rm_ncbi_overlap(ncbi_df, jgi_df, ncbi2jgi, true_ncbi, api = api)
        print('\t\t' + str(len(jgi2ncbi)) + ' overlapping genomes',
             flush = True)
        add_true_ncbi(true_ncbi, update_path + '../true_ncbi.tsv')
        add_ncbi2jgi(jgi2ncbi, update_path + '../ncbi2jgi.tsv')
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
            jgi_mtdb, jgi_failed1 = predb2mtdb(jgi_premtdb, refdbjgi, update_path,
                                            forbidden = forbid_omes, cpus = cpus, 
                                            remove = remove, spacer = '\t\t')
            jgi_failed.extend(jgi_failed1)
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
                df2db( jgi_db, jgi_db_path )
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

    print('\nAssimilating NCBI (10 download/minute w/API key, 3 w/o)', flush = True)
    new_db['version'] = new_db['version'].astype(str)
    if not os.path.isfile(update_path + date + '.ncbi.predb'):
#    if not os.path.isfile(update_path + date + '.ncbi.predb'):
        print('\tDownloading NCBI data', flush = True)
        ncbi_predb, new_db, ncbi_failed1 = ncbi2db( 
            update_path, ncbi_df, ref_db = new_db, 
            date = date, failed_dict = prev_failed, 
            rerun = rerun, duplicates = duplicates,
            check_MD5 = check_MD5
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

    print('\nAssimilating NCBI taxonomy data', flush = True)
    if kingdom.lower() == 'fungi':
        rank = 'kingdom'
    else:
        rank = 'superkingdom'
    tax_dicts = gather_taxonomy(new_db, api_key = ncbi_api, 
                                king=kingdom, rank = rank)
    new_db, genus_dicts = assimilate_tax(new_db, tax_dicts) 
    dupFiles = { 'fna': {}, 'faa': {}, 'gff3': {} }

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

    return new_db, update_mtdb 


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
                           for x in refDB['ome']}
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
    

def main():

    parser = argparse.ArgumentParser(description = "Initializes or updates " \
        + "MycotoolsDB (MTDB)")

    init_args = parser.add_argument_group('Initializiation')
    init_args.add_argument('-p', '--prokaryote', action = 'store_true',
        help = '[PROKARY]: experimental prokaryote MTDB')
    init_args.add_argument('-u', '--update', action = 'store_true')
    init_args.add_argument('-i', '--init', 
                           help = 'Initialize MTDB in dir; failed inits ' \
                                + 'must be unlinked via mtdb -u prior to resume')
    init_args.add_argument('-a', '--add', help = 'Curated .mtdb to add to database')
    init_args.add_argument('-r', '--reference', 
        help = '[-i] Initialize primary MTDB using a reference .mtdb')
    init_args.add_argument('--resume', type = int, help = 'Resume previous date (YYYYmmdd)')
#    init_args.add_argument('--reinit', action = 'store_true', help = 'Redownload all web data')
#    parser.add_argument('--rogue', action = 'store_true', 
 #       help = 'De novo MTDB') # currently required

    conf_args = parser.add_argument_group('Configuration')
    conf_args.add_argument('--nonpublished', action = 'store_true', 
        help = '[FUNGI]: Include MycoCosm restricted-use')
    conf_args.add_argument('--ncbi_only', help = '[FUNGI]: [-i] Forego MycoCosm', 
        action = 'store_true')
    conf_args.add_argument('--deviate', action = 'store_true', help = 'Deviate' \
        + ' from existing config without prompting')

    run_args = parser.add_argument_group('Runtime')
    run_args.add_argument('--save', action = 'store_true', 
        help = '[-u] Do not integrate/delete new data. -a to complete')
    run_args.add_argument('--failed', action = 'store_true', help = 'Rerun failed')
    run_args.add_argument('--forbidden', action = 'store_true', help = 'Rerun forbidden')
    run_args.add_argument('--no_md5', action = 'store_true', help = 'Skip NCBI MD5'
        + ' (expedite large reruns)')
#    run_args.add_argument('--overwrite', action = 'store_true',
#        help = 'Remove entries that violate new MTDB parameters')
    run_args.add_argument('-c', '--cpu', type = int, default = 1)
    args = parser.parse_args()

    if not args.init \
        and not args.update \
        and not args.reference \
        and not args.add:
        eprint('\nERROR: --update/--init/--reference/--add must be specified', flush = True)
        sys.exit(15)
    elif args.reference and not args.init:
        eprint('\nERROR: --reference requires a --init directory', flush = True)
        sys.exit(14)
    elif args.reference:
        if args.add:
            eprint('\nERROR: --add and --reference are incompatible')
            sys.exit(13)
        else:
            ref_db = mtdb(format_path(args.reference), add_paths = False)

#    if args.rogue:
    rogue_bool = True
    if args.ncbi_only:
        jgi = False
    else:
        jgi = True

    
    config = {'branch': None, 'nonpublished': None}
    if 'MYCODB' in os.environ:
        if not args.init: # is MYCODB initialized?
            config_path = format_path('$MYCODB/../config/mtdb.json')
            if os.path.isfile(config_path):
                config = read_json(format_path(config_path))
    #            rogue_bool = config['rogue']
#                args.nonpublished = config['nonpublished']
                if bool(args.nonpublished) and not bool(config['nonpublished']):
                    config['nonpublished'] = validate_t_and_c(config, discrepancy = True)
                    write_json(config, format_path(config_path))
                if bool(config['jgi']) and bool(args.ncbi_only): #and not args.overwrite:
                    eprint('\nERROR: --ncbi_only specified after initialization',
                           flush = True)
                    sys.exit(173)
        elif args.init:
            if format_path(args.init) != format_path(os.environ['MYCODB']):
                eprint('\nERROR: MTDB initialized elsewhere and linked. Unlink via `mtdb -u`')
                sys.exit(175)
    
    if args.prokaryote or config['branch'] == 'prokaryote':
        nonpublished = True
    elif args.nonpublished and rogue_bool: 
        nonpublished = validate_t_and_c(config)
    elif args.nonpublished and not rogue_bool:
        eprint('\nERROR: Specify --rogue to assimilate use-restricted data.' , flush = True)
        sys.exit(5)
    else:
        nonpublished = False

#    branch = 'stable'

    db_path = primaryDB()
    args_dict = { 
        'Database': db_path, 'Update': args.update, 'Initialize': args.init,
        'Add': format_path(args.add), #'Rogue': rogue_bool, 
        'Include Restricted': bool(nonpublished), 'Resume': args.resume,
        'Retry failed': args.failed, 'Retry forbidden': args.forbidden,
        'Save raw data': args.save
        }

    start_time = intro('Update MycotoolsDB', args_dict)
    if not args.resume or args.add:
        date = start_time.strftime('%Y%m%d')
    else:
        date = str(args.resume)

    ncbi_email, ncbi_api, jgi_email, jgi_pwd = loginCheck()
    Entrez.email = ncbi_email
    if ncbi_api:
        Entrez.api_key = ncbi_api

    if args.add: # add predb2mtdb 2 master database
        addDB = mtdb(format_path(args.add))
        addDB['aquisition_date'] = [date for x in addDB['ome']] 
        # make date the acquisition time
        orig_mtdb = mtdb(primaryDB())
        update_path = format_path('$MYCODB/../' + 'log/' + date + '/')
        if not os.path.isdir( update_path ):
            os.mkdir( update_path )
        shutil.copy(primaryDB(), update_path)

        if config['branch'] == 'prokaryote':
            jgi = False
            group = 'prokaryotes'
            king = 'bacteria' # NEED to make DB tools pull from this
            rank = 'superkingdom'
        else:
            jgi = not args.ncbi_only
            group = 'eukaryotes'
            king = 'fungi'
            rank = 'kingdom'
     
        tax_dicts = gather_taxonomy(addDB, api_key = ncbi_api, 
                                    king=king, rank = rank)
        addDB, genus_dicts = assimilate_tax(addDB, tax_dicts) 


        new_mtdb, update_omes = db2primary(addDB, orig_mtdb, save = args.save)
        new_db_path = format_path('$MYCODB/' + date + '.mtdb')


        new_mtdb.df2db(new_db_path)


#        if update_omes and args.clear_cache:
 #           for update_ome in update_omes:
  #              ome_gff3 = os.environ['MYCOGFF3'] + update_ome + '.gff3'
   #             ome_fna = os.environ['MYCOFNA'] + update_ome + '.fna'
    #            ome_faa = os.environ['MYCOFAA'] + update_ome + '.faa'
     #           for file_ in [ome_gff3, ome_fna, ome_faa]:
      #              if os.path.isfile(file_):
       #                 os.remove(file_)
        if new_db_path != db_path:
            os.remove(db_path)
        outro(start_time)

    if args.init:
        if args.prokaryote:
            dbtype = 'prokaryote'
        else:
            dbtype = 'fungi'
        init_dir = format_path(args.init)
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
            jgi = jgi, repo = format_path(args.reference)
            )
        for env in envs:
            os.environ[ env ] = envs[ env ]
        orig_db = db2df(mtdb()) # initialize a new database
        update_path = output + 'log/' + date + '/'
        if not os.path.isdir( update_path ):
            os.mkdir( update_path )
        mtdb_initialize(init_dir, init = True) #init_dir + 'config/mtdb.json', init = True)
    else:
        output = format_path('$MYCODB/..')
        update_path = output + 'log/' + date + '/'
        if not os.path.isdir( update_path ):
            os.mkdir( update_path )
        if not config['rogue']:            
            old_db = db2df( db_path )
            shutil.move( db_path, update_path + os.path.basename(db_path) )
            git_pull = subprocess.call( [
                'git', 'pull', '-C', output + 'mtdb', config['repository'], '-B', branch
                    ], stdout = subprocess.PIPE, stderr = subprocess.PIPE 
                )
            new_db = db2df(primaryDB())
            orig_db = pd.concat( [old_db, new_db.loc[~new_db['ome'].isin(old_db.index)]] )
        else:
            orig_db = db2df(primaryDB())

    orig_db = orig_db.dropna(subset = ['ome'])

    if config['branch'] == 'prokaryote':
        jgi = False
        group = 'prokaryotes'
        king = 'bacteria' # NEED to make DB tools pull from this
    else:
        jgi = not args.ncbi_only
        group = 'eukaryotes'
        king = 'fungi'
        
    # THIS IS WHERE WE INTEGRATE GIT LINKING
#    if not config['rogue']:
 #       missing_db = getMissingEntries(orig_db)
  #      new_db = stdUpdate(
   #         missing_db, ncbi_email, ncbi_api, update_path, 
    #        orig_db, jgi_email, jgi_pwd, date
     #       ) # NEED massive overhaul and git generation
#    else:

    if not args.reference:
        new_db, update_mtdb = rogue_update(
            orig_db, update_path, date, args.failed, jgi_email, jgi_pwd,
            config, ncbi_email, ncbi_api, cpus = args.cpu, 
            check_MD5 = not bool(args.no_md5), jgi = jgi, group = group,
            kingdom = king, remove = not args.save
            )
    else:
        if any(not x for x in ref_db['published']) and not args.nonpublished:
            eprint('\nWARNING: nonpublished data detected in reference and will be ignored', 
                   flush = True)
        new_db, update_mtdb = ref_update(
            ref_db, update_path, date, args.failed, jgi_email, jgi_pwd,
            config, ncbi_email, ncbi_api, cpus = args.cpu, check_MD5 = not bool(args.no_md5),
            jgi = jgi, group = group, kingdom = king,
            remove = not args.save
            )
    if not update_mtdb:
        eprint('\nNo new data acquired', flush = True)

    if not args.save: # add the predb2mtdb and remove files
#        df2db(db, format_path('$MYCODB/' + date + '.mtdb'))
        write_forbid_omes(set(new_db['ome']), format_path('$MYCODB/../log/relics.txt'))
        # output new database and new list of omes

        eprint('\nMoving data into database', flush = True)
        new_mtdb = mtdb.pd2mtdb(new_db)
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
        eprint('\nGathering assembly statistics', flush = True)
        assStats(primaryDB(), format_path('$MYCODB/../data/assemblyStats.tsv'),
                 args.cpu)
    
        eprint('\nGathering annotation statistics', flush = True)
        annStats(primaryDB(), format_path('$MYCODB/../data/annotationStats.tsv'),
                 args.cpu)
#        gen_algn_db(
 #           update_path, set(full_mtdb['ome'])
  #          )
    else:
        # NEED to: insert note aboutrunning updatedb on predb
        df2db(db, format_path(update_path + date + '.mtdb'))
        write_forbid_omes(set(db['ome']), format_path('$MYCODB/../log/relics.txt'))
        # output new database and new list of omes


    outro(start_time)

def cli():
    main()


if __name__ == '__main__':
    cli()
