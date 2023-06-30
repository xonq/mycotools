#! /usr/bin/env python3

# NEED a db check to ensure the log is relevant to the input

import os
import re
import sys
import gzip
import time
import urllib
import requests
import argparse
import subprocess
import numpy as np
import pandas as pd
from Bio import Entrez
from datetime import datetime
from mycotools.lib.kontools import intro, outro, format_path, prep_output, eprint, vprint, findExecs
from mycotools.lib.dbtools import log_editor, loginCheck, mtdb, read_tax

def db2df(data, stdin = False):
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

    if len(db_df.keys()) == 16: # legacy conversion TO BE DEPRECATED
        eprint('\tWARNING: Legacy MycotoolsDB format will be removed in the future.', flush = True)
        db_df.columns = [
            'ome', 'genus', 'species', 'strain', 'version',
            'biosample', 'fna', 'faa', 'gff3', 'taxonomy', 'ecology',
            'eco_conf',
            'source', 'published', 'assembly_acc', 'acquisition_date'
            ]
    else:
        db_df.columns = columns
    for i, row in db_df.iterrows():
        db_df.at[i, 'taxonomy'] = read_tax(
            row['taxonomy']
            )
        db_df.at[i, 'taxonomy']['genus'] = row['genus']
        db_df.at[i, 'taxonomy']['species'] = \
            row['genus'] + ' ' + row['species']
        # if malformatted due to decreased entries in some lines, this will raise an IndexError
        if row['fna'] or row['fna'] == row['ome'] + '.fna': # abbreviated line w/o file coordinates
           db_df.at[i, 'fna'] = os.environ['MYCOFNA'] + row['ome'] + '.fna'
           db_df.at[i, 'faa'] = os.environ['MYCOFAA'] + row['ome'] + '.faa'
           db_df.at[i, 'gff3'] = os.environ['MYCOGFF3'] + row['ome'] + '.gff3'
        else: # has file coordinates
           db_df.at[i, 'fna'] = format_path(row['fna'])
           db_df.at[i, 'faa'] = format_path(row['faa'])
           db_df.at[i, 'gff3'] = format_path(row['gff3'])

    if len(db_df) == 16: # LEGACY conversion to be deprecated
        del db_df['ecology']
        del db_df['eco_conf']

    return db_df


def prepare_folders( output_path, gff, prot, assem, transcript ):

    file_types = []
    if assem:
        if not os.path.exists(output_path + 'assembly'):
            os.mkdir(output_path + 'assembly')
        file_types.append( 'assembly' )
    if gff:
        if not os.path.exists(output_path + 'gff3'):
            os.mkdir(output_path + 'gff3')
        file_types.append( 'gff3' )
    if prot:
        if not os.path.exists(output_path + 'proteome'):
            os.mkdir(output_path + 'proteome')
        file_types.append( 'proteome' ) 
    if transcript:
        if not os.path.exists(output_path + 'transcript'):
            os.mkdir(output_path + 'transcript')
        file_types.append( 'transcript' )

    return file_types  


def compile_log( output_path, remove = False ):

    ass_prots = {}
    if not os.path.isfile( output_path + 'ncbiDwnld.log' ):
        with open( output_path + 'ncbiDwnld.log', 'w' ) as out:
            out.write('#ome\tassembly_acc\tassembly\tproteome\tgff3\ttranscript\t' + \
            'assemMD5\tprotMD5\tgff3MD5\ttransMD5\tgenome_id(s)')

        # too risky, too many things can go wrong and then users would be in a
        # loop, but necessary for huge downloads
    else:
        with open( output_path + 'ncbiDwnld.log', 'r' ) as raw:
            for line in raw:
                if not line.startswith('#'):
                    data = [x.rstrip() for x in line.split('\t')]
                    while len(data) < 10:
                        data.append('')
                    ass_prots[data[0]] = { 
                        'assembly_acc': str(data[1]), 'assembly': str(data[2]), 
                        'proteome': str(data[3]), 'gff3': str(data[4]), 
                        'transcript': str(data[5]), 'assembly_md5': str(data[6]),
                        'proteome_md5': str(data[7]), 'gff3_md5': str(data[8]),
                        'transcript_md5': str(data[9]), 'genome_id': data[10].split(',')
                        }

    return ass_prots

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
    search_term, esc_count = accession + '[' + column + ']', 0
    while esc_count < 10:
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
def collect_ftps(
    ncbi_df, ass_prots, api_key=0, column = 'assembly_acc',
    ncbi_column='Assembly Accession', database="assembly", output_path = '',
    verbose=True, remove = False, spacer = '\t\t'
    ):

    count, failed = 0, []

# for each row in the assembly, grab the accession number, form the search term for Entrez, use Entrez,
    if ncbi_column in {'assembly', 'genome', 'uid'}:
        out_df = ncbi_df[ncbi_df.index.isin(set(ass_prots.keys()))]
        ncbi_df = ncbi_df[~ncbi_df.index.isin(set(ass_prots.keys()))]
    for accession, row in ncbi_df.iterrows():
        if accession in ass_prots: # add all rows that have indices associated with this
            # query type
            out_df = pd.concat([out_df, row.to_frame().T])
            icount = 1
            test = str(accession) + '_' + str(icount)
            while test in ass_prots:
                count += 1
    #                row['assembly_acc'] = ass_prots[test]['genome_id']
                if 'ome' in row.keys():
                    row['ome'] = None # haven't assigned a mycotools ID yet
                out_df = pd.concat([out_df, row.to_frame().T])
                sys.exit()
                test = str(accession) + '_' + str(icount)
            continue

        elif pd.isnull(row[column]) or not row[column]: # ignore blank entries
            ass_prots[str(accession)] = {
                'accession': accession, 'assembly': '',
                'proteome': '', 'gff3': '', 'transcript': '',
                'assembly_md5': '', 'proteome_md5': '',
                'gff3_md5': '', 'transcript_md5': '', 'genome_id': ''
                }
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

        if ncbi_column in {'assembly', 'genome', 'uid'}: # be confident it is the most
        # recent assembly UID
            genome_id = [str(max([int(i) for i in genome_id]))]

        icount = 0
        for ID in genome_id:
            if count:
                new_acc = str(accession) + '$' + str(icount)
            else:
                new_acc = accession
# obtain the path from a summary of the ftp directory and create the standard paths for proteomes and assemblies
            ass_prots[str(new_acc)] = {
                'accession': accession, 'assembly': '', 'proteome': '', 
                'gff3': '', 'transcript': '', 'assembly_md5': '',
                'proteome_md5': '', 'gff3_md5': '', 
                'transcript_md5': '', 'genome_id': ID
                }

            record = esummary_ncbi(ID, database)
            assemblyID = record['DocumentSummarySet']['DocumentSummary'][0]['AssemblyAccession']
            ftp_path = str(record['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_GenBank'])

            if not ftp_path:
                eprint(spacer + '\t' + accession + ' failed to return any FTP path', flush = True)
                failed.append([accession, datetime.strftime(row['version'], '%Y%m%d')])
                continue

            esc_count = 0
            ass_md5, gff_md5, trans_md5, prot_md5, md5s = '', '', '', '', {}
            basename = os.path.basename(ftp_path)
            while True:
                try:
                    esc_count += 1
                    request = requests.get( 
                        (ftp_path + '/md5checksums.txt').replace('ftp://','https://'), timeout = 120
                        )
                    break
                except (requests.exceptions.ConnectionError, requests.exceptions.ReadTimeout) as e:
                    if esc_count == 20:
                        eprint(spacer + '\tERROR: Failed to fulfill ftp request!', flush = True)
                        sys.exit(70)
                    count = 0
                    time.sleep(1)
            count += 1
            if request.status_code != 200:
                eprint(spacer + '\t' + ome + '\tno ' + md5, flush = True)
                dwnld = 69
            else:
                dwnld = subprocess.call(
                    ['curl', ftp_path + '/md5checksums.txt', '-o', output_path + '.tmpmd5'], 
                    stdout = subprocess.PIPE, stderr = subprocess.PIPE
                    )
                count += 1
                with open( output_path + '.tmpmd5', 'r' ) as raw:
                    for line in raw:
                        data = line.rstrip().split('  ')
                        if data:
                            md5s[ftp_path + '/' + os.path.basename(data[1])] = data[0]
            tranname = os.path.basename(ftp_path.replace('/GCA','/GCF'))
            assembly = ftp_path + '/' + basename + '_genomic.fna.gz'
            if assembly in md5s:
                ass_md5 = md5s[assembly]
            else:
                assembly = ''
            proteome = ftp_path + '/' + basename + '_protein.faa.gz'
            if proteome in md5s:
                prot_md5 = md5s[proteome]
            else:
                proteome = ''
            gff3 = ftp_path + '/' + basename + '_genomic.gff3.gz'
            test_gff3 = re.sub( r'\.gff3\.gz$', '.gff.gz', gff3 )
            if gff3 in md5s:
                gff_md5 = md5s[gff3]
            elif test_gff3 in md5s:
                gff3 = test_gff3
                gff_md5 = md5s[gff3]
            else:
                gff3 = ''
    
            transcript = ftp_path.replace('/GCA', '/GCF') + '/' + tranname + '_rna.fna.gz'
            if transcript in md5s:
                trans_md5 = md5s[transcript]
    
            if (not assembly or not gff3) and remove:
                failed.append([accession, datetime.strftime(row['version'], '%Y%m%d')])
    
            log_editor( 
                output_path + 'ncbiDwnld.log', str(new_acc), 
                str(accession) + '\t' + accession + '\t' +  assembly + '\t' + \
                proteome + '\t' + gff3 + '\t' + transcript + '\t' + \
                ass_md5 + '\t' + prot_md5 + '\t' + gff_md5 + '\t' + trans_md5 + \
                '\t' + ID
                )
            ass_prots[str(new_acc)] = {
                'accession': accession, 'assembly': assembly, 'assembly_md5': ass_md5,
                'proteome': proteome, 'proteome_md5': prot_md5,
                'gff3': gff3, 'gff3_md5': gff_md5,
                'transcript': transcript, 'transcript_md5': trans_md5,
                'genome_id': ID
                }
            row['dwnld_id'] = ID
            if icount:
                if 'ome' in row:
                    row['ome'] = None
            out_df = pd.concat([out_df, row.to_frame().T])
            icount += 1
    
    # if no API key is used, we can only generate 3 queries per second, otherwise we can use 10
            count += 1
            count = wait_for_ncbi(count, api_key)
    
    return ass_prots, failed, out_df

# download the file depending on the type inputted
def download_files(acc_prots, acc, file_types, output_dir, count,
                   remove = False, spacer = '\t\t'):

    esc_count, dwnlds = 0, {}
    for file_type in file_types:
        ftp_link = acc_prots[file_type]
        dwnlds[file_type] = -1
        if file_type == 'assembly':
            file_path = output_dir + 'assembly/' + \
                os.path.basename(acc_prots[file_type])
        elif file_type == 'gff3':
            file_path = output_dir + 'gff3/' + \
                os.path.basename(acc_prots[file_type])
        elif file_type == 'proteome':
            file_path = output_dir + 'proteome/' + \
                os.path.basename(acc_prots[file_type])
        elif file_type == 'transcript':
            file_path = output_dir + 'transcript/' + \
                os.path.basename(acc_prots[file_type])

        if os.path.isfile( file_path ):
            count += 1
            md5_cmd = subprocess.run([
                'md5sum', file_path],
                stdout = subprocess.PIPE) # check the md5
            md5_res = md5_cmd.stdout.decode('utf-8')
            md5_find = re.search(r'\w+', md5_res)
            md5 = md5_find[0]
            if md5 == acc_prots[file_type + '_md5']:
                eprint(spacer + '\t' + file_type + ': ' + os.path.basename(file_path), flush = True)
                dwnlds[file_type] = 0
                continue
        elif os.path.isfile( file_path[:-3] ):
            eprint(spacer + '\t' + file_type + ': ' + os.path.basename(file_path), flush = True)
            dwnlds[file_type] = 0
            continue

        if ftp_link == '':
            dwnlds[file_type] = 15
            continue
        while True:
            try:
                esc_count += 1
                request = requests.get(ftp_link.replace('ftp://','https://'))
                break
            except requests.exceptions.ConnectionError:
                if esc_count == 20:
                    eprint(spacer + '\tERROR: Failed to fulfill ftp request!', flush = True)
                    sys.exit(71)
                time.sleep(1)
                count = 0
        count += 1
        if request.status_code != 200:
            eprint(spacer + '\t\t' + file_type + ': ERROR, no file exists', flush = True)
            dwnlds[file_type] = 69
            acc_prots[file_type] = ''
            log_editor(output_dir + 'ncbiDwnld.log', str(acc), str(acc) + \
                '\t' + str(acc_prots['accession']) + '\t' + str(acc_prots['assembly']) + \
                '\t' + str(acc_prots['proteome']) + '\t' + str(acc_prots['gff3']) + \
                '\t' + str(acc_prots['transcript']) + '\t' + \
                str(acc_prots['assembly_md5']) + '\t' + \
                str(acc_prots['proteome_md5']) + '\t' + \
                str(acc_prots['gff3_md5']) + '\t' + \
                str(acc_prots['transcript_md5']) + '\t' + \
                str(','.join([str(x) for x in acc_prots['genome_id']])))
            if remove and file_type in {'assembly', 'gff3'}:
                break
            continue

        count += 1
        dwnlds[file_type] = subprocess.call(
            ['curl', ftp_link, '-o', file_path], 
            stdout = subprocess.PIPE, stderr = subprocess.PIPE
            )

        if dwnlds[file_type] != 0:
            eprint(spacer + '\t' + file_type + ': ERROR, download failed', flush = True)
            if remove and file_type in {'assembly', 'gff3'}:
                break
        else:
            if os.stat(file_path).st_size < 150:
                eprint(spacer + '\t' + file_type + ': ERROR, file too small', flush = True)
#                       print('\t' + acc + '\t' + file_type + ' empty', flush = True)
                dwnlds[file_type] = 420
                if remove and file_type in {'assembly', 'gff3'}:
                    break
        eprint(spacer + '\t' + file_type + ': ' + os.path.basename(file_path), flush = True)

    return dwnlds, count


def dwnld_mngr(
    data, acc, file_types, output_path, count, remove, api, spacer
    ):
    fail = []
    exits, count = download_files( 
        data, acc, file_types, output_path, 
        count, remove = remove, spacer = spacer
        )
    count = wait_for_ncbi(count, api)
    t_acc = acc
    try:
        if exits['assembly'] != 0:
            if '$' in acc:
                t_acc = acc[:acc.find('$')]
            fail = [t_acc, ncbi_df['version'][t_acc]]
            return fail, count
    except KeyError:
        pass
    try:
        if exits['gff3'] != 0:
            fail = [t_acc, ncbi_df['version'][t_acc]]
            return fail, count
    except KeyError:
        pass
    return fail, count

def dwnld_mngr_no_MD5(
    data, acc, file_types, output_path, count, remove, api, spacer
    ):
    run, fail = False, []
    for file_type in file_types:
        file_path = output_path + file_type \
            + '/' + os.path.basename(data[file_type])
        if not os.path.isfile(file_path):
            run = True
            break
    if run:
        exits, count = download_files( 
            data, acc, file_types, output_path, 
            count, remove = remove,
            spacer = spacer
            )
        count = wait_for_ncbi(count, api)
    else:
        exits = {x: 0 for x in file_types}
        # we aren't checking md5s, so assume the exit is 0
        return fail, count
    try:
        if exits['assembly'] != 0:
            if '$' in acc:
                t_acc = acc[:acc.find('$')]
            else:
                t_acc = acc
            fail = [t_acc, ncbi_df['version'][t_acc]]
            return fail, count
    except KeyError:
        pass
    try:
        if exits['gff3'] != 0:
             fail = [t_acc, ncbi_df['version'][t_acc]]
             return fail, count
    except KeyError:
        pass
    return fail, count


def main( 
    api = None, 
    assembly = True, proteome = False, gff3 = True, transcript = False,
    ncbi_df = False, remove = False, output_path = os.getcwd(), verbose = False,
    column = 'assembly_acc', ncbi_column = 'Assembly', check_MD5 = True,
    spacer = '\t\t'
    ):

    # initialize run directory and information
    output_path = format_path(output_path)
    file_types = prepare_folders( 
        output_path, gff3, proteome, assembly, transcript
        )
    ass_prots = compile_log(output_path, remove )

    # check if ncbi_df is a dataframe, and import if not
    if not isinstance(ncbi_df, pd.DataFrame) and os.path.isfile(ncbi_df):
        ncbi_df = db2df(ncbi_df)
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

    ncbi_df = ncbi_df.set_index(pd.Index(list(ncbi_df[column])))
    # preserve the original column, but index ncbi_df on it as well
    vprint('\n' + spacer + 'Assembling NCBI ftp directories', v = verbose, flush = True)
    ass_prots, failed, ncbi_df = collect_ftps( 
            ncbi_df, ass_prots, remove = remove,
            ncbi_column = ncbi_column, column = column, api_key=api,
            output_path = output_path, verbose = verbose,
            spacer = spacer
            )

    if remove:
        ass_prots = { 
            o: ass_prots[o] for o in ass_prots \
            if all(ass_prots[o][p] for p in ['assembly', 'gff3'])
            }
    new_df = pd.DataFrame()

    vprint('\n' + spacer + 'Downloading NCBI files', v = verbose, flush = True)
    count = 0
    if check_MD5:
        for acc, data in ass_prots.items():
            eprint(spacer + '\t' + str(acc), flush = True)
            if data:
                fail, count = dwnld_mngr(
                    data, acc, file_types, output_path,
                    count, remove, api, spacer
                    )
                if fail:
                    failed.append(fail)
                else:
                    ncbi_df.at[acc, 'assemblyPath'] = output_path + 'assembly/' + \
                        os.path.basename(ass_prots[acc]['assembly'])
                    ncbi_df.at[acc, 'faa'] = output_path + 'proteome/' + \
                        os.path.basename(ass_prots[acc]['proteome'])
                    ncbi_df.at[acc, 'gffPath'] = output_path + 'gff3/' + \
                        os.path.basename(ass_prots[acc]['gff3'])
                    new_df = pd.concat([new_df, ncbi_df.loc[acc].to_frame().T])
            else:
                check = ncbi_df[ncbi_df[column] == acc[:acc.find('$')]]
                # check for entries in the inputted table that match the accession
                # provided without version modification
                db_vers = datetime.strftime(row['version'], '%Y%m%d')
                failed.append([acc, db_vers])
    else: # there is no checking md5, this is for efficient, so make it
    # efficient by avoiding conditional expressions
        for acc, data in ass_prots.items():
            eprint(spacer + '\t' + str(acc), flush = True)
            fail, count = dwnld_mngr_no_MD5(
                data, acc, file_types, output_path, count, remove, api, spacer
                )
            if fail:
                failed.append(fail)
            else:
                for file_type in file_types:
                    ncbi_df.at[acc, file_type] = output_path + file_type \
                        + '/' + os.path.basename(data[file_type])
                new_df = pd.concat([new_df, ncbi_df.loc[acc].to_frame().T])

    if 'assembly' in new_df.keys():
        new_df = new_df.rename(columns = {'assembly': 'assemblyPath'})
    if 'gff3' in new_df.keys():
        new_df = new_df.rename(columns = {'gff3': 'gffPath'})
    new_df = new_df.reset_index()
    return new_df, failed


def getSRA(assembly_acc, fastqdump = 'fastq-dump', pe = True):

    handle = Entrez.esearch(db='SRA', term=assembly_acc)
    ids = Entrez.read(handle)['IdList']
    for id in ids:
        handle = Entrez.esummary(db='SRA', id = id, report='full')
        records = Entrez.read(handle, validate = False)
        for record in records:
            srr = re.search(r'Run acc="(S\w+\d+)"', record['Runs'])[1]
            print('\t\t' + srr, flush = True)
            if pe:
                subprocess.call([
                    fastqdump, '--gzip', '--split-3', srr], 
                    stdout = subprocess.PIPE)
                if os.path.isfile(srr + '_1.fastq.gz'):
                    shutil.move(srr + '_1.fastq.gz', assembly_acc + '_' + srr + '_1.fq.gz')
                    shutil.move(srr + '_2.fastq.gz', assembly_acc + '_' + srr + '_2.fq.gz')
                else:
                    print('\t\t\tERROR: file failed', flush = True)
            else:
                subprocess.call([
                    fastqdump, '--gzip', srr],
                    stdout = subprocess.PIPE)
                if os.path.isfile(srr + '.fastq.gz'):
                    shutil.move(srr + '.fastq.gz', assembly_acc + '_' + srr + '.fq.gz')
                else:
                    print('\t\t\tERROR: file failed', flush = True)


def goSRA(df, output = os.getcwd() + '/', pe = True):

    print()
    sra_dir = output + 'sra/'
    if not os.path.isdir(sra_dir):
        os.mkdir(sra_dir)
    os.chdir(sra_dir)
    fastqdump = findExecs('fastq-dump', exit = set('fastq-dump'))
    count = 0

    if 'sra' in df.keys():
        row_key = 'sra'
    else:
        row_key = 'assembly_acc'

    for i, row in df.iterrows():
        print('\t' + row[row_key], flush = True)
        getSRA(row[row_key], fastqdump[0])
        count +=1
        if count >= 10:
            time.sleep(1)
            count = 0


def cli():
    parser = argparse.ArgumentParser()
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
    parser.add_argument( '-o', '--output', help = 'Output directory' )
    args = parser.parse_args()

    ncbi_email, ncbi_api, jgi_email, jgi_pwd = loginCheck(jgi = False) 

    Entrez.email = ncbi_email
    if ncbi_api:
        Entrez.api_key = ncbi_api

    if not args.output:
        output = os.getcwd() + '/'
    else:
        output = format_path(args.output)

    args_dict = {
        'NCBI Table': args.input,
        'email': ncbi_email,
        'Assemblies': args.assembly,
        'Proteomes': args.proteome,
        ".gff3's": args.gff3,
        'Transcripts': args.transcript,
        'SRA': args.sra
    }

    start_time = intro('Download NCBI files',args_dict)
    if not args.assembly and not args.proteome and not args.gff3 and not args.sra and not args.transcript:
        eprint('\nERROR: You must choose at least one download option\nExit code 37', flush = True)
        sys.exit( 37 )

    if args.sra:
        if os.path.isfile(format_path(args.input)):
            goSRA(pd.read_csv(format_path(args.input), sep = '\t'), output, pe = args.paired)
        else:
            goSRA(pd.DataFrame({'sra': [args.input.rstrip()]}), output, pe = args.paired)
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
            ncbi_df = pd.DataFrame({'assembly_acc': args.input.replace('"','').replace("'",'').split()})
            column = 'assembly_acc'
            ncbi_column = 'assembly'
    
        new_df, failed = main( 
            assembly = args.assembly, column = column, 
            ncbi_column = ncbi_column, proteome = args.proteome, 
            gff3 = args.gff3, transcript = args.transcript, ncbi_df = ncbi_df,
            output_path = output, verbose = True, spacer = ''
            )
        new_df.to_csv( args.input + '_dwnld', sep = '\t' )

    outro(start_time)


if __name__ == '__main__':
    cli()
