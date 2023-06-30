#! /usr/bin/env python3

import os
import re
import sys
import copy
import time
import getpass
import datetime
import argparse
import subprocess
import pandas as pd
import numpy as np
from io import StringIO
from mycotools.predb2mtdb import main as predb2mtdb
from mycotools.lib.kontools import intro, outro
from mycotools.lib.dbtools import db2df, df2db, readLog, log_editor
from mycotools.jgiDwnld import jgi_login as jgi_login
from mycotools.jgiDwnld import retrieveXML as retrieveXML
from mycotools.jgiDwnld import jgi_dwnld as jgi_dwnld


def compileLog( log_path ):

    log = {}
    if not os.path.isfile(log_path):
        with open(log_path, 'w') as out:
            out.write('#assembly_acc\tfna\tgff3\tfaa')
    else:
        log = readLog(log_path)

    return log



def break_name( jgi_df, name_col = 'name' ):

    for i, row in jgi_df.iterrows():
        taxa = row[name_col].split(' ')
        jgi_df.at[i, 'genus'] = taxa[0]
        if len(taxa) > 1:
            jgi_df.at[i, 'species'] = taxa[1]
            if len( taxa ) > 2:
                jgi_df.at[i, 'strain'] = ''.join(taxa[2:])
                vers_search = re.search(r'v(\d+\.\d+$)', jgi_df['strain'][i])
                if vers_search is not None:
                    version = vers_search[1]
                    jgi_df.at[i, 'version'] = version
                    jgi_df.at[i, 'strain'] = jgi_df['strain'][i][:-len(vers_search[0])]
                else:
                    jgi_df.at[i, 'version'] = 0
                jgi_df.at[i, 'strain'] = re.sub(r'[^a-zA-Z0-9]', '', jgi_df.at[i, 'strain'])
        else:
            jgi_df.at[i, 'species'] = 'sp.'

    jgi_df = jgi_df.sort_values(by = 'version', ascending = False)
    jgi_df = jgi_df.drop_duplicates('portal')

    return jgi_df


def JGIredundancyCheck(db, jgi_df, duplicates = {}, ome_col = 'portal',
                       jgi2ncbi = {}):
    '''everything in jgi_df (pd.DataFrame()) should be have a valid biosample 
       at this point.
       db will have the index set at "assembly_acc"'''

    db['version'].fillna(0.0)
    jgi_df['version'] = jgi_df.loc[:, 'version'].replace( np.nan, '0.0' )
    preOmes, biosamples = set( copy.deepcopy(db.index)), set(db['biosample'])
    updates, old_omes, todel_db, todel_jgi = {}, {}, [], []

    for i, row in jgi_df.iterrows():
#        version = float(row['version'].replace('v',''))
        version = float(row['version'])
        ome, biosample = copy.deepcopy(row[ome_col]), row['biosample'] # implied to exist
        if ome in duplicates:
            todel_jgi.append(i)
        elif ome in preOmes:
            db.at[ome, 'published'] = copy.deepcopy(row['publication(s)']) 
            if pd.isnull(db['version'][ome]):
                db.at[ome, 'version'] = 0
            # update most current pub status
            if not db['version'][ome] or pd.isnull(db['version'][ome]):
                db.at[ome, 'version'] = 0
                db_version = 0
            else:
                try:
                    db_version = float(db['version'][ome])
                except (ValueError, AttributeError) as e:
                    db_version = float(db['version'][ome].replace('v',''))
                    db.at[ome, 'version'] = db_version
            if version > db_version and db['source'][ome] == 'jgi':
                db_organism = db.loc[ome, 'genus'] + '_' \
                    + db.loc[ome, 'species'] + '_' + db.loc[ome, 'strain']
                organism = jgi_df['name'].replace(' ', '_')
                updates[ome] = [db.loc[ome, 'ome'], db_organism, organism]
                old_omes[db.loc[ome, 'ome']] = copy.deepcopy(db.loc[ome])
                # retain row information
#                jgi_df.at[i, 'ome'] = db.at[ome, 'ome'] # keep the same ome code.
                # I think this ought to be updated to a new code
                todel_db.append(ome)
            else:
                todel_jgi.append(i)
        elif biosample in biosamples: # genome_code isn't here, but the biosample is
            check_db = db[db['biosample'] == biosample]
            checkOme = list(check_db['ome'])[0]
            if len(check_db) > 1: # more than one entry to the biosample, keep as it is
                todel_jgi.append(i)
            else: # this is an NCBI accession that can be updated to JGI
                index = list(check_db.index)[0]
                db_organism = db.loc[index, 'genus'] + '_' \
                    + db.loc[index, 'species'] + '_' \
                    + db.loc[index, 'strain']
                organism = jgi_df['name'].replace(' ', '_')
                updates[ome] = [checkOme, db_organism, organism]
                todel_db.append(index)
        elif ome.lower() in jgi2ncbi: # NCBI can be updated to JGI
            low_ome = ome.lower()
            ncbi_acc = jgi2ncbi[low_ome]
            if ncbi_acc in set(db['assembly_acc']):
                print(low_ome, ncbi_acc)
                sys.exit(1)
                check_db = db[db['assembly_acc'] == ncbi_acc]
                checkOme = list(check_db['ome'])[0]
                index = list(check_db.index)[0]
                db_organism = db.loc[index, 'genus'] + '_' \
                    + db.loc[index, 'species'] + '_' \
                    + db.loc[index, 'strain']
                organism = jgi_df['name'].replace(' ', '_')
                updates[ome] = [checkOme, db_organism, organism]
                todel_db.append(index)
            else:
                organism = jgi_df['name'].replace(' ', '_')
                updates[ome] = [None, None, organism]
        else:
            organism = jgi_df['name'].replace(' ', '_')
            updates[ome] = [None, None, organism]

    todel_db.sort(reverse = True)
    todel_jgi.sort(reverse = True)
    for i in todel_db:
        db = db.drop(i)
    for i in todel_jgi:
        jgi_df = jgi_df.drop(i)

    db = db.reset_index()

    return jgi_df, db, updates, old_omes, duplicates


def runjgi_dwnld(
    jgi_df, i, user, pwd, ome_set,
    ome_col, output, log, log_path,
    dwnlds, failed, rerun, masked,
    spacer = '\t\t'
    ):

    ome = jgi_df[ome_col][i]
    jgi_login(user, pwd)

    ran_dwnld = False       
    if ome not in ome_set:
        print(spacer + '\t' + ome + ': ' + jgi_df['name'][i], flush = True)
        for typ in dwnlds:
            if log[ome][typ] == 'error':
                if not rerun:
                    print(spacer + '\t\t' + typ + ': ERROR', flush = True)
                    continue
            check, preexisting, new_typ, ran_dwnld = jgi_dwnld( 
                ome, typ, output, masked = masked
            )
            if type(check) != int:
                jgi_df.at[i, new_typ + '_path'] = output + '/' + new_typ + \
                    '/' + check
                check = os.path.basename( os.path.abspath( check ) )
                print(spacer + '\t\t' + new_typ + ': ' + str(check), flush = True)
                log[ome][typ] = check
            else:
                print(spacer + '\t\t' + new_typ + ': ERROR', flush = True)
                log[ome][typ] = 'error'
                if typ in {'gff3', 'fna'}:
                    log_editor(
                        log_path, ome, ome + '\t' + log[ome]['fna']
                        + '\t' + log[ome]['gff3'] + '\t' +
                        log[ome]['faa']
                        )
                    failed.append([ome, jgi_df['version'][i]])
                    jgi_df = jgi_df.drop(i)
                    if ran_dwnld:
                        time.sleep(60)
                    break
            if ran_dwnld:
                time.sleep(60)
        log_editor(
            log_path, ome, ome + '\t' + log[ome]['fna'] + '\t'
            + log[ome]['gff3'] + '\t' + log[ome]['faa']
            )

    return jgi_df, log, failed


def log2df(jgi_df, log, output):

    typ2col = {'fna': 'assemblyPath', 'gff': 'gffPath',
               'gff3': 'gffPath', 'faa': 'proteomePath'}
    for ome in list(jgi_df.index):
        for typ in log[ome]:
#            if log[ome][typ].endswith('.gff.gz'):
 #               out_file = output + '/gff/' + log[ome][typ]
  #              col = 'jgi_gff2_path'
   #         else:
            out_file = output + '/' + typ + '/' + log[ome][typ]
            jgi_df.at[ome, typ2col[typ]] = out_file

    return jgi_df


def main( 
    jgi_df, ref_db, output, user, pwd, date = datetime.datetime.now().strftime('%Y%m%d'), 
    assembly = True, proteome = False, gff3 = True, update = True,
    repeatmasked = True, nonpublished = False, rerun = False, failed_dict = {},
    duplicates = {}, spacer = '\t\t', jgi2ncbi = {}
    ):

    if not nonpublished:
        jgi_df = jgi_df[jgi_df['is published'] == 'Y']
    if 'assembly_acc' in jgi_df.columns:
        ome_col = 'assembly_acc'
    elif 'portal' in jgi_df.columns:
        ome_col = 'portal'
        name_col = 'name'
    else:
        print(spacer + 'ERROR: invalid headers for ' + mycocosm_output, flush = True)
        sys.exit( 3 )

    toDel = []
    jgi_df = break_name( jgi_df )
   
    if not rerun:
        jgi_df = jgi_df.set_index(ome_col)
        jgi_omes = set(jgi_df.index)
        for failure in failed_dict:
            if failure in jgi_omes:
                if 'v' in str(jgi_df['version'][failure]):
                    version = float(jgi_df['version'][failure].replace('v',''))
                    if failed_dict[failure]['version'] != '':
                        old_vers = float(failed_dict[failure]['version'].replace('v',''))
                        if not version > old_vers:
                            toDel.append(failure)
                else:
                    toDel.append(failure)
        for failed in toDel:
            jgi_df = jgi_df.drop(failed)
        jgi_df = jgi_df.reset_index()

    print(spacer + 'Redundancy check', flush = True)
    if isinstance(ref_db, pd.DataFrame):
        ref_db['index'] = ref_db['assembly_acc'].copy()
        ref_db = ref_db.set_index( 'index' )
        old_len = len(jgi_df)
        jgi_df, new_ref_db, updates, old_rows, duplicates = JGIredundancyCheck( 
            ref_db, jgi_df, duplicates = duplicates, ome_col = ome_col,
            jgi2ncbi = jgi2ncbi
            )
        output_str = 'ome\tdb_organism\tjgi_organism\tassembly_acc\n'
        if isinstance(updates, dict):
            for ome, update in updates.items():
                output_str += '\t'.join([str(x) for x in update]) + '\t' + ome + '\n'
            with open( output + '/jgiUpdates.tsv', 'w' ) as out:
                out.write( output_str )
        else:
            updates.to_csv(f'{output}/jgiUpdates.tsv', sep = '\t')
#        update_check = {i[-1]: i[0:3] for i in updates if i[0]}
        print(spacer + '\t' + str(len(jgi_df)) + ' genomes to assimilate' , flush = True)
    else:
        new_ref_db, updates = None, {}


    print(spacer + 'Logging into JGI', flush = True)
    jgi_login( user, pwd )

    if not os.path.exists( output + '/xml' ):
        os.mkdir( output + '/xml' )

    print(spacer + 'Retrieving `xml` directories' , flush = True)
    ome_set, failed, count = set(), [], 0
    for i, row in jgi_df.iterrows():
        error_check = retrieveXML( row[ ome_col ], output + '/xml' )
        if error_check > 0:
            ome_set.add( row[ ome_col ] )
        if error_check != -1:
            time.sleep( 0.1 )


    log_path = output + '/jgi2db.log'
    log = compileLog( log_path )
    if not rerun:
        prev_omes = set(jgi_df[ome_col])
        for ome in log:
            if ome not in prev_omes:
                continue
            if log[ome]['fna'] == 'error' or log[ome]['gff3'] == 'error':
                drop_index = list(jgi_df[jgi_df[ome_col] == ome].index)[0]
                failed.append([ome, jgi_df['version'][drop_index]])
                jgi_df = jgi_df.drop(drop_index)
                ome_set.add(ome)


    print(spacer + 'Downloading JGI data' , flush = True)
    dwnlds = []
    if assembly:
        dwnlds.append( 'fna' )
    if gff3:
        dwnlds.append( 'gff3' )
    if proteome:
        dwnlds.append( 'faa' )

    for typ in dwnlds:
        if not os.path.isdir( output + '/' + typ ):
            os.mkdir( output + '/' + typ )
        if typ == 'gff3':
            if not os.path.isdir( output + '/gff3' ):
                os.mkdir( output + '/gff3')

    if all(x in log for x in list(jgi_df[ome_col])) and not rerun:
        print(spacer + '\tAll downloaded, rerun off', flush = True)
        jgi_df = jgi_df.set_index(ome_col)
        jgi_df = log2df(jgi_df, log, output)
        jgi_df = jgi_df.reset_index()
    else:
        for i, row in jgi_df.iterrows():
            ome = row[ome_col]
            if ome not in log:
                log[ome] = {'fna': 'na', 'gff3': 'na', 'faa': 'na'}
            elif ome not in set(jgi_df[ome_col]):
                continue
            jgi_df, log, failed = runjgi_dwnld(
                jgi_df, i, user, pwd, ome_set,
                ome_col, output, log, log_path,
                dwnlds, failed, rerun, repeatmasked,
                spacer
                )

    if os.path.exists( 'cookies' ):
        os.remove( 'cookies' )
    if os.path.exists( os.path.expanduser( '~/.nullJGIdwnld' )):
        os.remove( os.path.expanduser( '~/.nullJGIdwnld' ))

    jgi_df = jgi_df.rename( columns = { 
        'published(s)': 'published',
        ome_col: 'assembly_acc',
        'fna': 'assemblyPath',
        'gff3': 'gffPath'
        } )
    jgi_df['source'] = 'jgi'

    for ome_d in failed: # add back the failed entries that were attempted updates
        ome = ome_d[0]
        if ome in old_rows: # if there is an old row to add back
            new_ref_db = new_ref_db.append(old_rows[ome])

    old_assembly_accs = set(new_ref_db['assembly_acc'])
    jgi_df = jgi_df.set_index('assembly_acc')
    new_ref_db = new_ref_db.set_index('assembly_acc')
    for ome in updates: # apply the old ome code to update the version
        if ome in old_assembly_accs:
            jgi_df.at[ome, 'ome'] = new_ref_db.loc[ome, 'ome']

    jgi_df = jgi_df.astype(str)
    jgi_premtdb_df = jgi_df.reset_index()
    for i, row in jgi_premtdb_df.iterrows():
        if not pd.isnull(row['publication(s)']) and row['publication(s)']:
            jgi_premtdb_df.at[i, 'published'] = row['publication(s)']
        elif not pd.isnull(row['is published']) and row['is published']:
            jgi_premtdb_df.at[i, 'published'] = 1
        elif not pd.isnull(row['is public']) and row['is public']:
            jgi_premtdb_df.at[i, 'published'] = 1

    return jgi_premtdb_df, new_ref_db.reset_index(), failed, duplicates


def cli():

    mycoCosmURL = 'https://mycocosm.jgi.doe.gov/ext-api/mycocosm/catalog/' + \
            'download-group?flt=&seq=all&pub=all&grp=fungi&srt=released&ord=desc'
    parser = argparse.ArgumentParser( description = ' Downloads or imports MycoCosm table.' \
        + ' Downloads files and outputs a mycotools `.db` of novel JGI downloads relative to' \
        + ' a reference `.db`. Can also begin a new `.db`'
    )
    parser.add_argument( '-l', '--login', 
        help = r'Login file: "<JGI username>\t<JGI Password>\n<NCBI email>\t<NCBI API key>"'
        )
    parser.add_argument( '-d', '--database', help = 'Existing myctools `.db` to reference' )
    parser.add_argument( 
        '-m', '--mycocosm', default = mycoCosmURL, 
        help = 'MycoCosm master table path / URL. This may change from time to time ' \
        + 'DEFAULT: ' + mycoCosmURL
        )
    parser.add_argument(
        '-u', '--update', help = 'Download potential JGI updates for existing omes ' \
        '- NOT FUNCTIONAL FOR DATABASE SAFETY',
        default = False, action = 'store_true'
        )
    parser.add_argument(
        '-o', '--output', help = 'Output directory. DEFAULT: Current date' 
        )
    parser.add_argument( '-a', '--assembly', default = True, action = 'store_false', \
        help = 'Do not download assemblies.' )
    parser.add_argument( '-p', '--proteome', default = True, action = 'store_false', \
        help = 'Do not download proteomes.' )
    parser.add_argument( '-g', '--gff', default = True, action = 'store_false', \
        help = 'Do not download gff3s.' )
    parser.add_argument( '-r', '--repeatmasked', default = True, action = 'store_false', \
        help = "Download nonmasked assemblies." )

    args = parser.parse_args()

    args_dict = {
            'Preexisting db': args.database,
            'MycoCosm Table': args.mycocosm,
            'Assemblies': args.assembly,
            'RepeatMasked': args.repeatmasked,
            'Proteomes': args.proteome,
            ".gff3's": args.gff,
            }

    start_time = intro( 'jgi2db', args_dict )
    if args.output:
        output = os.path.abspath( args.output )
    else:
        output = start_time.strftime('%Y%m%d') + '_jgi2db' 
    if not os.path.isdir( output ):
        os.mkdir( output )

    if args.login:
        with open( args.login, 'r' ) as raw:
            prep = raw.read()
        data = [ x.split('\t') for x in prep.split('\n') ]
        user = data[0][0]
        pwd = data[0][1]
        email = data[1][0]
        apikey = None
        if len(data[1]) > 1:
            if data[1][1] != '':
                apikey = data[1][1]

    ref_db = db2df( format_path(args.database) )
    jgi_df = main( args.mycocosm, refdb, output )

    df2db( jgi_df, output + '/new.db' )
    print('\nSuccess! ' + str(len(jgi_df), flush = True) + ' added to database\n \
            Run updateDB to confirm and finish update.')

    outro( start_time )


if __name__ == '__main__':
    cli()
