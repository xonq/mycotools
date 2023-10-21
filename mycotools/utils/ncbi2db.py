#! /usr/bin/env python3


# NEED to refine MTDB to consider GCA and the version in the assembly accession
#   which is important for version parsing and keeping them uniform with
#   biosample formatting in MTDB
#   NEED TO EDIT REDUNDANCY CHECK TO REFERENCE QUERIED ASSEMBLY ACCESSIONS FROM
#   BIOSAMPLES

import os
import re
import sys
import copy
import argparse
import subprocess
import numpy as np
import pandas as pd
from datetime import datetime
from mycotools.lib.kontools import intro, outro, eprint
from mycotools.lib.dbtools import db2df, df2db, primaryDB
from mycotools.ncbiDwnld import main as ncbi_dwnld
from mycotools.predb2mtdb import main as predb2mtdb

def redundancy_check(db, ncbi_df, ass_acc, duplicates = {}):
    '''should expect that ncbi_df (pd.DataFrame()) is only comprised of biosamples that are
    not overlapping with JGI'''
    db['index'], updates, old_rows = db.loc[:, 'ome'], [], {}
    db['versionless'] = [x[:x.find('.')] for x in list(db['assembly_acc'])]
    assembly_accs = set(db['versionless'])
    # prepare a set of assembly accessions without the version tag
    todel_db, todel_ncbi = [], []
    for i, row in ncbi_df.iterrows(): # totally based on assembly accession and version
        organism = '_'.join([row['genus'], row['species'], str(row['strain'])])
        assembly_acc = row[ass_acc]
        assembly_acc_check = \
            row[ass_acc][:row[ass_acc].find('.')]
        if assembly_acc_check in duplicates: # special case duplicate
            todel_ncbi.append(i)
        elif assembly_acc_check in assembly_accs: # assembly accession is in the
        # reference database
            check = db[db['versionless'] == assembly_acc_check]
            if len(check.index) > 1: # multiple versions of the same accession
            # in reference
                versions = [
                    int(re.search(r'\.(\d+)$', row['assembly_acc'])[1]) \
                    for i, row in check.iterrows()
                    ]
                indices = list(check.index)
                max_v = max(versions)
                if versions.count(max_v) > 1: # multiple max possible version
                    max_indices = [i for i,v in enumerate(indices) \
                                   if versions[i] == max_v]
                    new_check = check.iloc[max_indices]
                    dates = list(new_check['version'])
                    date_indices = list(new_check.index)
                    max_date = max(dates)
                    max_i = date_indices[dates.index(max_date)] # will remove
                    # duplicate versions with duplicate modification dates
                else:
                    max_i = indices[versions.index(max_v)]
                more_dels = [check.loc[v, 'ome'] \
                            for i, v in enumerate(indices) \
                            if v != max_i]
                todel_db.extend(more_dels) # delete the overlap
                db_row = copy.deepcopy(check.loc[max_i])
            else:
                db_row = copy.deepcopy(check.iloc[0])
            db_ome = db_row['ome']
            db_organism = '_'.join([db_row['genus'], db_row['species'], str(db_row['strain'])])
            version = row['version']
            try: # check if the version is new (by datetime)
                if version > datetime.strptime(db_row['version'].replace(' 00:00:00','').replace('-',''), '%Y%m%d'):
                    updates.append([
                        db_row['assembly_acc'], db_organism,
                        organism, assembly_acc, db_ome
                        ])
                    todel_db.append(db_ome)
                    old_rows[assembly_acc] = db_row
                else:
                    todel_ncbi.append(i)
            except ValueError:
                todel_ncbi.append(i)
        else: # new entry
            updates.append([None, None, organism, assembly_acc, None])

    db = db.set_index('index')
    todel_db.sort(reverse = True)
    todel_ncbi.sort(reverse = True)
    for i in set(todel_db):
        db = db.drop(i)
    for i in todel_ncbi:
        ncbi_df = ncbi_df.drop(i)
 
    return ncbi_df, db, updates, old_rows

   
def main( 
    out_dir, ncbi_df = None, date = datetime.today().strftime('%Y%m%d'),
    ref_db = None, assem = True, prot = False, gff = True, transcript = False,
    rerun = False, failed_dict = {}, duplicates = {}, check_MD5 = True,
    spacer = '\t\t'
    ):

    os.chdir( out_dir )
    todel = []

    if 'assembly_acc' in ncbi_df.columns:
        ass_acc = 'assembly_acc'
    else:
        ass_acc = 'Assembly Accession'
    if not rerun:
        ncbi_df = ncbi_df.set_index(ass_acc)
        if '-' in ncbi_df.keys():
            ncbi_df = ncbi_df.drop('-')
        prev_omes = set(ncbi_df.index)
        for failure in failed_dict:
            if failure in prev_omes:
                if failed_dict[failure]['version']:
                    version = ncbi_df['version'][failure]
                    prev_version = datetime.strptime(
                        failed_dict[failure]['version'], '%Y%m%d'
                        )
                    if not version > prev_version:
                        todel.append(failure)
        for i in todel:
            ncbi_df = ncbi_df.drop(i)
        ncbi_df = ncbi_df.reset_index()

    print(spacer + 'Redundancy check', flush = True)
    update_check = {}
    if ref_db is not None:
        if len(ref_db) > 0:
            ncbi_df, ref_db, updates, old_rows = redundancy_check( 
                ref_db, ncbi_df, ass_acc
                )
            output_str = 'old_assembly_acc\tdb_organism\tncbi_organism\tnew_assembly_acc\n'
            for update in updates:
                output_str += '\t'.join([str(x) for x in update]) + '\n'
            with open(out_dir + '/ncbiUpdates.tsv', 'w') as out:
                out.write(output_str)
            update_check = {i[-2]: i for i in updates if i[0]}
            # dict(update_check) = {assembly_accNEW: [organism, ref organism,
            # old_assembly_acc]}
            print(spacer + '\t' + str(len(ncbi_df)) + ' genomes to assimilate', flush = True)

    if len(ncbi_df) > 0:
        print(spacer + 'Downloading NCBI data', flush = True)
        ncbi_df, failed = ncbi_dwnld(
            assembly = assem, proteome = prot,
            gff3 = gff, ncbi_df = ncbi_df, remove = True, output_path = out_dir,
            column = ass_acc, ncbi_column = 'genome',
            check_MD5 = check_MD5
            )
    
        print(spacer + '\t' + str(len(ncbi_df)) + ' entries with assemblies and gffs', flush = True)
        ncbi_df = ncbi_df.rename(columns = {
            'Release Date': 'published', 'Assembly Accession': 'assembly_acc',
            'BioSample Accession': 'biosample'
            })
        ncbi_df['source'] = 'ncbi'
        for assembly_acc_info in failed:
            if assembly_acc_info[0] in update_check:
                ref_db = ref_db.append(old_rows[assembly_acc_info[0]])
                del update_check[assembly_acc_info[0]] # remove it from potential updates
    
        ref_db = ref_db.set_index('assembly_acc') # update ome codes
      #  try:
        ncbi_df = ncbi_df.set_index('assembly_acc') 
    #    except KeyError: # no entries
     #       return ncbi_df, ref_db.reset_index(), failed, duplicates
        for assembly_acc, update_d in update_check.items():
             old_ome = update_d[-1]
             ncbi_df.at[assembly_acc, 'ome'] = old_ome
    
        return ncbi_df.reset_index(), ref_db.reset_index(), failed
    else:
        return ncbi_df, ref_db, [], duplicates


if __name__ == '__main__':
    cli()
