#! /usr/bin/env python3

import argparse, sys, os, re, pandas as pd, datetime, numpy as np, subprocess
from mycotools.lib.dbtools import db2df, df2db, gen_omes, masterDB
from ncbiDwnld import main as ncbiDwnld
from predb2db import main as predb2db
from mycotools.lib.kontools import intro, outro, eprint


def dwnldNCBItable( ncbi_url, out_dir ):

    date = datetime.datetime.today().strftime('%Y%m%d')
    ncbi_file = out_dir + '/' + date + '_ncbi.tsv'
    getTbl = subprocess.call( 
        'curl ' + ncbi_url \
        + ' > ' + ncbi_file, shell = True )
    if getTbl != 0:
        eprint('\nERROR downloading NCBI table')
        sys.exit(13)
    ncbi_df = pd.read_csv( ncbi_file, sep = '\t' )

    return ncbi_df


def redundancyCheck( db, ncbi_df, overwrite = False ):

    db['strain'] = db.apply(lambda x: re.sub(r'_?v\d+\.\d+$', '', str(x['strain'])), \
        axis = 1)
    db['strain_check'] = db.apply(lambda x: str(x['genus']) + '_' + str(x['species']) + \
        '_' + x['strain'], axis = 1)
    db['strain_check1'] = db.apply(lambda x: str(x['genus']) + '_' + str(x['species']) + \
        '_' + x['strain'].replace('_','').replace('-',''), axis = 1)
    db['strain_check2'] = db.apply(lambda x: str(x['genus']) + '_' + str(x['species']) + \
        '_' + x['strain'].replace('_',''), axis = 1)
    db['strain_check3'] = db.apply(lambda x: str(x['genus']) + '_' + str(x['species']) + \
        '_' + x['strain'].replace('-',''), axis = 1)
    db['species_check'] = db.apply(lambda x: str(x['genus']) + '_' + str(x['species']), \
        axis = 1)
    strain_set, species_set = set(db['strain_check']), set(db['species_check'])
    strain_set = strain_set.union( set(db['strain_check1']) )
    strain_set = strain_set.union( set(db['strain_check2']) )
    strain_set = strain_set.union( set(db['strain_check3']) )
    toDel = []
    for i, row in ncbi_df.iterrows():
        organism = row['#Organism/Name'].split(' ')
        organism = [ x.replace('[','').replace(']','').replace('(','').replace(')','') for x in organism ]
        ncbi_df.at[i, 'genus'] = organism[0]
        if len(organism) > 1:
            ncbi_df.at[i, 'species'] = organism[1]
            if len(organism) > 2:
                ncbi_df.at[i, 'strain'] = '_'.join(organism[2:])
        else:
            ncbi_df.at[i, 'species'] = 'sp.'
        name = ncbi_df['genus'][i] + '_' + ncbi_df['species'][i]
        if not pd.isnull(row['strain']) and not row['strain'] == '':
            name += '_' + row['strain'].replace(' ','_')
            if name in strain_set:
                toDel.append( i )
            else:
                print( name )
        elif name in species_set:
            toDel.append( i )

    for i in toDel:
        ncbi_df = ncbi_df.drop( i )


    return ncbi_df, toDel
    

    
def main( 
    out_dir, email = '', api = None,
    ncbi_df = 'https://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/eukaryotes.txt', 
    ref_db = None, assem = True, prot = True, gff = True, transcript = False 
    ):

    if not type(ncbi_df) == pd.DataFrame:
        if os.path.isfile( ncbi_df ):
            ncbi_df = pd.read_csv( ncbi_df, sep = '\t', encoding = 'utf-8' )
        else:
            ncbi_df = dwnldNCBItable( ncbi_df, out_dir )

    ncbi_df['strain'] = ''
    if 'Group' in ncbi_df.columns:
        ncbi_df = ncbi_df[ncbi_df['Group'] == 'Fungi']
    if 'Genes' in ncbi_df.columns:
        ncbi_df['Genes'] = str(ncbi_df['Genes'])
        ncbi_df['Proteins'] = str(ncbi_df['Proteins'])
        ncbi_df = ncbi_df[ncbi_df['Genes'] != '-'] 
        ncbi_df = ncbi_df[ncbi_df['Genes'] != '0']
        ncbi_df = ncbi_df[ncbi_df['Proteins'] != '-']
        ncbi_df = ncbi_df[ncbi_df['Proteins'] != '0']
        ncbi_df = ncbi_df[ncbi_df['Genes'] != 'NaN']
        ncbi_df = ncbi_df[ncbi_df['Genes'] != '']
        ncbi_df = ncbi_df[ncbi_df['Proteins'] != 'NaN'] 
        ncbi_df = ncbi_df[ncbi_df['Proteins'] != '']
    ncbi_df = ncbi_df[ncbi_df['Assembly Accession'] != np.nan]
    ncbi_df = ncbi_df[ncbi_df['Assembly Accession'] != '']
    err_df = pd.DataFrame(ncbi_df['BioSample Accession'])
    os.chdir( out_dir )

    print('\nRedundancy check\n\t' + str(len(ncbi_df)) + ' initial entries')
    if ref_db is not None:
        ncbi_df, redundErrors = redundancyCheck( ref_db, ncbi_df )
        print('\t' + str(len(ncbi_df)) + ' following redundancy check')
        ncbi_df = gen_omes( ncbi_df, reference = ref_db )
    else:
        ncbi_df = gen_omes( ncbi_df )

    print('\nDownloading data')
    ncbi_df = ncbiDwnld(
        email = email, api = api, assembly = assembly, proteome = proteome,
        gff3 = gff3, ncbi_df = ncbi_df, remove = True, output_path = out_dir 
    )

    print('\t' + str(len(ncbi_df)) + ' entries with proteomes and gffs')
    ncbi_df['assembly_path'] = ncbi_df.apply(
        lambda x: out_dir + '/assembly/' + x['internal_ome'] + '_assembly.fasta.gz', \
        axis = 1 )
    ncbi_df['proteome_path'] = ncbi_df.apply(
        lambda x: out_dir + '/proteome/' + x['internal_ome'] + '_proteome.aa.fasta.gz', \
        axis = 1 )
    ncbi_df['gff3_path'] = ncbi_df.apply(
        lambda x: out_dir + '/gff3/' + x['internal_ome'] + '.gff3.gz', axis = 1
    )
    ncbi_df = ncbi_df.rename(columns = {'Release Date': 'publication', 'BioSample Accession': 'biosample'})
    ncbi_df['genome_code'] = pd.Series(ncbi_df['internal_ome'])
    ncbi_df['source'] = 'ncbi'

    ncbi_df = predb2db( ncbi_df, ref_db )

    return ncbi_df


if __name__ == '__main__':

    parser = argparse.ArgumentParser( description = 'Downloads and assimilates NCBI ' + \
        'tables into the database. Recommended to run following `jgi2db`.' )
    parser.add_argument( '-i', '--input', help = 'NCBI table (.tsv). DEFAULT: ' + \
        'https://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/eukaryotes.txt' )
    parser.add_argument( '-d', '--database', default = masterDB(), help = 'Existing DB to reference' ) 
    parser.add_argument( '-o', '--output', default = './ncbiDwnld', \
        help = 'Working directory output' )
    parser.add_argument( '-e', '--email', help = 'NCBI email for taxonomy' )
    parser.add_argument( '-k', '--key', help = 'NCBI API key to expedite queries/min' )
    args = parser.parse_args()

    args_dict = { 'Input': args.input, 'Database': args.database, 'Output': args.output,
        'Email': args.email, 'API Key': args.key }

    cur_dir = os.path.abspath( './' )
    out_dir = os.path.abspath( args.output )
    if not os.path.isdir( out_dir ):
        os.mkdir( out_dir )
    if args.database:
        ref_db = db2df( args.database )
    else:
        ref_db = None

    start_time = intro( 'NCBI to Database', args_dict )
    ncbi_df = main( out_dir, ncbi_df = args.input, email = args.email, api = args.api, ref_db = ref_db )

    df2db( ncbi_df, out_dir + '/new.db' )
    print('\nSuccess! ' + str(len(ncbi_df)) + ' added to database\n \
        Run updateDB to confirm and finish update.')

    outro( start_time )
