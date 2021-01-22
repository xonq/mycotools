#! /usr/bin/env python3

import pandas as pd, pandas
import numpy as np
import sys, os, time, mycotools.lib.fastatools, json, re, datetime
from mycotools.lib.kontools import collect_files, eprint, formatPath
from Bio import Entrez
from urllib.error import HTTPError


def genConfig( branch = 'stable' ):

    config = { 
        'forbidden': '$MYCODB/log/forbidden.tsv',
        'repository': 'https://gitlab.com/xonq/mycodb', 
        'branch': branch,
        'nonpublished': False,
        'rogue': False,
    }

    return config

# opens a `log` file path to read, searches for the `ome` code followed by a whitespace character, and edits the line with `edit` 
def log_editor( log, ome, edit ):

    with open( log, 'r' ) as raw:
        data = raw.read()

    if re.search( ome + r'\s', data ):
        new_data = re.sub( ome + r'\s.*', edit, data )
    else:
        data = data.rstrip()
        data += '\n'
        new_data = data + edit

    with open( log, 'w' ) as towrite:
        towrite.write(new_data)

def readLog( log, headers = '', sep = '\t' ):

    log_dict = {}
    with open( log, 'r' ) as raw:
        data = raw.read()
    dataLines = [ x.split( sep ) for x in data.split( '\n' ) ]

    if type( headers ) is str:
        if dataLines[0][0].startswith( '#' ):
            dataLines[0][0] = re.sub( r'^\#+', '', dataLines[0][0] )
        headers = [ head for head in dataLines[0] ]
        del dataLines[0]
    elif type( headers ) is int:
        x = 0
        headers = [ x + 1 for y in dataLines[0] ]

    for line in dataLines:
        log_dict[ line[0] ] = { headers[ x ]: line[ x ] for x in range( 1, len( line ) ) }

    return log_dict


def masterDB( path = '$MYCODB' ):

    path = path.replace( '$', '' )
    try:
        full_path = os.environ[ path ]
        files = collect_files( full_path, 'db' )
        basenames = [ os.path.basename( x ) for x in files ]
        dates = [ x.replace('.db','') for x in basenames if re.search( r'^\d+\.db$', x ) ]
        master = '19991231'
        for date in dates:
            dtDate = datetime.datetime.strptime( date, "%Y%m%d" )
            dtMast = datetime.datetime.strptime( master, "%Y%m%d" )
            if dtDate > dtMast:
                master = date
        if master == '19991231':
            eprint('\nERROR: master database not found in ' + path)
        master_path = formatPath( '$' + path + '/' + master + '.db' )
    except KeyError:
        eprint('\nERROR: ' + path + ' not in path' )
        master_path = None

    return master_path
 
# imports database, converts into df
# returns database dataframe
def db2df(db_path):

    db_path = formatPath( db_path )
    db_df = pd.read_csv( db_path, sep='\t' )
    if 'internal_ome' not in set( db_df.columns ) and 'genome_code' not in set( db_df.columns ):
        headers = [ 'internal_ome', 'genus', 'species', 'strain', 'biosample',
            'assembly', 'proteome', 'gff3',
            'taxonomy', 'ecology', 'eco_conf', 'source', 'published',
            'blastdb', 'genome_code' ]
        db_df = pd.read_csv( db_path, sep = '\t', names = headers )

    return db_df

def df2std( df ):
    '''
    Standardized organization of database columns
    '''
    # temporary check, will remove when all scripts account for new biosample column
    if 'biosample' in df.columns:
        trans_df = df[[
        'internal_ome', 'genus', 'species',
        'strain', 'biosample', 'assembly', 
        'proteome', 'gff3', 'taxonomy',
        'ecology', 'eco_conf', 'source', 'published', 
        'blastdb', 'genome_code'
    ]]
    else:
        trans_df = df[[
        'internal_ome', 'genus', 'species',
        'strain', 'assembly', 
        'proteome', 'gff', 'gff3', 'taxonomy',
        'ecology', 'eco_conf', 'source', 'published', 
        'blastdb', 'genome_code'
        ]]

    return trans_df

# imports dataframe and organizes it in accord with the `std df` format
# exports as database file into `db_path`. If rescue is set to 1, save in home folder if output doesn't exist
# if rescue is set to 0, do not output database if output dir does not exit
def df2db(df, db_path, header = True, overwrite = False, std_col = True, rescue = True):

    df = df.set_index('internal_ome')
    df = df.sort_index()
    df = df.reset_index()

    if db_path != sys.stdout:
        db_path = formatPath( db_path )
    elif overwrite:
        number = 0
        while os.path.exists(db_path):
            number += 1
            db_path = os.path.normpath(db_path) + '_' + str(number)

    if std_col:
        df = df2std( df )

    while True:
        try:
            df.to_csv(db_path, sep = '\t', index = None, header = header)
            break
        except FileNotFoundError:
            if rescue:
                eprint('\nOutput directory does not exist. Attempting to save in home folder.')
                db_path = '~/' + os.path.basename(os.path.normpath(db_path))
                df.to_csv(db_path, sep = '\t', index = None)
                raise FileNotFoundError
                break
            else:
                eprint('\nOutput directory does not exist. Rescue not enabled.')
                raise FileNotFoundError
                break


# collect fastas from a database dataframe. NEEDS to be deprecated
def collect_fastas(db_df,assembly=False,ome_index='internal_ome'):

    if not assembly:
        print('\n\nCompiling fasta files from database, using proteome if available ...')
    elif assembly == 1:
        print('\n\nCompiling assembly fasta files from database ...')
    ome_fasta_dict = {}
    for index,row in db_df.iterrows():
        ome = row[ome_index]
        print(ome)
        ome_fasta_dict[ome] = {}
        if assembly:
            if type(row['assembly']) != float:
                ome_fasta_dict[ome]['fasta'] = os.environ['MYCOFNA'] + '/' + row['assembly']
                ome_fasta_dict[ome]['type'] = 'nt'
            else:
                print('Assembly file does not exist for',ome,'\nExit code 3')
                sys.exit(3)

        elif not assembly:
            if type(row['proteome']) != float:
                print('\tproteome')
                ome_fasta_dict[ome]['fasta'] = os.environ['MYCOFAA'] + '/' + row['proteome']
                ome_fasta_dict[ome]['type'] = 'prot'
            elif type(row['assembly']) != float:
                print('\tassembly')
                ome_fasta_dict[ome]['fasta'] = os.environ['MYCOFNA'] + '/' + row['assembly']
                ome_fasta_dict[ome]['type'] = 'nt'
            else:
                print('Assembly file does not exist for',ome,'\nExit code 3')
                sys.exit(3)

    return ome_fasta_dict


def hit2taxonomy( 
    taxid, 
    rank = 'kingdom', 
    lineage = 'fungi', 
    skip = False,
    email = None,
    api = None
    ):
    ''' 
    Takes a searchTerm string, queries NCBI via Entrez, obtains TaxIDs,
    if there are multiple hits return an error an query the TaxID for 
    taxonomy to find a fungal hit. If there is a fungal hit, create a
    tax_dict and return it.
    tax_dict = { taxonomicRank: rankName }
    '''

    tax_dict, count, sleep = {}, 0, False
    while True:
        try:
            tax_handle = Entrez.efetch( db = "Taxonomy", id = taxid, retmode = "xml" )
            break
        except HTTPError as goon:
            count += 1
            if count == 5 and skip:
                print( '\n5 failed HTTP queries. Is NCBI down?' )
                sys.exit( 100 )
            if 500 <= goon.code <= 599:
                time.sleep( 1 )
            sleep = True

    count = 0
    while True:
        try:
            records = Entrez.read( tax_handle, validate = False )
            break
        except Entrez.Parser.NotXMLError:
            time.sleep(1)
            count += 1
            if count == 5:
                eprint( '\nERROR: 5 failed taxonomy queries. ' + str(taxid) )
                eprint( tax_handle )
                for line in tax_handle:
                    print( line )
                if skip:
                    sys.exit( 1 )
                records = False
                break
            tax_handle = Entrez.efetch( db = "Taxonomy", id = taxid, retmode = "xml" )
            if 'latin-1' in str(tax_handle):
                count = 0
            
                eprint( '\tERROR: latin-1 encoding' )
                for line in tax_handle:
                    print( line )
                time.sleep( 300 )
                Entrez.email = email
                Entrez.api_key = api
 #               wrapper = io.TextIOWrapper(
  #                  io.BytesIO(),
   #                 encoding = 'utf-8',
    #                line_buffering = True,
     #               )
      #          for line in tax_handle:
       #             wrapper.write(line + '\n')
        #        tax_handle = wrapper
            sleep = True

    if records:
        if 'LineageEx' in records[0]:
            lineages = records[0][ 'LineageEx' ]
            for tax in lineages:
                if tax['Rank'].lower() == rank.lower():
                    tax_dict = { str(x['Rank']): str(x['ScientificName']) for x in lineages }
                    continue

    return tax_dict, sleep


# gather taxonomy by querying NCBI
# if `api_key` is set to `1`, it assumes the `Entrez.api` method has been called already
def gather_taxonomy(df, api_key = None, king='fungi', ome_index = 'internal_ome'):

    print('\nGathering taxonomy information')
    tax_dicts = []
    count = 0

    check = re.compile( r'[Tt][Aa][Xx][Ii][Dd]' )
    genus = re.compile( r'[Gg][Ee][Nn][Uu][Ss]' )
    columns = set( x.lower() for x in df.columns )
    searchTax = check.search( '#$'.join(list(df.columns)) )
    searchGen = genus.search( '#$'.join(list(df.columns)) )
    if searchTax:
        taxColumn = searchTax[0]
    else:
        check = re.compile( r'NCBI Taxon' )
        searchTax = check.search( '#$'.join(list(df.columns)) )
        if searchTax:
            taxColumn = searchTax[0]

# print each ome code
    for i,row in df.iterrows():
        ome = row[ome_index]

# if there is no api key, sleep for a second after 3 queries
# if there is an api key, sleep for a second after 10 queries
        if not api_key and count == 2:
            time.sleep(1)
            count = 0
        elif api_key and count == 9:
            time.sleep(1)
            count = 0

# gather taxonomy from information present in the genus column and query NCBI using Entrez
        ids = None
        if searchTax:
            ids = [row[taxColumn]]
        else:
            genus = row[searchGen[0]]
            search_term = str(genus) + ' [GENUS]'
            handle = Entrez.esearch(db='Taxonomy', term=search_term)
            ids = Entrez.read(handle)['IdList']
            count += 1
            if len(ids) == 0:
                print('\t' + ome)
                print('\t\tNo taxonomy information')
                continue

# for each taxID acquired, fetch the actual taxonomy information
        taxid = 0
        for tax in ids:
            if not api_key and count == 2:
                time.sleep(1)
                count = 0
            elif api_key and count == 9:
                    time.sleep(1)
                    count = 0
            count += 1
            handle = Entrez.efetch(db="Taxonomy", id=tax, remode="xml")
            records = Entrez.read(handle)
            lineages = records[0]['LineageEx']

# for each lineage, if it is a part of the kingdom, `king`, use that TaxID
# if there are multiple TaxIDs, use the first one found
            for lineage in lineages:
                if lineage['Rank'] == 'kingdom':
                    if lineage['ScientificName'].lower() == king:
                        taxid = tax
                        if len(ids) > 1:
                            print('\t' + ome)
                            print('\t\tTax ID(s): ' + str(ids))
                            print('\t\t\tUsing first "' + king + '" ID: ' + str(taxid))
                        break
            if taxid != 0:
                break
       
        if taxid == 0:
            print('\t' + ome)
            print('\t\tNo taxonomy information')
            continue

# for each taxonomic classification, add it to the taxonomy dictionary string
# append each taxonomy dictionary string to a list of dicts
        tax_dicts.append({})
        tax_dicts[-1][ome_index] = ome
        for tax in lineages:
             tax_dicts[-1][tax['Rank']] = tax['ScientificName']
    
        count += 1
        if count == 3 or count == 10:
            if not api_key:
                time.sleep( 1 )
                count = 0
            elif api_key:
                if count == 10:
                    time.sleep( 1 )
                    count = 0

# to read, use `read_tax` below
    return tax_dicts


# read taxonomy by conterting the string into a dictionary using `json.loads`
def read_tax(taxonomy_string):
    
    dict_string = taxonomy_string.replace("'",'"')
    tax_dict = json.loads(dict_string)

    return tax_dict


# converts db data from row or full db into a dictionary of dictionaries
# the first dictionary is indexed via `ome_index` and the second via db columns
# specific columns, like `taxonomy`, `gff`, etc. have specific functions to acquire
# their path. Can also report values from a list `empty_list` as a `None` type
def acquire( db_or_row, ome_index = 'internal_ome', empty_list = [''], datatype = 'row' ):

    key2env = {
        'assembly': 'MYCOFNA',
        'proteome': 'MYCOFAA',
        'gff3': 'MYCOGFF3'
        }
    if ome_index != 'internal_ome':
        alt_ome = 'internal_ome'
    else:
        alt_ome = 'genome_code'

    info_dict = {}
    if datatype == 'row':
        row = db_or_row
        ome = row[ome_index]

        info_dict[ome] = {}
        for key in row.keys():
            if not pd.isnull(row[key]) and not row[key] in empty_list:
                if key == alt_ome:
                    info_dict[ome]['alt_ome'] = row[alt_ome]
                elif key == 'taxonomy':
                    info_dict[ome][key] = read_tax( row[key] )
                elif key == 'gff3' or key == 'proteome' or key == 'assembly':
                    new = key2env[key]
                    info_dict[ome][new] = os.environ[new] + '/' + str(row[key])
                else:
                    info_dict[ome][key] = row[key]
            else:
                info_dict[ome][key] = None

    else:
        for i,row in db_or_row.iterrows():

            info_dict[ome] = {}
            for key in row.keys():
                if not pd.isnull( row[key] ) or not row[key] in empty_list:
                    if key == alt_ome:
                        info_dict[ome]['alt_ome'] == row[alt_ome]
                    elif key == 'taxonomy':
                        info_dict[ome][key] == read_tax( row[key] )
                    elif key == 'gff3' or key == 'proteome' or key == 'assembly':
                        new = key2env[key]
                        info_dict[ome][key] == os.environ[new] + '/' + row[key]
                    else:
                        info_dict[ome][key] == row[key]
                else:
                    info_dict[ome][key] == None

    return info_dict


# assimilate taxonomy dictionary strings and append the resulting taxonomy string dicts to an inputted database
# forbid a list of taxonomic classifications you are not interested in and return a new database
def assimilate_tax(db, tax_dicts, ome_index = 'internal_ome', forbid=['no rank', 'superkingdom', 'subkingdom', 'genus', 'species', 'species group', 'varietas', 'forma']):

    if 'taxonomy' not in db.keys():
        db['taxonomy'] = ''

    tax_df = pd.DataFrame(tax_dicts)
    for forbidden in forbid:
        if forbidden in tax_df.keys():
            del tax_df[forbidden]

    output = {}
    keys = tax_df.keys()
    for i,row in tax_df.iterrows():
        output[row[ome_index]] = {}
        for key in keys:
            if key != ome_index:
                if pd.isnull(row[key]):
                    output[row[ome_index]][key] = ''
                else:
                    output[row[ome_index]][key] = row[key]

    for index in range(len(db)):
        if db[ome_index][index] in output.keys():
            db.at[index, 'taxonomy'] = str(output[db[ome_index][index]])

    return db


# abstract taxonomy and return a database with just the taxonomy you are interested in
# NEED TO SIMPLIFY THE PUBLISHED STATEMENT TO SIMPLY DELETE ALL VALUES THAT ARENT 1 WHEN PUBLISHED = 1
def abstract_tax(db, taxonomy, classification, inverse = False ):

    if not classification or not taxonomy:
        eprint('\nERROR: Abstracting via taxonomy requires both taxonomy name and classification.')

    new_db = pd.DataFrame()

# `genus` is a special case because it has its own column, so if it is not a genus, we need to read the taxonomy
# once it is read, we can parse the taxonomy dictionary via the inputted `classification` (taxonomic classification)
    if classification != 'genus':
        for i,row in db.iterrows():
            row_taxonomy = read_tax(row['taxonomy'])
            if row_taxonomy[classification].lower() == taxonomy.lower() and not inverse:
                new_db = new_db.append(row)
            elif row_taxonomy[classification].lower() != taxonomy.lower() and inverse:
                new_db = new_db.append(row)

# if `classification` is `genus`, simply check if that fits the criteria in `taxonomy`
    elif classification == 'genus':
        for i,row in db.iterrows():
            if row['genus'].lower() == taxonomy.lower() and not inverse:
                new_db = new_db.append(row)
            elif row['genus'].lower() != taxonomy.lower() and inverse:
                new_db = new_db.append(row)

    return new_db


# abstracts all omes from an `ome_list` or the inverse and returns the abstracted database
def abstract_omes(db, ome_list, index = 'internal_ome', inverse = False):

    new_db = pd.DataFrame()
    ref_set = set()
    ome_set = set(ome_list)
    if not inverse:
        for ome in ome_set:
            ref_set.add(ome.lower())
    else:
        for i, row in db.iterrows():
            if row[index] not in ome_set:
                ref_set.add( row[index].lower() )

    db = db.set_index(index, drop = False)
    for i,row in db.iterrows():
        if row[index].lower() in ref_set:
            new_db = new_db.append(row)

    return new_db


# creates a fastadict for blast2fasta.py - NEED to remove and restructure that script
def ome_fastas2dict(ome_fasta_dict):

    ome_fasta_seq_dict = {}
    for ome in ome_fasta_dict:
        temp_fasta = fasta2dict(ome['fasta'])
        ome_fasta_seq_dict['ome'] = temp_fasta

    return ome_fasta_seq_dict

# this can lead to duplicates, ie Lenrap1_155 and Lenrap1
def column_pop(db_df,column,directory,ome_index='genome_code',file_type='fasta'):

    db_df = db_df.set_index(ome_index)
    files = collect_files(directory,file_type)
    for ome,row in db_df.iterrows():
        for file_ in files:
            file_ = os.path.basename(file_)
        if ome + '_' in temp_list:
            index = temp_list.index(ome + '_')
            db_df.at[ome, column] = files[index]

    db_df.reset_index(level = 0, inplace = True)

    return db_df


# imports a database reference (if provided) and a `df` and adds a `new_col` with new ome codes that do not overlap current ones
def gen_omes(df, reference = 0, new_col = 'internal_ome', tag = None):

    print('\nGenerating omes')
    newdf = df
    newdf[new_col] = ''
    tax_set = set()
    access = '-'
    if type(reference) == pd.core.frame.DataFrame:
        tax_set = set(row['internal_ome'])
        print('\t' + str(len(tax_set)) + ' omes from reference dataset')

    for i, row in df.iterrows():
        if pd.isnull(row['species']):
            newdf.at[i, 'species'] = 'sp.'
            row['species'] = 'sp.'
        name = row['genus'][:3].lower() + row['species'][:3].lower()
        name = re.sub(r'\(|\)|\[|\]|\$|\#|\@| |\||\+|\=|\%|\^|\&|\*|\'|\"|\!|\~|\`|\,|\<|\>|\?|\;|\:|\\|\{|\}', '', name)
        name.replace('(', '')
        name.replace(')', '')
        while len(name) < 6:
            name += '.'
        numbers = [
            int(re.search(r'.*?(\d+)', x)[1]) for x in list(tax_set) if x.startswith(name)
        ]
        numbers.sort( reverse = True )
        number = numbers[-1] + 1
        ome = name + str(number)
        tax_set.add(ome)
        newdf.at[i, new_col] = ome
    
    return newdf

# queries database either by ome or column or ome and column. print the query, and return an exit_code (0 if fine, 1 if query1 empty)
def query(db, query_type, query1, query2 = None, ome_index='internal_ome'):

    exit_code = 0
    if query_type in ['column', 'c', 'ome', 'o']:

        db = db.set_index(ome_index)

        if query_type in ['ome', 'o']:
            if query1 in db.index.values:
                if query2 != None and query2 in db.columns.values:
                    print(query1 + '\t' + str(query2) + '\n' + str(db[query2][query1]))
                elif query2 != None:
                    print('Query ' + str(query2) + ' does not exist.')
                else:
                    print(query1 + '\n' + str(db.loc[query1]))
            else:
                print('Query ' + str(query1) + ' does not exist.')

        elif query_type in ['column', 'c']:
            if query1 in db.columns.values:
                if query2 != None and query2 in db.index.values:
                    print(query2 + '\t' + str(query1) +  '\n' + str(db[query1][query2]))
                elif query2 != None:
                    print('Query ' + str(query2) + ' does not exist.')
                else:
                    print(str(db[query1]))
            else:
                print('Query ' + str(query1) + ' does not exist.')

    if query1 == '':
        exit_code == 1

    return exit_code


# edits an entry by querying the column and ome code and editing the entry with `edit`
# returns the database and an exit code
def db_edit(query1, query2, edit, db_df, query_type = '', ome_index='internal_ome',new_col=0):

    exit_code = 0
    if query_type in ['column', 'c']:
        ome = query2
        column = query1
    else:
        ome = query1
        column = query2

    if column in db_df or new_col==1:
        db_df = db_df.set_index(ome_index)
        db_df.at[ome, column] = edit

    db_df.reset_index(level = 0, inplace = True)

    if ome == '' and column == '':
        exit_code = 1

    return db_df, exit_code


def add_jgi_data(db_df,jgi_df,db_column,jgi_column,binary='0'):

    jgi_df_ome = jgi_df.set_index('genome_code')
    db_df.set_index('genome_code',inplace=True)

    jgi_col_dict = {}
    for ome in jgi_df['genome_code']:
        jgi_col_dict[ome] = jgi_df_ome[jgi_column][ome]

    if binary == 1:
        for ome in jgi_col_dict:
            if ome in db_df.index:
                if pd.isnull(jgi_col_dict[ome]):    
                    db_df.at[ome, db_column] = 0
                else:
                    db_df[ome, db_column] = 1
            else:
                print(ome + '\tNot in database')

    else:
        for ome in jgi_col_dict:
            if ome in db_df.index:
                db_df[ome, db_column] = jgi_col_dict[ome]
            else:
                print(ome + '\tNot in database')

    db_df.reset_index(inplace=True)

    return db_df


def report_omes( db_df, ome_index = 'internal_ome' ):

    ome_list = db_df[ome_index].tolist()
    print('\n\n' + str(ome_list))



