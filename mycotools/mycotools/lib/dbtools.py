#! /usr/bin/env python3

import pandas as pd, pandas
import numpy as np
import sys, os, time, mycotools.lib.fastatools, json, re, datetime, hashlib, getpass, base64
from mycotools.lib.kontools import collect_files, eprint, formatPath
from Bio import Entrez
from urllib.error import HTTPError
from io import StringIO

def getLogin( ncbi, jgi ):

    ncbi_email, ncbi_api, jgi_email, jgi_pwd = None, None, None, None
    print(flush = True)
    if ncbi:
        ncbi_email = input( 'NCBI email: ' )
        ncbi_api = getpass.getpass( prompt = 'NCBI api key (blank if none): ' )
    if jgi:
        jgi_email = input( 'JGI email: ' )
        jgi_pwd = getpass.getpass( prompt = 'JGI password (required): ' )
    print(flush = True)

    return ncbi_email, ncbi_api, jgi_email, jgi_pwd


def loginCheck( info_path = '~/.mycodb', ncbi = True, jgi = True ):

    salt = b'D9\x82\xbfSibW(\xb1q\xeb\xd1\x84\x118'
    #NEED to make this store a password
    if os.path.isfile( formatPath(info_path) ):
        from cryptography.fernet import Fernet
        from cryptography.hazmat.primitives.kdf.pbkdf2 import PBKDF2HMAC
        from cryptography.hazmat.backends import default_backend
        kdf = PBKDF2HMAC(
            algorithm=hashlib.sha256(),
            length=32,
            salt=salt,
            iterations=100000,
            backend=default_backend()
        )
        print(flush = True)
        if sys.stdin.isatty():
            hash_pwd = getpass.getpass( prompt = 'MycoDB login password (stdin allowed): ' )
        else:
            hash_pwd = sys.stdin.readline().rstrip()
        key = base64.urlsafe_b64encode(kdf.derive(hash_pwd.encode('utf-8')))
        fernet = Fernet( key )
#        with open(formatPath(info_path) + '/.key', 'rb') as raw_key:
 #           fernet = Fernet(raw_key)
        with open(formatPath(info_path), 'rb') as raw_file:
            data = raw_file.read()
        decrypted = fernet.decrypt(data)
        data = decrypted.decode('UTF-8').split('\n')
        if len(data) != 4:
            eprint('BAD PASSWORD FILE. Delete ~/.mycodb to reset.', flush = True)
            sys.exit(8)
        ncbi_email = data[0].rstrip()
        ncbi_api = data[1].rstrip()
        jgi_email = data[2].rstrip()
        jgi_pwd = data[3].rstrip()
    else:
        ncbi_email, ncbi_api, jgi_email, jgi_pwd = getLogin(ncbi, jgi)
        if ncbi and jgi:
            hash_check = input( 'Would you like to encrypt your login ' + \
                'information to ' + info_path + ' [Y/n]: ' )
        else:
            hash_check = 'n'
        if hash_check not in {'n', 'N'}:
            from cryptography.fernet import Fernet
            from cryptography.hazmat.primitives.kdf.pbkdf2 import PBKDF2HMAC
            from cryptography.hazmat.backends import default_backend
            kdf = PBKDF2HMAC(
                algorithm=hashlib.sha256(),
                length=32,
                salt=salt,
                iterations=100000,
                backend=default_backend()
            )
            print(flush = True)
            hash_pwd, hash_check = True, False
            while hash_pwd != hash_check:
                if hash_pwd != True:
                    eprint('ERROR: passwords do not match', flush = True)
                hash_pwd = getpass.getpass( prompt = 'New MycoDB login password: ' )
                hash_check = getpass.getpass( prompt = 'Confirm password: ' )

            key = base64.urlsafe_b64encode(kdf.derive(hash_pwd.encode('utf-8')))
            fernet = Fernet( key )
            out_data = ncbi_email + '\t' + ncbi_api + '\t' + jgi_email + '\t' + jgi_pwd
            encrypt_data = fernet.encrypt(out_data)
            with open( formatPath(info_path), 'wb' ) as out:
                out.write(encrypt_data)
#
    return ncbi_email, ncbi_api, jgi_email, jgi_pwd


def genConfig( 
    branch = 'stable', forbidden = '$MYCODB/log/forbidden.tsv',
    repo = "https://gitlab.com/xonq/mycodb", 
    rogue = False, nonpublished = False
    ):

    config = { 
        'forbidden': forbidden,
        'repository': repo, 
        'branch': branch,
        'nonpublished': nonpublished,
        'rogue': rogue,
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
            eprint('\nERROR: master db not found in ' + path + '. Have you initialized MycoDB?', flush = True)
        master_path = formatPath( '$' + path + '/' + master + '.db' )
    except KeyError:
        eprint('\nERROR: ' + path + ' not in path' , flush = True)
        master_path = None

    return master_path
 
# imports database, converts into df
# returns database dataframe
def db2df(data, stdin = False):

    headers = [ 'internal_ome', 'genus', 'species', 'strain', 'version', 
        'biosample', 'assembly', 'proteome', 'gff3',
        'taxonomy', 'ecology', 'eco_conf', 'source', 'published',
        'genome_code', 'acquisition_date' ]
    if not stdin:
        data = formatPath( data )
        db_df = pd.read_csv( data, sep='\t' )
        if 'internal_ome' not in set( db_df.columns ) and 'genome_code' not in set( db_df.columns ):
            db_df = pd.read_csv( data, sep = '\t', names = headers )
    else:
        db_df = pd.read_csv( StringIO(data), sep='\t' )
        if 'internal_ome' not in set( db_df.columns ) and 'genome_code' not in set( db_df.columns ):
            db_df = pd.read_csv( StringIO(data), sep = '\t', names = headers )

    return db_df


def df2std( df ):
    '''
    Standardized organization of database columns
    '''
    # temporary check, will remove when all scripts account for new biosample column
    if 'biosample' in df.columns:
        trans_df = df[[
        'internal_ome', 'genus', 'species',
        'strain', 'version', 'biosample', 'assembly', 
        'proteome', 'gff3', 'taxonomy',
        'ecology', 'eco_conf', 'source', 'published', 
        'genome_code', 'acquisition_date'
    ]]
    else:
        trans_df = df[[
        'internal_ome', 'genus', 'species',
        'strain', 'assembly', 
        'proteome', 'gff', 'gff3', 'taxonomy',
        'ecology', 'eco_conf', 'source', 'published', 
        'genome_code'
        ]]

    return trans_df


# imports dataframe and organizes it in accord with the `std df` format
# exports as database file into `db_path`. If rescue is set to 1, save in home folder if output doesn't exist
# if rescue is set to 0, do not output database if output dir does not exit
def df2db(df, db_path, header = False, overwrite = False, std_col = True, rescue = True):

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
                eprint('\nOutput directory does not exist. Attempting to save in home folder.', flush = True)
                db_path = '~/' + os.path.basename(os.path.normpath(db_path))
                df.to_csv(db_path, sep = '\t', index = None)
                raise FileNotFoundError
                break
            else:
                eprint('\nOutput directory does not exist. Rescue not enabled.', flush = True)
                raise FileNotFoundError
                break


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
                print( '\n5 failed HTTP queries. Is NCBI down?' , flush = True)
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
                eprint( '\nERROR: 5 failed taxonomy queries. ' + str(taxid) , flush = True)
                eprint( tax_handle , flush = True)
                for line in tax_handle:
                    print( line , flush = True)
                if skip:
                    sys.exit( 1 )
                records = False
                break
            tax_handle = Entrez.efetch( db = "Taxonomy", id = taxid, retmode = "xml" )
            if 'latin-1' in str(tax_handle):
                count = 0
            
                eprint( '\tERROR: latin-1 encoding' , flush = True)
                for line in tax_handle:
                    print( line , flush = True)
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

    print('\nGathering taxonomy information', flush = True)
    df['taxonomy'] = df['taxonomy'].replace( np.nan, None )
    count = 0
    tax_dicts = {x['genus']: read_tax(x['taxonomy']) for i,x in df.iterrows()}

    for genus in tax_dicts:
        if any(tax_dicts[genus][x] for x in tax_dicts[genus]):
            continue
        print('\t' + genus, flush = True)
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
        search_term = str(genus) + ' [GENUS]'
        handle = Entrez.esearch(db='Taxonomy', term=search_term)
        ids = Entrez.read(handle)['IdList']
        count += 1
        if len(ids) == 0:
            print('\t\tNo taxonomy information', flush = True)
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
        for attempt in range(3):
            try:
                handle = Entrez.efetch(db="Taxonomy", id=tax, remode="xml")
                records = Entrez.read(handle)
                lineages = records[0]['LineageEx']
                break
            except IndexError:
                lineages = []

# for each lineage, if it is a part of the kingdom, `king`, use that TaxID
# if there are multiple TaxIDs, use the first one found
        if king:
            for lineage in lineages:
                if lineage['Rank'] == 'kingdom':
                    if lineage['ScientificName'].lower() == king:
                        taxid = tax
                        if len(ids) > 1:
                            print('\t\tTax ID(s): ' + str(ids), flush = True)
                        break
        else:
            for lineage in lineages:
                taxid = tax 
       
        if taxid == 0:
            print('\t\tNo taxonomy information', flush = True)
            continue

# for each taxonomic classification, add it to the taxonomy dictionary string
# append each taxonomy dictionary string to a list of dicts
        tax_dicts[genus] = {}
        for tax in lineages:
             tax_dicts[genus][tax['Rank'].lower()] = tax['ScientificName']
    
        count += 1
        if count == 3 or count == 10:
            if not api_key:
                time.sleep( 1 )
                count = 0
            elif api_key:
                if count == 10:
                    time.sleep( 1 )
                    count = 0

    return tax_dicts


# read taxonomy by conterting the string into a dictionary using `json.loads`
def read_tax(taxonomy_string):
   
    tax_strs = ['kingdom', 'phylum', 'subphylum', 'class', 'order', 'family', 'subfamily']
    if taxonomy_string: 
        dict_string = taxonomy_string.replace("'",'"')
        tax_dict = json.loads(dict_string)
        tax_dict = {**tax_dict, **{x: '' for x in tax_strs if x not in tax_dict}}
        return tax_dict
    else:
        return {}


# assimilate taxonomy dictionary strings and append the resulting taxonomy string dicts to an inputted database
# forbid a list of taxonomic classifications you are not interested in and return a new database
def assimilate_tax(db, tax_dicts, ome_index = 'internal_ome', forbid={'no rank', 'superkingdom', 'subkingdom', 'genus', 'species', 'species group', 'varietas', 'forma'}):

    genera = set(db['genus'])
    tax_dicts = {x: tax_dicts[x] for x in tax_dicts if tax_dicts[x]}
    for genus in tax_dicts:
        tax_dicts[genus] = {x: tax_dicts[genus][x] for x in tax_dicts[genus] if x not in forbid}
    missing = list(genera.difference(set(tax_dicts.keys())))
    
    for miss in missing:
        tax_dicts[miss] = {}
    db['taxonomy'] = db.apply(
        lambda x: tax_dicts[x['genus']], axis = 1
        )

    return db


# extract taxonomy and return a database with just the taxonomy you are interested in
def extract_tax(db, taxonomy, classification, inverse = False ):

    if type(taxonomy) == str:
        taxonomy = [taxonomy]

    taxonomy = set(x.lower() for x in taxonomy)
    new_db = pd.DataFrame()

# `genus` is a special case because it has its own column, so if it is not a genus, we need to read the taxonomy
# once it is read, we can parse the taxonomy dictionary via the inputted `classification` (taxonomic classification)
    if classification != 'genus':
        for i,row in db.iterrows():
            row_taxonomy = read_tax(row['taxonomy'])
            if classification in row_taxonomy:
                if row_taxonomy[classification].lower() in taxonomy and not inverse:
                    new_db = new_db.append(row)
                elif inverse and row_taxonomy[classification].lower() not in taxonomy:
                    new_db = new_db.append(row)
            elif inverse:
               new_db = new_db.append(row)

# if `classification` is `genus`, simply check if that fits the criteria in `taxonomy`
    elif classification == 'genus':
        taxonomy = set(x[0].upper() + x[1:].lower() for x in taxonomy)
        if not inverse:
            new_db = db[db['genus'].isin(taxonomy)]
        else:
            new_db = db[~db['genus'].isin(taxonomy)]

    return new_db


# extracts all omes from an `ome_list` or the inverse and returns the extracted database
def extract_omes(db, ome_list, index = 'internal_ome', inverse = False):

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


# imports a database reference (if provided) and a `df` and adds a `new_col` with new ome codes that do not overlap current ones
def gen_omes(df, reference = 0, new_col = 'internal_ome', tag = None):

#    print('\nGenerating omes', flush = True)
    newdf = df
    if not new_col in set(df.keys()):
        newdf[new_col] = ''
    tax_set = set()
    access = '-'
    if isinstance(reference, pd.core.frame.DataFrame):
        tax_set = set(reference['internal_ome'])
#        print('\t' + str(len(tax_set)) + ' omes from reference dataset', flush = True)

    for i, row in df.iterrows():
        if not pd.isnull(row['internal_ome']) and row['internal_ome']:
            continue
        if pd.isnull(row['species']) or not row['species']:
            newdf.at[i, 'species'] = 'sp.'
            row['species'] = 'sp.'
        name = row['genus'][:3].lower() + row['species'][:3].lower()
        name = re.sub(r'\(|\)|\[|\]|\$|\#|\@| |\||\+|\=|\%|\^|\&|\*|\'|\"|\!|\~|\`|\,|\<|\>|\?|\;|\:|\\|\{|\}', '', name)
        name.replace('(', '')
        name.replace(')', '')
        while len(name) < 6:
            name += '.'
        numbers = [
            int(re.search(r'.*?(\d+$)', x)[1]) for x in list(tax_set) if x.startswith(name)
            ]
        numbers.sort(reverse = True)
        if numbers:
            number = numbers[0] + 1
        else:
            number = 1
        ome = name + str(number)
        tax_set.add(ome)
        newdf.at[i, new_col] = ome
    
    return newdf
