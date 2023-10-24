#! /usr/bin/env python3

# NEED to make login check more intuitive and easier for NCBI only

import os
import re
import sys
import copy
import json
import time
import base64
import urllib
import getpass
import hashlib
import datetime
from Bio import Entrez
from io import StringIO
from collections import defaultdict
from mycotools.lib.kontools import collect_files, eprint, format_path, \
    read_json, write_json


class mtdb(dict):
    '''
    MycotoolsDB (mtdb) class. Designed to support high throughput "pandas-like"
    interface without the overhead of pandas.
    '''
    columns = [
       'ome', 'genus', 'species', 'strain', 'taxonomy',
       'version', 'source', 'biosample', 'assembly_acc',
       'acquisition_date', 'published', 'fna', 'faa', 'gff3'
       ]

    # NEED a detect index feature for adding dicts in with alternative indices
    def __init__(self, db = None, index = None, add_paths = True):
        self.columns = [
            'ome', 'genus', 'species', 'strain', 'taxonomy',
            'version', 'source', 'biosample', 'assembly_acc',
            'acquisition_date', 'published', 'fna', 'faa', 'gff3'
            ]
        if not db:
            if not index:
                super().__init__({x: [] for x in self.columns})
            else:
                super().__init__()
        elif isinstance(db, dict):
            if not index:
                if all(isinstance(db[x], list) for x in self.columns):
                    super().__init__(db)
                    self.index = None
            else:
                try:
                    if set(db[list(db.keys())[0]].keys()).union({index}) \
                        == set(self.columns):
                        super().__init__(db)
                except AttributeError:
                    if set(db[list(db.keys())[0]][0]).union({index}) \
                        == set(self.columns):
                        super().__init__(db)
        else:
            super().__init__(mtdb.db2df(self, db, add_paths = add_paths))
        self.index = index

    def mtdb2pd(self):
        import pandas as pd
        copy_mtdb = copy.deepcopy(self)
        copy_mtdb = copy_mtdb.reset_index()
        return pd.DataFrame(copy_mtdb) # assume pd is imported

    def pd2mtdb(df): # legacy integration
        df = df.fillna('')
        for i, row in df.iterrows():
            if not row['gff3']:
                row['gff3'] = os.environ['MYCOGFF3'] + row['ome'] + '.gff3'
            if not row['faa']:
                row['faa'] = os.environ['MYCOFAA'] + row['ome'] + '.faa'
            if not row['fna']:
                row['fna'] = os.environ['MYCOFNA'] + row['ome'] + '.fna'
        db = mtdb({
            x: list(df[x]) for x in mtdb.columns
            })
        return db
        
    def db2df(self, db_path, add_paths = True):
        df = defaultdict(list)
        if os.stat(db_path).st_size == 0:
            return {x: [] for x in mtdb.columns}
        with open(format_path(db_path), 'r') as raw:
            data = [x.rstrip().split('\t') for x in raw if not x.startswith('#')]
        columns = self.columns
        for entry in data:
            [df[c].append('') for c in columns] # add a blank entry to each
            # column
            for i, d in enumerate(entry):
                df[columns[i]][-1] = d
            try:
                df['taxonomy'][-1] = read_tax(
                    df['taxonomy'][-1]
                    )
            except json.decoder.JSONDecodeError:
                print(df['taxonomy'][-1])
                sys.exit()
            df['taxonomy'][-1]['genus'] = df['genus'][-1]
            df['taxonomy'][-1]['species'] = \
                df['genus'][-1] + ' ' + df['species'][-1]
            df['taxonomy'][-1]['strain'] = \
                df['strain'][-1]
        
        if not add_paths:
            return df
        try:
            for i, ome in enumerate(df['ome']): 
                if not df['fna'][i]:
                    if {'MYCOFNA', 'MYCOFAA', 'MYCOGFF3'}.difference(set(os.environ.keys())):
                        raise FileNotFoundError('You are not connected to a primary MTDB. ' \
                                              + 'Standalone databases need absolute paths')
                    df['fna'][i] = os.environ['MYCOFNA'] + ome + '.fna'
                    df['faa'][i] = os.environ['MYCOFAA'] + ome + '.faa'
                    df['gff3'][i] = os.environ['MYCOGFF3'] + ome + '.gff3'
                elif df['fna'][i] == ome + '.fna':
                    if {'MYCOFNA', 'MYCOFAA', 'MYCOGFF3'}.difference(set(os.environ.keys())):
                        raise FileNotFoundError('You are not connected to a primary MTDB. ' \
                                              + 'Standalone databases need absolute paths')
                    df['fna'][i] = os.environ['MYCOFNA'] + ome + '.fna'
                    df['faa'][i] = os.environ['MYCOFAA'] + ome + '.faa'
                    df['gff3'][i] = os.environ['MYCOGFF3'] + ome + '.gff3'
        except KeyError:
            eprint('ERROR: MycotoolsDB not in path, cannot delineate biofile paths', flush = True)
    
        return df

    def df2db(self, db_path = None, headers = False):
        df = copy.copy(self)
        df = df.reset_index()
        output = mtdb(
            {k: v for k,v in sorted(self.set_index('ome').items(), key = lambda x: x[0])},
            index = 'ome'
            ) 
        # does this work if its not an inplace change
        paths = {
            'faa': [os.environ['MYCOFAA'], '.faa'],
            'fna': [os.environ['MYCOFNA'], '.fna'],
            'gff3': [os.environ['MYCOGFF3'], '.gff3']
            }
        if db_path:
            with open(db_path, 'w') as out:
                if headers:
                    out.write('#' + '\t'.join(self.columns)+ '\n')
                for ome in output:
                    for file_type in ['fna', 'faa', 'gff3']:
                        output[ome][file_type] = output[ome][file_type].replace(
                            paths[file_type][0] + ome + paths[file_type][1], ''
                            ) # abbreviate when possible
                    if output[ome]['taxonomy']:
                        output[ome]['taxonomy'] = json.dumps(output[ome]['taxonomy'])
                    else:
                        output[ome]['taxonomy'] = '{}'
                    if not output[ome]['published']:
                        output[ome]['published'] = ''
                    out.write(
                        ome + '\t' + \
                        '\t'.join([str(output[ome][x]) for x in output[ome]]) + '\n'
                        )
        else:
            if headers:
                print(
                    '#' + '\t'.join(self.columns), flush = True
                    )

            for ome in output:
                for file_type in ['fna', 'faa', 'gff3']:
                    output[ome][file_type] = output[ome][file_type].replace(
                        paths[file_type][0] + ome + paths[file_type][1], ''
                        ) # abbreviate when possible
                for rank in ['species', 'genus', 'strain']:
                    try:
                        del output[ome]['taxonomy'][rank]
                    except KeyError:
                        pass
                output[ome]['taxonomy'] = json.dumps(output[ome]['taxonomy'])
                print(ome + '\t' + \
                    '\t'.join([str(output[ome][x]) for x in output[ome]]), flush = True
                    )

    def set_index(self, column = 'ome', inplace = False):
        data, retry, error, df, columns = {}, bool(column), False, copy.copy(self), copy.copy(self.columns)
        if not column:
            return df.reset_index()
        while retry:
            oldCol = set()
            try:
                columns.pop(columns.index(column))
            except (ValueError, IndexError) as e:
                df = df.reset_index() #will this actually reset the index
                columns.pop(columns.index(column))
            try:
                df[column]
            except KeyError:
                df = df.reset_index()
            if len(df[column]) == 0:
                return mtdb({}, index = column)
            if column in \
                {'assembly_acc', 'ome', 'fna', 'gff3', 'faa'}: # if it's unique
                for i, v in enumerate(df[column]):
                    data[v] = {}
                    for head in columns:
                        data[v][head] = df[head][i]
            else:
                for i, v in enumerate(df[column]):
                    if v not in data:
                        data[v] = []
                    data[v].append({})
                    for head in columns:
                        try:
                            data[v][-1][head] = df[head][i]
                        except KeyError: # if an index exists
                            if error:
                                eprint('\nERROR: invalid column', flush = True)
                                return self
                            df = df.reset_index()
                            error = True
                            break

            retry = False
        return mtdb(data, column)

    def reset_index(self):
        df = copy.copy(self)
        if df.index:
            data = {x: [] for x in mtdb().columns}
            data[df.index] = list(df.keys())
            if df.index in {'assembly_acc', 'ome', 'biosample', 'fna', 'gff3', 'faa'}:
                for key in df:
                    for otherCol in df[key]:
                        data[otherCol].append(df[key][otherCol])
            else:
                for key in df:
                    for v in df[key]:
                        for otherCol in v:
                            data[otherCol].append(v[otherCol])
            return mtdb(data, index = None)
        else:
            return df

    def append(self, info = {}):
#        if any(x not in set(self.columns) for x in info):
 #           raise KeyError('Invalid keys: ' + str(set(info.keys()).difference(set(self.columns))))
        index = self.index
        df = copy.copy(self)
        df = df.reset_index()
        info = {**info, **{k: None for k in set(self.columns).difference(set(info.keys()))}}
        for key in self.columns:
            df[key].append(info[key])
        return df.set_index(index)
                
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

def loginCheck(info_path = '~/.mycotools/mtdb_key', ncbi = True, jgi = True):
    salt = b'D9\x82\xbfSibW(\xb1q\xeb\xd1\x84\x118'
    #NEED to make this store a password
    if os.path.isfile( format_path(info_path) ):
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
            hash_pwd = getpass.getpass( prompt = 'MTDB login password (stdin allowed): ' )
        else:
            hash_pwd = sys.stdin.readline().rstrip()
        key = base64.urlsafe_b64encode(kdf.derive(hash_pwd.encode('utf-8')))
        fernet = Fernet( key )
#        with open(format_path(info_path) + '/.key', 'rb') as raw_key:
 #           fernet = Fernet(raw_key)
        with open(format_path(info_path), 'rb') as raw_file:
            data = raw_file.read()
        decrypted = fernet.decrypt(data)
        data = decrypted.decode('UTF-8').split('\n')
        if len(data) != 4:
            eprint('BAD PASSWORD FILE. Delete ~/.mycotools/mtdb_key to reset.', flush = True)
            sys.exit(8)
        ncbi_email = data[0].rstrip()
        ncbi_api = data[1].rstrip()
        jgi_email = data[2].rstrip()
        jgi_pwd = data[3].rstrip()
    else:
        ncbi_email, ncbi_api, jgi_email, jgi_pwd = getLogin(ncbi, jgi)
        # CURRENTLY THE REST DOESNT WORK, SO SKIP FOR NOW
        return ncbi_email, ncbi_api, jgi_email, jgi_pwd
        if ncbi and jgi:
            hash_check = input('Would you like to encrypt your login ' + \
                'information to ' + info_path + ' [Y/n]: ')
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
            with open( format_path(info_path), 'wb' ) as out:
                out.write(encrypt_data)
#
    return ncbi_email, ncbi_api, jgi_email, jgi_pwd


# opens a `log` file path to read, searches for the `ome` code followed by a whitespace character, and edits the line with `edit` 
def log_editor( log, ome, edit ):

    try:
        with open( log, 'r' ) as raw:
            data = raw.read()
    except FileNotFoundError:
        data = ''

    if re.search( ome + r'\s', data ):
        new_data = re.sub( ome + r'\s.*', edit, data )
    else:
        data = data.rstrip()
        data += '\n'
        new_data = data + edit

    with open( log, 'w' ) as towrite:
        towrite.write(new_data)

def readLog( log, columns = '', sep = '\t' ):

    log_dict = {}
    with open( log, 'r' ) as raw:
        data = raw.read()
    dataLines = [ x.split( sep ) for x in data.split( '\n' ) ]

    if type( columns ) is str:
        if dataLines[0][0].startswith( '#' ):
            dataLines[0][0] = re.sub( r'^\#+', '', dataLines[0][0] )
        columns = [ head for head in dataLines[0] ]
        del dataLines[0]
    elif type( columns ) is int:
        x = 0
        columns = [ x + 1 for y in dataLines[0] ]

    for line in dataLines:
        log_dict[ line[0] ] = { columns[ x ]: line[ x ] for x in range( 1, len( line ) ) }

    return log_dict

def primaryDB(path = '$MYCODB'):
    """Acquire the path of the master database by searching $MYCODB for a file
    with a basename that starts with a date string %Y%m%d."""

    path = path.replace('$', '')
    try:
        full_path = os.environ[path]
    except KeyError: # $MYCODB not initialized
        return None

    files = collect_files(full_path, 'mtdb')
    basenames = [os.path.basename(x) for x in files]
    dates = [x.replace('.mtdb','') \
            for x in basenames \
            if re.search(r'^\d+\.mtdb$', x)]
    master = '19991231' # arbitrary master for sorting
    for date in dates:
        dtDate = datetime.datetime.strptime(date, "%Y%m%d")
        dtMast = datetime.datetime.strptime(master, "%Y%m%d")
        if dtDate > dtMast:
            master = date
    if master == '19991231': # if it is the arbitrary start
        eprint('\nERROR: master db not found in ' + path \
              + '. Have you initialized MycoDB?', flush = True)
        return None
    master_path = format_path( '$' + path + '/' + master + '.mtdb')

    return master_path
 
# imports database, converts into df
# returns database dataframe
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

    return db_df

def df2std( df ):
    '''
    Standardized organization of database columns. DEPRECATED pandas MTDB
    '''
    # temporary check, will remove when all scripts account for new biosample column
    trans_df = df[
        mtdb.columns
        ]
    return trans_df

# imports dataframe and organizes it in accord with the `std df` format
# exports as database file into `db_path`. If rescue is set to 1, save in home folder if output doesn't exist
# if rescue is set to 0, do not output database if output dir does not exit
def df2db(df, db_path, header = False, overwrite = False, std_col = True, rescue = True):
    import pandas as pd, pandas

    df = df.set_index('ome')
    df = df.sort_index()
    df = df.reset_index()

    if db_path != sys.stdout:
        db_path = format_path( db_path )
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
        except urllib.error.HTTPError as goon:
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
def gather_taxonomy(df, api_key = None, king='fungi', ome_index = 'ome',
                    rank = 'kingdom'):

    if isinstance(df, mtdb) or isinstance(df['taxonomy'], list):
        df = df.set_index('ome')
        tax_dicts = {v['genus']: read_tax(v['taxonomy']) for k, v in df.items()}
    else:
        df['taxonomy'] = df['taxonomy'].fillna({})
        tax_dicts = {x['genus']: read_tax(x['taxonomy']) for i,x in df.iterrows()}
            

    count = 0

    for genus in tax_dicts:
        if any(tax_dicts[genus][x] for x in tax_dicts[genus] \
            if x not in {'genus', 'species', 'strain'}):
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
        fails = 0
        while True:
            try:
                search_term = str(genus) + ' [GENUS]'
                handle = Entrez.esearch(db='Taxonomy', term=search_term)
                ids = Entrez.read(handle)['IdList']
                break
            except RuntimeError:
                fails += 1
                if fails > 3:
                    break
                time.sleep(1)
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
            except:
                time.sleep(1)
                handle = Entrez.efetch(db="Taxonomy", id=tax, remode = "xml")
                records = Entrez.read(handle)
                lineages = records[0]['LineageEx']
                break

# for each lineage, if it is a part of the kingdom, `king`, use that TaxID
# if there are multiple TaxIDs, use the first one found
        if king:
            for lineage in lineages:
                if lineage['Rank'].lower() == rank.lower():
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
        if isinstance(taxonomy_string, str):
            dict_string = taxonomy_string.replace("'",'"')
            try:
                tax_dict = json.loads(dict_string)
            except TypeError:
                tax_dict = {}
        else:
            tax_dict = taxonomy_string
        try:
            tax_dict = {**tax_dict, **{x: '' for x in tax_strs if x not in tax_dict}}
        except TypeError: # inappropriate tax_dict in the column
            tax_dict = {x: '' for x in tax_strs}
        return tax_dict
    else:
        return {}

# assimilate taxonomy dictionary strings and append the resulting taxonomy string dicts to an inputted database
# forbid a list of taxonomic classifications you are not interested in and return a new database
def assimilate_tax(db, tax_dicts, ome_index = 'ome', 
                   forbid={'no rank', 'superkingdom', 'subkingdom', 'genus', 
                           'species', 'species group', 'varietas', 'forma'}):

    genera = set(db['genus'])
    tax_dicts = {x: tax_dicts[x] for x in tax_dicts if tax_dicts[x]}
    for genus in tax_dicts:
        tax_dicts[genus] = {x: tax_dicts[genus][x] for x in tax_dicts[genus] if x not in forbid}
    missing = list(genera.difference(set(tax_dicts.keys())))
    
    for miss in missing:
        tax_dicts[miss] = {}
    if isinstance(db, mtdb):
        for i, genus in enumerate(db['genus']):
            db['taxonomy'][i] = tax_dicts[genus]
        return mtdb(db)
    else:
        db['taxonomy'] = db.apply(
            lambda x: tax_dicts[x['genus']], axis = 1
            )

    return db, tax_dicts

def mtdb_disconnect(config, mtdb_config_file = format_path('~/.mycotools/config.json')):
    config['active'] = False
    write_json(config, mtdb_config_file)

def mtdb_connect(config, dbtype, 
                 mtdb_config_file= format_path('~/.mycotools/config.json')):
    config['active'] = dbtype
    write_json(config, mtdb_config_file)
    for var, env in config[config['active']].items():
        os.environ[var] = env

def mtdb_initialize(mycodb_loc, 
                    mtdb_config_file= format_path('~/.mycotools/config.json'),
                    init = False):
    mtdb_config = read_json(mycodb_loc + 'config/mtdb.json')
    dbtype = mtdb_config['branch']
    eprint('Establishing ' + dbtype + ' connection', flush = True)
    config = {}

    if not os.path.isdir(mycodb_loc + 'mtdb') and not init:
        raise FileNotFoundError('invalid MycotoolsDB path')
    dPath = mycodb_loc + 'data/'
    config[dbtype] = {
        'MYCODB': mycodb_loc + 'mtdb/',
        'MYCOFNA': dPath + 'fna/',
        'MYCOFAA': dPath + 'faa/',
        'MYCOGFF3': dPath + 'gff3/'
        }
    mtdb_connect(config, dbtype, mtdb_config_file)


interface = format_path('~/.mycotools/config.json')
if os.path.isfile(interface):
    envs_info = read_json(interface)
    if envs_info['active']:
        for var, env in envs_info[envs_info['active']].items():
            os.environ[var] = env

#if not primaryDB():
#    eprint('WARNING: Primary MycotoolsDB not connected; setup using `mtdb u/-i/-p/-f`', flush = True)

