#! /usr/bin/env python3

import subprocess, argparse, sys, os, re, time, getpass
import pandas as pd
from io import StringIO
from mycotools.predb2db import main as predb2db
from mycotools.lib.kontools import intro, outro
from mycotools.lib.dbtools import db2df, df2db, gen_omes
from mycotools.jgiDwnld import jgi_login as jgi_login
from mycotools.jgiDwnld import retrieveXML as retrieveXML
from mycotools.jgiDwnld import JGIdwnld as JGIdwnld


def dwnldMycoCosm( mycocosm_url, output ):

    print('\nDownloading MycoCosm URL')
    subprocess.call(
        ['curl', mycocosm_url, '-o', output + '/mycoCosm.csv'],
        stdout = subprocess.PIPE
        )
    mycocosm = output + '/' + start_time.strftime('%Y%m%d') + '_mycoCosm.csv'

    try:
        jgi_df = pd.read_csv( mycocosm, encoding = 'cp1252' )
    except UnicodeDecodeError:
        jgi_df = pd.read_csv( mycocosm, encoding = 'utf-8' )

    return jgi_df


def breakName( jgi_df, name_col = 'name' ):

    for i, row in jgi_df.iterrows():
        taxa = row['name'].split(' ')
        jgi_df.at[i, 'genus'] = taxa[0]
        if len( taxa ) > 1:
            jgi_df.at[i, 'species'] = taxa[1]
            if len( taxa ) > 2:
                jgi_df.at[i, 'strain'] = '_'.join(taxa[2:])

    return jgi_df


def JGIredundancyCheck( 
        db, jgi_df, ome_col = 'portal', name_col = 'name', update = False 
        ):

    preOmes = set( db.index )
    strain_check = db.apply(lambda x: str(x['genus']) + '_' + str(x['species']) + \
        '_' + str(x['strain']), axis = 1)
    strain_check1 = db.apply(lambda x: str(x['genus']) + '_' + str(x['species']) + \
        '_' + str(x['strain']).replace('_','').replace('-',''), axis = 1)
    strain_check2 = db.apply(lambda x: str(x['genus']) + '_' + str(x['species']) + \
        '_' + str(x['strain']).replace('_',''), axis = 1)
    strain_check3 = db.apply(lambda x: str(x['genus']) + '_' + str(x['species']) + \
        '_' + str(x['strain']).replace('-',''), axis = 1)
    species_check = db.apply(lambda x: str(x['genus']) + '_' + str(x['species']), \
        axis = 1)
    strain_set, species_set = set(strain_check), set(species_check)
    strain_set = strain_set.union( set(strain_check1) )
    strain_set = strain_set.union( set(strain_check2) )
    strain_set = strain_set.union( set(strain_check3) )
    strain_set = strain_set.union( species_set )
    toDel, updates = [], []

    for i, row in jgi_df.iterrows():
        ome = row[ome_col]
        if row[ ome_col ] in preOmes:
            organism = row[ name_col ].replace(' ', '_')
            db.at[ ome, 'published'] = row['publication(s)']
            if organism not in strain_set:
                db_organism = db['genus'][ome] + '_' + db['species'][ome] + '_'  +\
                    str(db['strain'][ome])
                updates.append([ db['internal_ome'][ome], db_organism, organism, ome ])
#                if update:
 #                   db = db.drop( ome )
  #                  continue
            toDel.append( i )
           
    for i in toDel:
        jgi_df = jgi_df.drop( i )

    db = db.reset_index()

    return jgi_df, db, toDel, updates


def main( 
    jgi_df, ref_db, output, assembly = True, 
    proteome = True, gff3 = True, update = True,
    repeatmasked = True
    ):

    if type( jgi_df ) is not pd.DataFrame:
        jgi_df = 'https://mycocosm.jgi.doe.gov/ext-api/mycocosm/catalog/' + \
            'download-group?flt=&seq=all&pub=all&grp=fungi&srt=released&ord=desc'
        jgi_df = dwnldMycoCosm( jgi_df, output )

    if 'genome_code' in jgi_df.columns:
        ome_col = 'genome_code'
    elif 'portal' in jgi_df.columns:
        ome_col = 'portal'
        name_col = 'name'
    else:
        print('\nERROR: no valid MycoCosm header for `ome`')
        sys.exit( 3 )

    print('\nChecking for redundancy')
    toDel = None
    jgi_df = breakName( jgi_df )
    if type(ref_db) is pd.DataFrame:
        db = ref_db.set_index( 'genome_code' )
        old_len = len(jgi_df)
        jgi_df, db, toDel, updates = JGIredundancyCheck( 
            db, jgi_df, ome_col, update = update 
            )
        output_str = 'internal_ome\tdb_organism\tjgi_organism\tgenome_code\n'
        for update in updates:
            output_str += '\t'.join( update ) + '\n'
        with open( output + '/jgiUpdates.tsv', 'w' ) as out:
            out.write( output_str )
        print( '\t' + str(len(updates)) + ' omes can be updated' )
        print('\t' + str(len(toDel)) + ' / ' + str(old_len) + ' omes are present w/o update')
        jgi_df = gen_omes( jgi_df, reference = db )
   #     df2db( 
  #          db, 
 #           output + '/' + os.path.basename(os.path.abspath(args.database)) + '_rmOldies' 
#        )
    else:
        db = None
        jgi_df = gen_omes( jgi_df )

    print('\nLogging into JGI')
    if not os.path.exists( output + '/xml' ):
        os.mkdir( output + '/xml' )

    print( '\nRetrieving `xml` directories' )
    ome_set = set()
    count = 0
    for i, row in jgi_df.iterrows():
        error_check = retrieveXML( row[ ome_col ], output + '/xml' )
        if error_check > 0:
            ome_set.add( row[ ome_col ] )
        if error_check != -1:
            time.sleep( 0.1 )

    print( '\nDownloading' )
    dwnlds = []
    if assembly:
        dwnlds.append( 'assembly' )
    if proteome:
        dwnlds.append( 'proteome' )
    if gff3:
        dwnlds.append( 'gff3' )

    for typ in dwnlds:
        if not os.path.isdir( output + '/' + typ ):
            os.mkdir( output + '/' + typ )
        if typ == 'gff3':
            if not os.path.isdir( output + '/gff' ):
                os.mkdir( output + '/gff')

    for i, row in jgi_df.iterrows():
        ome = row[ome_col]
        species_strain = re.search(
            r'(.*?)_(.*?)_(.*)',
            jgi_df['name'][i].replace(' ', '_')
            )
        if species_strain is not None:
            jgi_df.at[i, 'genus'] = species_strain[1]
            jgi_df.at[i, 'species'] = species_strain[2]
            jgi_df.at[i, 'strain'] = species_strain[3]
        elif '_' in row['name'].replace(' ', '_'):
            genus_species = row['name'].split('_')
            jgi_df.at[i, 'genus'] = genus_species[0]
            jgi_df.at[i, 'species'] = genus_species[1]
        else:
            jgi_df.at[i, 'genus'] = row['name']
            jgi_df.at[i, 'species'] = 'sp.'
            
        if ome not in ome_set:
            jgi_login( user, pwd )
            print( '\t' + ome + ': ' + row['name'] )
            for typ in dwnlds:
                check, preexisting, new_typ = JGIdwnld( 
                    ome, typ, output, masked = repeatmasked
                )
                if type(check) != int:
                    jgi_df.at[i, new_typ + '_path'] = output + '/' + new_typ + \
                        '/' + check
                    check = os.path.basename( os.path.abspath( check ) )
                print('\t\t' + new_typ + ': ' + str(check))
                if not preexisting:
                    time.sleep( 60 )


    print( '\nPreparing predb' )
    if os.path.exists( 'cookies' ):
        os.remove( 'cookies' )
    if os.path.exists( os.path.expanduser( '~/.nullJGIdwnld' )):
        os.remove( os.path.expanduser( '~/.nullJGIdwnld' ))

    jgi_df = jgi_df.rename( columns = { 
        'publication(s)': 'publication',
        ome_col: 'genome_code'    
        } )
    jgi_df['source'] = 'jgi'

#    if args.update:
 #       jgi_df = jgi_df.set_index( 'genome_code' )
  #      updates_df = pd.DataFrame()
   #     for up in updates:
    #        updates_df = updates_df.append( jgi_df.loc[up[3]] )
     #       jgi_df = jgi_df.drop( up[3] )
      #  jgi_df = jgi_df.reset_index()
       # print('\t' + str(len(updates_df)) + ' omes ready for update using `predb2db.py`')
	#updates_df['genome_code'] = updates_df.index
        #updates_df.to_csv( output + '/omes2update.predb', sep = '\t' )
            
    jgi_df = predb2db( jgi_df, db )


if __name__ == '__main__':

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
    parser.add_argument( '-g3', '--gff3', default = True, action = 'store_false', \
        help = 'Do not download gff3s.' )
    parser.add_argument( '-r', '--repeatmasked', default = True, action = 'store_false', \
        help = "Download nonmasked assemblies." )
    parser.add_argument( '-g', '--gff', default = False, action = 'store_true', \
        help = 'Download gffs.' )

    args = parser.parse_args()

    args_dict = {
            'Preexisting db': args.database,
            'MycoCosm Table': args.mycocosm,
            'Assemblies': args.assembly,
            'RepeatMasked': args.repeatmasked,
            'Proteomes': args.proteome,
            ".gff's": args.gff,
            ".gff3's": args.gff3,
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

    ref_db = db2df( formatPath(args.database) )
    jgi_df = main( args.mycocosm, refdb, output )

    df2db( jgi_df, output + '/new.db' )
    print('\nSuccess! ' + str(len(jgi_df)) + ' added to database\n \
            Run updateDB to confirm and finish update.')

    outro( start_time )
