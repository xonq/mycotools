#! /usr/bin/env python3

import argparse, sys, os, re, pandas as pd
from mycotools.lib.dbtools import db2df, df2db, masterDB
from mycotools.lib.kontools import formatPath, prep_output
from shutil import copy


def main( filetypes, envtypes, db, output_path, link = False ):

    filetypes = {x: filetypes[x] for x in filetypes if filetypes[x]}
    if not link:
        for ftype in filetypes:
            if not os.path.isdir( output_path + ftype ):
                os.mkdir( output_path + ftype )

        for i, row in db.iterrows():
            for ftype in filetypes:
                if not pd.isnull( row[ ftype ] ):
                    fpath = formatPath( '$' + envtypes[ftype] + '/' + row[ftype] )
                    copy( fpath, output_path + ftype + '/' + os.path.basename( fpath ) ) 
    else:
        for i, row in db.iterrows():
            for ftype in filetypes:
                if not pd.isnull( row[ ftype ] ):
                    fpath = formatPath( '$' + envtypes[ftype] + '/' + row[ftype] )
                    print( fpath , flush = True)


if __name__ == '__main__':

    parser = argparse.ArgumentParser( description = 'Copies selected files from database' )
    parser.add_argument( '-a', '--assembly', action = 'store_true', help = 'Grab assemblies' )
    parser.add_argument( '-p', '--proteome', action = 'store_true', help = 'Grab proteomes' )
    parser.add_argument( '-g', '--gff', action = 'store_true', help = 'Grab gff`s' )
    parser.add_argument( '-l', '--link', action = 'store_true', help = 'Print paths, no copy' )
    parser.add_argument( '-d', '--database', default = masterDB(), help = 'DEFAULT: masterDB' )
    parser.add_argument( '-o', '--output', default = os.getcwd() )
    parser.add_argument( '--yes', action = 'store_true', help = 'Assume yes' )
    args = parser.parse_args()

    db_path = formatPath( args.database )
    output_path = prep_output( formatPath( args.output ), cd = False )
    args_dict = {
        'DATABASE': db_path, 'OUTPUT': output_path,
        'ASSEMBLY': args.assembly, 'PROTEOME': args.proteome,
        'GFF3': args.gff, 'Links': args.link
        }

    filetypes = { 
        'assembly': args.assembly, 'proteome': args.proteome,
        'gff3': args.gff
        }
    envtypes = { 
        'assembly': 'MYCOFNA', 'proteome': 'MYCOFAA',
        'gff3': 'MYCOGFF3'
        }

    db = db2df( db_path )
    max_files = len(db) * len([x for x in filetypes if filetypes[x]])
    if max_files > 100:
        if not args.yes and not args.link:
            print('\nWARNING: the maximum size pull for your cmd is ' + \
                str( max_files ) + ' files' )
            check = input('Proceed? [y/N]: ')
            if check.lower() not in { 'y', 'yes' }:
                print(flush = True)
                sys.exit( 1 )
    main( filetypes, envtypes, db, output_path, link = args.link )
    sys.exit( 0 )
