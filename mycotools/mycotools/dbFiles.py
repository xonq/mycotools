#! /usr/bin/env python3

import argparse, sys, os, re, pandas as pd
from mycotools.lib.dbtools import db2df, df2db, masterDB
from mycotools.lib.kontools import formatPath, prep_output
from shutil import copy


def softMain( filetypes, envtypes, db, output_path, print_link = False ):

    filetypes = {x: filetypes[x] for x in filetypes if filetypes[x]}
    if not print_link:
        for ftype in filetypes:
            if not os.path.isdir(output_path + ftype):
                os.mkdir(output_path + ftype)
        for ftype in filetypes:
            check_db = db[~db[ftype].isnull()]
            check_db.apply(
                lambda x: os.symlink(envtypes[ftype] + '/' + x[ftype],
                output_path + ftype + '/' + x[ftype]), axis = 1
                )
        
    else:
        for i, row in db.iterrows():
            for ftype in filetypes:
                if not pd.isnull( row[ ftype ] ):
                    fpath = formatPath( '$' + envtypes[ftype] + '/' + row[ftype] )
                    print( fpath , flush = True)


def hardMain( filetypes, envtypes, db, output_path ):

    filetypes = {x: filetypes[x] for x in filetypes if filetypes[x]}
    for ftype in filetypes:
        if not os.path.isdir( output_path + ftype ):
            os.mkdir( output_path + ftype )

    for i, row in db.iterrows():
        for ftype in filetypes:
            if not pd.isnull( row[ ftype ] ):
                fpath = formatPath( '$' + envtypes[ftype] + '/' + row[ftype] )
                copy( fpath, output_path + ftype + '/' + os.path.basename( fpath ) ) 

if __name__ == '__main__':

    parser = argparse.ArgumentParser( description = 'Symlinks/copies selected files from database' )
    parser.add_argument( '-d', '--database', default = masterDB(), help = 'DEFAULT: masterDB' )
    parser.add_argument( '-a', '--assembly', action = 'store_true', help = 'Grab assemblies' )
    parser.add_argument( '-p', '--proteome', action = 'store_true', help = 'Grab proteomes' )
    parser.add_argument( '-g', '--gff', action = 'store_true', help = 'Grab gff`s' )
    parser.add_argument( '--print', action = 'store_true', help = 'Print paths, no copy' )
    parser.add_argument( '--hard', action = 'store_true', help = 'Hard copy files' )
    parser.add_argument( '-o', '--output', default = os.getcwd() )
    args = parser.parse_args()

    if not args.assembly and not args.proteome and not args.gff:
        print('\nERROR: no file type (-a, -g, -p) selected', flush = True)
        sys.exit(4)

    db_path = formatPath( args.database )
    output_path = prep_output( formatPath( args.output ), cd = False )
    args_dict = {
        'DATABASE': db_path, 'OUTPUT': output_path,
        'ASSEMBLY': args.assembly, 'PROTEOME': args.proteome,
        'GFF3': args.gff, 'Print links': args.print, 'Hard copy': args.hard
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
    if args.print or not args.hard:
        softMain(filetypes, envtypes, db, output_path, print_link = args.print)
    else:
        hardMain( filetypes, envtypes, db, output_path )

    sys.exit( 0 )
