#! /usr/bin/env python3

import pandas as pd, numpy as np
import argparse, os, sys
from mycotools.lib.kontools import file2list, intro, outro, formatPath
from mycotools.lib.dbtools import db2df, df2db, abstract_tax, abstract_omes, masterDB

def main( args_dict, db ):

    # if a taxonomy list is specified then open it, store each entry in a list
    # and abstract each taxonomic entry based on the classification specified
    if args.taxonomy_list:
        taxonomy_list = file2list( args.taxonomy_list )
        if not taxonomy_list:
            sys.exit( 15 )

        if args.classification:
            new_db = pd.DataFrame()
            new_db = new_db.append( 
                abstract_tax( 
                    db, taxonomy_list, 
                    args.classification, inverse = args.inverse 
                    ) 
                )

    # if there is one taxonomic group specified, abstract it 
    elif args.taxonomy:
        new_db = abstract_tax(db, args.taxonomy, args.classification, inverse = args.inverse )

    # if an ome list is specified then open it, store each entry in a list and pull each ome
    elif args.ome:
        ome_list = file2list( args.ome )
        new_db = abstract_omes( db, ome_list, inverse = args.inverse )

    # if none of these are specified then create a `new_db` variable to work for later
    else:
        new_db = db

    # if there is a source specified, abstract it or the opposite if inverse is specified
    if args.source and not args.inverse:
        new_db = new_db[new_db['source'] == args.source]
    elif args.source and args.inverse:
        new_db = new_db[new_db['source'] != args.source]

    # if you want publications, then just pull those out 
    if not args.nonpublished:
        for i, row in new_db.iterrows():
            if pd.isnull(row['published']) or not row['published'] or str(row['published']) == "0.0":
                new_db.at[i, 'published'] = 0
            else:
                new_db.at[i, 'published'] = row['published']

    

    return new_db


if __name__ == '__main__':

    parser = argparse.ArgumentParser( description = 'Imports database and file of ' + \
       'omes or taxonomy. Abstracts database from parameters. E.g.\t`abstractDB.py ' + \
       '-i $DATABASE -t Atheliaceae -c family`' )

    parser.add_argument( '-t', '--taxonomy', help = 'Taxonomic group to abstract. ' + \
        'Requires `-c`' )
    parser.add_argument( '-c', '--classification', help = 'Taxonomic classification to ' + \
        'abstract from. Requires `-t` or `-tl`' )
    parser.add_argument( '-s', '--source', help = 'Genome source to abstract' )
    parser.add_argument( '-n', '--nonpublished', action = 'store_true', help = 'Abstract ' + \
        'nonpublished rows.' )
    parser.add_argument( '--ome', help = "New line separated list of omes to include." )
    parser.add_argument( '-tl', '--taxonomy_list', help = 'New line separated list of ' + \
        'taxonomic groups - must be same classification. Requires `-c`' )
    parser.add_argument( '--inverse', action = 'store_true', help = 'Inverse non-binary arguments' )
    parser.add_argument( '--headers',  default = False, action = 'store_true', 
        help = 'Output db headers' )
    parser.add_argument( '-i', '--input', default = masterDB(), help = "Kontools database. " + \
        'DEFAULT: master database' )
    parser.add_argument( '-o', '--output', help = "Output path instead of print to stdout. " + \
        'Includes column headers (stdout does not)' )
    args = parser.parse_args()

    db_path = formatPath( args.input )

    output = 'stdout'
    if args.output:
        output = formatPath( args.output )
        if not output.endswith( '/' ):
            tag = ''
            if args.taxonomy:
                tag += '_' + args.taxonomy
            if args.taxonomy_list:
                tag += '_taxonomy'
            if args.source:
                tag += args.source.lower()
            if not args.nonpublished:
                tag += '_pub'
            output += '/' + os.path.basename( db_path ) + tag
        
    args_dict = {
        'Database': db_path,
        'Output': output,
        'Ome list': args.ome,
        'Taxonomy': args.taxonomy,
        'Taxonomy List': args.taxonomy_list,
        'Classification': args.classification,
        'Source': args.source,
        'Nonpublished': args.nonpublished,
        'Inverse': args.inverse,
        'Headers': args.headers
    }
    start_time = intro( 'Abstract Database', args_dict, stdout = args.output )

    db = db2df( db_path )
    new_db = main( args_dict, db )
    if args.output:
        df2db( new_db, output )
    else:
        df2db( new_db, sys.stdout, header = args.headers )

    outro( start_time, stdout = args.output )
