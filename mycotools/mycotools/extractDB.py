#! /usr/bin/env python3

import pandas as pd, numpy as np
import argparse, os, sys
from mycotools.lib.kontools import file2list, intro, outro, formatPath, eprint
from mycotools.lib.dbtools import mtdb, extract_tax, extract_omes, masterDB

def main( 
    db, rank = False, unique_species = False, lineage_list = False,
    inverse = False, lineage = False, ome_list = False,
    source = False, nonpublished = False, unique_strains = False
    ):

    if unique_species:
        db['full_name'] = db.apply(
            lambda x: str(x['genus']) + '_' + str(x['species']), axis = 1
            )
        db = db.drop_duplicates('full_name')
    elif unique_strains:
        db['full_name'] = db.apply(
            lambda x: str(x['genus']) + '_' + str(x['species']) + str(x['strain']), axis = 1
            )
        db = db.drop_duplicates('full_name')

    # if a taxonomy list is specified then open it, store each entry in a list
    # and extract each taxonomic entry based on the classification specified
    if lineage_list:
        taxonomy_list = file2list( formatPath(lineage_list) )
        if not taxonomy_list:
            sys.exit( 15 )

        if rank:
            new_db = mtdb()
            new_db = new_db.append( 
                extract_tax( 
                    db, taxonomy_list, 
                    rank, inverse = inverse 
                    ) 
                )

    # if there is one taxonomic group specified, extract it 
    elif lineage:
        new_db = extract_tax(
            db, lineage, rank, inverse = inverse )

    # if an ome list is specified then open it, store each entry in a list and pull each ome
    elif ome_list:
        omes = set(file2list( formatPath(ome_list) ))
        db = db.set_index('internal_ome')
        new_db = mtdb()
        new_db = new_db.set_index()
        for i in db:
            if i in omes:
                new_db[i] = db[i]
        new_db = new_db.reset_index() 
#        new_db = extract_omes( db, omes, inverse = inverse )

    # if none of these are specified then create a `new_db` variable to work for later
    else:
        new_db = db

    # if there is a source specified, extract it or the opposite if inverse is specified
    if source and not inverse:
        new_db = new_db[new_db['source'] == source]
    elif source and inverse:
        new_db = new_db[new_db['source'] != source]

    # if you want publications, then just pull those out 
    new_db = new_db.set_index()
    if not nonpublished:
        todel = []
        for ome in new_db:
            if not new_db[ome]['published']:
#                new_db.at[i, 'published'] = 0
#                new_db = new_db.drop(i)
                todel.append(ome)
 #           else:
#                new_db['published'][i] = .at[i, 'published'] = row['published']
        for v in todel:
            del new_db[v]

    return new_db


if __name__ == '__main__':

    parser = argparse.ArgumentParser( description = \
       'Extracts a MycotoolsDB from arguments. E.g.\t`extractDB.py ' + \
       '-l Atheliaceae -r family`' )
    parser.add_argument( '-d', '--database', default = masterDB(), help = 'DEFAULT: masterdb' )
    parser.add_argument( '-l', '--lineage' )
    parser.add_argument( '-r', '--rank', help = "Taxonomy rank" )
    parser.add_argument( '-s', '--source', help = 'Data source' )
    parser.add_argument( '-n', '--nonpublished', action = 'store_true', help = 'Include ' + \
        'restricted-use' )
    parser.add_argument( '-u', '--unique', action = 'store_true', help = 'Unique species' )
    parser.add_argument( '--unique_strains', action = 'store_true' )
    parser.add_argument( '-i', '--inverse', action = 'store_true', help = 'Inverse arguments' )
    parser.add_argument( '--headers', action = 'store_true', help = 'Include header')
    parser.add_argument( '-o', '--output' )
    parser.add_argument( '-', '--stdin', action = 'store_true', help = "Pipe MycotoolsDB from stdin" )
    parser.add_argument( '-ol', '--ome', help = "File w/list of omes" )
    parser.add_argument( '-ll', '--lineages', help = 'File w/list of lineages (same rank)' )
    args = parser.parse_args()
    db_path = formatPath( args.database )


    if (args.lineage or args.lineages) and not args.rank:
        eprint('\nERROR: need rank for lineage(s)', flush = True)
        sys.exit(5)
    elif args.lineage or args.lineages:
        eprint("\nWARNING: extracting taxonomy is subject to taxonomic error and errors in NCBI's hiearchy\n")

    output = 'stdout'
    if args.output:
        output = formatPath( args.output )
        if not output.endswith( '/' ):
            tag = ''
                       
            if args.lineage:
                tag += '_' + args.lineage
            if args.lineages:
                tag += '_taxonomy'
            if args.source:
                tag += args.source.lower()
            if not args.nonpublished:
                tag += '_pub'
            output += '/' + os.path.basename( db_path ) + tag
        
    args_dict = {
        'database': db_path,
        'output': output,
        'lineage': args.lineage,
        'lineages list': args.lineages,
        'rank': args.rank,
        'ome_list': args.ome,
        'Source': args.source,
        'unique': args.unique,
        'nonpublished': args.nonpublished,
        'inverse': args.inverse,
        'headers': bool(args.headers)
    	}
#    start_time = intro( 'Abstract Database', args_dict, stdout = args.output )

    if args.stdin:
        data = ''
        for line in sys.stdin:
            data += line.rstrip() + '\n'
        data = data.rstrip()
        db = mtdb(data, stdin=True)
    else:
       db = mtdb( db_path )
    new_db = main( 
        db, rank = args.rank, lineage = args.lineage, lineage_list = args.lineages,
        ome_list = args.ome, source = args.source, unique_species = args.unique, 
        nonpublished = args.nonpublished, inverse = args.inverse, unique_strains = args.unique_strains
        )
    if args.output:
        new_db.df2db( output )
    else:
        new_db.df2db( headers = bool(args.headers) )

 #   outro( start_time, stdout = args.output )
    sys.exit(0)
