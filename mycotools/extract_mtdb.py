#! /usr/bin/env python3

# NEED multiple lineages from command line
# NEED stdin acceptance for most of these arguments

import os
import re
import sys
import copy
import random
import argparse
from collections import defaultdict
from mycotools.lib.kontools import file2list, intro, outro, format_path, \
    eprint, mkOutput
from mycotools.lib.dbtools import mtdb, primaryDB
from mycotools.db2files import mtdb_main as gen_full_mtdb

# NEED to fix when same lineage multiple ranks, e.g. Tremellales sp. will be listed
# as an order and as a genus

def infer_rank(db, lineage):
    """Identify the taxonomic rank associated with an inputted lineage of
    interest"""
    linlow, rank = lineage.lower(), None
    for ome, row in db.items():
        if linlow in set([x.lower() for x in row['taxonomy'].values()]):
            rev_dict = {k.lower(): v \
                        for k,v in zip(row['taxonomy'].values(), 
                                       row['taxonomy'].keys())}
            rank = rev_dict[lineage]
            

    if not rank:
        raise KeyError(f'no entry for {lineage}')

    return rank
                            

def extract_unique(db, allowed = 1, rank = 'species'):
    """Extract unique rank from an MTDB"""
    keys = copy.deepcopy(list(db.keys()))
    random.shuffle(keys)
    prep_db0 = {x: db[x] for x in keys}
    prep_db1 = mtdb().set_index('ome')
    if rank == 'strain':
        found = set()
        for ome, row in prep_db0.items():
            name = row['taxonomy']['species'] + ' ' + row['strain']
            if name not in found:
                prep_db1[ome] = row
            found_prep = list(found)
            found_prep.append(name)
            found = set(found_prep)
    else:
        found = defaultdict(int)
        for ome, row in prep_db0.items():
            name = row['taxonomy'][rank]
            found[name] += 1
            if found[name] <= allowed:
                prep_db1[ome] = row
            
    return prep_db1

def extract_tax(db, lineages):
    """Extract specific taxonomic lineages of interest based on their rank"""
    if isinstance(lineages, str):
        lineages = [lineages]
    lineages = set(x.lower() for x in lineages)
    rank_dict = {k: infer_rank(db, k) for k in list(lineages)}
    ranks = list(set(rank_dict.values()))

    new_db = mtdb().set_index()
    for ome in db:
        for rank in ranks:
            try:
                if db[ome]['taxonomy'][rank].lower() in lineages:
                    new_db[ome] = db[ome]
            except KeyError: # invalid rank key for row
                pass # probably should standardize tax jsons period

    return new_db

def extract_ome(db, omes):
    """Extract a list of genome codes (omes) of interest"""
    new_db = mtdb().set_index()
    for i in db:
        if i in list(omes):
            new_db[i] = db[i]
    return new_db

def extract_source(db, source):
    """Extract an MTDB with genomes from a particular source"""
    return mtdb({ome: row for ome, row in db.items() \
                   if row['source'].lower() == source.lower()}, 
                   index = 'ome')

def extract_pub(db):
    """Extract only published and usable genomes"""
    new_db = mtdb().set_index()
    for ome, row in db.items():
        if row['published']:
            new_db[ome] = row
    return new_db

def main( 
    db, rank = None, x_number = 0, 
    lineage_list = [], omes_set = set(), by_rank = False,
    source = None, nonpublished = False, inverse = False
    ):
    """Python entry point for extract_mtdb"""

    db = db.set_index('ome')
    if x_number > 0:
        db = extract_unique(db, x_number, rank = rank)

    # extract each taxonomic entry based on the classification specified
    if lineage_list:
        new_db = extract_tax(db, lineage_list)
    # if an ome list is specified then open it, store each entry in a list and pull each ome
    elif omes_set:
        new_db = extract_ome(db, omes_set)
    # if none of these are specified then create a `new_db` variable to work for later
    else:
        new_db = db

    # if there is a source specified, extract it or the opposite
    if source:
        new_db = extract_source(new_db, source)

    # if you want publisheds, then just pull those out 
    if not nonpublished:
        new_db = extract_pub(new_db)

    if inverse:
        new_omes = set(new_db.keys())
        inv_db = mtdb().set_index()
        for ome, row in db.items():
            if ome not in new_omes:
                inv_db[ome] = row
        new_db = inv_db

    if by_rank:
        lineages = set(v['taxonomy'][rank] for k,v in new_db.items())
        dbs = {}
        for lineage in lineages:
            if lineage:
                dbs[lineage.lower()] = extract_tax(new_db, [lineage]).reset_index()
            else:
                dbs['unclassified'] = extract_tax(new_db, ['']).reset_index()
        return dbs
    else:
        return new_db.reset_index()


def cli():
    ranks = ['kingdom', 'phylum', 'subphylum', 'class', 'order',
             'family', 'genus', 'species', 'strain']
    parser = argparse.ArgumentParser( description = \
       'Extracts a MycotoolsDB from arguments. E.g.\t`mtdb extract ' + \
       '-l Atheliaceae`' )
    parser.add_argument('-l', '--lineage', default = '')
    parser.add_argument('-s', '--source')
    parser.add_argument('-n', '--nonpublished', action = 'store_true', 
        help = 'Include restricted')
    parser.add_argument('-r', '--rank', 
        help = f'[-a|-b] Taxonomic rank: {ranks}')
    parser.add_argument('-a', '--allowed_rank', default = 0, type = int, 
        help = '[-r] Number of randomly sampled --rank allowed' )
    parser.add_argument('-b', '--by_rank', action = 'store_true',
        help = '[-r] Output MTDBs for each lineage in --rank')
    parser.add_argument('-i', '--inverse', action = 'store_true', 
        help = 'Inverse [source|lineage(s)|nonpublished]')
    parser.add_argument('-ol', '--ome', help = "File w/list of omes" )
    parser.add_argument('-ll', '--lineages', help = 'File w/list of lineages' )
    parser.add_argument('-m', '--new_mtdb', action = 'store_true',
                        help = 'Create MTDB directory hierarchy')
    parser.add_argument('--headers', action = 'store_true')
    parser.add_argument('-d', '--mtdb', help = '- for stdin', default = primaryDB())
    parser.add_argument('-o', '--output')
    args = parser.parse_args()
    db_path = format_path(args.mtdb)


    if args.lineage or args.lineages:
        eprint("\nWARNING: extracting taxonomy is subject to " \
             + "errors in NCBI's hierarchy\n")

    # these arguments require one another
    if args.by_rank and not args.rank:
        eprint('\nERROR: --by_rank requires --rank', flush = True)
        sys.exit(10)
    elif args.allowed_rank and not args.rank:
        eprint('\nERROR: --allowed_rank requres --rank', flush = True)
        sys.exit(11)
    elif args.rank and not args.allowed_rank and not args.by_rank:
        eprint('\nERROR: --rank requires --allowed_rank or --by_rank', 
               flush = True)
    elif args.rank:
        args.rank = args.rank.lower()
        if args.rank not in set(ranks):
            eprint(f'\nERROR: --rank not in {ranks}', flush = True)
            sys.exit(12)

    args.lineage = args.lineage.replace('"','').replace("'",'')

    output = ''
    if args.output:
        output = format_path(args.output)
        if not output.endswith('/'):
            tag = ''
                       
            if args.lineage:
                tag += '_' + args.lineage
            if args.lineages:
                tag += '_taxonomy'
            if args.source:
                tag += args.source.lower()
            if not args.nonpublished:
                tag += '_pub'
            output += '/' + os.path.basename(db_path) + tag
        
    if args.mtdb == '-':
        data = ''
        for line in sys.stdin:
            data += line.rstrip() + '\n'
        data = data.rstrip()
        db = mtdb(data, stdin=True)
    else:
       db = mtdb(db_path)

    if args.ome:
        omes = set(file2list(format_path(args.ome)))
    else:
        omes = set()

    if args.lineages:
        lineage_list = file2list(format_path(args.lineages))
    elif args.lineage:
        lineage_list = [args.lineage]
    else:
        lineage_list = []

    new_db = main( 
        db, lineage_list = lineage_list,
        omes_set = omes, source = args.source, rank = args.rank, 
        x_number = args.allowed_rank, by_rank = args.by_rank,
        nonpublished = args.nonpublished, inverse = args.inverse
        )
    if args.new_mtdb:
        gen_full_mtdb(new_db, format_path(output),
            format_path(os.environ['MYCODB'] + '/../'))
    elif args.output or args.by_rank:
        if isinstance(new_db, mtdb):
            new_db.df2db(output)
        else:
            out_dir = mkOutput(output, 'extract_mtdb')
            prefix = re.sub(r'\.mtdb$', '', os.path.basename(db_path))
            for lineage, db in new_db.items():
                out_f = f'{out_dir}{prefix}.{lineage}.mtdb'
                db.df2db(out_f)
    else:
        new_db.df2db(headers = bool(args.headers))

    sys.exit(0)


if __name__ == '__main__':
    cli()
