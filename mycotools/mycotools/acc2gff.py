#! /usr/bin/env python3

# NEED to remove pandas requirement

import multiprocessing as mp
import os, sys, argparse, re
from mycotools.lib.biotools import gff2list, list2gff
from mycotools.lib.dbtools import mtdb, masterDB
from mycotools.lib.kontools import formatPath

def grabGffAcc( gff_list, acc ):
    '''grab acc from alias'''

    if ';Alias=' in gff_list[0]['attributes']:
        alias = 'Alias=' + acc
        alias_on = alias + ';'
    elif ' alias "' in gff_list[0]['attributes']:
        alias = 'alias "' + acc + '"'
        alias_on = 'alias "' + acc + '"'
    out_list = [ 
        x for x in gff_list \
        if x['attributes'].endswith(alias) or \
            alias_on in x['attributes'] 
        ]
    return out_list


def grabGffAccs(gff_list, acc_list, ome = None):

    aliases, ends = set(), set()
    for acc in acc_list:
        aliases.add('Alias=' + acc)
        ends.add('Alias=' + acc + ';')
    out_list = []
    for i in gff_list:
        if any(i['attributes'].endswith(x) for x in aliases):
            out_list.append(i)
        elif any(x in i['attributes'] for x in ends):
            out_list.append(i)

    if ome:
        return out_list, ome
    else:
        return out_list


def gffMain(gffData, accs):
    accGffs = {}
    if isinstance(gffData, str):
        gff = gff2list(gffData)
    else:
        gff = gffData
    for acc in accs:
        accGffs[acc] = grabGffAcc(gff, acc)
    return accGffs

def dbMain(db, accs, cpu = 1):

    omes = set([x[:x.find('_')] for x in accs])
    
    grabAcc_cmds = []
    for ome in list(omes):
        omeAccs = [acc for acc in accs if acc.startswith(ome + '_')]
        gffPath = formatPath('$MYCOGFF3/' + ome + '.gff3')
        gff_list = gff2list(gff)
        grabAcc_cmds.append([omeAccs, gff_list, ome])

    with mp.Pool(processes = cpus) as pool:
        grabAcc_res = pool.starmap(grabAccs, grabAcc_cmds)

    gffs = {}
    for res in grabAcc_res:
        gffs[res[1]] = res[0]

    return gffs


if __name__ == '__main__':

    mp.set_start_method('spawn')
    parser = argparse.ArgumentParser( 
        description = 'Inputs gff or database and new line delimitted file of accessions.'
        )
    parser.add_argument( '-a', '--accession' )
    parser.add_argument( '-i', '--input', help = 'File with accessions' )
    parser.add_argument( '-g', '--gff', help = 'GFF3 input' )
    parser.add_argument( '-c', '--column', default = 1, 
        help = 'Accssions column for -i (1 indexed). DEFAULT: 1' )
    parser.add_argument( '-o', '--ome', action = 'store_true', help = 'Output files by ome code' )
    parser.add_argument( '-d', '--database', default = masterDB(), help = 'mycodb DEFAULT: master' )
    parser.add_argument( '--cpu', type = int, default = mp.cpu_count(), \
        help = 'CPUs to parallelize accession searches. DEFAULT: all' )
    args = parser.parse_args()

    if args.input:
        input_file = formatPath( args.input )
        with open(input_file, 'r') as raw:
            accs = [x.rstrip().split('\t')[args.column-1] for x in raw]
    elif not args.accession:
        print('\nERROR: need input file or accession', flush = True)
        sys.exit( 1 )
    else:
        accs = [args.accession]

    if args.cpu < mp.cpu_count():
        cpu = args.cpu
    else:
        cpu = mp.cpu_count()

    db_path = formatPath( args.database )
    if not args.gff:
        db = mtdb( formatPath(args.database) )
        gff_lists = dbMain( db, accs, cpu = args.cpu )
    else:
        gff_path = formatPath( args.gff )
        gff_lists = gffMain(gffData, accs)

    if args.accession:
        print( gff2list(gff_lists[list(gff_lists.keys())[0]].rstrip()) , flush = True)
    elif args.ome:
        output = mkOutput(os.getcwd() + '/', 'acc2gff')
        for ome in gff_strs:
            with open( output + ome + '.accs.gff3', 'w' ) as out:
                out.write( list2gff(gff_lists[ome]) )
    else:
        out_str = ''
        for ome in gff_lists:
            out_str += list2gff(gff_lists[ome]) + '\n'
        print(out_str)

    sys.exit( 0 )
