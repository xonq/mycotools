#! /usr/bin/env python3

import os
import re
import sys
import argparse
import multiprocessing as mp
from mycotools.lib.biotools import gff2list, list2gff
from mycotools.lib.dbtools import mtdb, primaryDB
from mycotools.lib.kontools import format_path, stdin2str

def grabGffAcc( gff_list, acc, term = 'Alias=' ):
    '''grab acc from alias'''

    alias = term + acc
    alias_on = alias + ';'
#    elif ' alias "' in gff_list[0]['attributes']:
 #       alias = 'alias "' + acc + '"'
  #      alias_on = 'alias "' + acc + '"'
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


def dbMain(db, accs, cpus = 1):

    omes = set([x[:x.find('_')] for x in accs])
    db = db.set_index('ome')
    grabAcc_cmds = []
    for ome in list(omes):
        omeAccs = [acc for acc in accs if acc.startswith(ome + '_')]
        gff_list = gff2list(db[ome]['gff3'])
#        grabAcc_res = [grabGffAccs(gff_list, omeAccs, ome)]
        grabAcc_cmds.append([gff_list, omeAccs, ome])

    with mp.Pool(processes = cpus) as pool:
        grabAcc_res = pool.starmap(grabGffAccs, grabAcc_cmds)

    gffs = {}
    for res in grabAcc_res:
        gffs[res[1]] = res[0]

    return gffs


def cli():

#    mp.set_start_method('spawn')
    parser = argparse.ArgumentParser( 
        description = 'Inputs gff or database and new line delimitted file of accessions.'
        )
    parser.add_argument('-a', '--accession', help = '"-" for stdin')
    parser.add_argument('-i', '--input', help = 'File with accessions')
    parser.add_argument('-g', '--gff', help = 'GFF3 input')
    parser.add_argument('-c', '--column', default = 1, 
        help = 'Accssions column for -i (1 indexed). DEFAULT: 1')
    parser.add_argument('-o', '--ome', action = 'store_true', help = 'Output files by ome code')
    parser.add_argument('-d', '--mtdb', default = primaryDB(), help = 'mycodb DEFAULT: master')
    parser.add_argument('--cpu', type = int, default = mp.cpu_count())
    args = parser.parse_args()

    if args.input:
        input_file = format_path( args.input )
        with open(input_file, 'r') as raw:
            accs = [x.rstrip().split('\t')[args.column-1] for x in raw]
    elif not args.accession:
        raise ValueError('need input file or accession')
    else:
        if '-' == args.accession:
            data = stdin2str()
            accs = data.split()
        else:
            if {'"', "'"}.intersection(set(args.accession)):
                args.accession = args.accession.replace('"','').replace("'",'')
            if ',' in args.accession:
                accs = args.accession.split(',')
            elif re.search(r'\s', args.accession):
                accs = args.accession.split()
            else:
                accs = [args.accession]

    if args.cpu < mp.cpu_count():
        cpu = args.cpu
    else:
        cpu = mp.cpu_count()

    db_path = format_path( args.mtdb )
    if not args.gff:
        db = mtdb( format_path(args.mtdb) )
        gff_lists = dbMain( db, accs, cpus = args.cpu )
    else:
        gff_path = format_path( args.gff )
        gff_lists = gffMain(gffData, accs)

    if args.accession:
        print( list2gff(gff_lists[list(gff_lists.keys())[0]]).rstrip() , flush = True)
    elif args.ome:
        output = mkOutput(os.getcwd() + '/', 'acc2gff')
        for ome in gff_strs:
            if gff_lists[ome]:
                with open( output + ome + '.accs.gff3', 'w' ) as out:
                    out.write( list2gff(gff_lists[ome]) )
            else:
                eprint('ERROR: ' + ome + ' failed, no accessions retrieved', flush = True)
    else:
        out_str = ''
        for ome in gff_lists:
            if gff_lists[ome]:
                out_str += list2gff(gff_lists[ome]) + '\n'
            else:
                eprint('ERROR: ' + ome + ' does not have accession', flush = True)
        print(out_str)

    sys.exit( 0 )

if __name__ == '__main__':
    cli()
