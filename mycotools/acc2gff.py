#! /usr/bin/env python3

import os
import re
import sys
import argparse
import multiprocessing as mp
from mycotools.lib.biotools import gff2list, list2gff
from mycotools.lib.dbtools import mtdb, primaryDB
from mycotools.lib.kontools import format_path, stdin2str

def grab_gff_acc( gff_list, acc, term = 'Alias=' ):
    '''grab acc from alias'''

    alias = term + acc
    alias_on = alias + ';'
    out_list = [ 
        x for x in gff_list \
        if x['attributes'].endswith(alias) or \
            alias_on in x['attributes'] 
        ]
    return out_list


def grab_gff_accs(gff_list, acc_list, ome = None):
    """grab entries with any of a list of aliases"""
   # aliases, ends = set(), set()
#    for acc in acc_list:
 #       aliases.add('Alias=' + acc)
  #      ends.add('Alias=' + acc + ';')
    acc_set = set(acc_list)
    out_list = []
    for i in gff_list:
        alia = re.search(r'Alias=([^;]+)', i['attributes'])
        try:
            aliases = set(alia[1].split('|'))
        except TypeError:
            continue
        if acc_set.intersection(aliases):
            out_list.append(i)
#        if any(i['attributes'].endswith(x) for x in aliases):
 #           out_list.append(i)
  #      elif any(x in i['attributes'] for x in ends):
   #         out_list.append(i)

    if ome:
        return out_list, ome
    else:
        return out_list


def gff_main(gff_data, accs):
    """acquire a dictionary of gffs for each accession as a key"""
    acc_gffs = {}
    if isinstance(gff_data, str):
        gff = gff2list(gff_data)
    else:
        gff = gff_data
    for acc in accs:
        acc_gffs[acc] = grab_gff_acc(gff, acc)
    return acc_gffs


def db_main(db, accs, cpus = 1):
    """Grab a gff for accessions that may have multiple ome codes"""
    omes = set([x[:x.find('_')] for x in accs])
    db = db.set_index('ome')
    grab_acc_cmds = []
    # for each ome code, prepare a command to acquire the accessions
    for ome in list(omes):
        ome_accs = [acc for acc in accs if acc.startswith(ome + '_')]
        gff_list = gff2list(db[ome]['gff3'])
#        grab_acc_res = [grab_gff_accs(gff_list, ome_accs, ome)]
        grab_acc_cmds.append([gff_list, ome_accs, ome])

    with mp.Pool(processes = cpus) as pool:
        grab_acc_res = pool.starmap(grab_gff_accs, grab_acc_cmds)

    gffs = {}
    for res in grab_acc_res:
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

    # if there is an input file, extract the accessions from that
    if args.input:
        input_file = format_path(args.input)
        with open(input_file, 'r') as raw:
            accs = [x.rstrip().split('\t')[args.column-1] for x in raw]
    # no input and accession
    elif not args.accession:
        raise ValueError('need input file or accession')
    # parse accession input into a list of acc or accs
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

    # if no gff is provided, then acquire it from the primary database
    db_path = format_path(args.mtdb)
    if not args.gff:
        db = mtdb(format_path(args.mtdb))
        gff_lists = db_main(db, accs, cpus = args.cpu)
    # otherwise just use what is available
    else:
        gff_path = format_path(args.gff)
        gff_lists = gff_main(gff_data, accs)

    # if there is an inputted accession, then print the output to stdout
    if args.accession:
        print(list2gff(gff_lists[list(gff_lists.keys())[0]]).rstrip(), flush = True)
    # if it is specified output, open a folder for it
    elif args.ome:
        output = mkOutput(os.getcwd() + '/', 'acc2gff')
        for ome in gff_strs:
            if gff_lists[ome]:
                with open(output + ome + '.accs.gff3', 'w') as out:
                    out.write(list2gff(gff_lists[ome]))
            else:
                eprint('ERROR: ' + ome + ' failed, no accessions retrieved', flush = True)
    # print to stdout each gff_list
    else:
        out_str = ''
        for ome in gff_lists:
            if gff_lists[ome]:
                out_str += list2gff(gff_lists[ome]) + '\n'
            else:
                eprint('ERROR: ' + ome + ' does not have accession', flush = True)
        print(out_str)

    sys.exit(0)

if __name__ == '__main__':
    cli()
