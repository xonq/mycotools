#! /usr/bin/env python3

import os
import re
import sys
import argparse
from mycotools.lib.biotools import fa2dict, dict2fa, reverse_complement
from mycotools.lib.dbtools import mtdb, masterDB
from mycotools.lib.kontools import format_path, eprint, stdin2str

def extract_mtdb_accs(fa_dict, accs):
    """extract MTDB accessions from fa_dict using a list `accs`; can 
    extract coordinates of an accession if labeled 
    `acc[START:END] | acc[START-END]`. If start > end, extract the reverse
    complement of the sequence. Returns a fa_dict of the results"""

    out_fa = {}
    for acc in accs:
        coords = acc[acc.find('['):acc.find(']')]
        if coords: # if there are coordinates
            try: # try splitting by colon
                start, end = [int(x) for x in coords.split(':')]
            except ValueError: # try splitting by hyphen
                try:
                    start, end = [int(x) for x in coords.split('-')]
                except ValueError: # invalid coordinates, try whole acc
                    out_fa[acc] = fa_dict[acc]
                    continue
            acc_name = acc[:start]
            if start < end:
                out_fa[acc] = {
                    'sequence': fa_dict[acc_name]['sequence'][start:end+1],
                    'description': fa_dict[acc_name]['description']
                    }
            else:
                out_fa[acc] = {
                    'sequence': reverse_complement(
                        fa_dict[acc_name]['sequence'][end-1:start+1]
                        ),
                    'description': fa_dict[acc_name]['description']
                    }
        else: # no coordinates
            out_fa[acc] = fa_dict[acc] # extract the whole accession

    return out_fa
            

def extractHeaders(fasta_file, accessions, ome = None):
    '''searches headers for "[]", which indicate coordinate-based extraction.
    otherwise, just retrieves the accession from the fasta dictionary'''

    fasta = fa2dict(fasta_file)
    out_fasta = {}
    acc_comp = re.compile( r'([^/[]*)\[(\d+)\-(\d+)\]$' )

    for header in accessions:
        if '[' in header: # if coordinate based
            coord_search = acc_comp.search( header )
            if coord_search is not None:
                start = int(coord_search[2]) - 1
                end = int(coord_search[3])
                acc = coord_search[1]
                if ome:
                    header = ome + '_' + header
                if end > start:
                    out_fasta[header] = {
                        'sequence': fasta[acc]['sequence'][start:end],
                        'description': ''
                        }
                else:
                    out_fasta[header] = {
                        'sequence': reverse_complement(
                            fasta[acc]['sequence'][end-1:start+1]
                            ),
                        'description': ''
                        }
            else:
                if ome:
                    out_fasta[ome + '_' + header] = fasta[header]
                else:
                    out_fasta[header] = fasta[header]
        else:
            if ome:
                out_fasta[ome + '_' + header] = fasta[header]
            else:
                out_fasta[header] = fasta[header]

    return out_fasta
   

def dbmain(db, accs):
    '''takes in mtdb, takes accessions ome by ome from accs'''

    fa_dict = {}
    db = db.set_index( 'ome' )
    omes = set([x[:x.find('_')] for x in accs])
    for ome in omes:
        if ome:
            ome_accs = [x for x in accs if x.startswith(ome + '_')]
            ome_fasta = fa2dict(db[ome]['faa'])
            fa_dict = {**fa_dict, **extract_mtdb_accs(ome_fasta, ome_accs)}
    
    return fa_dict


def famain( accs, fa, ome = None ):
    '''takes in accessions, fasta, and retrieves accessions'''

    fa_dict = {}
    fa_dict = {**fa_dict, **extractHeaders(fa, accs, ome)}
        
    return fa_dict


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description = 'Inputs accession, extracts fasta')
    parser.add_argument('-a', '--accession', help = '"-" for stdin. For coordinates ' + \
        'append [$START-$END] - reverse coordinates for antisense')
    parser.add_argument('-i', '--input', help = 'File with accessions')
    parser.add_argument('-f', '--fasta', help = 'Fasta input')
    parser.add_argument('-c', '--column', default = 1, 
        help = 'Accessions column for -i (1 indexed). DEFAULT: 1', type = int)
    parser.add_argument('-s', '--start', help = 'Start index column (1 indexed)', type = int)
    parser.add_argument('-e', '--end', help = 'End index column (1 indexed)', type = int)
    parser.add_argument('-d', '--database', default = masterDB(), help = 'MTDB DEFAULT: master')
    args = parser.parse_args()

    if args.input: # input file
        input_file = format_path(args.input)
        if args.start: # if using columns for coordinates
            with open(input_file, 'r') as raw:
                accs = []
                for line in raw:
                    if not line.startswith('#'):
                        d = line.rstrip().split('\t')
                        accs.append([
                            d[args.column-1] + '[' + d[args.start-1] + '-' + d[args.end-1] + ']'
                            ])
        else:
            with open(input_file, 'r') as raw:
                accs = [x.rstrip().split('\t')[args.column-1] for x in raw]
    else: # assume we are using accessions
        if '-' == args.accession: # stdin
            data = stdin2str()
            accs = data.split()
        else:
            accs = [args.accession]
       
    db_path = format_path(args.database)
    if not args.fasta: # MTDB run
        db = mtdb(db_path)
        fa_dict = dbmain(db, accs)
        fasta_str = dict2fa(fa_dict)
    else: # non MTDB
        fa_path = format_path(args.fasta)
        fa_dict = famain(accs, fa_path)
        fasta_str = dict2fa(fa_dict)

    print(fasta_str.rstrip() , flush = True)
    sys.exit(0)
