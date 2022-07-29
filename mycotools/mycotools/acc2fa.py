#! /usr/bin/env python3

import os
import re
import sys
import argparse
from mycotools.lib.biotools import fa2dict, dict2fa, reverse_complement
from mycotools.lib.dbtools import mtdb, masterDB
from mycotools.lib.kontools import format_path, eprint, stdin2str


def extractHeaders( fasta_file, accessions, ome = None ):
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
   

def dbmain( db, accs ):
    '''takes in mtdb, takes accessions ome by ome from accs'''

    fa_dict = {}
    db = db.set_index( 'ome' )
    omes = set([x[:x.find('_')] for x in accs])
    for ome in omes:
        if ome:
            accessions = [x for x in accs if x.startswith(ome + '_')]
            fasta = db[ome]['faa']
            fa_dict = {**fa_dict, **extractHeaders(fasta, accessions)}
    
    return fa_dict


def famain( accs, fa, ome = None ):
    '''takes in accessions, fasta, and retrieves accessions'''

    fa_dict = {}
    fa_dict = {**fa_dict, **extractHeaders( fa, accs, ome )}
        
    return fa_dict


if __name__ == '__main__':

    parser = argparse.ArgumentParser( description = 'Inputs fasta or database and new line delimitted file of headers. ' + \
        'Outputs fasta with headers' )
    parser.add_argument( '-a', '--accession', help = '"-" for stdin. For coordinates ' + \
        'append [$START-$END] - reverse coordinates for antisense' )
    parser.add_argument( '-i', '--input', help = 'File with accessions' )
    parser.add_argument( '-f', '--fasta', help = 'Fasta input' )
    parser.add_argument( '-c', '--column', default = 1, 
        help = 'Accessions column for -i (1 indexed). DEFAULT: 1', type = int )
    parser.add_argument( '-s', '--start', help = 'Start index column (1 indexed)', type = int )
    parser.add_argument( '-e', '--end', help = 'End index column (1 indexed)', type = int )
    parser.add_argument( '-d', '--database', default = masterDB(), help = 'mycodb DEFAULT: master' )
    args = parser.parse_args()

    if args.input:
        input_file = format_path( args.input )
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
    else:
        if '-' == args.accession:
            data = stdin2str()
            accs = data.split()
        else:
            accs = [args.accession]
       
    db_path = format_path( args.database )
    if not args.fasta:
        db = mtdb( format_path(args.database) )
        fa_dict = dbmain( db, accs )
        fasta_str = dict2fa(fa_dict)
    else:
        fa_path = format_path( args.fasta )
        fa_dict = famain( accs, fa_path )
        fasta_str = dict2fa(fa_dict)

    print( fasta_str.rstrip() , flush = True)

    sys.exit( 0 )
