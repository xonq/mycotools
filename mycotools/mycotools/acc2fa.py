#! /usr/bin/env python3

import os
import re
import sys
import argparse
from collections import defaultdict
from mycotools.lib.biotools import fa2dict, dict2fa, reverse_complement
from mycotools.lib.dbtools import mtdb, primaryDB
from mycotools.lib.kontools import format_path, eprint, stdin2str

def extract_mtdb_accs_exp(fa_dict, accs):
    return {acc: fa_dict[acc] for acc in accs}

def extract_mtdb_accs(fa_dict, accs, spacer = ''):
    """extract MTDB accessions from fa_dict using a list `accs`; can 
    extract coordinates of an accession if labeled 
    `acc[START:END] | acc[START-END]`. If start > end, extract the reverse
    complement of the sequence. Returns a fa_dict of the results"""

    out_fa = {}
    for acc in accs:
        coords = acc[acc.find('[')+1:acc.find(']')]
        if coords: # if there are coordinates
            try: # try splitting by colon
                start, end = [int(x) for x in coords.split(':')]
            except ValueError: # try splitting by hyphen
                try:
                    start, end = [int(x) for x in coords.split('-')]
                except ValueError: # invalid coordinates, try whole acc
                    out_fa[acc] = fa_dict[acc]
                    continue
            acc_name = acc[:acc.find('[')]
            if start < end:
                try:
                    out_fa[acc] = {
                        'sequence': fa_dict[acc_name]['sequence'][start:end+1],
                        'description': fa_dict[acc_name]['description']
                        }
                except KeyError:
                    eprint(spacer + 'WARNING: invalid accession ' + acc_name, 
                           flush = True)
            else:
                try:
                    out_fa[acc] = {
                        'sequence': reverse_complement(
                        fa_dict[acc_name]['sequence'][end-1:start+1]
                        ),
                        'description': fa_dict[acc_name]['description']
                        }
                except KeyError:
                    eprint(spacer + 'WARNING: invalid accession ' + acc_name,
                           flush = True)
        else: # no coordinates
            try:
                out_fa[acc] = fa_dict[acc] # extract the whole accession
            except KeyError:
                eprint(spacer + 'WARNING: invalid accession ' + acc_name,
                       flush = True)

    return out_fa

def extractHeaders(fasta_file, accessions, ome = None):
    """searches headers for "[]", which indicate coordinate-based extraction.
    otherwise, just retrieves the accession from the fasta dictionary"""

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

def dbmain(db, accs, error = True, 
           spacer = '\t\t\t', coord_check = True):
    """takes in mtdb, takes accessions ome by ome from accs;
    error will exit if set to true, report to stderr if not;
    coord_check will check entries for coordinate signatures in accession names
    but will bypass if false (throughput increase). Returns a fa_dict of the
    accessions and coordinates if applicable"""

    # set the db index
    db = db.set_index( 'ome' )
    # create a dict to populate for each ome
    ome_data = defaultdict(list)
    # organize accession by ome
    for acc in accs:
        ome = acc[:acc.find('_')]
        ome_data[ome].append(acc)

    # populate the fasta dict with the accession information by referencing the
    # mycotools proteome entry for each ome
    fa_dict = {}
    if coord_check:
        for ome, ome_accs in ome_data.items():
            try:
                ome_fasta = fa2dict(db[ome]['faa'])
            except KeyError: # if there is a missing ome
                if error:
                    raise KeyError('invalid ome: ' + ome)
                else:
                    eprint(spacer + ome + ' not in database', flush = True)
            fa_dict = {**fa_dict, **extract_mtdb_accs(ome_fasta, ome_accs)}
    else:
        for ome, ome_accs in ome_data.items():
            try:
                ome_fasta = fa2dict(db[ome]['faa'])
            except KeyError: # if there is a missing ome
                if error:
                    raise KeyError
                else:
                    eprint(spacer + ome + ' not in database', flush = True)
            fa_dict = {**fa_dict, **extract_mtdb_accs_exp(ome_fasta, ome_accs)}

    return fa_dict


def famain(accs, fa, ome = None):
    """takes in accessions, fasta, and retrieves accessions"""

    fa_dict = {}
    fa_dict = {**fa_dict, **extractHeaders(fa, accs, ome)}
        
    return fa_dict


def cli():

    parser = argparse.ArgumentParser(description = 'Inputs accession, extracts fasta')
    parser.add_argument('-a', '--accession', help = '"-" for stdin. For coordinates ' + \
        'append [$START-$END] - reverse coordinates for antisense')
    parser.add_argument('-i', '--input', help = 'File with accessions')
    parser.add_argument('-f', '--fasta', help = 'Fasta input')
    parser.add_argument('-c', '--column', default = 1, 
        help = 'Accessions column for -i (1 indexed). DEFAULT: 1', type = int)
    parser.add_argument('-s', '--start', help = 'Start index column (1 indexed)', type = int)
    parser.add_argument('-e', '--end', help = 'End index column (1 indexed)', type = int)
    parser.add_argument('-d', '--mtdb', default = primaryDB())
    args = parser.parse_args()

    if args.input: # input file
        input_file = format_path(args.input)
        if args.start: # if using columns for coordinates
            with open(input_file, 'r') as raw:
                accs = []
                for line in raw:
                    if not line.startswith('#') and line.rstrip():
                        d = line.rstrip().split('\t')
                        accs.append([
                            d[args.column-1] + '[' + d[args.start-1] + '-' + d[args.end-1] + ']'
                            ])
        else:
            with open(input_file, 'r') as raw:
                accs = [x.rstrip().split('\t')[args.column-1] \
                        for x in raw if x.rstrip()]
    else: # assume we are using accessions
        if '-' == args.accession: # stdin
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
       
    db_path = format_path(args.mtdb)
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

if __name__ == '__main__':
    cli()
