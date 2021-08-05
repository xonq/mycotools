#! /usr/bin/env python3

import pandas as pd
import os, sys, argparse, re
from mycotools.lib.fastatools import fasta2dict, dict2fasta, reverse_complement
from mycotools.lib.dbtools import db2df, masterDB
from mycotools.lib.kontools import formatPath, eprint


def extractHeaders( fasta_file, accessions, ome = None ):

    fasta = fasta2dict(fasta_file)
    out_fasta = {}
    acc_comp = re.compile( r'([^/[]*)\[(\d+)\-(\d+)\]$' )

    for header in accessions:
        if '[' in header:
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
   

def dbmain( db, accs, column = None, start = None, end = None ):

    fasta_str = ''
    db = db.set_index( 'internal_ome' )
    if not isinstance(accs, str):
        omes_prep = re.findall( r'^.*?_', '\n'.join(list(accs[ column ])), re.M )
        omes = set( x[:-1] for x in omes_prep )
        accs = accs.set_index( column )
        for ome in list(omes):
            accessions = [ x for x in accs.index if x.startswith( ome + '_' ) ]
            if start:
                accessions = [
                    x + '[' + accs[start][x] + '-' + accs[end][x] + ']' \
                    for x in accessions
                    ]
            fasta = os.environ['MYCOFAA'] + '/' + db['proteome'][ome]
            fasta_str += dict2fasta(extractHeaders( fasta, accessions )) + '\n'
    else:
        omes_prep = re.search( r'(^.*?)_', accs )
        if omes_prep is None:
            eprint('\nERROR: invalid accession for database search', flush = True)
            sys.exit( 5 )
        fasta = os.environ['MYCOFAA'] + '/' + db['proteome'][omes_prep[1]]
        fasta_str += dict2fasta(extractHeaders( fasta, [accs] ))

    return fasta_str


def famain( accs, fa, column = None, ome = None, start = None, end = None ):

    fasta_str = ''
    if type(accs) is not str:
        accs = accs.set_index( column )
        accessions = list(accs.index)
        if start:
            accessions = [
                x + '[' + accs[start][x] + '-' + accs[end][x] + ']' \
                for x in accessions
                ]
        fasta_str += dict2fasta(extractHeaders( fa, accessions, ome )) + '\n'
    else:
        fasta_str += dict2fasta(extractHeaders( fa, [accs], ome ))
        
    return fasta_str.rstrip()


def main( accs_str_or_df, column = None, db = None, fa = None, start = None, end = None ):

    if db and not fa:
        return dbmain(db, accs_str_or_df, column = column, start = start, end = end)
    elif not db and fa:
        return famain(accs_str_or_df, fa, column, start, end)


if __name__ == '__main__':

    parser = argparse.ArgumentParser( description = 'Inputs fasta or database and new line delimitted file of headers. ' + \
        'Outputs fasta with headers' )
    parser.add_argument( '-i', '--input', help = 'Input file with accessions' )
    parser.add_argument( '-a', '--accession', help = 'Input accession. For coordinates ' + \
        'append [$START-$END] - accepts reverse coordinates for nucleotide accessions' )
    parser.add_argument( '-f', '--fasta', help = 'Fasta input' )
    parser.add_argument( '-c', '--column', default = 0, help = 'Accessions column. Must be specified if there ' \
        + 'are headers. Otherwise, the default is the first tab separated column' )
    parser.add_argument( '-s', '--start', help = 'Start index column. Will subtract 1 from this (python is ' \
        + 'a 0 based language)' )
    parser.add_argument( '-e', '--end', help = 'End index column' )
    parser.add_argument( '-d', '--database', default = masterDB(), help = 'mycodb DEFAULT: master' )
    args = parser.parse_args()

    if args.input:
        input_file = formatPath( args.input )
        if args.column == 0:
            df = pd.read_csv(input_file, sep = '\t', header = None)
        else:
            df = pd.read_csv(input_file, sep = '\t')

    db_path = formatPath( args.database )

    if not args.fasta:
        db = db2df( formatPath(args.database) )
        if args.accession:
            fasta_str = dbmain( db, args.accession )
        else:
            fasta_str = dbmain( db, df, column = args.column, start = args.start, end = args.end )

    else:
        fa_path = formatPath( args.fasta )
        if args.accession:
            fasta_str = famain( args.accession, fa_path )
        else:
            fasta_str = famain( df, fa_path, column = args.column, start = args.start, end = args.end )

    if not args.accession:
        with open( input_file + '.fa', 'w' ) as out:
            out.write( fasta_str )
    else:
        print( fasta_str.rstrip() , flush = True)

    sys.exit( 0 )
