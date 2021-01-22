#! /usr/bin/env python3

import pandas as pd
import os, sys, argparse, re
from mycotools.lib.fastatools import fasta2dict, dict2fasta
from mycotools.lib.dbtools import db2df, masterDB
from mycotools.lib.kontools import formatPath

def extractHeaders( fasta_file, accessions ):

    fasta = fasta2dict(fasta_file)
    out_fasta = {}

    if type(accessions) == dict:
        for header in accessions:
            start = int(accessions[ header ][ 0 ])
            end = int(accessions[ header ][ 1 ])
            header = re.sub( r'\$\$\d+$', '', header )
            if start - 1 < 0:
                start = 1
            new_header = header + '[' + str(start) + '-' + str(end) + ']' 
            out_fasta[ new_header ] = { 'sequence': fasta[ header ]['sequence'][ start - 1 : end ], 
                'description': '' }
    else:
        for header in accessions:
            out_fasta[ header ] = fasta[ header ]

    out_fasta_str = dict2fasta( out_fasta )

    return out_fasta_str
    

def main( df, column = None, db = None, fa = None, start = None, end = None ):

    acc_comp = re.compile( r'(.*?)\[(\d+)\-(\d+)\]$' )
    if start or end:
        if not start and not end:
            print( '\n`Start` column or `stop` column required together.\n' )
            sys.exit( 4 )
        new_df = pd.DataFrame()
        for i, row in df.iterrows():
            count = 0
            name = str(row[ column ])
            if len(new_df) > 0:
                while row[column] in set(new_df[ column ]):
                    row[ column ] = name + '$$' + str(count)
                    count += 1
            new_df = new_df.append( row )
        new_df = new_df.set_index( column )

    if fa: 
        if start:
            accessions_prep = list( new_df.index )
            accessions = {}
            for accession in accessions_prep:
                accessions[ accession ] = [ new_df[start][accession], \
                    new_df[end][accession], accession ]
        else:
            if type(df) is not str:
                accessions = list(df[column])
            else:
                accessions = [ df ]
        fasta_str = extractHeaders( fa, accessions )
    else:
        fasta_str = ''
        db = db.set_index( 'internal_ome' )
        if type(df) is not str:
            omes_prep = re.findall( r'^.*?_', '\n'.join(list(df[ column ])), re.M )
            omes = set( x[:-1] for x in omes_prep )
        else:
            omes_prep = re.search( r'(^.*?)_', df )
            if omes_prep is None:
                print('\nERROR: no ome found in accession and no fasta specified')
                sys.exit( 5 )
            omes = [ omes_prep[1] ]

        if start:
            for ome in list(omes):
                accessions_prep = [ x for x in new_df.index if x.startswith( ome + '_' ) ]
                accessions = {}
                for accession in accessions_prep:
                    accessions[ accession ] = [ int(new_df[args.start][accession]), 
                        int(new_df[args.end][accession]), accession ]
                fasta = os.environ['MYCOFAA'] + '/' + db['proteome'][ome]
                fasta_str += extractHeaders( fasta, accessions ) + '\n'
        else:
            for ome in list(omes):
                if type(df) is not str:
                    accessions = [ x for x in df[ column ] if x.startswith( ome + '_' ) ]
                else:
                    accessions = [ df ]
#                acc_search = acc_comp.search( accessions_prep[0] )
 #               if acc_search:
  #                  accessions = {}
   #                 for accession in accessions_prep:
    #                    accessions[ acc_search[1] ] = [ 
     #                       int(acc_search[2]),
      #                      int(acc_search[3]),
       #                     acc_search[1]
        #                    ]
         #       else:
          #          accessions = accessions_prep
                fasta = os.environ['MYCOFAA'] + '/' + db['proteome'][ome]
                fasta_str += extractHeaders( fasta, accessions ) + '\n'

    return fasta_str


if __name__ == '__main__':

    parser = argparse.ArgumentParser( description = 'Inputs fasta or database and new line delimitted file of headers. ' + \
        'Outputs fasta with headers' )
    parser.add_argument( '-i', '--input', help = 'Input file with accessions' )
    parser.add_argument( '-a', '--accession', help = 'Input accession' )
    parser.add_argument( '-f', '--fasta', help = 'Fasta input' )
    parser.add_argument( '-c', '--column', default = 0, help = 'Column to abstract from. Must be specified if there ' \
        + 'are headers. Otherwise, the default is the first tab separated column' )
    parser.add_argument( '-s', '--start', help = 'Optional start index column. Will subtract 1 from this (python is ' \
        + 'a 0 based language)' )
    parser.add_argument( '-e', '--end', help = 'Optional end index column' )
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
            fasta_str = main( args.accession, db = db )
        else:
            fasta_str = main( df, args.column, db = db, start = args.start, end = args.end )
    else:
        fa_path = formatPath( args.fasta )
        if args.accession:
            fasta_str = main( args.accession, fa = fa_path )
        else:
            fasta_str = main( df, args.column, fa = fa_path, start = args.start, end = args.end )

    if not args.accession:
        with open( input_file + '.fa', 'w' ) as out:
            out.write( fasta_str )
    else:
        print( fasta_str.rstrip() )

    sys.exit( 0 )
