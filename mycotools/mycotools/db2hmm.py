#! /usr/bin/env python3

import os
import re
import sys
import argparse
import subprocess
import pandas as pd
import multiprocessing as mp
from io import StringIO
from mycotools.lib.kontools import intro, outro, collect_files, multisub, findExecs
from mycotools.lib.dbtools import db2df
from mycotools.lib.biotools import dict2fa
from mycotools.extractHmmsearch import main as exHmm
from mycotools.acc2fa import dbmain as acc2fa
from mycotools.extractHmmAcc import grabAccs, main as absHmm

def compileHmmCmd( db, hmmdb_path, output, ome_set = set() ):
    '''
    Inputs: mycotools db, hmmdb_path, output directory, set of omes to ignore.
    Outputs: tuples of arguments for hmmsearches
    '''

    cmd_tuples = [ ]
    for i, row in db.iterrows():
        if pd.isnull( row['faa'] ) or row['ome'] in ome_set:
            continue
        proteome_path = os.environ['MYCOFAA'] + '/' + row['faa']
        output_path = output + '/' + row['ome'] + '.out'
        cmd = ( 'hmmsearch', '-o', output_path, hmmdb_path, proteome_path ) 
        cmd_tuples.append( cmd )

    return cmd_tuples
    

def compileextractHmmCmd( db, args, output ):
    '''
    Inputs: mycotools db, argparse arguments, and output path
    Outputs: tuples of arguments for `runExHmm`
    '''

    cmd_tuples = [ ]
    for i, row in db.iterrows():
        if pd.isnull( row['faa'] ):
            continue
        hmmsearch_out = output + '/' + row['ome'] + '.out'
        cmd_tuples.append( ( args, hmmsearch_out, output ) )

    return cmd_tuples


def runExHmm( args, hmmsearch_out, output ):

    with open( hmmsearch_out, 'r' ) as raw:
        data = raw.read()
    ome = os.path.basename(hmmsearch_out).replace( '.out', '' )
    out_dict = None
    if len(data) > 100:
        out_dict = exHmm(data, args[0], args[1], args[2], args[3], header = False)
        out_dict = { x: [ome, out_dict[x][1]] for x in out_dict }
    else:
        print('\t' + ome + ' empty hmmsearch', flush = True)

    return out_dict
    

def compileacc2fa( db, output, q_dict ):

    cmd_tuples = []
    if not os.path.isdir(output):
        os.mkdir( output )
    fa_dict = { row['ome']: os.environ['MYCOFAA'] + '/' + row['faa'] \
        for i, row in db.iterrows() if not pd.isnull( row['faa'] ) }
    for q in q_dict:
        cmd_tuples.append( ( q_dict[q], 0, fa_dict, output + '/' + q + '.aa.fa' ) )

    return cmd_tuples


def runacc2fa( q_files, header, fa_dict, output ):

    fa_str = ''
    for ome in q_files:
        q = q_files[ome]
        df = pd.read_csv( StringIO(q), sep = '\t', header = None )
        fa_str += dict2fa(acc2fa( df, header, fa = fa_dict[ ome ] )) + '\n'
    
    with open( output, 'w' ) as out:
        out.write( fa_str )


def compileHmmalignCmd( output, accessions ):

    cmd_tuples = [ [], [] ]
    for acc in accessions:
        fa = output + '/fastas/' + acc + '.aa.fa'
        hmm = output + '/hmms/' + acc + '.hmm'
        align = output + '/aligns/' + acc + '.a2m'
        conv = output + '/aligns/' + acc + '.phylip'
        trim = output + '/trimmed/' + acc + '.trimal'
        if os.path.isfile( conv ):
            with open( conv, 'r' ) as raw:
                data = raw.read()
            if len( data ) > 10:
                continue
        hmmalign = ( 'hmmalign', '--outformat', 'A2m', '-o', align, hmm, fa )
        esl = ( 'esl-reformat', '-o', conv, 'phylip', align )
        cmd_tuples[0].append( hmmalign )
        cmd_tuples[1].append( esl )

    return cmd_tuples


def compileTrimalCmd( output, trimal = '', trimmed = None):

    trimal_args = trimal.split( ' ' )
    cmd_tuples = []
    aligns = collect_files( output + '/aligns', 'phylip' )
    if trimmed:
        trimmed = set( os.path.basename(x).replace('.trimal','') \
            for x in trimmed )
        aligns = [x for x in aligns if os.path.basename(x).replace('.phylip','') \
             not in trimmed]
    for align in aligns:
        trim = output + '/trimmed/' + \
            os.path.basename( align ).replace('.phylip','.trimal')
        args = ['trimal', '-in', align, '-out', trim]
        args.extend( trimal_args ) 
        cmd_tuples.append( tuple( args ) )

    return cmd_tuples


if __name__ == '__main__':

    mp.set_start_method('spawn')
    parser = argparse.ArgumentParser( description = 'Runs `hmmsearch` v each proteome in a' \
        + ' `.db`. For each query, extracts results and optionally outputs a compiled fasta' \
        + ', hmmalignment and/or trimmed alignment.' )
    parser.add_argument( '-i', '--input', required = True, help = 'Input database `.db`' )
    parser.add_argument( '-d', '--hmmdb', required = True, help = 'Hmm database `.hmm`' )
    parser.add_argument( '-o', '--output', help = 'User-specified output directory' )
    parser.add_argument( '-p', '--previous', help = 'Previous db2hmmsearch dir. ' + \
        'Be wary of incomplete outputs from interrupted runs' )
    parser.add_argument( '-c', '--cpu', type = int, \
        help = 'Processors to use. Default = All' )
    parser.add_argument( '-f', '--fasta', action = 'store_true', \
        help = 'Compile fastas for each query' )
    parser.add_argument( '-l', '--align', default = False, action = 'store_true', \
        help = 'Align fastas to original hmm and trim via `trimal`. Calls `-f`' )
    parser.add_argument( '-a', '--accession', action = 'store_true', default = False, \
        help = 'Extract accessions instead of queries (Pfam, etc). Requires `-f`' )
    parser.add_argument( '-b', '--best', type = int, help = '# top hits for each organism' \
        + ' Requires `-f`' )
    parser.add_argument( '-t', '--threshold', help = 'Query percent threshold (+/-). Requires `-f`' )
    parser.add_argument( '-e', '--evalue', help = 'E value threshold, e.g. ' \
        + '10^(-x) where x is the input. Requires `-f`' )
    parser.add_argument( '--trimal', default = '', help = \
        'User-specified trimAl commands in "", e.g.: "-strictplus -fasta"' )

    args = parser.parse_args()
    deps = ['hmmsearch']

    if args.previous:
        args.output = args.previous

    if not args.output:
        output = os.getcwd() + '/db2hmmsearch'
    else:
        output = os.path.abspath( args.output )
    if not os.path.isdir( output ):
        os.mkdir( output )

    if args.cpu and args.cpu < mp.cpu_count():
        cpu = args.cpu
    else:
        cpu = mp.cpu_count()

    trimal = args.trimal.replace('"','')

    if args.align:
        args.fasta = True
        deps.extend( ['hmmalign', 'trimal'] )
    else:
        trimal = None

    args_dict = { 
        'Database': args.input, 'Hmm database': args.hmmdb, 'Output': output,
        'Previous run': args.previous, 'CPUs': cpu, 'Fasta': args.fasta, 
        'Alignment': args.align, 'Accession': args.accession, 'Top hits': args.best, 
        'Threshold': args.threshold, 'E-value': args.evalue, 'trimAl args': trimal 
    }
    start_time = intro( 'db2hmmsearch', args_dict )
    findExecs( deps, exit = {deps})
    db = db2df( args.input )

#    main( args, db, output, cpu, previous = args.previous )
#def main( args, db, output, cpu, previous = False ):

    previous = args.previous
    ome_set, skip1, skip2 = set(), False, False
    if previous:
        print( '\nCompiling previous run' , flush = True)
        if os.path.isdir( output + '/fastas' ):
            fas = collect_files( output + '/fastas', 'aa.fa' )
            ranQueries = [ os.path.basename(fa).replace('.aa.fa','') for fa in fas ]
            with open( args.hmmdb, 'r' ) as hmmdb_raw:
                queries = grabAccs( hmmdb_raw.read() )
            if len( set( queries ).intersection( set( ranQueries ) ) ) == len( set( ranQueries ) ):
                skip1 = True
                for fa in fas:
                    with open( fa, 'r' ) as raw:
                        data = raw.read()
                    if len(data) < 10:
                        skip1 = False
                        break
            if skip1:
                if os.path.isdir( output + '/aligns' ):
                    phylips = collect_files( output + '/aligns', 'phylip' )
                    ranQueries = [os.path.basename(x).replace('.phylip','') \
                        for x in phylips ]
                    if len( set( queries ).intersection( set ( ranQueries ) \
                        ) ) == len( set( ranQueries ) ):
                        skip2 = True
                        for phy in phylips:
                            with open( phy, 'r' ) as raw:
                                data = raw.read()
                            if len(data) < 10:
                                skip2 = False
                                break
                
        if not skip1:
            if os.path.isfile( output + '/omes.tar.gz' ):
                if not os.path.isdir( output + '/omes' ):
                    tar_code = subprocess.call( 'tar -xzf ' + output + \
                        '/omes.tar.gz', stdout = subprocess.DEVNULL, \
                        stderr = subprocess.DEVNULL, shell = True )
            omes = collect_files( os.path.abspath( previous ) + '/omes', 'out' )
            ome_set = set( os.path.basename(x).replace('.out', '') for x in omes )
        else:
            print('\thmmsearch -> hits.fa DONE', flush = True)
            if skip2:
                print('\thits.fa -> phylip DONE', flush = True)

    if not skip1:
        print( '\nSubmitting omes to hmmsearch' , flush = True)
        if not os.path.isdir( output + '/omes' ):
            os.mkdir( output + '/omes' )
        hmmsearch_tuples = compileHmmCmd( db, args.hmmdb, output + '/omes', ome_set = ome_set )
        hmmsearch_codes = multisub( hmmsearch_tuples, processes = cpu - 1 )
        for code in hmmsearch_codes:
            if code['code'] != 0:
                print( '\tExit ' + str(code['code'], flush = True) + ' ' + \
                     os.path.basename( code['stdin'][-1] ) )
    
        print( '\nExtracting hmmsearch output' , flush = True)
        exHmm_args = [args.accession, args.best, args.threshold, args.evalue]
        exHmm_tuples = compileextractHmmCmd( db, exHmm_args, output + '/omes' )
        with mp.get_context( 'spawn' ).Pool(processes = cpu) as pool:
            hmmAligns = pool.starmap( runExHmm, exHmm_tuples )

        if args.fasta:
            subprocess.Popen(['tar', '-czf', output + "/omes.tar.gz", \
                output + '/omes', '--remove-files'], stdout = subprocess.DEVNULL, \
                stderr = subprocess.DEVNULL )
            print( '\nCompiling fastas' , flush = True)
            q_dict = {}
            for hit in hmmAligns:
                if hit:
                    for query in hit:
                        if query not in q_dict:
                            q_dict[query] = {}
                        ome = hit[query][0]
                        align = hit[query][1]
                        q_dict[query][ ome ] = align

            acc2fa_tuples = compileacc2fa( db, output + '/fastas', q_dict )
            with mp.get_context( 'spawn').Pool( processes = cpu ) as pool:
                pool.starmap( runacc2fa, acc2fa_tuples )        

    if args.align:
        if not skip2:
            print( '\nAligning fastas to hmms' , flush = True)
            if not os.path.isdir( output + '/hmms' ):
                os.mkdir( output + '/hmms' )
            with open( args.hmmdb, 'r' ) as hmmdb:
                hmm_dict = absHmm( hmmdb.read() )
            for accession in hmm_dict:
                with open( output + '/hmms/' + accession + '.hmm', 'w' ) as out:
                    out.write( hmm_dict[ accession ] )

            if not os.path.isdir( output + '/aligns' ):
                os.mkdir( output + '/aligns' )

            
            hmmAlign_tuples = compileHmmalignCmd( output, list(hmm_dict.keys() ) )
            hmmAlign_codes = multisub( hmmAlign_tuples[0], processes = cpu - 1 )
            esl_codes = multisub( hmmAlign_tuples[1], processes = cpu - 1 )
            for index in range( len (hmmAlign_codes) ):
                code = hmmAlign_codes[ index ]
                if code['code'] != 0:
                    print( '\thmmalign exit ' + str(code['code'], flush = True) + ' ' + \
                         os.path.basename( code['stdin'][-1] ) )
                else:
                    code = esl_codes[ index ]
                    if code['code'] != 0:
                        print( '\tesl-reformat exit ' + str(code['code'], flush = True) + ' ' + \
                            os.path.basename( code['stdin'][-1] ) )
                    os.remove( code['stdin'][-1] )
        
        print( '\nTrimming alignments' , flush = True)
        if not os.path.isdir( output + '/trimmed' ):
            os.mkdir( output + '/trimmed' )
        trimmed = None
        if args.previous:
            trimmed = collect_files( output + '/trimmed', 'trimal' )
        trimal_tuples = compileTrimalCmd( output, trimal = trimal, trimmed = trimmed )

        trimal_codes = multisub( trimal_tuples, processes = cpu - 1 )
        for code in trimal_codes:
            if code['code'] != 0:
                print( '\tExit ' + str(code['code'], flush = True) + ' ' + \
                     os.path.basename( code['stdin'][2] ) )

    outro( start_time )
