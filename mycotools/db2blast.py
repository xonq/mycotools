#! /usr/bin/env python3

import argparse, subprocess, os, datetime, sys, multiprocessing as mp, pandas as pd, numpy as np
from mycotools.lib.kontools import eprint, intro, outro, formatPath, multisub, collect_files
from mycotools.lib.dbtools import db2df, masterDB
from mycotools.acc2fa import famain as acc2fa


def search4db( ome, blast_dir, extension, seq_dir ):

    subject = blast_dir + ome
    blast_file = blast_dir + ome + extension
    if os.path.isfile( blast_file ):
        return subject
    else:
        return False


def compBlastDB( db, seq_type ):

    vsearch4db = np.vectorize(search4db)
    if seq_type == 'prot':
        blast_dir = formatPath( '$MYCOFAA/blastdb' )
        seq_dir = os.environ['MYCOFAA'] + '/'
        db = db.assign( **{seq_type: vsearch4db( db['internal_ome'], blast_dir, '.psd', seq_dir ) } )
#        db[seq_type] = db.apply(search4db, args(blast_dir, '.psd', seq_dir))
    else:
        db = db.assign( **{seq_type: vsearch4db( db['internal_ome'], blast_dir, '.psd', seq_dir ) } )
        blast_dir = formatPath( '$MYCOFNA/blastdb' )
        seq_dir = os.environ['MYCOFNA'] + '/'
#        db[seq_type] = db.apply(search4db, args(blast_dir, '.nsd', seq_dir))

    blast_db = db.loc[db[seq_type] != False]
    seq_db = db.loc[db[seq_type] == False]

    return blast_db, seq_db


def compileBlastDBCmd( ome, out_dir, env_dir, blast_scaf ):
    return blast_scaf + ['-out', out_dir + ome + '.tsv', \
        '-db', env_dir + 'blastdb/' + ome]
def compileBlastCmd( ome, biofile, env_dir, out_dir, blast_scaf ):
    return blast_scaf + ['-out', out_dir + ome + '.tsv', \
        '-subject', env_dir + biofile]

def compBlastTups( 
    blast_db, seq_db, blast_type, seq_type, out_dir, 
    biotype, env_dir, query, hsps = 1, max_hits = None,
    evalue = 10
    ):

    blast_scaf = [
        blast_type, '-query', query, '-max_hsps',
        str(hsps), '-outfmt', '6', '-evalue', str(evalue)
        ]
    if max_hits:
        blast_scaf.extend([ '-max_target_seqs', str(max_hits) ])
#    vcomp_blast_db = np.vectorize(compileBlastDBCmd)
 #   vcomp_seq_db = np.vectorize(compileBlastCmd)
  #  blast_db['cmd'] = vcomp_blast_db(
   #     blast_db['internal_ome'], out_dir, env_dir, blast_scaf
    #    )
#    seq_db['cmd'] = vcomp_seq_db(
 #       seq_db['internal_ome'], seq_db['biotype'], env_dir,
  #      out_dir, blast_scaf
   #     )

    blast_cmds = []
    for i, row in blast_db.iterrows():
        blast_cmds.append( compileBlastDBCmd(
            row['internal_ome'], out_dir, env_dir, blast_scaf
            ) )
    for i, row in seq_db.iterrows():
        if not pd.isnull(row[biotype]):
            blast_cmds.append( compileBlastCmd(
                row['internal_ome'], row[biotype], env_dir,
                out_dir, blast_scaf
                ) )

    return blast_cmds


def parseOutput( ome, file_, bitscore = 0, pident = 0 ):

    ome_results = [ome, []]
    if os.path.exists( file_ ):
        with open( file_, 'r' ) as raw:
            raw_data = raw.read()
        data = [x.split('\t') for x in raw_data.split('\n') if x]
        for i in data:
            if int(float(i[-1])) < bitscore and int(1000*float(i[2])) > 100000 * pident:
                ome_results[-1].append(i)

    return ome_results


def compileResults( res_dict ):

    output_res = {}
    for i in res_dict:
        for hit in res_dict[i]:
            query = hit[0]
            subject = hit[1]
            pident = hit[2]
            start = hit[-4]
            end = hit[-3]
            if query not in output_res:
                output_res[query] = {}
            if i not in output_res[query]:
                output_res[query][i] = []

            output_res[query][i].append( [ subject, start, end ] )

    return output_res


def compAcc2fa( db, biotype, env_dir, output_res, subhit = False ):

    acc2fa_cmds, db, subs = {}, db.set_index('internal_ome'), set()
    if subhit:
        for query in output_res:
            acc2fa_cmds[query] = []
            for i in output_res[query]:
                temp_df, accs = pd.DataFrame(), []
                for hit in output_res[query][i]:
                    accs.append(hit[0] + '[' + hit[1] + '-' + x[2] + ']')
                temp_df[0] = accs
                acc2fa_cmds[query].append(
                    [temp_df, env_dir + db[biotype][i], 0]
                    )
    else:
        for query in output_res:
            acc2fa_cmds[query] = []
            for i in output_res[query]:
                temp_df, accs = pd.DataFrame(), []
                for hit in output_res[query][i]:
                    subject = hit[0]
                    if subject not in subs:
                        subs.add( subject )
                        accs.append(subject)
                temp_df[0] = accs
                if accs:
                    acc2fa_cmds[query].append( [
                        temp_df, env_dir + db[biotype][i], 0
                        ] )

    return acc2fa_cmds


def main( 
    db, blast, query, out_dir, hsps = 1, 
    max_hits = None, evalue = 10, bitscore = 0, 
    pident = 0, #coverage = 0, 
    prev = False, cpus = 1 
    ):

    out_dir, query = formatPath(out_dir), formatPath( query )
    if not out_dir.endswith('/'):
        out_dir += '/'
    if not os.path.isdir( out_dir ):
        os.mkdir( out_dir )
    finished, rundb = set(), db
    if not os.path.isdir( out_dir + 'reports/' ):
        os.mkdir( out_dir + 'reports/' )
    elif prev:
        reports = collect_files(out_dir + 'reports/', 'tsv')
        finished = set(os.path.basename(x)[:-4] for x in reports)
        rundb = db.loc[~db['internal_ome'].isin(finished)]

    if blast in {'tblastn', 'blastp'}:
        seq_type = 'prot'
        biotype = 'proteome'
        env_dir = os.environ['MYCOFAA'] + '/'
    elif blast in {'blastx', 'blastn'}:
        seq_type = 'nucl'
        biotype = 'assembly'
        env_dir = os.environ['MYCOFNA'] + '/'
    else:
        eprint('\nERROR: invalid blast type: ' + blast)

# NEED to add other threshold options
    print('\nCompiling blast databases')
    blast_db, seq_db = compBlastDB( rundb, seq_type )
    blast_tups = compBlastTups(
        blast_db, seq_db, blast, seq_type, 
        out_dir + 'reports/', biotype, env_dir, query, 
        hsps = hsps, max_hits = max_hits, evalue = evalue
        )

    print('\nBlasting')
    blast_outs = multisub( blast_tups, processes = cpus )
   
    omes = set(blast_db['internal_ome'])
    parse_tups = []
    for i, row in db.iterrows():
        if row['internal_ome'] in omes or not pd.isnull(row[biotype]):
            parse_tups.append( [row['internal_ome'],
                out_dir + 'reports/' + row['internal_ome'] + '.tsv',
                bitscore, pident ] )

    print('\nParsing reports')
    with mp.get_context( 'spawn' ).Pool( processes = cpus ) as pool:
        results = pool.starmap( parseOutput, tuple(parse_tups) )
    results_dict = {x[0]: x[1] for x in results}

    print('\nCompiling fastas')
    output_res = compileResults( results_dict )
    acc2fa_cmds = compAcc2fa( db, biotype, env_dir, output_res )
    output_fas = {}
    for query in acc2fa_cmds:
        with mp.get_context('spawn').Pool(processes = cpus) as pool:
            results = pool.starmap( acc2fa, acc2fa_cmds[query] )
        output_fas[query] = '\n'.join(results)

    return output_fas


if __name__ == '__main__':
# NEED abstract covered portion
# NEED max hits after blast compiled
# NEED hsp option
    mp.set_start_method('spawn')
    parser = argparse.ArgumentParser( 
        description = 'Blasts a query against a db and compiles output fastas ' + \
            'for each query.'
        )
    parser.add_argument( '-b', '--blast', required = True, 
        help = 'Blast type { blastn, blastp, tblastn, blastx }' )
    parser.add_argument( '-q', '--query', help = 'Query fasta' )
    parser.add_argument( '-e', '--evalue', default = -1, type = int,
        help = 'Negative e-value order max, e.g. 2 for 10^-2' )
    parser.add_argument( '-s', '--bitscore', default = 0,
        type = int, help = 'Bit score minimum' )
    parser.add_argument( '-i', '--identity', default = 0,
        type = float, help = 'Identity minimum, e.g. 0.6' )
    #parser.add_argument( '-c', '--coverage', type = float, help = 'Query coverage +/-, e.g. 0.5' )
    parser.add_argument( '-m', '--maxhits', type = int, help = 'Max hits for each organism' )
    parser.add_argument( '-d', '--database', default = masterDB(), help = 'mycodb, DEFAULT: masterdb' )
    parser.add_argument( '-o', '--output', default = os.getcwd() )
    parser.add_argument( '-p', '--previous', help = 'Previous run directory' )
    parser.add_argument( '--cpu', default = mp.cpu_count(), type = int, help = 'DEFAULT: all' )
    args = parser.parse_args()

    db_path = formatPath( args.database )
    if args.cpu > mp.cpu_count():
        cpus = mp.cpu_count()
    else:
        cpus = args.cpu
   
    if args.previous:
        output_dir = formatPath( args.previous )
        if not output_dir.endswith('/'):
            output_dir += '/'
        
    evalue = 10**( -args.evalue)
    start_args = {
        'Blast': args.blast, 'query': formatPath(args.query), 
        'database': db_path, 'output': args.output, 'evalue': evalue,
        'bitscore': args.bitscore, 'max hits/organism': args.maxhits, #'coverage': args.coverage, 
        'identity': args.identity, 'cores': cpus, 'previous': args.previous,
        }

    start_time = intro( 'db2blast', start_args )
    date = start_time.strftime('%Y%m%d')
   
    if not args.previous:
        output_dir = formatPath( args.output ) + date + '_db2blast/'

    db = db2df(db_path)

    output_fas = main( 
        db, args.blast, args.query, output_dir, 
        bitscore = args.bitscore, pident = args.identity,
        max_hits = args.maxhits, prev = args.previous, cpus = cpus 
        )
    if not os.path.isdir( output_dir + 'fastas/' ):
        os.mkdir( output_dir + 'fastas/' )

    for query in output_fas:
        with open( output_dir + 'fastas/' + query + '.fa', 'w' ) as out:
            out.write( output_fas[query] )

    outro( start_time )
