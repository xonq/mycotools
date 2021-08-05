#! /usr/bin/env python3

import argparse, subprocess, os, datetime, sys, re
import multiprocessing as mp, pandas as pd, numpy as np
from mycotools.lib.kontools import eprint, intro, outro, formatPath, multisub, collect_files, findExecs
from mycotools.lib.dbtools import db2df, masterDB
from mycotools.acc2fa import famain as acc2fa


def compileBlastCmd( ome, biofile, env_dir, out_dir, blast_scaf ):
    return blast_scaf + ['-out', out_dir + ome + '.tsv', \
        '-subject', env_dir + biofile]

def compBlastTups( 
    seq_db, blast_type, seq_type, out_dir, 
    biotype, env_dir, query, hsps = None, max_hits = None,
    evalue = None
    ):

    blast_scaf = [
        blast_type, '-query', query,
        '-outfmt', '6'
        ]
    if evalue:
        blast_scaf.extend(['-evalue', str(evalue)])
    if hsps:
        blast_scaf.extend(['-max_hsps', str(hsps)])
    if max_hits:
        blast_scaf.extend([ '-max_target_seqs', str(max_hits) ])

    blast_cmds = []
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
            if int(float(i[-1])) > bitscore and int(1000*float(i[2])) > 100000 * pident:
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
    if subhit or biotype == 'assembly':
        for query in output_res:
            acc2fa_cmds[query] = []
            for i in output_res[query]:
                temp_df, accs = pd.DataFrame(), []
                for hit in output_res[query][i]:
                    accs.append(hit[0] + '[' + hit[1] + '-' + hit[2] + ']')
                temp_df[0] = accs
                if biotype == 'assembly':
                    acc2fa_cmds[query].append(
                        [temp_df, env_dir + db[biotype][i], 0, i]
                        )
                else:
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
                    try:
                        acc2fa_cmds[query].append( [
                            temp_df, env_dir + db[biotype][i], 0
                            ] )
                    except KeyError:
                        eprint('\t' + i + ' not in db')

    return acc2fa_cmds


def prepOutput(out_dir, query):

    out_dir, query = formatPath(out_dir), formatPath( query )
    if not out_dir.endswith('/'):
        out_dir += '/'
    if not os.path.isdir( out_dir ):
        os.mkdir( out_dir )
    report_dir = out_dir + 'reports/'
    if not os.path.isdir( report_dir ):
        os.mkdir( report_dir )

    return report_dir


def ObyOprep(
    db, report_dir, blast, query, 
    max_hits, evalue, out_dir
    ):

    prev, finished, rundb = False, set(), db
    log_str = report_dir + '\n' + blast + '\n' + query + '\n' + \
        str(max_hits) + '\n' + str(evalue)
    log_name = out_dir + '.' + os.path.basename(out_dir[:-1]) + '.log'
    if not os.path.isfile( log_name ):
        with open( log_name, 'w' ) as out:
            out.write( log_str )
    else:
        with open( log_name, 'r' ) as raw:
            log = raw.read()
        if log.rstrip() == log_str.rstrip():
            eprint('\nPicking up previous run', flush = True)
            prev = True
        else:
            eprint('\nInconsistent run parameters in log, rerunning', flush = True)
            count = 0
            while os.path.isdir( report_dir ):
                count += 1
                report_dir = report_dir[:-1] + str(count) + '/'
            os.mkdir( report_dir )
            log_str = report_dir + '\n' + blast + '\n' + query + '\n' + \
                str(max_hits) + '\n' + str(evalue)  
            with open( log_name, 'w' ) as out:
                out.write( log_str )
    if prev:
        reports = collect_files(report_dir, 'tsv')
        finished = {
            os.path.basename(x)[:-4] for x in reports \
            if os.path.getsize(x) > 0
            }
        rundb = db.loc[~db['internal_ome'].isin(finished)]

    return rundb


def ObyOblast(
    rundb, blast, seq_type, report_dir, biotype,
    env_dir, query, hsps, max_hits, evalue, cpus,
    bitscore, pident
    ):

# NEED to add other threshold options
    if len(rundb) > 0:
        print('\nBlasting', flush = True)
        blast_tups = compBlastTups(
            rundb, blast, seq_type, 
            report_dir, biotype, env_dir, query, 
            hsps = hsps, max_hits = max_hits, evalue = evalue
            )

        blast_outs = multisub( blast_tups, processes = cpus )
   
    parse_tups = []
    for i, row in db.iterrows():
        if not pd.isnull(row[biotype]):
            parse_tups.append( [row['internal_ome'],
                report_dir + row['internal_ome'] + '.tsv',
                bitscore, pident ] )

    print('\nParsing reports', flush = True)
    with mp.get_context( 'spawn' ).Pool( processes = cpus ) as pool:
        results = pool.starmap( parseOutput, tuple(parse_tups) )
    results_dict = {x[0]: x[1] for x in results}

    return results_dict


def checkBlastDB():

    db_date = os.path.basename(masterDB())
    blast_db = formatPath('$MYCOFAA/blastdb/' + db_date + '.00.psd')
    if os.path.isfile(blast_db):
        return blast_db[:-7]
    elif os.path.isfile(blast_db[:-7] + '.psd'):
        return blast_db[:-7]


def dbBlast(
    db_path, blast_type, query, 
    evalue, hsps, cpus,
    report_dir
    ):

    out_file = report_dir + os.path.basename(db_path)[:-3] + '.out'
    blast_scaf = [
        blast_type, '-query', query, '-outfmt', '6',
        '-db', db_path, '-num_threads', str(cpus), '-out',
        out_file, '-num_descriptions', str(135904025),
        '-num_alignments', str(135904025)
        ]

    if evalue:
        blast_scaf.extend(['-evalue', str(evalue)])
    if hsps:
        blast_scaf.extend(['-max_hsps', str(hsps)])

    blast_call = subprocess.call(
        blast_scaf#, stdout = subprocess.PIPE,
    #    stderr = subprocess.PIPE
        )

    return blast_call, out_file


def parseDBout( db, file_, bitscore = 0, pident = 0, max_hits = None ):

    ome_results = {}
    with open( file_, 'r' ) as raw:
        for line in raw:
            data = [x for x in line.rstrip().split('\t')]
            ome = re.search(r'(^[^_]*)', data[1])[1]
            if int(float(data[-1])) > bitscore and int(1000*float(data[2])) > 100000 * pident:
                if ome not in ome_results:
                    ome_results[ome] = []
                ome_results[ome].append(data)

    x_omes = set(db['internal_ome'])
    omes_results = {
        x: ome_results[x] for x in ome_results if x in x_omes
        }

    if max_hits:
        out_results = {}
        for ome in ome_results:
            ome_results[ome].sort(key = lambda x: float(x[-1]), reverse = True)
            out_results[ome] = ome_results[ome][:max_hits]
        return out_results
    else:
        for ome in ome_results:
            ome_results[ome].sort(key = lambda x: float(x[-1]), reverse = True)
            return ome_results


def main( 
    db, blast, query, out_dir, hsps = 1, 
    max_hits = None, evalue = 10, bitscore = 0, 
    pident = 0, #coverage = 0, 
    cpus = 1, force = False
    ):

    if blast in {'tblastn', 'blastp'}:
        seq_type = 'prot'
        biotype = 'proteome'
        env_dir = os.environ['MYCOFAA'] + '/'
        blast_db = checkBlastDB()
    elif blast in {'blastx', 'blastn'}:
        seq_type = 'nucl'
        biotype = 'assembly'
        env_dir = os.environ['MYCOFNA'] + '/'
        blast_db = None
    else:
        eprint('\nERROR: invalid blast type: ' + blast, flush = True)

    report_dir = prepOutput(out_dir, query)
    if blast_db and not force:
        print('\nBlasting using MycotoolsDB blastdb')
        blast_exit, blast_output = dbBlast(
            blast_db, blast, query, evalue, 
            hsps, cpus, report_dir
            )
        if blast_exit:
            eprint('\nERROR: blast failed: ' + str(blast_exit))
            sys.exit(10)
        results_dict = parseDBout(
            db, blast_output, bitscore = bitscore, 
            pident = pident, max_hits = max_hits
            )
    else:
        print('\nBlasting on an ome-by-ome basis') 
        rundb = ObyOprep(
            db, report_dir, blast, query, 
            max_hits, evalue, out_dir
            )
        results_dict = ObyOblast(
            rundb, blast, seq_type, report_dir, biotype,
            env_dir, query, hsps, max_hits, evalue, cpus,
            bitscore, pident
            )

    print('\nCompiling fastas', flush = True)
    output_res = compileResults( results_dict )
    acc2fa_cmds = compAcc2fa( db, biotype, env_dir, output_res )
    output_fas = {}
    for query in acc2fa_cmds:
        print('\t' + query, flush = True)
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
            'for each query. --evalue and --maxhit thresholds are applied at ' + \
            'the BLAST step and are therefore incompatible with --previous. Be ' +
            'aware that db2blast.py will preferably blast against a blastdb of ' +
            'Mycotools -omes; however, if the blastdb was not generated from ' +
            'the current masterDB (' + os.path.basename(masterDB()) + ') then ' +
            'db2blast.py will blast on an ome-by-ome basis. In the latter case, ' +
            'e-value is dependent on database size, so it will threshold ' +
            'variably for each ome. If no blastdb exists for the current masterDB ' +
            'in ' + os.environ['MYCOFAA'] + '/blastdb/ then it is recommended to ' +
            'use the bit score threshold as it is not dependent on database size.'
        )
    parser.add_argument( '-b', '--blast', required = True, 
        help = 'Blast type { blastn, blastp, tblastn, blastx }' )
    parser.add_argument( '-q', '--query', required = True, help = 'Query fasta' )
    parser.add_argument( '-e', '--evalue', type = int,
        help = 'Negative e-value order max, e.g. 2 for 10^-2' )
    parser.add_argument( '-s', '--bitscore', default = 0,
        type = int, help = 'Bit score minimum' )
    parser.add_argument( '-i', '--identity', default = 0,
        type = float, help = 'Identity minimum, e.g. 0.6' )
    #parser.add_argument( '-c', '--coverage', type = float, help = 'Query coverage +/-, e.g. 0.5' )
    parser.add_argument( '-m', '--maxhits', type = int, help = 'Max hits for each organism' )
    parser.add_argument( '-d', '--database', default = masterDB(), 
        help = 'DEFAULT: masterdb' )
    parser.add_argument( '-o', '--output', default = os.getcwd() )
    parser.add_argument( '--cpu', default = mp.cpu_count(), type = int, help = 'DEFAULT: all' )
    parser.add_argument( '-f', '--force', action = 'store_true', help = 'Force ome-by-ome blast' )
    args = parser.parse_args()

    db_path = formatPath( args.database )
    if args.cpu > mp.cpu_count():
        cpus = mp.cpu_count()
    else:
        cpus = args.cpu
  
    evalue = None 
    if args.evalue:
        evalue = 10**( -args.evalue)
        
    start_args = {
        'Blast': args.blast, 'query': formatPath(args.query), 
        'database': db_path, 'output': args.output, 'max evalue': evalue,
        'min bitscore': args.bitscore, 'max hits/organism': args.maxhits, #'coverage': args.coverage, 
        'min identity': args.identity, 'cores': cpus, 'force ome-by-ome': args.force
        }

    start_time = intro( 'db2blast', start_args )
    date = start_time.strftime('%Y%m%d')
    findExecs( [args.blast], exit = {args.blast} )
   
    output_dir = formatPath( args.output ) + date + '_db2blast/'

    db = db2df(db_path)

    output_fas = main( 
        db, args.blast, args.query, output_dir, 
        bitscore = args.bitscore, pident = args.identity,
        max_hits = args.maxhits, cpus = cpus, force = args.force 
        )
    if not os.path.isdir( output_dir + 'fastas/' ):
        os.mkdir( output_dir + 'fastas/' )

    for query in output_fas:
        with open( output_dir + 'fastas/' + query + '.blast.fa', 'w' ) as out:
            out.write( output_fas[query] )

    outro( start_time )
