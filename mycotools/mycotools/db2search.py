#! /usr/bin/env python3

import argparse, subprocess, os, datetime, sys, re, copy
import multiprocessing as mp, numpy as np
from mycotools.lib.kontools import eprint, intro, outro, formatPath, multisub, collect_files, findExecs
from mycotools.lib.dbtools import mtdb, masterDB
from mycotools.lib.biotools import dict2fa, fa2dict
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
    for i, ome in enumerate(seq_db['internal_ome']):
        if seq_db[biotype][i]:
            blast_cmds.append( compileBlastCmd(
                ome, seq_db[biotype][i], env_dir,
                out_dir, blast_scaf
                ) )

    return blast_cmds


def compMMseqTups( 
    seq_db, out_dir, 
    biotype, env_dir, query, mmseqs = 'mmseqs',
    evalue = None
    ):

    cmd_scaf = [
        mmseqs, 'easy-search', query, 'tmp'#,
#        '--format-mode', '0'
        ]
    if evalue:
        cmd_scaf.extend(['-e', str(evalue)])
#    if max_hits:
#        cmd_scaf.extend([ '-max_target_seqs', str(max_hits) ])

    search_cmds = []
    for i, ome in enumerate(seq_db['internal_ome']):
        if seq_db[biotype][i]:
            new_cmd = copy.deepcopy(cmd_scaf)
            inputFile = formatPath(env_dir + seq_db[biotype][i])
            outputFile = formatPath(out_dir + ome + '.out')
            new_cmd.insert(3, inputFile)
            new_cmd.insert(4, outputFile)
#            search_cmds.append([
 #               mmseqs, query, formatPath(env_dir + row[biotype]),
  #              out_dir + ome + '.out', 'tmp'
   #             ])
            search_cmds.append(new_cmd)


    return search_cmds



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


def compAcc2fa( db, biotype, env_dir, output_res, subhit = False, skip = None ):

    acc2fa_cmds, db = {}, db.set_index('internal_ome')
    if subhit or biotype == 'assembly':
        for query in output_res:
            if skip:
                continue
            acc2fa_cmds[query] = []
            for i in output_res[query]:
                accs = []
                for hit in output_res[query][i]:
                    accs.append(hit[0] + '[' + hit[1] + '-' + hit[2] + ']')
                if biotype == 'assembly':
                    acc2fa_cmds[query].append(
                        [accs, env_dir + db[i][biotype]]
                        )
                else:
                    acc2fa_cmds[query].append(
                        [accs, env_dir + db[i][biotype]]
                        )
    else:
        for query in output_res:
            if skip:
                continue
            acc2fa_cmds[query] = []
            for i in output_res[query]:
                accs = []
                for hit in output_res[query][i]:
                    subject = hit[0]
                    accs.append(subject)
                if accs:
                    try:
                        acc2fa_cmds[query].append( [
                            list(set(accs)), env_dir + db[i][biotype]
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

def db2searchLog(report_dir, blast, query, max_hits, evalue, bit, pident, out_dir):
    log_list0 = [
        report_dir, blast, query, str(max_hits), str(evalue),
        str(bit), str(pident)
        ]
    
    prev = False    
    log_name = out_dir + '.' + os.path.basename(out_dir[:-1]) + '.log'
    log_list1 = None
    if not os.path.isfile( log_name ):
        with open( log_name, 'w' ) as out:
            out.write( '\n'.join(log_list0) )
    else:
        with open( log_name, 'r' ) as raw:
            log_list1 = raw.read().split('\n')
        if log_list1[:5] == log_list0[:5]:
            eprint('\nPicking up previous run', flush = True)
            prev = True
        else:
            eprint('\nInconsistent run parameters in log, rerunning', flush = True)
            count = 0
            while os.path.isdir( report_dir ):
                count += 1
                report_dir = report_dir[:-1] + str(count) + '/'
            os.mkdir( report_dir )
            with open( log_name, 'w' ) as out:
                out.write( '\n'.join(log_list0) )

    return log_list0, log_list1, prev


def ObyOprep(
    db, report_dir, blast, query, 
    max_hits, evalue, out_dir, bit, pident
    ):

    prev, finished, rundb = False, set(), db

    log_list0, log_list1, prev = db2searchLog(report_dir, blast, query, max_hits, evalue, bit, pident, out_dir)

    reparse = True
    if prev:
        reparse = False
        reports = collect_files(report_dir, 'tsv')
        finished = {
            os.path.basename(x)[:-4] for x in reports \
            if os.path.getsize(x) > 0
            }
        rundb, checkdb = mtdb({}).set_index('internal_ome'), db.set_index('internal_ome')
        for ome, val in checkdb.items():
            if ome not in finished:
                rundb[ome] = val
        rundb = rundb.reset_index()
        if log_list0[-2:] != log_list1[-2:]:
            reparse = True
#        print(rundb, flush = True)

    return rundb


def ObyOsearch(
    db, rundb, blast, seq_type, report_dir, biotype,
    env_dir, query, hsps, max_hits, evalue, cpus,
    bitscore, pident
    ):

# NEED to add other threshold options
    if len(rundb) > 0:
        print('\nSearching', flush = True)
        if 'blast' in blast:
            search_tups = compBlastTups(
                rundb, blast, seq_type, 
                report_dir, biotype, env_dir, query, 
                hsps = hsps, max_hits = max_hits, evalue = evalue
                )
        else:
            search_tups = compMMseqTups(
                rundb, report_dir, 
                biotype, env_dir, query, mmseqs = 'mmseqs',
                evalue = evalue
                )

        search_outs = multisub( search_tups, processes = cpus )
    
   
    parse_tups = []
    for i, ome in enumerate(db['internal_ome']):
        if db[biotype][i]:
            parse_tups.append( [ome,
                report_dir + ome + '.tsv',
                bitscore, pident ] )

    print('\nParsing reports', flush = True)
    with mp.get_context( 'spawn' ).Pool( processes = cpus ) as pool:
        results = pool.starmap( parseOutput, tuple(parse_tups) )
    results_dict = {x[0]: x[1] for x in results}

    return results_dict


def checkSearchDB(binary = 'blast'):

    db_date = os.path.basename(masterDB())
    if 'blast' in binary:
        search_db = formatPath('$MYCOFAA/blastdb/' + db_date + '.00.psd')
        if os.path.isfile(search_db):
            return search_db[:-7]
        elif os.path.isfile(search_db[:-7] + '.psd'):
            return search_db[:-7]
    else:
        search_db = formatPath('$MYCOFAA/blastdb/' + db_date.replace('.db','') + '.mmseqs.db')
        if os.path.isfile(search_db):
            return search_db


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


def dbmmseq(
    db_path, query, evalue, cpus, report_dir, mmseqs = "mmseqs", mem = None
    ):

    out_file = report_dir + os.path.basename(db_path)[:-3].replace('.mmseqs','') + '.out'
#    output_str = '"query,target,pident,alen,mismatch,gapopen,qstart,qend,sstart,send,evalue,bits"'
    cmd_scaf = [
        mmseqs, 'easy-search', formatPath(query), db_path, out_file, 'tmp',
#        '--format-output', "0", 
        '--threads', str(cpus*2),
        ]
    if mem:
        cmd_scaf.extend(['--split-memory-limit', str(mem)])

    if evalue:
        cmd_scaf.extend(['-e', str(evalue)])

    cmd_call = subprocess.call(
        cmd_scaf#, stdout = subprocess.PIPE,
    #    stderr = subprocess.PIPE
        )

    return cmd_call, out_file


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
    pident = 0, mem = None, #coverage = 0,
    cpus = 1, force = False, biotype = None,
    skip = None
    ):

    if blast in {'tblastn', 'blastp'}:
        seq_type = 'prot'
        biotype = 'proteome'
        env_dir = os.environ['MYCOFAA'] + '/'
        search_db = checkSearchDB()
    elif blast in {'blastx', 'blastn'}:
        seq_type = 'nucl'
        biotype = 'assembly'
        env_dir = os.environ['MYCOFNA'] + '/'
        search_db = None
    elif blast == 'mmseqs':
        search_db = checkSearchDB('mmseqs')
        if biotype == 'proteome':
            env_dir = os.environ['MYCOFAA'] + '/'
        elif biotype == 'assembly':
            env_dir = os.environ['MYCOFNA'] + '/'
    else:
        eprint('\nERROR: invalid search binary: ' + blast, flush = True)

    report_dir = prepOutput(out_dir, query)
    if search_db and not force:
        print('\nSearching using MycotoolsDB searchdb', flush = True)
        if 'blast' in blast:
            search_exit, search_output = dbBlast(
                search_db, blast, query, evalue, 
                hsps, cpus, report_dir
                )
        else:
            search_exit, search_output = dbmmseq(
                search_db, query, evalue, cpus, report_dir, mmseqs = "mmseqs", mem = mem
                )
        if search_exit:
            eprint('\nERROR: search failed: ' + str(search_exit))
            sys.exit(10)
        results_dict = parseDBout(
            db, search_output, bitscore = bitscore, 
            pident = pident, max_hits = max_hits
            )
    else:
        print('\nSearching on an ome-by-ome basis', flush = True) 
        rundb = ObyOprep(
            db, report_dir, blast, query, 
            max_hits, evalue, out_dir, bitscore, pident
            )
        results_dict = ObyOsearch(
            db, rundb, blast, None, report_dir, biotype,
            env_dir, query, hsps, max_hits, evalue, cpus,
            bitscore, pident
            )

    print('\nCompiling fastas', flush = True)
    output_res = compileResults( results_dict )
    acc2fa_cmds = compAcc2fa( db, biotype, env_dir, output_res, skip = None )
    output_fas = {}
    queryfa = fa2dict(query)
    for query1, cmd in acc2fa_cmds.items():
        output_fas[query1] = {}
        print('\t' + query1, flush = True)
        with mp.get_context('spawn').Pool(processes = cpus) as pool:
            results = pool.starmap( acc2fa, acc2fa_cmds[query1] )
        for x in results:
            output_fas[query1].update(x)
#        output_fas[query1][query1] = queryfa[query1]
        #'\n'.join([dict2fa(x) for x in results])

    return output_fas


if __name__ == '__main__':
# NEED abstract covered portion
# NEED max hits after blast compiled
# NEED hsp option
    mp.set_start_method('spawn')
    parser = argparse.ArgumentParser( 
        description = 'Searches a query against a db and compiles output fastas ' + \
            'for each query. --evalue and --maxhit thresholds are applied at ' + \
            'the search step and are therefore incompatible with --previous. Be ' +
            'aware that db2search.py will preferably search against a pre-indexed database of ' +
            'Mycotools -omes; however, if the database was not generated from ' +
            'the current masterDB (' + os.path.basename(masterDB()) + ') then ' +
            'db2search.py will search on an ome-by-ome basis. In the latter case, ' +
            'e-value is dependent on database size, so it will threshold ' +
            'variably for each ome. If no search database exists for the current masterDB ' +
            'in ' + os.environ['MYCOFAA'] + '/blastdb/ then it is recommended to ' +
            'use the bit score threshold as it is not dependent on database size.'
        )
    parser.add_argument( '-s', '--search', required = True, 
        help = 'Search binary { mmseqs, blastn, blastp, tblastn, blastx }' )
    parser.add_argument('-st', '--sequencetype', help = 'Subject sequence type {aa, nt} for mmseqs')
    parser.add_argument( '-q', '--query', required = True, help = 'Query fasta' )
    parser.add_argument( '-e', '--evalue', type = int,
        help = 'Negative e-value order max, e.g. 2 for 10^-2' )
    parser.add_argument( '-b', '--bitscore', default = 0,
        type = int, help = 'Bit score minimum' )
    parser.add_argument( '-i', '--identity', default = 0,
        type = float, help = 'Identity minimum, e.g. 0.6' )
    #parser.add_argument( '-c', '--coverage', type = float, help = 'Query coverage +/-, e.g. 0.5' )
    parser.add_argument( '-m', '--maxhits', type = int, help = 'Max hits for each organism' )
    parser.add_argument( '-d', '--database', default = masterDB(), 
        help = 'DEFAULT: masterdb' )
    parser.add_argument( '-o', '--output', default = os.getcwd() )
    parser.add_argument( '--cpu', default = mp.cpu_count(), type = int, help = 'DEFAULT: all' )
    parser.add_argument( '--ram', help = 'Useful for mmseqs: e.g. 10M or 5G' )
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

    if args.search == 'mmseqs':
        if not args.sequencetype or args.sequencetype not in {'aa', 'nt'}:
            eprint('\nERROR: -st required for mmseqs', flush = True)
            sys.exit(3)
        if args.sequencetype == 'aa':
            biotype = 'proteome'
        else:
            biotype = 'assembly'
    else:
        biotype = None
        
    start_args = {
        'Search binary': args.search, 'query': formatPath(args.query), 
        'database': db_path, 'output': args.output, 'max evalue': evalue,
        'min bitscore': args.bitscore, 'max hits/organism': args.maxhits, #'coverage': args.coverage, 
        'min identity': args.identity, 'cores': cpus, 'force ome-by-ome': args.force
        }

    start_time = intro( 'db2search', start_args )
    date = start_time.strftime('%Y%m%d')
    findExecs( [args.search], exit = {args.search} )
   
    output_dirPrep = formatPath( args.output ) 
    if not output_dirPrep.endswith('/'):
        output_dirPrep = re.sub(r'/+$', '', output_dirPrep)
    output_dir = output_dirPrep + date + '_db2search/'

    db = mtdb(db_path)

    output_fas = main( 
        db, args.search, args.query, output_dir, 
        bitscore = args.bitscore, pident = args.identity, mem = args.ram,
        max_hits = args.maxhits, cpus = cpus, force = args.force, biotype = biotype
        )
    if not os.path.isdir( output_dir + 'fastas/' ):
        os.mkdir( output_dir + 'fastas/' )

    for query in output_fas:
        with open( output_dir + 'fastas/' + query + '.search.fa', 'w' ) as out:
            out.write( dict2fa(output_fas[query]) )

    outro( start_time )
