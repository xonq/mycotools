#! /usr/bin/env python3

# NEED distinguish nt and aa mmseqs dbs
# NEED blastdb
# NEED streamline mmseqs parsing
# NEED mmseqs mtdb check
# NEED mmseqs save db option
# NEED profile mmseqs
# NEED concatenate mmseqs query dbs
# NEED option to fail upon any failures
# NEED log (hmmer)
# NEED nhmmer option
# NEED to .tmp and move files

import os
import re
import sys
import copy
import datetime
import argparse
import subprocess
import multiprocessing as mp
from io import StringIO
from collections import defaultdict
from mycotools.db2files import soft_main as db2files
from mycotools.lib.kontools import intro, outro, collect_files, multisub, \
    findExecs, untardir, eprint, format_path, mkOutput, tardir, inject_args, stdin2str
from mycotools.lib.dbtools import primaryDB, mtdb
from mycotools.lib.biotools import dict2fa, fa2dict
#from mycotools.extractHmmsearch import main as exHmm
from mycotools.acc2fa import dbmain as acc2fa_db, famain as acc2fa_fa
from mycotools.utils.extractHmmsearch import main as exHmm
from mycotools.utils.extractHmmAcc import grabAccs, main as absHmm


def compile_hmm_cmd(db, hmm_path, output, ome_set = set(), cpu = 1):
    """
    Inputs: mycotools db, hmm_path, output directory, set of omes to ignore.
    Outputs: tuples of arguments for hmmsearches
    """

    cmd_tuples = [ ]
    for ome, row in db.items():
        if ome not in ome_set:
            output_path = output + ome + '.out'
            cmd = (('hmmsearch', '-o', output_path + '.tmp', 
                   '--cpu', str(cpu),
                   hmm_path, row['faa'],
                   '&&',), ('mv', output_path + '.tmp', output_path,),) 
            cmd_tuples.append(cmd)

    return cmd_tuples

def compileextractHmmCmd(db, args, output):
    """
    Inputs: mycotools db, argparse arguments, and output path
    Outputs: tuples of arguments for `run_ex_hmm`
    """

    cmd_tuples = []
    for ome, row in db.items():
        hmmsearch_out = output + ome + '.out'
        cmd_tuples.append((args, hmmsearch_out, output,))

    return tuple(cmd_tuples)


def run_ex_hmm(args, hmmsearch_out, output):

    ome = os.path.basename(hmmsearch_out).replace('.out', '')

    try:
        with open(hmmsearch_out, 'r') as raw:
            data = raw.read()
    except FileNotFoundError:
        eprint('\tWARNING: ' + ome + ' failed', flush = True)
        return ome, False
    if len(data) > 100: # check for data # check for data
        hmm_data = exHmm(data, args[0], args[1], args[2], args[3], args[4],
                         header = False)
        out_dict = {}
        for q, data in hmm_data.items(): # parse alignment info
            hit_info = [x.split('\t') for x in data[1].rstrip().split('\n')]
            out_dict[q] = tuple([(x[0], int(x[9]), int(x[10]),) for x in hit_info])
            # hit_accession, coordinate start, coordinate end
        # does this overwrite other hits?
        return ome, out_dict
    else:
        eprint('\tWARNING: ' + ome + ' empty results', flush = True)
        return ome, False

    
def comp_hmm_acc2fa(db, q_dict, coords = True):

    cmd_tuples = []
    fa_dict = {ome: row['faa'] for ome, row in db.items()}
    for q in q_dict:
        cmd_tuples.append((db, q_dict[q], 
                           q, coords,))

    return cmd_tuples


def hmm_acc2fa(db, q_files, q, coords = False):

    accs = []
    if coords:
        for ome, hit_data in q_files.items():
            t_accs = [f'{x[0]}[{x[1]}:{x[2]}]' for x in hit_data]
            accs.extend(t_accs)
    else:
        for ome, hit_data in q_files.items():
            t_accs = [x[0] for x in hit_data]
            accs.extend(t_accs)

    return q, acc2fa_db(db, accs)
#    fa_str = dict2fa(acc2fa_db(db, accs))
    
 #   with open( output, 'w' ) as out:
  #      out.write( fa_str )


def compile_mafft_cmds(output, faa_dir):
    fas = collect_files(faa_dir, 'faa') # grab completed fastas
    cmds = []
    for fa in fas:
        acc = os.path.basename(fa)[:-4]
        align = output + '/aligns/' + acc + '.mafft.fasta'
        if os.path.isfile(align): # check for data
            with open(align, 'r') as raw:
                data = raw.read()
            if len(data) > 10:
                continue
        cmds.append(f'mafft {fa} > {align}')
    return tuple(cmds)


def compile_hmmalign_cmds( output, accessions ):

#    cmd_tuples = [ [], [] ]
    cmds = []
    for acc in accessions:
        fa = output + '/fastas/' + acc + '.faa'
        hmm = output + '/hmms/' + acc + '.hmm'
        align = output + '/aligns/' + acc + '.stockholm'
        conv = output + '/aligns/' + acc + '.phylip'
        trim = output + '/trimmed/' + acc + '.clipkit.fa'
        if os.path.isfile(conv):
            with open(conv, 'r') as raw:
                data = raw.read()
            if len(data) > 10:
                continue
        cmds.append((('hmmalign', '--trim', '-o', 
                     align, hmm, fa, '&&',),
                    ('esl-reformat', '--informat', 'stockholm',
                     '-o', conv, 'phylip', align),),)

    return tuple(cmds)


def compile_trim_cmd(output, mod = '', trimmed = None, ex = 'phylip'):

    mod_args = mod.split(' ')
    cmd_tuples = []
    aligns = collect_files(output + 'aligns/', ex)
    if trimmed:
        trimmed = set(os.path.basename(x).replace('.clipkit','') \
            for x in trimmed )
        aligns = [x for x in aligns if os.path.basename(x).replace(ex,'') \
             not in trimmed]
    for align in aligns:
        trim = f'{output}trimmed/' \
             + os.path.basename(align).replace(ex,'clipkit') \
             + '.' + ex
        args = ['clipkit', align, '-o', trim]
        if mod_args[0]:
            args.extend(mod_args) 
        cmd_tuples.append( tuple( args ) )

    return cmd_tuples


def compile_hmm_queries(hmm_paths, hmm_out):
    queries, complete_hmm = [], ''
    for hmm_path in hmm_paths:
        with open(hmm_path, 'r') as hmm_raw:
            hmm_data = hmm_raw.read()
            queries.extend(grabAccs(hmm_data))
            complete_hmm += hmm_data + '\n'
    with open(hmm_out, 'w') as hmm_oh:
        hmm_oh.write(complete_hmm.rstrip())
    return queries


def hmmer_main(db, hmm_paths, output, accessions, max_hits, query_cov,
               coords = False, evalue = 0.01, bitscore = 0,
               binary = 'hmmsearch', verbose = False, cpu = 1):

    ome_dir = output + 'omes/'
    faa_dir = output + 'fastas/'

    db = db.set_index('ome')
    hmm_out = output + 'query.hmm'

    queries = compile_hmm_queries(hmm_paths, hmm_out)

    ome_set, skip1 = set(), False
    if os.path.isdir(faa_dir): # is there a previous run?
        print('\nCompiling previous run', flush = True)
        # check if all fastas are made
        fas = collect_files(faa_dir, 'faa') # grab completed fastas
        ranQueries = [os.path.basename(fa).replace('.faa','') for fa in fas]
        if not set(queries).difference(set(ranQueries)):
            skip1 = True

    if not skip1: # do not skip the first step
        if os.path.isfile(output + 'omes.tar.gz'):
            if not os.path.isdir(ome_dir):
                untardir(output + 'omes.tar.gz')
        # check what reports have been generated
        omes = collect_files(ome_dir, 'out')
        ome_set = set(os.path.basename(x).replace('.out', '') for x in omes)
    else:
        print('\thmmsearch -> hits.faa DONE', flush = True)

    if not skip1:
        print('\nRunning ' + os.path.basename(binary), flush = True)
        if not os.path.isdir(ome_dir):
            os.mkdir(ome_dir)

        # run hmmer
        par_runs = round(((cpu - 1)/2) - 0.5)
        if par_runs < 1:
            par_cpus = 1
            par_runs = 1
        else:
            par_cpus = 2
        hmmsearch_tuples = compile_hmm_cmd(db, hmm_out, 
                                         ome_dir, 
                                         ome_set = ome_set,
                                         cpu = par_cpus)
        hmmsearch_codes = multisub(hmmsearch_tuples, processes = par_runs,
                                   injectable = True, verbose = verbose)
        for i, code in enumerate(hmmsearch_codes):
            if code:
                eprint('\tERROR: ' \
                        + str(hmmsearch_tuples[i]),
                        flush = True)

        # extract results
        print('\nExtracting hmmsearch output', flush = True)
        exHmm_args = [accessions, max_hits, query_cov, evalue, bitscore]
        exHmm_tuples = compileextractHmmCmd(db, exHmm_args, ome_dir)
        with mp.get_context('spawn').Pool(processes = cpu) as pool:
            hmmAligns = pool.starmap(run_ex_hmm, exHmm_tuples)

        mp.Process(target=tardir, args=[ome_dir])
        print( '\nCompiling fastas' , flush = True)
        # q_dict = {query: ome: alignment}
        q_dict = defaultdict(dict)
        for ome, hits in hmmAligns:
            if hits:
                for query, hit_info in hits.items():
#                    align = hit[query][1]
                    q_dict[query][ome] = hit_info

        acc2fa_tuples = comp_hmm_acc2fa(db, q_dict, coords = coords)
        with mp.get_context('spawn').Pool(processes = cpu) as pool:
            fa_info = pool.starmap(hmm_acc2fa, acc2fa_tuples)        
        fa_dicts = {q: fd for q, fd in fa_info}

        return fa_dicts


def compileBlastCmd( ome, biofile, out_dir, blast_scaf ):
    return blast_scaf + ['-out', out_dir + ome + '.tsv', \
        '-subject', biofile]

def comp_blast_tups( 
    seq_db, blast_type, seq_type, out_dir, 
    biotype, query, hsps = None, max_hits = None,
    evalue = None, coverage = None, search_args = None
    ):

    blast_scaf = [
        blast_type, '-query', query
        ]
    if evalue:
        blast_scaf.extend(['-evalue', str(evalue)])
    if hsps:
        blast_scaf.extend(['-max_hsps', str(hsps)])
    if max_hits:
        blast_scaf.extend(['-max_target_seqs', str(max_hits)])
    if coverage:
        blast_scaf.extend(['-qcov_hsp_perc', str(coverage)])
    if search_args:
        blast_scaf.extend(search_args)
    blast_scaf.extend(['-outfmt', 
                      '"6 qseqid sseqid pident ppos sstart send evalue bitscore"'])

    blast_cmds = []
    for i, ome in enumerate(seq_db['ome']):
        if seq_db[biotype][i]:
            blast_cmds.append( ' '.join(compileBlastCmd(
                ome, seq_db[biotype][i],
                out_dir, blast_scaf
                )))
    return blast_cmds

def compileDiamondCmd(ome, dmnd_db, out_dir, blast_scaf):
    return blast_scaf + ['--out', out_dir + ome + '.tsv', \
        '--subject', dmnd_db]

def comp_diamond_tups( 
    seq_db, diamond, blast_type, seq_type, out_dir, 
    biotype, query, hsps = None, max_hits = None,
    evalue = None, coverage = None, search_args = []
    ):

    if not os.path.isdir(out_dir + 'dmnd/'):
        os.mkdir(out_dir + 'dmnd/')
    blast_scaf = [
        diamond, blast_type, '-query', query,
        ]
    if evalue:
        blast_scaf.extend(['--evalue', str(evalue)])
    if hsps:
        blast_scaf.extend(['--max_hsps', str(hsps)])
    if max_hits:
        blast_scaf.extend(['--max_target_seqs', str(max_hits)])
    if coverage:
        blast_scaf.extend(['--query-cover', str(coverage)])
    if search_args:
        blast_scaf.extend(search_args)
    blast_scaf.extend('--outfmt', '6', 'qseqid', 'sseqid', 'pident',
        'ppos', 'sstart', 'send', 'evalue', 'bitscore')


    blast_cmds, db_cmds = [], []
    for i, ome in enumerate(seq_db['ome']):
        if seq_db[biotype][i]:
            db_cmds.append([
                diamond, 'makedb', '--in', seq_db[biotype][i],
                '--db', out_dir + 'dmnd/' + ome, '-p', '2'
                ])
            blast_cmds.append(compileDiamondCmd(
                ome, out_dir + 'dmnd/' + ome,
                out_dir, blast_scaf
                ))
    return db_cmds, blast_cmds



# concat, search, parse
def run_mmseq(
    seq_db, out_dir, biotype, query, mmseqs = 'mmseqs',
    coverage = None, search_args = [], cpus = 1,
    iterations = 3
    ):
    db_dir = format_path('$MYCOGFF3/../db')
    
    if not os.path.isdir(f'{out_dir}db/'):
        os.mkdir(f'{out_dir}db/')
    db_dir = format_path('$MYCOGFF3/../db/')
    createdb_cmds = []
    db_path = f'{out_dir}db/searchdb'
    
    # make mmseqs dbs for those that dont exist
    for i, ome in enumerate(seq_db['ome']):
        db_path = f'{db_dir}{ome}_{biotype}'
        out_file = out_dir + ome + '.tsv'
        if not os.path.isfile(db_path + '.dbtype'):
            # will fail at fastas that dont have sequences on one line
            createdb_cmds.append((mmseqs, 'createdb', seq_db[biotype][i], db_path,
                                  '--createdb-mode', '1', '--shuffle', '0'))

    if createdb_cmds:
        print(f'\nCreating {len(createdb_cmds)} mmseqs search dbs', flush = True)
        createdb_outs = multisub(createdb_cmds, processes = cpus, verbose = 2) 

    #if len(query) > 1:
   #     if not os.path.isfile(f'{out_dir}db/query.dbtype'):
  #          mergedbs_cmd = ['mmseqs', 'mergedbs']
 #           mergedbs_cmd.extend([q for q in query])
#            mergedbs_cmd.append('--prefixes')
           # mergedbs_cmd.extend([os.path.basename(q) for q in query])
     #       mergedbs_cmd.insert(3, f'{out_dir}db/querydb')
      #      mergedbs_out = subprocess.call(mergedbs_cmd)
        #query = [f'{out_dir}db/query'] # need to adjust check
        

    # create a concatenated mmseqs db for the search target
    if not os.path.isfile(f'{out_dir}db/searchdb.dbtype'):
        mergedbs_cmd = ['mmseqs', 'mergedbs']
        mergedbs_cmd.extend([f'{db_dir}{ome}_{biotype}' for ome in seq_db['ome']])
        mergedbs_cmd.insert(3, f'{out_dir}db/searchdb')
        print('\nMerging search dbs', flush = True)
        mergedbs_out = subprocess.call(mergedbs_cmd) #, stderr = subprocess.DEVNULL,
#                                        stdout = subprocess.DEVNULL)
    
    print('\nSearching', flush = True)
    for i, q in enumerate(query):
        out_file = f'{out_dir}{q}.tsv'
        if os.path.isfile(out_file):
            continue
        print('\t' + q, flush = True)
        search_cmd = [mmseqs, 'search', q, f'{out_dir}db/searchdb', 
                      out_file + '.tmp', f'{out_dir}tmp/', 
                      '--remove-tmp-files', '1',
                      '--threads', str(cpus*2), '--num-iterations', str(iterations)]
        search_cmd.extend(search_args)
        if coverage:
            search_cmd.extend(['-c', str(coverage)])

        search_out = subprocess.call(search_cmd) #, stderr = subprocess.DEVNULL,
#                                     stdout = subprocess.DEVNULL)
        results_cmd = [mmseqs, 'convertalis', q, f'{out_dir}db/searchdb',
                       out_file + '.tmp', out_file,                        
                       '--format-output',
                       'qset,target,pident,tstart,tend,evalue,bits']
        results_out = subprocess.call(results_cmd, stderr = subprocess.DEVNULL,
                                      stdout = subprocess.DEVNULL)

def parseOutput(algorithm, ome, file_, bitscore = 0, pident = 0, 
                evalue = 0, max_hits = None, ppos = 0, scale = 1000):

    ome_results = [ome, []]
    if os.path.exists( file_ ):
        with open( file_, 'r' ) as raw:
            data = [x.rstrip().split('\t') for x in raw if x.rstrip()]
        byq = defaultdict(list)
        for d in data:
            byq[d[0]].append(d)
        # sort by percent identity (we denote everything before here as a homolog)
        for q, data in byq.items():
            d = sorted(data, key = lambda x: float(x[2]), reverse = True)
            if max_hits:
                d = d[:max_hits]
            for i in d:
                if int(float(i[-1])) > bitscore \
                  and int(1000*float(i[2])) > scale * pident \
                  and float(i[-2]) <= evalue and int(1000*float(i[3])) > scale * ppos:
                    ome_results[-1].append(i)
    return ome_results

def parseOutput_mmseqs(algorithm, ome, file_, bitscore = 0, pident = 0, 
                evalue = 0, max_hits = None, ppos = None):

    ome_results = [ome, []]
    if os.path.exists( file_ ):
        with open( file_, 'r' ) as raw:
            data = [x.rstrip().split('\t') for x in raw if x.rstrip()]
        byq = defaultdict(list)
        for d in data:
            byq[d[0]].append(d)
        # sort by percent identity (we denote everything before here as a homolog)
        for q, data in byq.items():
            d = sorted(data, key = lambda x: float(x[2]), reverse = True)
            if max_hits:
                d = d[:max_hits]
            for i in d:
                if int(float(i[-1])) > bitscore \
                  and int(1000*float(i[2])) > 100000 * pident \
                  and float(i[-2]) <= evalue:
                    ome_results[-1].append(i)
    return ome_results


def compileResults( res_dict, skip = [] ):

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

    for query in skip:
        del output_res[query]

    return output_res

def comp_mmseq_acc2fa(db, biotype, output_res, coords = False, skip = None):
 
    cmd_tuples = []
    fa_dict = {ome: row['faa'] for ome, row in db.set_index().items()}
    for q, ome_dict in output_res.items():
        q_accs = []
        if coords:
            for ome, hit_list in ome_dict.items():
                for hit, start, end in hit_list:
                    q_accs.append(f'{hit}[{start}:{end}]')
        else:
            for ome, hit_list in ome_dict.items():
                for hit, null0, null1 in hit_list:
                    q_accs.append(f'{hit}')
            
        cmd_tuples.append((db, q, q_accs,))
  
    return cmd_tuples 

def mmseq_acc2fa(db, q, q_accs):
    return q, acc2fa_db(db, q_accs)

def comp_blast_acc2fa( db, biotype, output_res, coords = False, skip = None ):

    acc2fa_cmds, db = {}, db.set_index('ome')
    if coords or biotype == 'assembly':
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
                        [accs, db[i][biotype]]
                        )
                else:
                    acc2fa_cmds[query].append(
                        [accs, db[i][biotype]]
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
                            list(set(accs)), db[i][biotype]
                            ] )
                    except KeyError:
                        eprint('\t' + i + ' not in db')

    return acc2fa_cmds


def prepOutput(out_dir):

    out_dir = format_path(out_dir)
    if not out_dir.endswith('/'):
        out_dir += '/'
    if not os.path.isdir( out_dir ):
        os.mkdir( out_dir )
    report_dir = out_dir + 'reports/'
    if not os.path.isdir( report_dir ):
        os.mkdir( report_dir )

    return report_dir


def run_denovo(report_dir, log_list0, log_name):
    if os.path.isdir(report_dir):
        count = 0
        while os.path.isdir(report_dir[:-1] + str(count)):
            count += 1
        report_dir = report_dir[:-1] + str(count) + '/'
    log_list0[0] = report_dir
    os.mkdir(report_dir)
    with open(log_name, 'w') as out:
        out.write('\n'.join(log_list0))
    return log_list0, report_dir

def db2searchLog(report_dir, blast, query, max_hits,
                 evalue, bit, pident, coverage, out_dir,
                 ppos):
    log_list0 = [
        report_dir, blast, query, str(max_hits), str(evalue),
        str(bit), str(pident), str(coverage), str(ppos)
        ]
    
    prev, reparse = False, False
    log_name = out_dir + '.' + os.path.basename(out_dir[:-1]) + '.log'
    log_list1 = None
    if not os.path.isfile(log_name): # generate a new log
        with open(log_name, 'w') as out:
            out.write('\n'.join(log_list0))
    else: # check the old one
        with open(log_name, 'r') as raw:
            log_list1 = [x.rstrip() for x in raw if x]
        if blast != log_list1[1]:
            eprint('\tInconsistent search algorithm, rerunning', flush = True)
            log_list0, report_dir = run_denovo(report_dir, log_list0, log_name)
        elif blast == 'mmseqs':
            if log_list1[-1] != log_list0[-1]: # coverage is off, need a rerun
                eprint('\tCoverage changed, rerunning', flush = True)
                log_list0, report_dir = run_denovo(report_dir, log_list0, log_name)
            # need to reparse if anything is different
            elif any(log_list0[i] != log_list1[i] for i in range(len(log_list0))):
                eprint('\tDeleting old report compilations', flush = True)
                reports = collect_files(report_dir, 'tsv')
                for r in reports:
                    os.remove(r)
                reparse = True
            prev = True
        elif log_list1[1] != log_list0[1] and log_list1[3:] != log_list0[3:]:
            eprint('\tInconsistent thresholds, rerunning', flush = True)
            log_list0, report_dir = run_denovo(report_dir, log_list0, log_name)
        else:
            prev = True

    return log_list0, log_list1, prev, reparse


def prepare_search_run(
    db, report_dir, blast, query, 
    max_hits, evalue, out_dir, bit, pident,
    coverage, ppos
    ):


    prev, finished, rundb = False, set(), db
    if isinstance(query, list):
        query = ','.join(query)
    log_list0, log_list1, prev, reparse = db2searchLog(report_dir, blast, query, max_hits, 
                                              evalue, bit, pident, coverage, out_dir,
                                              ppos)
    report_dir = log_list0[0]
    if prev:
  #      reparse = False
        reports = collect_files(report_dir, 'tsv')
        finished = {
            os.path.basename(x)[:-4] for x in reports \
            if os.path.getsize(x) > 0
            }
        rundb, checkdb = mtdb({}).set_index('ome'), db.set_index('ome')
        for ome, val in checkdb.items():
            if ome not in finished:
                rundb[ome] = val
        rundb = rundb.reset_index()
 #       if log_list0[-3:] != log_list1[-3:]: # bitscore discrepancy
#            reparse = True
    return rundb, reparse, report_dir


def prep_mmseq_output(rundb, report_dir, queries, convert = False):
    ome_res = defaultdict(str)
    for i, q in enumerate(queries):
        if not os.path.isfile(f'{report_dir}{q}.tsv'):
            continue
        with open(f'{report_dir}{q}.tsv', 'r') as raw:
            base = os.path.basename(q)
            if not convert:
                for line in raw:
                    line_d = line.split()
                    ome = line_d[1][:line_d[1].find('_')]
                    ome_res[ome] += line.rstrip() + '\n'     
            else:
                for line in raw:
                    line_d = line.rstrip().split()
                    line_d[0] = base
                    ome = line_d[1][:line_d[1].find('_')]
                    ome_res[ome] += '\t'.join(line_d) + '\n'
    
    for ome, out_str in ome_res.items():
        with open(f'{report_dir}{ome}.tsv', 'r') as out:
            out.write(out_str.rstrip())

    
             
    

def comp_mmseq_res(rundb, report_dir, queries, convert = False):
    for ome in rundb['ome']:
        out_file = report_dir + ome + '.tsv'
        if os.path.isfile(out_file):
            continue
        out_str, todel = '', []
        for i, q in enumerate(queries):
            with open(out_file + '.tmp' + str(i), 'r') as raw:
                if not convert:
                    out_str += raw.read().rstrip() + '\n'
                else:
                    base = os.path.basename(q)
                    for line in raw:
                        line_d = line.split()
                        line_d[0] = base
                        out_str += '\t'.join(line_d) + '\n'
                todel.append(out_file + '.tmp' + str(i))
        with open(out_file, 'w') as out:
            out.write(out_str.rstrip())

        for todel_file in todel:
            os.remove(todel_file)

def ObyOsearch(
    db, rundb, blast, seq_type, report_dir, biotype,
    query, hsps, max_hits, evalue, cpus,
    bitscore, pident, coverage, diamond, search_arg,
    reparse = False, ppos = 0
    ):
    if len(rundb) > 0:
        print('\nSearching on an ome-by-ome basis', flush = True)
        if diamond:
            db_tups, search_tups = comp_diamond_tups(
                rundb, diamond, blast, seq_type, report_dir, biotype, query,
                hsps = hsps, evalue = evalue,
                coverage = coverage*100, search_args = search_arg,
                )
            db_outs = multisub(db_tups, processes = cpus)
            print(f'\t{len(search_tups)} searches to run', flush = True)
            search_outs = multisub(search_tups, processes = cpus, 
                                    verbose = 2, injectable = True)
            scale = 100000
        else:
            search_tups = comp_blast_tups(
                rundb, blast, seq_type, 
                report_dir, biotype, query, 
                hsps = hsps, evalue = evalue,
                coverage = coverage*100, search_args = search_arg
                )
            search_outs = multisub(search_tups, processes = cpus, 
                                    verbose = 2, shell = True, injectable = True)
            scale = 100000
   
    # prepare report parsing commands for multiprocessing
    parse_tups = []
    for i, ome in enumerate(db['ome']):
        if db[biotype][i]:
            parse_tups.append([blast, ome,
                report_dir + ome + '.tsv',
                bitscore, pident, evalue, max_hits,
                ppos, scale])

    print('\nParsing reports', flush = True)
    with mp.get_context('spawn').Pool(processes = cpus) as pool:
        results = pool.starmap(parseOutput, tuple(parse_tups))
    results_dict = {x[0]: x[1] for x in results}

    return results_dict


def mmseqs_mngr(
    db, rundb, mmseqs, report_dir, biotype,
    query, max_hits, evalue, cpus,
    bitscore, pident, coverage, search_arg,
    convert = False, reparse = False
    ):

    run_mmseq(
        rundb, report_dir, biotype, query, mmseqs = 'mmseqs',
        coverage = None, search_args = [], cpus = cpus
        )    

    print('\nExtracting ome reports', flush = True)
    prep_mmseq_output(rundb, report_dir, query, convert = convert)

    print('\nParsing output', flush = True)
    # prepare report parsing commands for multiprocessing
    parse_tups = []
    for i, ome in enumerate(db['ome']):
        if db[biotype][i]:
            parse_tups.append([mmseqs, ome,
                report_dir + ome + '.tsv',
                bitscore, pident, evalue, max_hits])

    with mp.get_context('spawn').Pool(processes = cpus) as pool:
        results = pool.starmap(parseOutput_mmseqs, tuple(parse_tups))
    results_dict = {x[0]: x[1] for x in results}

    return results_dict




def checkSearchDB(binary = 'blast'):

    db_date = os.path.basename(primaryDB())
    if 'blast' in binary:
        search_db = format_path('$MYCOFAA/blastdb/' + db_date + '.00.psd')
        if os.path.isfile(search_db):
            return search_db[:-7]
        elif os.path.isfile(search_db[:-7] + '.psd'):
            return search_db[:-7]
    else:
        search_db = format_path('$MYCOFAA/blastdb/' + db_date.replace('.db','') + '.mmseqs.db')
        if os.path.isfile(search_db):
            return search_db


def dbBlast(
    db_path, blast_type, query, 
    evalue, hsps, cpus,
    report_dir, diamond = False
    ):

    out_file = report_dir + os.path.basename(db_path)[:-3] + '.out'
    if not diamond:
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
    else:
        blast_scaf = [
            diamond, blast_type, '--query', query, '--outfmt', '6',
            '--db', db_path, '-p', str(cpus * 2), '--out',
            out_file, '--num_descriptions', '135904025',
            '--num_alignments', '135904025'
            ]
        if evalue:
            blast_scaf.extend(['--evalue', str(evalue)])
        if hsps:
            blast_scaf.extend(['--max_hsps', str(hsps)])
    
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
        mmseqs, 'easy-search', format_path(query), db_path, out_file, 'tmp',
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


def parseDBout(db, file_, bitscore = 0, pident = 0, 
               ppos = 0, max_hits = None):

    ome_results = {}
    with open( file_, 'r' ) as raw:
        for line in raw:
            data = [x for x in line.rstrip().split('\t')]
            ome = re.search(r'(^[^_]*)', data[1])[1]
            if int(float(data[-1])) > bitscore \
                and int(1000*float(data[2])) > 100000 * pident \
                and int(1000*float(data[3])) > 100000 * ppos:
                if ome not in ome_results:
                    ome_results[ome] = []
                ome_results[ome].append(data)

    x_omes = set(db['ome'])
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


def mmseqs_main( 
    db, mmseqs, query, out_dir, 
    max_hits = None, evalue = 10, bitscore = 0, 
    pident = 0, mem = None, coverage = None,
    cpus = 1, biotype = None,
    skip = [], coordinate = False,
    search_arg = [], convert = False, iterations = 3,
    ):

    report_dir = prepOutput(out_dir)
    if isinstance(query, str):
        query = [query]
    elif isinstance(query, dict):
        with open(out_dir + 'query.fa', 'w') as out:
            out.write(dict2fa(query))
        query = [out_dir + 'query.fa'] 



    rundb, reparse = prepare_search_run(
        db, report_dir, mmseqs, query, 
        max_hits, evalue, out_dir, bitscore, pident,
        coverage
        )
    results_dict = mmseqs_mngr(
        db, rundb, mmseqs, report_dir, biotype,
        query, max_hits, evalue, cpus,
        bitscore, pident, coverage, search_arg,
        convert = convert, reparse = reparse
        )

    print('\nCompiling fastas', flush = True)
    output_res = compileResults(results_dict, skip)
    output_fas = {}
    acc2fa_cmds = comp_mmseq_acc2fa(db, biotype, output_res, 
                                    coords = coordinate, skip = None)
    with mp.get_context('spawn').Pool(processes = cpus) as pool:
        results = pool.starmap(mmseq_acc2fa, acc2fa_cmds)
    output_fas = {q: v for q, v in results}
    
    return output_fas




def blast_main( 
    db, blast, query, out_dir, hsps = 1, 
    max_hits = None, evalue = 10, bitscore = 0, 
    pident = 0, coverage = None,
    cpus = 1, force = False,
    skip = [], diamond = None, coordinate = False,
    search_arg = [], ppos = 0
    ):

       
    if blast in {'tblastn', 'blastp'}:
        seq_type = 'prot'
        biotype = 'faa'
        search_db = checkSearchDB()
    elif blast in {'blastx', 'blastn'}:
        seq_type = 'nucl'
        biotype = 'fna'
        search_db = None
    else:
        eprint('\nERROR: invalid search binary: ' + blast, flush = True)

    report_dir = prepOutput(out_dir)
    if isinstance(query, str):
        query = [query]
    elif isinstance(query, dict):
        with open(out_dir + 'query.fa', 'w') as out:
            out.write(dict2fa(query))
        query = [out_dir + 'query.fa'] 


    query_dict = {}
    for q in query:
        query_dict = {**query_dict, **fa2dict(q)}
    with open(out_dir + 'query.fa', 'w') as out:
        out.write(dict2fa(query_dict))
    query = out_dir + 'query.fa'

    if search_db and not force:
        print('\nSearching using MycotoolsDB searchdb', flush = True)
        search_exit, search_output = dbBlast(
            search_db, blast, query, evalue, 
            hsps, cpus, report_dir, diamond = diamond
            )
        if search_exit:
            eprint('\nERROR: search failed: ' + str(search_exit))
            sys.exit(10)
        results_dict = parseDBout(
            db, search_output, bitscore = bitscore, 
            pident = pident, max_hits = max_hits,
            ppos = ppos
            )
    else:
        rundb, reparse, report_dir = prepare_search_run(
            db, report_dir, blast, query, 
            max_hits, evalue, out_dir, bitscore, pident,
            coverage, ppos
            )
        results_dict = ObyOsearch(
            db, rundb, blast, None, report_dir, biotype,
            query, hsps, max_hits, evalue, cpus,
            bitscore, pident, coverage, diamond, search_arg,
            reparse = reparse, ppos = ppos
            )

    print('\nCompiling fastas', flush = True)
    output_res = compileResults(results_dict, skip)
    output_fas = {}
    acc2fa_cmds = comp_blast_acc2fa(db, biotype, output_res, 
                             coords = coordinate, skip = None)
    queryfa = fa2dict(query)
    for query1, cmd in acc2fa_cmds.items():
        output_fas[query1] = {}
        print('\t' + query1, flush = True)
        with mp.get_context('spawn').Pool(processes = cpus) as pool:
            results = pool.starmap(acc2fa_fa, acc2fa_cmds[query1])
        for x in results:
            output_fas[query1].update(x)

    return output_fas


def cli():
# NEED abstract covered portion
# NEED max hits after blast compiled
# NEED hsp option

    sys.argv, manual_cmds = inject_args(sys.argv, ['--manual'])
    manual_cmd = manual_cmds[0]
    algorithms = {'mmseqs', 'blastn', 'blastp', 'tblastn', 'blastx', 'hmmsearch'}
    mp.set_start_method('spawn')
    parser = argparse.ArgumentParser(description = 'Searches a query sequence ' \
               + 'or profile database against an mtdb. ' \
               + 'Compiles output fastas for each query.'
               )

    i_arg = parser.add_argument_group('Inputs')
    i_arg.add_argument('-d', '--mtdb', default = primaryDB())
    i_arg.add_argument('-q', '--query', help = 'Profile database, sequence, or "-" for stdin fasta')
    i_arg.add_argument('-qd', '--query_dir', help = 'Dir of queries')
    i_arg.add_argument('-qf', '--query_file', help = 'File of query paths')

    sa_arg = parser.add_argument_group('Search algorithm')
    sa_arg.add_argument('-a', '--algorithm',
        help = f'Search binary {sorted(algorithms)}')
    sa_arg.add_argument('-s', '--seqtype', 
        help = '[mmseqs] Subject sequence type {aa, nt}')
    sa_arg.add_argument('--diamond', action = 'store_true',
        help = '[blast] Use diamond. Not recommended for ome-by-ome')

    p_arg = parser.add_argument_group('Search parameters')
    p_arg.add_argument('-e', '--evalue', help = 'E value threshold, e.g. ' \
        + '10^(-x) where x is the input', type = int, default = 0)
    p_arg.add_argument('-bit', '--bitscore', default = 0,
        type = int, help = 'Bit score minimum')
    p_arg.add_argument('-i', '--identity', default = 0,
        type = float, help = '[mmseqs|blast] Identity minimum, e.g. 0.6')
    p_arg.add_argument('-p', '--positives', default = 0,
        type = float, help = '[blast] Positives minimum, e.g. 0.6')
    p_arg.add_argument('-m', '--max_hits', type = int, help = 'Max hits for each ome')
    p_arg.add_argument('-qt', '--query_thresh', type = float, default = 0,
        help = 'Query percent hit threshold (+/-), e.g. 0.5')
    p_arg.add_argument('--coordinate', help = 'Extract alignment; DEFAULT: full protein ',
        action = 'store_true')
    p_arg.add_argument('--convert', action = 'store_true', 
        help = '[mmseqs] Convert query IDs to respective database name (profile search)')
    p_arg.add_argument('--iterations', default = 3,
        help = '[mmseqs] --num-iterations in search; DEFAULT: 3')
    p_arg.add_argument('--manual', 
        help = '[mmseqs] Additional search arguments in quotes')
    p_arg.add_argument('--acc', action = 'store_true', default = False, \
        help = '[hmmer] Extract accessions instead of queries (Pfam)' )

    r_arg = parser.add_argument_group('Runtime options')
    r_arg.add_argument('-v', '--verbose', action = 'store_true')
    r_arg.add_argument('-o', '--output')
    r_arg.add_argument('-c', '--cpu', type = int)
    r_arg.add_argument('--ram', help = 'Useful for mmseqs: e.g. 10M or 5G')

    #parser.add_argument( '-c', '--coverage', type = float, help = 'Query coverage +/-, e.g. 0.5' )
#    parser.add_argument('-f', '--force', action = 'store_true', help = 'Force ome-by-ome blast')
    args = parser.parse_args()

    if args.algorithm not in algorithms:
        eprint('\nERROR: {args.algorithm} not implemented', flush = True)
        sys.exit(19)

    # query parsing
    if not args.query and not args.query_dir and not args.query_file:
        eprint('\nERROR: --query, --query_dir, or --query_file not specified', flush = True)
        sys.exit(18)
    elif [args.query, args.query_dir, args.query_file].count(True) > 1:
        eprint('\nERROR: multiple query types specified', flush = True)
        sys.exit(20)
    else:
        if args.query:
            if args.query == '-':
                queries = fa2dict(stdin2str(), file_ = False)
            else:
                queries = [format_path(args.query)]
        elif args.query_dir:
            queries = [format_path(args.query_dir) + x \
                       for x in os.listdir(format_path(args.query_dir))]
        else:
            with open(format_path(args.query_file), 'r') as raw:
                queries = [format_path(x.rstrip()) for x in raw.read().split()]
    deps = [args.algorithm]

    # diamond specific
    if args.diamond:
        if 'blast' not in args.algorithm:
            eprint(f'\nERROR: --diamond incompatible with {args.algorithm}',
                   flush = True)
            sys.exit(4)
        del deps[0]
        deps.append('diamond')

    # mmseqs-specific 
    if args.algorithm == 'mmseqs':
        if not args.seqtype or args.seqtype not in {'aa', 'nt'}:
            eprint('\nERROR: -st required for mmseqs', flush = True)
            sys.exit(3)
        if args.seqtype == 'aa':
            biotype = 'faa'
        else:
            biotype = 'fna'
#        query_set = set(queries)
 #       for query in queries:
  #          if query + '.index' in query_set:
   #             for ext in ['_h', '_h.index', '.lookup', '.source',
    #                        '.dbtype', '_h.dbtype', '.index']:
     #               if query + ext in query_set:
      #                  query_set.remove(query + ext)
       # queries = sorted(query_set)
    else:
        biotype = None
    findExecs(deps, exit = set(deps))

    # identity
    if args.identity:
        if args.algorithm == 'hmmsearch':
            eprint(f'\nERROR: --identity incompatible with hmmsearch', flush = True)
            sys.exit(21)

    if not args.output:
        base = os.getcwd() + '/'
    else:
        base = format_path(args.output, force_dir = True)
    if base:
        output = base
        if not os.path.isdir(output):
            os.mkdir(output)
#            output = mkOutput(base, 'db2search')
    else:
        output = mkOutput(base, 'db2search')

    if args.cpu and args.cpu < mp.cpu_count():
        cpu = args.cpu
    else:
        cpu = mp.cpu_count()

    if args.evalue:
        evalue = 10**( -args.evalue)
    else:
        evalue = 1

    args_dict = { 
        'Algorithm': args.algorithm, 'MycotoolsDB': args.mtdb, 
        'Queries': ','.join(list(queries)), 
        'Output': output, 'Hits/ome': args.max_hits, 
        'Coverage threshold': args.query_thresh, 'Max E-value': evalue, 
        'Min bitscore': args.bitscore, 'Min identity': args.identity,
        'Min positives': args.positives,
        'Coordinate extract': args.coordinate, 'Iterations': args.iterations,
        'HMM accession': args.acc, 'Diamond': args.diamond,
        'Output': output, 'CPUs': cpu
        }
    start_time = intro('db2search', args_dict)
    db = mtdb(format_path(args.mtdb))


    # run hmm mode
    if args.algorithm in {'hmmsearch'}:
        output_fas = hmmer_main(db, queries, output, 
             args.acc, args.max_hits, args.query_thresh,
             evalue = evalue, bitscore = args.bitscore, binary = args.algorithm,
             coords = args.coordinate, verbose = args.verbose, cpu = cpu)    
    # run blast mode
    elif 'blast' in args.algorithm.lower():
        output_fas = blast_main( 
            db, args.algorithm, queries, output, evalue = evalue,
            bitscore = args.bitscore, pident = args.identity,
            max_hits = args.max_hits, cpus = cpu, force = True,
            search_arg = manual_cmd, coordinate = args.coordinate, 
            coverage = args.query_thresh, ppos = args.positives
            )
    # run mmseqs
    else:
        output_fas = mmseqs_main( 
            db, args.algorithm, queries, output, 
            max_hits = args.max_hits, evalue = evalue, 
            bitscore = args.bitscore, pident = args.identity, 
            mem = args.ram, coverage = args.query_thresh,
            cpus = cpu, biotype = biotype,
            skip = [], coordinate = args.coordinate,
            search_arg = manual_cmd, convert = args.convert, iterations = args.iterations
            )
    if not os.path.isdir(output + 'fastas/'):
        os.mkdir(output + 'fastas/')

    for query in output_fas:
        with open(output + 'fastas/' + query + '.search.fa', 'w' ) as out:
            out.write( dict2fa(output_fas[query]) )

    outro(start_time)


if __name__ == '__main__':
    cli()
