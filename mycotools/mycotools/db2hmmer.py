#! /usr/bin/env python3

# NEED option to fail upon any failures
# NEED log
# NEED nhmmer option
# NEED to .tmp and move files
# NEED to add alignment extraction only option

import os
import re
import sys
import argparse
import subprocess
import multiprocessing as mp
from io import StringIO
from collections import defaultdict
from mycotools.lib.kontools import intro, outro, collect_files, multisub, \
    findExecs, untardir, eprint, format_path, mkOutput, tardir
from mycotools.lib.dbtools import masterDB, mtdb
from mycotools.lib.biotools import dict2fa
#from mycotools.extractHmmsearch import main as exHmm
from mycotools.acc2fa import dbmain as acc2fa
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
        hmm_data = exHmm(data, args[0], args[1], args[2], args[3], header = False)
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

    
def compileacc2fa(db, output, q_dict, coords = True):

    cmd_tuples = []
    if not os.path.isdir(output):
        os.mkdir(output)
    fa_dict = {ome: row['faa'] for ome, row in db.items()}
    for q in q_dict:
        cmd_tuples.append((db, q_dict[q], 
                           output + '/' + q + '.faa', coords,))

    return cmd_tuples


def runacc2fa(db, q_files, output, coords = False):

    accs = []
    if coords:
        for ome, hit_data in q_files.items():
            t_accs = [f'{x[0]}[{x[1]}:{x[2]}]' for x in hit_data]
            accs.extend(t_accs)
    else:
        for ome, hit_data in q_files.items():
            t_accs = [x[0] for x in hit_data]
            accs.extend(t_accs)

    fa_str = dict2fa(acc2fa(db, accs))
    
    with open( output, 'w' ) as out:
        out.write( fa_str )


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


def main(db, hmm_path, output, accessions, max_hits, query_cov,
         coords = False, align = False, evalue = 0.01, 
         binary = 'hmmsearch', trim_mod = '', verbose = False):

    ome_dir = output + 'omes/'
    aln_dir = output + 'aligns/'
    faa_dir = output + 'fastas/'
    hmm_dir = output + 'hmms/'
    trm_dir = output + 'trimmed/'

    if align == 'hmmalign':
        aln_ext = 'phylip'
    elif align == 'mafft':
        aln_ext = 'fasta'

    db = db.set_index('ome')
    ome_set, skip1, skip2 = set(), False, False
    if os.path.isdir(faa_dir): # is there a previous run?
        print('\nCompiling previous run', flush = True)
        # check if all fastas are made
        fas = collect_files(faa_dir, 'faa') # grab completed fastas
        ranQueries = [os.path.basename(fa).replace('.faa','') for fa in fas]
        with open(hmm_path, 'r') as hmm_raw:
            queries = grabAccs(hmm_raw.read())
        if not set(queries).difference(set(ranQueries)):
            skip1 = True

        if skip1:
            if os.path.isdir(aln_dir):
                # check if all alignments are made
                aligns = collect_files(aln_dir, aln_ext)
                ranQueries = [os.path.basename(x).replace(f'.{aln_ext}','') \
                              for x in aligns]
                if len(set(queries).intersection(set(ranQueries) \
                    )) == len(set(queries)):
                    skip2 = True
            
    if not skip1: # do not skip the first step
        if os.path.isfile(output + 'omes.tar.gz'):
            if not os.path.isdir(ome_dir):
                untardir(output + 'omes.tar.gz')
        # check what reports have been generated
        omes = collect_files(ome_dir, 'out')
        ome_set = set(os.path.basename(x).replace('.out', '') for x in omes)
    else:
        print('\thmmsearch -> hits.faa DONE', flush = True)
        if skip2:
            print('\thits.faa -> phylip DONE', flush = True)

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
        hmmsearch_tuples = compile_hmm_cmd(db, hmm_path, 
                                         ome_dir, 
                                         ome_set = ome_set,
                                         cpu = par_cpus)
        hmmsearch_codes = multisub(hmmsearch_tuples, processes = par_runs,
                                   verbose = verbose)
        for code in hmmsearch_codes:
            if code['code']:
                eprint('\tERROR: ' \
                        + str(code['stdin']),
                        flush = True)

        # extract results
        print('\nExtracting hmmsearch output', flush = True)
        exHmm_args = [accessions, max_hits, query_cov, evalue]
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

        acc2fa_tuples = compileacc2fa(db, faa_dir, q_dict, coords = coords)
        with mp.get_context('spawn').Pool(processes = cpu) as pool:
            pool.starmap(runacc2fa, acc2fa_tuples)        

    if align:
        if not skip2:
            print('\nAligning' , flush = True)

            if not os.path.isdir(aln_dir):
                os.mkdir(aln_dir)

            if align == 'hmmalign':
                if not os.path.isdir(hmm_dir):
                    os.mkdir(hmm_dir)
                with open(hmm_path, 'r' ) as hmm:
                    hmm_dict = absHmm(hmm.read())
                for accession in hmm_dict:
                    with open(hmm_dir + accession + '.hmm', 'w') as out:
                        out.write(hmm_dict[accession])
                hmmAlign_tuples = compile_hmmalign_cmds(output,
                                                    list(hmm_dict.keys()))
                align_codes = multisub(hmmAlign_tuples, processes = cpu - 1,
                                          verbose = verbose)
            elif align == 'mafft':
                mafft_tuples = compile_mafft_cmds(output, faa_dir)
                align_codes = multisub(mafft_tuples, processes = cpu - 1,
                                       shell = True, verbose = verbose)
        
#        print( '\nTrimming alignments' , flush = True)
 #       trimmed = None
  #      if os.path.isdir(trm_dir):
   #         trimmed = collect_files(trm_dir, aln_ext)
    #    else:
     #       os.mkdir(trm_dir)
      #  trim_tuples = compile_trim_cmd(output, mod = trim_mod, 
       #                                trimmed = trimmed, ex = aln_ext)

        #trim_codes = multisub(trim_tuples, processes = cpu - 1,
         #                     verbose = verbose)


if __name__ == '__main__':

    mp.set_start_method('spawn')
    parser = argparse.ArgumentParser( description = 'Runs `hmmsearch` v each proteome in a' \
        + ' `.db`. For each query, extracts results and optionally outputs a compiled fasta' \
        + ', hmmalignment and/or trimmed alignment.' )

    i_arg = parser.add_argument_group('Inputs')
    i_arg.add_argument('-d', '--mtdb', default = masterDB())
    i_arg.add_argument('--hmm', required = True, help = '.hmm database' )
#    parser.add_argument('-p', '--previous', help = 'Previous db2hmmsearch dir. ' + \
 #       'Be wary of incomplete outputs from interrupted runs' )
#    parser.add_argument('-f', '--fasta', action = 'store_true', \
 #       help = 'Compile fastas for each query')

    p_arg = parser.add_argument_group('Search parameters')
    p_arg.add_argument('-b', '--binary', default = 'hmmsearch',
                       help = 'Search binary {"hmmsearch", \
                               "nhmmer"} NONFUNCTIONAL')
    p_arg.add_argument('-m', '--max_hits', type = int, help = 'Max hits for each organism')
    p_arg.add_argument('-q', '--query_thresh', type = float,
        help = 'Query percent hit threshold (+/-)')
    p_arg.add_argument('-e', '--evalue', help = 'E value threshold, e.g. ' \
        + '10^(-x) where x is the input', type = int)
    p_arg.add_argument('--coordinate', help = 'Extract alignment coordinates',
        action = 'store_true')
    p_arg.add_argument('-a', '--accession', action = 'store_true', default = False, \
        help = 'Extract accessions instead of queries (Pfam, etc)' )

    a_arg = parser.add_argument_group('Alignment parameters')
    a_arg.add_argument('--hmmalign', default = False, action = 'store_true', \
        help = 'Align via hmmalign')
    a_arg.add_argument('--mafft', action = 'store_true', 
        help = 'Align via mafft')

    r_arg = parser.add_argument_group('Runtime options')
    r_arg.add_argument('-v', '--verbose', action = 'store_true')
    r_arg.add_argument('-o', '--output')
    r_arg.add_argument('-c', '--cpu', type = int)

    args = parser.parse_args()
    deps = ['hmmsearch']
    if args.hmmalign and args.mafft:
        eprint('\nERROR: --hmmalign OR --mafft permitted', flush = True)
    elif args.mafft:
         deps.extend(['mafft'])
         align = 'mafft'
    elif args.hmmalign:
         deps.extend(['hmmalign', 'esl-reformat'])
         align = 'hmmalign'
    else:
        align = False
    findExecs(deps, exit = set(deps))

    if not args.output:
        output = mkOutput(os.getcwd() + '/', 'db2hmmer')
    else:
        output = format_path(args.output)

    if args.cpu and args.cpu < mp.cpu_count():
        cpu = args.cpu
    else:
        cpu = mp.cpu_count()

    args_dict = { 
        'MycotoolsDB': args.mtdb, 'Hmm database': args.hmm, 'Output': output,
        'Align': align, 'Use accession': args.accession, 'Hits/ome': args.max_hits, 
        'Coverage threshold': args.query_thresh, 'E-value': args.evalue, 
        'Coordinate extract': args.coordinate,
        'Output': output, 'CPUs': cpu
        }
    start_time = intro('db2hmmer', args_dict)
    db = mtdb(args.mtdb)

    main(db, format_path(args.hmm), output, 
         args.accession, args.max_hits, args.query_thresh, align = align,
         evalue = args.evalue, binary = args.binary,
         coords = args.coordinate, verbose = args.verbose)

    outro(start_time)
