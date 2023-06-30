#!/usr/bin/env python3

# NEED to ditch extracthmm and move to simplified output parsing

import os
import re
import sys
import argparse
import datetime
import subprocess
import multiprocessing as mp
from mycotools.extractHmmsearch import main as exHmm, grabNames
from mycotools.extractHmmAcc import main as extr_hmm
from mycotools.db2search import compAcc2fa
from mycotools.acc2fa import famain as acc2fa
from mycotools.lib.kontools import intro, outro, findExecs, eprint, format_path
from mycotools.lib.dbtools import mtdb, primaryDB
from mycotools.lib.biotools import dict2fa

def runextractHmmAcc(hmm, accession, output):

    with open(hmm, 'r') as raw:
        hmm_data = raw.read()
    hmm_str = extr_hmm(hmm_data, accession)
    with open(output, 'w') as out:
        out.write(hmm_str)

    return output


def runHmmer(fasta, hmm, output, cpu = 1, binary = 'hmmsearch'):

    hmm_status = subprocess.call([
        binary, '-o', output, '--cpu', str(cpu), hmm, fasta
        ], stdout = subprocess.PIPE,
        stderr = subprocess.PIPE
        )

    return hmm_status


def run_extract_hmm(
    hmm_out, top_hits, cov_threshold, evalue, query = True, acc = None
    ):

    with open(hmm_out, 'r') as raw:
        data = raw.read()
    accs = grabNames( data, query = query )
    if len(accs) > 1:
        hmm_data = exHmm(
            data, True, top_hits, cov_threshold, evalue, query = query
            )
    else:
        hmm_data = exHmm(
            data, list(accs)[0], top_hits, cov_threshold, evalue, query =
            query
            )
    
    return hmm_data


def parse_hmm_data(hmm_data):

    output_res = {}
    for query in hmm_data:
        output_res[query] = {}
        for v in hmm_data[query][1].split('\n'):
            if v.startswith('#') or not v.rstrip():
                continue
            data = [x.rstrip() for x in v.split('\t')]
            seq = data[0]
            start = data[-5]
            end = data[-4]
            ome = re.search(r'(.*?)\_', seq)[1]
            if ome not in output_res[query]:
                output_res[query][ome] = []
            output_res[query][ome].append([seq, start, end])

    return output_res


def run_acc2fa(db, biotype, output_res, subhit = True, cpu = 1):

    acc2fa_cmds = compAcc2fa(db, biotype, output_res, subhit)
    output_fas = {}
    for query in acc2fa_cmds:
        print('\t' + query, flush = True)
        with mp.get_context('spawn').Pool(processes = cpu) as pool:
            results = pool.starmap(acc2fa, acc2fa_cmds[query])
        output_fas[query] = '\n'.join([dict2fa(x) for x in results])

    return output_fas


def outputFas(output_fas, output_dir, fastaname):

    for query in output_fas:
        with open(output_dir + fastaname + '_' + query + '.fa', 'w') as out:
            out.write(output_fas[query])


def main(
    db, binary, fasta_path, hmm_path, out_dir, accession, 
    top_hits = None, cov_threshold = None, evalue = None,
    cpu = 1, accession_search = False, subhit = True 
    ):

    if binary == 'nhmmer':
        biotype = 'fna'
    else:
        biotype = 'faa'

    if accession:
        print('\nExtracting ' + accession, flush = True)
        hmm_path = runextractHmmAcc(
            hmm_path, accession, out_dir + accession + '.hmm'
            )
        if os.path.isfile(accession):
            accession = []
    else:
        print('\nExtracting accessions', flush = True)
        with open(hmm_path, 'r') as raw:
            hmm_data = raw.read()
        accession = []
   
    if cpu > 1:
        hmm_cpu = cpu - 1
    else:
        hmm_cpu = cpu
    hmmer_out = out_dir + 'hmmer.out'
    print('\nRunning ' + binary, flush = True)
    if runHmmer(fasta_path, hmm_path, hmmer_out, cpu = hmm_cpu, binary = binary):
        eprint('\tERROR: ' + binary + ' failed', flush = True)
        sys.exit(2)
    
    print('\nParsing output', flush = True)
    hmm_data = run_extract_hmm(
        hmmer_out, top_hits, cov_threshold, evalue, 
        not accession_search, accession
        )
    output_res = parse_hmm_data(hmm_data)

    print('\nCompiling fastas', flush = True)
    output_fas = run_acc2fa(db, biotype, output_res, subhit = subhit, cpu = cpu)

    return output_fas


def cli():

    parser = argparse.ArgumentParser( 
        description = 'Runs hmmsearch or nhmmer on a fasta and returns a fasta of hits'
        )
    parser.add_argument('-f', '--fasta', required = True, help = 'Input .fasta')
    parser.add_argument('--hmm', required = True, help = 'Input .hmm')
    parser.add_argument('-b', '--binary', required = True, help = "{'hmmsearch', 'nhmmer'}")
    parser.add_argument('-d', '--mtdb', default = primaryDB(), help = 'MycoDB. DEFAULT: master')
    parser.add_argument('-q', '--query', help = 'Query [acc if -a] from .hmm, or new line delimited file')
    parser.add_argument('-c', '--coverage', type = float, help = 'Decimal minimum percent hmm coverage')
    parser.add_argument('-w', '--whole', action = 'store_true', help = 'Extract entire hit sequence' )
    parser.add_argument('-e', '--evalue', type = int, help = 'Maximum evalue threshold 10^(-x)')
    parser.add_argument(
        '-a', '--accession', action = 'store_true', help = 'Extract accession, not query, from .hmm'
        )
    parser.add_argument('-o', '--output', help = 'Output directory')
    parser.add_argument('--cpu', default = 1, type = int)
    args = parser.parse_args()

    if args.binary not in {'hmmsearch', 'nhmmer'}:
        eprint('\nERROR: invalid --binary', flush = True)
        sys.exit(1)
    findExecs([args.binary], exit = {args.binary})

    if args.output:
        out_dir = args.output
    else:
        date = datetime.datetime.today().strftime('%Y%m%d')
        out_dir = os.getcwd() + '/' + date + '_fa2hmm2fa/'
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)
    out_dir = format_path(out_dir)


    if args.evalue:
        evalue = 10**-(args.evalue)
    else:
        evalue = None

    args_dict = {
        'Fasta': format_path(args.fasta), 'Hmm': format_path(args.hmm), 'Binary': args.binary,
        'MycoDB': format_path(args.mtdb), 'Hmm Acc': args.query, 'Min Coverage': args.coverage,
        'Max E-value': evalue, 'Accessions': args.accession, 'Full Seq': args.whole, 'Output': out_dir, 'CPU': args.cpu
        }
    start_time = intro('fa2hmmer2fa', args_dict)

#    if args.query:
 #       if os.path.isfile(format_path(args.query)):
  #          accs = file2list(args.query)
   #     else:
    #        accs = args.query

    output_fas = main(
        mtdb(format_path(args.mtdb)), args.binary, args.fasta, args.hmm, out_dir, args.query,
        cov_threshold = args.coverage, evalue = args.evalue, cpu = args.cpu,
        accession_search = args.accession, subhit = not args.whole
        )
    fastaname = re.sub(r'\.fa[^\.]*$', '', os.path.basename(os.path.abspath(args.fasta)))
    outputFas(output_fas, out_dir, fastaname)

    outro(start_time)


if __name__ == '__main__':
    cli()
