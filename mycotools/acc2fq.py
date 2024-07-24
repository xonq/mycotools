#! /usr/bin/env python3

import os
import re
import sys
import gzip
import argparse
from Bio import SeqIO
from collections import defaultdict
#from mycotools.lib.biotools import dict2fq 
from mycotools.lib.kontools import format_path, eprint, stdin2str

def dict2fq(fastq_dict, description = True):
    """Convert a fastq dictionary to a fastq string"""
    fastq_string = ''
    if description:
        for seq in fastq_dict:
            fastq_string += '@' + seq.rstrip()
            if 'description' in fastq_dict[seq]:
                fastq_string += (' ' + fastq_dict[seq]['description']).rstrip() + '\n'
            else:
                fastq_string += '\n'
            fastq_string += fastq_dict[seq]['sequence'].rstrip() + '\n+\n' \
                          + fastq_dict[seq]['score'].rstrip() + '\n'

    else:
        for seq in fastq_dict2:
            fastq_string += '@' + seq.rstrip() + '\n' \
                + fastq_dict[seq]['sequence'].rstrip() + '\n+\n' \
                + fastq_dict[seq]['score'].rstrip() + '\n'

    return fastq_string


def acc2fq(fq_path, accs):
    fq_dict, parse, header, indices, count = {}, False, False, [], 0
    if fq_path.endswith(('.gz', '.gzip')):
        with gzip.open(fq_path, 'rt') as raw:
            for line in raw:
                data = line.rstrip()
                if not parse:
                    if not accs:
                        break
                    if data.startswith('@'):
                        header = data[1:].split(' ')
                        seq_name = header[0]
                        if seq_name in accs:
                            fq_dict[seq_name] = {'sequence': '', 
                                              'description': ' '.join(header[1:]),
                                              'score': ''}
                            parse = True
                            append_seq = 'sequence'
                            accs.remove(seq_name)
                elif parse:
                    if data == '+' and append_seq != 'score':
                        append_seq = 'score'
                    else:
                        fq_dict[seq_name][append_seq] += data
                        seq_len = len(fq_dict[seq_name]['sequence'])
                        scr_len = len(fq_dict[seq_name]['score']) 
                        # the FASTQ is a terribly inefficient format
                        if seq_len <= scr_len:
                            if seq_len == scr_len:
                                parse = False
                            else:
                                raise ValueError('discrepancy between ' \
                                               + 'sequence and score ' \
                                               + f'lengths: {seq_name}')

    else:
        with open(fq_path, 'r') as raw:
            for line in raw:
                data = line.rstrip()
                if not parse:
                    if not accs:
                        break

                    if data.startswith('@'):
                        parse = False
                        header = data[1:].split(' ')
                        seq_name = header[0]
                        if seq_name in accs:
                            fq_dict[seq_name] = {'sequence': '', 
                                              'description': ' '.join(header[1:]),
                                              'score': ''}
                            parse = True
                            append_seq = 'sequence'
                            accs.remove(seq_name)
                elif parse:
                    if data == '+' and append_seq != 'score':
                        append_seq = 'score'
                    else:
                        fq_dict[seq_name][append_seq] += data
                        seq_len = len(fq_dict[seq_name]['sequence'])
                        scr_len = len(fq_dict[seq_name]['score']) 
                        if seq_len <= scr_len:
                            if seq_len == scr_len:
                                parse = False
                            else:
                                raise ValueError('discrepancy between ' \
                                               + 'sequence and score ' \
                                               + f'lengths: {seq_name}')

    return fq_dict

def acc2fq_pe(fq_path, indices):
    fq_dict, parse, count = {}, False, []
    if fq_path.endswith(('.gz', '.gzip')):
        with gzip.open(fq_path, 'rt') as raw:
            for line in raw:
                data = line.rstrip()
                if data.startswith('@'):
                    if not accs:
                        break
                    parse = False
                    header = data[1:].split(' ')
                    seq_name = header[0]
                    if seq_name in accs:
                        fq_dict[seq_name] = {'sequence': '', 
                                          'description': ' '.join(header[1:]),
                                          'score': ''}
                        parse = True
                        append_seq = 'sequence'
                        accs.remove(seq_name)
                elif parse:
                    if data == '+':
                        append_seq = 'score'
                    else:
                        fq_dict[seq_name][append_seq] += data
    else:
        with open(fq_path, 'r') as raw:
            for line in raw:
                data = line.rstrip()
                if data.startswith('@'):
                    if not accs:
                        break
                    parse = False
                    header = data[1:].split(' ')
                    seq_name = header[0]
                    if seq_name in accs:
                        fq_dict[seq_name] = {'sequence': '', 
                                          'description': ' '.join(header[1:]),
                                          'score': ''}
                        parse = True
                        append_seq = 'sequence'
                        accs.remove(seq_name)
                elif parse:
                    if data == '+':
                        append_seq = 'score'
                    else:
                        fq_dict[seq_name][append_seq] += data
    return fq_dict
 
 
            
def cli():

    parser = argparse.ArgumentParser(description = 'Inputs accession, extracts fastq')
    parser.add_argument('-a', '--accession', help = '"-" for stdin. For coordinates ' + \
        'append [$START-$END] - reverse coordinates for antisense')
    parser.add_argument('-i', '--input', help = 'File with accessions')
    parser.add_argument('-f', '--fastq', help = 'FASTQ input', required = True)
    args = parser.parse_args()

    if args.input: # input file
        input_file = format_path(args.input)
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
       
    fq_path = format_path(args.fastq)
    fq_dict = acc2fq(fq_path, set(accs))
    fastq_str = dict2fq(fq_dict)

    print(fastq_str.rstrip() , flush = True)
    sys.exit(0)

if __name__ == '__main__':
    cli()
