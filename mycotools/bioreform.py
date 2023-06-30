#! /usr/bin/env python3

import sys
from mycotools.lib.kontools import sys_start, format_path
from Bio import SeqIO

def main(in_file, in_fmt, out_fmt):
    SeqIO.convert(in_file, in_fmt, sys.stdout, out_fmt)

def cli():
    usage = 'bioreform.py <INPUT_FILE> <OUT_FORMAT> [IN_FORMAT]' \
          + '\nIN_FORMAT is required for ambiguous extensions'
    args = sys_start(sys.argv[1:], usage, 2, files = [sys.argv[1]])
    in_file, out_fmt = format_path(args[0]), args[1]
    if len(args) > 2:
        in_fmt = args[2]
    else:
        if in_file.endswith(('.nexus', '.nex', '.nxs')):
            in_fmt = 'nexus'
        elif in_file.endswith(('.phylip', '.phy', '.ph')):
            in_fmt = 'phylip'
        elif in_file.endswith(('.fa', '.fasta', '.fna', '.faa', '.fsa')):
            in_fmt = 'fasta'
        elif in_file.endswith(('.clustal', '.clus')):
            in_fmt = 'clustal'
        elif in_file.endswith(('.stockholm', '.sto', '.stk')):
            in_fmt = 'stockholm'
    
    main(in_file, in_fmt, out_fmt)
    sys.exit(0)


if __name__ == '__main__':
    cli()
