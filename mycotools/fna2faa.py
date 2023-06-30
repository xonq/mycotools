#! /usr/bin/env python3

import sys
from Bio.Seq import Seq
from mycotools.lib.kontools import format_path, stdin2str, sys_start
from mycotools.lib.biotools import fa2dict, dict2fa

def cli():
    usage = 'Input nucleotide fasta ("-" for stdin), translate to protein fasta'
    args = sys_start(sys.argv[1:], usage, 1)
    if args[0] == '-':
        data = stdin2str()
    else:
        with open(format_path(args[0]), 'r') as raw:
            data = raw.read()
    if data.startswith('>'):
        fna = fa2dict(data, file_ = False)
    else:
        fna = {'input': {'sequence': data, 'description': ''}}

    faa = {k: {'sequence': Seq(v['sequence']).translate(), 
               'description': v['description']} \
           for k, v in fna.items()}
    print(dict2fa(faa))
    sys.exit(0)


if __name__ == '__main__':
    cli()
