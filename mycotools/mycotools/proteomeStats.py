#! /usr/bin/env python3

from mycotools.lib.biotools import fa2dict
import os, sys

def protStats( fa_path, output ):

    fa = fa2dict(fa_path)
    prot_lens = []
    total = 0

    for prot in fa:
        length = len(fa[prot]['sequence'].replace('*',''))
        total += int(length)
        prot_lens.append(length)

    prot_lens.sort()
    if len(prot_lens) % 2 == 0:
        median = prot_lens[ int(len(prot_lens)/2 - 1) ]
    else:
        median = (prot_lens[ round(len(prot_lens)/2) - 2 ] + prot_lens[ round(len(prot_lens)/2) - 1 ])/2

    if not output:
        print( '{:<30}'.format('PROTEIN LENGTH:') + str(total) , flush = True)
        print( '{:<30}'.format('PROTEINS:') + str(len(fa)), flush = True)
        print( '{:<30}'.format('MEAN PROTEIN LENGTH:') + str(total/len(fa)), flush = True)
        print( '{:<30}'.format('MEDIAN PROTEIN LENGTH:') + str(median), flush = True)
        protStats = None
    else:
        protStats = [ total, len(fa), total/len(fa), median ]

    return protStats


if __name__ == '__main__':

    output = False
    usage = '\nUSAGE: protein `fasta`, optional output file\n'
    if '-h ' in sys.argv or '--help' in sys.argv or '-h' == sys.argv[-1]:
        print(usage, flush = True)
        sys.exit(1)
    elif len( sys.argv ) < 2:
        print(usage, flush = True)
        sys.exit(1)
    elif not os.path.isfile( sys.argv[1] ):
        print(usage, flush = True)
        sys.exit(1)
    elif len( sys.argv ) > 3:
        output = True

    protStats = protStats( sys.argv[1], output )
    if output:
        out_str = 'total_len\tproteins\tmean_len\tmedian_len\n' + str(protStats[0]) + '\t' + str(protStats[1]) + '\t' \
            + str(protStats[2]) + '\t' + str(protStats[3])

    sys.exit(0)
