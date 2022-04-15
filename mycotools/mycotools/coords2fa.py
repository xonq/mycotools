#! /usr/bin/env python3

from mycotools.lib.biotools import fa2dict, dict2fa
from mycotools.lib.kontools import sysStart, eprint, formatPath
from Bio.Seq import Seq
import argparse, sys, os

def extractCoords(fa_dict, seqid, coord_start = 0, coord_end = -1, sense = '+', fa_name = ''):

    new_fa, error = {}, ''
    name = seqid + '_' + str(coord_start) + '-' + str(coord_end) + '_' + sense
    if sense == '+':
        new_fa[name] = {
            'sequence': fa_dict[seqid]['sequence'][coord_start:coord_end],
            'description': ''
            }
    else:
        new_fa[name] = {
            'sequence': str(Seq(fa_dict[seqid]['sequence'][coord_start:coord_end]).reverse_complement()),
            'description': ''
            }
    if coord_end > len(fa_dict[seqid]['sequence']):
        error = '\nNOTICE: end coordinate beyond contig edge: ' + fa_name + ' ' + seqid + ' ' + str(coord_end)

    return new_fa, error

if __name__ == '__main__':

    usage = 'Input nucleotide fasta/tsv input, extract coordinates\n' \
        + 'coords2fa.py <FA> <SEQID> <START_COORD> <END_COORD> <STRAND_SENSE>' \
        + '\nExtract full sequence from sense strand: coords2fa.py test.fna scaffold_20 0 -1' \
        + '\nExtract coordinates from antisense strand: coords2fa.py test.fna scaffold_20 69 420 -' \
        + '\n\n\nBulk extraction tab delimitted row format:' \
        + '\nfasta\tsequence_id\tstart_coordinate\tend_coordinate\tstrand_sense' \
        +  '\n\ncoords2fa.py coords.tsv'

    args = sysStart(sys.argv[1:], usage, 1, files = [sys.argv[1]])

    out_fa, error = {}, ''
    try:
        fa_file = formatPath(args[0])
        if fa_file.endswith('/'):
            fa_name = os.path.basename(fa_file[:-1])
        else:
            fa_name = os.path.basename(fa_file)

        fa = fa2dict(fa_file)
        if not fa:
            raise IndexError
        if len(args) < 3:
            args.extend([0, -1, '+'])
        elif len(args) < 6:
            args.append('+')
        out_fa, error = extractCoords(
            fa, args[1], min([int(args[2]), int(args[3])]), max([int(args[2]), int(args[3])]), args[4], fa_name
            )

        print(dict2fa(out_fa))
        eprint(error)
        sys.exit(0)

    except IndexError:
        try: 
            with open(args[0], 'r') as raw:
                data = [x.rstrip().split('\t') for x in raw]
            files_data = {}
            for x in data:
                if x[0] not in files_data:
                    files_data[x[0]] = []
                files_data[x[0]].append((x[1], min([int(x[2]), int(x[3])]), max([int(x[2]), int(x[3])]), x[4]))
        except:
            eprint('\nERROR: incorrectly formatted input', flush = True)

    for fa_file, rows in files_data.items():
        fa = fa2dict(formatPath(fa_file))
        if fa_file.endswith('/'):
            fa_name = os.path.basename(fa_file[:-1])
        else:
            fa_name = os.path.basename(fa_file)
        for x in rows:
            new_fa, error_t = extractCoords(fa, x[0], min([int(x[1]), int(x[2])]), max([int(x[1]), int(x[2])]), x[3], fa_name)
            error += error_t
            out_fa = {**out_fa, **{fa_name + '_' + x: new_fa[x] for x in new_fa}}

    print(dict2fa(out_fa))
    eprint(error)
    sys.exit(0)
