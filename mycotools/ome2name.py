#! /usr/bin/env python3

import os
import re
import sys
from mycotools.lib.kontools import format_path, sys_start, eprint
from mycotools.lib.dbtools import primaryDB, mtdb


def parse_args(args):
    """Parse the command-line arguments of the script"""

    # metadata is true until otherwise
    genus, species, strain = True, True, True
    ome_code, go_on, alternative = True, True, True,
    allowable = {'_', '-', '!',
                '`', ',', '.',
                '~', "'", '"'}

    arg_index = len(args) - 1
    valid_ranks = {'kingdom', 'subphylum', 'phylum', 'class',
                   'order', 'family'}
    tax, rank = False, None
    for arg in args[2:]:
        if len(arg) > 1:
            # replace beginning and ending quotations
            if arg[0] in {'"', "'"}:
                arg = arg[1:]
                if arg[-1] in {'"', "'"}:
                    arg = arg[:-1]
        if tax:
            rank = arg.lower()
            if rank not in valid_ranks:
                eprint('\nERROR: --taxonomy not in ' \
                     + str(valid_ranks), flush = True)
                sys.exit(10)
            tax = False
        elif '-' in arg:
            if arg in {'-t', '--t', '--taxonomy'}:
                tax = True
                continue
            else:
                eprint('\nERROR: invalid argument `-`', flush = True)
                sys.exit(11)
        elif os.path.isfile(format_path(arg)):
            if not go_on:
                eprint('\nERROR: multiple files' , flush = True)
            db = mtdb(arg)
            go_on = False
            continue
                 
        for let in arg:
            if let == 'g':
                genus = False
            elif let == 's':
                species = False
            elif let == 'v':
                strain = False
            elif let == 'o':
                ome_code = False
            elif let == 'a':
                alternative = False
            elif let in allowable:
                forbidden.append( let )

    # import primary MTDB
    if go_on:
        db = mtdb(primaryDB())

    # read input data
    with open(args[1], 'r') as raw_input:
        data_input = raw_input.read()

    return db, data_input, genus, species, strain, ome_code, alternative, rank


def main(db, data_input, genus = True, species = True, strain = True, 
         ome_code = None, alternative = None, rank = None):
    """Python entry point for converting MTDB genome accessions to taxonomic
    metadata"""

    # forbidden characters 
    forbidden = [
        '/', '\\', '[', ']', '|', '+', '=', '(', ')',
        '{', '}', ':', ';', '<', '>', '?', '&', '*', '^',
        '%', '@', '#', '$', '~'
        ]

    db = db.set_index()
    ome2name = {}
    for ome, row in db.items():
        name = ''
        if genus:
            name += str(row['genus']) + '_'
        if species:
            if row['species']:
                name += str(row['species']) + '_'
            else:
                name += 'sp._'
        if strain:
            if row['strain']:
                name += str(row['strain']) + '_'
        if alternative:
            name += row['assembly_acc'] + '_'
        if ome_code:
            name += ome + '_'
        if rank:
            name += row['taxonomy'][rank]
#        name = name[:-1] + '_'
        # remove forbidden characters
        for i in forbidden:
            name = name.replace(i, '')
        ome2name[ome] = name

    # sort so the longest ome codes are replaced first and thus smaller
    # versions do not replace sub versions of larger
    ome2name = {k: v for k, v in sorted(ome2name.items(),
                                  key = lambda x: len(x[0]),
                                  reverse = True)}
    for k, v in ome2name.items():
        # change the ome code with the name, with the explicit caveat that a digit at the end will not be modified
#        data_input = data_input.replace(k, v)
        data_input = re.sub(k + r'([^\d])', v + r'\1', data_input)
        data_input = re.sub(k + r'$', v, data_input)

    return data_input.rstrip()


def cli():
    """Command line entry point"""
    usage = 'USAGE: ome2name <INPUTFILE> | ome2name.py <INPUTFILE>' \
        + ' [.mtdb] asvg*&\nDEFAULTS: master db, see script for default' \
        + ' forbidden characters' + \
        '\nInput file to regex sub omes with their name.\n' + \
        'optional mycotools db, string of forbidden characters\n' + \
        '"o" no ome | "g" no genus | "s" no species | "v" no strain' + \
        ' | "a" no source accession'
    args = sys_start(sys.argv, usage, 2, files = [sys.argv[1]])
    db, data_input, g, sp, st, ome_code, alt, rank = parse_args(args)
    data_output = main(db, data_input, g, sp, st, ome_code, alt,
                       rank = rank)
    print(data_output, flush = True)
    sys.exit(0)


if __name__ == '__main__':
    cli()
