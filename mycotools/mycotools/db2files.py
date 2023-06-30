#! /usr/bin/env python3

import os
import re
import sys
import argparse
from shutil import copy as cp
from mycotools.lib.dbtools import primaryDB, mtdb
from mycotools.lib.kontools import format_path, prep_output, eprint, vprint


def soft_main(filetypes, db, output_path, print_link = False, verbose = False):
    """Symlink or print files from each file_type"""

    db = db.set_index('ome')
    if not print_link:
        for ftype in filetypes:
            if not os.path.isdir(output_path + ftype):
                os.mkdir(output_path + ftype)
        for ome, row in db.items():
            for ftype in filetypes:
                if os.path.isfile(row[ftype]):
                    try:
                        os.symlink(row[ftype], output_path + ftype \
                                 + '/' + ome + '.' + ftype)
                    except FileExistsError:
                        vprint('\t' + ome + ' ' + ftype + ' exists',
                               v = verbose, flush = True)
                else:
                    vprint('\tERROR: ' + ome + ' ' + ftype, flush = True,
                           v = verbose)
    else:
        for ome, row in db.items():
            for ftype in filetypes:
                print(row[ftype], flush = True)


def hard_main(filetypes, db, output_path):
    """Hard copy files from filetypes to their filetype output directory"""

    db = db.set_index('ome')
    for ftype in filetypes:
        if not os.path.isdir(output_path + ftype):
            os.mkdir(output_path + ftype)

    for ome, row in db.items():
        for ftype in filetypes:
            try:
                cp(row[ftype], output_path + ftype + '/' \
                + os.path.basename(row[ftype])) 
            except FileNotFoundError:
                eprint('\tERROR: ' + ome + ' ' + ftype, flush = True)

def cli():

    parser = argparse.ArgumentParser( description = 'Symlinks/copies selected files from database' )
    parser.add_argument( '-d', '--mtdb', default = primaryDB(), help = 'DEFAULT: primaryDB' )
    parser.add_argument( '-a', '--assembly', action = 'store_true', help = 'Grab assemblies' )
    parser.add_argument( '-p', '--proteome', action = 'store_true', help = 'Grab proteomes' )
    parser.add_argument( '-g', '--gff', action = 'store_true', help = 'Grab gff`s' )
    parser.add_argument( '--print', action = 'store_true', help = 'Print paths, no copy' )
    parser.add_argument( '--hard', action = 'store_true', help = 'Hard copy files' )
    parser.add_argument( '-o', '--output', default = os.getcwd() )
    args = parser.parse_args()

    if not args.assembly and not args.proteome and not args.gff:
        print('\nERROR: no file type (-a, -g, -p) selected', flush = True)
        sys.exit(4)

    db_path = format_path(args.mtdb)
    args.output = format_path(args.output, force_dir = True)
    output_path = prep_output(args.output, cd = False)
    args_dict = {
        'DATABASE': db_path, 'OUTPUT': output_path,
        'ASSEMBLY': args.assembly, 'PROTEOME': args.proteome,
        'GFF3': args.gff, 'Print links': args.print, 'Hard copy': args.hard
        }

    filetypes = []
    if args.proteome:
        filetypes.append('faa')
    if args.gff:
        filetypes.append('gff3')
    if args.assembly:
        filetypes.append('fna')

    db = mtdb(db_path)
    if args.print or not args.hard:
        soft_main(filetypes, db, output_path, print_link = args.print)
    else:
        hard_main(filetypes, db, output_path)

    sys.exit(0)


if __name__ == '__main__':
    cli()
