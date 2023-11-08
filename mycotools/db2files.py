#! /usr/bin/env python3

import os
import re
import sys
import argparse
from datetime import datetime
from shutil import copy as cp
from mycotools.lib.dbtools import primaryDB, mtdb
from mycotools.lib.kontools import format_path, prep_output, eprint, vprint


def soft_main(filetypes, db, output_path, print_link = False, verbose = False):
    """Symlink or print files from each file_type"""

    db = db.set_index('ome')
    # symlink files
    if not print_link:
        # make the directories for each requested file type
        for ftype in filetypes:
            if not os.path.isdir(output_path + ftype):
                os.mkdir(output_path + ftype)
        # grab the files for each genome code
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
    # simply print the link for each file
    else:
        for ome, row in db.items():
            for ftype in filetypes:
                print(row[ftype], flush = True)


def hard_main(filetypes, db, output_path):
    """Hard copy files from filetypes to their filetype output directory"""

    db = db.set_index('ome')
    # create the directories to output each file type
    for ftype in filetypes:
        if not os.path.isdir(output_path + ftype):
            os.mkdir(output_path + ftype)

    # copy each file by genome
    for ome, row in db.items():
        for ftype in filetypes:
            try:
                cp(row[ftype], output_path + ftype + '/' \
                + os.path.basename(row[ftype])) 
            except FileNotFoundError:
                eprint('\tERROR: ' + ome + ' ' + ftype, flush = True)


def mtdb_main(db, output_path, og_mtdb_path):
    """Create a MycotoolsDB directory with the files wanted for copy"""
    
    # generate the base directory for output
    if not output_path:
        output_path = os.getcwd() + '/'
    if not os.path.isdir(output_path):
        os.mkdir(output_path)
    mtdb_dir = output_path + 'mycotoolsdb/'
    if not os.path.isdir(mtdb_dir):
        os.mkdir(mtdb_dir)

    # generate the MTDB hierarchy subdirectories
    sub_dirs = [f'{mtdb_dir}log/', f'{mtdb_dir}config/',
                f'{mtdb_dir}mtdb/', f'{mtdb_dir}data/']
    for dir_ in sub_dirs:
        if not os.path.isdir(dir_):
            os.mkdir(dir_)

    # copy the og_mtdb configuration
    cp(og_mtdb_path + 'config/mtdb.json', f'{mtdb_dir}config/mtdb.json')

    # output the database
    cdate = datetime.now().strftime('%Y%m%d')
    db.df2db(f'{mtdb_dir}mtdb/{cdate}.mtdb')

    # output the files
    hard_main(['gff3', 'faa', 'fna'], db, f'{mtdb_dir}data/')



def cli():

    parser = argparse.ArgumentParser(description = 'Symlinks/copies selected files from database')
    parser.add_argument('-d', '--mtdb', default = primaryDB(), help = 'DEFAULT: primaryDB')
    parser.add_argument('-a', '--assembly', action = 'store_true', help = 'Grab assemblies')
    parser.add_argument('-p', '--proteome', action = 'store_true', 
                        help = 'Grab proteomes')
    parser.add_argument('-g', '--gff', action = 'store_true', help = 'Grab gff`s')
    parser.add_argument('--print', action = 'store_true', help = 'Print paths, no copy')
    parser.add_argument('--hard', action = 'store_true', help = 'Hard copy files')
    parser.add_argument('-n', '--new_mtdb', action = 'store_true', 
        help = 'Create MTDB directory hierarchy')
    parser.add_argument('-o', '--output', default = os.getcwd())
    args = parser.parse_args()

    if not args.assembly and not args.proteome and not args.gff \
        and not args.new_mtdb:
        print('\nERROR: --assembly/--proteome/--gff/--new_mtdb required', 
              flush = True)
        sys.exit(4)
    if args.new_mtdb:
        args.hard = False
        args.print = False

    db_path = format_path(args.mtdb)
    args.output = format_path(args.output, force_dir = True)
    output_path = prep_output(args.output, cd = False)
    args_dict = {
        'DATABASE': db_path, 'OUTPUT': output_path,
        'ASSEMBLY': args.assembly, 'PROTEOME': args.proteome,
        'GFF3': args.gff, 'Print links': args.print, 'Hard copy': args.hard,
        'New MTDB': args.new_mtdb
        }

    filetypes = []
    if args.proteome:
        filetypes.append('faa')
    if args.gff:
        filetypes.append('gff3')
    if args.assembly:
        filetypes.append('fna')

    db = mtdb(db_path)
    if args.new_mtdb:
        mtdb_main(db, output_path, format_path(os.environ['MYCODB'] + '/../'))
    elif args.print or not args.hard:
        soft_main(filetypes, db, output_path, print_link = args.print)
    else:
        hard_main(filetypes, db, output_path)

    sys.exit(0)


if __name__ == '__main__':
    cli()
