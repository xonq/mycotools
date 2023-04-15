#! /usr/bin/env python3

import os
import sys
import argparse
from mycotools.lib.dbtools import loginCheck, masterDB, mtdb
from mycotools.lib.kontools import format_path, read_json

# NEED delete database feature

def restrictions(db, restr_list, mtdb_config = format_path('~/.mycotools/config.json'),
                 yes = False):
    mtdb_config = read_json(mtdb_config)
    restr_path = mtdb_config[mtdb_config['active']]['MYCODB'] + '../log/restricted.tsv'

    try:    
        with open(restr_path, 'r') as raw:
            restricted = [x.rstrip().split() for x in raw]
    except FileNotFoundError:
        restricted = []

    accs = set(x[0] for x in restricted)
    for r, s, reason in restr_list:
        if s.lower() in {'ncbi', 'jgi'} and r not in accs:
            restricted.append([r, s.lower(), str(reason)])
            print(r, s, flush = True)

    in_db = [x[0] for x in restricted if x[0] in db]
    while in_db:
        if not yes:
            check = input('Some restrictions are in the MTDB. Delete them? [y/N]: ')
            if check.lower() in {'yes', 'y'}:
                break
            else:
                sys.exit(1)
    

    with open(restr_path, 'w') as out:
        out.write('\n'.join(['\t'.join(x) for x in restricted]))
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'Primary MycotoolsDB management utility')
    parser.add_argument('-p', '--password', action = 'store_true',
                        help = 'Encrypt NCBI/JGI passwords to expedite access')
    parser.add_argument('-r', '--restrict', 
                        help = 'Restrict assembly accessions file, formatted: ' \
                             + '<ACCESSION>\t<SOURCE>\t[REASON]')
    parser.add_argument('-y', '--yes', help = 'Answer yes', action = 'store_true')
    args = parser.parse_args()

    db = mtdb(masterDB()).set_index('assembly_acc')

    if args.password:
        loginCheck()
    if args.restrict:
        restrict_path = format_path(args.restrict)
        with open(restrict_path, 'r') as raw:
            restricted = [x.rstrip().split('\t') for x in raw]
        for v in restricted:
            if len(v) < 3:
                v = v + [None]
        restrictions(db, restricted, yes = args.yes)

    sys.exit(0)
