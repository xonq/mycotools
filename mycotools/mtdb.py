#! /usr/bin/env python3

# NEED to add a log option
    #list of DBs, ability to change names, quickly change between
    # list update dates
    #report storage information
    #report taxonomy data
# NEED to remove standalone scripts from PATH and just reference mtdb (legacy)
# NEED to add option to export NCBI/JGI credentials
# NEED to pay attention to old ome versions

import os
import re
import sys
import subprocess
from mycotools.lib.kontools import format_path, read_json, write_json, eprint
from mycotools.lib.dbtools import primaryDB, mtdb_connect, mtdb_disconnect, \
    mtdb_initialize, mtdb, loginCheck


def parse_config(mtdb_config_file = format_path('~/.mycotools/config.json')):
    config_dir = format_path('~/.mycotools/')
    if not os.path.isdir(config_dir):
        os.mkdir(config_dir)
        config_dir += '/'

    if os.path.isfile(mtdb_config_file):
        config = read_json(mtdb_config_file)
    else:
        config = {}
    return config

def main(argv = sys.argv):
    description = 'MycotoolsDB (MTDB) utility usage' \
        + '\n\nmtdb: print master MTDB path' \
        + '\n\nmtdb <SUBSCRIPT>:' \
        + '\n[e]xtract\t\textract sub .mtdb file' \
        + '\n[u]pdate\t\tupdate/initialize primary MTDB' \
        + '\n[p]redb2mtdb\t\tadd local genomes to the primary MTDB' \
        + '\n[m]anage\t\tMTDB management utility' \
        + '\n\nmtdb <OME>[.gff3|.fna|.faa]: [PATH] print ome/ome code path' \
        + '\n\nmtdb <ARG>:' \
        + '\n[-i DBPATH]\t[--interface]\tInitialize MTDB connection' \
        + '\n[-f]\t\t[--fungi]     \tConnect to fungal MTDB' \
        + '\n[-p]\t\t[--prokaryote]\tConnect to prokaryote MTDB' \
        + '\n[-u]\t\t[--unlink]\tUnlink from MTDB' \
        + '\n[-d]\t\t[--dependency]\tInstall/update dependencies'

    config = parse_config()
    if any([not x.startswith('-') for x in argv[1:]]):
        script = [x for x in argv[1:] if not x.startswith('-')]
        script = script[0]
        ab2mt = {'extract': 'extract_mtdb', 'update': 'update_mtdb', 'predb2mtdb': 'predb2mtdb',
                 'e': 'extract_mtdb', 'u': 'update_mtdb', 'p': 'predb2mtdb', 'm': 'manage_mtdb',
                 'manage': 'mange_mtdb'}
        if script in ab2mt:
            exit_code = subprocess.call([ab2mt[script]] + [x for x in argv[1:] \
                                         if x != script])
            sys.exit(exit_code)

    if len([x for x in argv if x.startswith('-')]) > 1:
        eprint('\nERROR: one argument allowed.\n' + description, flush = True)
        sys.exit(2)

    set_argv = set(argv)
    if {'-h', '--help'}.intersection(set_argv):
        print('\n' + description + '\n', flush = True)
        sys.exit(0)
       
    if {'-i', '--interface'}.intersection(set_argv):
        if '-i' in set_argv:
            coord = '-i'
        else:
            coord = '--interface'
        try:
            mycodb_loc_prep = argv[argv.index(coord) + 1]
            if isinstance(mycodb_loc_prep, str):
                if mycodb_loc_prep.startswith('-'):
                    raise IndexError
                else:
                    mycodb_loc = format_path(mycodb_loc_prep)
                    if not os.path.isdir(mycodb_loc):
                        raise FileNotFoundError('invalid MTDB path')
                    elif not os.path.isdir(mycodb_loc + 'mtdb'):
                        raise FileNotFoundError('invalid MTDB path')
        except IndexError:
            raise ValueError('--interface requires a path')
        mtdb_initialize(mycodb_loc)
    elif {'-d', '--dependencies'}.intersection(set_argv):
        pip_deps = ['dna_features_viewer', 'mycotools']
        dep_cmds = [['conda', 'install', '-y', '-c', 'jlsteenwyk', 'clipkit'],
                    ['python3', '-m', 'pip', 'install'] + pip_deps + ['--upgrade']]
        for dep_cmd in dep_cmds:
            cmd = subprocess.call(dep_cmd)
            if cmd:
                print('\nUPDATE failed. Exit ' + str(dep_cmd))
                sys.exit(cmd)
    elif {'-f', '--fungi'}.intersection(set_argv):
        if len(set_argv) > 2:
            eprint('\nERROR: -f does not accept additional arguments. Did you mean -i?')
            sys.exit(14)
        if '-f' in set_argv:
            coord = '-f'
        else:
            coord = '--fungi'
        if 'fungi' not in config:
            raise ValueError('fungal MTDB not connected')
        else:
            mtdb_connect(config, 'fungi')
            if not os.path.isfile(format_path('~/.mycotools/mtdb_key')):
                loginCheck()
    elif {'-p', '--prokaryote'}.intersection(set_argv):
        if len(set_argv) > 2:
            eprint('\nERROR: -p does not accept additional arguments. Did you mean -i?')
            sys.exit(15)

        if '-p' in set_argv:
            coord = '-p'
        else:
            coord = '--prokaryote'
        if 'prokaryote' not in config:
            raise ValueError('prokaryote MTDB not connected')
        else:
            mtdb_connect(config, 'prokaryote')
            if not os.path.isfile(format_path('~/.mycotools/mtdb_key')):
                loginCheck()
    elif {'-u', '--unlink'}.intersection(set_argv):
        if '-u' in set_argv:
            coord = '-u'
        else:
            coord = '--unlink'
        mtdb_disconnect(config)
        sys.exit(0)
    elif len(argv) > 1:
        omes = argv[1].replace('"','').replace("'",'').split()
        for ome_prep in omes:
            db = mtdb(primaryDB()).set_index()
            if ome_prep in db:
                print(ome_prep + '\t' \
                    + '\t'.join([str(db[ome_prep][x]) \
                      for x in db[ome_prep]]))
                sys.exit(0)
            ome = re.sub(r'\.\w+[\w\d]$', '', ome_prep)
            extension_srch = re.search(r'^\d+\.?\d*\.(.*$)', ome_prep[6:])
            if extension_srch is not None:
                extension = extension_srch[1]
            else:
                extension = None
            if ome in db:
                try:
                    if extension:
                        print(db[ome][extension])
                    else:
                        print({**{'ome': ome}, **db[ome]})
                except KeyError:
                    raise KeyError('Invalid extension ' + extension)
            else:
                for ref_ome, row in db.items():
                    if ref_ome.startswith(ome + '.'):
                        if extension:
                            print(row[extension])
                        else:
                            print({**{'ome': ome}, **row}) 
                        break
                else:
                    raise KeyError('Invalid ome ' + ome)
        sys.exit(0)


    path = primaryDB()
    if path:
        print(path)
        sys.exit(0)
    else:
        eprint('Link a MycotoolsDB via `mtdb -i <MTDB_DIR>`')
        sys.exit(1)


def cli():
    main(sys.argv)


if __name__ == '__main__':
    cli()
