#! /usr/bin/env python3

import os
import re
import sys
from mycotools.lib.kontools import format_path, sys_start, eprint
from mycotools.lib.dbtools import primaryDB, mtdb

forbidden = [
    '/', '\\', '[', ']', '|', '+', '=', '(', ')',
    '{', '}', ':', ';', '<', '>', '?', '&', '*', '^',
    '%', '@', '#', '$', '~'
    ]

def parseArgs(args):

    genus, species, strain, ome_code, goOn, alternative = True, True, True, \
        True, True, True
    allowable = {
        '_', '-', '!',
        '`', ',', '.',
        '~', "'", '"'
        }

    for arg in args[2:]:
        if len( arg ) > 1:
            if arg[0] in { '"', "'" }:
                arg = arg[1:]
                if arg[-1] in {'"', "'"}:
                    arg = arg[:-1]
        if os.path.isfile( format_path(arg) ):
            if not goOn:
                eprint( '\nERROR: multiple files' , flush = True)
            db = mtdb( arg )
            goOn = False
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

    if goOn:
        db = mtdb( primaryDB() )

    with open( args[1], 'r' ) as raw_input:
        data_input = raw_input.read()

    return db, data_input, genus, species, strain, ome_code, alternative

def main( db, data_input, genus, species, strain, ome_code, alternative ):

    db = db.set_index()
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
        name = name[:-1]
        for i in forbidden:
            name = name.replace( i, '' )
        data_input = re.sub( ome, name, data_input )

    return data_input.rstrip()


def cli():

    usage = 'USAGE: ome2name.py <INPUTFILE> | ome2name.py <INPUTFILE>' \
        + ' [MYCODB] asvg*&\nDEFAULTS: master db, see script for default' \
        + ' forbidden characters' + \
        '\nInput file to regex sub omes with their name.\n' + \
        'optional mycotools db, string of forbidden characters\n' + \
        '"o" no ome | "g" no genus | "s" no species | "v" no strain' + \
        ' | "a" no source accession'
    args = sys_start( sys.argv, usage, 2, files = [sys.argv[1]] )
    db, data_input, g, sp, st, ome_code, alt = parseArgs(args)
    data_output = main(db, data_input, g, sp, st, ome_code, alt)
    print( data_output , flush = True)
    sys.exit( 0 )


if __name__ == '__main__':
    cli()
