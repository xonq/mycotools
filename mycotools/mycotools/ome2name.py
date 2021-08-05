#! /usr/bin/env python3

from mycotools.lib.kontools import formatPath, sysStart, eprint
from mycotools.lib.dbtools import db2df, masterDB
import pandas as pd, re, sys, os



def main( args ):

    genus, species, strain, ome_code, goOn, alternative = True, True, True, \
        True, True, True
    allowable = {
        '_', '-', '!',
        '`', ',', '.',
        '~', "'", '"'
        }
    forbidden = [
        '.', '/', '\\', '[', ']', '|', '+', '=', '(', ')',
        '{', '}', ':', ';', '<', '>', '?', '&', '*', '^',
        '%', '@', '#', '$', '~'
        ]
    for arg in args[2:]:
        if len( arg ) > 1:
            if arg[0] in { '"', "'" }:
                arg = arg[1:]
                if arg[-1] in {'"', "'"}:
                    arg = arg[:-1]
        if os.path.isfile( formatPath(arg) ):
            if not goOn:
                eprint( '\nERROR: multiple files' , flush = True)
            db = db2df( arg )
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
        db = db2df( masterDB() )

    with open( args[1], 'r' ) as raw_input:
        data_input = raw_input.read()

    for i, row in db.iterrows():
        ome = row['internal_ome']
        name = ''
        if ome_code:
            name += ome + '_'
        if genus:
            if not pd.isnull( row['genus'] ):
                name += str(row['genus']) + '_'
        if species:
            if not pd.isnull( row['species'] ):
                name += str(row['species']) + '_'
            else:
                name += 'sp._'
        if strain:
            if not pd.isnull( row['strain'] ):
                name += str(row['strain']) + '_'
        if alternative:
            name += row['genome_code'] + '_'
        name = name[:-1]
        for i in forbidden:
            name = name.replace( i, '' )
        data_input = re.sub( ome, name, data_input )

    return data_input.rstrip()


if __name__ == '__main__':

    usage = 'USAGE: ome2name.py <INPUTFILE> | ome2name.py <INPUTFILE>' \
        + ' [MYCODB] asvg*&\nDEFAULTS: master db, see script for default' \
        + ' forbidden characters' + \
        '\nInput file to regex sub omes with their name.\n' + \
        'optional mycotools db, string of forbidden characters\n' + \
        '"o" no ome | "g" no genus | "s" no species | "v" no strain' + \
        ' | "a" no alternative ome'
    args = sysStart( sys.argv, usage, 2, files = [sys.argv[1]] )
    data_output = main( args )
    print( data_output , flush = True)
    sys.exit( 0 )
    
