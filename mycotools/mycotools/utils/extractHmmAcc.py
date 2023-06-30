#! /usr/bin/env python3

import os
import re
import sys
from mycotools.lib.kontools import eprint, file2list


def grabAccs(db_str):

    accessions = []
    hmm_search = re.compile( r'HMMER\d\/f \[.*?\]\nNAME(.*?)\nACC +(.*?)\n[^\/]*?\/\/')
    if not hmm_search.search(db_str):
        hmm_all = re.findall( r'HMMER\d\/f \[.*?\]\nNAME +(.*?)' + \
            '\n[^\/]*?\/\/', db_str )
    else:
        hmm_all = hmm_search.findall( db_str )

    for hit in hmm_all:
        accessions.append( hit )

    return accessions


def hmmExtract(accession, db_str):

    hmm_search = re.search( 
        r'HMMER\d\/f \[.*?\].*?\n^NAME.*?\n^ACC   ' + accession + r'[\s\S]*?^\/\/', db_str, re.M
        )
    if hmm_search is None:
        hmm_search = re.search( r'HMMER\d\/f \[.*?\]\nNAME +' + accession + \
            r'[\s\S]*?\/\/', db_str )
        if not hmm_search:
            eprint('\nERROR: ' + accession + ' does not exist or unexpected error\n', flush = True)
            sys.exit( 2 )

    hmm = hmm_search[0] + '\n'

    return hmm


def main(hmm_db, accessions = False):

    if accessions:
        if os.path.isfile( accessions ):
            accessions = file2list( accessions )
            hmm_str = ''
            for accession in accessions:
                hmm_str += hmmExtract( accession, hmm_db )
        else:
            hmm_str = hmmExtract( accessions, hmm_db )

    else:
        accessions = grabAccs( hmm_db )
        hmm_str = {}
        if len(accessions[0]) == 3:
            accessions = [x[2] for x in accessions]
        elif len(accessions[0]) == 2:
            accessions = [x[1] for x in accessions]
        for accession in accessions: 
            hmm_str[ accession ] = hmmExtract( accession, hmm_db )

    return hmm_str


def cli():

    usage = '\nInputs <`.hmm`> database, optional <accession | new line delimitted accessions file>, returns `hmm` accession. If no accession is specified, a folder of all accessions will be created.\n'
    if len( sys.argv ) < 2:
        print( usage , flush = True)
        sys.exit( 1 )
    if '-h' in sys.argv or '--help' in sys.argv:
        print( usage , flush = True)
        sys.exit( 1 )

    if not os.path.isfile( sys.argv[1] ):
        eprint( '\nERROR: Invalid `.hmm` database path' , flush = True)
        sys.exit( 2 )

    print('\nReading hmm database ...', flush = True)
    with open( sys.argv[1], 'r' ) as raw_hmm_db:
        hmm_db = raw_hmm_db.read()

    if len( sys.argv ) > 2:
        hmm_str = main( hmm_db, accessions = sys.argv[2] )
        print('\nWriting abstracted hmms ...', flush = True)
        with open( sys.argv[2] + '.hmm', 'w' ) as out:
            out.write( hmm_str )
    else:
        if not os.path.isdir( 'hmm' ):
            os.mkdir( 'hmm' )
        hmm_strs = main( hmm_db )
        for accession in hmm_strs:
            with open( 'hmm/' + accession + '.hmm', 'w' ) as out:
                out.write( hmm_strs[ accession ] )

    sys.exit( 0 )


if __name__ == '__main__':
    cli()
