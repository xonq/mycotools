#! /usr/bin/env python3

import os
import sys
import time
import getpass
from Bio import Entrez
from mycotools.lib.kontools import file2list, eprint


def cli():

    usage = "\nInput NCBI accession or new line delimitted file of accessions and optionally the column name.\n"
    usage += 'ncbiAccs2fa.py <ACC> <COLNAME>\n'
    if "-h" in sys.argv or "--help" in sys.argv:
        eprint(usage, flush = True)
        sys.exit( 1 )
    elif len( sys.argv ) < 2:
        eprint(usage, flush = True)
        sys.exit( 1 )
    elif len( sys.argv ) <= 3:
        if os.path.isfile( sys.argv[1] ):
            if len( sys.argv ) == 3:
                accs = file2list( sys.argv[1], sep = '\t', col = sys.argv[2] )
            else:
                accs = file2list( sys.argv[1] )
        else:
            accs = [ sys.argv[1] ]
    else:
        eprint(usage, flush = True)
        sys.exit( 1 )

    email = input( "\nInput NCBI login email: " )
    limit = 3
    Entrez.email = email
    if len( accs ) > 3:
        api = getpass.getpass( prompt = "NCBI API key (leave blank if none): " )
        if api != "":
            Entrez.api_key = api
            limit = 10

    eprint(flush = True)

    count, amount, out_str = 0, 1, ''
    for acc in accs:
        count += 1
        if count >= limit:
            count = 0
            time.sleep( 1 )
        eprint( acc , flush = True)
        while True:
            try:
                handle = Entrez.efetch( db = "protein", id=acc, retmode = "xml" )
                records = Entrez.read(handle)
                out_str += ">" + acc + ' ' + records[0]['GBSeq_organism'].replace(' ','_') + '\n' + records[0]["GBSeq_sequence"].upper() + '\n'
                break
            except:
                time.sleep( 1)

    eprint(flush = True)
    with open( sys.argv[1] + '.retr.fa', 'w' ) as out:
        out.write( out_str )

    sys.exit( 0 )


if __name__ == '__main__':
    cli()
