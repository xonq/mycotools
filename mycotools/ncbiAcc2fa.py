#! /usr/bin/env python3

import os
import sys
import time
import getpass
from Bio import Entrez
from mycotools.lib.kontools import file2list, eprint, sys_start


def entrez_login():
    """Login to Entrez from user input"""
    email = input("\nInput NCBI login email: ")
    limit = 3
    Entrez.email = email
    if len(accs) > 3:
        api = getpass.getpass(prompt = "NCBI API key (leave blank if none): ")
        if api != "":
            Entrez.api_key = api
            limit = 10

    eprint(flush = True)
    return limit


def grab_accs(accs, limit):
    """Grab FASTAs of NCBI accessions"""
    count, amount, out_str = 0, 1, ''
    for acc in accs:
        count += 1
        # do not overwhelm the server
        if count >= limit:
            count = 0
            time.sleep(1)
        eprint(acc, flush = True)

        # iteratively query until successful
        attempt = 0
        while attempt < 3:
            attempt += 1
            try:
                handle = Entrez.efetch(db = "protein", id=acc, retmode = "xml")
                records = Entrez.read(handle)
                out_str += f">{acc} " \
                        + f"{records[0]['GBSeq_organism'].replace('','_')}\n" \
                        + f"{records[0]['GBSeq_sequence'].upper()}\n"
                break
            except:
                time.sleep(1)

    return out_str


def cli():
    """Command line entrance"""
    usage = "Input NCBI accession or new line delimitted " \
          + "file of accessions and optionally the column name." \
          + "\nncbiAccs2fa.py <ACC> <COLNAME>\n"

    # parse the arguments
    args = sys_start(sys.argv, usage, 2)

    if len(args) <= 3:
        # import a file of accessions
        if os.path.isfile(args[1]):
            if len(args) == 3:
                accs = file2list(args[1], sep = '\t', col = args[2])
            else:
                accs = file2list(args[1])
        # import the command line accessions
        else:
            accs = [args[1]]

    limit = entrez_login()
    out_str = grab_accs(accs, limit)

    eprint(flush = True)
    with open(args + '.retr.fa', 'w') as out:
        out.write(out_str)

    sys.exit(0)


if __name__ == '__main__':
    cli()
