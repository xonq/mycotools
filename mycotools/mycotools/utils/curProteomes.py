#! /usr/bin/env python3
import re, os, sys, subprocess
from mycotools.lib.kontools import sysStart, formatPath
from mycotools.lib.biotools import dict2fa

# converts JGI proteome hedaers for file from `>jgi|ome|accession|description` to `>ome_accession`
# converts all non letter/number characters in ome code to `.`
def curate(file_, internal_ome):
    '''Import file_, check for jgi or ncbi headers and curate, else return None'''

    jgi, ncbi = False, False
    with open(file_,'r') as fasta_raw:
        fasta = fasta_raw.read()
# are there jgi headers?
        if re.match('^>jgi',fasta):
            headers = re.compile(r'>jgi\|')
            data = headers.sub(r'>',fasta)
            jgi = True
# are there ncbi headers?
        elif re.match('^>\w+\.\d ',fasta):
            ncbi_ome_prep = os.path.basename(file_)
            ncbi_ome_prep2 = re.match(r'(.+)?_proteome.aa.fasta',ncbi_ome_prep)

            headers = re.compile(r'>(\w+\.\d).+')
            data = headers.sub(r'>' + internal_ome + r'_\1',fasta)
            ncbi = True
        else:
            data = None

    if jgi:
        ome_start = re.match(r'^>(.+?)\|',data)[1]
        while not re.match(r'(^>[A-Za-z0-9\.]+)\|',data):
            ome_match = re.compile(r'>([A-Za-z0-9\.]+).(.*?)\|')
            data = ome_match.sub(r'>\1.\2|',data)
        headers = re.compile(r'>([A-Za-z0-9\-.]+)\|([A-Za-z0-9\.]+).*')
        data = headers.sub(r'>' + internal_ome + r'_\2',data)

    return data.rstrip()


if __name__ == '__main__':

    usage = 'Imports jgi/ncbi proteome, ome, and curates headers'
    sysStart(sys.argv, usage, 3, files = [sys.argv[1]])
    new_fa = curate(formatPath(sys.argv[1]), sys.argv[2])
    print(new_fa, flush = True)
