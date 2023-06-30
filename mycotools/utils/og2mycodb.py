#! /usr/bin/env python3

import os
import re
import sys
import multiprocessing as mp
from mycotools.lib.biotools import gff2list, list2gff, gff3Comps
from mycotools.lib.kontools import format_path, sys_start

def og2dict(orthogroup_file):

    ome_ogs = {}
    with open(orthogroup_file, 'r') as raw:
        for line in raw:
            data = line.rstrip().split(' ')
            og = int(data[0].replace(':', '').replace('OG',''))
            for gene in data[1:]:
                ome = gene[:gene.index('_')]
                if ome not in ome_ogs:
                    ome_ogs[ome] = {}
                ome_ogs[ome][gene] = og

    return ome_ogs

def sortOGtag(ogtag_dict):
    sort_list = ['K', 'P', 'U', 'C', 'O', 'F', 'G', 'S']

    sorted_dict = {}
    for i in sort_list:
        if i in ogtag_dict:
            sorted_dict[i] = ogtag_dict[i]
            del ogtag_dict[i]
    for i in ogtag_dict:
        sorted_dict[i] = ogtag_dict[i]

    return sorted_dict

def readOGtag(ogtagData):

    ogs = ogtagData.split('|')
    ogtag_dict = {}
    for x in ogs:
        ogtag_info = x.split(':')
        ogtag_dict[ogtag_info[0]] = ogtag_info[1]

    return ogtag_dict

def writeOGtag(ogtag_dict):
    og_str = ''
    for i in ogtag_dict:
        og_str += i + ':' + str(ogtag_dict[i]) + '|'
    return og_str[:-1]

def editOGtag(ogtag_dict, ogtag, og):

    ogtag_dict[ogtag] = og

    new_oginfo = 'OG='
    ogtag_dict = sortOGtag(ogtag_dict)
    for i in ogtag_dict:
        new_oginfo += i + ':' + str(ogtag_dict[i]) + '|'
    new_oginfo = new_oginfo[:-1]

    return new_oginfo

def mycodbOGs(file_path = format_path('$MYCOGFF3/../ogs.tsv'), omes = set()):

    ogInfo_dict = {}
    with open(file_path, 'r') as raw:
        if omes:
            for line in raw:
                og_info = line.rstrip().split('\t')
                gene, ogtag_info = og_info[0], og_info[1]
                if gene[:gene.find('_')] in omes:
                    ogtag_dict = readOGtag(ogtag_info)
                ogInfo_dict[gene] = ogtag_dict
        else:
            for line in raw:
                og_info = line.rstrip().split('\t')
                gene, ogtag_info = og_info[0], og_info[1]
                ogtag_dict = readOGtag(ogtag_info)
                ogInfo_dict[gene] = ogtag_dict

    return ogInfo_dict

def extract_ogs(ogInfo_dict, ogtag):

    og2gene, gene2og = {}, {}
    for gene in ogInfo_dict:
        if ogtag in ogInfo_dict[gene]:
            og = ogInfo_dict[gene][ogtag]
            if og not in og2gene:
                og2gene[og] = []
            og2gene[og].append(gene)
            gene2og[gene] = og

    return og2gene, gene2og

def og2mycoDB(ogInfo_dict, omes = set(), file_path = format_path('$MYCOGFF3/../ogs.tsv')):

    out_list = []
    if os.path.isfile(file_path):
        with open(file_path, 'r') as raw:
            for line in raw:
                data = line.rstrip().split('\t')
                ome = data[0][:data[0].find('_')]
                if ome not in out_list:
                    out_list.append(line.split('\t'))

    for gene in ogInfo_dict:
        out_list.append([gene, writeOGtag(ogInfo_dict[gene])])

    sorted_list = ['\t'.join([str(x) for x in y]) for y in sorted(out_list, key = lambda x: x[0])]
    with open(file_path, 'w') as out:
        out.write('\n'.join(sorted_list))


def og2gff(ogs_dict, gff_path, ogtag):

    gff = gff2list(gff_path)

    for entry in gff:
        if entry['type'] == 'gene':
            gene = re.search(gff3Comps()['Alias'], entry['attributes'])[1]
            if gene in ogs_dict:
                ogSearch = re.search(gff3Comps()['OG'], entry['attributes'])
                if ogSearch:
                    ogtag_dict = readOGtag(ogSearch[1])
                    new_oginfo = editOGtag(ogtag_dict, ogtag, ogs_dict[gene])
                    entry['attributes'] = re.sub(gff3Comps()['OG'], new_oginfo, entry['attributes'])
                else:
                    if not entry['attributes'].endswith(';'):
                        entry['attributes'] += ';'
                    entry['attributes'] += 'OG=' + ogtag + ':' + str(ogs_dict[gene])
                    
    with open(gff_path, 'w') as out:
        out.write(list2gff(gff)) 


def dbMain(og_file, ogtag):
    ome_ogs = og2dict(og_file)
    omes = set(ome_ogs.keys())
    try:
        ogInfo_dict = mycodbOGs(file_path = format_path('$MYCOGFF3/../ogs.tsv'), omes = omes)
    except FileNotFoundError:
        ogInfo_dict = {}
    for ome in ome_ogs:
        for gene in ome_ogs[ome]:
            if gene in ogInfo_dict:
                ogInfo_dict[gene][ogtag] = ome_ogs[ome][gene]
            else:
                ogInfo_dict[gene] = {ogtag: ome_ogs[ome][gene]}

    og2mycoDB(ogInfo_dict, omes = omes, file_path = format_path('$MYCOGFF3/../ogs.tsv'))



def dbMainGff(og_file, ogtag, cpus = 1):

    ome_ogs = og2dict(og_file)

    og2gff_cmds = []
    for ome in ome_ogs:
        og2gff_cmds.append([ome_ogs[ome], format_path('$MYCOGFF3/' + ome + '.gff3'), ogtag])
    with mp.Pool(processes = cpus) as pool:
        pool.starmap(og2gff, og2gff_cmds)


def cli():

    usage = 'Adds orthogroup tag ("OG=") to mycotoolsDB GFFs from an OrthoFinder Orthogroups.txt file\n' \
        + 'og2mycodb.py <Orthogroups.txt> <TAG>\n\n' \
        + 'Default tags: kingdom [K], phylum [P], subphylum [U], class [C], order [O], family [F], genus [G],' \
        + 'species [S]'
    args = sys_start(sys.argv[1:], usage, 2, files = [sys.argv[1]])
    if len(args) > 2:
        cpus = int(args[2])
    else:
        cpus = 1
    dbMain(format_path(args[0]), args[1])
    sys.exit(0)


if __name__ == '__main__':
    cli()
