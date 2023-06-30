#! /usr/bin/env python3

import os
import re
import sys
import argparse
import multiprocessing as mp
from collections import defaultdict
from mycotools.lib.kontools import eprint, format_path, file2list, stdin2str
from mycotools.lib.dbtools import primaryDB, mtdb
from mycotools.lib.biotools import gff2list, fa2dict, dict2fa, list2gff, gff3Comps
from mycotools.acc2gff import grabGffAcc

def prepGffOutput(hit_list, gff_path, cpu = 1):

    gff_list = gff2list(gff_path)
    mp_cmds = []
    for hit in hit_list:
        mp_cmds.append([gff_list, hit])
    with mp.Pool(processes = cpu) as pool:
        gff_list_strs = pool.starmap(
            grabGffAcc, mp_cmds
            )
    gff_strs = [list2gff(x) for x in gff_list_strs]
    gff_str = '##gff-version 3\n'
    for x in gff_strs:
        for line in x.split('\n'):
            if not line.startswith('#'):
                gff_str += line + '\n'

    return gff_str

def prepFaaOutput(hit_list, proteome_path):

    proteome_dict = fa2dict(proteome_path)
    clus_fa = {}
    for hit in hit_list:
        clus_fa[hit] = proteome_dict[hit]

    return dict2fa(clus_fa)

def compileCDS(gff_list, accs_list = []):

    accs_set = set(accs_list)
    gene = False
    for entry in gff_list:
        if entry['type'] == 'CDS':
            break
        else:
            entry = None

    # need to make this predict in a better way
    prot_comps = [ 
    r'Alias=([^;]*)', r';Name=(.*?);', r'alias ([^;]*)', r'protein_id\=(.*?;|.*?$)', 
    r'proteinId (\d+)', r'gene_id "(.*?)"', r'Parent=(.*?)$'
    ]
    for tag in prot_comps:
        if re.search(tag, entry['attributes']):
            prot_comp = re.compile(tag)
            if tag == prot_comps[-1]:
                gene = True
            break

    if gene:
        mRNA_dict = {}
        for entry in [x for x in gff_list if x['type'] == 'mRNA']:
            gene = re.search(r'Alias=([^;]*)', entry['attributes'])[1]
            mRNA = re.search(r'ID=([^;]*)', entry['attributes'])[1]
            mRNA_dict[mRNA] = gene

    cds_dict = {}
    for entry in gff_list:
        if entry['type'] == 'CDS':
            if entry['seqid'] not in cds_dict:
                cds_dict[entry['seqid']] = defaultdict(list)
            prot = prot_comp.search(entry['attributes'])[1]
            if gene:
                prot = mRNA_dict[prot]
            cds_dict[entry['seqid']][prot].extend([entry['start'], entry['end']])
            
    acc2seqid = {}
    for seqid in cds_dict:
        for prot in cds_dict[seqid]:
            cds_dict[seqid][prot].sort()
            if prot in accs_set:
                acc2seqid[prot] = seqid
    #    cds_dict[seqid] = sorted(cds_dict[seqid].keys(), key = lambda k: cds_dict[seqid][k][0]])
        cds_dict[seqid] = {
            k: v for k, v in sorted(
                cds_dict[seqid].items(), key = lambda x: x[1][0]
                )
            }

    return cds_dict, acc2seqid


def compileCDS_mycotools(gff_list, accs_list):
    '''
    Inputs the gff_list and organism ome code. Compiles CDS entries for loci parsing and outputs:
    cds_dict = {contig: {protein: [[CDSstart, CDSstop]] } }
    Sorts the cds_dict for each protein from smallest to largest
    '''

    accs_set = set(accs_list)
    cds_dict, fail = {}, False
    for entry in gff_list:
        if entry['type'] == 'CDS':
            if entry['seqid'] not in cds_dict:
                cds_dict[entry['seqid']] = defaultdict(list)
            prot_prep_i0 = entry['attributes'].index(';Alias=')
            try:
                prot_prep_i1 = entry['attributes'][prot_prep_i0+1:].index(';') + prot_prep_i0+1
            except ValueError:
                prot_prep_i1 = len(entry['attributes'])
            prot = entry['attributes'][prot_prep_i0:prot_prep_i1].replace(';Alias=','')
#            if prot not in cds_dict[entry['seqid']] and prot:
 #               cds_dict[entry['seqid']][prot] = []
            if not prot:
                print(entry['attributes'], prot_prep_i0, prot_prep_i1)
                if not fail:
                    print('\tWARNING: proteins in gff with no Alias', flush = True)
                fail = True
                continue
            cds_dict[entry['seqid']][prot].extend(
                [entry['start'], entry['end']]
                )

    acc2seqid = {}
    for seqid in cds_dict:
        for prot in cds_dict[seqid]:
            if prot in accs_set:
                acc2seqid[prot] = seqid
            cds_dict[seqid][prot].sort()
#        cds_dict[seqid] = list(sorted(cds_dict[seqid].keys(), key = lambda k: cds_dict[seqid][k][0]))
        cds_dict[seqid] = {
            k: v for k, v in sorted(
                cds_dict[seqid].items(), key = lambda x: x[1][0]
                )
            }
    # dict(cds_dict) = {seqid: {prot: [sorted_coords]}}

    return cds_dict, acc2seqid


def prep_outputXgene(prot_dict, acc, plusminus):

    prot_list = list(prot_dict.keys())
    index = prot_list.index(acc)
    if index - plusminus < 0:
        lower = 0
    else:
        lower = int(index - plusminus)
    upper = int(index + plusminus) + 1
    out_index = prot_list[lower:upper]

    return out_index

def prep_outputXbase(prot_dict, acc, plusminus):

    prot_list = list(prot_dict.keys())
    index = prot_list.index(acc)
    start, end = prot_dict[acc][0], prot_dict[acc][1]
    low_bound, high_bound = start - plusminus, end + plusminus
    for i1, prot in enumerate(prot_list):
        coords = prot_dict[prot]
        high = coords[-1]
        if high >= low_bound:
            break
    for rev_i2, prot in enumerate(prot_list[::-1]):
        coords = prot_dict[prot]
        low = coords[0]
        if low <= high_bound:
            break
    i2 = len(prot_list) - rev_i2
    out_index = prot_list[i1:i2]

    return out_index

def grab_between(prot_dict, accs):
    prot_list = list(prot_dict.keys())
    a0i = prot_list.index(accs[0])
    a1i = prot_list.index(accs[1])
    mini = min(a0i, a1i)
    maxi = max(a0i, a1i)
    out_index = prot_list[mini:maxi+1]
    return out_index

def main(gff_list, accs, plusminus = 10, mycotools = False,
         geneGff = False, nt = False, between = False):

    out_indices = {}
    if mycotools:
        cds_dict, acc2seqid = compileCDS_mycotools(gff_list, accs)
    else:
        cds_dict, acc2seqid = compileCDS(gff_list, accs)
    if between:
        seqid = list(acc2seqid.values())[0]
        out_indices[accs[0]] = grab_between(cds_dict[seqid], accs)
    elif nt:
        for acc, seqid in acc2seqid.items():
            out_indices[acc] = prep_outputXbase(cds_dict[seqid], acc,
                                                plusminus)
    else:
        for acc, seqid in acc2seqid.items():
            out_indices[acc] = prep_outputXgene(cds_dict[seqid], acc, 
                                                plusminus)

    out_indices = {k: v for k, v in out_indices.items() if v}
    if geneGff:
        geneGffs_prep = {acc: {} for acc in out_indices}
        gene_sets = {acc: set(genes) for acc, genes in out_indices.items()}
        for entry in gff_list:
            if 'RNA' in entry['type']:
                try:
                    gene = re.search(gff3Comps()['Alias'],
                                     entry['attributes'])[1]
                    for acc, genes in gene_sets.items():
                        if gene in genes:
                            geneGffs_prep[acc][gene] = entry
#                            break # if one gene is in multiple loci it needs to show up
                except TypeError: # no alias
                    pass

        geneGffs, todel = {}, []
        for acc in out_indices:
            geneGffs[acc] = []
            for gene in out_indices[acc]:
                try:
                    geneGffs[acc].append(geneGffs_prep[acc][gene])
                except KeyError: # gene without rna
                    raise KeyError('gene without RNA: ' + gene)
#                    todel.append((acc, gene))
        for acc, gene in todel:
            del out_indices[acc][gene]
        return out_indices, geneGffs
    return out_indices

def mycotools_main(db, accs, plusminus = 10, cpus = 1, nt = False,
                   between = False):

    acc_dict = {}
    for acc in accs:
        ome = acc[:acc.find('_')]
        if ome not in acc_dict:
            acc_dict[ome] = []
        acc_dict[ome].append(acc)

    db = db.set_index('ome')
    cmds = [
        [gff2list(db[ome]['gff3']), accs, plusminus, True, False, nt, between] \
        for ome, accs in acc_dict.items()
        ]
    with mp.Pool(processes = cpus) as pool:
        acc_res = pool.starmap(main, cmds)

    out_indices = {}
    for res in acc_res:
        out_indices = {**out_indices, **res}

    return out_indices

def cli():

    parser = argparse.ArgumentParser(description = 'Extracts loci from acc(s)')
    parser.add_argument('-a', '--acc', help = '"-" for stdin')
    parser.add_argument('-i', '--input', help = 'File of accs')
    parser.add_argument('-b', '--between', help = 'Between two input accs',
                        action = 'store_true')
    parser.add_argument('-n', '--nucleotide', action = 'store_true',
        help = '+/- by base')
    parser.add_argument('-p', '--plusminus', default = 10, type = int,
        help = '+/- from acc; DEFAULT: 10 genes')
    parser.add_argument('-o', '--output', action = 'store_true',
        help = 'Output locus fasta(s) and gff(s)')
    parser.add_argument('-g', '--gff', help = 'Input GFF file')
    parser.add_argument('-f', '--faa', help = 'Input protein fasta file')
    parser.add_argument('-s', '--sep', help = 'Separator for input file.', default = '\n')
    parser.add_argument('-d', '--mtdb', default = primaryDB(), help = 'MTDB; DEFAULT: master')
    parser.add_argument('--cpu', type = int, default = 1)
    args = parser.parse_args()

    if args.cpu < mp.cpu_count():
        cpu = args.cpu
    else:
        cpu = mp.cpu_count()
    args.sep = args.sep.replace("'",'').replace('"','')
 
    if args.input:
        accs = file2list(format_path(args.input), sep = args.sep)
    elif args.acc:
        if args.acc == "-":
            accs = stdin2str().split()
        else:
            if {'"', "'"}.intersection(set(args.acc)):
                args.acc = args.acc.replace('"','').replace("'",'')
            if ',' in args.acc:
                accs = args.acc.split(',')
            elif re.search(r'\s', args.acc):
                accs = args.acc.split()
            else:
                accs = [args.acc]
    else:
        eprint('\nERROR: requires input or acc', flush = True)
        sys.exit(1)

    if args.between:
        if len(accs) > 2:
            eprint('\nERROR: -b needs 2 accessions', flush = True)
            sys.exit(2)
        if args.nucleotide:
            eprint('\nERROR: -b and -n are incompatible', flush = True)
            sys.exit(3)

    db = None
    out_indices = {}
    if args.gff: 
        gff = gff2list(format_path(args.gff))
        out_indices = main(gff, accs, args.plusminus, between = args.between,
                           nt = args.nucleotide)
    else:
        db = mtdb(format_path(args.mtdb)).set_index('ome')
        out_indices = mycotools_main(db, accs, plusminus = args.plusminus, between = args.between,
                                     cpus = args.cpu, nt = args.nucleotide)

    if args.output:
        if not db:
            db = mtdb(format_path(args.mtdb)).set_index('ome')
        for acc in out_indices:
            if args.gff:
                gff = format_path(args.gff)
                if args.faa:
                    prot = format_path(args.faa)
                else:
                    prot = None
            else:
                ome = acc[:acc.find('_')]
                gff, prot = db[ome]['gff3'], db[ome]['faa']
            if len(out_indices[acc]) > 0:
                gff_str = prepGffOutput( 
                    out_indices[acc], gff, cpu = args.cpu
                )
                with open(acc + '.locus.gff3', 'w') as out:
                    out.write(gff_str)
                if prot:
                    prot_str = prepFaaOutput(
                        out_indices[acc], prot
                        )
                    with open(acc + '.locus.faa', 'w') as out:
                        out.write(prot_str)
    else:
        for acc in out_indices:
            for hit in out_indices[acc]:
                print(hit, flush = True)
            print(flush = True)
    sys.exit(0)

if __name__ == '__main__':
    cli()
