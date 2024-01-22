#! /usr/bin/env python3

import os
import re
import sys
import argparse
import multiprocessing as mp
from itertools import chain
from collections import defaultdict
from mycotools.lib.kontools import eprint, format_path, file2list, stdin2str
from mycotools.lib.dbtools import primaryDB, mtdb
from mycotools.lib.biotools import gff2list, fa2dict, dict2fa, list2gff, gff3Comps
from mycotools.acc2gff import grab_gff_acc

def prep_gff_output(hit_list, gff_path, cpu = 1):
    """Prepare an output file for gffs"""
    gff_list = gff2list(gff_path)
    mp_cmds = []
    for hit in hit_list:
        mp_cmds.append([gff_list, hit])
    with mp.Pool(processes = cpu) as pool:
        gff_list_strs = pool.starmap(
            grab_gff_acc, mp_cmds
            )
    gff_strs = [list2gff(x) for x in gff_list_strs]
    gff_str = '##gff-version 3\n'
    for x in gff_strs:
        for line in x.split('\n'):
            if not line.startswith('#'):
                gff_str += line + '\n'

    return gff_str

def prep_faa_output(hit_list, proteome_path):
    """Prepare and output file for proteomes"""
    proteome_dict = fa2dict(proteome_path)
    clus_fa = {}
    for hit in hit_list:
        clus_fa[hit] = proteome_dict[hit]

    return dict2fa(clus_fa)

def compile_alias_coords(gff_list, accs_list = []):
    """Compile the coordinates of all genes and RNAs of a given gff into a
    dictionary that is accessed through the sequence ID, followed by the
    alias of each sequence"""
    accs_set = set(accs_list)
    alias_comp = re.compile(gff3Comps()['Alias'])

    # gather the coordinates for each RNA and genes without RNAs 
    coord_dict = defaultdict(lambda: defaultdict(list))
    gene2alia = defaultdict(lambda: defaultdict(list))
    alia2seqid = {}
    for entry in gff_list:
        seqid = entry['seqid']
        if 'gene' in entry['type']:
            # account for alternately spliced aliases
            alias = alias_comp.search(entry['attributes'])[1]
            for al in alias.split('|'):
                gene2alia[seqid][al].extend([entry['start'], entry['end']])
                alia2seqid[al] = seqid
        elif 'RNA' in entry['type']:
            alias = alias_comp.search(entry['attributes'])[1]
            coord_dict[seqid][alias].extend([entry['start'], entry['end']])
           
    # are there genes without RNA?
    gene_keys = set(k for k in chain(*list(gene2alia.values())))
    rna_keys = set(k for k in chain(*list(coord_dict.values())))
    missing = gene_keys.difference(rna_keys)
    for key in list(missing):
        seqid = alia2seqid[key]
        coord_dict[seqid][key] = gene2alia[seqid][key]

    # compile output structures for the queried sequences
    acc2seqid = {}
    for seqid in coord_dict:
        for alias in coord_dict[seqid]:
            coord_dict[seqid][alias].sort()
            if alias in accs_set:
                acc2seqid[alias] = seqid
    #    coord_dict[seqid] = sorted(coord_dict[seqid].keys(), key = lambda k: coord_dict[seqid][k][0]])
        coord_dict[seqid] = {
            k: v for k, v in sorted(
                coord_dict[seqid].items(), key = lambda x: x[1][0]
                )
            }

    return coord_dict, acc2seqid


def prep_outputXgene(coords_dict, acc, plusminus):
    """Prep the output for each gene accession using the coordinates
    dictionary as a sorting mechanism, and return the list of accession names"""
    alias_list = list(coords_dict.keys())
    index = alias_list.index(acc)
    if index - plusminus < 0:
        lower = 0
    else:
        lower = int(index - plusminus)
    upper = int(index + plusminus) + 1
    out_index = alias_list[lower:upper]

    return out_index

def prep_outputXbase(coords_dict, acc, plusminus):
    """Prep the output based on the coordinates of a list of accessions if they
    are within the range of the bases alotted, provided by plusminus"""
    alias_list = list(coords_dict.keys())
    index = alias_list.index(acc)
    start, end = coords_dict[acc][0], coords_dict[acc][1]
    low_bound, high_bound = start - plusminus, end + plusminus
    # acquire the indices of the accessions that fit the boundary
    for i1, prot in enumerate(alias_list):
        coords = coords_dict[prot]
        high = coords[-1]
        if high >= low_bound:
            break
    for rev_i2, prot in enumerate(alias_list[::-1]):
        coords = coords_dict[prot]
        low = coords[0]
        if low <= high_bound:
            break
    i2 = len(alias_list) - rev_i2
    out_index = alias_list[i1:i2]

    return out_index

def grab_between(coords_dict, accs):
    """Grab the accessions that are between two particular accessions by
    accessing their indices"""
    alias_list = list(coords_dict.keys())
    a0i = alias_list.index(accs[0])
    a1i = alias_list.index(accs[1])
    mini = min(a0i, a1i)
    maxi = max(a0i, a1i)
    out_index = alias_list[mini:maxi+1]
    return out_index

def main(gff_list, accs, plusminus = 10, mycotools = False,
         geneGff = False, nt = False, between = False):
    """Input a GFF data structure, and accessions, then return a dictionary set
    for each sequence ID to key the accessions that meet the coordinate
    extraction parameters"""
    out_indices = {}
    # compile the coordinates of all genes and RNAs that hit an accession
    coords_dict, acc2seqid = compile_alias_coords(gff_list, accs)
    if between: # if looking for accessions between a set of accessions
        seqid = list(acc2seqid.values())[0]
        out_indices[accs[0]] = grab_between(coords_dict[seqid], accs)
    elif nt: # if looking for accessions that are +/- a number of nucleotides
        for acc, seqid in acc2seqid.items():
            out_indices[acc] = prep_outputXbase(coords_dict[seqid], acc,
                                                plusminus)
    else: # if looking for accessions that are +/- a number of accessions
        for acc, seqid in acc2seqid.items():
            out_indices[acc] = prep_outputXgene(coords_dict[seqid], acc, 
                                                plusminus)

    out_indices = {k: v for k, v in out_indices.items() if v}
    if geneGff: # if a gff of the RNA entries is desired
        geneGffs_prep = {acc: {} for acc in out_indices}
        alt_geneGffs_prep = {acc: {} for acc in out_indices}
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
            elif 'gene' in entry['type']:
                try:
                    gene = re.search(gff3Comps()['Alias'],
                                     entry['attributes'])[1]
                    for acc, genes in gene_sets.items():
                        if gene in genes:
                            alt_geneGffs_prep[acc][gene] = entry
                except TypeError: # no alias
                    pass

        geneGffs, todel = {}, []
        for acc in out_indices:
            geneGffs[acc] = []
            for gene in out_indices[acc]:
                try:
                    geneGffs[acc].append(geneGffs_prep[acc][gene])
                except KeyError: # gene without rna
                    try:
                        geneGffs[acc].append(alt_geneGffs_prep[acc][gene])
                    except KeyError: # try to grab a gene
                        raise KeyError('gene without RNA/gene entry: ' + gene)
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
    parser.add_argument('-d', '--mtdb', default = primaryDB(), 
        help = 'MTDB; DEFAULT: primary')
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
                gff_str = prep_gff_output( 
                    out_indices[acc], gff, cpu = args.cpu
                )
                with open(acc + '.locus.gff3', 'w') as out:
                    out.write(gff_str)
                if prot:
                    prot_str = prep_faa_output(
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
