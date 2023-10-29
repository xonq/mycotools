#! /usr/bin/env python3

# NEED to add info on pseudogenes, tRNAs, account for those in output, base proteins off CDS

import os
import re
import sys
from itertools import chain
from collections import defaultdict
from mycotools.lib.biotools import gff2list, gff3Comps
from mycotools.lib.kontools import format_path, eprint


def compile_alia(gff_path, output, ome = None):
    # Prepare the data structures for population
    gff = gff2list(gff_path)
    prot_dict, exon_dict, mrna_dict, trna_dict, orna_dict, pseudogene_dict, gene_dict = \
        defaultdict(list), defaultdict(list), defaultdict(list), \
        defaultdict(list), defaultdict(list), defaultdict(list), \
        defaultdict(list)

    # iterate through the GFF entries and add them to their respective type
    # dictionary
    for entry in gff:
        try:
            # extract the alias from the attributes field
            alias = re.search(gff3Comps()['Alias'], entry['attributes'])[1]
        except TypeError:
            raise TypeError(f'entry without MTDB alias: {entry}')
        if entry['type'].lower() == 'gene':
            gene_dict[alias].append(sorted((entry['start'], entry['end'])))
        elif entry['type'].lower() == 'cds':
            prot_dict[alias].append(sorted((entry['start'], entry['end'])))
        elif entry['type'].lower() == 'exon':
            exon_dict[alias].append(sorted((entry['start'], entry['end'])))
        elif entry['type'].lower() == 'mrna':
            mrna_dict[alias].append(sorted((entry['start'], entry['end'])))
        elif entry['type'].lower() == 'trna':
            trna_dict[alias].append(sorted((entry['start'], entry['end'])))
        elif 'RNA' in entry['type']:
            orna_dict[alias].append(sorted((entry['start'], entry['end'])))
        elif entry['type'] == 'pseudogene':
            pseudogene_dict[alias].append(sorted((entry['start'], entry['end'])))

    # calculate the lengths of each specific type
    gene_lens = sorted([v[1] - v[0] for v in chain(*list(gene_dict.values()))])
    prot_lens = sorted([v[1] - v[0] for v in chain(*list(prot_dict.values()))])
    exon_lens = sorted([v[1] - v[0] for v in chain(*list(exon_dict.values()))])
    mrna_lens = sorted([v[1] - v[0] for v in chain(*list(mrna_dict.values()))])
    trna_lens = sorted([v[1] - v[0] for v in chain(*list(trna_dict.values()))])
    orna_lens = sorted([v[1] - v[0] for v in chain(*list(orna_dict.values()))])
    pseudogene_lens = sorted([v[1] - v[0] \
                             for v in chain(*list(pseudogene_dict.values()))])

    # calculate the number of each specific type
    gene_len = len(gene_lens)
    prot_len = len(prot_lens)
    exon_len = len(exon_lens)
    mrna_len = len(mrna_lens)
    trna_len = len(trna_lens)
    orna_len = len(orna_lens)
    pseu_len = len(pseudogene_lens)
    med_genes = gene_lens[round(gene_len/2) - 1]

    # if no output, then print data to stdout
    if not output:
        mean_genes = sum(gene_lens)/gene_len
        print('{:<25}'.format('GENES:') + str(gene_len), flush = True)
        print('{:<25}'.format('GENE LENGTH:') + str(sum(gene_lens)) , flush = True)
        print('{:<25}'.format('MEAN GENE LENGTH:') + str(mean_genes), flush = True)
        print('{:<25}'.format('MEDIAN GENE LENGTH:') + str(med_genes), flush = True)
        if prot_lens:
            mean_prots = sum(prot_lens)/prot_len
            print('{:<25}'.format('CDS LENGTH (bases):') + str(sum(prot_lens)) , flush = True)
            print('{:<25}'.format('MEAN CDS LENGTH:') + str(mean_prots), 
                flush = True)
            med_prots = prot_lens[round(prot_len/2) - 1]
            print('{:<25}'.format('MEDIAN CDS LENGTH:') + str(med_prots), flush = True)
        if mrna_lens:
            mean_mrnas = sum(mrna_lens)/mrna_len
            print('{:<25}'.format('mRNA LENGTH:') + str(sum(mrna_lens)) , flush = True)
            print('{:<25}'.format('MEAN mRNA LENGTH:') + str(mean_mrnas), 
                flush = True)
            med_mrnas = mrna_lens[round(mrna_len/2) - 1]
            print('{:<25}'.format('MEDIAN mRNA LENGTH:') + str(med_mrnas), flush = True)
        if trna_lens:
            mean_trnas = sum(trna_lens)/trna_len
            print('{:<25}'.format('tRNA LENGTH:') + str(sum(trna_lens)) , flush = True)
            print('{:<25}'.format('MEAN tRNA LENGTH:') + str(mean_trnas), 
                flush = True)
            med_trnas = trna_lens[round(trna_len/2) - 1]
            print('{:<25}'.format('MEDIAN tRNA LENGTH:') + str(med_trnas), flush = True)
        if orna_lens:
            mean_ornas = sum(orna_lens)/orna_len
            print('{:<25}'.format('OTHER RNA LENGTH:') + str(sum(orna_lens)) , flush = True)
            print('{:<25}'.format('MEAN OTHER RNA LENGTH:') + str(mean_ornas), 
                flush = True)
            med_ornas = orna_lens[round(orna_len/2) - 1]
            print('{:<25}'.format('MEDIAN OTHER RNA LENGTH:') + str(med_ornas), flush = True)
        if pseudogene_lens:
            mean_pseudogenes = sum(pseudogene_lens)/pseu_len
            print('{:<25}'.format('PSEUDOGENE LENGTH:') + str(sum(pseudogene_lens)) , flush = True)
            print('{:<25}'.format('MEAN PSEUDOGENE LENGTH:') + str(mean_pseudogenes), 
                flush = True)
            med_pseudogenes = pseudogene_lens[round(pseu_len/2) - 1]
            print('{:<25}'.format('MEDIAN PSEUDOGENE LENGTH:') + str(med_pseudogenes), flush = True)
    else:
        mean_genes = sum(gene_lens)/gene_len
        if prot_lens:
            mean_prots = sum(prot_lens)/prot_len
            med_prots = prot_lens[round(prot_len/2) - 1]
        if mrna_lens:
            mean_mrnas = sum(mrna_lens)/mrna_len
            med_mrnas = mrna_lens[round(mrna_len/2) - 1]
        if trna_lens:
            mean_trnas = sum(trna_lens)/trna_len
            med_trnas = trna_lens[round(trna_len/2) - 1]
        if orna_lens:
            mean_ornas = sum(orna_lens)/orna_len
            med_ornas = orna_lens[round(orna_len/2) - 1]
        if pseudogene_lens:
            mean_pseudogenes = sum(pseudogene_lens)/pseu_len
            med_pseudogenes = pseudogene_lens[round(pseu_len/2) - 1]

    # populate a dictionary of the statistics for each type
    geneStats = {
        'gene_len': sum(gene_lens), 'genes': gene_len,
        'mean_gene': mean_genes, 'median_gene': med_genes
        }
    if prot_lens:
        geneStats = {**geneStats, **{'cds_len_(bases)': sum(prot_lens), 'cdss': prot_len,
                                   'mean_cds': mean_prots, 'median_cds': med_prots}}
    else:
        geneStats = {**geneStats, **{'cds_len_(bases)': '', 'cdss': '',
                                   'mean_cds': '', 'median_cds': ''}}
    if mrna_lens:
        geneStats = {**geneStats, **{'mrna_len': sum(mrna_lens), 'mrnas': mrna_len,
                                   'mean_mrna': mean_mrnas, 'median_mrna': med_mrnas}}
    else:
        geneStats = {**geneStats, **{'mrna_len': '', 'mrnas': '',
                                   'mean_mrna': '', 'median_mrna': ''}}

    if trna_lens:
        geneStats = {**geneStats, **{'trna_len': sum(trna_lens), 'trnas': trna_len,
                                   'mean_trna': mean_trnas, 'median_trna': med_trnas}}
    else:
        geneStats = {**geneStats, **{'trna_len': '', 'trnas': '',
                                   'mean_trna': '', 'median_trna': ''}}
    if orna_lens:
        geneStats = {**geneStats, **{'other_rna_len': sum(orna_lens), 'ornas': orna_len,
                                   'mean_orna': mean_ornas, 'median_orna': med_ornas}}
    else:
        geneStats = {**geneStats, **{'other_rna_len': '', 'ornas': '',
                                   'mean_orna': '', 'median_orna': ''}}
    if pseudogene_lens:
        geneStats = {**geneStats, **{'pseudogene_len': sum(pseudogene_lens), 
                                   'pseudogenes': pseu_len,
                                   'mean_pseudogene': mean_pseudogenes, 
                                   'median_pseudogene': med_pseudogenes}}
    else:
        geneStats = {**geneStats, **{'pseudogene_len': '', 'pseudogenes': '',
                                   'mean_pseudogene': '', 'median_pseudogene': ''}}


    return ome, geneStats


def main(in_path, log_path = None, cpus = 1, db = None):

    if in_path[-4:] not in {'.gtf', '.gff', 'gff3'} or db:
        from mycotools.lib.dbtools import mtdb
        import multiprocessing as mp
        if not db:
            db = mtdb(in_path).set_index()

        prevOmes = {}
        if log_path and os.path.isfile(log_path):
             with open(log_path, 'r') as raw:
                 for line in raw:
                     if not line.startswith('#'):
                         omeI = line.index('\t')
                         ome = line[:omeI]
                         prevOmes[ome] = line[omeI+1:].rstrip()

        cmds = []
        for ome in db:
            if ome not in prevOmes:
                cmds.append((db[ome]['gff3'], True, ome,))
        with mp.Pool(processes = cpus) as pool:
            res = pool.starmap(compile_alia, cmds)

        outPrep = {}
        for result in res:
            keys = list(result[1].keys())
            if result:
                outPrep[result[0]] = '\t'.join([str(x) for x in result[1].values()])
        out = {
            k: v for k, v in \
            sorted({**outPrep, **prevOmes}.items(), key = lambda x: x[0])
            }


        if not log_path:
#            output_file = os.path.basename(format_path(sys.argv[1])) + '.annStats.tsv'
            print('#ome\t' + "\t".join(keys) + '\n', flush = True)
            for ome in out:
                print(ome + '\t' + out[ome], flush = True)
        else:
            with open(log_path, 'w') as write:
                write.write('#ome\t' + "\t".join(keys) + '\n')
                for ome in out:
                    write.write(ome + '\t' + out[ome] + '\n')
    else:
        ome, geneStats = compile_alia(in_path, log_path)
#        if log_path:
        


def cli():

    output = False
    usage = '\nUSAGE: `gff`/`gtf`/`gff3` OR mycotoolsDB, optional output file\n'
    if '-h ' in sys.argv or '--help' in sys.argv or '-h' == sys.argv[-1]:
        print(usage, flush = True)
        sys.exit(1)
    elif len( sys.argv ) < 2:
        print(usage, flush = True)
        sys.exit(1)
    elif not os.path.isfile( format_path(sys.argv[1]) ):
        print(usage, flush = True)
        sys.exit(1)
    elif len( sys.argv ) > 2:
        log_path = format_path(sys.argv[2])
    else:
        log_path = None

    in_path = format_path(sys.argv[1])
    main(in_path, log_path, os.cpu_count())

    sys.exit(0)


if __name__ == '__main__':
    cli()
