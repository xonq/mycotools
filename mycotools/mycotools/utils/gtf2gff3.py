#! /usr/bin/env python3

# NEED TO DEFAULT CHEKC FOR INTRONS FOR PREDB2DB
# NEED to arrive at a consensus for protein IDs
# NEED to name tRNAs, pseudogenes, etc. similar to curGFF.py

import os
import re
import sys
import copy
import argparse
from collections import defaultdict
from mycotools.lib.biotools import gff2list, list2gff, fa2dict, dict2fa, \
    gtfComps, gff3Comps, gff2Comps
from mycotools.lib.kontools import collect_files, eprint, format_path
from mycotools.gff2seq import aamain as gff2proteome
from mycotools.utils.curGFF3 import rename_and_organize


def grabOutput(output_pref):
    '''
    Inputs: orthofiller `output_path` for results
    Outputs: gff_dict and fasta_dict of results
    Collect the files in `output_path`, find the proteome and gtf, and make
    sure there is only one possible entry. Return the dictionary versions of
    both.
    '''

    output_path, pref = os.path.dirname(output_pref), os.path.basename(output_pref)
    output_files = collect_files(output_path, '*')
    hits = [x for x in output_files if os.path.basename(x).startswith(pref)]
    proteome = [x for x in hits if x.endswith('.results.aa.fasta')]
    gtf = [x for x in hits if x.endswith('.results.gtf')]
    if len(gtf) != 1 or len(proteome) != 1:
        if len(gtf) < 1 or len(proteome) < 1:
            print('\nERROR: complete output files not detected' , flush = True)
        else:
            print('\nERROR: multiple eligible complete output detected' , flush = True)
        sys.exit(3)

    return gff2list(gtf[0]), fa2dict(proteome[0])


def intron2exon(gff, gene_comp = re.compile(r'gene_id \"(.*?)\"')):
    '''
    Inputs: gff_dict
    Outputs: gff_dict with introns converted to exons
    For each entry in the gff_dict, search for the gene ID. If the gene is not
    in `gene_info` add the gene as a key and populate a blank intron list, 
    start codon list, stop codon list, strand string, and raw list. If the
    entry type is within the `gene_info` dict for the gene, then append the
    list of start and stop coordinates. Append the entire entry to the raw data
    key. 
    Create a dictionary `intron_genes` for each gene in gene_info if there is
    an intron entry. For each gene in `intron_genes` create a blank list for
    `exon_coords` dict under the key `gene`. If the intron coordinates' start
    codon end coordinate is greater than the start coordinate, then change the
    start codon entry in `intron_genes[gene]` to have the greater value first.
    Repeat for the stop codon. Then sort the intron coordinates of that gene.
    Append the appropriate exon coordinates for the intron based upon strand
    sense.
    For each gene in the gff, if it is not in `exon_coords` then simply append
    to the `new_gff`. Then add genes with new exons.
    '''

    comps = gtfComps()
    gff1, gene_info = [], {}
    for entry in gff:
        try:
            gene = gene_comp.search(entry['attributes'])[1]
        except TypeError:
            if entry['type'] in {'transcript', 'gene'} and entry['source'] == 'AUGUSTUS':
                continue
            gene_comp = re.compile(r'name \"(.*?)\"')
            gene = gene_comp.search(entry['attributes'])[1]
            comps = gff2Comps()
        if int(entry['start']) > int(entry['end']):
            start, end = copy.deepcopy(entry['start']), copy.deepcopy(entry['end'])
            entry['start'], entry['end'] = end, start
        if gene not in gene_info:
            gene_info[gene] = {
                'intron': [], 
                'start_codon': [], 
                'stop_codon': [], 
                'strand': str(entry['strand']), 
                'raw': [] 
            }
        if entry['type'] in gene_info[gene]:
            gene_info[gene][entry['type']].append(
                [int(entry['start']), int(entry['end'])] 
            )
        gene_info[gene]['raw'].append(entry)
        gff1.append(entry)


    intron_genes = {
        gene: gene_info[gene] for gene in gene_info if gene_info[gene]['intron']
        }
    exon_coords = {}
    for gene in intron_genes:
#        if intron_genes[gene]['intron']:
        exon_coords[gene] = []
        intron_genes[gene]['intron'].sort(key = lambda x: x[0])
    
        if intron_genes[gene]['strand'] == '+':
            exon_coords[gene].append([intron_genes[gene]['start_codon'][0][0]])
            for intron in intron_genes[gene]['intron']:
                exon_coords[gene][-1].append(intron[0] - 1)
                exon_coords[gene].append([intron[1] + 1])
            exon_coords[gene][-1].append(intron_genes[gene]['stop_codon'][0][1])

        else:
            exon_coords[gene].append([intron_genes[gene]['stop_codon'][0][0]])
            for intron in intron_genes[gene]['intron']:
                exon_coords[gene][-1].append(intron[0] - 1)
                exon_coords[gene].append([intron[1] + 1])
            exon_coords[gene][-1].append(intron_genes[gene]['start_codon'][0][1])

    gff2 = []
    for entry in gff1:
        gene = gene_comp.search(entry['attributes'])[1]
        if gene not in exon_coords:
            gff2.append(entry)
    for gene in exon_coords:
        add_data = [
            entry for entry in intron_genes[gene]['raw'] if entry['type'] != 'intron' 
        ]
        scaffold = list(add_data)[0]
        for exon in exon_coords[gene]:
            new_entry = dict(scaffold)
            new_entry['type'] = 'exon'
            new_entry['start'] = str(exon[0])
            new_entry['end'] = str(exon[1])
            new_entry['score'] = '.'
            new_entry['phase'] = '.'
            add_data.append(new_entry)
        gff2.extend(add_data)


    return gff2, comps


def curCDS(gff, gene_compile = re.compile(r'gene_id \"(.*?)\"')):

    new_gff, info_dict = [], {}
    for entry in gff:
        if not gene_compile.search(entry['attributes']):
            gene_compile = re.compile(r'ID=(.*?);')
        try:
            gene = gene_compile.search(entry['attributes'])[1]
        except TypeError:
            gene = gene_compile.search(entry['attributes'])[1]
        if gene not in info_dict:
            info_dict[gene] = {'start_codon': [], 'stop_codon': [], 'raw': []}
        if entry['type'] == 'gene':
#            info_dict[gene]['start_codon'].extend(sorted([int(entry['end']) - 2, int(entry['end'])]))
            info_dict[gene]['start_codon'].extend(sorted([int(entry['end']), int(entry['start'])]))
            info_dict[gene]['stop_codon'].extend(sorted([int(entry['end']), int(entry['start'])]))
        elif entry['type'] in {'start_codon', 'stop_codon'}:
            info_dict[gene][entry['type']].extend(sorted([int(entry['start']), int(entry['end'])]))
        info_dict[gene]['raw'].append(copy.deepcopy(entry))

    for gene in info_dict:
        if not info_dict[gene]['start_codon']:
            new_gff.extend(info_dict[gene]['raw'])
            continue
        cop = False
        try:
            start0 = int(info_dict[gene]['start_codon'][0])
            end1 = int(info_dict[gene]['stop_codon'][1])
        except IndexError: # need to flag here
            continue
        if start0 > end1:
            start0 = int(info_dict[gene]['stop_codon'][0])
            end1 = int(info_dict[gene]['start_codon'][1])
        for index, entry in enumerate(info_dict[gene]['raw']):
            if entry['type'] == 'exon':
                cop = False
                break
            elif entry['type'] == 'CDS': # or entry['type'] == 'exon':
                if (int(entry['end']) == end1 and int(entry['start']) == start0):
                    cop = copy.deepcopy(entry)
                    cop['type'] = 'exon'
                    cop['attributes'] = cop['attributes'].replace('cds', 'exon')
        new_gff.extend(info_dict[gene]['raw'])
        if cop:
            new_gff.append(cop)

    return new_gff


def liberalRemoval(gene_dict_prep, contigs):

    check_contigs, failed, gene_dict = {}, [], {}
    for i in contigs:
        check_contigs[i] = [['null', 10000000000000000000], ['null', 0]]
        for gene in contigs[i]:
            if contigs[i][gene]:
                contigs[i][gene].sort()
                if contigs[i][gene][0] < check_contigs[i][0][1]:
                    check_contigs[i][0][1] = contigs[i][gene][0]
                    check_contigs[i][0][0] = gene
                elif contigs[i][gene][-1] > check_contigs[i][1][1]:
                    check_contigs[i][1][1] = contigs[i][gene][-1]
                    check_contigs[i][1][0] = gene


    for gene in gene_dict_prep:
        if gene_dict_prep[gene]['start_codon'] and gene_dict_prep[gene]['stop_codon']:
            gene_dict[gene] = gene_dict_prep[gene]
        else:

            contig = gene_dict_prep[gene]['contig']
            if gene_dict_prep[gene]['strand'] == '+':
                if not gene_dict_prep[gene]['start_codon']:
                    if gene == check_contigs[contig][0][0]:
                        gene_dict_prep[gene]['start_codon'] = [
                            check_contigs[contig][0][1],
                            check_contigs[contig][0][1] + 2
                            ]
                if not gene_dict_prep[gene]['stop_codon']:
                    if gene == check_contigs[contig][1][0]:
                        gene_dict_prep[gene]['stop_codon'] = [
                            check_contigs[contig][1][1] - 2,
                            check_contigs[contig][1][1]
                            ]
            else:
                 if not gene_dict_prep[gene]['start_codon']:
                    if gene == check_contigs[contig][1][0]:
                        gene_dict_prep[gene]['start_codon'] = [
                            check_contigs[contig][1][1] - 2,
                            check_contigs[contig][1][1]
                            ]
                 if not gene_dict_prep[gene]['stop_codon']:
                    if gene == check_contigs[contig][0][0]:
                        gene_dict_prep[gene]['stop_codon'] = [
                            check_contigs[contig][0][1],
                            check_contigs[contig][0][1] + 2
                            ]
            if not gene_dict_prep[gene]['start_codon'] or not gene_dict_prep[gene]['stop_codon']:            
                failed.append([gene, 'no start and/or stop'])
            else:
                gene_dict[gene] = gene_dict_prep[gene]

    return gene_dict, failed


def conservativeRemoval(gene_dict_prep):

    gene_dict, flagged, failed = {}, [], []
    for gene, temp in gene_dict_prep.items():
        rna_keys = []
        for typ in ['exon', 'cds', 'rna', 'start_codon', 'stop_codon']:
            rna_keys.extend(set(temp[typ].keys()).difference({None}))
        rnas = list(set(rna_keys))
        try:
            for tran in rnas:
                temp['exon'][tran].sort()
                if not temp['start_codon'][tran]:
                    if temp['strand'] == '+':
                        temp['start_codon'][tran] = [
                            temp['exon'][tran][0],
                            temp['exon'][tran][0] + 2
                            ]
                    else:
                        temp['start_codon'][tran] = [
                            temp['exon'][tran][-1] - 2,
                            temp['exon'][tran][-1]
                            ]
                if not temp['stop_codon'][tran]:
                    if temp['strand'] == '+':
                        temp['stop_codon'][tran] = [
                            temp['exon'][tran][-1] - 2,
                            temp['exon'][tran][-1]
                            ]
                    else:
                        temp['stop_codon'][tran] = [
                            temp['exon'][tran][0],
                            temp['exon'][tran][0] + 3
                            ]
                flagged.append(gene)
        except IndexError:
            eprint(gene + ' cannot create gene coordinates' , flush = True)
            continue

        gene_dict[gene] = temp
    
    return gene_dict, flagged, failed
        
def fill_transcripts(gene_dict_prep):
    # currently doesnt handle alternate splicing
    for gene, typ_dict in gene_dict_prep.items():
        if 'cds' in typ_dict and 'exon' in typ_dict:
            if {None} == set(typ_dict['cds'].keys()):
                if not None in typ_dict['exon']:
                    if len(typ_dict['exon']) != 1:
                        raise ValueError('cannot reconcile JGI gff2 alternate splicing')
                    else:
                        tran = list(typ_dict['exon'].keys())[0]
                    typ_dict['cds'][tran] = typ_dict['cds'][None]
                    del typ_dict['cds'][None]
    return gene_dict_prep

def add_genes(gtf, safe = True, 
              comps = gtfComps(),
              gene_prefix = 'gene_id'):

    contigs = defaultdict(dict)
    tran_compile = re.compile(comps['transcript'])
    gene_compile = re.compile(comps['id'])
    gene_dict_prep, gene_dict = {}, {}
    for entry in gtf:
        gene = gene_compile.search(entry['attributes'])[1]
        try:
            tran = tran_compile.search(entry['attributes'])[1]
        except TypeError:
            tran = None
        contig = entry['seqid']
        # prepare a dictionary that contains the hierarchy 
        # for each gene
        if gene not in gene_dict_prep:
            gene_dict_prep[gene] = {
                'start_codon': defaultdict(list), 
                'stop_codon': defaultdict(list), 
                'exon': defaultdict(list),
                'cds': defaultdict(list),
                'strand': str(entry['strand']), 'contig': contig,
                'rna': {}
                }

        if entry['type'].lower() in {'start_codon', 'stop_codon'}:
            gene_dict_prep[gene][entry['type'].lower()][tran].extend(
                sorted([int(entry['start']), int(entry['end'])])
               )
        elif entry['type'].lower() in {'cds', 'exon'}:
            gene_dict_prep[gene][entry['type'].lower()][tran].extend(
                sorted([int(entry['start']), int(entry['end'])])
               )
        elif 'RNA' in entry['type']:
            gene_dict_prep[gene]['rna'][tran] = entry

    gene_dict_prep = fill_transcripts(gene_dict_prep)
    gene_dict, flagged, failed = conservativeRemoval(gene_dict_prep)

    check, insert_list = set(), []
    for index, entry in enumerate(gtf):
        gene = gene_compile.search(entry['attributes'])[1]
        if gene not in check and gene in gene_dict:
            # make a new entry for this gene
            new_entry = dict(entry)
            new_entry['type'] = 'gene'
            new_entry['phase'], new_entry['score'] = '.', '.'
            new_entry['attributes'] = f'{gene_prefix} "' + gene + '";'

            # grab all start and stop coordinates for entries
            start_stop_coords = []
            for rna, coords in gene_dict[gene]['start_codon'].items():
                start_stop_coords.extend(coords)
            for rna, coords in gene_dict[gene]['stop_codon'].items():
                start_stop_coords.extend(coords)
            for rna, coords in gene_dict[gene]['exon'].items():
                start_stop_coords.extend(coords)
 
            # use the minimum and maximum start stops to define the gene
            new_entry['start'] = min(start_stop_coords)
            new_entry['end'] = max(start_stop_coords)


            # for each existing entry in the rnas of gene_dict, append a new
            if gene_dict[gene]['rna']:
                for rna, rna_entry in gene_dict[gene]['rna']:
                    insert_list.extend([[index, rna_entry], [index, new_entry]])
            # we must infer the coordinates from the exons
            else:
                # for each transcript id in the exon information
                for rna, coords in gene_dict[gene]['exon'].items():
                    # if there aren't CDS then its an ambiguous RNA (tRNA/rRNA)
                    if not gene_dict[gene]['cds'][rna]:
                        rna_type = 'RNA' # usually tRNA, but its not a safe assumption
                        # from the data we are using
                    else:
                        rna_type = 'mRNA'
                    new2 = copy.deepcopy(new_entry)
                    coords.sort()
                    new2['type'] = rna_type
                    # define the RNA from the minimum and maximum coordintes
                    new2['start'] = min(coords)
                    new2['end'] = max(coords)
                    # give it a new entry
                    if new2['attributes'].endswith(';'):
                        new2['attributes'] += f' transcript_id "{rna}";'
                    else:
                        new2['attributes'] += f'; transcript_id "{rna}";'
                    insert_list.append([index, new2])
                # add one gene entry
                insert_list.append([index, new_entry])
            check.add(gene)

    insert_list.sort(key = lambda x: x[0], reverse = True)

    for insert in insert_list:
        gtf.insert(insert[0], insert[1])

    return gtf, failed, flagged


def remove_start_stop(gtf):
    return [x for x in gtf if x['type'] not in {'start_codon', 'stop_codon'}]


def curate(gff, prefix = None, failed = set(), comps = gtfComps()):

    failed = set(failed)
    gene_comp = re.compile(comps['id'])
    tran_comp = re.compile(comps['transcript'])

    #acquire the start and end coordinates for all entries with the same gene attribute
    temp_exon_dict, rna_info = defaultdict(dict), {}
    tran2gene, gene2rnas = {}, defaultdict(list)
     
    for line in gff:
        if line['type'] in {'exon', 'start_codon', 'stop_codon', 'CDS'} \
            or 'RNA' in line['type']:
            try:
                gene_id = gene_comp.search(line['attributes'])[1]
            except TypeError:
                gene_id = re.search(gtfComps()['id'], line['attributes'])[1]
            try:
                tran_id = tran_comp.search(line['attributes'])[1]
            except TypeError:
                print(line)
                raise TypeError
            if tran_id not in tran2gene:
                tran2gene[tran_id] = gene_id
            elif not tran2gene[tran_id] == gene_id: # transcript points to
            # multiple genes?
                failed = failed.add(tran_id)
                failed = failed.add(gene_id)
            gene2rnas[gene_id].append(tran_id)
            if tran_id not in temp_exon_dict[line['seqid']]:
                temp_exon_dict[line['seqid']][tran_id] = []
            if 'RNA' in line['type']:
                rna_info[tran_id] = line['type']
            temp_exon_dict[line['seqid']][tran_id].extend([int(line['start']), int(line['end'])])

    exon_dict = {k: v for k, v in sorted(temp_exon_dict.items())}

    count = 1
    alias_dict = {}
    for seqid, tran_hits in exon_dict.items():
        for trans_exons in tran_hits.values():
            trans_exons.sort()
        exon_dict[seqid] = dict(sorted(tran_hits.items(), key = lambda x: x[1][0]))
        if prefix:
            for trans in exon_dict[seqid]:
                if trans not in alias_dict:
                    alias_dict[trans] = prefix + '_' + str(count)
                count += 1

    crudesortGff = preSortGFF(gff, gene_comp)

    newGff, exon_check, cds_check = [], defaultdict(int), defaultdict(int)
    for entry in crudesortGff:
        gene = gene_comp.search(entry['attributes'])[1]
        if gene in failed or trans in failed:
            continue
        if 'RNA' in entry['type']:
            trans = tran_comp.search(entry['attributes'])[1]
            count = 1
            if entry['type'] == 'mRNA':
                entry['attributes'] = 'ID=' + trans + ';Parent=' + gene
                if prefix:
                    entry['attributes'] += f';Alias={alias_dict[trans]}'
            else:
                entry['attributes'] = 'ID=' + trans + ';Parent=' + gene
                if prefix:
                    entry['attributes'] += f';Alias={alias_dict[trans]}'
        elif entry['type'] == 'exon':
            trans = tran_comp.search(entry['attributes'])[1]
            exon = exon_check[trans] + 1
            entry['attributes'] = f'ID={trans}.exon{exon};Parent=' \
                                + f'{trans}'
            if prefix:
                entry['attributes'] += f';Alias={alias_dict[trans]}'
            exon_check[trans] += 1
        elif entry['type'] == 'CDS':
            alias = alias_dict[trans]
            trans = tran_comp.search(entry['attributes'])[1]
            cds = cds_check[trans] + 1
            entry['attributes'] = f'ID={trans}.cds{cds};Parent='
            if prefix:
                entry['attributes'] += f'{trans};Alias={alias_dict[trans]}'
            cds_check[trans] += 1
        elif entry['type'] == 'gene' and prefix: # what about other types of entries?
            alia = sorted({alias_dict[trans] for trans in gene2rnas[gene]})
            entry['attributes'] = 'ID=' + gene + ';Alias=' + '|'.join(alia)

        if int(entry['end']) < int(entry['start']):
            start = copy.deepcopy(entry['end'])
            end = copy.deepcopy(entry['start'])
            entry['start'], entry['end'] = start, end
        newGff.append(entry)

    translation_str = ''
    for entry, translation in alias_dict.items():
        translation_str += f'{entry}\t{translation}\n'

    return newGff, translation_str


def sortGene(sorting_group):

    out_group = []
    for entryType in ['gene', 'mrna', 'trna', 'rrna', 'exon', 'cds']:
        todel = []
        for i, entry in enumerate(sorting_group):
            if entry['type'].lower() == entryType:
                out_group.append(entry)
                todel.append(i)
        for i in reversed(todel):
            sorting_group.pop(i)

    for entry in sorting_group:
        out_group.append(entry)

    return out_group 


def sortContig(contigData):

    coordinates = {}
    for gene in contigData:
        t_coords = []
        for entry in contigData[gene]:
            t_coords.extend([int(entry['start']), int(entry['end'])])
        coordinates[gene] = min(t_coords)
    coordinates = {
        k: v for k, v in sorted(coordinates.items(), key = lambda item: item[1])
        }
    outContig = []
    for gene in coordinates:
        outContig.extend(contigData[gene])

    return outContig

def preSortGFF(unsorted_gff, idComp):

    sorting_groups, oldGene = {}, None
    for i, entry in enumerate(unsorted_gff):
        seqid = entry['seqid']
        if seqid not in sorting_groups:
            sorting_groups[seqid] = {}
        gene = idComp.search(entry['attributes'])[1]
#        gene = re.sub(r'(.*?_\d+).*', r'\1', gene)
        if gene not in sorting_groups[seqid]:
            sorting_groups[seqid][gene] = []
        sorting_groups[seqid][gene].append(entry)

    sortedGff = []
    for seqid in sorting_groups:
        contigData = {}
        for gene in sorting_groups[seqid]:
            contigData[gene] = sortGene(sorting_groups[seqid][gene])
        sortedGff.extend(sortContig(contigData))

    return sortedGff


def sortGFF(unsorted_gff, idComp):

    sorting_groups, oldGene = {}, None
    for i, entry in enumerate(unsorted_gff):
        seqid = entry['seqid']
        if seqid not in sorting_groups:
            sorting_groups[seqid] = {}
        gene = idComp.search(entry['attributes'])[1]
        gene = re.sub(r'(.*?_\d+).*', r'\1', gene)
        if gene not in sorting_groups[seqid]:
            sorting_groups[seqid][gene] = []
        sorting_groups[seqid][gene].append(entry)

    sortedGff = []
    for seqid in sorting_groups:
        contigData = {}
        for gene in sorting_groups[seqid]:
            contigData[gene] = sortGene(sorting_groups[seqid][gene])
        sortedGff.extend(sortContig(contigData))

    return sortedGff
        

def addExons(gff):

    exon_check = {}
    for i in range(len(gff)):
        entry = gff[i]
        if entry['source'] == 'AUGUSTUS':
            if entry['type'] in {'CDS', 'exon'}:
                prot = re.search(r'ID=(.*?_\d+)', entry['attributes'])[1]
                if prot not in exon_check:
                    exon_check[prot] = [False, i, entry]
                if entry['type'] == 'exon':
                    exon_check[prot] = [True, i, None]
           
    prots = reversed(list(exon_check.keys())) 
    for prot in prots:
        if not exon_check[prot][0]:
            newEntry = copy.deepcopy(exon_check[prot][2])
            newEntry['attributes'] = newEntry['attributes'].replace('.cds;', '.exon1;')
            newEntry['type'] = 'exon'
            gff.insert(exon_check[prot][1] + 1, newEntry)

    return gff

def main(gff_path, prefix, fail = True):

    if isinstance(gff_path, str):
        gff = gff2list(gff_path)
    elif isinstance(gff_path, list):
        gff = gff_path

    exonGtf, comps = intron2exon(gff)
    exonGtfCur = curCDS(exonGtf, re.compile(comps['id']))
    exonGtfCurGenes, failed, flagged = add_genes(exonGtfCur, safe = fail, 
                                                 comps = comps)

    preGff = remove_start_stop(exonGtfCurGenes)
    unsortedGff, trans_str = curate(preGff, prefix, failed, comps) 
#    unsortedGff = addExons(gffUncur)

    if prefix:
        out_gff = sortGFF(unsortedGff, re.compile(gff3Comps()['Alias']))
        gff = rename_and_organize(out_gff)
    else:
        gff = unsortedGff

    return gff, trans_str, failed, flagged


def sortMain(gff, prefix):

    id_comp = re.compile(gff3Comps()['id'])
    crude_sort = sorted(gff, key = lambda x: \
        int(re.search(r'ID=' + prefix + '_(\d+)', x['attributes'])[1]))
    gff = preSortGFF(crude_sort, id_comp)

    return gff


def cli():

    parser = argparse.ArgumentParser(description = 'Curates Funannotate or ' + \
        'post OrthoFiller output by naming accessions as <PREFIX>_####, where ' + \
        '#### is the is the index from ordered accessions. If you are planning' + \
        ' on using OrthoFiller, you may want to wait to not mix accessions.')
    parser.add_argument('-g', '--gff', required = True, help = '.gtf, .gff/.gff3')
    parser.add_argument(
        '-a', '--assembly', 
        help = 'Assembly fasta for proteome export. Requires -p' 
       )
    parser.add_argument('-p', '--prefix', help = 'Ome code for alias')
    parser.add_argument('-s', '--sort', action = 'store_true',
        help = 'Only sort compatible gff3')
    parser.add_argument('-o', '--output', help = 'Output directory')
#    parser.add_argument('--fail', default = True, action = 'store_false',
#        help = 'Fail genes without start or stop codons.')
    args = parser.parse_args()

    if args.output:
        output = format_path(args.output)
        if not os.path.isdir(output):
            os.mkdir(output)
            output += '/'
        output += args.prefix
    else:
        output = args.prefix

    if args.sort:
        gff = sortMain(gff2list(format_path(args.gff)), args.prefix)
        with open(output + '.gff3', 'w') as out:
            out.write(list2gff(gff) + '\n')
        sys.exit(0)
#    elif not args.assembly:
 #       eprint('\nERROR: assembly (`-f`) required for curation\n', flush = True)
  #      sys.exit(2)

    if args.prefix:
        if '_' in args.prefix:
            eprint('\nERROR: "_" not allowed in prefix\n', flush = True)
            sys.exit(1)

    gff, trans_str, failed, flagged = main(
        format_path(args.gff), args.prefix #, args.fail
        )
    if args.assembly and args.prefix:
        fna = fa2dict(format_path(args.assembly))
        faa = gff2proteome(gff, fna)
        with open(output + '.faa', 'w') as out:
            out.write(dict2fa(faa))

    with open(output + '.gff3', 'w') as out:
        out.write(list2gff(gff) + '\n')
    with open(output + '.transitions', 'w') as out:
        out.write(trans_str)
    if failed:
        print('\n' + str(len(failed)) + ' failures\n', flush = True)
        with open(output + '.failed', 'w') as out:
            out.write('\n'.join(['\t'.join(x) for x in failed]))
    if flagged:
        print('\n' + str(len(flagged)) + ' flagged\n', flush = True)
        with open(output + '.flagged', 'w') as out:
            out.write('\n'.join(flagged))

    sys.exit(0)


if __name__ == '__main__':
    cli()
