#! /usr/bin/env python3

# NEED to remove introns
# NEED to make compatible with rerunning
# NEED to curate CDS parents to ensure they relate to the RNA
# NEED to curate protein ids to be a uniform title
# NEED to fix mobile_genetic_element conversion from NCBI
    # entries contain a gene, mRNA, exon
# NEED to fix standalone gene alias
    # an alias is applied like "Alias" with nothing else
    # only sometimes
        # if followed by unreviewed type
# NEED validation module to check common errors
    # lacking alias
    # overlapping genes
    # RNAs with no parent
    # CDS/exons with no parent

import os
import re
import sys
import copy
from collections import defaultdict
from mycotools.lib.kontools import format_path, sys_start, eprint
from mycotools.lib.biotools import gff2list, list2gff, gff3Comps

class RNAError(Exception):
    pass

class GeneError(Exception):
    pass


def add_missing(gff_list, intron, comps, ome):
    accepted_types = {'pseudogene', 'pseudogenic_trna', 'trna',
                      'rrna', 'mrna', 'rna', 'transcript', 'gene',
                      'cds', 'exon', 'cds', 'three_prime_utr', 'intron',
                      'five_prime_utr', '5_prime_utr', '3_prime_utr'}
    out_genes, t_list, rnas, introns = {}, [], {}, {}
    mtdb_count, pseudocount, alt_alias = 1, 1, {}
    rna_changes = {} # a dictionary for changing ambigious rna id names for
    # explicit RNA type references
    for entry in gff_list:
        if entry['type'].lower() not in accepted_types:
            if 'RNA' in entry['type']: # converts snRNAs
                entry['type'] = 'RNA'
            else:
                continue
        entry['attributes'] = \
            entry['attributes'].replace('proteinId', 'protein_id')
        entry['attributes'] = \
            entry['attributes'].replace('transcriptId', 'transcript_id')
        # remove old alias if present
        entry['attributes'] = re.sub(r';?Alias=[^;]+', '', entry['attributes'])
        addEntry = None
        id_ = re.search(comps['id'], entry['attributes'])[1]
        if 'gene' in entry['type']: # includes pseudogenes
            out_genes[id_] = {
                'gene': [entry], 'tmrna': [], 'rna': [], 
                'cds': [], 'exon': [], 'texon': [], 
                'etc': [], 'pseudo': entry['type'] == 'pseudogene'
                } # create an entry for the genes
            if 'gene_biotype=protein_coding' in entry['attributes']: # if it is
            # NCBI protein coding gene (or conforming format)
                addEntry = copy.deepcopy(entry)
                addEntry['type'] = 'mRNA'
                addEntry['attributes'] = addEntry['attributes'].replace('ID=gene-', 'ID=mrna-')
                addEntry['attributes'] += ';Parent=' + id_
                addEntry['attributes'] = addEntry['attributes'].replace('gbkey=Gene', 'gbkey=mRNA')
                addEntry['attributes'] = addEntry['attributes'].replace('gene_biotype=protein_coding;','')
                out_genes[id_]['tmrna'].append(addEntry) # add a temporary RNA
                # in case there isn't one
        elif entry['type'] == 'CDS':
            try: # attempt to acquire the parent
                par = re.search(comps['par'], entry['attributes'])[1]
            except TypeError:
                continue # if there isn't a parent then skip it
            if par in rna_changes:
                entry['attributes'] = re.sub(comps['par'], 'Parent=' \
                                             + rna_changes[par],
                                             entry['attributes'])
                par = rna_changes[par]
            if par in rnas: # if the parent is an RNA, acquire the gene
                rna_par = par
                par = rnas[par] # change par to gene
            else:
                rna_par = None
            search = re.search(comps['prot'], entry['attributes'])
            if search: # search for protein
                prot = search.groups()[0] # try the first search
                if not prot:
                    prot = search.groups()[1] # try the second
                if re.search(comps['Alias'], entry['attributes'] ) is None:
                # is there a mycotools alias or associated alias?
                    if not entry['attributes'].endswith(';'):
                        entry['attributes'] += ';'
                    alias = ome + '_' + prot # create an alias with the found
                    # protein
                    entry['attributes'] += 'Alias=' + ome + '_' + prot
                    alt_alias[rna_par] = ome + '_' + prot
            else: # doesnt have a protein id
                if rna_par: # if there is an RNA
                    if rna_par not in alt_alias: # add it in
                        alias = ome + '_mtdb' + str(mtdb_count)
                        mtdb_count += 1
                        alt_alias[rna_par] = alias
                    else: # acquire the alias that is used
                        alias = alt_alias[rna_par]
                else: # no RNA parent
                    if par not in alt_alias:
                        if out_genes[par]['pseudo']: # if it is a pseudogene
                            alias = ome + '_pseudogene' + str(pseudocount)
                            pseudocount += 1
                        else:
                            alias = ome + '_mtdb' + str(mtdb_count)
                            mtdb_count += 1
                        alt_alias[par] = alias
                    else:
                        alias = alt_alias[par]
                if not entry['attributes'].endswith(';'):
                    entry['attributes'] += ';'
                entry['attributes'] += 'Alias=' + alias

            if not out_genes[par]['pseudo']: # if it isnt a defined pesudogene
                entry['attributes'] = entry['attributes'].replace('Parent=gene-', 'Parent=mrna-')
                out_genes[par]['cds'].append(entry)
                if not intron:
                    addEntry = copy.deepcopy(entry)
                    addEntry['type'] = 'exon'
                    addEntry['attributes'] = addEntry['attributes'].replace('ID=cds-', 'ID=exon-')
                    addEntry['attributes'] = addEntry['attributes'].replace('gbkey=CDS', 'gbkey=mRNA')
                    out_genes[par]['texon'].append(addEntry)
            else: # else just add the unmodified CDS
                out_genes[par]['cds'].append(entry)
        elif 'RNA' in entry['type'] or entry['type'] == 'transcript':
            entry['type'] = entry['type'].replace('transcript','RNA')
            if id_.startswith('rna'): # make RNA ID explicit
                new_id = entry['type'].lower() + id_[3:]
                entry['attributes'] = re.sub(comps['id'], 'ID=' + new_id,
                                             entry['attributes'])
                rna_changes[id_] = new_id
                id_ = new_id
            try:
                par = re.search(comps['par'], entry['attributes'])[1]
            except TypeError: # no gene in the ids, initialize new entry
                addEntry = copy.deepcopy(entry)
                addEntry['type'] = 'gene'
                geneID = id_.replace('rna-', 'gene-')
                if geneID == id_:
                    if id_.startswith('trna'):
                        geneID = 'gene-' + id_
                    else:
                        raise RNAError(f'RNA ID annotated as gene: {entry}')
                addEntry['attributes'] = re.sub(r'ID=' + id_, 'ID=' + geneID, addEntry['attributes'])
                #addEntry['attributes'].replace('ID=rna-', 'ID=gene-')
                addEntry['attributes'] = addEntry['attributes'].replace(entry['type'], 'gene')
                addEntry['attributes'] = re.sub(r'gbkey=[^;]+', 'gbkey=gene',
                                                addEntry['attributes'])
                addEntry['attributes'] = re.sub(comps['par'] + ';', '',
                                                addEntry['attributes'])
                                                # remove parent
                entry['attributes'] = re.sub(comps['par'] + ';', '',
                                             entry['attributes'])
                out_genes[geneID] = {
                    'gene': [addEntry], 'tmrna': [], 'rna': [], 'cds': [], 
                    'exon': [], 'texon': [], 'etc': [], 'pseudo': False 
                    }
#                out_genes[geneID]['gene'].append(addEntry)
                entry['attributes'] += ';Parent=' + geneID
                par = re.search(comps['par'], entry['attributes'])[1]

            search = re.search(comps['prot'], entry['attributes'])
            if search: # search for protein
                prot = search.groups()[0] # try the first search
                if not prot:
                    prot = search.groups()[1] # try the second
                if re.search(comps['Alias'], entry['attributes'] ) is None:
                # is there a mycotools alias or associated alias?
                    if not entry['attributes'].endswith(';'):
                        entry['attributes'] += ';'
                    alias = ome + '_' + prot # create an alias with the found
                    # protein
                    entry['attributes'] += 'Alias=' + ome + '_' + prot
                    alt_alias[id_] = ome + '_' + prot


            rnas[id_] = par
            out_genes[par]['rna'].append(entry)
        elif entry['type'] == 'intron':
            try:
                par = re.search(comps['par'], entry['attributes'])[1]
            except TypeError:
                continue
            if par not in introns:
                introns[par] = []
            introns[par].append(entry)
        else:
            try:
                par = re.search(comps['par'], entry['attributes'])[1]
            except TypeError: # skip no parents, should be a flag
                continue
            if entry['type'] == 'exon':
                if par in rna_changes:
                    entry['attributes'] = re.sub(comps['par'], 'Parent=' \
                                                 + rna_changes[par],
                                                 entry['attributes'])
                    par = rna_changes[par]
                try:
                    out_genes[rnas[par]]['exon'].append(entry)
                except KeyError:
                    entry['attributes'] = entry['attributes'].replace('Parent=gene-', 'Parent=mrna-')
                    # assumes it is derived from an mRNA
                    out_genes[par]['exon'].append(entry)
            else:
                if entry['type'].lower() in {'3_prime_utr', 'three_prime_utr'}:
                    entry['type'] = 'three_prime_UTR'
                elif entry['type'].lower() in {'5_prime_utr', 'five_prime_utr'}:
                    entry['type'] = 'five_prime_UTR'
                try:
                    out_genes[par]['etc'].append(entry)
                except KeyError: # derived from an RNA
                    out_genes[rnas[par]]['etc'].append(entry)

    if introns:
        exons = {}
        for gene, intronList in introns.items():
            if out_genes[gene]['exon']:
                continue
            geneCoords = sorted([int(out_genes[gene]['gene'][0]['start']), int(out_genes[gene]['gene'][0]['end'])])
            intronList = sorted(intronList, key = lambda x: int(x['start']))
            intronCoords = [sorted([int(x['start']), int(x['end'])]) for x in intronList]
            if len(intronCoords) == 1:
                exonCoords0 = [geneCoords[0], intronCoords[0][1] - 1]
                exonCoords1 = [intronCoords[0][1] + 1, geneCoords[1]]
                addEntry = copy.deepcopy(intronList[0])
                addEntry['attributes'] = addEntry['attributes'].replace('gbkey=intron', 'gbkey=mRNA')
                addEntry['attributes'] = addEntry['attributes'].replace('ID=intron-', 'ID=exon-')
                addEntry['type'] = 'exon'
                addEntry['start'], addEntry['end'] = str(exonCoords0[0]), str(exonCoords0[1])
                out_genes[gene]['exon'].append(addEntry)
                addEntry1 = copy.deepcopy(addEntry)
                addEntry1['start'], addEntry1['end'] = str(exonCoords1[0]), str(exonCoords1[1])
                out_genes[gene]['exon'].append(addEntry1)
            else:
                for i, entryCoords in enumerate(intronCoords):
                    if i == 0:
                        exonCoords = [geneCoords[0], entryCoords[0] - 1]
                    elif i == len(intronCoords) - 1:
                        exonCoords = [entryCoords[1] + 1, geneCoords[1]]
                    else:
                        exonCoords = [intronCoords[i-1][1]+1, intronCoords[i][0] - 1]
                    addEntry = copy.deepcopy(intronList[i])
                    addEntry['attributes'] = addEntry['attributes'].replace('gbkey=intron', 'gbkey=mRNA')
                    addEntry['attributes'] = addEntry['attributes'].replace('ID=intron-', 'ID=exon-')
                    addEntry['type'] = 'exon'
                    addEntry['start'], addEntry['end'] = str(exonCoords[0]), str(exonCoords[1])
                    out_genes[gene]['exon'].append(addEntry)

    out_list = []
    for geneID, geneInfo in out_genes.items():
        multiRNA = False
        if not any(x['type'] in {'RNA', 'mRNA'} for x in geneInfo['rna']) \
            and not geneInfo['tmrna'] and geneInfo['cds']:
            # assume there is an mrna involved - I don't like having to do this
            # but I don't like having to do most all of this because  the files
            # are so inconsistently formatted
            geneInfo['tmrna'] = copy.deepcopy(geneInfo['gene'])
            id_ = re.search(gff3Comps()['id'], geneInfo['tmrna'][0]['attributes'])[1]
            new_id = 'mrna' + id_[4:]
            geneInfo['tmrna'][0]['attributes'] = f'ID={new_id};Parent={geneID};gbkey=mRNA'
            geneInfo['tmrna'][0]['type'] = 'mRNA'
            for cds in geneInfo['cds']:
                cds['attributes'] = re.sub(gff3Comps()['par'], 'Parent=' + new_id, cds['attributes'])

            
        if geneInfo['rna']:
#            if geneInfo['rna'][0]['type'] != 'mRNA' and not geneInfo['cds']:
 #               del geneInfo['tmrna']
            if len(geneInfo['rna']) > 1:
                multiRNA = True
            if any(x['type'] == 'mRNA' for x in geneInfo['rna']):
                del geneInfo['tmrna']
        elif geneInfo['tmrna'] and not geneInfo['cds']:
            del geneInfo['tmrna']
        elif geneInfo['tmrna']:
            # check for difficult-to-detect alternate splicing
            prots = []
            for cds in geneInfo['cds']:
                search = re.search(comps['prot'], cds['attributes'])
                if search:
                    prot = search.groups()[0] # try the first search
                    if not prot:
                        prot = search.groups()[1] # try the second
                    prots.append(prot)
            prots = set(prots)
            if len(prots) > 1: # alternate splicing detected
                prot_coords, prot_dir, prot2ali, prot2i = {}, {}, {}, {}
                count = 0
                for cds in geneInfo['cds']:
                    search = re.search(comps['prot'], cds['attributes'])
                    if search:
                        prot = search.groups()[0] # try the first search
                        if not prot:
                            prot = search.groups()[1] # try the second
                        alias = re.search(comps['Alias'], cds['attributes'])[1]
                        prot_dir[prot] = cds['strand']
                        prot2ali[prot] = alias
                        if prot not in prot_coords:
                            prot_coords[prot] = []
                            prot2i[prot] = count
                            count += 1
                        parent = re.search(comps['par'], cds['attributes'])[1].replace('gene-', 'mrna=') + \
                                 '-T' + str(prot2i[prot])
                        cds['attributes'] = re.sub(comps['par'], 'Parent=' + parent, cds['attributes'])
                        prot_coords[prot].extend([cds['start'], cds['end']])

                # reinitialize the trmnas
                geneInfo['tmrna'] = []
                for prot, coords in prot_coords.items():
                    start, stop = min(coords), max(coords)
                    geneInfo['tmrna'].append(copy.deepcopy(geneInfo['gene'][0]))
                    geneInfo['tmrna'][-1]['start'] = start
                    geneInfo['tmrna'][-1]['end'] = stop
                    geneInfo['tmrna'][-1]['strand'] = prot_dir[prot]
                    geneInfo['tmrna'][-1]['type'] = 'mRNA'

                    parent = re.search(comps['id'], geneInfo['gene'][0]['attributes'])[1]
                    new_id = 'mrna-' + parent[4:] + '-T' + str(prot2i[prot])
                    new_id = new_id.replace('--','-')
                    geneInfo['tmrna'][-1]['attributes'] = re.sub(comps['id'], 'ID=' + new_id + ';' \
                                                            + 'Parent=' + parent,
                                                 geneInfo['tmrna'][-1]['attributes'])
                    geneInfo['tmrna'][-1]['attributes'] += ';Alias=' + prot2ali[prot]
                    count += 1

                geneInfo['exon'] = copy.deepcopy(geneInfo['cds'])
                for exon in geneInfo['exon']:
                    id_ = re.search(comps['id'], exon['attributes'])[1]
                    par = re.search(comps['par'], exon['attributes'])[1]
                    ali = re.search(comps['Alias'], exon['attributes'])[1]
                    new_id = 'exon' + id_[3:]
                    exon['attributes'] = f'ID={new_id};Parent={par};Alias={ali}'
                    exon['type'] = 'exon'

                       
        if geneInfo['exon']:
            del geneInfo['texon']
#        if geneInfo['texon'] and geneInfo['exon']:
    #        if multiRNA:
#                exonCoords0 = set(tuple(sorted([int(x['start']), int(x['end'])])) for x in geneInfo['exon'])
 #               exonCoords1 = set(tuple(sorted([int(x['start']), int(x['end'])])) for x in geneInfo['texon'])
  #              if all(x in exonCoords0 for x in exonCoords1):
   #                 del geneInfo['texon']
#            else:
 #               del geneInfo['texon']
        if geneInfo['pseudo']:
            for cds in geneInfo['cds']:
                cds['attributes'] = re.sub(gff3Comps()['Alias'], '', cds['attributes'])

        for seqType, seqEntry in geneInfo.items():
            if seqType != 'pseudo':
                out_list.extend(seqEntry)

    for entry in out_list: # this will relate modified RNA IDs from
    # nonsequential gffs post hoc. not sure if there are more downstream
    # implications to changing RNAs after the above curation
        try:
            par = re.search(comps['par'], entry['attributes'])[1]
        except TypeError:
            continue
        if par in rna_changes:
            new_par = rna_changes[par]
            entry['attributes'] = re.sub(comps['par'], 'Parent=' + new_par,
                                         entry['attributes'])
    return out_list, pseudocount

def acquireFormat( gff_list ):

    prot_comp = re.compile(gff3Comps()['id'])
    gene = False
    for line in gff_list:
        if line['type'] == 'gene':
            gene = True
            break

    if gene:
        if prot_comp.search( line['attributes'] ):
            if 'Name=jgi.p' in line['attributes']:
                return 'jgi_gff3'
            else:
                return 'misc_gff3'

    return None

def compile_genes(cur_list, ome, pseudocount = 0, comps = gff3Comps(),
                 cur_seqids = False):

    genes, pseudogenes, rnas = [], [], {}
    mrnas = defaultdict(list) 
    trnas = defaultdict(list)
    rrnas = defaultdict(list)
    ptrnas = defaultdict(list)
    transcripts = defaultdict(list)
    rna_etcs = defaultdict(list)
    gen_etcs = defaultdict(list)
    utrs = False
    for i, v in enumerate(cur_list):
        id_ = re.search(comps['id'], v['attributes'])[1]
        if 'gene' in v['type']:
            genes.append((i, id_,))
            if v['type'] == 'pseudogene':
                pseudogenes.append(id_)
        elif 'RNA' in v['type']:
            par = re.search(comps['par'], v['attributes'])[1]
            rnas[id_] = par # rnaID2geneID
            if v['type'] == 'mRNA':
                mrnas[par].append(id_)
            elif v['type'] == 'tRNA':
                trnas[par].append(id_)
            elif v['type'] == 'rRNA':
                rrnas[par].append(id_)
            elif v['type'] == 'pseudogenic_tRNA':
                ptrnas[par].append(id_)
            elif v['type'] == 'RNA':
                transcripts[par].append(id_)
        elif 'UTR' in v['type']:
            utrs = True

    if utrs:
        for i, v in enumerate(cur_list):
            if 'UTR' in v['type']:
                id_ = re.search(comps['id'], v['attributes'])[1]
                par = re.search(comps['par'], v['attributes'])[1]
                if par in rnas:
                    rna_etcs[id_] = par
                else:
                    gen_etcs[id_] = par
    

    todel = []
    for par, trna in trnas.items():
        if par in pseudogenes:
            ptrnas[par].extend(trna)
            todel.append(par)
    for i in todel:
        del trnas[i]

    id_dict, pseudogenes, estranged_cds = {}, set(pseudogenes), []
    for i, entry in enumerate(cur_list):
        if 'gene' in entry['type']:
            id_ = re.search(comps['id'], entry['attributes'])[1]
            id_dict[id_] = {}
        elif entry['type'] == 'CDS':
            try:
                par = re.search(comps['par'], entry['attributes'])[1]
            except TypeError:
                continue
            try:
                id_dict[rnas[par]][par] = re.search(comps['Alias'], entry['attributes'])[1]
            except KeyError: # pseudogene or CDS parent is related to gene, not
            # RNA
                id_dict[par][par] = re.search(comps['Alias'],
                                              entry['attributes'])[1]
            except TypeError: # missing alias, might be stripped from before
                estranged_cds.append((id_, par))

    trna_count, rrna_count, ptrna_count, rna_count, etc_count = 1, 1, 1, 1, 1
    gene2alias = {}
    estranged, gene_estranged = [], [] 
    for i, entry in enumerate(cur_list):
        alias_tag, rna_check = 'Alias=', False
        if 'gene' in entry['type']:
            id_ = re.search(comps['id'], entry['attributes'])[1]
            if id_ in mrnas:
                rna_check = True
                if entry['type'] != 'pseudogene':
                    if id_dict[id_]:
                        for rna,alias in id_dict[id_].items():
                            alias_tag += alias + '|'
                    else: # no CDS associated with mRNA
                        alias = ome + '_etc' + str(etc_count)
                        alias_tag += alias + '|'
                        etc_count += 1
                        for rna in mrnas[id_]:
                            id_dict[id_][rna] = alias
                            
                else:
                    for rna in mrnas[id_]:
                        alias = ome + '_pseudogene' + str(pseudocount)
                        id_dict[id_][rna] = alias
                        alias_tag += alias + '|'
                        pseudocount += 1
            if id_ in trnas:
                rna_check = True
                for rna in trnas[id_]: 
                    alias = ome + '_tRNA' + str(trna_count)
                    id_dict[id_][rna] = alias 
                    # only works if the gff list is sequential
                    alias_tag += alias + '|'
                    trna_count += 1
            if id_ in rrnas:
                rna_check = True
                for rna in rrnas[id_]:
                    alias = ome + '_rRNA' + str(rrna_count)
                    id_dict[id_][rna] = alias
                    # only works if gff list is sequential
                    alias_tag += alias + '|'
                    rrna_count += 1
            if id_ in ptrnas:
                rna_check = True
                for rna in ptrnas[id_]:
                    alias = ome + '_pseudo-tRNA' + str(ptrna_count)
                    id_dict[id_][rna] = alias
                    alias_tag += alias + '|'
                    ptrna_count += 1
            if id_ in transcripts:
                rna_check = True
                for rna in transcripts[id_]:
                    alias = ome + '_RNA' + str(rna_count)
                    id_dict[id_][rna] = alias
                    alias_tag += alias + '|'
                    rna_count += 1
            if not rna_check and id_ not in pseudogenes: # no association with RNA
                alias = ome + '_etc' + str(etc_count)
                alias_tag += alias + '|'
                etc_count += 1
            alias_tag = alias_tag[:-1]
            if entry['attributes'].rstrip().endswith(';'):
                entry['attributes'] += alias_tag
            else:
                entry['attributes'] += ';' + alias_tag
            gene2alias[id_] = alias_tag
        else:
            if 'Alias=' not in entry['attributes']:
                if 'RNA' in entry['type']:
                    id_ = re.search(comps['id'], entry['attributes'])[1]
                    try:
                        alias_tag += id_dict[rnas[id_]][id_]
                    except KeyError:
                        if id_ in rnas:
                            estranged.append((i, id_,))
                            continue
#                            raise KeyError
                        else:
                            raise RNAError('unknown RNA type: ' + entry['type'])
                else:
                    id_ = re.search(comps['id'], entry['attributes'])[1]
                    try:
                        par = re.search(comps['par'], entry['attributes'])[1]
                    except TypeError: # no parent
                        continue
                    if id_ in gen_etcs:
                        if par in gene2alias:
                            alias_tag = gene2alias[par]
                        elif par in rnas:
                            alias_tag = id_dict[rnas[id_]][par]
                        else:
                            gene_estranged.append((i, par,))
                    else:
                        try:
                            alias_tag += id_dict[rnas[par]][par]
                        except KeyError: # no reference gene/reference transcript
#                           eprint(re.search(comps['id'], entry['attributes'])[1], entry['type'])
                            estranged.append((i, id_,))
                            continue
                if entry['attributes'].rstrip().endswith(';'):
                    entry['attributes'] += alias_tag
                else:
                    entry['attributes'] += ';' + alias_tag

    for i, id_ in estranged:
        if not 'Alias=' in entry['attributes']:
            gene = rnas[id_]
            try: # acquire from the RNA retroactively
                if cur_list[i]['attributes'].rstrip().endswith(';'):
                    cur_list[i]['attributes'] += 'Alias=' + id_dict[gene][id_]
                else:
                    cur_list[i]['attributes'] += ';Alias=' + id_dict[gene][id_]
            except KeyError: # acquire from the gene
                if cur_list[i]['attributes'].rstrip().endswith(';'):
                    cur_list[i]['attributes'] += 'Alias=' + id_dict[gene][gene]
                else:
                    cur_list[i]['attributes'] += ';Alias=' \
                                              +  id_dict[gene][gene]
                id_dict[gene][id_] = id_dict[gene][gene]
                # this assumes that for the case where CDSs' parents are genes,
                # and mRNAs thus aren't related to the CDS alias directly, then
                # the RNA can be related to the alias through just the gene;
                # however, this will overlook alternately spliced genes in this
                # rare instance. Ideally, CDS parents should be curated to the
                # RNA somehow
    for i, id_ in gene_estranged:
        if not 'Alias=' in entry['attributes']:
            if id_ in gene2alias:
                if cur_list[i]['attributes'].rstrip().endswith(';'):
                    cur_list[i]['attributes'] += 'Alias=' + gene2alias[id_]
                else:
                    cur_list[i]['attributes'] += ';Alias=' + gene2alias[id_]

    return cur_list


def rename_and_organize(gff_list):
    alias2geneid, geneid2alias = {}, {}
    scaf2gene2entries = defaultdict(lambda: defaultdict(lambda: {'gene': None, 'rna': {}}))
    todel = []
    for i, entry in enumerate(gff_list):
        if entry['type'] in {'pseudogene', 'gene'}:
            seqid = entry['seqid']
            alias = re.search(gff3Comps()['Alias'], entry['attributes'])[1]
            alias_list = alias.split('|')
            if len(alias_list) > 1:
                gene_id = 'gene_' + alias_list[0] + '-' \
                        + alias_list[-1][alias_list[-1].find('_') + 1:]
            elif len(alias_list) == 1:
                gene_id = 'gene_' + alias_list[0]
            else:
                raise GeneError('gene with no alias')
            for a in alias_list:
                alias2geneid[a] = gene_id
                geneid2alias[gene_id] = a
            entry['attributes'] = re.sub(gff3Comps()['id'], 'ID=' + gene_id, 
                                         entry['attributes'])
            scaf2gene2entries[seqid][gene_id]['gene'] = [entry]
            todel.append(i)

    for i in reversed(todel):
        del gff_list[i]

    todel = []
    alias2rnaid, rnaid2alias = {}, {}
    for i, entry in enumerate(gff_list):
        if 'RNA' in entry['type']:
            alias = re.search(gff3Comps()['Alias'], entry['attributes'])[1]
            rna_id = entry['type'].lower() + '_' + alias
            entry['attributes'] = re.sub(gff3Comps()['id'], 'ID=' + rna_id,
                                         entry['attributes'])
            try:
                gene_id = alias2geneid[alias]
            except KeyError:
                raise KeyError('RNA alias that does not link with gene')
            entry['attributes'] = re.sub(gff3Comps()['par'], 'Parent=' + gene_id,
                                     entry['attributes'])
            seqid = entry['seqid']
            scaf2gene2entries[seqid][gene_id]['rna'][alias] = [entry]
            alias2rnaid[alias] = rna_id
            rnaid2alias[rna_id] = alias
            todel.append(i)

    for i in reversed(todel):
        del gff_list[i]

    id_dict = defaultdict(dict)
    for entry in gff_list:
        try:
            alias = re.search(gff3Comps()['Alias'], entry['attributes'])[1]
        except TypeError:
            raise TypeError(entry)
        typ = entry['type'].lower()
        if typ not in id_dict[alias]:
            id_dict[alias][typ] = 0
        else:
            id_dict[alias][typ] += 1
        oth_id = typ + str(id_dict[alias][typ]) + '_' + alias
        entry['attributes'] = re.sub(gff3Comps()['id'], 'ID=' + oth_id,
                                     entry['attributes'])
        if alias in alias2rnaid:
            par_id = alias2rnaid[alias]
            gene_id = alias2geneid[alias]
            entry['attributes'] = re.sub(gff3Comps()['par'], 'Parent=' + par_id,
                                         entry['attributes'])
            scaf2gene2entries[entry['seqid']][gene_id]['rna'][alias].append(entry)
        elif alias in alias2geneid:
            par_id = alias2geneid[alias]
            entry['attributes'] = re.sub(gff3Comps()['par'], 'Parent=' + par_id,
                                         entry['attributes'])
            scaf2gene2entries[entry['seqid']][par_id]['gene'][alias].append(entry)
        else:
            par_id = re.search(gff3Comps()['par'], entry['attributes'])[1]
            if par_id in rnaid2alias:
                new_alias = rnaid2alias[par_id]
                gene_id = alias2geneid[new_alias]
            elif par_id in geneid2alias:
                new_alias = geneid2alias[par_id]
                gene_id = par_id
            else:
                raise KeyError(f'missing parent: {entry}')
            entry['attributes'] = re.sub(gff3Comps()['par'], 'Parent=' + par_id,
                                         entry['attributes'])
            entry['attributes'] = re.sub(gff3Comps()['Alias'], 'Alias=' + new_alias,
                                         entry['attributes'])
            scaf2gene2entries[entry['seqid']][gene_id]['gene'][new_alias].append(entry)

    out_gff = []
    for seqid in sorted(scaf2gene2entries.keys()):
        gene_data = scaf2gene2entries[seqid]
        sorted_genes = {k: v for k, v in sorted(gene_data.items(), 
                                                key = lambda x: x[1]['gene'][0]['start'])}
        for gene_id, entry_dict in sorted_genes.items():
            out_gff.extend(entry_dict['gene'])
            sorted_rnas = {k: v for k, v in sorted(entry_dict['rna'].items(),
                                                   key = lambda x: x[1][0]['start'])}
            for rna, entries in sorted_rnas.items():
                out_gff.append(entries[0])
                out_gff.extend(sorted(entries[1:], key = lambda x: x['start']))

    return out_gff
        
  


def curGff3(gff_list, ome, cur_seqids = False):

    cur_list, intron = [], False
    for line in gff_list:
        if line['type'] == 'intron':
            intron = True

    cur_list, pseudocount = add_missing(gff_list, intron, gff3Comps(), ome)
    final_list = compile_genes(cur_list, ome, pseudocount,
                              cur_seqids = cur_seqids)

    return final_list


def main(gff_path, ome, cur_seqids = False):

    if isinstance(gff_path, str):
        gff = gff2list( format_path(gff_path) )
    elif isinstance(gff_path, list):
        gff = gff_path
    typ = acquireFormat( gff )

    if not typ:
        eprint('\tERROR: type unknown ', flush = True)
        return None

    new_gff = curGff3(gff, ome, cur_seqids)
    clean_gff = rename_and_organize(new_gff)

    return clean_gff


def cli():
    usage = 'Imports gene coordinates file gff3, ome and curates headers'
    sys_start( sys.argv, usage, 3, files = [sys.argv[1]] )
    cur_gff = main( format_path(sys.argv[1]), sys.argv[2] )
    print( list2gff( cur_gff ) , flush = True)
    sys.exit( 0 )


if __name__ == '__main__':
    cli()
