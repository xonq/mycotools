#! /usr/bin/env python3

# NEED to curate CDS parents to ensure they relate to the RNA
# NEED to curate protein ids to be a uniform title

import os
import re
import sys
import copy
from collections import defaultdict
from mycotools.lib.kontools import format_path, sys_start, eprint
from mycotools.lib.biotools import gff2list, list2gff, gff3Comps

class RNAError(Exception):
    pass

def addMissing(gff_list, intron, comps, ome):
    out_genes, t_list, rnas, introns = {}, [], {}, {}
    mtdb_count, pseudocount, alt_alias = 1, 1, {}
    rna_changes = {} # a dictionary for changing ambigious rna id names for
    # explicit RNA type references
    for entry in gff_list:
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
                    if not entry['attributes'].endswith( ';' ):
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
                    raise RNAError('RNA ID annotated as gene')
                addEntry['attributes'] = addEntry['attributes'].replace('ID=rna-', 'ID=gene-')
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
                    'exon': [], 'texon': [], 'etc': [] 
                    }
#                out_genes[geneID]['gene'].append(addEntry)
                entry['attributes'] += ';Parent=' + geneID
                par = re.search(comps['par'], entry['attributes'])[1]

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
        if geneInfo['rna']:
#            if geneInfo['rna'][0]['type'] != 'mRNA' and not geneInfo['cds']:
 #               del geneInfo['tmrna']
            if len(geneInfo['rna']) > 1:
                multiRNA = True
            if any(x['type'] == 'mRNA' for x in geneInfo['rna']):
                del geneInfo['tmrna']
        elif geneInfo['tmrna'] and not geneInfo['cds']:
            del geneInfo['tmrna']
        if geneInfo['texon'] and geneInfo['exon']:
            if multiRNA:
                exonCoords0 = set(tuple(sorted([int(x['start']), int(x['end'])])) for x in geneInfo['exon'])
                exonCoords1 = set(tuple(sorted([int(x['start']), int(x['end'])])) for x in geneInfo['texon'])
                if all(x in exonCoords0 for x in exonCoords1):
                    del geneInfo['texon']
            else:
                del geneInfo['texon']
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

    prot_comp = re.compile( r'ID=(.*?)[;|$]' )
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

def compileGenes(cur_list, ome, pseudocount = 0, comps = gff3Comps(),
                 cur_seqids = False):

    genes, pseudogenes, rnas = [], [], {}
    mrnas = defaultdict(list) 
    trnas = defaultdict(list)
    rrnas = defaultdict(list)
    ptrnas = defaultdict(list)
    transcripts = defaultdict(list)
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

    id_dict, pseudogenes = {}, set(pseudogenes)
    if cur_seqids:
        for i, entry in enumerate(cur_list):
            if not entry['seqid'].startswith(ome + '_'):
                entry['seqid'] = ome + '_' + entry['seqid']
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
    else:
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

    trna_count, rrna_count, ptrna_count, rna_count, etc_count = 1, 1, 1, 1, 1
    estranged = []
    for i, entry in enumerate(cur_list):
        alias_tag, rna_check = ';Alias=', False
        if 'gene' in entry['type']:
            id_ = re.search(comps['id'], entry['attributes'])[1]
            if id_ in mrnas:
                rna_check = True
                if entry['type'] != 'pseudogene':
                    for rna,alias in id_dict[id_].items():
                        alias_tag += alias + '|'
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
            entry['attributes'] += alias_tag
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
                    try:
                        par = re.search(comps['par'], entry['attributes'])[1]
                    except TypeError: # no parent
                        continue
                    try:
                        alias_tag += id_dict[rnas[par]][par]
                    except KeyError: # no reference gene/reference transcript
#                       eprint(re.search(comps['id'], entry['attributes'])[1], entry['type'])
                        estranged.append((i, id_,))
                        continue
                entry['attributes'] += alias_tag

    for i, id_ in estranged:
        if not 'Alias=' in entry['attributes']:
            gene = rnas[id_]
            try: # acquire from the RNA retroactively
                cur_list[i]['attributes'] += ';Alias=' + id_dict[gene][id_]
            except KeyError: # acquire from the gene
                cur_list[i]['attributes'] += ';Alias=' + id_dict[gene][gene]
                id_dict[gene][id_] = id_dict[gene][gene]
                # this assumes that for the case where CDSs' parents are genes,
                # and mRNAs thus aren't related to the CDS alias directly, then
                # the RNA can be related to the alias through just the gene;
                # however, this will overlook alternately spliced genes in this
                # rare instance. Ideally, CDS parents should be curated to the
                # RNA somehow


    return cur_list


def curGff3(gff_list, ome, cur_seqids = False):

    cur_list, intron = [], False
    for line in gff_list:
        if line['type'] == 'intron':
            intron = True

    cur_list, pseudocount = addMissing(gff_list, intron, gff3Comps(), ome)
    final_list = compileGenes(cur_list, ome, pseudocount,
                              cur_seqids = cur_seqids)

    return final_list


def main(gff_path, ome, cur_seqids = False):

    if isinstance(gff_path, str):
        gff = gff2list( format_path(gff_path) )
    elif isinstance(gff_path, list):
        gff = gff_path
    typ = acquireFormat( gff )

    if not typ:
        eprint('\tERROR: type unknown ' + gff_path, flush = True)
        return None

    new_gff = curGff3(gff, ome, cur_seqids)

    return new_gff


if __name__ == '__main__':
    usage = 'Imports gene coordinates file gff3, ome and curates headers'
    sys_start( sys.argv, usage, 3, files = [sys.argv[1]] )
    cur_gff = main( format_path(sys.argv[1]), sys.argv[2] )
    print( list2gff( cur_gff ) , flush = True)
    sys.exit( 0 )
