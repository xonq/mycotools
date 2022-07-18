#! /usr/bin/env python3

import os
import re
import sys
import copy
from collections import defaultdict
from mycotools.lib.kontools import formatPath, sysStart, eprint
from mycotools.lib.biotools import gff2list, list2gff, gff3Comps

#class GffError(exception):
 #   pass

def addMissing(gff_list, intron, comps, ome):
    out_genes, t_list, rnas, introns = {}, [], {}, {}
    mtdb_count, pseudocount, alt_alias = 1, 1, {}
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
        elif 'RNA' in entry['type']:
            try:
                par = re.search(comps['par'], entry['attributes'])[1]
            except TypeError: # no gene in the ids
                addEntry = copy.deepcopy(entry)
                addEntry['type'] = 'gene'
                geneID = id_.replace('rna-', 'gene-')
                if geneID == id_:
                    raise TypeError
                addEntry['attributes'] = addEntry['attributes'].replace('ID=rna-', 'ID=gene-')
                addEntry['attributes'] = addEntry['attributes'].replace(entry['type'], 'gene')
                out_genes[geneID] = {'gene': [entry], 'tmrna': [], 'rna': [], 'cds': [], 'exon': [], 'texon': [], 'etc': [] }
                out_genes[geneID]['gene'].append(addEntry)
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
            except TypeError:
                continue
            if entry['type'] == 'exon':
                try:
                    out_genes[rnas[par]]['exon'].append(entry)
                except KeyError:
                    entry['attributes'] = entry['attributes'].replace('Parent=gene-', 'Parent=mrna-')
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
    return out_list

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

def compileGenes(cur_list, ome, comps = gff3Comps()):

    genes, pseudogenes, rnas = [], [], {}
    mrnas = defaultdict(list) 
    trnas = defaultdict(list)
    rrnas = defaultdict(list)
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


    id_dict, pseudogenes = {}, set(pseudogenes)
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
            except KeyError:
                if par in pseudogenes:
                    id_dict[par][par] = re.search(comps['Alias'],
                                                  entry['attributes'])[1]

    etc, trna_count, rrna_count = set(), 1, 1
    for i, entry in enumerate(cur_list):
        alias_tag = ';Alias='
        if 'gene' in entry['type']:
            id_ = re.search(comps['id'], entry['attributes'])[1]
            if id_ in mrnas:
                for rna,alias in id_dict[id_].items():
                    alias_tag += alias + '|'
            if id_ in trnas:
                for rna in trnas[id_]: 
                    alias = ome + '_tRNA' + str(trna_count)
                    id_dict[id_][rna] = alias 
                    # only works if the gff list is sequential
                    alias_tag += alias + '|'
                    trna_count += 1
            if id_ in rrnas:
                for rna in rrnas[id_]:
                    alias = ome + '_rRNA' + str(rrna_count)
                    id_dict[id_][rna] = alias
                    # only works if gff list is sequential
                    alias_tag += alias + '|'
                    rrna_count += 1
            alias_tag = alias_tag[:-1]
            entry['attributes'] += alias_tag

                   
        else:
           if 'Alias=' not in entry['attributes']:
               if 'RNA' in entry['type']:
                   id_ = re.search(comps['id'], entry['attributes'])[1]
        #           try:
                   alias_tag += id_dict[rnas[id_]][id_]
#                   except KeyError: # no CDS
 #                      count = 1
  #                     new_tag = entry['type']
   #                    while new_tag + str(count) in etc:
    #                       count += 1
     #                  alias_tag += new_tag + str(count)
      #                 etc.add(new_tag + str(count))
       #                id_dict[rnas[id_]][id_] = new_tag + str(count)
               else:
                   try:
                       par = re.search(comps['par'], entry['attributes'])[1]
                   except TypeError: # no parent
                       continue
                   try:
                       alias_tag += id_dict[rnas[par]][par]
                   except KeyError: # no reference gene
#                       eprint(re.search(comps['id'], entry['attributes'])[1], entry['type'])
                       continue
               entry['attributes'] += alias_tag

    return cur_list


def curGff3( gff_list, ome ):

    cur_list, intron = [], False
    for line in gff_list:
        if line['type'] == 'intron':
            intron = True

    cur_list = addMissing(gff_list, intron, gff3Comps(), ome)
    final_list = compileGenes(cur_list, ome) #, par2id)

    return final_list


def main( gff_path, ome):

    if isinstance(gff_path, str):
        gff = gff2list( formatPath(gff_path) )
    elif isinstance(gff_path, list):
        gff = gff_path
    typ = acquireFormat( gff )

    if not typ:
        eprint('\tERROR: type unknown ' + gff_path, flush = True)
        return None

    new_gff = curGff3( gff, ome )

    return new_gff


if __name__ == '__main__':
    usage = 'Imports gene coordinates file gff3, ome, and curates headers'
    sysStart( sys.argv, usage, 3, files = [sys.argv[1]] )
    cur_gff = main( formatPath(sys.argv[1]), sys.argv[2] )
    print( list2gff( cur_gff ) , flush = True)
    sys.exit( 0 )
