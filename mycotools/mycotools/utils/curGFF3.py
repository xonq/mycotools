#! /usr/bin/env python3

from mycotools.lib.kontools import formatPath, sysStart, eprint
from mycotools.lib.biotools import gff2list, list2gff, gff3Comps
import re, sys, os, copy

#class GffError(exception):
 #   pass

def addMissing(gff_list, intron, comps, ome):
    out_genes, t_list, rnas, introns = {}, [], {}, {}
    for entry in gff_list:
        addEntry = None
        id_ = re.search(comps['id'], entry['attributes'])[1]
        if 'gene' in entry['type']:
            out_genes[id_] = {'gene': [entry], 'tmrna': [], 'rna': [], 'cds': [], 'exon': [], 'texon': [], 'etc': [] }
            if 'gene_biotype=protein_coding' in entry['attributes']:
                addEntry = copy.deepcopy(entry)
                addEntry['type'] = 'mRNA'
                addEntry['attributes'] = addEntry['attributes'].replace('ID=gene-', 'ID=mrna-')
                addEntry['attributes'] += ';Parent=' + id_
                addEntry['attributes'] = addEntry['attributes'].replace('gbkey=Gene', 'gbkey=mRNA')
                addEntry['attributes'] = addEntry['attributes'].replace('gene_biotype=protein_coding;','')
                out_genes[id_]['tmrna'].append(addEntry)
        elif entry['type'] == 'CDS':
            try:
                par = re.search(comps['par'], entry['attributes'])[1]
            except TypeError:
                continue
            if par in rnas:
                par = rnas[par]
            search = re.search(comps['prot'], entry['attributes'] )
            if search:
                prot = search.groups()[0]
                if not prot:
                    prot = search.groups()[1]
                if re.search( r'Alias=[^;]$', entry['attributes'] ) is None:
                    if not entry['attributes'].endswith( ';' ):
                        entry['attributes'] += ';'
                    entry['attributes'] += 'Alias=' + ome + '_' + prot
    
            entry['attributes'] = entry['attributes'].replace('Parent=gene-', 'Parent=mrna-')
            out_genes[par]['cds'].append(entry)
            if not intron:
                addEntry = copy.deepcopy(entry)
                addEntry['type'] = 'exon'
                addEntry['attributes'] = addEntry['attributes'].replace('ID=cds-', 'ID=exon-')
                addEntry['attributes'] = addEntry['attributes'].replace('gbkey=CDS', 'gbkey=mRNA')
                out_genes[par]['texon'].append(addEntry)
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
                out_genes[par]['etc'].append(entry)

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
            if geneInfo['rna'][0]['type'] != 'mRNA' and not geneInfo['cds']:
                del geneInfo['tmrna']
            else:
                multiRNA = True
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
            out_list.extend(seqEntry)

    return out_list

def hashGff3IDs( gff_list, comps, tag = 'Alias' ):

    id_hash = {}
    for index, entry in enumerate(gff_list):
        hit_id = re.search( comps['id'], entry['attributes'] )[1]
        par_id = re.search( comps['par'], entry['attributes'] )
        alias = re.search( comps[tag], entry['attributes'] )
        if par_id is not None:
            if par_id[1] not in id_hash:
                id_hash[par_id[1]] = [ [], None, [] ]
            if alias is not None:
                id_hash[par_id[1]][1] = alias[1]
            id_hash[par_id[1]][0].append( index )
            id_hash[par_id[1]][2].append( hit_id )

    return id_hash


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
#    else:
 #       prot_comp = re.compile( r'proteinId (.*?);' )
  #      for line in gff_list:
   #         if line['type'] == 'CDS':
    #            break
     #   if prot_comp.search( line['attributes'] ):
      #      return 'jgi_gff'

    return None


#def indices2gff3( gff_list, indices ): 
 #   return [ gff_list[x] for x in indices ]


def checkGff3Alias( alias_dict, cur_list ):

    out_list = []
    for alias in alias_dict:
        temp_list = [cur_list[x] for x in alias_dict[alias]]
        for line in temp_list:
            if not line['attributes'].endswith(';Alias=' + alias):
                line['attributes'] += ';Alias=' + alias
        out_list.extend( temp_list )

    return out_list


def compileGenes( cur_list, comps = gff3Comps() ):

    genes = [
        (i, re.search(comps['id'], v['attributes'])[1]) \
        for i,v in enumerate(cur_list) \
        if 'gene' in v['type']
        ]
    rnas = {
        re.search(comps['id'], v['attributes'])[1]: \
        re.search(comps['par'], v['attributes'])[1]
        for i,v in enumerate(cur_list) \
        if 'RNA' in v['type']
        }

    id_dict = {}
    for i, entry in enumerate(cur_list):
        if 'gene' in entry['type']:
            id_ = re.search(comps['id'], entry['attributes'])[1]
            id_dict[id_] = {}
        elif entry['type'] == 'CDS':
            try:
                par = re.search(comps['par'], entry['attributes'])[1]
            except TypeError:
                continue
            id_dict[rnas[par]][par] = re.search(comps['Alias'], entry['attributes'])[1]

    etc = set()
    for i, entry in enumerate(cur_list):
        alias_tag = ';Alias='
        if 'gene' in entry['type']:
            id_ = re.search(comps['id'], entry['attributes'])[1]
            for rna,alias in id_dict[id_].items():
                alias_tag += alias + '|'
            alias_tag = alias_tag[:-1]
            entry['attributes'] += alias_tag
        else:
           if 'Alias=' not in entry['attributes']:
               if 'RNA' in entry['type']:
                   id_ = re.search(comps['id'], entry['attributes'])[1]
                   try:
                       alias_tag += id_dict[rnas[id_]][id_]
                   except KeyError: # no CDS
                       count = 1
                       new_tag = entry['type']
                       while new_tag + str(count) in etc:
                           count += 1
                       alias_tag += new_tag + str(count)
                       etc.add(new_tag + str(count))
                       id_dict[rnas[id_]][id_] = new_tag + str(count)
               else:
                   try:
                       par = re.search(comps['par'], entry['attributes'])[1]
                       alias_tag += id_dict[rnas[par]][par]
                   except (KeyError, TypeError) as e:
#                       eprint(re.search(comps['id'], entry['attributes'])[1], entry['type'])
                       continue
               entry['attributes'] += alias_tag

    '''
    gene_dict = {}
    for info in genes:
        gene = info[1]
        gene_dict[gene] = {
            'children': [info[0]],
            'alias': [],
            'tocheck': { gene }
            }
        while gene_dict[gene]['tocheck']:
            hit = list(gene_dict[gene]['tocheck'])[0]
            if hit in par2id:
                gene_dict[gene]['children'].extend(
                    par2id[hit][0]
                    )
                gene_dict[gene]['tocheck'] = gene_dict[gene]['tocheck'] \
                    | set( par2id[hit][2] )
                if par2id[hit][1]:
                    gene_dict[gene]['alias'].append(par2id[hit][1])
            gene_dict[gene]['tocheck'].remove(
                hit
                )

    alias_dict = { 
        gene_dict[x]['alias']: sorted(gene_dict[x]['children']) for x in gene_dict 
        }
    if None in alias_dict:
        del alias_dict[None]
    etc = [ gene_dict[x]['children'] for x in gene_dict if not gene_dict[x]['alias']]
    etc_list = []
    for hit in etc:
        for i in hit:
            etc_list.append( cur_list[i] )

    return alias_dict, etc_list  
    '''
    return cur_list


def curGff3pass1( gff_list, prot_comp, ome ):

    cur_list = []
    for line in gff_list:
        search = prot_comp.search( line['attributes'] )
        if search:
            prot = search[1]
            if re.search( r'Alias=[^;]$', line['attributes'] ) is None:
                if not line['attributes'].endswith( ';' ):
                    line['attributes'] += ';'
                line['attributes'] += 'Alias=' + ome + '_' + prot
        cur_list.append( line )

    return cur_list
                

def curGff3( gff_list, ome ):

    cur_list, intron = [], False
    for line in gff_list:
    # this order is sketchy, really don't know if it's conserved
        if line['type'] == 'intron':
            intron = True

#    if intron:
    cur_list = addMissing(gff_list, intron, gff3Comps(), ome)

 #   if not cur_list:
#        cur_list = curGff3pass1( gff_list, prot_comp, ome )

#    par2id = hashGff3IDs( cur_list, gff3Comps(), 'Alias' )
    final_list = compileGenes(cur_list) #, par2id)
#    alias_dict, etc_list = compileGenes( cur_list, par2id )
 #   final_list = checkGff3Alias( alias_dict, cur_list )
  #  final_list.extend( etc_list )

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
