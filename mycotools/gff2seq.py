#! /usr/bin/env python3

import re
import sys
import argparse
from Bio.Seq import Seq
from mycotools.lib.dbtools import mtdb, primaryDB
from mycotools.lib.biotools import fa2dict, gff2list, gff3Comps, dict2fa
from mycotools.lib.kontools import format_path, sys_start, eprint, stdin2str

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
 
def sortMain(gff):

    id_comp = re.compile(gff3Comps()['id'])
#    crude_sort = sorted( gff, key = lambda x: \
 #       int( re.search(r'ID=([^;]+)', x['attributes'])[1] ))
    gff = sortGFF(gff, id_comp)

    return gff

def grabCDS(gff_dicts, spacer = '\t'):
    """Grab CDSs that are associated with genes. gff_dicts is a
    mycotools.lib.biotools gff2list() list"""

    mrnas, genes, warning, ome = [], [], False, None
    for entry in gff_dicts:
        if entry['type'] == 'mRNA':
            try:
                alias = re.search(gff3Comps()['Alias'], entry['attributes'])[1]
                if not ome:
                    ome = alias[:alias.find('_')]
            except TypeError:
                raise TypeError(str(entry))
            mrnas.extend(alias.split('|')) # account for posttranslational mods
        elif 'gene' in entry['type']:
            try:
                alias = re.search(gff3Comps()['Alias'], entry['attributes'])[1]
                genes.extend(alias.split('|'))
            except TypeError:
                if not warning and ome:
                    eprint(spacer + 'WARNING: ' + str(ome) + ' missing aliases', flush = True)
                    warning = True
                raise TypeError(str(entry))

    mrna_set = set([x for x in mrnas if x in set(genes)])
    out_cds = []
    for entry in gff_dicts:
        if entry['type'] == 'CDS':
            alias = re.search(gff3Comps()['Alias'], entry['attributes'])[1]
            if alias in mrna_set:
                out_cds.append(entry)

    return out_cds

def order_neg_dict( neg_dict ):

    for seqid in neg_dict:
        for gene in neg_dict[seqid]:
            for cds in neg_dict[seqid][gene]:
                cds.sort( reverse = True )
            neg_dict[seqid][gene] = sorted( 
                neg_dict[seqid][gene], key = lambda x: x[0], reverse = True
                )
    
    return neg_dict

def grabCoords( cdss ):
    
    pos_dict, neg_dict = {}, {}
    for entry in cdss:
        gene = re.search( gff3Comps()['Alias'], entry['attributes'] )[1]
        seqid = entry['seqid']
        if entry['strand'].rstrip() == '+':
            if seqid not in pos_dict:
                pos_dict[seqid], pos_dict[seqid][gene] = {}, []
            elif gene not in pos_dict[seqid]:
                pos_dict[seqid][gene] = []
            pos_dict[seqid][gene].append( [
                int(entry['start']) - 1, int(entry['end'])
                ] )
        elif entry['strand'].rstrip() == '-':
            negseq = seqid
            if seqid not in neg_dict:
                neg_dict[seqid], neg_dict[seqid][gene] = {}, []
            elif gene not in neg_dict[seqid]:
                neg_dict[seqid][gene] = []
            neg_dict[seqid][gene].append( [
                int(entry['start']), int(entry['end'])
                ] )

    return pos_dict, order_neg_dict(neg_dict)

def translatePos( contig_seq, postig_dict ):

    genes_fa_dict = {}
    for gene in postig_dict:
        nt_seq = ''.join([
            contig_seq[x[0]:x[1]] for x in postig_dict[gene]
            ])
        genes_fa_dict[gene] = {
            'sequence': str(Seq(nt_seq).translate()), 'description': ''
            }

    return genes_fa_dict

def translateNeg( rev_seq, negtig_dict ):

    genes_fa_dict = {}
    seqlen = len(rev_seq.rstrip())
    for gene in negtig_dict:
        nt_seq = ''.join([
            rev_seq[seqlen-x[0]:seqlen-x[1]+1] for x in negtig_dict[gene]
            ])
        genes_fa_dict[gene] = {
            'sequence': str(Seq(nt_seq).translate()), 'description': ''
            }

    return genes_fa_dict 

def ntPos( contig_seq, postig_dict ):

    genes_fa_dict = {}
    for gene in postig_dict:
        genes_fa_dict[gene] = {
            'sequence': ''.join([
                contig_seq[x[0]:x[1]] for x in postig_dict[gene]
                ]),
            'description': ''
            }

    return genes_fa_dict


def posPlusMinusCode(genePostig_dict, plusminus, contig_seq, plus = True, minus = True):
    geneDict = {'sequence': '', 'description': ''}
    minimumList = []
    for x in genePostig_dict:
        minimumList.extend(x)
        geneDict['sequence'] += contig_seq[x[0]:x[1]]
    minimum = min(minimumList)
    maximum = max(minimumList)
    if minus:
        if minimum:
            if minimum - plusminus > 0:
                geneDict['sequence'] = contig_seq[minimum-plusminus:minimum] + \
                    geneDict['sequence']
                geneDict['description'] = '-' + str(plusminus)
            else:
                geneDict['sequence'] = contig_seq[0:minimum] + \
                    geneDict['sequence']
                geneDict['description'] = '-' + str(minimum)
    if plus:
        geneDict['sequence'] += contig_seq[maximum+1:maximum+plusminus+1]
        if len(contig_seq) < maximum + plusminus:
            geneDict['description'] += '+' + str(len(contig_seq) - maximum)
        else:
            geneDict['description'] += '+' + str(plusminus)
    return geneDict


def posPlusMinus(genePostig_dict, plusminus, contig_seq, plus = True, minus = True):
    minimumList = []
    for x in genePostig_dict:
        minimumList.extend(x)
    minimum = min(minimumList)
    maximum = max(minimumList)

    geneDict = {'sequence': '', 'description': ''}
    if minus:
        if minimum - plusminus > 0:
            geneDict['sequence'] = contig_seq[minimum-plusminus:maximum+1]
            if plusminus:
                geneDict['description'] = '-' + str(plusminus)
        else:
            geneDict['sequence'] = contig_seq[:maximum+1]
            if plusminus:
                geneDict['description'] = '-' + str(minimum)
    if plus:
        geneDict['sequence'] += contig_seq[maximum+1:maximum+1+plusminus]
        if plusminus:
            if maximum + plusminus > len(contig_seq):
                if geneDict['description']:
                    geneDict['description'] += ';'
                geneDict['description'] += '+' + str(len(contig_seq) - maximum)
            else:
                if geneDict['description']:
                    geneDict['description'] += ';'
                geneDict['description'] += '+' + str(plusminus)
    return geneDict

def ntPosNoncode(contig_seq, postig_dict, plusminus = 0):

    genes_fa_dict = {}
    for i, gene in enumerate(list(postig_dict.keys())):
        genes_fa_dict[gene] = {'sequence': '', 'description': ''}
        if plusminus:
            genes_fa_dict[gene] = posPlusMinus(postig_dict[gene], plusminus, contig_seq)
        else:
            minimumList = []
            for x in postig_dict[gene]:
                minimumList.extend(x)
            minimum = min(minimumList)
            maximum = max(minimumList)
    
            genes_fa_dict[gene]['sequence'] = contig_seq[minimum:maximum]
   
    return genes_fa_dict

def negPlusMinusCode(geneNegtig_dict, seqlen, plusminus, rev_seq, plus = True, minus = True):
    minimumList, geneDict = [], {'sequence': '', 'description': ''}
    for x in geneNegtig_dict:
        minimumList.extend(x)
        geneDict['sequence'] += rev_seq[seqlen-x[0]:seqlen-x[1]+1]
    minimum = min(minimumList)
    maximum = max(minimumList)
    if minus:
        if maximum + plusminus < seqlen:
            geneDict['sequence'] = \
                rev_seq[seqlen-maximum-plusminus:seqlen-maximum] + geneDict['sequence']
            geneDict['description'] = '-' + str(plusminus) + ';+'
        else:
            geneDict['sequence'] = \
                rev_seq[:seqlen-maximum] + geneDict['sequence']
            geneDict['description'] = '-' + str(seqlen-maximum) + ';+'
    if plus:
        if minimum - plusminus > 0:
            geneDict['sequence'] += rev_seq[seqlen-minimum+1:seqlen-minimum+plusminus+1]
            geneDict['description'] += str(plusminus)
        else:
            geneDict['sequence'] += rev_seq[seqlen-minimum+1:]
            geneDict['description'] += str(seqlen-minimum)   
    return geneDict

def negPlusMinus(neg_gene, seqlen, plusminus, rev_seq, plus = True, minus = True):
    minimumList = []
    for x in neg_gene:
        minimumList.extend(x)
    minimum = min(minimumList)
    maximum = max(minimumList)

    geneDict = {'sequence': '', 'description': ''}
    if minus:
        if maximum + plusminus < seqlen:
            geneDict['sequence'] = \
                rev_seq[seqlen-maximum-plusminus:seqlen-maximum]
            if plusminus:
                geneDict['description'] = '-' + str(plusminus)
        else:
            geneDict['sequence'] = \
                rev_seq[:seqlen-maximum]
            if plusminus:
                geneDict['description'] = '-' + str(seqlen-maximum)
    if plus:
        if minimum - plusminus > 0:
            geneDict['sequence'] += rev_seq[seqlen-maximum:seqlen-minimum+plusminus+1]
            if plusminus:
                if geneDict['description']:
                    geneDict['description'] += ';'
                geneDict['description'] += '+' + str(plusminus)
        else:
            geneDict['sequence'] += rev_seq[seqlen-maximum+1:]
            if plusminus:
                if geneDict['description']:
                    geneDict['description'] += ';'
                geneDict['description'] += '+' + str(seqlen-minimum)   
    return geneDict


def ntNeg( rev_seq, negtig_dict, plusminus = 0 ):

    genes_fa_dict = {}
    seqlen = len(rev_seq.rstrip())
    for gene in negtig_dict:
        genes_fa_dict[gene] = {
            'sequence': ''.join([
                rev_seq[seqlen-x[0]:seqlen-x[1]+1] for x in negtig_dict[gene]
                ]),
            'description': ''
            }

    return genes_fa_dict 

def ntNegNoncode(rev_seq, negtig_dict, plusminus = 0):

    genes_fa_dict = {}
    seqlen = len(rev_seq.rstrip())
    for gene in negtig_dict:
        if plusminus:
            genes_fa_dict[gene] = negPlusMinus(negtig_dict[gene], plusminus, seqlen, rev_seq)
        else:
            minimumList = []
            for x in negtig_dict[gene]:
                minimumList.extend(x)
            minimum = min(minimumList)
            maximum = max(minimumList)

            genes_fa_dict[gene] = {
                'sequence': rev_seq[seqlen - maximum:seqlen - minimum + 1],
                'description': ''
                }

    return genes_fa_dict


def ntmain(gff_dicts, assem_dict, coding = True, 
           flanks = True, fullRegion = False, plusminus = 0,
           spacer = '\t'):
    
    if fullRegion:
        flanks, coding = True, False
    contig_dict, contig_info, genes_fa_dict, geneOrder = {}, {}, {}, {}
    cdss = grabCDS( sortMain(gff_dicts), spacer = spacer )
    for cds in cdss:
        seqid = cds['seqid']
        if seqid not in contig_dict:
            contig_dict[seqid] = []
            geneOrder[seqid] = []
        contig_dict[seqid].append(cds)
        gene = re.search(gff3Comps()['Alias'], cds['attributes'])[1]
        if gene not in set(geneOrder[seqid]):
            geneOrder[seqid].append(gene)
    for seqid in contig_dict:
        contig_info[seqid] = [(
            re.search(gff3Comps()['Alias'], contig_dict[seqid][0]['attributes'])[1],
            contig_dict[seqid][0]['strand']
            ),(
            re.search(gff3Comps()['Alias'], contig_dict[seqid][-1]['attributes'])[1],
            contig_dict[seqid][-1]['strand']
            )]

    pos_dict, neg_dict = grabCoords( cdss )

    contig_fa_dict, startFlanks, endFlanks = {seqid: {} for seqid in contig_info}, {}, {}
    if flanks:
        if coding:
            for seqid in contig_info:
                startGene, startStrand = contig_info[seqid][0][0], contig_info[seqid][0][1]
                endGene, endStrand = contig_info[seqid][1][0], contig_info[seqid][1][1]
                if startGene == endGene:
                    if startStrand == '+':
                        startFlanks[seqid] = posPlusMinusCode(
                            pos_dict[seqid][startGene],
                            plusminus,
                            assem_dict[seqid]['sequence'],
                            )
                        del pos_dict[seqid][startGene]
                    else:
                        rev_comp = str(Seq(assem_dict[seqid]['sequence']).reverse_complement())
                        startFlanks[seqid] = negPlusMinusCode(
                            neg_dict[seqid][startGene], len(rev_comp), plusminus, rev_comp
                            )
                        del neg_dict[seqid][startGene]
                    continue
                if startStrand == '+':
                    startFlanks[seqid] = posPlusMinusCode(
                        pos_dict[seqid][startGene],
                        plusminus,
                        assem_dict[seqid]['sequence'],
                        plus = False
                        )
                    del pos_dict[seqid][startGene]
                else:
                    rev_comp = str(Seq(assem_dict[seqid]['sequence']).reverse_complement())
                    startFlanks[seqid] = negPlusMinusCode(
                        neg_dict[seqid][startGene], len(rev_comp), plusminus, rev_comp, plus = False
                        )
                    del neg_dict[seqid][startGene]
                if endStrand == '+':
                    endFlanks[seqid] = posPlusMinusCode(
                        pos_dict[seqid][endGene],
                        plusminus,
                        assem_dict[seqid]['sequence'],
                        minus = False
                        )
                    del pos_dict[seqid][endGene]
                else:
                    rev_comp = str(Seq(assem_dict[seqid]['sequence']).reverse_complement())
                    endFlanks[seqid] = negPlusMinusCode(
                        neg_dict[seqid][endGene], len(rev_comp), plusminus, rev_comp, minus = False
                        )
                    del neg_dict[seqid][endGene]
        else:
            for seqid in contig_info:
                startGene, startStrand = contig_info[seqid][0][0], contig_info[seqid][0][1],
                endGene, endStrand = contig_info[seqid][1][0], contig_info[seqid][1][1]
                      
                if startGene == endGene:
                    if startStrand == '+':
                        startFlanks[seqid] = posPlusMinus(
                            pos_dict[seqid][startGene],
                            plusminus,
                            assem_dict[seqid]['sequence'],
                            )
                        del pos_dict[seqid][startGene]
                    else:
                        rev_comp = str(Seq(assem_dict[seqid]['sequence']).reverse_complement())
                        startFlanks[seqid] = negPlusMinus(
                            neg_dict[seqid][startGene], len(rev_comp), plusminus, rev_comp
                            )
                        del neg_dict[seqid][startGene]
                    continue

                elif fullRegion:
                    region = []
                    if startStrand == '+':
                        region.append(min(pos_dict[seqid][startGene]))
                    else:
                        region.append(min(neg_dict[seqid][startGene]))
                    if endStrand == '+':
                        region.append(max(pos_dict[seqid][endGene]))
                    else:
                        region.append(max(neg_dict[seqid][endGene]))
                    name = startGene + '-' + endGene
                    genes_fa_dict[name + '_sense'] = posPlusMinus(
                        region, plusminus, assem_dict[seqid]['sequence']
                        )
                    rev_comp = str(Seq(assem_dict[seqid]['sequence']).reverse_complement())
                    genes_fa_dict[name + '_antisense'] = negPlusMinus(
                        region, len(rev_comp), plusminus, rev_comp
                        )

                
                if startStrand == '+':
                    startFlanks[seqid] = posPlusMinus(
                        pos_dict[seqid][startGene],
                        plusminus,
                        assem_dict[seqid]['sequence'],
                        plus = False
                        )
                    del pos_dict[seqid][startGene]
                else:
                    rev_comp = str(Seq(assem_dict[contig]['sequence']).reverse_complement())
                    startFlanks[seqid] = negPlusMinus(
                        neg_dict[seqid][startGene], len(rev_comp), plusminus, rev_comp, plus = False
                        )
                    del neg_dict[seqid][startGene]
                if endStrand == '+':
                    endFlanks[seqid] = posPlusMinus(
                        pos_dict[seqid][endGene],
                        plusminus,
                        assem_dict[seqid]['sequence'],
                        minus = False
                        )
                    del pos_dict[seqid][endGene]
                else:
                    rev_comp = str(Seq(assem_dict[seqid]['sequence']).reverse_complement())
                    endFlanks[seqid] = negPlusMinus(
                        neg_dict[seqid][endGene], len(rev_comp), plusminus, rev_comp, minus = False
                        )
                    del neg_dict[seqid][endGene]

        if coding:
            for contig in contig_fa_dict:
                if contig in pos_dict:
                    contig_fa_dict[contig] = { 
                        **ntPos(assem_dict[contig]['sequence'], pos_dict[contig]), 
                        **contig_fa_dict[contig]
                        }
                if contig in neg_dict:
                    rev_comp = str(Seq(assem_dict[contig]['sequence']).reverse_complement())
                    contig_fa_dict[contig] = {
                        **ntNeg(rev_comp, neg_dict[contig]),
                        **contig_fa_dict[contig]
                        }
        elif not fullRegion:
            for contig in contig_fa_dict:
                if contig in pos_dict:
                    contig_fa_dict[contig] = { 
                        **ntPosNoncode(assem_dict[contig]['sequence'], pos_dict[contig]), 
                        **contig_fa_dict[contig]
                        }
                if contig in neg_dict:
                    rev_comp = str(Seq(assem_dict[contig]['sequence']).reverse_complement())
                    contig_fa_dict[contig] = {
                        **ntNegNoncode(rev_comp, neg_dict[contig]),
                        **contig_fa_dict[contig]
                        }
    if not fullRegion:
        for contig in geneOrder:
            if contig in startFlanks:
                genes_fa_dict[geneOrder[contig][0]] = startFlanks[contig]
            for gene in geneOrder[contig][1:-1]:
                genes_fa_dict[gene] = contig_fa_dict[contig][gene]
            if contig in endFlanks:
                genes_fa_dict[geneOrder[contig][-1]] = endFlanks[contig]
    
    return genes_fa_dict

def aamain( gff_dicts, assem_dict, spacer = '\t' ):
    
    cdss = grabCDS(gff_dicts, spacer)
    pos_dict, neg_dict = grabCoords(cdss)


    genes_fa_dict = {}
    for contig in pos_dict:
        genes_fa_dict = { 
            **translatePos(assem_dict[contig]['sequence'], pos_dict[contig]), 
            **genes_fa_dict 
            }
    for contig in neg_dict:
        rev_comp = str(Seq(assem_dict[contig]['sequence']).reverse_complement())
        genes_fa_dict = {
            **translateNeg(rev_comp, neg_dict[contig]),
            **genes_fa_dict
            }

    return genes_fa_dict

            
def cli():

    parser = argparse.ArgumentParser(
       description = 'Inputs MycoDB compatible gff3, assembly (optional), ' + \
           'and outputs nucleotides/proteins. Use an abstracted gene gff ' + \
           '(acc2gff.py) to only output a smaller set of gene(s).'
       )
    parser.add_argument('-g', '--gff', help = '"-" for stdin', required = True)
    parser.add_argument('-n', '--nucleotide', action = 'store_true')
    parser.add_argument('-p', '--protein', action = 'store_true')
    parser.add_argument('-a', '--assembly')
    parser.add_argument('-i', '--intergenic', action = 'store_true', help = '-n only')
    parser.add_argument('-nc', '--noncoding', action = 'store_true', help = '-n only')
    parser.add_argument('-pm', '--plusminus', default = 0, type = int, help = '-n only')
    parser.add_argument('-af', '--all_flanks', action = 'store_true', help = '-n and -nc only')
    args = parser.parse_args()

    if args.gff == '-':
        data = stdin2str()
        input_gff = gff2list(data, path = False)
    else:
        input_gff = gff2list(format_path(args.gff))
    if args.assembly:
        assembly_dicts = {'input': fa2dict( format_path(args.assembly) )}
        gff_dicts = {'input': input_gff}
    else:
        db = mtdb(primaryDB()).set_index('ome')
        gff_dicts, assembly_dicts = {}, {}
        try:
            for line in input_gff:
                gene = re.search(gff3Comps()['Alias'], line['attributes'])[1]
                ome = re.search( r'(.*?)_', gene )[1]
                if ome not in gff_dicts:
                    gff_dicts[ome] = []
                    assembly_dicts[ome] = fa2dict(db[ome]['fna'])
                gff_dicts[ome].append(line)
        except IndexError:
            eprint('\nERROR: ' + args.gff + ' is incompatible with MycotoolsDB', flush = True)
            sys.exit(1)
        
    if args.protein:
        for ome in gff_dicts:
            print(dict2fa(aamain(gff_dicts[ome], assembly_dicts[ome])), flush = True)
    else:
        for ome in assembly_dicts:
            print(dict2fa(
                ntmain(gff_dicts[ome], assembly_dicts[ome], not args.noncoding, not args.all_flanks, args.intergenic, args.plusminus)
                ), flush = True)

    sys.exit(0)


if __name__ == '__main__':
    cli()
