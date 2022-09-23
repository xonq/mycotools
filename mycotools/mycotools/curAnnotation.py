#! /usr/bin/env python3

# NEED TO DEFAULT CHEKC FOR INTRONS FOR PREDB2DB
# NEED to arrive at a consensus for protein IDs
# NEED to name tRNAs, pseudogenes, etc. similar to curGFF.py

import os
import re
import sys
import copy
import argparse
from mycotools.lib.biotools import gff2list, list2gff, fa2dict, dict2fa, \
    gtfComps, gff3Comps, gff2Comps
from mycotools.lib.kontools import collect_files, eprint, format_path
from mycotools.gff2seq import aamain as gff2proteome


def grabOutput( output_pref ):
    '''
    Inputs: orthofiller `output_path` for results
    Outputs: gff_dict and fasta_dict of results
    Collect the files in `output_path`, find the proteome and gtf, and make
    sure there is only one possible entry. Return the dictionary versions of
    both.
    '''

    output_path, pref = os.path.dirname(output_pref), os.path.basename(output_pref)
    output_files = collect_files( output_path, '*' )
    hits = [x for x in output_files if os.path.basename(x).startswith(pref)]
    proteome = [x for x in hits if x.endswith('.results.aa.fasta')]
    gtf = [x for x in hits if x.endswith( '.results.gtf' )]
    if len( gtf ) != 1 or len( proteome ) != 1:
        if len( gtf ) < 1 or len( proteome ) < 1:
            print( '\nERROR: complete output files not detected' , flush = True)
        else:
            print( '\nERROR: multiple eligible complete output detected' , flush = True)
        sys.exit( 3 )

    return gff2list( gtf[0] ), fa2dict( proteome[0] )


def intron2exon( gff, gene_comp = re.compile( r'gene_id \"(.*?)\"' ) ):
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

    gff1, gene_info = [], {}
    for entry in gff:
        if int(entry['start']) > int(entry['end']):
            start, end = copy.deepcopy(entry['start']), copy.deepcopy(entry['end'])
            entry['start'], entry['end'] = end, start
        gene = gene_comp.search( entry['attributes'] )[1]
        if gene not in gene_info:
            gene_info[ gene ] = { 
                'intron': [], 
                'start_codon': [], 
                'stop_codon': [], 
                'strand': str(entry['strand']), 
                'raw': [] 
            }
        if entry[ 'type' ] in gene_info[ gene ]:
            gene_info[ gene ][ entry['type'] ].append( 
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
        exon_coords[ gene ] = []
        intron_genes[gene]['intron'].sort( key = lambda x: x[0] )
    
        if intron_genes[gene]['strand'] == '+':
            exon_coords[gene].append( [intron_genes[gene]['start_codon'][0][0]] )
            for intron in intron_genes[gene]['intron']:
                exon_coords[gene][-1].append( intron[0] - 1 )
                exon_coords[gene].append([intron[1] + 1])
            exon_coords[gene][-1].append( intron_genes[gene]['stop_codon'][0][1] )

        else:
            exon_coords[gene].append( [intron_genes[gene]['stop_codon'][0][0]] )
            for intron in intron_genes[gene]['intron']:
                exon_coords[gene][-1].append( intron[0] - 1 )
                exon_coords[gene].append([intron[1] + 1])
            exon_coords[gene][-1].append( intron_genes[gene]['start_codon'][0][1] )
#        elif intron_genes[gene]['start_codon']:
 #           if intron_genes[gene]['strand'] == '+':
  #              exon_coords[gene] = [[
   #                 intron_genes[gene]['start_codon'][0][0], intron_genes[gene]['stop_codon'][0][1]
    #                ]]
     #       else:
      #          exon_coords[gene] = [[
       #             intron_genes[gene]['stop_codon'][0][0], intron_genes[gene]['start_codon'][0][1]
        #            ]]

    gff2 = []
    for entry in gff1:
        gene = gene_comp.search( entry['attributes'] )[1]
        if gene not in exon_coords:
            gff2.append(entry)
    for gene in exon_coords:
        add_data = [ 
            entry for entry in intron_genes[gene]['raw'] if entry['type'] != 'intron' 
        ]
        scaffold = list(add_data)[0]
        for exon in exon_coords[gene]:
            new_entry = dict( scaffold )
            new_entry['type'] = 'exon'
            new_entry['start'] = str( exon[0] )
            new_entry['end'] = str( exon[1] )
            new_entry['score'] = '.'
            new_entry['phase'] = '.'
            add_data.append( new_entry )
        gff2.extend( add_data )

    return gff2


def curCDS( gff ):

    new_gff, info_dict = [], {}
    gene_compile = re.compile( r'gene_id \"(.*?)\"' )
    for entry in gff:
#        if entry['source'] != 'AUGUSTUS':
#or entry['strand'] == '+':
 #           new_gff.append(entry)
  #          continue
        if not gene_compile.search(entry['attributes']):
            gene_compile = re.compile( r'ID=(.*?);' )
        gene = gene_compile.search(entry['attributes'])[1]
        gene = re.sub( r'\-T\d+.*$', '', gene )

        if gene not in info_dict:
            info_dict[gene] = { 'start_codon': [], 'stop_codon': [], 'raw': [] }
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
#        ready, ready1 = False, False
        cop = False
#        info_dict[gene]['start_codon'].sort()
 #       info_dict[gene]['stop_codon'].sort()
#        end0 = int(info_dict[gene]['start_codon'][1])
        try:
            start0 = int(info_dict[gene]['start_codon'][0])
            end1 = int(info_dict[gene]['stop_codon'][1])
        except IndexError: # need to flag here
            continue
        if start0 > end1:
            start0 = int(info_dict[gene]['stop_codon'][0])
            end1 = int(info_dict[gene]['start_codon'][1])
 #       start1 = int(info_dict[gene]['stop_codon'][0])
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

                    
#                    info_dict[gene]['raw'][index]['end'] = new
 #                   if ready:
  #                      ready1 = True
   #                     break
    #                ready = True
#don't need once fixed also don't need the readies
#                    cop = copy.deepcopy(info_dict[gene]['raw'][index])
# don't need once fixed
 #               elif int(entry['end']) == new and entry['strand'] == '+':
  #                  if ready:
   #                     ready1 = True
    #                    break
     #               ready = True
      #              cop = copy.deepcopy(info_dict[gene]['raw'][index])

# fixes 1 exon from intron2exon script - needs to be removed once fixed
       # if not ready1:
        #    if cop and cop['type'] == 'CDS':
         #       cop['type'], cop['score'], cop['phase'] = 'exon', '.', '.'
          #      cop['attributes'] = cop['attributes'].replace('cds', 'exon')
           #     info_dict[gene]['raw'].append(cop)


    return new_gff


def liberalRemoval( gene_dict_prep, contigs ):

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
            gene_dict[ gene ] = gene_dict_prep[gene]
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


def conservativeRemoval( gene_dict_prep ):

    gene_dict, flagged, failed = {}, [], []
    for gene in gene_dict_prep:
        temp = gene_dict_prep[gene]
        try:
            if not temp['start_codon'] \
                or not temp['stop_codon']:
                temp['exon'].sort()
                if not temp['start_codon']:
                    if temp['strand'] == '+':
                        temp['start_codon'] = [
                            temp['exon'][0],
                            temp['exon'][0] + 2
                            ]
                    else:
                        temp['start_codon'] = [
                            temp['exon'][-1] - 2,
                            temp['exon'][-1]
                            ]
                if not temp['stop_codon']:
                    if temp['strand'] == '+':
                        temp['stop_codon'] = [
                            temp['exon'][-1] - 2,
                            temp['exon'][-1]
                            ]
                    else:
                        temp['stop_codon'] = [
                            temp['exon'][0],
                            temp['exon'][0] + 3
                            ]
                flagged.append( gene )

        except IndexError:
            eprint( gene + ' cannot create gene coordinates' , flush = True)
            continue

        gene_dict[gene] = temp

    
    return gene_dict, flagged, failed
        


def addGenes( gtf, safe = True ):

    comps = gtfComps()
    gene_compile = re.compile( comps['id'] )
    if not gene_compile.search( gtf[0]['attributes'] ):
        comps = gff2Comps()
        gene_compile = re.compile( comps['id'] )
    gene_dict_prep, gene_dict, contigs = {}, {}, {}
#    gene_compile = re.compile( r'gene_id "(.*?)";')
    for entry in gtf:
        gene = gene_compile.search( entry['attributes'] )[1]

        contig = entry['seqid']
        if gene not in gene_dict_prep:
            gene_dict_prep[ gene ] = { 
                'start_codon': [], 'stop_codon': [], 'exon': [],
                'strand': str(entry['strand']), 'contig': contig,
                'rna': None, 'cds': []
                }
            if contig not in contigs:
                contigs[contig] = {}
            contigs[contig][gene] = []

        if entry['type'].lower() in gene_dict_prep[gene]:
            contigs[contig][gene].extend( 
                [int(entry['start']), int(entry['end'])]
                )
            gene_dict_prep[gene][entry['type'].lower()].extend(
                [int(entry['start']), int(entry['end'])]
                )
        elif entry['type'] in {'mRNA', 'tRNA', 'rRNA'}:
            gene_dict_prep[gene]['rna'] = entry['type']
   
    flagged, failed = [], []
    if safe:
        gene_dict, flagged, failed = conservativeRemoval( gene_dict_prep )
    else:
        gene_dict, failed = liberalRemoval( gene_dict_prep, contigs )

    check, insert_list = set(), []
    for index, entry in enumerate(gtf):
        gene = gene_compile.search( entry['attributes'] )[1]
        if gene not in check and gene in gene_dict:
            new_entry = dict(entry)
            new_entry['type'] = 'gene'
            new_entry['phase'], new_entry['score'] = '.', '.'
            if comps['ver'] == 'gff2':
                new_entry['attributes'] = 'name "' + gene + '";'
            else:
                new_entry['attributes'] = 'gene_id "' + gene + '";'
            if gene_dict[ gene ][ 'start_codon' ][0] > gene_dict[gene]['start_codon'][1]:
                gene_dict[gene]['start_codon'] = [
                    gene_dict[gene]['start_codon'][1], 
                    gene_dict[gene]['start_codon'][0]
                    ]
                gene_dict[gene]['stop_codon'] = [
                    gene_dict[gene]['stop_codon'][1], 
                    gene_dict[gene]['stop_codon'][1]
                    ]
            if gene_dict[gene]['strand'] == '+':
                new_entry['start'] = gene_dict[gene]['start_codon'][0]
                new_entry['end'] = gene_dict[gene]['stop_codon'][1]
            else:
                new_entry['start'] = gene_dict[gene]['stop_codon'][0]
                new_entry['end'] = gene_dict[gene]['start_codon'][1]
            if gene_dict_prep[gene]['rna']:
                new2 = copy.deepcopy( new_entry )
                new2['type'] = gene_dict_prep[gene]['rna']
                insert_list.extend( [[index, new2], [index, new_entry]] )
            elif not gene_dict_prep[gene]['cds']:
                new2 = copy.deepcopy(new_entry)
                new2['type'] = 'tRNA' # assume tRNA, but not valid to discredit rRNA
                insert_list.extend([[index, new2], [index, new_entry]])
            else:
                new2 = copy.deepcopy(new_entry)
                new2['type'] = 'mRNA'
                insert_list.extend([[index, new2], [index, dict(new_entry)]])
            check.add( gene )

    insert_list.sort( key = lambda x: x[0], reverse = True )
    for insert in insert_list:
        gtf.insert( insert[0], insert[1] )


    return gtf, failed, flagged


def removeStartStop( gtf ):

    out_gtf = []
    for entry in gtf:
        if entry['type'] not in { 'start_codon', 'stop_codon' }:
            out_gtf.append( entry )

    return out_gtf


def curate( gff, prefix, failed = set() ):

    for index in range(len(gff)):
        if gff[index]['type'].lower() == 'exon':
            break

    comps = gtfComps()
    geneComp = re.compile(comps['id'])
    if geneComp.search( gff[0]['attributes'] ) is None:
        comps = gff3Comps()
        geneComp = re.compile(comps['id'])

    # geneComp = comps['id']
    #acquire the start and end coordinates for all entries with the same gene attribute
    temp_exon_dict, rna_info = {}, {}
    for line in gff:
        if line['type'] in {'exon', 'start_codon', 'stop_codon', 'mRNA', 'CDS', 'tRNA', 'rRNA'}:
            if line['seqid'] not in temp_exon_dict:
                temp_exon_dict[ line['seqid'] ] = {}
            prot = geneComp.search( line['attributes'] )[1]
            if comps == gff3Comps():
                protsub = re.compile( r'[\.|-].*' ) # funannotate remove transcript tag
                prot = protsub.sub( '', prot )
            if prot not in temp_exon_dict[ line['seqid'] ]:
                temp_exon_dict[ line['seqid'] ][ prot ] = []
            if line['type'] in {'mRNA', 'tRNA', 'rRNA'}:
                rna_info[prot] = line['type']
            temp_exon_dict[ line['seqid'] ][ prot ].extend( [ int( line['start'] ), int( line['end'] ) ] )
    exon_dict = {}
    for key, value in sorted( temp_exon_dict.items() ):
        exon_dict[ key ] = value

    count = 1
    change_dict = {}
    for seq in exon_dict:
        for prot in exon_dict[ seq ]:
            exon_dict[ seq ][ prot ].sort()
        exon_dict[ seq ] = dict(sorted( exon_dict[ seq ].items(), key = lambda e: e[1][0] ))
        for prot in exon_dict[ seq ]:
            new_prot = re.sub(r'-T\d+$', '', prot)
            if new_prot not in change_dict:
                change_dict[ new_prot ] = prefix + '_' + str(count)
            count += 1

    crudesortGff = preSortGFF(gff, re.compile(comps['id']))

    newGff, trans_set, exonCheck, cdsCheck, failed = [], set(), {}, {}, set( x[0] for x in failed )
    transComp = re.compile( r'transcript_id "(.*?)"\;' )
    # transComp = comps['transcript']
    for entry in crudesortGff:
        prot = geneComp.search( entry['attributes'] )
        trans = transComp.search( entry['attributes'] )
        prot = prot[1]
        if prot in failed:
            continue
        if comps['ver'] == 'gff3':
            prot = protsub.sub('', prot)
        prot = change_dict[ prot ]
        if entry['type'] in {'mRNA', 'tRNA', 'rRNA'}:
            count = 1
            transProt = prot + '-T' + str( count )
            while transProt in trans_set:
                count += 1
                transProt = prot + '-T' + str( count )
            if entry['type'] == 'mRNA':
                entry['attributes'] = 'ID=' + transProt + ';Parent=' + prot + ';product=' + \
                    'hypothetical protein;Alias=' + prot
            else:
                entry['attributes'] = 'ID=' + transProt + ';Parent=' + prot + ';Alias=' + prot
            trans_set.add( transProt )
        elif entry['type'] == 'exon':
            if transProt not in exonCheck:
                exonCheck[ transProt ] = 0
            exon = exonCheck[ transProt ] + 1
            entry['attributes'] = 'ID=' + transProt + '.exon' + str(exon) + ';Parent=' + \
                transProt + ';Alias=' + prot
            exonCheck[transProt] += 1
        elif entry['type'] == 'CDS':
            if transProt not in cdsCheck:
                cdsCheck[ transProt ] = 0
            cds = cdsCheck[transProt] + 1
            entry['attributes'] = 'ID=' + transProt + '.cds' + str(cds) + \
                ';Parent=' + transProt + ';' + 'Alias=' + prot 
            cdsCheck[transProt] += 1
        else: #if entry['type'] == 'gene':
            entry['attributes'] = 'ID=' + prot + ';Alias=' + prot

        if int(entry['end']) < int(entry['start']):
            start = copy.deepcopy(entry['end'])
            end = copy.deepcopy(entry['start'])
            entry['start'], entry['end'] = start, end
        newGff.append( entry )

    translation_str = ''
    for entry in change_dict:
        translation_str += entry + '\t' + change_dict[ entry ] + '\n'

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
        

def addExons( gff ):

    exonCheck = {}
    for i in range(len(gff)):
        entry = gff[i]
        if entry['source'] == 'AUGUSTUS':
            if entry['type'] in { 'CDS', 'exon' }:
                prot = re.search( r'ID=(.*?_\d+)', entry['attributes'] )[1]
                if prot not in exonCheck:
                    exonCheck[ prot ] = [False, i, entry]
                if entry['type'] == 'exon':
                    exonCheck[prot] = [True, i, None]
           
    prots = reversed(list(exonCheck.keys())) 
    for prot in prots:
        if not exonCheck[prot][0]:
            newEntry = copy.deepcopy(exonCheck[prot][2])
            newEntry['attributes'] = newEntry['attributes'].replace('.cds;', '.exon1;')
            newEntry['type'] = 'exon'
            gff.insert( exonCheck[prot][1] + 1, newEntry )

    return gff


def main(gff_path, prefix, fail = True):

    if isinstance(gff_path, str):
        gff = gff2list(gff_path)
    elif isinstance(gff_path, list):
        gff = gff_path

    if re.search(gtfComps()['id'], gff[0]['attributes']) is not None:
        exonGtf = intron2exon( gff )
        exonGtfCur = curCDS( exonGtf )
        exonGtfCurGenes, failed, flagged = addGenes( exonGtfCur, safe = fail )
        preGff = removeStartStop( exonGtfCurGenes )
        gffUncur, trans_str = curate( preGff, prefix, failed ) 
        unsortedGff = addExons( gffUncur )
    else:
        unsortedGff, trans_str = curate( 
            gff, prefix
            ) 
        failed, flagged = None, None

#    id_comp = re.compile(gff3Comps()['id'])
    crude_sort = sorted( unsortedGff, key = lambda x: \
        int( re.search(r'ID=' + prefix + '_(\d+)', x['attributes'])[1] ))
    gff = sortGFF(crude_sort, re.compile(gff3Comps()['Alias']))

    return gff, trans_str, failed, flagged


def sortMain(gff, prefix):

    id_comp = re.compile(gff3Comps()['id'])
    crude_sort = sorted( gff, key = lambda x: \
        int( re.search(r'ID=' + prefix + '_(\d+)', x['attributes'])[1] ))
    gff = preSortGFF(crude_sort, id_comp)

    return gff


if __name__ == '__main__':

    parser = argparse.ArgumentParser( description = 'Curates Funannotate or ' + \
        'post OrthoFiller output by naming accessions as <PREFIX>_####, where ' + \
        '#### is the is the index from ordered accessions. If you are planning' + \
        ' on using OrthoFiller, you may want to wait to not mix accessions.' )
    parser.add_argument( '-g', '--gff', required = True, help = '.gtf, .gff/.gff3' )
    parser.add_argument( 
        '-a', '--assembly', 
        help = 'Assembly fasta for proteome export' 
        )
    parser.add_argument( '-p', '--prefix', required = True,
        help = 'Accession prefix' )
    parser.add_argument( '-s', '--sort', action = 'store_true',
        help = 'Only sort compatible gff3' )
    parser.add_argument( '-o', '--output', help = 'Output directory' )
    parser.add_argument( '--fail', default = True, action = 'store_false',
        help = 'Fail genes without start or stop codons.' )
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
        with open( output + '.gff3', 'w' ) as out:
            out.write( list2gff( gff ) + '\n' )
        sys.exit(0)
    elif not args.assembly:
        eprint('\nERROR: assembly (`-f`) required for curation\n', flush = True)
        sys.exit(2)

    if '_' in args.prefix:
        eprint('\nERROR: "_" not allowed in prefix\n', flush = True)
        sys.exit(1)

    gff, trans_str, failed, flagged = main(
        format_path(args.gff), args.prefix, args.fail
        )
    if args.assembly:
        fna = fa2dict(format_path(args.assembly))
        faa = gff2proteome( gff, fna )
        with open( output + '.faa', 'w' ) as out:
            out.write( dict2fa( faa ) )

    with open( output + '.gff3', 'w' ) as out:
        out.write( list2gff( gff ) + '\n' )
    with open( output + '.transitions', 'w' ) as out:
        out.write( trans_str )
    if failed:
        print('\n' + str(len(failed)) + ' failures\n', flush = True)
        with open( output + '.failed', 'w' ) as out:
            out.write('\n'.join(['\t'.join(x) for x in failed]))
    if flagged:
        print('\n' + str(len(flagged)) + ' flagged\n', flush = True)
        with open( output + '.flagged', 'w' ) as out:
            out.write('\n'.join(flagged))

    sys.exit(0)
