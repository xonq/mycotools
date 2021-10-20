#! /usr/bin/env python3

import sys, re, os, copy, argparse
from mycotools.lib.fastatools import gff2dict, dict2gff, fasta2dict, dict2fasta, \
    gtfComps, gff3Comps, gff2Comps
from mycotools.lib.kontools import collect_files, eprint, formatPath
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

    return gff2dict( gtf[0] ), fasta2dict( proteome[0] )


def intron2exon( gff ):
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

    new_gff, gene_info = [], {}
    gene_comp = re.compile( r'gene_id \"(.*?)\"' )
    for entry in gff:
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
        gene_info[gene]['raw'].append( dict(entry) )

    intron_genes = { 
        gene: gene_info[gene] for gene in gene_info if len(gene_info[gene]['intron']) > 0 
        }
    exon_coords = {}
    for gene in intron_genes:
        exon_coords[ gene ] = []
        if intron_genes[gene]['start_codon'][0][0] > intron_genes[gene]['start_codon'][0][1]:
            intron_genes[gene]['start_codon'] = [
                [
                    intron_genes[gene]['start_codon'][0][1], 
                    intron_genes[gene]['start_codon'][0][0]
                ]
            ] 
            intron_genes[gene]['stop_codon'] = [
                [
                    intron_genes[gene]['stop_codon'][0][1], 
                    intron_genes[gene]['stop_codon'][0][0]
                ]
            ]
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
            exon_coords[gene][-1].append( intron_genes[gene]['start_codon'][0][0] )

    for entry in gff:
        gene = gene_comp.search( entry['attributes'] )[1]
        if gene not in exon_coords:
            new_gff.append( dict(entry) )
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
        new_gff.extend( add_data )

    return new_gff


def curCDS( gff ):

    new_gff, info_dict = [], {}
    gene_compile = re.compile( r'gene_id \"(.*?)\"' )
    for entry in gff:
        if entry['source'] != 'AUGUSTUS':
#or entry['strand'] == '+':
            new_gff.append(dict(entry))
            continue
        if not gene_compile.search(entry['attributes']):
            gene_compile = re.compile( r'ID=(.*?);' )
        gene = gene_compile.search(entry['attributes'])[1]
        gene = re.sub( r'\-T\d+.*$', '', gene )
        if gene not in info_dict:
            info_dict[gene] = { 'start_codon': [], 'raw': [] }
        if entry['type'] == 'gene':
            info_dict[gene]['start_codon'].extend([int(entry['end']) - 2, int(entry['end'])])
        elif entry['type'] == 'start_codon':
            info_dict[gene][entry['type']].extend([int(entry['start']), int(entry['end'])])
        info_dict[gene]['raw'].append(copy.deepcopy(entry))

    for gene in info_dict:
        ready, ready1 = False, False
        cop = False
        info_dict[gene]['start_codon'].sort()
        new = int(info_dict[gene]['start_codon'][1])
        old = int(info_dict[gene]['start_codon'][0])
        for index in range(len(info_dict[gene]['raw'])):
            entry = info_dict[gene]['raw'][index]
            if entry['type'] == 'CDS' or entry['type'] == 'exon':
                if int(entry['end']) == old:
                    info_dict[gene]['raw'][index]['end'] = new
                    if ready:
                        ready1 = True
                        break
                    ready = True
#don't need once fixed also don't need the readies
                    cop = copy.deepcopy(info_dict[gene]['raw'][index])
# don't need once fixed
                elif int(entry['end']) == new and entry['strand'] == '+':
                    if ready:
                        ready1 = True
                        break
                    ready = True
                    cop = copy.deepcopy(info_dict[gene]['raw'][index])

# fixes 1 exon from intron2exon script - needs to be removed once fixed
        if not ready1:
            if cop and cop['type'] == 'CDS':
                cop['type'], cop['score'], cop['phase'] = 'exon', '.', '.'
                cop['attributes'] = cop['attributes'].replace('cds', 'exon')
                info_dict[gene]['raw'].append(cop)

        new_gff.extend( info_dict[gene]['raw'] )

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
                'strand': str(entry['strand']), 'contig': contig
                }
            if contig not in contigs:
                contigs[contig] = {}
            contigs[contig][gene] = []

        if entry['type'] in gene_dict_prep[gene]:
            contigs[contig][gene].extend( 
                [int(entry['start']), int(entry['end'])]
                )
            gene_dict_prep[gene][entry['type']].extend(
                [int(entry['start']), int(entry['end'])]
                )
   
    flagged, failed = [], []
    if safe:
        gene_dict, flagged, failed = conservativeRemoval( gene_dict_prep )
    else:
        gene_dict, failed = liberalRemoval( gene_dict_prep, contigs )


    check, insert_list = set(), []
    for index in range(len(gtf)):
        entry = dict(gtf[index])
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
            new2 = copy.deepcopy( new_entry )
            new2['type'] = 'mRNA'
            insert_list.extend( [[index, dict(new2)], [index, dict(new_entry)]] )
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
    temp_exon_dict = {}
    for line in gff:
        if line['type'] in { 'exon', 'start_codon', 'stop_codon', 'mRNA', 'CDS' }:
            if line['seqid'] not in temp_exon_dict:
                temp_exon_dict[ line['seqid'] ] = {}
            prot = geneComp.search( line['attributes'] )[1]
            if comps == gff3Comps():
                protsub = re.compile( r'[\.|-].*' )
                prot = protsub.sub( '', prot )
            if prot not in temp_exon_dict[ line['seqid'] ]:
                temp_exon_dict[ line['seqid'] ][ prot ] = []
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
            new_prot = re.sub( r'-T\d+$', '', prot )
            if new_prot not in change_dict:
                change_dict[ new_prot ] = prefix + '_' + str(count)
            count += 1

    newGff, trans_set, exonCheck, cdsCheck, failed = [], set(), {}, {}, set( x[0] for x in failed )
    transComp = re.compile( r'transcript_id "(.*?)"\;' )
    # transComp = comps['transcript']
    for entry in gff:
        prot = geneComp.search( entry['attributes'] )
        trans = transComp.search( entry['attributes'] )
        prot = prot[1]
        if prot in failed:
            continue
        if comps['ver'] == 'gff3':
            prot = protsub.sub('', prot)
        prot = change_dict[ prot ]
        if entry['type'] == 'gene':
            entry['attributes'] = 'ID=' + prot + ';Alias=' + prot
        elif entry['type'] == 'mRNA':
            count = 1
            transProt = prot + '-T' + str( count )
            while transProt in trans_set:
                count += 1
                transProt = prot + '-T' + str( count )
            entry['attributes'] = 'ID=' + transProt + ';Parent=' + prot + ';product=' + \
                'hypothetical protein;Alias=' + prot
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
        newGff.append( entry )

    translation_str = ''
    for entry in change_dict:
        translation_str += entry + '\t' + change_dict[ entry ] + '\n'

    newGff = sorted( newGff, key = lambda x: \
        int( re.search( r'ID=.*?\_(\d+)', x['attributes'])[1] ))

    return newGff, translation_str


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
            
    for prot in exonCheck:
        if not exonCheck[prot][0]:
            newEntry = copy.deepcopy(exonCheck[prot][2])
            newEntry['attributes'] = newEntry['attributes'].replace('.cds;', '.exon1;')
            newEntry['type'] = 'exon'
            gff.insert( exonCheck[prot][1] + 1, newEntry )

    return gff


def main( gff_path, fasta_path, prefix, fail = True ):

    if isinstance(gff_path, str):
        gff = gff2dict( gff_path )
    elif isinstance(gff_path, list):
        gff = gff_path
    if isinstance(fasta_path, str):
        assembly = fasta2dict( fasta_path )
    elif isinstance(fasta_path, dict):
        assembly = fasta_path

    if gff_path.endswith( 'gtf' ) or re.search(gtfComps()['id'], gff[0]['attributes']) is not None:
        exonGtf = intron2exon( gff )
        exonGtfCur = curCDS( exonGtf )
        exonGtfCurGenes, failed, flagged = addGenes( exonGtfCur, safe = fail )
        preGff = removeStartStop( exonGtfCurGenes )
        gffUncur, trans_str = curate( preGff, prefix, failed ) 
        gff = addExons( gffUncur )
    else:
        gff, trans_str = curate( 
            gff, prefix
            ) 
        failed, flagged = None, None
  
    fa = gff2proteome( gff, assembly )
 
    return gff, fa, trans_str, failed, flagged


if __name__ == '__main__':

    parser = argparse.ArgumentParser( description = 'Curates Funannotate or ' + \
        'post OrthoFiller output by naming accessions as <PREFIX>_####, where ' + \
        '#### is the is the index from ordered accessions. If you are planning' + \
        ' on using OrthoFiller, you may want to wait to not mix accessions.' )
    parser.add_argument( '-g', '--gff', required = True, help = '.gtf, .gff/.gff3' )
    parser.add_argument( '-f', '--fasta', required = True, help = 'Assembly fasta' )
    parser.add_argument( '-p', '--prefix', required = True,
        help = 'Accession prefix' )
    parser.add_argument( '-o', '--output', help = 'Output directory' )
    parser.add_argument( '--fail', default = True, action = 'store_false',
        help = 'Fail genes without start or stop codons.' )
    args = parser.parse_args()

    if args.output:
        output = formatPath(args.output, isdir = True)
        if not os.path.isdir(output):
            os.mkdir(output)
        output += args.prefix
    else:
        output = args.prefix

    gff, fa, trans_str, failed, flagged = main(
        formatPath(args.gff), formatPath(args.fasta), args.prefix, args.fail
        )

    with open( output + '.aa.fa', 'w' ) as out:
        out.write( dict2fasta( fa ) )
    with open( output + '.gff3', 'w' ) as out:
        out.write( dict2gff( gff ) + '\n' )
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
