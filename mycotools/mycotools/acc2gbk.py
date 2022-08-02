#! /usr/bin/env python3

import os
import re
import sys
import argparse
import multiprocessing as mp
from collections import defaultdict
from Bio.Seq import Seq
from mycotools.lib.kontools import eprint, format_path, stdin2str
from mycotools.lib.dbtools import mtdb, masterDB
from mycotools.lib.biotools import fa2dict, gff2list, gff3Comps
from mycotools.acc2gff import dbMain as acc2gff

def col_CDS(gff_list, 
            types = {'gene', 'CDS', 'exon', 'mRNA', 
                     'tRNA', 'rRNA', 'RNA', 'pseudogene'}):
    '''Collect all CDS entries from a `gff` and store them into cds_dict. If the gene found in the CDS
    is in the set of `ids` then add that type to the CDS dict for that protein
    header.'''

    cds_dict = defaultdict(dict)
    for entry in gff_list:
        if entry['type'] in types:
            contig = entry['seqid']
            try:
                alias = re.search(gff3Comps()['Alias'], entry['attributes'])[1]
            except TypeError:
                eprint( '\n\tERROR: could not extract Alias ID from ' + gff , flush = True)
                continue
            aliases = alias.split('|') # to address alternate splicing in gene
            # aliases
            for alias in aliases:
                ome = alias[:alias.find('_')]
                if contig not in cds_dict[ome]:
                    cds_dict[ome][contig] = defaultdict(list)
                cds_dict[ome][contig][alias].append(entry)
            
    return cds_dict

def contig2gbk(ome, row, contig, contig_dict, 
               contig_seq, faa, product_search = r'product=([^;]+)'):

    seq_coords = []
    for prot, prot_list in contig_dict.items():
        for entry in prot_list:
            t_start, t_end = int(entry['start']), int(entry['end'])
            seq_coords.append(sorted([t_start, t_end]))

    seq_coords.sort(key = lambda x: x[0])
    init_startTest = int(seq_coords[0][0]) - 1
    init_endTest = int(seq_coords[-1][1])
    endTest = init_endTest - init_startTest
    startTest = 1

    # specify the contig edge - may cause parsers (BiG-SCAPE/antiSMASH) to fail
    edge = '/contig_edge=' + str(init_startTest) + '-;' \
         + str(len(contig_seq) - init_endTest) + '+'
    name = ome + '_' + contig
    relativeEnd = seq_coords[-1][1] - seq_coords[0][0]
    gbk = "LOCUS       " + name + '\n' + \
        "DEFINITION  " + name + '\n' + \
        "ACCESSION   " + contig + '\n' + \
        "VERSION     " + name + '\n' + \
        "KEYWORDS    .\nSOURCE    " + row['source'] + "\n  ORGANISM  " + \
        str(row['genus']) + '_' + str(row['species']) + '_' + \
        str(row['strain']) + '_' + ome + '\n' + \
        "COMMENT     coordinates are relative to adjustment: " + str(init_startTest) + \
        '\n            .\nFEATURES             Location/Qualifiers\n' + \
        '     region          ' + str(startTest) + '..' + str(endTest) + '\n' + \
        '                     ' + str(edge) + '\n'
#        '     region          1..' + str(relativeEnd) + '\n' + \
        # this should be
        # related to the start of the contig

#        '                     /product="' + product + '"\n'

    used_aliases = set() # to address alternate splicing and 1 entry per gene
    for prot, prot_list in contig_dict.items():
        # make products unique
#        prot_dict['products'] = set(prot_dict['products'])
 #       if len(prot_dict['products']) > 1:
  #          product = '|'.join(sorted(prot_dict['products']))
   #     else:
    #        product = list(prot_dict['products'])[0]
        final_coords = ''
        entries = defaultdict(list)
        for entry in prot_list:
            if entry['type'].lower() == 'cds':
                entries['CDS'].append(entry)
            elif entry['type'].lower() == 'exon':
                entries['exon'].append(entry)
            elif entry['type'].lower() in {'mrna', 'rrna', 'trna', 'rna'}:
                entries['RNA'].append(entry)
            elif entry['type'].lower() in {'gene', 'pseudogene'}:
                entries['gene'].append(entry)

        if entries['CDS']:
            entry = entries['CDS'][0]
            strand = entry['strand']
            cds_coords = ''
            if strand == '+':
                for cds in entries['CDS']:
                    coords = [cds['start'], cds['end']]
                    max_c, min_c = max(coords), min(coords)
                    cds_coords += str(min_c) + '..' + str(max_c) + ','
            else:
                for cds in entries['CDS']:
                    coords = [cds['start'] - init_startTest, 
                              cds['end'] - init_startTest]
                    max_c, min_c = max(coords), min(coords)
                    cds_coords += 'complement(' + str(min_c) + '..' \
                               + str(max_c) + '),'
            cds_coords = cds_coords[:-1]
            cds_coords_list = cds_coords.split(',')
            final_coords = '' 
            total_len = 0
            for coord in cds_coords_list[:-1]:
                if total_len + len(coord) < 57 or total_len == 0:
                    final_coords += coord + ','
                    total_len += len(coord) + 1
                else:
                    total_len = len(coord) + 1
                    final_coords += '\n                     ' \
                                 + coord + ','
            if total_len == 0 or total_len + len(cds_coords_list[-1]) < 58:
                final_coords += cds_coords_list[-1]
            else:
                final_coords += '\n                     ' \
                             + cds_coords_list[-1]
            if len(cds_coords_list) > 1:
                final_coords = 'join(' + final_coords + ')'
            alias = re.search(gff3Comps()['Alias'], entry['attributes'])[1]
            id_ = re.search(gff3Comps()['id'], entry['attributes'])[1]

            try:
                product = re.search(product_search, entry['attributes'])[1]
            except TypeError: # no product
                product = ''
            gbk += '     ' + entry['type'] + \
                '                '[:-len(entry['type'])]
            final_coords += '\n                     '
            gbk += final_coords
            gbk += '/locus_tag="' + alias + '"\n                     ' 

            if product:
                gbk += '/product="' + product + '"\n                     '
            if entry['phase'] != '.':
                gbk += '/phase="' + entry['phase'] + '"\n                     '
            gbk += '/source="' + entry['source'] + '"\n'
            gbk += '                     '
            gbk += '/protein_id="' + alias + '"\n'
            gbk += '                     ' + \
                '/transl_table=1\n                     /translation="'
            inputSeq = faa[alias]['sequence']
            lines = [inputSeq[:45]]
            if len(inputSeq) > 45:
                restSeq = inputSeq[45:]
                n = 59
                toAdd = [restSeq[i:i+n] for i in range(0, len(restSeq), n)]
                lines.extend(toAdd)
            gbk += lines[0] + '\n'
            if len(lines) > 1:
                for line in lines[1:len(lines)-1]:
                    gbk += '                     ' + line + '\n'
                gbk += '                     ' + lines[-1] + '"\n'
            else:
                gbk += '"\n'

        for entry in entries['gene']:
            alias = re.search(gff3Comps()['Alias'], entry['attributes'])[1]
            if '|' in alias: # alternately spliced gene
                if alias in used_aliases:
                    continue
                else:
                    used_aliases.add(alias)
            id_ = re.search(gff3Comps()['id'], entry['attributes'])[1]

            try:
                product = re.search(product_search, entry['attributes'])[1]
            except TypeError: # no product
                product = ''
            start, end = sorted([int(entry['start']), int(entry['end'])])
            gbk += '     ' + entry['type'] + \
                '                '[:-len(entry['type'])]
            if entry['strand'] == '+':
                gene_coords = str(start-init_startTest) + '..' \
                    + str(end-init_startTest) + '\n                     '
            else:
                gene_coords = 'complement(' + str(start-init_startTest) \
                    + '..' + str(end-init_startTest) + ')\n                     '
            gbk += gene_coords
            gbk += '/locus_tag="' + alias + '"\n                     ' 

            if product:
                gbk += '/product="' + product + '"\n                     '
            if entry['phase'] != '.':
                gbk += '/phase="' + entry['phase'] + '"\n                     '
            gbk += '/source="' + entry['source'] + '"\n'

        for entry in entries['RNA']:
            gbk += '     ' + entry['type'] + \
                '                '[:-len(entry['type'])]
            if final_coords:
                gbk += final_coords
            else:
                gbk += gene_coords
            gbk += '/locus_tag="' + alias + '"\n                     ' 

            if product:
                gbk += '/product="' + product + '"\n                     '
            if entry['phase'] != '.':
                gbk += '/phase="' + entry['phase'] + '"\n                     '
            gbk += '/source="' + entry['source'] + '"\n'
            gbk += '                     /transcript_id="' + id_ + '"\n'


    total_coords = []
    for coord in seq_coords:
        total_coords.extend(coord)
    min_c, max_c = min(total_coords), max(total_coords)
#    for coordinates in seq_coords:
    assSeq = contig_seq[min_c:max_c]
    seqLines = [assSeq[i:i+60] for i in range(0, len(assSeq), 60)]
    count = -59

    gbk += 'ORIGIN\n'
    for line in seqLines:
        count += 60
        gbk += '{:>9}'.format(str(count)) + ' '
        seps = [line[i:i+10] for i in range(0, len(line), 10)]
        for sep in seps:
            gbk += sep.lower() + ' '
        gbk = gbk.rstrip() + '\n'
    
    gbk += '//'
    return gbk

def genGBK(ome_dict, row, ome, product_search):

    assembly = fa2dict(row['fna'])
    faa = fa2dict(row['faa'])
    gbk = {}
    for contig, contig_dict in ome_dict.items():
        gbk[contig] = contig2gbk(ome, row, contig, contig_dict, 
                   assembly[contig]['sequence'], faa, 
                   product_search = product_search)

    return gbk
        

def main(gff_list, db, product_search = r'product=([^;]+)'):
    cds_dict = col_CDS(gff_list, 
            types = {'gene', 'CDS', 'exon', 'mRNA', 
                     'tRNA', 'rRNA', 'RNA', 'pseudogene'})
    ome_gbks = {}
    for ome, ome_dict in cds_dict.items():
        ome_gbks[ome] = genGBK(ome_dict, db[ome], ome, product_search)

    return ome_gbks

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description = "Inputs MTDB gff or accessions, outputs GenBank file" )
    parser.add_argument('-a', '--accession', help = '"-" for stdin')
    parser.add_argument('-i', '--input', help = 'File with accessions')
    parser.add_argument('-g', '--gff', help = 'GFF3 input')
    parser.add_argument('-p', '--product', 
        help = 'Product regular expression: DEFAULT: "Product=([^;]+)"')
    parser.add_argument('-o', '--ome', action = 'store_true', 
        help = 'Output files by ome code')
    parser.add_argument('-s', '--seqid', action = 'store_true',
        help = 'Output files by sequence ID (contig/scaffold/chromosome)'
        )
    parser.add_argument('-d', '--database',
        help = "MTDB; DEFAULT: master", default = masterDB())
    parser.add_argument('-c', '--cpu', type = int, default = 1)
    args = parser.parse_args()

    if args.input:
        input_file = format_path(args.input)
        with open(input_file, 'r') as raw:
            accs_prep = [x.rstrip().split() for x in raw]
        accs = []
        for acc in accs_prep:
            accs.extend(acc)
    elif args.accession:
        if args.accession == '-':
            data = stdin2str()
            accs = data.split()
        else: 
            accs = [args.accession]
    elif not args.gff:
        raise FileNotFoundError('need input file or accession')

    if not args.product:
        regex = r'Product=([^;])+'
    else:
        regex = r''
        if args.regex.startswith(("'",'"')):
            args.regex = args.regex[1:]
        if args.regex.endswith(('"',"'")):
            args.regex = args.regex[:-1]
        for char in args.regex:
            regex += char

    # import database and set index
    db = mtdb(format_path(args.database))
    db = db.set_index() 

    gbk_dict = {}
    if args.gff:
        gbk_dict = main(gff2list(format_path(args.gff)), db, regex)
    else:
        gffs = acc2gff(db, accs, args.cpu)
        for ome, gff in gffs.items():
            gbk_dict[ome] = main(gff, db, regex)[ome]
    
    if args.ome:
        for ome, gbks in gbk_dict.items():
            with open(ome + '.acc2gbk.gbk', 'w') as out:
                out.write('\n'.join([gbk for gbk in gbks]))
    elif args.seqid:
        for ome, gbks in gbk_dict.items():
            for contig, gbk in gbks.items():
                with open(ome + '_' + contig + '.gbk', 'w') as out:
                    out.write(gbk)
    else:
        for ome, gbks in gbk_dict.items():
            for contig, gbk in gbks.items():
                print(gbk, flush = True)

    sys.exit(0)
