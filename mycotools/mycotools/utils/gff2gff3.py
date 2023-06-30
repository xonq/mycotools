#! /usr/bin/env python3

# NEED to arrive at a consensus for protein and transcript IDs

import re
import sys
import copy
import argparse
from collections import defaultdict
from mycotools.lib.biotools import gff2list, list2gff, gff2Comps, gff3Comps
from mycotools.lib.kontools import format_path, eprint, vprint
from mycotools.utils.gtf2gff3 import add_genes, remove_start_stop
from mycotools.utils.curGFF3 import rename_and_organize

def gff2gff3(gff_list, ome, jgi_ome):

    comps2, exon_dict, cds_dict, out_list, gene_dict = gff2Comps(), {}, {}, [], {}
    for entry in gff_list:
        if entry['type'] == 'exon':
            name = re.search( comps2['id'], entry['attributes'] )[1]
  #          try:
            transcript = re.search( comps2['transcript'], entry['attributes'] )[1]
  
  #except TypeError:
 #               print(entry['attributes'], comps2['transcript'], flush = True)
#                sys.exit()
            if name not in exon_dict:
                exon_dict[name] = 0
            if name not in gene_dict:
                gene_dict[name] = { 
                    'protein': '', 'transcript': '', 'product': ''
                    }
            gene_dict[name]['transcript'] = transcript
            exon_dict[name] += 1
            exon_id = 'ID=exon_' + transcript + '_' + str(exon_dict[name])
            par_id = 'mRNA_' + transcript
            entry['attributes'] = exon_id + ';Parent=' + par_id + ';Alias=' + name
        elif entry['type'] == 'CDS':
            name = re.search( comps2['id'], entry['attributes'] )[1]
            protein = re.search( comps2['prot'], entry['attributes'] )[1]
            product = re.search( comps2['product'], entry['attributes'] )
            if name not in cds_dict:
                cds_dict[name] = 0
            if name not in gene_dict:
                gene_dict[name] = { 
                    'protein': '', 'transcript': '', 'product': ''
                    }
            if product is not None:
                gene_dict[name]['product'] = product[1]
            gene_dict[name]['protein'] = protein
            cds_dict[name] += 1
            cds_id = 'CDS_$_' + str( cds_dict[name] )
            entry['attributes'] = 'ID=' + cds_id + ';Alias=' + name

    comps3 = gff3Comps()
    for entry in gff_list:
        if entry['type'] not in {'start_codon', 'stop_codon'}:
            if entry['type'] == 'exon':
                name = re.search( comps3['Alias'], entry['attributes'] )[1]
                trans = gene_dict[name]['transcript']
                prot = ome + '_' + gene_dict[name]['protein']
                entry['attributes'] = re.sub(
                    comps3['Alias'], 
                    'Parent=mRNA_' + trans + ';Alias=' + prot, 
                    entry['attributes']
                    )
            elif entry['type'] == 'CDS':
                name = re.search( comps3['Alias'], entry['attributes'] )[1]
                trans = gene_dict[name]['transcript']
                prot = ome + '_' + gene_dict[name]['protein']
                entry['attributes'] = re.sub(
                    r'ID=CDS_\$_(\d+)', 
                    'ID=CDS_' + trans + r'_\1;' + 'Parent=mRNA_' + trans,
                    entry['attributes']
                    )
                entry['attributes'] = re.sub(
                    comps3['Alias'], 'Alias=' + prot, entry['attributes']
                    )
            elif entry['type'] == 'gene':
                name = re.search( comps2['id'], entry['attributes'] )[1]
                gene = 'gene_' + gene_dict[name]['transcript']
                prot = gene_dict[name]['protein']
                jgi = 'jgi.p|' + jgi_ome + '|' + prot
                prod = gene_dict[name]['product']
                trans = gene_dict[name]['transcript']
                alias = ome + '_' + prot
                entry['attributes'] = 'ID=' + gene + ';Name=' + jgi + \
                    ';portal_id=' + jgi_ome + ';product_name=' + prod + \
                    ';proteinId=' + prot + ';transcriptId=' + trans + \
                    ';Alias=' + alias
            elif 'RNA' in entry['type']:
                name = re.search( comps2['id'], entry['attributes'] )[1]
                gene = 'gene_' + gene_dict[name]['transcript']
                mrna = entry['type'] + '_' + gene_dict[name]['transcript']
                prot = gene_dict[name]['protein']
                jgi = 'jgi.p|' + jgi_ome + '|' + prot
                trans = gene_dict[name]['transcript']
                alias = ome + '_' + prot
                entry['attributes'] = 'ID=' + mrna + ';Name=' + jgi + \
                    ';Parent=' + gene + \
                    ';proteinId=' + prot + ';transcriptId=' + trans + \
                    ';Alias=' + alias
            out_list.append( entry )

    return out_list


def resolve_alternate_splicing(gff):
    # this is a hack job and should be done during add_genes
    contig2gene, a2z, a2gi = defaultdict(dict), defaultdict(list), {}
    for i, entry in enumerate(gff):
        alias = re.search(gff3Comps()['Alias'], entry['attributes'])[1]
        if entry['type'] == 'gene':
            a2gi[alias] = i
            contig = entry['seqid']
            start, end = entry['start'], entry['end']
            gene = re.search(gff3Comps()['id'], entry['attributes'])[1]
            contig2gene[contig][gene] = (start, end, alias)
        a2z[alias].append(entry)

    gene2gene = defaultdict(list)
    for contig, gene_dict in contig2gene.items():
        genes = list(gene_dict.keys())
        for i0, gene0 in enumerate(genes):
            s0, e0, a0 = gene_dict[gene0]
            for i1, gene1 in enumerate(genes[i0+1:]):
                s1, e1, a1 = gene_dict[gene1]
                if max(s0, s1) < min(e0, e1):
                    gene2gene[a0].append(a1)
            if not gene2gene[a0]:
                del gene2gene[a0]

    final = []
    while gene2gene:
        a0 = list(gene2gene.keys())[0]
        az = gene2gene[a0]
        all_hits = set(az)
        for a1 in az:
            for a2 in gene2gene[a1]:
                if a2 not in all_hits:
                    all_hits.add(a2)
                    az.append(a2)
        full = sorted({a0}.union(set(az)))
        final.append(full)
        for a in full:
            del gene2gene[a]

    todel = []
    for accs in final:
        if len(accs) == 1:
            continue
        a0 = accs[0]
        max_iz, min_iz = [], []
        for g_entry in a2z[a0]:
            if g_entry['type'] == 'gene':
                gid = re.search(gff3Comps()['id'], g_entry['attributes'])[1]
                g_entry['attributes'] += '|' + '|'.join(accs[1:])
                max_iz.append(max(g_entry['start'],g_entry['end']))
                min_iz.append(min(g_entry['start'], g_entry['end']))
                break
        for a1 in accs[1:]:
            todel.append(a2gi[a1])
            for entry in a2z[a1]:
                if 'RNA' in entry['type']:
                    entry['attributes'] = re.sub(gff3Comps()['par'], 
                                         f'Parent={gid}', entry['attributes'])
                    max_iz.append(max(entry['start'], entry['end']))
                    min_iz.append(min(entry['start'], entry['end']))
        g_entry['start'] = min(min_iz)
        g_entry['end'] = max(max_iz)
    for i in sorted(todel, reverse = True):
        del gff[i]
        
    return gff

def find_jgi_problems(gff3):
    # check for out of bounds CDS and exons
    rna2gene = {}
    gene_coords, other_coords, comps = {}, defaultdict(list), gff3Comps()
    for entry in gff3:
        if entry['type'] == 'gene':
            gene = re.search(comps['id'], entry['attributes'])[1]
            gene_coords[gene] = sorted([entry['start'], entry['end']])
        elif entry['type'] in {'CDS', 'exon'}:
            rna = re.search(comps['par'], entry['attributes'])[1]
            other_coords[rna].extend([entry['start'], entry['end']])
        elif 'RNA' in entry['type']:
            rna = re.search(comps['id'], entry['attributes'])[1]
            gene = re.search(comps['par'], entry['attributes'])[1]
            rna2gene[rna] = gene

    errors = {'ob':[], 'nr': []}
    for rna, coords in other_coords.items():
        if rna in rna2gene:
            gene = rna2gene[rna]
        else:
            errors['nr'].append(rna)
            continue
        max_other, min_other = max(coords), min(coords)
        max_gene, min_gene = max(gene_coords[gene]), min(gene_coords[gene])
        if max_other > max_gene or min_other < min_gene:
            errors['ob'].append(gene)
    return errors
            
    

def main(gff_list, ome, jgi_ome, safe = True, verbose = True):

    if gff_list[0]['attributes'].startswith('gene_id'):
        comps = gtfComps()
        gene_prefix = 'gene_id'
    else:
        comps = gff2Comps()
        gene_prefix = 'name'
    gff_prep, failed, flagged = add_genes(gff_list, safe = safe, 
                                          comps = comps,
                                          gene_prefix = gene_prefix)
    if failed:
        vprint( str(len(failed)) + '\tgenes failed', v = verbose , e = True, flush = True)
    if flagged:
        vprint( str(len(flagged)) + '\tgene coordinates from exons', v = verbose, e = True, flush = True)
    pregff3 = gff2gff3(gff_prep, ome, jgi_ome)
    gff3 = resolve_alternate_splicing(pregff3)
    gff3 = rename_and_organize(gff3)
    errors, err_name = find_jgi_problems(gff3), []
    for err, err_list in errors.items():
        if err_list:
            err_name.append(err.upper())
            if verbose:
                if err == 'ob':
                    eprint('ERROR: genes with out of bounds coordinates', flush = True)
                    eprint(",".join(err_list), flush = True)
                elif err == 'nr':
                    eprint('ERROR: missing RNA', flush = True)
                    eprint(",".join(err_list), flush = True)

    return gff3, errors

def cli():
    
    parser = argparse.ArgumentParser( description = 'Converts jgi gff2 to gff3' )
    parser.add_argument( '-i', '--input', required = True, help = 'JGI gff2' )
    parser.add_argument( '-o', '--ome', required = True, help = 'Internal ome' )
    parser.add_argument( '-j', '--jgi', required = True, help = 'JGI ome' )
    parser.add_argument( '--fail', default = True, action = 'store_false', \
        help = 'Fail genes w/o CDS sequences that lack start or stop codons' )
    args = parser.parse_args()

    gff_list = gff2list(format_path( args.input ))
    eprint( args.ome + '\t' + args.input , flush = True)
    gff3, errors = main( gff_list, args.ome, args.jgi, args.fail )
    print( list2gff( gff3 ) , flush = True)

    sys.exit( 0 )


if __name__ == '__main__':
    cli()
