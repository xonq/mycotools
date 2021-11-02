#! /usr/bin/env python3

from mycotools.lib.dbtools import db2df, masterDB
from mycotools.lib.biotools import fa2dict, gff2list, gff3Comps, dict2fa
from mycotools.lib.kontools import formatPath, sysStart, eprint
from Bio.Seq import Seq
import re, sys, argparse


def grabCDS( gff_dicts ):
    return [x for x in gff_dicts if x['type'].lower() == 'cds']

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

def ntNeg( rev_seq, negtig_dict ):

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

def ntmain( gff_dicts, assem_dict ):
    
    cdss = grabCDS( gff_dicts )
    pos_dict, neg_dict = grabCoords( cdss )

    genes_fa_dict = {}
    for contig in pos_dict:
        genes_fa_dict = { 
            **ntPos(assem_dict[contig]['sequence'], pos_dict[contig]), 
            **genes_fa_dict 
            }
    for contig in neg_dict:
        rev_comp = str(Seq(assem_dict[contig]['sequence']).reverse_complement())
        genes_fa_dict = {
            **ntNeg(rev_comp, neg_dict[contig]),
            **genes_fa_dict
            }

    return genes_fa_dict

def aamain( gff_dicts, assem_dict ):
    
    cdss = grabCDS( gff_dicts )
    pos_dict, neg_dict = grabCoords( cdss )

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

            
if __name__ == '__main__':

    parser = argparse.ArgumentParser(
       description = 'Inputs MycoDB compatible gff3, assembly (optional), ' + \
           'and outputs nucleotides/proteins. Use an abstracted gene gff ' + \
           '(acc2gff.py) to only output a smaller set of gene(s).'
       )
    parser.add_argument('-g', '--gff', required = True)
    parser.add_argument('-n', '--nucleotide', action = 'store_true')
    parser.add_argument('-p', '--protein', action = 'store_true')
    parser.add_argument('-a', '--assembly')
    args = parser.parse_args()

    input_gff = gff2list(formatPath(args.gff))
    if args.assembly:
        assembly_dicts = {'input': fa2dict( formatPath(args.assembly) )}
        gff_dicts = {'input': input_gff}
    else:
        db = db2df(masterDB()).set_index('internal_ome')
        gff_dicts, assembly_dicts = {}, {}
        try:
            for line in input_gff:
                gene = re.search(gff3Comps()['Alias'], line['attributes'])[1]
                ome = re.search( r'(.*?)_', gene )[1]
                if ome not in gff_dicts:
                    gff_dicts[ome] = []
                    assembly_dicts[ome] = fa2dict(
                        formatPath('$MYCOFNA/' + db.loc[ome]['assembly'])
                        )
                gff_dicts[ome].append(line)
        except IndexError:
            eprint('\nERROR: ' + args.gff + ' is incompatible with MycotoolsDB', flush = True)
            sys.exit(1)
        
    if args.protein:
        for ome in gff_dicts:
            print(dict2fa(aamain(gff_dicts[ome], assembly_dicts[ome])), flush = True)
    else:
        for ome in assembly_dicts:
            print(dict2fa(ntmain(gff_dicts[ome], assembly_dicts[ome])), flush = True)

    sys.exit(0)
