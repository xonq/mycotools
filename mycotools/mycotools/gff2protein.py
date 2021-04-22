#! /usr/bin/env python3

from mycotools.lib.dbtools import db2df, masterDB
from mycotools.lib.fastatools import fasta2dict, gff2dict, gff3Comps, dict2fasta
from mycotools.lib.kontools import formatPath, sysStart, eprint
from Bio.Seq import Seq
import re, sys


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

def main( gff_dicts, assem_dict ):
    
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

    usage = 'Inputs MycoDB compatible gff3, assembly (optional), and outputs proteins.' + \
       '\nUse an abstracted gene gff (acc2gff.py) to only output a ' + \
       'smaller set of gene(s).'
    args = sysStart( sys.argv, usage, 2, files = [sys.argv[1]] )

    gff_dicts = gff2dict(formatPath(args[1]))
    if len(sys.argv) > 2:
        assembly_dict = fasta2dict( formatPath(args[2]) )
    else:
        try:
            gene = re.search(gff3Comps()['Alias'], gff_dicts[0]['attributes'])[1]
            ome = re.search( r'(.*?)_', gene )[1]
        except:
            eprint( '\nERROR: ' + args[1] + ' is incompatible with MycoDB' , flush = True)
            sys.exit(1)
        db = db2df( masterDB() ).set_index('internal_ome')
        assembly_dict = fasta2dict(formatPath('$MYCOFNA/' + db.loc[ome]['assembly']))
        
    print(dict2fasta(main(gff_dicts, assembly_dict)), flush = True)
