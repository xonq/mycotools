#! /usr/bin/env python3

# NEED to convert gff list to appropriate types

import re
import sys
from mycotools.lib.kontools import eprint

aa_weights = {
    'A': 89.1, 'R': 174.2, 'N': 132.1, 'D': 133.1, 'C': 121.2, 'E': 147.1,
    'Q': 146.2, 'G': 75.1, 'H': 155.2, 'I': 131.2, 'L': 131.2, 'K': 146.2,
    'M': 149.2, 'F': 165.2, 'P': 115.1, 'S': 105.1, 'T': 119.1, 'W': 204.2,
    'Y': 181.2, 'V': 117.1
    }
aa_weights['X'] = sum(aa_weights.values())/len(aa_weights)

# codon translation table
translation_table = {
    'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
    'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R', 'AGA':'R', 'AGG':'R',
    'AAT':'N', 'AAC':'N',
    'GAT':'D', 'GAC':'D',
    'TGT':'C', 'TGC':'C',
    'CAA':'Q', 'CAG':'Q',
    'GAA':'E', 'GAG':'E',
    'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G',
    'CAT':'H', 'CAC':'H',
    'ATT':'I', 'ATC':'I', 'ATA':'I',
    'TTA':'L', 'TTG':'L', 'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
    'AAA':'K', 'AAG':'K',
    'ATG':'M',
    'TTT':'F', 'TTC':'F',
    'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
    'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S', 'AGT':'S', 'AGC':'S',
    'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
    'TGG':'W',
    'TAT':'Y', 'TAC':'Y',
    'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',
    'TAA':'*', 'TGA':'*', 'TAG':'*'
}

rev_comp_table = {
    'T': 'A', 'G': 'C', 'C': 'G', 'A': 'T',
    'N': 'N', 'Y': 'R', 'R': 'Y', 'K': 'M',
    'M': 'K', 'B': 'V', 'D': 'H', 'H': 'D',
    'V': 'B', 't': 'a', 'g': 'c', 'c': 'g',
    'a': 't', 'n': 'n', 'y': 'r', 'r': 'y',
    'k': 'm', 'm': 'k', 'b': 'v', 'd': 'h',
    'h': 'd', 'v': 'b', 'W': 'W', 'S': 'S',
    'w': 'w', 's': 's'
    }


def calc_weight(seq):
    seq = seq.replace('*','')
    weight = aa_weights[seq[0]]
    for char in seq[1:]:
        weight += aa_weights[char] - 18.01
    return weight



def reverse_complement(seq):

    new_seq = ''
    for i in seq[::-1]:
        new_seq += rev_comp_table[i]

    return new_seq


def fa2dict(fasta_input, file_ = True):
    fasta_dict = {}
    if file_:
        with open(fasta_input, 'r') as fasta:
            str_fasta = ''
            for line in fasta:
                data = line.rstrip()
                if data.startswith('>'):
                    str_fasta += '\n' + data.rstrip() + '\n'
                elif data:
                    str_fasta += data
                # concatenates the prep string with the prepped line from the fasta file
        str_fasta = str_fasta.lstrip()
    else:
        str_fasta = fasta_input

    # extracts 0) seq ID and description and 1) sequence
    extracted = re.findall( r'(^>[^\n]*)\n([^>]*)', str_fasta, re.M)

    # adds a new dictionary for each gene
    for val in extracted:
        step1 = val[0]
        step2 = re.search(r'^>([^ ]*)', step1)
        gene = step2[1]
        step3 = re.search(r' (.*)', step1)
        if step3:
            descrip = step3[1]
        else:
            descrip = ''

        seq = val[1]
        if seq:
            if seq[-1] == '\n' or seq[-1] == '\r':
                seq = seq.rstrip()
        else:
            continue

        # prepares dictionaries for each gene with description, seq, rvcmpl_seq, and codons
        fasta_dict[gene] = {}
        if descrip != '\n':
            fasta_dict[gene]['description'] = descrip
        fasta_dict[gene]['sequence'] = seq

    return fasta_dict


# truncates sequences based on inputted lenght
def dnatrunc(fasta_dict,trunc_length):
    for gene in fasta_dict:
        if len(fasta_dict[gene]['sequence']) >= trunc_length:
            fasta_dict[gene]['sequence'] = fasta_dict[gene]['sequence'][:(trunc_length - 1)]
            # I should be able to do this in one line
            fasta_dict[gene]['reverse_complement'] = fasta_dict[gene]['reverse_complement'][-1:(0-trunc_length):-1]
            fasta_dict[gene]['reverse_complement'] = fasta_dict[gene]['reverse_complement'][::-1]
    return(fasta_dict)

# takes output dict from fasta2Dict and extracts codons from all 6 reading frames
def dict2codon(fasta_dict):

    fasta_dict[gene]['codons'] = {}

    # outputs a dictionary of reading frames and the possible codons for each
    for index in range(3):
        fasta_dict[gene]['codons']['reading_frame_' + str(index)] = []
        codon = ''
        for nt in fasta_dict[gene]['sequence'][index:]:
            if len(codon) < 3:
                codon = codon + nt
            if len(codon) == 3:
                fasta_dict[gene]['codons']['reading_frame_' + str(index)].append(codon)
                codon = ''

    # let's extract codons from all reading frames for the reverse sequences too
    for index in range(3):
        fasta_dict[gene]['codons']['reverse_reading_frame_' + str(index)] = []
        codon = ''
        for nt in fasta_dict[gene]['reverse_complement'][index:]:
            if len(codon) < 3:
                codon = codon + nt
            if len(codon) == 3:
                fasta_dict[gene]['codons']['reverse_reading_frame_' + str(index)].append(codon)
                codon = ''

    return codondict


def dict2fa(fasta_dict, description = True):

    fasta_string = ''
    if description:
        for gene in fasta_dict:
            fasta_string += '>' + gene.rstrip()
            if 'description' in fasta_dict[gene]:
                fasta_string += (' ' + fasta_dict[gene]['description']).rstrip() + '\n'
            else:
                fasta_string += '\n'
            fasta_string += fasta_dict[gene]['sequence'].rstrip() + '\n'

    else:
        for gene in fasta_dict:
            fasta_string += '>' + gene.rstrip() + '\n' + \
            fasta_dict[ gene ][ 'sequence' ].rstrip() + '\n'

    return fasta_string

# need to turn these into classes and class functions
def calc_gc(gene):

    G = gene['sequence'].count('G')
    C = gene['sequence'].count('C')
    GC = (G + C)/len(gene['sequence'])
    GC_con = '{:.2%}'.format(GC)

    return(GC_con)

# need to change into a class
def gff2list(gff_info, path = True, error = True):

    gff_list_dict = []
    if path:
        with open( gff_info, 'r' ) as raw_gff:
            data = [x.split('\t') for x in raw_gff.read().split('\n') \
                    if x and not x.startswith('#')]
    else:
        data = [x.split('\t') for x in gff_info.split('\n') \
                if x and not x.startswith('#')]
    try:
        for col_list in data:
            gff_list_dict.append({
                'seqid': col_list[0], 'source': col_list[1], 'type': col_list[2],
                'start': int(col_list[3]), 'end': int(col_list[4]), 'score': col_list[5],
                'strand': col_list[6], 'phase': col_list[7],
                'attributes': col_list[8]
                })
    except IndexError:
        raise IndexError(str(len(col_list)) + '/9 expected tab-' \
                        + 'delimitted fields: ' + str(col_list))
    except ValueError:
        raise ValueError(str(col_list[4:6]) + ' invalid integer ' \
                        + 'conversion: ' + str(col_list))
            
    return gff_list_dict


def list2gff(gff_list, ver = 3):

    if ver:
        gff_str = '##gff-version ' + str(ver) + '\n'
    else:
        gff_str = ''
    for line in gff_list:
        add_str = '\t'.join(str(val) for val in list(line.values()))
        gff_str += add_str + '\n'

    return gff_str.rstrip()

def gff3Comps( source = None ):

    comps = {}
    comps['par'] = '(?:^|(?<=;))' + r'Parent=["\']?([^;\'"]+)'
    comps['id'] = '(?:^|(?<=;))' + r'ID=["\']?([^;\'"]+)'
    comps['Alias'] = '(?:^|(?<=;))' + r'Alias=["\']?([^;\'"]+)'
    comps['product'] = '(?:^|(?<=;))' + r"""product=["\']?([^;"']+)"""
    comps['OG'] = '(?:^|(?<=;))' + r'OG=["\']?([\w+:\d+\|]+)'
    comps['ver'] = 'gff3'
    comps['prot'] = '(?:^|(?<=;))' + r'protein_id=["\']?([^;\'"]+)|' \
                  + '(?:^|(?<=;))' + r'proteinId=["\']?([^"\';]+)'
    comps['transcript'] = '(?:^|(?<=;))' + r'transcriptId=["\']?([^;"\'])|' \
                  + '(?:^|(?<=;))' + r'transcript_id=["\']?([^;\'"]+)'
   
    return comps

def gff2Comps():

    comps = {}
    comps['id'] = r'name "([^"]+)"'
    # this unfortunately includes "gene_name" from gtf naming
    comps['prot'] = r'proteinId ([^;]+)'
    comps['transcript'] = r'transcriptId ([^;]+)'
    comps['alias'] = r'alias "([^"]+)"'
    comps['product'] = r'product_name "([^"]+)'
    comps['ver'] = 'gff2'

    return comps

def gtfComps():

    comps = {}
    comps['id'] = r'gene_id "?([^"]+)"?'
    comps['transcript'] = r'transcript_id "([^"]+)"'
    comps['alias'] = r'alias "([^"]+)"'
    comps['ver'] = 'gtf'

    return comps

def compileExon( gff ):
	
    exon_dict = {}

    for index in range(len(gff)):
        if gff[index]['type'].lower() == 'exon':
            break

    protComp = re.compile( r';Parent\=([^;]*)' )
    if not protComp.search( gff[index]['attributes'] ):
        protComp = re.compile( r'gene_id "(.*?)"' )
        if not protComp.search( gff[index]['attributes'] ):
            protComp = re.compile( r'name "(.*?)"\;' )
            if not protComp.search( gff[index]['attributes'] ):
                protComp = re.compile( r'ID=(.*?);' )

    for line in gff:
        if line['type'].lower() == 'exon':
            prot = protComp.search( line['attributes'] )[1]
            if prot not in exon_dict:
                exon_dict[ prot ] = []
            if line['strand'] == '+':
                exon_dict[prot].append( [int(line['start'] ) - 1, int( line['end'] )] )
            else:
                if int(line['start']) > int(line['end']):
                    exon_dict[prot].append( [int (line['end'] ) - 1, int( line['start'])])
                else:
                    exon_dict[prot].append( [int (line['start'] ) - 1, int( line['end'])])

    for prot in exon_dict:
        exon_dict[ prot ] = sorted(exon_dict[ prot ], key = lambda i: i[0])

    return exon_dict
