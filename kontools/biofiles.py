#! /usr/bin/env python3
# parses all fasta files inputted as arguments and stores them as {filename: {genes: {description: }{sequence: }}}

import re
import sys
import os

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

class fa:

    def fa2dict( self, data, datatype = 'path'):
        fasta_dict = {}

        if datatype == 'path':
            with open( os.path.normpath( os.path.expanduser( str( data ))), 'r') as fasta:
                str_fasta = ''
                for line in fasta:
                    # rstrips if it's not a sequence line and isn't a line with only \n
                    if line[0] != '>' and line[0] != '\n':
                        line = line.rstrip()
                    # concatenates a new line character if it is the first sequence line
                    if line[0] == '>' and str_fasta != '':
                        str_fasta = str_fasta + '\n'
                    # concatenates the prep string with the prepped line from the fasta file
                    str_fasta = str_fasta + line

        elif datatype == str:
            str_fasta = data

        # extracts 0) seq ID 1) description (or the \n if there is non) 2) sequence
        extracted = re.findall(r'(>\S+)([^\n]\s\w+.+\n|\n)(\S+[^>])', str_fasta)

        # adds a new dictionary for each gene
        for index in range(len(extracted)):
            gene = extracted[index][0]
            gene = gene[1:]

            # remove description if it is not actually present
            descrip = extracted[index][1]
            descrip = descrip[1:]
            seq = extracted[index][2]
            if seq[-1] is '\n' or '\r':
                seq = seq.rstrip()

            # prepares dictionaries for each gene with description, seq, rvcmpl_seq, and codons
            fasta_dict[gene] = {}
            if descrip != '\n':
                fasta_dict[gene]['description'] = descrip
            fasta_dict[gene]['sequence'] = seq.upper()

        self = fasta_dict

        return self

    def __init__( self, data, datatype ):
        if datatype == 'path':
            self = self.fa2dict( data )
        elif datatype == str or datatype == 'str':
            self = self.fa2dict( data, str )
        else:
            raise ValueError('datatype is not specified as "path" or str')


    def tostr(fasta_dict):

        fasta_string = ''
        for gene in fasta_dict:
            fasta_string = fasta_string + '>' + gene
            if fasta_dict[gene].get('description'):
                fasta_string = fasta_string + '' + fasta_dict[gene]['description']
            else:
                fasta_string = fasta_string + '\n'
            fasta_string = fasta_string + fasta_dict[gene]['sequence'] + '\n'

        return fasta_string

##NEED a reverse complement function and then integrate dict2codon



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

# need to turn these into classes and class functions
def calc_gc(gene):

    G = gene['sequence'].count('G')
    C = gene['sequence'].count('C')
    GC = (G + C)/len(gene['sequence'])
    GC_con = '{:.2%}'.format(GC)

    return(GC_con)



# need to change into a class
def gff2dict( gff_path ):

    gff_list_dict = []
    with open( gff_path, 'r' ) as raw_gff:
        for line in raw_gff:
            if not line.startswith('##'):
                if line.endswith('\n'):
                    line = line.rstrip()
                col_list = line.split(sep='\t')
                gff_list_dict.append({
                    'seqid': col_list[0], 'source': col_list[1], 'type': col_list[2],
                    'start': col_list[3], 'end': col_list[4], 'score': col_list[5],
                    'strand': col_list[6], 'phase': col_list[7], 'attributes': col_list[8]
                    })
                
    return gff_list_dict
