#! /usr/bin/env python3

from mycotools.lib.fastatools import gff2dict
from mycotools.lib.kontools import formatPath
import sys, re, os


def compileExon( gff_path, output ):
    '''
    Inputs: `gff_path` and `output` (converted to boolean)
    Outputs: summary annotation statistics
    Import the gff, find the first exon, and test which protein regular
    expression pattern to use. For each line in the gff, if it is an exon grab
    the protein and include it in `exon_dict`. Extend the start and stop
    coordinates for that exon. For each protein in the `exon_dict`, sort the
    entry, add the exon length to the total by subtracting the smallest entry
    from the largest for that protein. Append this value to the length list for
    median evaluation. Sort the `len_list` and calculate/output annotation
    statistics.
    '''

    gff = gff2dict( gff_path )
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
            exon_dict[ prot ].extend( [ int( line['start'] ), int( line['end'] ) ] )

    total, len_list = 0, []
    for prot in exon_dict:
        exon_dict[ prot ].sort()
        total += exon_dict[prot][-1] - exon_dict[prot][0]
        len_list.append( exon_dict[prot][-1] - exon_dict[prot][0] )

    len_list.sort()
    if len( len_list ) % 2 == 0:
        median = len_list[int(len(len_list)/2 - 1)]
    else:
        median = (len_list[round(len(len_list)/2 - 1)] + len_list[round(len(len_list)/2 - 2)])/2


    if not output:
        print( '{:<25}'.format('GENE LENGTH:') + str(total) , flush = True)
        print( '{:<25}'.format('GENES:') + str(len(exon_dict)), flush = True)
        print( '{:<25}'.format('MEAN GENE LENGTH:') + str(total/len(exon_dict)), flush = True)
        print( '{:<25}'.format('MEDIAN GENE LENGTH:') + str(median), flush = True)
        geneStats = None
    else:
        geneStats = [ total, len(exon_dict), total/len(exon_dict), median ]

    if any( line for line in gff if line['type'] == 'intron' ):
        print('WARNING: OrthoFiller outputs do not have exons and will be skipped.\n' + \
            'Curate OrthoFiller with `curAnnotation.py`')


    return exon_dict, geneStats


if __name__ == '__main__':

    output = False
    usage = '\nUSAGE: `gff`/`gtf`/`gff3`, optional output file\n'
    if '-h ' in sys.argv or '--help' in sys.argv or '-h' == sys.argv[-1]:
        print(usage, flush = True)
        sys.exit(1)
    elif len( sys.argv ) < 2:
        print(usage, flush = True)
        sys.exit(1)
    elif not os.path.isfile( formatPath(sys.argv[1]) ):
        print(usage, flush = True)
        sys.exit(1)
    elif len( sys.argv ) > 3:
        output = True       

    exon_dict, geneStats = compileExon( formatPath(sys.argv[1]), output )
    if output:
        out_str = 'total_len\tgenes\tmean_len\tmedian_len\n' + str(geneStats[0]) + '\t' + str(geneStats[1]) + '\t' \
            + str(geneStats[2]) + '\t' + str(geneStats[3])

    sys.exit(0)    
