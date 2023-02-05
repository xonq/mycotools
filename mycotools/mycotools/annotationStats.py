#! /usr/bin/env python3

import os
import re
import sys
from mycotools.lib.biotools import gff2list
from mycotools.lib.kontools import format_path, eprint


def compileExon( gff_path, output, ome = None ):
    '''
    Inputs: `gff_path` or mycotoolsDB file
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

    gff = gff2list( gff_path )
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
    try:
        if len( len_list ) % 2 == 0:
            median = len_list[int(len(len_list)/2 - 1)]
        else:
            median = (len_list[round(len(len_list)/2 - 1)] + len_list[round(len(len_list)/2 - 2)])/2
    except IndexError:
        eprint('ERROR: ' + os.path.basename(gff_path) + ' - no exons detected. Skipping', flush = True)
        return None

    if any( line for line in gff if line['type'] == 'intron' ):
        eprint('ERROR: ' + os.path.basename(gff_path) + ' - introns detected. Exons only considered', flush = True)

    if not output:
        print( '{:<25}'.format('GENE LENGTH:') + str(total) , flush = True)
        print( '{:<25}'.format('GENES:') + str(len(exon_dict)), flush = True)
        print( '{:<25}'.format('MEAN GENE LENGTH:') + str(total/len(exon_dict)), flush = True)
        print( '{:<25}'.format('MEDIAN GENE LENGTH:') + str(median), flush = True)
    geneStats = {
        'total_len': total, 'genes': len(exon_dict), 
        'mean_len': total/len(exon_dict), 'median_len': median 
        }

    return ome, geneStats


def main(in_path, log_path = None, cpus = 1):

    if in_path[-4:] not in {'.gtf', '.gff', 'gff3'}:
        from mycotools.lib.dbtools import mtdb
        import multiprocessing as mp
        db = mtdb(in_path).set_index()

        prevOmes = {}
        if log_path and os.path.isfile(log_path):
             with open(log_path, 'r') as raw:
                 for line in raw:
                     if not line.startswith('#'):
                         omeI = line.index('\t')
                         ome = line[:omeI]
                         prevOmes[ome] = line[omeI+1:].rstrip()

        cmds = []
        for ome in db:
            if ome not in prevOmes:
                cmds.append((db[ome]['gff3'], True, ome,))
        with mp.Pool(processes = cpus) as pool:
            res = pool.starmap(compileExon, cmds)

        outPrep = {}
        for result in res:
            if result:
                outPrep[result[0]] = '\t'.join([str(x) for x in result[1].values()])
        out = {
            k: v for k, v in \
            sorted({**outPrep, **prevOmes}.items(), key = lambda x: x[0])
            }


        if not log_path:
#            output_file = os.path.basename(format_path(sys.argv[1])) + '.annStats.tsv'
            print('#ome\ttotal_length\tgenes\tmean_length\tmedian_length', flush = True)
            for ome in out:
                print(ome + '\t' + out[ome], flush = True)
        else:
            with open(log_path, 'w') as write:
                write.write('#ome\ttotal_length\tgenes\tmean_length\tmedian_length\n')
                for ome in out:
                    write.write(ome + '\t' + out[ome] + '\n')
    else:
        ome, geneStats = compileExon(in_path, log_path)
        if log_path:
            out_str = 'total_length\tgenes\tmean_length\tmedian_length\n' + str(geneStats['total_len']) + \
                '\t' + str(geneStats['genes']) + '\t' \
                + str(geneStats['mean_len']) + '\t' + str(geneStats['median_len'])


if __name__ == '__main__':

    output = False
    usage = '\nUSAGE: `gff`/`gtf`/`gff3` OR mycotoolsDB, optional output file\n'
    if '-h ' in sys.argv or '--help' in sys.argv or '-h' == sys.argv[-1]:
        print(usage, flush = True)
        sys.exit(1)
    elif len( sys.argv ) < 2:
        print(usage, flush = True)
        sys.exit(1)
    elif not os.path.isfile( format_path(sys.argv[1]) ):
        print(usage, flush = True)
        sys.exit(1)
    elif len( sys.argv ) > 2:
        log_path = format_path(sys.argv[2])
    else:
        log_path = None

    in_path = format_path(sys.argv[1])
    main(in_path, log_path, os.cpu_count())

    sys.exit(0)
