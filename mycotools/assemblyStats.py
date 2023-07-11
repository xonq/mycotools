#! /usr/bin/env python3

# NEED to come up with more clever way of detecting mtdb

'''
Takes a database as argument 1 and an output for the .tsv as argument 2.
Calculates basic genome statistics.
'''

import os
import sys
import copy
import multiprocessing as mp
from mycotools.lib.dbtools import mtdb
from mycotools.lib.biotools import fa2dict
from mycotools.lib.kontools import format_path, eprint


def calcMask( contig_list ):

    seq = ''.join([x['sequence'] for x in contig_list])
    mask = seq.count('a')
    mask += seq.count('t')
    mask += seq.count('g')
    mask += seq.count('c')

    return mask


def sortContigs( assembly_path ):
    '''Imports fasta, creates a list of dicts for each contig length and its name.
       Sorts the list in descending order by length'''

    contigList = []
    assembly = fa2dict( assembly_path )
    for contig in assembly:
        contigList.append( {
            'len': len(assembly[contig]['sequence']), 
            'name': contig, 
            'sequence': str(assembly[contig]['sequence'])
        } )

    sortedList = sorted(contigList, key = lambda i: i['len'], reverse = True)

    return sortedList


def n50l50( sortedContigs ):
    '''Calculates contigs greater than 1000 bp, the total length excluding 1000 bp contigs.
       Then caculates l50, n50, largest contig, total contigs, shortest contig, and l50%.
       Returns a dictionary of the results.'''

    pass_fa = []
    total, total1000, bp1000, gc, gc1000, gctot, gctot1000 = 0, 0, 0, 0, 0, 0, 0
    for contig in sortedContigs:
        total += contig['len']
        gc += contig['sequence'].lower().count( "g" ) + contig['sequence'].lower().count( 'c' )
        gctot += contig['sequence'].lower().count( "g" ) + contig['sequence'].lower().count( 'c' ) + \
           contig['sequence'].lower().count( 'a' ) + contig['sequence'].lower().count( 't' ) 
        if contig['len'] >= 1000:
            gc1000 += contig['sequence'].lower().count( "g" ) + contig['sequence'].lower().count( 'c' )
            gctot1000 += contig['sequence'].lower().count( "g" ) + contig['sequence'].lower().count( 'c' ) + \
                contig['sequence'].lower().count( 'a' ) + contig['sequence'].lower().count( 't' ) 
            total1000 += contig['len']
            bp1000 += 1
            pass_fa.append( contig )


    out = {}
    count50, count, check = 0, 0, 0
    for contig in sortedContigs:
        count50 += contig['len']
        count += 1
        if count50 >= total1000/2:
            if contig['len'] >= 1000 and check != 1:
                out['n50-1000bp'] = int(contig['len'])
                out['l50-1000bp'] = int(count)
                out['l50%-1000bp'] = out['l50-1000bp']/int(bp1000)
                check = 1
            if count50 >= total/2:
                out['n50'] = int(contig['len'])
                out['l50'] = int(count)
                break

    try:
        out['l50%'] = out['l50']/len(sortedContigs)
        out['largest_contig'] = sortedContigs[0]['len']
        out['shortest_contig'] = sortedContigs[-1]['len']
        out['contigs'] = len(sortedContigs)
        out['contigs-1000bp'] = bp1000
        out['assembly_len'] = int(total)
        out['assembly_len-1000bp'] = int(total1000)
        out['gc'] = float(gc/gctot)
        out['gc-1000bp'] = float(gc1000/gctot1000)
        if 'n50-1000bp' not in out:
            out['n50-1000bp'] = 'na'
            out['l50-1000bp'] = 'na'
        maskCount = calcMask( sortedContigs )
        maskCount1000 = calcMask( pass_fa )
        out['mask%'] = maskCount / int( total ) * 100
        out['mask%-1000bp'] = maskCount1000 / int( total1000) * 100
        
    except KeyError:
        out = {}


    return out

def mngr(assembly_path, ome):
    sortedContigs = sortContigs(assembly_path)
    calcs = n50l50(sortedContigs)
    return ome, tuple([(x, calcs[x]) for x in calcs])


def main(in_path, log_path = None, cpus = 1, db = None):

    stats = {}

    if in_path.endswith('db') or db:
        head = '#ome\tn50-1000bp\tl50-1000bp\tl50%-1000bp\tn50\tl50\tl50%\tlargest_contig\tshortest_contig\tcontigs' + \
            '\tcontigs-1000bp\tassembly_len\tassembly_len-1000bp\tgc\tgc-1000bp\tmask%\tmask%-1000bp'

        prevOmes = {}
        if log_path:
            if not os.path.isfile( log_path ):
                with open( log_path, 'w' ) as log_open:
                    log_open.write( head )
            else:
                with open(log_path, 'r') as raw:
                    for line in raw:
                        if not line.startswith('#'):
                            omeI = line.index('\t')
                            ome = line[:omeI]
                            prevOmes[ome] = line[omeI+1:].rstrip()

        if not db:
            db = mtdb(in_path).set_index()

        cmds = []
        for ome in db:
            row = db[ome]
            if row['fna'] and ome not in prevOmes:
                cmds.append((row['fna'], ome,))
        with mp.Pool(processes=cpus) as pool:
            results = pool.starmap(mngr, cmds)

        calcs = {}
        for res in results:
            if res[1]:
                calcs[res[0]] = '\t'.join([str(x[1]) for x in res[1]])
            else:
                eprint('\t\tERROR:\t' + ome, flush = True)

        calcs = {
            k: v for k,v in \
            sorted(
                {**calcs, **prevOmes}.items(), 
                key = lambda x: x[0]
                )
            }

        if log_path:
            with open(log_path, 'w') as logWrite:
                logWrite.write(head + '\n')
                for ome in calcs:
                    data = calcs[ome]
                    logWrite.write(
                        ome + '\t' + data + '\n'
                        )
        else:
            print(head, flush = True)
            for ome in calcs:
                data = calcs[ome]
                print(
                    ome + '\t' + data, flush = True
                    )

    else:
        sortedContigs = sortContigs(in_path)
        calculations = n50l50( sortedContigs )
        if calculations:
            stats[ os.path.basename( os.path.abspath(in_path)) ] = n50l50( sortedContigs )
        else:
            eprint('\tERROR:\t' + in_path, flush = True)

        for stat in stats:
            if stats[stat]['shortest_contig'] >= 1000:
                stats[stat] = { 
                    info: stats[stat][info]  for info in stats[stat] if '1000bp' not in info 
                }
            for info in stats[stat]:
                print( '{:<25}'.format( info.upper() + ':' ) + str( stats[stat][info] ) , flush = True)
            

def cli():
    usage = '\nUSAGE: assembly statistics\nAssembly `fasta` or mycotoolsDB, optional output file if using database\n'
    if {'-h', '--help'}.intersection(set(sys.argv)):
        print(usage, flush = True)
        sys.exit(1)
    if len(sys.argv) < 2:
        print( usage , flush = True)
        sys.exit( 1 )

    in_path = format_path(sys.argv[1])
    if len(sys.argv) > 2:
        log_path = format_path(sys.argv[2])
    else:
        log_path = None

    main(in_path, log_path, os.cpu_count())

    sys.exit(0)



if __name__ == '__main__':
    cli()
