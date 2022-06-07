#! /usr/bin/env python3
# NEED to make outgroup detection within the bounds of maximum sequences, optional to do so
# NEED to refine further
# extract a fasta of the cluster of interest
# fix output

from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform
#from sklearn.cluster import AgglomerativeClustering
from mycotools.lib.kontools import multisub, findExecs, formatPath, eprint, vprint, readJson, writeJson, mkOutput
from mycotools.lib.biotools import fa2dict, dict2fa
import string, argparse, os, sys, itertools, tempfile, re, random, copy, pandas as pd, subprocess

sys.setrecursionlimit(1000000)

def splitFasta( fa_path, output ):

    fastas = []
    fa = fa2dict( fa_path )
    for header in fa:
        with open( output + '/' + header + '.fa', 'w' ) as out:
            out.write( '>' + header + '\n' + fa[header]['sequence'] + '\n' )
        fastas.append( output + '/' + header + '.fa' )

    return fastas 

 
def makeDmndDB(diamond, queryFile, output_dir, cpus = 1):
    outputDB = output_dir + re.sub(r'\.[^\.]+$', '', os.path.basename(queryFile))
    cmd = [
        diamond, 'makedb', '--in', queryFile,
        '--db', outputDB, '--threads', str(cpus)
        ]
    dmndDBcode = subprocess.call(
        cmd, stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL
        )
    return outputDB, dmndDBcode


def runDmnd(
    diamond, queryFile, queryDB, outputFile, pid = True,
    cpus = 1, verbose = False, blast = 'blastp'
    ):

    if pid:
        distanceVal = 'pident'
    else:
        distanceVal = 'bitscore'

    cmd = [
        diamond, blast, '-d', queryDB, '-q', queryFile,
        '--out', outputFile, '--threads', str(cpus),
        '--outfmt', '6', 'qseqid', 'sseqid', distanceVal,
        ]

    if not verbose:
        dmndCode = subprocess.call(
            cmd, stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL
            )
    else:
        dmndCode = subprocess.call(
            cmd
            )
    return outputFile, dmndCode


def readDmndDist( outputFile, minVal, pid = True ):
    # should convert this to numpy and sparse

    dist_dict = {}
    with open(outputFile, 'r') as raw:
        for line in raw:
            q, s, v0 = line.rstrip().split('\t')
            v = round(float(v0))
            if q == s:
                continue
            if v > minVal:
  #              if q not in out_dict:
 #                   out_dict[q] = {}
                if q not in dist_dict:
                    dist_dict[q] = {}
                    dist_dict[q][q] = 0.0
                if s not in dist_dict:
                    dist_dict[s] = {}
                    dist_dict[s][s] = 0.0
                dist_dict[q][s] = 100 - v
                dist_dict[s][q] = 100 - v
#                out_dict[q][s] = v

    distanceMatrix = pd.DataFrame(dist_dict).sort_index(1).sort_index(0).fillna(100)
#    outMatrix = pd.DataFrame(out_dict).sort_index(1).sort_index(0).fillna(1)

    if pid:
        distanceMatrix /= 100
 #       outMatrix /= 100
    else:
        minV, maxV = min(distanceMatrix), max(distanceMatrix)
        denom = maxV - minV
        distanceMatrix -= minV
#        outMatrix -= minV
        distanceMatrix /= denom
 #       outMatrix /= denom
        distanceMatrix = 1 - distanceMatrix
  #      outMatrix = 1 - outMatrix

    return distanceMatrix #, outMatrix


def writeDist(outputFile, outMatrix):
    dist_out = ''
    for i, row in outMatrix.iterrows():
        for column in outMatrix.columns:
            if str(column) != str(i) and outMatrix.at[i, column] < 1:
                dist_out += str(i) + '\t' + str(column) + '\t' + \
                    str(outMatrix.at[i, column]) + '\n'
    with open( output, 'w' ) as out:
        out.write( dist_out )


def runUsearch( fasta, output, maxdist, cpus = 1, verbose = False ):

    if verbose:
        subprocess.call([
            'usearch', '-calc_distmx', fasta,
            '-tabbedout', output, '-maxdist',
            maxdist, '-threads', str(cpus)
            ] #stdout = subprocess.PIPE,
            #stderr = subprocess.PIPE
            )
    else:
        subprocess.call([
            'usearch', '-calc_distmx', fasta,
            '-tabbedout', output, '-maxdist',
            maxdist, '-threads', str(cpus)
            ], stdout = subprocess.DEVNULL,
            stderr = subprocess.DEVNULL
            )



def importDist( dis_path, sep = '\t' ):
    '''Imports a distance matrix with each line formatted as `organism $SEP organism $SEP distance`.
    - this is equivalent to the `-tabbedout` argument in `usearch -calc_distmx`. The function 
    compiles a dictionary with the information and reciprocal information for each organism in each
    line, then converts this dictionary of dictionaries into a pandas dataframe. As a distance matrix,
    NA values are converted to maximum distance (1).'''

    distanceMatrix = pd.DataFrame()
    with open( dis_path, 'r' ) as raw:
        data = raw.read()
    prepData = [ x.rstrip() for x in data.split('\n') if x != '' ]

    dist_dict = {}
    for line in prepData:
        vals = line.split( sep )
        if vals[0] not in dist_dict:
            dist_dict[ vals[0] ] = {}
        if vals[1] not in dist_dict:
            dist_dict[ vals[1] ] = {}
        dist_dict[ vals[0] ][ vals[1] ] = float(vals[2])
        dist_dict[ vals[1] ][ vals[0] ] = float(vals[2])

    distanceMatrix = pd.DataFrame( dist_dict ).sort_index(0).sort_index(1)

    return distanceMatrix.fillna( 1.0 )


def scikitaggd( distanceMatrix, maxDist = 0.6, linkage = 'single' ):
    '''Performs agglomerative clustering and extracts the cluster labels, then sorts according to
    cluster number. In the future, this will also extract a Newick tree.'''

    clustering = AgglomerativeClustering( 
        affinity = 'precomputed',
        distance_threshold = maxDist,
        n_clusters = None,
        linkage = linkage
        ).fit( distanceMatrix )

    i = iter( distanceMatrix.columns )
    clusters = { next( i ): x for x in clustering.labels_ }
    clusters = { k: v for k, v in sorted( clusters.items(), key = lambda item: item[1] ) }

    return clusters


def getNewick(node, newick, parentdist, leaf_names):
    '''Code from https://stackoverflow.com/questions/28222179/save-dendrogram-to-newick-format
    adopted from @jfn'''
    if node.is_leaf():
        return "%s:%.2f%s" % (leaf_names[node.id], parentdist - node.dist, newick)
    else:
        if len(newick) > 0:
            newick = "):%.2f%s" % (parentdist - node.dist, newick)
        else:
            newick = ");"
        newick = getNewick(node.get_left(), newick, node.dist, leaf_names)
        newick = getNewick(node.get_right(), ",%s" % (newick), node.dist, leaf_names)
        newick = "(%s" % (newick)
        return newick


def getClusterLabels( labels, clusters ):

    i = iter(labels)
    protoclusters = {next(i): x for x in clusters}
    clusters = {k: v for k, v in sorted(protoclusters.items(), key = lambda item: item[1])}

    return clusters


def scipyaggd( distMat, maxDist, method = 'single' ):
    '''Performs agglomerative clustering using SciPy and extracts the cluster labels, then sorts
    according to cluster number.'''

    squareform_matrix = squareform(distMat.values)
    linkage_matrix = hierarchy.linkage(squareform_matrix, method)
    tree = hierarchy.to_tree(linkage_matrix)
    fcluster = hierarchy.fcluster(
        linkage_matrix, maxDist, 
        criterion = 'distance'
        )
    clusters = getClusterLabels(distMat.index, fcluster)

    return clusters, tree

def dmndMain(fastaPath, minVal, output_dir, distFile, pid = True, verbose = False, cpus = 1):
    queryDB, makeDBcode = makeDmndDB('diamond', fastaPath, output_dir, cpus = cpus)
    dmndOut, dmndCode = runDmnd(
        'diamond', fastaPath, queryDB, distFile, pid = pid,
        cpus = cpus, verbose = verbose, blast = 'blastp'
        )
    distanceMatrix = readDmndDist( dmndOut, minVal, pid = pid )
#    writeDist(distFile, outMatrix)
    return distanceMatrix
    
def usearchMain(fasta, min_id, output, cpus = 1, verbose = False):
    vprint('\nusearch aligning', flush = True, v = verbose)
    runUsearch(fasta, output + '.dist', str(1 - min_id), cpus, verbose)
    distanceMatrix = importDist( output + '.dist' )
    return distanceMatrix

def readLog(log_path, newLog):
    oldLog = readJson(log_path)
    if oldLog['fasta'] == newLog['fasta'] and \
        oldLog['minimum_id'] == newLog['minimum_id'] and \
        oldLog['search_program'] == newLog['search_program']:
        newLog['iterations'] = oldLog['iterations']
    elif os.path.isfile(newLog['distance_matrix']):
        os.remove(newLog['distance_matrix'])
    return newLog

def iterativeRun(
    minseq, maxseq, maxdist, minid, distanceMatrix, linkage, focalGene, interval = 0.1,
    log_dict = None, log_path = None, minval = 0, maxval = 1, verbose = False, spacer = '\t'
    ):

    exitCode = 1
    attempt, direction, clusters, tree = 0, None, None, None
    oldDirection, oldClusters, oldTree = None, None, None
    while maxdist >= minval and maxdist <= maxval and 1 - maxdist >= minid:
        attempt += 1
        oldClusters, oldTree = clusters, tree
        clusters, tree = scipyaggd(distanceMatrix, float(maxdist), linkage)
        cluster_dict = {} # could be more efficient by grabbing in getClusterLabels
        for gene, index in clusters.items():
            if index not in cluster_dict:
                cluster_dict[index] = []
            cluster_dict[index].append(gene)
        focalLen = len(cluster_dict[clusters[focalGene]])
        vprint('\nITERATION ' + str(attempt) + ': ' + focalGene + ' cluster size: ' + str(focalLen), flush = True, v = verbose)
        vprint('\tMinimum identity: ' + str(minid) + '; Maximum Distance ' + str(maxdist), flush = True, v = verbose)
        if log_path:
            iteration_dict = {'size': focalLen, 'maximum_distance': maxdist}
            log_dict['iterations'].append(iteration_dict)
            writeJson(log_path, log_dict)
        if focalLen >= minseq:
            if maxseq:
                if focalLen <= maxseq:
                    vprint(spacer + '\tSUCCESS!', flush = True, v = verbose)
                    exitCode = 0
                    break
                else:
                    direction = 1
            else:
                vprint(spacer + '\tSUCCESS!', flush = True, v = verbose)
                exitCode = 0
                break
        else:
            direction = -1

        if oldDirection:
            if direction != oldDirection:
                eprint(spacer + 'WARNING: Overshot. Could not find parameters using current interval.', flush = True)
                vprint(spacer + 'Outputting Iterations ' + str(attempt) + ' & ' + str(attempt - 1), flush = True, v = verbose)
                exitCode = 2
                break
        else:
            oldDirection = direction
#         minid += (interval*direction)
        maxdist -= (interval*direction)

    return clusters, tree, oldClusters, oldTree, log_dict, exitCode

def reiterativeRun(
    minseq, maxseq, maxdist, minid, distanceMatrix, linkage, focalGene, interval = 0.01,
    log_dict = None, log_path = None, minval = 0, maxval = 1, verbose = False
    ):

    attempt, direction, clusters, tree = 0, None, None, None
    oldDirection, refines = None, []
#    print(maxdist, minval, maxval)
    while maxdist >= minval and maxdist <= maxval:
        attempt += 1
        clusters, tree = scipyaggd(distanceMatrix, float(maxdist), linkage)
        cluster_dict = {} # could be more efficient by grabbing in getClusterLabels
        for gene, index in clusters.items():
            if index not in cluster_dict:
                cluster_dict[index] = []
            cluster_dict[index].append(gene)
        focalLen = len(cluster_dict[clusters[focalGene]])
        vprint('\nITERATION ' + str(attempt) + ': ' + focalGene + ' cluster size: ' + str(focalLen), flush = True, v = verbose)
        vprint('\tMinimum identity: ' + str(minid) + '; Maximum Distance ' + str(maxdist), flush = True, v = verbose)
        if log_path:
            iteration_dict = {'size': focalLen, 'maximum_distance': maxdist}
            log_dict['iterations'].append(iteration_dict)
            writeJson(log_path, log_dict)
        if focalLen >= minseq:
            if maxseq:
                if focalLen <= maxseq:
                    refines.append([clusters, tree, True])
                    direction = -1
                else:
                    refines.append([clusters, tree, False])
                    direction = 1
            else:
                return clusters, tree, None, None, log_dict
        else:
            refines.append([clusters, tree, False])
            direction = -1

        if oldDirection:
            if direction != oldDirection: # if the directions changed it was
            # met/overshot
                if refines[-2][2]: # if the second to last run was successful
                    return refines[-2][0], refines[-2][1], None, None, log_dict
                else: # return the two closest to success
                    return refines[-1][0], refines[-1][1], refines[-2][0], refines[-2][1], log_dict
        else:
            oldDirection = direction
        maxdist -= (interval*direction)

#    if refines[-1][2]: # if the final run was successful
    return refines[-1][0], refines[1][1], None, None, log_dict


def main(
    fastaPath, minid, maxdist, minseq, maxseq, 
    searchProg = 'usearch', linkage = 'single',
    iterative = False, interval = None, output= None, cpus = 1,
    verbose = False, spacer = '\t', log_path = None,
    refine = False, pid = True, dmnd_dir = None
    ):

    log_dict = {
        'algorithm': searchProg, 'fasta': fastaPath, 'minimum_id': minid, 'maximum_distance': maxdist,
        'minimum_sequences': minseq, 'maximum_sequences': maxseq, 'search_program': searchProg,
        'distance_matrix': output + '.dist', 'linkage': linkage, 'focal_gene': iterative,
        'iterations': []
        }
    if log_path:
        if os.path.isfile(log_path):
            log_dict = readLog(log_path, log_dict)
        writeJson(log_path, log_dict)

    if os.path.isfile(log_dict['distance_matrix']):
        if searchProg == 'usearch':
            distanceMatrix = importDist(log_dict['distance_matrix'])
        else:
            distanceMatrix = readDmndDist(log_dict['distance_matrix'], minid, pid = pid)
    else:
        if searchProg == 'diamond': #elif if above lines not highlighted
            if not dmnd_dir:
                dmnd_dir = os.path.dirname(log_dict['distance_matrix']) + '/'
            distanceMatrix = dmndMain(
                fastaPath, minid, dmnd_dir, 
                log_dict['distance_matrix'], pid = pid, 
                verbose = verbose, cpus = cpus
                )
        elif searchProg == 'usearch':
            distanceMatrix = usearchMain(fastaPath, minid, output, cpus, verbose)

    vprint('\nClustering\n', flush = True, v = verbose)
    if iterative and interval:
        clusters, tree, oldClusters, oldTree, log_dict, exitCode = iterativeRun(
            minseq, maxseq, maxdist, minid, distanceMatrix, linkage, 
            iterative , interval = 0.1, log_dict = log_dict, log_path = log_path, minval = 0, maxval = 1,
            verbose = verbose, spacer = spacer
            )
        if refine:
            if exitCode == 0 or exitCode == 2: # parameters were met/overshot
                if len(log_dict['iterations']) > 1: # if there's room to work
                # with
                    finalIter = log_dict['iterations'][-1]
                    compIter = log_dict['iterations'][-2]
                    interval = 0.005
                    minval = min([
                            compIter['maximum_distance'], finalIter['maximum_distance']
                            ]) 
                    maxval = max([
                            compIter['maximum_distance'],
                            finalIter['maximum_distance']
                            ])
                    clusters, tree, oldClusters, oldTree, log_dict = reiterativeRun(
                        minseq, maxseq, minval, minid, distanceMatrix, linkage, iterative, 
                        interval = interval, log_dict = log_dict, log_path = log_path, minval = minval,
                        maxval = maxval, verbose = verbose
                        )
                    if oldTree:
                        vprint(
                            spacer + 'WARNING: Cluster does not exist within parameters', 
                            v = verbose, flush = True
                            )
                    else:
                        vprint(
                            spacer + 'SUCCESS!', v = verbose, flush = True
                            )
            else: # parameters were not met
                vprint(
                    spacer + 'WARNING: Could not find a cluster that meets parameters', 
                    v = verbose, flush = True
                    )
        elif exitCode == 0:
            vprint(spacer + 'SUCCESS!', v = verbose, flush = True)
        elif exitCode == 1:
            vprint(
                spacer + 'WARNING: Could not find a cluster that meets parameters', 
                v = verbose, flush = True
                )
        else:
            vprint(
                spacer + 'WARNING: Overshot parameters',
                v = verbose, flush = True
                )
    else:
        oldClusters, oldTree = None, None
        clusters, tree = scipyaggd(distanceMatrix, float(maxdist), linkage)

    return distanceMatrix, clusters, tree, oldClusters, oldTree

def writeData(tree, distanceMatrix, clusters, output):

    #clusters = scikitaggd( distanceMatrix, float(args.max_dist), args.linkage )

    out_str = '\n'.join( [i + '\t' + str(clusters[i]) for i in clusters] )
    with open( output + '.clus', 'w' ) as out:
        out.write( out_str )
    if tree:
        newick = getNewick(tree, "", tree.dist, list(distanceMatrix.index))

        with open( output + '.newick', 'w' ) as out:
            out.write(newick)


if __name__ == '__main__':

    parser = argparse.ArgumentParser( 
        description = "Hiearchical agglomerative clustering. " 
        + "Outputs `.clus` cluster assignments and `.newick` dendrogram"
        )
    parser.add_argument( '-f', '--fasta', required = True )
    parser.add_argument( 
        '-a', '--alignment', help = 'Alignment software '
        + '{"diamond", "usearch"}', default = 'diamond' 
        )
    parser.add_argument( '-o', '--output', help = 'e.g. "~/name" will output name.clus, etc' )
    parser.add_argument( '-m', '--minid', help = 'Minimum identity'
        + ' for alignment. DEFAULT 0.2', type = float, default = 0.2 )
    parser.add_argument( '-x', '--maxdist', help = 'Maximum distance'
        + ' for clustering (1 - identity). DEFAULT: 0.75', type =
        float, default = 0.75)
    parser.add_argument( '-l', '--linkage', default = 'single', help = "Linkage criterion:" \
        + " 'complete' (maximum distance), 'average', 'weighted', 'centroid', or 'single'" \
        + " (minimum distance) DEFAULT: 'single' (minimum)" )
    parser.add_argument('--iterative', help = 'Gene to iteratively cluster for. Requires ' \
        + ' --minseq, eliminates --maxdist')
    parser.add_argument('--minseq', type = int, help = 'Iteratively adjust -m and -x until a cluster ' \
        + ' of minimum sequences is obtained. Requires --iterative. DEFAULT: 100', default = 100)
    parser.add_argument('--maxseq', type = int, help = 'Maximum sequences for iterative clustering. ' \
        + 'Requires --iterative. DEFAULT 1000', default = 1000)
#    parser.add_argument('--interval', type = float, help = 'Distance adjustment for each iteration. ' \
 #       + 'Requires --iterative. DEFAULT: 0.05', default = 0.05)
    parser.add_argument('-r', '--refine', action = 'store_true', 
        help = 'Refine cluster to maximize sequences within parameters.')
    parser.add_argument( '-c', '--cpus', default = 1, type = int, \
        help = 'Cores for alignment during distance matrix construction' )
    args = parser.parse_args()


    if args.linkage not in { 'complete', 'average', 'weighted', 'centroid', 'single' }:
        eprint('\nERROR: Invalid linkage criterium', flush = True)
        sys.exit( 1 )
    elif args.alignment not in {'diamond', 'usearch'}:
        eprint('\nERROR: Invalid alignment method', flush = True)
        sys.exit( 2 )
    else:
        findExecs( [args.alignment], exit = set(args.alignment) )

    if args.iterative:
        args.maxdist = 1 - args.minid
    elif 1 - args.maxdist <= minid:
        eprint('\nWARNING: 1 - maximum distance exceeds minimum identity, clustering is ineffective', flush = True)

    fastaPath = formatPath(args.fasta)
    fa = fa2dict(fastaPath)
    if args.iterative:
        focalGene = args.iterative
        if len(fa) < args.minseq:
            eprint('\nERROR: minimum sequences is greater than fasta input', flush = True)
            sys.exit(5)
        elif focalGene not in fa:
            eprint('\nERROR: ' + focalGene + ' not in ' + fastaPath, flush = True)
            sys.exit(6)

    if args.output:
        if not os.path.isdir(formatPath(args.output)):
            os.mkdir(formatPath(args.output))
        dmnd_dir = formatPath(args.output)
        output = out_dir + os.path.basename(fastaPath)
    else:
        dmnd_dir = mkOutput(os.getcwd() + '/', 'aggClus')
        output = os.getcwd() + '/' + os.path.basename(fastaPath)

    distanceMatrix, clusters, tree, ocl, ot = main(
        fastaPath, args.minid, args.maxdist, args.minseq, args.maxseq, 
        searchProg = args.alignment, linkage = args.linkage, verbose = True,
        iterative = args.iterative, interval = 0.10, output = output, spacer = '\n',
        log_path = os.path.dirname(output) + '.' + os.path.basename(fastaPath) + '.aggClus.json',
        refine = args.refine, pid = True, dmnd_dir = dmnd_dir, cpus = args.cpus
        )

    if ocl:
        writeData(tree, distanceMatrix, clusters, output + '0')
        writeData(ot, odm, ocl, output + '1')
    else:
        writeData(tree, distanceMatrix, clusters, output)

    sys.exit( 0 )
