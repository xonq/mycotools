#! /usr/bin/env python3

# NEED an option to run iterative for larger-than sets on the subset generated
    # by the larger than clustering step
# NEED to try MCL using the binary
# NEED to make outgroup detection within the bounds of maximum sequences, optional to do so
    # NEEDs heavy optimization in general
# extract a fasta of the cluster of interest
# fix output
# NEED to avoid rerunning already ran values, only when necessary
# NEED to archive old runs'
# NEED to implement linclust

import os
import re
import sys
import copy
import shutil
import string
import random
import tempfile
import argparse
import itertools
import subprocess
import pandas as pd
from collections import defaultdict
from mycotools.lib.kontools import multisub, findExecs, format_path, \
    eprint, vprint, read_json, write_json, mkOutput, fmt_float
from mycotools.lib.biotools import fa2dict, dict2fa
sys.setrecursionlimit(1000000)

class ClusteringError(Exception):
    pass

class RefineError(Exception):
    pass

def run_mmseqs(fa_path, res_base, wrk_dir, algorithm = 'mmseqs easy-linclust',
                 min_id = 0.3, min_cov = 0.5, cpus = 1, verbose = False):
    res_path = res_base + '_cluster.tsv'
    tmp_dir = os.path.dirname(res_path) + '/tmp/'
    if verbose:
        stdout, stderr = None, None
    else:
        stdout = subprocess.DEVNULL
        stderr = subprocess.DEVNULL
    mmseqs_cmd = algorithm.split(' ')
    mmseqs_cmd.extend([fa_path, res_base, tmp_dir,
                       '--min-seq-id', fmt_float(min_id), '--threads',
                       str(cpus), '--compressed', '1',
                       '--cov-mode', '0', '-c', str(min_cov),
                       '-e', '0.1'])
    mmseqs_exit = subprocess.call(mmseqs_cmd,
                                  stdout = stdout,
                                  stderr = stderr)
    if mmseqs_exit:
        raise ClusteringError('Clustering failed: ' + str(mmseqs_exit) \
                            + ' ' + str(mmseqs_cmd))
    if os.path.isfile(res_base + '_all_seqs.fasta'):
        os.remove(res_base + '_all_seqs.fasta')
    if os.path.isfile(res_base + '_rep_seq.fasta'):
        os.remove(res_base + '_rep_seq.fasta')
    if os.path.isdir(tmp_dir):
        shutil.rmtree(tmp_dir)
    return res_path

def parse_mmseqs_clus(res_path):
    derivations = defaultdict(list)
    with open(res_path, 'r') as raw:
        for line in raw:
            k, v = line.rstrip().split() # all white space
            derivations[k].append(v)

    hg_list = []
    for k, v in derivations.items():
        v.append(k)
        hg_list.append(sorted(set(v)))

    hg2gene = {i: v for i, v in enumerate(
                    sorted(hg_list, key = lambda x: len(x),
                    reverse = True)
                    )}
    gene2hg = {}
    for hg, genes in hg2gene.items():
        for gene in genes:
            gene2hg[gene] = hg
    return hg2gene, gene2hg

def makeDmndDB(diamond, queryFile, output_dir, cpus = 1):
    """create a diamond database:
    diamond: binary path, queryFile: query_path"""
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
    cpus = 1, verbose = False, blast = 'blastp', minid = 15
    ):
    """run diamond, no self-hits, more-sensitive"""

    if pid:
        distanceVal = 'pident'
    else:
        distanceVal = 'bitscore'

    cmd = [
        diamond, blast, '-d', queryDB, '-q', queryFile,
        '--out', outputFile, '--threads', str(cpus),
        '--outfmt', '6', 'qseqid', 'sseqid', distanceVal,
        '--id', str(minid), '--no-self-hits', '--more-sensitive'
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

def rd_dmnd_distmtx( outputFile, minVal, pid = True ):
    """Read a diamond output as a distance matrix"""

    dist_dict = {}
    with open(outputFile, 'r') as raw:
        if pid:
            for line in raw:
                q, s, v0 = line.rstrip().split('\t')
                v = round(float(v0))
                if q == s:
                    continue
                if v > minVal:
                    if q not in dist_dict:
                        dist_dict[q] = {}
                        dist_dict[q][q] = 0.0
                    if s not in dist_dict:
                        dist_dict[s] = {}
                        dist_dict[s][s] = 0.0
                    dist_dict[q][s] = 100 - v
                    dist_dict[s][q] = 100 - v
            distance_matrix = pd.DataFrame(dist_dict).sort_index(1).sort_index(0).fillna(100)
        else:
             for line in raw:
                q, s, v0 = line.rstrip().split('\t')
                v = float(v0)
                if q == s:
                    continue
                if v > minVal:
                    if q not in dist_dict:
                        dist_dict[q] = {}
                        dist_dict[q][q] = 0.0
                    if s not in dist_dict:
                        dist_dict[s] = {}
                        dist_dict[s][s] = 0.0
             distance_matrix = pd.DataFrame(dist_dict).sort_index(1).sort_index(0).fillna(0)


    if pid:
        distance_matrix /= 100
    else:
        minV, maxV = min(distance_matrix), max(distance_matrix)
        denom = maxV - minV
        distance_matrix -= minV
        distance_matrix /= denom
        distance_matrix = 1 - distance_matrix

    return distance_matrix #, outMatrix

def runUsearch( fasta, output, clus_parameter, cpus = 1, verbose = False ):

    if verbose:
        subprocess.call([
            'usearch', '-calc_distmx', fasta,
            '-tabbedout', output, '-clus_parameter',
            clus_parameter, '-threads', str(cpus)
            ] #stdout = subprocess.PIPE,
            #stderr = subprocess.PIPE
            )
    else:
        subprocess.call([
            'usearch', '-calc_distmx', fasta,
            '-tabbedout', output, '-clus_parameter',
            clus_parameter, '-threads', str(cpus)
            ], stdout = subprocess.DEVNULL,
            stderr = subprocess.DEVNULL
            )

def rd_usrch_distmtx( dis_path, sep = '\t' ):
    '''Imports a distance matrix with each line formatted as `organism $SEP organism $SEP distance`.
    - this is equivalent to the `-tabbedout` argument in `usearch -calc_distmx`. The function 
    compiles a dictionary with the information and reciprocal information for each organism in each
    line, then converts this dictionary of dictionaries into a pandas dataframe. As a distance matrix,
    NA values are converted to maximum distance (1).'''

    distance_matrix = pd.DataFrame()
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

    distance_matrix = pd.DataFrame( dist_dict ).sort_index(0).sort_index(1)

    return distance_matrix.fillna( 1.0 )

def scikitaggd( distance_matrix, maxDist = 0.6, linkage = 'single' ):
    '''Performs agglomerative clustering and extracts the cluster labels, then sorts according to
    cluster number. In the future, this will also extract a Newick tree.'''

    clustering = AgglomerativeClustering( 
        affinity = 'precomputed',
        distance_threshold = maxDist,
        n_clusters = None,
        linkage = linkage
        ).fit( distance_matrix )

    i = iter( distance_matrix.columns )
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

def runMCL(distMat, inflation):

    key2gene = {i: v for i, v in enumerate(list(distMat.keys()))}
    npMat = distMat.to_numpy()
    result = mc.run_mcl(npMat, inflation = inflation)
    protoClusters = mc.get_clusters(result)

    clusters = {}
    for i, keys in enumerate(protoClusters):
        clusters = {**clusters, **{key2gene[key]: i for key in keys}}

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

def dmnd_main(fasta_path, minVal, output_dir, distFile, pid = True, verbose = False, cpus = 1):
    queryDB, makeDBcode = makeDmndDB('diamond', fasta_path, output_dir, cpus = cpus)
    dmndOut, dmndCode = runDmnd(
        'diamond', fasta_path, queryDB, distFile + '.tmp', pid = pid,
        cpus = cpus, verbose = verbose, blast = 'blastp'
        )
    shutil.move(distFile + '.tmp', distFile)
    distance_matrix = rd_dmnd_distmtx( distFile, minVal, pid = pid )
    return distance_matrix
    
def usrch_main(fasta, min_id, output, cpus = 1, verbose = False):
    vprint('\nusearch aligning', flush = True, v = verbose)
    runUsearch(fasta, output + '.dist', str(1 - min_id), cpus, verbose)
    distance_matrix = rd_usrch_distmtx( output + '.dist' )
    return distance_matrix

def readLog(log_path, newLog):
    oldLog = read_json(log_path)
    if oldLog['fasta'] == newLog['fasta'] and \
        oldLog['minimum_id'] == newLog['minimum_id'] and \
        oldLog['search_program'] == newLog['search_program']:
        newLog['iterations'] = oldLog['iterations']
    elif newLog['distance_matrix']:
        if os.path.isfile(newLog['distance_matrix']):
            os.remove(newLog['distance_matrix'])
    return newLog

def run_iterative(
    minseq, maxseq, clus_parameter, mincon, params, focal_gene, 
    interval = 0.1, log_dict = None, log_path = None, minval = 0, 
    maxval = 1, verbose = False, spacer = '\t', cpus = 1
    ):
    """
    Iteratively clusters until input parameters are met.

    minseq: int(minimum sequences in a cluster)
    maxseq: int(maximum sequences in a cluster)
    clus_parameter: float(maximum_distance for agglomerative), float(minimum
        coverage for linclust)
    mincon: minimum connection value
    distance_matrix: distance_matrix generated from fa2clus functions
    focal_gene: str(gene to cluster around)
    linkage: agglomerative clustering linkage parameter
    interval: interval to adjust clus_parameter per iteration
    log_dict: log
    log_path: log_file path
    minval: int(minimum clus_parameter value)
    maxval: int(maximum clus_parameter value)
    verbose: bool(verbosity)
    spacer: str(stderr/stdout spacer)
    """

    exit_code = 1 # default exit code is failure
    attempt, direction, clusters, tree = 0, None, None, None
    oldDirection, oldClusters, oldTree = None, None, None
    while clus_parameter >= minval and clus_parameter <= maxval:
        attempt += 1 # keep track of iterations for output details
        if 'dist' in params: # hierarchical clustering
            if 1 - clus_parameter < mincon:
                break
            oldClusters, oldTree = clusters, tree
            clusters, tree = scipyaggd(params['distance_matrix'],
                                       float(clus_parameter), 
                                       params['linkage']) # cluster
            cluster_dict = defaultdict(list) # could be more efficient by grabbing in getClusterLabels
            for gene, index in clusters.items(): # create a dictionary clusID:
            # [genes]
                cluster_dict[index].append(gene)
        else:
            res_base = params['dir'] + focal_gene
            res_path = run_mmseqs(params['fa'], res_base, params['dir'], 
                                    algorithm = params['bin'],
                                    min_id = clus_parameter, verbose = verbose,
                                    min_cov = params['min_cov'], cpus = cpus)
            cluster_dict, clusters = parse_mmseqs_clus(res_path)
        focalLen = len(cluster_dict[clusters[focal_gene]]) 
        # size of cluster with focal_gene

        vprint('\nITERATION ' + str(attempt) + ': ' + focal_gene \
             + ' cluster size: ' + str(focalLen), flush = True, v = verbose)
        vprint('Cluster parameter: ' + str(clus_parameter), flush = True, v = verbose)
        if log_path: # output the iteration details to the log
            iteration_dict = {'size': focalLen, 'cluster_parameter': clus_parameter}
            log_dict['iterations'].append(iteration_dict)
            write_json(log_dict, log_path)
        if focalLen >= minseq: # if greater than minimum sequences
            if maxseq: # if there is a max set of sequences
                if focalLen <= maxseq: # if less than max sequences
                    vprint(spacer + '\tSUCCESS!', flush = True, v = verbose)
                    exit_code = 0
                    break
                else: # descend
                    direction = 1
            else: # no max sequences, minimum is met, this was successful
                vprint(spacer + '\tSUCCESS!', flush = True, v = verbose)
                exit_code = 0
                break
        else: # need to increase cluster size
            direction = -1

        if oldDirection: # if there is an older iteration
            if direction != oldDirection: # if the direction switched
                eprint(spacer + 'WARNING: Overshot - ' \
                     + 'could not find parameters using current interval', 
                       flush = True)
                vprint(spacer + 'Outputting Iterations ' + str(attempt) \
                     + ' & ' + str(attempt - 1), flush = True, v = verbose)
                exit_code = 2
                break
        else: # set the old iteration direction detail
            oldDirection = direction
#         mincon += (interval*direction)

        oclus_parameter = copy.copy(clus_parameter)
        if 'dist' in params:
            clus_parameter -= (interval*direction) 
        else:
            clus_parameter += (interval*direction)

        # clus_parameter should be to 3 significant figures
        clus_parameter = round(clus_parameter * 100)/100
        # adjust the clus_parametereter based on the direction it should change
        if oclus_parameter != minval and oclus_parameter != maxval:
            if clus_parameter < minval: # set as minimum value if it descends below
                clus_parameter = minval
            elif clus_parameter > maxval: # set as maximum if it ascends above
                clus_parameter = maxval

    return clusters, tree, oldClusters, oldTree, log_dict, exit_code

def run_reiterative(
    minseq, maxseq, clus_parameter, mincon, params, focal_gene, interval = 0.01,
    log_dict = None, log_path = None, minval = 0, maxval = 1, verbose = False,
    cpus = 1
    ):
    """Goal is to refine cluster size as close to maxseq as much as possible"""

    attempt, direction, clusters, tree = 0, None, None, None
    oldDirection, refines = None, []
    while clus_parameter >= minval and clus_parameter <= maxval: # while a valid cluster
    # parameter
        attempt += 1 # log attempts
        if 'dist' in params: # hierarchical clustering
            oldClusters, oldTree = clusters, tree
            clusters, tree = scipyaggd(params['distance_matrix'],
                                       float(clus_parameter), 
                                       params['linkage']) # cluster
            cluster_dict = defaultdict(list) # could be more efficient by grabbing in getClusterLabels
            for gene, index in clusters.items(): # create a dictionary clusID:
            # [genes]
                cluster_dict[index].append(gene)
            focalLen = len(cluster_dict[clusters[focal_gene]])
            vprint('\nITERATION ' + str(attempt) + ': ' + focal_gene \
                 + ' cluster size: ' + str(focalLen), flush = True, v = verbose)
            vprint('\tMinimum connection: ' + str(mincon) \
                 + '; Maximum cluster parameter' + str(clus_parameter), flush = True, v = verbose)
        else: # mmseqs
            res_base = params['dir'] + focal_gene
            res_path = run_mmseqs(params['fa'], res_base, params['dir'], 
                                    algorithm = params['bin'],
                                    min_id = clus_parameter, verbose = verbose,
                                    min_cov = params['min_cov'], cpus = cpus)
            cluster_dict, clusters = parse_mmseqs_clus(res_path)
            focalLen = len(cluster_dict[clusters[focal_gene]])
            vprint('\nITERATION ' + str(attempt) + ': ' + focal_gene \
                 + ' cluster size: ' + str(focalLen), flush = True, v = verbose)
            vprint('\tMinimum coverage: ' + str(mincon) \
                 + '; Minimum identity: ' + str(clus_parameter), flush = True, v = verbose)
        # find the length of the cluster with the gene
        if log_path: # output the log
            iteration_dict = {'size': focalLen, 'cluster_parameter': clus_parameter}
            log_dict['iterations'].append(iteration_dict)
            write_json(log_dict, log_path)
        if focalLen >= minseq: # if the focal cluster is > minimum size
            if maxseq: # check if there is a maximum size
                if focalLen <= maxseq: # if less than the maximum
                    refines.append([clusters, tree, True]) # append the info to
                    # refine info
                    direction = -1 # we can climb
                else:
                    refines.append([clusters, tree, False]) # else we need to
                    # go down
                    direction = 1
            else:
                # if there is no maximum size, then we met the goal
                return clusters, tree, None, None, log_dict
        else: # if less than the minimum size, then we need to climb
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

        oclus_parameter = copy.copy(clus_parameter)
        if 'dist' in params: # if hierachical agglomerative clustering
            clus_parameter -= (interval*direction)
        else:
            clus_parameter += (interval*direction)

        # clus_parameter should be to 3 significant figures
        clus_parameter = round(clus_parameter * 100)/100

        if oclus_parameter != minval and oclus_parameter != maxval:
            if clus_parameter < minval:
                clus_parameter = minval
            elif clus_parameter > maxval:
                clus_parameter = maxval

    return refines[-1][0], refines[-1][1], None, None, log_dict


def main(
    fasta_path, mincon, clus_parameter, minseq, maxseq, 
    search_program = 'diamond', linkage = 'single',
    iterative = False, interval = 0.1, output= None, cpus = 1,
    verbose = False, spacer = '\t\t', log_path = None,
    refine = False, pid = True, dmnd_dir = None
    ):

    if search_program in {'usearch', 'diamond'}:
        algorithm = 'hierarchical'
        dist_path = output + '.dist'
    else:
        algorithm = search_program
        dist_path = None

    minval, maxval = 0, 1
    overshot = False
    log_dict = {
        'algorithm': algorithm, 'fasta': fasta_path, 
        'minimum_id': mincon, 'cluster_parameter': clus_parameter,
        'minimum_sequences': minseq, 'maximum_sequences': maxseq, 
        'search_program': search_program, 'distance_matrix': dist_path, 
        'linkage': linkage, 'focal_gene': iterative,
        'iterations': []
        }
    if log_path:
        if os.path.isfile(log_path):
            log_dict = readLog(log_path, log_dict)
        write_json(log_dict, log_path)

    if algorithm == 'hierarchical':
        global hierarchy
        global squareform
        from scipy.cluster import hierarchy
        from scipy.spatial.distance import squareform
        param_dict = {'dist': None, 'link': linkage}
        if os.path.isfile(log_dict['distance_matrix']):
            if search_program == 'usearch':
                param_dict['dist'] = rd_usrch_distmtx(log_dict['distance_matrix'])
            else:
                param_dict['dist'] = rd_dmnd_distmtx(log_dict['distance_matrix'], mincon, pid = pid)
        else:
            if search_program == 'diamond': #elif if above lines not highlighted
                if not dmnd_dir:
                    dmnd_dir = os.path.dirname(log_dict['distance_matrix']) + '/'
                param_dict['dist'] = dmnd_main(
                    fasta_path, mincon, dmnd_dir, 
                    log_dict['distance_matrix'], pid = pid, 
                    verbose = verbose, cpus = cpus
                    )
            elif search_program == 'usearch':
                param_dict['dist'] = usrch_main(fasta_path, mincon, output, cpus, verbose)
    else:
        param_dict = {'fa': fasta_path, 'dir': os.path.dirname(output) + '/',
                      'bin': search_program, 'min_cov': mincon}

    vprint('\nClustering', flush = True, v = verbose)
    if iterative:
        clusters, tree, oldClusters, oldTree, log_dict, exit_code = run_iterative(
            minseq, maxseq, clus_parameter, mincon, param_dict,
            iterative, interval = interval, 
            log_dict = log_dict, log_path = log_path, minval = minval, maxval = maxval,
            verbose = verbose, spacer = spacer, cpus = cpus
            )
        if refine:
            if exit_code == 0 or exit_code == 2: # parameters were met/overshot
                if len(log_dict['iterations']) > 1: # if there's room to work
                    vprint('\nMaximizing cluster size', v = verbose, flush = True)
                    beyond_iters = [x for x in log_dict['iterations'] \
                                    if int(x['size']) > maxseq]
                                    # iterations with too many sequences
                    sorted_beyonds = sorted(beyond_iters, key = lambda x: \
                                            (x['size'], x['cluster_parameter']))
                                            # sort smallest to largest
                    beyonds = [x for x in sorted_beyonds \
                               if x['size'] == sorted_beyonds[0]['size']]
                               # grab all iterations of the fewest of too
                               # many
                    final_iteration = log_dict['iterations'][-1]
                    if algorithm == 'hierarchical':
                        minval = final_iteration['cluster_parameter']
                        new_cp = minval
                        if beyonds:
                            interval = 0.005
                            maxval = float(beyonds[0]['cluster_parameter'])
                        else:
                            interval = 0.05
                            if pid:
                                maxval = 1 - mincon # max is relative to % ID cutoff
                            else:
                                maxval = 1
                    else: # mmseqs
                        if beyonds:
                            interval = 0.01
                            minval = float(beyonds[-1]['cluster_parameter'])
                            # smallest value surpassing max sequences
                        else:
                            interval = 0.05
                            minval = 0.0
                        maxval = final_iteration['cluster_parameter']
                        # the working value
                        if minval > maxval:
                            raise RefineError('Cannot refine mmseqs cluster')
                        new_cp = maxval
                    clusters, tree, oldClusters, oldTree, log_dict = run_reiterative(
                        minseq, maxseq, new_cp, mincon, param_dict, iterative, 
                        interval = interval, log_dict = log_dict, 
                        log_path = log_path, minval = minval,
                        maxval = maxval, verbose = verbose, cpus = cpus
                        )
                    if oldTree:
                        vprint(
                            spacer + 'WARNING: Cluster does not exist within parameters', 
                            v = verbose, flush = True
                            )
            else: # parameters were not met
                vprint(
                    spacer + 'WARNING: Could not find a cluster that meets parameters', 
                    v = verbose, flush = True
                    )
        elif exit_code == 0:
            vprint(spacer + 'SUCCESS!', v = verbose, flush = True)
        elif exit_code == 1:
            vprint(
                spacer + 'WARNING: Could not find a cluster that meets parameters', 
                v = verbose, flush = True
                )
        else:
            vprint(
                spacer + 'WARNING: Overshot parameters',
                v = verbose, flush = True
                )
            overshot = True
    else:
        oldClusters, oldTree = None, None
        if algorithm == 'hierarchical':
            clusters, tree = scipyaggd(distance_matrix, clus_parameter, linkage)
        else:
            res_base = output
            res_path = run_mmseqs(fasta_path, output, 
                                    os.path.dirname(output) + '/', 
                                    verbose = verbose,
                                    algorithm = algorithm, min_id = clus_parameter, 
                                    min_cov = mincon, cpus = cpus)
            null_dict, clusters = parse_mmseqs_clus(res_path)

    if algorithm == 'hierarchical':
        return distance_matrix, clusters, tree, oldClusters, oldTree, overshot
    else:
        return None, clusters, None, oldClusters, None, overshot

def write_data(tree, distance_matrix, clusters, output):

    #clusters = scikitaggd( distance_matrix, float(args.max_dist), args.linkage )

    out_str = '\n'.join([i + '\t' + str(clusters[i]) for i in clusters])
    with open(output + '.clus', 'w') as out:
        out.write( out_str )
    if tree:
        newick = getNewick(tree, "", tree.dist, list(distance_matrix.index))
        with open(output + '.newick', 'w') as out:
            out.write(newick)


if __name__ == '__main__':

    parser = argparse.ArgumentParser( 
        description = "Sequence clustering with iterative search option." 
        )
    parser.add_argument('-f', '--fasta', required = True)
    parser.add_argument( 
        '-a', '--alignment', 
        help = 'Alignment software {"mmseqs", "diamond", "usearch"}; ' \
             + 'diamond and usearch call hierarchical clustering. ' \
             + 'DEFAULT: mmseqs', 
        default = 'mmseqs' 
        )
    parser.add_argument('-l', '--linclust', action = 'store_true',
        help = 'Use linclust instead of mmseqs cluster (faster, less sensitive)')
    parser.add_argument('-m', '--mincon', help = 'Minimum connection' \
        + ' for alignment; DEFAULT: [mmseqs]: 0.1 percent query coverage, '  \
        + '[aggclus]: 0.2 percent identity, [aggclus]: 30 bitscore', 
        type = float)
    parser.add_argument('-x', '--cluster_parameter', 
        help = '[mmseqs]: Minimum identity; DEFAULT: 0.3; ' \
            +  '[aggclus]: Maximum distance; DEFAULT: 0.75', 
            type = float)
    parser.add_argument('-i', '--iterative', 
        help = '[--minseq]: Gene to iteratively cluster for; '  \
            + 'eliminates -x')
    parser.add_argument('--minseq', type = int, help = '[-i]: Iteratively adjust -m and -x until a cluster ' \
        + ' of minimum sequences is obtained; DEFAULT: 2', default = 2)
    parser.add_argument('--maxseq', type = int, 
        help = '[-i]: Maximum sequences for iterative clustering; ' \
            + 'DEFAULT 100', default = 100)
    parser.add_argument('-r', '--refine', action = 'store_true', 
        help = '[--maxseq]: Refine cluster to maximize sequences within parameters.')
#    parser.add_argument('-i', '--inflation', type = float, default = 1.5,
 #       help = 'Inflation value for MCL; DEFAULT: 1.5')
    parser.add_argument(
        '-d', '--distance_type', help = '[aggclus]: {"identity", "bitscore"} ' + \
        'bitscore is normalized 0-1 following -m and -x should be 0-1; DEFAULT: identity',
        default = 'identity'
        )

    parser.add_argument('--linkage', default = 'single', 
        help = "[aggclus]: Linkage criterion: " \
        + "'complete' (maximum distance), 'average', 'weighted', " \
        + "'centroid', or 'single' (minimum distance) DEFAULT: 'single' (minimum)")
    parser.add_argument('-o', '--output')
    parser.add_argument('-c', '--cpus', default = 1, type = int)
    parser.add_argument('-v', '--verbose', action = 'store_true')
    args = parser.parse_args()

    if args.refine and not args.maxseq:
        eprint('\nERROR: --maxseq required for refinement', flush = True)
    if args.alignment in {'diamond', 'usearch'}:
        if args.linkage not in { 'complete', 'average', 'weighted', 'centroid', 'single' }:
            eprint('\nERROR: Invalid linkage criterium', flush = True)
            sys.exit( 1 )
        elif args.distance_type not in {'identity', 'bitscore'}:
            eprint('\nERROR: Invalid distance type', flush = True)
            sys.exit( 3 )
        if args.distance_type == 'identity':
            pid = True
            if not args.mincon:
                args.mincon = 0.2
        else:
            pid = False
            if not args.mincon:
                args.mincon = 30
        if not args.cluster_parameter:
            args.cluster_parameter = 0.75
    elif args.algnment == 'mmseqs':
        if args.linclust:
            args.alignment = 'mmseqs easy-linclust'
        else:
            args.alignment = 'mmseqs easy-cluster'
        if not args.mincon:
            args.mincon = 0.1
            pid = None
        if not args.cluster_parameter:
            args.cluster_parameter = 0.3
    else:
        eprint('\nERROR: Invalid alignment software', flush = True)
        sys.exit( 2 )

    findExecs([args.alignment], exit = set(args.alignment))
    interval = 0.1

    if args.iterative and args.alignment in {'diamond', 'usearch'}:
        clus_parameter = 1 - args.mincon
    elif 1 - args.cluster_parameter <= args.mincon \
        and args.alignment != 'mmseqs':
        eprint('\nWARNING: 1 - maximum distance exceeds minimum connection, clustering is ineffective', flush = True)
    else:
        clus_parameter = args.cluster_parameter
 #   else:
#        clus_parameter = args.inflation

    fasta_path = format_path(args.fasta)
    fa = fa2dict(fasta_path)
    if args.iterative:
        focal_gene = args.iterative
        if len(fa) < args.minseq:
            eprint('\nERROR: minimum sequences is greater than fasta input', flush = True)
            sys.exit(5)
        elif focal_gene not in fa:
            eprint('\nERROR: ' + focal_gene + ' not in ' + fasta_path, flush = True)
            sys.exit(6)

    if args.output:
        if not os.path.isdir(format_path(args.output)):
            os.mkdir(format_path(args.output))
        dmnd_dir = format_path(args.output)
        output = dmnd_dir + os.path.basename(fasta_path)
    else:
        dmnd_dir = mkOutput(os.getcwd() + '/', 'fa2clus')
        output = dmnd_dir + os.path.basename(fasta_path)

    distance_matrix, clusters, tree, ocl, ot, overshot = main(
        fasta_path, args.mincon, clus_parameter, args.minseq, args.maxseq, 
        search_program = args.alignment, linkage = args.linkage, verbose = args.verbose,
        iterative = args.iterative, interval = interval, output = output,
        log_path = dmnd_dir + '.' + os.path.basename(fasta_path) + '.fa2clus.json',
        refine = args.refine, pid = pid, dmnd_dir = dmnd_dir, cpus = args.cpus,
        spacer = '\n'
        )

    if ocl:
        write_data(tree, distance_matrix, clusters, output + '.0')
        write_data(ot, distance_matrix, ocl, output + '.1')
    else:
        write_data(tree, distance_matrix, clusters, output)

    if args.iterative:
        cluster_dict = {}
        for gene, index in clusters.items():
            if index not in cluster_dict:
                cluster_dict[index] = []
            cluster_dict[index].append(gene)
    
        gene_module = cluster_dict[clusters[focal_gene]]
        input_fa = fa2dict(fasta_path)
        output_fa = {x: input_fa[x] for x in gene_module}
   
        with open(output + '.fa', 'w') as out:
            out.write(dict2fa(output_fa)) 

    sys.exit(0)
