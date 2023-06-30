#! /usr/bin/env python3

# NEED log sorted by default, and only unique run parameters
# NEED percent positive mode
# NEED to try MCL using the binary
# NEED to make rerunning aggclus not use old data

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

class ClusterParameterError(Exception):
    pass

def run_mmseqs(fa_path, res_base, wrk_dir, algorithm = 'mmseqs easy-linclust',
                 min_id = 0.3, clus_const = 0.5, cpus = 1, verbose = False):
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
                       '--cov-mode', '0', '-c', str(clus_const),
                       '-e', '0.1', '-s', '7.5', '--createdb-mode', '0'])
                       # createdb-mode is necessary for how I format fastas due
                       # to some optimization in mmseqs, they will address in a
                       # future update
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

def runUsearch( fasta, output, clus_var, cpus = 1, verbose = False ):

    if verbose:
        subprocess.call([
            'usearch', '-calc_distmx', fasta,
            '-tabbedout', output, '-clus_var',
            clus_var, '-threads', str(cpus)
            ] #stdout = subprocess.PIPE,
            #stderr = subprocess.PIPE
            )
    else:
        subprocess.call([
            'usearch', '-calc_distmx', fasta,
            '-tabbedout', output, '-clus_var',
            clus_var, '-threads', str(cpus)
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

def dmnd_main(fa_path, minVal, output_dir, distFile, pid = True, verbose = False, cpus = 1):
    queryDB, makeDBcode = makeDmndDB('diamond', fa_path, output_dir, cpus = cpus)
    dmndOut, dmndCode = runDmnd(
        'diamond', fa_path, queryDB, distFile + '.tmp', pid = pid,
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
        newLog['successes'] = oldLog['successes']
        for i in newLog['iterations']:
            i['cluster'] = tuple(i['cluster'])
        for i in newLog['successes']:
            i['cluster'] = tuple(i['cluster'])
    elif newLog['distance_matrix']:
        if os.path.isfile(newLog['distance_matrix']):
            os.remove(newLog['distance_matrix'])
    return newLog

def extract_closest_cluster(iterations, min_seq, max_seq):
    for i in iterations:
        if i['size'] < min_seq:
            i['discrep'] = min_seq - i['size']
        else:
            i['discrep'] = i['size'] - max_seq
    sorts = sorted(iterations, key = lambda x: x['discrep'])
    return sorts[0]

def sort_iterations(iterations, reverse = True):
    return sorted([dict(t) for t in {
                    tuple(d.items()) for d in iterations
                    }], key = lambda x: (x['size'], x['cluster_variable']),
                    reverse = reverse)

def cluster_iter_mmseqs(params, min_seq, max_seq, clus_const, clus_var,
               interval, min_var, max_var, log_dict, log_path, 
               focal_gene = None, verbose = False, spacer = '\t', 
               cpus = 1):

    oldDirection, attempt = None, 0
    if focal_gene:
        res_base = params['dir'] + 'working/' + focal_gene
    else:
        res_base = params['dir'] + 'working/' \
                 + re.sub(r'\.[^\.]+$', '', params['fa']) 
    while clus_var >= min_var and clus_var <= max_var:
        attempt += 1
        name = '_c' + str(round(clus_const*1000) / 1000) \
             + '_v' + str(round(clus_var*1000) / 1000)
        res_info = res_base + name
        res_path = res_info + '_cluster.tsv'
        if not os.path.isfile(res_path):
            run_mmseqs(params['fa'], res_info, 
                       params['dir'] + 'working/',
                       algorithm = params['bin'],
                       min_id = clus_var, verbose = verbose,
                       clus_const = params['clus_const'], cpus = cpus)
        cluster_dict, clusters = parse_mmseqs_clus(res_path)
        cluster = cluster_dict[clusters[focal_gene]]
        focal_len = len(cluster) 

        vprint('\nITERATION ' + str(attempt) + ': ' + focal_gene \
             + ' cluster size: ' + str(focal_len), flush = True, v = verbose)
        vprint('Cluster parameter: ' + str(clus_var), flush = True, v = verbose)
        iteration_dict = {'size': focal_len, 'cluster_variable': clus_var,
                          'cluster': tuple(cluster)}
        log_dict['iterations'].append(iteration_dict)
        log_dict['iterations'] = sort_iterations(log_dict['iterations'])
        if focal_len >= min_seq: # if greater than minimum sequences
            if max_seq: # if there is a max set of sequences
                if focal_len <= max_seq: # if less than max sequences
                    vprint(spacer + '\tSUCCESS!', flush = True, v = verbose)
                    exit_code = 0
                    direction = -1
                    log_dict['successes'].append(iteration_dict)
                    log_dict['successes'] = \
                        sort_iterations(log_dict['successes'])
                else: # descend
                    direction = 1
            else: # no max sequences, minimum is met, this was successful
                vprint(spacer + '\tSUCCESS!', flush = True, v = verbose)
                exit_code = 0
                log_dict['successes'].append(iteration_dict)
                log_dict['successes'] = sort_iterations(log_dict['successes'])
                break
        else: # need to increase cluster size
            direction = -1

        if oldDirection: # if there is an older iteration
            if direction != oldDirection: # if the direction switched
                diff_len = ofocal_len - focal_len
                if oldDirection * -1 == diff_len/abs(diff_len):
                # percent ID change is congruent with the size change, which should
                # not be the case, i.e. decreasing % ID shouldn't decrease
                # module size
                    pass
                elif log_dict['successes']:
                    # grab the largest success
                    cluster = log_dict['successes'][0]['cluster']
                    exit_code = 0
                    break
                else:
                    eprint(spacer + 'WARNING: Overshot - ' \
                         + 'could not find parameters using current interval', 
                           flush = True)
                    iteration = extract_closest_cluster(log_dict['iterations'],
                                                        min_seq, max_seq)
                    cluster = iteration['cluster']
                    exit_code = 2
                    break
        else: # set the old iteration direction detail
            oldDirection = direction
            exit_code = 1
#         clus_const += (interval*direction)

        ofocal_len = copy.copy(focal_len)
        oclus_var = copy.copy(clus_var)
        clus_var += (interval*direction)

        # clus_var should be to 3 significant figures
        clus_var = round(clus_var * 100)/100
        # adjust the clus_vareter based on the direction it should change
        if oclus_var != min_var and oclus_var != max_var:
            if clus_var < min_var: # set as minimum value if it descends below
                clus_var = min_var
            elif clus_var > max_var: # set as maximum if it ascends above
                clus_var = max_var

    write_json(log_dict, log_path)
    try:
        return cluster, log_dict, exit_code
    except UnboundLocalError: # cluster not defined, didn't enter while loop
        raise ClusterParameterError('failed to run with given parameters')

def cluster_iter_aggclus(params, min_seq, max_seq, clus_const, clus_var,
                interval, min_var, max_var, log_dict, log_path, 
                focal_gene = None, verbose = False, spacer = '\t', 
                cpus = 1):

    oldDirection, attempt = None, 0
    if focal_gene:
        res_base = params['dir'] + 'working/' + focal_gene
    else:
        res_base = params['dir'] + 'working/' \
                 + re.sub(r'\.[^\.]+$', '', params['fa']) 
    while clus_var >= min_var and clus_var <= max_var:
        attempt += 1
        if 100 - round(100*clus_var) < round(100*clus_const):
            break
        name = '_c' + str(round(clus_const*1000) / 1000) \
             + '_v' + str(round(clus_var*1000) / 1000)
        clusters, tree = scipyaggd(params['dist'],
                                   float(clus_var), 
                                   params['link']) # cluster
        newick = getNewick(tree, "", tree.dist, list(params['dist'].index))
        write_data(newick, clusters, res_base + name)

        cluster_dict = defaultdict(list) # could be more efficient by grabbing in getClusterLabels
        for gene, index in clusters.items(): # create a dictionary clusID:
            # [genes]
            cluster_dict[index].append(gene)

        try:
            cluster = cluster_dict[clusters[focal_gene]]
            focal_len = len(cluster) 
        except KeyError: # focal gene not in output
            focal_len = 0
            cluster = None
            newick = ''


        vprint('\nITERATION ' + str(attempt) + ': ' + focal_gene \
             + ' cluster size: ' + str(focal_len), flush = True, v = verbose)
        vprint('Cluster parameter: ' + str(clus_var), flush = True, v = verbose)
        iteration_dict = {'size': focal_len, 'cluster_variable': clus_var,
                          'cluster': tuple(cluster), 'tree': newick}
        log_dict['iterations'].append(iteration_dict)
        log_dict['iterations'] = sort_iterations(log_dict['iterations'])
        if focal_len >= min_seq: # if greater than minimum sequences
            if max_seq: # if there is a max set of sequences
                if focal_len <= max_seq: # if less than max sequences
                    vprint(spacer + '\tSUCCESS!', flush = True, v = verbose)
                    exit_code = 0
                    direction = -1
                    log_dict['successes'].append(iteration_dict)
                    log_dict['successes'] = \
                        sort_iterations(log_dict['successes'])
                else: # descend
                    direction = 1
            else: # no max sequences, minimum is met, this was successful
                vprint(spacer + '\tSUCCESS!', flush = True, v = verbose)
                exit_code = 0
                log_dict['successes'].append(iteration_dict)
                log_dict['successes'] = \
                    sort_iterations(log_dict['successes'])
                break
        else: # need to increase cluster size
            exit_code = 1
            direction = -1

        if oldDirection: # if there is an older iteration
            if direction != oldDirection: # if the direction switched
                if log_dict['successes']:
                    # sort by the largest success
                    cluster = log_dict['successes'][0]['cluster']
                    newick = log_dict['successes'][0]['tree']
                    exit_code = 0
                else:
                    eprint(spacer + 'WARNING: Overshot - ' \
                         + 'could not find parameters using current interval', 
                           flush = True)
                    iteration = extract_closest_cluster(log_dict['iterations'],
                                                        min_seq, max_seq)
                    cluster = iteration['cluster']
                    newick = iteration['tree']
                    exit_code = 2
                break
        else: # set the old iteration direction detail
            oldDirection = direction
#         clus_const += (interval*direction)

        oclus_var = copy.copy(clus_var)
        clus_var -= (interval*direction)

        # clus_var should be to 3 significant figures
        clus_var = round(clus_var * 100)/100
        # adjust the clus_vareter based on the direction it should change
        if oclus_var != min_var and oclus_var != max_var:
            if clus_var < min_var: # set as minimum value if it descends below
                clus_var = min_var
            elif clus_var > max_var: # set as maximum if it ascends above
                clus_var = max_var

    write_json(log_dict, log_path)
    try:
        return cluster, newick, log_dict, exit_code
    except UnboundLocalError: # cluster not defined, didn't enter while loop
        raise ClusterParameterError('failed to run with given parameters')


def main(
    fa_path, clus_const, clus_var, min_seq, max_seq, 
    search_program = 'diamond', linkage = 'single',
    focal_gene = False, interval = 0.1, output= None, cpus = 1,
    verbose = False, spacer = '\t\t', log_path = None,
    pid = True, dmnd_dir = None, min_var = 0.05,
    max_var = 1
    ):

    if not os.path.isdir(os.path.dirname(output) + '/working/'):
        os.mkdir(os.path.dirname(output) + '/working/')

    if search_program in {'usearch', 'diamond'}:
        algorithm = 'hierarchical'
        dist_path = output + '.dist'
    else:
        algorithm = search_program
        dist_path = None

    overshot = False
    log_dict = {
        'algorithm': algorithm, 'fasta': fa_path, 
        'minimum_id': clus_const, 'cluster_variable': clus_var,
        'minimum_sequences': min_seq, 'maximum_sequences': max_seq, 
        'search_program': search_program, 'distance_matrix': dist_path, 
        'linkage': linkage, 'focal_gene': focal_gene,
        'iterations': [], 'successes': []
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
        param_dict = {'dist': None, 'link': linkage, 'dir': os.path.dirname(output) + '/'}
        if os.path.isfile(log_dict['distance_matrix']):
            if search_program == 'usearch':
                param_dict['dist'] = rd_usrch_distmtx(log_dict['distance_matrix'])
            else:
                param_dict['dist'] = rd_dmnd_distmtx(log_dict['distance_matrix'], clus_const, pid = pid)
        else:
            if search_program == 'diamond': #elif if above lines not highlighted
                if not dmnd_dir:
                    dmnd_dir = os.path.dirname(log_dict['distance_matrix']) + '/'
                param_dict['dist'] = dmnd_main(
                    fa_path, clus_const, dmnd_dir, 
                    log_dict['distance_matrix'], pid = pid, 
                    verbose = verbose, cpus = cpus
                    )
            elif search_program == 'usearch':
                param_dict['dist'] = usrch_main(fa_path, clus_const, output, cpus, verbose)
    else:
        param_dict = {'fa': fa_path, 'bin': search_program, 
                      'clus_const': clus_const, 
                      'dir': os.path.dirname(output) + '/'}

    if focal_gene:
        vprint('\nClustering', flush = True, v = verbose)
        if algorithm == 'hierarchical':
            cluster, newick, log_dict, error = cluster_iter_aggclus(
                   param_dict, min_seq, max_seq, clus_const, clus_var,
                   interval, min_var, max_var, log_dict, log_path, 
                   focal_gene, verbose, spacer, cpus
                   )
            if error == 2:
                overshot = True
            else:
                overshot = False
            return cluster, newick, overshot, log_dict
        else:
            cluster, log_dict, error = cluster_iter_mmseqs(
                   param_dict, min_seq, max_seq, clus_const, clus_var,
                   interval, min_var, max_var, log_dict, log_path, 
                   focal_gene, verbose, spacer, cpus
                   )
            if error == 2:
                overshot = True
            else:
                overshot = False
            return cluster, None, overshot, log_dict
    else:
        if focal_gene:
            res_base = param_dict['dir'] + focal_gene
        else:
            res_base = param_dict['dir'] \
                     + re.sub(r'\.[^\.]+$', '', os.path.basename(fa_path))
        if algorithm == 'hierarchical':
            clusters, tree = scipyaggd(param_dict['dist'],
                                       float(clus_var), 
                                       param_dict['link']) # cluster
            newick = getNewick(tree, "", tree.dist, 
                               list(param_dict['dist'].index))
            write_data(newick, clusters, res_base)
            return clusters, None, None, log_dict
        else:
            res_path = run_mmseqs(param_dict['fa'], res_base, 
                                    param_dict['dir'] + 'working/',
                                    algorithm = param_dict['bin'],
                                    min_id = clus_var, verbose = verbose,
                                    clus_const = param_dict['clus_const'], cpus = cpus)
            cluster_dict, clusters = parse_mmseqs_clus(res_path)
            write_data(None, clusters, param_dict['dir'], res_base)
            return clusters, None, None, log_dict
    
def write_data(newick, clusters, output):

    #clusters = scikitaggd( distance_matrix, float(args.max_dist), args.linkage )
    out_str = '\n'.join([i + '\t' + str(clusters[i]) for i in clusters])
    with open(output + '.clus', 'w') as out:
        out.write(out_str)
    if newick:
        with open(output + '.newick', 'w') as out:
            out.write(newick)


def cli():

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
    parser.add_argument('-m', '--cluster_constant', 
        help = 'DEFAULT: [mmseqs]: 0.4 percent query coverage, '  \
             + '[aggclus]: 0.2 percent identity, [aggclus]: 30 bitscore', 
        type = float)
    parser.add_argument('-x', '--cluster_variable', 
        help = '[mmseqs]: Minimum identity; DEFAULT: 0.3; ' \
            +  '[aggclus]: Maximum distance; DEFAULT: 0.75', 
            type = float)
    parser.add_argument('-i', '--iterative', 
        help = '[--min_seq]: Gene to iteratively cluster for; '  \
            + 'eliminates -x')
    parser.add_argument('--min_seq', type = int, help = '[-i]: Iteratively adjust -m and -x until a cluster ' \
        + ' of minimum sequences is obtained; DEFAULT: 2', default = 2)
    parser.add_argument('--max_seq', type = int, 
        help = '[-i]: Maximum sequences for iterative clustering; ' \
            + 'DEFAULT 100', default = 100)
#    parser.add_argument('-r', '--refine', action = 'store_true', 
 #       help = '[--max_seq]: Refine cluster to maximize sequences within parameters.')
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

#    if args.refine and not args.max_seq:
 #       eprint('\nERROR: --max_seq required for refinement', flush = True)
    if args.alignment in {'diamond', 'usearch'}:
        if args.linkage not in {'complete', 'average', 'weighted', 'centroid', 'single'}:
            eprint('\nERROR: Invalid linkage criterium', flush = True)
            sys.exit(1)
        elif args.distance_type not in {'identity', 'bitscore'}:
            eprint('\nERROR: Invalid distance type', flush = True)
            sys.exit(3)
        if args.distance_type == 'identity':
            pid = True
            if not args.cluster_constant:
                args.cluster_constant = 0.2
        else:
            pid = False
            if not args.cluster_constant:
                args.cluster_constant = 30
        if not args.cluster_variable:
            args.cluster_variable = 0.75
    elif args.alignment == 'mmseqs':
        if args.linclust:
            args.alignment = 'mmseqs easy-linclust'
        else:
            args.alignment = 'mmseqs easy-cluster'
        pid = None
        if not args.cluster_constant:
            args.cluster_constant = 0.4
        if not args.cluster_variable:
            args.cluster_variable = 0.3
    else:
        eprint('\nERROR: Invalid alignment software', flush = True)
        sys.exit(2)

    findExecs([args.alignment.split()[0]], 
               exit = set(args.alignment.split()[0]))
    interval = 0.1

#    if args.iterative and args.alignment in {'diamond', 'usearch'}:
 #       clus_var = 1 - args.cluster_constant
    if 1 - args.cluster_variable <= args.cluster_constant \
        and args.alignment != 'mmseqs':
        eprint('\nWARNING: 1 - maximum distance exceeds minimum connection, clustering is ineffective', flush = True)
        sys.exit(3)
    else:
        clus_var = args.cluster_variable
 #   else:
#        clus_var = args.inflation

    fa_path = format_path(args.fasta)
    fa = fa2dict(fa_path)
    if args.iterative:
        if len(fa) < args.min_seq:
            eprint('\nERROR: minimum sequences is greater than fasta input', flush = True)
            sys.exit(5)
        elif args.iterative not in fa:
            eprint('\nERROR: ' + args.iterative + ' not in ' + fa_path, flush = True)
            sys.exit(6)

    if args.output:
        if not os.path.isdir(format_path(args.output)):
            os.mkdir(format_path(args.output))
        dmnd_dir = format_path(args.output)
        output = dmnd_dir \
            + re.sub(r'\.[^\.]+$', '', os.path.basename(fa_path))
    else:
        dmnd_dir = mkOutput(os.getcwd() + '/', 'fa2clus')
        output = dmnd_dir \
            + re.sub(r'\.[^\.]+$', '', os.path.basename(fa_path))

    cluster, tree, overshot, log_dict = main(
        fa_path, args.cluster_constant, clus_var, args.min_seq, args.max_seq, 
        search_program = args.alignment, linkage = args.linkage, verbose = args.verbose,
        focal_gene = args.iterative, interval = interval, output = output,
        log_path = dmnd_dir + '.' + os.path.basename(fa_path) + '.fa2clus.json',
        pid = pid, dmnd_dir = dmnd_dir, cpus = args.cpus,
        spacer = '\n'
        )

    if tree:
        with open(output + '.newick', 'w') as out:
            out.write(tree)
    if args.iterative:
        input_fa = fa2dict(fa_path)
        try:
            output_fa = {x: input_fa[x] for x in cluster}
        except TypeError:
            eprint('\nERROR: empty cluster', flush = True)
            sys.exit(10)
        with open(output + '.fa', 'w') as out:
            out.write(dict2fa(output_fa)) 

    sys.exit(0)


if __name__ == '__main__':
    cli()
