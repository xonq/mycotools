#! /usr/bin/env python3

from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform
#from sklearn.cluster import AgglomerativeClustering
from mycotools.lib.kontools import multisub, findExecs, formatPath, eprint, vprint, readJson, writeJson
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


def prepAlign( fas, output ):

    fas.sort()
    outprep = (os.path.basename(fas[0]) + '.v.' + os.path.basename(fas[1])).replace('.fa','')
    outname = output + '/' + outprep
    cmd = [
        'needle', '-outfile', outname, 
        '-asequence', fas[0], '-bsequence', fas[1],
        '-gapopen=10', '-gapextend=0.5'
        ]

    return (outname, cmd)

 
def grabIdentity( align ):

    with open( align, 'r' ) as raw:
        data = raw.read()
    idPrep = re.search( r'Identity:\W+\d+\/\d+\W+\((.*?)\%\)', data, re.M )
    identity = float( idPrep[1] ) / 100

    return identity


def createDist( alignments, minIdentity ):

    dist_dict, out_dict = {}, {}
    out_df = pd.DataFrame()
    for align in alignments:
        names = os.path.basename( align ).split( '.v.' )
        identity = grabIdentity( align )
        if identity > minIdentity:
            if names[0] not in out_dict:
                out_dict[ names[0] ] = {}
            if names[0] not in dist_dict:
                dist_dict[ names[0] ] = {}
                dist_dict[names[0]][names[0]] = 0.0
            if names[1] not in dist_dict:
                dist_dict[ names[1] ] = {}
                dist_dict[names[1]][names[1]] = 0.0
            dist_dict[ names[0] ][ names[1] ] = float(identity)
            dist_dict[ names[1] ][ names[0] ] = float(identity)
            out_dict[ names[0] ][ names[1] ] = float(identity)

    distanceMatrix = pd.DataFrame( dist_dict ).sort_index(1).sort_index(0)
    outMatrix = pd.DataFrame( out_dict ).sort_index(1).sort_index(0)

    return distanceMatrix.fillna( 1 ), outMatrix.fillna( 1 )


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
        criterion = 'distance')
    clusters = getClusterLabels(distMat.index, fcluster)

    return clusters, tree


def needleMain(fasta, min_id, cpus = 1):
    print( '\nCalculating needle alignments' , flush = True)
    aligns, alignCmds = set(), list()
    tmpdir = tempfile.gettempdir() + '/aggTmpZK' 
    tmpdir += ''.join(random.choices(string.ascii_uppercase + string.digits, k=6))
    if not os.path.isdir( tmpdir ):
        os.mkdir( tmpdir )
    if not os.path.isdir( tmpdir + '/fastas' ):
        os.mkdir( tmpdir + '/fastas' )
    if not os.path.isdir( tmpdir + '/aligns' ):
        os.mkdir( tmpdir + '/aligns' )
    print('\tSTEP 1/3: Preparing data. ' + tmpdir + ' will be removed upon completion', flush = True)
    fastas = splitFasta( fasta, tmpdir + '/fastas' )
    for fa1 in fastas:
        for fa2 in fastas:
            if fa1 != fa2: 
                out = prepAlign( [ fa1, fa2 ], tmpdir + '/aligns' )
                if out[0] not in aligns:
                    aligns.add( out[0] )
                    alignCmds.append( out[1] )
    print('\tSTEP 2/3: Aligning sequences. Using ' + str(cpus) + ' cores', flush = True)
    codes = multisub( list(alignCmds), processes = cpus )
    print('\tSTEP 3/3: Generating distance matrix', flush = True)
    distanceMatrix, outMatrix = createDist( aligns, min_id )
    dist_out = ''
    for i, row in outMatrix.iterrows():
        for column in outMatrix.columns:
            if str(column) != str(i) and outMatrix.at[i, column] < 1:
                dist_out += str(i) + '\t' + str(column) + '\t' + \
                    str(outMatrix.at[i, column]) + '\n'
    with open( output + '.dist', 'w' ) as out:
        out.write( dist_out )
    subprocess.call( ['rm', '-rf', tmpdir], stdout = subprocess.PIPE, stderr = \
        subprocess.PIPE )

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
    return newLog

def iterativeRun(maxdist, minid, distanceMatrix, linkage, focalGene, log_dict = None, log_path = None, minval = 0, maxval = 1):

    attempt, direction, clusters, distanceMatrix, tree = 0, None, None, None, None
    oldDirection, oldClusters, oldTree = None, None, None, None
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
        vprint('\nITERATION ' + str(attempt) + ': ' + focalGene + ' cluster size: ' + str(focalLen), flush = True)
        vprint('\tMinimum identity: ' + str(minid) + '; Maximum Distance ' + str(maxdist), flush = True, v = verbose)
        if log_path:
            iteration_dict = {'size': focalLen, 'maximum_distance': maxdist}
            log_dict['iterations'].append(iteration_dict)
            writeJson(log_path, log_dict)
        if focalLen >= minseq:
            if maxseq:
                if focalLen <= maxseq:
                    vprint('\t\tSUCCESS!', flush = True, v = verbose)
                    break
                else:
                    direction = 1
            else:
                vprint('\t\tSUCCESS!', flush = True, v = verbose)
                break
        else:
            direction = -1

        if oldDirection:
            if direction != oldDirection:
                eprint(spacer + 'WARNING: Overshot. Could not find parameters using current interval.', flush = True)
                vprint('Outputting Iterations ' + str(attempt) + ' & ' + str(attempt - 1), flush = True, v = verbose)
                break
        else:
            oldDirection = direction
#         minid += (interval*direction)
        maxdist -= (interval*direction)
    else:
        eprint(spacer + 'WARNING: Could not find a maximum distance that satisfies parameters', flush = True)

    return clusters, tree, oldClusters, oldTree, log_dict    

def main(
    fastaPath, minid, maxdist, minseq, maxseq, 
    searchProg = 'usearch', linkage = 'single', 
    iterative = False, interval = 0.05, output= None, cpus = 1,
    verbose = False, spacer = '\t', log_path = None,
    refine = False
    ):

    log_dict = {
        'fasta': fastaPath, 'minimum_id': minid, 'maximum_distance': maxdist,
        'minimum_sequences': minseq, 'maximum_sequences': maxseq, 'search_program': searchProg,
        'linkage': linkage, 'interval': interval, 'focal_gene': iterative,
        'iterations': []
        }
    if log_path:
        if os.path.isfile(log_path):
            log_dict = readLog(log_path, log_dict)
        writeJson(log_path, log_dict)

    if searchProg == 'needle': #elif if above lines not highlighted
        distanceMatrix = needleMain(fastaPath, minid, cpus = 1)
    elif searchProg == 'usearch':
        distanceMatrix = usearchMain(fastaPath, minid, output, cpus, verbose)

    vprint('\nClustering\n', flush = True, v = verbose)
    if iterative:
        clusters, tree, oldClusters, oldTree, log_dict = iterativeRun(
            maxdist, minid, distanceMatrix, linkage, 
            focalGene, log_dict = log_dict, log_path = log_path, minval = 0, maxval = 1
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
        + '{"needle", "usearch"}', default = 'usearch' 
        )
    parser.add_argument( '-o', '--output', help = 'e.g. "~/name" will output name.clus, etc' )
    parser.add_argument( '-m', '--minid', help = 'Minimum identity'
        + ' for alignment. DEFAULT 0.3', type = float, default = 0.3 )
    parser.add_argument( '-x', '--maxdist', help = 'Maximum distance'
        + ' for clustering (1 - minimum identity). DEFAULT: 0.6', type = float, default = 0.6)
    parser.add_argument( '-l', '--linkage', default = 'single', help = "Linkage criterion:" \
        + " 'complete' (maximum distance), 'average', 'weighted', 'centroid', or 'single'" \
        + " (minimum distance) DEFAULT: 'single' (minimum)" )
    parser.add_argument('--iterative', help = 'Gene to iteratively cluster for. Requires ' \
        + ' --minseq.')
    parser.add_argument('--minseq', type = int, help = 'Iteratively adjust -m and -x until a cluster ' \
        + ' of minimum sequences is obtained. Requires --iterative. DEFAULT: 100', default = 100)
    parser.add_argument('--maxseq', type = int, help = 'Maximum sequences for iterative clustering. ' \
        + 'Requires --iterative. DEFAULT 1000', default = 1000)
    parser.add_argument('--interval', type = float, help = 'Distance adjustment for each iteration. ' \
        + 'Requires --iterative. DEFAULT: 0.05', default = 0.05)
    parser.add_argument('-r', '--refine', action = 'store_true', help = 'Refine --interval when ' \
        + 'parameters are met.')
    parser.add_argument( '-c', '--cpus', default = 1, type = int, \
        help = 'Cores for alignment during distance matrix construction' )
    args = parser.parse_args()


    if args.linkage not in { 'complete', 'average', 'weighted', 'centroid', 'single' }:
        eprint('\nERROR: Invalid linkage criterium', flush = True)
        sys.exit( 1 )
    elif args.alignment not in {'needle', 'usearch'}:
        eprint('\nERROR: Invalid alignment method', flush = True)
        sys.exit( 2 )
    else:
        findExecs( [args.alignment], exit = set(args.alignment) )

    if 1 - args.maxdist <= minid:
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
        output = formatPath(args.output)
    else:
        output = os.getcwd() + '/' + os.path.basename(fastaPath)

    distanceMatrix, clusters, tree, ocl, ot = main(
        fastaPath, args.minid, args.maxdist, args.minseq, args.maxseq, 
        searchProg = args.alignment, linkage = args.linkage, verbose = True,
        iterative = args.iterative, interval = args.interval, output = output, spacer = '\n',
        log_path = os.dirname(output) + '.' + os.path.basename(fastaPath) + '.aggClus.json',
        refine = args.refine
        )

    if ocl:
        writeData(tree, distanceMatrix, clusters, output + '0')
        writeData(ot, odm, ocl, output + '1')
    else:
        writeData(tree, distanceMatrix, clusters, output)

    sys.exit( 0 )
