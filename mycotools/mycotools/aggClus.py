#! /usr/bin/env python3

from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform
#from sklearn.cluster import AgglomerativeClustering
from mycotools.lib.kontools import multisub, findExecs, formatPath
from mycotools.lib.fastatools import fasta2dict, dict2fasta
import string, argparse, os, sys, itertools, tempfile, re, random, pandas as pd, subprocess

sys.setrecursionlimit(1000000)

def splitFasta( fa_path, output ):

    fastas = []
    fa = fasta2dict( fa_path )
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


def runUsearch( fasta, output, maxdist ):

    subprocess.call([
        'usearch', '-calc_distmx', fasta,
        '-tabbedout', output, '-maxdist',
        maxdist
        ] #stdout = subprocess.PIPE,
        #stderr = subprocess.PIPE
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


if __name__ == '__main__':

    parser = argparse.ArgumentParser( 
        description = "Hiearchical agglomerative clustering. " 
        + "Outputs `.clus` cluster assignments and `.newick` dendrogram"
        )
    parser.add_argument( '-f', '--fasta', required = True )
    parser.add_argument( 
        '-d', '--distance', help = 'Alignment software '
        + '{"needle", "usearch"}', default = 'usearch' 
        )
    parser.add_argument( '-o', '--output', help = 'e.g. "~/name" will output name.clus, etc' )
    parser.add_argument( '-m', '--min_id', help = 'Minimum identity'
        + ' for alignment, written as decimal', type = float )
    parser.add_argument( '-x', '--max_dist', required = True, help = 'Maximum distance'
        + ' for clustering (1 - minimum identity), written as decimal', type = float)
    parser.add_argument( '-l', '--linkage', default = 'single', help = "Linkage criterion:" \
        + " 'complete' (maximum distance), 'average', 'weighted', 'centroid', or 'single'" \
        + " (minimum distance) DEFAULT: 'single' (minimum)" )
    parser.add_argument( '-c', '--cpus', default = 1, type = int, \
        help = 'Cores for alignment during distance matrix construction' )
    args = parser.parse_args()


    if args.linkage not in { 'complete', 'average', 'weighted', 'centroid', 'single' }:
        eprint('\nERROR: Invalid linkage criterium', flush = True)
        sys.exit( 1 )
    elif args.distance not in {'needle', 'usearch'}:
        eprint('\nERROR: Invalid alignment method', flush = True)
        sys.exit( 2 )

    if args.output:
        output = args.output
    else:
        output = args.fasta

#    if os.path.isfile(output + '.dist'):
 #       print( '\nDistance matrix found! Ignoring "--min_id". Specify new output to rerun.' , flush = True)
  #      distanceMatrix = importDist( output + '.dist' )
    if args.distance == 'needle': #elif if above lines not highlighted
        findExecs( ['needle'], exit = set('needle') )
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
        fastas = splitFasta( args.fasta, tmpdir + '/fastas' )
        for fa1 in fastas:
            for fa2 in fastas:
                if fa1 != fa2: 
                    out = prepAlign( [ fa1, fa2 ], tmpdir + '/aligns' )
                    if out[0] not in aligns:
                        aligns.add( out[0] )
                        alignCmds.append( out[1] )
        print('\tSTEP 2/3: Aligning sequences. Using ' + str(args.cpus) + ' cores', flush = True)
        codes = multisub( list(alignCmds), processes = args.cpus )
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
    elif args.distance == 'usearch':
        findExecs( ['usearch'], exit = set('usearch') )
        print('\nCalculating usearch alignments', flush = True)
        runUsearch(args.fasta, output + '.dist', str(1 - args.min_id))
        distanceMatrix = importDist( output + '.dist' )

    print('\nClustering\n', flush = True)
    clusters, tree = scipyaggd(distanceMatrix, float(args.max_dist), args.linkage)
    newick = getNewick(tree, "", tree.dist, list(distanceMatrix.index))

    #clusters = scikitaggd( distanceMatrix, float(args.max_dist), args.linkage )

    out_str = '\n'.join( [i + '\t' + str(clusters[i]) for i in clusters] )
    with open( output + '.clus', 'w' ) as out:
        out.write( out_str )
    with open( output + '.newick', 'w' ) as out:
        out.write(newick)

    sys.exit( 0 )
