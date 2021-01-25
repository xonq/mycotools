#! /usr/bin/env python3

from sklearn.cluster import AgglomerativeClustering
from mycotools.lib.kontools import multisub, findExecs
from mycotools.lib.fastatools import fasta2dict, dict2fasta
import string, argparse, os, sys, itertools, tempfile, re, random, pandas as pd, subprocess


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
            if names[1] not in dist_dict:
                dist_dict[ names[1] ] = {}
            dist_dict[ names[0] ][ names[1] ] = identity
            dist_dict[ names[1] ][ names[0] ] = identity
            out_dict[ names[0] ][ names[1] ] = identity

    distanceMatrix = pd.DataFrame( dist_dict )
    outMatrix = pd.DataFrame( out_dict )

    return distanceMatrix.fillna( 1 ), outMatrix.fillna( 1 )


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
        dist_dict[ vals[0] ][ vals[1] ] = vals[2]
        dist_dict[ vals[1] ][ vals[0] ] = vals[2]

    distanceMatrix = pd.DataFrame( dist_dict )

    return distanceMatrix.fillna( 1 )


def aggd( distanceMatrix, maxDist = 0.6, linkage = 'single' ):
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


if __name__ == '__main__':

    parser = argparse.ArgumentParser( description = "Hiearchical agglomerative clustering" )
    parser.add_argument( '-f', '--fasta',  \
        help = 'Start from fasta (generate distance matrix)' )
    parser.add_argument( '-d', '--distance', help = 'Start from distance matrix' )
    parser.add_argument( '-o', '--output', help = 'Output prefix other than input name' )
    parser.add_argument( '-m', '--max_dist', required = True, help = 'Distance or identity ' \
        + 'threshold - in an identity-based distance matrix (usearch -calc_distmx), this ' \
        + 'value is equal to 1 - minimum_identity.' )
    parser.add_argument( '-l', '--linkage', default = 'single', help = "Linkage criterion:" \
        + " 'complete' (maximum distance), 'average', or 'single' (minimum distance) " \
        + "DEFAULT: 'single' (minimum)" )
    parser.add_argument( '-c', '--cpus', default = 1, type = int, \
        help = 'Cores for alignment during distance matrix construction' )
    args = parser.parse_args()


    # `ward` can be implemented if an option is added to create the distance matrix
    if args.linkage not in { 'complete', 'average', 'single' }:
        eprint('\nERROR: Invalid linkage criterium')
        sys.exit( 1 )

    if args.output:
        output = args.output + '.clus'
    else:
        output = None

    if args.fasta:
        findExecs( ['needle'], exit = set('needle') )
        print( '\nCalculating global alignment distance matrix' )
        aligns, alignCmds = set(), list()
        if not output:
            output = args.fasta + '.clus'
        tmpdir = tempfile.gettempdir() + '/aggTmpZK' 
        tmpdir += ''.join(random.choices(string.ascii_uppercase + string.digits, k=6))
        if not os.path.isdir( tmpdir ):
            os.mkdir( tmpdir )
        if not os.path.isdir( tmpdir + '/fastas' ):
            os.mkdir( tmpdir + '/fastas' )
        if not os.path.isdir( tmpdir + '/aligns' ):
            os.mkdir( tmpdir + '/aligns' )
        print('\tSTEP 1/3: Preparing data. ' + tmpdir + ' will be removed upon completion')
        fastas = splitFasta( args.fasta, tmpdir + '/fastas' )
        for fa1 in fastas:
            for fa2 in fastas:
                if fa1 != fa2: 
                    out = prepAlign( [ fa1, fa2 ], tmpdir + '/aligns' )
                    if out[0] not in aligns:
                        aligns.add( out[0] )
                        alignCmds.append( out[1] )
        print('\tSTEP 2/3: Aligning sequences. Using ' + str(args.cpus) + ' cores')
        codes = multisub( list(alignCmds), processes = args.cpus )
        print('\tSTEP 3/3: Generating distance matrix')
        distanceMatrix, outMatrix = createDist( aligns, 1 - float(args.max_dist) )
        dist_out = ''
        for i, row in outMatrix.iterrows():
            for column in outMatrix.columns:
                if str(column) != str(i) and outMatrix.at[i, column] < 1:
                    dist_out += str(i) + '\t' + str(column) + '\t' + \
                        str(outMatrix.at[i, column]) + '\n'
        with open( re.sub( r'.clus$', '.dist', output ), 'w' ) as out:
            out.write( dist_out )
        subprocess.call( ['rm', '-rf', tmpdir], stdout = subprocess.PIPE, stderr = \
            subprocess.PIPE )

    elif args.distance:
        if not output:
            output = args.distance + '.clus'
        print('\nImporting distance matrix')
        distanceMatrix = importDist( dist )

    print('\nClustering')
    clusters = aggd( distanceMatrix, float(args.max_dist), args.linkage )

    out_str = '\n'.join( [i + '\t' + str(clusters[i]) for i in clusters] )
    with open( output, 'w' ) as out:
        out.write( out_str )

    sys.exit( 0 )
