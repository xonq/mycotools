#! /usr/bin/env python3

# NEED end of contig option in gff2svg
# NEED an intelligent taxonomy checker to validate orthogroup tag is within DB scope
# NEED to figure out min/maxseq for aggclus on the fly
# NEED logoutput and intelligent resume
# NEED to check for multiple instances of a single orthogroup/query, and use all genes as focal
    # color code differently if not part of the same cluster eventually
        # or a merged query e.g. query1 | query2
    # preliminary output svgs, then change if necessary
# NEED a relational query to new query name table

from mycotools.lib.dbtools import mtdb, masterDB
from mycotools.lib.kontools import eprint, formatPath, findExecs, intro, outro
from mycotools.lib.biotools import fa2dict, dict2fa, gff2list, list2gff, gff3Comps
from mycotools.acc2fa import dbmain as acc2fa
from mycotools.aggClus import writeData, main as aggClus
from mycotools.fa2tree import main as fa2tree, PhyloError
from mycotools.acc2loci import main as acc2loci
from mycotools.gff2svg import main as gff2svg
from mycotools.db2search import main as db2search
from mycotools.utils.og2mycodb import mycodbOGs, extractOGs
try:
    from ete3 import Tree, faces, TreeStyle
except ImportError:
    raise ImportError('Install ete3 into your conda environment via `conda install ete3`')
import argparse, re, sys, os, datetime, multiprocessing as mp

os.environ['QT_QPA_PLATFORM'] = 'offscreen'

def prepGff(gff, prots, comps, hits = set(), parDict = {}):
    RNA = False
    for entry in gff:
        if ';SearchQuery=' in entry['attributes']:
            continue
        if entry['type'] == 'gene':
            geneID = re.search(comps['id'], entry['attributes'])[1]
            if geneID in prots:
                entry['attributes'] += ';SearchQuery=' + geneID
                hits.add(geneID)
            elif geneID in parDict:
                entry['attributes'] += ';SearchQuery=' + parDict[geneID]
                del parDict[geneID]
            else:
                protID = re.search(comps['prot'], entry['attributes'])[1]
                if protID in prots:
                    entry['attributes'] += ';SearchQuery=' + protID
                    hits.add(protID)
        elif 'RNA' in entry['type']:
            RNA = True
            rnaID = re.search(comps['id'], entry['attributes'])[1]
            parID = re.search(comps['par'], entry['attributes'])[1]
            try:
                protID = re.search(comps['prot'], entry['attributes']).groups()[1]
                if protID in prots:
                    entry['attributes'] += ';SearchQuery=' + protID
                    parDict[parID] = protID
                    hits.add(protID)
                    continue
            except TypeError:
                pass
            if parID in prots:
                entry['attributes'] += ';SearchQuery=' + parID
                hits.add(parID)
            elif rnaID in parDict:
                entry['attributes'] += ';SearchQuery=' + parDict[rnaID]
                parDict[parID] = parDict[rnaID]
                del parDict[rnaID]
            else:
                protID = re.search(comps['prot'], entry['attributes'])[1]
                if protID in prots:
                    entry['attributes'] += ';SearchQuery=' + protID
                    hits.add(protID)
        elif entry['type'] == 'CDS':
            try:
                protID = re.search(comps['prot'], entry['attributes'])[1]
            except TypeError: # protein field not available
                continue
            if protID in prots:
                entry['attributes'] += ';SearchQuery=' + protID
                parID = re.search(comps['par'], entry['attributes'])[1]
                parDict[parID] = protID
                hits.add(protID)
    return parDict, hits, RNA, gff

def rogueLocus(locusID, rnaGFF, wrk_dir, query2color):
    with open(wrk_dir + 'genes/' + locusID + '.locus.genes', 'w') as out:
            out.write(list2gff(rnaGFF))
    svg_path = wrk_dir + 'svg/' + locusID + '.locus.svg'
    gff2svg(rnaGFF, svg_path, product_dict = query2color, prod_comp = r';SearchQuery=([^;]+$)', width = 10, null = 'na')


def inputGenes2inputOGs(inputGenes, gene2og, ogtag):
    inputOGs = {}
    for gene in inputGenes:
        try:
            inputOGs[gene] = gene2og[gene]
        except KeyError:
            eprint('ERROR: ' + gene + ' no valid OG under tag ' + ogtag, flush = True)
            eprint('\t' + gene + ' will be ignored. Future updates will place the gene in an OG', flush = True)
    return inputOGs

def compileOGfa(db, gene_list, og):
    fa_dict = acc2fa(db, gene_list)
    return og, fa_dict

def checkFaSize(fas, max_size):
    fas4clus, fas4trees = {}, {}
    for name, fa in fas.items():
        print('\t' + name + '\t' + str(len(fa)) + ' genes', flush = True)
        if len(fa) > max_size:
            fas4clus[name] = fa
        else:
            fas4trees[name] = fa

    fas4clus = {k: v for k, v in sorted(fas4clus.items(), key = lambda x: len(x[1]))}
    fas4trees = {k: v for k, v in sorted(fas4trees.items(), key = lambda x: len(x[1]))}

    return fas4clus, fas4trees

def runAggclus(
    fa_path, db, focalGene, minseq, maxseq, output, out_name,
    minid = 0.3, maxdist = 0.65, direction = -1, cpus = 1,
    verbose = False,
    ):

    output_path = output + str(focalGene)    
    tree, clusters, distanceMatrix, ot, ocl, odm  = aggClus(
        fa_path, minid, maxdist, minseq, maxseq,
        searchProg = 'usearch', linkage = 'single', cpus = cpus,
        iterative = focalGene, interval = 0.05, output = output_path,
        verbose = verbose
        )

    cluster_dict = {}
    for gene, index in clusters.items():
        if index not in cluster_dict:
            cluster_dict[index] = []
        cluster_dict[index].append(gene)

    geneModule = cluster_dict[clusters[focalGene]]
    fa_dict = acc2fa(db, geneModule)
# need to implement some method to choose if the max and min parameters couldn't be met
#    if ocl:
 #       cd = {}
  #      for gene, index in ocl.items():
   #         if index not in cd:
    #            cd[index] = []
     #       cd[index].append(gene)
      #  focalLen0 = len(cluster_dict[clusters[focalGene]])
       # focalLen1 = len(cd[ocl[focalGene]])
        #if not focalLen0 

 
#        if direction < 0: # if preferring the minimal size
 #           if focalLen0 > focalLen1 and not :
  #              geneModule = cd[ocl[focalGene]]
   #     else:
    #        if focalLen0 < focalLen1:
    print('\t\t\t' + str(len(fa_dict)) + ' genes in cluster', flush = True)
    with open(output + '../' + str(out_name) + '.fa', 'w') as out:
        out.write(dict2fa(fa_dict))
    writeData(None, distanceMatrix, clusters, output_path)


def makeOutput(base_dir):

    if not base_dir:
#        if not os.path.exists(base_dir):
 #           eprint('\nERROR: base output directory missing: ' + base_dir, flush = True)
  #          sys.exit(2)
        curdate = datetime.datetime.now().strftime('%Y%m%d')
        output_dir = os.getcwd() + '/crap_' + curdate + '/'
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)
    else:
        if not os.path.isdir(base_dir):
            os.mkdir(base_dir)
        output_dir = formatPath(base_dir)

    wrk_dir = output_dir + 'working/'
    global svg_dir # needs to be global for etetree

    svg_dir = wrk_dir + 'svg/'
    gff_dir = wrk_dir + 'genes/'
    tre_dir = wrk_dir + 'trees/'
    if not os.path.isdir(wrk_dir):
        os.mkdir(wrk_dir)
    if not os.path.isdir(svg_dir):
        os.mkdir(svg_dir)
    if not os.path.isdir(gff_dir):
        os.mkdir(gff_dir)


    return output_dir, wrk_dir, gff_dir, tre_dir
    

def gene2og2ome2og(gene2og):
    omeGene2og = {}
    for gene, og in gene2og.items():
        ome = gene[:gene.find('_')]
        if ome not in omeGene2og:
            omeGene2og[ome] = {}
        omeGene2og[ome][gene] = og
    return omeGene2og


def compileGenesByOme(inputs, wrk_dir):
    locusIDs, omeGenes = [], {}
    for query in inputs:
        fa = fa2dict(wrk_dir + str(query) + '.fa')
        locusIDs.extend(list(fa.keys()))

        for locusID in locusIDs:
            ome = locusID[:locusID.find('_')]
            if ome not in omeGenes:
                omeGenes[ome] = []
            omeGenes[ome].append(locusID)
    omeGenes = {ome: list(set(genes)) for ome, genes in omeGenes.items()}

    return omeGenes

def compileGenesByOme4queries(queries, wrk_dir):
    # need to account for hits that appear in multiple files
    queryGenes = {}
    for query in queries:
        if os.path.isfile(wrk_dir + str(query) + '.fa'):
            fa = fa2dict(wrk_dir + str(query) + '.fa')
        elif os.path.isfile(wrk_dir + 'clus/' + str(query) + '.fa'):
            fa = fa2dict(wrk_dir + 'clus/' + str(query) + '.fa')
        else:
            eprint('\tERROR: no fasta for ' + query, flush = True)
            continue
        locusIDs = list(fa.keys())

        for locusID in locusIDs:
            ome = locusID[:locusID.find('_')]
            if ome not in queryGenes:
                queryGenes[ome] = {}
            queryGenes[ome][locusID] = query

    return queryGenes



def extractLocusOG(ome, genesTograb, ogs, omeGene2og, plusminus, og2color, wrk_dir):
    gff_list = gff2list(formatPath('$MYCOGFF3/' + ome + '.gff3'))
    out_indices, geneGffs = acc2loci(gff_list, genesTograb, plusminus, mycotools = True, geneGff = True)

    extractedGenes = {}    
    for locusID, genes in out_indices.items():
        startI, endI = None, None
        for i, gene in enumerate(genes):
            try:
                geneGffs[locusID][i]['attributes'] += ';OG=' + str(omeGene2og[gene])
            except KeyError: #gene not in gene2og
                geneGffs[locusID][i]['attributes'] += ';OG=na'
                continue
            if gene == locusID or omeGene2og[gene] in ogs:
                if startI is None:
                    startI = i
                else:
                    endI = i
        if not endI:
            endI = startI
        extractedGenes[locusID] = geneGffs[locusID][startI:endI+1]

    for locusID, geneGff in extractedGenes.items():
        with open(wrk_dir + 'genes/' + locusID + '.locus.genes', 'w') as out:
            out.write(list2gff(geneGff))

    for locusID, geneGff in extractedGenes.items():
        svg_path = wrk_dir + 'svg/' + locusID + '.locus.svg'
        gff2svg(geneGff, svg_path, product_dict = og2color, prod_comp = r';OG=([^;]+$)', width = 10, null = 'na')
   
    return extractedGenes


def extractLocusGene(ome, accs, gene2query, plusminus, query2color, wrk_dir):
    gff_list = gff2list(formatPath('$MYCOGFF3/' + ome + '.gff3'))

    try:
        out_indices, geneGffs = acc2loci(gff_list, accs, plusminus, mycotools = True, geneGff = True)
    except KeyError:
        eprint('\t\t' + ome + ' incorrectly formatted GFF3', flush = True)
        return

    extractedGenes = {}
    for locusID, genes in out_indices.items():
        if os.path.isfile(wrk_dir + 'svg/' + locusID + '.locus.svg'):
            continue
        startI, endI = None, None
        for i, gene in enumerate(genes):
            if gene in gene2query:
                geneGffs[locusID][i]['attributes'] += ';SearchQuery=' + gene2query[gene]
                if startI is None:
                    startI = i
                else:
                    endI = i
            else:
                geneGffs[locusID][i]['attributes'] += ';SearchQuery=na'
                continue
        if not endI:
            endI = startI
        extractedGenes[locusID] = geneGffs[locusID][startI:endI+1]

    for locusID, geneGff in extractedGenes.items():
        with open(wrk_dir + 'genes/' + locusID + '.locus.genes', 'w') as out:
            out.write(list2gff(geneGff))

    for locusID, geneGff in extractedGenes.items():
        svg_path = wrk_dir + 'svg/' + locusID + '.locus.svg'
        gff2svg(geneGff, svg_path, product_dict = query2color, prod_comp = r';SearchQuery=([^;]+$)', width = 10, null = 'na')
   
    return extractedGenes


def svg2node(node):
    if node.is_leaf():
        svg_path = svg_dir + node.name + '.locus.svg'
        nodeFace = faces.ImgFace(svg_path)
        nodeFace.margin_top = 10
        nodeFace.margin_bottom = 10
        nodeFace.border.margin = 1
        faces.add_face_to_node(nodeFace, node, column = 0)
 #   else:
#        node.add_features(confidence = node.support) #redundant

def svgs2tree(inputGene, og, tree_file, out_dir):#svg_dir, out_dir):

    with open(tree_file, 'r') as raw:
        tree = Tree(raw.read())
    ts = TreeStyle()
    ts.layout_fn = svg2node
    ts.show_branch_support = True
    if og is not None:
        tree.render(out_dir + inputGene + '_OG' + str(og) + '.svg', w=1000, tree_style = ts)
    else:
        tree.render(out_dir + inputGene + '.svg', w=1000, tree_style = ts)


def makeColorPalette(inputs):
    colors = [
        "#000000","#004949","#009292","#ff6db6","#ffb6db",
        "#490092","#006ddb","#b66dff","#6db6ff","#b6dbff",
        "#920000","#924900","#db6d00","#24ff24","#ffff6d"
        ]
    extColors = [
        "F0A3FF", "#0075DC", "#993F00", "#4C005C", "#191919",
        "#005C31", "#2BCE48", "#FFCC99", "#808080", "#94FFB5",
        "#8F7C00", "#9DCC00", "#C20088", "#003380", "#FFA405",
        "#FFA8BB", "#426600", "#FF0010", "#5EF1F2", "#00998F",
        "#E0FF66", "#740AFF", "#990000", "#FFFF80", "#FFFF00",
        "#FF5005"
        ]

    color_dict = {'na': '#ffffff'}
    if len(inputs) <= len(colors):
        for i, v in enumerate(inputs):
            color_dict[str(v)] = colors[i]
    elif len(inputs) <= len(extColors):
        for i, v in enumerate(inputs):
            color_dict[str(v)] = extColors[i]
    else:
        eprint('\nWARNING: input too large for gene arrow colors', flush = True)
        try:
            for i, v in enumerate(inputs):
                color_dict[str(v)] = extColors[i]
        except IndexError:
            for i1, v in enumerate(inputs[i:]):
                color_dict[str(v)] = extColors[i1-i]

    return color_dict


def treeMngr(
    query, out_dir, wrk_dir, tre_dir, fast, treeSuffix, cpus = 1, verbose = False
    ):

    queryFa_path = wrk_dir + query + '.fa'
    if os.path.isfile(tre_dir + str(query) + treeSuffix):
        return
    else:
        try:
            fa2tree(queryFa_path, output_dir = tre_dir, fast = fast, cpus = cpus, verbose = verbose)
        except PhyloError:
            return query


def initLog(
    log_file, log_dict
    ):
    with open(log_file, 'w') as out:
        out.write(
            'og_file\t' + log_dict['og_file'] + '\n' + \
            'plusminus\t' + str(log_dict['plusminus']) + '\n' + \
            'pair_percentile\t' + str(log_dict['pair_percentile']) + '\n' + \
            'ogx_percentile\t' + str(log_dict['ogx_percentile']) + '\n' + \
            'border_percentile\t' + str(log_dict['border_percentile']) + '\n' + \
            'inflation\t' + str(log_dict['inflation']) + '\n' + \
            'null_samples\t' + str(log_dict['null_samples']) + '\n' + \
            'n50\t' + str(log_dict['n50'])
            )

def readLog(
    log_file, log_dict
    ):
    log_res = {}
    with open(log_file, 'r') as raw:
        for line in raw:
            key = line[:line.find('\t')]
            res = line[line.find('\t') + 1:].rstrip()
            if res != str(log_dict[key]):
                log_res[key] = False
            else:
                log_res[key] = True
    try:
        if not log_res['n50']:
            log_res['plusminus'] = False
        if not log_res['plusminus']:
            log_res['null_samples'] = False
        if not log_res['null_samples']:
            log_res['pair_percentile'] = False
        if not log_res['pair_percentile']:
            log_res['border_percentile'] = False
        if not log_res['border_percentile']:
            log_res['inflation'] = False
        if not log_res['inflation']:
            log_res['ogx_percentile'] = False
    except KeyError:
        print('\nERROR: corrupted log.txt.' + \
            '\nIf not rectified, future runs will completely overwrite the current\n')
        sys.exit(149)

    return log_res


def rmOldData(
    log_res, out_dir, wrk_dir
    ):
    if not log_res['null_samples']:
        nulls = collect_files(wrk_dir, 'null.txt')
        for null in nulls:
            os.remove(null)
    if not log_res['pair_percentile']:
        seed_file = out_dir + 'seed_scores.tsv.gz'
        seed_arr = wrk_dir + 'ogpair.arr.npy'
        if os.path.isfile(seed_file):
            os.remove(seed_file)
        if os.path.isfile(seed_arr):
            os.remove(seed_arr)
        clus_pickle = wrk_dir + 'clus_loc.pickle'
        ogx_pickle = wrk_dir + 'ogx_scores.pickle'
        ome_pickle = wrk_dir + 'ogx_omes.pickle'
        if os.path.isfile(clus_pickle):
            os.remove(clus_pickle)
        if os.path.isfile(ogx_pickle):
            os.remove(ogx_pickle)
        if os.path.isfile(ome_pickle):
            os.remove(ome_pickle)
    if not log_res['border_percentile']:
        row_file = wrk_dir + 'mtx/mcl.prep.rows'
        prep_file = wrk_dir + 'mtx/mcl.prep.gz'
        if os.path.isfile(row_file):
            os.remove(row_file)    
        if os.path.isfile(prep_file):
            os.remove(prep_file)
    if not log_res['inflation']:
        if os.path.isfile(wrk_dir + 'mtx.tar.gz'):
            os.remove(wrk_dir + 'mtx.tar.gz')
        if os.path.isfile(wrk_dir + 'mcl.res.gz'):
            os.remove(wrk_dir + 'mcl.res.gz')
        if os.path.isdir(wrk_dir + 'mtx/'):
            shutil.rmtree(wrk_dir + 'mtx/')
    if not log_res['ogx_percentile']:
        kern_file = out_dir + 'kernel_scores.tsv.gz'
        clus_file = out_dir + 'ogx_scores.tsv.gz'
        patch_pickle = wrk_dir + 'patchiness.nulls.pickle'
        ome_dir = wrk_dir + 'ome/'
        ogx_dir = wrk_dir + 'ogx/'
        hmm_dir = wrk_dir + 'hmm/'
        if os.path.isfile(kern_file):
            os.remove(kern_file)
        if os.path.isfile(clus_file):
            os.remove(clus_file)
        if os.path.isdir(ome_dir):
            shutil.rmtree(ome_dir)

def crapMngr(
    query, queryHits, out_dir, wrk_dir, tre_dir, fast, 
    treeSuffix, genes2query, plusminus, query2color, 
    cpus = 1, verbose = False, og = False, ogs = None
    ):

    info = treeMngr(
        query, out_dir, wrk_dir, tre_dir, fast, treeSuffix, cpus, verbose
        )

    if info:
        return

    print('\t\tExtracting loci and generating synteny diagrams', flush = True)
    extractLoci_cmds = []
    if og:
        for ome, ome_genes2og in genes2query.items():
            omeHits = [x for x in queryHits if x.startswith(ome + '_')]
            extractLoci_cmds.append([ome, omeHits, ogs, ome_genes2og, plusminus, query2color, wrk_dir])
        with mp.Pool(processes = cpus) as pool:
            pool.starmap(extractLocusOG, extractLoci_cmds)
    else:
        for ome, ome_genes2query in genes2query.items():
            omeHits = [x for x in queryHits if x.startswith(ome + '_')]
            extractLoci_cmds.append([ome, omeHits, ome_genes2query, plusminus, query2color, wrk_dir])
        with mp.Pool(processes = cpus) as pool:
            pool.starmap(extractLocusGene, extractLoci_cmds)
    
    print('\t\tMapping synteny diagrams on phylogeny', flush = True)
    tree_file = tre_dir + query + treeSuffix
    svgs2tree(
        query, None, tree_file, out_dir #svg_dir, out_dir
        )



def OGmain(
    db, inputGenes, ogtag, fast = True, output = None,
    minid = 0.3, maxdist = 0.65, minseq = 20, max_size = 250, cpus = 1,
    plusminus = 5, verbose = False
    ):
    '''inputGenes is a list of genes within an inputted cluster'''

    print('\nPreparing output directory', flush = True)
    out_dir, wrk_dir, gff_dir, tre_dir = makeOutput(output)


    print('\nCompiling orthogroup data', flush = True)
    print('\tCompiling orthogroups', flush = True)
    ogInfo_dict = mycodbOGs(omes = set(db['internal_ome']))
    og2gene, gene2og = extractOGs(ogInfo_dict, ogtag)
    inputOGs = inputGenes2inputOGs(inputGenes, gene2og, ogtag)
    og2color = makeColorPalette(inputOGs)


    # in the future, genes without OGs will be placed into OGs via RBH
    if not inputOGs:
        eprint('\nERROR: no OGs for any inputted genes', flush = True)
        sys.exit(3)

    print('\tPreparing orthogroup fastas', flush = True)
    og_fas = {}
    compileOGfa_cmds = [
        [db, og2gene[og], og] for og in list(inputOGs.values()) if not os.path.isfile(wrk_dir + str(og) + '.fa')
        ]
    with mp.Pool(processes = cpus) as pool:
        og_fas = {x[0]: x[1] for x in pool.starmap(compileOGfa, compileOGfa_cmds)}
    for og in inputOGs.values():
        if os.path.isfile(wrk_dir + str(og) + '.fa'):
            og_fas = {**og_fas, **fa2dict(wrk_dir + str(og) + '.fa')}

    print('\nChecking fasta sizes', flush = True)
    fas4clus, fas4trees = checkFaSize(og_fas, max_size)
    for query, hitFa in fas4trees.items():
        hitFa_path = wrk_dir + query + '.fa'
        with open(hitFa_path, 'w') as out:
            out.write(dict2fa(hitFa))

    if fas4clus:
        clus_dir = wrk_dir + 'clus/'
        if not os.path.isdir(clus_dir):
            os.mkdir(clus_dir)
        for query, fa in fas4clus.items():
            clusFa_path = clus_dir + query + '.fa'
            if not os.path.isfile(clusFa_path):
                with open(clusFa_path, 'w') as out:
                    out.write(dict2fa(fa))
        print('\nRunning hierarchical agglomerative clustering on ' + str(len(fas4clus)) + ' fastas', flush = True)

    print('\nCRAP', flush = True)
    omeGene2og = gene2og2ome2og(gene2og)
    if fast:
        treeSuffix = '.fa.clipkit.treefile'
    else:
        treeSuffix = '.fa.clipkit.contree'
    fas4trees = {k: v for k, v in sorted(fas4trees.items(), key = lambda x: len(x[1]))}
    for query, queryFa in fas4trees.items():
        queryHits = list(queryFa.keys())
        print('\tQuery: ' + str(query), flush = True)
        crapMngr(
            query, queryHits, out_dir, wrk_dir, tre_dir, fast, treeSuffix, 
            omeGene2og, plusminus, og2color, cpus, verbose, inputOGs[query], list(inputOGs.values())
            )

    for query in fas4clus:
        print('\tQuery: ' + str(query), flush = True)
        print('\t\tHierarchical agglomerative clustering', flush = True)
        runAggclus(
            clus_dir + str(query) + '.fa', db, query, minseq, max_size, clus_dir,
            query, minid, maxdist, cpus = cpus, verbose = verbose
            )
        queryFa = fa2dict(wrk_dir + query + '.fa')
        queryHits = list(queryFa.keys())
        crapMngr(
            query, queryHits, out_dir, wrk_dir, tre_dir, fast, treeSuffix,
            omeGene2og, plusminus, og2color, cpus, verbose, inputOGs[query], list(inputOGs.values())
            )


def SearchMain(
    db, inputGenes, queryFa, queryGff, binary = 'mmseqs', fast = True, output = None,
    minid = 0.3, maxdist = 0.65, minseq = 20, max_size = 250, cpus = 1,
    plusminus = 5, evalue = None, bitscore = 40, pident = 0, mem = None, verbose = False
    ):
    '''inputGenes is a list of genes within an inputted cluster'''

    print('\nPreparing run', flush = True)
    out_dir, wrk_dir, gff_dir, tre_dir = makeOutput(output)
    clus_dir = wrk_dir + 'clus/'

    query2color = makeColorPalette(inputGenes)

    if queryGff:
        print('\tCleaning input GFF', flush = True)
        parDict, protHits, RNA, queryGff = prepGff(queryGff, set(inputGenes), gff3Comps())
        count = 0
        while parDict and count < 4:
            count += 1
            parDict, protHits, RNA, queryGff = prepGff(queryGff, set(inputGenes), gff3Comps(), protHits, parDict)
        if parDict:
            eprint('\nIncorrectly formatted GFF', flush = True)
            sys.exit(5)
        elif set(inputGenes).difference(protHits):
            eprint('\nProteins missing from GFF', flush = True)
            eprint('\t' + ','.join([str(x) for x in list(set(inputGenes).difference(protHits))]))
            sys.exit(6)
        cleanGff = [x for x in queryGff if ';SearchQuery=' in x['attributes'] and 'RNA' in x['type']]
        for query in inputGenes:
            rogueLocus(query, cleanGff, wrk_dir, query2color)

    query_path = wrk_dir + 'query.crap.fa'
    if not queryFa:
        queryFa = acc2fa(db, inputGenes)
        with open(query_path, 'w') as out:
            out.write(dict2fa(queryFa))
    elif not os.path.isfile(query_path):
        with open(query_path, 'w') as out:
            out.write(dict2fa(queryFa))

    search_fas = {}
    for query in queryFa:
        if os.path.isfile(wrk_dir + query + '.fa'):
            search_fas[query] = fa2dict(wrk_dir + query + '.fa')
        elif os.path.isfile(clus_dir + query + '.fa'):
            search_fas[query] = fa2dict(clus_dir + query + '.fa')

    omes = set(db['internal_ome'])
    if not len(search_fas) == len(queryFa) and False is True: # need more intelligent check
        print('\nRunning search against database', flush = True)
        search_fas = {**search_fas, **db2search(
            db, binary, query_path, wrk_dir, evalue = evalue, bitscore = bitscore,
            pident = pident, mem = mem, biotype = 'proteome', force = True
            )}
        if queryGff:
            for query in search_fas:
                search_fas[query][query] = queryFa[query]
    
    print('\nChecking hit fasta sizes', flush = True)
    fas4clus, fas4trees = checkFaSize(search_fas, max_size)
    for query, hitFa in fas4trees.items():
        hitFa_path = wrk_dir + query + '.fa'
        with open(hitFa_path, 'w') as out:
            out.write(dict2fa(hitFa))

    if fas4clus:
        if not os.path.isdir(clus_dir):
            os.mkdir(clus_dir)
        for query, fa in fas4clus.items():
            clusFa_path = clus_dir + query + '.fa'
            if not os.path.isfile(clusFa_path):
                with open(clusFa_path, 'w') as out:
                    out.write(dict2fa(fa))
        print('\tRunning hierarchical agglomerative clustering on ' + str(len(fas4clus)) + ' fastas', flush = True)

    print('\nCRAP', flush = True)
    if fast:
        treeSuffix = '.fa.clipkit.treefile'
    else:
        treeSuffix = '.fa.clipkit.contree'
    fas4trees = {k: v for k, v in sorted(fas4trees.items(), key = lambda x: len(x[1]))}
    genes2query = compileGenesByOme4queries(inputGenes, wrk_dir)
    for query, queryFa in fas4trees.items():
        queryHits = list(queryFa.keys())
        print('\tQuery: ' + str(query), flush = True)
        crapMngr(
            query, queryHits, out_dir, wrk_dir, tre_dir, 
            fast, treeSuffix, genes2query, plusminus, query2color, cpus, verbose
            )
        
    for query in fas4clus:
        print('\tQuery: ' + str(query), flush = True)
        print('\t\tHierarchical agglomerative clustering', flush = True)
        runAggclus(
            clus_dir + str(query) + '.fa', 
            db, query, minseq, max_size, 
            clus_dir, query, minid, maxdist, cpus = cpus, verbose = verbose
            )
        queryFa = fa2dict(wrk_dir + query + '.fa')
        queryHits = list(queryFa.keys())
        crapMngr(
            query, queryHits, out_dir, wrk_dir, tre_dir, 
            fast, treeSuffix, genes2query, plusminus, query2color, cpus, verbose
            )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = 'Mycotools integrated Cluster Reconstruction and Phylogeny (CRAP) pipeline'
        )
    parser.add_argument(
        '-i', '--input', 
        help = 'Fasta or white-space delimited file/string of cluster genes',
        required = True
        )
    parser.add_argument('-d', '--database', default = masterDB())
    parser.add_argument(
        '-s', '--search',
        help = 'Search binary {mmseqs, blastp} for search-based CRAP'
        )
    parser.add_argument(
        '-og', '--orthogroups', 
        help = 'MycotoolsDB Orthogroup tag for OG searches. DEFAULT: P for phylum',
        default = 'P'
        )
    parser.add_argument(
        '-p', '--plusminus',
        help = 'Genes up-/downstream to analyze from loci. DEFAULT: 7',
        default = 7, type = int
        )
    parser.add_argument(
        '-m', '--maximum', 
        help = 'Max sequences for trees/min for aggClus.py. DEFAULT: 250', 
        default = 250, type = int
        )
    parser.add_argument('-f', '--fast', action = 'store_true', help = 'Fast IQtree construction.')
    parser.add_argument('--interval', help = 'Agglomerative clustering identity/distance interval. ' \
        + 'DEFAULT: 0.05', default = 0.05, type = float)
    parser.add_argument(
        '-b', '--bitscore', type = float,
        default = 30, help = 'Bitscore minimum for search algorithm. DEFAULT: 30'
        )
    parser.add_argument(
        '-g', '--gff',
        help = 'GFF for non-mycotools input. Requires -s and a fasta for -i'
        )
    parser.add_argument('-o', '--output', help = 'Output base dir')
    parser.add_argument('-c', '--cpu', default = 1, type = int)
    parser.add_argument('-v', '--verbose', default = False, action = 'store_true')
    args = parser.parse_args()

    execs = ['usearch', 'clipkit', 'mafft', 'iqtree']
    if args.search:
        if args.search not in {'mmseqs', 'blastp'}:
            eprint('\nERROR: invalid -b', flush = True)
            sys.exit(3)
        else:
            execs.append(args.search)
            args.orthogroups = None
    findExecs(execs, exit = set(execs))

    inputFa, inputGFF = False, False
    if os.path.isfile(args.input):
        if args.input.lower().endswith((
            '.fasta', '.fa', '.faa', '.fna', '.fsa'
            )):
            inputFa = fa2dict(args.input)
            inputGenes = list(inputFa.keys())
            if args.gff:
                inputGFF = gff2list(formatPath(args.gff))
        else:
            with open(args.input, 'r') as raw:
                data = raw.read()
            inputGenes = data.rstrip().split()
    elif "'" in args.input or '"' in args.input:
        if args.gff:
            eprint('\nERROR: GFF input requires fasta input', flush = True)
            sys.exit(3)
        inputGenes = args.input.replace('"','').replace("'",'').split()
    else:
        eprint('\nERROR: invalid input', flush = True)
        sys.exit(2)


    if round(len(inputGenes)/2) > args.plusminus:
        eprint('\nERROR: -p is less than input cluster size', flush = True)
        sys.exit(1)

    output = formatPath(args.output)
    args_dict = {
        'Input': ','.join(inputGenes), 
        'Database': args.database,
        'Orthogroup Tag': args.orthogroups,
        'Search binary': args.search,
        'Locus +/-': args.plusminus,
        'Fast tree': args.fast,
        'Maximum seq': args.maximum,
        'Bitscore': args.bitscore,
        'GFF': args.gff,
        'CPU': args.cpu,
        'Output directory': args.output,
        'Verbose': args.verbose
        }

    start_time = intro('Mycotools OHCrap', args_dict, 'Jason Slot, Zachary Konkel')

    db = mtdb(args.database)
    gene0 = inputGenes[0]
    ome = inputGenes[0][:inputGenes[0].find('_')]
    if not ome in set(db['internal_ome']):
        print('\nDetected non-mycotools input', flush = True)
        if not inputGFF or not inputFa or not args.search:
            eprint('\tnon-mycotools input requires -s, -i as a fasta & -g', flush = True)
            sys.exit(4)


    if args.orthogroups:    
        OGmain(
            db, inputGenes, args.orthogroups, fast = args.fast, output = output, 
            minid = 0.3, maxdist = 0.65, minseq = 20, max_size = args.maximum, cpus = args.cpu,
            verbose = args.verbose, plusminus = args.plusminus
            )
    else:
        SearchMain(
            db, inputGenes, inputFa, inputGFF, binary = args.search, fast = args.fast, output = output,
            minid = 0.3, maxdist = 0.65, minseq = 20, max_size = 250, cpus = 1,
            plusminus = args.plusminus, bitscore = args.bitscore, pident = 0, mem = None, verbose = args.verbose
            )
    outro(start_time)
