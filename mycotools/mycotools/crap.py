#! /usr/bin/env python3

# NEED intelligent resume
# NEED to check for multiple instances of a single orthogroup/query, and use all genes as focal
# NEED to update OGmain with conversion function
# NEED root to be based off the furthest in outgroup % ID from focal gene
    # dont know if possible since diamond doesnt calc all pairwise if low enough
# NEED assembly reference method
# NEED to pull distance matrix from blast results

import os
import re
import sys
import argparse
import datetime
import hashlib
import random
import multiprocessing as mp
import shutil
try:
    from ete3 import Tree, faces, TreeStyle, NodeStyle, AttrFace
except ImportError:
    raise ImportError('Install ete3 into your conda environment via `conda install ete3`')
from mycotools.lib.dbtools import mtdb, masterDB
from mycotools.lib.kontools import eprint, format_path, findExecs, intro, outro, \
    read_json, write_json
from mycotools.lib.biotools import fa2dict, dict2fa, gff2list, list2gff, gff3Comps
from mycotools.acc2fa import dbmain as acc2fa
from mycotools.fa2clus import writeData, main as fa2clus
from mycotools.fa2tree import main as fa2tree, PhyloError
from mycotools.acc2locus import main as acc2locus
from mycotools.gff2svg import main as gff2svg
from mycotools.db2search import main as db2search
from mycotools.ome2name import main as ome2name
from mycotools.utils.og2mycodb import mycodbOGs, extractOGs
os.environ['QT_QPA_PLATFORM'] = 'offscreen'

def parseConversion(conversion_file):
    with open(conversion_file, 'r') as raw:
        data = [x.rstrip().split('\t') for x in raw]
    return {x[0]: x[1] for x in data}

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

def rogueLocus(locusID, rnaGFF, wrk_dir, query2color, labels = True):
    with open(wrk_dir + 'genes/' + locusID + '.locus.genes', 'w') as out:
            out.write(list2gff(rnaGFF))
    svg_path = wrk_dir + 'svg/' + locusID + '.locus.svg'
    gff2svg(
        rnaGFF, svg_path, product_dict = query2color, labels = labels,
        prod_comp = r';SearchQuery=([^;]+$)', width = 10, null = 'na'
        )


def inputGenes2inputOGs(inputGenes, gene2og, ogtag):
    inputOGs = {}
    for gene in inputGenes:
        try:
            inputOGs[gene] = gene2og[gene]
        except KeyError:
            eprint('ERROR: ' + gene + ' no valid OG under tag ' + ogtag, flush = True)
            eprint('\t' + gene + ' will be ignored. Future updates will place the gene in an OG', flush = True)
    return inputOGs

def compileOGfa(db, gene_list, query):
    fa_dict = acc2fa(db, gene_list)
    return query, fa_dict


def checkFaSize(fas, max_size):
    fas4clus, fas4trees = {}, {}
    for query, fa in fas.items():
        print('\t' + query + '\t' + str(len(fa)) + ' genes', flush = True)
        if len(fa) > max_size:
            fas4clus[query] = fa
        else:
            fas4trees[query] = fa

    fas4clus = {k: v for k, v in sorted(fas4clus.items(), key = lambda x: len(x[1]))}
    fas4trees = {k: v for k, v in sorted(fas4trees.items(), key = lambda x: len(x[1]))}

    return fas4clus, fas4trees

def runFA2clus(
    fa_path, db, focalGene, minseq, maxseq, output, out_name,
    minid = 0.05, clusParam = 0.65, direction = -1, cpus = 1,
    verbose = False, interval = 0.1
    ):

    dmnd_dir = output + 'dmnd/'
    if not os.path.isdir(dmnd_dir):
        os.mkdir(dmnd_dir)
    output_path = output + str(focalGene) 
    log_path = output + '.' + str(focalGene) + '.log'   
    try:
        tree, clusters, distanceMatrix, ot, ocl = fa2clus(
            fa_path, minid, clusParam, minseq, maxseq,
            searchProg = 'diamond', linkage = 'single', cpus = cpus,
            iterative = focalGene, interval = interval, output = output_path,
            verbose = verbose, log_path = log_path, refine = True, dmnd_dir = dmnd_dir
            )
    except KeyError: # query not in adjacency matrix
        return False
    
    cluster_dict = {}
    for gene, index in clusters.items():
        if index not in cluster_dict:
            cluster_dict[index] = []
        cluster_dict[index].append(gene)

    geneModule = cluster_dict[clusters[focalGene]]
    fa_dict = acc2fa(db, geneModule)
# need to implement some method to choose if the max and min parameters couldn't be met
    print('\t\t\t' + str(len(fa_dict)) + ' genes in group', flush = True)
    with open(output + '../' + str(out_name) + '.fa', 'w') as out:
        out.write(dict2fa(fa_dict))
    writeData(None, distanceMatrix, clusters, output_path)
    return True

def outgroupMngr(
    db, focalGene, minseq, maxseq, clus_dir,
    minid = 0.05, direction = -1, cpus = 1, interval = 0.1,
    verbose = False, spacer = '\t\t\t'
    ):

    fa_path = clus_dir + focalGene + '.fa'
    out_name = str(focalGene) + '.outgroup'
    log_path = clus_dir + '.' + str(focalGene) + '.log'
    clusLog = read_json(log_path)
    iterations = sorted(
        clusLog['iterations'], key = lambda x: x['size'], reverse = True
        )
    prevSize = len(fa2dict(clus_dir + '../' + str(focalGene) + '.fa'))

    prevIndex = [int(v['size']) for i,v in enumerate(iterations)].index(prevSize)
    newInfo = iterations[prevIndex-1]

    if prevIndex > 0:
        res = runFA2clus(
            fa_path, db, focalGene, minseq, None, clus_dir, out_name,
            minid, float(newInfo['cluster_parameter']), direction,
            cpus, verbose, None
            )
    elif 1-float(iterations[prevIndex-1]['cluster_parameter']) - minid < 0:
        eprint(
            spacer + 'WARNING: cannot find outgroups with minimum ID', flush = True
            )
        return False # can't go any further in the clusters
    else: # can cluster further than what's been done, need to run iterative
        clusParam = float(newInfo['cluster_parameter'])
        minseq = int(newInfo['size'] + 1)
        res = runFA2clus(
            fa_path, db, focalGene, minseq, None, clus_dir, out_name,
            minid, clusParam, direction, cpus, verbose, 0.01
            )

    return res
    

def makeOutput(base_dir, newLog):

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
        output_dir = format_path(base_dir)


    logPath = output_dir + '.craplog.json'
    parseLog(logPath, newLog, output_dir)
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

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

def compileGenesByOme4queries(search_fas, conversion_dict, omes = set()):
    queryGenes = {}
    for query, fa in search_fas.items():
        locusIDs = list(fa.keys())
        for locusID in locusIDs:
            ome = locusID[:locusID.find('_')]
            if ome not in omes:
                continue
            if ome not in queryGenes:
                queryGenes[ome] = {}
            if locusID not in queryGenes[ome]:
                queryGenes[ome][locusID] = []
            queryGenes[ome][locusID].append(conversion_dict[query])

    merges = []
    for ome in queryGenes:
        for locusID in queryGenes[ome]:
            if len(queryGenes[ome][locusID]) > 1:
                queryGenes[ome][locusID] = sorted(queryGenes[ome][locusID])
            merges.append(tuple(queryGenes[ome][locusID]))
    merges = list(set(merges))

    return queryGenes, merges



def extractLocusOG(gff3, ome, genesTograb, ogs, omeGene2og, plusminus, og2color, wrk_dir, labels = True):
    gff_list = gff2list(gff3)
    genesTograb = [x for x in genesTograb if not os.path.isfile(wrk_dir + 'svg/' + x + '.locus.svg')]
    try:
        out_indices, geneGffs = acc2locus(gff_list, genesTograb, plusminus, mycotools = True, geneGff = True)
    except KeyError:
        eprint('\t\t\t' + ome + ' incorrectly formatted GFF3', flush = True)
        return

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
        gff2svg(
            geneGff, svg_path, product_dict = og2color, prod_comp = r';OG=([^;]+$)', 
            width = 10, null = 'na', labels = labels
            )
   
    return extractedGenes


def extractLocusGene(gff3, ome, accs, gene2query, plusminus, query2color, wrk_dir, labels = True):
    gff_list = gff2list(gff3)
    accs = [x for x in accs if not os.path.isfile(wrk_dir + 'svg/' + x + '.locus.svg')]
    try:
        out_indices, geneGffs = acc2locus(gff_list, accs, plusminus, mycotools = True, geneGff = True)
    except KeyError:
        eprint('\t\t\t' + ome + ' incorrectly formatted GFF3', flush = True)
        return

    extractedGenes = {}
    for locusID, genes in out_indices.items():
        if os.path.isfile(wrk_dir + 'svg/' + locusID + '.locus.svg'):
            continue
        startI, endI = None, None
        for i, gene in enumerate(genes):
            if gene in gene2query:
                geneGffs[locusID][i]['attributes'] += ';SearchQuery=' + '|'.join(gene2query[gene])
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
        gff2svg(
            geneGff, svg_path, product_dict = query2color, labels = labels,
            prod_comp = r';SearchQuery=([^;]+$)', width = 8, null = 'na'
            )
   
    return extractedGenes


def svg2node(node):
    if node.is_leaf():
        try:
            accName = re.search(r'_([^_]+_[^_]+$)', node.name)[1]
            svg_path = svg_dir + accName + '.locus.svg'
            nodeFace = faces.ImgFace(svg_path)
            nodeFace.margin_top = 5
            nodeFace.margin_bottom = 5
            nodeFace.border.margin = 1
            faces.add_face_to_node(nodeFace, node, column = 0)
        except TypeError:
            pass

def svgs2tree(inputGene, og, tree_data, out_dir, rootKey = None):#svg_dir, out_dir):
    tree = Tree(tree_data)
    if rootKey:
        tree.set_outgroup(rootKey)
    ts = TreeStyle()
    ts.layout_fn = svg2node
    ts.show_branch_support = True
    ts.branch_vertical_margin = 20
    ts.scale = 400
    nstyle = NodeStyle()
    nstyle["size"] = 0
    for n in tree.traverse():
        n.set_style(nstyle)
    if og is not None:
        tree.render(out_dir + inputGene + '_OG' + str(og) + '.svg', w=800, tree_style = ts)
    else:
        tree.render(out_dir + inputGene + '.svg', w=800, tree_style = ts)

def mergeColorPalette(merges, query2color):
    for merge in merges:
        query2color['|'.join(merge)] = query2color[merge[0]]
    return query2color

def makeColorPalette(inputs, conversion_dict = {}):
    if isinstance(inputs, list): # genes list
        for i in inputs:
            if i not in conversion_dict:
                conversion_dict[i] = str(i)
    else:
        for i, v in inputs.items(): # OGs
            if v not in conversion_dict:
                conversion_dict[v] = str(v)

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
            color_dict[conversion_dict[v]] = colors[i]
    elif len(inputs) <= len(extColors):
        for i, v in enumerate(inputs):
            color_dict[conversion_dict[v]] = extColors[i]
    else:
        eprint('\nWARNING: input too large for gene arrow colors', flush = True)
        try:
            for i, v in enumerate(inputs):
                color_dict[conversion_dict[v]] = extColors[i]
        except IndexError:
            for i1, v in enumerate(inputs[i:]):
                color_dict[conversion_dict[v]] = extColors[i1-i]

    return color_dict, conversion_dict


def treeMngr(
    query, out_dir, wrk_dir, tre_dir, fast, treeSuffix, 
    cpus = 1, verbose = False, reoutput = True
    ):

    queryFa_path = wrk_dir + query + '.fa'
    if os.path.isfile(tre_dir + str(query) + treeSuffix) and reoutput:
        return
    else:
        try:
            fa2tree(queryFa_path, output_dir = tre_dir, fast = fast, cpus = cpus, verbose = verbose)
        except PhyloError:
            return query


def initLog(
    db_path, queries, searchMech, bitscore,
    maximum, plusminus, labels
    ):
    with open(format_path(db_path), 'rb') as raw:
        db5 = hashlib.md5(raw.read()).hexdigest()
    log_dict = {
        'db_path':  format_path(db_path),
        'db':   db5,
        'queries':  ','.join(sorted(queries)),
        'search':   searchMech,
        'bitscore': bitscore,
        'clusMax':  maximum,
        'plusminus':    plusminus,
        'labels': labels
        }
    return log_dict

def parseLog(logPath, newLog, out_dir):
    wrk_dir = out_dir + 'working/'
    try:
        oldLog = read_json(logPath)
    except FileNotFoundError:
        oldLog = None

    if oldLog:
        try:
            with open(oldLog['db_path'], 'rb') as raw:
                db5 = hashlib.md5(raw.read()).hexdigest()
            if db5 != oldLog['db']:
                if os.path.isdir(out_dir):
                    shutil.rmtree(out_dir)
                return
            elif oldLog['search'] != newLog['search']:
                if os.path.isdir(out_dir):
                    shutil.rmtree(out_dir)
                return
            elif oldLog['bitscore'] != newLog['bitscore']:
                fas = collect_files(wrk_dir, 'fa')
                for fa in fas:
                    os.remove(fa)
                fas = collect_files(clus_dir, 'fa')
                for fa in fas:
                    os.remove(fa)
                if os.path.isdir(wrk_dir + 'tree/'):
                    shutil.rmtree(wrk_dir + 'tree/')
            elif oldLog['plusminus'] != newLog['plusminus']:
                if os.path.isdir(wrk_dir + 'genes/'):
                    shutil.rmtree(wrk_dir + 'genes/')
                if os.path.isdir(wrk_dir + 'svg/'):
                    shutil.rmtree(wrk_dir + 'svg/')
            elif oldLog['labels'] != newLog['labels']:
                if os.path.isdir(wrk_dir + 'svg/'):
                    shutil.rmtree(wrk_dir + 'svg/')
        except KeyError:
            eprint(
                '\tERROR: log file corrupted. Hoping for the best.', 
                flush = True
                )
    write_json(newLog, logPath)


def crapMngr(
    db, query, queryHits, out_dir, wrk_dir, tre_dir, fast, 
    treeSuffix, genes2query, plusminus, query2color, 
    cpus = 1, verbose = False, og = False, ogs = None, reoutput = True,
    outKeys = [], labels = True
    ):

    db = db.set_index('ome')
    info = treeMngr(
        query, out_dir, wrk_dir, tre_dir, fast, treeSuffix, cpus, verbose, reoutput
        )

    if info:
        return

    print('\t\tExtracting loci and generating synteny diagrams', flush = True)
    extractLoci_cmds = []
    if og:
        for ome, ome_genes2og in genes2query.items():
            omeHits = [x for x in queryHits if x.startswith(ome + '_')]
            if omeHits:
                extractLoci_cmds.append([db[ome]['gff3'], ome, omeHits, ogs, ome_genes2og, plusminus, query2color, wrk_dir, labels])
        with mp.Pool(processes = cpus) as pool:
            pool.starmap(extractLocusOG, extractLoci_cmds)
    else:
        for ome, ome_genes2query in genes2query.items():
            omeHits = [x for x in queryHits if x.startswith(ome + '_')]
            if omeHits:
                extractLoci_cmds.append([db[ome]['gff3'], ome, omeHits, ome_genes2query, plusminus, query2color, wrk_dir, labels])
        with mp.Pool(processes = cpus) as pool:
            pool.starmap(extractLocusGene, extractLoci_cmds)
    
    print('\t\tMapping synteny diagrams on phylogeny', flush = True)
    tree_file = tre_dir + query + treeSuffix
    with open(tree_file, 'r') as raw:
        raw_tree = raw.read()
    name_tree = ome2name(db, raw_tree, True, True, True, True, False)

    if outKeys:
        rootKey = ome2name(
            db, random.choice(outKeys), True, True, True, True, False
            )
        svgs2tree(
            query, None, name_tree, out_dir, rootKey #svg_dir, out_dir
            )
    else:
        svgs2tree(
            query, None, name_tree, out_dir #svg_dir, out_dir
            )

def OGmain(
    db, inputGenes, ogtag, fast = True, out_dir = None,
    minid = 0.05, clusParam = 0.65, minseq = 3, max_size = 250, cpus = 1,
    plusminus = 5, verbose = False, reoutput = True, interval = 0.1,
    outgroups = True, labels = True
    ):
    '''inputGenes is a list of genes within an inputted cluster'''

    wrk_dir = out_dir + 'working/'
    gff_dir, tre_dir = wrk_dir + 'genes/', wrk_dir + 'trees/'

    print('\nCompiling orthogroup data', flush = True)
    print('\tCompiling orthogroups', flush = True)
    ogInfo_dict = mycodbOGs(omes = set(db['ome']))
    og2gene, gene2og = extractOGs(ogInfo_dict, ogtag)
    inputOGs = inputGenes2inputOGs(inputGenes, gene2og, ogtag)
    inputOG2gene = {v: k for k, v in inputOGs.items()} # create hashes for transitioning

    todel, hits = [], set()
    for i, og in enumerate(inputOGs):
        if og in hits:
            todel.append(i)
        else:
            hits.add(og)
    for i in reversed(todel):
        del inputOGs[i]
        del inputGenes[i]

    og2color, conversion_dict = makeColorPalette(inputOGs)

    # in the future, genes without OGs will be placed into OGs via RBH
    if not inputOGs:
        eprint('\nERROR: no OGs for any inputted genes', flush = True)
        sys.exit(3)

    print('\tPreparing orthogroup fastas', flush = True)
    og_fas = {}
    compileOGfa_cmds = [
        [db, og2gene[og], gene] for gene, og in inputOGs.items() if not os.path.isfile(wrk_dir + gene + '.fa')
        ]
    with mp.Pool(processes = cpus) as pool:
        og_fas = {x[0]: x[1] for x in pool.starmap(compileOGfa, compileOGfa_cmds)}
    for gene, og in inputOGs.items():
        if os.path.isfile(wrk_dir + gene + '.fa'): # add finished in working directory back
            og_fas = {**og_fas, **{gene: fa2dict(wrk_dir + gene + '.fa')}}


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
        outKeys = None
        queryHits = list(queryFa.keys())
        print('\tQuery: ' + str(query), flush = True)
        if outgroups:
            print('\t\tOutgroup detection', flush = True)
            if os.path.isfile(clus_dir + query + '.fa'):
                if not os.path.isfile(wrk_dir + query + '.outgroup.fa'):
                    res = outgroupMngr(
                        db, query, minseq, max_size, clus_dir,
                        minid = minid, cpus = cpus, interval = interval,
                        verbose = False
                        )
                    if not res:
                        print('\t\t\tQuery failed, query not in ' \
                             + 'diamond output.', flush = True)
                        continue
                inKeys = set(queryHits)
                query = query + '.outgroup'
                queryHits = list(fa2dict(wrk_dir + query + '.fa'))
                outKeys = list(set(queryHits).difference(inKeys))
            else:
                eprint(
                    '\t\t\tWARNING: search results fewer than --maxseq; no outgroup detection',
                    flush = True
                    )            
        OG = inputOGs[re.sub(r'\.outgroup$','',query)] # bulletproof against outgroups
        crapMngr(
            db, query, queryHits, out_dir, wrk_dir, tre_dir, fast, treeSuffix, 
            omeGene2og, plusminus, og2color, cpus, verbose, OG, list(inputOGs.values()),
            reoutput = reoutput, outKeys = outKeys, labels = labels
            )

    for query in fas4clus:
        print('\tQuery: ' + str(query), flush = True)
        print('\t\tHierarchical agglomerative clustering', flush = True)
        outKeys = None
        res = runFA2clus(
            clus_dir + query + '.fa', db, query, minseq, max_size, clus_dir,
            query, minid, clusParam, cpus = cpus, verbose = verbose, interval = interval
            )
        if not res:
            print('\t\t\tQuery failed, sequence not in diamond matrix',
                  flush = True)
            continue
        if outgroups:
            print('\t\tOutgroup detection', flush = True)
            if not os.path.isfile(wrk_dir + query + '.outgroup.fa'):
                res = outgroupMngr(
                    db, query, minseq, max_size, clus_dir,
                    minid = minid, cpus = cpus,
                    interval = interval, verbose = False
                    )
                if not res:
                    print('\t\t\tQuery failed, query not in ' \
                         + 'diamond output.', flush = True)
                    continue
            query = query + '.outgroup'
            queryFa = fa2dict(wrk_dir + query + '.fa')
            queryHits = list(queryFa.keys())
            inKeys = set(fa2dict(wrk_dir + query + '.fa').keys())
            outKeys = list(set(queryHits).difference(inKeys))
        else:
            queryFa = fa2dict(wrk_dir + query + '.fa')
            queryHits = list(queryFa.keys())
        OG = inputOGs[re.sub(r'\.outgroup$','',query)] # bulletproof against outgroups
        crapMngr(
            db, query, queryHits, out_dir, wrk_dir, tre_dir, fast, treeSuffix,
            omeGene2og, plusminus, og2color, cpus, verbose, OG, 
            list(inputOGs.values()), reoutput = reoutput, labels = labels,
            outKeys = outKeys
            )


def SearchMain(
    db, inputGenes, queryFa, queryGff, binary = 'mmseqs', fast = True, out_dir = None,
    minid = 0.05, clusParam = 0.65, minseq = 3, max_size = 250, cpus = 1, reoutput = True,
    plusminus = 5, evalue = None, bitscore = 40, pident = 0, mem = None, verbose = False,
    interval = 0.01, outgroups = False, conversion_dict = {}, labels = True
    ):
    '''inputGenes is a list of genes within an inputted cluster'''

    print('\nPreparing run', flush = True)
    wrk_dir = out_dir + 'working/'
    gff_dir, tre_dir = wrk_dir + 'genes/', wrk_dir + 'trees/'
    clus_dir = wrk_dir + 'clus/'

    query2color, conversion_dict = makeColorPalette(inputGenes, conversion_dict)

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
            rogueLocus(query, cleanGff, wrk_dir, query2color, labels = labels)

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
        if os.path.isfile(clus_dir + query + '.fa'):
            search_fas[query] = fa2dict(clus_dir + query + '.fa')
            search_fas[query][query] = queryFa[query]
        elif os.path.isfile(wrk_dir + query + '.fa'):
            search_fas[query] = fa2dict(wrk_dir + query + '.fa')
            search_fas[query][query] = queryFa[query]

    skips = list(search_fas.keys())
    omes = set(db['ome'])
    if not len(search_fas) == len(queryFa):
        search_fas = {**search_fas, **db2search(
            db, binary, query_path, wrk_dir, evalue = evalue, bitscore = bitscore,
            pident = pident, mem = mem, biotype = 'proteome', force = True,
            skip = skips
            )}
        for query in search_fas:
            search_fas[query][query] = queryFa[query]
    
    print('\nChecking hit fasta sizes', flush = True)
    genes2query, merges = compileGenesByOme4queries(
        search_fas,
        conversion_dict, 
        set(db['ome'])
        )
    query2color = mergeColorPalette(merges, query2color)

    for query in search_fas: # revert back to other fas
        if os.path.isfile(wrk_dir + query + '.fa'):
            search_fas[query] = fa2dict(wrk_dir + query + '.fa')
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
    for query, queryFa in fas4trees.items():
        outKeys = None
        queryHits = list(queryFa.keys())
        print('\tQuery: ' + str(query), flush = True)
        if outgroups:
            print('\t\tOutgroup detection', flush = True)
            if os.path.isfile(clus_dir + query + '.fa'):
                if not os.path.isfile(wrk_dir + query + '.outgroup.fa'):
                    res = outgroupMngr(
                        db, query, minseq, max_size, clus_dir,
                        minid = minid, cpus = cpus, 
                        interval = interval, verbose = False
                        )
                    if not res:
                        print('\t\t\tQuery failed, query not in ' \
                             + 'diamond output.', flush = True)
                        continue
                inKeys = set(queryHits)
                query = query + '.outgroup'
                queryHits = list(fa2dict(wrk_dir + query + '.fa'))
                outKeys = list(set(queryHits).difference(inKeys))
            else:
                eprint(
                    '\t\t\tWARNING: search results fewer than --maxseq; no outgroup detection', 
                    flush = True
                    )
        crapMngr(
            db, query, queryHits, out_dir, wrk_dir, tre_dir, 
            fast, treeSuffix, genes2query, plusminus, query2color, 
            cpus = cpus, verbose = verbose, reoutput = reoutput,
            outKeys = outKeys, labels = labels
            )

    if fas4clus:
        if cpus > 5:
            clusCpus = 5
        else:
            clusCpus = cpus

    for query in fas4clus:
        print('\tQuery: ' + str(query), flush = True)
        print('\t\tHierarchical agglomerative clustering', flush = True)
        outKeys = None
        res = runFA2clus(
            clus_dir + str(query) + '.fa', 
            db, query, minseq, max_size, 
            clus_dir, query, minid, clusParam, cpus = clusCpus, verbose = verbose, interval = interval
            )
        if not res:
             print('\t\t\tQuery failed, sequence not in diamond matrix',
                  flush = True)
             continue
        if outgroups:
            print('\t\tOutgroup detection', flush = True)
            if not os.path.isfile(wrk_dir + query + '.outgroup.fa'):
                res = outgroupMngr(
                    db, query, minseq, max_size, clus_dir,
                    minid = minid, cpus = cpus, 
                    interval = interval, verbose = False
                    )
                if not res:
                    print('\t\t\tQuery failed, query not in ' \
                         + 'diamond output.', flush = True)
                    continue

            query = query + '.outgroup'
            queryFa = fa2dict(wrk_dir + query + '.fa')
            queryHits = list(queryFa.keys())
            inKeys = set(fa2dict(wrk_dir + query + '.fa').keys())
            outKeys = list(set(queryHits).difference(inKeys))
        else:
            queryFa = fa2dict(wrk_dir + query + '.fa')
            queryHits = list(queryFa.keys())
        crapMngr(
            db, query, queryHits, out_dir, wrk_dir, tre_dir, 
            fast, treeSuffix, genes2query, plusminus, query2color, 
            cpus = cpus, verbose = verbose, reoutput = reoutput,
            outKeys = outKeys, labels = labels
            )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = 'Mycotools integrated Cluster Reconstruction and Phylogeny (CRAP) pipeline'
        )
    parser.add_argument(
        '-q', '--query', 
        help = 'Fasta or white-space delimited file/string of cluster genes',
        required = True
        )
    parser.add_argument('-d', '--database', default = masterDB())
    parser.add_argument(
        '-s', '--search',
        help = 'Search binary {mmseqs, blastp} for search-based CRAP'
        )
    parser.add_argument(
        '-b', '--bitscore', type = float,
        default = 30, help = 'Bitscore minimum for search algorithm. DEFAULT: 30'
        )
    parser.add_argument(
        '-og', '--orthogroups', 
        help = 'MycotoolsDB Orthogroup tag for OG-based CRAP. DEFAULT: P for phylum',
        default = 'P'
        )
    parser.add_argument(
        '-p', '--plusminus',
        help = 'Genes up-/downstream to analyze from loci. DEFAULT: 10',
        default = 10, type = int
        )
    parser.add_argument(
        '--maxseq', 
        help = 'Max sequences for trees/min for fa2clus.py. ' + \
        'Outgroup detection may exceed this number. DEFAULT: 250', 
        default = 250, type = int
        )
#    parser.add_argument(
 #       '--minid', 
  #      help = 'Minimum identity for fa2clus.py. DEFAULT: 0.2', 
   #     default = 0.2, type = float
    #    )
    parser.add_argument('-i', '--iqtree', action = 'store_true', 
        help = 'IQTree2 1000 bootstrap iterations. DEFAULT: fasttree')
    parser.add_argument(
        '--conversion',
        help = 'Tab delimited conversion file for annotations: Query\tConversion'
        )
#    parser.add_argument(
 #       '-r', '--reoutput', action = 'store_true',
  #      help = 'Reoutput - permits replacing tree file, e.g. w/a rooted version'
   #     )
#    parser.add_argument(
 #       '--outgroup',
  #      help = 'Stepback agglomerative clustering on query to include outgroups'
   #     )
    parser.add_argument(
        '--ingroup',
        help = 'Do not detect outgroups, do not root trees', 
        action = 'store_true'
        )
    parser.add_argument(
        '--no_label',
        help = 'Do not label synteny diagrams',
        action = 'store_true'
        )
    parser.add_argument(
        '-g', '--gff',
        help = 'GFF for non-mycotools locus diagram. Requires -s and a fasta for -i'
        )
    parser.add_argument(
        '-o', '--output', 
        help = 'Output base dir. Will rerun if previous director exists.'
        )
    parser.add_argument('-c', '--cpu', default = 1, type = int)
    parser.add_argument('-v', '--verbose', default = False, action = 'store_true')
    args = parser.parse_args()

    execs = ['diamond', 'clipkit', 'mafft', 'iqtree']
    if args.search:
        if args.search not in {'mmseqs', 'blastp'}:
            eprint('\nERROR: invalid -b', flush = True)
            sys.exit(3)
        else:
            execs.append(args.search)
            args.orthogroups = None
    findExecs(execs, exit = set(execs))

    inputFa, inputGFF = False, False
    if os.path.isfile(args.query):
        if args.query.lower().endswith((
            '.fasta', '.fa', '.faa', '.fna', '.fsa'
            )):
            inputFa = fa2dict(args.query)
            inputGenes = list(inputFa.keys())
            if args.gff:
                inputGFF = gff2list(format_path(args.gff))
        else:
            with open(args.query, 'r') as raw:
                data = raw.read()
            inputGenes = data.rstrip().split()
    elif "'" in args.query or '"' in args.query:
        if args.gff:
            eprint('\nERROR: GFF input requires fasta input', flush = True)
            sys.exit(3)
        inputGenes = args.query.replace('"','').replace("'",'').split()
    else:
        eprint('\nERROR: invalid input', flush = True)
        sys.exit(2)


    if round(len(inputGenes)/2) > args.plusminus:
        eprint('\nERROR: -p is less than input cluster size', flush = True)
        sys.exit(1)

    if args.iqtree:
        tree = 'iqtree'
        fast = False
    else:
        tree = 'fasttree'
        fast = True
    output = format_path(args.output)
    args_dict = {
        'Input': ','.join(inputGenes), 
        'Database': args.database,
        'Orthogroup Tag': args.orthogroups,
        'Search binary': args.search,
        'Locus +/-': args.plusminus,
        'Tree': tree,
        'Maximum seq': args.maxseq,
        'Bitscore': args.bitscore,
        'GFF': args.gff,
        'Conversion file': args.conversion,
        'Labels': not args.no_label,
        'CPU': args.cpu,
        'Output directory': args.output,
        'Verbose': args.verbose
        }

    start_time = intro('Mycotools OHCrap', args_dict, 'Jason Slot, Zachary Konkel')

    db = mtdb(args.database)
    gene0 = inputGenes[0]
    ome = inputGenes[0][:inputGenes[0].find('_')]
    if not ome in set(db['ome']):
        print('\nDetected non-mycotools input', flush = True)
        if not inputFa or not args.search:
            eprint('\tnon-mycotools input requires -s, -i as a fasta, optionally -g', flush = True)
            sys.exit(4)

    if args.conversion:
        conversion_dict = parseConversion(format_path(args.conversion))
    else:
        conversion_dict = {}

    if args.orthogroups:
        newLog = initLog(
            args.database, inputGenes, 'orthogroups', args.bitscore,
            args.maxseq, args.plusminus, not args.no_label
            )
        print('\nPreparing output directory', flush = True)
        out_dir, wrk_dir, gff_dir, tre_dir = makeOutput(output, newLog)
        OGmain(
            db, inputGenes, args.orthogroups, fast = fast, 
            out_dir = out_dir, 
            minid = 0.05, clusParam = 0.65, minseq = 3, max_size = args.maxseq, cpus = args.cpu,
            verbose = args.verbose, plusminus = args.plusminus, interval = 0.1, labels = not args.no_label,
            outgroups = not args.ingroup
            )
    else:
        newLog = initLog(
            args.database, inputGenes, args.search, args.bitscore,
            args.maxseq, args.plusminus, not args.no_label
            )
        print('\nPreparing output directory', flush = True)
        out_dir, wrk_dir, gff_dir, tre_dir = makeOutput(output, newLog)
        SearchMain(
            db, inputGenes, inputFa, inputGFF, binary = args.search, fast =
            fast, out_dir = out_dir,
            minid = 0.05, clusParam = 0.65, minseq = 3, max_size = args.maxseq, cpus = 1,
            plusminus = args.plusminus, bitscore = args.bitscore, pident = 0,
            mem = None, verbose = args.verbose, interval = 0.1, 
            outgroups = not args.ingroup, conversion_dict = conversion_dict, labels = not args.no_label
            )
    outro(start_time)
