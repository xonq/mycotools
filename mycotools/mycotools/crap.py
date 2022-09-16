#! /usr/bin/env python3


# NEED to revert to aggclus when failing linclust
    # provide skip option in arguments
    # provide run aggclus in arguments
        # provide minimum connectivity fall-back option?
# NEED run input order option
# NEED to run all-v-all blast on failed diamonds
# NEED intelligent resume
    # edit gffs if conversion dict is added
# NEED to check for multiple instances of a single orthogroup/query, and use all genes as focal
# NEED to update og_main with conversion function
# NEED root to be based off the furthest in outgroup % ID from focal gene
    # dont know if possible since diamond doesnt calc all pairwise if low enough
# NEED assembly reference method, i.e. tblastn
# NEED aggclus option
# NEED to fix why rerunning can switch outgroup to midpoint

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
from mycotools.fa2clus import write_data, ClusteringError, RefineError, main as fa2clus
from mycotools.fa2tree import main as fa2tree, PhyloError
from mycotools.acc2locus import main as acc2locus
from mycotools.gff2svg import main as gff2svg
from mycotools.db2search import main as db2search
from mycotools.ome2name import main as ome2name
from mycotools.utils.og2mycodb import mycodbOGs, extract_ogs
os.environ['QT_QPA_PLATFORM'] = 'offscreen'

def parse_conversion_file(conversion_file):
    with open(conversion_file, 'r') as raw:
        data = [x.rstrip().split('\t') for x in raw]
    return {x[0]: x[1] for x in data}

def prep_gff(gff, prots, comps, hits = set(), par_dict = {}):
    RNA = False
    for entry in gff:
        if ';SearchQuery=' in entry['attributes']:
            continue
        if entry['type'] == 'gene':
            geneID = re.search(comps['id'], entry['attributes'])[1]
            if geneID in prots:
                entry['attributes'] += ';SearchQuery=' + geneID
                hits.add(geneID)
            elif geneID in par_dict:
                entry['attributes'] += ';SearchQuery=' + par_dict[geneID]
                del par_dict[geneID]
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
                    par_dict[parID] = protID
                    hits.add(protID)
                    continue
            except TypeError:
                pass
            if parID in prots:
                entry['attributes'] += ';SearchQuery=' + parID
                hits.add(parID)
            elif rnaID in par_dict:
                entry['attributes'] += ';SearchQuery=' + par_dict[rnaID]
                par_dict[parID] = par_dict[rnaID]
                del par_dict[rnaID]
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
                par_dict[parID] = protID
                hits.add(protID)
    return par_dict, hits, RNA, gff

def rogue_locus(locusID, rnaGFF, wrk_dir, query2color, labels = True):
    with open(wrk_dir + 'genes/' + locusID + '.locus.genes', 'w') as out:
            out.write(list2gff(rnaGFF))
    svg_path = wrk_dir + 'svg/' + locusID + '.locus.svg'
    gff2svg(
        rnaGFF, svg_path, product_dict = query2color, labels = labels,
        prod_comp = r';SearchQuery=([^;]+$)', width = 10, null = 'na'
        )


def input_genes2input_ogs(input_genes, gene2og, ogtag):
    input_ogs = {}
    for gene in input_genes:
        try:
            input_ogs[gene] = gene2og[gene]
        except KeyError:
            eprint('ERROR: ' + gene + ' no valid OG under tag ' + ogtag, flush = True)
            eprint('\t' + gene + ' will be ignored. Future updates will place the gene in an OG', flush = True)
    return input_ogs

def compile_og_fa(db, gene_list, query):
    fa_dict = acc2fa(db, gene_list)
    return query, fa_dict


def check_fa_size(fas, max_size):
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

def write_seq_clus(clusters, focal_gene, db, output_path, out_fa):
    cluster_dict = {}
    for gene, index in clusters.items():
        if index not in cluster_dict:
            cluster_dict[index] = []
        cluster_dict[index].append(gene)

    geneModule = cluster_dict[clusters[focal_gene]]
    try:
        fa_dict = acc2fa(db, geneModule)
    except KeyError:
        eprint('\t\t\tWARNING: some hits not in database', flush = True)
        db_omes = set(db['ome'])
        geneModule = [x for x in geneModule if x[:x.find('_')] in db_omes]
        fa_dict = acc2fa(db, geneModule)

# need to implement some method to choose if the max and min parameters couldn't be met
    print('\t\t\t' + str(len(fa_dict)) + ' genes in group', flush = True)
    with open(out_fa, 'w') as out:
        out.write(dict2fa(fa_dict))
#    write_data(None, distanceMatrix, clusters, output_path)
    write_data(None, None, clusters, output_path)


def run_fa2clus(
    fa_path, db, focal_gene, minseq, maxseq, output, out_name,
    mincon = 0.05, clus_param = 0.3, direction = -1, cpus = 1,
    verbose = False, interval = 0.1, spacer = '\t\t\t', error = False,
    algorithm = 'linclust'
    ):

    dmnd_dir = output + 'dmnd/'
    if not os.path.isdir(dmnd_dir):
        os.mkdir(dmnd_dir)
    output_path = output + str(focal_gene) 
    log_path = output + '.' + str(focal_gene) + '.log'
    try:
#        tree, clusters, distanceMatrix, ot, ocl, overshot = fa2clus(
 #           fa_path, mincon, clus_param, minseq, maxseq,
  #          search_program = 'diamond', linkage = 'single', cpus = cpus,
   #         iterative = focal_gene, interval = interval, output = output_path,
    #        verbose = verbose, log_path = log_path, refine = True, 
     #       dmnd_dir = dmnd_dir, spacer = spacer
      #      )

        null, clusters, null, null, ocl, overshot = fa2clus(
            fa_path, mincon, clus_param, minseq, maxseq,
            search_program = algorithm, iterative = focal_gene,
            interval = interval, output = output_path,
            verbose = verbose, log_path = log_path, refine = True,
            spacer = spacer, cpus = cpus
            ) # mmseqs
    except KeyError: # query not in adjacency matrix
        return False, False
    except ClusteringError as le: # if cluster error, try aggclus
        if not error:
            eprint(spacer + 'mmseqs failed, attempting hierarchical ' + \
                   'clustering', flush = True)
            try:
                null, clusters, null, null, ocl, overshot = fa2clus(
                    fa_path, 0.2, 0.7, minseq, maxseq, search_program = 'diamond',
                    iterative = focal_gene, interval = interval,
                    output = output_path, verbose = verbose,
                    log_path = log_path, refine = True, spacer = spacer,
                    cpus = cpus
                    ) # aggclus will have a higher minimum connectivity
            except KeyError:
                return False, False
        else:
            raise le

    write_seq_clus(clusters, focal_gene, db, 
                   output_path, output + '../' + str(out_name) + '.fa')

    return True, overshot

def outgroup_mngr(
    db, focal_gene, minseq, maxseq, clus_dir,
    mincon = 0.05, direction = -1, cpus = 1, interval = 0.1,
    verbose = False, spacer = '\t\t\t\t', error = False
    ):

    fa_path = clus_dir + focal_gene + '.fa'
    out_name = str(focal_gene) + '.outgroup'
    output_path = output + str(focal_gene) 
    log_path = clus_dir + '.' + str(focal_gene) + '.log'
    clusLog = read_json(log_path)
    algorithm = clusLog['search_program'] # use previous search algorithm
    iterations = sorted(
        clusLog['iterations'], key = lambda x: x['size'], reverse = True
        )
    prev_size = len(fa2dict(clus_dir + '../' + str(focal_gene) + '.fa'))
    try:
        prev_index = [int(v['size']) for i,v in enumerate(iterations)].index(prev_size)
    except ValueError: # some hits were not in the database (re: run_fa2clus)
        prev_sizes = {int(v['size']): i for i,v in enumerate(iterations)}
        hit_sizes = {k: v for k, v in prev_sizes.items() \
                     if k > minseq and k < maxseq}
        if not hit_sizes: # nothing in bound
            hit_sizes = {k: v for k, v in prev_sizes.items() \
                         if k > minseq}
            if not hit_sizes: # nothing over minimum, skip
                return False, False
        prev_index_1 = max(hit_sizes.keys())
        prev_index = hit_sizes[prev_index_1]

    new_info = iterations[prev_index-1]

    # we are trying to make it smaller than the previous size because the
    # previous size is refined to the maximum number of sequences
    # so we need to find the minimum cluster_parameter that yields the maximum
    # sequence so that we can decrease efficiently
    if prev_size > maxseq:
        prev_size = maxseq
    if prev_index > 0: # an index exists
        try:
            res, overshot = run_fa2clus(
                fa_path, db, focal_gene, minseq, prev_size - 1, clus_dir, 
                out_name, mincon, float(new_info['cluster_parameter']), 
                direction, cpus, verbose, 0.01, spacer = spacer, 
                error = error, algorithm = algorithm
                )
        except RefineError: # adjusting percent ID doesnt have reliable
        # direction in terms of cluster size
            passing_iters = [x for x in iterations \
                             if x['size'] < maxseq and x['size'] > minseq]
            sizes = sorted(set([x['size'] for x in passing_iters]),
                           reverse = True)[:2]
            if len(sizes) > 1:
                out_iter = passing_iters[0]['cluster_parameter']
                for i in passing_iters[1:]:
                    if i['size'] == sizes[1]:
                        in_iter = i['cluster_parameter']
                        break
                null, in_clusters, null, null, ocl, overshot = fa2clus(
                    fa_path, mincon, in_iter, minseq, maxseq,
                    search_program = algorithm,
                    output = output_path,
                    verbose = verbose, log_path = log_path,
                    spacer = spacer, cpus = cpus
                    )
                write_seq_clus(in_clusters, focal_gene, db, 
                               clus_dir + str(focal_gene), 
                               clus_dir + '../' + focal_gene + '.fa')
                null, out_clusters, null, null, ocl, overshot = fa2clus(
                    fa_path, mincon, out_iter, minseq, maxseq,
                    search_program = algorithm,
                    output = output_path,
                    verbose = verbose, log_path = log_path,
                    spacer = spacer, cpus = cpus
                    )
                write_seq_clus(out_clusters, focal_gene, db, 
                               clus_dir + str(focal_gene) + '.outgroup', 
                               clus_dir + '../' + out_name + '.fa')
                res, overshot = True, False # mitigate downstream error check
            else:
                eprint(
                    spacer + 'WARNING: cannot find outgroups with minimum \
                    ID', flush = True
                    )
                return False, False        
    elif 1-float(iterations[prev_index-1]['cluster_parameter']) - mincon < 0:
        eprint(
            spacer + 'WARNING: cannot find outgroups with minimum ID', flush = True
            )
        return False, False # can't go any further in the clusters
    else: # can cluster further than what's been done, need to run iterative
        clus_param = float(new_info['cluster_parameter'])
        minseq = int(new_info['size'] + 1)
        res, overshot = run_fa2clus(
            fa_path, db, focal_gene, minseq, maxseq, clus_dir, out_name,
            mincon, clus_param, direction, cpus, verbose, 0.01, 
            spacer = spacer, algorithm = algorithm
            )

    return res, overshot
    

def make_output(base_dir, new_log):

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
    parseLog(logPath, new_log, output_dir)
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
    ome_gene2og = {}
    for gene, og in gene2og.items():
        ome = gene[:gene.find('_')]
        if ome not in ome_gene2og:
            ome_gene2og[ome] = {}
        ome_gene2og[ome][gene] = og
    return ome_gene2og


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

def compile_genesXome4queries(search_fas, conversion_dict, omes = set()):
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



def extract_locus_og(gff3, ome, genesTograb, ogs, ome_gene2og, plusminus, og2color, wrk_dir, labels = True):
    try:
        gff_list = gff2list(gff3)
    except FileNotFoundError:
        eprint('\t\t\tWARNING: ' + ome + ' mycotoolsdb entry without GFF3', flush = True)
        return
    genesTograb = [x for x in genesTograb if not os.path.isfile(wrk_dir + 'svg/' + x + '.locus.svg')]
    try:
        out_indices, geneGffs = acc2locus(gff_list, genesTograb, 
                                          plusminus, mycotools = True, 
                                          geneGff = True, nt = True)
    except KeyError:
        eprint('\t\t\tWARNING: ' + ome + ' incorrectly formatted GFF3', flush = True)
        return

    extractedGenes = {}    
    for locusID, genes in out_indices.items():
        startI, endI = None, None
        for i, gene in enumerate(genes):
            try:
                geneGffs[locusID][i]['attributes'] += ';OG=' + str(ome_gene2og[gene])
            except KeyError: #gene not in gene2og
                geneGffs[locusID][i]['attributes'] += ';OG=na'
                continue
            if gene == locusID or ome_gene2og[gene] in ogs:
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


def extract_locus_gene(gff3, ome, accs, gene2query, plusminus, query2color, wrk_dir, labels = True):
    try:
        gff_list = gff2list(gff3)
    except FileNotFoundError:
        eprint('\t\t\tWARNING: ' + ome + ' mycotoolsdb entry without GFF3', flush = True)
        return
    accs = [x for x in accs if not os.path.isfile(wrk_dir + 'svg/' + x + '.locus.svg')]
    try:
        out_indices, geneGffs = acc2locus(gff_list, accs, 
                                          plusminus, mycotools = True, 
                                          geneGff = True, nt = True)
    except KeyError:
        eprint('\t\t\tWARNING: ' + ome + ' incorrectly formatted GFF3', flush = True)
        return

    extractedGenes = {}
    for locusID, genes in out_indices.items():
        if os.path.isfile(wrk_dir + 'svg/' + locusID + '.locus.svg'):
            # NEED to rerun if there's a change in output parameters
            continue
        startI, endI = None, None
        for i, gene in enumerate(genes):
            if gene in gene2query:
                geneGffs[locusID][i]['attributes'] += ';SearchQuery=' + '|'.join(gene2query[gene])
                print(gene2query[gene])
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

def svgs2tree(input_gene, og, tree_data, db, tree_path,
              out_dir, root_key = None, midpoint = True):#svg_dir, out_dir):
    tree = Tree(tree_data)
    if root_key:
        tree.set_outgroup(root_key)
        adj = 'outgroup'
    elif midpoint:
        mid = tree.get_midpoint_outgroup()
        tree.set_outgroup(mid)
        adj = 'midpoint'
    else:
        adj = 'nonrooted'
    new_tree = ome2name(db, tree.write(), True, True, True, True, False)
    tree_info = re.search(r'(.*\.)([^\.]+)$', tree_path)
    with open(tree_info[1] + adj + '.' + tree_info[2], 'w') as out:
        out.write(new_tree)

    tree = Tree(new_tree)
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
        tree.render(out_dir + input_gene + '_OG' + str(og) + '.' + adj + '.svg', 
                    w=800, tree_style = ts)
    else:
        tree.render(out_dir + input_gene + '.' + adj + '.svg', 
                    w=800, tree_style = ts)

def merge_color_palette(merges, query2color):
    for merge in merges:
        query2color['|'.join(merge)] = query2color[merge[0]]
    return query2color

def make_color_palette(inputs, conversion_dict = {}):
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


def tree_mngr(
    query, out_dir, wrk_dir, tre_dir, fast, tree_suffix, 
    cpus = 1, verbose = False, reoutput = True
    ):

    query_fa_path = wrk_dir + query + '.fa'
    if os.path.isfile(tre_dir + str(query) + tree_suffix) and reoutput:
        return
    else:
        try:
            fa2tree(query_fa_path, output_dir = tre_dir, fast = fast, cpus = cpus, verbose = verbose)
        except PhyloError:
            return query


def init_log(
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

def parseLog(logPath, new_log, out_dir):
    wrk_dir = out_dir + 'working/'
    try:
        oldLog = read_json(logPath)
    except FileNotFoundError:
        oldLog = None

    rereun_search = False
    if oldLog:
        try:
            with open(oldLog['db_path'], 'rb') as raw:
                db5 = hashlib.md5(raw.read()).hexdigest()
            if db5 != oldLog['db']:
                rerun_search = True
            elif oldLog['search'] != new_log['search']:
                if os.path.isdir(out_dir):
                    shutil.rmtree(out_dir)
                return
            elif oldLog['bitscore'] != new_log['bitscore']:
                fas = collect_files(wrk_dir, 'fa')
                for fa in fas:
                    os.remove(fa)
                fas = collect_files(clus_dir, 'fa')
                for fa in fas:
                    os.remove(fa)
                if os.path.isdir(wrk_dir + 'tree/'):
                    shutil.rmtree(wrk_dir + 'tree/')
            elif oldLog['plusminus'] != new_log['plusminus']:
                if os.path.isdir(wrk_dir + 'genes/'):
                    shutil.rmtree(wrk_dir + 'genes/')
                if os.path.isdir(wrk_dir + 'svg/'):
                    shutil.rmtree(wrk_dir + 'svg/')
            elif oldLog['labels'] != new_log['labels']:
                if os.path.isdir(wrk_dir + 'svg/'):
                    shutil.rmtree(wrk_dir + 'svg/')
        except KeyError:
            eprint(
                '\tERROR: log file corrupted. Hoping for the best.', 
                flush = True
                )
    write_json(new_log, logPath)


def crap_mngr(
    db, query, query_hits, out_dir, wrk_dir, tre_dir, fast, 
    tree_suffix, genes2query, plusminus, query2color, 
    cpus = 1, verbose = False, og = False, ogs = None, reoutput = True,
    out_keys = [], labels = True, midpoint = True
    ):

    db = db.set_index('ome')
    info = tree_mngr(
        query, out_dir, wrk_dir, tre_dir, fast, 
        tree_suffix, cpus, verbose, reoutput
        )

    if info:
        return

    print('\t\tExtracting loci and generating synteny diagrams', flush = True)
    extract_loci_cmds = []
    if og:
        for ome, ome_genes2og in genes2query.items():
            ome_hits = [x for x in query_hits if x.startswith(ome + '_')]
            if ome_hits:
                extract_loci_cmds.append([db[ome]['gff3'], ome, ome_hits, ogs, ome_genes2og, 
                                          plusminus, query2color, wrk_dir, labels])
        with mp.Pool(processes = cpus) as pool:
            pool.starmap(extract_locus_og, extract_loci_cmds)
    else:
        for ome, ome_genes2query in genes2query.items():
            ome_hits = [x for x in query_hits if x.startswith(ome + '_')]
            if ome_hits:
                extract_loci_cmds.append([db[ome]['gff3'], ome, ome_hits, ome_genes2query, 
                                          plusminus, query2color, wrk_dir, labels])
        with mp.Pool(processes = cpus) as pool:
            pool.starmap(extract_locus_gene, extract_loci_cmds)
    
    print('\t\tMapping synteny diagrams on phylogeny', flush = True)
    tree_file = tre_dir + query + tree_suffix
    with open(tree_file, 'r') as raw:
        raw_tree = raw.read()
#    name_tree = ome2name(db, raw_tree, True, True, True, True, False)

    if out_keys:
#        root_key = ome2name(
 #           db, random.choice(out_keys), True, True, True, True, False
  #          )
        root_key = random.choice(out_keys)
        svgs2tree(
            query, None, raw_tree, db, tree_file,
            out_dir, root_key, midpoint = midpoint #svg_dir, out_dir
            )
    else:
        svgs2tree(
            query, None, raw_tree, db, tree_file,
            out_dir, midpoint = midpoint #svg_dir, out_dir
            )

def og_main(
    db, input_genes, ogtag, fast = True, out_dir = None,
    mincon = 0.05, clus_param = 0.65, minseq = 3, max_size = 250, cpus = 1,
    plusminus = 10000, verbose = False, reoutput = True, interval = 0.1,
    outgroups = True, labels = True, midpoint = True, 
    clus_meth = 'mmseqs easy-cluster'
    ):
    """input_genes is a list of genes within an inputted cluster"""

    wrk_dir = out_dir + 'working/'
    gff_dir, tre_dir = wrk_dir + 'genes/', wrk_dir + 'trees/'

    print('\nCompiling orthogroup data', flush = True)
    print('\tCompiling orthogroups', flush = True)
    og_info_dict = mycodbOGs(omes = set(db['ome']))
    og2gene, gene2og = extract_ogs(og_info_dict, ogtag)
    input_ogs = input_genes2input_ogs(input_genes, gene2og, ogtag)
    input_og2gene = {v: k for k, v in input_ogs.items()} # create hashes for transitioning

    todel, hits = [], set()
    for i, og in enumerate(input_ogs):
        if og in hits:
            todel.append(i)
        else:
            hits.add(og)
    for i in reversed(todel):
        del input_ogs[i]
        del input_genes[i]

    og2color, conversion_dict = make_color_palette(input_ogs)

    # in the future, genes without OGs will be placed into OGs via RBH
    if not input_ogs:
        eprint('\nERROR: no OGs for any inputted genes', flush = True)
        sys.exit(3)

    print('\tPreparing orthogroup fastas', flush = True)
    og_fas = {}
    compile_og_fa_cmds = [
        [db, og2gene[og], gene] for gene, og in input_ogs.items() if not os.path.isfile(wrk_dir + gene + '.fa')
        ]
    with mp.Pool(processes = cpus) as pool:
        og_fas = {x[0]: x[1] for x in pool.starmap(compile_og_fa, compile_og_fa_cmds)}
    for gene, og in input_ogs.items():
        if os.path.isfile(wrk_dir + gene + '.fa'): # add finished in working directory back
            og_fas = {**og_fas, **{gene: fa2dict(wrk_dir + gene + '.fa')}}


    print('\nChecking fasta sizes', flush = True)
    fas4clus, fas4trees = check_fa_size(og_fas, max_size)
    for query, hit_fa in fas4trees.items():
        hit_fa_path = wrk_dir + query + '.fa'
        with open(hit_fa_path, 'w') as out:
            out.write(dict2fa(hit_fa))

    if fas4clus:
        clus_dir = wrk_dir + 'clus/'
        if not os.path.isdir(clus_dir):
            os.mkdir(clus_dir)
        for query, fa in fas4clus.items():
            clus_fa_path = clus_dir + query + '.fa'
            if not os.path.isfile(clus_fa_path):
                with open(clus_fa_path, 'w') as out:
                    out.write(dict2fa(fa))
        print('\nRunning clustering on ' + str(len(fas4clus)) + ' fastas', flush = True)

    print('\nCRAP', flush = True)
    ome_gene2og = gene2og2ome2og(gene2og)
    if fast:
        tree_suffix = '.fa.clipkit.treefile'
    else:
        tree_suffix = '.fa.clipkit.contree'
    fas4trees = {k: v for k, v in sorted(fas4trees.items(), key = lambda x: len(x[1]))}
    for query, query_fa in fas4trees.items():
        out_keys = None
        query_hits = list(query_fa.keys())
        print('\tQuery: ' + str(query), flush = True)
        if outgroups:
            print('\t\tOutgroup detection', flush = True)
            if os.path.isfile(clus_dir + query + '.fa'):
                if not os.path.isfile(wrk_dir + query + '.outgroup.fa'):
                    res, overshot = outgroup_mngr(
                        db, query, minseq, max_size, clus_dir,
                        mincon = mincon, cpus = cpus, interval = interval,
                        verbose = False
                        )
                    if not res:
                        print('\t\t\tERROR: query not in ' \
                             + 'output/clusters are granular', flush = True)
                        continue
                in_keys = set(query_hits)
                out_query = query + '.outgroup'
                out_query_hits = list(fa2dict(wrk_dir + out_query + '.fa'))
                out_keys = list(in_keys.difference(set(out_query_hits)))
            else:
                eprint(
                    '\t\t\tWARNING: search results fewer than --maxseq; no outgroup detection',
                    flush = True
                    )            
        OG = input_ogs[re.sub(r'\.outgroup$','',query)] # bulletproof against outgroups
        crap_mngr(
            db, query, query_hits, out_dir, wrk_dir, tre_dir, fast, tree_suffix, 
            ome_gene2og, plusminus, og2color, cpus, verbose, OG, list(input_ogs.values()),
            reoutput = reoutput, out_keys = out_keys, labels = labels,
            midpoint = midpoint
            )

    for query in fas4clus:
        print('\tQuery: ' + str(query), flush = True)
        print('\t\tSequence clustering', flush = True)
        out_keys = None
        res, overshot = run_fa2clus(
            clus_dir + query + '.fa', db, query, minseq, max_size, clus_dir,
            query, mincon, clus_param, cpus = cpus, verbose = verbose, 
            interval = interval, algorithm = clus_meth
            )
        if not res:
            print('\t\t\tERROR: query had no significant hits',
                  flush = True)
            continue
        if outgroups and not overshot:
            print('\t\tOutgroup detection', flush = True)
            if not os.path.isfile(wrk_dir + query + '.outgroup.fa'):
                res, overshot = outgroup_mngr(
                    db, query, minseq, max_size, clus_dir,
                    mincon = mincon, cpus = cpus,
                    interval = interval, verbose = False
                    )
                if not res:
                    print('\t\t\tERROR: query not in ' \
                         + 'diamond output/clusters are granular', flush = True)
                    continue
            out_query = query + '.outgroup'
            out_query_fa = fa2dict(wrk_dir + out_query + '.fa')
            query_hits = list(query_fa.keys())
            in_keys = set(out_query_fa.keys())
            out_keys = list(set(query_hits).difference(in_keys))
        else:
            query_fa = fa2dict(wrk_dir + query + '.fa')
            query_hits = list(query_fa.keys())
        OG = input_ogs[re.sub(r'\.outgroup$','',query)] # bulletproof against outgroups
        crap_mngr(
            db, query, query_hits, out_dir, wrk_dir, tre_dir, fast, tree_suffix,
            ome_gene2og, plusminus, og2color, cpus, verbose, OG, 
            list(input_ogs.values()), reoutput = reoutput, labels = labels,
            out_keys = out_keys, midpoint = midpoint
            )


def search_main(
    db, input_genes, query_fa, query_gff, binary = 'mmseqs', fast = True, 
    out_dir = None, mincon = 0.05, clus_param = 0.65, minseq = 3, 
    max_size = 250, cpus = 1, reoutput = True, plusminus = 10000, 
    evalue = None, bitscore = 40, pident = 0, mem = None, verbose = False,
    interval = 0.01, outgroups = False, conversion_dict = {}, labels = True,
    midpoint = True, clus_meth = 'mmseqs easy-linclust'
    ):
    """input_genes is a list of genes within an inputted cluster"""

    print('\nPreparing run', flush = True)
    wrk_dir = out_dir + 'working/'
    gff_dir, tre_dir = wrk_dir + 'genes/', wrk_dir + 'trees/'
    clus_dir = wrk_dir + 'clus/'

    query2color, conversion_dict = make_color_palette(input_genes, conversion_dict)

    if query_gff:
        print('\tCleaning input GFF', flush = True)
        par_dict, prot_hits, RNA, query_gff = prep_gff(query_gff, set(input_genes), 
                                                       gff3Comps())
        count = 0
        while par_dict and count < 4:
            count += 1
            par_dict, prot_hits, RNA, query_gff = prep_gff(query_gff, set(input_genes),
                                                           gff3Comps(), prot_hits, 
                                                           par_dict)
        if par_dict:
            eprint('\nIncorrectly formatted GFF', flush = True)
            sys.exit(5)
        elif set(input_genes).difference(prot_hits):
            eprint('\nProteins missing from GFF', flush = True)
            eprint('\t' + ','.join([
                                str(x) for x in list(
                                   set(input_genes).difference(prot_hits)
                                   )
                              ])
                    )
            sys.exit(6)
        clean_gff = [x for x in query_gff if ';SearchQuery=' in x['attributes'] and 'RNA' in x['type']]
        for query in input_genes:
            rogue_locus(query, clean_gff, wrk_dir, query2color, labels = labels)

    query_path = wrk_dir + 'query.crap.fa'
    if not query_fa:
        query_fa = acc2fa(db, input_genes)
        with open(query_path, 'w') as out:
            out.write(dict2fa(query_fa))
    elif not os.path.isfile(query_path):
        with open(query_path, 'w') as out:
            out.write(dict2fa(query_fa))

    search_fas = {}
    for query in query_fa:
        if os.path.isfile(clus_dir + query + '.fa'):
            search_fas[query] = fa2dict(clus_dir + query + '.fa')
            search_fas[query][query] = query_fa[query]
        elif os.path.isfile(wrk_dir + query + '.fa'):
            search_fas[query] = fa2dict(wrk_dir + query + '.fa')
            search_fas[query][query] = query_fa[query]

    skips = list(search_fas.keys())
    omes = set(db['ome'])
    if not len(search_fas) == len(query_fa):
        if binary == 'diamond':
            binary = 'blastp'
            diamond = 'diamond'
        else:
            diamond = False
        search_fas = {**search_fas, **db2search(
            db, binary, query_path, wrk_dir, evalue = evalue, bitscore = bitscore,
            pident = pident, mem = mem, biotype = 'proteome', force = True,
            skip = skips, diamond = diamond
            )}
        for query in search_fas:
            search_fas[query][query] = query_fa[query]
    
    print('\nChecking hit fasta sizes', flush = True)
    genes2query, merges = compile_genesXome4queries(
        search_fas,
        conversion_dict, 
        set(db['ome'])
        )
    query2color = merge_color_palette(merges, query2color)

    for query in search_fas: # revert back to other fas
        if os.path.isfile(wrk_dir + query + '.fa'):
            search_fas[query] = fa2dict(wrk_dir + query + '.fa')
    fas4clus, fas4trees = check_fa_size(search_fas, max_size)
    for query, hit_fa in fas4trees.items():
        hit_fa_path = wrk_dir + query + '.fa'
        with open(hit_fa_path, 'w') as out:
            out.write(dict2fa(hit_fa))

    if fas4clus:
        if not os.path.isdir(clus_dir):
            os.mkdir(clus_dir)
        for query, fa in fas4clus.items():
            clus_fa_path = clus_dir + query + '.fa'
            if not os.path.isfile(clus_fa_path):
                with open(clus_fa_path, 'w') as out:
                    out.write(dict2fa(fa))
        print('\tRunning clustering on ' + str(len(fas4clus)) + ' fastas', flush = True)

    print('\nCRAP', flush = True)
    if fast:
        tree_suffix = '.fa.clipkit.treefile'
    else:
        tree_suffix = '.fa.clipkit.contree'
    fas4trees = {k: v for k, v in sorted(fas4trees.items(), key = lambda x: len(x[1]))}
    for query, query_fa in fas4trees.items():
        out_keys = None
        query_hits = list(query_fa.keys())
        print('\tQuery: ' + str(query), flush = True)
        if outgroups:
            print('\t\tOutgroup detection', flush = True)
            if os.path.isfile(clus_dir + query + '.fa'):
                if not os.path.isfile(wrk_dir + query + '.outgroup.fa'):
                    res, overshot = outgroup_mngr(
                        db, query, minseq, max_size, clus_dir,
                        mincon = mincon, cpus = cpus, 
                        interval = interval, verbose = False
                        )
                    if not res:
                        print('\t\t\tERROR: query not in ' \
                             + 'diamond output/clusters are granular', flush = True)
                        continue
                in_keys = set(query_hits)
                out_query = query + '.outgroup'
                out_query_hits = list(fa2dict(wrk_dir + out_query + '.fa'))
                out_keys = list(in_keys.difference(set(out_query_hits)))
            else:
                eprint(
                    '\t\t\tWARNING: search results fewer than --maxseq; no outgroup detection', 
                    flush = True
                    )
        crap_mngr(
            db, query, query_hits, out_dir, wrk_dir, tre_dir, 
            fast, tree_suffix, genes2query, plusminus, query2color, 
            cpus = cpus, verbose = verbose, reoutput = reoutput,
            out_keys = out_keys, labels = labels, midpoint = midpoint
            )

    if fas4clus:
        if cpus > 5:
            clus_cpus = 5
        else:
            clus_cpus = cpus

    for query in fas4clus:
        print('\tQuery: ' + str(query), flush = True)
        print('\t\tSequence clustering', flush = True)
        out_keys = None
        res, overshot = run_fa2clus(
            clus_dir + str(query) + '.fa', 
            db, query, minseq, max_size, 
            clus_dir, query, mincon, clus_param, cpus = clus_cpus, 
            verbose = verbose, interval = interval, algorithm = clus_meth
            )
        if not res:
             print('\t\t\tFAILED: query had no significant hits',
                  flush = True)
             continue
        if outgroups and not overshot:
            print('\t\tOutgroup detection', flush = True)
            if not os.path.isfile(wrk_dir + query + '.outgroup.fa'):
                res, overshot = outgroup_mngr(
                    db, query, minseq, max_size, clus_dir,
                    mincon = mincon, cpus = cpus, 
                    interval = interval, verbose = False
                    )
                if not res:
                    print('\t\t\tFAILED: query not in ' \
                         + 'diamond output/clusters are granular', flush = True)
                    continue

            query_hits = list(fa2dict(wrk_dir + query + '.fa').keys())
            in_keys = set(query_hits)
            out_query = query + '.outgroup'
            out_query_fa = fa2dict(wrk_dir + out_query + '.fa')
            out_query_hits = set(query_fa.keys())
            out_keys = list(in_keys.difference(out_query_hits))
        else:
            query_fa = fa2dict(wrk_dir + query + '.fa')
            query_hits = list(query_fa.keys())
        crap_mngr(
            db, query, query_hits, out_dir, wrk_dir, tre_dir, 
            fast, tree_suffix, genes2query, plusminus, query2color, 
            cpus = cpus, verbose = verbose, reoutput = reoutput,
            out_keys = out_keys, labels = labels, midpoint = midpoint
            )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = 'Mycotools integrated Cluster Reconstruction and Phylogeny (CRAP) pipeline'
        )

    i_opt = parser.add_argument_group('Inputs')
    i_opt.add_argument(
        '-q', '--query', 
        help = 'Fasta, white-space delimited file/string of cluster genes, \
                "-" for stdin',
        required = True
        )
    i_opt.add_argument('-d', '--mtdb', default = masterDB())
    i_opt.add_argument(
        '-g', '--gff',
        help = 'GFF for non-mycotools locus diagram. Requires -q fasta input'
        )


    hg_opt = parser.add_argument_group('Homolog inference')
    hg_opt.add_argument(
        '-s', '--search',
        help = 'Search binary {mmseqs, diamond, blastp} for search-based CRAP'
        )
    hg_opt.add_argument(
        '-b', '--bitscore', type = float,
        default = 30, help = 'Bitscore minimum for search algorithm. DEFAULT: 30'
        )
    hg_opt.add_argument(
        '-og', '--orthogroups', 
        help = 'MycotoolsDB Orthogroup tag for OG-based CRAP. DEFAULT: P for phylum',
        default = 'P'
        )

    clus_opt = parser.add_argument_group('Phylogeny size management')
    clus_opt.add_argument(
        '--minseq',
        help = 'Min sequences for trees/cluster size; ' \
             + 'DEFAULT: 3',
        default = 3, type = int
        )
    clus_opt.add_argument(
        '--maxseq', 
        help = 'Max sequences for trees/min for fa2clus.py; ' \
             + 'DEFAULT: 250', 
        default = 250, type = int
        )
#    clus_opt.add_argument(
 #       '--mincon', 
  #      help = 'Minimum identity for fa2clus.py. DEFAULT: 0.2', 
   #     default = 0.2, type = float
    #    )
    clus_opt.add_argument(
        '-l', '--linclust',
        help = 'Cluster large gene sets via mmseqs linclust; DEFAULT: \
                easy-cluster',
        action = 'store_true'
        )
    clus_opt.add_argument(
        '-a', '--agg_clus',
        help = 'Cluster large gene sets via hierarchical clustering ' \
             + 'NONFUNCTIONAL',
        action = 'store_true'
        ) # NOT IMPLEMENTED
    clus_opt.add_argument(
        '-f', '--fail',
        help = 'Do not fallback to alternative clustering method upon ' \
             + 'failure NONFUNCTIONAL',
        action = 'store_true'
        ) # NOT IMPLEMENTED


    phy_opt = parser.add_argument_group('Phylogeny annotation')
    phy_opt.add_argument(
        '-p', '--plusminus',
        help = 'Bases up-/downstream query hits for homolog search; DEFAULT: 20,000',
        default = 20000, type = int
        )
    phy_opt.add_argument(
        '--conversion',
        help = 'Tab delimited alternative annotation file: "query,name"'
        )
    phy_opt.add_argument('-i', '--iqtree', action = 'store_true', 
        help = 'IQTree2 1000 bootstrap iterations. DEFAULT: fasttree')
    phy_opt.add_argument(
        '--no_outgroup',
        help = 'Do not detect outgroups, do not root trees', 
        action = 'store_true'
        )
    phy_opt.add_argument(
        '--no_midpoint',
        help = 'Do not infer midpoint roots for outgroup failures',
        action = 'store_true'
        )
    phy_opt.add_argument(
        '--no_label',
        help = 'Do not label synteny diagrams',
        action = 'store_true'
        )

    run_opt = parser.add_argument_group('Runtime options')
    run_opt.add_argument(
        '-o', '--output', 
        help = 'Output base dir - will rerun if previous directory exists'
        )
    run_opt.add_argument('-c', '--cpu', default = 1, type = int)
    run_opt.add_argument('-v', '--verbose', default = False, action = 'store_true')
    args = parser.parse_args()

    execs = ['diamond', 'clipkit', 'mafft', 'iqtree']
    if args.search:
        if args.search not in {'mmseqs', 'blastp', 'diamond'}:
            eprint('\nERROR: invalid -b', flush = True)
            sys.exit(3)
        else:
            execs.append(args.search)
            args.orthogroups = None
    findExecs(execs, exit = set(execs))

    input_fa, input_GFF = False, False
    if args.query == '-':
        if args.gff:
            eprint('\nERROR: GFF input requires fasta input', flush = True)
            sys.exit(3)
        input_genes = stdin2str().replace('"','').replace("'",'').split()
    elif os.path.isfile(args.query):
        if args.query.lower().endswith((
            '.fasta', '.fa', '.faa', '.fna', '.fsa'
            )):
            input_fa = fa2dict(args.query)
            input_genes = list(input_fa.keys())
            if args.gff:
                input_GFF = gff2list(format_path(args.gff))
        else:
            with open(args.query, 'r') as raw:
                data = raw.read()
            input_genes = data.rstrip().split()
    elif "'" in args.query or '"' in args.query:
        if args.gff:
            eprint('\nERROR: GFF input requires fasta input', flush = True)
            sys.exit(3)
        input_genes = args.query.replace('"','').replace("'",'').split()
    else:
        eprint('\nERROR: invalid input', flush = True)
        sys.exit(2)

# would be nice to make this work from a base standpoint now
#    if round(len(input_genes)/2) > args.plusminus:
 #       eprint('\nERROR: -p is less than input cluster size', flush = True)
 #       sys.exit(1) 

    if args.iqtree:
        tree = 'iqtree'
        fast = False
    else:
        tree = 'fasttree'
        fast = True

    if args.agg_clus and args.linclust:
        eprint('\nERROR: multiple clustering methods specified', flush = True)
        sys.exit(5)
    elif args.agg_clus:
        clus_meth = 'hierarchical'
    elif args.linclust:
        clus_meth = 'mmseqs easy-linclust'
    else:
        clus_meth = 'mmseqs easy-cluster'
    output = format_path(args.output)
    args_dict = {
        'Input': ','.join(input_genes), 
        'Database': args.mtdb,
        'Orthogroup Tag': args.orthogroups,
        'Search binary': args.search,
        'Locus +/-': args.plusminus,
        'Cluster method': clus_meth,
        'Tree': tree,
        'Minimum seq': args.minseq,
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

    db = mtdb(args.mtdb)
    gene0 = input_genes[0]
    ome = input_genes[0][:input_genes[0].find('_')]
    if not ome in set(db['ome']):
        print('\nDetected non-mycotools input', flush = True)
        if not input_fa or not args.search:
            eprint('\tnon-mycotools input requires -s, -i as a fasta, optionally -g', flush = True)
            sys.exit(4)

    if args.conversion:
        conversion_dict = parse_conversion_file(format_path(args.conversion))
    else:
        conversion_dict = {}

    if args.orthogroups:
        new_log = init_log(
            args.mtdb, input_genes, 'orthogroups', args.bitscore,
            args.maxseq, args.plusminus, not args.no_label
            )
        print('\nPreparing output directory', flush = True)
        out_dir, wrk_dir, gff_dir, tre_dir = make_output(output, new_log)
        og_main(
            db, input_genes, args.orthogroups, fast = fast, 
            out_dir = out_dir, 
            mincon = 0.05, clus_param = 0.65, minseq = args.minseq, 
            max_size = args.maxseq, cpus = args.cpu,
            verbose = args.verbose, plusminus = args.plusminus, 
            interval = 0.1, labels = not args.no_label,
            outgroups = not args.no_outgroup, clus_meth = clus_meth,
            midpoint = not bool(args.no_midpoint)
            )
    else:
        new_log = init_log(
            args.mtdb, input_genes, args.search, args.bitscore,
            args.maxseq, args.plusminus, not args.no_label
            )
        print('\nPreparing output directory', flush = True)
        out_dir, wrk_dir, gff_dir, tre_dir = make_output(output, new_log)
        search_main(
            db, input_genes, input_fa, input_GFF, 
            binary = args.search, fast = fast, out_dir = out_dir,
            mincon = 0.05, clus_param = 0.65, minseq = args.minseq, 
            max_size = args.maxseq, cpus = 1, plusminus = args.plusminus, 
            bitscore = args.bitscore, pident = 0, mem = None, 
            verbose = args.verbose, interval = 0.1, 
            midpoint = not bool(args.no_midpoint),
            outgroups = not args.no_outgroup, conversion_dict = conversion_dict, 
            clus_meth = clus_meth, labels = not args.no_label
            )
    outro(start_time)
