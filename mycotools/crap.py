#! /usr/bin/env python3

# NEED outgroup manager for clusters that fit within min and max seq
    # crude identify the largest cluster below the total sequences,
    # then refine up
# NEED percent positive filter
# NEED provide run aggclus in arguments
        # provide minimum connectivity fall-back option?
# NEED run input order option
# NEED intelligent resume
    # edit gffs if conversion dict is added
# NEED to check for multiple instances of a single homolog/query, and use all genes as focal
# NEED to update hg_main with conversion function
# NEED more intelligent root inference from outgroup
# NEED assembly reference method, i.e. tblastn
# NEED cluster variable input default

import os
import re
import sys
import argparse
import datetime
import hashlib
import random
import multiprocessing as mp
import shutil
from collections import Counter
try:
    from ete3 import Tree, faces, TreeStyle, NodeStyle, AttrFace
    from ete3.parser.newick import NewickError
except ImportError:
    raise ImportError('Install ete3 into your conda environment via `conda install ete3`')
from mycotools.lib.dbtools import mtdb, primaryDB
from mycotools.lib.kontools import eprint, format_path, findExecs, intro, outro, \
    read_json, write_json, stdin2str, getColors, collect_files
from mycotools.lib.biotools import fa2dict, dict2fa, gff2list, list2gff, gff3Comps
from mycotools.acc2fa import dbmain as acc2fa
from mycotools.fa2clus import write_data, ClusteringError, \
    ClusterParameterError, main as fa2clus, sort_iterations
from mycotools.fa2tree import main as fa2tree, PhyloError
from mycotools.acc2locus import main as acc2locus
from mycotools.gff2svg import main as gff2svg
from mycotools.db2search import blast_main as db2search
from mycotools.ome2name import main as ome2name
#from mycotools.utils.og2mycodb import mycodbHGs, extract_ogs
from mycotools.db2microsyntree import compile_homolog_groups
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


def input_genes2input_hgs(input_genes, gene2hg):
    input_hgs = {}
    for gene in input_genes:
        try:
            input_hgs[gene] = gene2hg[gene]
        except KeyError:
            eprint('ERROR: ' + gene + ' query with no HG', flush = True)
            eprint('\t' + gene + ' will be ignored.', flush = True)
    return input_hgs

def compile_hg_fa(db, gene_list, query):
    fa_dict = acc2fa(db, gene_list)
    return query, fa_dict


def check_fa_size(fas, max_size):
    fas4clus, fas4trees = {}, {}
    for query, fa in fas.items():
        print('\t' + query + '\t' + str(len(fa)) + ' genes', flush = True)
        if len(fa) > max_size:
            fas4clus[query] = fa
        elif len(fa) < 3:
            eprint(f'\t\tWARNING: too few hits ({len(fa)})', flush = True)
        else:
            fas4trees[query] = fa

    fas4clus = {k: v for k, v in sorted(fas4clus.items(), key = lambda x: len(x[1]))}
    fas4trees = {k: v for k, v in sorted(fas4trees.items(), key = lambda x: len(x[1]))}

    return fas4clus, fas4trees

def write_seq_clus(gene_module, focal_gene, db, output_path, out_fa):
#    gene_module = cluster_dict[clusters[focal_gene]]
    try:
        fa_dict = acc2fa(db, gene_module)
    except KeyError:
        eprint('\t\t\tWARNING: some hits not in database', flush = True)
        db_omes = set(db.keys())
        gene_module = [x for x in gene_module if x[:x.find('_')] in db_omes]
        fa_dict = acc2fa(db, gene_module)

# need to implement some method to choose if the max and min parameters couldn't be met
    print('\t\t\t' + str(len(fa_dict)) + ' genes in group', flush = True)
    with open(out_fa, 'w') as out:
        out.write(dict2fa(fa_dict))
#    write_data(None, distanceMatrix, clusters, output_path)
#    write_data(None, None, clusters, output_path)


def run_fa2clus(
    fa_path, db, focal_gene, min_seq, max_seq, output, out_name,
    clus_cons = 0.4, clus_var = 0.3, min_var = 0.1, max_var = 1,
    cpus = 1, verbose = False, 
    interval = 0.1, spacer = '\t\t\t', error = False, 
    algorithm = 'linclust'
    ):

    dmnd_dir = output + 'dmnd/'
    if not os.path.isdir(dmnd_dir):
        os.mkdir(dmnd_dir)
    output_path = output + str(focal_gene) 
    log_path = output + '.' + str(focal_gene) + '.log'
    if os.path.isfile(log_path):
        fa2clus_log = read_json(log_path)
    else:
        fa2clus_log = {'algorithm': 'null'}
    try:
        if fa2clus_log['algorithm'] == 'hierarchical':
            raise ClusteringError # go straight to hierarchical
        cluster, null, overshot, fa2clus_log = fa2clus(
            fa_path, clus_cons, clus_var, min_seq, max_seq,
            search_program = algorithm, focal_gene = focal_gene,
            interval = interval, output = output_path,
            verbose = verbose, log_path = log_path,
            spacer = spacer, cpus = cpus
            ) # mmseqs
    except KeyError: # query not in adjacency matrix
        return False, False, fa2clus_log
    except ClusteringError as le: # if cluster error, try aggclus
        if not error:
            eprint(spacer + 'mmseqs failed, attempting hierarchical ' + \
                   'clustering', flush = True)
            try:
                cluster, newick, overshot, fa2clus_log = fa2clus(
                    fa_path, 0.2, 0.7, min_seq, max_seq, 
                    search_program = 'diamond',
                    focal_gene = focal_gene, interval = interval,
                    output = output_path, verbose = verbose,
                    log_path = log_path, spacer = spacer,
                    cpus = cpus, dmnd_dir = dmnd_dir
                    ) # aggclus will have a higher minimum connectivity
            except KeyError:
                return False, False, fa2clus_log
        else:
            raise le

    write_seq_clus(cluster, focal_gene, db, 
                   output_path, output + '../' + str(out_name) + '.fa')

    return True, overshot, fa2clus_log


def outgroup_mngr(
    db, focal_gene, min_seq, max_seq, clus_dir, wrk_dir,
    clus_cons = 0.4, direction = -1, cpus = 1, interval = 0.1,
    verbose = False, spacer = '\t\t\t\t', error = False
    ):
    """There needs to be a run in the fa2clus_log['successes'] for this
    function to operate, or else an IndexError will result"""

    fa_path = clus_dir + focal_gene + '.fa'
    out_name = str(focal_gene) + '.outgroup'
    output_path = clus_dir + str(focal_gene) 
    dmnd_dir = clus_dir + 'dmnd/'
    log_path = clus_dir + '.' + str(focal_gene) + '.log'
    fa2clus_log = read_json(log_path)
    algorithm = fa2clus_log['algorithm'] # use previous search algorithm
    successes = fa2clus_log['successes']
    iterations = fa2clus_log['iterations']
    prev_size = len(fa2dict(clus_dir + '../' + str(focal_gene) + '.fa'))

    max_success = successes[0]
    # it is best to have the largest cluster, so attempt to refine upward
    if max_seq - max_success['size'] > 1: # if there is room to refine upward
        run_min_seq = max_success['size']
        if algorithm == 'hierarchical':
            try:
                min_var, max_var = max_success['cluster_variable'], 0.7
                clus_var = max_success['cluster_variable'] + 0.01
                cluster, newick, overshot, fa2clus_log = fa2clus(
                    fa_path, 0.2, clus_var, run_min_seq, max_seq, 
                    search_program = 'diamond', dmnd_dir = dmnd_dir,
                    focal_gene = focal_gene, interval = interval,
                    output = output_path, verbose = verbose,
                    log_path = log_path, spacer = spacer,
                    cpus = cpus, min_var = min_var, max_var = max_var
                    ) # aggclus will have a higher minimum connectivity
            except ClusterParameterError as e: # has run previously refined all it
                pass

        else:
            min_var, max_var = 0.1, max_success['cluster_variable']
            clus_var = max_success['cluster_variable'] - 0.01
            cluster, null, overshot, fa2clus_log = fa2clus(
                fa_path, clus_cons, clus_var, run_min_seq, max_seq,
                search_program = algorithm, focal_gene = focal_gene,
                interval = interval, output = output_path,
                verbose = verbose, log_path = log_path,
                spacer = spacer, cpus = cpus, min_var = min_var,
                max_var = max_var
                ) # mmseqs
    if not any(x['size'] > max_success['size'] \
        for x in fa2clus_log['successes']):
        # was not able to refine upward
        run_max_seq = max_success['size'] - 1
        if any(x['size'] < max_success['size'] \
            for x in fa2clus_log['successes']):
            # are there two successes to refine between?
            for i, suc in enumerate(successes[1:]):
                if suc['size'] < run_max_seq:
                    run_min_seq = suc['size']
                    run_min_seq_var = suc['cluster_variable']
                    break
        else:
            run_min_seq = min_seq
            run_min_seq_var = None
        if algorithm == 'hierarchical':
            max_var = max_success['cluster_variable']
            if not run_min_seq_var:
                min_var = 0.05
            else:
                min_var = run_min_seq_var
            clus_var = max_var - 0.01
            cluster, newick, overshot, fa2clus_log = fa2clus(
                fa_path, 0.2, clus_var, run_min_seq, run_max_seq, 
                search_program = 'diamond', dmnd_dir = dmnd_dir,
                focal_gene = focal_gene, interval = interval,
                output = output_path, verbose = verbose,
                log_path = log_path, spacer = spacer,
                cpus = cpus, min_var = min_var, max_var = max_var
                ) # aggclus  (clus_cons is held at 0.2)
        else:
            min_var = max_success['cluster_variable']
            try:
                if not run_min_seq_var:
                    max_var = 1
                else:
                    max_var = run_min_seq_var
            except UnboundLocalError: # run_min_seq_var not here
                max_var = 1
            clus_var = max_var - 0.01
            try:
                cluster, null, overshot, fa2clus_log = fa2clus(
                    fa_path, clus_cons, clus_var, run_min_seq, run_max_seq,
                    search_program = algorithm, focal_gene = focal_gene,
                    interval = interval, output = output_path,
                    verbose = verbose, log_path = log_path,
                    spacer = spacer, cpus = cpus, min_var = min_var,
                    max_var = max_var
                    ) # mmseqs
            except ClusterParameterError as e:
                pass

    [fa2clus_log['successes'].append(x) for x in fa2clus_log['iterations'] \
     if x['size'] > min_seq and x['size'] < max_seq]
     # account for successes that were only acquired in the outgroup detection
    fa2clus_log['successes'] = sort_iterations(fa2clus_log['successes'])

    # set the out_group to the largest module
    out_group = fa2clus_log['successes'][0]['cluster']


    if any(x['size'] < fa2clus_log['successes'][0]['size'] \
           for x in fa2clus_log['successes'][1:]): # if two different size
           # modules exist
        for i in fa2clus_log['successes'][1:]:
            if i['size'] < len(out_group):
                in_group = i['cluster']
                break
        with open(wrk_dir + focal_gene + '.fa', 'w') as out:
            out.write(dict2fa(acc2fa(db, out_group)))
        return in_group, out_group
    else:
        return out_group, set()
    

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
    

def gene2hg2ome2hg(gene2hg):
    ome_gene2hg = {}
    for gene, og in gene2hg.items():
        ome = gene[:gene.find('_')]
        if ome not in ome_gene2hg:
            ome_gene2hg[ome] = {}
        ome_gene2hg[ome][gene] = og
    return ome_gene2hg


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



def extract_locus_hg(gff3, ome, genesTograb, ogs, ome_gene2hg, 
                     plusminus, wrk_dir, labels = True):
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
        eprint('\t\t\tWARNING: ' + ome + ' could not parse gff', flush = True)
        return

    extractedGenes, new_hgs = {}, []
    for locusID, genes in out_indices.items():
        startI, endI = None, None
        for i, gene in enumerate(genes):
            try:
                geneGffs[locusID][i]['attributes'] += ';HG=' + str(ome_gene2hg[gene])
            except KeyError: #gene not in gene2hg
                geneGffs[locusID][i]['attributes'] += ';HG=na'
                continue
            if ome_gene2hg[gene] not in ogs:
                new_hgs.append(ome_gene2hg[gene])
            if gene == locusID or ome_gene2hg[gene] in ogs:
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

    return ome, extractedGenes, new_hgs

def extract_locus_svg_hg(extractedGenes, wrk_dir, hg2color, labels):
    for locusID, geneGff in extractedGenes.items():
        svg_path = wrk_dir + 'svg/' + locusID + '.locus.svg'
        gff2svg(
            geneGff, svg_path, product_dict = hg2color, prod_comp = r';HG=([^;]+$)', 
            width = 10, null = 'na', labels = labels, gen_new_colors = False
            )


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
        eprint('\t\t\tWARNING: ' + ome + ' could not parse gff', flush = True)
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
              out_dir, root_key = None, midpoint = True,
              ext = '.svg'):#svg_dir, out_dir):
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
    with open(tree_info[1] + adj + '.omes.' + tree_info[2], 'w') as out:
        out.write(tree.write())

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
    try:
        if og is not None:
            tree.render(out_dir + input_gene + '_HG' + str(og) + '.' + adj + ext, 
                        w=800, tree_style = ts)
        else:
            tree.render(out_dir + input_gene + '.' + adj + ext, 
                        w=800, tree_style = ts)
    except TypeError: # QStandardPaths doesn't have permissions
        # this does not work
        if not os.path.isdir(out_dir + '.XDG/'):
             os.mkdir(out_dir + '.XDG/')
        os.environ['XDG_RUNTIME_DIR'] = out_dir + '.XDG/'
        if og is not None:
            tree.render(out_dir + input_gene + '_HG' + str(og) + '.' + adj + ext, 
                        w=800, tree_style = ts)
        else:
            tree.render(out_dir + input_gene + '.' + adj + ext, 
                        w=800, tree_style = ts)

def merge_color_palette(merges, query2color):
    for merge in merges:
        query2color['|'.join(merge)] = query2color[merge[0]]
    return query2color

def extend_color_palette(hgs, color_dict):
    new_hgs = set(hgs).difference(set(color_dict.keys()))
    colors = getColors(len(hgs))
    new_colors = list(set(colors).difference(set(color_dict.values())))
    cor_i = 0
    for i, v in enumerate(list(new_hgs)):
        try:
            color_dict[str(v)] = new_colors[i-cor_i]
        except IndexError:
            eprint('\nWARNING: input too large for discrete arrow colors', flush = True)
            cor_i = i
            new_colors = colors
            color_dict[str(v)] = new_colors[i-cor_i]
    return color_dict
            

def make_color_palette(inputs, conversion_dict = {}):
    if isinstance(inputs, list): # genes list
        for i in inputs:
            if i not in conversion_dict:
                conversion_dict[i] = str(i)
    else:
        for i, v in inputs.items(): # HGs
            if str(v) in conversion_dict:
                conversion_dict[int(v)] = conversion_dict[str(v)]
            elif v not in conversion_dict:
                conversion_dict[v] = str(v)
            
        inputs = list(inputs.values())


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
        eprint('\nWARNING: input too large for discrete arrow colors', flush = True)
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
            fa2tree(query_fa_path, output_dir = tre_dir, fast = fast, 
                    cpus = cpus, verbose = verbose)
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
    clus_dir = wrk_dir + 'clus/'
    try:
        oldLog = read_json(logPath)
    except FileNotFoundError:
        oldLog = None

    rereun_search = False
    if oldLog:
        try:
            with open(oldLog['db_path'], 'rb') as raw:
                db5 = hashlib.md5(raw.read()).hexdigest()
 #           if db5 != oldLog['db']: if md5 changes do somethign
#                rerun_search = True
            if oldLog['search'] != new_log['search']:
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
    cpus = 1, verbose = False, hg = None, hgs = None, reoutput = True,
    out_keys = [], labels = True, midpoint = True, ext = '.svg'
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
    if hg:
        for ome, ome_genes2hg in genes2query.items():
            ome_hits = [x for x in query_hits if x.startswith(ome + '_')]
            if ome_hits:
                extract_loci_cmds.append([db[ome]['gff3'], ome, ome_hits, hgs, ome_genes2hg, 
                                          plusminus, wrk_dir, labels])
        with mp.Pool(processes = cpus) as pool:
            ex_locs_res = pool.starmap(extract_locus_hg, extract_loci_cmds)
        ex_locs_res = [x for x in ex_locs_res if x]
        etc_hgs, ex_locs = [], []
        for ome, ex_loc, new_hgs in ex_locs_res:
            etc_hgs.extend(new_hgs)
            ex_locs.append(ex_loc)
        new_hgs = [k for k,v in Counter(etc_hgs).items() if v > 1]
        query2color = extend_color_palette(hgs + new_hgs, query2color)
        with mp.Pool(processes = cpus) as pool:
            pool.starmap(extract_locus_svg_hg, ((ex_loc, wrk_dir, query2color, labels) \
                                                for ex_loc in ex_locs))
        
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
        try:
            svgs2tree(
                query, hg, raw_tree, db, tree_file,
                out_dir, root_key, midpoint = midpoint,
                ext = ext #svg_dir, out_dir
                )
        except NewickError:
            eprint('\t\t\tERROR: newick malformatted', flush = True)
    else:
        try:
            svgs2tree(
                query, hg, raw_tree, db, tree_file,
                out_dir, midpoint = midpoint, ext = ext #svg_dir, out_dir
                )
        except NewickError:
            eprint('\t\t\tERROR: newick malformatted', flush = True)

    return query2color
    
def hg_main(
    db, input_genes, hg_file, fast = True, out_dir = None,
    clus_cons = 0.05, clus_var = 0.65, min_seq = 3, max_size = 250, cpus = 1,
    plusminus = 10000, verbose = False, reoutput = True, interval = 0.1,
    outgroups = True, labels = True, midpoint = True, faa_dir = None,
    clus_meth = 'mmseqs easy-cluster', ext = '.svg', conversion_dict = {}
    ):
    """input_genes is a list of genes within an inputted cluster"""

    wrk_dir = out_dir + 'working/'
    gff_dir, tre_dir = wrk_dir + 'genes/', wrk_dir + 'trees/'
    clus_dir = wrk_dir + 'clus/'
    db = db.set_index()

    print('\nCompiling homolog data', flush = True)
    print('\tCompiling homologs', flush = True)
#    og_info_dict = mycodbHGs(omes = set(db['ome']))
 #   hg2gene, gene2hg = extract_ogs(og_info_dict, ogtag)
    ome2i, gene2hg, i2ome, hg2gene = compile_homolog_groups(hg_file, wrk_dir,
                                                            useableOmes = set(db.keys())) 
    input_hgs = input_genes2input_hgs(input_genes, gene2hg)
    input_hg2gene = {v: k for k, v in input_hgs.items()} # create hashes for transitioning

    todel, hits = [], set()
    for i, hg in enumerate(input_hgs):
        if hg in hits:
            todel.append(i)
        else:
            hits.add(hg)
    for i in reversed(todel):
        del input_hgs[i]
        del input_genes[i]

    hg2color, conversion_dict = make_color_palette(input_hgs, conversion_dict)

    # in the future, genes without HGs will be placed into HGs via RBH
    if not input_hgs:
        eprint('\nERROR: no HGs for any inputted genes', flush = True)
        sys.exit(3)

    hg_fas = {}
    if not all(os.path.isfile(f'{faa_dir}{hg}.faa') for gene, hg in input_hgs.items()):
        print('\tPreparing homolog fastas', flush = True)
        compile_hg_fa_cmds = [
            [db, hg2gene[hg], gene] for gene, hg in input_hgs.items() \
             if not os.path.isfile(wrk_dir + gene + '.fa')
            ]
        with mp.Pool(processes = cpus) as pool:
            hg_fas = {x[0]: x[1] for x in pool.starmap(compile_hg_fa, compile_hg_fa_cmds)}
        for gene, hg in input_hgs.items():
            if os.path.isfile(wrk_dir + gene + '.fa'): # add finished in working directory back
                hg_fas = {**hg_fas, **{gene: fa2dict(wrk_dir + gene + '.fa')}}
    else:
        for gene, hg in input_hgs.items():
            hg_fas = {**hg_fas, **{gene: fa2dict(f'{faa_dir}{hg}.faa')}}
        for gene in input_hgs:
            if gene not in hg_fas[gene]:
                eprint(f'\t\tERROR: {gene} not in homologs. Incorrect input?', flush = True)
                sys.exit(3) 


    print('\nChecking fasta sizes', flush = True)
    fas4clus, fas4trees = check_fa_size(hg_fas, max_size)
    for query, hit_fa in fas4trees.items():
        hit_fa_path = wrk_dir + query + '.fa'
        with open(hit_fa_path, 'w') as out:
            out.write(dict2fa(hit_fa))

    if not os.path.isdir(clus_dir):
        os.mkdir(clus_dir)
    if fas4clus:
        for query, fa in fas4clus.items():
            clus_fa_path = clus_dir + query + '.fa'
            if not os.path.isfile(clus_fa_path):
                with open(clus_fa_path, 'w') as out:
                    out.write(dict2fa(fa))
        print('\nRunning clustering on ' + str(len(fas4clus)) + ' fastas', flush = True)

    print('\nCRAP', flush = True)
    ome_gene2hg = gene2hg2ome2hg(gene2hg)
    if fast:
        tree_suffix = '.fa.mafft.clipkit.treefile'
    else:
        tree_suffix = '.fa.mafft.clipkit.contree'
    fas4trees = {k: v for k, v in sorted(fas4trees.items(), key = lambda x: len(x[1]))}
    for query, query_fa in fas4trees.items():
        out_keys = None
        query_hits = list(query_fa.keys())
        print('\tQuery: ' + str(query), flush = True)
        if outgroups:
            print('\t\tOutgroup detection', flush = True)
            if os.path.isfile(clus_dir + query + '.fa'):
                if not os.path.isfile(wrk_dir + query + '.outgroup.fa'):
                    in_keys, all_keys = outgroup_mngr(
                        db, query, min_seq, max_size, clus_dir, wrk_dir,
                        clus_cons = clus_cons, cpus = cpus, interval = interval,
                        verbose = False
                        )
                out_query = query + '.outgroup'
                query_hits = all_keys
                out_keys = list(set(all_keys).difference(set(in_keys)))
                print('\t\t\t' + str(len(in_keys)) + ' gene ingroup', 
                      flush = True)
                if out_keys:
                    print('\t\t\t' + str(len(out_keys)) + ' gene outgroup', 
                          flush = True)
            else:
                eprint(
                    '\t\t\tWARNING: Could not detect outgroup for root',
                    flush = True
                    )            
        HG = input_hgs[re.sub(r'\.outgroup$','',query)] # bulletproof against outgroups
        hg2color = crap_mngr(
            db, query, query_hits, out_dir, wrk_dir, tre_dir, fast, tree_suffix, 
            ome_gene2hg, plusminus, hg2color, cpus, verbose, HG, list(input_hgs.values()),
            reoutput = reoutput, out_keys = out_keys, labels = labels,
            midpoint = midpoint, ext = ext
            )

    for query in fas4clus:
        print('\tQuery: ' + str(query), flush = True)
        print('\t\tSequence clustering', flush = True)
        out_keys = None
        res, overshot, fa2clus_log = run_fa2clus(
            clus_dir + query + '.fa', db, query, min_seq, max_size, clus_dir,
            query, clus_cons, clus_var, cpus = cpus, verbose = verbose, 
            interval = interval, algorithm = clus_meth
            )
        if not res:
            print('\t\t\tERROR: query had no significant hits',
                  flush = True)
            continue
        if outgroups and not overshot:
            print('\t\tOutgroup detection', flush = True)
            if not os.path.isfile(wrk_dir + query + '.outgroup.fa'):
                in_keys, all_keys = outgroup_mngr(
                    db, query, min_seq, max_size, clus_dir, wrk_dir,
                    clus_cons = clus_cons, cpus = cpus,
                    interval = interval, verbose = False
                    )
            out_query = query + '.outgroup'
            query_hits = all_keys
            out_keys = list(set(all_keys).difference(in_keys))
            print('\t\t\t' + str(len(in_keys)) + ' gene ingroup', 
                  flush = True)
            if out_keys:
                print('\t\t\t' + str(len(out_keys)) + ' gene outgroup', 
                      flush = True)
        else:
            query_fa = fa2dict(wrk_dir + query + '.fa')
            query_hits = list(query_fa.keys())
        HG = input_hgs[re.sub(r'\.outgroup$','',query)] # bulletproof against outgroups
        hg2color = crap_mngr(
            db, query, query_hits, out_dir, wrk_dir, tre_dir, fast, tree_suffix,
            ome_gene2hg, plusminus, hg2color, cpus, verbose, HG, 
            list(input_hgs.values()), reoutput = reoutput, labels = labels,
            out_keys = out_keys, midpoint = midpoint, ext = ext
            )


def search_main(
    db, input_genes, query_fa, query_gff, binary = 'mmseqs', fast = True, 
    out_dir = None, clus_cons = 0.4, clus_var = 0.65, min_seq = 3, 
    max_size = 250, cpus = 1, reoutput = True, plusminus = 10000, 
    evalue = 1, bitscore = 40, pident = 0, mem = None, verbose = False,
    interval = 0.01, outgroups = False, conversion_dict = {}, labels = True,
    midpoint = True, clus_meth = 'mmseqs easy-linclust', ppos = 0,
    max_hits = 1000, ext = '.svg'
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
        # max hits is arbitrarily set here, can be higher
        search_fas = {**search_fas, **db2search(
            db, binary, query_path, wrk_dir, evalue = evalue, bitscore = bitscore,
            pident = pident, force = True, coverage = clus_cons,
            skip = skips, diamond = diamond, ppos = ppos, max_hits = max_hits
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
        tree_suffix = '.fa.mafft.clipkit.treefile'
    else:
        tree_suffix = '.fa.mafft.clipkit.contree'
    fas4trees = {k: v for k, v in sorted(fas4trees.items(), key = lambda x: len(x[1]))}
    for query, query_fa in fas4trees.items():
        out_keys = None
        query_hits = list(query_fa.keys())
        print('\tQuery: ' + str(query), flush = True)
        if outgroups:
            print('\t\tOutgroup detection', flush = True)
            if os.path.isfile(clus_dir + query + '.fa'):
                if not os.path.isfile(wrk_dir + query + '.outgroup.fa'):
                    in_keys, all_keys = outgroup_mngr(
                        db, query, min_seq, max_size, clus_dir, wrk_dir,
                        clus_cons = clus_cons, cpus = cpus, 
                        interval = interval, verbose = False
                        )
                out_query = query + '.outgroup'
                query_hits = all_keys
                out_keys = list(set(all_keys).difference(set(in_keys)))
                print('\t\t\t' + str(len(in_keys)) + ' gene ingroup', 
                      flush = True)
                if out_keys:
                    print('\t\t\t' + str(len(out_keys)) + ' gene outgroup', 
                          flush = True)
            else:
                eprint(
                    '\t\t\tWARNING: could not detect outgroup for root', 
                    flush = True
                    )
        crap_mngr(
            db, query, query_hits, out_dir, wrk_dir, tre_dir, 
            fast, tree_suffix, genes2query, plusminus, query2color, 
            cpus = cpus, verbose = verbose, reoutput = reoutput,
            out_keys = out_keys, labels = labels, midpoint = midpoint,
            ext = ext
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
        query_hits = list(query_fa.keys())

        res, overshot, fa2clus_log = run_fa2clus(
            clus_dir + str(query) + '.fa', 
            db, query, min_seq, max_size, 
            clus_dir, query, clus_cons, clus_var, cpus = clus_cpus, 
            verbose = verbose, interval = interval, algorithm = clus_meth,
            )
        if not res:
             print('\t\t\tFAILED: query had no significant hits',
                  flush = True)
             continue

        if outgroups and not overshot:
            print('\t\tOutgroup detection', flush = True)
            if not os.path.isfile(wrk_dir + query + '.outgroup.fa'):
            # WILL RAISE AN ERROR IF THIS EXISTS BECAUSE ALL_KEYS DOESNT
                in_keys, all_keys = outgroup_mngr(
                    db, query, min_seq, max_size, clus_dir, wrk_dir,
                    clus_cons = clus_cons, cpus = cpus, 
                    interval = interval, verbose = False
                    )
            query_hits = all_keys
            out_query = query + '.outgroup'
            out_keys = list(set(all_keys).difference(set(in_keys)))
            print('\t\t\t' + str(len(in_keys)) + ' gene ingroup', 
                  flush = True)
            if out_keys:
                print('\t\t\t' + str(len(out_keys)) + ' gene outgroup', 
                          flush = True)
        else:
            query_fa = fa2dict(wrk_dir + query + '.fa')
        crap_mngr(
            db, query, query_hits, out_dir, wrk_dir, tre_dir, 
            fast, tree_suffix, genes2query, plusminus, query2color, 
            cpus = cpus, verbose = verbose, reoutput = reoutput,
            out_keys = out_keys, labels = labels, midpoint = midpoint,
            ext = ext
            )


def cli():
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
    i_opt.add_argument('-d', '--mtdb', default = primaryDB())
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
        '-id', '--identity', type = float, default = 0,
        help = '[0 < -i < 1] Alignment identity minimum'
        )
    hg_opt.add_argument(
        '-pos', '--positives', type = float, default = 0,
        help = '[0 < -p < 1] Alignment positives minimum'
        )
    hg_opt.add_argument(
        '-b', '--bitscore', type = float,
        default = 25, help = 'Bitscore minimum for search algorithm. DEFAULT: 25'
        )
    hg_opt.add_argument(
        '-mts', '--max_target_seq', type = int, default = 1000,
        help = 'Max alignment target sequences; DEFAULT: 1000'
        )
    hg_opt.add_argument(
        '-hg', '--homologs', 
        help = 'OrthoFinder orthogroups.tsv or mmseqs cluster output'
        )
    hg_opt.add_argument(
        '-hf', '--faa',
        help = 'Directory of -hg fastas'
        )

    clus_opt = parser.add_argument_group('Phylogeny size management')
    clus_opt.add_argument(
        '--min_seq',
        help = 'Min sequences for trees/cluster size; ' \
             + 'DEFAULT: 3',
        default = 3, type = int
        )
    clus_opt.add_argument(
        '--max_seq', 
        help = 'Max sequences for trees/min for fa2clus.py; ' \
             + 'DEFAULT: 500', 
        default = 500, type = int
        )
    clus_opt.add_argument(
        '--min_cov', 
        help = 'Minimum query coverage; DEFAULT: 0.3', 
        default = 0.3, type = float
        )
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
    run_opt.add_argument('-of', '--out_format', default = 'svg',
                         help = 'Output format: ["svg", "pdf", "png"]')
    args = parser.parse_args()

    execs = ['diamond', 'clipkit', 'mafft', 'iqtree']
    if args.search:
        if args.search not in {'mmseqs', 'blastp', 'diamond'}:
            eprint('\nERROR: invalid -s', flush = True)
            sys.exit(3)
        else:
            execs.append(args.search)
            args.homologs = None
    findExecs(execs, exit = set(execs))

    if args.out_format.lower() not in {'svg', 'pdf', 'png'}:
        eprint('\nERROR: invalid -of', flush = True)
        sys.exit(10)
    else:
        out_ext = '.' + args.out_format.lower()

    if not args.homologs and args.faa:
        eprint('\nERROR: -hf requires -hg', flush = True)
        sys.exit(12)

    input_fa, input_GFF = False, False


    if args.query == '-':
        if args.gff:
            eprint('\nERROR: GFF input requires fasta input', flush = True)
            sys.exit(3)
        input_genes = stdin2str().replace('"','').replace("'",'').split()
    elif "'" in args.query or '"' in args.query:
        if args.gff:
            eprint('\nERROR: GFF input requires fasta input', flush = True)
            sys.exit(3)
        input_genes = args.query.replace('"','').replace("'",'').replace(',', ' ').split()
    else:
        input_genes = args.query.replace('"','').replace("'",'').replace(',', ' ').split()

    db = mtdb(args.mtdb)
    gene0 = input_genes[0]
    ome = input_genes[0][:input_genes[0].find('_')]
    if not ome in set(db['ome']):
        if os.path.isfile(args.query):
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
        else:
            print('\nDetected non-mycotools input', flush = True)
            if not input_fa or not args.search:
                eprint('\tnon-mycotools input requires -s, -i as a fasta, optionally -g', flush = True)
                sys.exit(4)

#        eprint('\nERROR: invalid input', flush = True)
 #       sys.exit(2)

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
        clus_meth = 'diamond'
    elif args.linclust:
        clus_meth = 'mmseqs easy-linclust'
    else:
        clus_meth = 'mmseqs easy-cluster'
    output = format_path(args.output, force_dir = True)
    args_dict = {
        'Input': ','.join(input_genes), 
        'Database': args.mtdb,
        'Orthogroup Tag': args.homologs,
        'Search binary': args.search,
        'Maximum alignment seq': args.max_target_seq,
        'Locus +/-': args.plusminus,
        'Cluster method': clus_meth,
        'Tree': tree,
        'Minimum seq': args.min_seq,
        'Maximum seq': args.max_seq,
        'Minimum coverage': args.min_cov,
        'Minimum identity': args.identity,
        'Minimum positives': args.positives,
        'Minimum Bitscore': args.bitscore,
        'GFF': args.gff,
        'Conversion file': args.conversion,
        'Labels': not args.no_label,
        'CPU': args.cpu,
        'Output directory': args.output,
        'Output extension': out_ext,
        'Verbose': args.verbose
        }

    start_time = intro('Mycotools OHCrap', args_dict, 'Jason Slot, Zachary Konkel')

    if args.conversion:
        conversion_dict = parse_conversion_file(format_path(args.conversion))
    else:
        conversion_dict = {}

    if args.homologs:
        new_log = init_log(
            args.mtdb, input_genes, 'homologs', args.bitscore,
            args.max_seq, args.plusminus, not args.no_label
            )
        print('\nPreparing output directory', flush = True)
        out_dir, wrk_dir, gff_dir, tre_dir = make_output(output, new_log)
        hg_main(
            db, input_genes, format_path(args.homologs), fast = fast, 
            out_dir = out_dir, faa_dir = format_path(args.faa),
            clus_cons = args.min_cov, clus_var = 0.65, min_seq = args.min_seq, 
            max_size = args.max_seq, cpus = args.cpu,
            verbose = args.verbose, plusminus = args.plusminus, 
            interval = 0.1, labels = not args.no_label,
            outgroups = not args.no_outgroup, clus_meth = clus_meth,
            midpoint = not bool(args.no_midpoint), ext = out_ext,
            conversion_dict = conversion_dict
            )
    else:
        new_log = init_log(
            args.mtdb, input_genes, args.search, args.bitscore,
            args.max_seq, args.plusminus, not args.no_label
            )
        print('\nPreparing output directory', flush = True)
        out_dir, wrk_dir, gff_dir, tre_dir = make_output(output, new_log)
        search_main(
            db, input_genes, input_fa, input_GFF, 
            binary = args.search, fast = fast, out_dir = out_dir,
            clus_cons = args.min_cov, clus_var = 0.65, min_seq = args.min_seq, 
            max_size = args.max_seq, cpus = 1, plusminus = args.plusminus, 
            bitscore = args.bitscore, pident = args.identity, mem = None, 
            verbose = args.verbose, interval = 0.1, ppos = args.positives,
            midpoint = not bool(args.no_midpoint),
            outgroups = not args.no_outgroup, conversion_dict = conversion_dict, 
            clus_meth = clus_meth, labels = not args.no_label,
            max_hits = args.max_target_seq, ext = out_ext
            )
    outro(start_time)


if __name__ == '__main__':
    cli()
