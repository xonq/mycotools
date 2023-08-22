#!/usr/bin/env python3

import os
import sys
import shutil
import argparse
import subprocess
import numpy as np
import multiprocessing as mp
from tqdm import tqdm
from itertools import combinations
from collections import defaultdict, Counter
from mycotools.db2files import soft_main as symlink_files
from mycotools.lib.kontools import format_path, mkOutput, \
                                   findExecs, intro, outro, \
                                   eprint
from mycotools.lib.dbtools import mtdb, primaryDB
from mycotools.lib.biotools import gff2list


def run_mmseqs(db, wrk_dir, algorithm = 'mmseqs easy-cluster', 
               min_id = 0.3, min_cov = 0.3, sensitivity = 7.5,
               cpus = 1):
    symlink_files(['faa'], db, wrk_dir, verbose = False) # symlink proteomes
    cluster_res_file = wrk_dir + 'homolog_groups.tsv'
    if not os.path.isfile(cluster_res_file): # NEED to add to log removal
        # be cautious about shell injection because we need to glob
        int(cpus)
        float(min_id)
        float(min_cov)
        if not os.path.isdir(wrk_dir):
            raise OSError('invalid working directory')
        elif not algorithm in {'mmseqs easy-linclust', 'mmseqs easy-cluster'}:
            raise OSError('invalid mmseqs binary')
        mmseqs_cmd = subprocess.call(algorithm + ' ' + ' '.join([wrk_dir + 'faa/*faa',
                                      wrk_dir + 'cluster', wrk_dir + 'tmp/',
                                      '--min-seq-id', str(min_id), '--threads',
                                      str(cpus), '--compressed', '1', '-s', 
                                      str(sensitivity),
                                      '--cov-mode', '0', '-c', str(min_cov)]),
                                      shell = True, stdout = subprocess.DEVNULL)
#                                      stderr = subprocess.DEVNULL)
        shutil.move(wrk_dir + 'cluster_cluster.tsv', cluster_res_file)
    elif os.path.getsize(cluster_res_file):
        mmseqs_cmd = 0
    else:
        mmseqs_cmd = 1
    if mmseqs_cmd:
        eprint('\tERROR: cluster failed')
        sys.exit(1)
    if os.path.isfile(wrk_dir + 'cluster_all_seqs.fasta'):
        os.remove(wrk_dir + 'cluster_all_seqs.fasta')
    if os.path.isfile(wrk_dir + 'cluster_rep_seq.fasta'):
        os.remove(wrk_dir + 'cluster_rep_seq.fasta')
    if os.path.isdir(wrk_dir + 'tmp/'):
        shutil.rmtree(wrk_dir + 'tmp/')
    return cluster_res_file

        
def parse_orthofinder(hg_file, useableOmes = set()):
    """
    imports orthofinder Orthogroups.txt "hg_file". outputs several data structures:
    ome_num = {ome: number}, gene2hg = {gene: og}, i2ome = [ome0, ome1, ome2]
    """

    gene2hg, ome_num, i2ome, hg2gene = \
        {}, {}, [], {}
    with open(hg_file, 'r') as raw:
        for line in raw:
            data = line.rstrip().split() # OrthoFinder input
            og = int(data[0].replace(':','').replace('OG',''))
            hits = [x for x in data[1:] if x[:x.find('_')] in useableOmes]
            omes = [x[:x.find('_')] for x in hits]
            if len(set(omes)) < 2: # skip singleton organisms
                continue
            hg2gene[og] = hits
            for i, gene in enumerate(hits):
                ome = omes[i]
#                ome = gene[:gene.find('_')] # wtf, parsing omes raises and index error
                if ome not in ome_num:
                    i2ome.append(ome)
                    ome_num[ome] = len(i2ome) - 1
                gene2hg[gene] = og

    return ome_num, gene2hg, i2ome, hg2gene

def parse_1to1(hg_file, useableOmes = set()):
    derivations = defaultdict(list)
    with open(hg_file, 'r') as raw:
        for line in raw:
            k, v = line.rstrip().split() # all white space
            derivations[k].append(v)
  
    hgs_list = []
    for k, v in derivations.items():
        v.append(k) 
        hgs_list.append(sorted(set(v)))
        
    hg2gene = {i: v for i, v in enumerate(
                    sorted(hgs_list, key = lambda x: len(x),
                    reverse = True)
                    )}
    
    todel = []
    for og, genes in hg2gene.items():
        genes = [x for x in genes if x[:x.find('_')] in useableOmes]
        omes = set([x[:x.find('_')] for x in genes])
        if len(omes) < 2: # skip singleton omes
            todel.append(og)
    for og in todel:
        del hg2gene[og]    
    hg2gene = {k: v for k, v in sorted(hg2gene.items(), key = lambda x:
                                       len(x[1]), reverse = True)}
    
    gene2hg = {}
    for i, og_list in hg2gene.items():
        for gene in og_list:
            gene2hg[gene] = i
        
    i2ome = sorted(set([x[:x.find('_')] for x in list(gene2hg.keys())]))
    ome_num = {v: i for i, v in enumerate(i2ome)}
        
    return ome_num, gene2hg, i2ome, hg2gene


def compile_homolog_groups(hg_file, wrk_dir = None,
                           useableOmes = set()):
    try:
        hg_info = parse_orthofinder(hg_file, useableOmes)
    except ValueError:
        hg_info = parse_1to1(hg_file, useableOmes)
        hg2genes = hg_info[-1]
        with open(wrk_dir + 'homolog_groups.tsv', 'w') as out:
            for hg, genes in hg2genes.items():
                out.write(str(hg) + '\t' + ' '.join(genes) + '\n')

    return hg_info


def compile_cds(gff_list, ome, gene2hg):
    """ 
    Inputs the gff_list and organism ome code. Compiles CDS entries for loci parsing and outputs:
    cds_dict = {contig: {protein: [[CDSstart, CDSstop]] } }
    Sorts the cds_dict for each protein from smallest to largest
    """

    cds_dict, fail = defaultdict(dict), False
    for entry in gff_list:
        if entry['type'] == 'CDS':
            prot_prep_i0 = entry['attributes'].index(';Alias=') # grab the
            # mycotools accession index
            try: # grab the end of the Alias by grabbing the index for a ';',
            # that is after the accession tag's occurrence.
            # then add that to the start of the accession tag and add one to
            # remove that semicolon. If there is a value error, the accession 
            # is at the end of the attributes (gff col 8 (?))
                prot_prep_i1 = entry['attributes'][prot_prep_i0+1:].index(';') + prot_prep_i0+1
            except ValueError:
                prot_prep_i1 = len(entry['attributes'])
            prot = entry['attributes'][prot_prep_i0:prot_prep_i1].replace(';Alias=','')
            # obtain the protein accession
            if prot not in cds_dict[entry['seqid']] and prot: # add the protein
            # to the cds_dict entry
                cds_dict[entry['seqid']][prot] = []
            elif not prot: # if there isn't a valid accession it may mean the
            # mycotools curation did not work or the user did not curate
            # correctly 
                print(entry['attributes'], prot_prep_i0, prot_prep_i1)
                if not fail:
                    print('\tWARNING: ' + ome + ' has proteins in gff with no Alias', flush = True)
                fail = True
                continue
            cds_dict[entry['seqid']][prot].extend([int(entry['start']),
                int(entry['end'])]) # add the coordinates

    hg_dict = {}
    for contig in cds_dict:
        for prot in cds_dict[contig]:
            cds_dict[contig][prot].sort() # sort the coordinates of the proteins
            # lowest to highest
        cds_dict[contig] = list(sorted(cds_dict[contig].keys(), key = lambda k:
            cds_dict[contig][k][0])) # sort the proteins in the contigs lowest
            # to highest coordinates
        hg_dict[contig] = []
        for prot in cds_dict[contig]:
            try:
                hg_dict[contig].append(gene2hg[prot])
            except KeyError:
                hg_dict[contig].append(None)

    return hg_dict


def parse_loci(
    gff_path, ome, gene2hg, window = 6
    ):
    """obtain a set of tuples of HG pairs {(OG0, OG1)...}"""

    gff_list = gff2list(gff_path) # open here to improve pickling
    hg_dict = compile_cds(gff_list, os.path.basename(gff_path).replace('.gff3',''),
                          gene2hg)
    pairs = []
    for scaf, hgs in hg_dict.items(): # for each contig
        windows = [sorted(set([x for x in hgs[i:i+window+1] \
                               if x is not None])) \
                   for i in range(len(hgs) - window + 1)]
        for w in windows:
            pairs.extend([x for x in combinations(w, 2)])

    out_pairs = set(pairs) # unique pairs of OGs
    return ome, out_pairs


def compile_loci(
    db, ome2i, gene2hg, window, cpus = 1
    ):

    loci_hash_cmds = [
        [v['gff3'], ome,
        gene2hg, window]
        for ome, v in db.items() \
        if ome in ome2i
        ]
    with mp.get_context('fork').Pool(processes = cpus) as pool:
        loci_hashes = pool.starmap(parse_loci, tqdm(loci_hash_cmds, total = len(ome2i)))
    pool.join()
    pairs = {}
    pairs = {x[0]: x[1] for x in loci_hashes if x[1]}

    return pairs


def form_cooccur_array(cooccur_dict, ome2i):

    count, hgx2i, size_dict, cooccur_arrays, i2hgx = 0, {}, {}, {}, {}
    cooccur_dict = {k: tuple(sorted(v)) for k, v in sorted(cooccur_dict.items(), key = lambda x: len(x[0]))}
    cooccur_arrays = np.zeros([len(ome2i), len(cooccur_dict)], dtype = np.int32)

    old_len = len(list(cooccur_dict.keys())[0])
    for i, hgx in enumerate(list(cooccur_dict.keys())):
        i2hgx[i] = hgx
        hgx2i[hgx] = i
        for ome in cooccur_dict[hgx]:
#            cooccur_arrays[ome, i] = hgx_dict[ome][hgx]
            cooccur_arrays[ome2i[ome], i] = 1

    return i2hgx, hgx2i, cooccur_dict, cooccur_arrays



def form_cooccur_structures(pairs, min_omes, ome2i, cc_arr_path = None):
    """
    Imports the out_dicts from parse_loci and creates index vectors that bear
    the ome_num's with a given  cooccurence. e.g. cooccur_dict[(og1, og2)] = [1, 2, 3]
    """

    cooccur_dict = defaultdict(list)
    for ome in pairs:
        for key in pairs[ome]:
            new_key = tuple(sorted(key))
            cooccur_dict[new_key].append(ome)

    cooccur_dict = {
        x: tuple(cooccur_dict[x]) for x in cooccur_dict \
        if len(cooccur_dict[x]) > min_omes
        }
            
    i2hgpair, hgpair2i, cooccur_dict, cooccur_array = form_cooccur_array(
        cooccur_dict, ome2i
        )   
            
    return cooccur_array, dict(cooccur_dict), hgpair2i, i2hgpair


def id_near_schgs(hg2gene, omes, max_hgs = 10000, max_median = 2, max_mean = 2):
    """Identify near single copy homology groups. Criteria:
       1) has all omes; 2) minimal overall size of homology group;
       3) lowest median"""
    near_schgs = []
    min_hg2gene = {k: v for k, v in sorted(hg2gene.items(),
                    key = lambda x: len(x[1])) \
                   if len(v) >= len(omes)}

    max_len = max_mean * len(omes) # dont need to compute means iteratively
    median_i = round((len(omes) / 2) - 0.5)

    for hg, genes in min_hg2gene.items():
        hg_omes = [x[:x.find('_')] for x in genes]
        if not omes.difference(hg_omes): # all omes are present
            count_omes = list(Counter(hg_omes).values()) # prepare for median calculation
            count_omes.sort()
            median = count_omes[median_i]
#            print(median, max_median, len(genes), max_len)
            if median <= max_median and len(genes) <= max_len:
 #               print('\t', median, len(genes))
                near_schgs.append(hg)
            elif len(genes) >= max_len: # everything else will be higher
                break
            if len(near_schgs) == max_hgs:
                break

    return near_schgs


def extract_nschg_pairs(nschgs, hgpair2i, m_arr):
    nschg_set = set(nschgs)
    nschg_pairs = [i for hgpair, i in hgpair2i.items() \
                   if any(x in nschg_set for x in hgpair)]
    nschgs_arr = np.array(nschg_pairs)
    valid_arr = m_arr[:, nschgs_arr]
    return valid_arr


def align_microsynt_np(m_arr, i2ome, hg2gene, hgpair2i, wrk_dir, nschgs = None):
    if not nschgs:
        nschgs = id_near_schgs(hg2gene, set(i2ome), max_hgs = 100,
                           max_median = 4, max_mean = 3)
    trm_arr = extract_nschg_pairs(nschgs, hgpair2i, m_arr)
    with open(wrk_dir + 'microsynt.align.phy', 'w') as out:
        out.write(f'{trm_arr.shape[0]} {trm_arr.shape[1]}\n')
        [out.write(f'{i2ome[i]} {"".join([str(x) for x in trm_arr[i]])}\n') \
                   for i in range(trm_arr.shape[0])]
    return wrk_dir + 'microsynt.align.phy'


def remove_nulls(cc_arr):
    sum_arr = np.sum(cc_arr, axis = 1) # sum all ome rows
    null_i_list = list(np.where(sum_arr == 0)[0])
    del_list = sorted(null_i_list, reverse = True)
    for i in sorted(null_i_list, reverse = True):
        cc_arr = np.delete(cc_arr, i, axis = 0)
    return cc_arr, del_list


def run_tree(alignment, wrk_dir, constraint = False, iqtree = 'iqtree',
             model = 'GTR2+FO+ASC+R5', verbose = False, cpus = 1):
    
    tree_dir = wrk_dir + 'tree/'
    if not os.path.isdir(tree_dir):
        os.mkdir(tree_dir)
    
    prefix = tree_dir + 'microsynt'
    tree_cmd = [iqtree, '-s', alignment, '-m', model,
                '--prefix', prefix] 
    if constraint:
        tree_cmd.extend(['-te', constraint])
        comp_file = prefix + '.treefile'
    else: # use bootstrap
        tree_cmd.extend(['-B', '1000'])     
        comp_file = prefix + '.contree'
        
    if verbose:                             
        v = None                            
    else:
        v = subprocess.DEVNULL
    tree_res = subprocess.call(tree_cmd, stdout = v)
    shutil.copy(comp_file, wrk_dir + '../microsynt.newick')
        
    return tree_res


def main(db, hg_file, out_dir, wrk_dir, algorithm, 
         tree_path, plusminus = 3, min_cov = 0, min_id = 0.3, 
         n50thresh = None, near_single_copy_genes = [], 
         constraint = None, verbose = False, return_post_compile = False,
         cpus = 1):

    # obtain useable omes 
    useableOmes, dbOmes = set(), set(db.keys())
    print('\nI. Inputting data', flush = True)
    if n50thresh: # optional n50 threshold via mycotoolsDB
        assemblyPath = format_path('$MYCODB/../data/assemblyStats.tsv')
        if os.path.isfile(assemblyPath): 
            with open(assemblyPath, 'r') as raw:
                for line in raw:
                    d = line.rstrip().split('\t')
                    ome, omeN50 = d[1], float(d[2])
                    if omeN50 > n50thresh and ome in dbOmes:
                        useableOmes.add(ome)
        else:
            raise FileNotFoundError(assemblyPath + ' not found. N50 threshold not applied')
    else:
        useableOmes = dbOmes
    
    # initialize orthogroup data structures            
    if not hg_file and not os.path.isfile(wrk_dir + 'homology_groups.tsv'):
        hg_file = run_mmseqs(db, wrk_dir, algorithm = algorithm,
                             min_id = min_id, cpus = cpus)
    print('\tParsing homology groups (HGs)', flush = True)
    ome2i, gene2hg, i2ome, hg2gene = compile_homolog_groups(hg_file, wrk_dir, 
                                                            useableOmes)

# remove and reactivate other once otefa is finished
#    if return_post_compile:
 #       return ome2i, gene2hg, i2ome, hg2gene, None, None


    
    missing_from_db = set(ome2i.keys()).difference(set(db.keys()))
    print('\t\tOmes:', len(ome2i), flush = True)
    if missing_from_db:
        print(f'\t\t\t{len(missing_from_db)} omes in HGs but not mtdb')
        for ome in list(missing_from_db):
            del ome2i[ome]
        todel = []
        for gene in gene2hg:
            ome = gene[:gene.find('_')]
            if ome in missing_from_db:
                todel.append(gene)
        for gene in todel:
            hg = gene2hg[gene]
            del gene2hg[gene]
            hg2gene[hg].pop(hg2gene[hg].index(gene))
            if not hg2gene[hg]:
                del hg2gene[hg]
        i2ome = [ome for ome in ome2i]
        ome2i = {ome: i for i, ome in enumerate(i2ome)} 
        with open(wrk_dir + 'ome2i.tsv', 'w') as out:
            out.write(
                '\n'.join([k + '\t' + str(v) for k, v in ome2i.items()])
                )

    print('\t\tHGs:', len(hg2gene), flush = True)
    print('\t\tGenes:', len(gene2hg), flush = True)



     # compile cooccuring pairs of homogroups in each genome
    print('\tCompiling all loci', flush = True)
    cc_arr_path = wrk_dir + 'microsynt'
    ome2pairs = compile_loci(
        db, ome2i, gene2hg, plusminus*2+1,
        cpus = cpus
        )

    cooccur_dict = None
#    if not os.path.isfile(out_dir + 'hgps.tsv.gz'):     
    # assimilate cooccurrences across omes
    print('\tIdentifying cooccurences', flush = True)
    
    seed_len = sum([len(ome2pairs[x]) for x in ome2pairs])
    print('\t\t' + str(seed_len) + ' initial HG-pairs', flush = True)
    cooccur_array, cooccur_dict, hgpair2i, i2hgpair = \
        form_cooccur_structures(ome2pairs, 2, ome2i, cc_arr_path)
    max_omes = max([len(cooccur_dict[x]) for x in cooccur_dict])
    print('\t\t' + str(max_omes) + ' maximum organisms with HG-pair', flush = True)
    cooccur_array[cooccur_array > 0] = 1
    cooccur_array.astype(np.uint8)
    print('\t\t' + str(sys.getsizeof(cooccur_array)/1000000) + ' MB', flush = True)
    cooccur_array, del_omes = remove_nulls(cooccur_array)
    for i in del_omes:
        print(f'\t\t\t{i2ome[i]} removed for lacking overlap')
        del i2ome[i]
    
    ome2i = {v: i for i, v in enumerate(i2ome)}
    with open(wrk_dir + 'ome2i.tsv', 'w') as out:
        out.write(
            '\n'.join([k + '\t' + str(v) for k, v in ome2i.items()])
            )
    cooccur_dict = {k: tuple(sorted([ome2i[x] for x in v])) \
                    for k, v in cooccur_dict.items()}
    if return_post_compile:
        return ome2i, gene2hg, i2ome, hg2gene, None, None


#        np.save(cc_arr_path, cooccur_array)
#    elif not os.path.isfile(tree_path):
#    elif not os.path.isfile(tree_path):
 #       cooccur_array = np.load(cc_arr_path + '.npy')

    ome2pairs = {ome2i[ome]: v for ome, v in ome2pairs.items() \
                 if ome in ome2i}
    microsynt_dict = {}
    print('\nII. Microsynteny tree', flush = True)
    if not os.path.isfile(tree_path):
        nschgs = []
        if near_single_copy_genes:
            try:
                for entry in near_single_copy_genes:
                    nschgs.append(int(entry))
                nschgs = sorted(set(nschgs))
            except ValueError:
                nschgs = sorted({gene2hg[gene] \
                             for gene in near_single_copy_genes \
                             if gene in gene2hg})
        # create microsynteny distance matrix and make tree
        print('\tPreparing microsynteny alignment', flush = True)
        align_file = align_microsynt_np(cooccur_array, i2ome, hg2gene,
                                        hgpair2i, wrk_dir, nschgs)
        print('\tBuilding microsynteny tree', flush = True)
        run_tree(align_file, wrk_dir, constraint = constraint, iqtree = 'iqtree',
                 model = 'GTR2+FO+ASC+R5', verbose = False, cpus = cpus)
        # too bulky to justify keeping
        os.remove(cc_arr_path + '.npy')

    return ome2i, gene2hg, i2ome, hg2gene, ome2pairs, cooccur_dict

def cli():
    parser = argparse.ArgumentParser(description = 'Generate a microsynteny ' \
      + 'tree that recapitulates the divergence of gene order within small ' \
      + 'loci.')
    parser.add_argument('-d', '--db', default = primaryDB(), 
        help = 'MycotoolsDB. DEFAULT: masterdb')
    parser.add_argument('-f', '--focal_genes',
        help = 'File of genes for neighborhood extraction, e.g. SCOs')
    parser.add_argument('-of', '--orthofinder', 
        help = 'Precomputed OrthoFinder output directory')
    parser.add_argument('-i', '--input', 
        help = 'Precomputed cluster results file "gene\tcluster#"')
    parser.add_argument('-l', '--linclust', action = 'store_true',
        help = 'Run linclust instead of easy-cluster if lacking precomputed ' \
             + 'groups')
    parser.add_argument('-w', '--window', default = 5, type = int,
        help = 'Max genes +/- for locus window. DEFAULT: 5 (11 gene window)')
    parser.add_argument('-t', '--topology_constraint',
        help = 'Constrain microsynteny topology to species tree')
    parser.add_argument('-c', '--cpus', default = 1, type = int)
    parser.add_argument('-o', '--output')
    args = parser.parse_args()

    execs = ['iqtree']
    if args.orthofinder:
        of_out = format_path(args.orthofinder)
        if os.path.isdir(of_out):
            homogroups = of_out + '/Orthogroups/Orthogroups.txt'
            hg_dir = of_out + '/Orthogroup_Sequences/'
        else:
            homogroups = of_out
            hg_dir = os.path.dirname(of_out) + '../Orthogroup_Sequences/'
        if not os.path.isfile(hg_dir + 'OG0000000.fa'):
            hg_dir = None 
        method = 'orthofinder' 
    elif args.input:
        homogroups = format_path(args.input)
        hg_dir = None
        method = 'mmseqs easy-cluster'
    elif args.linclust:
        method = 'mmseqs easy-linclust'
        hg_dir = None
        homogroups = None
        execs.append('mmseqs')
    else:
        method = 'mmseqs easy-cluster'
        hg_dir = None
        homogroups = None
        execs.append('mmseqs')
    findExecs(execs, exit = set(execs))

    if not args.output:
        out_dir = mkOutput(os.getcwd() + '/', 'db2microsyntree')
    else:
        out_dir = mkOutput(format_path(args.output), 'db2microsyntree')

    args_dict = {'Database': args.db, 'Focal genes': args.focal_genes,
                 'Preclassified homologs': homogroups,
                 'Window size': args.window,
                 'Constraint': args.topology_constraint,
                 'CPUS': args.cpus, 'output': out_dir}
    intro('db2microsyntree', args_dict, 'Zachary Konkel')

    wrk_dir = out_dir + 'working/'
    if not os.path.isdir(wrk_dir):
        os.mkdir(wrk_dir)

    if args.focal_genes:
        with open(format_path(args.focal_genes), 'r') as raw:
            focal_genes = [x.rstrip() for x in raw.read().split() \
                           if x.rstrip()]
    else:
        focal_genes = []

    tree_path = out_dir + 'microsynteny.newick'

    db = mtdb(format_path(args.db)).set_index()

    main(db, homogroups, out_dir, wrk_dir, method,
         tree_path, plusminus = args.window, min_cov = 0, min_id = 0.3,
         n50thresh = None, near_single_copy_genes = focal_genes,
         constraint = format_path(args.topology_constraint), verbose = False, 
         cpus = args.cpus)


if __name__ == '__main__':
    cli()
