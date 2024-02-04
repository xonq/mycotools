#! /usr/bin/env python3

import os
import sys
import argparse
from itertools import chain
from cogent3 import PhyloNode, load_tree
from mycotools.lib.kontools import format_path, split_input, eprint, vprint
from mycotools.lib.dbtools import mtdb

def compile_tree(tree_path, root = [], verbose = False):
    """Compile a phylogeny from a path and conver the tip names to an index"""
    with open(tree_path, 'r') as raw:
        nwk = raw.read()

    phylo = load_tree(tree_path)
    tips = set(phylo.get_tip_names())
    vprint(f'{len(tips)} tips on input', v = verbose, e = True)

    if len(root) == 1:
        phylo = phylo.rooted_with_tip(root[0])
        phylo.write(tree_path, with_distances = True)
    elif len(root) > 1:
        nodes = {k: (v, len(v.get_tip_names())) \
                 for k, v in phylo.get_nodes_dict().items() \
                 if set(root).issubset(set(v.get_tip_names()))}
        mrca_tip_len = min([v[1] for v in list(nodes.values())])
        mrca_edge = [k for k, v in nodes.items() if v[1] == mrca_tip_len]
        phylo = phylo.rooted_at(mrca_edge)

    return phylo


def prune_to_new(phylo, tips = [], spacer = ''):
    """Prune missing tips from a rooted phylogeny"""
    missing = set(phylo.get_tip_names()).difference(tips)
    if missing:
        eprint(f'{spacer}Removing {len(missing)} tips', flush = True)
        todel = []
        for n in phylo.tips():
            if n.name in missing:
                todel.append(n)
        for n in todel:
            n.parent.remove(n)
            phylo.prune()
    return phylo


def main(phylo_path, db = None, 
         tips = [], root = [], trim = False,
         phylo = None, verbose = True):

    db = db.set_index()
    if phylo_path:
        phylo = compile_tree(phylo_path, root, verbose = verbose)

    if trim:
        phylo = prune_to_new(phylo, tips)

    return phylo


def cli():
    parser = argparse.ArgumentParser(description = \
        'Phylogeny manipulation tools')
    parser.add_argument('-i', '--input', help = 'Input phylogeny',
                        required = True)
    parser.add_argument('-t', '--tips', help = 'Tips to reference')
    parser.add_argument('-d', '--mtdb', help = 'MTDB to reference')
    parser.add_argument('-p', '--prune', action = 'store_true',
                        help = 'Remove omes absent from input')
    parser.add_argument('-r', '--root', action = 'store_true',
                        help = 'Root on tip(s)')
    args = parser.parse_args()

    root = split_input(args.root)

    if args.tips:
        tips_path = format_path(args.tips)
        if os.path.isfile(tips_path):
            eprint('\nDetecting path input', flush = True)
            with open(tips_path, 'r') as raw:
                tips = list(chain(*[line.rstrip().split() for line in raw]))
        else:
            tips = split_input(args.tips)
    else:
        tips = []

    if args.mtdb:
        db = mtdb(format_path(args.mtdb))
        if not tips:
            tips = sorted(db['ome'])
    else:
        db = None

    phylo = main(format_path(args.input), db, tips, trim = args.prune)
    eprint(f'{len(phylo.get_tip_names())} tips on output', flush = True)
    print(phylo.get_newick(with_distances = True), flush = True)


if __name__ == '__main__':
    cli()
    sys.exit(0)
