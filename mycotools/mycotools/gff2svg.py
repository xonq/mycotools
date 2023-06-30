#! /usr/bin/env python3

#NEED TO PARSE FOR IN GENE COORDINATES AND ANNOTATIONS
#NEED to create a single file output option for multiple inputs

import os
import re
import sys
import random
import argparse
from mycotools.lib.kontools import sys_start, format_path, file2list, getColors
from dna_features_viewer import GraphicFeature, GraphicRecord
from mycotools.lib.biotools import gff2list, gff3Comps

def compileProducts(gff, prod_comp, types = {'tRNA', 'mRNA', 'rRNA'}):
    # find all the product attributes for color pallette selection
    products = []
    for entry in gff:
        if entry['type'] in types:
            prod = re.search(prod_comp, entry['attributes'])
            if prod is not None:
                product = prod[1].lower()
                products.append(product)
    out_products = list(set([
        x.replace('"', '') for x in products
        ]))
    return out_products

def gff2svg( 
    gff, svg_path, product_dict, colors,
    prod_comp = gff3Comps()['product'], width = 10,
    null = 'hypothetical protein', types = {'tRNA', 'mRNA', 'rRNA'},
    max_size = 100, labels = True, gen_new_colors = True
    ):
    # needs to end with .svg
    if not svg_path.endswith('.svg'):
        svg_path += '.svg'

    seqs, lens = {}, {}
    count = len(product_dict) - 1
    for i in gff: # for each line in the gff_list
        if i['seqid'] not in seqs: # add the contig
            seqs[i['seqid']], lens[i['seqid']] = [], []
        if i['type'] in types: # if it is a sequence type of interest
            prod = re.search(prod_comp, i['attributes']) # find the product
            if prod is not None:
                product = prod[1].lower() # lower case for uniformity
                product = product.replace('"', '')
                if product in product_dict:
                    pass
                elif product == null or not product: # if the product is a null
                # value it should be annotated as such
                    product = null
                elif not gen_new_colors:
                    product = null
                elif product not in product_dict:
                    print(product)
                    try: # try to use the next color
                        product_dict[product] = colors[count]
                        count += 1
                    except IndexError: # if out of colors, must restart
                        count = 0
                        product_dict[product] = colors[count]
                        count += 1
            else:
                product = null
            # append the start, end, strand sense, color, and max name size
            seqs[i['seqid']].append([
                int(i['start']), int(i['end']),
                int(i['strand'] + '1'), product_dict[product],
                product[:max_size+1]])
            # append the lengths independently to lens
            lens[i['seqid']].append( int(i['start']) )
            lens[i['seqid']].append( int(i['end']) )

    for seq in seqs:
        # sort each contig length set, start becomes the smallest
        lens[seq].sort()
        start, features = lens[seq][0], []
        for i in seqs[seq]:
#            i[0] -= start
 #           i[1] -= start
            if i[-1] == null or not labels:
                features.append(GraphicFeature(
                    start = i[0], end = i[1], strand = i[2],
                    color = i[3])
                    )
            else:
                features.append(GraphicFeature(
                    start = i[0], end = i[1], strand = i[2],
                    color = i[3], label = i[4])
                    )
        length = lens[seq][-1] - lens[seq][0]
        record = GraphicRecord(
            sequence_length=length, features = features, first_index = start
            )
        ax, _ = record.plot(figure_width=width, strand_in_label_threshold=7)
        ax.figure.savefig( svg_path, bbox_inches = 'tight'  )

#    CLOSE(ax)

    return product_dict 


def main(
    gff_list, svg_path, product_dict = {}, 
    width = 10, prod_comp = gff3Comps()['product'], 
    null = 'hypothetical protein', types = {'tRNA', 'mRNA', 'rRNA'},
    labels = True, wheel = None, shuffle = False, gen_new_colors = True
    ):

    product_dict = {str(k.lower()): v for k, v in product_dict.items()}
    if not wheel and not product_dict:
        if not set(product_dict.keys()).difference({null}): 
        # if no keys or null is the only product key
            products = compileProducts(gff_list, prod_comp, types = types)
        else:
            products = list(product_dict.keys())
        colors = getColors(len(products))
    elif wheel == 1: # spoof function to get wheel
        colors = getColors(1)
    elif wheel == 2:
        colors = getColors(17)
    elif wheel == 3:
        colors = getColors(28)
    else:
        colors = None

    if shuffle:
        random.shuffle(colors)

    if null not in product_dict:
        product_dict[null] = '#ffffff'

    product_dict = gff2svg(
        gff_list, svg_path, product_dict, colors,
        prod_comp = prod_comp, width = width, null = null,
        labels = labels, types = types, gen_new_colors = gen_new_colors
        )

    return product_dict


def cli():

    parser = argparse.ArgumentParser(description = 'Converts .gff3(s) to .svg')
    parser.add_argument('-g', '--gff')
    parser.add_argument('-i', '--input', help = 'New line delimited list of gffs')
    parser.add_argument('-w', '--width', type = int, default = 10, 
        help = '.svg width; DEFAULT: 10')
    parser.add_argument(
        '-r', '--regex',
        help = 'Expression to gather annotation. DEFAULT: "product=([^;]*)"'
        )
    parser.add_argument(
        '-t', '--type', default = '"tRNA mRNA rRNA"',
        help = 'Type to extract. DEFAULT: "tRNA mRNA rRNA"'
        )
    parser.add_argument(
        '-n', '--nolabel', action = 'store_true', default = False,
        help = "Don't output labels on svg"
        )
    parser.add_argument(
        '-c', '--color', default = None, type = int,
        help = 'Color wheel {1 smallest, 2, 3 largest}; DEFAULT: auto'
        )
    parser.add_argument(
        '-s', '--shuffle', help = 'Randomly shuffle color wheel', action = 'store_true'
        )
    parser.add_argument('-o', '--output', help = 'Optional output directory')
    args = parser.parse_args()

    if not args.regex:
        regex = gff3Comps()['product']
    else:
        regex = r''
        if args.regex.startswith(("'",'"')):
            args.regex = args.regex[1:]
        if args.regex.endswith(("'",'"')):
            args.regex = args.regex[:-1]
        for char in args.regex:
            regex += char

    if args.color:
        if args.color not in {1,2,3}:
            raise ValueError('--color must be in {1, 2, 3}')

    args.type = args.type.replace('"','').replace("'",'')
    types = set(args.type.split())

    if args.input:
        if args.output:
            out_dir = format_path(args.output)
        else:
            out_dir = format_path(os.path.dirname(args.input))
        gffs = file2list(args.input)
        for entry in gffs:
            if entry.endswith('/'):
                entry = re.sub(r'/+$', '', entry)
        svg_path = os.path.dirname(gffs[0]) + re.sub(r'\.gf[^\.]+$', '.svg', os.path.basename(gffs[0]))
        product_dict = main(
            gff2list(gffs[0]), svg_path, prod_comp = regex, types = types, labels = not args.nolabel,
            wheel = args.color, shuffle = args.shuffle
            )
        for gff in gffs[1:]:
            svg_path = os.path.dirname(gff) + re.sub(r'\.gf[^\.]+$', '.svg', os.path.basename(gff))
            product_dict = main(
                gff2list(gff), svg_path, product_dict, args.width, prod_comp = regex, 
                types = types, labels = not args.nolabel, wheel = args.color, shuffle = args.shuffle
                )
    else:
        if args.gff.endswith('/'):
            args.gff = re.sub(r'/+$', '', args.gff)
        if args.output:
            svg_path = format_path(args.output) + re.sub(
                r'\.gf[^\.]+$', '.svg', os.path.basename(args.gff)
                )
        else:
            svg_path = os.path.dirname(args.gff) + re.sub(
                r'\.gf[^\.]+$', '.svg', os.path.basename(args.gff)
                )
        main(
            gff2list(args.gff), svg_path, width = args.width, prod_comp = regex, 
            types = types, labels = not args.nolabel, wheel = args.color,
            shuffle = args.shuffle
            )

    sys.exit(0)


if __name__ == '__main__':
    cli()
