#! /usr/bin/env python3
#NEED TO PARSE FOR IN GENE COORDINATES AND ANNOTATIONS

import os, sys, re, argparse
from mycotools.lib.kontools import sysStart, formatPath, file2list
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
    return list(set(products))

def gff2svg( 
    gff, svg_path, product_dict, colors,
    prod_comp = gff3Comps()['product'], width = 10,
    null = 'hypothetical protein', types = {'tRNA', 'mRNA', 'rRNA'},
    max_size = 100, labels = True
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
                if product == null or not product: # if the product is a null
                # value it should be annotated as such
                    product = null
                elif product not in product_dict:
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
    labels = True
    ):

    if not set(product_dict.keys()).difference({null}): 
    # if no keys or null is the only product key
        products = compileProducts(gff_list, prod_comp, types = types)
    else:
        products = list(product_dict.keys())
    if len(products) < 16:
        colors = [
            "#000000","#004949","#009292","#ff6db6","#ffb6db",
            "#490092","#006ddb","#b66dff","#6db6ff","#b6dbff",
            "#920000","#924900","#db6d00","#24ff24","#ffff6d"
            ]
    elif len(products) < 27:
        colors = [
            "F0A3FF", "#0075DC", "#993F00", "#4C005C", "#191919",
            "#005C31", "#2BCE48", "#FFCC99", "#808080", "#94FFB5",
            "#8F7C00", "#9DCC00", "#C20088", "#003380", "#FFA405",
            "#FFA8BB", "#426600", "#FF0010", "#5EF1F2", "#00998F",
            "#E0FF66", "#740AFF", "#990000", "#FFFF80", "#FFFF00",
            "#FF5005"
            ]
    else:
        colors = [
            '#000000', '#010067', '#d5ff00', '#ff0056', '#9e008e', 
            '#0e4ca1', '#ffe502', '#005f39', '#00ff00', '#95003a', 
            '#ff937e', '#a42400', '#001544', '#91d0cb', '#620e00', 
            '#6b6882', '#0000ff', '#007db5', '#6a826c', '#00ae7e', 
            '#c28c9f', '#be9970', '#008f9c', '#5fad4e', '#ff0000', 
            '#ff00f6', '#ff029d', '#683d3b', '#ff74a3', '#968ae8', 
            '#98ff52', '#a75740', '#01fffe', '#ffeee8', '#fe8900', 
            '#bdc6ff', '#01d0ff', '#bb8800', '#7544b1', '#a5ffd2', 
            '#ffa6fe', '#774d00', '#7a4782', '#263400', '#004754', 
            '#43002c', '#b500ff', '#ffb167', '#ffdb66', '#90fb92', 
            '#7e2dd2', '#bdd393', '#e56ffe', '#deff74', '#00ff78', 
            '#009bff', '#006401', '#0076ff', '#85a900', '#00b917', 
            '#788231', '#00ffc6', '#ff6e41', '#e85ebe'
            ]
    if null not in product_dict:
        product_dict[null] = '#ffffff'

    product_dict = gff2svg(
        gff_list, svg_path, product_dict, colors,
        prod_comp = prod_comp, width = width, null = null,
        labels = labels
        )

    return product_dict


if __name__ == "__main__":

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
    parser.add_argument('-o', '--output', help = 'Optional output directory')
    args = parser.parse_args()

    if not args.regex:
        regex = gff3Comps()['product']
    else:
        regex = r''
        if args.regex.startswith(("'",'"')):
            args.regex = args.regex[1:]
        if args.regex.startswith(("'",'"')):
            args.regex = args.regex[:-1]
        for char in args.regex:
            regex += char

    args.type = args.type.replace('"','').replace("'",'')
    types = set(args.type.split())

    if args.input:
        if args.output:
            out_dir = formatPath(args.output)
        else:
            out_dir = formatPath(os.path.dirname(args.input))
        gffs = file2list(args.input)
        for entry in gffs:
            if entry.endswith('/'):
                entry = re.sub(r'/+$', '', entry)
        svg_path = os.path.dirname(gffs[0]) + re.sub(r'\.gf[^\.]+$', '.svg', os.path.basename(gffs[0]))
        product_dict = main(gff2list(gffs[0]), svg_path, prod_comp = regex, types = types)
        for gff in gffs[1:]:
            svg_path = os.path.dirname(gff) + re.sub(r'\.gf[^\.]+$', '.svg', os.path.basename(gff))
            product_dict = main(gff2list(gff), svg_path, product_dict, args.width, prod_comp = regex, types = types)
    else:
        if args.gff.endswith('/'):
            args.gff = re.sub(r'/+$', '', args.gff)
        if args.output:
            svg_path = formatPath(args.output) + re.sub(
                r'\.gf[^\.]+$', '.svg', os.path.basename(args.gff)
                )
        else:
            svg_path = os.path.dirname(args.gff) + re.sub(
                r'\.gf[^\.]+$', '.svg', os.path.basename(args.gff)
                )
        main( gff2list(args.gff), svg_path, width = args.width, prod_comp = regex, types = types )

    sys.exit(0)
