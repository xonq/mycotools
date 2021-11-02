#! /usr/bin/env python3
#NEED TO PARSE FOR IN GENE COORDINATES AND ANNOTATIONS

import os, sys, re, argparse
from mycotools.lib.kontools import sysStart, formatPath, file2list
from dna_features_viewer import GraphicFeature, GraphicRecord
from mycotools.lib.biotools import gff2list, gff3Comps



def gff2svg( 
    gff_path, svg_dir, product_dict,
    prod_comp = gff3Comps()['product'], width = 10
    ):

    colors = ['#000000', '#010067', '#d5ff00', '#ff0056', '#9e008e', '#0e4ca1', '#ffe502', '#005f39', '#00ff00', '#95003a', '#ff937e', '#a42400', '#001544', '#91d0cb', '#620e00', '#6b6882', '#0000ff', '#007db5', '#6a826c', '#00ae7e', '#c28c9f', '#be9970', '#008f9c', '#5fad4e', '#ff0000', '#ff00f6', '#ff029d', '#683d3b', '#ff74a3', '#968ae8', '#98ff52', '#a75740', '#01fffe', '#ffeee8', '#fe8900', '#bdc6ff', '#01d0ff', '#bb8800', '#7544b1', '#a5ffd2', '#ffa6fe', '#774d00', '#7a4782', '#263400', '#004754', '#43002c', '#b500ff', '#ffb167', '#ffdb66', '#90fb92', '#7e2dd2', '#bdd393', '#e56ffe', '#deff74', '#00ff78', '#009bff', '#006401', '#0076ff', '#85a900', '#00b917', '#788231', '#00ffc6', '#ff6e41', '#e85ebe']


    svg_path = svg_dir + re.sub(r'\.gf[^\.]+$', '.svg', os.path.basename(gff_path))
    if not svg_path.endswith('.svg'):
        svg_path += '.svg'

    gff = gff2list(gff_path)
    seqs, lens = {}, {}
    count = len(product_dict) - 1
    for i in gff:
        if i['seqid'] not in seqs:
            seqs[i['seqid']], lens[i['seqid']] = [], []
        if i['type'] in {'gene'}:
            prod = re.search(prod_comp, i['attributes'])
            if prod is not None:
                product = prod[1].lower()
                if product == 'hypothetical protein' or not product:
                    product = 'hypothetical protein'
                elif product not in product_dict:
                    try:
                        product_dict[product] = colors[count]
                        count += 1
                    except IndexError:
                        count = 0
                        product_dict[product] = colors[count]
                        count += 1
            else:
                product = 'hypothetical protein'
#            print(i['start'],i['end'],i['strand']+'1',product_dict[product],product[:100], flush = True)
            seqs[i['seqid']].append([
                int(i['start']), int(i['end']),
                int(i['strand'] + '1'), product_dict[product],
                product[:100]])
            lens[i['seqid']].append( int(i['start']) )
            lens[i['seqid']].append( int(i['end']) )

    for seq in seqs:
        lens[seq].sort()
        start, features = lens[seq][0], []
        for i in seqs[seq]:
            i[0] -= start
            i[1] -= start
            if i[-1] == 'hypothetical protein':
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
        record = GraphicRecord(sequence_length=length, features = features)
        ax, _ = record.plot(figure_width=width, strand_in_label_threshold=7)
        ax.figure.savefig( svg_path, bbox_inches = 'tight'  )
        

    return product_dict 


def main(
    gff_path, out_dir, product_dict = {'hypothetical protein': '#ffffcc'}, 
    width = 10, prod_comp = gff3Comps()['product']
    ):

    gff_path, out_dir = formatPath( gff_path ), formatPath( out_dir )
    product_dict = gff2svg( gff_path, out_dir, product_dict, width = width )
    
    return product_dict


if __name__ == "__main__":

    parser = argparse.ArgumentParser( description = 'Converts .gff3(s) to .svg' )
    parser.add_argument( '-g', '--gff' )
    parser.add_argument( '-i', '--input', help = 'New line delimited list of gffs' )
    parser.add_argument( '-w', '--width', type = int, default = 10, 
        help = '.svg width; DEFAULT: 10' )
    parser.add_argument( '-o', '--output', help = 'Optional output directory' )
    args = parser.parse_args()

    if args.input:
        if args.output:
            out_dir = formatPath(args.output)
        else:
            out_dir = formatPath(os.path.dirname(args.input))
        gffs = file2list(args.input)
        product_dict = main(gffs[0], out_dir )
        for gff in gffs[1:]:
            product_dict = main(gff, out_dir, product_dict, args.width)
    else:
        if args.output:
            out_dir = formatPath(args.output)
        else:
            out_dir = os.path.dirname(formatPath(args.gff))
        main( args.gff, out_dir, width = args.width )

    sys.exit(0)
