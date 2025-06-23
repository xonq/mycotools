#! /usr/bin/env python3

import os
import re
import sys
import argparse
import multiprocessing as mp
from itertools import chain
from collections import defaultdict
from mycotools.lib.kontools import eprint, format_path, stdin2str
from mycotools.lib.dbtools import mtdb, primaryDB
from mycotools.lib.biotools import fa2dict, gff2list, gff3Comps
from mycotools.acc2gff import db_main as acc2gff


def col_CDS(
    gff_list, types={"gene", "CDS", "exon", "mRNA", "tRNA", "rRNA", "RNA", "pseudogene"}
):
    """Collect all CDS entries from a `gff` and store them into cds_dict.
    If the gene found in the CDS is in the set of `ids` then add that type
    to the CDS dict for that protein header."""

    cds_dict = defaultdict(dict)
    # parse the gff to acquire the aliases
    for entry in gff_list:
        if entry["type"] in types:
            contig = entry["seqid"]
            try:
                alias = re.search(gff3Comps()["Alias"], entry["attributes"])[1]
            except TypeError:
                eprint("\n\tERROR: could not extract Alias ID from " + gff, flush=True)
                continue
            aliases = alias.split("|")  # to address alternate splicing in gene
            # aliases
            # compile all CDS for each alias
            for alias in aliases:
                ome = alias[: alias.find("_")]
                if contig not in cds_dict[ome]:
                    cds_dict[ome][contig] = defaultdict(list)
                cds_dict[ome][contig][alias].append(entry)

    return cds_dict


def contig2gbk(
    ome,
    row,
    contig,
    contig_dict,
    contig_seq,
    faa,
    product_searches={"product": r"product=([^;]+)"},
    full=False,
):
    """Generate a genbank string for a contig_dict, which contains a dictionary
    of protein keys and list of each gff entry associated with the protein.
    ome: ome_code
    row: mtdb row
    contig: contig name
    contig_dict: {prot: [gff_entries]}
    contig_seq: contig sequence string
    faa: proteome
    product_searches: {prod_id: product regular expression}"""

    seq_coords = []  # a list of all the coordinates from the contig
    for prot, prot_list in contig_dict.items():
        for entry in prot_list:
            t_start, t_end = int(entry["start"]), int(entry["end"])
            seq_coords.append(sorted([t_start, t_end]))

    # sort the coordinates
    seq_coords.sort(key=lambda x: x[0])
    # the starting index
    if full:
        seq_coords.insert(0, [1, 1])
        seq_coords.append([len(contig_seq), len(contig_seq)])

    init_start_test = int(seq_coords[0][0]) - 1
    init_end_test = int(seq_coords[-1][1])
    end_test = init_end_test - init_start_test  # relative end
    start_test = 1  # relative start

    # specify the contig edge - may cause parsers (BiG-SCAPE/antiSMASH) to fail
    edge = (
        "/contig_edge="
        + str(init_start_test)
        + "-;"
        + str(len(contig_seq) - init_end_test)
        + "+"
    )
    name = ome + "_" + contig
    relative_end = seq_coords[-1][1] - seq_coords[0][0]
    gbk = (
        "LOCUS       "
        + name
        + "\n"
        + "DEFINITION  "
        + name
        + "\n"
        + "ACCESSION   "
        + contig
        + "\n"
        + "VERSION     "
        + name
        + "\n"
        + "KEYWORDS    .\nSOURCE    "
        + row["source"]
        + "\n  ORGANISM  "
        + str(row["genus"])
        + "_"
        + str(row["species"])
        + "_"
        + str(row["strain"])
        + "_"
        + ome
        + "\n"
        + "COMMENT     coordinates are relative to adjustment: "
        + str(init_start_test)
        + "\n            .\nFEATURES             Location/Qualifiers\n"
        + "     region          "
        + str(start_test)
        + ".."
        + str(end_test)
        + "\n"
        + "                     "
        + str(edge)
        + "\n"
    )
    #        '     region          1..' + str(relative_end) + '\n' + \
    # this should be
    # related to the start of the contig

    used_aliases = set()  # to address alternate splicing and 1 entry per gene
    for prot, prot_list in contig_dict.items():
        final_coords = ""
        entries = defaultdict(list)
        for entry in prot_list:
            if entry["type"].lower() == "cds":
                entries["CDS"].append(entry)
            elif entry["type"].lower() == "exon":
                entries["exon"].append(entry)
            elif entry["type"].lower() in {"mrna", "rrna", "trna", "rna"}:
                entries["RNA"].append(entry)
            elif entry["type"].lower() in {"gene", "pseudogene"}:
                entries["gene"].append(entry)

        # for each gene, let it be the parent entry
        for entry in entries["gene"]:
            alias = re.search(gff3Comps()["Alias"], entry["attributes"])[1]
            #            eprint(alias)
            if "|" in alias:  # alternately spliced gene
                if alias in used_aliases:
                    continue
                else:
                    used_aliases.add(alias)
            id_ = re.search(gff3Comps()["id"], entry["attributes"])[1]

            products = {}
            # try to acquire the product name
            for prod_id, product_search in product_searches.items():
                try:
                    product = re.search(product_search, entry["attributes"])[1]
                    products[prod_id] = product
                except TypeError:  # no product
                    pass

            # begin process for assigning gene coordinates
            start, end = sorted([int(entry["start"]), int(entry["end"])])
            gbk += "     " + entry["type"] + "                "[: -len(entry["type"])]
            if entry["strand"] == "+":
                gene_coords = (
                    str(start - init_start_test)
                    + ".."
                    + str(end - init_start_test)
                    + "\n                     "
                )
            else:
                gene_coords = (
                    "complement("
                    + str(start - init_start_test)
                    + ".."
                    + str(end - init_start_test)
                    + ")\n                     "
                )
            gbk += gene_coords

            # assign etc information to entry
            gbk += '/locus_tag="' + alias + '"\n                     '
            for prod_id, product in products.items():
                if product:
                    gbk += f'/{prod_id}="' + product + '"\n                     '
            if entry["phase"] != ".":
                gbk += '/phase="' + entry["phase"] + '"\n                     '
            gbk += '/source="' + entry["source"] + '"\n'

        # for each RNA entry
        for entry in entries["RNA"]:
            gbk += "     " + entry["type"] + "                "[: -len(entry["type"])]
            # acquire the products
            products = {}
            for prod_id, product_search in product_searches.items():
                try:
                    product = re.search(product_search, entry["attributes"])[1]
                    products[prod_id] = product
                except TypeError:  # no product
                    pass
            alias = re.search(gff3Comps()["Alias"], entry["attributes"])[1]
            id_ = re.search(gff3Comps()["par"], entry["attributes"])[1]

            # append to the final gene coordinates
            if final_coords:
                gbk += final_coords
            else:  # if there is not alternative splicing
                start, end = sorted([int(entry["start"]), int(entry["end"])])
                if entry["strand"] == "+":
                    gene_coords = (
                        str(start - init_start_test)
                        + ".."
                        + str(end - init_start_test)
                        + "\n                     "
                    )
                else:
                    gene_coords = (
                        "complement("
                        + str(start - init_start_test)
                        + ".."
                        + str(end - init_start_test)
                        + ")\n                     "
                    )
                gbk += gene_coords

            # add the etcs
            gbk += '/locus_tag="' + alias + '"\n                     '

            for prod_id, product in products.items():
                if product:
                    gbk += f'/{prod_id}="' + product + '"\n                     '
            if entry["phase"] != ".":
                gbk += '/phase="' + entry["phase"] + '"\n                     '
            gbk += '/source="' + entry["source"] + '"\n'
            gbk += '                     /transcript_id="' + id_ + '"\n'

        # for each CDS entry
        if entries["CDS"]:
            entry = entries["CDS"][0]
            strand = entry["strand"]
            cds_coords = ""
            if strand == "+":  # sense strand, easy set-up
                for cds in entries["CDS"]:
                    coords = [
                        cds["start"] - init_start_test,
                        cds["end"] - init_start_test,
                    ]
                    max_c, min_c = max(coords), min(coords)
                    cds_coords += str(min_c) + ".." + str(max_c) + ","
            else:  # antisense, prepare the negative coordinate list
                for cds in entries["CDS"]:
                    coords = [
                        cds["start"] - init_start_test,
                        cds["end"] - init_start_test,
                    ]
                    max_c, min_c = max(coords), min(coords)
                    cds_coords += "complement(" + str(min_c) + ".." + str(max_c) + "),"
            cds_coords = cds_coords[:-1]
            cds_coords_list = cds_coords.split(",")
            final_coords = ""
            total_len = 0
            # prepare the coordinates for each set of sequence
            for coord in cds_coords_list[:-1]:
                if total_len + len(coord) < 57 or total_len == 0:
                    final_coords += coord + ","
                    total_len += len(coord) + 1
                else:
                    total_len = len(coord) + 1
                    final_coords += "\n                     " + coord + ","
            if total_len == 0 or total_len + len(cds_coords_list[-1]) < 58:
                final_coords += cds_coords_list[-1]
            else:
                final_coords += "\n                     " + cds_coords_list[-1]
            if len(cds_coords_list) > 1:
                final_coords = "join(" + final_coords + ")"
            alias = re.search(gff3Comps()["Alias"], entry["attributes"])[1]
            id_ = re.search(gff3Comps()["id"], entry["attributes"])[1]

            products = {}
            for prod_id, product_search in product_searches.items():
                try:
                    product = re.search(product_search, entry["attributes"])[1]
                    products[prod_id] = product
                except TypeError:  # no product
                    pass
            gbk += "     " + entry["type"] + "                "[: -len(entry["type"])]
            final_coords += "\n                     "
            gbk += final_coords
            gbk += '/locus_tag="' + alias + '"\n                     '

            for prod_id, product in products.items():
                if product:
                    gbk += f'/{prod_id}="' + product + '"\n                     '
            if entry["phase"] != ".":
                gbk += '/phase="' + entry["phase"] + '"\n                     '
            gbk += '/source="' + entry["source"] + '"\n'

            # acquire the sequence for the protein
            try:
                inputSeq = faa[alias]["sequence"]
            except KeyError:
                continue
            gbk += "                     "
            gbk += '/protein_id="' + alias + '"\n'
            gbk += (
                "                     "
                + '/transl_table=1\n                     /translation="'
            )

            lines = [inputSeq[:45]]
            if len(inputSeq) > 45:
                restSeq = inputSeq[45:]
                n = 59
                toAdd = [restSeq[i : i + n] for i in range(0, len(restSeq), n)]
                lines.extend(toAdd)
            gbk += lines[0] + "\n"
            if len(lines) > 1:
                for line in lines[1 : len(lines) - 1]:
                    gbk += "                     " + line + "\n"
                gbk += "                     " + lines[-1] + '"\n'
            else:
                gbk += '"\n'
            if gbk.endswith('\n"\n'):
                gbk = gbk[:-2] + '                     "\n'
    total_coords = []
    for coord in seq_coords:
        total_coords.extend(coord)
    min_c, max_c = min(total_coords), max(total_coords)
    #    for coordinates in seq_coords:
    assSeq = contig_seq[min_c:max_c]
    seqLines = [assSeq[i : i + 60] for i in range(0, len(assSeq), 60)]
    count = -59

    # format the header
    gbk += "ORIGIN\n"
    for line in seqLines:
        count += 60
        gbk += "{:>9}".format(str(count)) + " "
        seps = [line[i : i + 10] for i in range(0, len(line), 10)]
        for sep in seps:
            gbk += sep.lower() + " "
        gbk = gbk.rstrip() + "\n"

    gbk += "//"
    return gbk


def gen_gbk(
    ome_dict, row, fna, faa, ome, product_searches, full=False, break_contigs=False
):
    """manage generating the GBK for each contiguous sequence"""
    gbk = {}
    if break_contigs:
        for contig, contig_dict in ome_dict.items():
            seq_coords = {}  # a dict of all the coordinates from the contig
            for prot, prot_list in contig_dict.items():
                coords = []
                for entry in prot_list:
                    t_start, t_end = int(entry["start"]), int(entry["end"])
                    coords.append([t_start, t_end])
                new_coords = [min(chain(*coords)), max(chain(*coords))]
                seq_coords[prot] = new_coords

            seq_coords = {
                k: v for k, v in sorted(seq_coords.items(), key=lambda x: x[1][0])
            }
            contig_dict = {k: contig_dict[k] for k in seq_coords}

            breaks, count0 = set(), 0
            for k1, v1 in seq_coords.items():
                if count0:
                    if v1[0] - v0[1] > break_contigs:
                        breaks.add(k0)
                k0 = k1
                v0 = v1
                count0 += 1

            count1 = 0
            if breaks:
                new_contigs = {}
                for k, v in contig_dict.items():
                    new_contigs[k] = v
                    if k in breaks:
                        nc_name = f"{contig}_{count1}"
                        count1 += 1
                        gbk[nc_name] = contig2gbk(
                            ome,
                            row,
                            nc_name,
                            new_contigs,
                            fna[contig]["sequence"],
                            faa,
                            full=full,
                            product_searches=product_searches,
                        )
                        new_contigs = {}
                nc_name = f"{contig}_{count1}"
                gbk[nc_name] = contig2gbk(
                    ome,
                    row,
                    nc_name,
                    new_contigs,
                    fna[contig]["sequence"],
                    faa,
                    full=full,
                    product_searches=product_searches,
                )

            else:
                gbk[contig] = contig2gbk(
                    ome,
                    row,
                    contig,
                    contig_dict,
                    fna[contig]["sequence"],
                    faa,
                    full=full,
                    product_searches=product_searches,
                )
    else:
        for contig, contig_dict in ome_dict.items():
            gbk[contig] = contig2gbk(
                ome,
                row,
                contig,
                contig_dict,
                fna[contig]["sequence"],
                faa,
                full=full,
                product_searches=product_searches,
            )

    return gbk


def ome_main(
    ome,
    gff_lists,
    row,
    product_searches={"product": r"product=([^;]+)"},
    break_contigs=False,
):
    """Create a genbank for the whole genome"""
    gbks = defaultdict(list)
    faa, fna = fa2dict(row["faa"]), fa2dict(row["fna"])
    for key, gffs in gff_lists.items():
        for gff_list in gffs:
            cds_dict = col_CDS(
                gff_list,
                types={
                    "gene",
                    "CDS",
                    "exon",
                    "mRNA",
                    "tRNA",
                    "rRNA",
                    "RNA",
                    "pseudogene",
                },
            )
            gbks[key].append(
                gen_gbk(
                    cds_dict[ome], row, fna, faa, ome, product_searches, break_contigs
                )
            )
    return gbks


def main(
    gff_list,
    db,
    product_searches={"product": r"product=([^;]+)"},
    full=False,
    break_contigs=False,
):
    """Generate a genbank for each inputed gff, its associated fna, and ome"""
    cds_dict = col_CDS(
        gff_list,
        types={"gene", "CDS", "exon", "mRNA", "tRNA", "rRNA", "RNA", "pseudogene"},
    )
    ome_gbks = {}
    for ome, ome_dict in cds_dict.items():
        fna = fa2dict(db[ome]["fna"])
        faa = fa2dict(db[ome]["faa"])
        ome_gbks[ome] = gen_gbk(
            ome_dict,
            db[ome],
            fna,
            faa,
            ome,
            product_searches,
            full=full,
            break_contigs=break_contigs,
        )

    return ome_gbks


def cli():

    parser = argparse.ArgumentParser(
        description="Inputs MTDB gff or accessions, outputs GenBank file"
    )
    parser.add_argument("-a", "--accession", help='"-" for stdin')
    parser.add_argument("-i", "--input", help="File with accessions")
    parser.add_argument(
        "-f", "--full", action="store_true", help="Whole genome; -a must be an ome code"
    )
    parser.add_argument("-g", "--gff", help="GFF3 input")
    parser.add_argument(
        "-p", "--product", help='Product regular expression: DEFAULT: "Product=([^;]+)"'
    )
    parser.add_argument(
        "-o", "--ome", action="store_true", help="Output files by ome code"
    )
    parser.add_argument(
        "-s",
        "--seqid",
        action="store_true",
        help="Output files by sequence ID (contig/scaffold/chromosome)",
    )
    parser.add_argument(
        "-b",
        "--split",
        default=0,
        type=int,
        help="Base pairs to split long breaks between genes",
    )
    parser.add_argument("-d", "--mtdb", help="DEFAULT: master", default=primaryDB())
    parser.add_argument("-c", "--cpu", type=int, default=1)
    args = parser.parse_args()

    # gather accessions from various input styles and split into a list
    if args.input:
        input_file = format_path(args.input)
        with open(input_file, "r") as raw:
            accs_prep = [x.rstrip().split() for x in raw]
        accs = []
        for acc in accs_prep:
            accs.extend(acc)
    # if it is a single accession or coming from standard input, then acquire
    # that
    elif args.accession:
        if args.accession == "-":
            data = stdin2str()
            accs = data.split()
        else:
            if {'"', "'"}.intersection(set(args.accession)):
                args.accession = args.accession.replace('"', "").replace("'", "")
            if "," in args.accession:
                accs = args.accession.split(",")
            elif re.search(r"\s", args.accession):
                accs = args.accession.split()
            else:
                accs = [args.accession]
    # if there is no gff then no input was received
    elif not args.gff:
        raise FileNotFoundError("need input file or accession")

    # only full genomes for args.full, no accessions
    if args.full:
        if any("_" in x for x in accs):  # _ is forbidden from ome codes
            eprint("\nERROR: -f requires omes, not accessions", flush=True)
            sys.exit(3)

    # create a default regular expression for the product name
    if not args.product:
        regex = r"Product=([^;])+"
    else:
        regex = r""
        if args.regex.startswith(("'", '"')):
            args.regex = args.regex[1:]
        if args.regex.endswith(('"', "'")):
            args.regex = args.regex[:-1]
        for char in args.regex:
            regex += char

    # import database and set index
    db = mtdb(format_path(args.mtdb))
    db = db.set_index()

    # various output formats
    gbk_dict = {}
    if args.gff:  # use an inputted gff
        gbk_dict = main(
            gff2list(format_path(args.gff)),
            db,
            {"product": regex},
            full=args.full,
            break_contigs=args.split,
        )
    elif not args.full:  # acquire an inputted gff from the accessions
        gffs = acc2gff(db, accs, args.cpu)
        for ome, gff in gffs.items():
            gbk_dict[ome] = main(gff, db, {"product": regex}, break_contigs=args.split)[
                ome
            ]
    else:  # acquire a gff for each inputted ome
        gffs = {ome: db[ome]["gff3"] for ome in accs}
        for ome, gff_p in gffs.items():
            gbk_dict[ome] = main(
                gff2list(gff_p),
                db,
                {"product": regex},
                args.full,
                break_contigs=args.split,
            )[ome]

    if args.ome:  # output a file for each genome inputted
        for ome, gbks in gbk_dict.items():
            with open(ome + ".acc2gbk.gbk", "w") as out:
                out.write("\n".join([gbk for gbk in gbks]))
    #    elif args.accession: # output a file for each contig inputted
    #       for ome, gbks in gbk_dict.items():
    #          for contig, gbk in gbks.items():
    #             with open(contig + '.gbk', 'w') as out:
    #                out.write(gbk)
    else:  # otherwise just print to stdout
        for ome, gbks in gbk_dict.items():
            for contig, gbk in gbks.items():
                print(gbk, flush=True)


if __name__ == "__main__":
    cli()
    sys.exit(0)
