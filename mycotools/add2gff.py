#! /usr/bin/env python3
"""
Add a gff entry to an existing gff
Zachary Konkel
2022 / 11 / 01
"""

# NEED a protein_id import option

import os
import re
import sys
import argparse
from collections import defaultdict
from mycotools.lib.kontools import sys_start, format_path, eprint, mkOutput
from mycotools.lib.biotools import gff2list, list2gff, gff3Comps, \
    gff2Comps, gtfComps
from mycotools.lib.dbtools import mtdb, primaryDB
from mycotools.utils.curGFF3 import rename_and_organize


def determine_version(toadd_gff, ome = None):
    """determine gff version for regex compilations"""
    for entry in toadd_gff:
        if entry['type'] == 'gene':
            if re.search(gff3Comps()['id'], entry['attributes']):
                return toadd_gff, gff3Comps()
            elif re.search(gtfComps()['id'], entry['attributes']):
                return toadd_gff, gtfComps()
            elif re.search(gff2Comps()['id'], entry['attributes']):
                return toadd_gff, gff2Comps()
        elif entry['type'] == 'start_codon': # get this shit out
            from mycotools.utils.gtf2gff3 import main as curAnn
            if not ome:
                eprint('\nOme required for gtf input', flush = True)
                sys.exit(1)
            new_gff = curAnn(toadd_gff, ome)[0]
            return new_gff, gff3Comps()
    else: # finished the for loop and no version detected
        eprint('\nCould not determine gff version', flush = True)
        sys.exit(4)


def id_mtdb_accs(gff):
    """identify mtdb accessions (<OME>_mtdb###) to ID initial index for new
    additions"""
    mtdb_accs = []
    for entry in gff:
        try:
            mtdb_acc = re.search(gff3Comps()['Alias'], entry['attributes'])[1]
        except TypeError: # no alias
            continue
        if '_manual' in mtdb_acc: # mtdb accession explicit
            mtdb_accs.append(int(mtdb_acc[mtdb_acc.find('_manual')+7:]))
    if mtdb_accs:
        max_mtdb_acc = max(mtdb_accs)
    else:
        max_mtdb_acc = 0
    ome = mtdb_acc[:mtdb_acc.find('_')]
    return max_mtdb_acc, ome


def parse_toadd(toadd_gff, comps, ome, mtdb_acc = 0):
    """Compile genes and hierarchy of toadd gff"""

    genes = {}
    for entry in toadd_gff:
        if entry['type'] == 'gene':
            id_ = re.search(comps['id'], entry['attributes'])[1]
            genes[id_] = {'gene': entry, 'rna': [], 'cds': [], 'exon': []}
        elif entry['type'].lower() in {'cds', 'exon'}:
            try: # attempt to grab parent, if not assume it is organized
                if 'par' in comps:
                    par = re.search(comps['par'], entry['attributes'])[1]
                else:
                    par = id_
            except IndexError:
                par = id_
            start, stop = sorted([int(entry['start']),
                                  int(entry['end'])])
            genes[id_][entry['type'].lower()].append([start, stop])
        elif 'rna' in entry['type'].lower():
            start, stop = sorted([int(entry['start']),
                                  int(entry['end'])])
            try: # attempt to grab parent, if not assume it is organized
                par = re.search(comps['par'], entry['attributes'])[1]
            except IndexError:
                par = id_
            genes[id_]['rna'].append([entry['type'], start, stop])
    
    new_gff = []
    for gene_key, gene_dict in genes.items():
        mtdb_acc += 1
        gene = gene_dict['gene']
        start, stop = sorted([int(gene_dict['gene']['start']), 
                              int(gene_dict['gene']['end'])])
        scaf, source, strand = gene['seqid'], gene['source'], gene['strand']
        score, phase = gene['score'], gene['phase']
        gene_acc = ome + '_manual' + str(mtdb_acc)
        if len(gene_dict['rna']) > 1:
            eprint('\nAlternately spliced loci currently not supported')
            sys.exit(2)

        gene_dict['cds'] = sorted(gene_dict['cds'], key = lambda x: x[0])
        gene_dict['exon'] = sorted(gene_dict['exon'], key = lambda x: x[0])
        if not gene_dict['rna']: # if no mRNA
            if gene_dict['cds']: # but there is a protein
                gene_dict['rna'] = [['mRNA', start, stop]]
        if not gene_dict['exon']: # if no exon
            if gene_dict['cds']: # but there is a protein
                gene_dict['exon'] = [x for x in gene_dict['cds']]
            elif gene_dict['rna']: # but there is an RNA
                gene_dict['exon'] = [[gene_dict['rna'][0][1],
                                     gene_dict['rna'][0][2]]]

        gene_entry = {'seqid': scaf, 'source': source, 'type': 'gene',
                      'start': start, 'end': stop, 'score': score,
                      'strand': strand, 'phase': phase,
                      'attributes': f'ID=gene_{gene_acc}'}
        if gene_dict['cds']:
            gene_entry['attributes'] += f';protein_id={gene_acc}'
        if gene_dict['rna']:
            gene_entry['attributes'] += f';transcript_id={gene_acc}'
        gene_entry['attributes'] += f';Alias={gene_acc}'
        new_gff.append(gene_entry)

        for rna in gene_dict['rna']:
            rna_entry = {'seqid': scaf, 'source': source, 'type': rna[0],
                         'start': start, 'end': stop, 'score': '.',
                         'strand': strand, 'phase': '.',
                         'attributes': f'ID={rna[0]}_{gene_acc};' \
                                     + f'Parent=gene_{gene_acc};' \
                                     + f'transcript_id={gene_acc}'}
            if gene_dict['cds']:
                rna_entry['attributes'] += f';protein_id={gene_acc}'
            rna_entry['attributes'] += f';Alias={gene_acc}'
            new_gff.append(rna_entry)

        for i, cds in enumerate(gene_dict['cds']):
            rna_type = rna[0]
            cds_entry = {'seqid': scaf, 'source': source, 'type': 'CDS',
                         'start': cds[0], 'end': cds[1], 'score': '.',
                         'strand': strand, 'phase': '.', 
                         'attributes': f'ID=cds_{gene_acc}.{i+1};' \
                                     + f'Parent={rna_type}_{gene_acc};' \
                                     + f'protein_id={gene_acc};' \
                                     + f'transcript_id={gene_acc};' \
                                     + f'Alias={gene_acc}'}
            new_gff.append(cds_entry)

        for i, exon in enumerate(gene_dict['exon']):
            rna_type = rna[0]
            exon_entry = {'seqid': scaf, 'source': source, 'type': 'exon',
                          'start': exon[0], 'end': exon[1], 'score': '.',
                          'strand': strand, 'phase': '.',
                          'attributes': f'ID=exon_{gene_acc}.{i+1};' \
                                      + f'Parent={rna_type}_{gene_acc};' \
                                      + f'protein_id={gene_acc};' \
                                      + f'transcript_id={gene_acc};' \
                                      + f'Alias={gene_acc}'}
            new_gff.append(exon_entry)

    return new_gff


def compile_mtdb_scaf(scafs, gff):
    old_coords, rev_coords = defaultdict(dict), defaultdict(dict)
    for entry in gff:
        if entry['type'] == 'gene':
            if entry['seqid'] in scafs:
                mtdb_acc = re.search(gff3Comps()['Alias'], 
                                     entry['attributes'])[1]
                coord_tup = tuple(sorted([entry['start'], entry['end']]))
                old_coords[entry['seqid']][mtdb_acc] = \
                    sorted([entry['start'], entry['end']])
                rev_coords[entry['seqid']][coord_tup] = mtdb_acc

    return rev_coords, old_coords 


def add_to_mtdb_gff(curadd_gff, addto_gff, ome, replace = False):
    update_scafs = set([x['seqid'] for x in curadd_gff])
    if addto_gff:
        ref_coords, old_coords = compile_mtdb_scaf(update_scafs, 
                                                   addto_gff)

    update = defaultdict(dict)
    if replace:
        for entry in curadd_gff:
            if entry['type'] == 'gene':
                start, stop = entry['start'], entry['end']
                seqid = entry['seqid']
                new_acc = re.search(gff3Comps()['Alias'], entry['attributes'])[1]
                if seqid in ref_coords:
                    for ref_coord, ref_acc in ref_coords[seqid].items():
                        ref_start, ref_stop = ref_coord
                        if (start >= ref_start and start <= ref_stop) or \
                           (stop >= ref_start and stop <= ref_stop): # if overlap
                            update[entry['seqid']][ref_acc] = new_acc # what about multiple gene
                            # updates? e.g. fusions? this would delete
                            eprint(ref_acc + '\t->\t' + new_acc, flush = True)

        out_gff = []
        for entry in addto_gff:
            if entry['seqid'] in update:
                mtdb_acc = re.search(gff3Comps()['Alias'],
                                     entry['attributes'])[1]
                if mtdb_acc in update[entry['seqid']]:
                    continue
            out_gff.append(entry)
    else:
        out_gff = addto_gff

    out_gff.extend(curadd_gff)
    return out_gff

def main(toadd_gff, addto_gff = [], ome = None, replace = False):
    if addto_gff:
        max_mtdb_acc, ome = id_mtdb_accs(addto_gff)
    else:
        max_mtdb_acc = 0

    toadd_gff, comps = determine_version(toadd_gff, ome = ome)
    curadd_gff = parse_toadd(toadd_gff, comps, ome, max_mtdb_acc)
    out_gff = add_to_mtdb_gff(curadd_gff, addto_gff, ome, replace)
    final_gff = rename_and_organize(out_gff)
    return final_gff, ome


def prep_mtdb_update(new_gff, ome, db):
    from mycotools.predb2mtdb import main as predb2mtdb
    from mycotools.lib.dbtools import mtdb, primaryDB
    out_dir = mkOutput(format_path(os.getcwd()), 'add2gff')
    wrk_dir = out_dir + 'working/'
    if not os.path.isdir(wrk_dir):
        os.mkdir(out_dir + 'working/')
    if not os.path.isdir(wrk_dir + 'fna/'):
        os.mkdir(out_dir + 'working/fna/')
    if not os.path.isdir(wrk_dir + 'gff3/'):
        os.mkdir(out_dir + 'working/gff3/')
    if not os.path.isdir(wrk_dir + 'faa/'):
        os.mkdir(out_dir + 'working/faa/')

    new_gff_path = out_dir + 'working/gff3/' + ome + '.new.gff3'
    with open(new_gff_path, 'w') as out:
        out.write(list2gff(new_gff))

    # make predb output
    predb = {'assembly_acc': [db[ome]['assembly_acc']],
             'previous_ome': [ome], 'genus': [db[ome]['genus']],
             'species': [db[ome]['species']], 
             'strain': [db[ome]['strain']],
             'version': [db[ome]['version']], 
             'biosample': [db[ome]['biosample']],
             'assemblyPath': [db[ome]['fna']], 
             'gffPath': [new_gff_path],
             'source': [db[ome]['source']],
             'useRestriction (yes/no)': [db[ome]['published']],
             'published': [db[ome]['published']]}
    new_db, failed = predb2mtdb(predb, db.reset_index(), out_dir + 'working/', 
                              exit = True, remove = False)

    update_db_path = out_dir + 'add2gff.mtdb'
    count = 1
    while os.path.isfile(update_db_path):
        update_db_path = re.sub(r'_\d+$', '', update_db_path)
        update_db_path += f'_{count}'
        count += 1
    new_db.df2db(update_db_path)


def import_toadd_gff(toadd_file):
    try:
        toadd_gff = gff2list(toadd_file)
    except IndexError: # not completely a gff
        toadd_gff_str = ''
        with open(toadd_file, 'r') as raw:
            for line in raw:
                if len(line.split('\t')) == 9:
                    toadd_gff_str += line 
        toadd_gff = gff2list(toadd_gff_str, path = False)
    return toadd_gff

def cli():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help = 'To add gff', required = True)
    parser.add_argument('-a', '--addto', help = 'Add to gff')
    parser.add_argument('-o', '--ome', help = 'Input ome if no [-a]')
    parser.add_argument('-r', '--replace', action = 'store_true',
        help = 'Replace overlapping accession(s)')
    parser.add_argument('-u', '--update', action = 'store_true',
        help = '[-a] Prepare output for mtdb update')
    parser.add_argument('-d', '--mtdb', default = primaryDB())
    args = parser.parse_args()

#    usage = 'Add gff to an existing mtdb gff.\n' \
 #         + 'add2gff.py <TOADD_GFF> <ADDTO_GFF> <REPLACE ACCS [y|N]>'
  #  args = sys_start(sys.argv[1:], usage, 2,
   #                  files = sys.argv[2:4])

    if not args.addto:
        if args.update:
            eprint('\nERROR: -u requires -a', flush = True)
            sys.exit(7)
        if not args.ome:
            eprint('\nERROR: need -a or -o', flush = True)
            sys.exit(5)
        else:
            ome = args.ome
            addto_gff = [] # empty gff list
    else:
        ome = None
        addto_file = format_path(args.addto)
        if not os.path.isfile(addto_file):
            eprint('\nERROR: -a does not exist', flush = True)
            sys.exit(6)
        addto_gff = gff2list(addto_file)

    toadd_file = format_path(args.input)
    if not os.path.isfile(toadd_file):
        eprint('\nERROR: -i does not exist', flush = True)
        sys.exit(7)

    toadd_gff = import_toadd_gff(toadd_file)

    replace = bool(args.replace)
    gff_list, ome = main(toadd_gff, addto_gff, ome, replace)
    if args.update:
        db = mtdb(format_path(args.mtdb)).set_index()
        prep_mtdb_update(gff_list, ome, db)
    else:
        print(list2gff(gff_list))
    sys.exit(0)


if __name__ == '__main__':
    cli()
