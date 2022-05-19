#! /usr/bin/env python3

from mycotools.lib.kontools import eprint, formatPath, file2list
from mycotools.lib.dbtools import masterDB, db2df
from mycotools.lib.biotools import gff2list, fa2dict, dict2fa, list2gff, gff3Comps
from mycotools.acc2gff import grabGffAcc
import sys, os, re, argparse, pandas as pd, multiprocessing as mp


def grabFiles( ome ):
    gff = formatPath( '$MYCOGFF3/' + ome + '.gff3' )
    if not os.path.isfile(gff):
        print( 'ERROR: ' + ome + ' gff not detected' , flush = True)
    prot = formatPath( '$MYCOFAA/' + ome + '.aa.fa' )
    return gff, prot


def prepGffOutput( hit_list, gff_path, cpu = 1):

    gff_list = gff2list( gff_path )
    mp_cmds = []
    for hit in hit_list:
        mp_cmds.append( [gff_list, hit] )
    with mp.get_context('spawn').Pool( processes = cpu ) as pool:
        gff_list_strs = pool.starmap(
            grabGffAcc, mp_cmds
        )
    gff_strs = [ list2gff(x) for x in gff_list_strs ]

    return '\n'.join( gff_strs )


def prepFaaOutput( hit_list, proteome_path ):

    proteome_dict = fa2dict( proteome_path )
    clus_fa = {}
    for hit in hit_list:
        clus_fa[hit] = proteome_dict[hit]

    return dict2fa( clus_fa )


def compileCDS( gff_list ):

    gene = False
    for entry in gff_list:
        if entry['type'] == 'CDS':
            break
        else:
            entry = None


    # need to make this predict in a better way
    prot_comps = [ 
    r'Alias=([^;]*)', r';Name=(.*?);', r'alias ([^;]*)', r'protein_id\=(.*?;|.*?$)', 
    r'proteinId (\d+)', r'gene_id "(.*?)"', r'Parent=(.*?)$'
    ]
    for tag in prot_comps:
        if re.search( tag, entry['attributes'] ):
            prot_comp = re.compile( tag )
            if tag == prot_comps[-1]:
                gene = True
            break

    if gene:
        mRNA_dict = {}
        for entry in [ x for x in gff_list if x['type'] == 'mRNA' ]:
            gene = re.search( r'Alias=([^;]*)', entry['attributes'] )[1]
            mRNA = re.search( r'ID=([^;]*)', entry['attributes'] )[1]
            mRNA_dict[mRNA] = gene


    cds_dict = {}
    for entry in gff_list:
        if entry['type'] == 'CDS':
            if entry['seqid'] not in cds_dict:
                cds_dict[ entry['seqid'] ] = {}
            prot = prot_comp.search( entry['attributes'] )[1]
            if gene:
                prot = mRNA_dict[ prot ]
            if prot not in cds_dict[ entry['seqid'] ]:
                cds_dict[ entry['seqid'] ][ prot ] = []
            cds_dict[ entry['seqid'] ][ prot ].extend( [int( entry['start']), int( entry['end'] )] )
            
    for seqid in cds_dict:
        for prot in cds_dict[ seqid ]:
            cds_dict[ seqid ][ prot ].sort()
        cds_dict[seqid] = sorted( cds_dict[seqid].keys(), key = lambda k: cds_dict[seqid][k][0] )

    return cds_dict


def compileCDS_mycotools(gff_list):
    '''
    Inputs the gff_list and organism ome code. Compiles CDS entries for loci parsing and outputs:
    cds_dict = {contig: {protein: [[CDSstart, CDSstop]] } }
    Sorts the cds_dict for each protein from smallest to largest
    '''

    cds_dict, fail = {}, False
    for entry in gff_list:
        if entry['type'] == 'CDS':
            if entry['seqid'] not in cds_dict:
                cds_dict[entry['seqid']] = {}
            prot_prep_i0 = entry['attributes'].index(';Alias=')
            try:
                prot_prep_i1 = entry['attributes'][prot_prep_i0+1:].index(';') + prot_prep_i0+1
            except ValueError:
                prot_prep_i1 = len(entry['attributes'])
            prot = entry['attributes'][prot_prep_i0:prot_prep_i1].replace(';Alias=','')
            if prot not in cds_dict[entry['seqid']] and prot:
                cds_dict[entry['seqid']][prot] = []
            elif not prot:
                print(entry['attributes'], prot_prep_i0, prot_prep_i1)
                if not fail:
                    print('\tWARNING: proteins in gff with no Alias', flush = True)
                fail = True
                continue
            cds_dict[entry['seqid']][prot].extend([int(entry['start']), int(entry['end'])])

    for seqid in cds_dict:
        for prot in cds_dict[seqid]:
            cds_dict[seqid][prot].sort()
        cds_dict[seqid] = list(sorted(cds_dict[seqid].keys(), key = lambda k: cds_dict[seqid][k][0]))

    return cds_dict


def prepOutput( cds_dict, accession, plusminus ):

    out_index = []
    for seqid in cds_dict:
        if accession in set(cds_dict[ seqid ]):
            index = cds_dict[ seqid ].index( accession )
            if index - plusminus < 0:
                lower = 0
            else:
                lower = int(index - plusminus)
            upper = int( index + plusminus ) + 1
            for index in cds_dict[ seqid ][ lower:upper ]:
                out_index.append( index )
            break    

    return out_index


def main( gff_list, accessions, plusminus = 10, mycotools = False, geneGff = False ):

    out_indices = {}
    if mycotools:
        cds_dict = compileCDS_mycotools(gff_list)
    else:
        cds_dict = compileCDS( gff_list )
    for accession in accessions:
        out_indices[accession] = prepOutput( cds_dict, accession, plusminus )

    out_indices = {k: v for k, v in out_indices.items() if v}
    if geneGff:
        geneGffs_prep = {acc: {} for acc in out_indices}
        gene_sets = {acc: set(genes) for acc, genes in out_indices.items()}
        for entry in gff_list:
            if 'RNA' in entry['type']:
                try:
                    gene = re.search(gff3Comps()['Alias'], entry['attributes'])[1]
                    for acc, genes in gene_sets.items():
                        if gene in genes:
                            geneGffs_prep[acc][gene] = entry
#                            break # if one gene is in multiple loci it needs to show up
                except TypeError: # no alias
                    pass

        geneGffs = {}
        for acc in out_indices:
            geneGffs[acc] = []
            for gene in out_indices[acc]:
                geneGffs[acc].append(geneGffs_prep[acc][gene])
        return out_indices, geneGffs
    return out_indices


def mycotools_main(db, accessions, plusminus = 10, cpus = 1):

    acc_dict = {}
    for acc in accessions:
        ome = acc[:acc.find('_')]
        if ome not in acc_dict:
            acc_dict[ome] = []
        acc_dict[ome].append(acc)

    cmds = [
        [gff2list(grabFiles(ome)[0]), accs, args.plusminus, True] for ome, accs in acc_dict.items()
        ]
    with mp.get_context('spawn').Pool(processes = cpu) as pool:
        acc_res = pool.starmap(main, cmds)

    out_indices = {}
    for res in acc_res:
        out_indices = {**out_indices, **res}

    return out_indices


if __name__ == '__main__':

    mp.set_start_method('spawn')
    parser = argparse.ArgumentParser( description = 'Extracts loci from accession(s)' )
    parser.add_argument( '-i', '--input', help = 'File of accessions' )
    parser.add_argument( '-a', '--accession' )
    parser.add_argument( '-p', '--plusminus', default = 10, type = int,
        help = 'Genes +/- from accession. DEFAULT: 10' )
    parser.add_argument( '-o', '--output', action = 'store_true',
        help = 'Output locus fasta(s) and gff(s)' )
    parser.add_argument( '-g', '--gff', help = 'Input GFF file' )
    parser.add_argument( '-f', '--faa', help = 'Input protein fasta file' )
    parser.add_argument( '-s', '--sep', help = 'Separator for input file.', default = '\n')
    parser.add_argument( '--cpu', type = int, default = 1, help = 'CPUs, DEFAULT: 1' )
    args = parser.parse_args()

    if args.cpu < mp.cpu_count():
        cpu = args.cpu
    else:
        cpu = mp.cpu_count()
    args.sep = args.sep.replace("'",'').replace('"','')
 
    if args.input:
        accessions = file2list( formatPath(args.input), sep = args.sep )
    elif args.accession:
        accessions = [ args.accession ]
    else:
        print('\nERROR: requires input or accession', flush = True)
        sys.exit( 1 )

    out_indices = {}
    if args.gff: 
        gff = gff2list( formatPath(args.gff) )
#        if not args.faa:
 #           print('\nERROR: no proteome', flush = True)
  #          sys.exit( 3 )
        out_indices = main( gff, accession, args.plusminus )
    else:
        db = db2df( formatPath( masterDB() ) ).set_index('internal_ome')
        out_indices = mycotools_main(db, accessions, plusminus = 10, cpus = args.cpus)

#        for i in acc_res:
 #           out_indices[i[1]] = i[0] 

    if args.output:
        for accession in out_indices:
            if args.gff:
                gff = formatPath(args.gff)
                if args.faa:
                    prot = formatPath(args.faa)
                else:
                    prot = None
            else:
                gff, prot = grabFiles( accession, db )
            for hit in out_indices[accession]:
                if len( out_indices[accession][hit] ) > 0:
                    gff_str = prepGffOutput( 
                        out_indices[accession][hit], gff, cpu = args.cpu
                    )
                    with open( accession + '.locus.gff3', 'w' ) as out:
                        out.write( gff_str )
                    if prot:
                        prot_str = prepFaaOutput(
                            out_indices[accession][hit], prot
                            )
                        with open( accession + '.locus.aa.fa', 'w' ) as out:
                            out.write( prot_str )
    else:
        for accession in out_indices:
#            print('\n' + accession + ' locus +/- ' + str(args.plusminus), flush = True)
            for hit in out_indices[accession]:
                for index in out_indices[accession][ hit ]:
                    print( index , flush = True)
            print(flush = True)
 
    sys.exit( 0 )
