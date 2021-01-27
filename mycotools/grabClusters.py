#! /usr/bin/env python3

from mycotools.lib.kontools import eprint, formatPath, file2list
from mycotools.lib.dbtools import masterDB, db2df
from mycotools.lib.fastatools import gff2dict, fasta2dict, dict2fasta, dict2gff, grabGffAcc
import sys, os, re, argparse, pandas as pd, multiprocessing as mp


def grabFiles( accession, db ):

    ome = re.search( r'(.*?)_', accession )
    if ome is not None:
        ome = ome[1]
        if not pd.isnull( db.at[ome, 'gff3'] ):
            gff = formatPath( '$MYCOGFF3/' + db.at[ome, 'gff3'] )
        else:
            print( 'ERROR: ' + ome + ' gff not detected' )
        if not pd.isnull( db.at[ome, 'proteome'] ):
            prot = formatPath( '$MYCOFAA/' + db.at[ome, 'proteome'] )
        else:
            prot = None
    else:
        print( 'ERROR: ' + accession + ' ome not detected' )
        gff, prot = None, None
        
    return gff, prot


def prepGffOutput( hit_list, gff_path, cpu = 1):

    gff_list = gff2dict( gff_path )
    mp_cmds = []
    for hit in hit_list:
        mp_cmds.append( [gff_list, hit] )
    with mp.get_context('spawn').Pool( processes = cpu ) as pool:
        gff_list_strs = pool.starmap(
            grabGffAcc, mp_cmds
        )
    gff_strs = [ dict2gff(x) for x in gff_list_strs ]

    return '\n'.join( gff_strs )


def prepFaaOutput( hit_list, proteome_path ):

    proteome_dict = fasta2dict( proteome_path )
    clus_fa = {}
    for hit in hit_list:
        clus_fa[hit] = proteome_dict[hit]

    return dict2fasta( clus_fa )


def compileCDS( gff_dict ):

    gene = False
    for entry in gff_dict:
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
        for entry in [ x for x in gff_dict if x['type'] == 'mRNA' ]:
            gene = re.search( r'Alias=([^;]*)', entry['attributes'] )[1]
            mRNA = re.search( r'ID=([^;]*)', entry['attributes'] )[1]
            mRNA_dict[mRNA] = gene


    cds_dict = {}
    for entry in gff_dict:
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


def prepOutput( cds_dict, accession, plusminus ):

    out_index = {}
    for seqid in cds_dict:
        out_index[ seqid ] = []
        if accession in set(cds_dict[ seqid ]):
            index = cds_dict[ seqid ].index( accession )
            if index - plusminus < 0:
                lower = 0
            else:
                lower = int(index - plusminus)
            upper = int( index + plusminus ) + 1
            for index in cds_dict[ seqid ][ lower:upper ]:
                out_index[ seqid ].append( index )
            break    

    return out_index


def main( gff_path, accession, plusminus = 10 ):

    gff = gff2dict( gff_path )
    cds_dict = compileCDS( gff )
    out_index = prepOutput( cds_dict, accession, plusminus )
 
    return out_index


if __name__ == '__main__':

    mp.set_start_method('spawn')
    parser = argparse.ArgumentParser( description = 'Abstracts clusters from accession(s)' )
    parser.add_argument( '-i', '--input', help = 'New line delimitted file of accessions' )
    parser.add_argument( '-a', '--accession', help = 'Input accession' )
    parser.add_argument( '-p', '--plusminus', default = 10, type = int,
        help = 'Genes +/- from accession. DEFAULT: 10' )
    parser.add_argument( '-o', '--output', action = 'store_true',
        help = 'Output cluster fasta(s) and gff(s)' )
    parser.add_argument( '-g', '--gff', help = 'Input GFF file' )
    parser.add_argument( '-f', '--faa', help = 'Input protein fasta file' )
    parser.add_argument( '--cpu', type = int, default = 1, help = 'CPUs, DEFAULT: 1' )
    args = parser.parse_args()

    if args.cpu < mp.cpu_count():
        cpu = args.cpu
    else:
        cpu = mp.cpu_count()
 
    if args.input:
        accessions = file2list( formatPath(args.input) )
    elif args.accession:
        accessions = [ args.accession ]
    else:
        print('\nERROR: requires input or accession')
        sys.exit( 1 )

    out_indices = {}
    if args.gff: 
#        if not args.faa:
 #           print('\nERROR: no proteome')
  #          sys.exit( 3 )
        for accession in accessions:
            out_indices[accession] = main( formatPath(args.gff), accession, args.plusminus )
    else:
        db = db2df( formatPath( masterDB() ) ).set_index('internal_ome')
        for accession in accessions:
            gff, prot = grabFiles( accession, db )
            if gff:
                out_indices[ accession ] = main( gff, accession, args.plusminus )

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
                    with open( accession + '_clus.gff3', 'w' ) as out:
                        out.write( gff_str )
                    if prot:
                        prot_str = prepFaaOutput(
                            out_indices[accession][hit], prot
                            )
                        with open( accession + '_clus.aa.fa', 'w' ) as out:
                            out.write( prot_str )
    else:
        for accession in out_indices:
            print('\n' + accession + ' cluster +/- ' + str(args.plusminus))
            for hit in out_indices[accession]:
                for index in out_indices[accession][ hit ]:
                    print( index )
        print()
 
    sys.exit( 0 )
