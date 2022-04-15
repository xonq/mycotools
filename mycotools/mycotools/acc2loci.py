#! /usr/bin/env python3

from mycotools.lib.kontools import eprint, formatPath, file2list
from mycotools.lib.dbtools import masterDB, db2df
from mycotools.lib.biotools import gff2list, fa2dict, dict2fa, list2gff, grabGffAcc
import sys, os, re, argparse, pandas as pd, multiprocessing as mp


def grabFiles( accession, db ):

    ome = re.search( r'(.*?)_', accession )
    if ome is not None:
        ome = ome[1]
        if not pd.isnull( db.at[ome, 'gff3'] ):
            gff = formatPath( '$MYCOGFF3/' + db.at[ome, 'gff3'] )
        else:
            print( 'ERROR: ' + ome + ' gff not detected' , flush = True)
        if not pd.isnull( db.at[ome, 'proteome'] ):
            prot = formatPath( '$MYCOFAA/' + db.at[ome, 'proteome'] )
        else:
            prot = None
    else:
        print( 'ERROR: ' + accession + ' ome not detected' , flush = True)
        gff, prot = None, None
        
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


def compileCDS_mycotools(gff_dict, ome):
    '''
    Inputs the gff_dict and organism ome code. Compiles CDS entries for loci parsing and outputs:
    cds_dict = {contig: {protein: [[CDSstart, CDSstop]] } }
    Sorts the cds_dict for each protein from smallest to largest
    '''

    cds_dict, fail = {}, False
    for entry in gff_dict:
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
                    print('\tWARNING: ' + ome + ' has proteins in gff with no Alias', flush = True)
                fail = True
                continue
            cds_dict[entry['seqid']][prot].extend([int(entry['start']), int(entry['end'])])

    for seqid in cds_dict:
        for prot in cds_dict[seqid]:
            cds_dict[seqid][prot].sort()
        cds_dict[seqid] = list(sorted(cds_dict[seqid].keys(), key = lambda k: cds_dict[seqid][k][0]))

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


def main( gff_dict, accession, plusminus = 10 ):

    cds_dict = compileCDS( gff_dict )
    out_index = prepOutput( cds_dict, accession, plusminus )
 
    return out_index

def mycotools_main(gff_path, accession, plusminus = 10):

    if gff_path:
        gff = gff2list( gff_path )
        ome = os.path.basename(gff_path).replace('.gff3','')
        cds_dict = compileCDS_mycotools(gff, ome)
        accession_list = prepOutput(cds_dict, accession, plusminus)

    return accession_list, accession


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
        for accession in accessions:
            out_indices[accession] = main( gff, accession, args.plusminus )
    else:
        db = db2df( formatPath( masterDB() ) ).set_index('internal_ome')
        cmds = []
        for accession in accessions:
            cmds.append([grabFiles(accession,db)[0], accession, args.plusminus])
        with mp.get_context('spawn').Pool(processes = cpu) as pool:
            acc_res = pool.starmap(mycotools_main, cmds)
        for i in acc_res:
            out_indices[i[1]] = i[0] 

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
            print('\n' + accession + ' locus +/- ' + str(args.plusminus), flush = True)
            for hit in out_indices[accession]:
                for index in out_indices[accession][ hit ]:
                    print( index , flush = True)
        print(flush = True)
 
    sys.exit( 0 )
