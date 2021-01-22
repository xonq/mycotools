#! /usr/bin/env python3

from Bio import Entrez
import pandas as pd, numpy as np
import argparse, os, time, subprocess, sys, requests, re
from mycotools.lib.kontools import intro, outro, formatPath, prep_output, eprint
from mycotools.lib.dbtools import log_editor


def prepareFolders( output_path, gff, prot, assem, transcript ):

    file_types = []
    if assem:
        if not os.path.exists(output_path + 'assembly'):
            os.mkdir(output_path + 'assembly')
        file_types.append( 'assembly' )
    if prot:
        if not os.path.exists(output_path + 'proteome'):
            os.mkdir(output_path + 'proteome')
        file_types.append( 'proteome' ) 
    if gff:
        if not os.path.exists(output_path + 'gff3'):
            os.mkdir(output_path + 'gff3')
        file_types.append( 'gff3' )
    if transcript:
        if not os.path.exists(output_path + 'transcript'):
            os.mkdir(output_path + 'transcript')
        file_types.append( 'transcript' )

    return file_types  


def compileLog( output_path, remove = False ):

    ass_prots = {}
    if not os.path.isfile( output_path + 'ncbiDwnld.log' ):
        with open( output_path + 'ncbiDwnld.log', 'w' ) as out:
            out.write('#ome\tbiosample\tassembly\tproteome\tgff3\ttranscript\t' + \
            'assemMD5\tprotMD5\tgff3MD5\ttransMD5')

        # too risky, too many things can go wrong and then users would be in a
        # loop, but necessary for huge downloads
    else:
        with open( output_path + 'ncbiDwnld.log', 'r' ) as raw:
            for line in raw:
                if not line.startswith('#'):
                    data = [x.rstrip() for x in line.split('\t')]
                    while len(data) < 10:
                        data.append('')
                    ass_prots[data[0]] = { 
                        'biosample': str(data[1]), 'assembly': str(data[2]), 
                        'proteome': str(data[3]), 'gff3': str(data[4]), 
                        'transcript': str(data[5]), 'assembly_md5': str(data[6]),
                        'proteome_md5': str(data[7]), 'gff3_md5': str(data[8]),
                        'transcript_md5': str(data[9])
                        }

    return ass_prots


# collects paths to download proteomes and assemblies
def collect_ftps(
    ncbi_df, ass_prots, api_key=0,
    column='Assembly Accession', database="assembly", output_path = ''
    ):

    print('\nAssembling ftp directories ...')
    count = 0

# for each row in the assembly, grab the accession number, form the search term for Entrez, use Entrez,
    for index, row in ncbi_df.iterrows():
        accession = row[column]
#        ome = row['internal_ome']
        if index in ass_prots:
            continue

        ass_prots[str(index)] = {}

        search_term = accession + '[' + column + ']'
        handle = Entrez.esearch(db=database, term=search_term)
        genome_id = Entrez.read(handle)['IdList']
        ID = genome_id[0]

# if more than one genome id is found, this needs to be reported because we're only using the first (shouldn't happen)
        if len(genome_id) > 1:
            print('\t' + accession, 'yields more than one genome ID:')
            for ID in genome_id:
                print('\t\t' + str(ID))
            ID = genome_id[0]
            print('\t\tUsing', str(ID))

# obtain the path from a summary of the ftp directory and create the standard paths for proteomes and assemblies
        esc_count = 0
        while True:
            try:
                esc_count += 1
                handle = Entrez.esummary(db=database, id=ID, report="full")
                record = Entrez.read(handle, validate = False)
                ftp_path = str(record['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_GenBank'])
                break
            except IndexError:
                if esc_count == 20:
                    print('\t\tFailed to retrieve.')
                    sys.exit(69)
                time.sleep(0.1)

        ass_md5, gff_md5, trans_md5, prot_md5, md5s = '', '', '', '', {}
        basename = os.path.basename(ftp_path)
        request = requests.get( 
            (ftp_path + '/md5checksums.txt').replace('ftp://','https://'), timeout = 10
            )
        if request.status_code != 200:
            eprint('\t' + ome + '\tno ' + md5)
            dwnld = 69
        else:
            dwnld = subprocess.call(
                ['curl', ftp_path + '/md5checksums.txt', '-o', output_path + '.tmpmd5'], 
                stdout = subprocess.PIPE, stderr = subprocess.PIPE
                )
            count += 1
            with open( output_path + '.tmpmd5', 'r' ) as raw:
                for line in raw:
                    data = line.rstrip().split('  ')
                    if data:
                        md5s[ftp_path + '/' + os.path.basename(data[1])] = data[0]
        tranname = os.path.basename(ftp_path.replace('/GCA','/GCF'))
        assembly = ftp_path + '/' + basename + '_genomic.fna.gz'
        if assembly in md5s:
            ass_md5 = md5s[assembly]
        proteome = ftp_path + '/' + basename + '_protein.faa.gz'
        if proteome in md5s:
            prot_md5 = md5s[proteome]
        gff3 = ftp_path + '/' + basename + '_genomic.gff3.gz'
        if gff3 in md5s:
            gff_md5 = md5s[gff3]
        else:
            test_gff3 = re.sub( r'\.gff3\.gz$', '.gff.gz', gff3 )
            if test_gff3 in md5s:
               gff3 = test_gff3
               gff_md5 = md5s[gff3]

        transcript = ftp_path.replace('/GCA', '/GCF') + '/' + tranname + '_rna.fna.gz'
        if transcript in md5s:
            trans_md5 = md5s[transcript]

        log_editor( 
            output_path + 'ncbiDwnld.log', str(index), 
            str(index) + '\t' + accession + '\t' +  assembly + '\t' + \
            proteome + '\t' + gff3 + '\t' + transcript + '\t' + \
            ass_md5 + '\t' + prot_md5 + '\t' + gff_md5 + '\t' + trans_md5
            )
        ass_prots[str(index)] = {
            'assembly': assembly, 'assembly_md5': ass_md5,
            'proteome': proteome, 'proteome_md5': prot_md5,
            'gff3': gff3, 'gff3_md5': gff_md5,
            'transcript': transcript, 'transcript_md5': trans_md5,
            'biosample': accession
            }

        time.sleep(0.1)
# if no API key is used, we can only generate 3 queries per second, otherwise we can use 10
        count += 1
        if not api_key:
            if count >= 2:
                time.sleep( 1 )
                count = 0
        else:
            if count >= 9:
                time.sleep( 1 )
                count = 0

    return ass_prots

# download the file depending on the type inputted
def download_files( ome_prots, ome, file_types, output_dir, remove = False ):

        for file_type in file_types:
            ftp_link = ome_prots[file_type]
            if file_type == 'assembly':
                file_path = output_dir + 'assembly/' + \
                    os.path.basename(ome_prots[file_type])
            elif file_type == 'proteome':
                file_path = output_dir + 'proteome/' + \
                    os.path.basename(ome_prots[file_type])
            elif file_type == 'gff3':
                file_path = output_dir + 'gff3/' + \
                    os.path.basename(ome_prots[file_type])
            elif file_type == 'transcript':
                file_path = output_dir + 'transcript/' + \
                    os.path.basename(ome_prots[file_type])

            if os.path.isfile( file_path ):
                md5_cmd = subprocess.run( [
                        'md5sum', file_path],
                        stdout = subprocess.PIPE )
                md5_res = md5_cmd.stdout.decode( 'utf-8' )
                md5_find = re.search(r'\w+', md5_res)
                md5 = md5_find[0]
                print(md5 + ' ' + ome_prots[file_type + '_md5'])
                if md5 == ome_prots[file_type + '_md5']:
                    print('yes')
                    continue

            if ftp_link == '':
                dwnld = 1
                continue
            request = requests.get( ftp_link.replace('ftp://','https://' ))
            if request.status_code != 200:
                eprint('\t' + ome + '\tno ' + file_type)
                dwnld = 69
                ome_prots[file_type] = ''
                log_editor( output_dir + 'ncbiDwnld.log', str(ome), str(ome) + \
                    '\t' + ome_prots['biosample'] + '\t' + ome_prots['assembly'] + \
                    '\t' + ome_prots['proteome'] + '\t' + ome_prots['gff3'] + \
                    '\t' + ome_prots['transcript'] + '\t' + \
                    ome_prots['assembly_md5'] + '\t' + \
                    ome_prots['proteome_md5'] + '\t' + \
                    ome_prots['gff3_md5'] + '\t' +
                    ome_prots['transcript_md5'] )
                if remove:
                    break
                continue

            dwnld = subprocess.call(
                ['curl', ftp_link, '-o', file_path], 
                stdout = subprocess.PIPE, stderr = subprocess.PIPE
                )

            if dwnld != 0:
                if remove:
                    break
            else:
                with open( file_path, 'r' ) as data:
                    try:
                        if len(data.read()) < 500:
#                           print('\t' + ome + '\t' + file_type + ' empty')
                            dwnld = 420
                            if remove:
                                break
                    except UnicodeDecodeError:
                        continue
        return dwnld

# generates ome codes from biosample accession (default) and first 3 letters of genus and species
def gen_omes(df,column='#Organism/Name',tag='BioSample Accession'):

    newdf = df
    newdf['internal_ome'] = np.nan
    tax_set = set()
    for index in range(len(df)):
        taxonomy = df[column][index]
        taxonomy = taxonomy.split(' ')
        genus = taxonomy[0]
        species = taxonomy[1]
        name = genus[:3] + species[:3]
        access = df[tag][index]
        number = 1

        while name + str(number) in tax_set:
            number += 1

        tax_set.add(name + str(number))

        if access != '-':
            ome = name + str(number) + '.' + access
        else:
            ome = name + str(number)

        newdf['internal_ome'][index] = ome
    
    return newdf


def callDwnldFTPs( ncbi_df, ass_prots, api, output_path ):

    if 'internal_ome' in ncbi_df.keys() and 'biosample' in ncbi_df.keys():
        ncbi_df = ncbi_df.set_index('internal_ome')
        if not all( x in ass_prots for x in ncbi_df.index ):
            ass_prots = collect_ftps( 
                ncbi_df, ass_prots, 
                column = 'biosample', api_key=api,
                output_path = output_path
                )
    elif 'biosample' in ncbi_df.keys():
        ncbi_df['index'] = ncbi_df['biosample']
        ncbi_df = ncbi_df.set_index( 'index' )
        ncbi_df.index = ncbi_df.index.astype(str)
        if not all( x in ass_prots for x in ncbi_df.index ):
            ass_prots = collect_ftps( 
            ncbi_df, ass_prots, 
            column = 'biosample', api_key=api, 
            output_path = output_path 
            )
    elif 'internal_ome' in ncbi_df.keys():
        ncbi_df = ncbi_df.set_index('internal_ome')
        if not all( x in ass_prots for x in ncbi_df.index ):
            ass_prots = collect_ftps( 
                ncbi_df, ass_prots, 
                api_key=api, output_path = output_path
                )
    else:
        if len( ncbi_df.columns ) > 1:
            ncbi_df = ncbi_df.set_index( indices[1] )
        ncbi_df.index = ncbi_df.index.astype(str)
        if not all( x in ass_prots for x in ncbi_df.index ):
            ass_prots = collect_ftps( 
                ncbi_df, ass_prots, 
                column = ncbi_df.columns[0], api_key=api,
                output_path = output_path
            )

    return ncbi_df, ass_prots


def main( 
    email = '', api = None, 
    assembly = True, proteome = True, gff3 = True, transcript = False,
    ncbi_df = False, remove = False, output_path = os.getcwd()
    ):

    output_path = formatPath(output_path)
    file_types = prepareFolders( 
        output_path, gff3, proteome, assembly, transcript
    )

    Entrez.email = email
    if api:
        Entrez.api_key = api

    ass_prots = compileLog(output_path, remove )

    if type(ncbi_df) != pd.DataFrame:
        ncbi_df = pd.read_csv( ncbi_df, sep = '\t' )
    if 'BioSample Accession' in ncbi_df:
        ncbi_df['biosample'] = ncbi_df['BioSample Accession']
        del ncbi_df['BioSample Accession']

    ncbi_df, ass_prots = callDwnldFTPs( ncbi_df, ass_prots, api, output_path )
    if remove:
        ass_prots = { 
    o: ass_prots[o] for o in ass_prots if all(ass_prots[o][p] != '' for p in ass_prots[o])
    }
    new_df = pd.DataFrame()
    
    print('\nDownloading NCBI files')
    for ome in ass_prots:
        print('\t' + str(ome))
        if ome not in set( ncbi_df.index ):
            continue
        if type(ome) is int:
            new_ome = ass_prots[ome]['biosample']
            exit = download_files( ass_prots[ome], new_ome, file_types, output_path, remove = remove )
        else:
            exit = download_files( ass_prots[ome], ome, file_types, output_path, remove = remove )
        if exit != 0 and remove:
            continue
        if ome in set(ncbi_df.index):
            ncbi_df.at[ome, 'assembly_path'] = os.path.basename( ass_prots[ome]['assembly'] )
            ncbi_df.at[ome, 'proteome_path'] = os.path.basename( ass_prots[ome]['proteome'] )
            ncbi_df.at[ome, 'gff3_path'] = os.path.basename( ass_prots[ome]['gff3'] )
        new_df = new_df.append( ncbi_df.loc[ome] )
 
    new_df = new_df.reset_index()
# won't work nice if ome isn't in the original dataframe
    new_df = new_df.rename(columns = { 'index': 'internal_ome' } )
    del new_df['gff3']
    del new_df['proteome']
    del new_df['assembly']

    return new_df


if __name__ == "__main__":

    parser = argparse.ArgumentParser( description = \
    'Downloads assemblies, proteomes, and/or gff3s in current directory.' 
    )
    parser.add_argument( '-i', '--input', required = True, \
    help = 'Imports tab delimited file with BioSample accessions under ' + \
         'column "biosample" or an NCBI master tsv with curated codenames' )
    parser.add_argument( '-a', '--assembly', action = 'store_true', \
    help = 'Download assemblies')
    parser.add_argument( '-p', '--proteome', action = 'store_true', \
    help = 'Download proteomes')
    parser.add_argument( '-g', '--gff3', action = 'store_true', \
    help = 'Download gff3s')
    parser.add_argument( '-t', '--transcript', action = 'store_true', \
    help = 'Download transcripts' )
    parser.add_argument( '-e', '--email', help = 'NCBI login email' )
    parser.add_argument( '-k', '--key', \
    help = 'API key to expedite interfacing with NCBI' )
    parser.add_argument( '-o', '--output', help = 'Output directory' )
    args = parser.parse_args()

    if not args.email:
        args.email = input( 'NCBI login email: ' )
        if not args.key:
            args.key = input( 'NCBI api key (leave blank if none: ' )

    if not args.output:
        output = os.getcwd()
    else:
        output = args.output

    args_dict = {
        'NCBI Table': args.input,
        'email': args.email,
        'Assemblies': args.assembly,
        'Proteomes': args.proteome,
        ".gff3's": args.gff3,
        'Transcripts': args.transcript
    }

    start_time = intro('Download NCBI files',args_dict)
    if not args.assembly and not args.proteome and not args.gff3:
        eprint('\nERROR: You must choose at least one download option\nExit code 37')
        sys.exit( 37 )

    new_df = main( 
        email = args.email, api = args.key, assembly = args.assembly,
        proteome = args.proteome, gff3 = args.gff3, transcript = args.transcript,
        ncbi_df = args.input, output_path = output
        )
    new_df.to_csv( args.input + '_dwnld', sep = '\t' )

    outro(start_time)
