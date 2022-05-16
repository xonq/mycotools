#! /usr/bin/env python3

from Bio import Entrez
import pandas as pd, numpy as np
import argparse, os, time, subprocess, sys, requests, re, gzip
from mycotools.lib.kontools import intro, outro, formatPath, prep_output, eprint, vprint, findExecs
from mycotools.lib.dbtools import log_editor, loginCheck, db2df
from datetime import datetime


def prepareFolders( output_path, gff, prot, assem, transcript ):

    file_types = []
    if assem:
        if not os.path.exists(output_path + 'assembly'):
            os.mkdir(output_path + 'assembly')
        file_types.append( 'assembly' )
    if gff:
        if not os.path.exists(output_path + 'gff3'):
            os.mkdir(output_path + 'gff3')
        file_types.append( 'gff3' )
    if prot:
        if not os.path.exists(output_path + 'proteome'):
            os.mkdir(output_path + 'proteome')
        file_types.append( 'proteome' ) 
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
            'assemMD5\tprotMD5\tgff3MD5\ttransMD5\tgenome_id(s)')

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
                        'transcript_md5': str(data[9]), 'genome_id': data[10].split(',')
                        }

    return ass_prots


# collects paths to download proteomes and assemblies
def collect_ftps(
    ncbi_df, ass_prots, api_key=0,
    column='Assembly Accession', database="assembly", output_path = '',
    verbose=True, remove = False
    ):

    eprint('\nAssembling NCBI ftp directories', flush = True)
    count, failed = 0, []

# for each row in the assembly, grab the accession number, form the search term for Entrez, use Entrez,
    for index, row in ncbi_df.iterrows():
        accession = row[column]
        if 'biosample' in set(row.keys()):
            failure_acc = row['biosample']
        else:
            failure_acc = row[column]
        if pd.isnull(row[column]):
            ass_prots[str(index)] = {
                'biosample': failure_acc, 'assembly': '',
                'proteome': '', 'gff3': '', 'transcript': '',
                'assembly_md5': '', 'proteome_md5': '',
                'gff3_md5': '', 'transcript_md5': '', 'genome_id': ''
                }
            failed.append([failure_acc, datetime.strftime(row['version'], '%Y%m%d')])
            continue

#        ome = row['internal_ome']
        if index in ass_prots:
            continue

        search_term, esc_count = accession + '[' + str(column) + ']', 0
        while esc_count < 10:
            try:
                handle = Entrez.esearch(db=database, term=search_term)
                genome_id = Entrez.read(handle)['IdList']
                break
            except RuntimeError:
                time.sleep(1)
                esc_count += 1
        else:
            print('\tERROR:', index, 'failed to search NCBI')
        try:
            ID = genome_id[0]
        except IndexError:
            if 'internal_ome' in row.keys():
                accession = row['internal_ome']
            eprint('\t' + accession + ' failed to find genome ID', flush = True)
            try:
                failed.append( [failure_acc, datetime.strftime(row['version'], '%Y%m%d')] )
            except TypeError:
                failed.append( [failure_acc, row['version']] )
            continue
        ass_prots[str(index)] = {
            'biosample': '', 'assembly': '', 'proteome': '',
            'gff3': '', 'transcript': '', 'assembly_md5': '',
            'proteome_md5': '', 'gff3_md5': '', 
            'transcript_md5': '', 'genome_id': genome_id
            }


# if more than one genome id is found, this needs to be reported because we're only using the first (shouldn't happen)
        if len(genome_id) > 1:
            vprint('\t' + accession + ' yields more than one genome ID:', v = verbose, flush = True)
            for ID in genome_id:
                vprint('\t\t' + str(ID), v = verbose, flush = True)
            ID = genome_id[0]
            vprint('\t\tUsing ' + str(ID), v = verbose, flush = True)

# obtain the path from a summary of the ftp directory and create the standard paths for proteomes and assemblies
        esc_count = 0
        while esc_count < 10:
            try:
                esc_count += 1
                handle = Entrez.esummary(db=database, id=ID, report="full")
                record = Entrez.read(handle, validate = False)
                ftp_path = str(record['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_GenBank'])
                esc_count = 0
                break
            except IndexError:
                time.sleep(0.1)
        else:
            if esc_count == 20:
                eprint('\tERROR: Failed to fulfill ftp request!', flush = True)
                sys.exit(69)


        if not ftp_path:
            if 'internal_ome' in row.keys():
                accession = row['internal_ome']
            eprint('\t' + accession + ' failed to return any FTP path', flush = True)
            failed.append( [failure_acc, datetime.strftime(row['version'], '%Y%m%d')] )
            continue

        ass_md5, gff_md5, trans_md5, prot_md5, md5s = '', '', '', '', {}
        basename = os.path.basename(ftp_path)
        while True:
            try:
                esc_count += 1
                request = requests.get( 
                    (ftp_path + '/md5checksums.txt').replace('ftp://','https://'), timeout = 120
                    )
                break
            except requests.exceptions.ConnectionError:
                if esc_count == 20:
                    eprint('\n\tERROR: Failed to fulfill ftp request!', flush = True)
                    sys.exit(70)
                count = 0
                time.sleep(1)
        count += 1
        if request.status_code != 200:
            eprint('\t' + ome + '\tno ' + md5, flush = True)
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
        else:
            assembly = ''
        proteome = ftp_path + '/' + basename + '_protein.faa.gz'
        if proteome in md5s:
            prot_md5 = md5s[proteome]
        else:
            proteome = ''
        gff3 = ftp_path + '/' + basename + '_genomic.gff3.gz'
        test_gff3 = re.sub( r'\.gff3\.gz$', '.gff.gz', gff3 )
        if gff3 in md5s:
            gff_md5 = md5s[gff3]
        elif test_gff3 in md5s:
            gff3 = test_gff3
            gff_md5 = md5s[gff3]
        else:
            gff3 = ''

        transcript = ftp_path.replace('/GCA', '/GCF') + '/' + tranname + '_rna.fna.gz'
        if transcript in md5s:
            trans_md5 = md5s[transcript]

        if (not assembly or not gff3) and remove:
            failed.append([row['biosample'], datetime.strftime(row['version'], '%Y%m%d')])

        log_editor( 
            output_path + 'ncbiDwnld.log', str(index), 
            str(index) + '\t' + accession + '\t' +  assembly + '\t' + \
            proteome + '\t' + gff3 + '\t' + transcript + '\t' + \
            ass_md5 + '\t' + prot_md5 + '\t' + gff_md5 + '\t' + trans_md5 + \
            '\t' + ','.join(genome_id)
            )
        ass_prots[str(index)] = {
            'assembly': assembly, 'assembly_md5': ass_md5,
            'proteome': proteome, 'proteome_md5': prot_md5,
            'gff3': gff3, 'gff3_md5': gff_md5,
            'transcript': transcript, 'transcript_md5': trans_md5,
            'biosample': accession, 'genome_id': genome_id
            }

# if no API key is used, we can only generate 3 queries per second, otherwise we can use 10
        count += 1
        if not api_key:
            if count >= 2:
                time.sleep( 1 )
                count = 0
        else:
            if count >= 7:
                time.sleep( 1 )
                count = 0

    return ass_prots, failed

# download the file depending on the type inputted
def download_files( ome_prots, ome, file_types, output_dir, count, remove = False ):

    esc_count, dwnlds = 0, {}
    for file_type in file_types:
        ftp_link = ome_prots[file_type]
        dwnlds[file_type] = -1
        if file_type == 'assembly':
            file_path = output_dir + 'assembly/' + \
                os.path.basename(ome_prots[file_type])
        elif file_type == 'gff3':
            file_path = output_dir + 'gff3/' + \
                os.path.basename(ome_prots[file_type])
        elif file_type == 'proteome':
            file_path = output_dir + 'proteome/' + \
                os.path.basename(ome_prots[file_type])
        elif file_type == 'transcript':
            file_path = output_dir + 'transcript/' + \
                os.path.basename(ome_prots[file_type])

        if os.path.isfile( file_path ):
            count += 1
            md5_cmd = subprocess.run( [
                'md5sum', file_path],
                stdout = subprocess.PIPE )
            md5_res = md5_cmd.stdout.decode( 'utf-8' )
            md5_find = re.search(r'\w+', md5_res)
            md5 = md5_find[0]
            if md5 == ome_prots[file_type + '_md5']:
                eprint('\t\t' + file_type + ': ' + os.path.basename(file_path), flush = True)
                dwnlds[file_type] = 0
                continue
        elif os.path.isfile( file_path[:-3] ):
            eprint('\t\t' + file_type + ': ' + os.path.basename(file_path), flush = True)
            dwnlds[file_type] = 0
            continue

        if ftp_link == '':
            dwnlds[file_type] = 15
            continue
        while True:
            try:
                esc_count += 1
                request = requests.get( ftp_link.replace('ftp://','https://' ))
                break
            except requests.exceptions.ConnectionError:
                if esc_count == 20:
                    eprint('\n\tERROR: Failed to fulfill ftp request!', flush = True)
                    sys.exit(71)
                time.sleep(1)
                count = 0
        count += 1
        if request.status_code != 200:
            eprint('\t\t' + file_type + ': ERROR, no file exists', flush = True)
            dwnlds[file_type] = 69
            ome_prots[file_type] = ''
            log_editor( output_dir + 'ncbiDwnld.log', str(ome), str(ome) + \
                '\t' + str(ome_prots['biosample']) + '\t' + str(ome_prots['assembly']) + \
                '\t' + str(ome_prots['proteome']) + '\t' + str(ome_prots['gff3']) + \
                '\t' + str(ome_prots['transcript']) + '\t' + \
                str(ome_prots['assembly_md5']) + '\t' + \
                str(ome_prots['proteome_md5']) + '\t' + \
                str(ome_prots['gff3_md5']) + '\t' + \
                str(ome_prots['transcript_md5']) + '\t' + \
                str(','.join([str(x) for x in ome_prots['genome_id']])))
            if remove and file_type in {'assembly', 'gff3'}:
                break
            continue

        count += 1
        dwnlds[file_type] = subprocess.call(
            ['curl', ftp_link, '-o', file_path], 
            stdout = subprocess.PIPE, stderr = subprocess.PIPE
            )

        if dwnlds[file_type] != 0:
            eprint('\t\t' + file_type + ': ERROR, download failed', flush = True)
            if remove and file_type in {'assembly', 'gff3'}:
                break
        else:
            if os.stat(file_path).st_size < 150:
                eprint('\t\t' + file_type + ': ERROR, file too small', flush = True)
#                       print('\t' + ome + '\t' + file_type + ' empty', flush = True)
                dwnlds[file_type] = 420
                if remove and file_type in {'assembly', 'gff3'}:
                    break
        eprint('\t\t' + file_type + ': ' + os.path.basename(file_path), flush = True)

    return dwnlds, count


def callDwnldFTPs( ncbi_df, ass_prots, api, output_path, verbose = False, remove = False ):

    failed = []
    if 'biosample' in ncbi_df.keys():
        ncbi_df['index'] = ncbi_df['biosample']
        ncbi_df = ncbi_df.set_index('index')
        if not all( x in ass_prots for x in list(ncbi_df.index) ):
            ass_prots, failed = collect_ftps( 
                ncbi_df, ass_prots, 
                column = 'biosample', api_key=api,
                output_path = output_path, verbose = verbose,
                remove = remove
                )
    elif 'internal_ome' in ncbi_df.keys():
        ncbi_df = ncbi_df.set_index('internal_ome')
        if not all( x in ass_prots for x in list(ncbi_df.index) ):
            ass_prots, failed = collect_ftps( 
                ncbi_df, ass_prots, remove = remove,
                api_key=api, output_path = output_path, verbose = verbose
                )
    else:
  #  elif len(ncbi_df.index) > 1:
        ncbi_df['index'] = ncbi_df[ncbi_df.columns[0]]
        ncbi_df['biosample'] = ncbi_df[ncbi_df.columns[0]]
        ncbi_df = ncbi_df.set_index('index')
        ncbi_df.index = ncbi_df.index.astype(str)
        if not all( x in ass_prots for x in list(ncbi_df.index) ):
            ass_prots, failed = collect_ftps( 
                ncbi_df, ass_prots, remove = remove,
                column = ncbi_df.columns[0], api_key=api,
                output_path = output_path, verbose = verbose
                )
#    else:
 #       print(ncbi_df)

    return ncbi_df, ass_prots, failed


def main( 
    api = None, 
    assembly = True, proteome = False, gff3 = True, transcript = False,
    ncbi_df = False, remove = False, output_path = os.getcwd(), verbose = False
    ):

    output_path = formatPath(output_path)
    file_types = prepareFolders( 
        output_path, gff3, proteome, assembly, transcript
    )

    ass_prots = compileLog(output_path, remove )

    if not isinstance(ncbi_df, pd.DataFrame) and os.path.isfile(ncbi_df):
  #      try:
        ncbi_df = db2df(ncbi_df)
#        except:
 #           ncbi_df = pd.read_csv( ncbi_df, sep = '\t' )
    elif not isinstance(ncbi_df, pd.DataFrame):
        ncbi_df = pd.DataFrame({'biosample': [ncbi_df]})
    if len(ncbi_df.index) == 0:
        ncbi_df = pd.DataFrame(
            {i: [v] for i, v in enumerate(list(ncbi_df.keys()))}
            )
    if 'BioSample Accession' in ncbi_df.keys() and not 'biosample' in ncbi_df.keys():
        ncbi_df['biosample'] = ncbi_df['BioSample Accession']
    if 'Modify Date' in ncbi_df.keys() and not 'version' in ncbi_df.keys():
        ncbi_df['version'] = pd.to_datetime(ncbi_df['Modify Date'])
    elif 'version' not in ncbi_df.keys():
        ncbi_df['version'] = ''

    failed = []
    ncbi_df, ass_prots, failed = callDwnldFTPs( 
        ncbi_df, ass_prots, api, output_path, verbose, remove
        )
    if remove:
        ass_prots = { 
            o: ass_prots[o] for o in ass_prots \
            if all(ass_prots[o][p] for p in ['assembly', 'gff3'])
            }
    new_df = pd.DataFrame()

    eprint('\nDownloading NCBI files', flush = True)
    count = 0
    for ome in ass_prots:
        if count >= 2:
            if not api:
                time.sleep(1)
                count = 0
            elif count >= 7:
                time.sleep(1)
                count = 0
        if ome not in set( ncbi_df['biosample'] ):
            continue
        eprint('\t' + str(ome), flush = True)
        new_ome = ass_prots[ome]['biosample']
        if ass_prots[ome]:
            exits, count = download_files( 
                ass_prots[ome], new_ome, file_types, output_path, 
                count, remove = remove 
                )
        else:
            check = ncbi_df[ncbi_df['biosample'] == ome]
            db_vers = datetime.strftime(check.iloc[0]['version'], '%Y%m%d')
            failed.append([ome, db_vers])
            continue
        if 'assembly' in exits:
            if exits['assembly'] != 0:
#        if any(exits[x] != 0 for x in {'assembly', 'gff3'}):
                failed.append([ome, ncbi_df['version'][ome]])
                continue
        if 'gff3' in exits:
            if exits['gff3'] != 0:
                 failed.append([ome, ncbi_df['version'][ome]])
                 continue
        ncbi_df.loc[ome, 'assemblyPath'] = output_path + 'assembly/' + \
            os.path.basename(ass_prots[ome]['assembly'])
        ncbi_df.loc[ome, 'proteome'] = output_path + 'proteome/' + \
            os.path.basename(ass_prots[ome]['proteome'])
        ncbi_df.loc[ome, 'gffPath'] = output_path + 'gff3/' + \
            os.path.basename(ass_prots[ome]['gff3'])
        new_df = new_df.append( ncbi_df.loc[ome] )
 
    new_df = new_df.reset_index()
# won't work nice if ome isn't in the original dataframe
#    new_df = new_df.rename(columns = { 
 #       'index': 'internal_ome'
  #      } )

    return new_df, failed


def getSRA(biosample, fastqdump = 'fastq-dump', pe = True):

    handle = Entrez.esearch(db='SRA', term=biosample)
    ids = Entrez.read(handle)['IdList']
    for id in ids:
        handle = Entrez.esummary(db='SRA', id = id, report='full')
        records = Entrez.read(handle, validate = False)
        for record in records:
            srr = re.search(r'Run acc="(S\w+\d+)"', record['Runs'])[1]
            print('\t\t' + srr, flush = True)
            if pe:
                subprocess.call([
                    fastqdump, '--gzip', '--split-3', srr], 
                    stdout = subprocess.PIPE)
                if os.path.isfile(srr + '_1.fastq.gz'):
                    os.rename(srr + '_1.fastq.gz', biosample + '_' + srr + '_1.fq.gz')
                    os.rename(srr + '_2.fastq.gz', biosample + '_' + srr + '_2.fq.gz')
                else:
                    print('\t\t\tERROR: file failed', flush = True)
            else:
                subprocess.call([
                    fastqdump, '--gzip', srr],
                    stdout = subprocess.PIPE)
                if os.path.isfile(srr + '.fastq.gz'):
                    os.rename(srr + '.fastq.gz', biosample + '_' + srr + '.fq.gz')
                else:
                    print('\t\t\tERROR: file failed', flush = True)


def goSRA(df, output = os.getcwd() + '/', pe = True):

    print()
    sra_dir = output + 'sra/'
    if not os.path.isdir(sra_dir):
        os.mkdir(sra_dir)
    os.chdir(sra_dir)
    fastqdump = findExecs('fastq-dump', exit = set('fastq-dump'))
    count = 0

    if 'sra' in df.keys():
        row_key = 'sra'
    else:
        row_key = 'biosample'

    for i, row in df.iterrows():
        print('\t' + row[row_key], flush = True)
        getSRA(row[row_key], fastqdump[0])
        count +=1
        if count >= 10:
            time.sleep(1)
            count = 0


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument( '-i', '--input', required = True, \
    help = 'Standalone accession; tab delimited file with accessions under ' + \
         'column "biosample" or "sra"; or tsv directly downloaded from NCBI' )
    parser.add_argument( '-a', '--assembly', action = 'store_true')
    parser.add_argument( '-p', '--proteome', action = 'store_true')
    parser.add_argument( '-g', '--gff3', action = 'store_true')
    parser.add_argument( '-t', '--transcript', action = 'store_true')
    parser.add_argument( '-s', '--sra', action = 'store_true', \
    help = 'Download SRAs only' )
    parser.add_argument( '-pe', '--paired', action = 'store_true', \
    help = 'Download paired-end SRAs. (REQUIRES -s)')
    parser.add_argument( '-o', '--output', help = 'Output directory' )
    args = parser.parse_args()

    ncbi_email, ncbi_api, jgi_email, jgi_pwd = loginCheck(jgi = False) 

    Entrez.email = ncbi_email
    if ncbi_api:
        Entrez.api_key = ncbi_api

    if not args.output:
        output = os.getcwd() + '/'
    else:
        output = formatPath(args.output)

    args_dict = {
        'NCBI Table': args.input,
        'email': ncbi_email,
        'Assemblies': args.assembly,
        'Proteomes': args.proteome,
        ".gff3's": args.gff3,
        'Transcripts': args.transcript,
        'SRA': args.sra
    }

    start_time = intro('Download NCBI files',args_dict)
    if not args.assembly and not args.proteome and not args.gff3 and not args.sra and not args.transcript:
        eprint('\nERROR: You must choose at least one download option\nExit code 37', flush = True)
        sys.exit( 37 )

    if args.sra:
        if os.path.isfile(args.input):
            goSRA(pd.read_csv(args.input, sep = '\t'), output, pe = args.paired)
        else:
            goSRA(pd.DataFrame({'sra': [args.input.rstrip()]}), output, pe = args.paired)
    else:
        new_df, failed = main( 
            assembly = args.assembly,
            proteome = args.proteome, gff3 = args.gff3, transcript = args.transcript,
            ncbi_df = args.input, output_path = output, verbose = True
            )
        new_df.to_csv( args.input + '_dwnld', sep = '\t' )

    outro(start_time)
