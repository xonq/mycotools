#! /usr/bin/env python3

import subprocess, argparse, sys, os, re, time, getpass
import pandas as pd
from mycotools.lib.kontools import eprint, formatPath, outro, intro

def jgi_login( user, pwd ):
    '''Login via JGI's prescribed method by creating a cookie cache and 
    downloading JGI's sign-in file.'''

    null = os.path.expanduser( '~/.nullJGIdwnld' )
        
    login_cmd = subprocess.call(
        ['curl', 
        'https://signon.jgi.doe.gov/signon/create', '--data-urlencode', 'login=' + str(user), 
        '--data-urlencode', 'password=' + str(pwd), '-c', 'cookies', '-o', null],
        stdout = subprocess.PIPE, stderr = subprocess.PIPE
        )

    return login_cmd


def retrieveXML( ome, output ):
    '''Retrieve JGI xml file tree. First check if it already exists, if not then 
    download it using JGI's prescribed method. Then open the xml and check for 
    the common 'Portal does not exist' error. If so, report.'''

    if os.path.exists( output + '/' + str(ome) + '.xml' ):
        with open( output + '/' + str(ome) + '.xml', 'r' ) as xml_raw:
            xml_data = xml_raw.read()
        if xml_data == 'Portal does not exist':
            print( '\tERROR: `' + ome + ' not in JGIs `organism` database.' )
            xml_cmd = 1
            os.remove( output + '/' + ome + '.xml' )
        else:
            xml_cmd = -1

    else:
        xml_cmd = subprocess.call( [
            'curl', 
            "https://genome.jgi.doe.gov/portal/ext-api/downloads/get-directory?organism=" +
            str(ome), '-b', 'cookies', '-o', output + '/' + str(ome) + ".xml" ],
            stdout = subprocess.PIPE, stderr = subprocess.PIPE )
        if xml_cmd != 0:
            print('\tERROR: `xml` directory ' + ome + ' - `curl` error: ' + str(xml_cmd))

    if xml_cmd == 0:
        with open( output + '/' + ome + '.xml', 'r' ) as xml_raw:
            xml_data = xml_raw.read()
        if xml_data == 'Portal does not exist':
            print( '\tERROR: `' + ome + ' not in JGIs `organism` database.' )
            xml_cmd = 1
            os.remove( output + '/' + ome + '.xml' )

    return xml_cmd


def checkMatches( match1, match2, xml_data, file_type):

    findall, preexisting = False, None
    matches = match1.search(xml_data)
    if matches is None or 'get_tape_file' in matches[0]:
        matches = match2.search(xml_data)
        if matches is None or 'get_tape_file' in matches[0]:
            matches = match1.findall(xml_data)
            matches = [ x for x in matches if 'get_tape_file' not in x ]
            findall = True
            if matches:
                matches = match2.findall(xml_data)
            if len(matches) > 1:
                matches = matches[1]
            else:
                print('\t\tERROR: no ' + file_type + ' detected.')
                matches = None
                preexisting = True

    return matches, preexisting, findall


def JGIdwnld( ome, file_type, output, masked = True ):
    '''Download JGI files. For each type of file, use regular expressions to 
    gather the URL from the file. Grab the md5checksum if possible from the 
    xml as well. Create arbitrary values for md5 and curl_cmd. If the download file
    already exists, then run an md5 checksum if an md5 value exists in the xml.

    If the md5 does not match then open and check for typical errors. If those 
    errors exist, obtain the alternative download URL from the xml. If the md5 
    matches, pass through the rest of the function.

    If there is no download md5, then check to see if the file is greater than 10 
    lines as a proxy to make sure that the file isn't empty/blatantly wrong. If it 
    passes this test, change the values to not enter the while loop at the end of 
    the function.

    The while loop following the file exists error allows for 3 attempts. If a file 
    is downloaded it will check its md5 using the xml reference. If it fails, it will 
    proceed via the error checking above, wait a minute, and reattempt. If there is no 
    download md5 from the xml it will proceed via the error checking above as well. If it 
    passes the line count check, it will exit the loop - otherwise, it will attempt to 
    gather a new URL and restart after another minute wait. This is an unfortunate
    circumnavigation of JGI's lack of a consistent maximum ping / time before 
    booting.'''


    prefix = 'https://genome.jgi.doe.gov'

    with open( output + '/xml/' + ome + '.xml', 'r') as xml:
        data = xml.read()


    gffMatch1 = re.compile(r'filename\="([^ ]+\_GeneCatalog_genes_\d+\.gff\.gz)".*?url\="([^ ]+)".*?md5\="([^ ]+)"')
    gffMatch2 = re.compile(r'filename\="([^ ]+\_GeneCatalog_genes_\d+\.gff\.gz)".*?url\="([^ ]+)"')
    if file_type == 'gff':
        match1 = gffMatch1
        match2 = gffMatch2
    elif file_type == 'assembly':
        if masked:
            match1 = re.compile(r'filename\="(' + ome + r'\_AssemblyScaffolds_Repeatmasked\.fasta\.gz)".*?url\="([^ ]+)".*?md5\="([^ ]+)"')
            match2 = re.compile(r'filename\="(' + ome + r'\_AssemblyScaffolds_Repeatmasked\.fasta\.gz)".*?url\="([^ ]+)"')
        else:
            match1 = re.compile(r'filename\="(' + ome + r'\_AssemblyScaffolds\.fasta\.gz)".*?url\="([^ ]+)".*?md5\="([^ ]+)"')
            match2 = re.compile(r'filename\="(' + ome + r'\_AssemblyScaffolds\.fasta\.gz)".*?url\="([^ ]+)"')

    elif file_type == 'proteome':
        match1 = re.compile(r'filename\="([^ ]+\_GeneCatalog_proteins_\d+\.aa\.fasta\.gz)".*?url\="([^ ]+)".*?md5\="([^ ]+)"')
        match2 = re.compile(r'filename\="([^ ]+\_GeneCatalog_proteins_\d+\.aa\.fasta\.gz)".*?url\="([^ ]+)"')

    elif file_type == 'gff3':
        match1 = re.compile(r'filename\="([^ ]+\_GeneCatalog\_\d+\.gff3\.gz)".*?url\="([^ ]+)".*?md5\="([^ ]+)"')
        match2 = re.compile(r'filename\="([^ ]+\_GeneCatalog\_\d+\.gff3\.gz)".*?url\="([^ ]+)"')

    elif file_type == 'transcript':
        match1 = re.compile(r'filename\="([^ ]+_GeneCatalog_transcripts_\d+\.nt\.fasta\.gz)".*?url\="([^ ]+)".*?md5\="([^ ]+)"')
        match2 = re.compile(r'filename\="([^ ]+\_GeneCatalog_transcripts_\d+\.nt\.fasta\.gz)".*?url\="([^ ]+)"')

    elif file_type == 'AllTranscript':
        match1 = re.compile(r'filename\="([^ ]+\_all_transcripts_\d+\.nt\.fasta\.gz)".*?url\="([^ ]+)".*?md5\="([^ ]+)"')
        match2 = re.compile(r'filename\="([^ ]+\_all_transcripts_\d+\.nt\.fasta\.gz)".*?url\="([^ ]+)"')

    elif file_type == 'est':
        match1 = re.compile(r'filename\="([^ ]+\_EST_\d+_cluster_consensi\.fasta\.gz)".*?url\="([^ ]+)".*?md5\="([^ ]+)"')
        match2 = re.compile(r'filename\="([^ ]+\_EST_\d+_cluster_consensi\.fasta\.gz)".*?url\="([^ ]+)"')

    check = 1
    matches, preexisting, findall = checkMatches( match1, match2, data, file_type )
    if not matches and file_type == 'gff3':
        file_type = 'gff'
        matches, preexisting, findall = checkMatches( gffMatch1, gffMatch2, data, file_type )
    if matches:

        md5 = 420
        attempt = 0
        curl_cmd = 420

        if findall:
            dwnld_url = prefix + matches[1].replace( '&amp;', '&' )
        else:
            dwnld_url = prefix + matches[2].replace( '&amp;', '&' )

#        print('\t\t' + dwnld_url)
        if match1.search(data):
            dwnld_md5 = match1.search(data)[3]
        else:
            dwnld_md5 = None
        dwnld = os.path.basename( dwnld_url )

        if os.path.exists( output + '/' + file_type + '/' + dwnld):
            if dwnld_md5:
                md5_cmd = subprocess.run( [
                        'md5sum', output + '/' + file_type + '/' + dwnld], 
                        stdout = subprocess.PIPE )
                md5_res = md5_cmd.stdout.decode( 'utf-8' )
                md5_find = re.search(r'\w+', md5_res)
                md5 = md5_find[0]

            if md5 == dwnld_md5:
                curl_cmd = 0
                check = dwnld
                preexisting = True
            else:
                while True:
                    try:
                        with open(output + '/' + file_type + '/' + dwnld, 'r') \
                            as dwnld_data_raw:
                            dwnld_data = dwnld_data_raw.read()
                        if re.search( '307 Temporary Redirect', dwnld_data ):
                            print('\t\t' + dwnld + ' link has moved. Trying a different link.')
                            matches1 = re.search( 
                                    r'url\="(\/portal\/' + ome + r'\/download\/' + \
                                    dwnld + r')', data)
                            if matches1:
                                dwnld_url = prefix + matches1[1].replace( '&amp;', '&' )
                                print('\t\t' + dwnld_url)
                            break
                        else:
                            if not dwnld_md5:
                                check_size = subprocess.run( 
                                        ['wc', '-l', output + '/' + file_type + '/' + dwnld], 
                                        stdout = subprocess.PIPE 
                                        )
                                check_size_res = check_size.stdout.decode( 'utf-8' )
                                print('\t\tFile exists - no md5 to check.')
                                check_size_find = re.search(r'\d+', check_size_res)
                                size = check_size_find[0]
                                if int(size) < 10:
                                    print('\t\tInvalid file size.')
                                else:
                                    md5 = None
                            else:
                                print('\t\tFile exists, but md5 does not match.')
                            break
                    except UnicodeDecodeError:
                        if not dwnld_md5:
                            print('\t\tFile exists - no md5 to check.')
                            check_size = subprocess.run( 
                                    ['wc', '-l', output + '/' + file_type + '/' + dwnld], 
                                    stdout = subprocess.PIPE )
                            check_size_res = check_size.stdout.decode( 'utf-8' )
                            check_size_find = re.search(r'\d+', check_size_res)
                            size = check_size_find[0]
                            if int(size) < 10:
                                print('\t\tInvalid file size.')
                            else:
                                md5 = None
                                preexisting = True
                                check = dwnld
                        else:
                            print('\t\tFile exists, but md5 does not match.')
                        break

        while md5 != dwnld_md5 and curl_cmd != 0 and attempt < 3:
            attempt += 1
            curl_cmd = subprocess.call( 
                    ['curl', dwnld_url, '-b', 'cookies', '-o', 
                    output + '/' + file_type + '/' + dwnld], 
                    stdout = subprocess.PIPE, stderr = subprocess.PIPE 
                    )

            if curl_cmd == 0:
                check = dwnld
                
                if dwnld_md5:
                    md5_cmd = subprocess.run( 
                            ["md5sum", output + '/' + file_type + '/' + dwnld], 
                            stdout = subprocess.PIPE 
                            )
                    md5_res = md5_cmd.stdout.decode( 'utf-8' )
                    md5_find = re.search(r'\w+', md5_res)
                    md5 = md5_find[0]

                if not dwnld_md5:
                        try:
                            with open(output + '/' + file_type + '/' + dwnld, 'r') \
                                as dwnld_data_raw:
                                dwnld_data = dwnld_data_raw.read()
                            if re.search( '307 Temporary Redirect', dwnld_data ):
                                print(
                                    '\t\t\t' + dwnld + ' link has moved. ' + \
                                    'Trying a different link.'
                                    )
                                matches1 = re.search( 
                                        r'url\="(\/portal\/' + ome + r'\/download\/' + \
                                        dwnld + r')', data
                                        )
                                if matches1:
                                    dwnld_url = prefix + matches1[1].replace( '&amp;', '&' )
                                    print('\t\t' + dwnld_url)
                                else:
                                    print('\t\t\t\tNo valid alternative link. Skip.')
                                    attempt = 4
                                    break
                            else:
                                print('\t\tFile exists - no md5 to check.')
                                check_size = subprocess.run( 
                                        "wc -l " + output + '/' + file_type + '/' + dwnld, 
                                        stdout = subprocess.PIPE 
                                        )
                                check_size_res = check_size.stdout.decode( 'utf-8' )
                                check_size_find = re.search(r'\d+', check_size_res)
                                size = check_size_find[0]
                                if int(size) < 10:
                                    print('\t\tInvalid file size. Attempt ' + str(attempt))
                                else:
                                    md5 = None
                                    break
                        except UnicodeDecodeError:
                            pass
                        time.sleep( 60 )

                elif md5 != dwnld_md5 and attempt == 1:
                    print('\t\tERROR: md5 does not match JGI. Attempt ' + str(attempt))
                    curl_cmd = -1
                    check = 2
                    while True:
                        try:
                            with open(output + '/' + file_type + '/' + dwnld, 'r') \
                                as dwnld_data_raw:
                                dwnld_data = dwnld_data_raw.read()
                            if re.search( '307 Temporary Redirect', dwnld_data ):
                                print(
                                    '\t\t\t' + dwnld + ' link has moved. ' + \
                                    'Trying a different link.'
                                )
                                matches1 = re.search( r'url\="(\/portal\/' + \
                                        ome + r'\/download\/' + dwnld + r')', data)
                                if matches1:
                                    dwnld_url = prefix + matches1[1].replace( '&amp;', '&' )
                                    print('\t\t' + dwnld_url)
                                else:
                                    print('\t\t\t\tNo valid alternative link. Skip.')
                                    attempt = 4
                            break
                        except UnicodeDecodeError:
                            pass
                            break
                        time.sleep( 60 )
                elif md5 != dwnld_md5 and attempt == 2:
                    print('\t\tERROR: md5 does not match JGI. Attempt ' + str(attempt))
                    curl_cmd = -1
                    matches1 = re.search( 
                        r'url\="(\/portal\/' + ome + r'\/download\/' + dwnld + r')', data
                        )
                    if matches1:
                        dwnld_url = prefix + matches1[1].replace( '&amp;', '&' )
                        print('\t\t' + dwnld_url)
                    time.sleep( 60 )
                    check = 2
                elif md5 != dwnld_md5:
                    print('\t\tERROR: md5 does not match JGI. Attempt ' + str(attempt))
                    check = 2
            else:
                print(
                    '\t\tERROR: Failed to retrieve ' + file_type + '. `curl` error: ' \
                    + str(curl_cmd) + '\nAttempt' + str(attempt)
                    )
                check = 2

        if attempt == 3:
            if md5 != dwnld_md5:
                print('\t\tExcluding from database - potential failure')
                curl_cmd = 0
            if curl_cmd != 0:
                print('\t\tFile failed to download')

#    else:
 #       print('\t\tfailed')

    return check, preexisting, file_type



def main( 
    df, output, user, pwd, assembly = True, proteome = True, gff3 = True, 
    gff = False, transcript = False, est = False, masked = True
    ):

    if not 'genome_code' in df.columns:
        if len( df.columns ) != 1:
            print( '\nInvalid input. No genome_code column and more than one column.' )
        else:
            ome_col = list(df.columns)[0]
    else:
        ome_col = 'genome_code'

    print('\nLogging into JGI ...')
    login_attempt = 0
    while jgi_login( user, pwd ) != 0 and login_attempt < 5:
        eprint('\nJGI Login Failed. Attempt: ' + str(login_attempt))
        jgi_login( user, pwd )
        time.sleep( 5 )
        login_attempt += 1
        if login_attempt == 5:
            eprint('\nERROR: Failed login 5 attempts.')
            sys.exit( 100 )

    if not os.path.exists( output + '/xml' ):
        os.mkdir( output + '/xml' )
# perhaps add a counter here, but one that checks if it is actually querying jgi
    print( '\nRetrieving `xml` directories ...' )
    ome_set, count = set(), 0
    for i,row in df.iterrows():
        error_check = retrieveXML( row[ome_col], output + '/xml' )
        if error_check > 0:
            ome_set.add( row[ome_col] )
        if error_check != -1:
            time.sleep( 0.1 )

    print( '\nDownloading' )
    dwnlds = []
    if assembly:
        dwnlds.append( 'assembly' )
    if proteome:
        dwnlds.append( 'proteome' )
    if gff3:
        dwnlds.append( 'gff3' )
    if gff:
        dwnlds.append( 'gff' )
    if transcript:
        dwnlds.append( 'transcript' )
#        dwnlds.append( 'AllTranscript' )
    if est:
        dwnlds.append( 'est' )

    for typ in dwnlds:
        if not os.path.isdir( output + '/' + typ ):
            os.mkdir( output + '/' + typ )
        if typ == 'gff3':
            if not os.path.isdir( output + '/gff' ):
                os.mkdir( output + '/gff')

    preexisting = True
    for i, row in df.iterrows():
        ome = row[ome_col]
        if not preexisting:
            time.sleep( 60 )
        if ome not in ome_set:
            jgi_login( user, pwd )
            print( '\t' + ome )
            for typ in dwnlds:
                check, preexisting, new_typ = JGIdwnld(
                    ome, typ, output, masked = masked
                )
                if type(check) != int:
                    if new_typ == 'gff':
                        df.at[i, 'jgi_gff_path'] = output + '/' + new_typ \
                            + '/' + check
                    else:
                        df.at[i, new_typ + '_path'] = output + '/' + new_typ + \
                            '/' + check
                    check = os.path.basename( os.path.abspath( check ) )
                print('\t\t' + new_typ + ': exit status ' + str(check))
        else:
            eprint( '\t' + ome + ' failed.' )

    if os.path.exists( 'cookies' ):
        os.remove( 'cookies' )
    if os.path.exists( os.path.expanduser( '~/.null' )):
        os.remove( os.path.expanduser( '~/.null' ))

    if 'gff3' in df.columns:
        del df['gff3']
    if 'proteome' in df.columns:
        del df['proteome']
    if 'assembly' in df.columns:
        del df['assembly']

    return df


if __name__ == '__main__':

    parser = argparse.ArgumentParser( description = \
        "Imports table/database with JGI `genome_code` column and downloads assembly, proteome, gff, " + \
        'and/or gff3. This script supports rerunning/continuing previous runs in the same directory. ' + \
        'JGI has stringent, nondefined ping limits, so file downloads are limited to 1 per minute.' )
    parser.add_argument( '-i', '--input_table', required = True, help = 'Database/table with `genome_code` column of JGI ome codes' )
    parser.add_argument( '-a', '--assembly', default = False, action = 'store_true', \
        help = 'Download assemblies.' )
    parser.add_argument( '-p', '--proteome', default = False, action = 'store_true', \
        help = 'Download proteomes.' )
    parser.add_argument( '-g', '--gff', default = False, action = 'store_true', \
        help = 'Download gffs.' )
    parser.add_argument( '-g3', '--gff3', default = False, action = 'store_true', \
        help = 'Download gff3s.' )
    parser.add_argument( '-t', '--transcript', default = False, action = 'store_true', \
        help = 'Download transcripts.' )
    parser.add_argument( '-e', '--est', default = False, action = 'store_true', \
        help = 'Download EST data.' )
    parser.add_argument( '--nonmasked', default = False, action = 'store_true', \
        help = "Download nonmasked assemblies. Requires -a" )
    parser.add_argument( '-o', '--output', default = os.getcwd(), help = 'Output dir' )
    args = parser.parse_args()

    if args.repeatmasked:
        args.assembly = True

    if not args.assembly and not args.proteome \
        and not args.gff and not args.transcript \
        and not args.est and not args.gff3:
        eprint( '\nERROR: You must choose at least one download option.' )

    user = input( 'JGI username: ' )
    pwd = getpass.getpass( prompt='JGI Login Password: ' )

    args_dict = {
            'JGI Table': args.input_table,
            'Assemblies': args.assembly,
            'RepeatMasked': not args.nonmasked,
            'Proteomes': args.proteome,
            ".gff's": args.gff,
            ".gff3's": args.gff3,
            'Transcripts': args.transcript,
            'EST': args.est
            }

    start_time = intro( 'Download JGI files', args_dict )
    df = pd.read_csv( args.input_table, sep = '\t' )
    output = formatPath(args.output)

    jgi_df = main( 
        df, output, user, pwd, args.assembly, args.proteome, 
        args.gff3, args.gff2, args.transcript, args.est, not args.nonmasked
        )
    df.to_csv( os.path.normpath(args.input_table) + '_jgiDwnld', sep = '\t', index = False )

    outro( start_time )
