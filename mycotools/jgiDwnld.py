#! /usr/bin/env python3
"""
PLEASE respect JGI's ping time limits. I've tuned it to respect their
unannounced limit.


NEED to remove gff v gff3 option
"""

import os
import re
import sys
import time
import getpass
import argparse
import subprocess
import pandas as pd
import xml.etree.ElementTree as ET
from tqdm import tqdm
from mycotools.lib.kontools import eprint, format_path, outro, intro
from mycotools.lib.dbtools import loginCheck

def jgi_login(user, pwd):
    '''Login via JGI's prescribed method by creating a cookie cache and 
    downloading JGI's sign-in file.'''

    null = os.path.expanduser('~/.nulljgi_dwnld')
        
    login_cmd = subprocess.call(
        ['curl', 
        'https://signon.jgi.doe.gov/signon/create', '--data-urlencode', 'login=' + str(user), 
        '--data-urlencode', 'password=' + str(pwd), '-c', 'cookies', '-o', null],
        stdout = subprocess.PIPE, stderr = subprocess.PIPE
        )

    return login_cmd

def dwnld_xml(output, ome, max_tempts = 2):
    attempts = 0
    while not os.path.isfile(f'{output}/{ome}.xml') and attempts < max_tempts:
        attempts += 1
        xml_cmd = subprocess.call([
            'curl', 
            "https://genome.jgi.doe.gov/portal/ext-api/downloads/get-directory?organism=" \
             + str(ome), '-b', 'cookies', '-o', f'{output}/{ome}.xml'],
            stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        if xml_cmd != 0:
            print(f'\tERROR: {ome} xml curl error: {xml_cmd}', flush = True)
    if not os.path.isfile(f'{output}/{ome}.xml'):
        return -1
    else:
        return xml_cmd


def retrieve_xml(ome, output):
    '''Retrieve JGI xml file tree. First check if it already exists, if not then 
    download it using JGI's prescribed method. Then open the xml and check for 
    the common 'Portal does not exist' error. If so, report.'''

    if os.path.exists(output + '/' + str(ome) + '.xml'):
        with open(output + '/' + str(ome) + '.xml', 'r') as xml_raw:
            xml_data = xml_raw.read()
        if xml_data == 'Portal does not exist':
            print('\tERROR: `' + ome + ' not in JGIs `organism` database', flush = True)
            xml_cmd = 1
            os.remove(output + '/' + ome + '.xml')
        elif not xml_data:
            xml_cmd = None
            os.remove(f'{output}/{ome}.xml')
        else:
            xml_cmd = -1
    else:
        xml_cmd = dwnld_xml(output, ome)

    if xml_cmd == 0:
        with open(output + '/' + ome + '.xml', 'r') as xml_raw:
            xml_data = xml_raw.read()
        if xml_data == 'Portal does not exist':
            print('\tERROR: `' + ome + ' not in JGIs `organism` database', flush = True)
            xml_cmd = 1
            os.remove(output + '/' + ome + '.xml')
        elif not xml_data:
            xml_cmd = None
            os.remove(f'{output}/{ome}.xml')

    return xml_cmd


def parse_xml(ft, xml_file, masked = False, forbidden = {}, filtered = True):
    """Parse the XML data to obtain the file types of interest based on
    predefined hashes that contain the known subdirectories associated with JGI
    organization"""

    if ft == 'fna':
        if masked:
            ft += '$masked'
        else:
            ft += '$unmasked'

    # set the initial hashes for the XML hierarchy - relate file types to their
    # hierarchy structure
    ft2xt = {'fna$masked': {'assembly'}, 'fna$unmasked': {'assembly'},
             'gff': {'annotation'}, 'gff3': {'annotation'},
             'transcripts': {'annotation'},
             'est': {'ests and est clusters', 'transcriptome'}}
    ft2fh = {'fna$masked': ['genome assembly (masked)', 
                            'assembled scaffolds (masked)'],
             'fna$unmasked': ['assembled scaffolds (unmasked)',
                            'genome assembly (unmasked)'],
             'gff': ['genes'], 'gff3': ['genes'], 
             'transcripts': ['transcripts'],
             'est': ['ests', 'transcriptome assembly']}

    ft2fn = {'fna$masked': ['masked', 'Genome Assembly (masked)'], 
             'fna$unmasked': ['AssembledScaffolds', 'scaffolds',
                              'Genome Assembly (unmasked)',
                              'AssemblyScaffolds'],
             'gff3': ['GeneCatalog', 'FilteredModels'],
             'transcripts': ['transcripts'],
             'est': ['EST']}
    ft2fe = {'fna$masked': {'fasta', 'fa', 'fna', 'fsa'},
             'fna$unmasked': {'fasta', 'fa', 'fna', 'fsa'},
             'gff': {'gff', 'gff3'}, 'gff3': {'gff', 'gff3'},
             'transcripts': {'fa', 'fasta', 'fna', 'fsa'},
             'est': {'fa', 'fasta', 'fna', 'fsa'}}

    url, md5, filename = None, False, None

    # parse the XML file
    tree = ET.parse(xml_file)
    root = tree.getroot()
    flip = True
    org_name = None
    has_flipped = False
    attempt = 0

    # flip is a way to rerun the loop if the file type changes (e.g. from
    # masked to unmasked); parse through the XML hiearchy in accord with the
    # hashes established above
    while flip and attempt < 10:
        attempt += 1
        for child in root:
            # conserved subdirectory we need
            if 'Files' == child.attrib['name']:
                for chil1 in child:
                    # does this subdirectory match what we need for our
                    # filetype?
                    if chil1.attrib['name'].lower() in ft2xt[ft]:
                        # we only want the filtered models for annotations/RNA
                        if ft in {'gff3', 'gff', 'transcripts', 'est'}:
                            for chil2 in chil1:
                                if chil2.attrib['name'].lower().startswith('filtered models ('):
                                    chil1 = chil2
                                    break
                        # continue parsing toward the files of interest
                        for chil2 in chil1:
                            if any(x == chil2.attrib['name'].lower() \
                                    for x in ft2fh[ft]):
                                for chil3 in chil2:
                                    t_url = chil3.attrib['url']
                                    try:
                                        org_name = chil3.attrib['label']
                                    except KeyError:
                                        pass
                                    # we want to avoid tape files as we cannot
                                    # download them readily
                                    if 'get_tape_file' not in t_url \
                                        and t_url not in forbidden:
                                        if all(x not in chil3.attrib['filename'] \
                                               for x in ft2fn[ft]):
                                            continue
                                        file_ext_srch = \
                                            re.search(r'\.([^\.]+)$', 
                                                chil3.attrib['filename'])
                                        if file_ext_srch is not None:
                                            file_ext = file_ext_srch[1]
                                            if file_ext == 'gz':
                                                file_ext_srch = \
                                                    re.search(r'\.([^\.]+)\.gz$',
                                                       chil3.attrib['filename'])
                                                if file_ext_srch is not None:
                                                    file_ext = file_ext_srch[1]
                                            if file_ext not in ft2fe[ft]:
                                                continue
                                        else:
                                            continue

                                        url = chil3.attrib['url']
                                        filename = chil3.attrib['filename']
                                        # all requirements satisfied
                                        try:
                                            md5 = chil3.attrib['md5']
                                            break
                                        # continue on to find an md5, or omit
                                        # if exhaustively searched
                                        except KeyError:
                                            pass

        # if the unmasked genome is not present, then query for the masked and
        # vice versa
        if not url and ft == 'fna$masked':
            ft = 'fna$unmasked'
            if not has_flipped:
                flip = True 
            else:
                flip = False
        elif not url and ft == 'fna$unmasked':
            ft = 'fna$masked'
            if not has_flipped:
                flip = True
            else:
                flip = False
        else:
            break

    return filename, url, md5, org_name


def handle_redirect_307(dwnld_data, dwnld, dwnld_url, file_type, xml_file,
                       masked, url, urls, spacer):
    """Handle a redirection error by identifying a new file URL to download
    from, or return the original if none exist"""
    print(spacer +'\t' + dwnld + ' link has moved. ' \
        + 'Trying a different link.', flush = True)
    filename, n_url, dwnld_md5, t_org_name = parse_xml(file_type, xml_file, 
                                         masked = masked,
                                         forbidden = {url}.union(urls))

    if n_url:
        url = n_url
        dwnld_url = prefix + url.replace('&amp;', '&')
        dwnld = f'{output}{file_type}/{os.path.basename(dwnld_url)}'
    return url, dwnld_url, dwnld, {url}.union(urls), t_org_name


def no_md5_checks(dwnld, md5, spacer):
    """If there is no MD5, simply check the file has content in it"""
    check_size = subprocess.run( 
            ['wc', '-l', dwnld], 
            stdout = subprocess.PIPE 
            )
    check_size_res = check_size.stdout.decode('utf-8')
    print(spacer + '\t\tFile exists - no md5 to check.', flush = True)
    check_size_find = re.search(r'\d+', check_size_res)
    size = check_size_find[0]
    if int(size) < 10:
        print(spacer + '\tInvalid file size.', flush = True)
    else:
        md5 = None
    return md5


def jgi_dwnld(ome, file_type, output, masked = True, spacer = '\t'):
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
    circumnavigation of JGI's cryptic maximum ping / time before booting.'''

    # prepare data structures
    prefix = 'https://genome.jgi.doe.gov'
    xml_file = f'{output}xml/{ome}.xml'
    preexisting, check = False, 1

    # stop gap for legacy input
    if file_type == 'gff':
        file_type += '3'

    ran_dwnld = False
    org_name = None

    # acquire the filename, URL, and MD5 from the xml for the file type of
    # interest
    filename, url, dwnld_md5, t_org_name = parse_xml(file_type, xml_file, 
                                         masked = masked)
    if t_org_name:
        org_name = t_org_name
    if not dwnld_md5:
        dwnld_md5 = None

    # if there is a URL present, begin the downloading process
    if url:
        f_urls = {url}
        md5 = False 
        attempt = 0
        curl_cmd = 420

        dwnld_url = prefix + url.replace('&amp;', '&')
   
        dwnld = output + file_type + '/' + os.path.basename(dwnld_url)
        # if the file currently exists, then check its MD5
        if os.path.exists(dwnld):
            if dwnld_md5:
                md5_cmd = subprocess.run([
                        'md5sum', dwnld], 
                        stdout = subprocess.PIPE, stderr = subprocess.PIPE)
                md5_res = md5_cmd.stdout.decode('utf-8')
                md5_find = re.search(r'\w+', md5_res)
                md5 = md5_find[0]

            if md5 == dwnld_md5:
                curl_cmd = 0
                check = dwnld
                preexisting = True
            # if the MD5 does not equal the download MD5 then check the file
            else:
                while True:
                    try:
                        with open(dwnld, 'r') \
                            as dwnld_data_raw:
                            dwnld_data = dwnld_data_raw.read()
                        if re.search('307 Temporary Redirect', dwnld_data):    
                            url, dwnld_url, dwnld, f_urls, t_org_name \
                                = handle_redirect_307(dwnld_data, 
                                                      dwnld, dwnld_url,
                                                      file_type, xml_file,
                                                      masked, url, f_urls, spacer)
                            if t_org_name:
                                org_name = t_org_name
                        break
                    except FileNotFoundError:
                        md5 = False
                        break
                    except UnicodeDecodeError:
                        if not dwnld_md5:
                            md5 = no_md5_checks(dwnld, md5, spacer)
                            if not md5:
                                preexisting = True
                                check = dwnld
                        else:
                            print(spacer + '\tmd5 does not match.', flush = True)
                        break

        # while the MD5 doesn't match, or there is a curl error, try up to 3
        # times to download the file
        while md5 != dwnld_md5 and curl_cmd != 0 and attempt < 3:
            attempt += 1
            curl_cmd = subprocess.call( 
                    ['curl', dwnld_url, '-b', 'cookies', '-o', 
                    dwnld], 
                    stdout = subprocess.PIPE, stderr = subprocess.PIPE 
                    )
            ran_dwnld = True

            if curl_cmd == 0:
                check = dwnld

                # acquire the MD5
                if dwnld_md5:
                    md5_cmd = subprocess.run( 
                            ["md5sum", dwnld], 
                            stdout = subprocess.PIPE, stderr = subprocess.PIPE
                            )
                    md5_res = md5_cmd.stdout.decode('utf-8')
                    md5_find = re.search(r'\w+', md5_res)
                    try:
                        md5 = md5_find[0]
                    except TypeError:
                        md5 = False
                        if not os.path.isfile(dwnld):
                            attempt += 1
                            continue
                # if there is no MD5 attempt the crude file check
                if not dwnld_md5:
                    try:
                        with open(dwnld, 'r') \
                            as dwnld_data_raw:
                            dwnld_data = dwnld_data_raw.read()
                        if re.search('307 Temporary Redirect', dwnld_data):
                            t_url, dwnld_url, dwnld, f_urls, t_org_name = \
                                handle_redirect_307(dwnld_data, 
                                                    dwnld, 
                                                    dwnld_url,
                                                    file_type, xml_file,
                                                    masked, url, f_urls,
                                                    spacer)
                            if t_org_name:
                                org_name = t_org_name
                            
                            if t_url == url:
                                print(spacer + '\t\tNo valid alternative', flush = True)
                                attempt = 4
                                break
                            else:
                                url = t_url
                        else:
                            md5 = no_md5_checks(dwnld, md5, spacer)
                            if not md5:
                                break
                    except FileNotFoundError:
                        pass
                    except UnicodeDecodeError:
                        pass
                    # this is slow, and was arbitrarily set to not overping JGI
                    time.sleep(60)

                # the download may have failed, so prepare to retry
                elif md5 != dwnld_md5 and attempt == 1:
                    print(f'{spacer}\tERROR: md5 does not match JGI. ' \
                        + f'Attempt {attempt}', flush = True)
                    curl_cmd = -1
                    check = 2
                    while True:
                        try:
                            with open(dwnld, 'r') \
                                as dwnld_data_raw:
                                dwnld_data = dwnld_data_raw.read()
                            if re.search('307 Temporary Redirect', dwnld_data ):
                                t_url, dwnld_url, dwnld, f_urls, t_org_name = \
                                    handle_redirect_307(dwnld_data, dwnld, 
                                                        dwnld_url,
                                                        file_type, xml_file,
                                                        masked, url, f_urls,
                                                        spacer)
                                if t_org_name:
                                    org_name = t_org_name
                                if t_url == url:
                                    print(spacer + '\t\tNo valid alternative', flush = True)
                                    attempt = 4
                                    break
                            break
                        except FileNotFoundError:
                            pass
                            break
                        except UnicodeDecodeError:
                            pass
                            break
                        time.sleep(60)
                # if there are two fails, attempt a new URL
                elif md5 != dwnld_md5 and attempt == 2:
                    print(f'{spacer}\tERROR: md5 does not match JGI. Attempt {attempt}', 
                          flush = True)
                    curl_cmd = -1
                    filename, n_url, dwnld_md5, t_org_name \
                        = parse_xml(file_type, xml_file,
                                    masked = masked,
                                    forbidden = f_urls)
                    if t_org_name:
                        org_name = t_org_name
                
                    if n_url:
                        url = n_url
                        dwnld_url = prefix + url.replace('&amp;', '&')
                        f_ulrs = {url}.union(f_urls)
                        dwnld = f'{output}{file_type}/{os.path.basename(dwnld_url)}'
                    time.sleep(60)
                    check = 2
                elif md5 != dwnld_md5:
                    print(f'{spacer}\tERROR: md5 does not match JGI. Attempt {attempt}', 
                          flush = True)
                    check = 2
            else:
                print(f'{spacer}\tERROR: Failed to retrieve {file_type}. `curl` error: ' \
                    + f'{curl_cmd}\n{spacer}\tAttempt {attempt}', flush = True)
                check = 2

        # three strikes and the file is out
        if attempt == 3:
            if md5 != dwnld_md5:
                print(spacer + '\tExcluding from database - potential failure', flush = True)
                curl_cmd = 0
            if curl_cmd != 0:
                print(spacer + '\tFile failed to download', flush = True)

    return check, preexisting, file_type, ran_dwnld, org_name


def main( 
    df, output, user, pwd, assembly = True, proteome = False, gff3 = True, 
    transcript = False, est = False, masked = True, spacer = '\t'
    ):
#    pd.options.mode.chained_assignment = None  # default='warn'
    if not 'assembly_acc' in df.columns:
        if len(df.columns) != 1:
            eprint('\nInvalid input. No assembly_acc column and more than one column.', flush = True)
        else:
            ome_col = list(df.columns)[0]
    else:
        ome_col = 'assembly_acc'

    eprint(spacer + 'Logging into JGI', flush = True)
    login_attempt = 0
    while jgi_login(user, pwd) != 0 and login_attempt < 5:
        eprint(spacer + '\tJGI Login Failed. Attempt: ' + str(login_attempt), flush = True)
        time.sleep(5)
        login_attempt += 1
        if login_attempt == 3:
            eprint(spacer + '\tERROR: Failed 3 login attempts.', flush = True)
            sys.exit(100)

    if not os.path.exists(output + '/xml'):
        os.mkdir(output + '/xml')
# perhaps add a counter here, but one that checks if it is actually querying jgi
    print('\nRetrieving `xml` directories', flush = True)
    ome_set, count = set(), 0
    for i, row in tqdm(df.iterrows(), total = len(df)):
        error_check, attempt = True, 0
        while error_check != -1 and attempt < 3:
            attempt += 1
            error_check = retrieve_xml(row[ome_col], output + '/xml')
            if error_check is None:
                time.sleep(1)
                continue
#            elif error_check > 0:
 #               ome_set.add(row[ome_col])
            elif error_check != -1:
                time.sleep(0.3)
        if error_check != -1:
            eprint(f'{spacer}\t{row[ome_col]} failed to retrieve XML', 
                flush = True)
            ome_set.add(row[ome_col])
    
    eprint(f'{spacer}Downloading {len(df)} JGI files\n\t' \
          + 'Maximum rate: 1 file/min' , flush = True)
    
    dwnlds = []
    if assembly:
        dwnlds.append('fna')
    if proteome:
        dwnlds.append('faa')
    if gff3:
        dwnlds.append('gff3')
    if transcript:
        dwnlds.append('transcript')
#        dwnlds.append( 'AllTranscript' )
    if est:
        dwnlds.append('est')

    for typ in dwnlds:
        if not os.path.isdir(output + '/' + typ):
            os.mkdir(output + '/' + typ)

    preexisting, ran_dwnld = True, False
    for i, row in df.iterrows():
        ome = row[ome_col]
        if ran_dwnld:
            time.sleep(60)
        if ome not in ome_set:
            jgi_login(user, pwd)
            if 'ome' in row.keys():
                eprint(spacer + row['ome'] + '\t' + ome , flush = True)
            else:
                eprint(spacer + ome , flush = True)     
            for typ in dwnlds:
                check, preexisting, new_typ, ran_dwnld, org_name = jgi_dwnld(
                    ome, typ, output, masked = masked, spacer = spacer
                )
                if type(check) != int:
                    df.at[i, new_typ + '_path'] = output + '/' + new_typ + \
                        '/' + os.path.basename(check)
                    check = os.path.basename(os.path.abspath(check))
                    if org_name:
                        org_d = org_name.split()
                        genus = org_d[0]
                        if len(org_d) > 1:
                            sp = org_d[1]
                        else:
                            sp = 'sp.'
                        if len(org_d) > 2:
                            strain = ''.join(org_d[2:])
                        else:
                            strain = ''
                    else:
                        genus, sp, strain = '', '', ''
                    df.at[i, 'genus'] = genus
                    df.at[i, 'species'] = sp
                    df.at[i, 'strain'] = strain
                elif type(check) == int:
                    ome_set.add(row[ome_col])
                eprint(spacer + '\t' + new_typ + ': exit status ' + str(check), flush = True)
        else:
            eprint(spacer + ome + ' failed.' , flush = True)

    if os.path.exists('cookies'):
        os.remove('cookies')
    if os.path.exists(os.path.expanduser('~/.null')):
        os.remove(os.path.expanduser('~/.null'))

    if 'gff3' in df.columns:
        del df['gff3']
    if 'faa' in df.columns:
        del df['faa']
    if 'fna' in df.columns:
        del df['fna']

    return df, ome_set


def cli():

    parser = argparse.ArgumentParser( description = \
        "Imports table/database with JGI `assembly_acc` column and downloads assembly, proteome, gff, " + \
        'and/or gff3. This script supports rerunning/continuing previous runs in the same directory. ' + \
        'JGI has stringent, nondefined ping limits, so file downloads are limited to 1 per minute.' )
    parser.add_argument('-i', '--input', required = True, help = 'Genome code or table with `assembly_acc` column of JGI ome codes')
    parser.add_argument('-a', '--assembly', default = False, action = 'store_true', \
        help = 'Download assembly fastas')
    parser.add_argument('-p', '--proteome', default = False, action = 'store_true', \
        help = 'Download proteome fastas')
#    parser.add_argument( '-g', '--gff', default = False, action = 'store_true', \
 #       help = 'Download gffs.' )
    parser.add_argument('-g', '--gff', default = False, action = 'store_true', \
        help = 'Download gff3s')
    parser.add_argument('-t', '--transcript', default = False, action = 'store_true', \
        help = 'Download transcripts fastas')
    parser.add_argument('-e', '--est', default = False, action = 'store_true', \
        help = 'Download EST fastas')
    parser.add_argument('--nonmasked', default = False, action = 'store_true', \
        help = "[-a] Download nonmasked assemblies")
    parser.add_argument('-o', '--output', default = os.getcwd(), help = 'Output dir')
    args = parser.parse_args()

    if args.nonmasked:
        args.assembly = True

    if not args.assembly and not args.proteome \
        and not args.transcript \
        and not args.est and not args.gff:
        eprint('\nERROR: You must choose at least one download option.' , flush = True)

    ncbi_email, ncbi_api, user, pwd = loginCheck( ncbi = False )
#    user = input( 'JGI username: ' )
 #   pwd = getpass.getpass( prompt='JGI Login Password: ' )

    args_dict = {
            'JGI Table': args.input,
            'Assemblies': args.assembly,
            'RepeatMasked': not args.nonmasked,
            'Proteomes': args.proteome,
            ".gff3's": args.gff,
            'Transcripts': args.transcript,
            'EST': args.est
            }

    start_time = intro('Download JGI files', args_dict)
    eprint('\nWARNING: This script does NOT account for use-restricted data. ' \
         + 'It is user responsibility to determine use restriction status ' \
         + 'in accord with the MycoCosm terms and conditions: ' \
         + 'https://jgi.doe.gov/user-programs/pmo-overview/policies/legacy-data-policies/',
         flush = True)
    eprint(flush = True)
   
    if os.path.isfile(args.input): 
        with open(args.input, 'r') as raw:
            for line in raw:
                if 'assembly_acc' in line.rstrip().split('\t'):
                    df = pd.read_csv(args.input, sep = '\t', index_col = None)
                else:
                    df = pd.read_csv(args.input, sep='\t', header = None)
                break
    else:
        in_data = args.input.replace('"','').replace("'",'').replace(',',' ').split()
        df = pd.DataFrame({'assembly_acc': in_data})

    output = format_path(args.output)

    jgi_df, ome_set = main( 
        df, output, user, pwd, args.assembly, args.proteome, 
        args.gff, args.transcript, args.est, not args.nonmasked,
        spacer = ''
        )
    jgi_df = jgi_df.rename(columns = {'assembly_acc': '#assembly_acc'})
    jgi_df['source'] = 'jgi'
    jgi_df['restriction'] = 'no'
    jgi_df.to_csv(os.path.normpath(args.input) + '.predb.tsv', sep = '\t', index = False)

    outro(start_time)


if __name__ == '__main__':
    cli()
