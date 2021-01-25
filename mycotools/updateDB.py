#! /usr/bin/env python3

import argparse, pandas as pd, sys, os, re, subprocess, json, getpass, shutil
from mycotools.lib.dbtools import db2df, df2db, gather_taxonomy, assimilate_tax, masterDB, genConfig
from mycotools.lib.kontools import intro, outro, formatPath, eprint, prep_output
from mycotools.lib.fastatools import fasta2dict
from Bio import Entrez
from mycotools.ncbiDwnld import main as ncbiDwnld
from mycotools.jgiDwnld import main as jgiDwnld
from mycotools.utils.ncbi2db import main as ncbi2db
from mycotools.utils.jgi2db import main as jgi2db
from mycotools.predb2db import main as predb2db


def initDB( init_dir, branch, envs, rogue = False ):
    '''Initialize database in `init_dir`'''

    new_dirs = [
        init_dir + 'data/', init_dir + 'config/', 
        init_dir + 'data/fna/', init_dir + 'data/faa/',
        init_dir + 'data/fna/blastdb/',
        init_dir + 'data/faa/blastdb/', init_dir + 'data/gff3/', 
        init_dir + 'log/'
    ]
    output = prep_output( init_dir, cd = True )
    if not output.endswith('/'):
        output += '/'
    for new_dir in new_dirs:
        if not os.path.isdir( new_dir ):
            os.mkdir( new_dir )

    config = genConfig( branch = branch )
    with open( init_dir + 'config/mycodb.conf', 'w' ) as out:
        json.dump( config, out, indent = 1 )    

    if not rogue:
        if not os.path.isdir( init_dir + 'mycodb' ):
            git_exit = subprocess.call( [
                'git', 'clone', config['repository'], init_dir + 'mycodb'
#                '-b', branch
                ] )
            if git_exit != 0:
                eprint('\nERROR: git clone failed.')
                sys.exit( 2 )
        else:
            print('\nmycodb directory already exists')
# NEED TO ADD GITIGNORE TO GIT
        if not masterDB():
            eprint('\nERROR: no YYYYmmdd.db in ' + formatPath( envs['MYCODB'] ))
            sys.exit( 3 )

    if any( x not in os.environ for x in envs ):

        # don't know how this will perform on non-bash or non-profile/non-bash_profile setups
        # need to change BLAST env var, maybe just insinuate that it is in proteome
        # need to update scripts to reference mycotoolsdb envvar
        env_vars = 'export MYCOFNA=' + envs['MYCOFNA'] + '\n' + \
            'export MYCOFAA=' + envs['MYCOFAA'] + '\n' + \
            'export MYCOGFF3=' + envs['MYCOGFF3'] + '\n' + \
            'export MYCODB=' + envs['MYCODB'] + '\n' + \
            '# end mycotools database\n'

        if os.path.isfile( formatPath( '~/.bash_profile' ) ):
            print('\nExporting environment variables to ~/.bash_profile')
            with open( formatPath( '~/.bash_profile' ), 'r' ) as raw:
                data = raw.read()
            data += '\n\n# myctools database\n' + env_vars
            with open( formatPath( '~/.bash_profile' ), 'w' ) as out:
                out.write( data )
        elif os.path.isfile( formatPath( '~/.profile' ) ):
            print('\nExporting environment variables to ~/.profile')
            with open( formatPath( '~/.profile' ), 'r' ) as raw:
                data = raw.read()
            data += '\n\n# mycotools database\n' + env_vars
            with open( formatPath( '~/.profile' ), 'w' ) as out:
                out.write( data )
        else:
            eprint('\nERROR: environment variables not in PATH and no existing ' + \
                '~/.bash_profile or ~/.profile to append')
            sys.exit( 4 )
    else:
        envs = { envs[x]: os.environ[x] for x in envs }

    return output, config
    

def getMissingEntries( db ):

    missing_db = pd.DataFrame()
    for i, row in db.iterrows():
        if not pd.isnull(row['proteome']):
            if not os.path.isfile(os.environ['MYCOFAA'] + '/' + str(row['proteome'])):
                missing_db = missing_db.append(row)
                continue
        if not pd.isnull(row['gff3']):
            if not os.path.isfile(os.environ['MYCOGFF3'] + '/' + str(row['gff3'])):
                missing_db = missing_db.append(row)
                continue
        if not pd.isnull(row['assembly']):
            if not os.path.isfile(os.environ['MYCOFNA'] + '/' + str(row['assembly'])):
                missing_db = missing_db.append(row)

    return missing_db


def checkAntiSMASH( row, smash_dict, update = True, overwrite = False ):
    '''Check that an antiSMASH dir exists. If not and if the user wants to update, then prepare `antiSMASH` bash script.'''

    gff = None
    if not pd.isnull(row['gff3']):
        gff = os.environ['MYCOGFF3'] + '/' +row['gff3']

    if overwrite and gff and not pd.isnull(row['assembly']):
        smash_dict[ 'omes' ][ row['internal_ome'] ] = ( gff, row['assembly'], True )
        smash_dict[ 'walltime' ] += 30

    elif not os.path.isdir( os.environ['MYCOSMASH'] + '/' + row['internal_ome'] ):
        if update and gff and not pd.isnull(row['assembly']):
            smash_dict[ 'omes' ][ row['internal_ome'] ] = ( gff, row['assembly'], False )           
            smash_dict[ 'walltime' ] += 30

    return smash_dict


def smashScript( smash_dict ):

    hours = round( smash_dict['walltime'] / 60 ) + 1
    smash_str = '#PBS -l walltime=' + str(hours) + ':00:00\n#PBS -l nodes=1:ppn=1\n#PBS -A $INSERTPROJECTNUMBER\n' + \
        '#PBS -o updateSMASH\n\nsource activate antismash\n'
    for ome in smash_dict[ 'omes' ]:
        if smash_dict['omes'][ome][2]:
            smash_str += '\nrm -rf $MYCOSMASH/' + ome


        if not smash_dict['omes'][ome][0].endswith('3'):
            smash_str += '\n'
            smash_str += r"""grep -v \"$(printf '\t')exon$(printf '\t')\" """ + smash_dict['omes'][ome][0] + ' > ' + \
                smash_dict['omes'][ome][0] + '_noexon\nantismash $MYCOFNA/' + smash_dict['omes'][ome][1] + ' --taxon' + \
                ' fungi --cb-knownclusters --clusterhmmer --cpus 1 --skip-zip-file --verbose --genefinding-tool glimmerhmm ' + \
                '--genefinding-gff ' + smash_dict['omes'][ome][0] + '_noexon --logfile $MYCOSMASH/log/' + ome + '.out --output-dir $MYCOSMASH/' + \
                ome + '\nrm ' + smash_dict['omes'][ome][1] + '_noexon'

        else:
            smash_str += '\nantismash $MYCOFNA/' + smash_dict['omes'][ome][1] + ' --taxon fungi --cb-knownclusters ' + \
                '--clusterhmmer --cpus 1 --skip-zip-file --verbose --genefinding-tool glimmerhmm --genefinding-gff ' + \
                smash_dict['omes'][ome][0] + ' --logfile $MYCOSMASH/log/' + ome + '.out --output-dir $MYCOSMASH/' + ome

    return smash_str


def checkEcology( row, funguild, update = True ):
    '''Check that ecology data is present and if not query the funguild database.'''

    code = (0, 0)
    ome = row['internal_ome']
    taxon = row['genus'] + ' ' + row['species']
    if pd.isnull(row['ecology']):
        if update and taxon in funguild.index:
            code = (funguild['guild'][taxon], funguild['confidenceRanking'][taxon])

    return code 


def checkBiofile( row, biotype ):
    '''Check whether or not the proteome has valid headers. First, if the proteome is specified
    as `uncur` or `gzipped` then remove `uncur` or `gunzip`. If the proteome does not exist,
    return code 1. Otherwise, use RE to make sure all headers are `$ome_$accession`. If not, 
    return code 2.'''

    hit = ''
    code = 0
    ome = row['internal_ome']
    if not pd.isnull( row[biotype] ):
        hit = os.environ[biotype.upper()] + '/' + row[biotype]
        if hit.endswith('uncur') or hit.endswith('gz'):
            if os.path.isfile(hit):
                if hit.endswith('.gz.uncur'):
                    gunzipPrep = subprocess.call( 'mv ' + hit + ' ' + hit.replace('.gz.uncur', '.uncur.gz'), shell = True )
                    hit = hit.replace('.gz.uncur', '.uncur.gz')
                if hit.endswith('.gz'):
                    gunzip = subprocess.call( 'gunzip ' + hit, shell = True )
                    hit = re.sub(r'\.gz$', '', hit)
            elif os.path.isfile(hit.replace('.gz','')):
                hit = hit.replace('.gz', '')
            elif os.path.isfile(hit.replace('.uncur','')):
                hit = hit.replace('.uncur','')
        if not os.path.isfile(hit):
            print('\t' + ome + '\n\t\tInvalid ' + biotype + ' dir')
            code = 1
        else:
            if biotype.lower() == 'proteome':
                fastadict = fasta2dict(hit)
                for header in fastadict:
                    if not re.search(r'^' + ome + r'\_.+$', header):
                        code = 2
                        break
            else:
                gffdict = gff2dict( hit )
                for line in gffdict:
                    if re.search( r' alias ' + ome + '_', line['attributes'] ):
                        code = 0
                    elif re.search( r';Alias=' + ome, line['attributes'] ):
                        code = 0
                if code == 2:
                    print('\t' + ome + '\n\t\tInvalid ' + biotype + ' headers')
    else:
        if not os.path.isfile((os.environ[biotype.upper()] + '/' + str(row[biotype])).replace('.gz', '').replace('.uncur', '')):
            print('\t' + ome + '\n\t\tNo ' + biotype)
            code = 3
        else:
            hit = row[ biotype ].replace('.gz','').replace('.uncur','')

    hit = os.path.basename( hit )
    return code, hit


def dupCheck( row, types ):
    '''Check for duplicate entries in the database by their assembly, proteome, and gff3 references.
    For each ome, append its data to types, and if that data is already present, report the omes that have
    the same reference. Does not report omes with no entry for a certain column.'''

    for file_type in types:
        if row[file_type] not in types[file_type]:
            types[file_type][row[file_type]] = [row['internal_ome']]
        else:
            if not pd.isnull(row[file_type]):
                types[file_type][row[file_type]].append(row['internal_ome'])
                print('\t' + str(types[file_type][row[file_type]]) + ' have the same ' + file_type)

    return types




def main():

    parser = argparse.ArgumentParser( description = "Initializes or updates myctools database." )
    parser.add_argument( '-u', '--update', action = 'store_true', help = 'Update' )
    parser.add_argument( '-i', '--init', help = 'New directory to initialize database' )
    parser.add_argument( '-d', '--database', default = masterDB(), help = 'Input database DEFAULT: gitlab_stable' )
    parser.add_argument( '--unstable', action = 'store_true', help = 'Attach to unstable db branch' )
    parser.add_argument( '--rogue', action = 'store_true', help = 'Detach from gitlab database' )
    parser.add_argument( '--forbidden', help = 'File of ignored accessions ("acc\ttaxid\\tsource"); DEFAULT: gitlab_stable' )
    parser.add_argument( '--funguild', help = 'Path to funguild' )
    parser.add_argument( '--blastdb', action = 'store_true', help = 'Create protein blastdbs' ) 
    parser.add_argument( '--nonpublished', action = 'store_true', help = 'Agree to MycoCosm terms and ' + \
        'download nonpublished genomes. Requires `--rogue`' )
    parser.add_argument( '-t', '--taxonomy', action = 'store_true', help = 'Query taxonomy' )
#    parser.add_argument( '-e', '--email', help = 'NCBI email for taxonomy' )
#    parser.add_argument( '-a', '--api', help = 'NCBI API key to expedite queries for taxonomy' )
    parser.add_argument( '-y', '--yes', help = 'Agree to all (config changes necessary for nonpublished)' )

    args = parser.parse_args()

    if args.nonpublished and args.rogue:
        print( '\nIMPORTANT NOTE: This setting exists for the ease of access to nonpublished sequences.\n' + \
            'Though many sequence data are not published, all recent JGI data is obligated to be released\n' + \
            'within 2 years of upload. It is inappropriate and against the terms of JGI MycoCosm to\n' + \
            'disseminate, announce, or publish data from nonpublished genomes without explicit consent\n' + \
            "from the genomes' corresponding author(s). Failure to abide by these guidelines in good\n" + \
            "faith can negatively impact the data producers' profession by scooping analyses from data\n" + \
            "they invested resources (time, materials, personnel, training, and funding) to procure.\n" + \
            'Regardless of your personal persuasion, it is the terms agreed upon when both posting\n' + \
            'and using MycoCosm. Disobeying this directive will jeopardize availability for everyone.\n' + \
            'Please respect this responsibility.' )
        acknowledge = \
            "I will not disseminate, announce, or publish nonpublished data in accord with JGI" + \
            "MycoCosm terms and conditions."
        print( '\n' + acknowledge )
        check = input( '\nDo you agree to honor the terms and conditions of using nonpublished' + \
            'JGI data? If so, rewrite the final statement:\n\n' )
        if check.rstrip().upper() != acknowledge.replace( '\n', ' ' ).rstrip().upper():
            if acknowledge.replace( '\n', ' ' ).startswith( check ):
                eprint( '\nCopy and paste prohibited.' )
            else:
                eprint('\nDoes not match. Agreement declined')
            sys.exit( 6 )
    elif args.nonpublished and not args.rogue:
        eprint( '\nERROR: --rogue not called' )
        sys.exit( 5 )

    print()
    ncbi_email = input( 'NCBI email: ' )
    ncbi_api = getpass.getpass( prompt = 'NCBI api key (blank if none): ' )
    jgi_email = input( 'JGI email: ' )
    jgi_pwd = getpass.getpass( prompt = 'JGI password (required): ' )


    if args.unstable:
        branch = 'unstable'
    else:
        branch = 'stable'

    db_path = formatPath( args.database )
    args_dict = { 
        'Database': db_path, 'Update': args.update, 'Initialize': args.init,
        'Branch': branch , 'Rogue': args.rogue, 'Published': not args.nonpublished
    }

    start_time = intro('Update database', args_dict)
    date = start_time.strftime('%Y%m%d')

    if args.init:
        init_dir = formatPath( args.init )
        if not init_dir.endswith( '/' ):
             init_dir += '/'
        envs = { 
            'MYCOFNA': init_dir + 'data/fna', 
            'MYCOFAA': init_dir + 'data/faa', 
            'MYCOGFF3': init_dir + 'data/gff3', 
            'MYCODB': init_dir + 'mycodb'
            }
        os.environ['MYCODB'] = init_dir + 'mycodb'
        output, config = initDB( init_dir, branch, envs, rogue = args.rogue )
        for env in envs:
            os.environ[ env ] = envs[ env ]


        args.database = masterDB()      
    else:
# NEED TO MAKE readDBconfig
        output, config = formatPath('$MYCODB/..'), readDBconfig()

    db_path = formatPath( args.database )
    if db_path:
        db = db2df( db_path )
        missing_db = getMissingEntries( db )
    else:
        db = None

    update_path = output + 'log/' + date + '/'
    if not os.path.isdir( update_path ):
        os.mkdir( update_path )

# NEED TO OUTPUT TO CONFIGURED DB IF NOT THERE
# NEED TO CHECK IF ROGUE HERE VIA CONFIG
# NEED TO PULL EARLIER TO GRAB THE NEWEST DB
    if not config['rogue']:
        git_pull = subprocess.call( [
            'git', 'pull', '-C', output + 'mycodb', config['repository'], '-B', branch
                ], stdout = subprocess.PIPE, stderr = subprocess.PIPE 
            )
        missing_db_ncbi = missing_db[missing_db['source'] == 'ncbi']
        missing_db_jgi = missing_db[missing_db['source'] == 'jgi']
        ncbi_predb = ncbiDwnld( 
            ncbi_df = missing_db_ncbi, email = ncbi_email, api = ncbi_api,
            output_path = update_path
            )
        ncbi_db = predb2db( ncbi_predb, db )
        df2db( ncbi_db, update_path + date + '.ncbi.db' )
        jgi_predb = jgiDwnld(
            missing_db_jgi, update_path, jgi_email, jgi_pwd
            )
        jgi_db = predb2db( jgi_predb, db )
        df2db( jgi_db, update_path + date + '.jgi.db' )
        for i in ['assembly', 'proteome', 'gff3', 'gff']:
            if os.path.isdir( update_path + i ):
                shutil.rmtree( update_path + i )
        
    elif config['rogue']:
# NEED preference check (ncbi or jgi or allow all)
# NEED to mark none for new databases' refdb
        jgi_db = jgi2db( None, db, update_path )
        if type(db) is not None:
            db = pd.concat([jgi_db, db])
        else:
            db = jgi_db
        ncbi_db = ncbi2db( 
            update_path, email = email, api = api,
            ref_db = db
        )
        db = pd.concat([db, ncbi_db])
        missing_db = pd.DataFrame()
        for i, row in db.iterrows():
             if type(row['taxonomy']) != dict:
                 missing_db = missing_db.append( row )
        tax_dicts = gather_taxonomy( missing_db, api_key = api )
        db = assimilate_tax( db, tax_dicts )
        if args.funguild:
            funguild = pd.read_csv(args.funguild, sep ='\t')
            funguild = funguild.set_index('taxon')

        dupFiles = { 'assembly': {}, 'proteome': {}, 'gff3': {} }
        print('\nChecking database ...')
        for i, row in db.iterrows():
            eco, conf = checkEcology( row, funguild, update = args.update )
            if eco != 0:
                db.at[i, 'ecology'] = eco
                db.at[i, 'eco_conf'] = conf
            dupFiles = dupCheck( row, dupFiles )
            for biotype in [ 'gff3', 'proteome']:
                hit, file_ = checkBiofile( row, biotype )
                if file_ != '' and file_ != row[ biotype ]:
                    db.at[i, biotype] = file_
                if hit == 2 and args.update:
                    hiteome = os.environ[biotype.upper()] + '/' + db[biotype][i]
                    if biotype == 'proteome':
                        new_name = curate(hiteome, row['internal_ome'])
                    else:
                        if new_name != hiteome:
                            os.rename( hiteome, new_name )
                        with open( new_name, 'w' ) as out:
                            out.write( dict2gff(new_gff) ) 
                    if new_name:
                        db.at[i, biotype] = os.path.basename( os.path.abspath( new_name ))

    outro(start_time)


if __name__ == '__main__':
    main()
