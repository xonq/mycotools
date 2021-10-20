#! /usr/bin/env python3

import argparse, os, sys, subprocess, re, shutil, datetime, numpy as np
from mycotools.lib.kontools import outro, intro, eprint, gunzip, formatPath
from mycotools.lib.fastatools import dict2gff, gff2dict, fasta2dict, dict2fasta
from mycotools.lib.dbtools import db2df, df2db, gen_omes, masterDB, df2std, loginCheck, gather_taxonomy, assimilate_tax
from mycotools.utils.curGFF3 import main as curGFF3
from mycotools.utils.gff2gff3 import main as gff2gff3
from mycotools.utils.curProteomes import curate as curProteome
from mycotools.gff2seq import aamain as gff2prot
from mycotools.curAnnotation import main as curRogue
from Bio import Entrez
import pandas as pd


def fixVerDate(db):

    db['version'] = db['version'].astype(str)
    db['version'] = db['version'].replace('-','').replace(' 00:00:00','')

    return db


def checkGff3(gff_dict, ome):

    mrna = [x for x in gff_dict if x['type'].lower() == 'mrna']
    if mrna:
        try:
            aliases = [re.search(r'Alias=([^;]+)', x['attributes'])[1] for x in mrna]
            if not all(x.startswith(ome) for x in aliases):
                old_ome = re.search(r'([^_]*)', aliases[0])[1]
                return set(aliases), old_ome
            else:
                return set(aliases), ome
        except IndexError:
            pass


def checkProt(fa_dict, aliases):

    fa_accs = list(fa_dict.keys())
    return set(fa_accs).difference(aliases)

    
def predb2db( pre_db ):

    start_time = datetime.datetime.now()
    date = start_time.strftime('%Y%m%d')
    keys = [
        'genome_code', 'genus', 'biosample', 'strain',
        'proteome_path', 'jgi_gff2_path', 'publication',
        'gff3_path', 'assembly_path', 'source', 'internal_ome',
        'version', 'taxonomy'
        ]
    data_dict_list = []
    for key in keys:
        if key not in pre_db.columns:
            pre_db[key] = ''
    for i, row in pre_db.iterrows():
        if pd.isnull( row['species'] ):
            row['species'] = 'sp.'
        data_dict_list.append( {
            'internal_ome': row['internal_ome'],
            'genome_code': row['genome_code'],
            'genus': row['genus'],
            'biosample': row['biosample'],
            'strain': row['strain'],
            'version': row['version'],
            'ecology': '',
            'eco_conf': '',
            'species': row['species'],
            'proteome': row['proteome_path'],
            'gff': row['jgi_gff2_path'],
            'gff3': row['gff3_path'],
            'assembly': row['assembly_path'],
            'taxonomy': row['taxonomy'],
            'source': row['source'],
            'published': row['publication'],
            'acquisition_date': int(date)
        } )

    new_db = pd.DataFrame( data_dict_list )

    return new_db


def copyFile( old_path, new_path, ome, typ ):

    try:
        shutil.copy( old_path, new_path )
        return True
    except:
        eprint('\t' + ome + '\tcopy ' + typ + ' failed' , flush = True)
        return False


def moveBioFile( old_path, ome, typ, env, uncur = '' ):

    if old_path.endswith('.gz'):
        if not os.path.isfile(old_path[:-3]):
            temp_path = gunzip(old_path)
            if temp_path:
                new_path = os.environ[env] + '/' + ome + '.' + typ +  uncur
            else:
                return False
        else:
            new_path = os.environ[env] + '/' + ome + '.' + typ + uncur
            temp_path = old_path[:-3]
        if not copyFile(formatPath(temp_path), new_path, ome, typ):
            return False
    else:
        new_path = os.environ[env] + '/' + ome + '.' + typ +  uncur
        if not os.path.isfile(new_path) and os.path.isfile(old_path):
            if not copyFile(formatPath(old_path), new_path, ome, typ):
                return False
        elif not os.path.isfile(new_path):
            return False
        else: 
            if not copyFile(formatPath(old_path), new_path, ome, typ):
                return False

    return os.path.basename(new_path)


def main( prepdb, refdb ):

    failed, to_del = [], []
    predb = predb2db( prepdb )
   # if 'internal_ome' not in predb.columns:
    predb_omes = gen_omes( predb, reference = refdb )
# this is a copout. there will be rare occassions of duplicate updates at the same time
# in those cases, the same internal_ome will be upgraded twice. therefore, to fix this
# I will have to identify those dual updates. At present, this will be mitigated through
# a second update
    predb_omes = predb_omes.drop_duplicates('internal_ome')
#    else:
 #       predb_omes = predb
  #      predb_omes['internal_ome'] = predb['internal_ome']

    for key in ['assembly', 'proteome', 'gff', 'gff3']:
        if key not in predb_omes.columns:
            predb_omes[key] = np.nan

    ## need to multiprocess here
    print('\nCopying to database', flush = True)
    for i, row in predb_omes.iterrows():
        gff3, proteome = False, False
        if not pd.isnull( row['assembly'] ) and row['assembly']:
            new_path = moveBioFile( row['assembly'], row['internal_ome'], 'fa', 'MYCOFNA' )
            if new_path:
                predb_omes.at[i, 'assembly'] = new_path
            else:
                if row['source'] == 'ncbi':
                    failed.append([row['biosample'], row['version']])
                else:
                    failed.append([row['genome_code'], row['version']])
                to_del.append(i)
                continue
        else:
            eprint('\t' + row['internal_ome'] + ' no assembly, removing entry', flush = True)
            if row['source'] == 'ncbi':
                failed.append([row['biosample'], row['version']])
            else:
                failed.append([row['genome_code'], row['version']])
            to_del.append(i)
            continue

        if not pd.isnull(row['gff3']) and row['gff3']:
            new_path = moveBioFile( row['gff3'], row['internal_ome'], 'gff3', 'MYCOGFF3', uncur = '.uncur' )
            predb_omes.at[i, 'gff3'] = new_path
            if new_path:
                if row['source'].lower() in {'ncbi', 'jgi'}:
                    try:
                        new_gff = curGFF3(formatPath('$MYCOGFF3/' + new_path), row['internal_ome'])
                        cur_gff = re.sub(r'\.uncur$', '', new_path)
                        with open( formatPath('$MYCOGFF3/' + cur_gff), 'w' ) as out:
                            out.write( dict2gff(new_gff) )
                        predb_omes.loc[i, 'gff3'] = cur_gff
                        os.remove(formatPath('$MYCOGFF3/' + new_path))
                    except KeyboardInterrupt:
                        sys.exit(-1)
                    except:
                        eprint('\t' + row['internal_ome'] + ' gff3 failed curation', flush = True)
                        if row['source'].lower() == 'ncbi':
                            failed.append([row['biosample'], row['version']])
                        else:
                            failed.append([row['genome_code'], row['version']])
                        to_del.append(i)
                        continue
                else:
                    gff, fa = gff2dict(formatPath('$MYCOGFF3/' + new_path)), None
                    aliases, old_ome = checkGff3(gff, row['internal_ome'])
                    if old_ome != row['internal_ome']:
                        with open(formatPath('$MYCOGFF3/' + new_path), 'r') as raw:
                            data = raw.read()
                        with open(formatPath('$MYCOGFF3/' + new_path), 'w') as out:
                            out.write(re.sub(old_ome, row['internal_ome'], data))
                        gff = gff2dict(formatPath('$MYCOGFF3/' + new_path))
                    if aliases:
                        if not pd.isnull(row['proteome']) and row['proteome']:
                            fa = fasta2dict(formatPath(row['proteome']))
                            acc_check = checkProt(fa, aliases)
                        else:
                            assembly = fasta2dict(formatPath('$MYCOFNA/' + row['assembly']))
                            fa = gff2prot(gff, assembly)
                            acc_check = checkProt(fa, aliases)
                        if acc_check:
                            eprint('\t' + row['internal_ome'] + ' discrepant curation', flush = True)
                            failed.append([row['genome_code'], row['version']])
                            os.remove(formatPath('$MYCOGFF3/' + new_path))
                            predb_omes.at[i, 'gff3'] = None
                            predb_omes.at[i, 'proteome'] = None
                        else:
                            os.rename(
                                formatPath('$MYCOGFF3/' + new_path), 
                                formatPath('$MYCOGFF3/' + row['internal_ome'] + '.gff3')
                                )
                            predb_omes.at[i, 'gff3'] = row['internal_ome'] + '.gff3'
                            if old_ome == row['internal_ome']:
                                with open(formatPath('$MYCOFAA/' + row['internal_ome']) + '.aa.fa', 'w') as out:
                                    out.write(dict2fasta(fa))
                            else:
                                 with open(formatPath('$MYCOFAA/' + row['internal_ome']) + '.aa.fa', 'w') as out:
                                     out.write(re.sub(old_ome, row['internal_ome'], dict2fasta(fa)))
                            predb_omes.at[i, 'proteome'] = row['internal_ome'] + '.aa.fa'
                            proteome = True
                    else:
                        try:
                            gff, fa, trans_str, rog_failed, flagged = curRogue(
                                os.environ['MYCOGFF3'] + '/' + row['gff3'],
                                os.environ['MYCOFNA'] + '/' + row['assembly'],
                                row['internal_ome']
                                )
                            cur_gff = re.sub(r'\.uncur$', '', row['gff3'])
                            with open(formatPath('$MYCOGFF3/' + cur_gff), 'w') as out:
                                out.write(dict2gff(gff))
                            os.remove( formatPath('$MYCOGFF3/' + row['gff3'] ))
                            predb_omes.at[i, 'gff3'] = cur_gff
                            cur_prot = re.sub(r'\.uncur$', '', row['proteome'])
                            with open(formatPath('$MYCOFAA/' + cur_prot), 'w') as out:
                                out.write(dict2fasta(fa))
                            os.remove( formatPath('$MYCOFAA/' + row['proteome'] ) )
                            predb_omes.at[i, 'proteome'] = cur_prot    
                            proteome = True
                        except KeyboardInterrupt:
                           sys.exit(-1)
                        except:
                            eprint('\t' + row['internal_ome'] + ' rogue curation failed. ' +
                                'Manually curate', flush = True)
                            failed.append([row['genome_code'], row['version']])
                            to_del.append(i)
                            continue
            else:
                if row['source'] == 'ncbi':
                    failed.append([row['biosample'], row['version']])
                else:
                    failed.append([row['genome_code'], row['version']])
                to_del.append(i)
                continue
        else:
            if row['source'] == 'ncbi':
                failed.append([row['biosample'], row['version']])
            else:
                failed.append([row['genome_code'], row['version']])
            to_del.append(i)
            continue       

        if not pd.isnull(row['gff']) and row['gff']:
            if row['gff'].endswith('.gz'):
                if os.path.isfile(row['gff'][:-3]):
                    gff_path = row['gff'][:-3]
                else:
                    gff_path = gunzip(row['gff'])
                    if not gff_path:
                        eprint('\t' + row['internal_ome'] + ' gff2 gunzip failed', flush = True)
                        continue
            else:
                gff_path = formatPath(row['gff'])
            new_path = os.environ['MYCOGFF3'] + '/' + row['internal_ome'] + '.gff3'
            try:
                new_gff3 = gff2gff3(
                    gff2dict(gff_path), 
                    row['internal_ome'], row['genome_code'], verbose = False
                    )
                with open( new_path, 'w' ) as out:
                    out.write( dict2gff(new_gff3) )
                predb_omes.at[i, 'gff3'] = os.path.basename(new_path)
            except KeyboardInterrupt:
                sys.exit(-1)
            except:
                predb_omes.at[i, 'gff3'] = None

        if not pd.isnull( row['proteome'] ) and row['proteome'] and not proteome:
            new_path = moveBioFile( row['proteome'], row['internal_ome'], 'aa.fa', 'MYCOFAA', uncur = '.uncur' )
            if new_path:
                predb_omes.at[i, 'proteome'] = new_path
                if row['source'].lower() in {'ncbi', 'jgi'}:
                    try:
                        new_prot = curProteome(formatPath('$MYCOFAA/' + new_path), row['internal_ome'])
                        cur_prot = re.sub(r'\.uncur$', '', new_path)
                        with open( formatPath('$MYCOFAA/' + cur_prot), 'w' ) as out:
                            out.write( new_prot )
                        predb_omes.loc[i, 'proteome'] = cur_prot
                        os.remove(formatPath('$MYCOFAA/' + new_path))
                    except KeyboardInterrupt:
                        sys.exit(-1)
                    except:
                        eprint('\t' + row['internal_ome'] + ' proteome failed curation', flush = True)
            else:
                predb_omes.at[i, 'proteome'] = None
        elif not proteome:
            prot_path = formatPath('$MYCOFAA/' + row['internal_ome'])
            if not pd.isnull(row['gff3']) and row['gff3']:
                gff = gff2dict(formatPath('$MYCOGFF3/' + row['internal_ome']))
                assembly = fasta2dict(formatPath('$MYCOFNA/' + row['internal_ome']))
                fa = gff2prot(gff, assembly)
                with open(formatPath('$MYCOFAA/' + row['internal_ome'] + '.aa.fa', 'w'), 'w') as out:
                    out.write(dict2fasta(fa))
                predb_omes.at[i, 'proteome'] = row['internal_ome'] + '.aa.fa'
   
 
    to_del.sort(reverse = True)        
    for i in to_del:
        predb_omes = predb_omes.drop(i)
    del predb_omes['gff']
   
    db = df2std(predb_omes)
    db = fixVerDate(db)

    return db, failed


if __name__ == '__main__':

    parser = argparse.ArgumentParser( 
        description = 'Takes a .predb file (or creates one and exits), ' + 
            'curates, moves files to database path, and updates masterDB. ' +  
             'NOTE: "ncbi" & "jgi" should only be specified as the genome source ' + 
             'ONLY if the gff retains the original source format. Mycotools can only ' +
             'curate these gff formats and Funannotate and/or OrthoFiller outputs.'
        )
    parser.add_argument( '-p', '--predb' )
    parser.add_argument( 
        '-g', '--generate', action = 'store_true', 
        help = 'Generate .predb' 
        )
    parser.add_argument( 
        '-d', '--database', default = masterDB(),
        help = 'DEFAULT: masterDB' )
    args = parser.parse_args()

    args_dict = { 
        'Database': args.database, 
        'Predatabase': args.predb, 
        'Generate .predb': args.generate,
    }



    if args.generate:
        data = 'genome_code\tgenus\tspecies\tstrain\tproteome_path' + \
                    '\tgff3_path\tjgi_gff2_path\tassembly_path\tsource\tpublication\n'
        with open('template.predb', 'w') as file2write:
            file2write.write(data)
        sys.exit(0)
    elif not args.predb or not args.database:
        eprint('\nERROR: Need a `predb` and reference `db` file.\nExit code 3.', flush = True)
        sys.exit( 3 )

    refdb = db2df( formatPath(args.database) )
    predb = pd.read_csv(formatPath(args.predb), sep = '\t')

    start_time = intro( '`.predb` to `.db`', args_dict )

    ncbi_email, ncbi_api, jgi_email, jgi_pwd = loginCheck(jgi = False)
    Entrez.email = ncbi_email
    if ncbi_api:
        Entrez.api_key = ncbi_api

    predb_omes_tax, failed = main( predb, refdb ) 
    date = start_time.strftime('%Y%m%d')
    if failed:
        eprint('\nERROR: failures in predb2db - updates not transferred to masterDB')
        sys.exit(79)

    tax_dicts = gather_taxonomy(predb_omes_tax, api_key = ncbi_api)
    new_db = assimilate_tax(predb_omes_tax, tax_dicts)

    if formatPath(args.database) != masterDB():
        df2db( new_db, 'new.db' )
        eprint('\nUpdate outputted to new.db due to -d, manually update')
    else:
        update_path = os.environ['MYCODB'] + '/../log/' + date + '/'
        if not os.path.isdir(update_path):
            os.mkdir(update_path)
        os.rename(masterDB(), update_path + os.path.basename(masterDB()))
        df2db(pd.concat([new_db, refdb]), os.environ['MYCODB'] + '/' + date + '.db')

    outro(start_time)
