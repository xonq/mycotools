#! /usr/bin/env python3

import argparse, os, sys, subprocess, re, shutil, numpy as np
from mycotools.lib.kontools import outro, intro, eprint, gunzip, formatPath
from mycotools.lib.fastatools import dict2gff
from mycotools.lib.dbtools import db2df, df2db, gen_omes, gather_taxonomy, assimilate_tax, masterDB, df2std
from mycotools.db.curGFF3 import main as curGFF3
from mycotools.db.gff2gff3 import main as gff2gff3
from mycotools.db.curProteomes import curate as curProteome
from mycotools.curAnnotation import main as curRogue
from Bio import Entrez
import pandas as pd


def predb2db( pre_db ):

    keys = [
        'genome_code', 'genus', 'biosample', 'strain',
        'proteome_path', 'jgi_gff2_path', 'publication',
        'gff3_path', 'assembly_path', 'source', 'internal_ome'
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
            'ecology': '',
            'eco_conf': '',
            'species': row['species'],
            'proteome': row['proteome_path'],
            'gff': row['jgi_gff2_path'],
            'gff3': row['gff3_path'],
            'assembly': row['assembly_path'],
            'blastdb': 0,
            'taxonomy': row['taxonomy'],
            'source': row['source'],
            'published': row['publication']
        } )

    new_db = df2std(pd.DataFrame( data_dict_list ))

    return new_db


def copyFile( old_path, new_path, ome, typ ):

    try:
        shutil.copy( old_path, new_path )
        return True
    except:
        eprint('\t' + ome + '\tcopy ' + typ + ' failed' )
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
        if not os.path.isfile(new_path):
            if not copyFile(formatPath(temp_path), new_path, ome, typ):
                return False
    else:
        new_path = os.environ[env] + '/' + ome + '.' + typ +  uncur
        if not os.path.isfile(new_path) and os.path.isfile(old_path):
            if not copyFile(formatPath(old_path), new_path, ome, typ):
                return False
        elif not os.path.isfile(new_path):
            return False

    return os.path.basename(new_path)


def main( prepdb, refdb, rogue = False ):

    predb = predb2db( prepdb )
    if 'internal_ome' not in predb.columns:
        predb_omes = gen_omes( predb, reference = refdb )
    else:
        predb_omes = predb
        predb_omes['internal_ome'] = predb['internal_ome']

    for key in ['assembly', 'proteome', 'gff', 'gff3']:
        if key not in predb_omes.columns:
            predb_omes[key] = np.nan

    ## need to multiprocess here
    print('\nCopying to database')
    to_del = []
    for i, row in predb_omes.iterrows():
        gff3 = False
        if not pd.isnull( row['assembly'] ):
            new_path = moveBioFile( row['assembly'], row['internal_ome'], 'fa', 'MYCOFNA' )
            if new_path:
                predb_omes.at[i, 'assembly'] = new_path
            else:
                predb_omes.at[i, 'assembly'] = None
        if not row['assembly'] or pd.isnull(row['assembly']):
            eprint('\t' + row['internal_ome'] + ' no assembly, removing entry')
            predb_omes = predb_omes.drop(i)
            continue

        if not pd.isnull(row['gff3']) and row['gff3']:
            new_path = moveBioFile( row['gff3'], row['internal_ome'], 'gff3', 'MYCOGFF3', uncur = '.uncur' )
            if new_path:
                predb_omes.at[i, 'gff3'] = new_path
                
                if row['source'].lower() in {'ncbi', 'jgi'}:
                    try:
                        new_gff = curGFF3(formatPath('$MYCOGFF3/' + new_path), row['internal_ome'])
                        cur_gff = re.sub(r'\.uncur$', '', new_path)
                        to_del.append(formatPath('$MYCOGFF3/' + new_path))
                        with open( formatPath('$MYCOGFF3/' + cur_gff), 'w' ) as out:
                            out.write( dict2gff(new_gff) )
                        predb_omes.loc[i, 'gff3'] = cur_gff
                    except:
                        eprint('\t' + row['internal_ome'] + ' gff3 failed curation')
                else:
                    gff3 = True    
            else:
                predb_omes.at[i, 'gff3'] = None

        if not pd.isnull(row['gff']) and row['gff']:
            if row['gff'].endswith('.gz'):
                gff_path = gunzip(row['gff'])
                if not gff_path:
                    eprint('\t' + row['internal_ome'] + ' gff2 gunzip failed')
                    continue
            else:
                gff_path = formatPath(row['gff'])
            try:
                new_gff3 = gff2gff3(
                    gff2dict(gff_path), 
                    row['internal_ome'], row['genome_code']
                    )
                new_path = os.environ['MYCOGFF3'] + '/' + row['internal_ome'] + '.gff3'
                with open( new_path, 'w' ) as out:
                    out.write( gff2dict(new_gff3) )
                predb_omes.at[i, 'gff3'] = os.path.basename(new_path)
            except:
                predb_omes.at[i, 'gff3'] = None
            
        if not pd.isnull( row['proteome'] ) and row['proteome']:
            new_path = moveBioFile( row['proteome'], row['internal_ome'], 'aa.fa', 'MYCOFAA', uncur = '.uncur' )
            if new_path:
                predb_omes.at[i, 'proteome'] = new_path
                if row['source'].lower() in {'ncbi', 'jgi'}:
                    try:
                        new_prot = curProteome(formatPath('$MYCOFAA/' + new_path), row['internal_ome'])
                        cur_prot = re.sub(r'\.uncur$', '', new_path)
                        to_del.append(formatPath('$MYCOFAA/' + new_path))
                        with open( formatPath('$MYCOFAA/' + cur_prot), 'w' ) as out:
                            out.write( new_prot )
                        predb_omes.loc[i, 'proteome'] = cur_prot
                    except:
                        eprint('\t' + row['internal_ome'] + ' proteome failed curation')
                elif rogue:
                    if not gff3:
                        eprint('\t' + row['internal_ome'] + ' no gff3, cannot curate headers')
                    else:
                        try:
                            gff, fa, trans_str, failed, flagged = curRogue(
                                os.environ['MYCOGFF3'] + '/' + row['gff3'],
                                os.environ['MYCOFAA'] + '/' + row['proteome'],
                                row['internal_ome']
                                )
                            to_del.append( formatPath('$MYCOGFF3/' + row['gff3'] ) )
                            cur_gff = re.sub(r'\.uncur$', '', row['gff3'])
                            with open(formatPath('$MYCOGFF3/' + cur_gff), 'w') as out:
                                out.write(dict2gff(gff))
                            predb_omes.at[i, 'gff3'] = cur_prot
                            to_del.append( formatPath('$MYCOFAA/' + row['proteome'] ) )
                            cur_prot = re.sub(r'\.uncur$', '', row['proteome'])
                            with open(formatPath('$MYCOFAA/' + cur_prot), 'w') as out:
                                out.write(dict2fasta(fa))
                            predb_omes.at[i, 'proteome'] = cur_prot    
            # NEED TO MAKE ROGUE OMES A THING AND UPDATE THE CONFIG
                        except:
                            eprint('\t' + row['internal_ome'] + ' rogue curation failed. Manually curate')
                            predb_omes = predb_omes.drop(i)
            else:
                predb_omes.at[i, 'proteome'] = None

    for i in to_del:
        os.remove(i)
    del predb_omes['gff']
    return df2std(predb_omes)



if __name__ == '__main__':

    parser = argparse.ArgumentParser( 
        description = 'Takes a `predb` file (or creates one and exits), ' + \
            'curates, moves files to database path, and exports a mycotools db' )
    parser.add_argument( '-p', '--predb', help = '`.predb` file to add data from' )
    parser.add_argument( 
        '-g', '--generate', action = 'store_true', 
        help = 'Specify to generate `.predb` file.' 
        )
    parser.add_argument(
        '-r', '--rogue', action = 'store_true',
        help = 'Curate non-JGI/NCBI entries. MUST BE FUNANNOTATE/ORTHOFILLER outputs'
        )
    parser.add_argument( 
        '-d', '--database', default = masterDB(),
        help = 'Existing database to reference. DEFAULT: masterdb' )
    parser.add_argument( '-e', '--email', help = 'NCBI email to gather taxonomy information' )
    parser.add_argument( '-k', '--key', help = 'NCBI API key' )
    args = parser.parse_args()

    args_dict = { 
        'Database': args.database, 
        '`.predb`': args.predb, 
        'Generate `.predb`': args.generate,
        'Email': args.email, 
        'API': args.key
    }


    if args.generate:
        data = 'genome_code\tgenus\tspecies\tstrain\tproteome_path' + \
                    '\tgff3_path\tjgi_gff2_path\tassembly_path\tsource\tpublication\n'
        with open('template.predb', 'w') as file2write:
            file2write.write(data)
        sys.exit(0)
    elif not args.predb or not args.database:
        eprint('\nERROR: Need a `predb` and reference `db` file.\nExit code 3.')
        sys.exit( 3 )

    refdb = db2df( formatPath(args.database) )
    predb = pd.read_csv(formatPath(args.predb), sep = '\t')

    start_time = intro( '`.predb` to `.db`', args_dict )
    predb_omes_tax = main( args_dict, predb, refdb )    
    df2db( predb_omes_tax, 'new.db' )
