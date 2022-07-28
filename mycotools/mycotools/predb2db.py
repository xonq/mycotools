#!/usr/bin/env python3

# NEED to curate assembly and gff references to ome_contig

import os
import re
import sys
import copy
import shutil
import multiprocessing as mp
from collections import Counter
from mycotools.lib.kontools import gunzip, mkOutput, format_path, eprint, vprint
from mycotools.lib.biotools import gff2list, list2gff, fa2dict, dict2fa, \
    gff3Comps, gff2Comps
from mycotools.lib.dbtools import mtdb, masterDB, loginCheck
from mycotools.curAnnotation import main as curRogue
from mycotools.utils.curGFF3 import main as curGFF3
from mycotools.utils.gff2gff3 import main as gff2gff3
from mycotools.gff2seq import aamain as gff2seq

predb_headers = [
    'assembly_accession', 'genus', 'species', 'strain',
    'version', 'biosample',
    'assemblyPath', 'gffPath', 'genomeSource (ncbi/jgi/new)', 
    'useRestriction (yes/no)', 'published'
    ]

def prep_output(base_dir):
    out_dir = mkOutput(base_dir, 'predb2db')
    wrk_dir = out_dir + 'working/'
    dirs = [
        out_dir, wrk_dir, wrk_dir + 'gff3/', wrk_dir + 'fna/', wrk_dir + 'faa/'
        ]
    for dir_ in dirs:
        if not os.path.isdir(dir_):
            os.mkdir(dir_)
    return dirs[:2]

def copy_file( old_path, new_path ):

    try:
        shutil.copy( old_path, new_path )
        return True
    except:
        raise IOError


def move_biofile(old_path, ome, typ, wrk_dir, suffix = '' ):

    if old_path.endswith('.gz'):
        if not os.path.isfile(old_path[:-3]):
            temp_path = gunzip(old_path)
            new_path = wrk_dir + ome + '.' + typ +  suffix
        else:
            new_path = wrk_dir + ome + '.' + typ + suffix
            temp_path = old_path[:-3]
        copy_file(format_path(temp_path), new_path)
    else:
        new_path = wrk_dir + ome + '.' + typ +  suffix
        if not os.path.isfile(new_path) and os.path.isfile(old_path):
            copy_file(format_path(old_path), new_path)
        elif not os.path.isfile(new_path):
            raise IOError
        else:
            copy_file(format_path(old_path), new_path)

    return new_path


def gen_predb():
    example = [
        'Fibsp1', 'Fibularhizoctonia', 'psychrophila', 'CBS',
        '1.0', 'n', '<PATH/TO/ASSEMBLY>', '<PATH/TO/GFF3>', 'jgi',
        'no', '2018'
        ]
    eprint('INSTRUCTIONS: fill in each column with the relevant information and ' + \
        'separate each column by a tab. The predb can be filled in ' + \
        'via spreadsheet software and exported as a tab delimited `.tsv`. ' + \
        'Alternatively, use a plain text editor and separate by tabs. ' + \
        'De novo annotations produced by Funannotate/Orthofiller must be ' + \
        'filled in as "new" for the genomeSource column. The row below is an ' + \
        'example', flush = True)
    outputStr = '#' + '\t'.join(predb_headers)
    outputStr += '\n#' + '\t'.join(example) + '\n'

    return outputStr

def read_predb(predb_path, spacer = '\t'):

    predb, headers = {}, None
    with open(predb_path, 'r') as raw:
        for line in raw:
            if line.startswith('#'):
                if not headers:
                    predb = {x: [] for x in line.rstrip()[1:].split('\t')}
                    headers = list(predb.keys())
            elif line.rstrip():
                entry = line.split('\t')
                if len(entry) != len(headers):
                    eprint(spacer + 'ERROR: all columns must have an entry.', flush = True)
                    eprint(headers, '\n', entry, flush = True)
                    sys.exit(3)
                for i, v in enumerate(entry):
                    predb[headers[i]].append(v.rstrip())
#                (predb[headers[i]].append(v.rstrip()) for i,v in enumerate(entry)) 

    try:
        predb['assembly_acc'] = predb['assembly_accession']
        del predb['assembly_accession']
    except KeyError:
        if not 'assembly_acc' in predb and not 'assembly_accession' in predb:
            predb['assembly_acc'] = ['' for x in predb['genus']]
    try:
        predb['source'] = predb['genomeSource (ncbi/jgi/new)']
        del predb['genomeSource (ncbi/jgi/new)']
    except KeyError:
        if not 'genomeSource' in predb and not 'source' in predb:
            raise KeyError
    try:
        predb['restriction'] = predb['useRestriction (yes/no)']
        del predb['useRestriction (yes/no)']
        if any(x.lower() not in {'y', 'n', 'yes', 'no', ''} for x in predb['restriction']):
            eprint(spacer + 'ERROR: useRestriction entries must be in {y, n, yes, no}', flush = True)
            sys.exit(4)
    except KeyError:
        if not 'restriction' in predb and 'published' not in predb:
            raise KeyError
        elif 'published' in predb:
            predb['restriction'] = [bool(x) for x in predb['published']]

    if any(x.lower() not in {'jgi', 'ncbi', 'new'} \
        for x in predb['source']):
        eprint(spacer + 'ERROR: genoemSource entries must be in {jgi, ncbi, new}', flush = True)
        sys.exit(5)
    for i, path in enumerate(predb['assemblyPath']):
        predb['assemblyPath'][i] = format_path(predb['assemblyPath'][i])
        predb['gffPath'][i] = format_path(predb['gffPath'][i])
        if predb['restriction'][i].lower() in {'y', 'yes'}:
            predb['restriction'][i] = True
        if not predb['species'][i]:
            predb['species'][i] = 'sp.'
        for col in predb.keys():
            if predb[col][i].lower() in {'n', 'no'}:
                predb[col][i] = ''
    return predb

def predb2mtdb(predb):
    infdb = mtdb()
    for i, code in enumerate(predb['genus']):
        if predb['published'][i]:
            toAdd = {
                'assembly_acc': predb['assembly_acc'][i],
                'genus': predb['genus'][i],
                'species': predb['species'][i],
                'strain': predb['strain'][i],
                'version': predb['version'][i],
                'biosample': predb['biosample'][i],
                'fna': predb['assemblyPath'][i],
                'gff3': predb['gffPath'][i],
                'source': predb['source'][i],
                'published': predb['published'][i]
                }
        else:
            toAdd = {
                'assembly_acc': predb['assembly_acc'][i],
                'genus': predb['genus'][i],
                'species': predb['species'][i],
                'strain': predb['strain'][i],
                'version': predb['version'][i],
                'biosample': predb['biosample'][i],
                'fna': predb['assemblyPath'][i],
                'gff3': predb['gffPath'][i],
                'source': predb['source'][i],
                'published': not predb['restriction'][i]
                }
#        if predb['published'][i]:
 #           toAdd['published'] = predb['published'][i]
        infdb = infdb.append(toAdd)
    return infdb

def gen_omes(
    newdb, refdb = None, ome_col = 'ome', forbidden = set(),
    spacer = '\t'
    ):

    t_failed = []
    ref_ome_check = [k for k,v in Counter(refdb['ome']).items() if v > 1]
    new_ome_check = [k for k,v in Counter(newdb['ome']).items() if v > 1 and k]
    if ref_ome_check:
        raise ValueError('corrupted reference database with non-unique omes: '
                        + str(ref_ome_check))
    elif new_ome_check:
        raise ValueError('corrupted pre-database with non-unique omes: '
                        + str(new_ome_check))
    tax_list = list(set(refdb['ome']).union(forbidden))
    tax_count = {}
    for tax in tax_list:
        abb = tax[:6]
        num = int(tax[6:])
        if abb in tax_count:
            if tax_count[abb] < num:
                tax_count[abb] = num
        else:
            tax_count[abb] = num
    tax_count = Counter(tax_count)
    todel, tax_count = [], Counter()
    for i, ome in enumerate(newdb['ome']):
        if not ome: # if there isn't an ome for this entry yet
            try:
                name = newdb['genus'][i][:3].lower() + newdb['species'][i][:3].lower()
            except TypeError:
                todel.append(i)
                if not isinstance(newdb['assembly_acc'][i], float):
                    eprint(spacer + newdb['assembly_acc'][i] + ' no metadata - ' \
                         + 'failed', flush = True)
                elif 'index' in newdb: # for updateDB
                    if not isinstance(newdb, float):
                        eprint(spacer + newdb['index'][i] + ' no metadata - ' \
                            + 'failed', flush = True)
                    continue 
                else: # no use appending failed when there's no identifiable
                # info
                    continue
                row = {key: newdb[key][i] for key in mtdb.columns \
                       if key != 'ome'}
                t_failed.append(add2failed(row))
                continue
            name = re.sub(r'\(|\)|\[|\]|\$|\#|\@| |\||\+|\=|\%|\^|\&|\*|\'|\"|\!|\~|\`|\,|\<|\>|\?|\;|\:|\\|\{|\}', '', name)
            name.replace('(', '')
            name.replace(')', '')
            while len(name) < 6:
                name += '.'
            tax_count[name] += 1
            new_ome = name + str(tax_count[name])
            newdb['ome'][i] = new_ome
        elif ome:
            format_search = re.search(r'^[^\d_,\'";:\\\|\[\]\{\}\=\+\!@#\$\%\^' 
                                      + r'&\*\(\)]{6}\d+[^-_+=\\\|\{\[\}\]\:;'
                                      + r'\'\"\,\<\>\?/\`\~\!\@\#\$\%\^\&\*\('
                                      + r'\)]+', ome) # crude format check
            if not format_search:
                raise ValueError('invalid ome ' + ome)
            if ome in tax_set: # it's an update
                v_search = re.search(r'\.(\d+)$', ome)
                if v_search:
                    v = int(v_search[1]) + 1 # new version
                    new_ome = re.sub(r'\.\d+$', '.' + str(v), ome)
                else:
                    new_ome = ome + '.1' # first modified version
                eprint(spacer + ome + ' update -> ' + new_ome, flush = True)
                newdb['ome'][i] = new_ome

    for i in reversed(todel):
        for key in mtdb.columns:
            del newdb[key][i]

    return newdb, t_failed

def cur_fna(cur_raw_fna_path, uncur_raw_fna_path, ome):
    with open(cur_raw_fna_path + '.tmp', 'w') as out:
        with open(uncur_raw_fna_path, 'r') as in_:
            for line in in_:
                if line.startswith('>'):
                    if not line.startswith('>' + ome + '_'):
                        out.write('>' + ome + '_' + line[1:])
                    else: # already curated
                        out.write(line)
                else:
                    out.write(line)
    os.rename(cur_raw_fna_path + '.tmp', cur_raw_fna_path)

def cur_mngr(ome, raw_fna_path, raw_gff_path, wrk_dir, 
            source, assembly_accession, exit = False,
            remove = False, spacer = '\t'):

    # assembly FNAs
    uncur_fna_path = wrk_dir + 'fna/' + ome + '.fna.uncur'
    cur_fna_path = wrk_dir + 'fna/' + ome + '.fna'
    if not os.path.isfile(cur_fna_path):
        try:
            uncur_fna_path = move_biofile(raw_fna_path, ome, 'fa', wrk_dir + 'fna/',
                                      suffix = '.uncur')
        except IOError:
            eprint(spacer + ome + '|' + assembly_accession \
                 + ' failed FNA parsing', flush = True)
            return ome, False, 'fna'
        cur_fna(cur_fna_path, uncur_fna_path, ome)

    # gene coordinate GFF3s
    uncur_gff_path = wrk_dir + 'gff3/' + ome + '.gff3.uncur'
    cur_gff_path = wrk_dir + 'gff3/' + ome + '.gff3'
    if not os.path.isfile(cur_gff_path):
        try:
            uncur_gff_path = move_biofile(raw_gff_path, ome, 'gff3', 
                                          wrk_dir + 'gff3/', suffix = '.uncur')
        except IOError:
            eprint(spacer + ome + '|' + assembly_accession \
                 + ' failed GFF3 parsing', flush = True)
            return ome, False, 'gff3'
    
        gff = gff2list(uncur_gff_path)
        try:
            gff_mngr(ome, gff, cur_gff_path, source, assembly_accession)
        except: # catch all errors to continue script
            eprint(spacer + ome + '|' + assembly_accession \
                + ' failed GFF3 curation', flush = True)
            if exit:
                sys.exit(17)
            return ome, False, 'gff3'

    # proteome FAAs
    faa_path = wrk_dir + 'faa/' + ome + '.faa'
    if not os.path.isfile(faa_path):
        try:
            faa = gff2seq(gff2list(cur_gff_path), fa2dict(uncur_fna_path),
                          spacer = spacer)
        except: # catch all errors
            eprint(spacer + ome + '|' + assembly_accession \
                 + ' failed proteome generation', flush  = True)
            if exit:
                sys.exit(18)
            return ome, False, 'faa'
        with open(faa_path + '.tmp', 'w') as out:
            out.write(dict2fa(faa))
        os.rename(faa_path + '.tmp', faa_path)

    if remove:
        if os.path.isfile(uncur_gff_path):
            os.remove(uncur_gff_path)
        if os.path.isfile(raw_gff_path):
            os.remove(raw_gff_path)
        if os.path.isfile(re.sub(r'\.gz$', '', raw_gff_path)):
            os.remove(re.sub(r'\.gz$', '', raw_gff_path))
        if os.path.isfile(uncur_fna_path):
            os.remove(uncur_fna_path)
        if os.path.isfile(raw_fna_path):
            os.remove(raw_fna_path)
        if os.path.isfile(re.sub(r'\.gz$', '', raw_fna_path)):
            os.remove(re.sub(r'\.gz$', '', raw_fna_path))

    return ome, cur_fna_path, cur_gff_path, faa_path


def gff_mngr(ome, gff, cur_path, source, assembly_accession):

    gffVer, alias = None, False
    for entry in gff:
        if re.search(gff3Comps()['id'], entry['attributes']):
            gffVer = 3
            alias = re.search(gff3Comps()['Alias'], entry['attributes'])
            break
        elif re.search(gff2Comps()['id'], entry['attributes']):
            gffVer = 2
            break

    if gffVer == 3:
        if source == 'new':
            if alias: #already curated
                aliasOme = alias[1]
                for entry in gff:
                    alias0 = re.search(gff3Comps()['Alias'], entry['attributes'])[1]
                    alias_num = alias0[alias0.find('_') + 1:]
                    new_alias = ome + '_' + alias_num
                    entry['attributes'] = re.sub(
                        gff3Comps()['Alias'], 'Alias='+ new_alias,
                        entry['attributes']
                        )
            else:
                gff, trans_str, failed, flagged = curRogue(gff, ome,
                                                           cur_seqids = True)
        else:
            gff = curGFF3(gff, ome, cur_seqids = True)
    else:
        gff = gff2gff3(gff, ome, assembly_accession, verbose = False,
                       cur_seqids = True)

    with open(cur_path + '.tmp', 'w') as out:
        out.write(list2gff(gff))
    os.rename(cur_path + '.tmp', cur_path)

def add2failed(row):
    if isinstance(row['assembly_acc'], float) or not row['assembly_acc']:
        return False
    else:
        return [row['assembly_acc'], row['version']]

def main(
    predb, refdb, wrk_dir, 
    verbose = False, spacer = '\t', forbidden = set(), 
    cpus = 1, exit = False, remove = False
    ):

    infdb = predb2mtdb(predb)
    omedb, failed = gen_omes(infdb, refdb, ome_col = 'ome', forbidden = forbidden,
                     spacer = spacer)
    
    cur_cmds = []
    omedb = omedb.set_index('ome')
    for ome, row in omedb.items():
        cur_cmds.append([
            ome, row['fna'], row['gff3'], 
            wrk_dir, row['source'], row['assembly_acc'],
            exit, remove, spacer
            ])
    with mp.Pool(processes = cpus) as pool:
        cur_data = pool.starmap(cur_mngr, cur_cmds)

    for data in cur_data:
        if not data[1]:
            failed.append(add2failed(omedb[data[0]]))
            del omedb[data[0]]
        else:
            ome, fna_path, gff3_path, faa_path = data
            omedb[ome]['fna'] = fna_path
            omedb[ome]['gff3'] = gff3_path
            omedb[ome]['faa'] = faa_path

    return omedb.reset_index(), failed


if __name__ == '__main__':
    usage = 'Generate a predb file:\npredb2db.py\n\nCreate a mycotoolsdb ' + \
    'from a predb file:\npredb2db.py <PREDBFILE>\n\nCreate a mycotoolsdb ' + \
    'referencing an alternative master database:\npredb2db.py <PREDBFILE> ' + \
    '<REFERENCEDB>'

    if any(x in {'-h', '--help', '-help'} for x in sys.argv):
        eprint('\n' + usage + '\n', flush = True)
        sys.exit(1)
    elif len(sys.argv) == 1:
        print(gen_predb())
        sys.exit(0)
    elif len(sys.argv) == 3:
        refDB = mtdb(format_path(sys.argv[2]))
    else:
        refDB = mtdb(masterDB())

    from Bio import Entrez
    ncbi_email, ncbi_api, jgi_email, jgi_pwd = loginCheck(jgi = False)
    Entrez.email = ncbi_email
    if ncbi_api:
        Entrez.api_key = ncbi_api

    predb = read_predb(format_path(sys.argv[1]), spacer = '\t')
    out_dir, wrk_dir = prep_output(os.path.dirname(format_path(sys.argv[1])))
    omedb, failed = main(predb, refDB, wrk_dir, exit = True)

    from mycotools.lib.dbtools import gather_taxonomy, assimilate_tax
    tax_dicts = gather_taxonomy(omedb, api_key = ncbi_api)
    outdb = assimilate_tax(omedb, tax_dicts)
    outdb.df2db(out_dir + 'predb2db.db')
    sys.exit(0)
