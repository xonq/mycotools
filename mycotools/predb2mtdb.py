#!/usr/bin/env python3

# NEED source to reference the annotation source
# NEED to error check FAA generation simply by file size

import os
import re
import sys
import copy
import shutil
import multiprocessing as mp
from tqdm import tqdm
from collections import Counter, defaultdict
from mycotools.lib.kontools import gunzip, mkOutput, format_path, eprint, vprint
from mycotools.lib.biotools import gff2list, list2gff, fa2dict, dict2fa, \
    gff3Comps, gff2Comps, gtfComps
from mycotools.lib.dbtools import mtdb, primaryDB, loginCheck
from mycotools.utils.gtf2gff3 import main as gtf2gff3
from mycotools.utils.curGFF3 import main as curGFF3
from mycotools.utils.gff2gff3 import main as gff2gff3
from mycotools.utils.curGFF3 import rename_and_organize as rename_and_organize
from mycotools.gff2seq import aamain as gff2seq

predb_headers = [
    'assembly_accession', 'previous_ome', 
    'genus', 'species', 'strain', 'version', 'biosample',
    'assemblyPath', 'gffPath', 'genomeSource (ncbi/jgi/new)', 
    'useRestriction (yes/no)', 'published'
    ]

def prep_output(base_dir):
    out_dir = mkOutput(base_dir, 'predb2mtdb')
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
            raise IOError(old_path, new_path)
        else:
            copy_file(format_path(old_path), new_path)

    return new_path


def gen_predb():
    example = [
        'Fibsp1', 'fibpsy1', 'Fibularhizoctonia', 'psychrophila', 'CBS',
        '1.0', 'n', '<PATH/TO/ASSEMBLY>', '<PATH/TO/GFF3>', 'jgi',
        'no', '2018'
        ]
    eprint('INSTRUCTIONS: fill in each column with the relevant information and \
        separate each column by a tab. The predb can be filled in \
        via spreadsheet software and exported as a tab delimited `.tsv`. \
        ASSEMBLY ACCESSIONS and PREVIOUS_OME fields must be unique to the \
        genome; otherwise predb2mtdb will update the corresponding database entry. \
        Novel data must be filled in as "new" for the genomeSource column.', flush = True)
    outputStr = '#' + '\t'.join(predb_headers)
    outputStr += '\n#' + '\t'.join(example) + '\n'

    return outputStr

def read_predb(predb_path, spacer = '\t'):

#    predb, headers = {}, None
    predb = defaultdict(list)
    with open(predb_path, 'r') as raw:
        for i, line in enumerate(raw):
#            if line.startswith('#'):
 #               if not headers:
  #                  predb = {x: [] for x in line.rstrip()[1:].split('\t')}
   #                 headers = list(predb.keys())
            if not line.startswith('#') and line.rstrip():
#            if line.rstrip():
                line = line.replace('\n','')
                entry = line.split('\t')
                if len(entry) != len(predb_headers):
                    eprint(spacer + 'ERROR: Incorrect columns, line ' + str(i),
                           flush = True)
                    eprint(predb_headers, '\n', entry, flush = True)
                    sys.exit(3)
                for i, v in enumerate(entry):
                    predb[predb_headers[i]].append(v.rstrip())

    predb = dict(predb)
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
        if any(x.lower() not in {'y', 'n', 'yes', 'no', '', 'true', 'false'} for x in predb['restriction']):
            eprint(spacer + 'ERROR: useRestriction entries must be in {y, n, yes, no}', flush = True)
            sys.exit(4)
    except KeyError:
        if not 'restriction' in predb and 'published' not in predb:
            raise KeyError
        elif 'published' in predb:
            predb['restriction'] = [bool(x) for x in predb['published']]

    if any(x.lower() not in {'jgi', 'ncbi', 'new'} \
        for x in predb['source']):
        eprint(spacer + 'ERROR: genomeSource entries must be in {jgi, ncbi, new}', flush = True)
        sys.exit(5)
    for i, path in enumerate(predb['assemblyPath']):
        predb['assemblyPath'][i] = format_path(predb['assemblyPath'][i])
        predb['gffPath'][i] = format_path(predb['gffPath'][i])
        if predb['restriction'][i].lower() in {'y', 'yes', 'true'}:
            predb['restriction'][i] = True
        elif predb['restriction'][i].lower() in {'n', 'no', 'false'}:
            predb['restriction'][i] = ''
        if not predb['species'][i]:
            predb['species'][i] = 'sp.'
    return predb

def sub_disallowed(data, disallowed = r"""[^\w\d]"""):
    if data:
        return re.sub(disallowed, '', data)
    else:
        return data

def predb2mtdb(predb):
    infdb = mtdb()
    if not 'assemblyPath' in predb and 'fna' in predb:
        predb['assemblyPath'] = predb['fna']
    if not 'gffPath' in predb and 'gff3' in predb:
        predb['gffPath'] = predb['gff3']
    if not 'previous_ome' in predb and 'ome' in predb:
        predb['previous_ome'] = predb['ome']
    elif not 'previous_ome' in predb:
        predb['previous_ome'] = \
            [None for x in predb[list(predb.keys())[0]]]
    for i, code in enumerate(predb['genus']):
        toAdd = {
            'assembly_acc': predb['assembly_acc'][i],
            'ome': predb['previous_ome'][i].lower(),
            'genus': sub_disallowed(predb['genus'][i]),
            'species': sub_disallowed(predb['species'][i]),
            'strain': re.sub(r'[^a-zA-Z0-9]', '', predb['strain'][i]),
            'version': sub_disallowed(predb['version'][i]),
            'biosample': sub_disallowed(predb['biosample'][i]),
            'fna': predb['assemblyPath'][i],
            'gff3': predb['gffPath'][i],
            'source': predb['source'][i]
            }
        if predb['published'][i]:
            toAdd['published'] = predb['published'][i]
        else:
            toAdd['published'] = not predb['restriction'][i]
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
#    new_ome_check = [k for k,v in Counter(newdb['ome']).items() if v > 1 and k]
    if ref_ome_check:
        raise ValueError('corrupted reference database with non-unique omes: '
                        + str(ref_ome_check))
 #   elif new_ome_check:
  #      raise ValueError('corrupted pre-database with non-unique omes: '
   #                     + str(new_ome_check))
    tax_list = list(set(refdb['ome']).union(forbidden))
    tax_count = {}
    for tax in tax_list:
        abb = tax[:6]
        try:
            num = int(tax[6:])
        except ValueError: # version included in number
            num = int(re.search(r'(^\d+)', tax[6:])[1])
        if abb in tax_count:
            if tax_count[abb] < num:
                tax_count[abb] = num
        else:
            tax_count[abb] = num
    tax_count = Counter(tax_count)

    todel = []
    refdb_aas = refdb.set_index('assembly_acc')
    refdb_accs = set([x for x in list(refdb['assembly_acc']) if x])
    refdb_omes = set([x for x in list(refdb['ome'])])
    refdb_nover = {re.sub(r'(^.{6}\d+)\.\d+', r'\1', x): x \
                   for x in list(refdb['ome'])}
    for i, ome in enumerate(newdb['ome']):
        if not ome: # if there isn't an ome for this entry yet
            if newdb['assembly_acc'][i] in refdb_accs: # if this is an
            # established assembly accession; maybe make sure this works for
            # changed MycoCosm or NCBI assembly accs
                ome = refdb_aas[newdb['assembly_acc'][i]]['ome']
                v_search = re.search(r'\.(\d+)$', ome[6:])
                if v_search:
                    v = int(v_search[1]) + 1 # new version
                    new_ome = re.sub(r'\.\d+$', '.' + str(v), ome)
                else:
                    new_ome = ome + '.1' # first modified version
                eprint(spacer + ome + ' update -> ' + new_ome, flush = True)
                newdb['ome'][i] = new_ome
                continue
               
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
            format_search = re.search(r'^[^\d_,\'";:\\\|\[\]\{\}\=\+\!@#\$\%\^' \
                                    + r'&\*\(\)]{6}\d+[^-_+=\\\|\{\[\}\]\:;' \
                                    + r'\'\"\,\<\>\?/\`\~\!\@\#\$\%\^\&\*\(' \
                                    + r'\)\w\W]*\.{0,1}\d*$', ome) # crude format check
            if not format_search:
                raise ValueError('invalid ome ' + ome)
            if ome in refdb_omes: # it's an update
                v_search = re.search(r'\.(\d+)$', ome[6:])
                if v_search:
                    v = int(v_search[1]) + 1 # new version
                    new_ome = re.sub(r'\.\d+$', '.' + str(v), ome)
                else:
                    new_ome = ome + '.1' # first modified version
                eprint(spacer + ome + ' update -> ' + new_ome, flush = True)
                newdb['ome'][i] = new_ome
            elif ome in refdb_nover: # has a version, wasn't given in predb
                version_ome = refdb_nover[ome]
                v_search = re.search(r'\.(\d+)$', version_ome[6:])
                if v_search:
                    v = int(v_search[1]) + 1 # new version
                    new_ome = ome + '.' + str(v)
                else:
                    raise TypeError('unknown error ' + ome)
                eprint(spacer + ome + ' update -> ' + new_ome, flush = True)
                newdb['ome'][i] = new_ome
                    

    for i in reversed(todel):
        for key in mtdb.columns:
            del newdb[key][i]

    return newdb, t_failed

def cur_fna(cur_raw_fna_path, uncur_raw_fna_path, ome):
    ome_ver = re.search(r'(.{6}\d+).(\d+)$', ome)
    if ome_ver:
        less_ome = ome_ver[1]
        ver_num = ome_ver[2]
    else:
        less_ome = ome
        ver_num = 0
    with open(cur_raw_fna_path + '.tmp', 'w') as out:
        with open(uncur_raw_fna_path, 'r') as in_:
            for line in in_:
                if line.startswith('>'):
                    if not line.startswith('>' + ome + '_'):
                        if re.search(r'^>' + less_ome + '_', line) is not None:
                            new_line = line.replace('>' + less_ome + '_', '>' + ome + '_')
                        elif re.search(r'^>' + less_ome + r'\.\d+_', line) is not None:
                            new_line = re.sub(r'^>' + less_ome + r'\.\d+_', '>' + ome + '_',
                                              line)
                        else:
                            new_line = '>' + ome + '_' + line[1:]
                        out.write(new_line)
                    else: # already curated
                        out.write(line)
                else:
                    out.write(line)
    shutil.move(cur_raw_fna_path + '.tmp', cur_raw_fna_path)

def cur_mngr(ome, raw_fna_path, raw_gff_path, wrk_dir, 
            source, assembly_accession, exit = False,
            remove = False, spacer = '\t\t\t', verbose = False):

    predb_dir = os.path.basename(os.path.dirname(wrk_dir[:-1])) + '/working/'

    # assembly FNAs
    vprint('\t' + ome, v = verbose, flush = True)
    uncur_fna_path = wrk_dir + 'fna/' + ome + '.fna.uncur'
    cur_fna_path = wrk_dir + 'fna/' + ome + '.fna'
    vprint('\t\t' + predb_dir + 'fna/' + ome + '.fna', v = verbose, flush = True)
    if not os.path.isfile(cur_fna_path):
        try:
            uncur_fna_path = move_biofile(raw_fna_path, ome, 'fa', wrk_dir + 'fna/',
                                      suffix = '.uncur')
        except IOError as ie:
            eprint(spacer + ome + '|' + assembly_accession \
                 + ' failed FNA parsing', flush = True)
            if exit:
                raise ie from None
            return ome, False, 'fna'
        cur_fna(cur_fna_path, uncur_fna_path, ome)

    # gene coordinate GFF3s
    uncur_gff_path = wrk_dir + 'gff3/' + ome + '.gff3.uncur'
    cur_gff_path = wrk_dir + 'gff3/' + ome + '.gff3'
    vprint('\t\t' + predb_dir + 'gff3/' + ome + '.gff3', v = verbose, flush = True)
    if not os.path.isfile(cur_gff_path):
        try:
            uncur_gff_path = move_biofile(raw_gff_path, ome, 'gff3', 
                                          wrk_dir + 'gff3/', suffix = '.uncur')
        except IOError as ie:
            eprint(spacer + ome + '|' + assembly_accession \
                 + ' failed GFF3 parsing', flush = True)
            if exit:
                raise ie from None
            return ome, False, 'gff3'
        try:   
            gff = gff2list(uncur_gff_path)
        # malformatted
        except IndexError:
            return ome, False, 'gff3'
        try:
            gff_mngr(ome, gff, cur_gff_path, source, assembly_accession)
        except Exception as e: # catch all errors to continue script
            eprint(spacer + ome + '|' + assembly_accession \
                + ' failed GFF3 curation', flush = True)
            if exit:
                raise e from None
            return ome, False, 'gff3'

    # proteome FAAs
    faa_path = wrk_dir + 'faa/' + ome + '.faa'
    vprint('\t\t' + predb_dir + 'faa/' + ome + '.faa', v = verbose, flush = True)
    if not os.path.isfile(faa_path):
        try:
            faa = gff2seq(gff2list(cur_gff_path), fa2dict(cur_fna_path),
                          spacer = spacer)
        except Exception as e: # catch all errors
            eprint(spacer + ome + '|' + assembly_accession \
                 + ' failed proteome generation', flush  = True)
            if exit:
                raise e
            return ome, False, 'faa'
        with open(faa_path + '.tmp', 'w') as out:
            out.write(dict2fa(faa))
        shutil.move(faa_path + '.tmp', faa_path)

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
            if alias is not None:
#                if entry['seqid'].startswith(ome + '_'):
                alias = True
                #else:
                 #   alias = False # remove the old aliases
                  #  entry['attributes'] = re.sub(r';?Alias=[^;]+', '',
                  #                              entry['attributes'])
            else:
                break
        elif re.search(gtfComps()['id'], entry['attributes']):
            gffVer = 2.5
            break
        elif re.search(gff2Comps()['id'], entry['attributes']):
            gffVer = 2
            break

    if gffVer == 3:
#        if source == 'new':
        if alias: #already curated
            try:
                new_gff = copy.deepcopy(gff)
                old_ome_p = re.search(gff3Comps()['Alias'], new_gff[0]['attributes'])[1]
                old_ome = old_ome_p[:old_ome_p.find('_')]
                for entry in new_gff:
#                    alias0 = re.search(gff3Comps()['Alias'], entry['attributes'])[1]
 #                   alias_num = alias0[alias0.find('_') + 1:]
#                    new_alias = ome + '_' + alias_num
 #                   entry['attributes'] = re.sub(
  #                      gff3Comps()['Alias'], 'Alias='+ new_alias,
   #                     entry['attributes']
    #                    )
                    entry['attributes'] = entry['attributes'].replace(old_ome, ome)
                new_gff = rename_and_organize(new_gff)
                gff = new_gff
            except:
                gff = curGFF3(gff, ome, cur_seqids = True)
#        else:
 #           gff = curGFF3(gff, ome)
        else:
            gff = curGFF3(gff, ome, cur_seqids = True)
    elif gffVer == 2.5:
        gff, trans_str, failed, flagged = gtf2gff3(gff, ome)
    else:
        gff, errors = gff2gff3(gff, ome, assembly_accession, verbose = False)

    ver_search = re.search(r'(.{6}\d+)\.(\d+)', ome)
    if ver_search is not None:
        less_ome = ver_search[1]
        ome_ver = ver_search[2]
    else:
        less_ome = ome
    for line in gff:
        seqid = line['seqid']
        if seqid.startswith(less_ome + '_'):
            line['seqid'] = re.sub(r'^' + less_ome + '_', 
                                 ome + '_', seqid)
        elif re.search(r'^' + less_ome + r'\.\d+_', seqid):
            line['seqid'] = re.sub(r'^' + less_ome + r'[^_]+_', 
                                 ome + '_', seqid)
        else:
            line['seqid'] = ome + '_' + seqid

    with open(cur_path + '.tmp', 'w') as out:
        out.write(list2gff(gff))
    shutil.move(cur_path + '.tmp', cur_path)

def add2failed(row):
    if isinstance(row['assembly_acc'], float) or not row['assembly_acc']:
        return False
    else:
        return [row['assembly_acc'], row['version']]

def main(
    predb, refdb, wrk_dir, 
    verbose = False, spacer = '\t\t\t', forbidden = set(), 
    cpus = 1, exit = False, remove = False
    ):

    infdb = predb2mtdb(predb)
    vprint('\nGenerating omes', v = verbose, flush = True)
    omedb, failed = gen_omes(infdb, refdb, ome_col = 'ome', forbidden = forbidden,
                     spacer = spacer)
    
    cur_cmds = []
    omedb = omedb.set_index('ome')
    for ome, row in omedb.items():
        cur_cmds.append([
            ome, row['fna'], row['gff3'], 
            wrk_dir, row['source'], row['assembly_acc'],
            exit, remove, spacer, verbose
            ])

    vprint('\nCurating data', v = verbose, flush = True)
    with mp.Pool(processes = cpus) as pool:
        cur_data = pool.starmap(cur_mngr, tqdm(cur_cmds, total = len(cur_cmds)))

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


def cli():
    usage = 'Generate a predb file:\npredb2mtdb.py\n\nCreate a mycotoolsdb ' + \
    'from a predb file:\npredb2mtdb.py <PREDBFILE>\n\nCreate a mycotoolsdb ' + \
    'referencing an alternative master database:\npredb2mtdb.py <PREDBFILE> ' + \
    '<REFERENCEDB>\nSkip failing genomes:\npredb2mtdb.py <PREDBFILE> -s'

    if any(x in {'-h', '--help', '-help'} for x in sys.argv):
        eprint('\n' + usage + '\n', flush = True)
        sys.exit(1)
    elif len(sys.argv) == 1:
        print(gen_predb())
        sys.exit(0)
    elif len(sys.argv) >= 3:
        if sys.argv[2] not in {'-s', '--skip'}:
            refDB = mtdb(format_path(sys.argv[2]))
        elif len(sys.argv) > 3:
            refDB = mtdb(format_path(sys.argv[3]))
        else:
            refDB = mtdb(primaryDB())
    else:
        refDB = mtdb(primaryDB())

    if set(sys.argv).intersection({'-s', '--skip'}):
        exit = False
    else:
        exit = True

    from Bio import Entrez
    ncbi_email, ncbi_api, jgi_email, jgi_pwd = loginCheck(jgi = False)
    Entrez.email = ncbi_email
    if ncbi_api:
        Entrez.api_key = ncbi_api

    eprint('\nPreparing run', flush = True)
    predb = read_predb(format_path(sys.argv[1]), spacer = '\t')
    out_dir, wrk_dir = prep_output(os.path.dirname(format_path(sys.argv[1])))

    omedb, failed = main(predb, refDB, wrk_dir, exit = exit, verbose = True)

    from mycotools.lib.dbtools import gather_taxonomy, assimilate_tax
    tax_dicts = gather_taxonomy(omedb, api_key = ncbi_api)
    outdb = assimilate_tax(omedb, tax_dicts)
    outdb.df2db(out_dir + 'predb2mtdb.mtdb')
    sys.exit(0)


if __name__ == '__main__':
    cli()
