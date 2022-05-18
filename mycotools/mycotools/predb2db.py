#!/usr/bin/env python3

import sys, os, re, shutil, multiprocessing as mp
from mycotools.lib.kontools import gunzip, mkOutput, formatPath, eprint, vprint
from mycotools.lib.biotools import gff2list, list2gff, fa2dict, dict2fa, \
    gff3Comps, gff2Comps
from mycotools.lib.dbtools import mtdb, masterDB, loginCheck
from mycotools.curAnnotation import main as curRogue
from mycotools.utils.curGFF3 import main as curGFF3
from mycotools.utils.gff2gff3 import main as gff2gff3
from mycotools.gff2seq import aamain as gff2seq

predbHeaders = [
    'genomeCode', 'genus', 'species', 'strain',
    'version', 'biosample',
    'assemblyPath', 'gffPath', 'genomeSource (ncbi/jgi/new)', 
    'useRestriction (yes/no)', 'publication'
    ]

def prepOutput(base_dir):
    out_dir = mkOutput(base_dir, 'predb2db')
    wrk_dir = out_dir + 'working/'
    dirs = [
        out_dir, wrk_dir, wrk_dir + 'gff/', wrk_dir + 'fna/', wrk_dir + 'faa/'
        ]
    for dir_ in dirs:
        if not os.path.isdir(dir_):
            os.mkdir(dir_)
    return dirs[:2]

def copyFile( old_path, new_path ):

    try:
        shutil.copy( old_path, new_path )
        return True
    except:
        raise IOError


def moveBioFile( old_path, ome, typ, wrk_dir, suffix = '' ):

    if old_path.endswith('.gz'):
        if not os.path.isfile(old_path[:-3]):
            temp_path = gunzip(old_path)
            if temp_path:
                new_path = wrk_dir + ome + '.' + typ +  suffix
            else:
                raise IOError
        else:
            new_path = wrk_dir + ome + '.' + typ + suffix
            temp_path = old_path[:-3]
        copyFile(formatPath(temp_path), new_path)
    else:
        new_path = wrk_dir + ome + '.' + typ +  suffix
        if not os.path.isfile(new_path) and os.path.isfile(old_path):
            copyFile(formatPath(old_path), new_path)
        elif not os.path.isfile(new_path):
            raise IOError
        else:
            copyFile(formatPath(old_path), new_path)

    return new_path


def genPredb():
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
    outputStr = '#' + '\t'.join(predbHeaders)
    outputStr += '\n#' + '\t'.join(example) + '\n'

    return outputStr

def readPredb(predb_path, spacer = '\t'):

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
        predb['genome_code'] = predb['genomeCode']
        del predb['genomeCode']
    except KeyError:
        if not 'genome_code' in predb and not 'genomeCode' in predb:
            predb['genome_code'] = ['' for x in predb['genus']]
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
        if not 'restriction' in predb and 'publication' not in predb:
            raise KeyError
        elif 'publication' in predb:
            predb['restriction'] = [bool(x) for x in predb['publication']]

    if any(x.lower() not in {'jgi', 'ncbi', 'new'} \
        for x in predb['source']):
        eprint(spacer + 'ERROR: genoemSource entries must be in {jgi, ncbi, new}', flush = True)
        sys.exit(5)
    for i, path in enumerate(predb['assemblyPath']):
        predb['assemblyPath'][i] = formatPath(predb['assemblyPath'][i])
        predb['gffPath'][i] = formatPath(predb['gffPath'][i])
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
        if predb['publication'][i]:
            toAdd = {
                'genome_code': predb['genome_code'][i],
                'genus': predb['genus'][i],
                'species': predb['species'][i],
                'strain': predb['strain'][i],
                'version': predb['version'][i],
                'biosample': predb['biosample'][i],
                'assembly': predb['assemblyPath'][i],
                'gff3': predb['gffPath'][i],
                'source': predb['source'][i],
                'published': predb['publication'][i]
                }
        else:
            toAdd = {
                'genome_code': predb['genome_code'][i],
                'genus': predb['genus'][i],
                'species': predb['species'][i],
                'strain': predb['strain'][i],
                'version': predb['version'][i],
                'biosample': predb['biosample'][i],
                'assembly': predb['assemblyPath'][i],
                'gff3': predb['gffPath'][i],
                'source': predb['source'][i],
                'published': not predb['restriction'][i]
                }
#        if predb['publication'][i]:
 #           toAdd['published'] = predb['publication'][i]
        infdb = infdb.append(toAdd)
    return infdb

def genOmes(newdb, refdb = None, ome_col = 'internal_ome'):

    tax_set = set(refdb['internal_ome'])
    for i, ome in enumerate(newdb['internal_ome']):
        if not ome:
            name = newdb['genus'][i][:3].lower() + newdb['species'][i][:3].lower()
            name = re.sub(r'\(|\)|\[|\]|\$|\#|\@| |\||\+|\=|\%|\^|\&|\*|\'|\"|\!|\~|\`|\,|\<|\>|\?|\;|\:|\\|\{|\}', '', name)
            name.replace('(', '')
            name.replace(')', '')
            while len(name) < 6:
                name += '.'
            numbers = [
                int(re.search(r'.*?(\d+$)', x)[1]) for x in list(tax_set) if x.startswith(name)
                ]
            numbers.sort(reverse = True)
            if numbers:
                number = numbers[0] + 1
            else:
                number = 1
            newOme = name + str(number)
            tax_set.add(newOme)
            newdb['internal_ome'][i] = newOme

    return newdb

def curMngr(ome, fna_path, gff_path, wrk_dir, source, genomeCode):

    try:
        newFNA_path = moveBioFile( fna_path, ome, 'fa', wrk_dir + 'fna/' )
    except IOError:
        return ome, False, 'assembly'
    try:
        newGFF_path = moveBioFile(gff_path, ome, 'gff3', wrk_dir + 'gff/', suffix = '.uncur')
    except IOError:
        return ome, False, 'gff3'

    gff = gff2list(newGFF_path)
    cur_path = re.sub(r'\.uncur$', '', newGFF_path)

    curGFF_path = gffMngr(ome, gff, cur_path, source, genomeCode)
    faa = gff2seq(gff, fa2dict(newFNA_path))
    faa_path = wrk_dir + 'faa/' + ome + '.aa.fa'
    with open(faa_path, 'w') as out:
        out.write(dict2fa(faa))

    return ome, newFNA_path, curGFF_path, faa_path


def gffMngr(ome, gff, cur_path, source, genomeCode):

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
                    aliasNum = alias0[alias0.find('_') + 1:]
                    newAlias = ome + '_' + aliasNum
                    entry['attributes'] = re.sub(
                        gff3Comps()['Alias'], 'Alias='+ newAlias,
                        entry['attributes']
                        )
            else:
                gff, trans_str, failed, flagged = curRogue(gff, ome)
        else:
            gff = curGFF3(gff, ome)
    else:
        gff = gff2gff3(gff, ome, genomeCode, verbose = False)

    with open(cur_path, 'w') as out:
        out.write(list2gff(gff))

    return cur_path 

def add2failed(row):
    if row['source'] == 'ncbi':
        return [row['biosample'], row['version']]
    else:
        return [row['genome_code'], row['version']]

def main(
    predb, refdb, wrk_dir, 
    verbose = False, spacer = '\t', cpus = 1
    ):

    failed = []

    infdb = predb2mtdb(predb)
    omedb = genOmes(infdb, refdb, ome_col = 'internal_ome')
    
    curCmds = []
    omedb = omedb.set_index('internal_ome')
    for ome, row in omedb.items():
        curCmds.append([
            ome, row['assembly'], row['gff3'], 
            wrk_dir, row['source'], row['genome_code']
            ])
    with mp.Pool(processes = cpus) as pool:
        curData = pool.starmap(curMngr, curCmds)

    for data in curData:
        if not data[1]:
            eprint(spacer + data[0] + ' failed ' + data[2], flush = True)
            failed.append(add2failed(omedb[data[0]]))
            del omedb[data[0]]
        else:
            ome, fnaPath, gffPath, faaPath = data
            omedb[ome]['assembly'] = fnaPath
            omedb[ome]['gff3'] = gffPath
            omedb[ome]['proteome'] = faaPath

    return omedb.reset_index(), failed


if __name__ == '__main__':
    usage = 'Generate a predb file:\npredb2db.py\n\nCreate a mycotoolsdb ' + \
    'from a predb file:\npredb2db.py <PREDBFILE>\n\nCreate a mycotoolsdb ' + \
    'referencing an alternative master database:\npredb2db.py <PREDBFILE> ' + \
    '<REFERENCEDB>'

    if any(x in {'-h', '--help', '-help'} for x in sys.argv):
        print('\n' + usage + '\n', flush = True)
        sys.exit(1)
    elif len(sys.argv) == 1:
        print(genPredb())
        sys.exit(0)
    elif len(sys.argv) == 3:
        refDB = mtdb(formatPath(sys.argv[2]))
    else:
        refDB = mtdb(masterDB())

    from Bio import Entrez
    ncbi_email, ncbi_api, jgi_email, jgi_pwd = loginCheck(jgi = False)
    Entrez.email = ncbi_email
    if ncbi_api:
        Entrez.api_key = ncbi_api

    predb = readPredb(formatPath(sys.argv[1]), spacer = '\t')
    out_dir, wrk_dir = prepOutput(os.path.dirname(formatPath(sys.argv[1])))
    omedb, failed = main(predb, refDB, wrk_dir)

    from mycotools.lib.dbtools import gather_taxonomy, assimilate_tax
    tax_dicts = gather_taxonomy(omedb, api_key = ncbi_api)
    outdb = assimilate_tax(omedb, tax_dicts)
    outdb.df2db(out_dir + 'predb2db.db')
    sys.exit(0)
