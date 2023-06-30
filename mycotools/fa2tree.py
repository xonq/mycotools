#! /usr/bin/env python3

# NEED to implement fa2clus?
# NEED to ignore non fasta inputs

import os
import re
import sys
import shutil
import argparse
import subprocess
from collections import defaultdict
from mycotools.lib.kontools import eprint, vprint, collect_files, \
    format_path, intro, outro, findExecs, mkOutput, multisub
from mycotools.lib.biotools import fa2dict, dict2fa

class PhyloError(Exception):
    pass

def mafftRun(name, fasta, out_dir, hpc, verbose = True, cpus = 1, spacer = '\t'):

    cmd = [
        'mafft', '--auto', '--thread', '-' + str(cpus), fasta
        ]

    if not hpc:
        print(spacer + 'Aligning', flush = True)
        with open(out_dir + name + '.mafft', 'w') as out_file:
            if verbose:
                run_mafft = subprocess.call(cmd, stdout = out_file)
            else:
                run_mafft = subprocess.call(cmd, stdout = out_file, stderr = subprocess.DEVNULL)
         
        if run_mafft != 0:
            eprint(spacer + '\tERROR: mafft failed: ' + str(run_mafft), flush = True)
            if os.path.isfile(out_dir + name + '.mafft'):
                os.remove(out_dir + name + '.mafft')	
            raise PhyloError
    else:
        with open(out_dir + name + '_mafft.sh', 'w') as out:
            if '#PBS' in hpc:
                out.write( hpc + '\n\n' + ' '.join([str(x) for x in cmd]) + \
                ' > ' + out_dir + '/' + name + '.mafft' + \
                '\ncd ' + out_dir + '\nqsub ' + name + '_clipkit.sh')
            else:
                 out.write( hpc + '\n\n' + ' '.join([str(x) for x in cmd]) + \
		' > ' + out_dir + '/' + name + '.mafft' + \
                '\ncd ' + out_dir + '\nsbatch ' + name + '_clipkit.sh')
        
    return out_dir + '/' + name + '.mafft'


def trimRun(name, mafft, out_dir, hpc, param, 
            verbose, cpus = 1, spacer = '\t'):

    name2 = out_dir + os.path.basename(mafft) + '.clipkit'
    cmd = ['clipkit', mafft, '--output', name2]
    cmd.extend(param)
    if not hpc:
        print(spacer + 'Trimming', flush = True)
        if verbose:
            run_clipkit = subprocess.call(cmd, stdout = subprocess.PIPE)
        else:
            run_clipkit = subprocess.call( 
                cmd, stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL
                )
        if run_clipkit != 0:
            eprint(spacer + '\tERROR: `clipkit` failed: ' + str(run_clipkit) , flush = True)  
            if os.path.isfile(out_dir + name2):
                os.remove(out_dir + name2)
            raise PhyloError

    else:
        with open(out_dir + name + '_clipkit.sh', 'w') as out:
            if '#PBS' in hpc:
                out.write( hpc + '\n\n' + #cmd1 + '\n' + 
                    ' '.join([str(x) for x in cmd]) + '\n\ncd ' + out_dir \
                    + '\nqsub ' + name + '_tree.sh' )
            else:
                out.write( hpc + '\n\n' + #cmd1 + '\n' + 
                    ' '.join([str(x) for x in cmd]) + '\n\ncd ' + out_dir \
                    + '\nsbatch ' + name + '_tree.sh' )
    
    return name2


def run_mf(clipkit_files, out_dir, constraint = False, 
           verbose = False, cpus = 1, spacer = '\t', scripts = False):
    cpus_per_cmd = 3
    concurrent_cmds = round((cpus - 1)/4 - 0.5) # round down
    cmds = [['iqtree', '-m', 'MFP+MERGE', '-nt', 'AUTO',
             '-B', '1000', '-s', f_, '--prefix', 
             out_dir + os.path.basename(f_)] \
            for f_ in clipkit_files]
    if constraint:
        [x.extend(['-g', constraint]) for x in cmds]
    if scripts:
        return {os.path.basename(v): cmds[i] for i, v in enumerate(clipkit_files)}
    multisub(cmds, verbose = verbose, processes = concurrent_cmds)

    models = {}
    for f_ in clipkit_files:
        with open(out_dir + os.path.basename(f_) + '.log', 'r') as raw:
            for line in raw:
                if line.startswith('Best-fit model:'):
                    model = re.search(r'Best-fit model: (.+) chosen', line)[1]
                    models[f_] = model
                    break

    return models

def prepare_nexus(concat_fa, models, spacer = '\t'):

    nex_data = '#nexus\nbegin sets;'
    index = 1
    for i, trimmed_f in enumerate(list(models.keys())):
        trim_fa = fa2dict(trimmed_f)
        len0 = len(trim_fa[list(trim_fa.keys())[0]]['sequence'])
        if not all(len(x['sequence']) == len0 for x in trim_fa.values()):
            eprint(spacer + '\tERROR: alignment sequences are not same length ' \
                   + trimmed_f, flush = True)
            sys.exit(17)
        new_index = index + len0
        nex_data += '\n    charset part' + str(i + 1) + ' = ' \
                  + str(index) + '-' + str(new_index - 1) + ';'
        index = new_index

        for seq, data in trim_fa.items():
            if '_' in seq:
                ome = seq[:seq.find('_')]
            else:
                ome = seq
            if ome in concat_fa:
                concat_fa[ome]['sequence'] += data['sequence']

    nex_data += '\n    charpartition mine = '
    for i, model in enumerate(models.values()):
        nex_data += model + ':part' + str(i+1) + ','
    nex_data = nex_data[:-1] + ';\nend;'

    return nex_data, concat_fa 

def run_partition_tree(fa_file, nex_file, constraint,
                       verbose, cpus, spacer = '\t'):
    cmd = ['iqtree', '-s', fa_file, '-p', nex_file, '-nt', str(cpus),
           '-B', '1000', '--sampling', 'GENESITE'] # this is not using the
           # specified CPUs, but instead the most efficient
    if constraint:
        cmd.extend(['-g', constraint])
    print(spacer + 'Tree building', flush = True)
    if verbose:
        run_tree = subprocess.call(cmd)
    else:
        run_tree = subprocess.call(cmd, stdout = subprocess.DEVNULL,
                                   stderr = subprocess.DEVNULL)
    return run_tree



def treeRun(name, clipkit_file, out_dir, hpc, constraint,
            verbose, fast = False, cpus = 1, spacer = '\t'):

    tree_file = out_dir + os.path.basename(clipkit_file)
    if fast:
        cmd = ['fasttree', '-out', tree_file + '.treefile', clipkit_file]
        if constraint:
            cmd.extend(['-constraints', constraint])
    else:
        # the following will identify the optimum number of threads
        cmd = ['iqtree', '-s', clipkit_file, '-B', '1000', '-nt', str(cpus),
               '--prefix', tree_file]
        if constraint:
            cmd.extend(['-g', constraint])
    if not hpc:
        print(spacer + 'Tree building', flush = True)
        if fast:
            cmd[2] += '.tmp'
        if verbose:
            run_tree = subprocess.call(cmd)
        else:
            run_tree = subprocess.call(
                cmd, stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL
                )
        if fast:
            shutil.move(cmd[2], cmd[2][:-4])
        if run_tree != 0:
            eprint(spacer + '\tERROR: tree failed: ' + str(run_tree), flush = True)
            raise PhyloError
    else: 
        vprint('\nOutputting bash script `tree.sh`.\n', v = verbose, flush = True)
        with open(out_dir + name + '_tree.sh', 'w') as out:
            out.write(hpc + '\n\n' + ' '.join([str(x) for x in cmd]))


def main( 
    fasta_path, slurm = False, torque = False, fast = False, project = '', 
    output_dir = None, verbose = True, alignment = False, constraint = False,
    ome_conv = False, mem = '60GB', cpus = 1, spacer = '\t\t', align_stop = False, 
    tree_stop = False, flag_incomplete = True, start = 0, partition = False,
    clip_list = []
    ):


    hpc, hpc_prep = False, False
    if slurm:
        vprint('\nHPC mode, preparing submission scripts (Slurm)\n', v = verbose, flush = True)
        hpc_prep = '#!/bin/bash\n#SBATCH --time=24:00:00\n#SBATCH --nodes=1\n' + \
            '#SBATCH --ntasks-per-node=' + str(cpus) + '\n#SBATCH -A ' + project + '\n' + \
            '#SBATCH --mem=' + str(mem)
#            '\n\nsource activate ' + source
    elif torque:
        vprint('\nHPC mode, preparing submission scripts (PBS)\n', v = verbose, flush = True)
        hpc_prep = '#PBS -l walltime=10:00:00\n#PBS -l nodes=1:ppn=4\n#PBS -A ' + project
    if not partition and hpc_prep:
        hpc = hpc_prep

    output_dir_prep = os.getcwd() + '/'
    if isinstance(fasta_path, str):
        if os.path.isfile(fasta_path):
            if not output_dir:
                dir_name = re.sub(r'\..*?$', '_tree', os.path.basename( os.path.abspath(fasta_path)))
                out_dir = output_dir + '/' + dir_name
                if not os.path.isdir(out_dir):
                    os.mkdir(out_dir)
                out_dir = format_path(out_dir)
            else:
                if not os.path.isdir(output_dir):
                    os.mkdir(output_dir)
                out_dir = format_path(output_dir)

            wrk_dir = out_dir + 'working/'
            if not os.path.isdir(wrk_dir):
                os.mkdir(wrk_dir)

            files = [fasta_path]
        elif os.path.isdir(fasta_path):
            if not output_dir:
                out_dir = mkOutput(output_dir_prep, 'fa2tree')
            else:
                if not os.path.isdir(output_dir):
                    os.mkdir(output_dir)
                out_dir = format_path(output_dir)

            wrk_dir = out_dir + 'working/'
            if not os.path.isdir(wrk_dir):
                os.mkdir(wrk_dir)

            files = collect_files(fasta_path, '*')
            check_fas = set(files)
            new_fas = []
            for f in files:
                if f + '.mafft.clipkit' not in check_fas:
                    if f + '.mafft' not in check_fas:
                        new_fas.append(f)
                    else:
                        shutil.copy(f + '.mafft', wrk_dir + os.path.basename(f) \
                                                          + '.mafft')
                else:
                    shutil.copy(f + '.mafft.clipkit',
                                wrk_dir + os.path.basename(f) + '.mafft.clipkit')
            files = new_fas
    else:
        if not output_dir:
            out_dir = mkOutput(output_dir_prep, 'fa2tree')
        else:
            if not os.path.isdir(output_dir):
                os.mkdir(output_dir)
            out_dir = format_path(output_dir)
        wrk_dir = out_dir + 'working/'
        if not os.path.isdir(wrk_dir):
            os.mkdir(wrk_dir)
        files = []
        for path in fasta_path:
            if os.path.isdir(path):
                files.extend([format_path(x) for x in collect_files(path, '*')])
            else:
                files.append(path)

    conv_dir = wrk_dir + 'conv/'
    if not os.path.isdir(conv_dir):
        os.mkdir(conv_dir)

    # check for non fasta inputs
    if any(f_.endswith(('.nexus', '.nex', '.nxs', 
                        '.phylip', '.phy', '.ph',
                        '.clustal', '.clus',)) \
           for f_ in files):
        print(spacer + 'Converting to fastas', flush = True)
        from Bio import SeqIO
        for i, f_ in enumerate(files):
            out_name = conv_dir + os.path.basename(
                re.search(r'(.*)\.[^.]+$', f_)[1] + '.fa'
                )
            if not os.path.isfile(out_name):
                if f_.endswith(('.nexus', '.nex', '.nxs',)):
                    records = SeqIO.parse(f_, 'nexus')
                    SeqIO.write(records, out_name, 'fasta')
                elif f_.endswith(('.phylip', '.phy', '.ph',)):
                    records = SeqIO.parse(f_, 'phylip')
                    SeqIO.write(records, out_name, 'fasta')
                elif f_.endswith(('.clustal', '.clus',)):
                    records = SeqIO.parse(f_, 'clustal')
                    SeqIO.write(records, out_name, 'fasta')
                else:
                    continue
            files[i] = out_name

    files = [x for x in files \
             if x.endswith(('fasta', 'fa', 'fna', 'faa', 'fsa'))]

    if partition: #multigene partition mode
        ome2fa2gene = defaultdict(dict)
        incomp_omes = {}
        for f in files:
            if f.endswith('/'):
                f = f[:-1]
            fa_dict = fa2dict(f)
            incomp_omes[os.path.basename(f)] = []
            for seq in fa_dict:
                ome = seq[:seq.find('_')]
                if f not in ome2fa2gene[ome]:
                    incomp_omes[os.path.basename(f)].append(ome)
                    ome2fa2gene[ome][f] = seq
                else:
                    eprint('\nERROR: multiple sequences for a ' \
                         + 'single ome in ' + f, flush = True)
                    sys.exit(5)

        incomp_files, comp_files = {}, []
        for f, omes in incomp_omes.items():
            if len(omes) < len(ome2fa2gene):
                incomp_files[f] = sorted(
                    set(ome2fa2gene.keys()).difference(set(omes))
                    )
            else:
                comp_files.append(f)

        if incomp_files:
            incomp_files = {k: v for k,v in sorted(incomp_files.items(), 
                                                   key = lambda x: len(x[1]))}
            if flag_incomplete:
                eprint('\nERROR: omes without sequences in all fastas: ', flush = True)
            else:
                eprint('\nWARNING: omes removed without sequences in all \
                        fastas: ', flush = True)
            for f, omes in incomp_files.items():
                eprint('\t' + f + ': ' + ','.join([x for x in omes]), flush = True)
            if comp_files:
                eprint('\tComplete files: ' + ','.join(comp_files), flush = True)
            else:
                eprint('\tNo complete files', flush = True)
            if flag_incomplete:
                sys.exit(6)

    trimmed_files = []
    for fasta in files:
        name = re.search(r'.*?/*([^/]+)/*$', fasta)[1]
        fasta_name = os.path.basename(os.path.abspath(fasta))
        clipkit = wrk_dir + fasta_name + '.clipkit'
        if start == 0:
            mafft = wrk_dir + fasta_name + '.mafft'
            clipkit = mafft + '.clipkit'
            try:
                if os.stat(mafft):
                    vprint('\nAlignment exists', v = verbose, flush = True)
                else:
                    raise ValueError
            except (FileNotFoundError, ValueError) as e:
                if not alignment:
                    mafft = mafftRun(name, os.path.abspath(fasta), wrk_dir, hpc, 
                                     verbose, cpus, spacer = spacer)
                else:
                    mafft = fasta
        else:
            mafft = fasta
            clipkit = mafft + '.clipkit'


        if start == 0 or start == 1: 
            try:
                if os.stat(clipkit):
                    vprint('\nTrim exists', v = verbose, flush = True)
                    clipkit_out = clipkit
                else:
                    raise ValueError
            except (FileNotFoundError, ValueError) as e:
                clipkit_out = trimRun(name, mafft, wrk_dir, hpc, clip_list,
                                      verbose, spacer = spacer)
            trimmed_files.append(clipkit_out)
        else:
            trimmed_files.append(fasta)


    if ome_conv:
        new_trimmed_files = []
        for trimmed_f in trimmed_files:
            new_f = conv_dir + os.path.basename(trimmed_f)
            in_fa = fa2dict(trimmed_f)
            out_fa = {}
            for seq, data in in_fa.items():
                ome = seq[:seq.find('_')]
                out_fa[ome] = data
            with open(new_f, 'w') as out:
                out.write(dict2fa(out_fa))
            new_trimmed_files.append(new_f)
        trimmed_files = new_trimmed_files



    if not partition:
        for trimmed_f in trimmed_files:
            name_prep = re.search(r'.*?/*([^/]+)/*$', trimmed_f)[1]
            name = re.sub(r'\.mafft\.clipkit$', '', name_prep)
            treeRun(name, trimmed_f, out_dir, hpc, constraint,
                    verbose, fast, cpus, spacer = spacer)
    else:
        if align_stop:
            script_dict = run_mf(trimmed_files, wrk_dir, constraint,
                            verbose, cpus, spacer = spacer + '\t',
                            scripts = True)
            for f, cmd in script_dict.items():
                with open(out_dir + f + '.sh', 'w') as out:
                    if hpc_prep:
                        out.write(f'{hpc_prep}\n#SBATCH --output={f}_tree\n' \
                                + f'#SBATCH --error={f}_tree\n' \
                                + f'#SBATCH --job-name={f}_tree\n')
                    out.write(' '.join(cmd))
            print(spacer + 'Alignments outputed', flush = True)
            sys.exit(0)

        print(spacer + 'Model finding', flush = True)
        models = run_mf(trimmed_files, wrk_dir, constraint,
                        verbose, cpus, spacer = spacer + '\t')

        print(spacer + 'Concatenating', flush = True)
        empty_concat_fa = {ome: {'description': '', 'sequence': ''}
                     for ome in ome2fa2gene if ome not in incomp_omes}
        nex_data, concat_fa = prepare_nexus(empty_concat_fa, models)
        with open(out_dir + 'concatenated.fa', 'w') as out:
            out.write(dict2fa(concat_fa))
        with open(out_dir + 'concatenated.nex', 'w') as out:
            out.write(nex_data)

        if tree_stop:
            print(spacer + 'Concatenated nexus and fasta outputted', 
                  flush = True)
            sys.exit(0)

        run_partition_tree(out_dir + 'concatenated.fa', 
                           out_dir + 'concatenated.nex', constraint,
                           verbose, cpus, spacer = spacer)

    if hpc:
        if slurm:
            print('\nStart pipeline via `sbatch <STARTSTEP>.sh` in ' \
                + out_dir + '\n')
        else:
            vprint('\nStart pipeline via `qsub <STARTSTEP>.sh` in ' + out_dir + '\n',
                v = verbose)


def cli():

    parser = argparse.ArgumentParser( 
        description = 'Takes in a multifasta, or directory/list of fastas [partition]  ' + \
        ' aligns (mafft), trims (clipkit), and infers phylogeny \
        (iqtree/fastree). When using non-Mycotools fastas for multigene \
        partition analysis, all accessions must be formatted as \
        "<GENOME>_<ACCESSION>" and each <GENOME> needs to be in ALL fastas'
        )

    io_opt = parser.add_argument_group('Input options')
    io_opt.add_argument('-i', '--input', required = True,
        help = 'Fasta or directory of fastas/alignments (fasta, phylip, \
                nexus, or clustal). "-" for stdin')
    io_opt.add_argument('-a', '--align', help = 'Start from alignments',
                        action = 'store_true')
    io_opt.add_argument('-t', '--trim', help = 'Start from trimmed alignments',
                        action = 'store_true')

    gt_opt = parser.add_argument_group('General options')
    gt_opt.add_argument('-g', '--gappy', 
        help = 'ClipKIT gappy threshold')
    gt_opt.add_argument('--constraint', help = 'Topological constraint')
    gt_opt.add_argument('--ome', help = 'Convert input sequence to mtdb omes',
                        action = 'store_true')

    s_opt = parser.add_argument_group('Single tree options')
    s_opt.add_argument('-f', '--fast', action = 'store_true', \
        help = 'Fasttree' )
    s_opt.add_argument('-s', '--slurm', action = 'store_true', \
        help = 'Submit ALL steps to Slurm' )
    s_opt.add_argument('--torque', action = 'store_true', \
        help = 'Submit ALL steps to Torque.' )
    s_opt.add_argument('-A', '--project', help = 'HPC project')

    m_opt = parser.add_argument_group('Partition tree options')
    m_opt.add_argument('-p', '--partition', action = 'store_true',
        help = 'Multigene partition mode')
    m_opt.add_argument('-sa', '--stop_align', action = 'store_true',
        help = 'Stop after aligning and trimming, output tree scripts')
    m_opt.add_argument('-sc', '--stop_cat', action = 'store_true',
        help = 'Stop after concatenating sequences')
    m_opt.add_argument('-m', '--missing', action = 'store_true',
        help = 'Remove omes with missing sequences')

    r_opt = parser.add_argument_group('Runtime options')
    r_opt.add_argument('-v', '--verbose', action = 'store_true')
    r_opt.add_argument('-o', '--output')
    r_opt.add_argument('-c', '--cpus', default = 1, type = int)

    args = parser.parse_args()

    execs = ['mafft', 'clipkit']
    if args.fast:
        execs.append('fasttree')
    else:
        execs.append('iqtree')
    if args.torque or args.slurm:
        findExecs(execs, execs)
    else:
        findExecs(execs, execs)

    output = format_path(args.output)
    if output:
        if not output.endswith('/'):
            output += '/' # bring to mycotools path expectations
    if args.gappy:
        clip_list = ['-m', 'gappy', '-g', args.gappy]
    else:
        clip_list = []


    args_dict = {
        'Fasta': args.input, 'Fast': args.fast, 
        'Multigene partition': args.partition, 'Input alignments': bool(args.align),
        'Input trimmed': bool(args.trim), 'Constraint': args.constraint,
        'Ome conversion': args.ome, 'ClipKIT parameters': ' '.join(clip_list),
        'Alignment stop': args.stop_align, 'Concatenation stop': args.stop_cat,
        'Ignore missing': args.missing, 'Torque': args.torque, 'Slurm': args.slurm, 
        'HPC project': args.project, 'Output': output, 'CPUs': args.cpus
        }
    start_time = intro('Fasta2Tree', args_dict)


    if args.input == '-':
        args.input = \
            stdin2str().replace('"','').replace("'",'').split()
    else:
        args.input = args.input.replace('"','').replace("'",'').split()


    start = 0
    if args.trim:
        start = 2
    elif args.align:
        start = 1
    main( 
        args.input, slurm = args.slurm, fast = args.fast,
        torque = args.torque, project = args.project,
        output_dir = output, verbose = bool(args.verbose), 
        alignment = args.align, constraint = format_path(args.constraint),
        ome_conv = args.ome, clip_list = clip_list,
        align_stop = args.stop_align, tree_stop = args.stop_cat,
        flag_incomplete = bool(not args.missing), start = start,
        partition = bool(args.partition), cpus = args.cpus
        )

    outro(start_time)

if __name__ == '__main__':
    cli()
