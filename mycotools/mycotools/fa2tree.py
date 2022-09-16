#! /usr/bin/env python3

# NEED to implement multigene phylogeny
# NEED to implement fa2clus?

import os
import re
import sys
import argparse
import subprocess
from collections import defaultdict
from mycotools.lib.kontools import eprint, vprint, collect_files, \
    format_path, intro, outro, findExecs
from mycotools.lib.biotools import fa2dict

class PhyloError(Exception):
    pass

def mafftRun(fasta, out_dir, hpc, verbose = True, cpus = 1, spacer = '\t'):

    name = fasta
    if '/' in fasta:
        name = os.path.basename(os.path.abspath(fasta)) 
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
        with open(out_dir + '/mafft.sh', 'w') as out:
            if '#PBS' in hpc:
                out.write( hpc + '\n\n' + ' '.join([str(x) for x in cmd]) + \
                ' > ' + out_dir + '/' + name + '.mafft' + \
                '\ncd ' + out_dir + '\nqsub clipkit.sh')
            else:
                 out.write( hpc + '\n\n' + ' '.join([str(x) for x in cmd]) + \
		' > ' + out_dir + '/' + name + '.mafft' + \
                '\ncd ' + out_dir + '\nsbatch clipkit.sh')
        
    return out_dir + '/' + name + '.mafft'


def trimRun(mafft, out_dir, hpc, verbose, cpus = 1, spacer = '\t'):

    name2 = os.path.basename(os.path.abspath(re.sub( r'\.mafft$', '.clipkit', mafft)))
    cmd = ['clipkit', mafft, '--output', out_dir + '/' + name2]
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
        with open(out_dir + '/clipkit.sh', 'w') as out:
            if '#PBS' in hpc:
                out.write( hpc + '\n\n' + #cmd1 + '\n' + 
                    ' '.join([str(x) for x in cmd]) + '\n\ncd ' + out_dir + '\nqsub tree.sh' )
            else:
                out.write( hpc + '\n\n' + #cmd1 + '\n' + 
                    ' '.join([str(x) for x in cmd]) + '\n\ncd ' + out_dir + '\nsbatch tree.sh' )
    
    return out_dir + '/' + name2


def run_mf(clipkit_file, out_dir, hpc, verbose, cpus = 1, spacer = '\t'):
    cmd = ['iqtree', '-m', 'MFP+MERGE', '-T', '1', '-B', '1000', '-s',
           clipkit_file]
    vprint(spacer + clipkit_file, v = verbose, flush = True)
    if verbose:
        run_mf = subprocess.call(cmd)
    else:
        run_mf = subprocess.call(cmd, stdout = subprocess.DEVNULL,
                                 stderr = subprocess.DEVNULL)
    with open(clipkit_file + '.log', 'r') as raw:
        for line in raw:
            if line.startswith('Best-fit model:'):
                model = re.search(r'Best-fit model: ([^\W]+)', line)[1]
                break

    return model

def prepare_nexus(concat_fa, models):

    nex_data = '#nexus\nbegin sets;'
    index = 1
    for i, trimmed_f in enumerate(list(models.keys())):
        trim_fa = fa2dict(trimmed_f)
        len0 = trim_fa[list(trim_fa.keys())[0]]['sequence']
        if not all(len(x['sequence']) == len0):
            eprint(spacer + 'ERROR: alignment sequences are not same length ' \
                   + trimmed_f, flush = True)
            sys.exit(17)
        new_index = index + len0
        nex_data += '\n    charset part' + str(i + 1) + ' = ' \
                  + str(index) + '-' + str(new_index) + ';'
        index = newindex

        for seq, data in trim_fa.items():
            ome = seq[:seq.find('_')]
            if ome in concat_fa:
                concatfa_data[ome]['sequence'] += data['sequence']

    nex_data += '\n    charpartition mine = '
    for i, model in models.values():
        nex_data += model + ':part' + str(i +1) + ','
    nex_data = nex_data[:-1] + ';\nend;'

    return nex_data, concat_fa 

def run_partition_tree(fa_file, nex_file, verbose, cpus, spacer = '\t'):
    cmd = ['iqtree', '-s', fa_file, '-p', nex_file, '-T', str(cpus),
           '-B', '1000', '--sampling', 'GENESITE']
    print(spacer + 'Tree building', flush = True)
    if verbose:
        run_tree = subprocess.call(cmd)
    else:
        run_tree = subprocess.call(cmd, stdout = subprocess.DEVNULL,
                                   stderr = subprocess.DEVNULL)
    return run_tree



def treeRun(clipkit_file, out_dir, hpc, verbose, fast = False, cpus = 1, spacer = '\t'):

    if fast:
        cmd = ['fasttree', '-out', clipkit_file + '.treefile', clipkit_file]
    else:
        cmd = ['iqtree', '-s', clipkit_file, '-B', '1000', '-T', str(cpus)]

   
    if not hpc:
        print(spacer + 'Tree building', flush = True)
        if verbose:
            run_tree = subprocess.call(
                cmd
                )
        else:
            run_tree = subprocess.call(
                cmd, stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL
                )
        if run_tree != 0:
            eprint(spacer + '\tERROR: tree failed: ' + str(run_tree), flush = True)
            raise PhyloError
    else: 
        vprint('\nOutputting bash script `tree.sh`.\n', v = verbose, flush = True)

        with open(out_dir + '/tree.sh', 'w') as out:
            out.write(hpc + '\n\n' + ' '.join([str(x) for x in cmd]))


def main( 
    fasta_path, slurm = False, torque = False, fast = False, project = '', 
    output_dir = None, verbose = True, alignment = False, mem = '60GB', cpus = 1,
    spacer = '\t\t', align_stop = False, tree_stop = False, 
    flag_incomplete = True
    ):

    hpc = False
    if slurm:
        vprint('\nHPC mode, preparing submission scripts (Slurm)\n', v = verbose, flush = True)
        hpc = '#!/bin/bash\n#SBATCH --time=24:00:00\n#SBATCH --nodes=1\n' + \
            '#SBATCH --ntasks-per-node=' + str(cpus) + '\n#SBATCH -A ' + project + '\n' + \
            '#SBATCH --mem=' + str(mem)
#            '\n\nsource activate ' + source
    elif torque:
        vprint('\nHPC mode, preparing submission scripts (PBS)\n', v = verbose, flush = True)
        hpc = '#PBS -l walltime=10:00:00\n#PBS -l nodes=1:ppn=4\n#PBS -A ' + project

    partition = False
    output_dir_prep = os.getcwd() + '/'
    if len(fasta_path) == 1:
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
            fastas = [fasta_path]
        elif os.path.isdir(fasta_path):
            vprint('\nDirectory input, partition analysis mode', 
                   v = verbose, flush = True)
            partition = True
            if not output_dir:
                out_dir = mkOutput(output_dir_prep, 'fa2tree')
            else:
                if not os.path.isdir(output_dir):
                    os.mkdir(output_dir)
                out_dir = format_path(output_dir)
            fastas = collect_files(fasta_path, 'fa')
            fastas.extend(collect_files(fasta_path, 'faa'))
            fastas.extend(collect_files(fasta_path, 'fna'))
            fastas.extend(collect_files(fasta_path, 'fasta'))
            fastas.extend(collect_files(fasta_path, 'fsa'))
            if prealigned:
                fastas.extend(collect_files(fasta_path, 'mafft'))
                fastas.extend(collect_files(fasta_path, 'clipkit'))
                check_fas = set(fastas)
                new_fas = []
                for fa in fastas:
                    if fa + '.clipkit' not in check_fas:
                        if fa + '.mafft' not in check_fas:
                            new_fas.append(fa)
                fastas = new_fas
    elif len(fasta_path) > 1:
        vprint('\nMultiple fastas inputted, partition analysis mode',
               v = verbose, flush = True)
        partition = True
        if not output_dir:
            out_dir = mkOutput(output_dir_prep, 'fa2tree')
        else:
            if not os.path.isdir(output_dir):
                os.mkdir(output_dir)
            out_dir = format_path(output_dir)
        fastas = [format_path(x) for x in fasta_path]

    if partition: #multigene partition mode
        ome2fa2gene = defaultdict(dict)
        complete_omes = []
        for fasta in fastas:
            fa_dict = fa2dict(fasta)
            for seq in fa_dict:
                ome = seq[:seq.find('_')]
                if fasta not in ome2fa2gene[ome]:
                    ome2fa2gene[ome][fasta] = seq
                else:
                    eprint('\nERROR: multiple sequences for a \
                            single ome in ' + fasta, flush = True)
                    sys.exit(5)
            incomp_omes = set([x for x, fas in ome2fa2gene.items() \
                                 if len(fas) != len(fastas)])
            if incomp_omes:
                if flag_incomplete:
                    eprint('\nERROR: omes without sequences in all fastas: ' \
                         + ','.join([str(x) for x in list(incomp_omes)]))
                    sys.exit(6)
                else:
                    eprint('\nWARNING: omes removed without sequences in all \
                            fastas: ' \
                         + ','.join([str(x) for x in list(incomp_omes)]))

    trimmed_files = []
    for fasta in fastas:
        if not prealigned:
            fasta_name = os.path.basename(os.path.abspath(fasta))
            mafft = out_dir + '/' + fasta_name + '.mafft'
            clipkit = out_dir + '/' + fasta_name + '.clipkit'
            
            try:
                if os.stat(mafft):
                    vprint('\nAlignment exists', v = verbose, flush = True)
                else:
                    raise ValueError
            except (FileNotFoundError, ValueError) as e:
                if not alignment:
                    mafft = mafftRun(os.path.abspath(fasta_path), out_dir, hpc, 
                                     verbose, cpus, spacer = spacer)
                else:
                    mafft = fasta_path
        
            try:
                if os.stat(clipkit):
                    vprint('\nTrim exists', v = verbose, flush = True)
                    clipkit_out = clipkit
                else:
                    raise ValueError
            except (FileNotFoundError, ValueError) as e:
                clipkit_out = trimRun(mafft, out_dir, hpc, verbose, spacer = spacer)
            trimmed_files.append(clipkit_out)

    if align_stop:
        vprint(spacer + 'Alignments outputed', v = verbose, flush = True)
        sys.exit(0)

    if not partition:
        for trimmed_f in trimmed_files:
            treeRun(trimmed_f, out_dir, hpc, verbose, fast, cpus, spacer = spacer)
    else:
        print(spacer + 'Model finding', flush = True)
        models = {}
        for trimmed_f in trimmed_files:
            model = run_mf(trimmed_f, out_dir, hpc, verbose, cpus,
                   spacer = spacer + '\t')
            models[trimmed_f] = model

        print(spacer + 'Concatenating', flush = True)
        concat_fa = {ome: {'description': '', 'sequence': ''}
                     for ome in ome2fa2gene if ome not in incomp_omes}
        nex_data, concat_fa = prepare_nexus(concat_fa, models)
        with open(out_dir + 'concatenated.fa', 'w') as out:
            out.write(dict2fa(concat_fa))
        with open(out_dir + 'concatenated.nex', 'w') as out:
            out.write(nex_data)

        if tree_stop:
            vprint(spacer + 'Concatenated nexus and fasta outputted', 
                    v = verbose, flush = True)
            sys.exit(0)

        run_partition_tree(out_dir + 'concatenated.fa', 
                           out_dir + 'concatenated.nex',
                           verbose, cpus, spacer = spacer)

    if hpc:
        if slurm:
            vprint('\nStart pipeline via `sbatch <STARTSTEP>.sh` in ' + out_dir + '\n',
                v = verbose)
        else:
            vprint('\nStart pipeline via `qsub <STARTSTEP>.sh` in ' + out_dir + '\n',
                v = verbose)


if __name__ == '__main__':

    parser = argparse.ArgumentParser( 
        description = 'Takes in a multifasta, or directory/list of fastas [partition]  ' + \
        ' aligns (mafft), trims (clipkit), and infers phylogeny \
        (iqtree/fastree). When using non-Mycotools fastas for multigene \
        partition analysis, all accessions must be formatted as \
        "<GENOME>_<ACCESSION>" and each <GENOME> needs to be in ALL fastas'
        )

    io_opt = parser.add_argument_group('Input/Output')
    io_opt.add_argument('-i', '--input', required = True,
        help = 'Fasta or directory of fastas/alignments [partition]. \
                "-" for stdin')
    io_opt.add_argument('-o', '--output')

    s_opt = parser.add_argument_group('Single tree options')
    s_opt.add_argument('--fast', action = 'store_true', \
        help = 'Fasttree' )
    s_opt.add_argument('-s', '--slurm', action = 'store_true', \
        help = 'Submit ALL steps to Slurm' )
    s_opt.add_argument('-t', '--torque', action = 'store_true', \
        help = 'Submit ALL steps to Torque.' )
    s_opt.add_argument('-A', '--project', help = 'HPC project')

    m_opt = parser.add_argument_group('Partition tree options')
    m_opt.add_argument('-a', '--stop_align', action = 'store_true',
        help = 'Stop after aligning and trimming')
    m_opt.add_argument('-p', '--prealigned', action = 'store_true',
        help = 'Use alignments')
    m_opt.add_argument('-c', '--stop_cat', action = 'store_true',
        help = 'Stop after concatenating sequences')
    m_opt.add_argument('-m', '--missing', action = 'store_true',
        help = 'Remove omes with missing sequences')
    args = parser.parse_args()

    output = format_path(args.output)
    if output:
        if not output.endswith('/'):
            output += '/' # bring to mycotools path expectations
    args_dict = {
        'Fasta': args.input, 'Fast': args.fast, 'Prealigned': args.prealigned,
        'Alignment stop': args.stop_align, 'Concatenation stop': args.stop_cat,
        'Ignore missing': args.missing,
        'Torque': args.torque, 'Slurm': args.slurm, 'HPC project': args.project, 'Output': output
        }
    start_time = intro('Fasta2Tree', args_dict)

    execs = ['mafft', 'clipkit']
    if args.fast:
        execs.append('fasttree')
    else:
        execs.append('iqtree')
    if args.torque or args.slurm:
        findExecs(execs, execs)
    else:
        findExecs(execs, execs)

    if args.input == '-':
        args.input = \
            stdin2str().replace('"','').replace("'",'').split()
    else:
        args.input = args.input.replace('"','').replace("'",'').split()
    main( 
        args.input, slurm = args.slurm, fast = args.fast,
        torque = args.torque, project = args.project,
        output_dir = output, verbose = True, alignment = args.prealigned,
        align_stop = args.stop_align, tree_stop = args.stop_cat,
        flag_incomplete = bool(not args.missing)
        )
#    elif os.path.isdir(input_check):
 #       fas = collect_files(input_check, 'fa')
  #      fas.extend(collect_files( input_check, 'fasta' ))
#        fas.extend(collect_files( input_check, 'fna' ))
 #       fas.extend(collect_files( input_check, 'fsa' ))
  #      if not fas:
   #         eprint('\nERROR: no {.fa, .fasta, .fna, .fsa} detected', flush = True)
    #        sys.exit(2)
     #   print(flush = True)
      #  for fa in fas:
       #     print(fa, flush = True)
        #    main(
         #       fa, slurm = args.slurm, torque = args.torque, alignment = alignment, fast = args.fast,
          #      project = args.project, output_dir = output, verbose = False, spacer = '\n'
           #     )
#    else:
 #       eprint('\nERROR: invalid input', flush = True)
  #      sys.exit(3)

    outro(start_time)
