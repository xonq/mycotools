#! /usr/bin/env python3

# NEED to implement fa2clus?
# NEED to ignore non fasta inputs
# NEED to work as a standalone script

import os
import re
import sys
import shutil
import argparse
import subprocess
import contextlib
from collections import defaultdict
from mycotools.lib.kontools import eprint, vprint, collect_files, \
    format_path, intro, outro, findExecs, mkOutput, multisub
from mycotools.lib.biotools import fa2dict, dict2fa
from clipkit import clipkit
try:
    from ete3 import Tree
except ImportError:
    eprint('WARNING: ete3 not installed.\nInstall ete3 into your ' \
                + 'conda environment via `conda install ete3`')


# adopted from https://stackoverflow.com/a/2829036 for verbosity control
class DummyFile(object):
    def write(self, x): pass

@contextlib.contextmanager
def nostdout():
    save_stdout = sys.stdout
    sys.stdout = DummyFile()
    yield
    sys.stdout = save_stdout


class PhyloError(Exception):
    pass


def run_mafft(name, fasta, out_dir, hpc, verbose = True, 
              cpus = 1, spacer = '\t'):
    """Call Mafft to align an inputted fasta and either output commands for
    sequential execution on an HPC, or run directly through subprocess"""

    cmd = ['mafft', '--auto', '--thread', '-' + str(cpus), fasta]

    # run the command directly
    if not hpc:
        print(spacer + 'Aligning', flush = True)
        with open(out_dir + name + '.mafft', 'w') as out_file:
            if verbose:
                run_mafft = subprocess.call(cmd, stdout = out_file)
            else:
                run_mafft = subprocess.call(cmd, stdout = out_file, 
                                        stderr = subprocess.DEVNULL)
         
        if run_mafft != 0:
            eprint(spacer + '\tERROR: mafft failed: ' + str(run_mafft), flush = True)
            if os.path.isfile(out_dir + name + '.mafft'):
                os.remove(out_dir + name + '.mafft')	
            raise PhyloError
    # prepare an execution script for the HPC submission
    else:
        with open(out_dir + name + '_mafft.sh', 'w') as out:
            if '#PBS' in hpc:
                out.write(hpc + '\n\n' + ' '.join([str(x) for x in cmd]) + \
                ' > ' + out_dir + '/' + name + '.mafft' + \
                '\ncd ' + out_dir + '\nqsub ' + name + '_clipkit.sh')
            else:
                 out.write(hpc + '\n\n' + ' '.join([str(x) for x in cmd]) + \
		' > ' + out_dir + '/' + name + '.mafft' + \
                '\ncd ' + out_dir + '\nsbatch ' + name + '_clipkit.sh')
        
    return out_dir + '/' + name + '.mafft'


def run_clipkit(name, mafft_name, out_dir, hpc, gappy, 
            verbose, cpus = 1, spacer = '\t'):
    """Execute ClipKIT from a complete Mafft run"""

    clipkit_out_name = out_dir + os.path.basename(mafft_name) + '.clipkit'
    if gappy:
        mode = "gappy"
        cmd = ['clipkit', mafft_name, '--output', clipkit_out_name, 
               '-m', 'gappy', str(gappy)]
    else:
        cmd = ['clipkit', mafft_name, '--output', clipkit_out_name]
        mode = "smart-gap"
        gappy = None


    # execute immediately
    if not hpc:
        print(spacer + 'Trimming', flush = True)
        if not verbose:
            run_clipkit = subprocess.call(cmd, stdout = subprocess.PIPE)
        else:
            run_clipkit = subprocess.call( 
                cmd, stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL
                )
#            clipkit_out = clipkit.execute(
 #       	    input_file=mafft_name,
#	            output_file=clipkit_out_name,
 #               output_file_format='fasta',
  #              input_file_format='fasta',
   #             use_log = False,
    #            complement = False,
	 #           mode=mode, gaps=gappy
      #          )
       # else:
        #    with nostdout():
         #       clipkit_out = clipkit.execute(
        #	        input_file=mafft_name,
	     #           output_file=clipkit_out_name,
          #          output_file_format='fasta',
           #         input_file_format='fasta',
            #        use_log = False,
             #       complement = False,
	          #      mode=mode, gaps=gappy
               #     )

        # no output file, the run failed
        if not os.path.isfile(clipkit_out_name):
            eprint(spacer + '\tERROR: `clipkit` failed: ' + str(run_clipkit), 
                   flush = True)  
            raise PhyloError
    # prepare a command for HPC execution
    else:
        cmd = ['clipkit', mafft_name, '--output', clipkit_out_name]
        if gappy:
            cmd.extend(['-m', 'gappy', '-g', str(gappy)])
        with open(out_dir + name + '_clipkit.sh', 'w') as out:
            if '#PBS' in hpc:
                out.write( hpc + '\n\n' + #cmd1 + '\n' + 
                    ' '.join([str(x) for x in cmd]) + '\n\ncd ' + out_dir \
                    + '\nqsub ' + name + '_tree.sh' )
            else:
                out.write( hpc + '\n\n' + #cmd1 + '\n' + 
                    ' '.join([str(x) for x in cmd]) + '\n\ncd ' + out_dir \
                    + '\nsbatch ' + name + '_tree.sh' )
    
    return clipkit_out_name


def run_mf(clipkit_files, out_dir, constraint = False, 
           verbose = False, cpus = 1, spacer = '\t', scripts = False):
    """Run ModelFinder on its own in preparation for multigene phylogeny
    reconstruction"""
    # set the CPUs for each ModelFinder run arbitrarily
    cpus_per_cmd = 3
    # determine the concurrent ModelFinder runs that are possible
    concurrent_cmds = round((cpus - 1)/4 - 0.5) # round down

    # prepare a scaffold for each command
    cmds = [['iqtree', '-m', 'MFP+MERGE', '-nt', 'AUTO',
             '-B', '1000', '-s', f_, '--prefix', 
             out_dir + os.path.basename(f_)] \
            for f_ in clipkit_files]
    # append the topological constraint to each command if it is present
    if constraint:
        [x.extend(['-g', constraint]) for x in cmds]

    # if scripts of each ModelFinder run are desired to be run independently
    if scripts:
        return {os.path.basename(v): cmds[i] for i, v in enumerate(clipkit_files)}

    # otherwise parallelize and run ModelFinder directly
    multisub(cmds, verbose = verbose, processes = concurrent_cmds)

    # parse and identify the evolutionary models determined by ModelFinder
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
    """Prepare a NEXUS file of the concatenated sequence coordinates and their
    associated evolutionary models to feed into IQTree for multigene phylogeny
    reconstruction."""
    nex_data = '#nexus\nbegin sets;'
    index = 1
    for i, trimmed_f in enumerate(list(models.keys())):
        # determine the length of the fasta
        trim_fa = fa2dict(trimmed_f)
        len0 = len(trim_fa[list(trim_fa.keys())[0]]['sequence'])
        # if some alignments aren't the same length then there is some cryptic
        # issue, likely user-caused
        if not all(len(x['sequence']) == len0 for x in trim_fa.values()):
            eprint(spacer + '\tERROR: alignment sequences are not same length ' \
                   + trimmed_f, flush = True)
            sys.exit(17)
        # adjust the index of the coordinates of each sequence based on the
        # previous sequences' length
        new_index = index + len0
        nex_data += '\n    charset part' + str(i + 1) + ' = ' \
                  + str(index) + '-' + str(new_index - 1) + ';'
        index = new_index

        # concatenate the sequences
        for seq, data in trim_fa.items():
            if '_' in seq:
                ome = seq[:seq.find('_')]
            else:
                ome = seq
            if ome in concat_fa:
                concat_fa[ome]['sequence'] += data['sequence']

    # build the NEXUS string
    nex_data += '\n    charpartition mine = '
    for i, model in enumerate(models.values()):
        nex_data += model + ':part' + str(i+1) + ','
    nex_data = nex_data[:-1] + ';\nend;'

    return nex_data, concat_fa 


def run_partition_tree(fa_file, nex_file, constraint,
                       tree_stop, verbose, cpus, spacer = '\t'):
    """Prepare and execute the command for multigene phylogeny
    reconstruction"""
    cmd = ['iqtree', '-s', fa_file, '-p', nex_file, '-nt', str(cpus),
           '-B', '1000', '--sampling', 'GENESITE', '-safe'] # this is not using the
           # specified CPUs, but instead the most efficient
    # add a topological constraint if necessary
    if constraint:
        cmd.extend(['-g', constraint])
    if tree_stop:
        print(spacer + 'Concatenated nexus and fasta outputted', 
              flush = True)
        print(' '.join(cmd), flush = True)
        sys.exit(0)

    print(spacer + 'Tree building', flush = True)
    if verbose:
        run_tree = subprocess.call(cmd)
    else:
        run_tree = subprocess.call(cmd, stdout = subprocess.DEVNULL,
                                   stderr = subprocess.DEVNULL)
    return run_tree


def run_tree_reconstruction(prefix, clipkit_file, out_dir, hpc, constraint,
            verbose, fast = False, cpus = 1, spacer = '\t'):
    """Manage and execute phylogeny reconstruction from an inputted ClipKIT
    output trimmed alignment."""

    tree_file = out_dir + os.path.basename(clipkit_file)
    # prepare a fasttree ommand
    if fast:
        cmd = ['fasttree', '-out', tree_file + '.treefile', clipkit_file]
        if constraint:
            cmd.extend(['-constraints', constraint])
    # prepare an iqtree command
    else:
        # the following will identify the optimum number of threads
        cmd = ['iqtree', '-s', clipkit_file, '-B', '1000', '-nt', str(cpus),
               '--prefix', tree_file]
        if constraint:
            cmd.extend(['-g', constraint])
    
    # execute in the current terminal
    if not hpc:
        print(spacer + 'Tree building', flush = True)
        if fast:
            cmd[2] += '.tmp'
        if verbose:
            run_tree = subprocess.call(cmd)
        else:
            run_tree = subprocess.call(cmd, 
                stdout = subprocess.DEVNULL, 
                stderr = subprocess.DEVNULL)
        if fast:
            shutil.move(cmd[2], cmd[2][:-4])
        if run_tree != 0:
            eprint(spacer + '\tERROR: tree failed: ' + str(run_tree), 
            flush = True)
            raise PhyloError
    # prepare an .sh file for user execution, or sequential execution from the
    # previous scripts
    else: 
        vprint('\nOutputting bash script `tree.sh`.\n', v = verbose, flush = True)
        with open(f'{out_dir}../{prefix}_tree.sh', 'w') as out:
            out.write(hpc + '\n\n' + ' '.join([str(x) for x in cmd]))


def prep_fasta_path_input(fasta_path, output_dir):
    """Prepare the output directories and collect the files from an input that
    is a fasta path (file or directory)"""
    # start from a file input
    if os.path.isfile(fasta_path):
        if not output_dir:
            dir_name = re.sub(r'\..*?$', '_tree', 
                        os.path.basename(os.path.abspath(fasta_path)))
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
    # start from a directory of fastas
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
    return out_dir, wrk_dir, files


def prep_fasta_list_input(fastas, output_dir):
    if not output_dir:
        out_dir = mkOutput('./', 'fa2tree')
    else:
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)
        out_dir = format_path(output_dir)
    wrk_dir = out_dir + 'working/'
    if not os.path.isdir(wrk_dir):
        os.mkdir(wrk_dir)
    files = []
    for path in fastas:
        if os.path.isdir(path):
            files.extend([format_path(x) for x in collect_files(path, '*')])
        else:
            files.append(path)
    return out_dir, wrk_dir, files


def nonfasta2fasta(files, conv_dir):
    """Use BioPython to convert non-fastas into fastas"""
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
    return files


def identify_incomplete_files(files, flag_incomplete, wrk_dir):
    """Identify fastas that do not have every genome represented in any other
    file and handle this as an error or warning and genome deletion."""
    ome2fa2gene = defaultdict(dict)
    incomp_omes = {}
    # identify files with missing sequences
    for f in files:
        if f.endswith('/'):
            f = f[:-1]
        fa_dict = fa2dict(f)
        incomp_omes[os.path.basename(f)] = []
        # populate a hash with each gene for each ome associated with each
        # fasta
        for seq in fa_dict:
            ome = seq[:seq.find('_')]
            if f not in ome2fa2gene[ome]:
                incomp_omes[os.path.basename(f)].append(ome)
                ome2fa2gene[ome][f] = seq
            # a multigene partition cannot be reconstructed for a genome with
            # multiple genes in the same alignment
            else:
                eprint('\nERROR: multiple sequences for a ' \
                     + 'single ome in ' + f, flush = True)
                sys.exit(5)

    # identify the files with missing genomes and the missing omes themselves
    del_omes = set()
    incomp_files, comp_files = {}, []
    for f, omes in incomp_omes.items():
        if len(omes) < len(ome2fa2gene):
            missing_omes = set(ome2fa2gene.keys()).difference(set(omes))
            incomp_files[f] = sorted(missing_omes)
            del_omes = del_omes.union(missing_omes)
        else:
            comp_files.append(f)

    # if there are incomplete files, then run through the flagging process
    if incomp_files:
        incomp_files = {k: v for k,v in sorted(incomp_files.items(), 
                                               key = lambda x: len(x[1]))}
        # it is an error if genomes are missing and it isn't explicitly
        # permitted
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
        # exit if incomplete flags
        if flag_incomplete:
            sys.exit(6)
        # otherwise remove the genomes from each fasta that are not present in
        # all of them
        else:
            new_files = []
            for fasta in files:
                fasta_name = os.path.basename(fasta)
                fa = fa2dict(fasta)
                new_fa = {k: v for k, v in fa.items() \
                          if k[:k.find('_')] not in del_omes}
                with open(wrk_dir + fasta_name + '.tmp', 'w') as out:
                    out.write(dict2fa(new_fa))
                os.rename(wrk_dir + fasta_name + '.tmp', 
                          wrk_dir + fasta_name)
                new_files.append(wrk_dir + fasta_name)
            files = new_files

    return files, ome2fa2gene, del_omes


def algn_mngr(start, files, wrk_dir, hpc, gappy, alignment,
              verbose, cpus, spacer):
    """Manage the alignment and trimming executor for inputted fasta files"""

    trimmed_files = []
    for fasta in files:
        name = re.search(r'.*?/*([^/]+)/*$', fasta)[1]
        fasta_name = os.path.basename(os.path.abspath(fasta))
        clipkit = wrk_dir + fasta_name + '.clipkit'

        # if the starting point is a non-aligned fasta
        if start == 0:
            mafft = wrk_dir + fasta_name + '.mafft'
            clipkit = mafft + '.clipkit'
            try:
                # check if it exists
                if os.stat(mafft):
                    vprint('\nAlignment exists', v = verbose, flush = True)
                else:
                    raise ValueError
            except (FileNotFoundError, ValueError) as e:
                # otherwise run the alignment
                if not alignment:
                    mafft = run_mafft(name, os.path.abspath(fasta), wrk_dir, hpc, 
                                     verbose, cpus, spacer = spacer)
                else:
                    mafft = fasta
        # else the starting point is an alignment or trimmed alignment
        else:
            mafft = fasta
            clipkit = mafft + '.clipkit'

        # if the starting point is an alignment or trimmed alignment
        if start == 0 or start == 1: 
            try:
                # check for a trimmed fasta
                if os.stat(clipkit):
                    vprint('\nTrim exists', v = verbose, flush = True)
                    clipkit_out = clipkit
                else:
                    raise ValueError
            except (FileNotFoundError, ValueError) as e:
                clipkit_out = run_clipkit(name, mafft, wrk_dir, hpc, gappy,
                                          verbose, spacer = spacer)
            trimmed_files.append(clipkit_out)
        else:
            trimmed_files.append(fasta)

    return trimmed_files


def convert_seq_to_ome_name(trimmed_files, conv_dir):
    """Convert the inputted sequence name to its genome name"""
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

    return new_trimmed_files


def extract_supported(trimmed_files, min_mean_support, out_dir):
    """Extract trees that meet the minimum summary support value"""
    out_files = []
    for f_ in trimmed_files:
        t_path = f'{out_dir}working/{os.path.basename(f_)}.contree'
        supports = check_tree_support(t_path)
        mean_support = sum(supports)/len(supports)
        if mean_support >= min_mean_support:
            print(f'\t{os.path.basename(t_path)} {mean_support} passed', 
                  flush = True)
            out_files.append(f_)
        else:
            print(f'\t{os.path.basename(t_path)} {mean_support} failed', 
                  flush = True)
    print(f'\t{len(out_files)} ({len(out_files)/len(trimmed_files)*100}%)' \
         + ' passed', flush = True)
    return out_files


def multigene_mngr(align_stop, trimmed_files, wrk_dir, constraint, out_dir,
                  ome2fa2gene, del_omes, verbose, cpus, spacer, hpc_prep,
                  min_mean_support, tree_stop):
    """Manage the execution of the multigene phylogeny execution from
    ModelFinding forward"""

    # stop at the end of the alignment and output a command
    if align_stop:
        # run modelfinding script generator
        script_dict = run_mf(trimmed_files, wrk_dir, constraint,
                        verbose, cpus, spacer = spacer + '\t',
                        scripts = True)
        for f, cmd in script_dict.items():
            # output a script, SLURM formatted or otherwise
            with open(out_dir + f + '.sh', 'w') as out:
                if hpc_prep: # doesnt work with PBS/Torque
                    out.write(f'{hpc_prep}\n#SBATCH --output={f}_tree\n' \
                            + f'#SBATCH --error={f}_tree\n' \
                            + f'#SBATCH --job-name={f}_tree\n')
                out.write(' '.join(cmd))
        print(spacer + 'Alignments outputed', flush = True)
        sys.exit(0)

    print(spacer + 'Model finding', flush = True)
    models = run_mf(trimmed_files, wrk_dir, constraint,
                    verbose, cpus, spacer = spacer + '\t')
    if min_mean_support:
        print(spacer + 'Extracting passing trees', flush = True)
        passing = extract_supported(trimmed_files, min_mean_support, out_dir)
        models = {k: models[k] for k in passing}

    # build the concatenated NEXUS 
    print(spacer + 'Concatenating', flush = True)
    empty_concat_fa = {ome: {'description': '', 'sequence': ''}
                 for ome in ome2fa2gene if ome not in del_omes}
    nex_data, concat_fa = prepare_nexus(empty_concat_fa, models)
    with open(out_dir + 'concatenated.fa', 'w') as out:
        out.write(dict2fa(concat_fa))
    with open(out_dir + 'concatenated.nex', 'w') as out:
        out.write(nex_data)

    # stop and allow the user to build the tree

    # run the multigene phylogeny reconstruction
    run_partition_tree(out_dir + 'concatenated.fa', 
                       out_dir + 'concatenated.nex', constraint,
                       tree_stop, verbose, cpus, spacer = spacer)


def check_tree_support(tree_file):
    """Acquire the support values of a newick"""
    t = Tree(tree_file)
    supports = []
    # parse the tree
    for n in t.traverse():
        # skip leaves and roots with null support
        if not n.is_leaf():
            if not n.is_root():
                supports.append(n.support)

    # adjust supports if they are 0 - 1
    if any(x > 1 for x in supports): 
        supports = [x/100 for x in supports]

    return supports


def main( 
    fasta_path, slurm = False, torque = False, fast = False, project = '', 
    output_dir = None, verbose = True, alignment = False, constraint = False,
    ome_conv = False, mem = '60GB', cpus = 1, spacer = '\t\t', align_stop = False, 
    tree_stop = False, flag_incomplete = True, start = 0, partition = False,
    gappy = None, min_mean_support = 0
    ):
    """Python entry-point for fa2tree. Build a phylogeny from an inputted fasta
    of homologs. Align (mafft) -> trim (clipkit) -> ModelFinder (multigene
    mode) -> tree reconstruction (fasttree/iqtree)."""

    # prepare the scaffold of an HPC command if that is requested
    hpc, hpc_prep = False, False
    if slurm:
        vprint('\nHPC mode, preparing submission scripts (Slurm)\n', 
               v = verbose, flush = True)
        hpc_prep = '#!/bin/bash\n#SBATCH --time=24:00:00\n' \
                 + '#SBATCH --nodes=1\n#SBATCH --ntasks-per-node=' \
                 + f'{cpus}\n#SBATCH -A {project}\n' \
                 + f'#SBATCH --mem={mem}'
#            '\n\nsource activate ' + source
    elif torque:
        vprint('\nHPC mode, preparing submission scripts (PBS)\n', 
               v = verbose, flush = True)
        hpc_prep = '#PBS -l walltime=10:00:00\n#PBS -l nodes=1:ppn=4\n#PBS ' \
                 + '-A ' + project
    if not partition and hpc_prep:
        hpc = hpc_prep

    output_dir_prep = os.getcwd() + '/'
    # the fasta data should be a path if it is a string
    if isinstance(fasta_path, str):
        out_dir, wrk_dir, files = prep_fasta_path_input(fasta_path, output_dir)
    # otherwise the fasta should be a provided list of fastas
    else:
        out_dir, wrk_dir, files = prep_fasta_list_input(fasta_path, output_dir)

    # create a directory for converting files
    conv_dir = wrk_dir + 'conv/'
    if not os.path.isdir(conv_dir):
        os.mkdir(conv_dir)

    # check for non fasta inputs - if these files exist in a directory then
    # they will be used anyway
    if any(f_.endswith(('.nexus', '.nex', '.nxs', 
                        '.phylip', '.phy', '.ph',
                        '.clustal', '.clus',)) \
           for f_ in files):
        print(spacer + 'Converting to fastas', flush = True)
        files = nonfasta2fasta(files, conv_dir)

    # only proceed with fastas from the inputted files
    files = [x for x in files \
             if x.endswith(('fasta', 'fa', 'fna', 'faa', 'fsa'))]

    if partition: #multigene partition mode
        files, ome2fa2gene, del_omes = identify_incomplete_files(files,
                                                                 flag_incomplete,
                                                                 wrk_dir)

    # manage the alignment and trimming
    trimmed_files = algn_mngr(start, files, wrk_dir, hpc, gappy, 
                              alignment, verbose, cpus, spacer)

    # convert the sequences to genome codenames    
    if ome_conv:
        trimmed_files = convert_seq_to_ome_name(trimmed_files, conv_dir)

    # reconstruct a single gene phylogeny
    if not partition:
        for trimmed_f in trimmed_files:
            name_prep = re.search(r'.*?/*([^/]+)/*$', trimmed_f)[1]
            name = re.sub(r'\.mafft\.clipkit$', '', name_prep)
            run_tree_reconstruction(name, trimmed_f, out_dir, hpc, constraint,
                    verbose, fast, cpus, spacer = spacer)
    # run multigene partition model phylogeny mode
    else:
        multigene_mngr(align_stop, trimmed_files, wrk_dir, constraint, out_dir,
                      ome2fa2gene, del_omes, verbose, cpus, spacer, hpc_prep,
                      min_mean_support, tree_stop)

    if hpc:
        if slurm:
            print('\nStart pipeline via `sbatch <STARTSTEP>.sh` in ' \
                + out_dir + '\n')
        else:
            vprint('\nStart pipeline via `qsub <STARTSTEP>.sh` in ' + out_dir + '\n',
                v = verbose)


def cli():

    parser = argparse.ArgumentParser( 
        description = 'Takes in a multifasta, or directory/list of ' \
                    + 'fastas/alignments [partition] aligns (mafft) ' \
                    + 'trims (clipkit), and reconstructs phylogeny ' \
                    + '(iqtree/fastree). When using non-Mycotools ' \
                    + 'fastas for multigene partition analysis, all ' \
                    + 'accessions must be formatted as ' \
                    + '"<GENOME>_<ACCESSION>" and each <GENOME> needs ' \
                    + 'to be in ALL fastas')

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
        help = 'Fasttree')
    s_opt.add_argument('--slurm', action = 'store_true', \
        help = 'Submit ALL steps to Slurm')
    s_opt.add_argument('--torque', action = 'store_true', \
        help = 'Submit ALL steps to Torque')
    s_opt.add_argument('-A', '--project', help = 'HPC project')

    m_opt = parser.add_argument_group('Partition tree options')
    m_opt.add_argument('-p', '--partition', action = 'store_true',
        help = 'Multigene partition mode')
    m_opt.add_argument('-sa', '--stop_align', action = 'store_true',
        help = 'Stop after aligning and trimming, output tree scripts')
    m_opt.add_argument('-s', '--support', default = 0, type = float,
        help = '[0 < -s < 1] Minimum mean node support to consider in final phylogeny' \
             + '. Adopted from Gluck-Thaler et al. 2020 Scientific Advances')
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

    if args.support:
        if args.support > 1 or args.support < 0:
            eprint('\nERROR: --support must be between 0 and 1')
            sys.exit(3)
        
    output = format_path(args.output)
    if output:
        if not output.endswith('/'):
            output += '/' # bring to mycotools path expectations

    if args.gappy:
        if args.gappy > 1:
            eprint('\nERROR: gappy threshold must be less than 1')
            sys.exit(3)
    
    args_dict = {
        'Fasta': args.input, 'Fast': args.fast, 
        'Multigene partition': args.partition, 'Input alignments': bool(args.align),
        'Input trimmed': bool(args.trim), 'Constraint': args.constraint,
        'Ome conversion': args.ome, 'Gappy threshold': args.gappy,
        'Minimum mean support': args.support,
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
        ome_conv = args.ome, gappy = args.gappy,
        align_stop = args.stop_align, tree_stop = args.stop_cat,
        flag_incomplete = bool(not args.missing), start = start,
        partition = bool(args.partition), cpus = args.cpus,
        min_mean_support = args.support
        )

    outro(start_time)

if __name__ == '__main__':
    cli()
