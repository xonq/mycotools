#! /usr/bin/env python3

import os
import re
import sys
import argparse
import subprocess
from mycotools.lib.kontools import eprint, vprint, collect_files, format_path, intro, outro, findExecs
from mycotools.lib.biotools import fa2dict

class PhyloError(Exception):
    pass


def mafftRun( fasta, out_dir, hpc, verbose = True, cpus = 1, spacer = '\t' ):

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
                run_mafft = subprocess.call( cmd, stdout = out_file )
            else:
                run_mafft = subprocess.call( cmd, stdout = out_file, stderr = subprocess.DEVNULL )
         
        if run_mafft != 0:
            eprint(spacer + '\tERROR: mafft failed: ' + str(run_mafft), flush = True)
            if os.path.isfile(out_dir + name + '.mafft'):
                os.remove( out_dir + name + '.mafft' )	
            raise PhyloError
    else:
        with open( out_dir + '/mafft.sh', 'w' ) as out:
            if '#PBS' in hpc:
                out.write( hpc + '\n\n' + ' '.join([str(x) for x in cmd]) + \
                ' > ' + out_dir + '/' + name + '.mafft' + \
                '\ncd ' + out_dir + '\nqsub clipkit.sh')
            else:
                 out.write( hpc + '\n\n' + ' '.join([str(x) for x in cmd]) + \
		' > ' + out_dir + '/' + name + '.mafft' + \
                '\ncd ' + out_dir + '\nsbatch clipkit.sh')
        
    return out_dir + '/' + name + '.mafft'


def trimRun( mafft, out_dir, hpc, verbose, cpus = 1, spacer = '\t' ):

    name2 = os.path.basename(os.path.abspath(re.sub( r'\.mafft$', '.clipkit', mafft)))
    cmd = ['clipkit', mafft, '--output', out_dir + '/' + name2]
    if not hpc:

        print(spacer + 'Trimming', flush = True)
        if verbose:
            run_clipkit = subprocess.call( cmd, stdout = subprocess.PIPE )
        else:
            run_clipkit = subprocess.call( 
                cmd, stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL
                )
        if run_clipkit != 0:
            eprint(spacer + '\tERROR: `clipkit` failed: ' + str(run_clipkit) , flush = True)  
            if os.path.isfile(out_dir + name2):
                os.remove( out_dir + name2 )
            raise PhyloError

    else:
        with open( out_dir + '/clipkit.sh', 'w' ) as out:
            if '#PBS' in hpc:
                out.write( hpc + '\n\n' + #cmd1 + '\n' + 
                    ' '.join([str(x) for x in cmd]) + '\n\ncd ' + out_dir + '\nqsub tree.sh' )
            else:
                out.write( hpc + '\n\n' + #cmd1 + '\n' + 
                    ' '.join([str(x) for x in cmd]) + '\n\ncd ' + out_dir + '\nsbatch tree.sh' )
    
    return out_dir + '/' + name2


def treeRun( clipkit_file, out_dir, hpc, verbose, fast = False, cpus = 1, spacer = '\t' ):

    if fast:
        cmd = ['fasttree', '-out', clipkit_file + '.treefile', clipkit_file]
    else:
        cmd = ['iqtree', '-s', clipkit_file, '-B', '1000', '-T', str(cpus)]

    vprint('\nOutputting bash script `tree.sh`.\n', v = verbose, flush = True)
   
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
        with open( out_dir + '/tree.sh', 'w' ) as out:
            out.write( hpc + '\n\n' + ' '.join([str(x) for x in cmd]) )


def main( 
    fasta_path, slurm = False, pbs = False, fast = False, project = '', 
    output_dir = None, verbose = True, alignment = False, mem = '60GB', cpus = 1,
    spacer = '\t\t'
    ):

    hpc = False
    if slurm:
        vprint('\nHPC mode, preparing submission scripts (Slurm)\n', v = verbose, flush = True)
        hpc = '#!/bin/bash\n#SBATCH --time=24:00:00\n#SBATCH --nodes=1\n' + \
            '#SBATCH --ntasks-per-node=' + str(cpus) + '\n#SBATCH -A ' + project + '\n' + \
            '#SBATCH --mem=' + str(mem)
#            '\n\nsource activate ' + source
    elif pbs:
        vprint('\nHPC mode, preparing submission scripts (PBS)\n', v = verbose, flush = True)
        hpc = '#PBS -l walltime=10:00:00\n#PBS -l nodes=1:ppn=4\n#PBS -A ' + project

    if not output_dir:
        output_dir = os.path.dirname(os.path.abspath(fasta_path))
        dir_name = re.sub( r'\..*?$', '_tree', os.path.basename( os.path.abspath(fasta_path) ))
        out_dir = output_dir + '/' + dir_name
        if not os.path.isdir( out_dir ):
            os.mkdir( out_dir )
        out_dir = format_path(out_dir)
    else:
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)
        out_dir = format_path(output_dir)

    fasta_name = os.path.basename( os.path.abspath(fasta_path))
    mafft = out_dir + '/' + fasta_name + '.mafft'
    clipkit = out_dir + '/' + fasta_name + '.clipkit'

    try:
        if os.stat(mafft):
            vprint('\nAlignment exists ...', v = verbose, flush = True)
        else:
            raise ValueError
    except (FileNotFoundError, ValueError) as e:
        if not alignment:
            mafft = mafftRun( os.path.abspath(fasta_path), out_dir, hpc, verbose, cpus, spacer = spacer )
        else:
            mafft = fasta_path

    try:
        if os.stat(clipkit):
            vprint('\nTrim exists ...', v = verbose, flush = True)
            clipkitOut = clipkit
        else:
            raise ValueError
    except (FileNotFoundError, ValueError) as e:
        clipkitOut = trimRun( mafft, out_dir, hpc, verbose, spacer = spacer )

    treeRun( clipkitOut, out_dir, hpc, verbose, fast, cpus, spacer = spacer )

    if hpc:
        if slurm:
            vprint('\nStart pipeline via `sbatch <STARTSTEP>.sh` in ' + out_dir + '\n',
                v = verbose)
        else:
            vprint('\nStart pipeline via `qsub <STARTSTEP>.sh` in ' + out_dir + '\n',
                v = verbose)


if __name__ == '__main__':

    parser = argparse.ArgumentParser( 
        description = 'Takes in a multifasta or directory of multifastas, ' + \
        ' aligns (mafft), trims (clipkit), and infers phylogeny (iqtree).'
        )
    parser.add_argument( '-f', '--fasta', \
        help = 'Fasta or directory of fastas' )
    parser.add_argument( '-a', '--alignment', \
        help = 'Alignment or directory of alignments (need fasta file extension)' )
    parser.add_argument( '--fast', action = 'store_true', \
        help = 'Fasttree' )
    parser.add_argument( '-p', '--pbs', action = 'store_true', \
        help = 'Submit ALL steps to Torque.' )
    parser.add_argument( '-s', '--slurm', action = 'store_true', \
        help = 'Submit ALL steps to Slurm' )
    parser.add_argument( '-A', '--project', help = 'HPC project' )
    parser.add_argument( '-o', '--output' )
    args = parser.parse_args()

    output = format_path(args.output)
    if output:
        if not output.endswith('/'):
            output += '/' # bring to mycotools path expectations
    args_dict = {
        'Fasta': args.fasta, 'Alignment': args.alignment, 'Fast': args.fast,
        'Torque': args.pbs, 'Slurm': args.slurm, 'HPC project': args.project, 'Output': output
        }
    start_time = intro( 'Fasta2Tree', args_dict )

    if not args.fasta and not args.alignment:
        eprint('\nERROR: no input', flush = True)
        sys.exit(7)

    execs = ['mafft', 'clipkit']
    if args.fast:
        execs.append('fasttree')
    else:
        execs.append('iqtree')
    if args.pbs or args.slurm:
        findExecs(execs, execs)
    else:
        findExecs(execs, execs)

    if args.fasta:
        input_check = format_path(args.fasta)
        alignment = False
    else:
        input_check = format_path(args.alignment)
        alignment = True
    if os.path.isfile(input_check):
        main( 
            input_check, slurm = args.slurm, fast = args.fast,
            pbs = args.pbs, project = args.project,
            output_dir = output, verbose = True, alignment = alignment
            )
    elif os.path.isdir(input_check):
        fas = collect_files( input_check, 'fa' )
        fas.extend( collect_files( input_check, 'fasta' ) )
        fas.extend( collect_files( input_check, 'fna' ) )
        fas.extend( collect_files( input_check, 'fsa' ) )
        if not fas:
            eprint('\nERROR: no {.fa, .fasta, .fna, .fsa} detected', flush = True)
            sys.exit(2)
        print(flush = True)
        for fa in fas:
            print(fa, flush = True)
            main(
                fa, slurm = args.slurm, pbs = args.pbs, alignment = alignment, fast = args.fast,
                project = args.project, output_dir = output, verbose = False, spacer = '\n'
                )
    else:
        eprint('\nERROR: invalid input', flush = True)
        sys.exit(3)

    outro( start_time )
