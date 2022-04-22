#! /usr/bin/env python3

import os, re, sys, subprocess, argparse
from mycotools.lib.kontools import eprint, vprint, collect_files, formatPath, intro, outro, findExecs
from mycotools.lib.biotools import fa2dict


def mafftRun( fasta, out_dir, hpc, verbose = True ):

    name = fasta
    if '/' in fasta:
        name = os.path.basename(os.path.abspath(fasta)) 
    cmd = 'mafft --auto --thread -1 ' + fasta + ' > ' + out_dir + '/' + name + '.mafft 2>>/dev/null'

    if not hpc:
        vprint('\nAligning', v = verbose, flush = True)
        run_mafft = subprocess.call( cmd, shell = True, stdout = subprocess.PIPE )
        if run_mafft != 0:
            eprint( '\nERROR: mafft failed\n' , flush = True)
            os.remove( out_dir + '/' + name )
            sys.exit( 2 )
    else:
        with open( out_dir + '/mafft.sh', 'w' ) as out:
            if '#PBS' in hpc:
                out.write( hpc + '\n\n' + cmd + '\ncd ' + out_dir + '\nqsub clipkit.sh')
            else:
                 out.write( hpc + '\n\n' + cmd + '\ncd ' + out_dir + '\nsbatch clipkit.sh')
               
        
    return out_dir + '/' + name + '.mafft'


def trimRun( mafft, out_dir, hpc, verbose ):

    name2 = os.path.basename(os.path.abspath(re.sub( r'\.mafft$', '.clipkit', mafft)))
    cmd = ['clipkit', mafft, '--output', out_dir + '/' + name2]

    if not hpc:

        vprint('\nTrimming', v = verbose, flush = True)
        run_clipkit = subprocess.call( cmd, shell = True, stdout = subprocess.PIPE )
        if run_clipkit != 0:
            eprint( '\nERROR: `clipkit` failed\n' , flush = True)
            os.remove( out_dir + '/' + name2 )
            sys.exit( 4 )

    else:
        with open( out_dir + '/clipkit.sh', 'w' ) as out:
            if '#PBS' in hpc:
                out.write( hpc + '\n\n' + #cmd1 + '\n' + 
                    ' '.join([str(x) for x in cmd]) + '\n\ncd ' + out_dir + '\nqsub iqtree.sh' )
            else:
                out.write( hpc + '\n\n' + #cmd1 + '\n' + 
                    ' '.join([str(x) for x in cmd]) + '\n\ncd ' + out_dir + '\nsbatch iqtree.sh' )
    
    return out_dir + '/' + name2


def iqtreeRun( clipkit_file, out_dir, hpc, verbose, fast = False ):

    cmd = ['iqtree', '-s', clipkit_file, '-B', '1000', '-T', 'AUTO']
    if fast:
        del cmd[4]
        del cmd[3]
        cmd.append('--fast')

    vprint('\nOutputting bash script `iqtree.sh`.\n', v = verbose, flush = True)
   
    if not hpc:
        run_iqtree = subprocess.call(cmd.split(' '))
        if run_iqtree != 0:
            eprint('\nERROR: iqtree failed\n', flush = True)
            sys.exit(6)
    else: 
        with open( out_dir + '/iqtree.sh', 'w' ) as out:
            out.write( hpc + '\n\n' + ' '.join([str(x) for x in cmd]) )


def main( 
    fasta_path, slurm = False, pbs = False, fast = False, project = '', 
    output_dir = None, verbose = True, alignment = False 
    ):

    hpc = False
    if slurm:
        vprint('\nHPC mode, preparing submission scripts (Slurm)\n', v = verbose, flush = True)
        hpc = '#!/bin/bash\n#SBATCH --time=24:00:00\n#SBATCH --nodes=1\n' + \
            '#SBATCH --ntasks-per-node=4\n#SBATCH -A ' + project + '\n' + \
            '#SBATCH --mem=60000'
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

    fasta_name = os.path.basename( os.path.abspath(fasta_path))
    mafft = out_dir + '/' + fasta_name + '.mafft'
    clipkit = out_dir + '/' + fasta_name + '.clipkit'

    if not os.path.isfile( mafft ) and not alignment:
        mafft = mafftRun( os.path.abspath(fasta_path), out_dir, hpc, verbose )
    elif alignment:
        mafft = fasta_path
    else:
        vprint('\nAlignment exists ...', v = verbose, flush = True)
    if not os.path.isfile( clipkit ):
        clipkitOut = trimRun( mafft, out_dir, hpc, verbose )
    else:
        vprint('\nTrim exists ...', v = verbose , flush = True)

    iqtreeRun( clipkitOut, out_dir, hpc, verbose, args.fast )

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
        help = 'Run IQTree fast option' )
    parser.add_argument( '-p', '--pbs', action = 'store_true', \
        help = 'Submit ALL steps to Torque.' )
    parser.add_argument( '-s', '--slurm', action = 'store_true', \
        help = 'Submit ALL steps to Slurm' )
    parser.add_argument( '-A', '--project', help = 'HPC project' )
    parser.add_argument( '-o', '--output' )
    args = parser.parse_args()

    output = formatPath( args.output, isdir = True )
    args_dict = {
        'Fasta': args.fasta, 'Alignment': args.alignment, 'Fast': args.fast,
        'Torque': args.pbs, 'Slurm': args.slurm, 'HPC project': args.project, 'Output': output
        }
    start_time = intro( 'Fasta2Tree', args_dict )

    if not args.fasta and not args.alignment:
        eprint('\nERROR: no input', flush = True)
        sys.exit(7)

    findExecs( ['iqtree', 'mafft', 'clipkit'] )

    if args.fasta:
        input_check = formatPath(args.fasta)
        alignment = False
    else:
        input_check = formatPath(args.alignment)
        alignment = True
    if os.path.isfile(input_check):
        main( 
            input_check, slurm = args.slurm, 
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
                project = args.project, output_dir = output, verbose = False
                )
    else:
        eprint('\nERROR: invalid input', flush = True)
        sys.exit(3)

    outro( start_time )
