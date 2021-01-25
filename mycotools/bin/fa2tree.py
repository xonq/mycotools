#! /usr/bin/env python3

import os, re, sys, subprocess, argparse
from mycotools.lib.kontools import eprint, vprint, collect_files
from mycotools.lib.fastatools import fasta2dict


def mafftRun( fasta, out_dir, hpc ):

    name = fasta
    if '/' in fasta:
        name = os.path.basename(os.path.abspath(fasta)) 
    cmd = 'mafft --auto --thread -1 ' + fasta + ' > ' + out_dir + '/' + name + '.mafft 2>>/dev/null'

    if not hpc:
        vprint('\nAligning', v = verbose)
        run_mafft = subprocess.call( cmd, shell = True, stdout = subprocess.PIPE )
        if run_mafft != 0:
            eprint( '\nERROR: mafft failed\n' )
            os.remove( out_dir + '/' + name )
            sys.exit( 2 )
    else:
        with open( out_dir + '/mafft.sh', 'w' ) as out:
            out.write( hpc + '\n\n' + cmd + '\ncd ' + out_dir + '\nqsub trimal.sh')
        
    return out_dir + '/' + name + '.mafft'


def trimalRun( mafft, out_dir, hpc, tree, verbose ):

    name1 = os.path.basename( os.path.abspath(mafft.replace('mafft', 'prelim_trimal')))
    name2 = re.sub( r'\.prelim_trimal$', '.trimal', name1 )
    cmd1 = 'trimal -in ' + mafft + ' -automated1 -out ' + out_dir + '/' + name1
    cmd2 = 'trimal -in ' + out_dir + '/' + name1 + ' -out ' + out_dir + '/' + name2 + ' -gt 0.25'

    if not hpc:
        if not os.path.isfile( out_dir + '/' + name1 ):
            vprint('\nTrimming step 1/2', v = verbose)
            run_trimal1 = subprocess.call( cmd1, shell = True, stdout = subprocess.PIPE )
            if run_trimal1 != 0:
                os.remove( out_dir + '/' + name1 )
                eprint( '\nERROR: preliminary `trimal` failed\n' )
                sys.exit( 3 )

        vprint('\nTrimming step 2/2', v = verbose)
        run_trimal2 = subprocess.call( cmd2, shell = True, stdout = subprocess.PIPE )
        if run_trimal2 != 0:
            eprint( '\nERROR: `trimal` failed\n' )
            os.remove( out_dir + '/' + name2 )
            sys.exit( 4 )

    else:
        with open( out_dir + '/trimal.sh', 'w' ) as out:
            out.write( hpc + '\n\n' + cmd1 + '\n' + cmd2 + '\n\ncd ' + out_dir + '\nqsub ' + tree + '.sh' )
    
    return out_dir + '/' + name2


def fasttreeRun( trimal, out_dir, hpc, verbose ):


    name = os.path.basename(os.path.abspath(trimal + '.fast.tre'))
    cmd = 'fasttree -lg ' + trimal + ' > ' + out_dir + '/' + name

    if not hpc:
        vprint('\nMaking fasttree ...', v = verbose)
        run_fasttree = subprocess.call( cmd, shell = True, stdout = subprocess.PIPE )
        if run_fasttree != 0:
            eprint( '\nERROR: `fasttree` failed\n' )
            os.remove( out_dir + '/' + name )
            sys.exit( 5 )
    else:
        with open( out_dir + '/fasttree.sh', 'w' ) as out:
            out.write( hpc + '\n\n' + cmd )

    return name


def iqtreeRun( trimal, out_dir, hpc, verbose ):

    cmd = 'iqtree -s ' + trimal + ' -B 1000 -T AUTO'

    vprint('\nOutputting bash script `iqtree.sh`.\n', v = verbose)
    
    with open( out_dir + '/iqtree.sh', 'w' ) as out:
        out.write( hpc + '\n\n' + cmd )


def main( fasta_path, tree, slurm = False, torque = False, project = '', output_dir = None, verbose = True ):

    hpc = False
    if slurm:
        vprint('\nHPC mode, preparing submission scripts (Slurm)\n', v = verbose)
        hpc = '#SBATCH --time=10:00:00\n#SBATCH --nodes=1\n' + \
            '#SBATCH --ntaskspernode=4\n#SBATCH -A ' + project
#            '\n\nsource activate ' + source
    elif pbs:
        vprint('\nHPC mode, preparing submission scripts (PBS)\n', v = verbose)
        hpc = '#PBS -l walltime=10:00:00\n#PBS -l nodes=1:ppn=4\n#PBS -A ' + project

    if not output_dir:
        output_dir = os.path.dirname(os.path.abspath(fasta_path))
        
    dir_name = re.sub( r'\..*?$', '_tree', os.path.basename( os.path.abspath(fasta_path) ))
    out_dir = output_dir + '/' + dir_name
    if not os.path.isdir( out_dir ):
        os.mkdir( out_dir )

    fasta_name = os.path.basename( os.path.abspath(fasta_path))
    mafft = out_dir + '/' + fasta_name + '.mafft'
    trimal = out_dir + '/' + fasta_name + '.trimal'

    if not os.path.isfile( mafft ):
        mafft = mafftRun( os.path.abspath(fasta_path), out_dir, hpc, verbose )
    else:
        vprint('\nAlignment exists ...', v = verbose)
    if not os.path.isfile( trimal ):
        trimal = trimalRun( mafft, out_dir, hpc, tree, verbose )
    else:
        vprint('\nTrim exists ...', v = verbose )

    if tree == 'fasttree':
        fasttreeRun( trimal , out_dir, hpc, verbose )
    elif tree == 'iqtree':
        iqtreeRun( trimal, out_dir, hpc, verbose )

    if hpc:
        if slurm:
            vprint('\nStart pipeline via `sbatch <STARTSTEP>.sh` in ' + out_dir + '\n', v = verbose)
        else:
            vprint('\nStart pipeline via `qsub <STARTSTEP>.sh` in ' + out_dir + \
            '\n', v = verbose)

   


if __name__ == '__main__':

    parser = argparse.ArgumentParser( 
        description = 'Takes in a multifasta or directory of multifastas, ' + \
        ' aligns, trims, and runs treebuilding: {`fasttree`, `iqtree`, `raxml` (future)}.' )
    parser.add_argument( '-i', '--input', required = True, \
        help = 'Fasta or directory of fastas' )
    parser.add_argument( '-t', '--tree', default = 'fasttree', \
        help = 'Tree-building software. { fasttree, iqtree }' )
    parser.add_argument( '-p', '--pbs', action = 'store_true', \
        help = 'Submit ALL steps to Torque.' )
    parser.add_argument( '-s', '--slurm', action = 'store_true', \
        help = 'Submit ALL steps to Slurm' )
    parser.add_argument( '-A', '--project', help = 'HPC project' )
    parser.add_argument( '-o', '--output' )
    args = parser.parse_args()

    output = formatPath( args.output, isdir = True )
    args_dict = {
        'Input': args.input, 'Tree': args.tree, 'Torque': args.pbs,
        'Slurm': args.slurm, 'HPC project': args.project, 'Output': output
        }
    start_time = intro( 'Fasta2Tree', args_dict )

    if args.tree not in { 'fasttree', 'iqtree' }:
        eprint('\nERROR: invalid tree software: ' + args.tree )
        sys.exit(1)
    if os.path.isfile(args.input):
        main( 
            args.input, args.tree, slurm = args.slurm, 
            torque = args.torque, project = args.project,
            output_dir = output, verbose = True
            )
    elif os.path.isdir(args.input):
        fas = collect_files( args.input, 'fa' )
        fas.extend( collect_files( args.input, 'fasta' ) )
        fas.extend( collect_files( args.input, 'fna' ) )
        fas.extend( collect_files( args.input, 'fsa' ) )
        if not fas:
            eprint('\nERROR: no {.fa, .fasta, .fna, .fsa} detected')
            sys.exit(2)
        print()
        for fa in fas:
            print(fa)
            main(
                fa, args.tree, slurm = args.slurm, torque = args.torque,
                project = args.project, output_dir = output, verbose = False
                )
    else:
        eprint('\nERROR: invalid input')
        sys.exit(3)

    outro( start_time )
