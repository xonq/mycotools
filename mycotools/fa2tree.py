#! /usr/bin/env python3

import os, re, sys, subprocess, argparse
from mycotools.lib.kontools import eprint
from mycotools.lib.fastatools import fasta2dict

project = 'PAS1046'
source = 'python3'

def mafftRun( fasta, out_dir, hpc ):

    name = fasta
    if '/' in fasta:
        name = os.path.basename(os.path.abspath(fasta)) 
    cmd = 'mafft --auto --thread -1 ' + fasta + ' > ' + out_dir + '/' + name + '.mafft 2>>/dev/null'

    if not hpc:
        print('\nAligning ...')
        run_mafft = subprocess.call( cmd, shell = True, stdout = subprocess.PIPE )
        if run_mafft != 0:
            eprint( '\nERROR: mafft failed\n' )
            os.remove( out_dir + '/' + name )
            sys.exit( 2 )
    else:
        with open( out_dir + '/mafft.sh', 'w' ) as out:
            out.write( hpc + '\n\n' + cmd + '\ncd ' + out_dir + '\nqsub trimal.sh')
        
    return out_dir + '/' + name + '.mafft'


def trimalRun( mafft, out_dir, hpc, tree ):

    name1 = os.path.basename( os.path.abspath(mafft.replace('mafft', 'prelim_trimal')))
    name2 = re.sub( r'\.prelim_trimal$', '.trimal', name1 )
    cmd1 = 'trimal -in ' + mafft + ' -automated1 -out ' + out_dir + '/' + name1
    cmd2 = 'trimal -in ' + out_dir + '/' + name1 + ' -out ' + out_dir + '/' + name2 + ' -gt 0.25'

    if not hpc:
        if not os.path.isfile( out_dir + '/' + name1 ):
            print('\nTrimming step 1/2 ...')
            run_trimal1 = subprocess.call( cmd1, shell = True, stdout = subprocess.PIPE )
            if run_trimal1 != 0:
                os.remove( out_dir + '/' + name1 )
                eprint( '\nERROR: preliminary `trimal` failed\n' )
                sys.exit( 3 )

        print('\nTrimming step 2/2 ...')
        run_trimal2 = subprocess.call( cmd2, shell = True, stdout = subprocess.PIPE )
        if run_trimal2 != 0:
            eprint( '\nERROR: `trimal` failed\n' )
            os.remove( out_dir + '/' + name2 )
            sys.exit( 4 )

    else:
        with open( out_dir + '/trimal.sh', 'w' ) as out:
            out.write( hpc + '\n\n' + cmd1 + '\n' + cmd2 + '\n\ncd ' + out_dir + '\nqsub ' + tree + '.sh' )
    
    return out_dir + '/' + name2


def fasttreeRun( trimal, out_dir, hpc ):


    name = os.path.basename(os.path.abspath(trimal + '.fast.tre'))
    cmd = 'fasttree -lg ' + trimal + ' > ' + out_dir + '/' + name

    if not hpc:
        print('\nMaking fasttree ...')
        run_fasttree = subprocess.call( cmd, shell = True, stdout = subprocess.PIPE )
        if run_fasttree != 0:
            eprint( '\nERROR: `fasttree` failed\n' )
            os.remove( out_dir + '/' + name )
            sys.exit( 5 )
    else:
        with open( out_dir + '/fasttree.sh', 'w' ) as out:
            out.write( hpc + '\n\n' + cmd )

    return name


def iqtreeRun( trimal, out_dir, hpc ):

    cmd = 'iqtree -s ' + trimal + ' -B 1000 -T AUTO'

    print('\nOutputting bash script `iqtree.sh`.\n')
    
    with open( out_dir + '/iqtree.sh', 'w' ) as out:
        out.write( hpc + '\n\n' + cmd )



if __name__ == '__main__':

    parser = argparse.ArgumentParser( description = 'Takes in a multifasta, aligns, trims, and runs treebuilding: ' \
        + '`fasttree`, `iqtree`, `raxml` (future).' )
    parser.add_argument( '-i', '--input', required = True )
    parser.add_argument( '-t', '--tree', default = 'fasttree', \
        help = 'Tree-building software. DEFAULT: fasttree' )
    parser.add_argument( '-p', '--pbs', action = 'store_true', \
        help = 'Submit ALL steps to Torque.' )
    parser.add_argument( '-s', '--slurm', action = 'store_true', \
        help = 'Submit ALL steps to Slurm' )
    args = parser.parse_args()

    if not os.path.isfile(args.input):
        eprint('ERROR: input file invalid')
        sys.exit( 1 )
    elif args.tree not in { 'fasttree', 'iqtree' }:
        eprint('ERROR: invalid tree software: ' + args.tree )
        sys.exit( 1 )

    hpc = False
    if args.slurm:
        print('\nHPC mode, preparing submission scripts (Slurm)\n')
        hpc = '#SBATCH --time=10:00:00\n#SBATCH --nodes=1\n' + \
            '#SBATCH --ntaskspernode=4\n#SBATCH -A ' + project
#            '\n\nsource activate ' + source
    elif args.pbs:
        print('\nHPC mode, preparing submission scripts (PBS)\n')
        hpc = '#PBS -l walltime=10:00:00\n#PBS -l nodes=1:ppn=4\n#PBS -A ' + project



    output_dir = os.path.dirname(os.path.abspath(args.input))
        
    dir_name = re.sub( r'\..*?$', '_tree', os.path.basename( os.path.abspath( args.input ) ))
    out_dir = output_dir + '/' + dir_name
    if not os.path.isdir( out_dir ):
        os.mkdir( out_dir )

    fasta_name = os.path.basename( os.path.abspath( args.input ))
    mafft = out_dir + '/' + fasta_name + '.mafft'
    trimal = out_dir + '/' + fasta_name + '.trimal'

    if not os.path.isfile( mafft ):
        mafft = mafftRun( os.path.abspath(args.input), out_dir, hpc )
    else:
        print('\nAlignment exists ...')
    if not os.path.isfile( trimal ):
        trimal = trimalRun( mafft, out_dir, hpc, args.tree )
    else:
        print('\nTrim exists ...' )

    if args.tree == 'fasttree':
        fasttreeRun( trimal , out_dir, hpc )
    elif args.tree == 'iqtree':
        iqtreeRun( trimal, out_dir, hpc )

    if hpc:
        if args.slurm:
            print('\nStart pipeline via `sbatch <STARTSTEP>.sh` in ' + out_dir + '\n')
        else:
            print('\nStart pipeline via `qsub <STARTSTEP>.sh` in ' + out_dir + \
            '\n')

    sys.exit( 0 )
