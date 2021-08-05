#! /usr/bin/env python3

import os, re, sys, subprocess, argparse
from mycotools.lib.kontools import eprint, vprint, collect_files, formatPath, intro, outro, findExecs
from mycotools.lib.fastatools import fasta2dict


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
                out.write( hpc + '\n\n' + cmd + '\ncd ' + out_dir + '\nqsub trimal.sh')
            else:
                 out.write( hpc + '\n\n' + cmd + '\ncd ' + out_dir + '\nsbatch trimal.sh')
               
        
    return out_dir + '/' + name + '.mafft'


def trimalRun( mafft, out_dir, hpc, tree, verbose ):

#    name1 = os.path.basename( os.path.abspath(mafft.replace('mafft', 'prelim_trimal')))
    name2 = os.path.basename(os.path.abspath(re.sub( r'\.mafft$', '.trimal', mafft)))
#    cmd1 = 'trimal -in ' + mafft + ' -automated1 -out ' + out_dir + '/' + name1
    cmd2 = 'trimal -in ' + mafft + ' -out ' + out_dir + '/' + name2 + ' -gappyout'

    if not hpc:
#        if not os.path.isfile( out_dir + '/' + name1 ):
 #           vprint('\nTrimming step 1/2', v = verbose, flush = True)
  #          run_trimal1 = subprocess.call( cmd1, shell = True, stdout = subprocess.PIPE )
   #         if run_trimal1 != 0:
    #            os.remove( out_dir + '/' + name1 )
     #           eprint( '\nERROR: preliminary `trimal` failed\n' , flush = True)
      #          sys.exit( 3 )

        vprint('\nTrimming step 2/2', v = verbose, flush = True)
        run_trimal2 = subprocess.call( cmd2, shell = True, stdout = subprocess.PIPE )
        if run_trimal2 != 0:
            eprint( '\nERROR: `trimal` failed\n' , flush = True)
            os.remove( out_dir + '/' + name2 )
            sys.exit( 4 )

    else:
        with open( out_dir + '/trimal.sh', 'w' ) as out:
            if '#PBS' in hpc:
                out.write( hpc + '\n\n' + #cmd1 + '\n' + 
                    cmd2 + '\n\ncd ' + out_dir + '\nqsub ' + tree + '.sh' )
            else:
                out.write( hpc + '\n\n' + #cmd1 + '\n' + 
                    cmd2 + '\n\ncd ' + out_dir + '\nsbatch ' + tree + '.sh' )
    
    return out_dir + '/' + name2


def fasttreeRun( trimal, out_dir, hpc, verbose ):


    name = os.path.basename(os.path.abspath(trimal + '.fast.tre'))
    cmd = 'fasttree -lg ' + trimal + ' > ' + out_dir + '/' + name

    if not hpc:
        vprint('\nMaking fasttree ...', v = verbose, flush = True)
        run_fasttree = subprocess.call( cmd, shell = True, stdout = subprocess.PIPE )
        if run_fasttree != 0:
            eprint( '\nERROR: `fasttree` failed\n' , flush = True)
            os.remove( out_dir + '/' + name )
            sys.exit( 5 )
    else:
        with open( out_dir + '/fasttree.sh', 'w' ) as out:
            out.write( hpc + '\n\n' + cmd )

    return name


def iqtreeRun( trimal, out_dir, hpc, verbose ):

    cmd = 'iqtree -s ' + trimal + ' -B 1000 -T AUTO'

    vprint('\nOutputting bash script `iqtree.sh`.\n', v = verbose, flush = True)
   
    if not hpc:
        run_iqtree = subprocess.call(cmd.split(' '))
        if run_iqtree != 0:
            eprint('\nERROR: iqtree failed\n', flush = True)
            sys.exit(6)
    else: 
        with open( out_dir + '/iqtree.sh', 'w' ) as out:
            out.write( hpc + '\n\n' + cmd )


def main( fasta_path, tree, slurm = False, pbs = False, project = '', output_dir = None, verbose = True, alignment = False ):

    hpc = False
    if slurm:
        vprint('\nHPC mode, preparing submission scripts (Slurm)\n', v = verbose, flush = True)
        hpc = '#!/bin/bash\n#SBATCH --time=10:00:00\n#SBATCH --nodes=1\n' + \
            '#SBATCH --ntasks-per-node=4\n#SBATCH -A ' + project + '\n' + \
            '#SBATCH --mem=40000'
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
    trimal = out_dir + '/' + fasta_name + '.trimal'

    if not os.path.isfile( mafft ) and not alignment:
        mafft = mafftRun( os.path.abspath(fasta_path), out_dir, hpc, verbose )
    elif alignment:
        mafft = fasta_path
    else:
        vprint('\nAlignment exists ...', v = verbose, flush = True)
    if not os.path.isfile( trimal ):
        trimal = trimalRun( mafft, out_dir, hpc, tree, verbose )
    else:
        vprint('\nTrim exists ...', v = verbose , flush = True)

    if tree == 'fasttree':
        fasttreeRun( trimal , out_dir, hpc, verbose )
    elif tree == 'iqtree':
        iqtreeRun( trimal, out_dir, hpc, verbose )

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
        ' aligns, trims, and runs treebuilding: {`fasttree`, `iqtree`, `raxml` (future)}.' )
    parser.add_argument( '-f', '--fasta', \
        help = 'Fasta or directory of fastas' )
    parser.add_argument( '-a', '--alignment', \
        help = 'Alignment or directory of alignments (need fasta file extension)' )
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
        'Input': args.fasta, 'Tree': args.tree, 'Torque': args.pbs,
        'Slurm': args.slurm, 'HPC project': args.project, 'Output': output
        }
    start_time = intro( 'Fasta2Tree', args_dict )

    if not args.fasta and not args.alignment:
        eprint('\nERROR: no input', flush = True)
        sys.exit(7)

    if args.tree not in { 'fasttree', 'iqtree' }:
        eprint('\nERROR: invalid tree software: ' + args.tree , flush = True)
        sys.exit(1)

    findExecs( [args.tree, 'mafft', 'trimal'] )

    if args.fasta:
        input_check = formatPath(args.fasta)
        alignment = False
    else:
        input_check = formatPath(args.alignment)
        alignment = True
    if os.path.isfile(input_check):
        main( 
            input_check, args.tree, slurm = args.slurm, 
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
                fa, args.tree, slurm = args.slurm, pbs = args.pbs, alignment = alignment,
                project = args.project, output_dir = output, verbose = False
                )
    else:
        eprint('\nERROR: invalid input', flush = True)
        sys.exit(3)

    outro( start_time )
