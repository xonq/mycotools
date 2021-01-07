#! /usr/bin/env python3

import subprocess, sys, os, glob, re, argparse

def collect_files(directory = './', filetype = '*', recursive = True, verbose = False):

    if verbose == True:
        print('\n\nCompiling files ...')

    if type(filetype) == list:
        filetypes = filetype.split()
    else:
        filetypes = [filetype]

    directory = os.path.normpath(os.path.expanduser(directory)) + '/'
    filelist = []
    for filetype in filetypes:
        if recursive == True:
            filelist = filelist + (glob.glob(directory + "/**/*." + filetype, recursive = recursive))
        elif recursive == False:
            filelist = filelist + (glob.glob(directory + "/*." + filetype, recursive = recursive))
        else:
            raise TypeError

    return filelist


def changeVersion( setup_path, version ):

    with open( setup_path, 'r' ) as raw:
        data = raw.read()

    new_data = re.sub( r'version = .*', 'version = "' + version + '",', data, re.M )
    with open( setup_path, 'w' ) as out:
        out.write( new_data )


if __name__ == '__main__':

    parser = argparse.ArgumentParser( description = 'Prepares mycotools update for pip' )
    parser.add_argument( '-i', '--input', help = 'Input directory' )
    parser.add_argument( '-o', '--output', help = 'Output directory' )
    parser.add_argument( '-v', '--version', help = 'Version update setup.py' )
    args = parser.parse_args()

    files = collect_files( args.input, 'py' )
    if not os.path.isdir( args.output ):
        os.mkdir( args.output )

    for file in files:
        subprocess.Popen( ['cp', file, args.output], stdout = subprocess.PIPE, \
            stderr = subprocess.PIPE )

    prevOutput = os.path.dirname( os.path.abspath( args.output ) ) + '/'
    if not os.path.isdir( prevOutput + 'archive' ):
        os.mkdir( prevOutput + 'archive/' )
    os.chdir( prevOutput )
    delOld = subprocess.call( ['rm', '-rf', '*dist*', os.path.basename(os.path.abspath(args.output)) + \
        '.egg_info'], stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True )
    if args.version:
        changeVersion( prevOutput + 'setup.py', args.version )

    setup = subprocess.call(
        'python3 ' + prevOutput + 'setup.py sdist bdist_wheel',
        shell = True
    )
    if setup == 0: 
        upload = subprocess.call(
            'python3 -m twine upload --repository pypi dist/*',
            shell = True
        )

    sys.exit( 0 )
            

