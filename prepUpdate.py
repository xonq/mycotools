#! /usr/bin/env python3

import subprocess, sys, os, glob, re, argparse, shutil

def collect_files(directory = './', filetype = '*', recursive = False, verbose = False):

    if verbose == True:
        print('\nCompiling files ...')

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


def changeVersion( setup_path, version, execs ):

    with open( setup_path, 'r' ) as raw:
        data = raw.read()

    new_data = re.sub( r'version = .*', 'version = "' + version + '",', data, re.M )
    new_data = re.sub( r'    scripts .+', '    scripts = ' + str(execs) + ',', new_data )
    with open( setup_path, 'w' ) as out:
        out.write( new_data )


def addExecutables( files ):

    execs = []
    for file in files:
        if os.access( file, os.X_OK ):
            dir_ = os.path.basename( os.path.dirname( file ) )
            if dir_ in { 'lib', 'utils', 'db' }:
                execs.append( 'mycotools/' + dir_ + '/' + os.path.basename(file) )
            elif dir_ == 'mycotools':
                execs.append( 'mycotools/' + os.path.basename(file) )

    return execs


if __name__ == '__main__':

    parser = argparse.ArgumentParser( description = 'Prepares mycotools update for pip' )
    parser.add_argument( '-i', '--input', help = 'Input directory' )
    parser.add_argument( '-o', '--output', help = 'Output directory' )
    parser.add_argument( '-v', '--version', help = 'Version update setup.py' )
    args = parser.parse_args()

    if not args.output.endswith( '/' ):
        args.output += '/'
    prevOutput = os.path.dirname( os.path.abspath( args.output ) ) + '/'

    main_scripts = collect_files( args.input, '*' )
    db_scripts = collect_files( args.input + '/db/', '*' )
    tool_scripts = collect_files( args.input + '/lib/', '*' )
    if not os.path.isdir( args.output ):
        os.mkdir( args.output )
        for dir_ in ['lib', 'util', 'db']:
            os.mkdir( args.output + '/' + dir_ )

    with open( args.output + '../include.txt', 'r' ) as raw:
        include = {x.rstrip() for x in raw.read().split('\n')}


    files = []
    files.extend(main_scripts)
    files.extend(db_scripts)
    files.extend(tool_scripts)
    files = [x for x in files if os.path.basename(x) in include]
    for file_ in files:
        dir_ = os.path.basename(os.path.dirname(file_))
        if dir_ != os.path.basename(args.output[:-1]):
            shutil.copy( file_, args.output + dir_ + '/' + os.path.basename(file_))
        else:
            shutil.copy(file_, args.output + os.path.basename(file_))

    if not os.path.isdir( prevOutput + 'archive' ):
        os.mkdir( prevOutput + 'archive/' )
    os.chdir( prevOutput )
    delOld = subprocess.call( ['rm', '-rf', '*dist*', os.path.basename(os.path.abspath(args.output)) + \
        '.egg_info'], stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True )
    execs = addExecutables( files )
    changeVersion( prevOutput + 'setup.py', args.version, execs )

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
            

