#! /usr/bin/env python3

import subprocess, sys, os, glob, re

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



if __name__ == '__main__':

	usage = 'Input update directory of files and output for package, prepares update'
	if '-h' in sys.argv or '--help' in sys.argv:
		print( '\n' + usage + '\n' )
		sys.exit( 1 )
	elif len( sys.argv ) < 3:
		print( '\n' + usage + '\n' )
		sys.exit( 2 )

	files = collect_files( sys.argv[1], 'py' )
	if not os.path.isdir( sys.argv[2] ):
		os.mkdir( sys.argv[2] )

	names, new_files = [], []
	for file in files:
		subprocess.Popen( ['cp', file, sys.argv[2]], stdout = subprocess.PIPE, \
			stderr = subprocess.PIPE )
		name = os.path.basename( file )
		names.append( name.replace( '.py', '' ) )
		new_files.append( os.path.abspath( sys.argv[2]  ) + '/' + name )

	for file in new_files:
		with open( file, 'r' ) as raw:
			data = raw.read()
		for name in names:
			data1 = re.sub( r'^import(.*?)' + name + r'(.*?)', r'import\1kontools.' + name + r'\2', data )
			data2 = re.sub( r'^from(.*?)' + name + r'(.*?)', r'from\1kontools.' + name + r'\2', data1 )
			if data2 != data:
				with open( file, 'w' ) as out:
					out.write( data2 )

	sys.exit( 0 )
			

