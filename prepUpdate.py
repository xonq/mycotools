#! /usr/bin/env python3

from kontools import collect_files
import subprocess, sys


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

	for file in files:
		subprocess.Popen( ['cp', file, sys.argv[2]], stdout = subprocess.PIPE, \
			stderr = subprocess.PIPE )
