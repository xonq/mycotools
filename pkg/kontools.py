#! /usr/bin/env python3

import datetime, sys, glob, os, re, subprocess, numpy as np

# prints error to std error
def eprint( *args, **kwargs ):
    print(*args, file = sys.stderr, **kwargs )

def vprint( toPrint, v = False ):

    if v:
        print( toPrint )


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

def collect_folders( directory, path = True ):

	get_folders = subprocess.run( 'ls ' + directory, shell = True, stdout = subprocess.PIPE )
	folders = get_folders.stdout.decode( 'UTF-8' ).split( '\n' )
	if path:
		folders = [ os.path.abspath( directory ) + '/' + folder for folder in folders if folder != '' ]
	else:
		folders = [ folder for folder in folders if folder != '' ]

	folders = [ folder for folder in folders if os.path.isdir(folder) ]
	return folders


def dictSplit( Dict, factor ):

	list_dict = []
	keys_list = np.array_split( list(Dict.keys()), factor )
	
	for keys in keys_list:
		list_dict.append( {} )
		for key in keys:
			list_dict[-1][key] = Dict[key]

	return list_dict
	
	

def intro(script_name,args_dict,credit='', log = False):

    start_time = datetime.datetime.now()
    date = start_time.strftime( '%Y%m%d' )

    out_str = '\n' + script_name + '\n' + credit + '\nExecution began: ' + str(start_time)

    for arg in args_dict:
        out_str += '\n' + '{:<30}'.format(arg + ':') + str(args_dict[arg])

    print( out_str )
    if log:
        out_file = os.path.abspath( log ) + '/' + date + '.log'
        count = 0
        while os.path.exists( out_file ):
            count += 1
            out_file = os.path.abspath( log ) + '/' + date + '_' + str(count) + '.log'

        with open(out_file, 'w') as out:
            out.write( out_str )

    return start_time


def outro(start_time, log = False):

    end_time = datetime.datetime.now()
    duration = end_time - start_time
    dur_min = duration.seconds/60
    print('\nExecution finished:', str(end_time) + '\t' + '\n\t{:.2}'.format(dur_min), 'minutes\n')

#    out_str = '\n\nExecution finished:', str(end_time) + '\t' + '\n\t{:.2}'.format(dur_min), 'minutes'

    sys.exit(0)


def prep_output(output, mkdir=1, require_newdir=0, cd=1):

    output = os.path.realpath(output)
    
    if os.path.isdir(output):
        if require_newdir == 1:
            print('\n\nERROR: directory exists. Exit code 172')
            sys.exit( 172 )
    elif os.path.exists(output):
        output = os.path.dirname(output)
    else:
        if mkdir == 0:
            print('\n\nERROR: directory does not exist. Exit code 173')
            sys.exit( 173 )
        os.mkdir(output)

    if cd == 1:
        os.chdir(output)

    return output


def file2list(file_, sep = '\n', col = None):

    if sep not in ['\n', '\t', ',', ' ']:
        print('\n\nERROR: invalid separator.\nExit code 229')
        sys.exit( 229 )

    with open(file_, 'r') as raw_data:
        data = raw_data.read() + '\n'

    if col:
        check = data.split('\n')
        check_col = check[0].split(sep).index( col )
        data_list = [ x[check_col].rstrip() for x in [ y.split(sep) for y in check[1:] ] if len(x) >= check_col + 1 ]
    else:
        data_list = data.split( sep = sep )
        data_list = [ x.rstrip() for x in data_list if x != '' ]

    return data_list


def check_cmd( cmd ):

    exit_code = 0
    cmd_check = subprocess.call( 'which ' + str(cmd), shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE )
    if cmd_check != 0:
        exit_code = 1

    return exit_code 


def multisub( args_lists, processes = 1, timeout = None, shell = False ):

    running, outputs = [ ], []
    if processes <= 0:
        processes = 1
    while len( args_lists ) != 0:
        while len( running ) < processes and len( args_lists ) != 0:
            cmd = args_lists[ 0 ]
            run_temp = subprocess.Popen( cmd, stdout = subprocess.DEVNULL, \
                stderr = subprocess.DEVNULL, shell = shell )
            running.append( [run_temp, cmd ] )
            del args_lists[ 0 ]

        if len( args_lists ) > 0:
            for index in range( len( running ) ):
                handle = running[ index ]
                handle[0].poll()
                returncode = handle[0].returncode
                if returncode is not None:
                    outputs.append( {
                        'stdin': handle[1], 'code': returncode
                        } )
                    del running[ index ]
                    break

        else:
            while len( running ) > 0:
                handle = running[ 0 ]
                handle[0].poll()
                returncode = handle[0].returncode
                if returncode is not None:
                    outputs.append( {
                        'stdin': handle[1], 'code': returncode
                        } )
                    del running[ 0 ]

    return outputs  
