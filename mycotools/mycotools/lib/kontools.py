#! /usr/bin/env python3

from datetime import datetime
import sys, glob, os, re, subprocess, gzip, shutil, json, tarfile
import numpy as np


def tardir(dir_, rm = True):
    
    if not os.path.isdir(formatPath(dir_)):
        return False
    with tarfile.open(formatPath(dir_)[:-1] + '.tar.gz', 'w:gz') as tar:
        tar.add(dir_, arcname = os.path.basename(formatPath(dir_)[:-1]))
    shutil.rmtree(dir_)

def untardir(dir_, rm = False, to = None):
    if not to:
        to = os.path.dirname(dir_[:-1])
    tar = tarfile.TarFile.open(dir_)
    tar.extractall(path = to)
    tar.close()
    if rm:
        os.remove(dir_)

def checkdir(dir_, unzip = False, to = None, rm = False):
    if os.path.isdir(dir_):
        return True
    elif os.path.isfile(formatPath(dir_) + '.tar.gz'):
        if unzip:
            if dir_.endswith('/'):
                dir_ = dir_[:-1]
            untardir(dir_ + '.tar.gz', rm = rm, to = to)
            return True
    return False


def eprint( *args, **kwargs ):
    '''Prints to stderr'''

    print(*args, file = sys.stderr, **kwargs)


def fprint(out_str, log):
    with open(log, 'a') as out:
        out.write(args)

def zprint(out_str, log = None, flush = True):
    fprint(out_str, log)
    print(out_str, flush = flush)

def vprint( toPrint, v = False, e = False , flush = True):
    '''Boolean print option to stdout or stderr (e)'''

    if v:
        if e:
            eprint( toPrint, flush = True)
        else:
            print( toPrint, flush = True)


def readJson( config_path ):

    with open( config_path, 'r' ) as json_raw:
        json_dict = json.load( json_raw )

    return json_dict
   

def writeJson(json_path, obj, **kwargs):
    with open(json_path, 'w') as json_out:
        json.dump(obj, json_out, **kwargs)


def gunzip( gzip_file, remove = True, spacer = '\t' ):
    '''gunzips gzip_file and removes if successful'''

    new_file = re.sub( r'\.gz$', '', formatPath(gzip_file) )
    try:
        with gzip.open(formatPath(gzip_file), 'rb') as f_in:
            with open(new_file, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        if remove:
            os.remove( gzip_file )
        return new_file
    except:
        eprint('\n' + spacer + 'ERROR: gunzip ' + str(gzip_file) + ' failed.', flush = True)
        if os.path.isfile( new_file ):
            if os.path.isfile( gzip_file ):
                os.remove( new_file )
        return False


def fmt_float(val, sig_dig = None):

    val_str = str(val)
    try:
        e = val_str.index('e')
        dig = val_str[:e].replace('.','')
        num = val_str[e + 1:]
        if num.startswith('-'):
            val_op = '0.'
            for i in range(abs(int(num))-1):
                val_op += '0'
            val_str = val_op + dig.replace('.','')
        else:
            num = abs(int(num))
            if len(dig) > num:
                val_op = list(dig) 
                val_op.insert('.',num+1)
                val_str = ''.join(val_op)
            else:
                zeroes = num + 1 - len(dig)
                val_op = dig
                for i in range(zeroes):
                    val_op += '0'
                val_str = val_op
    except ValueError:
        val_str = str(val)

    if sig_dig:
        if len(val_str.replace('.','')) > sig_dig:
            try:
                per_i = val_str.index('.')
            except ValueError:
                per_i = None
            val_list = list(val_str)
            if per_i:
                if sig_dig > per_i:
                    i = val_list[sig_dig]
                    post = int(val_list[sig_dig+1])
                    if post >= 5:
                        val_list.insert(sig_dig, int(i) + 1)
                        val_list.pop(sig_dig + 1)
                    val_str = ''.join(str(v) for v in val_list[:sig_dig + 1])
                else:
                    i = val_list[sig_dig - 1]
                    try:
                        post = val_list[sig_dig]
                    except IndexError:
                        post = 0
                    if post == '.':
                        post = val_list[sig_dig + 1]
                    elif i == '.':
                        i = val_list[sig_dig]
                        try:
                            post = val_list[sig_dig + 1]
                        except IndexError:
                            post = 0
                    if int(post) >= 5:
                        val_list.insert(sig_dig, int(i) + 1)
                        val_list.pop(sig_dig + 1)
                    val_list = val_list[:sig_dig]
                    while len(val_str) > len(val_list):
                        val_list.append('0')
                    val_str = ''.join(str(v) for v in val_list[:sig_dig+3])
            else:
                i = val_list[sig_dig + 1]
                post = int(val_list[sig_dig + 2])
                if post == '.':
                    post = int(val_list[sig_dig + 3])
                elif i == '.':
                    i = val_list[sig_dig + 2]
                    try:
                        post = val_list[sig_dig + 3]
                    except IndexError:
                        post = 0
                if int(post) >= 5:
                    val_list.insert(sig_dig + 1, int(i) + 1)
                    val_list.pop(sig_dig + 2)
                val_list = val_list[:sig_dig + 2]
                while len(val_str) > len(val_list):
                    val_list.append('0')
                val_str = ''.join(str(v) for v in val_list[:sig_dig + 1])
                

    return val_str




def findExecs( deps, exit = set(), verbose = True ):
    '''
    Inputs list of dependencies, `dep`, to check path.
    If dependency is in exit and dependency is not in path,
    then exit.
    '''

    vprint('\nDependency check:', v = verbose, e = True, flush = True)
    checks = []
    if type(deps) is str:
        deps = [deps]
    for dep in deps:
        check = shutil.which( dep )
        vprint('{:<15}'.format(dep + ':', flush = True) + \
            str(check), v = verbose, e = True)
        if not check and dep in exit:
            eprint('\nERROR: ' + dep + ' not in PATH', flush = True)
            sys.exit(300)
        else:
            checks.append(check)

    return checks


def findEnvs( envs, exit = set(), verbose = True ):
    '''
    Inputs list of paths, `envs`, to check path.
    If env is not in path and it is in exit, exit.
    '''

    vprint('\nEnvironment check:', v = verbose, e = True, flush = True)
    if type(envs) is str:
        envs = [envs]
    eprint(flush = True)
    for env in envs:
        try:
            vprint('{:<15}'.format(env + ':', flush = True) + \
                str(os.environ[env]), v = verbose, e = True)
        except KeyError:
            vprint('{:<15}'.format(env + ':', flush = True) + \
                'None', v = verbose, e = True)
            if env in exit:
                eprint('\nERROR: ' + env + ' not in PATH', flush = True)
                sys.exit(301)


def expandEnvVar( path ):
    '''Expands environment variables by regex substitution'''

    if path.startswith( '/$' ):
        path = '/' + path
    env_comp = re.compile( r'/\$([^/]+)' )
    if not env_comp.search( path ):
        env_comp = re.compile( r'^\$([^/]+)' )
    var_search = env_comp.search( path )
    if var_search:
        var = var_search[1]
        pathChange = os.environ[ var ]
        path = env_comp.sub( pathChange, path )

    return path


def formatPath( path, isdir = None ):
    '''Goal is to convert all path types to absolute path with explicit dirs'''
   
    if path:
        path = os.path.expanduser( path )
        path = expandEnvVar( path )
        path = os.path.abspath( path )
        if isdir:
            if not path.endswith('/'):
                path += '/'
        else:
            if path.endswith( '/' ):
                if path.endswith( '//' ):
                    path = path[:-1]
                if not os.path.isdir( path ):
                    path = path[:-1]
            elif not path.endswith( '/' ):
                if os.path.isdir( path ):
                    path += '/'

    return path


# need to change recursive to false
def collect_files( directory = './', filetype = '*', recursive = False ):
    '''
    Inputs: directory path, file extension (no "."), recursivity bool
    Outputs: list of files with `filetype`
    If the filetype is a list, split it, else make a list of the one entry.
    Parse the environment variable if applicable. Then, obtain a clean, full 
    version of the input directory. Glob to obtain the filelist for each 
    filetype based on whether or not it is recursive.
    '''

    if type(filetype) == list:
        filetypes = filetype.split()
    else:
        filetypes = [filetype]

    directory = formatPath( directory )
    filelist = []
    for filetype in filetypes:
        if recursive:
            filelist.extend( 
                glob.glob( directory + "/**/*." + filetype, recursive = recursive )
            )
        else:
            filelist.extend(
                glob.glob( directory + "/*." + filetype, recursive = recursive ) 
            )

    return filelist


def collect_folders( directory, recursive = False ):
    '''
    Inputs: directory path, recursive search boolean
    Outputs: list of folders
    Get folders via os, recursively extend output list if True.
    '''

    directory = formatPath( directory )
    folders_prep = os.listdir( directory )
    folders = [ 
        directory + '/' + folder + '/' for folder in folders_prep \
        if os.path.isdir( directory + '/' + folder ) 
    ]

    if recursive:
        for folder in folders:
            folders.extend(collect_folders( folder, recursive ))

    return folders


def dictSplit( Dict, factor ):
    '''
    Inputs: a dictionary `Dict`, and an integer `factor` to split by
    Outputs: a list of split dictionaries `list_dict`
    Split dictionary keys into a numpy array via the factor. For each
    key, append a new dictionary to the output list and populate it 
    with the list of keys within that factor. The last entry will be
    whatever amount is leftover.
    '''
    

    list_dict = []
    keys_list = np.array_split( list(Dict.keys()), factor )

    try:    
        for keys in keys_list:
            list_dict.append( {} )
            for key in keys:
                list_dict[-1][key] = Dict[key]
    except TypeError:
        list_dict = []
        for keys in keys_list:
            list_dict.append( {} )
            for key in keys:
                list_dict[-1][tuple(key)] = Dict[tuple(key)]

    return list_dict


def sysStart( args, usage, min_len, dirs = [], files = [] ):

    if '-h' in args or '--help' in args:
        print( '\n' + usage + '\n' , flush = True)
        sys.exit( 1 )
    elif len( args ) < min_len:
        print( '\n' + usage + '\n' , flush = True)
        sys.exit( 2 )
    elif not all( os.path.isfile( formatPath(x) ) for x in files ):
        print( '\n' + usage , flush = True)
        eprint( 'ERROR: input file(s) do not exist\n' , flush = True)
        sys.exit( 3 )
    elif not all( os.path.isfile( formatPath(x) ) for x in dirs ):
        print( '\n' + usage , flush = True)
        eprint( 'ERROR: input directory does not exist\n' , flush = True)
        sys.exit( 4 )

    return args

 
def intro( script_name, args_dict, credit='', log = False, stdout = True):
    '''
    Inputs: script_name string, args_dict dictionary of arguments, 
    credit string bool / path for output log
    Outputs: prints an introduction, returns start_time in YYYYmmdd format
    Creates a string to populate and format for the introduction using 
    keys as the left-most descriptor and arguments (values) as the right
    most. Optionally outputs a log according to `log` path.
    '''

    start_time = datetime.now()
    date = start_time.strftime( '%Y%m%d' )

    out_str = '\n' + script_name + '\n' + credit + \
        '\nExecution began: ' + str(start_time)

    for arg in args_dict:
        out_str += '\n' + '{:<30}'.format(arg.upper() + ':') + \
            str(args_dict[ arg ])

    if log:
        zprint(out_str, log)
    elif stdout:
        print(out_str, flush = True)
    else:
        eprint( out_str, flush = True )

    return start_time


def outro( start_time, log = False, stdout = True ):
    '''
    Inputs: start time string formatted YYYYmmdd, log path
    Outputs: prints execution time and exits with 0 status
    '''

    
    end_time = datetime.now()
    duration = end_time - start_time
    dur_min = duration.seconds/60
    out_str = '\nExecution finished: ' + str(end_time) + '\t' + \
            '\n\t{:.2}'.format(dur_min) + ' minutes\n'

    if log:
        zprint(out_str, log)
    elif not stdout:
        eprint(out_str, flush = True)
    else:
        print(out_str, flush = True)

    sys.exit(0)


def prep_output(output, mkdir = True, require_newdir = False, cd = False):
    '''
    Inputs: output path, bool `mkdir`, bool `require_newdir`, bool cd.
    Outputs: creates directory, returns None if bool inputs fail, returns 
    formatted path
    '''

    output = formatPath( output )
    if os.path.isdir( output ):
        if require_newdir:
            eprint('\nERROR: directory exists.', flush = True)
            return None
    elif os.path.exists( output ):
        output = os.path.dirname(output)
    else:
        if not mkdir:
            eprint('\nERROR: directory does not exist.', flush = True)
            return None
        os.mkdir(output)

    if cd:
        os.chdir(output)

    return output

def mkOutput(base_dir, program, reuse = True, suffix = datetime.now().strftime('%Y%m%d')):
    if not os.path.isdir(formatPath(base_dir)):
        raise FileNotFoundError(base_dir + ' does not exist')
    out_dir = formatPath(base_dir) + program + '_' + suffix

    if not reuse:
        count, count_dir = 1, out_dir
        while os.path.isdir(count_dir):
            count_dir  += '_' + str(count)
            count += 1
        os.mkdir(count_dir)
        return count_dir + '/'
    else:
        if not os.path.isdir(out_dir):
            os.mkdir(out_dir)
        return out_dir + '/'


def checkDep( dep_list = [], var_list = [], exempt = set() ):
    '''Checks all dependencies in path from list, optional exemption set'''

    failedVars = []
    check = [ shutil.which( dep ) for dep in dep_list if dep not in exempt ]
    for var in var_list:
        if var not in exempt:
            try:
                os.environ[ var ]
                check.append( True )
            except KeyError:
                check.append( False )
                failedVars.append( var )
    if not all( check ):
        eprint('\nERROR: Dependencies not met:', flush = True)
        for dep in dep_list:
            if not shutil.which( dep ):
                eprint( dep + ' not in PATH' , flush = True)
        for failed in failedVars:
            eprint( failed + ' variable not set' , flush = True)
        sys.exit( 135 )


def file2list( file_, types = '', sep = '\n', col = None, compress = False):
    '''
    Inputs: `file_` path, separator string, column integer index
    Outputs: reads file and returns a list given arguments
    If the `sep` is not valid, return None, else read the file and
    split each line, if `col`, split each by the delimiter, and
    acquire the column index at each line. If not `col`, then 
    simply split by the separator and grab the only entry.
    '''

    if not compress:
        with open(file_, 'r') as raw_data:
            data = raw_data.read() + '\n'
    else:
        with gzip.open(file_, 'rt') as raw_data:
            data = raw_data.read() + '\n'

    if col:
        check = data.split('\n')
        check_col = check[0].split(sep).index( col )
        data1_list = [ 
            x[check_col].rstrip() for x in [ y.split(sep) for y in check[1:] ] if len(x) >= check_col + 1 
        ]
    elif types:
        data_list = data.split(sep = sep)
        if types == 'int':
            try:
                data1_list = [int(x.rstrip()) for x in data_list if x]
            except ValueError:
                data1_list = [int(float(x.rstrip())) for x in data_list if x]
        elif types == 'float':
            data_list = [float(x.rstrip()) for x in data_list if x]
    else:
        data_list = data.split( sep = sep )
        data1_list = [ x.rstrip() for x in data_list if x ]

    return data1_list


def multisub( args_lists, processes = 1, shell = False, stdout = None, rm = False ):
    '''
    Inputs: list of arguments, integer of processes, subprocess `shell` bool
    Outputs: launches and monitors subprocesses and returns list of exit 
    information
    If processes are too few, opt to default. While there are entries in the
    `args_lists`: 
        while there are less entries in `running` than max processes
    and arguments remaining in the `args_lists`: Popen the first argument in 
    the list, append the subprocess information and command to `running` and
    delete the argument from `args_lists`. 
        If there are still args in the args_lists, for each index in running
        grab the handle, poll if it is done, and if it is done output the 
        return code and stdin information to outputs and delete the index
        from `running`.
        IF there are not arguments in args_lists, wait until all processes
        are complete before exiting the function.
    '''

    running, outputs = [ ], []
    if processes <= 0:
        processes = 1

    if not stdout:
        while args_lists:
            while len( running ) < processes and args_lists:
                cmd = args_lists[ 0 ]
                run_temp = subprocess.Popen( cmd, stdout = subprocess.DEVNULL, \
                    stderr = subprocess.DEVNULL, shell = shell )
                running.append( [run_temp, cmd ] )
                del args_lists[ 0 ]

            if len( args_lists ) > 0:
                for index, handle in enumerate(running):
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

    else:
        count, outs = 0, []
        while args_lists:
            while len( running ) < processes and args_lists:
                cmd = args_lists[ 0 ][0]
                running.append([
                    cmd, args_lists[0][1], 
                    open(args_lists[0][1], 'wb')
                    ])
                run_temp = subprocess.Popen( cmd, stdout = running[-1][2], \
                    stderr = subprocess.DEVNULL, shell = shell )
                running[-1].append(run_temp)
                if rm:
                    running[-1].append(args_lists[0][2])
                del args_lists[ 0 ]

            if len( args_lists ) > 0:
                for index, handle in enumerate(running):
                    handle[3].poll()
                    returncode = handle[3].returncode
                    if returncode is not None:
#                        output = handle[0].stdout
 #                       with open(handle[2], 'wb') as out:
  #                          for line in output:
   #                             out.write(line)
                        if rm:
                            os.remove(handle[4])
                        handle[2].close()
                        outputs.append( {
                            'stdin': handle[0], 'code': returncode,
                            'output': handle[1]
                            } )
                        del running[ index ]
                        break

            else:
                while len( running ) > 0:
                    handle = running[ 0 ]
                    handle[3].poll()
                    returncode = handle[3].returncode
                    if returncode is not None:
#                        output = handle[0].stdout
 #                       with open(handle[2], 'wb') as out:
  #                          for line in output:
   #                             out.write(line)
                        handle[2].close()
                        if rm:
                            os.remove(handle[4])
    #                        handle[2] = handle[2][0]
                        outputs.append( {
                            'stdin': handle[0], 'code': returncode,
                            'output': handle[1]
                            } )
                        del running[ 0 ]

    return outputs  
