#! /usr/bin/env python3

from mycotools.kontools import eprint
import os, sys, re, subprocess, datetime, time, paramiko, signal
import pandas as pd

def timeout( ):
    
    raise Exception( )

# generates arguments file 
def gen_args(arg_list,program_name):

    print('\n\nGenerating', program_name + '.args file')
    date = datetime.datetime.now().strftime('%Y%m%d')
    output_str = ''

# for the first argument, populate a new string with arg
# add a tab for each and a new line for all other
    for arg in arg_list:
        if output_str == '':
            output_str = arg + '\t'
        else:
            output_str += '\n' + arg + '\t'

# if the argument ends with exe, run a which command to
# autopopulate if possible
        if arg[-3:] == 'exe':
            exe_loc_shell = subprocess.run('which ' + arg[:-4],shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            exe_loc_prep = exe_loc_shell.stdout.rstrip()
            exe_loc = exe_loc_prep.decode('utf-8')
            output_str += exe_loc

# output the arguments with the date and program name
    output_name = date + '_' + program_name + '.args'
    with open(output_name,'w') as output_file:
        output_file.write(output_str)

    print('Complete.')



def compile_args(arg_file):

# check for the argument files
    if not os.path.exists(arg_file):
        eprint('\n\nInvalid arguments path. Exit code 999')
        sys.exit(999)

# open argument file and transpose it with headers as arg
# titles
    args = pd.read_csv(arg_file,sep='\t',index_col=0,header=None)
    args_df = args.T

# populate a dictionary with the headers and their
# associated arguments
    args_dict = {}
    for header in args_df:
        args_dict[header] = args.T[header][1]

    return args_dict



def start_check(osc_number, server='owens.osc.edu', username='osu10393', pw='', spacer = '\t'):

# generate command to check job status and grep job number
    check_cmd = "/usr/bin/squeue -u " + username + " | grep " + osc_number
    count = 0
    while True:
        try:
            ssh = paramiko.SSHClient()
            ssh.load_system_host_keys()
            ssh.connect(server, username=username, password=pw)
            (stdin, stdout, stderr) = ssh.exec_command(check_cmd)
            exit_status = stdout.channel.recv_exit_status()
            break
        except OSError:
            time.sleep( 5 )
            count += 1
            if count == 5:
                print('\nERROR:\t5 failed SSH attempts. Sleep 15 minutes.')
                time.sleep( 900 )
                count = 0
# NEED to add error check
    try:
        status = stdout.readlines()[0]
    except IndexError:
        status = '?'
        start_code = '?'

    ssh.close()

# check status of qstat by seeing what letter code it has
# in the output 
    if status != '?':
        if re.search( username + r'\W+PD', status):
            start_code = 'Q'
        elif re.search( username + r'\W+R', status):
            start_code = 'R'
        else:
            start_code = '?'

    return start_code



# runs a simple command on a server and returns the error and out as a list of lines
def servcmd( cmd_str, pwd, user = 'osu10393', server = 'owens.osc.edu' ):

    count = 0
    while True:
        try:
            ssh = paramiko.SSHClient()
            ssh.load_system_host_keys()
            ssh.connect( server, username = user, password = pwd )
            ( stdin, stdout, stderr ) = ssh.exec_command( cmd_str )
            exit_status = stdout.channel.recv_exit_status()
            break
        except OSError:
            time.sleep( 5 )
            count += 1
            if count == 5:
                print('\nERROR:\t5 failed SSH attempts. Sleep 15 minutes.')
                time.sleep( 900 )
                count = 0
    out, err = [], []
    for line in stdout.readlines():
        out.append( line )
    for line in stderr.readlines():
        err.append( line )

    ssh.close()

    return out, err



# need to make user optional for start_check
def execution(
    prog_statement, osc_param_dict, server='owens.osc.edu', username='osu10393', 
    pw='', spacer='\t', start_cmd = 1, interval = 1, verbose = 0
    ):

# populate variables for output and log file
    output_var = osc_param_dict['output_dir']
    log_var = osc_param_dict['log_name']
    count = 0
# if the output directory doesn't exist, exit
# does this NEED to be the log file itself?
# DOESNT WORK BECAUASE IT IS NOT QUERYING THE REMOTE FS
#    if not os.path.exists(output_var):
#        print('\n\n' + spacer + 'Output directory does not exist.\nExit code 420.')
#        sys.exit(420)

# add the program statement to the osc param and make cmd
    osc_param_dict['prog_run'] = prog_statement
    exec_statement = 'echo -e "{prog_run}" | /usr/bin/sbatch --time={walltime} --nodes={nodes} --ntasks-per-node={threads} --output {output_dir}/{log_name} --error {output_dir}/{log_name} --job-name {log_name} -A {project}'.format(**osc_param_dict)
    if verbose == 1:
        print(spacer + exec_statement)


    while True:
        try:
            ssh = paramiko.SSHClient()
            ssh.load_system_host_keys()
            ssh.connect(server, username=username, password=pw)
            (stdin, stdout, stderr) = ssh.exec_command(exec_statement)
            exit_status = stdout.channel.recv_exit_status()
            stdout_read = stdout.readlines()
            ssh.close()
            break
        except OSError:
            time.sleep( 5 )
            count += 1
            if count == 5:
                print('\nERROR:\t5 failed SSH attempts. Sleep 15 minutes.')
                time.sleep( 900 )
                count = 0
# NEED to check for errors here
    for line in stdout_read:
        job_number = re.search(r'Submitted batch job (\d+)', line)
        if job_number:
            job_number = job_number[1]

    log_path = output_var + '/' + log_var

# wait for execution by checking start code
    #print('\t' + spacer + 'Waiting for OSC execution ...')
    time.sleep( interval )
    start_code = start_check(job_number, pw=pw, username=username, server=server, spacer=spacer + '\t')

    if start_cmd == 1:
        time.sleep( 15 )

# while it is queued, wait 60 s, check job again
        while start_code == 'Q':
            time.sleep( 60 )
            start_code = start_check(job_number, remote=remote, pw=pw, username=username, server=server, spacer=spacer + '\t')

# return log file as absolute path
    return log_path, job_number



def status_check(
    output_dir, log_path, pw='', server='owens.osc.edu', completeRE = '',
    username='osu10393', error_exit=1, spacer='\t', port = 22
    ):

    log_file = os.path.basename(log_path)
    exit_code = ''
    count = 0

    while True:
        try:
            signal.signal(signal.SIGALRM, timeout)
            signal.alarm(120)
            if count > 5:
                print( spacer + '5 failed attempts. Resting 10 minutes.' )
                time.sleep( 600 )

            ssh = paramiko.SSHClient()
            ssh.load_system_host_keys()
            ssh.connect(server, username=username, password=pw)
            (stdin, stdout, stderr) = ssh.exec_command('ls ' + os.path.dirname( os.path.abspath (log_path )))
            exit_status = stdout.channel.recv_exit_status()
            break
        except Exception:
            print(spacer + 'SSH timedout (120s)... reattempting')
            count += 1

# NEED to check for errors 
# populate a list of files from the ls cmd stdout
    ls_check = []
    for file_ in stdout.readlines():
        ls_check.append(file_.rstrip())

# if the log file doesnt exist, exit with error
    if log_file not in ls_check:
        eprint( spacer + 'ERROR:\tLog file `' + log_file + '` does not exist.' )
        ssh.close()
        exit_code = 'E'
        if error_exit == 1:
            sys.exit(710)

# grab the log file and save in current dir
    else:
        count = 0
        while True:
            try:
                if count > 5:
                    print( spacer + '5 failed attempts. Resting 10 minutes' )
                    time.sleep( 600 )

                signal.signal(signal.SIGALRM, timeout)
                signal.alarm(120)

                transport = paramiko.Transport(( server, port ))
                transport.connect( None, username, pw )
                sftp = paramiko.SFTPClient.from_transport( transport )
                sftp.get(log_path, './' + log_file)
                if sftp:
                    sftp.close()
                if transport:
                    transport.close()
                ssh.close()
                break
            except Exception:
                print( spacer + 'SFTP timeout (120s)... reattempting' )
                count += 1

    # open log file, save data in variable, remove log_file
        with open(log_file, 'r') as log2import:
            log_data = log2import.read()
        os.remove(log_file)

    # if the log file has resource units charged, it is done
        if re.search( completeRE, log_data):
            exit_code = 'C'

    # useless until updating for slurm
    # if job killed is in file, then exit code is E for error
            error_check = re.search(r'=>> PBS: job killed.*',log_file)
            if error_check:
                eprint(spacer + 'ERROR: OSC early exit see ' + log_path)
                exit_code = 'E'

    # if error_exit is on, exit the program, otherwise go on
                if error_exit == 1:
                    eprint(spacer + 'Exit on error active: code 310')
                    sys.exit(310)
                else:
                    eprint('\t' + spacer + 'Continuing')

    return exit_code
