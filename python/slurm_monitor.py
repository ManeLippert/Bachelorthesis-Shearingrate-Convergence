#!/usr/bin/env python3

# AUTHOR ===================================================================================================================
# Manuel Lippert (GitHub: ManeLippert (https://github.com/ManeLippert))
# ==========================================================================================================================

# MODULES ==================================================================================================================

import datetime, time, os, sys, subprocess, math, argparse

# PARSER ===================================================================================================================

description_text = '''

===================================== DESCRIPTION =====================================

FEATURES:
o Creats jobscript file with defined content (look into file itself for that)
o Start/Restarts job until criteria is suffused (default=0)
o Makes backup after each run before Restart amd Restore files after fail
o Sends mail at the beginning, end and restart (default=False) with status file 
  as attachment (mailx has to be defined)
o Creates status file with current status and progress bars
  Progress: Total progress of simulation
  Run X:    Progress of current run

START SCRIPT IN BACKGROUND:
>>> nohup python3 -u slurm_monitor.py &> status.txt &

With arguments (example for 30000 timesteps):
>>> nohup python3 -u slurm_monitor.py -n 30000 &> status.txt &

Output:
>>> [1] 10537
This will write every output in the status file 'status.txt'

LIST PROCESS:
>>> ps ax | grep slurm_monitor.py

Output:
>>> 10537 pts/1    S      0:00 python3 -u slurm_monitor.py
>>> 23426 pts/1    S+     0:00 grep --color=auto slurm_monitor.py

KILL PROCESS:
>>> kill 10537

OUTPUT STATUS:
>>> cd $DATA && find . -name status.txt -exec tail --lines=+2 {} \;

====================================== ARGUMENTS ======================================
'''

parser = argparse.ArgumentParser(description=description_text, formatter_class=argparse.RawTextHelpFormatter)
#parser._action_groups.pop()

required = parser.add_argument_group('required arguments')
additional = parser.add_argument_group('additional arguments')

required.add_argument('-n', dest='timesteps', nargs='?', type=int, required=True,
                    help='required timesteps                   (REQUIRED)')

additional.add_argument('--mail', dest='mail', nargs='?', type=str,
                    help='mail address (mail@server.de)        (default=None)')

additional.add_argument('--restart-mail', dest='bool', nargs='?', type=bool, default=False,
                    help='mail after every restart             (default=False)')

additional.add_argument('-b', '--backup', dest='backup', nargs='?', type=str,
                    help='backup location for files            (default=None)')

additional.add_argument('-m', '--monitorfile', dest='monitorFile', nargs='?', type=str, default='status.txt',
                    help='file with output from nohup command  (default=status.txt)')

additional.add_argument('-r', '--restartfile', dest='restartFile', nargs='?', type=str, default='FDS.dat',
                    help='restart file with data               (default=FDS.dat)')

additional.add_argument('-j', '--job', dest='jobscriptFile', nargs='?', type=str, default='jobscript',
                    help='jobscript to run SLURM job           (default=jobscript-create)')

additional.add_argument('--job-name', dest='jobname', nargs='?', type=str, default='Name',
                    help='job name not longer than 8 character (default=Name)')

additional.add_argument('--ntask-per-node', dest='tasks', nargs='?', type=str, default='32',
                    help='MPI task per node                    (default=32)')

additional.add_argument('--nodes', dest='nodes', nargs='?', type=str, default='3',
                    help='number of nodes                      (default=3)')

additional.add_argument('--time', dest='walltime', nargs='?', type=str, default='0-24:00:00',
                    help='walltime of server (d-hh:mm:ss)      (default=0-24:00:00)')

args = parser.parse_args()

# VARIABLES ================================================================================================================

outputCriteria = '0'

sleepTime = 5

runCounter = 0
currentTime = '00:00:00'

# PARSER VARIABLES =========================================================================================================

emailAddress = args.mail
restartMail = args.bool

backupLocation = args.backup

jobName = args.jobname
tasks = args.tasks
nodes = args.nodes
walltime = args.walltime

nTimestepsRequired = args.timesteps

jobscriptFilename = args.jobscriptFile
restartFilename = args.restartFile

monitorFilename = args.monitorFile

def outputFilename(info):
    return './slurm-' + info + '.out'

## JOBSCRIPT CONTENT =======================================================================================================

jobscriptContent = '''#!/bin/bash -l

# jobname
#SBATCH --job-name=''' + jobName + '''

# MPI tasks
#SBATCH --ntasks-per-node=''' + tasks + '''

# number of nodes
#SBATCH --nodes=''' + nodes + '''

# walltime
#              d-hh:mm:ss
#SBATCH --time=''' + walltime + '''

# execute the job
time mpirun -np $SLURM_NTASKS ./gkw.x > output.dat

# end
exit 0
'''

jobName = jobName[0:8]

## FLAGS ===================================================================================================================

startOutputFlag = 'Submitted batch job '
restartFlag = 'FILE_COUNT'

## COMMANDS ================================================================================================================

commandJobRunning = 'squeue --states=running -u '
commandJobPending = 'squeue --states=pending -u '

commandJobStarting = 'sbatch'

# SWITCHES =================================================================================================================

if emailAddress==None:
    emailNotification = False
else:
    emailNotification = True

if backupLocation==None:
    backup = False
else:
    backup = True

# FUNCTIONS ================================================================================================================

## PROGRESSBAR =============================================================================================================

def progressbar(required_value, current_value, barsize=42,
                prefix='', 
                progress_fill='=', progress_fill_top='', progress_fill_bot='',
                progress_unfill='.', 
                progress_bracket=['[',']']):
    
    x = int(barsize*current_value/required_value)
    percent = int(100*current_value/required_value)
    
    try:
        percent_format = '  (' + (int(math.log10(100)) - int(math.log10(percent)))*' ' + '{}%)'
    except ValueError:
        percent_format = '  (' + int(math.log10(100))*' ' + '{}%)'
        
    try:
        ratio_format = '  ' + (int(math.log10(required_value)) - int(math.log10(current_value)))*' ' + '{}/{}'
    except ValueError:
        ratio_format = '  ' + int(math.log10(required_value))*' ' + '{}/{}'

        
    bar_format =  '{}  ' + progress_bracket[0] + progress_fill_bot + '{}' + progress_fill_top + '{}' + progress_bracket[1] + percent_format + ratio_format
    bar = bar_format.format(prefix, progress_fill*x, progress_unfill*(barsize-x), percent, current_value, required_value)
        
    return bar

## OUTPUT TABLE ============================================================================================================

def print_table_row(content,
                    current_value, required_value,
                    run_conter, current_time, required_time = walltime,
                    delete_line_index = -5,
                    output_type = None, time_info = True,
                    table_outline = ['╭─', '─╮', '╰─', '─╯', '├─', '─┤', '│ ', ' │', '─'],):
    
    # for better format across plattforms
    # table_outline = ['o-', '-o', 'o-', '-o', 'o-', '-o', '| ', ' |', '-']
    
    delete_line_in_file(monitorFilename, end=delete_line_index)
    
    sep_top = table_outline[0] + 76*table_outline[8] + table_outline[1]
    sep_mid = table_outline[4] + 76*table_outline[8] + table_outline[5]
    sep_end = table_outline[2] + 76*table_outline[8] + table_outline[3]
    
    row_cols = [2, 10, 29, 13, 11, 13, 2]
    row_format = ''.join(['{:<' + str(col) + '}' for col in row_cols])
    
    progress_cols = [2, 76, 2]
    progress_format = ''.join(['{:<' + str(col) + '}' for col in progress_cols])
    
    progressbar_content = [progressbar(required_value, current_value, progress_fill_top='>',
                                       prefix='PROGRESS')]
    progressbar_content.insert(0, table_outline[6])
    progressbar_content.insert(len(progressbar_content), table_outline[7])
    
    required_time = get_time_in_seconds(required_time)
    current_time = get_time_in_seconds(current_time)
    
    try:
        progressbartime_content = [progressbar(required_time, current_time, progress_fill_top='>',
                                               prefix='RUN ' + str(run_conter) + (2 - int(math.log10(run_conter)))*' ' + ' ')]
    except ValueError:
        progressbartime_content = [progressbar(required_time, current_time, progress_fill_top='>',
                                               prefix='RUN ' + str(run_conter) + 2*' ' + ' ')]
    progressbartime_content.insert(0, table_outline[6])
    progressbartime_content.insert(len(progressbartime_content), table_outline[7])
    
    content.insert(0, table_outline[6])
    if time_info:
        if output_type == 'header':
            content.append('DATE')
            content.append('TIME')
            content.append('W:DD:HH:MM:SS')
        else:
            content.append(time_date())
            content.append(time_time())
            content.append(time_duration(startTime))
            
    content.insert(len(content), table_outline[7])
    
    if output_type == 'header':
        sys.stdout.write("\x1b[1A"*(1))
        
        write_line_to_monitorfile(sep_top)
        write_line_to_monitorfile(row_format.format(*content))
        write_line_to_monitorfile(sep_end)
        
    elif output_type == 'middle':
        sys.stdout.write("\x1b[1A"*(-delete_line_index))
        
        write_line_to_monitorfile(sep_mid)
        write_line_to_monitorfile(row_format.format(*content))
        write_line_to_monitorfile(sep_end)
    elif output_type == 'update':
        sys.stdout.write("\x1b[1A"*(-delete_line_index))
        
        write_line_to_monitorfile(sep_end)
    else:
        sys.stdout.write("\x1b[1A"*(-delete_line_index))
        
        write_line_to_monitorfile(row_format.format(*content))
        write_line_to_monitorfile(sep_end)
    
    write_line_to_monitorfile(sep_top)
    write_line_to_monitorfile(progress_format.format(*progressbar_content))
    write_line_to_monitorfile(progress_format.format(*progressbartime_content))
    write_line_to_monitorfile(sep_end)

## INFORMATIONS ============================================================================================================

def get_value_of_variable_from_file(file, file_index, relative_index, string):
    try:
        content = [i.strip().split() for i in open(file).readlines()]
        index = [idx for idx, s in enumerate(content) if string in s][file_index]
        value = content[index][relative_index]
        return value
    except IndexError:
        print_table_row(['ERROR', 'String not found in file'], 
                        nTimestepsRequired, nTimestepsRequired, runCounter, currentTime)
        quit()

## FILE ====================================================================================================================

def write_add_string_into_file(file, substring, add, comment = None):

    with open(file, 'r') as f:
        data = f.readlines()

    try:
        index = [idx for idx, s in enumerate(data) if substring in s][0]
        data[index] = substring + add + '\n'
        
    except IndexError:
        index = [idx for idx, s in enumerate(data) if '\n' in s][1]
        data.insert(index, '\n')
        data.insert(index + 1, comment)
        data.insert(index + 2, substring + add + '\n')
        
    with open(file, 'w') as f:
        f.writelines(data)
        
def write_file(filename, content):
    file = open(filename, 'w+')
    file.write(content)
    file.close()
        
def delete_line_in_file(file, start=None, end=None):
    f = open(file, )
    lines = f.readlines()[start:end]
    
    try:
        content = ''.join(lines)
        f = open(file, 'w+')
        f.write(content)
        f.close()
    except IndexError:
        pass
    
def write_line_to_monitorfile(content, mode = 'a'):
    with open(monitorFilename, mode) as f:
        f.write(content + '\n')
        f.close()
    print(content)

## TIME ====================================================================================================================

def format_num(time):
    if time < 10:
        return '0' + str(time)
    else:
        return str(time)
    
def get_time_as_string(sec):
    
    mins,  sec   = sec   // 60, sec   % 60
    hours, mins  = mins  // 60, mins  % 60
    days,  hours = hours // 24, hours % 24
    weeks, days  = days  // 7,  days  % 7
    
    timeConvertedString = (str(int(weeks)) + ':' + 
                     format_num(int(days)) + ':' + format_num(int(hours)) + ':' + 
                     format_num(int(mins)) + ':' + format_num(int(sec)))
    
    return timeConvertedString

def get_time_in_seconds(time):
    
    # Format D-HH:MM:SS or HH:MM:SS
    
    try:
        day_split = time.split('-')
        time_split = day_split[1].split(':')

        day_sec = int(day_split[0])*24*60*60
        hour_sec = int(time_split[0])*60*60
        min_sec = int(time_split[1])*60
        sec = int(time_split[2])

        time_sec = day_sec + hour_sec + min_sec + sec
        
    except IndexError:
        time_split = time.split(':')
        
        hour_sec = int(time_split[0])*60*60
        min_sec = int(time_split[1])*60
        sec = int(time_split[2])

        time_sec = hour_sec + min_sec + sec

    return time_sec

def time_date():
    e = datetime.datetime.now()
    return '%s-%s-%s' % (format_num(e.year), format_num(e.month), format_num(e.day))

def time_time():
    e = datetime.datetime.now()
    return '%s:%s:%s' % (format_num(e.hour), format_num(e.minute), format_num(e.second))

def time_duration(startTime):
    stop = time.time()
    return get_time_as_string(stop - startTime)

## STATUS ==================================================================================================================

def get_job_status(user):
    
    jobStatusRunning = subprocess.getoutput(commandJobRunning + user).strip().split()
    jobStatusPending = subprocess.getoutput(commandJobPending + user).strip().split()

    return jobStatusRunning, jobStatusPending

def set_output_type(user):
    jobStatusRunning, jobStatusPending = get_job_status(user)

    if jobName in jobStatusRunning:
        outputType = 'running'
    elif jobName in jobStatusPending:
        outputType = 'pending'
    else:
        outputType ='no Output'

    return outputType

## MAIL ====================================================================================================================

def send_mail(recipient, subject, body = None):

    recipient = recipient.encode('utf_8')
    
    subject = '"' + subject + '"'
    subject = subject.encode('utf_8')
    
    if body == None:
        body = 'For futher information open attachment'
    
    body = body.encode('utf_8')
    
    attachmentPath = folder + '/' + monitorFilename
    attachment = attachmentPath.encode('utf_8')

    process = subprocess.Popen(['ssh', 'master', '/usr/bin/mailx', '-s', subject, '-a',  attachment, recipient],
                               stdin=subprocess.PIPE)
    process.communicate(body)


# JOB INIT =================================================================================================================

startTime = time.time()
user = os.getlogin()
write_line_to_monitorfile('STATUS OF RUN', 'w+')

print_table_row(['OUTPUT', 'INFO'], 
                0, nTimestepsRequired, runCounter, currentTime, output_type='header')

folder = os.path.dirname(os.path.abspath(__file__))
path = os.path.dirname(os.path.abspath(__file__)).split(user + '/')[1]

## BACKUP PATH =============================================================================================================

if backup:
    if backupLocation[-1] != '/':
        backupLocation += '/'
    backupPath = backupLocation + path
    
    if not os.path.exists(backupPath):
        os.makedirs(backupPath)

# START/RESTART JOB ========================================================================================================

outputType = set_output_type(user)

## BEGIN ===================================================================================================================

while True:
    try:
        nTimestepsCurrent = int(get_value_of_variable_from_file('./' + restartFilename, 0, 2, restartFlag))
        
        if nTimestepsCurrent >= nTimestepsRequired:
            print_table_row(['CONTROL', 'Current Timesteps ' + str(nTimestepsCurrent)], 
                            nTimestepsCurrent, nTimestepsRequired, runCounter, currentTime)
            
            print_table_row(['SUCCESS', 'Stop monitoring ' + jobName], 
                            nTimestepsRequired, nTimestepsRequired, runCounter, walltime, output_type='middle')
            
            if emailNotification:
                send_mail(emailAddress, 'Ended Job ' + jobName)

            quit()
        
        # Continue
        else:
            print_table_row(['CONTINUE', 'Continue monitoring ' + jobName], 
                            nTimestepsCurrent, nTimestepsRequired, runCounter, currentTime, output_type='middle')
            
            if outputType == 'running':
                runCounter += 1
                
                jobStatusRunning, jobStatusPending = get_job_status(user)
                jobStatusRunningNameIndex = [idx for idx, s in enumerate(jobStatusRunning) if jobName in s][0]
                currentTime = jobStatusRunning[jobStatusRunningNameIndex + 3]
                
                print_table_row(['CONTROL', 'Current Timesteps ' + str(nTimestepsCurrent)], 
                                nTimestepsCurrent, nTimestepsRequired, runCounter, currentTime)
            
            if emailNotification:
                send_mail(emailAddress, 'Continued Job ' + jobName)
            break
    
    # Start    
    except FileNotFoundError:
        nTimestepsCurrent = 0
        
        print_table_row(['STARTING', 'Start monitoring ' + jobName], 
                        nTimestepsCurrent, nTimestepsRequired, runCounter, currentTime, output_type='middle')
        
        if backup:
            print_table_row(['BACKUP', backupLocation])
            subprocess.run(['rsync', '-a', '', backupPath])
            
        if emailNotification:
            send_mail(emailAddress, 'Started Job ' + jobName)
        break

## MONITOR ROUTINE =========================================================================================================

while True:
    
    jobStatusRunning, jobStatusPending = get_job_status(user)

    # Job running
    if jobName in jobStatusRunning:
        
        jobStatusRunningNameIndex = [idx for idx, s in enumerate(jobStatusRunning) if jobName in s][0]
        jobID = jobStatusRunning[jobStatusRunningNameIndex - 2]
        
        currentTime = jobStatusRunning[jobStatusRunningNameIndex + 3]
        
        if outputType == 'running':
            print_table_row(['RUNNING', 'Job is executed'], 
                            nTimestepsCurrent, nTimestepsRequired, runCounter, currentTime)
            outputType = 'no Output'
        else:
            print_table_row(['RUNNING', 'Job is executed'], 
                            nTimestepsCurrent, nTimestepsRequired, runCounter, currentTime, output_type='update')
            
        time.sleep(sleepTime)
        
    # Job pending
    elif jobName in jobStatusPending:

        if outputType == 'pending':
            print_table_row(['WAITING', 'Job is pending'], 
                            nTimestepsCurrent, nTimestepsRequired, runCounter, currentTime)
            outputType = 'running'
            
        time.sleep(sleepTime)

    # Job start/restart
    else:
        
        # Check error and making Backup
        while True:
            try:
                outputContent = open(outputFilename(jobID)).readlines()[0].replace('\n','')
                
                if outputContent == outputCriteria: 
                    
                    if backup:
                        print_table_row(['BACKUP', backupLocation])
                        subprocess.run(['rsync', '-a', '', backupPath], 
                                       nTimestepsCurrent, nTimestepsRequired, runCounter, currentTime)
                    
                    break
                else:
                    print_table_row(['ERROR', outputContent[0:27]], 
                                    nTimestepsCurrent, nTimestepsRequired, runCounter, currentTime)
                        
                    if backup:
                        print_table_row(['RESTORE', backupLocation])
                        subprocess.run(['rsync', '-a', '-I', '--exclude=status.txt', backupPath + '/', ''], 
                                       nTimestepsCurrent, nTimestepsRequired, runCounter, currentTime)
                    
                    break
                    
            except (IndexError, FileNotFoundError):
                time.sleep(30)
            except NameError:
                break
        
        try:
            nTimestepsCurrent = int(get_value_of_variable_from_file('./' + restartFilename, 0, 2, restartFlag))
            print_table_row(['CONTROL', 'Current Timesteps ' + str(nTimestepsCurrent)], 
                            nTimestepsCurrent, nTimestepsRequired, runCounter, currentTime)
            
            if nTimestepsCurrent >= nTimestepsRequired:
                print_table_row(['SUCCESS', 'Stop monitoring ' + jobName], 
                                nTimestepsRequired, nTimestepsRequired, runCounter, walltime)
                
                if emailNotification:
                    send_mail(emailAddress, 'Ended Job ' + jobName)
                break
                
        except FileNotFoundError:
            pass   
        
        # Create jobscript
        if jobscriptFilename == 'jobscript-create':
            jobscriptFilename = 'jobscript'
            write_file(jobscriptFilename, jobscriptContent)
        
        # Start Job
        startOutput = subprocess.check_output([commandJobStarting, jobscriptFilename]).decode('utf-8').replace('\n', '')
        jobID = startOutput.split(startOutputFlag)[1]
        
        runCounter += 1
        
        print_table_row(['STARTING', startOutput], 
                        nTimestepsCurrent, nTimestepsRequired, runCounter, currentTime, output_type='middle')
        
        try:
            if restartMail and emailNotification:
                send_mail(emailAddress, 'Restart Job ' + jobName)
        except NameError:
            continue
        
        time.sleep(30)
        
        outputType = set_output_type(user)
        
## RESTART =================================================================================================================
