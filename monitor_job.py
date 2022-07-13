#!/usr/bin/env python3

# AUTHOR: Manuel Lippert (GitHub: ManeLippert (https://github.com/ManeLippert))
#
# DESCRIPTION:
# Script that starts a given job wth GKW with slurm manager until the previous defined timesteps are completed. 
# Everything has to be defined in the jobscript with name "jobscript*" (Can bei modified under the variable section of this script). 
# It also sends an mail alert when the job gets started and ended.

# ADDITIONAL FLAGS IN JOBSCRIPT:
# 1) repeat=3 
#    
#    number how often GKW should restart its current timesteps until required timesteps are given. If flag is not given repeat=1 automatically.
#    (Example: 10000 current timesteps and repeat=3 than required timesteps are 30000)
#
# 2) email=mail@email.de
#
#    mail address to which email the notifaction shold be sended
#
# 3) backup-location=$PATH
#
#    path to the backup location. Has to be the absoute path (example ~/$PATH)!
#
# Every flag can be modified under the variable section of this script.

# START SCRIPT:
# To start the script in the background following command is needed:
#
# >>> nohup python3 -u monitor_job.py &> status.txt &
#
# Output:
#
# >>> [1] 10537
#
# This will write every output in the file jobstatus.txt that will be send as attachment to the defined mail address.

# LIST PRGRESS:
# To see which progress is in background running following command is needed:
#
# >>> ps ax | grep monitor_job.py
#
# Output:
#
# >>> 10537 pts/1    S      0:00 python3 -u monitor_job.py
# >>> 23426 pts/1    S+     0:00 grep --color=auto monitor_job.py
# 
# This will give you the ID to kill monitoring script with the command:
#
# >>> kill 10537
#
#  This will kill the monitor script.


import datetime
import time
from time import sleep 
import os
import subprocess


'''
#########################################################
#                                                       #
#                       VARIABLES                       #
#                                                       #
#########################################################
'''
# Filenames
jobscriptFilename = 'jobscript'
restartFilename = 'FDS.dat'
monitorFilename = 'status.txt'
inputFilename = 'input.dat'

# Flags of the specific informations
jobNameFlag = 'job-name='
jobRepetitionFlag = 'repeat='
emailAddressFlag = 'email='
walltimeFlag = 'time='
backupFlag = 'backup-location='

inputFlag = 'NTIME'
restartFlag = 'FILE_COUNT'

# Sleep time 5min (time code checks status)
sleepTime = 5*60

# Commands
commandJobRunning = 'squeue --states=running -u '
commandJobPending = 'squeue --states=pending -u '

commandJobStarting = 'sbatch'

'''
'#########################################################'
'#                                                       #'
'#                      FUNCTIONS                        #'
'#                                                       #'
'#########################################################'
'''

# Functions to get informations out of files
def read_file_to_string(file):
    content = ''.join(open(file).readlines())
    
    return content

def get_value_of_variable_from_input_file(file, string):
    try:
        content = [i.strip().split() for i in open(file).readlines()]
        index = [idx for idx, s in enumerate(content) if string in s][0]
        value = int(content[index][2])
        return value
    except IndexError:
        print('! String not in input file !')
        quit()
        
def get_job_information_from_jobscript_flag(content, flag):
    index = [idx for idx, s in enumerate(content) if flag in s][0]
    info = content[index].split(flag, 1)[1]
    
    return info

# Time functions
def get_time_in_seconds(time):
    # Check time format of time
    if len(time) < 5:
        ## d:hh:mm:ss
        if len(time) == 4:
            timeSeconds = int(time[0])*24*60*60 + int(time[1])*60*60 + int(time[2])*60 + int(time[3])
        ## hh:mm:ss
        elif len(time) == 3:
            timeSeconds = int(time[0])*60*60 + int(time[1])*60 + int(time[2])
        ## mm:ss
        elif len(time) == 2:
            timeSeconds = int(time[0])*60 + int(time[1])
        ## ss
        elif len(time) == 1:
            timeSeconds = int(time[0])
    else:
        print('! Time format is not supported !')
        quit()
    
    return timeSeconds

def format_num(time):
    if time < 10:
        return '0' + str(time)
    else:
        return str(time)

def get_time_converted(sec):
    mins, sec  = sec // 60, sec % 60
    hours, mins = mins // 60, mins % 60
    days, hours = hours // 24, hours % 24
    weeks, days = days // 7, days % 7
    
    converted = (str(int(weeks)) + ':' + 
                 format_num(int(days)) + ':' + format_num(int(hours)) + ':' + 
                 format_num(int(mins)) + ':' + format_num(int(sec)))
    
    return converted

def time_date():
    e = datetime.datetime.now()
    return '%s-%s-%s' % (format_num(e.year), format_num(e.month), format_num(e.day))

def time_time():
    e = datetime.datetime.now()
    return '%s:%s:%s' % (format_num(e.hour), format_num(e.minute), format_num(e.second))

def time_duration(startTime):
    stop = time.time()
    return get_time_converted(stop - startTime)

# Table functions
def table_row_format(content):
    if len(content) == 5:
        cols = [2, 10, 25, 41, 2]
    elif len(content) == 7:
        cols = [2, 10, 29, 13, 11, 13, 2]
    elif len(content) == 4:
        cols = [2, 10, 66, 2]
    
    row_format = "".join(["{:<" + str(col) + "}" for col in cols])
    
    return row_format

def print_table_row(content):
    
    if isinstance(content[0], list):
        row_format = table_row_format(content[0])
        for row in content:
            print(row_format.format(*row))
        
    else:
        row_format = table_row_format(content)
        print(row_format.format(*content))

# Output Strings as function
def job_init():
    content = [['o-', '----------', '------------------------------------------------------------------', '-o'],
               ['| ', 'OUTPUT', 'JOB INITIALIZE', ' |'],
               ['o-', '----------', '------------------------------------------------------------------', '-o']]
    return content

def job_informations():
    content = [['o-', '----------', '-------------------------', '-----------------------------------------', '-o'],
               ['| ', 'OUTPUT', 'JOB INFORMATIONS', 'VALUE', ' |'],
               ['o-', '----------', '-------------------------', '-----------------------------------------', '-o'],
               ['| ', 'INFO', 'Name', jobName, ' |'],
               ['| ', 'INFO', 'Repetitions', str(jobRepetition), ' |'],
               ['| ', 'INFO', 'E-Mail', emailAddress, ' |'],
               ['| ', 'INFO', 'Timesteps', str(nTimesteps), ' |'],
               ['| ', 'INFO', 'Required Timesteps', str(nTimestepsRequired), ' |'],
               ['| ', 'INFO', 'Walltime/s', str(walltimeSeconds), ' |'],
               ['| ', 'INFO', 'Sleep Time/s', str(sleepTime), ' |'],
               ['o-', '----------', '-------------------------', '-----------------------------------------', '-o']]
    return content

def job_monitor():
    content = [['o-', '----------', '-----------------------------', '-------------', '-----------', '-------------', '-o'],
               ['| ', 'OUTPUT', 'JOB MONITORING', 'DATE', 'TIME', 'W:DD:HH:MM:SS', ' |'],
               ['o-', '----------', '-----------------------------', '-------------', '-----------', '-------------', '-o']]
    return content

# Status functions for setting values for the job
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

# Mail command
def send_mail(recipient, subject, body):

    recipient = recipient.encode('utf_8')
    subject = '"' + subject + '"'
    subject = subject.encode('utf_8')
    body = body.encode('utf_8')
    attachment = monitorFilename.encode('utf_8')

    process = subprocess.Popen(['ssh', 'master', '/usr/bin/mailx', '-s', subject, '-a',  path + '/' + monitorFilename, recipient],
                               stdin=subprocess.PIPE)
    process.communicate(body)

'''
#########################################################
#                                                       #
#                  JOB INFORMATIONS                     #
#                                                       #
#########################################################
'''
# Output init header
print('\n')
print_table_row(job_init())

# User name
user = os.getlogin()

# Job script informations
try:
    jobscript = [filename for filename in os.listdir('.') if filename.startswith(jobscriptFilename)][0]
    print_table_row(['| ', 'SUCCESS', 'Found ' + jobscriptFilename, ' |'])
except IndexError:
    print_table_row([['| ', 'ERROR', 'No file "' + jobscriptFilename + '" found', ' |'],
                     ['o-', '----------', '------------------------------------------------------------------', '-o']])
    quit()
    
jobscriptContent = open(jobscript, 'r').read().splitlines()

# Job name
jobName = get_job_information_from_jobscript_flag(jobscriptContent, jobNameFlag)

# Job repetition
try:
    jobRepetition = float(get_job_information_from_jobscript_flag(jobscriptContent, jobRepetitionFlag))
except IndexError:
    jobRepetition = 1.0

# Timesteps GKW makes from input.dat
try:
    nTimesteps = get_value_of_variable_from_input_file('./' + inputFilename, inputFlag)
    print_table_row([['| ', 'SUCCESS', 'Found ' + inputFilename, ' |'],
                     ['o-', '----------', '------------------------------------------------------------------', '-o']])
except FileNotFoundError:
    print_table_row([['| ', 'ERROR', 'No file "' + inputFilename + '" found', ' |'],
                     ['o-', '----------', '------------------------------------------------------------------', '-o']])
    quit()

# Mail address
try:
    emailAddress = get_job_information_from_jobscript_flag(jobscriptContent, emailAddressFlag)
    
    # Email notification switch
    emailNotification = True
    
except IndexError:
    # Email notification switch
    emailNotification = False

# Number of required timesteps        
nTimestepsRequired = int(nTimesteps * jobRepetition)

## Walltime
walltime = get_job_information_from_jobscript_flag(jobscriptContent, walltimeFlag).replace('-', ':').split(':')
walltimeSeconds = get_time_in_seconds(walltime)

# Job status informations
#jobHeader = ['JOBID', 'PARTITION', 'NAME', 'USER', 'ST', 'TIME', 'NODES', 'NODELIST(REASON)']
#jobStatus = subprocess.getoutput("squeue -u bt712347").strip().split()

# Backup and data path
try:
    # Backup and data location
    backupLocation = get_job_information_from_jobscript_flag(jobscriptContent, backupFlag)
    
    path = os.path.dirname(os.path.abspath(__file__)).split(user + '/')[1]
    backupPath = backupLocation + path
    
    ## Create backup directory if do not exist
    if not os.path.exists(backupPath):
        os.makedirs(backupPath)
    
    # BackUp switch
    backup = True
except IndexError:
    # BackUp switch
    backup = False

# Output job informations
print('\n')
print_table_row(job_informations())

# Send start mail
if emailNotification:
    send_mail(emailAddress, 'Started Job ' + jobName, 'For futher information open attachment')

'''
#########################################################
#                                                       #
#                  START/RESTART JOB                    #
#                                                       #
#########################################################
'''

# Output job monitor header
print('\n')
print_table_row(job_monitor())

# Set start time for stop watch
startTime = time.time()

print_table_row(['| ', 'STARTING', 'Start monitoring', time_date(), time_time(), time_duration(startTime) , ' |'])

####################    BEGIN    ########################

# To limit repeating outputs
outputType = set_output_type(user)

# read FDS.dat (restart file) to a list of lists
## If gkw has run requiered timesteps stop already here
while True:
    try:
        nTimestepsCurrent = get_value_of_variable_from_input_file('./' + restartFilename, restartFlag)
        
        # Check if gkw has run requiered timesteps
        if nTimestepsCurrent >= nTimestepsRequired:
            print_table_row([['| ', 'SUCCESS', 'Stop monitoring', time_date(), time_time(), time_duration(startTime) , ' |'],
                             ['o-', '----------', '-----------------------------', '-------------', '-----------', '-------------', '-o']])
            
            # Send end email
            if emailNotification:
                send_mail(emailAddress, 'Ended Job ' + jobName, 'For futher information open attachment')

            quit()
        else:
            break
        
    except FileNotFoundError:
        break

################## MONITOR ROUTINE #######################

while True:
    
    # Check if no Job is running or pending and start/restart job
    ## Check running jobs status and get job running time
    ### If job is not running than wait 30min and check again
    while True:
        
        jobStatusRunning, jobStatusPending = get_job_status(user)

        # Check Status if job is running
        if jobName in jobStatusRunning:
            # Job ID Running
            jobStatusRunningNameIndex = [idx for idx, s in enumerate(jobStatusRunning) if jobName in s][0]
            jobID = jobStatusRunning[jobStatusRunningNameIndex - 2]

            # Job time
            jobTime = jobStatusRunning[jobStatusRunningNameIndex + 3].split(':')
            jobTimeSeconds = get_time_in_seconds(jobTime)

            # Set output type
            if outputType == 'running':
                print_table_row(['| ', 'RUNNING', 'Job is executed', time_date(), time_time(), time_duration(startTime) , ' |'])
                
                outputType = 'check ' + restartFilename

            break

        # Wait 30min if job is pending
        elif jobName in jobStatusPending:
            # Set output type
            if outputType == 'pending':
                print_table_row(['| ', 'WATING', 'Job is pending', time_date(), time_time(), time_duration(startTime) , ' |'])
                outputType = 'running'

            sleep(sleepTime)

        # Check if no Job is running or pending and start/restart job
        ## Check slurm output file for any errors 
        ### Print current timesteps 
        else:
            # Making backup of data
            if backup:
                print_table_row(['| ', 'BACKUP', backupLocation, time_date(), time_time(), time_duration(startTime) , ' |'])
                subprocess.run(['rsync', '-a', '', backupPath])
                
            # check slurm output for any errors
            while True:
                try:
                    slurmContent = open('./slurm-' + jobID + '.out').readlines()[0].replace('\n','')
                    
                    # Check for any errors
                    if slurmContent == '0':
                        break
                    else:
                        print_table_row([['| ', 'ERROR', 'GKW stopped job', time_date(), time_time(), time_duration(startTime) , ' |'],
                                         ['o-', '----------', '-----------------------------', '-------------', '-----------', '-------------', '-o']])
                        
                        # Send fail mail
                        if emailNotification:
                            send_mail(emailAddress, 'Failed Job ' + jobName, 'For futher information open attachment')
                        quit()
                # If jobID is not defieed
                except NameError:
                    break
                # If file is not generated
                except FileNotFoundError:
                    sleep(30)
                except IndexError:
                    sleep(30)
                
            # Timestep output
            try:
                nTimestepsCurrent = get_value_of_variable_from_input_file('./' + restartFilename, restartFlag)
                print_table_row(['| ', 'CONTROL', 'Current Timesteps' + str(nTimestepsCurrent), time_date(), time_time(), time_duration(startTime) , ' |'])
            except FileNotFoundError:
                pass    
            
            # Start Job
            startOutput = subprocess.check_output([commandJobStarting, jobscript]).decode('utf-8').replace('\n', '')
            print_table_row(['| ', 'STARTING', startOutput, time_date(), time_time(), time_duration(startTime) , ' |'])
            sleep(30)
            outputType = set_output_type(user)


    # Read FDS.dat (restart file) to a list of lists
    ## If job is not existing than wait 30min time
    while True:
        try:
            nTimestepsCurrent = get_value_of_variable_from_input_file('./' + restartFilename, restartFlag)
            break
        except FileNotFoundError:
            # Set output type
            if outputType == 'check ' + restartFilename:
                print_table_row(['| ', 'WAITING', restartFilename + ' not found', time_date(), time_time(), time_duration(startTime) , ' |'])
                outputType = 'no Output'

            sleep(sleepTime)

    # Check if gkw has run requiered timesteps
    if nTimestepsCurrent >= nTimestepsRequired:
        print_table_row([['| ', 'SUCCESS', 'Stop monitoring', time_date(), time_time(), time_duration(startTime) , ' |'],
                         ['o-', '----------', '-----------------------------', '-------------', '-----------', '-------------', '-o']])

        # Send end email
        if emailNotification:
            send_mail(emailAddress, 'Ended Job ' + jobName, 'For futher information open attachment')

        break
    else:
        sleep(sleepTime)