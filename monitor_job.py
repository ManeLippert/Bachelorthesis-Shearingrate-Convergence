#!/usr/bin/env python3

# AUTHOR: Manuel Lippert (GitHub: ManeLippert (https://github.com/ManeLippert))
#
# DESCRIPTION:
# Script that starts a given job wth GKW with slurm manager until the previous defined timesteps are completed. 
# It also sends an mail alert when the job gets started and ended.

# IMPORTANT:
# Take a look in the Variable Section

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
# This will kill the monitor script.

# OUTPUT STATUS:
# To see all status.txt files with one command follwing command is needed:
#
# >>> cd $DATA && find . -name status.txt -exec tail --lines=+10 {} \;


'''
#########################################################
#                                                       #
#                        Modules                        #
#                                                       #
#########################################################
'''

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

##################### ADDITIONAL ########################

emailAddress = 'Manuel.Lippert@uni-bayreuth.de'
backupLocation = '/scratch/bt712347/backup'
jobNamePathInfo = ['boxsize', '/Nsgrid']
jobRepetition = 4

##################### FILENAMES #########################

jobscriptFilename = 'jobscript'
restartFilename = 'FDS.dat'
monitorFilename = 'status.txt'
inputFilename = 'input.dat'

# Declared as function for dynamic changes
def outputFilename(info):
    return './slurm-' + info + '.out'

####################### FLAGS ###########################

jobNameComment = '# jobname\n'
jobNameFlag = '#SBATCH --job-name='
walltimeFlag = 'time='
startOutputFlag = 'Submitted batch job '

outputCriteria = '0'

inputFlag = 'NTIME'
restartFlag = 'FILE_COUNT'

#################### SLEEP TIME #########################

sleepTime = 5*60

##################### COMMANDS ##########################

commandJobRunning = 'squeue --states=running -u '
commandJobPending = 'squeue --states=pending -u '

commandJobStarting = 'sbatch'

'''
#########################################################
#                                                       #
#                      FUNCTIONS                        #
#                                                       #
#########################################################
'''

#################### INFORMATIONS #######################

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

######################## FILE ###########################

def write_add_string_into_file(file, substring, add, comment = None):
    # open file
    with open(file, 'r') as f:
        # read a list of lines into data
        data = f.readlines()

    try:
        index = [idx for idx, s in enumerate(data) if substring in s][0]
        data[index] = substring + add + '\n'
        
    except IndexError:
        index = [idx for idx, s in enumerate(data) if '\n' in s][1]
        data.insert(index, '\n')
        data.insert(index + 1, comment)
        data.insert(index + 2, substring + add + '\n')
        
    # replace line
    with open(file, 'w') as f:
        f.writelines(data)

######################## TIME ###########################

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
    
    timeConverted = (str(int(weeks)) + ':' + 
                 format_num(int(days)) + ':' + format_num(int(hours)) + ':' + 
                 format_num(int(mins)) + ':' + format_num(int(sec)))
    
    return timeConverted

def time_date():
    e = datetime.datetime.now()
    return '%s-%s-%s' % (format_num(e.year), format_num(e.month), format_num(e.day))

def time_time():
    e = datetime.datetime.now()
    return '%s:%s:%s' % (format_num(e.hour), format_num(e.minute), format_num(e.second))

def time_duration(startTime):
    stop = time.time()
    return get_time_converted(stop - startTime)

####################### TABLE ###########################

def table_row_format(content):
    if len(content) == 7:
        cols = [2, 10, 29, 13, 11, 13, 2]
    elif len(content) == 6:
        cols = [2, 10, 29, 13, 24, 2]
    elif len(content) == 5:
        cols = [2, 10, 25, 41, 2]
    elif len(content) == 4:
        cols = [2, 10, 66, 2]
    else:
        cols = [2, 76, 2]

    i, sep = 0, []
    
    while i < len(cols):
        if i == 0:
            sep.append('o' + (cols[i]-1)*'-')
        elif i == (len(cols)-1):
            sep.append((cols[i]-1)*'-' + 'o')
        else:
            sep.append(cols[i]*'-')
        
        i += 1

    row_format = "".join(["{:<" + str(col) + "}" for col in cols])
    
    return row_format, sep

def print_table_row(content, output_type = None, time_info = True):
    
    if isinstance(content[0], list):
        
        i = 0
        while i < len(content):
            content[i].insert(0, '| ')
            if time_info:
                if output_type == 'header':
                    content[i].append('DATE')
                    content[i].append('TIME')
                    content[i].append('W:DD:HH:MM:SS')
                else:
                    content[i].append(time_date())
                    content[i].append(time_time())
                    content[i].append(time_duration(startTime))
            content[i].insert(len(content[i]), ' |')
            
            i += 1
        
        row_format, sep = table_row_format(content[0])
        
        if output_type == 'header':
            print('\n')
            print(row_format.format(*sep))
            for row in content:
                print(row_format.format(*row))
            print(row_format.format(*sep))
            
        elif output_type == 'end':
            for row in content:
                print(row_format.format(*row))
            print(row_format.format(*sep))
            
        else:
            for row in content:
                print(row_format.format(*row))
    else:
        content.insert(0, '| ')
        if time_info:
            if output_type == 'header':
                content.append('DATE')
                content.append('TIME')
                content.append('W:DD:HH:MM:SS')
            else:
                content.append(time_date())
                content.append(time_time())
                content.append(time_duration(startTime))
        content.insert(len(content), ' |')
        
        row_format, sep = table_row_format(content)
        
        if output_type == 'header':
            print('\n')
            print(row_format.format(*sep))
            print(row_format.format(*content))
            print(row_format.format(*sep))
            
        elif output_type == 'end':
            print(row_format.format(*content))
            print(row_format.format(*sep))
            
        else:
            print(row_format.format(*content))

####################### STATUS ##########################

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

####################### MAIL ############################

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

'''
#########################################################
#                                                       #
#                  JOB INFORMATIONS                     #
#                                                       #
#########################################################
'''

################### OUTPUT STRING #######################

jobInformations = []

#################### START WATCH ########################

startTime = time.time()

#################### OUTPUT INIT ########################

print_table_row(['OUTPUT', 'JOB INITIALIZE'], output_type='header')

####################### USER ############################

user = os.getlogin()

#################### FILE PATH ##########################

folder = os.path.dirname(os.path.abspath(__file__))
path = os.path.dirname(os.path.abspath(__file__)).split(user + '/')[1]

################ JOB NAME FROM PATH #####################

try:
    i, jobName = 0, ''

    while i < (len(jobNamePathInfo)-1):
        jobName += path[path.find(jobNamePathInfo[i])+len(jobNamePathInfo[i]):path.rfind(jobNamePathInfo[i+1])]
        i += 1
    
except NameError:
    jobName = None

################## MAIL ADDRESS #########################

try:
    emailNotification = True
    jobInformations.append(['INFO', 'E-Mail', emailAddress])
except NameError:
    emailNotification = False
    
#################### JOB SCRIPT #########################

try:
    jobscript = [filename for filename in os.listdir('.') if filename.startswith(jobscriptFilename)][0]
    print_table_row(['SUCCESS', 'Found ' + jobscriptFilename])
    
except IndexError:
    print_table_row(['ERROR', 'No file "' + jobscriptFilename + '" found'])
    # Send error mail
    if emailNotification:
        send_mail(emailAddress, 'Failed Job ' + jobName)
    quit()

jobscriptContent = open(jobscript, 'r').read().splitlines()

####################### JOB NAME ########################

try:
    name = get_job_information_from_jobscript_flag(jobscriptContent, jobNameFlag)

    # Write job name from path into jobscript
    if name == '':
        write_add_string_into_file(jobscript, jobNameFlag, jobName, jobNameComment)
    else:
        jobName = name

except IndexError:
    write_add_string_into_file(jobscript, jobNameFlag, jobName, jobNameComment)

if len(jobName) > 8:
    jobName = jobName[0:8]

jobInformations.append(['INFO', 'Name', jobName])

################### JOB REPETITION ######################

try:
    jobRepetition = float(jobRepetition)
except NameError:
    jobRepetition = 1.0
    
jobInformations.append(['INFO', 'Repetition', jobRepetition])

##################### TIMESTEPS #########################

# nTimesteps
try:
    nTimesteps = get_value_of_variable_from_input_file('./' + inputFilename, inputFlag)
    print_table_row(['SUCCESS', 'Found ' + inputFilename], output_type='end', time_info=True)
except FileNotFoundError:
    print_table_row(['ERROR', 'No file "' + inputFilename + '" found'], output_type='end')
    # Send error mail
    if emailNotification:
        send_mail(emailAddress, 'Failed Job ' + jobName)
    quit()

jobInformations.append(['INFO', 'Timesteps', nTimesteps])

# Number of required timesteps        
nTimestepsRequired = int(nTimesteps * jobRepetition)

jobInformations.append(['INFO', 'Required Timesteps', nTimestepsRequired])

#################### WALLTIME ###########################

walltime = get_job_information_from_jobscript_flag(jobscriptContent, walltimeFlag).replace('-', ':').split(':')
walltimeSeconds = get_time_in_seconds(walltime)

jobInformations.append(['INFO', 'Walltime/s', walltimeSeconds])

################### BACKUP PATH #########################

try:
    # Backup location
    if backupLocation[-1] != '/':
        backupLocation += '/'
        
    backupPath = backupLocation + path
    
    # Create backup directory if do not exist
    if not os.path.exists(backupPath):
        os.makedirs(backupPath)
    
    # BackUp switch
    backup = True
    
    jobInformations.append(['INFO', 'Backup Path', backupLocation])
    
except NameError:
    # BackUp switch
    backup = False

################### OUTPUT INFO #########################

print_table_row(['OUTPUT', 'JOB INFORMATIONS', 'VALUE'], output_type='header', time_info=False)
print_table_row(jobInformations, output_type='end', time_info=False)


'''
#########################################################
#                                                       #
#                  START/RESTART JOB                    #
#                                                       #
#########################################################
'''

################## OUTPUT MONITOR #######################

print_table_row(['OUTPUT', 'JOB MONITORING'], output_type='header')

####################    BEGIN    ########################

outputType = set_output_type(user)

# read restart file 
## If gkw has run requiered timesteps stop already here
while True:
    try:
        nTimestepsCurrent = get_value_of_variable_from_input_file('./' + restartFilename, restartFlag)
        
        # Output current timesteps
        if outputType == 'running':
            print_table_row(['CONTROL', 'Current Timesteps ' + str(nTimestepsCurrent)])
        
        # Check if gkw has run requiered timesteps
        if nTimestepsCurrent >= nTimestepsRequired:
            print_table_row(['SUCCESS', 'Stop monitoring'], output_type='end')
            
            # Send end email
            if emailNotification:
                send_mail(emailAddress, 'Ended Job ' + jobName)

            quit()
        
        # Continue
        else:
            print_table_row(['CONTINUE', 'Continue monitoring'])
            
            # Send continue mail
            if emailNotification:
                send_mail(emailAddress, 'Continued Job ' + jobName)
            break
    
    # Start    
    except FileNotFoundError:
        print_table_row(['STARTING', 'Start monitoring'])
        
        # Making backup
        if backup:
            print_table_row(['BACKUP', backupLocation], output_type='end')
            subprocess.run(['rsync', '-a', '', backupPath])
        
        # Send start mail
        if emailNotification:
            send_mail(emailAddress, 'Started Job ' + jobName)
        break

################## MONITOR ROUTINE ######################

while True:
    
    jobStatusRunning, jobStatusPending = get_job_status(user)

    # Job running
    if jobName in jobStatusRunning:
        # Job ID Running
        jobStatusRunningNameIndex = [idx for idx, s in enumerate(jobStatusRunning) if jobName in s][0]
        jobID = jobStatusRunning[jobStatusRunningNameIndex - 2]
        
        # Set output type
        if outputType == 'running':
            print_table_row(['RUNNING', 'Job is executed'])
            outputType = 'no Output'
            
        sleep(sleepTime)
        
    # Job pending
    elif jobName in jobStatusPending:
        # Set output type
        if outputType == 'pending':
            print_table_row(['WAITING', 'Job is pending'])
            outputType = 'running'
            
        sleep(sleepTime)

    # Job start/restart
    else:
        
        # Check error and making Backup
        while True:
            try:
                outputContent = open(outputFilename(jobID)).readlines()[0].replace('\n','')
                
                if outputContent == outputCriteria: 
                    
                    # Making backup
                    if backup:
                        print_table_row(['BACKUP', backupLocation], output_type='end')
                        subprocess.run(['rsync', '-a', '', backupPath])
                    
                    break
                else:
                    print_table_row(['ERROR', outputContent[0:27]])
                        
                    # Restore Files from last Backup
                    if backup:
                        print_table_row(['RESTORE', backupLocation], output_type='end')
                        subprocess.run(['rsync', '-a', '-I', '--exclude=status.txt', backupPath + '/', ''])
                    
                    break
                    
            # If jobID is not defined or file is not generated
            except (IndexError, FileNotFoundError):
                sleep(30)
            except NameError:
                break
        
        # Check Timesteps
        try:
            nTimestepsCurrent = get_value_of_variable_from_input_file('./' + restartFilename, restartFlag)
            print_table_row(['CONTROL', 'Current Timesteps ' + str(nTimestepsCurrent)])
            
            # Check if gkw has run requiered timesteps
            if nTimestepsCurrent >= nTimestepsRequired:
                print_table_row(['SUCCESS', 'Stop monitoring'], output_type='end')
                
                # Send end email
                if emailNotification:
                        send_mail(emailAddress, 'Ended Job ' + jobName)
                break
                
        except FileNotFoundError:
            pass   
        
        # Start Job
        startOutput = subprocess.check_output([commandJobStarting, jobscript]).decode('utf-8').replace('\n', '')
        jobID = startOutput.split(startOutputFlag)[1]
        
        print_table_row(['STARTING', startOutput])
        
        sleep(30)
        
        outputType = set_output_type(user)
        
##################### RESTART ###########################
