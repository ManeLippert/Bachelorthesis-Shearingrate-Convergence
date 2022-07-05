#!/usr/bin/env python3

# AUTHOR: Manuel Lippert (GitHub: ManeLippert (https://github.com/ManeLippert))
#
# DESCRIPTION:
# Script that starts a given job wth GKW with slurm manager until the previous defined timesteps are completed. 
# Everything has to be defined in the jobscript with name "jobscript*". It also sends an mail alert when the job gets started and ended although the end alert is still in work.

# ADDITIONAL FLAGS IN JOBSCRIPT:
# 1) repeat=3 
#    
#    number how often GKW should restart its current timesteps until required timesteps are given. If flag is not given repeat=1 automatically.
#    (Example: 10000 current timesteps and repeat=3 than required timesteps are 30000)

# 2) email=mail@email.de
#
#    mail address to which email the notifaction shold be sended

# START SCRIPT:
# To start the script in the background following command is needed:
#
# >>> nohup python3 -u monitor_job.py &
#
# This will write every output in the file nohup.out that will be send as mail body to the defined mail address


from time import sleep 
import sys
import os
import subprocess
'''
#########################################################
#                                                       #
#                       VARIABLES                       #
#                                                       #
#########################################################
'''
# Filename of jobscript to get data
jobscriptFilename = 'jobscript'

# Flags of the specific informations
jobNameFlag = 'job-name='
jobRepetitionFlag = 'repeat='
emailAddressFlag = 'email='
walltimeFlag = 'time='

'''
'#########################################################'
'#                                                       #'
'#                      FUNCTIONS                        #'
'#                                                       #'
'#########################################################'
'''

def read_file_to_string(file):
    content = ''.join(open(file).readlines())
    
    return content

def get_value_of_variable_in_input_file(file, string):
    try:
        content = [i.strip().split() for i in open(file).readlines()]
        index = [idx for idx, s in enumerate(content) if string in s][0]
        value = int(content[index][2])
        return value
    except IndexError:
        print('! String not in input file !')

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

def job_informations():
    content = ('################' +                               '\n' +
               'JOB INFORMATIONS' +                               '\n' +
               '################' +                               '\n' +
                                                                  '\n' +
               'Name: ' + jobName +                               '\n' +
               'Repetitions: ' + str(jobRepetition) +             '\n' +
               'E-Mail: ' + emailAddress +                        '\n' +
               'Timesteps: ' + str(nTimesteps) +                  '\n' +
               'Required Timesteps: ' + str(nTimestepsRequired) + '\n' +
               'Walltime/s: ' + str(walltimeSeconds) +            '\n' +
               'Sleep Time/s: ' + str(sleepTime) +                '\n' )

    return content

def job_start():
    content = ('################' +                               '\n' +
               'RESTART JOB'      +                               '\n' +
               '################' +                               '\n' )
    
    return content

def get_job_status(user):
    jobStatusRunning = subprocess.getoutput('squeue --states=running -u ' + user).strip().split()
    jobStatusPending = subprocess.getoutput('squeue --states=pending -u ' + user).strip().split() 

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


def send_mail(recipient, subject, body):

    recipient = recipient.encode('utf_8')
    subject = '"' + subject + '"'
    subject = subject.encode('utf_8')
    body = body.encode('utf_8')

    process = subprocess.Popen(['ssh', 'master', '/usr/bin/mailx', '-s', subject, recipient],
                               stdin=subprocess.PIPE)
    process.communicate(body)

'''
#########################################################
#                                                       #
#                  JOB INFORMATIONS                     #
#                                                       #
#########################################################
'''

# User name
user = os.getlogin()

# Job script informations
try:
    jobscript = [filename for filename in os.listdir('.') if filename.startswith(jobscriptFilename)][0]
except IndexError:
    print('! No jobscript found !')
    quit()
jobscriptContent = open(jobscript, 'r').read().splitlines()

# Job name
jobNameIndex = [idx for idx, s in enumerate(jobscriptContent) if jobNameFlag in s][0]
jobName = jobscriptContent[jobNameIndex].split(jobNameFlag, 1)[1]

# Job repetition
try:
    jobRepetitionIndex = [idx for idx, s in enumerate(jobscriptContent) if jobRepetitionFlag in s][0]
    jobRepetition = float(jobscriptContent[jobRepetitionIndex].split(jobRepetitionFlag, 1)[1])
except IndexError:
    jobRepetition = 1.0

# Timesteps GKW makes from input.dat
try:
    nTimesteps = get_value_of_variable_in_input_file('./input.dat', 'NTIME')
except FileNotFoundError:
    print('! No input file found !')
    quit()

# Mail address
try:
    emailAddressIndex = [idx for idx, s in enumerate(jobscriptContent) if emailAddressFlag in s][0]
    emailAddress = jobscriptContent[emailAddressIndex].split(emailAddressFlag, 1)[1]
except IndexError:
    print('! No mail address given in jobscript !')
    quit()

# Number of required timesteps        
nTimestepsRequired = int(nTimesteps * jobRepetition)

## Walltime
walltimeIndex = [idx for idx, s in enumerate(jobscriptContent) if walltimeFlag in s][0]
walltime = jobscriptContent[walltimeIndex].split(walltimeFlag, 1)[1].replace('-', ':').split(':')
walltimeSeconds = get_time_in_seconds(walltime)

## Sleep time 5min (time code checks status)
sleepTime = 5*60

# Job status informations
#jobHeader = ['JOBID', 'PARTITION', 'NAME', 'USER', 'ST', 'TIME', 'NODES', 'NODELIST(REASON)']
#jobStatus = subprocess.getoutput("squeue -u bt712347").strip().split()

print(job_informations())

# Send start mail
send_mail(emailAddress, 'Started Job ' + jobName, read_file_to_string('./nohup.out'))

'''
#########################################################
#                                                       #
#                  START/RESTART JOB                    #
#                                                       #
#########################################################
'''

print(job_start())

# To limit repeating outputs
outputType = set_output_type(user)

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
            jobIDRunning = jobStatusRunning[jobStatusRunningNameIndex - 2]

            # Job time
            jobTime = jobStatusRunning[jobStatusRunningNameIndex + 3].split(':')
            jobTimeSeconds = get_time_in_seconds(jobTime)

            # Set output type
            if outputType == 'running':
                print('Job is running with ID: ' + jobIDRunning)
                outputType = 'check FDS.dat'

            break

        # Wait 30min if job is pending
        elif jobName in jobStatusPending:
            # Set output type
            if outputType == 'pending':
                print('Job is pending -> wait till job is running')
                outputType = 'running'

            sleep(sleepTime)

        # Check if no Job is running or pending and start/restart job
        else:
            subprocess.run(["sbatch", jobscript])
            sleep(30)
            outputType = set_output_type(user)


    # read FDS.dat (restart file) to a list of lists
    ## If job is not existing than wait 30min time
    while True:
        try:
            nTimestepsCurrent = get_value_of_variable_in_input_file('./FDS.dat', 'FILE_COUNT')

            # Set output type
            if outputType == 'check Timestep' or outputType == 'check FDS.dat':
                print('Current Timesteps:', nTimestepsCurrent)
                outputType = 'no Output'
                
            break
        except FileNotFoundError:
            # Set output type
            if outputType == 'check FDS.dat':
                print('FDS.dat not found -> wait until file gets generated')
                outputType = 'check Timestep'

            sleep(sleepTime)

    # Check if gkw has run requiered timesteps
    if nTimestepsCurrent >= nTimestepsRequired:
        print('Current Timesteps greater or equals Reqired Timesteps\n'+
              'Stop monitoring of GKW!')
        send_mail(emailAddress, 'Ended Job ' + jobName, read_file_to_string('./nohup.out'))

        break
    else:
        sleep(sleepTime)