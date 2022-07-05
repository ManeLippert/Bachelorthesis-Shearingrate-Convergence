#!/usr/bin/env python3

# Author: Manuel Lippert (GitHub: ManeLippert)
# Description:
# Script that starts a given job wth GKW with slurm manager until the previous defined timesteps are completed. 
# Everything has to be defined in the jobscript with name "jobscript*" and "jobscript_end*". 
# To get a notification at the end of the run add in "jobscript_end*" the necessary flag for sbatch this job gets
# executed after the timesteps are completed, so this script runs on iteration longer as needed 
# (it was easier and quicker to code). If you do not need a notification than do not include a file named
# "jobscript_end*"

from time import sleep 
import sys
import os
import subprocess

'''
        '#########################################################'
        '#                                                       #'
        '#                      FUNCTIONS                        #'
        '#                                                       #'
        '#########################################################'
'''

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

print(  '                                                         \n'+
        '#########################################################\n'+
        '#                                                       #\n'+
        '#                   JOB INFORMATIONS                    #\n'+
        '#                                                       #\n'+
        '#########################################################\n'+
        '                                                           ')

# User name
user = os.getlogin()

# Job script informations
try:
    jobscript = [filename for filename in os.listdir('.') if filename.startswith("jobscript")][0]
except IndexError:
    print('! No jobscript found !')
    quit()
jobscriptContent = open(jobscript, 'r').read().splitlines()

# Job name
jobNameIndex = [idx for idx, s in enumerate(jobscriptContent) if 'job-name=' in s][0]
jobName = jobscriptContent[jobNameIndex].split('job-name=', 1)[1]
print('Name: ', jobName)

# Job repetition
try:
    jobRepetitionIndex = [idx for idx, s in enumerate(jobscriptContent) if 'repeat=' in s][0]
    jobRepetition = float(jobscriptContent[jobRepetitionIndex].split('repeat=', 1)[1])
except IndexError:
    jobRepetition = 1.0
print('Repetitions:', jobRepetition)

# Timesteps GKW makes from input.dat
try:
    nTimesteps = get_value_of_variable_in_input_file("./input.dat", "NTIME")
except FileNotFoundError:
    print('! No input file found !')
    quit()

# Number of required timesteps        
nTimestepsRequired = int(nTimesteps * jobRepetition)
print('Timesteps:', nTimesteps)
print('Required Timesteps:', nTimestepsRequired)

## Walltime
walltimeIndex = [idx for idx, s in enumerate(jobscriptContent) if 'time=' in s][0]
walltime = jobscriptContent[walltimeIndex].split('time=', 1)[1].replace('-', ':').split(':')
walltimeSeconds = get_time_in_seconds(walltime)

print('Walltime/s:', walltimeSeconds)

## Sleep time until one job is completed = 24h + 5min
sleepTime = walltimeSeconds + 5*30
print('Sleep Time/s:', sleepTime)

# Job status informations
jobHeader = ['JOBID', 'PARTITION', 'NAME', 'USER', 'ST', 'TIME', 'NODES', 'NODELIST(REASON)']
#jobStatus = subprocess.getoutput("squeue -u bt712347").strip().split()

'''
        '#########################################################'
        '#                                                       #'
        '#                       E -MAIL                         #'
        '#                                                       #'
        '#########################################################'
'''

# mail address
try:
    emailAddressIndex = [idx for idx, s in enumerate(jobscriptContent) if 'email=' in s][0]
    emailAddress = jobscriptContent[emailAddressIndex].split('repeat=', 1)[1]
except IndexError:
    print('! No mail address given in jobscript !')
    quit()


'''
print(  '                                                         \n'+
        '#########################################################\n'+
        '#                                                       #\n'+
        '#                  START/RESTART JOB                    #\n'+
        '#                                                       #\n'+
        '#########################################################\n'+
        '                                                           ')

while True:
    
    # Check if no Job is running or pending and start/restart job
    ## Check running jobs status and get job running time
    ### If job is not running than wait 30min
    while True:
        
        jobStatusRunning = subprocess.getoutput("squeue --states=running -u " + user).strip().split()
        jobStatusPending = subprocess.getoutput("squeue --states=pending -u " + user).strip().split()

        # Check Status if job is running
        if jobName in jobStatusRunning:
            # Job ID Running
            jobStatusRunningNameIndex = [idx for idx, s in enumerate(jobStatusRunning) if jobName in s][0]
            jobIDRunning = jobStatusRunning[jobStatusRunningNameIndex - 2]

            # Job time
            jobTime = jobStatusRunning[jobStatusRunningNameIndex + 3].split(':')
            jobTimeSeconds = get_time_in_seconds(jobTime)

            break

        # Wait 30min if job is pending
        elif jobName in jobStatusPending:
            print('Job is pending -> sleeping for 30min')
            sleep(30*60)

        # Check if no Job is running or pending and start/restart job
        else:
            subprocess.run(["sbatch", jobscript])

    # Sleep time  = 24h - job time + 5min
    sleepTimeCurrent = sleepTime - jobTimeSeconds
    print('Current Sleep Time/s:', sleepTimeCurrent)
    print('Sleeping for', sleepTimeCurrent, 's')
    sleep(sleepTimeCurrent)
    
    # read FDS.dat (restart file) to a list of lists
    ## If job is not existing than wait 30min time
    while True:
        try:
            nTimestepsCurrent = get_value_of_variable_in_input_file("./FDS.dat", 'FILE_COUNT')
            print('Current Timesteps:', nTimestepsCurrent)
            break
        except FileNotFoundError:
            print('FDS.dat not found -> sleeping for 30min')
            sleep(30*60)
            
    # Check if gkw has run requiered timesteps
    if nTimestepsCurrent >= nTimestepsRequired:

        try:
            # Restart job a last time to get notification per email
            jobscriptEnd = [filename for filename in os.listdir('.') if filename.startswith("jobscript_end")][0]
            print('Current Timesteps greater or equals Reqired Timesteps\n'+
                  'Starting last run of GKW!')
            subprocess.run(["sbatch", jobscriptEnd])
            break
        except IndexError:
            # Finished monitoring
            print('Current Timesteps greater or equals Reqired Timesteps\n'+
                  'Stop monitoring of GKW!')
            break
'''