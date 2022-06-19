#!/usr/bin/env python

from time import sleep 
import sys
import subprocess

'''
#########################################################
#                                                       #
#                      FUNCTIONS                        #
#                                                       #
#########################################################
'''

def get_value_of_variable_in_input_file(file, string):
    try:
        content = [i.strip().split() for i in open(file).readlines()]
        index = [idx for idx, s in enumerate(content) if string in s][0]
        value = int(content[index][2])
        return value
    except IndexError:
        print('! String not in input file !')

'''
#########################################################
#                                                       #
#                   JOB INFORMATIONS                    #
#                                                       #
#########################################################
'''

# nTimesteps from input.dat
while True:
    try:
        nTimesteps = get_value_of_variable_in_input_file("./input.dat", "NTIME")
        break
    except FileNotFoundError:
        print('! No input file found !')

# Job script informations
jobscriptContent = open('jobscript_btrzx1', 'r').read().splitlines()

# Job name
jobNameIndex = [idx for idx, s in enumerate(jobscriptContent) if 'job-name=' in s][0]
jobName = jobscriptContent[jobNameIndex].split('job-name=', 1)[1]

# Job repetition
jobRepetitionIndex = [idx for idx, s in enumerate(jobscriptContent) if 'repeat=' in s][0]
jobRepetition = float(jobscriptContent[jobRepetitionIndex].split('repeat=', 1)[1])

# Number of required timesteps        
nTimestepsRequired = int(nTimesteps * jobRepetition)

## Walltime
walltimeIndex = [idx for idx, s in enumerate(jobscriptContent) if 'time=0-' in s][0]
walltime = jobscriptContent[walltimeIndex].split('time=0-', 1)[1].split(':')
walltimeSeconds = int(walltime[0])*60*60 + int(walltime[1])*60 + int(walltime[2])

## Sleep time until one job is completed = 24h + 5min
sleepTime = walltimeSeconds + 5*30

# Job status informations
jobHeader = ['JOBID', 'PARTITION', 'NAME', 'USER', 'ST', 'TIME', 'NODES', 'NODELIST(REASON)']
#jobStatus = subprocess.getoutput("squeue -u bt712347").strip().split() 


'''
#########################################################
#                                                       #
#                  START/RESTART JOB                    #
#                                                       #
#########################################################
'''

while True:
    
    jobStatusRunning = subprocess.getoutput("squeue --states=running -u bt712347").strip().split()
    jobStatusPending = subprocess.getoutput("squeue --states=pending -u bt712347").strip().split()
    
    # Check if no Job is running or pending and start/restart job
    if jobName not in jobStatusPending and jobName not in jobStatusRunning:
        subprocess.run(["sbatch", "jobscript_btrzx1"])
    
    # Check running jobs status and get job running time
    ## If job is not running than wait 30min
    while True:

        # Check Status if job is running
        if jobName in jobStatusRunning:

            # Job ID Running
            jobStatusRunningNameIndex = [idx for idx, s in enumerate(jobStatusRunning) if jobName in s][0]
            jobIDRunning = jobStatusRunning[jobStatusRunningNameIndex - 2]

            # Job time
            jobTime = jobStatusRunning[jobStatusRunningNameIndex + 3].split(':')
            jobTimeSeconds = int(jobTime[0])*60*60 + int(jobTime[1])*60 + int(jobTime[2])

            break
        # Wait 30min if job is pending
        elif jobName not in jobStatusRunning and jobName in jobStatusPending:
            sleep(30*60)
        
    # Sleep time  = 24h - job time + 5min
    sleepTimeCurrent = sleepTime - jobTimeSeconds
    sleep(sleepTimeCurrent)
    
    # read FDS.dat (restart file) to a list of lists
    ## If job is not existing than wait 30min time
    while True:
        try:
            nTimestepsCurrent = get_value_of_variable_in_input_file("./FDS.dat", 'FILE_COUNT')
            break
        except FileNotFoundError:
            sleep(30*60)
            
    # Check if gkw has run requiered timesteps
    if nTimestepsCurrent >= nTimestepsRequired:
        break
