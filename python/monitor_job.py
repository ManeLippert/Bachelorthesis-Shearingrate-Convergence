#!/usr/bin/env python3

# AUTHOR ===================================================================================================================
# Manuel Lippert (GitHub: ManeLippert (https://github.com/ManeLippert))
# ==========================================================================================================================

# FEATURES =================================================================================================================
# o Restarts job until criteria is suffused 
# o Makes backup after each run before Restart
# o Sends mail at the beginning, end and restart with status file as attachment
# ==========================================================================================================================

# START SCRIPT =============================================================================================================
# Command (Script runs in the background):
# >>> nohup python3 -u monitor_job.py &> status.txt &
#
# Output:
# >>> [1] 10537
#
# This will write every output in the status file 'status.txt'
# ==========================================================================================================================

# LIST PROCESS =============================================================================================================
# Command:
# >>> ps ax | grep monitor_job.py
#
# Output:
# >>> 10537 pts/1    S      0:00 python3 -u monitor_job.py
# >>> 23426 pts/1    S+     0:00 grep --color=auto monitor_job.py
# 
# KILL PROCESS =============================================================================================================
# Command:
# >>> kill 10537
# ==========================================================================================================================


# OUTPUT STATUS ============================================================================================================
# Command:
# >>> cd $DATA && find . -name status.txt -exec tail --lines=+10 {} \;
# ==========================================================================================================================



# MODULES ==================================================================================================================

import datetime, time, os, subprocess

# VARIABLES ================================================================================================================

emailAddress = 'Manuel.Lippert@uni-bayreuth.de'
backupLocation = '/scratch/bt712347/backup'
jobName = 'Name'

nTimestepsRequired = 20000
outputCriteria = '0'

sleepTime = 5*60

restartMail = False

## FILENAMES ===============================================================================================================

jobscriptFilename = 'jobscript'
restartFilename = 'FDS.dat'
monitorFilename = 'status.txt'
inputFilename = 'input.dat'

def outputFilename(info):
    return './slurm-' + info + '.out'

## FLAGS ===================================================================================================================

jobNameFlag = '#SBATCH --job-name='
startOutputFlag = 'Submitted batch job '
restartFlag = 'FILE_COUNT'

## COMMANDS ================================================================================================================

commandJobRunning = 'squeue --states=running -u '
commandJobPending = 'squeue --states=pending -u '

commandJobStarting = 'sbatch'

# FUNCTIONS ================================================================================================================

## OUTPUT TABLE ============================================================================================================

def print_table_row(content, output_type = None, time_info = True):
    
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

## INFORMATIONS ============================================================================================================

def get_value_of_variable_from_file(file, file_index, relative_index, string):
    try:
        content = [i.strip().split() for i in open(file).readlines()]
        index = [idx for idx, s in enumerate(content) if string in s][file_index]
        value = content[index][relative_index]
        return value
    except IndexError:
        print_table_row(['ERROR', 'String not found in file'], output_type='end')
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

def get_job_status(command, user):
    
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


# JOB INFORMATIONS =========================================================================================================

jobInformations = []
startTime = time.time()
user = os.getlogin()

print_table_row(['OUTPUT', 'JOB INITIALIZE'], output_type='header')


folder = os.path.dirname(os.path.abspath(__file__))
path = os.path.dirname(os.path.abspath(__file__)).split(user + '/')[1]

## MAIL SWITCH =============================================================================================================

try:
    emailNotification = True
    jobInformations.append(['INFO', 'E-Mail', emailAddress])
except NameError:
    emailNotification = False

## JOB SCRIPT ==============================================================================================================

if os.path.exists(jobscriptFilename):    
    print_table_row(['SUCCESS', 'Found ' + jobscriptFilename])
else:
    print_table_row(['ERROR', 'No file "' + jobscriptFilename + '" found'])
    
    if emailNotification:
        send_mail(emailAddress, 'Failed Job ' + jobName)
    quit()

jobscriptContent = open(jobscriptFilename, 'r').read().splitlines()

## JOB NAME ================================================================================================================

write_add_string_into_file(jobscriptFilename, jobNameFlag, jobName)

if len(jobName) > 8:
    jobName = jobName[0:8]

jobInformations.append(['INFO', 'Name', jobName])

## TIMESTEPS ===============================================================================================================

jobInformations.append(['INFO', 'Required Timesteps', nTimestepsRequired])

## BACKUP PATH/SWITCH ======================================================================================================

try:
    
    if backupLocation[-1] != '/':
        backupLocation += '/'
    backupPath = backupLocation + path
    
    if not os.path.exists(backupPath):
        os.makedirs(backupPath)
    
    backup = True
    jobInformations.append(['INFO', 'Backup Path', backupLocation])
    
except NameError:
    backup = False

print_table_row(['OUTPUT', 'JOB INFORMATIONS', 'VALUE'], output_type='header', time_info=False)
print_table_row(jobInformations, output_type='end', time_info=False)

# START/RESTART JOB ========================================================================================================

print_table_row(['OUTPUT', 'JOB MONITORING'], output_type='header')

outputType = set_output_type(user)

## BEGIN ===================================================================================================================

while True:
    try:
        nTimestepsCurrent = int(get_value_of_variable_from_file('./' + restartFilename, 0, 2, restartFlag))
        
        if outputType == 'running':
            print_table_row(['CONTROL', 'Current Timesteps ' + str(nTimestepsCurrent)])
        
        if nTimestepsCurrent >= nTimestepsRequired:
            print_table_row(['SUCCESS', 'Stop monitoring'], output_type='end')
            
            if emailNotification:
                send_mail(emailAddress, 'Ended Job ' + jobName)

            quit()
        
        # Continue
        else:
            print_table_row(['CONTINUE', 'Continue monitoring'])
            
            if emailNotification:
                send_mail(emailAddress, 'Continued Job ' + jobName)
            break
    
    # Start    
    except FileNotFoundError:
        print_table_row(['STARTING', 'Start monitoring'])
        
        if backup:
            print_table_row(['BACKUP', backupLocation], output_type='end')
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
        
        if outputType == 'running':
            print_table_row(['RUNNING', 'Job is executed'])
            outputType = 'no Output'
            
        time.sleep(sleepTime)
        
    # Job pending
    elif jobName in jobStatusPending:

        if outputType == 'pending':
            print_table_row(['WAITING', 'Job is pending'])
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
                        print_table_row(['BACKUP', backupLocation], output_type='end')
                        subprocess.run(['rsync', '-a', '', backupPath])
                    
                    break
                else:
                    print_table_row(['ERROR', outputContent[0:27]])
                        
                    if backup:
                        print_table_row(['RESTORE', backupLocation], output_type='end')
                        subprocess.run(['rsync', '-a', '-I', '--exclude=status.txt', backupPath + '/', ''])
                    
                    break
                    
            except (IndexError, FileNotFoundError):
                time.sleep(30)
            except NameError:
                break
        
        try:
            nTimestepsCurrent = int(get_value_of_variable_from_file('./' + restartFilename, 0, 2, restartFlag))
            print_table_row(['CONTROL', 'Current Timesteps ' + str(nTimestepsCurrent)])
            
            if nTimestepsCurrent >= nTimestepsRequired:
                print_table_row(['SUCCESS', 'Stop monitoring'], output_type='end')
                
                if emailNotification:
                    send_mail(emailAddress, 'Ended Job ' + jobName)
                break
                
        except FileNotFoundError:
            pass   
        
        # Start Job
        startOutput = subprocess.check_output([commandJobStarting, jobscriptFilename]).decode('utf-8').replace('\n', '')
        jobID = startOutput.split(startOutputFlag)[1]
        
        print_table_row(['STARTING', startOutput])
        
        try:
            if restartMail and emailNotification:
                send_mail(emailAddress, 'Restart Job ' + jobName)
        except NameError:
            continue
        
        time.sleep(30)
        
        outputType = set_output_type(user)
        
## RESTART =================================================================================================================