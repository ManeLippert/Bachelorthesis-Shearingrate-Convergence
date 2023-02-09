#!/usr/bin/env python3

# AUTHOR ===================================================================================================================
# Manuel Lippert (GitHub: ManeLippert (https://github.com/ManeLippert))
# ==========================================================================================================================

# MODULES ==================================================================================================================

import datetime, time, os, sys, subprocess, math, argparse

# PARSER ===================================================================================================================

description_text = """

===================================== DESCRIPTION =====================================

FEATURES:
o NO REQUIREMENTS, runs with standard python3 libary
o Creats jobscript file with defined content (look into file itself for the jobscript)
o Start/Restarts job until criteria is suffused (default=0)
o Makes backup after each run before Restart and Restore files after fail
o Sends mail at the beginning, end and restart (default=False) with status file 
  as attachment (mailx has to be installed and working look into file for more info)
o Creates status file with current status and progress bars and updates it dynamically
  Progress: Total progress of simulation
  Run X:    Progress of current run

START SCRIPT IN BACKGROUND:

  o WITH NOHUP:
    Command:
    >>> nohup python3 -u slurm_monitor.py &> /dev/null &

    With arguments (example for 30000 timesteps):
    >>> nohup python3 -u slurm_monitor.py -n 30000 &> /dev/null &

    Output:
    >>> [1] 10537

    List Process:
    >>> ps ax | grep python3

    Output:
    >>> 10537 pts/1    S      0:00 python3 -u slurm_monitor.py
    >>> 23426 pts/1    S+     0:00 grep --color=auto slurm_monitor.py

    Kill Process:
    grep --color=auto slurm_monitor.py gets generated automatically no need to kill!
    >>> kill 10537
    
    Output after Enter or next Command:
    [1]+ Terminated              nohup python3 slurm_monitor.py &> /dev/null &
    
  o WITH SCREEN (has to be installed):
    Create Screen:
    >>> screen -S $SESSION
    
    Command:
    >>> python3 -u slurm_monitor.py --screen

    With arguments (example for 30000 timesteps):
    >>> python3 -u slurm_monitor.py --screen -n 30000
    
    Leave Screen:
    >>> ((Strg + a) + d)
    
    Enter Screen:
    >>> screen -r 
    or
    >>> screen -r $SESSION
    
    Kill Process:
    >>> ^C (Strg + c)
    or
    >>> (Strg + d) (kill Screen itself)

OUTPUT STATUS:
Just run the file will create an dynamic output (recommended with using screen) or
>>> cd $DATA && find . -name status.txt -exec cat {} \;

====================================== ARGUMENTS ======================================
"""

parser = argparse.ArgumentParser(description=description_text, formatter_class=argparse.RawTextHelpFormatter)
#parser._action_groups.pop()

required = parser.add_argument_group("required arguments")
additional = parser.add_argument_group("additional arguments")

additional.add_argument("-n", dest="timesteps", nargs="?", type=int, default=10000,
                        help="required timesteps                   (default=10000)")

additional.add_argument("--mail", dest="mail", nargs="?", type=str,
                        help="mail address (mail@server.de)        (default=None)")

additional.add_argument("--restart-mail", dest="restart", action="store_true",
                        help="mail after every restart             (default=False)")

additional.add_argument("-b", "--backup", dest="backup", nargs="?", type=str,
                        help="backup location for files            (default=None)\n"+
                             "- local (creates backup in simulation folder)")

additional.add_argument("-s", "--statusfile", dest="statusFile", nargs="?", type=str, default="status.txt",
                        help="file with output from nohup command  (default=status.txt)")

additional.add_argument("-r", "--restartfile", dest="restartFile", nargs="?", type=str, default="FDS.dat",
                        help="restart file with data               (default=FDS.dat)")

additional.add_argument("-j", "--job", dest="jobscriptFile", nargs="?", type=str, default="jobscript-create",
                        help="jobscript to run SLURM job           (default=jobscript-create)")

additional.add_argument("--job-name", dest="jobname", nargs="?", type=str, default="Name",
                        help="job name not longer than 8 character (default=Name)")

additional.add_argument("--ntask-per-node", dest="tasks", nargs="?", type=str, default="32",
                        help="MPI task per node                    (default=32)")

additional.add_argument("--nodes", dest="nodes", nargs="?", type=str, default="3",
                        help="number of nodes                      (default=3)")

additional.add_argument("--time", dest="walltime", nargs="?", type=str, default="0-24:00:00",
                        help="walltime of server (d-hh:mm:ss)      (default=0-24:00:00)")

additional.add_argument("--format", dest="formattable", nargs="?", type=str, default="universal",
                        help="format of output table               (default=universal)\n"+
                             "- fancy (round box)\n"+
                             "- universal (crossplattform)")

additional.add_argument("--refresh-rate", dest="sleepTime", nargs="?", type=int, default=300,
                        help="time interval to check status in sec (default=60)")

additional.add_argument("--screen", dest="screen", action="store_true",
                        help="activate output of script for screen (default=False)")

args = parser.parse_args()

# VARIABLES ================================================================================================================

outputCriteria = ["0", "Run successfully completed"]

slurmFiles = [f for f in os.listdir() if "slurm-" in f]
runCounter = len(slurmFiles)

currentTime = "00:00:00"

# PARSER VARIABLES =========================================================================================================

emailAddress = args.mail
restartMail = args.restart

backupLocation = args.backup

jobName = args.jobname
tasks = args.tasks
nodes = args.nodes
walltime = args.walltime

nTimestepsRequired = args.timesteps

jobscriptFilename = args.jobscriptFile
restartFilename = args.restartFile

statusFilename = args.statusFile
#statusFile = open(statusFilename, "r+")

formatTable = args.formattable
screen = args.screen

# Changing this value can cause problems in writing status file
sleepTime = args.sleepTime

def outputFilename(info):
    return "./slurm-" + info + ".out"

## JOBSCRIPT CONTENT =======================================================================================================

jobscriptContent = """#!/bin/bash -l

# jobname
#SBATCH --job-name=""" + jobName + """

# MPI tasks
#SBATCH --ntasks-per-node=""" + tasks + """

# number of nodes
#SBATCH --nodes=""" + nodes + """

# walltime
#              d-hh:mm:ss
#SBATCH --time=""" + walltime + """

# execute the job
time mpirun -np $SLURM_NTASKS ./gkw.x > output.dat

# end
exit 0
"""

jobName = jobName[0:8]

## FLAGS ===================================================================================================================

startOutputFlag = "Submitted batch job "
restartFlag = "FILE_COUNT"

## COMMANDS ================================================================================================================

commandJobRunning = "squeue --states=running --name " + jobName
commandJobPending = "squeue --states=pending --name " + jobName

commandJobStarting = "sbatch"

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
                prefix="", 
                progress_fill="=", progress_fill_top="", progress_fill_bot="",
                progress_unfill=".", 
                progress_bracket=["[","]"]):
    
    x = int(barsize*current_value/required_value)
    percent = int(100*current_value/required_value)
    
    try:
        percent_format = "  (" + (int(math.log10(100)) - int(math.log10(percent)))*" " + "{}%)"
    except ValueError:
        percent_format = "  (" + int(math.log10(100))*" " + "{}%)"
        
    try:
        ratio_format = "  " + (int(math.log10(required_value)) - int(math.log10(current_value)))*" " + "{}/{}"
    except ValueError:
        ratio_format = "  " + int(math.log10(required_value))*" " + "{}/{}"

        
    bar_format =  "{}  " + progress_bracket[0] + progress_fill_bot + "{}" + progress_fill_top + "{}" + progress_bracket[1] + percent_format + ratio_format
    bar = bar_format.format(prefix, progress_fill*x, progress_unfill*(barsize-x), percent, current_value, required_value)
        
    return bar

## OUTPUT TABLE ============================================================================================================

def message(string, add):
    string = string + add + "\n"
    return string

def print_table_row(content,
                    current_value, required_value,
                    run_conter, current_time, required_time = walltime,
                    delete_line_index = -10, table_width = 80,
                    output_type = None, time_info = True):
    
    msg=""
    table_inner_width = table_width - 4
    
    if formatTable == "fancy":
        table_outline = ["╭─", "─╮", "╰─", "─╯", "├─", "─┤", "│ ", " │", "─"]
    
    if formatTable == "universal":
        table_outline = ["+-", "-+", "+-", "-+", "+-", "-+", "| ", " |", "-"]
    
    sep_top = table_outline[0] + table_inner_width*table_outline[8] + table_outline[1]
    sep_mid = table_outline[4] + table_inner_width*table_outline[8] + table_outline[5]
    sep_end = table_outline[2] + table_inner_width*table_outline[8] + table_outline[3]
    
    row_cols = [2, 10, table_inner_width - 47, 13, 11, 13, 2]
    row_format = "".join(["{:<" + str(col) + "}" for col in row_cols])
    
    progress_cols = [2, table_inner_width, 2]
    progress_format = "".join(["{:<" + str(col) + "}" for col in progress_cols])
    
    progressbar_content = [progressbar(required_value, current_value, progress_fill_top=">",
                                       prefix="PROGRESS")]
    progressbar_content.insert(0, table_outline[6])
    progressbar_content.insert(len(progressbar_content), table_outline[7])
    
    required_time = get_time_in_seconds(required_time)
    current_time = get_time_in_seconds(current_time)
    
    jobStatus = subprocess.getoutput("squeue --name " + jobName).strip().split("\n")
        
    jobStatusHeader = [" " + jobStatus[0]]
    jobStatusHeader.insert(0, table_outline[6])
    jobStatusHeader.insert(len(jobStatusHeader), table_outline[7])
        
    try:
        jobStatusInfo = [jobStatus[1][12:12+table_inner_width]]
    except IndexError:
        jobStatusInfo = [table_inner_width*" "]
    jobStatusInfo.insert(0, table_outline[6])
    jobStatusInfo.insert(len(jobStatusInfo), table_outline[7])
    
    try:
        progressbartime_content = [progressbar(required_time, current_time, progress_fill_top=">",
                                               prefix="RUN " + str(run_conter) + (2 - int(math.log10(run_conter)))*" " + " ")]
    except ValueError:
        progressbartime_content = [progressbar(required_time, current_time, progress_fill_top=">",
                                               prefix="RUN " + str(run_conter) + 2*" " + " ")]
    progressbartime_content.insert(0, table_outline[6])
    progressbartime_content.insert(len(progressbartime_content), table_outline[7])
    
    content.insert(0, table_outline[6])
    if time_info:
        if output_type == "header":
            content.append("DATE")
            content.append("TIME")
            content.append("W:DD:HH:MM:SS")
        else:
            content.append(time_date())
            content.append(time_time())
            content.append(time_duration(startTime, pastTime))
            
    content.insert(len(content), table_outline[7])
    
    if output_type == "header":
        delete_line_index = -1
        
        msg += sep_top + "\n"
        msg += row_format.format(*content) + "\n"
        msg += sep_end + "\n"
        
    elif output_type == "middle":
        if screen:
            sys.stdout.write("\x1b[1A"*(-delete_line_index + 1))
        
        msg += sep_mid + "\n"
        msg += row_format.format(*content) + "\n"
        msg += sep_end + "\n"
        
    elif output_type == "update":
        if screen:
            sys.stdout.write("\x1b[1A"*(-delete_line_index + 1))
        
        msg += sep_end + "\n"
        
    else:
        if screen:
            sys.stdout.write("\x1b[1A"*(-delete_line_index + 1))
        
        msg += row_format.format(*content) + "\n"
        msg += sep_end + "\n"
        
    msg += " "*table_width + "\n" + " "*table_width + "\n"
    
    msg += sep_top + "\n"
    msg += progress_format.format(*jobStatusHeader) + "\n"
    msg += progress_format.format(*jobStatusInfo) + "\n"   
    msg += sep_mid + "\n"
    msg += progress_format.format(*progressbar_content) + "\n"
    msg += progress_format.format(*progressbartime_content) + "\n"
    msg += sep_end + "\n"
    
    if screen:
        print(msg, flush=True)
    
    delete_write_line_to_file(statusFilename, msg, end=delete_line_index)

## INFORMATIONS ============================================================================================================

def get_value_of_variable_from_file(file, file_index, relative_index, string):
    try:
        content = [i.strip().split() for i in open(file).readlines()]
        index = [idx for idx, s in enumerate(content) if string in s][file_index]
        value = content[index][relative_index]
        return value
    except IndexError:
        print_table_row(["ERROR", "String not found in file"], 
                        nTimestepsRequired, nTimestepsRequired, runCounter, currentTime)
        quit()

def find_string_in_file(file, string):
    
    with open(file) as f:
        if string in f.read():
            return True
        else:
            return False
        
## FILE ====================================================================================================================

def write_add_string_into_file(filename, substring, add, comment = None):

    with open(filename, "r") as file:
        data = file.readlines()

    try:
        index = [idx for idx, s in enumerate(data) if substring in s][0]
        data[index] = substring + add + "\n"
        
    except IndexError:
        index = [idx for idx, s in enumerate(data) if "\n" in s][1]
        data.insert(index, "\n")
        data.insert(index + 1, comment)
        data.insert(index + 2, substring + add + "\n")
        
    with open(filename, "w") as file:
        file.writelines(data)
        file.flush()
        
def delete_write_line_to_file(filename, add = "", start=None, end=None):
    
    with open(filename, "r") as file:
    
        try:
            lines = file.readlines()[start:end]
            content = "".join(lines)
            
        except IndexError:
            pass
        
    content += add
    
    write_file(filename, content)
    
def write_file(filename, content):
    
    with open(filename, "w") as file:
        file.write(content)
        file.flush()

## TIME ====================================================================================================================

def format_num(time):
    if time < 10:
        return "0" + str(time)
    else:
        return str(time)
    
def get_time_as_string(sec):
    
    mins,  sec   = sec   // 60, sec   % 60
    hours, mins  = mins  // 60, mins  % 60
    days,  hours = hours // 24, hours % 24
    weeks, days  = days  // 7,  days  % 7
    
    timeConvertedString = (str(int(weeks)) + ":" + 
                     format_num(int(days)) + ":" + format_num(int(hours)) + ":" + 
                     format_num(int(mins)) + ":" + format_num(int(sec)))
    
    return timeConvertedString

def get_time_in_seconds(time):
    
    # Format D-HH:MM:SS or HH:MM:SS or MM:SS
    
    time = time.replace("-", ":")
    time_split = time.split(":")

    seconds = [24*60*60, 60*60, 60, 1]
    seconds = seconds[-len(time_split):]

    time_sec = sum([a*b for a,b in zip(seconds, map(int,time_split))])

    return time_sec

def time_date():
    e = datetime.datetime.now()
    return "%s-%s-%s" % (format_num(e.year), format_num(e.month), format_num(e.day))

def time_time():
    e = datetime.datetime.now()
    return "%s:%s:%s" % (format_num(e.hour), format_num(e.minute), format_num(e.second))

def time_duration(startTime, pastTime):
    stop = time.time()
    return get_time_as_string(stop - startTime + pastTime)

def get_time_from_statusfile(filename, line_index):
    with open(filename, "r") as file:
        line = file.readlines()[line_index]
        
        content = line.split(" ")
        time = content[-2]
        time_sec = get_time_in_seconds(time)
        
        return time_sec

## STATUS ==================================================================================================================

def get_job_status():
    
    jobStatusRunning = subprocess.getoutput(commandJobRunning).strip().split()
    jobStatusPending = subprocess.getoutput(commandJobPending).strip().split()

    return jobStatusRunning, jobStatusPending

def set_output_type():
    jobStatusRunning, jobStatusPending = get_job_status()

    if jobName in jobStatusRunning:
        outputType = "running"
    elif jobName in jobStatusPending:
        outputType = "pending"
    else:
        outputType ="no Output"

    return outputType

## MAIL ====================================================================================================================

def send_mail(recipient, subject, body = None):

    recipient = recipient.encode("utf_8")
    
    subject = subject.replace(" ", "_")
    subject = subject.encode("utf_8")
    
    if body == None:
        body = "For futher information open attachment"
    
    body = body.encode("utf_8")
    
    attachmentPath = folder + "/" + statusFilename
    attachment = attachmentPath.encode("utf_8")

    process = subprocess.Popen(["ssh", "master", "/usr/bin/mailx", "-s", subject, "-a",  attachment, recipient],
                               stdin=subprocess.PIPE)
    process.communicate(body)


# JOB INIT =================================================================================================================

user = os.getlogin()
startTime = time.time()

if not os.path.isfile(statusFilename):
    partTime = 0
    
    write_file(statusFilename,"")

    print_table_row(["OUTPUT", "INFO"], 
                    0, nTimestepsRequired, runCounter, currentTime, 
                    output_type="header")
else:
    partTime = get_time_from_statusfile(statusFilename, -11)

folder = os.getcwd()
path = folder.split(user + "/")[1]

## BACKUP PATH =============================================================================================================

if backup:
    # If local has been specified as backupLocation, then the backup is copied to the same directory as the simulation folder 
    # (simFolder) is located in. The backup is copied to a directory with name simFolder + "-backup". 
    if backupLocation == "local":
        simFolder = path.split("/")[-1]
        backupPath = folder + "/../" + simFolder + "-backup"
     
    # Otherwise, use the backupLocation parsed as argument.
    else:
        if backupLocation[-1] != "/":
            backupLocation += "/"
        backupPath = backupLocation + path
    
    if not os.path.exists(backupPath):
        os.makedirs(backupPath)

# START/RESTART JOB ========================================================================================================

outputType = set_output_type()

## BEGIN ===================================================================================================================

while True:
    try:
        nTimestepsCurrent = int(get_value_of_variable_from_file("./" + restartFilename, 0, 2, restartFlag))
        
        if nTimestepsCurrent >= nTimestepsRequired:
            print_table_row(["CONTROL", "Current Timesteps " + str(nTimestepsCurrent)], 
                            nTimestepsCurrent, nTimestepsRequired, runCounter, currentTime)
            
            print_table_row(["SUCCESS", "Stop monitoring " + jobName], 
                            nTimestepsRequired, nTimestepsRequired, runCounter, walltime, 
                            output_type="middle")
            
            if emailNotification:
                send_mail(emailAddress, "Ended Job " + jobName)

            quit()
        
        # Continue
        else:
            print_table_row(["CONTINUE", "Continue monitoring " + jobName], 
                            nTimestepsCurrent, nTimestepsRequired, runCounter, currentTime, 
                            output_type="middle")
            
            jobStatusRunning, jobStatusPending = get_job_status()
            
            if outputType == "running":
                #runCounter += 1
                
                jobStatusRunningNameIndex = [idx for idx, s in enumerate(jobStatusRunning) if jobName in s][0]
                currentTime = jobStatusRunning[jobStatusRunningNameIndex + 3]
                
                print_table_row(["CONTROL", "Current Timesteps " + str(nTimestepsCurrent)], 
                                nTimestepsCurrent, nTimestepsRequired, runCounter, currentTime)
            
            elif outputType == "pending":
                #runCounter += 1
                
                jobStatusPendingNameIndex = [idx for idx, s in enumerate(jobStatusPending) if jobName in s][0]
                currentTime = jobStatusPending[jobStatusPendingNameIndex + 3]
                
                print_table_row(["CONTROL", "Current Timesteps " + str(nTimestepsCurrent)], 
                                nTimestepsCurrent, nTimestepsRequired, runCounter, currentTime)
            
            else:
                if backup:
                    print_table_row(["BACKUP", backupLocation],
                                    nTimestepsCurrent, nTimestepsRequired, runCounter, currentTime)
                    subprocess.run(["rsync", "-a", "", backupPath])
            
            if emailNotification:
                send_mail(emailAddress, "Continued Job " + jobName)
            break
    
    # Start    
    except FileNotFoundError:
        nTimestepsCurrent = 0
        
        print_table_row(["STARTING", "Start monitoring " + jobName], 
                        nTimestepsCurrent, nTimestepsRequired, runCounter, currentTime, 
                        output_type="middle")
        
        if backup:
            print_table_row(["BACKUP", backupLocation],
                            nTimestepsCurrent, nTimestepsRequired, runCounter, currentTime)
            subprocess.run(["rsync", "-a", "", backupPath])
            
        if emailNotification:
            send_mail(emailAddress, "Started Job " + jobName)
        break

## MONITOR ROUTINE =========================================================================================================

while True:
    
    jobStatusRunning, jobStatusPending = get_job_status()

    # Job running
    if jobName in jobStatusRunning:
        
        jobStatusRunningNameIndex = [idx for idx, s in enumerate(jobStatusRunning) if jobName in s][0]
        jobID = jobStatusRunning[jobStatusRunningNameIndex - 2]
        
        currentTime = jobStatusRunning[jobStatusRunningNameIndex + 3]
        
        if outputType == "running":
            print_table_row(["RUNNING", "Job is executed"], 
                            nTimestepsCurrent, nTimestepsRequired, runCounter, currentTime)
            outputType = "no Output"
        else:
            print_table_row(["RUNNING", "Job is executed"], 
                            nTimestepsCurrent, nTimestepsRequired, runCounter, currentTime, 
                            output_type="update")
            
        time.sleep(sleepTime)
        
    # Job pending
    elif jobName in jobStatusPending:

        jobStatusPendingNameIndex = [idx for idx, s in enumerate(jobStatusPending) if jobName in s][0]
        jobID = jobStatusPending[jobStatusPendingNameIndex - 2]
        
        currentTime = jobStatusPending[jobStatusPendingNameIndex + 3]
        
        if outputType == "pending":
            print_table_row(["WAITING", "Job is pending"], 
                            nTimestepsCurrent, nTimestepsRequired, runCounter, currentTime)
            outputType = "running"
        else:
            print_table_row(["WAITING", "Job is pending"], 
                            nTimestepsCurrent, nTimestepsRequired, runCounter, currentTime, 
                            output_type="update")
            
        time.sleep(sleepTime)

    # Job start/restart
    else:
        
        # Check error and making Backup
        while True:
            try:
                outputContent = open(outputFilename(jobID)).readlines()[-5].replace("\n","")
                
                # If scan of output.dat is needed: Scans for string "Run Successful in output.dat and returns bool value"
                #runSuccess = find_string_in_file("output.dat", outputCriteria[1])
                #if outputCriteria[0] in outputContent and runSuccess: 
                
                if outputCriteria[0] in outputContent: 
                    
                    if backup:
                        print_table_row(["BACKUP", backupLocation],
                                        nTimestepsCurrent, nTimestepsRequired, runCounter, currentTime)
                        subprocess.run(["rsync", "-a", "", backupPath])
                    
                    break
                else:
                    print_table_row(["ERROR", "SLURM Job failed to execute"], 
                                    nTimestepsCurrent, nTimestepsRequired, runCounter, currentTime)
                        
                    if backup:
                        print_table_row(["RESTORE", backupLocation], 
                                        nTimestepsCurrent, nTimestepsRequired, runCounter, currentTime)
                        subprocess.run(["rsync", "-a", "-I", "--exclude=status.txt", backupPath + "/", ""])
                    
                    break
                    
            except (IndexError, FileNotFoundError):
                time.sleep(30)
            except NameError:
                break
        
        try:
            nTimestepsCurrent = int(get_value_of_variable_from_file("./" + restartFilename, 0, 2, restartFlag))
            print_table_row(["CONTROL", "Current Timesteps " + str(nTimestepsCurrent)], 
                            nTimestepsCurrent, nTimestepsRequired, runCounter, currentTime)
            
            if nTimestepsCurrent >= nTimestepsRequired:
                print_table_row(["SUCCESS", "Stop monitoring " + jobName], 
                                nTimestepsRequired, nTimestepsRequired, runCounter, walltime)
                
                if emailNotification:
                    send_mail(emailAddress, "Ended Job " + jobName)
                break
                
        except FileNotFoundError:
            pass   
        
        # Create jobscript
        if jobscriptFilename == "jobscript-create":
            jobscriptFilename = "jobscript"
            write_file(jobscriptFilename, jobscriptContent)
        
        # Start Job
        startOutput = subprocess.check_output([commandJobStarting, jobscriptFilename]).decode("utf-8").replace("\n", "")
        jobID = startOutput.split(startOutputFlag)[1]
        
        runCounter += 1
        
        print_table_row(["STARTING", startOutput], 
                        nTimestepsCurrent, nTimestepsRequired, runCounter, currentTime, 
                        output_type="middle")
        
        try:
            if restartMail and emailNotification:
                send_mail(emailAddress, "Restart Job " + jobName)
        except NameError:
            continue
        
        time.sleep(30)
        
        outputType = set_output_type()
        
## RESTART =================================================================================================================