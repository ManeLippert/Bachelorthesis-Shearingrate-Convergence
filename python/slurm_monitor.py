#!/usr/bin/env python3

# AUTHOR ===================================================================================================================
# Manuel Lippert (GitHub: ManeLippert (https://github.com/ManeLippert))
# ==========================================================================================================================

# MODULES ==================================================================================================================

import datetime, time, os, sys, subprocess, math, argparse, pkg_resources

# PARSER ===================================================================================================================

description_text = """

===================================== DESCRIPTION =====================================

FEATURES:
o NO REQUIREMENTS: Default script runs with standard python3 library
o Creates jobscript file with defined content 
  (look into variable "jobscriptContent" to add more option)
o Start/Restarts job until criteria is suffused (default=10000)
o Makes backup after each run before Restart and Restore files after failed run
  (uses rsync command line utility)
o Reset option after run fails and dump files were written
  (rely on "h5py" & "pandas" & "numpy" which get installed by the script)
o Creates status file with current status and progress bars and updates it dynamically
  Progress: Total progress of simulation
  Run X   : Progress of current run
o Sends mail at the beginning, end and restart (default=False) with status file 
  as attachment 
  (mailx has to be installed, look into "send_mail" function for more info)

START SCRIPT IN BACKGROUND:

  o WITH NOHUP:
    Command:
    >>> nohup python3 -u slurm_monitor.py --job-name $JOBNAME &> /dev/null &

    Kill Process:
    >>> python3 -u slurm_monitor.py --job-name $JOBNAME --kill
    
  o WITH SCREEN (has to be installed):
    Create Screen:
    >>> screen -S $SESSION
    
    Command:
    >>> python3 -u slurm_monitor.py --job-name $JOBNAME --verbose
    
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
>>> cd $DATA && find . -name $STATUSFILENAME -exec tail -8 {} \;

====================================== ARGUMENTS ======================================
"""

parser = argparse.ArgumentParser(description=description_text, formatter_class=argparse.RawTextHelpFormatter)
#parser._action_groups.pop()

required = parser.add_argument_group("required arguments")

required.add_argument("-j", "--job-name", dest="jobname", nargs="?", type=str, required= True,
                      help="job name not longer than 8 character")

additional = parser.add_argument_group("additional arguments")

additional.add_argument("-n", dest="timesteps", nargs="?", type=int, default=10000,
                        help="required timesteps                   (default=10000)")

additional.add_argument("-r", "--reset", dest="reset", action="store_true",
                        help="Uses Dumpfiles to reset Simulation   (default=False)")

additional.add_argument("-v", "--verbose", dest="verbose", action="store_true",
                        help="activate output of script            (default=False)")

additional.add_argument("-b", "--backup", dest="backup", nargs="?", type=str,
                        help="backup location for files            (default=None)\n"+
                             "- local (creates backup in simulation folder)\n"+
                             "- home  (creates backup in home folder)")

additional.add_argument("--jobscript", dest="jobscriptFile", nargs="?", type=str, default="jobscript-create",
                        help="jobscript to run SLURM job           (default=jobscript-create)")

additional.add_argument("--ntask-per-node", dest="tasks", nargs="?", type=str, default="32",
                        help="MPI task per node                    (default=32)")

additional.add_argument("--nodes", dest="nodes", nargs="?", type=str, default="3",
                        help="number of nodes                      (default=3)")

additional.add_argument("--walltime", dest="wallTime", nargs="?", type=str, default="0-24:00:00",
                        help="walltime of server (d-hh:mm:ss)      (default=0-24:00:00)")

additional.add_argument("--mail", dest="mail", nargs="?", type=str,
                        help="mail address (mail@server.de)        (default=None)")

additional.add_argument("--restart-mail", dest="restartmail", action="store_true",
                        help="mail after every restart             (default=False)")

additional.add_argument("--statusfile", dest="statusFile", nargs="?", type=str, default="status.txt",
                        help="file with output from nohup command  (default=status.txt)")

additional.add_argument("--restartfile", dest="restartFile", nargs="?", type=str, default="FDS.dat",
                        help="restart file with data               (default=FDS.dat)")

additional.add_argument("--format", dest="formattable", nargs="?", type=str,
                        help="format of output table               (default=None)\n"+
                             "- fancy (round box)\n"+
                             "- universal (crossplattform)\n"+
                             "- None (No frame)")

additional.add_argument("--refresh-rate", dest="sleepTime", nargs="?", type=int, default=300,
                        help="time interval to check status in sec (default=300)")

additional.add_argument("--kill", dest="kill", action="store_true",
                        help="kills monitor process                (default=False)")

args = parser.parse_args()

# VARIABLES ================================================================================================================

outputCriteria = ["0", "Run successfully completed"]

slurmFiles = [f for f in os.listdir() if "slurm-" in f]
runCounter = len(slurmFiles)

nTimestepsCurrent = 0

currentTime = "00:00:00"

dataFilename = "gkwdata.h5"

check1Filename, check2Filename = "DM1.dat", "DM2.dat"
check1Bin, check2Bin = check1Filename.replace(".dat", ""), check2Filename.replace(".dat", "")

# PARSER VARIABLES =========================================================================================================

emailAddress = args.mail
RESTARTMAIL = args.restartmail

backupLocation = args.backup

jobName = args.jobname
tasks = args.tasks
nodes = args.nodes
wallTime = args.wallTime

nTimestepsRequired = args.timesteps

jobscriptFilename = args.jobscriptFile
restartFilename = args.restartFile
restartBin = restartFilename.replace(".dat", "")

statusFilename = args.statusFile
#statusFile = open(statusFilename, "r+")

formatTable = args.formattable
VERBOSE = args.verbose
RESET = args.reset

# Kill process of monitoring
kill = args.kill

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
#SBATCH --time=""" + wallTime + """

# execute the job
time mpirun -np $SLURM_NTASKS ./gkw.x > output.dat

# end
exit 0
"""

jobName = jobName[0:8]

## FLAGS ===================================================================================================================

startOutputFlag = "Submitted batch job "
restartFlag = "FILE_COUNT"

# SWITCHES =================================================================================================================

if emailAddress==None:
    EMAIL = False
else:
    EMAIL = True

if backupLocation==None:
    BACKUP = False
else:
    BACKUP = True
    
## PATHS ===================================================================================================================

user = os.getlogin()
folder = os.getcwd()
path = folder.split(user + "/")[1]

# Set backup location
if BACKUP:
    # If local has been specified as backupLocation, then the backup is copied to the same directory as the simulation folder 
    # (simFolder) is located in. The backup is copied to a directory with name simFolder + "-backup". 
    if backupLocation == "local":
        simFolder = path.split("/")[-1]
        backupPath = folder + "/../" + simFolder + "-backup"
     
    # Create backup in home folder    
    elif backupLocation == "home":
        simFolder = path.split("/")[-1]
        backupPath = "~/" + simFolder + "-backup"
     
    # Otherwise, use the backupLocation parsed as argument.
    else:
        if backupLocation[-1] != "/":
            backupLocation += "/"
        backupPath = backupLocation + path
    
    if not os.path.exists(backupPath):
        os.makedirs(backupPath)
        
## COMMANDS ================================================================================================================

commandJobRunning = "squeue --states=running --name " + jobName
commandJobPending = "squeue --states=pending --name " + jobName

commandMonitorKill = "ps ax | grep " + jobName + " | awk '{print $1}'"

commandJobStarting = "sbatch"

if BACKUP:
    commandBackup  = ["rsync", "-a", "", backupPath]
    commandRestore = ["rsync", "-a", "-I", "--exclude=" + statusFilename, backupPath + "/", ""]

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
                    delete_line_index = -10, table_width = 80,
                    output_type = None, 
                    TIMEINFO = True, WRITEFILE=True):
    
    current_value, required_value = nTimestepsCurrent, nTimestepsRequired,
    run_conter = runCounter
    current_time, required_time = currentTime, wallTime
    
    msg=""
    table_inner_width = table_width - 4
    table_content_width = table_inner_width - 47
    
    if formatTable == "fancy":
        table_outline = ["╭─", "─╮", "╰─", "─╯", "├─", "─┤", "│ ", " │", "─"]
    
    if formatTable == "universal":
        table_outline = ["+-", "-+", "+-", "-+", "+-", "-+", "| ", " |", "-"]
        
    if formatTable == None:
        table_outline = ["--", "--", "--", "--", "--", "--", "  ", "  ", "-"]
    
    sep_top = table_outline[0] + table_inner_width*table_outline[8] + table_outline[1]
    sep_mid = table_outline[4] + table_inner_width*table_outline[8] + table_outline[5]
    sep_end = table_outline[2] + table_inner_width*table_outline[8] + table_outline[3]
    
    row_cols   = [2, 10, table_content_width, 13, 11, 13, 2]
    row_format = "".join(["{:<" + str(col) + "}" for col in row_cols])
    
    progress_cols   = [2, table_inner_width, 2]
    progress_format = "".join(["{:<" + str(col) + "}" for col in progress_cols])
    
    progressbar_content = [progressbar(required_value, current_value, progress_fill_top=">", 
                                       prefix="PROGRESS")]
    progressbar_content.insert(0, table_outline[6])
    progressbar_content.insert(len(progressbar_content), table_outline[7])
    
    required_time = get_time_in_seconds(required_time)
    current_time  = get_time_in_seconds(current_time)
    
    jobStatus = subprocess.getoutput("squeue --name " + jobName).strip().split("\n")
        
    try:
        jobStatusHeader = ["  " + jobStatus[0]]
        jobStatusInfo   = [jobStatus[1][11:11+table_inner_width]]
    except IndexError:
        jobStatus_cols   = [12, 12, 10, table_inner_width- 71, 13, 11, 13]
        jobStatus_format = "".join(["{:<" + str(col) + "}" for col in jobStatus_cols])
        
        jobStatusHeader = ["OUTPUT", "NAME", "USER", "" ,"DATE", "TIME", "W:DD:HH:MM:SS"]
        jobStatusInfo   = [content[0], jobName, user, "", time_date(), time_time(), time_duration(startTime, pastTime)]

        jobStatusHeader   = [jobStatus_format.format(*jobStatusHeader)]
        jobStatusInfo     = [jobStatus_format.format(*jobStatusInfo)]
        
    jobStatusHeader.insert(0, table_outline[6])
    jobStatusHeader.insert(len(jobStatusHeader), table_outline[7])
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
    
    content[1] = content[1][-(table_content_width-1):]
    content.insert(0, table_outline[6])
    if TIMEINFO:
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
        if VERBOSE:
            sys.stdout.write("\x1b[1A"*(-delete_line_index + 1))
        
        msg += sep_mid + "\n"
        msg += row_format.format(*content) + "\n"
        msg += sep_end + "\n"
        
    elif output_type == "update":
        if VERBOSE:
            sys.stdout.write("\x1b[1A"*(-delete_line_index + 1))
        
        msg += sep_end + "\n"
        
    else:
        if VERBOSE:
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
    
    if VERBOSE:
        print(msg, flush=True)
    if WRITEFILE:
        delete_write_line_to_file(statusFilename, msg, end=delete_line_index)

## PIP INSTALL =============================================================================================================

def pip_install(modules):
    required  = modules
    installed = {pkg.key for pkg in pkg_resources.working_set}
    missing   = required - installed

    if missing:
        subprocess.check_call([sys.executable, "-m", "pip", "install", "--user", *missing], 
                              stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        
        print_table_row(["INSTALL", "Required Modules installed"])
    else:
        print_table_row(["CHECK", "Modules already installed"])

## FILE ====================================================================================================================

def get_value_of_variable_from_file(filename, file_index, relative_index, string):
    try:
        content = [i.strip().split() for i in open(filename).readlines()]
        index = [idx for idx, s in enumerate(content) if string in s][file_index]
        value = content[index][relative_index]
        return value
    except IndexError:
        print_table_row(["ERROR", "String not found in file"])
        quit()

def find_string_in_file(filename, string):
    
    with open(filename) as f:
        if string in f.read():
            return True
        else:
            return False
        
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

# Reset gkwdata.h5 with the use of dump files DM1, DM2
# AUTHOR: Florian Rath
# IMPORT: gkw_reset_checkpiont.py (https://bitbucket.org/gkw/gkw/src/develop/python/gkw_reset_checkpoint.py)  
def reset_simulation(SIM_DIR, NTIME=None, use_ntime=False):
    
    # Function that delets all data in interval [nt_reset:nt_broke].
    # ncol considers, if data series is ordered by a multiple interger of
    # ntime.
    def reset_time_trace(indata, dim, ncol, nt_reset):
    
        # Shift dimension dim to 0.
        data_shifted = np.moveaxis(indata, dim, 0)
  
        # Reset data.
        data_shifted_trimmed = data_shifted[0:int(nt_reset*ncol),]
    
        # Shift dimension back.
        out = np.moveaxis(data_shifted_trimmed, 0, dim)
  
        # free memory
        del data_shifted
        del data_shifted_trimmed
  
        return out
    
    # Checks if file has binary format.
    def is_binary(filename):
        try:
            with open(filename, 'tr') as check_file:
                check_file.read()
                return False
        except:
            return True
    
    
    
    # Check for specific files that are no ordinary data files.
    def is_file_exception(filename):
    
        # substrings that have to be checked
        check_list = ['geom.dat', 'DM1.dat', 'DM2.dat', 'FDS.dat', '.o', 'FDS', 
                      'input.dat', 'perform_first.dat', 'perform.dat', 'output.dat', 
                      'gkwdata.meta', 'gkw_hdf5_errors.txt', 'kx_connect.dat', 
                      'jobscript', 'Poincare1.mat', 'perfloop_first.dat', 'par.dat', 
                      'input_init.dat', 'sgrid', 'gkw', '.out', 'status.txt']
        for key in check_list:
            if key in filename:
                return True
    
        return False
    
    
    # Check if given file is a PBS or SLURM jobscript.
    def is_jobscript(filename):
        with open(filename,'r') as file:                                            
            for line in file:
                if '#PBS -l' in line:
                    return True
                if '#SBATCH' in line:
                    return True
            
        return False
    
    def get_timestep(FILE):
        if(os.path.isfile(SIM_DIR+'/'+FILE+'.dat')):
            EXISTS = True
            with open(SIM_DIR+'/'+FILE+'.dat','r') as file:                                                                            
                for line in file:
                    if 'NT_REMAIN' in line:
                        expr = line.replace(' ','')
                        expr = expr.replace(',','')
                        expr = expr.replace('\n','')
                        REMAIN = int(expr.split('=')[1])
                    if 'NT_COMPLETE' in line:
                        expr = line.replace(' ','')
                        expr = expr.replace(',','')
                        expr = expr.replace('\n','')
                        COMPLETE = int(expr.split('=')[1])
        
        else:
            REMAIN, COMPLETE, EXISTS = None, None, False
            
        return REMAIN, COMPLETE, EXISTS
    
    HDF5_FILENAME = "gkwdata.h5"
    RESTARTFILE, DUMPFILE1, DUMPFILE2, = "FDS", "DM1", "DM2"
  
    # Change to simulation directory.
    if(not os.path.isdir(SIM_DIR)):
        return
    else:
        os.chdir(SIM_DIR)
  
    # Check if hdf5-file exists.
    if(not os.path.isfile(HDF5_FILENAME)):
        return
  
    # First, read hdf5 file and determine the number of big time steps NTIME, 
    # requested in the input.dat file.
    f = h5py.File(HDF5_FILENAME, "r+")
    
  
    # Get requested big time steps from the /control group in the hdf5-file.
    if(NTIME==None):
        NTIME = int(f['input/control/ntime'][:])
  
    # Get number of big time steps after which simulation broke.
    # If time.dat exists read this file to obtain number of time steps after
    # which simulation broke.
    if(os.path.isfile('time.dat')):
      tim = pd.read_csv(SIM_DIR+'time.dat', header=None, sep='\s+').values
      NT_BROKE = tim.shape[0]
    # Else, get time from hdf5-file.
    else:
      NT_BROKE = f['diagnostic/diagnos_growth_freq/time'].shape[1]
  
  
    # Set NT_BROKE for output files holding temporal derivates and therefore 
    # one timestep less.
    NT_BROKE_DERIV = NT_BROKE-1
  
    # Close the hdf5-file again.
    f.close()
    
    # Get the number of remaining big time steps NT_REMAIN from checkpoint 
    # files FDS.dat. This is used lateron to determine the most recent 
    # checkpoint file.
    NT_REMAIN, NT_COMPLETE, FDS_EXISTS = get_timestep(RESTARTFILE)
    
    # Get the number of remaining big time steps NT_REMAIN[1/2] from checkpoint 
    # files DM[1/2]. This is used lateron to determine the most recent 
    # checkpoint file.
    NT_REMAIN1, NT_COMPLETE1, DM1_EXISTS = get_timestep(DUMPFILE1)
    NT_REMAIN2, NT_COMPLETE2, DM2_EXISTS = get_timestep(DUMPFILE1)
          
          
    # Check if FDS is the most recent checkpoint file. In this case
    # resetting the simulation makes no sense.
    if(FDS_EXISTS):
        DM1_OLD, DM2_OLD = False, False
        if(DM1_EXISTS):
            if(NT_COMPLETE1 < NT_COMPLETE):
                DM1_OLD = True
        if(DM2_EXISTS):
            if(NT_COMPLETE2 < NT_COMPLETE):
                DM2_OLD = True
        if(DM1_OLD and DM2_OLD):
            return
          
          
    # Now determine which checkpoint file is the recent one.
    if(DM1_EXISTS and DM2_EXISTS):
        if(NT_REMAIN1 > NT_REMAIN2):
            NT_REMAIN, NT_COMPLETE, DUMPFILE = NT_REMAIN2, NT_COMPLETE2, DUMPFILE2
        else:
            NT_REMAIN, NT_COMPLETE, DUMPFILE = NT_REMAIN1, NT_COMPLETE1, DUMPFILE1
            
    elif(DM1_EXISTS and not DM2_EXISTS):
        NT_REMAIN, NT_COMPLETE, DUMPFILE = NT_REMAIN1, NT_COMPLETE1, DUMPFILE1
        
    elif(DM2_EXISTS and not DM1_EXISTS):
        NT_REMAIN, NT_COMPLETE, DUMPFILE = NT_REMAIN2, NT_COMPLETE2, DUMPFILE2 
        
    else:
        return
    
    
    if(not use_ntime):
        # Find the total number ob big time steps the simulation time trace should have, 
        # when considering the big time steps already completed as well as the big time 
        # steps that remain. Can be different from NTIME, since the simulation could 
        # have been restarted several times such that NTIME > NT_COMPLETE.
        NTOT = NT_COMPLETE + NT_REMAIN
        N_REQUEST = NTOT
    
        # Determine the time steps to which the time trace files have to be reset.
        # Use NTOT here, since NTIME could have been changed at some point, or NT_COMPLETE
        # could be larger than NTIME.
        NT_RESET = NTOT-NT_REMAIN
    else:
        # Determine the time steps to which the time trace files have to be reset.
        NT_RESET = NTIME-NT_REMAIN
        N_REQUEST = NTIME
    
    # Same for files holding time derivatives.
    NT_RESET_DERIV = NT_RESET-1
  
    # ----------------------------------------------------------------------
    # Cycle over all nodes of hdf5-file and reset time trace datasets.
  
    # Check if hdf5-file exists.
    if(os.path.isfile(HDF5_FILENAME)):
    
        # Find all possible keys items, i.e. both groups and datasets
        f = h5py.File(HDF5_FILENAME, "a")
        h5_keys = []
        f.visit(h5_keys.append)
  
        # Cycle over all keys items and check, if any dimension has size NT_BROKE, 
        # i.e., it is a time trace file.
        for item in enumerate(h5_keys):
      
            data = f.get(item)
          
            # Consider datasets only.
            if(isinstance(data, h5py.Dataset)):
        
                #Check if any dimension is an integer multiple of NT_BROKE, by checking
                # the residual of the division.
                res = [None]*len(data.shape)
                for i in range(len(data.shape)):
                    res[i] = data.shape[i]/NT_BROKE-np.floor(data.shape[i]/NT_BROKE)
          
                if 0.0 in res:
        
                    # Check which dimension is integer multiple of NT_BROKE
                    # and save dimension as well as integer.
                    new_shape = data.shape
                    for i in range(len(data.shape)):
                        ncol = data.shape[i]/NT_BROKE
                        res = ncol - np.floor(ncol)
                        if(res == 0):
                            dim = i
                            ncol = int(ncol)
                            # adjust new shape to ncol*NT_RESET
                            y = list(new_shape)
                            y[dim] = int(ncol*NT_RESET)
                            new_shape = tuple(y)
                            break
          
                    # Reset dataset (.resize discards data with indices larger than 
                    # ncol*NT_RESET along dimension dim).
                    dset = f[item]
                    dset.resize(int(ncol*NT_RESET),dim)
  
  
        # After having repaired all datasets, close the hdf5-file again.
        f.close()
    
    
    # ----------------------------------------------------------------------
    # Cycle over all csv-files and reset time trace.
    for filename in os.listdir(SIM_DIR):
      
        # First perform some checks on files; cycle if file is binary, an exception
        # or a jobscript.
        if(is_binary(filename)):
            continue
        if(is_file_exception(filename)):
            continue
        if(is_jobscript(filename)):
            continue
        
        # no file
        if(not os.path.isfile(filename)):
            continue
    
        # Load csv file.
        data =  pd.read_csv(SIM_DIR+'/'+filename, header=None, sep='\s+').values
    
        #Check if any dimension is an integer multiple of NT_BROKE, by checking
        #the residual of the division.
        res = [None]*len(data.shape)
    
        # residul for output that holds time derivatives and therefore nt-1 datapoints
        res_deriv = [None]*len(data.shape)
        for i in range(len(data.shape)):
            res[i] = data.shape[i]/NT_BROKE-np.floor(data.shape[i]/NT_BROKE)
            res[i] = data.shape[i]/(NT_BROKE_DERIV)-np.floor(data.shape[i]/(NT_BROKE_DERIV))
      
        # ordinary files
        if 0.0 in res:
    
            #Check which dimension is integer multiple of NT_BROKE
            #and save dimension as well as integer.
            for i in range(len(data.shape)):
                ncol = data.shape[i]/NT_BROKE
                res = ncol - np.floor(ncol)
                if(res == 0):
                    dim = i
                    ncol = int(ncol)
                    break
              
            # Load original dataset.
            original_data = data
            # print('\t Original shape: \t'+str(original_data.shape))
      
            # Reset time trace.
            reset_data = reset_time_trace(original_data,dim,ncol,NT_RESET)
            # print('\t Reset shape: \t' +str(reset_data.shape))
    
            # Save resetted data.
            pd.DataFrame(reset_data).to_csv(filename, sep='\t', header=None, index=None)
 
        # files holding time derivatives
        if 0.0 in res_deriv:
    
            # Check which dimension is integer multiple of NT_BROKE
            # and save dimension as well as integer.
            for i in range(len(data.shape)):
                ncol = data.shape[i]/NT_BROKE_DERIV
                res = ncol - np.floor(ncol)
                if(res == 0):
                    dim = i
                    ncol = int(ncol)
                    break
              
            # Load original dataset.
            original_data = data
      
            # Reset time trace.
            reset_data = reset_time_trace(original_data,dim,ncol,NT_RESET_DERIV)
    
            # Save resetted data.
            pd.DataFrame(reset_data).to_csv(filename, sep='\t', header=None, index=None)
      
  
    # ----------------------------------------------------------------------
    # Finally, copy most recent dump file to FDS[/.dat].
  
    # Copy the most recent dump file to FDS[/.dat].
    copyfile(SIM_DIR+'/'+DUMPFILE, SIM_DIR+'/'+'FDS')
    copyfile(SIM_DIR+'/'+DUMPFILE+'.dat', SIM_DIR+'/'+'FDS.dat')
  
    # First line of the so produced FDS.dat has to be modified. 
    old_text = '!Dump filename: '+DUMPFILE
    new_text = '!Dump filename: '+'FDS'
  
    # Replace first line in FDS.dat to set the correct file name.
    with fileinput.input(SIM_DIR+'/'+'FDS.dat',inplace=True) as f:
        for line in f:
            line.replace(old_text, new_text)
            
def check_and_delete_file(filenames):
    
    for filename in filenames:
        if(os.path.isfile(filename)):
            os.remove(filename)


def check_checkpoint_files():
    DM1, DM2 = False, False
    
    if(os.path.isfile(check1Filename) and os.path.isfile(check1Bin)):
        DM1 = True
    if(os.path.isfile(check2Filename) and os.path.isfile(check2Bin)):
        DM2 = True
    
    return DM1, DM2

def get_timestep_from_restartfile(filename, flag):
    return int(get_value_of_variable_from_file("./" + filename, 0, 2, flag).replace(",", ""))

def get_ntimestepCurrent(filenames):
    
    ntimestep = 0
    
    for f in filenames:
        if(os.path.isfile(f)):
            ntimestepFile = get_timestep_from_restartfile(f, restartFlag)
            
            if ntimestepFile > ntimestep:
                ntimestep = ntimestepFile
        
    return ntimestep

def get_time_from_statusfile(filename, line_index):
    with open(filename, "r") as file:
        line = file.readlines()[line_index]
        
        content = line.split(" ")
        time = content[-2]
        time_sec = get_time_in_seconds(time)
        
        return time_sec

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

    seconds = [7*24*60*60,24*60*60, 60*60, 60, 1]
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

def get_error_type(filename):
    
    slurm_errors = {"executable":["error on file ./gkw.x (No such file or directory)", "No executable found"],
                    "walltime":["process killed (SIGTERM)", "Exceeded wall time"],
                    "timeout":["DUE TO TIME LIMIT", "Exceeded time limit"],
                    "config":["couldn't open config directory", "Config not loading"],
                    "hdf5":["HDF5-DIAG", "Writing h5 file failed"]}
    
    for key in slurm_errors:
        if find_string_in_file(filename, slurm_errors[key][0]):
            return slurm_errors[key][1]
    
    return "Unknown error occurred"

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

## KILL JOB ================================================================================================================

PID = subprocess.getoutput(commandMonitorKill).split("\n")[0]
if kill:
    subprocess.run(["kill", PID])
    quit()

# START/RESTART JOB ========================================================================================================

startTime = time.time()
outputType = set_output_type()
jobStatusRunning, jobStatusPending = get_job_status()

# Set current Time for progress bar
if outputType == "running":
    jobStatusRunningNameIndex = [idx for idx, s in enumerate(jobStatusRunning) if jobName in s][0]
    currentTime = jobStatusRunning[jobStatusRunningNameIndex + 3]
elif outputType == "pending":
    jobStatusPendingNameIndex = [idx for idx, s in enumerate(jobStatusPending) if jobName in s][0]
    currentTime = jobStatusPending[jobStatusPendingNameIndex + 3]
else:
    currentTime = "00:00:00"

# Set pastTime and create status file if necessary. When status file exist append next lines
WRITEHEADER = True
if not os.path.isfile(statusFilename):
    pastTime = 0
    write_file(statusFilename,"")
else:
    try:
        pastTime = get_time_from_statusfile(statusFilename, -11)
        WRITEHEADER = False
    except IndexError:
        pastTime = 0  

print_table_row(["OUTPUT", "INFO"], output_type="header", WRITEFILE=WRITEHEADER)

## BEGIN ===================================================================================================================

# Check if timesteps criterion is satisfied, send mail and end monitoring 
# else continue monitoring and set output type accordingly
nTimestepsCurrent = get_ntimestepCurrent([restartFilename, check1Filename, check2Filename])
    
if nTimestepsCurrent >= nTimestepsRequired:
    print_table_row(["CONTROL", "Current Timesteps " + str(nTimestepsCurrent)], output_type="middle")
    print_table_row(["SUCCESS", "Stop monitoring " + jobName])
    
    if EMAIL:
        send_mail(emailAddress, "Ended Job " + jobName)
        
    quit()

# Continue monitoring and send mail
else:
    
    if outputType != "no Output":
        print_table_row(["CONTINUE", "Continue monitoring " + jobName], output_type="middle")
        print_table_row(["CONTROL", "Current Timesteps " + str(nTimestepsCurrent)])
        
        if EMAIL:
            send_mail(emailAddress, "Continued Job " + jobName)
            
    elif os.path.isfile(restartFilename):
        print_table_row(["CONTINUE", "Continue monitoring " + jobName], output_type="middle")
        
        if EMAIL:
            send_mail(emailAddress, "Continued Job " + jobName)
    
    else:
        print_table_row(["STARTING", "Start monitoring " + jobName], output_type="middle")
        
        if BACKUP:
            print_table_row(["BACKUP", backupPath])
            subprocess.run(commandBackup)
            
        if EMAIL:
            send_mail(emailAddress, "Started Job " + jobName)

if RESET:
    pip_install({"h5py", "pandas", "numpy"})
    
    import h5py, fileinput
    import pandas as pd
    import numpy as np
    from shutil import copyfile
    
    print_table_row(["IMPORT", "Load numpy, pandas, h5py"])

#else:
#    pip_install({"h5py"})
#    import h5py
#    print_table_row(["IMPORT", "Load module h5py"])

## MONITOR ROUTINE =========================================================================================================

while True:
    
    # Check job status to monitor current state 
    jobStatusRunning, jobStatusPending = get_job_status()

    # Job running
    if jobName in jobStatusRunning:
        
        jobStatusRunningNameIndex = [idx for idx, s in enumerate(jobStatusRunning) if jobName in s][0]
        jobID = jobStatusRunning[jobStatusRunningNameIndex - 2]
        
        currentTime = jobStatusRunning[jobStatusRunningNameIndex + 3]
        nTimestepsCurrent = get_ntimestepCurrent([restartFilename, check1Filename, check2Filename])
        
        if outputType == "running":
            print_table_row(["RUNNING", "Job is executed"])
            outputType = "no Output"
        else:
            print_table_row(["RUNNING", "Job is executed"], output_type="update")
            
        time.sleep(sleepTime)
        
    # Job pending
    elif jobName in jobStatusPending:

        jobStatusPendingNameIndex = [idx for idx, s in enumerate(jobStatusPending) if jobName in s][0]
        jobID = jobStatusPending[jobStatusPendingNameIndex - 2]
        
        currentTime = jobStatusPending[jobStatusPendingNameIndex + 3]
        nTimestepsCurrent = get_ntimestepCurrent([restartFilename, check1Filename, check2Filename])
        
        if outputType == "pending":
            print_table_row(["WAITING", "Job is pending"])
            outputType = "running"
        else:
            print_table_row(["WAITING", "Job is pending"], output_type="update")
            
        time.sleep(sleepTime)

    # Job start/restart
    else:
        
        # Check errors and making Backup
        while True:
            try:
                outputContent = open(outputFilename(jobID)).readlines()[-5].replace("\n","")
                
                # Create Backup if run is successful
                # If scan of output.dat is needed: Scans for string "Run Successful in output.dat and returns bool value"
                #runSuccess = find_string_in_file("output.dat", outputCriteria[1])
                #if outputCriteria[0] in outputContent and runSuccess: 
                
                if outputCriteria[0] in outputContent: 
                    
                    # Check if h5 file is closed before start/restart simulation
                    # than check if FDS/FDS.dat is updated
                    try:
                        f = open(dataFilename)
                        #f = h5py.File(dataFilename)
                        f.close()
                        
                        # Check if FDS/FDS.dat is updated after run and has equially time stamp as gkwdata.h5                    
                        timestamp_data    = int(os.path.getmtime(dataFilename))
                        timestamp_restart = int(os.path.getmtime(restartFilename))

                        wallTime_sec = get_time_in_seconds(wallTime)
                        timestamp_remain = timestamp_data - timestamp_restart

                        # FDS/FDS.dat does not get written at the same time as gkwdata.h5
                        # For that a time interval have to be considered 
                        # To be certain the half wall time is set aus time interval
                        if timestamp_remain > wallTime_sec/2:

                            print_table_row(["ERROR", "FDS/FDS.dat not updated"])

                            # Reset simulation and save as backup
                            if RESET:
                            
                                DM1, DM2 = check_checkpoint_files()

                                if (DM1 or DM2):
                                    print_table_row(["RESET", "Reset to last checkpoint."])
                                    reset_simulation(folder)

                                    # Update backup
                                    if BACKUP:
                                        print_table_row(["BACKUP", backupPath])
                                        subprocess.run(commandBackup)

                            # Restore backup to rerun simulation    
                            elif BACKUP:
                                print_table_row(["RESTORE", backupPath])
                                subprocess.run(commandRestore)

                        elif BACKUP:
                            print_table_row(["BACKUP", backupPath])
                            subprocess.run(commandBackup)
                    
                        break
                    except OSError:
                        time.sleep(sleepTime)
                        
                else:
                    print_table_row(["ERROR", get_error_type(outputFilename(jobID))])
                    
                    # Reset simulation and save as backup
                    if RESET:
                        
                        DM1, DM2 = check_checkpoint_files()
                    
                        if (DM1 or DM2):
                            print_table_row(["RESET", "Reset to last checkpoint."])
                            reset_simulation(folder)
                        
                            # Update backup
                            if BACKUP:
                                print_table_row(["BACKUP", backupPath])
                                subprocess.run(commandBackup)
                    
                    # Restore backup to rerun simulation    
                    elif BACKUP:
                        print_table_row(["RESTORE", backupPath])
                        subprocess.run(commandRestore)
                    
                    break
            
            # Wait sleepTime and check output file again
            except (IndexError, FileNotFoundError):
                time.sleep(sleepTime)
                
            # If jobID undefined break cycle
            except NameError:
                break

        # Check if timesteps criterion is satisfied, send mail and end monitoring
        nTimestepsCurrent = get_ntimestepCurrent([restartFilename, check1Filename, check2Filename])
        print_table_row(["CONTROL", "Current Timesteps " + str(nTimestepsCurrent)])
        
        if nTimestepsCurrent >= nTimestepsRequired:
            print_table_row(["SUCCESS", "Stop monitoring " + jobName])
            
            if EMAIL:
                send_mail(emailAddress, "Ended Job " + jobName)
            break
        
        # Delete checkpoint files
        check_and_delete_file([check1Bin, check1Filename, check2Bin, check2Filename])
        
        # Create jobscript
        if jobscriptFilename == "jobscript-create":
            jobscriptFilename = "jobscript"
            write_file(jobscriptFilename, jobscriptContent)
        
        # Start Job and send restart mail (if activated)
        startOutput = subprocess.check_output([commandJobStarting, jobscriptFilename]).decode("utf-8").replace("\n", "")
        jobID = startOutput.split(startOutputFlag)[1]
        
        runCounter += 1
        
        print_table_row(["STARTING", startOutput], output_type="middle")
        
        try:
            if RESTARTMAIL and EMAIL:
                send_mail(emailAddress, "Restart Job " + jobName)
        except NameError:
            continue
        
        time.sleep(30)
        
        # Set output type to running or pending
        outputType = set_output_type()
        
## RESTART =================================================================================================================