#!/usr/bin/env python

# This little script is useful if one has accidentally appended data
# to gkwdata.h5 by continuing a run.

# This script first looks for a jump in the time grid, and
# then it removes the data before that jump from the timetraces,
# making use of another python script.


import h5py
import bisect
import sys
import subprocess
import os

h5filename = "gkwdata.h5"


if __name__ == "__main__":
    f = h5py.File("gkwdata.h5", "r+")

    try:
        time_dset = "/grid/time"
        time = f[time_dset][0,:]
    except:
        # it is a (1,n) dimensional array
        time_dset = "/diagnostic/diagnos_growth_freq/time"
        time = f[time_dset][0,:]

    try:
        gkw_home = os.environ['GKW_HOME']
    except:
        print("GKW_HOME environment variable is not set.")
        exit(1)

    # find index where the time jumps to a lower value
    for i in range(1,len(time)):
        if(time[i]< time[i-1]):
            # ok, nailed it.
            i-=1
            print("timegrid:" + str(time[i-5:i+5]))
            print("position:" + str(i))
            # now copy the data to a backup file
            print(subprocess.call(["cp", "-v","gkwdata.h5","backup_gkwdata.h5"]))
            print("Data backup written to : backup_gkwdata.h5")
            # and remove the chunk of data in front of the time grid jump
            print("call:    " + " ".join([gkw_home+"/python/clip_gkwdatah5.py", "--remove-before","-n",str(i)]))
            print(subprocess.call([gkw_home+"/python/clip_gkwdatah5.py", "--remove-before","-n",str(i)]))
            exit(0)

    print("No jump in " + time_dset + " was found.")
    
    
    
