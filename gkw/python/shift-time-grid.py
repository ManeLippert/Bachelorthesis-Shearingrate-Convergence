#!/usr/bin/env python
import h5py
import bisect
import sys

h5filename = "gkwdata.h5"


def print_usage_and_exit():
    print("""
  Usage:
        %s [delta-time]

  This script acts on %s in the current directory.
  Without argument, it shifts the time grid so that
  it starts at
    
    0 + (time[1] - time[0])

  If the real number delta-time is given, then the shift
  is by delta-time instead.
     
  See also: defer-gkwdatah5-attrs-from-2nd-file.py
"""
    % ( sys.argv[0], h5filename))
    exit(1)
        
if __name__ == "__main__":
    if(len(sys.argv) >= 3):
        print_usage_and_exit()
    else:
        f = h5py.File("gkwdata.h5", "r+")

        try:
            time = f["/grid/time"][0,:]
        except:
            print("/grid/time does not exist in " + h5filename)
            exit(1)

        if(len(sys.argv) >= 2):
            delta_time = float(sys.argv[1])
        else:
            delta_time = time[0] - (time[1] - time[0])

        print(time[0:5])
        print("...")
        print("Shift time by " + str(delta_time))
        time = time - delta_time
        print(time[0:5])
        print("...")
        f["/grid/time"][0,:] = time

        
