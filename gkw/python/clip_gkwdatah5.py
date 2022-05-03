#!/usr/bin/env python3
import h5py
import bisect
import sys

h5filename = "gkwdata.h5"
remove_before = None
clip_nltransfer = False
f = None
dsets_to_delete = []
clip_time = None
time_attr_entry = None

#a function used to visit the objects in the hdf5 file.
def clip_timedep_dsets(name, object):
    """ discard some of the data from datasets which have the grid1=time attribute. """
    if("grid1" in object.attrs and object.attrs["grid1"] == time_attr_entry):
        clip_last_dim(object, clip_index)

def clip_last_dim(object, n):
    """ Clip the data in object, so that the last dimension has n elements."""
    old_shape = object.shape
    newsize = list(object.shape)
    if(remove_before):
        # -> remove the first n values, with respect to the last
        # -> dimension.
        # In order to do this, buffer the data we want to keep.
        tmp = object[..., (n+1):]
        # then resize the dataset
        newsize[-1] -= n+1
        object.resize(newsize)
        # and copy back all the data. probably there is a better syntax...
        object[...,0:newsize[-1]] = tmp
    else:
        # note that the first array dimension as seen from Fortran
        # corresponds here to the *last* dimension (as in C and Java)
        newsize[-1] = n
        object.resize(newsize)
    print(object.name+" is clipped from "+str(old_shape)+" to "+str(object.shape))

#a function used to visit the objects in the hdf5 file.
def rm_dsets_at_latest_timestamps(name, object):
    """delete datasets which have a "time" attribute and whose value is
larger than clip_time."""
    if("time" in object.attrs):
        if((not remove_before and object.attrs[time_attr_entry] <= clip_time) or
           (remove_before and object.attrs[time_attr_entry] > clip_time)):
            #then discard this dataset, because it was created after the
            #request timestamp.
            dsets_to_delete.append(name)

def print_usage_and_exit():
    print("""
  Usage:
        %s [--remove-before | --remove-after] [--nonlin-transfer-time] [ -t <time> | -n <timesteps> ]

  This script acts on %s in the current directory.
  It looks for datasets with the attribute
     grid1 = 'time'
  and clips them according to the value
  given as an argument.
     --remove-before  Remove everything before the specified time, including
                       the timestep itself
     --remove-after    Remove everything after the specified time, excluding
                       the timestep itself
     --nonlin-transfer-time
                       Work on the special time grid of the nonlin. transfer
                       diagnostic and clip only data of that diagnostic.
                       This leaves the regular time grid and other diagnostics
                       untouched.
        
     -t <time>         the timestamp given in code time units (in general,
                       this is a real number)
     -n <timestep>     the timestamp number (this is an integer), where the
                       starting index is 0

  NOTE: If GKW fails then %s does most likely not contain
  attributes. In that case this script here is useless.

  NOTE: The HDF5 filesize only decreases when the file is repacked afterwards.
  To repack, call the tool
      h5repack %s tmp-repacked.h5
      # and if this has worked correctly:
      mv tmp-repacked.h5 %s

  See also: defer-gkwdatah5-attrs-from-2nd-file.py
"""
    % ( sys.argv[0], h5filename, h5filename, h5filename, h5filename))
    exit(1)
        
if __name__ == "__main__":
    if(len(sys.argv) < 4 or len(sys.argv) > 5):
        print_usage_and_exit()
    else:
        before_option = '--remove-before'
        after_option = '--remove-after'
        i = 1
        before_or_after_arg = sys.argv[i]
        if(not before_or_after_arg in (before_option, after_option)):
            print_usage_and_exit()
        remove_before = (before_or_after_arg == before_option)
        i += 1
        
        nltransfer_option = '--nonlin-transfer-time'
        arg = sys.argv[i]
        
        if(arg == nltransfer_option):
            clip_nltransfer = True
            time_attr_entry = "nonlin_transfer_time"
            i += 1
        else:
            time_attr_entry = "time"
        
        timestamp_option = '-t'
        timestep_option = '-n'
        timetype_arg = sys.argv[i]
        if(not timetype_arg in (timestamp_option, timestep_option)):
            print_usage_and_exit()
        i += 1
        clip_time = timetype_arg == timestamp_option and float(sys.argv[i]) or int(sys.argv[i])
    
    f = h5py.File(h5filename, "r+")

    if(clip_nltransfer):
        try:
            dsetname = "/diagnostic/diagnos_nonlin_transfer/nonlin_transfer_time"
            time = f[dsetname][0,:]
        except:
            print("The dataset %s does not exist." % dsetname)
            exit(1)
    else:
        try:
            time = f["/grid/time"][0,:]
        except:
            # it is a (1,n) dimensional array
            time = f["/diagnostic/diagnos_growth_freq/time"][0,:]

    #find out how many chunks this is.
    if(timetype_arg == timestamp_option):
        #Find rightmost value less than clip_time
        clip_index = bisect.bisect_left(time, clip_time)
        print(clip_index)
    else:
        clip_index = clip_time
        clip_time = time[clip_index]

    # clip datasets in those groups:
    if(clip_nltransfer):
        f["diagnostic/diagnos_nonlin_transfer"].visititems(clip_timedep_dsets)
    else:
        f["diagnostic"].visititems(clip_timedep_dsets)
        f["grid"].visititems(clip_timedep_dsets)

    # finally clip the time grid itself
    if(clip_nltransfer):
        clip_last_dim(f["/diagnostic/diagnos_nonlin_transfer/nonlin_transfer_time"], clip_index)
    else:
        try:
            clip_last_dim(f["/grid/time"], clip_index)
        except:
            clip_last_dim(f["/diagnostic/diagnos_growth_freq/time"], clip_index)

            # and remove datasets that are created as a whole and have timestamps
            f.visititems(rm_dsets_at_latest_timestamps)
            for dset in dsets_to_delete:
                print("DELETED: "+dset)
                del f[dset]

    
