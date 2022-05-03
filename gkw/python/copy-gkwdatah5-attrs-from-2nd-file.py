#!/usr/bin/env python
import h5py
import bisect
import sys

h5filename = "gkwdata.h5"

def copy_attrs_if_obj_exists(name, object):
    """It is checked if the respective group
  (by name) exists also in <file-to-be-edited>. If this is the case,
  then the attributes are copied."""
    if name in f:
        print(name+" exists also in f!")
        if object.attrs:
            #print object.attrs.items()
            for key in object.attrs.keys():
                f[name].attrs[key] = object.attrs[key]

if __name__ == "__main__":
    if(len(sys.argv) != 3):
        print("""
  Usage:
        %s <HDF5-file-to-be-edited> <reference-HDF5-file>

  This script iterates through the datasets and groups in
  <reference-file>. For each one it is checked if the respective object
  (by name) exists also in <file-to-be-edited>. If this is the case,
  then the attributes are copied.

  Typical use case: If GKW fails then %s does most likely not contain
  attributes. In that case one may take %s from a successful run of
  the same scan, and copy the attributes.

  Then a few attributes may be wrong in <file-to-be-edited>, such as
  the input settings, because they show the setting of the
  reference. If this is a problem they may be edited, e.g. with
  HDFView.
  
  Now, one may use the script clip_gkwdatah5_for_checkpoint.py .
  This script makes use of some time-related information in the
  attributes of datasets.

  Finally one can try to continue the run at the last checkpoint. GKW
  will rewrite the attributes with the input settings, as they are read
  from the namelist file.
        
"""
    % ( sys.argv[0], h5filename, h5filename))
        exit(1)
    else:
        filename_to_edit = sys.argv[1]
        filename_ref = sys.argv[2]
    
    f = h5py.File(filename_to_edit, "r+")
    fref = h5py.File(filename_ref, "r")

    fref.visititems(copy_attrs_if_obj_exists)

    
