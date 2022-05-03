#!/usr/bin/env python3

def make_run_folder(dirname, orig_input_file=None, input_file_namelist=None, input_prof_data=None, remove_existing=True):
    import shutil
    import os
    
    try:
        shutil.rmtree(dirname)
    except:
        pass
    print("Writing the file %s/input.prof..." % dirname)
    os.mkdir(dirname)
    if(input_prof_data is not None):
        write_input_prof(os.path.join(dirname, "input.prof"), input_prof_data)
    if(input_file_namelist is not None):
        with open(os.path.join(dirname,"input.dat"),'w') as f:
            f.write(input_file_namelist.dump())
    elif(orig_input_file is not None):
        shutil.copy(orig_input_file, dirname)
    os.symlink(os.path.join(os.environ["HOME"],"bayreuth-scripts","btcluster","jobscript_bt"),
               os.path.join(dirname,"jobscript"))
    os.symlink(os.path.join(os.environ["HOME"],"bayreuth-scripts","marconi","jobscript"),
               os.path.join(dirname,"jobscript_marconi"))

def read_input_prof_data(input_prof_filename, blockname=None):
    import re
    import numpy as np
    ret = []
    with open(input_prof_filename, "r") as f:
        content = f.read()
        for block in re.split(r'^#', content, flags=re.MULTILINE):
            if(len(block) == 0):
                continue
            header, data = block.split(sep='\n', maxsplit=1)
            # we have to use strip, because notation with and without
            # blanks may occur mixed.
            data_parsed = np.array([[float(v) for v in line.split()] for line in data.split(sep='\n') if len(line)>0])
            if(blockname is not None and header.strip().startswith(blockname.strip())):
                #parse the whole thing
                return data_parsed
            elif(blockname is None):
                ret.append((header, data_parsed))
                
    if(blockname is not None):
        raise KeyError("A block '%s' could not be found in %s" % (blockname, input_prof_filename))
    elif(blockname is None):
        return ret

def get_input_prof_block_index(complete_input_prof_data, block_name):
    for i in range(len(complete_input_prof_data)):
        if(complete_input_prof_data[i][0].strip().startswith(block_name)):
            return i
    raise KeyError("There is no column %s" % column_name)

input_prof_background_prof_columns = ["xgr", "dens", "temp", "rln","rlt"]
input_prof_miller_columns = ['xgr',
                      'q',
                      'shat',
                      'kappax',
                      'deltax',
                      'squarex',
                      'skappax',
                      'sdeltax',
                      'ssquarex',
                      'Zmilx',
                      'dRmilx',
                      'dZmilx',
                      'gradpx',
                  ]

def write_input_prof(filename, complete_data):
    with open(filename,"w") as f:
        for header, data in complete_data:
            f.write("#"+header+"\n")
            for row in data:
                f.write(" ".join(str(elem) for elem in row)+"\n")
