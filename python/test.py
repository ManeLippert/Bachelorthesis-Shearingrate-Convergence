restartFlag = "NT_COMPLETE"


import os

def get_value_of_variable_from_file(filename, file_index, relative_index, string):
    try:
        content = [i.strip().split() for i in open(filename).readlines()]
        index = [idx for idx, s in enumerate(content) if string in s][file_index]
        value = content[index][relative_index]
        return value
    except IndexError:
        quit()

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
    
print(get_ntimestepCurrent(["FDS.dat", "DM1.dat", "DM2.dat"]))

check1Filename, check2Filename = "DM1.dat", "DM2.dat"
check1Bin, check2Bin = check1Filename.replace(".dat", ""), check2Filename.replace(".dat", "")

print(check1Bin)