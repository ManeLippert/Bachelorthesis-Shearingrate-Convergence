#!/usr/bin/env python3

import h5py
import os

# Returns a list with three levels, with all keys of gkw data
def get_keys(f):

    # Maximum level until datasets => 3
    data = []

    # First level
    i = 0
    for keys0 in list(f.keys()):
        data.append([])
        # Second level
        j = 0
        for keys1 in list(f[keys0].keys()):
            data[i].append([])
            # Third level
            if type(f[keys0 + '/' + keys1]) == h5py._hl.dataset.Dataset:
                data[i][j] = keys0 + '/' + keys1
                print('"' + keys0 + '/' + keys1 + '"')
            elif type(f[keys0 + '/' + keys1]) == h5py._hl.group.Group:
                for keys2 in list(f[keys0 + '/' + keys1].keys()):
                    data[i][j].append(keys0 + '/' + keys1 + '/' + keys2)
                    print('"' + keys0 + '/' + keys1 + '/' + keys2 + '"')
            j += 1
        i += 1
        
    return data


filename = [f for f in os.listdir() if f.endswith('.h5')]
for file in filename:
    
    print('\n')
    print('<----------------------- ' + file + ' ----------------------->')
    f = h5py.File(file, 'r+')
    get_keys(f)
    f.close()  
    print('\n')