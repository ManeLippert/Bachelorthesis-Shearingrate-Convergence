import h5py
import numpy as np

import derivative

import matplotlib.pyplot as plt
from matplotlib import rc


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
            elif type(f[keys0 + '/' + keys1]) == h5py._hl.group.Group:
                for keys2 in list(f[keys0 + '/' + keys1].keys()):
                    data[i][j].append(keys0 + '/' + keys1 + '/' + keys2)
            j += 1
        i += 1
        
    return data

# From Florian Rath
def get_eflux_from_hdf5_file(hdf5_file):
    
    # find out number of time steps
    tim = hdf5_file['diagnostic/diagnos_growth_freq/time'][()]
    nt = tim.shape[1]
    
    # name of the dataset
    node_name = 'diagnostic/diagnos_fluxes/eflux_species01'
    
    # load data into array
    data = np.transpose(hdf5_file[node_name][()])
    
    # reshape GKW flux ordering
    flux = np.reshape(data,(nt,2))[:,0]
    
    return flux, tim

# Returns key if the end of the level is equal the vairable search
def find_key(f, search):

    # i, j, k are indexes for the function get_data_keys()
    
    # First level
    #i = 0
    for keys0 in list(f.keys()):
        # Second level
        #j = 0
        for keys1 in list(f[keys0].keys()):
            # Third level
            #k = 0
            if type(f[keys0 + '/' + keys1]) == h5py._hl.dataset.Dataset:
                key = keys0 + '/' + keys1
                if search == keys1:
                    return key
            elif type(f[keys0 + '/' + keys1]) == h5py._hl.group.Group:
                for keys2 in list(f[keys0 + '/' + keys1].keys()):
                    key = keys0 + '/' + keys1 + '/' + keys2
                    if search == keys2:
                        return key
                    #k += 1
            #j += 1
        #i += 1


# Prints all possible keys with the variable search in it
def find_keys(f, search):
    
    key_list = []

    # i, j, k are indexes for the function get_data_keys()
    
    # First level
    #i = 0
    for keys0 in list(f.keys()):
        # Second level
        #j = 0
        for keys1 in list(f[keys0].keys()):
            # Third level
            #k = 0
            if type(f[keys0 + '/' + keys1]) == h5py._hl.dataset.Dataset:
                key = keys0 + '/' + keys1
                if search in key:
                    key_list.append(key)
            elif type(f[keys0 + '/' + keys1]) == h5py._hl.group.Group:
                for keys2 in list(f[keys0 + '/' + keys1].keys()):
                    key = keys0 + '/' + keys1 + '/' + keys2
                    if search in key:
                        key_list.append(key)
                    #k += 1
            #j += 1
        #i += 1
     
    for i in key_list:
        print(i)

        
def shearing_rate_wexb_interval(zonal_pot, dx, start, end):
    
    middle = int((end - start)/2)

    ddphi= derivative.finite_second_order(zonal_pot[:,start:end], dx, 'period')
    wexb = 0.5 * np.mean(ddphi,1)
    wexb_middle = 0.5*ddphi[:,middle]
    
    return wexb, wexb_middle
