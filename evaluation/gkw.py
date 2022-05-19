import h5py
import numpy as np

def get_data_keys(f):

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
    
    return flux