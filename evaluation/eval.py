import gkwdata
import h5py
import numpy as np
import matplotlib.pyplot as plt

filename = 'evaluation/data/S6_rtl6.0/gkwdata.h5'
f = h5py.File(filename,"r+")

data = gkwdata.gkw_data(f)
print(data)

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


# !Important! close h5 file after usage 
f.close()