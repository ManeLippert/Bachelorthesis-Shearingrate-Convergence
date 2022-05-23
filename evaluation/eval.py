import gkw
import derivative
import h5py
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

data = 'S6_rtl6.0'
filename = 'data/'+data+'/gkwdata.h5'
f = h5py.File(filename,"r+")
eflux_data, time = gkw.get_eflux_from_hdf5_file(f)

gkw.find_keys(f, 'phi')

# Elektrostatic potencial
phi = f[gkw.find_key(f, 'phi')][()]
nx = phi.shape[0]

# Mean over y to get a approximation for the zonal potenzial
zonal_pot = np.mean(phi,1)

# Finite Differnece for shearing rate omega_ExB

# Stepsize
rad_boxsize = f[gkw.find_key(f, 'lxn')][()][0]
dx = rad_boxsize/nx

# Plot all
#start, end = 0, 999 

gkw.plot_shearing_rate_wexb_all(zonal_pot, dx, 0, 999)
plt.show()
    
#plt.savefig('pictures/'+data+'/'+data+'wexb_all.pdf', bbox_inches='tight')

# !Important! close h5 file after usage 
f.close()