import gkw
import derivative
import h5py
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Flux rlt = 6.0
filename = 'evaluation/data/S6_rtl6.0/gkwdata.h5'
f = h5py.File(filename,"r+")
eflux_data, time = gkw.get_eflux_from_hdf5_file(f)

# Elektrostatic potencial
phi = f[gkw.find_key(f, 'phi')][()]
nx = phi.shape[0]

# Mean over y to get a approximation for the zonal potenzial
zonal_pot = np.mean(phi,1)

# Finite Differnece for shearing rate omega_ExB

# Stepsize
rad_boxsize = f[gkw.find_key(f, 'lxn')][()][0]
dx = rad_boxsize/nx

gkw.find_keys(f, 'phi')

plt.plot(eflux_data)
#plt.show()
#plt.clear()

i = 4999
while i <= zonal_pot.shape[1]:
    # Zonal potenzial for a fixed time
    y = zonal_pot[:,i]

    ddphi = derivative.finite_second_order(y, dx, 'period')
    wexb = 0.5 * ddphi

    plt.plot(np.arange(0, len(ddphi)), ddphi)
    
    i += 1000
    
#plt.show()


# !Important! close h5 file after usage 
f.close()