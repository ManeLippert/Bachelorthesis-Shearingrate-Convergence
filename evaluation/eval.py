import gkw
import derivative
import h5py
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Flux rlt = 6.0
filename = 'evaluation/data/S6_rtl6.0/gkwdata.h5'
f = h5py.File(filename,"r+")
data = gkw.get_data_keys(f)
eflux_data = gkw.get_eflux_from_hdf5_file(f)

# Elektrostatic potencial
phi = f[data[0][0][4]]

# Mean over n_y
# !TODO

# Finite Differnece for shearing rate omega_ExB
# !TODO

start = 0
stop = 200
# define grid
x = np.arange(start, stop) 
# compute function
y = phi[0][0][start:stop]

# compute vector of forward differences
forward_diff = derivative.finite_first_order(y, 0.1, 'forward')
backward_diff = derivative.finite_first_order(y, 0.1, 'backward')
central_diff = derivative.finite_first_order(y, 0.1, 'central')
# compute corresponding grid
x_diff = x[:-1:] 

#Data:

# Plotting
plt.plot(x, y)
plt.plot(x_diff, forward_diff, '--', label = 'Finite difference approximation')
plt.plot(x_diff, backward_diff, '--', label = 'Finite difference approximation')
plt.plot(x[1:-1], central_diff, '--', label = 'Finite difference approximation')
plt.legend()
plt.show()


# !Important! close h5 file after usage 
f.close()