# Import modules
import sys, os, h5py

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.transforms import Bbox

import zonalflow, h5tools, derivative, plot

# Plot parameter
plot.parameters(True, 22, (24,8), 300)

homepath = "/Users/manu/Bachelorthesis-Shearingrate-Wavelength/"
filepath = homepath + "/data/S6_rlt6.4/boxsize3x3/Ns16/Nvpar48/Nmu9"

# File import and Create picture folder
data = filepath[filepath.find('data/')+len('data/'):filepath.rfind('/boxsize')]
path = filepath.split(data + '/')[1]

resolution = path.replace('/', '_')
rlt = data.split("rlt",1)[1]
filename = filepath + '/' + 'data.h5'
f = h5py.File(filename,"r+")

picpath = homepath+'pictures/'+data+'/'+path+'/'
## Create target Directory if don't exist

if not os.path.exists(picpath):
    os.makedirs(picpath)

# linear growth rate
lineargrowth = pd.read_csv(homepath + 'data/linear_growthrate_ITG.dat', index_col=0)

#lineargrowth_rlt = lineargrowth['gamma_N'][float(data.split('rlt')[1])]
lineargrowth_rlt = 0.21367
lineargrowth_rlt_color = 'grey'

# Heat flux
eflux_data, time = zonalflow.get_eflux_time(f)

plot.eflux_time(time, eflux_data, (24,8))
plt.savefig(picpath+data+'_'+resolution+'_eflux.pdf', bbox_inches='tight')

# Finite Difference for shearing rate omega_ExB
wexb, rad_coord, rad_boxsize, ddphi, dx, zonal_pot = zonalflow.get_shearingrate_radialcoordinate_radialboxsize_ddphi_dx_zonalpot(f)

print('rad_boxsize :', rad_boxsize, '; stepsize :',dx)

# Fourier plot of time domain

wexb_max = zonalflow.get_max_shearingrate(f, wexb, time, 5)

plot.max_shearingrate_time(time, wexb_max, [1, 2, 3, 4], (24,8))
plt.savefig(picpath+data+'_'+resolution+'_wexb_max.pdf', bbox_inches='tight')

# Plot all shearing rate with mean overe intervall in 1000 steps
plot.all_shearingrate_radialcoordinate(rad_coord, wexb, (12,8), 7000)
plt.savefig(picpath+data+'_'+resolution+'_wexb_all.pdf', bbox_inches='tight')

# Time interval to display shearing rate
interval = np.array([[ 500, 1500, 4000],
                     [1000, 2000, 5500]])

# Dimension Subplot
for i in [3,2,1]:
    if interval.shape[1] % i == 0:
        xdim, ydim = i, int(interval.shape[1]/i)
        break

grid_x, grid_y = np.arange(xdim), np.arange(ydim)

fig, ax = plt.subplots(ydim, xdim,figsize=(12*xdim, 8*ydim), sharey=True, sharex=True, squeeze=True)

i = 0

for y in grid_y:
    for x in grid_x:
        start, end = zonalflow.get_index_from_value(time,interval[0][i]) , zonalflow.get_index_from_value(time,interval[1][i])
        start_time, end_time = interval[0][i], interval[1][i]
        
        # Shearing rate with mean over time
        wexb_rad_mean, wexb_rad_middle = zonalflow.get_mean_middle_shearingrate(start, end, wexb)
        # FT{shearing rate}
        wexb_rad_mean_amp, wexb_rad_mean_amp_max = zonalflow.get_fft_mean_max_shearingrate_amplitude(wexb_rad_mean)

        plot.mean_shearingrate_radialcoordinate_subplot(rad_coord, rad_boxsize, wexb_rad_mean, wexb_rad_middle, wexb_rad_mean_amp_max, 
                                                        ax, x, y, xdim, ydim, start_time, end_time)
        
        i += 1

fig.text(0.51, 0.05, r'$\psi~[\rho]$', ha='center')
fig.text(0.1, 0.5, r'$\omega_{\mathrm{E \times B}}$', va='center', rotation='vertical')

plt.subplots_adjust(wspace=0.03, hspace=0.05)
plt.savefig(picpath+data+'_'+resolution+'_wexb_evolution.pdf', bbox_inches='tight')

#Mean shearing rate for time interval in one plot
# Plot parameter
plot.parameters(True, 50, (24,8), 300)

fig, ax = plt.subplots(1, 1,figsize=(18,8))

interval = np.array([[4000],
                     [5000]])

colors = ['#029e73', '#d55e00']

i = 0

while i < len(interval[0]):
    start_time, end_time = interval[0][i], interval[1][i]
    label_time = r' $t \in$ [' + str(start_time) + r', ' + str(end_time) + r']'
    start, end = zonalflow.get_index_from_value(time, start_time) , zonalflow.get_index_from_value(time, end_time)

    # Shearing rate with mean over time
    wexb_rad_mean, wexb_rad_middle = zonalflow.get_mean_middle_shearingrate(start, end, wexb)
    # FT{shearing rate}
    wexb_rad_mean_amp, wexb_rad_mean_amp_max = zonalflow.get_fft_mean_max_shearingrate_amplitude(wexb_rad_mean)

    ax.plot(rad_coord, wexb_rad_mean, linewidth=4, label = label_time, color=colors[i])
    
    i += 1

# linear growth rate
ax.plot(rad_coord,  np.repeat(lineargrowth_rlt, len(rad_coord)), linewidth = 4, linestyle = 'dashed', color = lineargrowth_rlt_color)
ax.plot(rad_coord, -np.repeat(lineargrowth_rlt, len(rad_coord)), linewidth = 4, linestyle = 'dashed', color = lineargrowth_rlt_color)
        
ax.text(1.01, 0.735, r'\boldmath{$+\gamma$}', color = lineargrowth_rlt_color, transform=ax.transAxes)
ax.text(1.01, 0.205, r'\boldmath{$-\gamma$}', color = lineargrowth_rlt_color, transform=ax.transAxes)

ax.legend(loc='upper center', bbox_to_anchor=(1/2, 1.25), frameon=False, ncol = 4, handlelength=1)

ax.set_xlim(0, rad_boxsize)
ax.set_ylim(-0.4, 0.4)

ax.set_yticks(np.arange(-0.4, 0.5, 0.2))

ax.set_xlabel(r'$\psi~[\rho]$')
ax.set_ylabel(r'$\omega_{\mathrm{E \times B}}~[\nu_{\mathrm{th}}/R]$')

plt.savefig(picpath+data+'_'+resolution+'_wexb_selection.pdf', bbox_inches='tight')

# write evaluated data into h5 file

data_eval = [dx, ddphi, zonal_pot, wexb, wexb_max]
groupname = ['evaluation' + '/' + x for x in ['derivative_stepsize', 'second_derivative_phi', 'zonalflow_potential', 
                                              'shearing_rate', 'shearing_rate_maximum']]

h5tools.hdf5_write_data(f, data_eval, groupname)