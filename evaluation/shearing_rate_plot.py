
import sys, os, h5py, ast

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.transforms import Bbox

filepath = os.getcwd()
homepath = filepath.split('evaluation')[0]
sys.path.insert(1, homepath + 'python')

import zonalflow, h5tools, plot

def rotate(l, n):
    return np.concatenate((l[n:],l[:n]))

boxsize, Ns, Nvpar, Nmu = '1x1', 16, 64, 9

datainfopath = homepath + '/data/data.csv'
datainfo = zonalflow.get_data_info(datainfopath, boxsize, Ns = Ns, Nvpar = Nvpar, Nmu = Nmu)

filepath = homepath + datainfo['path'].values[0]

# File import and Create picture folder
data = filepath[filepath.find('data/')+len('data/'):filepath.rfind('/boxsize')]
path = filepath.split(data + '/')[1]

filename = filepath + '/' + 'data.h5'
f = h5py.File(filename,"r+")

picpath = homepath+'pictures/Theory'
## Create target Directory if don't exist

if not os.path.exists(picpath):
    os.makedirs(picpath)

# linear growth rate
lineargrowth = pd.read_csv(homepath + 'data/linear_growthrate_ITG.dat', index_col=0)
try:
    lineargrowth_rlt = lineargrowth['gamma_N'][float(data.split('rlt')[1])]
except KeyError:
    lineargrowth_rlt = 0.21
lineargrowth_rlt_color = 'grey'


eflux_data, time = zonalflow.get_eflux_time(f)
wexb, rad_coord, rad_boxsize, ddphi, dx, zonal_pot = zonalflow.get_shearingrate_radialcoordinate_radialboxsize_ddphi_dx_zonalpot(f)


# Plot parameter
plot.parameters(50, (12,12), 300)

fig, ax = plt.subplots(1, 1)

start_time, end_time = 3000, 4000
start, end = zonalflow.get_index_from_value(time, start_time) , zonalflow.get_index_from_value(time, end_time)

# Shearing rate with mean over time
wexb_rad_mean, wexb_rad_middle = zonalflow.get_mean_middle_shearingrate(start, end, wexb)
wexb_rad_mean = rotate(wexb_rad_mean, -4)
        
ax.plot(rad_coord, wexb_rad_mean, linewidth=4, label=r'fully')


start_time, end_time = 2200, 3000
start, end = zonalflow.get_index_from_value(time, start_time) , zonalflow.get_index_from_value(time, end_time)

# Shearing rate with mean over time
wexb_rad_mean, wexb_rad_middle = zonalflow.get_mean_middle_shearingrate(start, end, wexb)
wexb_rad_mean = rotate(wexb_rad_mean, -45)
        
ax.plot(rad_coord, wexb_rad_mean, linewidth=4, label=r'partially')

# linear growth rate
ax.plot(rad_coord,  np.repeat(lineargrowth_rlt, len(rad_coord)), linewidth = 4, linestyle = 'dashed', color = lineargrowth_rlt_color)
ax.plot(rad_coord, -np.repeat(lineargrowth_rlt, len(rad_coord)), linewidth = 4, linestyle = 'dashed', color = lineargrowth_rlt_color)

ax.text(1.01, 0.735, r'\boldmath{$+\gamma$}', color = lineargrowth_rlt_color, transform=ax.transAxes)
ax.text(1.01, 0.205, r'\boldmath{$-\gamma$}', color = lineargrowth_rlt_color, transform=ax.transAxes)

ax.set_xlim(0, rad_boxsize)
ax.set_ylim(-0.4, 0.4)

ax.set_xticks([x.round(1) for x in np.arange(0, rad_boxsize + rad_boxsize/4, rad_boxsize/4)])
ax.set_yticks(np.arange(-0.4, 0.4 + 0.2, 0.2)[1:-1])

ax.legend(frameon=False, ncol = 4, handlelength=1)

ax.set_xlabel(r'$x~[\rho]$')
ax.set_ylabel(r'$\langle\omega_{\mathrm{E \times B}}\rangle~[\nu_{\mathrm{th}}/R]$')

plt.savefig(picpath+'/'+'Shearing-Rate.pdf', bbox_inches='tight')