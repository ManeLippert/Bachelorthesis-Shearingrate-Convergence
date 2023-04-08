# MODULES ===========================================================================================================================================

import sys, os, h5py
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.transforms import Bbox

# PATH ==============================================================================================================================================

filepath = os.getcwd()
homepath = filepath.split('evaluation')[0]
sys.path.insert(1, homepath + 'python')

import zonalflow, h5tools, derivative, plot

# PARAMETER =========================================================================================================================================

plot.parameters(42, (24,8), 300, linewidth = 4, tickwidth = 2)

# FUNCTIONS =========================================================================================================================================

def eflux_compare_plot(files, label, xlim = (0, None), ylim = (0, None), invert_legend = True):
    
    fig, ax = plt.subplots()
    
    for f, n in zip(files, label):
        eflux_data, time = zonalflow.get_eflux_time(f)
        
        plot.eflux_time(time, eflux_data, axis = ax, xlim = xlim, ylim = ylim, create_plot = False, label = n)
    
    handles_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
    handles, labels = [sum(lol, []) for lol in zip(*handles_labels)]
    
    if invert_legend:
        ax.legend(handles[::-1], labels[::-1], ncol = len(files))
    else:
        ax.legend(handles, labels, ncol = len(files))

# R/L_T 6.0 & 6.3 ===================================================================================================================================

# File import and Create picture folder
path = ['S6_rlt6.0/boxsize1x1/Ns16/Nvpar64/Nmu9', 'S6_rlt6.3/boxsize1x1/Ns16/Nvpar64/Nmu9']

filename = [homepath + 'data/'+i+'/data.h5' for i in path]
files = [h5py.File(i,"r+") for i in filename]

picDir = homepath + 'pictures/Comparison/Gradient-Length/'
# Create target Directory if don't exist
if not os.path.exists(picDir):
    os.makedirs(picDir)

resolution = [r'$R/L_{\mathrm{T}}$ = 6.0', r'$R/L_{\mathrm{T}}$ = 6.3']

eflux_compare_plot(files, resolution, xlim=(0, 6000), ylim=(0,25))

plt.savefig(picDir + '/S6_rlt6.0-6.3_boxsize1x1_Ns16_Nvpar64_Nmu9_eflux_comparison.pdf', bbox_inches='tight')

# R/L_T 6.0 & 6.2 & 6.4  ============================================================================================================================

# File import and Create picture folder
path = ['S6_rlt6.0/boxsize3x3/Ns16/Nvpar48/Nmu9', 'S6_rlt6.2/boxsize3x3/Ns16/Nvpar48/Nmu9', 'S6_rlt6.4/boxsize3x3/Ns16/Nvpar48/Nmu9']

filename = [homepath + 'data/'+i+'/data.h5' for i in path]
files = [h5py.File(i,"r+") for i in filename]

picDir = homepath + 'pictures/Comparison/Gradient-Length/'
# Create target Directory if don't exist
if not os.path.exists(picDir):
    os.makedirs(picDir)

resolution = [r'$R/L_{\mathrm{T}}$ = 6.0', r'$R/L_{\mathrm{T}}$ = 6.2',  r'$R/L_{\mathrm{T}}$ = 6.4']

eflux_compare_plot(files, resolution, xlim=(0, 12000), ylim=(0,25), invert_legend=False)

plt.savefig(picDir + '/S6_rlt6.0-6.2-6.4_boxsize3x3_Ns16_Nvpar48_Nmu9_eflux_comparison.pdf', bbox_inches='tight')

# Ns 12, 16 =========================================================================================================================================

# File import and Create picture folder
path = ['S6_rlt6.0/boxsize1x1/Ns16/Nvpar64/Nmu9', 'S6_rlt6.0/boxsize1x1/Ns12/Nvpar64/Nmu9']

filename = [homepath + 'data/'+i+'/data.h5' for i in path]
files = [h5py.File(i,"r+") for i in filename]

picDir = homepath + 'pictures/Comparison/Resolution/Ns/'
# Create target Directory if don't exist
if not os.path.exists(picDir):
    os.makedirs(picDir)

resolution = [r'$N_{\mathrm{s}}$ = 16', r'$N_{\mathrm{s}}$ = 12']

eflux_compare_plot(files, resolution, xlim=(0, 6000), ylim=(0,25))

plt.savefig(picDir + '/S6_rlt6.0_boxsize1x1_Ns12-16_Nvpar64_Nmu9_eflux_comparison.pdf', bbox_inches='tight')

# Nvpar 16, 32, 48, 64 ==============================================================================================================================

# File import and Create picture folder
path = ['S6_rlt6.0/boxsize1x1/Ns16/Nvpar64/Nmu9', 'S6_rlt6.0/boxsize1x1/Ns16/Nvpar48/Nmu9',
        'S6_rlt6.0/boxsize1x1/Ns16/Nvpar32/Nmu9', 'S6_rlt6.0/boxsize1x1/Ns16/Nvpar16/Nmu9']

filename = [homepath + 'data/'+i+'/data.h5' for i in path]
files = [h5py.File(i,"r+") for i in filename]

picDir = homepath + 'pictures/Comparison/Resolution/Nvpar/'
# Create target Directory if don't exist
if not os.path.exists(picDir):
    os.makedirs(picDir)

resolution = [r'$N_{\nu_\parallel}$ = 64', r'$N_{\nu_\parallel}$ = 48',
              r'$N_{\nu_\parallel}$ = 32', r'$N_{\nu_\parallel}$ = 16']

eflux_compare_plot(files, resolution, xlim=(0, 6000), ylim=(0,25))

plt.savefig(picDir + '/S6_rlt6.0_boxsize1x1_Ns16_Nvpar16-32-48-64_Nmu9_eflux_comparison.pdf', bbox_inches='tight')

# Nmu 6, 9 ==========================================================================================================================================

# File import and Create picture folder
path = ['S6_rlt6.0/boxsize1x1/Ns16/Nvpar64/Nmu9', 'S6_rlt6.0/boxsize1x1/Ns16/Nvpar64/Nmu6']

filename = [homepath + 'data/'+i+'/data.h5' for i in path]
files = [h5py.File(i,"r+") for i in filename]

picDir = homepath + 'pictures/Comparison/Resolution/Nmu/'
# Create target Directory if don't exist
if not os.path.exists(picDir):
    os.makedirs(picDir)

resolution = [r'$N_{\mu}$ = 9', r'$N_{\mu}$ = 6']

eflux_compare_plot(files, resolution, xlim=(0, 6000), ylim=(0,25))

plt.savefig(picDir + '/S6_rlt6.0_boxsize1x1_Ns16_Nvpar64_Nmu6-9_eflux_comparison.pdf', bbox_inches='tight')

# File import and Create picture folder
path = ['S6_rlt6.0/boxsize1x1/Ns16/Nvpar48/Nmu9', 'S6_rlt6.0/boxsize1x1/Ns16/Nvpar48/Nmu6']

filename = [homepath + 'data/'+i+'/data.h5' for i in path]
files = [h5py.File(i,"r+") for i in filename]

picDir = homepath + 'pictures/Comparison/Resolution/Nmu/'
# Create target Directory if don't exist
if not os.path.exists(picDir):
    os.makedirs(picDir)

resolution = [r'$N_{\mu}$ = 9', r'$N_{\mu}$ = 6']

eflux_compare_plot(files, resolution, xlim=(0, 6000), ylim=(0,25))

plt.savefig(picDir + '/S6_rlt6.0_boxsize1x1_Ns16_Nvpar48_Nmu6-9_eflux_comparison.pdf', bbox_inches='tight')