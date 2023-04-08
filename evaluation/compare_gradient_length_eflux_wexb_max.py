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

plot.parameters(40, (24,16), 300)

'''
# COMPARE: 6.0 & 6.3 ================================================================================================================================

# File import and Create picture folder
path = ['S6_rlt6.0/boxsize1x1/Ns16/Nvpar64/Nmu9', 'S6_rlt6.3/boxsize1x1/Ns16/Nvpar64/Nmu9']

filename = [homepath + 'data/'+i+'/data.h5' for i in path]
f = [h5py.File(i,"r+") for i in filename]

picDir = homepath + 'pictures/Comparison/Gradient-Length/'
# Create target Directory if don't exist
if not os.path.exists(picDir):
    os.makedirs(picDir)
    
    
# Compare eflux and amplitude in time domain
fig, (ax_eflux, ax_wexb_max) = plt.subplots(2, 1) #, sharex=True)
#fig.suptitle('Comparsion of finite heat flux threshold of $R/L_{\mathrm{T}}$')

resolution = [r'$R/L_{\mathrm{T}}$ = 6.0', r'$R/L_{\mathrm{T}}$ = 6.3']

#ax_eflux.set_title(r'$N_{\mathrm{vpar}}$ = 64; $N_{\mathrm{\mu}}$ = 9', pad=20)
ax_eflux.set_xlabel(r'$t~[R/ \nu_{\mathrm{th}}]$')
ax_eflux.set_ylabel(r'$\chi~[\rho^2 \nu_{\mathrm{th}} / R]$')
#ax_eflux.yaxis.set_label_coords(-0.09,0.5)

#ax_wexb_max.set_title(r'$N_{\mathrm{vpar}}$ = 64; $N_{\mathrm{\mu}}$ = 9', pad=20)
ax_wexb_max.set_xlabel(r'$t~[R/ \nu_{\mathrm{th}}]$')
ax_wexb_max.set_ylabel(r'$|k_1^2 \phi|$')
#ax_wexb_max.yaxis.set_label_coords(-0.09,0.5)

max_index = 10000

for i, n in zip(f, resolution):
    #eflux
    eflux, time = zonalflow.get_eflux_time(i)
    
    eflux, time = eflux[:max_index], time[:max_index]

    ax_eflux.plot(time, eflux, label=n)
    
    plot.ax_ticks_subplot(ax_eflux)
    
    xmax = time[-2]
    
    ax_eflux.set_xlim(xmin=0, xmax=6000)
    ax_eflux.set_ylim(ymin=0, ymax=25)
    
    ax_eflux.legend(loc='upper center', ncol=4, frameon=False)
    
    #wexb_max
    wexb, rad_coord, rad_boxsize, ddphi, dx, zonal_pot = zonalflow.get_shearingrate_radialcoordinate_radialboxsize_ddphi_dx_zonalpot(i)
    wexb_max = zonalflow.get_max_shearingrate(i, wexb, time, 1)
    
    ax_wexb_max.plot(time, wexb_max[1][:max_index], label=n)
    
    plot.ax_ticks_subplot(ax_wexb_max)
    
    ax_wexb_max.set_xlim(xmin=0, xmax=xmax)
    ax_wexb_max.set_ylim(ymin=0, ymax=0.30)
    ax_wexb_max.yaxis.set_ticks(np.arange(0, 0.40, 0.1))
    
    #ax_wexb_max.legend(loc='lower right')

#plt.subplots_adjust(wspace=0, hspace=0)
plt.subplots_adjust(hspace=0.5)

plot.savefig_subplot(fig, ax_eflux   , picDir + '/S6_rlt6.0-6.3_boxsize1x1_Ns16_Nvpar64_Nmu9_eflux_comparison.pdf'   , pad=0.02)
plot.savefig_subplot(fig, ax_wexb_max, picDir + '/S6_rlt6.0-6.3_boxsize1x1_Ns16_Nvpar64_Nmu9_wexb_max_comparison.pdf', pad=0.02)

ax_eflux.legend(loc='upper center', ncol=4, frameon=False)


ax_eflux.text(-0.11, -0.17, r'\bf{(a)}', transform=ax_eflux.transAxes)
ax_wexb_max.text(-0.11, -0.17, r'\bf{(b)}', transform=ax_wexb_max.transAxes)

plt.savefig(picDir + '/S6_rlt6.0-6.3_boxsize1x1_Ns16_Nvpar64_Nmu9_comparison.pdf', bbox_inches='tight')
#'''

'''
# COMPARE: 6.0 & 6.2 ================================================================================================================================

# File import and Create picture folder
path = ['S6_rlt6.0/boxsize2x2/Ns16/Nvpar48/Nmu9', 'S6_rlt6.2/boxsize2x2/Ns16/Nvpar64/Nmu9']

filename = [homepath + 'data/'+i+'/data.h5' for i in path]
f = [h5py.File(i,"r+") for i in filename]

picDir = homepath + 'pictures/Comparison/Gradient-Length/'
# Create target Directory if don't exist
if not os.path.exists(picDir):
    os.makedirs(picDir)
    
# Compare eflux and amplitude in time domain
fig, (ax_eflux, ax_wexb_max) = plt.subplots(1, 2, figsize = (24,6)) #, sharex=True)
#fig.suptitle('Comparsion of finite heat flux threshold of $R/L_{\mathrm{T}}$')

resolution = [r'$R/L_{\mathrm{T}}$ = 6.0', r'$R/L_{\mathrm{T}}$ = 6.2']

#ax_eflux.set_title(r'$N_{\mathrm{vpar}}$ = 64; $N_{\mathrm{\mu}}$ = 9', pad=20)
ax_eflux.set_xlabel(r'$t~[R/ \nu_{\mathrm{th}}]$')
ax_eflux.set_ylabel(r'$\chi~[\rho^2 \nu_{\mathrm{th}} / R]$')
ax_eflux.yaxis.set_label_coords(-0.09,0.5)

#ax_wexb_max.set_title(r'$N_{\mathrm{vpar}}$ = 64; $N_{\mathrm{\mu}}$ = 9', pad=20)
ax_wexb_max.set_xlabel(r'$t~[R/ \nu_{\mathrm{th}}]$')
ax_wexb_max.set_ylabel(r'$|k_2^2 \phi|$')
ax_wexb_max.yaxis.set_label_coords(-0.09,0.5)

max_index = None

for i, n in zip(f, resolution):
    #eflux
    eflux, time = zonalflow.get_eflux_time(i)
    
    eflux, time = eflux[:max_index], time[:max_index]

    ax_eflux.plot(time, eflux, label=n, linewidth = 2.5)
    
    plot.ax_ticks_subplot(ax_eflux)
    
    xmax = time[-2]
    
    ax_eflux.set_xlim(xmin=0, xmax=xmax)
    ax_eflux.set_ylim(ymin=0, ymax=20)
    
    ax_eflux.legend(loc='upper center', bbox_to_anchor=(0.5, 1.2), ncol=4, frameon=False)
    
    #wexb_max
    wexb, rad_coord, rad_boxsize, ddphi, dx, zonal_pot = zonalflow.get_shearingrate_radialcoordinate_radialboxsize_ddphi_dx_zonalpot(i)
    wexb_max = zonalflow.get_max_shearingrate(i, wexb, time, 1)
    
    ax_wexb_max.plot(time, wexb_max[2][:max_index], label=n, linewidth = 2.5)
    
    plot.ax_ticks_subplot(ax_wexb_max)
    
    ax_wexb_max.set_xlim(xmin=0, xmax=xmax)
    ax_wexb_max.set_ylim(ymin=0, ymax=0.30)
    ax_wexb_max.yaxis.set_ticks(np.arange(0, 0.40, 0.1))
    
    #ax_wexb_max.legend(loc='lower right')

#plt.subplots_adjust(wspace=0, hspace=0)
plt.subplots_adjust(top=0.9, wspace=0.2, hspace=0.2)

plot.savefig_subplot(fig, ax_eflux   , picDir + '/S6_rlt6.0-6.2_boxsize2x2_Ns16_Nvpar48-64_Nmu9_eflux_comparison.pdf'   , pad=0.02)
plot.savefig_subplot(fig, ax_wexb_max, picDir + '/S6_rlt6.0-6.2_boxsize2x2_Ns16_Nvpar48-64_Nmu9_wexb_max_comparison.pdf', pad=0.02)

ax_eflux.legend(loc='upper center', bbox_to_anchor=(1.1, 1.2), ncol=4, frameon=False)


ax_eflux.text(-0.11, -0.17, r'\bf{(a)}', transform=ax_eflux.transAxes)
ax_wexb_max.text(-0.11, -0.17, r'\bf{(b)}', transform=ax_wexb_max.transAxes)

plt.savefig(picDir + '/S6_rlt6.0-6.2_boxsize2x2_Ns16_Nvpar48-64_Nmu9_comparison.pdf', bbox_inches='tight')
#'''