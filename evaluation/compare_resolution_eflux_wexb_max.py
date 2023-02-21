# MODULES ===========================================================================================================================================

# Import modules
import sys, os, h5py
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.transforms import Bbox

# PATH ==============================================================================================================================================

filepath = os.getcwd()
homepath = filepath.split('evaluation')[0]
sys.path.insert(1, homepath + 'python')

import zonalflow, h5tools, plot

# PARAMETER =========================================================================================================================================

plot.parameters(True, 22, (24,8), 300)

# NS GRID ===========================================================================================================================================

## NVPARr=64 & NMU=9 =================================================================================================================================

# File import and Create picture folder
data = 'S6_rlt6.0'
path = ['boxsize1x1/Ns12/Nvpar64/Nmu9', 'boxsize1x1/Ns16/Nvpar64/Nmu9']
compare = 'Ns'

filename = [homepath + 'data/'+data+'/'+i+'/data.h5' for i in path]
f = [h5py.File(i,"r+") for i in filename]

picDir = homepath + 'pictures/Comparison/Resolution/' + compare
# Create target Directory if don't exist
if not os.path.exists(picDir):
    os.makedirs(picDir)
    
# Compare eflux and amplitude in time domain
fig, (ax_eflux, ax_wexb_max) = plt.subplots(2, 1, figsize = (24,16)) #, sharex=True)
#fig.suptitle('Comparsion of resolution of $N_{\mathrm{s}}$ grid size')

resolution = [r'$N_{\mathrm{s}}$ = 12', r'$N_{\mathrm{s}}$ = 16']

#ax_eflux.set_title(r'$N_{\mathrm{vpar}}$ = 64; $N_{\mathrm{\mu}}$ = 9', pad=20)
ax_eflux.set_xlabel(r'$t~[R/ \nu_{\mathrm{th}}]$')
ax_eflux.set_ylabel(r'$\chi~[\rho^2 \nu_{\mathrm{th}} / R]$')
ax_eflux.yaxis.set_label_coords(-0.035,0.5)

#ax_wexb_max.set_title(r'$N_{\mathrm{vpar}}$ = 64; $N_{\mathrm{\mu}}$ = 9', pad=20)
ax_wexb_max.set_xlabel(r'$t~[R/ \nu_{\mathrm{th}}]$')
ax_wexb_max.set_ylabel(r'$|k_1^2 \phi|$')
ax_wexb_max.yaxis.set_label_coords(-0.035,0.5)

max_index = 10000

for i, n in zip(f, resolution):
    #eflux
    eflux, time = zonalflow.get_eflux_time(i)
    
    eflux, time = eflux[:max_index], time[:max_index]

    ax_eflux.plot(time, eflux, label=n)
    
    plot.ax_ticks_subplot(ax_eflux)
    
    xmax = time[-2]
    
    ax_eflux.set_xlim(xmin=0, xmax=xmax)
    ax_eflux.set_ylim(ymin=0, ymax=30)
    
    ax_eflux.legend(loc=(0.87, 0.8))
    
    #wexb_max
    wexb, rad_coord, rad_boxsize, ddphi, dx, zonal_pot = zonalflow.get_shearingrate_radialcoordinate_radialboxsize_ddphi_dx_zonalpot(i)
    wexb_max = zonalflow.get_max_shearingrate(i, wexb, time, 1)
    
    ax_wexb_max.plot(time, wexb_max[1][:max_index], label=n)
    
    plot.ax_ticks_subplot(ax_wexb_max)
    
    ax_wexb_max.set_xlim(xmin=0, xmax=xmax)
    ax_wexb_max.set_ylim(ymin=0, ymax=0.30)
    
    #ax_wexb_max.legend(loc='lower right')

#plt.subplots_adjust(wspace=0, hspace=0)
plt.subplots_adjust(top=0.9, wspace=0.4, hspace=0.3)

plot.savefig_subplot(fig, ax_eflux   , picDir + '/S6_rlt6.0_boxsize1x1_Ns12-16_Nvpar64_Nmu9_eflux_comparison.pdf'   , pad=0.02)
plot.savefig_subplot(fig, ax_wexb_max, picDir + '/S6_rlt6.0_boxsize1x1_Ns12-16_Nvpar64_Nmu9_wexb_max_comparison.pdf', pad=0.02)

ax_eflux.text(-0.05, -0.09, r'\bf{(a)}', transform=ax_eflux.transAxes)
ax_wexb_max.text(-0.05, -0.09, r'\bf{(b)}', transform=ax_wexb_max.transAxes)

plt.savefig(picDir + '/S6_rlt6.0_boxsize1x1_Ns12-16_Nvpar64_Nmu9_comparison.pdf', bbox_inches='tight')

## NVPAR=48 & NMU=9 =================================================================================================================================

# File import and Create picture folder
data = 'S6_rlt6.0'
path = ['boxsize1x1/Ns12/Nvpar48/Nmu9', 'boxsize1x1/Ns16/Nvpar48/Nmu9']
compare = 'Ns'

filename = [homepath + 'data/'+data+'/'+i+'/data.h5' for i in path]
f = [h5py.File(i,"r+") for i in filename]

picDir = homepath + 'pictures/Comparison/Resolution/' + compare
# Create target Directory if don't exist
if not os.path.exists(picDir):
    os.makedirs(picDir)
    
# Compare eflux and amplitude in time domain
fig, (ax_eflux, ax_wexb_max) = plt.subplots(2, 1, figsize = (24,16)) #, sharex=True)
#fig.suptitle('Comparsion of resolution of $N_{\mathrm{s}}$ grid size')

resolution = [r'$N_{\mathrm{s}}$ = 12', r'$N_{\mathrm{s}}$ = 16']

#ax_eflux.set_title(r'$N_{\mathrm{vpar}}$ = 64; $N_{\mathrm{\mu}}$ = 9', pad=20)
ax_eflux.set_xlabel(r'$t~[R/ \nu_{\mathrm{th}}]$')
ax_eflux.set_ylabel(r'$\chi~[\rho^2 \nu_{\mathrm{th}} / R]$')
ax_eflux.yaxis.set_label_coords(-0.035,0.5)

#ax_wexb_max.set_title(r'$N_{\mathrm{vpar}}$ = 64; $N_{\mathrm{\mu}}$ = 9', pad=20)
ax_wexb_max.set_xlabel(r'$t~[R/ \nu_{\mathrm{th}}]$')
ax_wexb_max.set_ylabel(r'$|k_1^2 \phi|$')
ax_wexb_max.yaxis.set_label_coords(-0.035,0.5)

max_index = 10000

for i, n in zip(f, resolution):
    #eflux
    eflux, time = zonalflow.get_eflux_time(i)
    
    eflux, time = eflux[:max_index], time[:max_index]

    ax_eflux.plot(time, eflux, label=n)
    
    plot.ax_ticks_subplot(ax_eflux)
    
    xmax = time[-2]
    
    ax_eflux.set_xlim(xmin=0, xmax=xmax)
    ax_eflux.set_ylim(ymin=0, ymax=30)
    
    ax_eflux.legend(loc=(0.87, 0.8))
    
    #wexb_max
    wexb, rad_coord, rad_boxsize, ddphi, dx, zonal_pot = zonalflow.get_shearingrate_radialcoordinate_radialboxsize_ddphi_dx_zonalpot(i)
    wexb_max = zonalflow.get_max_shearingrate(i, wexb, time, 1)
    
    ax_wexb_max.plot(time, wexb_max[1][:max_index], label=n)
    
    plot.ax_ticks_subplot(ax_wexb_max)
    
    ax_wexb_max.set_xlim(xmin=0, xmax=xmax)
    ax_wexb_max.set_ylim(ymin=0, ymax=0.30)
    
    #ax_wexb_max.legend(loc='lower right')

#plt.subplots_adjust(wspace=0, hspace=0)
plt.subplots_adjust(top=0.9, wspace=0.4, hspace=0.3)

plot.savefig_subplot(fig, ax_eflux   , picDir + '/S6_rlt6.0_boxsize1x1_Ns12-16_Nvpar48_Nmu9_eflux_comparison.pdf'   , pad=0.02)
plot.savefig_subplot(fig, ax_wexb_max, picDir + '/S6_rlt6.0_boxsize1x1_Ns12-16_Nvpar48_Nmu9_wexb_max_comparison.pdf', pad=0.02)

ax_eflux.text(-0.05, -0.09, r'\bf{(a)}', transform=ax_eflux.transAxes)
ax_wexb_max.text(-0.05, -0.09, r'\bf{(b)}', transform=ax_wexb_max.transAxes)

plt.savefig(picDir + '/S6_rlt6.0_boxsize1x1_Ns12-16_Nvpar48_Nmu9_comparison.pdf', bbox_inches='tight')

# NVPAR GRID ========================================================================================================================================

## NS=16 & NMU=9 ====================================================================================================================================


# File import and Create picture folder
data = 'S6_rlt6.0'
path = ['boxsize1x1/Ns16/Nvpar16/Nmu9', 'boxsize1x1/Ns16/Nvpar32/Nmu9', 
        'boxsize1x1/Ns16/Nvpar48/Nmu9', 'boxsize1x1/Ns16/Nvpar64/Nmu9']
compare = 'Nvpar'

filename = [homepath + 'data/'+data+'/'+i+'/data.h5' for i in path]
f = [h5py.File(i,"r+") for i in filename]

picDir = homepath + 'pictures/Comparison/Resolution/'+ compare
# Create target Directory if don't exist
if not os.path.exists(picDir):
    os.makedirs(picDir)
    
# Compare eflux and amplitude in time domain
fig, (ax_eflux, ax_wexb_max) = plt.subplots(2, 1, figsize = (24,16)) #, sharex=True)
#fig.suptitle('Comparsion of resolution of $N_{\mathrm{vpar}}$ grid size')

resolution = [r'$N_{\mathrm{vpar}}$ = 16', r'$N_{\mathrm{vpar}}$ = 32', r'$N_{\mathrm{vpar}}$ = 48', r'$N_{\mathrm{vpar}}$ = 64']

#ax_eflux.set_title(r'$N_{\mathrm{vpar}}$ = 64; $N_{\mathrm{\mu}}$ = 9', pad=20)
ax_eflux.set_xlabel(r'$t~[R/ \nu_{\mathrm{th}}]$')
ax_eflux.set_ylabel(r'$\chi~[\rho^2 \nu_{\mathrm{th}} / R]$')
ax_eflux.yaxis.set_label_coords(-0.035,0.5)

#ax_wexb_max.set_title(r'$N_{\mathrm{vpar}}$ = 64; $N_{\mathrm{\mu}}$ = 9', pad=20)
ax_wexb_max.set_xlabel(r'$t~[R/ \nu_{\mathrm{th}}]$')
ax_wexb_max.set_ylabel(r'$|k_1^2 \phi|$')
ax_wexb_max.yaxis.set_label_coords(-0.035,0.5)

max_index = 10000

for i, n in zip(f, resolution):
    #eflux
    eflux, time = zonalflow.get_eflux_time(i)
    
    eflux, time = eflux[:max_index], time[:max_index]

    ax_eflux.plot(time, eflux, label=n)
    
    plot.ax_ticks_subplot(ax_eflux)
    
    xmax = time[-2]
    
    ax_eflux.set_xlim(xmin=0, xmax=xmax)
    ax_eflux.set_ylim(ymin=0, ymax=30)
    
    ax_eflux.legend(loc=(0.855, 0.62))
    
    #wexb_max
    wexb, rad_coord, rad_boxsize, ddphi, dx, zonal_pot = zonalflow.get_shearingrate_radialcoordinate_radialboxsize_ddphi_dx_zonalpot(i)
    wexb_max = zonalflow.get_max_shearingrate(i, wexb, time, 1)
    
    ax_wexb_max.plot(time, wexb_max[1][:max_index], label=n)
    
    plot.ax_ticks_subplot(ax_wexb_max)
    
    ax_wexb_max.set_xlim(xmin=0, xmax=xmax)
    ax_wexb_max.set_ylim(ymin=0, ymax=0.30)
    
    #ax_wexb_max.legend(loc='lower right')

#plt.subplots_adjust(wspace=0, hspace=0)
plt.subplots_adjust(top=0.9, wspace=0.4, hspace=0.3)

plot.savefig_subplot(fig, ax_eflux   , picDir + '/S6_rlt6.0_boxsize1x1_Ns16_Nvpar16-32-48-64_Nmu9_eflux_comparison.pdf'   , pad=0.02)
plot.savefig_subplot(fig, ax_wexb_max, picDir + '/S6_rlt6.0_boxsize1x1_Ns16_Nvpar16-32-48-64_Nmu9_wexb_max_comparison.pdf', pad=0.02)

ax_eflux.text(-0.05, -0.09, r'\bf{(a)}', transform=ax_eflux.transAxes)
ax_wexb_max.text(-0.05, -0.09, r'\bf{(b)}', transform=ax_wexb_max.transAxes)

plt.savefig(picDir + '/S6_rlt6.0_boxsize1x1_Ns16_Nvpar16-32-48-64_Nmu9_comparison.pdf', bbox_inches='tight')

# NMU GRID ==========================================================================================================================================

## NS=16 & NVPAR=64 =================================================================================================================================


# File import and Create picture folder
data = 'S6_rlt6.0'
path = ['boxsize1x1/Ns16/Nvpar64/Nmu6', 'boxsize1x1/Ns16/Nvpar64/Nmu9']
compare = 'Nmu'

filename = [homepath + 'data/'+data+'/'+i+'/data.h5' for i in path]
f = [h5py.File(i,"r+") for i in filename]

picDir = homepath + 'pictures/Comparison/Resolution/'+ compare
# Create target Directory if don't exist
if not os.path.exists(picDir):
    os.makedirs(picDir)
    
# Compare eflux and amplitude in time domain
fig, (ax_eflux, ax_wexb_max) = plt.subplots(2, 1, figsize = (24,16)) #, sharex=True)
#fig.suptitle('Comparsion of resolution of $N_{\mathrm{\mu}}$ grid size')

resolution = [r'$N_{\mathrm{\mu}}$ = 6', r'$N_{\mathrm{\mu}}$ = 9']

#ax_eflux.set_title(r'$N_{\mathrm{vpar}}$ = 64; $N_{\mathrm{\mu}}$ = 9', pad=20)
ax_eflux.set_xlabel(r'$t~[R/ \nu_{\mathrm{th}}]$')
ax_eflux.set_ylabel(r'$\chi~[\rho^2 \nu_{\mathrm{th}} / R]$')
ax_eflux.yaxis.set_label_coords(-0.035,0.5)

#ax_wexb_max.set_title(r'$N_{\mathrm{vpar}}$ = 64; $N_{\mathrm{\mu}}$ = 9', pad=20)
ax_wexb_max.set_xlabel(r'$t~[R/ \nu_{\mathrm{th}}]$')
ax_wexb_max.set_ylabel(r'$|k_1^2 \phi|$')
ax_wexb_max.yaxis.set_label_coords(-0.035,0.5)

max_index = 10000

for i, n in zip(f, resolution):
    #eflux
    eflux, time = zonalflow.get_eflux_time(i)
    
    eflux, time = eflux[:max_index], time[:max_index]

    ax_eflux.plot(time, eflux, label=n)
    
    plot.ax_ticks_subplot(ax_eflux)
    
    xmax = time[-2]
    
    ax_eflux.set_xlim(xmin=0, xmax=xmax)
    ax_eflux.set_ylim(ymin=0, ymax=30)
    
    ax_eflux.legend(loc=(0.875, 0.79))
    
    #wexb_max
    wexb, rad_coord, rad_boxsize, ddphi, dx, zonal_pot = zonalflow.get_shearingrate_radialcoordinate_radialboxsize_ddphi_dx_zonalpot(i)
    wexb_max = zonalflow.get_max_shearingrate(i, wexb, time, 1)
    
    ax_wexb_max.plot(time, wexb_max[1][:max_index], label=n)
    
    plot.ax_ticks_subplot(ax_wexb_max)
    
    ax_wexb_max.set_xlim(xmin=0, xmax=xmax)
    ax_wexb_max.set_ylim(ymin=0, ymax=0.35)
    
    #ax_wexb_max.legend(loc='lower right')

#plt.subplots_adjust(wspace=0, hspace=0)
plt.subplots_adjust(top=0.9, wspace=0.4, hspace=0.3)

plot.savefig_subplot(fig, ax_eflux   , picDir + '/S6_rlt6.0_boxsize1x1_Ns16_Nvpar64_Nmu6-9_eflux_comparison.pdf'   , pad=0.02)
plot.savefig_subplot(fig, ax_wexb_max, picDir + '/S6_rlt6.0_boxsize1x1_Ns16_Nvpar64_Nmu6-9_wexb_max_comparison.pdf', pad=0.02)

ax_eflux.text(-0.05, -0.09, r'\bf{(a)}', transform=ax_eflux.transAxes)
ax_wexb_max.text(-0.05, -0.09, r'\bf{(b)}', transform=ax_wexb_max.transAxes)

plt.savefig(picDir + '/S6_rlt6.0_boxsize1x1_Ns16_Nvpar64_Nmu6-9_comparison.pdf', bbox_inches='tight')

## NS=16 & NVPAR=48 =================================================================================================================================

# File import and Create picture folder
data = 'S6_rlt6.0'
path = ['boxsize1x1/Ns16/Nvpar48/Nmu6', 'boxsize1x1/Ns16/Nvpar48/Nmu9']
compare = 'Nmu'

filename = [homepath + 'data/'+data+'/'+i+'/data.h5' for i in path]
f = [h5py.File(i,"r+") for i in filename]

picDir = homepath + 'pictures/Comparison/Resolution/'+ compare
# Create target Directory if don't exist
if not os.path.exists(picDir):
    os.makedirs(picDir)
    
# Compare eflux and amplitude in time domain
fig, (ax_eflux, ax_wexb_max) = plt.subplots(2, 1, figsize = (24,16)) #, sharex=True)
#fig.suptitle('Comparsion of resolution of $N_{\mathrm{\mu}}$ grid size')

resolution = [r'$N_{\mathrm{\mu}}$ = 6', r'$N_{\mathrm{\mu}}$ = 9']

#ax_eflux.set_title(r'$N_{\mathrm{vpar}}$ = 64; $N_{\mathrm{\mu}}$ = 9', pad=20)
ax_eflux.set_xlabel(r'$t~[R/ \nu_{\mathrm{th}}]$')
ax_eflux.set_ylabel(r'$\chi~[\rho^2 \nu_{\mathrm{th}} / R]$')
ax_eflux.yaxis.set_label_coords(-0.035,0.5)

#ax_wexb_max.set_title(r'$N_{\mathrm{vpar}}$ = 64; $N_{\mathrm{\mu}}$ = 9', pad=20)
ax_wexb_max.set_xlabel(r'$t~[R/ \nu_{\mathrm{th}}]$')
ax_wexb_max.set_ylabel(r'$|k_1^2 \phi|$')
ax_wexb_max.yaxis.set_label_coords(-0.035,0.5)

max_index = 10000

for i, n in zip(f, resolution):
    #eflux
    eflux, time = zonalflow.get_eflux_time(i)
    
    eflux, time = eflux[:max_index], time[:max_index]

    ax_eflux.plot(time, eflux, label=n)
    
    plot.ax_ticks_subplot(ax_eflux)
    
    xmax = time[-2]
    
    ax_eflux.set_xlim(xmin=0, xmax=xmax)
    ax_eflux.set_ylim(ymin=0, ymax=30)
    
    ax_eflux.legend(loc=(0.875, 0.79))
    
    #wexb_max
    wexb, rad_coord, rad_boxsize, ddphi, dx, zonal_pot = zonalflow.get_shearingrate_radialcoordinate_radialboxsize_ddphi_dx_zonalpot(i)
    wexb_max = zonalflow.get_max_shearingrate(i, wexb, time, 1)
    
    ax_wexb_max.plot(time, wexb_max[1][:max_index], label=n)
    
    plot.ax_ticks_subplot(ax_wexb_max)
    
    ax_wexb_max.set_xlim(xmin=0, xmax=xmax)
    ax_wexb_max.set_ylim(ymin=0, ymax=0.30)
    
    #ax_wexb_max.legend(loc='lower right')

#plt.subplots_adjust(wspace=0, hspace=0)
plt.subplots_adjust(top=0.9, wspace=0.4, hspace=0.3)

plot.savefig_subplot(fig, ax_eflux   , picDir + '/S6_rlt6.0_boxsize1x1_Ns16_Nvpar48_Nmu6-9_eflux_comparison.pdf'   , pad=0.02)
plot.savefig_subplot(fig, ax_wexb_max, picDir + '/S6_rlt6.0_boxsize1x1_Ns16_Nvpar48_Nmu6-9_wexb_max_comparison.pdf', pad=0.02)

ax_eflux.text(-0.05, -0.09, r'\bf{(a)}', transform=ax_eflux.transAxes)
ax_wexb_max.text(-0.05, -0.09, r'\bf{(b)}', transform=ax_wexb_max.transAxes)

plt.savefig(picDir + '/S6_rlt6.0_boxsize1x1_Ns16_Nvpar48_Nmu6-9_comparison.pdf', bbox_inches='tight')

# DTIM ==============================================================================================================================================

## NS=16 & NVPAR=48 & NMU=9 =========================================================================================================================

# File import and Create picture folder
data = 'S6_rlt6.0'
path = ['boxsize1x1/Ns16/Nvpar48/Nmu9', 'boxsize1x1/Ns16/Nvpar48/Nmu9/dtim0.025']
compare = 'dtim'

filename = [homepath + 'data/'+data+'/'+i+'/data.h5' for i in path]
f = [h5py.File(i,"r+") for i in filename]

picDir = homepath + 'pictures/Comparison/Resolution/'+ compare
# Create target Directory if don't exist
if not os.path.exists(picDir):
    os.makedirs(picDir)
    
# Compare eflux and amplitude in time domain
fig, (ax_eflux, ax_wexb_max) = plt.subplots(2, 1, figsize = (24,16)) #, sharex=True)
#fig.suptitle('Comparsion of resolution of $d_{\mathrm{tim}}$ timstep')

resolution = [r'$d_{\mathrm{tim}}$ = 0.020', r'$d_{\mathrm{tim}}$ = 0.025']

#ax_eflux.set_title(r'$N_{\mathrm{vpar}}$ = 64; $N_{\mathrm{\mu}}$ = 9', pad=20)
ax_eflux.set_xlabel(r'$t~[R/ \nu_{\mathrm{th}}]$')
ax_eflux.set_ylabel(r'$\chi~[\rho^2 \nu_{\mathrm{th}} / R]$')
ax_eflux.yaxis.set_label_coords(-0.035,0.5)

#ax_wexb_max.set_title(r'$N_{\mathrm{vpar}}$ = 64; $N_{\mathrm{\mu}}$ = 9', pad=20)
ax_wexb_max.set_xlabel(r'$t~[R/ \nu_{\mathrm{th}}]$')
ax_wexb_max.set_ylabel(r'$|k_1^2 \phi|$')
ax_wexb_max.yaxis.set_label_coords(-0.035,0.5)

max_index = 9000

for i, n in zip(f, resolution):
    #eflux
    eflux, time = zonalflow.get_eflux_time(i)
    
    eflux, time = eflux[:max_index], time[:max_index]

    ax_eflux.plot(time, eflux, label=n)
    
    plot.ax_ticks_subplot(ax_eflux)
    
    xmax = time[-2]
    
    ax_eflux.set_xlim(xmin=0, xmax=xmax)
    ax_eflux.set_ylim(ymin=0, ymax=30)
    
    ax_eflux.legend(loc=(0.845, 0.78))
    
    #wexb_max
    wexb, rad_coord, rad_boxsize, ddphi, dx, zonal_pot = zonalflow.get_shearingrate_radialcoordinate_radialboxsize_ddphi_dx_zonalpot(i)
    wexb_max = zonalflow.get_max_shearingrate(i, wexb, time, 1)
    
    ax_wexb_max.plot(time, wexb_max[1][:max_index], label=n)
    
    plot.ax_ticks_subplot(ax_wexb_max)
    
    ax_wexb_max.set_xlim(xmin=0, xmax=xmax)
    ax_wexb_max.set_ylim(ymin=0, ymax=0.30)
    
    #ax_wexb_max.legend(loc='lower right')

#plt.subplots_adjust(wspace=0, hspace=0)
plt.subplots_adjust(top=0.9, wspace=0.4, hspace=0.3)

plot.savefig_subplot(fig, ax_eflux   , picDir + '/S6_rlt6.0_boxsize1x1_Ns16_Nvpar48_Nmu9_dtim0.020-0.025_eflux_comparison.pdf'   , pad=0.02)
plot.savefig_subplot(fig, ax_wexb_max, picDir + '/S6_rlt6.0_boxsize1x1_Ns16_Nvpar48_Nmu9_dtim0.020-0.025_wexb_max_comparison.pdf', pad=0.02)

ax_eflux.text(-0.05, -0.09, r'\bf{(a)}', transform=ax_eflux.transAxes)
ax_wexb_max.text(-0.05, -0.09, r'\bf{(b)}', transform=ax_wexb_max.transAxes)

plt.savefig(picDir + '/S6_rlt6.0_boxsize1x1_Ns16_Nvpar48_Nmu9_dtim0.020-0.025_comparison.pdf', bbox_inches='tight')