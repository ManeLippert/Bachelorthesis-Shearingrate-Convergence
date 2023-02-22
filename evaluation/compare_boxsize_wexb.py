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

# RADIAL ============================================================================================================================================


# File import and Create picture folder
data = 'S6_rlt6.0'
path = ['boxsize4x1/Ns16/Nvpar48/Nmu9', 'boxsize3x1/Ns16/Nvpar48/Nmu9',
        'boxsize2x1/Ns16/Nvpar48/Nmu9', 'boxsize1x1/Ns16/Nvpar48/Nmu9']

filename = [homepath + 'data/'+data+'/'+i+'/data.h5' for i in path]
f = [h5py.File(i,"r+") for i in filename]

picDir = homepath + 'pictures/Comparison/Boxsize/'
# Create target Directory if don't exist
if not os.path.exists(picDir):
    os.makedirs(picDir)
    
# linear growth rate
lineargrowth = pd.read_csv(homepath + 'data/linear_growthrate_ITG.dat', index_col=0)

lineargrowth_rlt = lineargrowth['gamma_N'][float(data.split('rlt')[1])]
lineargrowth_rlt_color = 'grey'

plot.parameters(True, 40, (24,8), 300, linewidth=2)

'''
# parameter vor white plots
plt.rcParams['xtick.color']='white'
plt.rcParams['ytick.color']='white'
plt.rcParams['axes.labelcolor']='white'
plt.rcParams['axes.edgecolor']='white'
plt.rcParams['lines.color']='white'
plt.rcParams['text.color']='white'
'''

# Compare shearing rate in radial domain
fig = plt.figure(1, figsize = (24,6))

ax4 = fig.add_axes([0, 0, 4/4, 1])
ax3 = fig.add_axes([0, 0, 3/4, 1])
ax2 = fig.add_axes([0, 0, 2/4, 1])
ax1 = fig.add_axes([0, 0, 1/4, 1])

boxsize = [r'4 \times 1', r'3 \times 1', r'2 \times 1', r'1 \times 1']
colors = ['#0173b2', '#de8f05', '#029e73', '#d55e00']
#colors = ['#9a00ff','#ffad00', '#0cbf0f', '#fc0307']


axes = [ax4, ax3, ax2, ax1]
interval = np.array([[26000, 43000, 15000, 2000],
                     [28000, 45000, 18000, 5000]])

def rotate(l, n):
    return np.concatenate((l[n:],l[:n]))

i = 4

for ax, b, file in zip(axes, boxsize, f):

    if i < 4:
        ax.set_facecolor('none')        
        
        axR = ax.secondary_yaxis('right')
        axR.tick_params(direction = "out")
        axR.yaxis.set_ticklabels([])
        
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_xticks([])
    
    # Load shering rate and time
    eflux_data, time = zonalflow.get_eflux_time(file)
    wexb, rad_coord, rad_boxsize, ddphi, dx, zonal_pot = zonalflow.get_shearingrate_radialcoordinate_radialboxsize_ddphi_dx_zonalpot(file)
    
    start, end = zonalflow.get_index_from_value(time,interval[0][4-i]) , zonalflow.get_index_from_value(time,interval[1][4-i])
    start_time, end_time = interval[0][4-i], interval[1][4-i]
    
    # Shearing rate with mean over time
    wexb_rad_mean, wexb_rad_middle = zonalflow.get_mean_middle_shearingrate(start, end, wexb)
    # FT{shearing rate}
    #wexb_rad_mean_amp, wexb_rad_mean_amp_max = zonalflow.get_fft_mean_max_shearingrate_amplitude(wexb_rad_mean)

    #plot.mean_shearingrate_radialcoordinate_subplot(rad_coord, rad_boxsize, wexb_rad_mean, wexb_rad_middle, wexb_rad_mean_amp_max, 
    #                                                ax, x, y, xdim, ydim, start_time, end_time)
    
    
    
    lenght = len(rad_coord)
    
    if i == 4:
        wexb_rad_mean = rotate(wexb_rad_mean, 20)
        label_time = r' $t_\mathrm{' + b + r'} \in$ [' + str(start_time) + r', ' + str(end_time) + r']' # + '\n' + b
        
        #ax.plot(rad_coord[int(0/4*lenght):int(1/4*lenght)+1], wexb_rad_mean[int(0/4*lenght):int(1/4*lenght)+1], color=colors[i-1], linewidth = 4, linestyle='-')
        #ax.plot(rad_coord[int(1/4*lenght):int(2/4*lenght)+1], wexb_rad_mean[int(1/4*lenght):int(2/4*lenght)+1], color=colors[i-1], linewidth = 4, linestyle='-')
        #ax.plot(rad_coord[int(2/4*lenght):int(3/4*lenght)+1], wexb_rad_mean[int(2/4*lenght):int(3/4*lenght)+1], color=colors[i-1], linewidth = 4, linestyle='-')
        #ax.plot(rad_coord[int(3/4*lenght):int(4/4*lenght)+1], wexb_rad_mean[int(3/4*lenght):int(4/4*lenght)+1], color=colors[i-1], linewidth = 4, linestyle='-', label = label_time)
        
        ax.plot(rad_coord, wexb_rad_mean, color=colors[i-1], linewidth = 4, linestyle='-', label = label_time)
        
        # linear growth rate
        ax.plot(rad_coord,  np.repeat(lineargrowth_rlt, len(rad_coord)), linewidth = 4, linestyle = 'dashed', color = lineargrowth_rlt_color)
        ax.plot(rad_coord, -np.repeat(lineargrowth_rlt, len(rad_coord)), linewidth = 4, linestyle = 'dashed', color = lineargrowth_rlt_color)
        
        ax.text(1.01, 0.74, r'\boldmath{$+\gamma$}', color = lineargrowth_rlt_color, transform=ax.transAxes)
        ax.text(1.01, 0.21, r'\boldmath{$-\gamma$}', color = lineargrowth_rlt_color, transform=ax.transAxes)
        
        ax.legend(loc='upper center', bbox_to_anchor=(7/8, 1.45), frameon=False, handlelength=1)
        
        ax.set_xticks([0, 7/8*rad_boxsize, 8/8*rad_boxsize])
        ax.set_xticklabels([str(0), str(round(7/8*rad_boxsize, 1)), str(round(8/8*rad_boxsize, 1))])
        
        axT = ax.secondary_xaxis('top')
        axT.set_xticks([0, 7/8*rad_boxsize, rad_boxsize])
        axT.set_xticklabels(['0','3.5', '4.0'])
        
    elif i == 3:
        wexb_rad_mean = rotate(wexb_rad_mean, 2)
        label_time = r' $t_\mathrm{' + b + r'} \in$ [' + str(start_time) + r', ' + str(end_time) + r']' # + '\n' + b
        
        #ax.plot(rad_coord[int(0/3*lenght):int(1/3*lenght)+1], wexb_rad_mean[int(0/3*lenght):int(1/3*lenght)+1], color=colors[i-1], linewidth = 4, linestyle='-')
        #ax.plot(rad_coord[int(1/3*lenght):int(2/3*lenght)+1], wexb_rad_mean[int(1/3*lenght):int(2/3*lenght)+1], color=colors[i-1], linewidth = 4, linestyle='-')
        #ax.plot(rad_coord[int(2/3*lenght):int(3/3*lenght)+1], wexb_rad_mean[int(2/3*lenght):int(3/3*lenght)+1], color=colors[i-1], linewidth = 4, linestyle='-', label = label_time)
        
        ax.plot(rad_coord, wexb_rad_mean, color=colors[i-1], linewidth = 4, linestyle='-', label = label_time)

        ax.legend(loc='upper center', bbox_to_anchor=(5/6, 1.45), frameon=False, handlelength=1)
        
        ax.set_xticks([5/6*rad_boxsize, rad_boxsize])
        ax.set_xticklabels([str(round(5/6*rad_boxsize, 1)), str(round(rad_boxsize, 1))])
        
        axT = ax.secondary_xaxis('top')
        axT.set_xticks([5/6*rad_boxsize, rad_boxsize])
        axT.set_xticklabels(['2.5', '3.0'])
        
    elif i == 2:
        wexb_rad_mean = rotate(wexb_rad_mean, -5)
        label_time = r' $t_\mathrm{' + b + r'} \in$ [' + str(start_time) + r', ' + str(end_time) + r']' # + '\n' + b

        #ax.plot(rad_coord[int(0/2*lenght):int(1/2*lenght)+1], wexb_rad_mean[int(0/2*lenght):int(1/2*lenght)+1], color=colors[i-1], linewidth = 4, linestyle='-')
        #ax.plot(rad_coord[int(1/2*lenght):int(2/2*lenght)+1], wexb_rad_mean[int(1/2*lenght):int(2/2*lenght)+1], color=colors[i-1], linewidth = 4, linestyle='-', label = label_time)
        
        ax.plot(rad_coord, wexb_rad_mean, color=colors[i-1], linewidth = 4, linestyle='-', label = label_time)
        
        ax.legend(loc='upper center', bbox_to_anchor=(3/4, 1.45), frameon=False, handlelength=1)
        
        ax.set_xticks([3/4*rad_boxsize, rad_boxsize])
        ax.set_xticklabels([str(round(3/4*rad_boxsize, 1)), str(round(rad_boxsize, 1))])
        
        axT = ax.secondary_xaxis('top')
        axT.set_xticks([3/4*rad_boxsize, rad_boxsize])
        axT.set_xticklabels(['1.5', '2.0'])
        
    elif i == 1:
        wexb_rad_mean = rotate(wexb_rad_mean, 60)
        label_time = r' $t_\mathrm{' + b + r'} \in$ [' + str(start_time) + r', ' + str(end_time) + r']' # + '\n' + b
    
        #ax.plot(rad_coord[int(0/1*lenght):int(1/1*lenght)+1], wexb_rad_mean[int(0/1*lenght):int(1/1*lenght)+1], color=colors[i-1], linewidth = 4, label = label_time)
        
        ax.plot(rad_coord, wexb_rad_mean, color=colors[i-1], linewidth = 4, linestyle='-', label = label_time)
        
        ax.legend(loc='upper center', bbox_to_anchor=(1/2, 1.45), frameon=False, handlelength=1)
        
        ax.set_xticks([1/2*rad_boxsize, rad_boxsize])
        ax.set_xticklabels([str(round(1/2*rad_boxsize, 1)), str(round(rad_boxsize, 1))])
        
        axT = ax.secondary_xaxis('top')
        axT.set_xticks([1/2*rad_boxsize, rad_boxsize])
        axT.set_xticklabels(['0.5', '1.0'])
        
        ax.set_xlabel(r'$\psi~[\rho]$')
        ax.xaxis.set_label_coords(2.0, -0.14, transform=ax.transAxes)
        
        ax.set_ylabel(r'$\omega_{\mathrm{E \times B}}~[\nu_{\mathrm{th}}/R]$')
        ax.yaxis.set_label_coords(-0.25, 0.5, transform=ax.transAxes)

        axT.set_xlabel(r'box size')
        axT.xaxis.set_label_coords(2.0, 1.14, transform=ax.transAxes)
        
        
    
    ax.set_xlim(xmin=0, xmax = rad_boxsize)
    ax.set_ylim(ymin=-0.4, ymax = 0.4)
    
    i -= 1

plt.savefig(picDir + '/S6_rlt6.0_boxsize1-2-3-4x1_Ns16_Nvpar48_Nmu9_wexb_comparison.pdf', bbox_inches='tight')

#plt.savefig(picDir + '/S6_rlt6.0_boxsize1-2-3-4x1_Ns16_Nvpar48_Nmu9_wexb_comparison.png', bbox_inches='tight', transparent = True)

# BINORMAL ==========================================================================================================================================

# File import and Create picture folder
data = 'S6_rlt6.0'
path = ['boxsize3x3/Ns16/Nvpar48/Nmu9', 'boxsize2x2/Ns16/Nvpar48/Nmu9']

filename = [homepath + 'data/'+data+'/'+i+'/data.h5' for i in path]
f = [h5py.File(i,"r+") for i in filename]

picDir = homepath + 'pictures/Comparison/Boxsize/'
# Create target Directory if don't exist
if not os.path.exists(picDir):
    os.makedirs(picDir)
    
# linear growth rate
lineargrowth = pd.read_csv(homepath + 'data/linear_growthrate_ITG.dat', index_col=0)

lineargrowth_rlt = lineargrowth['gamma_N'][float(data.split('rlt')[1])]
lineargrowth_rlt_color = 'grey'

plot.parameters(True, 40, (24,8), 300, linewidth=2)

# Compare shearing rate in radial domain
fig = plt.figure(1, figsize = (24,6))

ax3 = fig.add_axes([0, 0, 3/3, 1])
ax2 = fig.add_axes([0, 0, 2/3, 1])

boxsize = [r'3 \times 3', r'2 \times 2']
colors = ['#de8f05', '#d55e00']

axes = [ax3, ax2]
interval = np.array([[2000,2000],
                     [3000,3000]])

def rotate(l, n):
    return np.concatenate((l[n:],l[:n]))

i = 2

for ax, b, file in zip(axes, boxsize, f):

    if i < 2:
        ax.set_facecolor('none')        
        
        axR = ax.secondary_yaxis('right')
        axR.tick_params(direction = "out")
        axR.yaxis.set_ticklabels([])
        
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_xticks([])
    
    # Load shering rate and time
    eflux_data, time = zonalflow.get_eflux_time(file)
    wexb, rad_coord, rad_boxsize, ddphi, dx, zonal_pot = zonalflow.get_shearingrate_radialcoordinate_radialboxsize_ddphi_dx_zonalpot(file)
    
    start, end = zonalflow.get_index_from_value(time,interval[0][2-i]) , zonalflow.get_index_from_value(time,interval[1][2-i])
    start_time, end_time = interval[0][2-i], interval[1][2-i]
    
    # Shearing rate with mean over time
    wexb_rad_mean, wexb_rad_middle = zonalflow.get_mean_middle_shearingrate(start, end, wexb)
    # FT{shearing rate}
    #wexb_rad_mean_amp, wexb_rad_mean_amp_max = zonalflow.get_fft_mean_max_shearingrate_amplitude(wexb_rad_mean)

    #plot.mean_shearingrate_radialcoordinate_subplot(rad_coord, rad_boxsize, wexb_rad_mean, wexb_rad_middle, wexb_rad_mean_amp_max, 
    #                                                ax, x, y, xdim, ydim, start_time, end_time)
    
    
    
    lenght = len(rad_coord)
    
    if i == 2:
        #wexb_rad_mean = rotate(wexb_rad_mean, 20)
        label_time = r' $t_\mathrm{' + b + r'} \in$ [' + str(start_time) + r', ' + str(end_time) + r']' # + '\n' + b
        
        #ax.plot(rad_coord[int(0/4*lenght):int(1/4*lenght)+1], wexb_rad_mean[int(0/4*lenght):int(1/4*lenght)+1], color=colors[i-1], linewidth = 4, linestyle='-')
        #ax.plot(rad_coord[int(1/4*lenght):int(2/4*lenght)+1], wexb_rad_mean[int(1/4*lenght):int(2/4*lenght)+1], color=colors[i-1], linewidth = 4, linestyle='-')
        #ax.plot(rad_coord[int(2/4*lenght):int(3/4*lenght)+1], wexb_rad_mean[int(2/4*lenght):int(3/4*lenght)+1], color=colors[i-1], linewidth = 4, linestyle='-')
        #ax.plot(rad_coord[int(3/4*lenght):int(4/4*lenght)+1], wexb_rad_mean[int(3/4*lenght):int(4/4*lenght)+1], color=colors[i-1], linewidth = 4, linestyle='-', label = label_time)
        
        ax.plot(rad_coord, wexb_rad_mean, color=colors[i-1], linewidth = 4, linestyle='-', label = label_time)
        
        # linear growth rate
        ax.plot(rad_coord,  np.repeat(lineargrowth_rlt, len(rad_coord)), linewidth = 4, linestyle = 'dashed', color = lineargrowth_rlt_color)
        ax.plot(rad_coord, -np.repeat(lineargrowth_rlt, len(rad_coord)), linewidth = 4, linestyle = 'dashed', color = lineargrowth_rlt_color)
        
        ax.text(1.01, 0.74, r'\boldmath{$+\gamma$}', color = lineargrowth_rlt_color, transform=ax.transAxes)
        ax.text(1.01, 0.21, r'\boldmath{$-\gamma$}', color = lineargrowth_rlt_color, transform=ax.transAxes)
        
        ax.legend(loc='upper center', bbox_to_anchor=(2/3, 1.45), frameon=False, handlelength=1)
        
        ax.set_xticks([0, 5/6*rad_boxsize, 6/6*rad_boxsize])
        ax.set_xticklabels([str(0), str(round(5/6*rad_boxsize, 1)), str(round(5/6*rad_boxsize, 1))])
        
        axT = ax.secondary_xaxis('top')
        axT.set_xticks([0, 5/6*rad_boxsize, rad_boxsize])
        axT.set_xticklabels(['0','2.5', '3.0'])
        
    elif i == 1:
        #wexb_rad_mean = rotate(wexb_rad_mean, 60)
        label_time = r' $t_\mathrm{' + b + r'} \in$ [' + str(start_time) + r', ' + str(end_time) + r']' # + '\n' + b
    
        #ax.plot(rad_coord[int(0/1*lenght):int(1/1*lenght)+1], wexb_rad_mean[int(0/1*lenght):int(1/1*lenght)+1], color=colors[i-1], linewidth = 4, label = label_time)
        
        ax.plot(rad_coord, wexb_rad_mean, color=colors[i-1], linewidth = 4, linestyle='-', label = label_time)
        
        ax.legend(loc='upper center', bbox_to_anchor=(1/2, 1.45), frameon=False, handlelength=1)
        
        ax.set_xticks([1/4*rad_boxsize, 2/4*rad_boxsize, 3/4*rad_boxsize, rad_boxsize])
        ax.set_xticklabels([str(round(1/4*rad_boxsize, 1)), str(round(2/4*rad_boxsize, 1)), str(round(3/4*rad_boxsize, 1)), str(round(rad_boxsize, 1))])
        
        axT = ax.secondary_xaxis('top')
        axT.set_xticks([1/4*rad_boxsize, 2/4*rad_boxsize, 3/4*rad_boxsize, rad_boxsize])
        axT.set_xticklabels(['0.5', '1.0', '1.5', '2.0'])
        
        ax.set_xlabel(r'$\psi~[\rho]$')
        ax.xaxis.set_label_coords(3/4, -0.14, transform=ax.transAxes)
        
        ax.set_ylabel(r'$\omega_{\mathrm{E \times B}}~[\nu_{\mathrm{th}}/R]$')
        ax.yaxis.set_label_coords(-0.1, 0.5, transform=ax.transAxes)

        axT.set_xlabel(r'box size')
        axT.xaxis.set_label_coords(3/4, 1.14, transform=ax.transAxes)
        
        
    
    ax.set_xlim(xmin=0, xmax = rad_boxsize)
    ax.set_ylim(ymin=-0.4, ymax = 0.4)
    
    i -= 1

plt.savefig(picDir + '/S6_rlt6.0_boxsize1x1-2x2-3x3_Ns16_Nvpar48_Nmu9_wexb_comparison.pdf', bbox_inches='tight')

#plt.savefig(picDir + '/S6_rlt6.0_boxsize1-2-3-4x1_Ns16_Nvpar48_Nmu9_wexb_comparison.png', bbox_inches='tight', transparent = True)

# BINORMAL SCAN =====================================================================================================================================

