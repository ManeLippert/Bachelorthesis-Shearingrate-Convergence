# MODULES ===========================================================================================================================================

# Import modules
import sys, os, h5py
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.transforms import Bbox

# PATH & IMPORT  ====================================================================================================================================

filepath = os.getcwd()
homepath = filepath.split('evaluation')[0]
sys.path.insert(1, homepath + 'python')

import zonalflow, h5tools, plot

picDir = homepath + 'pictures/Comparison/Boxsize/'
# Create target Directory if don't exist
if not os.path.exists(picDir):
    os.makedirs(picDir)

# linear growth rate
lineargrowth = pd.read_csv(homepath + 'data/linear_growthrate_ITG.dat', index_col=0)

lineargrowth_rlt = lineargrowth['gamma_N'][6.0]
lineargrowth_rlt_color = '#949494'

plot.parameters(True, 40, (24,8), 300, linewidth=2)

# FUNCTION ==========================================================================================================================================

def rotate(l, n):
    return np.concatenate((l[n:],l[:n]))

def box_plot(fig, time_label, plot_label,
             hdf5_file, xshift, interval, colors, 
             box_label, box_index, box_max, box_min, box_maximal,
             layer, layer_hspace = 0, legend_vpos = 0.5,
             WHITE = False, LEGEND = True, LEGEND_TOP = True,
             XLABEL = True, XLABEL_TOP = True,
             YLABEL = True,
             LABEL = False):
    
    if WHITE:
        plt.rcParams['xtick.color']='white'
        plt.rcParams['ytick.color']='white'
        plt.rcParams['axes.labelcolor']='white'
        plt.rcParams['axes.edgecolor']='white'
        plt.rcParams['lines.color']='white'
        plt.rcParams['text.color']='white'
    
    for b, n, c, f, s, t_start, t_end in zip(box_index, box_label, colors, hdf5_file, xshift, interval[0], interval[1]):
    
        ax = fig.add_axes([0, -(layer  - 1)*(1 + layer_hspace), b/box_maximal, 1])

        ax.set_facecolor('none')

        if b < box_max:        

            axR = ax.secondary_yaxis('right')
            axR.tick_params(direction = "out")
            axR.yaxis.set_ticklabels([])

            ax.set_xticklabels([])
            ax.set_yticklabels([])
            ax.set_xticks([])

        # Load shering rate and time
        eflux_data, time = zonalflow.get_eflux_time(f)
        wexb, rad_coord, rad_boxsize, ddphi, dx, zonal_pot = zonalflow.get_shearingrate_radialcoordinate_radialboxsize_ddphi_dx_zonalpot(f)

        start, end = zonalflow.get_index_from_value(time,t_start) , zonalflow.get_index_from_value(time,t_end)
        start_time, end_time = t_start, t_end
        # Shearing rate with mean over time
        wexb_rad_mean, wexb_rad_middle = zonalflow.get_mean_middle_shearingrate(start, end, wexb)
        # FT{shearing rate}
        #wexb_rad_mean_amp, wexb_rad_mean_amp_max = zonalflow.get_fft_mean_max_shearingrate_amplitude(wexb_rad_mean)

        lenght = len(rad_coord)

        wexb_rad_mean = rotate(wexb_rad_mean, s)

        label_time = r' $t_\mathrm{' + n + r'}^' + time_label + '\in$ [' + str(start_time) + r', ' + str(end_time) + r']' # + '\n' + n

        ax.plot(rad_coord, wexb_rad_mean, color=c, linewidth = 5, linestyle='-', label = label_time)

        if b == box_max:

            # linear growth rate
            ax.plot(rad_coord,  np.repeat(lineargrowth_rlt, len(rad_coord)), linewidth = 4, linestyle = 'dashed', color = lineargrowth_rlt_color)
            ax.plot(rad_coord, -np.repeat(lineargrowth_rlt, len(rad_coord)), linewidth = 4, linestyle = 'dashed', color = lineargrowth_rlt_color)

            ax.text(1.01, 0.74, r'\boldmath{$+\gamma$}', color = lineargrowth_rlt_color, transform=ax.transAxes)
            ax.text(1.01, 0.21, r'\boldmath{$-\gamma$}', color = lineargrowth_rlt_color, transform=ax.transAxes)
            
            if LABEL:
                ax.text(1.03, 0.9, plot_label, transform=ax.transAxes)
            
            ax.set_xticks(np.arange(0, rad_boxsize + 1, rad_boxsize/(2*box_max)))
            
            if XLABEL:
                ax_label = [str(round(x,1)) for x in np.arange(0, rad_boxsize + 1, rad_boxsize/(2*box_max))]
                ax_label[0] = '0'
                ax.set_xticklabels(ax_label)
                
                ax.set_xlabel(r'$\psi~[\rho]$')
                ax.xaxis.set_label_coords(0.5, -0.14, transform=ax.transAxes)
            else:
                ax.set_xticklabels([])
            
            axT = ax.secondary_xaxis('top')
            axT.set_xticks(np.arange(0, rad_boxsize + 1, rad_boxsize/(2*box_max)))
            
            if XLABEL_TOP:
                axT_label = [str(x) for x in np.arange(0, box_max + 0.5, box_max/(2*box_max))]
                axT_label[0] = '0'
                axT.set_xticklabels(axT_label)
                
                axT.set_xlabel(r'box size')
                axT.xaxis.set_label_coords(0.5, 1.14, transform=ax.transAxes)
            else:
                axT.set_xticklabels([])
            
            if YLABEL:
                ax.set_ylabel(r'$\omega_{\mathrm{E \times B}}~[\nu_{\mathrm{th}}/R]$')
                ax.yaxis.set_label_coords(-0.25/box_max, 0.5, transform=ax.transAxes)


        ax.set_xlim(xmin=0, xmax = rad_boxsize)
        ax.set_ylim(ymin=-0.5, ymax = 0.5)
    
    handles_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
    handles, labels = [sum(lol, []) for lol in zip(*handles_labels)]
    
    if LEGEND:
        if LEGEND_TOP:
            if XLABEL_TOP:  
                fig.legend(handles[::-1], labels[::-1],loc='upper center', bbox_to_anchor=(0.5, 1.45), frameon=False, ncol=len(box_index), handlelength=1)
            else:
                fig.legend(handles[::-1], labels[::-1],loc='upper center', bbox_to_anchor=(0.5, 1.2), frameon=False, ncol=len(box_index), handlelength=1)
        else:
            fig.legend(handles[::-1], labels[::-1],loc='center left', bbox_to_anchor=(1 + 0.15/box_max, legend_vpos), frameon=False, ncol=1, handlelength=1)
    
    if box_min > 1:
        
        while box_min > 1:
            box_min -= 1
            
            ax = fig.add_axes([0, -(layer - 1)*(1 + layer_hspace), box_min/box_maximal, 1])
            
            ax.set_facecolor('none')        

            axR = ax.secondary_yaxis('right')
            axR.tick_params(direction = "out")
            axR.yaxis.set_ticklabels([])

            ax.set_xticklabels([])
            ax.set_yticklabels([])
            ax.set_xticks([])

def box_twin_xspine(fig, hdf5_file, box_max, layer_max, linewidth = 2, layer_hspace = 0):
    
    ax = fig.add_axes([0, 0, 1, 1])
    ax.patch.set_visible(False)
    
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_xticks([])
    ax.set_yticks([])
    
    # Load shering rate and time
    wexb, rad_coord, rad_boxsize, ddphi, dx, zonal_pot = zonalflow.get_shearingrate_radialcoordinate_radialboxsize_ddphi_dx_zonalpot(hdf5_file)
    
    newax = ax.twiny()
                
    newax.set_xticks(np.arange(0, rad_boxsize + 1, rad_boxsize/(2*box_max)))
    newax_label = [str(round(x,1)) for x in np.arange(0, rad_boxsize + 1, rad_boxsize/(2*box_max))]
    newax_label[0] = '0'
    newax.set_xticklabels(newax_label)
    
    newax.set_xlabel(r'$\psi~[\rho]$')
    
    newax.tick_params(width=linewidth)
    
    newax.set_frame_on(True)
    newax.patch.set_visible(False)
    newax.xaxis.set_ticks_position('bottom')
    newax.xaxis.set_label_position('bottom')
    newax.spines['bottom'].set_position(('axes', -(layer_max - 1)*(1 + layer_hspace) - 0.1))
    newax.spines['bottom'].set_linewidth(linewidth)
    

# RADIAL ============================================================================================================================================

data = 'S6_rlt6.0'
path_rad = ['boxsize4x1/Ns16/Nvpar48/Nmu9', 'boxsize3x1/Ns16/Nvpar48/Nmu9', 'boxsize2x1/Ns16/Nvpar48/Nmu9', 'boxsize1x1/Ns16/Nvpar48/Nmu9']

filename_rad = [homepath + 'data/'+data+'/'+i+'/data.h5' for i in path_rad]
file_rad = [h5py.File(i,"r+") for i in filename_rad]

boxes_rad   = [4, 3, 2, 1]
shifts_rad  = [20, 2, -5, 60]
boxsize_rad = [r'4 \times 1\;\:\:', r'3 \times 1\;\:\:', r'2 \times 1\;\:\:', r'1 \times 1\;\:\:']
colors_rad  = ['#d55e00', '#029e73', '#de8f05', '#0173b2']

interval_rad = np.array([[26000, 43000, 15000, 2000],
                       [28000, 45000, 18000, 5000]])
                     
layer_max, box_max_rad, box_min_rad = 1, max(boxes_rad), min(boxes_rad)

time_label_rad = '{R}'

'''
fig_rad = plt.figure(1, figsize = (6*box_max_rad,6*layer_max))

box_plot(fig_rad, layer, time_label_rad, 
         file_rad, shifts_rad, interval_rad, colors_rad, 
         boxsize_rad, boxes_rad, box_max_rad, box_min_rad)

plt.savefig(picDir + '/S6_rlt6.0_boxsize1-2-3-4x1_Ns16_Nvpar48_Nmu9_wexb_comparison.pdf', bbox_inches='tight')

#plt.savefig(picDir + '/S6_rlt6.0_boxsize1-2-3-4x1_Ns16_Nvpar48_Nmu9_wexb_comparison.png', bbox_inches='tight', transparent = True)
#'''


# BINORMAL ==========================================================================================================================================

data = 'S6_rlt6.0'
path_rad_bi = ['boxsize3x3/Ns16/Nvpar48/Nmu9', 'boxsize2x2/Ns16/Nvpar48/Nmu9', 'boxsize1x1/Ns16/Nvpar48/Nmu9']

filename_rad_bi = [homepath + 'data/'+data+'/'+i+'/data.h5' for i in path_rad_bi]
file_rad_bi = [h5py.File(i,"r+") for i in filename_rad_bi]

plot.parameters(True, 40, (24,8), 300, linewidth=2)

boxes_rad_bi   = [3, 2, 1]
shifts_rad_bi  = [-5, -28, 60]
boxsize_rad_bi = [r'3 \times 3\;\:\:', r'2 \times 2\;\:\:', r'1 \times 1\;\:\:']
colors_rad_bi  = ['#57db5f', '#db57a4', '#0173b2']

interval_rad_bi = np.array([[2000, 2000, 2000],
                        [3000, 3000, 5000]])

layer_max, box_max_rad_bi, box_min_rad_bi = 1, max(boxes_rad_bi), min(boxes_rad_bi)

time_label_rad_bi = '{RB}'

'''
fig_rad_bi = plt.figure(1, figsize = (6*box_max_rad_bi,6*layer_max))

box_plot(fig_rad_bi, layer, time_label_rad_bi, 
         file_rad_bi, shifts_rad_bi, interval_rad_bi, colors_rad_bi, 
         boxsize_rad_bi, boxes_rad_bi, box_max_rad_bi, box_min_rad_bi, 
         LABEL_TOP=False)

plt.savefig(picDir + '/S6_rlt6.0_boxsize1x1-2x2-3x3_Ns16_Nvpar48_Nmu9_wexb_comparison.pdf', bbox_inches='tight')

#plt.savefig(picDir + '/S6_rlt6.0_boxsize1x1-2x2-3x3_Ns16_Nvpar48_Nmu9_wexb_comparison.png', bbox_inches='tight', transparent = True)
#'''


# BINORMAL SCAN =====================================================================================================================================

data = 'S6_rlt6.0'
#path_bi = ['boxsize3x5/Ns16/Nvpar48/Nmu9', 'boxsize3x3/Ns16/Nvpar48/Nmu9', 'boxsize3x2.5/Ns16/Nvpar48/Nmu9',
#          'boxsize3x1.5/Ns16/Nvpar48/Nmu9', 'boxsize3x1/Ns16/Nvpar48/Nmu9']
path_bi = ['boxsize3x5/Ns16/Nvpar48/Nmu9', 'boxsize3x2.5/Ns16/Nvpar48/Nmu9', 'boxsize3x1.5/Ns16/Nvpar48/Nmu9/Run2']

filename_bi = [homepath + 'data/'+data+'/'+i+'/data.h5' for i in path_bi]
file_bi = [h5py.File(i,"r+") for i in filename_bi]

plot.parameters(True, 40, (24,8), 300, linewidth=2)

#boxes_bi   = [3, 3, 3, 3, 3]
boxes_bi   = [3, 3, 3]
#shifts_bi  = [-10, -5, -10, -10, 2]
shifts_bi  = [-2, -2, -2]
#boxsize_bi = [r'3 \times 5\;\:\:', r'3 \times 3\;\:\:', r'3 \times 2.5', r'3 \times 1.5', r'3 \times 1\;\:\:']
boxsize_bi = [r'3 \times 5\;\:\:', r'3 \times 2.5', r'3 \times 1.5']
#colors_bi  = ['#cc78bc', '#029e73', '#d55e00', '#de8f05', '#0173b2']
colors_bi  = ['#cc78bc', '#56b4e9', '#ca9161']


#interval_bi = np.array([[2000, 2000, 2000, 11000, 43000],
#                       [ 4000, 3000, 3000, 12000, 45000]])

interval_bi = np.array([[2000, 2000, 11000],
                       [ 4000, 3000, 12000]])

layer_max, box_max_bi, box_min_bi = 1, max(boxes_bi), min(boxes_bi)

time_label_bi = '{B}'

'''
fig_bi = plt.figure(1, figsize = (6*box_max_bi,6*layer_max))

box_plot(fig_bi, layer, time_label_bi, 
         file_bi, shifts_bi, interval_bi, colors_bi, 
         boxsize_bi, boxes_bi, box_max_bi, box_min_bi, 
         LABEL_TOP=False)

plt.savefig(picDir + '/S6_rlt6.0_boxsize3x1-1.5-2.5-3-5_Ns16_Nvpar48_Nmu9_wexb_comparison.pdf', bbox_inches='tight')

#plt.savefig(picDir + '/S6_rlt6.0_boxsize1x1-2x2-3x3_Ns16_Nvpar48_Nmu9_wexb_comparison.png', bbox_inches='tight', transparent = True)
#'''



# ALL BOXES =========================================================================================================================================

file = [file_rad, file_rad_bi, file_bi]

plot.parameters(True, 40, (24,8), 300, linewidth=2)

boxes   = [boxes_rad, boxes_rad_bi, boxes_bi]
shifts  = [shifts_rad, shifts_rad_bi, shifts_bi]
boxsize = [boxsize_rad, boxsize_rad_bi, boxsize_bi]
colors  = [colors_rad, colors_rad_bi, colors_bi]

interval = [interval_rad, interval_rad_bi, interval_bi]

time_label = [time_label_rad, time_label_rad_bi, time_label_bi]

box_max, box_min = [box_max_rad, box_max_rad_bi, box_max_bi], [box_min_rad, box_min_rad_bi, box_min_bi]

xlabel = [False, False, False]
ylabel = [False, True, False]
xlabel_top = [True, False, False]

label = [ r'\bf{(a)}',  r'\bf{(b)}',  r'\bf{(c)}']


i = 0
layer, layer_hspace = 1, 0.1
layer_max, box_maximal, box_minimal = len(boxes), max(box_max), min(box_min)

for f in file:
    
    fig = plt.figure(1, figsize = (6*box_maximal,6*layer))
    
    box_plot(fig, time_label[i], label[i], f, shifts[i], interval[i], colors[i], boxsize[i], boxes[i], box_max[i], box_min[i], box_maximal, layer, layer_hspace,
             LEGEND_TOP=False, XLABEL=xlabel[i], XLABEL_TOP=xlabel_top[i], YLABEL=ylabel[i], LEGEND=False, LABEL=True)
    
    layer += 1
    i += 1
    
handles_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
handles, labels = [sum(lol, []) for lol in zip(*handles_labels)]
l = plt.Line2D([0],[0],color="w",alpha=0)

handles = handles[7:len(labels)] + [l] + handles[4:7] + [l] + handles[0:4]
labels = labels[7:len(labels)] + [''] + labels[4:7] + [''] + labels[0:4]

fig.legend(handles[::-1], labels[::-1],loc='center left', bbox_to_anchor=(3/box_maximal + 0.3/box_maximal, -1.1), frameon=False, ncol=1, handlelength=1)

box_twin_xspine(fig, file[0][0], box_maximal, layer_max, layer_hspace=layer_hspace)

plt.savefig(picDir + '/S6_rlt6.0_boxsize1-2-3-4x1-1.5-2-2.5-3-5_Ns16_Nvpar48_Nmu9_wexb_comparison.pdf', bbox_inches='tight')

#'''