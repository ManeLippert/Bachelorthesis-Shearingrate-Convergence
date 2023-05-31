# MODULES ===========================================================================================================================================

import sys, os, h5py, ast

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.transforms import Bbox

# PATH ==============================================================================================================================================

filepath = os.getcwd()
homepath = filepath.split('evaluation')[0]
sys.path.insert(1, homepath + 'python')

import zonalflow, h5tools, plot

print('\n<------------------------- Evaluation ------------------------->\n')

# PARAMETER =========================================================================================================================================

plot.parameters(40, (24,8), 300, linewidth = 4, tickwidth = 3, legendontop = True)

boxsize = '3x3'

ALL, EVOLUTION, SELECTION, PROFILE = False, False, True, True

#print('Input:', boxsize, ALL, EVOLUTION, SELECTION, '\n')

# DATA IMPORT =======================================================================================================================================

datainfopath = homepath + '/data/data.csv'
datainfo = zonalflow.get_data_info(datainfopath, boxsize, finit = "noise")

print(datainfo)

filepath = homepath + datainfo['path'].values[0]

if EVOLUTION:
    try:
        interval_evo = np.array([ast.literal_eval(datainfo['evo_start'].values[0]),
                                 ast.literal_eval(datainfo['evo_end'].values[0])])
    except ValueError:
        print('Evolution interval has to be defined in "data.csv"')
        EVOLUTION = False

if SELECTION or PROFILE:
    try:
        interval_sel = np.array([ast.literal_eval(datainfo['sel_start'].values[0]),
                                 ast.literal_eval(datainfo['sel_end'].values[0])])
    except ValueError:
        print('Selection interval has to be defined in "data.csv"')
        SELECTION, PROFILE = False, False

if ALL:
    stepsize_all = datainfo['stepsize_all'].values[0]
    

try:
    error_index = ast.literal_eval(datainfo['error_index'].values[0])
    
    try:
        if np.isnan(error_index):
            ERROR = False
        else:
            ERROR = True
            error_start, error_end = None, error_index
    except ValueError:
        ERROR = True
        error_index = np.array(error_index)
        error_start, error_end = error_index[0], error_index[1]
    
except ValueError:
    ERROR = False
    
if ERROR:
    print("Error Index : ", error_index)

# DATA ==============================================================================================================================================

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
try:
    lineargrowth_rlt = lineargrowth['gamma_N'][float(data.split('rlt')[1])]
except KeyError:
    lineargrowth_rlt = 0.21
lineargrowth_rlt_color = 'grey'

# HEAT FLUX =========================================================================================================================================

eflux_data, time = zonalflow.get_eflux_time(f)

if ERROR:
    if error_start == None:
        eflux_data, time = eflux_data[:error_end], time[:error_end]
    else:
        eflux_data, time = np.concatenate((eflux_data[:error_start], eflux_data[error_end:])), np.concatenate((time[:error_start], time[error_end:]))

plot.eflux_time(time, eflux_data, (24,8), xlim = (0, None), ylim = (0, 25))
plt.savefig(picpath+data+'_'+resolution+'_eflux.pdf', bbox_inches='tight')

# Omega_ExB =========================================================================================================================================

wexb, rad_coord, rad_boxsize, ddphi, dx, zonal_pot = zonalflow.get_shearingrate_radialcoordinate_radialboxsize_ddphi_dx_zonalpot(f)

print('rad_boxsize :', rad_boxsize, '; stepsize :',dx)

print(rad_coord)

# FOURIER ===========================================================================================================================================

wexb_max = zonalflow.get_max_shearingrate(f, wexb, time, 5)

if ERROR:
    wexb_max_list = []
    if error_start == None:
        for i in wexb_max:
            wexb_max_list.append(i[:error_end])
    else:
        for i in wexb_max:
            wexb_max_list.append(np.concatenate((i[:error_start], i[error_end:])))
        
    wexb_max = np.array(wexb_max_list)
    
#print(len(wexb_max))

plot.max_shearingrate_time(time, wexb_max, [1, 2, 3, 4], (24,8))
plt.savefig(picpath+data+'_'+resolution+'_wexb_max.pdf', bbox_inches='tight')


# ALL SHEARING RATE =================================================================================================================================

if ALL:
    
    plot.all_shearingrate_radialcoordinate(rad_coord, wexb, (12,8), stepsize_all)
    plt.savefig(picpath+data+'_'+resolution+'_wexb_all.pdf', bbox_inches='tight')

# TIME EVOLUTION ====================================================================================================================================

if EVOLUTION:
    
    try:
        # Dimension Subplot
        for i in [3,2,1]:
            if interval_evo.shape[1] % i == 0:
                xdim, ydim = i, int(interval_evo.shape[1]/i)
                break

        grid_x, grid_y = np.arange(xdim), np.arange(ydim)

        fig, ax = plt.subplots(ydim, xdim,figsize=(12*xdim, 8*ydim), sharey=True, sharex=True, squeeze=True)

        i = 0

        for y in grid_y:
            for x in grid_x:
                start, end = zonalflow.get_index_from_value(time,interval_evo[0][i]) , zonalflow.get_index_from_value(time,interval_evo[1][i])
                start_time, end_time = interval_evo[0][i], interval_evo[1][i]

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
    
    except TypeError:
        print('Evolution interval does not match time series')

# SELECTION =========================================================================================================================================

if SELECTION:
    
    try:
        # Plot parameter
        plot.parameters(50, (18,8), 300)

        fig, ax = plt.subplots(1, 1)

        colors = ['#029e73', '#d55e00']

        i = 0

        while i < len(interval_sel[0]):
            start_time, end_time = interval_sel[0][i], interval_sel[1][i]
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

        ax.legend(frameon=False, ncol = 4, handlelength=1)

        ax.set_xlim(0, rad_boxsize)
        ax.set_ylim(-0.4, 0.4)

        ax.set_yticks(np.arange(-0.4, 0.5, 0.2))

        ax.set_xlabel(r'$x~[\rho]$')
        ax.set_ylabel(r'$\langle\omega_{\mathrm{E \times B}}\rangle~[\nu_{\mathrm{th}}/R]$')

        plt.savefig(picpath+data+'_'+resolution+'_wexb_selection.pdf', bbox_inches='tight')
    
    except TypeError:
        print('Selection interval does not match time series')
        
# PROFILE ===========================================================================================================================================

if PROFILE:
    
    try:
        # Plot parameter
        plot.parameters(40, (18,16), 300)

        fig, ax = plt.subplots(2, 1, sharex=True)

        colors = ['#0173b2', '#de8f05', '#029e73','#a11a5b']

        i = 0

        while i < len(interval_sel[0]):
            
            start_time, end_time = interval_sel[0][i], interval_sel[1][i]
            label_time = r' $t \in$ [' + str(start_time) + r', ' + str(end_time) + r']'
            start, end = zonalflow.get_index_from_value(time, start_time) , zonalflow.get_index_from_value(time, end_time)

            wexb_rad_mean, wexb_rad_middle = zonalflow.get_mean_middle_shearingrate(start, end, wexb)
            dr_dens_mean, rad_coord = zonalflow.get_radial_density_profile(f, start, end)
            dr_ene_mean, rad_coord = zonalflow.get_radial_energy_profile(f, start, end)
            dr_zonal_pot_mean, rad_coord = zonalflow.get_radial_zonal_potential_profile(f, start, end)
            
            ax[0].plot(rad_coord, wexb_rad_mean, linewidth=4, color=colors[2])
            
            # linear growth rate
            ax[0].plot(rad_coord,  np.repeat(lineargrowth_rlt, len(rad_coord)), linewidth = 4, linestyle = 'dashed', color = lineargrowth_rlt_color)
            ax[0].plot(rad_coord, -np.repeat(lineargrowth_rlt, len(rad_coord)), linewidth = 4, linestyle = 'dashed', color = lineargrowth_rlt_color)

            ax[0].text(1.01, 0.735, r'\boldmath{$+\gamma$}', color = lineargrowth_rlt_color, transform=ax[0].transAxes)
            ax[0].text(1.01, 0.205, r'\boldmath{$-\gamma$}', color = lineargrowth_rlt_color, transform=ax[0].transAxes)
            
            ax[0].set_xlim(0, rad_boxsize)
            ax[0].set_ylim(-0.4, 0.4)

            ax[0].set_yticks(np.arange(-0.4, 0.5, 0.2))

            #ax[0].set_xlabel(r'$x~[\rho]$')
            ax[0].set_ylabel(r'$\langle\omega_{\mathrm{E \times B}}\rangle~[\nu_{\mathrm{th}}/R]$')
            

            #ax[1].plot(rad_coord, dr_zonal_pot_mean, linewidth=4, color = colors[0])
            
            #ax[1].set_xlim(0, rad_boxsize)
            #ax[1].set_ylim(-0.2, 0.2)
            #ax[1].set_xlabel(r'$x~[\rho]$')
            #ax[1].set_ylabel(r'$-\nabla_{\!\!x} \phi~[n_0T_0/R_0]$')
                        
            ax[1].plot(rad_coord, dr_ene_mean[0], linewidth=4, label = r"$E_{\perp, \mathrm{i}}$", color = colors[1])
            ax[1].plot(rad_coord, dr_ene_mean[1], linewidth=4, label = r"$E_{\parallel, \mathrm{i}}$", color = colors[3])

            ax[1].set_xlim(0, rad_boxsize)
            ax[1].set_ylim(-1.5, 1.5)
            ax[1].set_yticks(np.arange(-1.5, 1.6, 0.5))
            ax[1].set_xlabel(r'$x~[\rho]$')
            ax[1].set_ylabel(r'$-\nabla_{\!\!x} E~[n_0T_0/R_0]$')
            
            ax_right = ax[1].twinx()
            
            ax_right.plot(rad_coord, dr_zonal_pot_mean, linewidth=4, color = colors[0])
            ax_right.set_ylabel(r'$e n_0 \nabla_{\!\!x} \phi~[n_0T_0/R_0]$', color = colors[0])
            ax_right.set_ylim(-8, 8)
            ax_right.set_yticks(np.arange(-8, 9, 2))
            ax_right.tick_params(axis='y', colors = colors[0])
            
            i += 1
        
        fig.suptitle(label_time)
        
        handles_labels = [axis.get_legend_handles_labels() for axis in ax]
        handles, labels = [sum(lol, []) for lol in zip(*handles_labels)]
        
        fig.legend(handles, labels, ncol=3, bbox_to_anchor=(0.5, 0.88))
        
        plt.subplots_adjust(hspace=0.1)
        
        plt.savefig(picpath+data+'_'+resolution+'_wexb_dens_ene_selection.pdf', bbox_inches='tight')
    
    except TypeError:
        print('Selection interval does not match time series')
        

# SAVE DATA =========================================================================================================================================

data_eval = [dx, ddphi, zonal_pot, wexb, wexb_max]
data_group = ['derivative_stepsize', 'second_derivative_phi', 'zonalflow_potential', 'shearing_rate', 'shearing_rate_maximum']
    
groupname = ['evaluation' + '/' + x for x in data_group]
h5tools.hdf5_write_data(f, data_eval, groupname)

if PROFILE:
    f = h5py.File(filename,"r+")
    
    data_eval = [dr_dens_mean, dr_ene_mean[0], dr_ene_mean[1], dr_zonal_pot_mean, dr_zonal_pot_mean]
    data_group = ['derivative_dens', 'derivative_energy_perp', 'derivative_energy_par', 'derivative_zonalflow_potential', 'derivative_zonalflow_potential']

    groupname = ['evaluation' + '/' + x for x in data_group]
    h5tools.hdf5_write_data(f, data_eval, groupname)