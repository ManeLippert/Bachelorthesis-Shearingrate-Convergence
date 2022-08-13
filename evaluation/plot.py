from cProfile import label
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.transforms import Bbox

# Plot parameters
def parameters(usetex, fontsize, figsize, dpi):
    
    plt.rcParams['text.usetex'] = usetex
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.size'] = fontsize
    plt.rcParams['figure.figsize'] = figsize
    plt.rcParams['figure.dpi'] = dpi
    
def ax_ticks_subplot(ax):
    ax.tick_params(direction = "in")

    axT = ax.secondary_xaxis('top')
    axT.tick_params(direction = "in")
    axT.xaxis.set_ticklabels([])

    axR = ax.secondary_yaxis('right')
    axR.tick_params(direction = "in")
    axR.yaxis.set_ticklabels([])

# Bbox-Size for saving subplots
def full_extent(ax, pad):
    """Get the full extent of an axes, including axes labels, tick labels, and
    titles."""
    # For text objects, we need to draw the figure first, otherwise the extents
    # are undefined.
    ax.figure.canvas.draw()
    items = ax.get_xticklabels() + ax.get_yticklabels() 
    items += [ax, ax.title, ax.get_xaxis().get_label(), ax.get_yaxis().get_label()]
    bbox = Bbox.union([item.get_window_extent() for item in items])
    
    return bbox.expanded(1.0 + pad, 1.0 + pad)

def savefig_subplot(fig, ax, path, pad):
    #extent = full_extent(ax, pad).transformed(fig.dpi_scale_trans.inverted())
    bbox = ax.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.dpi_scale_trans.inverted())
    bbox = bbox.expanded(1.0 + pad, 1.0 + pad)
    fig.savefig(path, bbox_inches=bbox)

def eflux_time(time, eflux, figuresize):
    fig, ax = plt.subplots(figsize=figuresize)
    
    ax.plot(time, eflux)
    ax.set_xlabel(r'$t~[R/ \nu_{\mathrm{th}}]$')
    ax.set_ylabel(r'$\chi~[\rho^2 \nu_{\mathrm{th}} / R]$')
    
    ax_ticks_subplot(ax)
    
def max_shearingrate_time(time, wexb_max, fourier_index, figuresize):
    fig, ax = plt.subplots(figsize=figuresize)
    
    if type(fourier_index) == int:
        ax.plot(time,wexb_max[fourier_index], label = 'fourier mode' + str(fourier_index))
    else:
        for i in fourier_index:
            ax.plot(time,wexb_max[i], label = 'fourier mode ' + str(i))

    ax.set_xlabel(r'$t~[R/ \nu_{\mathrm{th}}]$')
    ax.set_ylabel(r'$|k_x^2 \phi|$')
    
    ax_ticks_subplot(ax)
    
    plt.legend()
    
def all_shearingrate_radialcoordinate(rad_coord, wexb, figuresize, stepsize):
    fig, ax = plt.subplots(figsize=figuresize)

    start, end = 0, stepsize - 1 
    
    while end <= wexb.shape[1]:
        
        wexb_mean = np.mean(wexb[:,start:end],1)
    
        ax.plot(rad_coord, wexb_mean)
    
        start += stepsize
        end += stepsize
    
    
    #ax.set_title(r'$R/L_T =$ ' + rlt + ', time interval [0 '+str(wexb.shape[1])+']', pad=20)
    ax.set_xlabel(r'$x[\rho]$')
    ax.set_ylabel(r'$\omega_{\mathrm{E \times B}}$')
    
    ax_ticks_subplot(ax)
    
def mean_shearingrate_radialcoordinate_amplitude(rad_coord, wexb_rad_mean, wexb_rad_middle, wexb_rad_mean_amp, wexb_rad_mean_amp_max, 
                                                 figuresize):
    
    if figuresize == (24,8):
        fig, ax = plt.subplots(1, 2, figsize=figuresize)

        # Plot shearing rate
        ax[0].plot(rad_coord, wexb_rad_mean)
        ax[0].plot(rad_coord, wexb_rad_middle, 'black', linestyle='--', linewidth=1)
        ax[0].plot(rad_coord, np.repeat(wexb_rad_mean_amp_max, len(rad_coord)), 'r', linestyle='--', linewidth=1)
        ax[0].plot(rad_coord, -np.repeat(wexb_rad_mean_amp_max, len(rad_coord)), 'r', linestyle='--', linewidth=1)
        #ax[0].set_title(r'$R/L_T =$ ' + rlt + ', time interval [' + str(start) + ' ' + str(end) + ']', pad=20)
        ax[0].set_xlabel(r'$x[\rho]$')
        ax[0].set_ylabel(r'$\omega_{\mathrm{E \times B}}$')

        ax_ticks_subplot(ax[0])

        #savefig_subplot(fig, ax[0], '../pictures/'+data+'/'+path+'/'+data+'_'+resolution+'_wexb_'+str(start)+'_'+str(end)+'.pdf', 0.02)

        # FT{shearing rate}
        ax[1].plot(rad_coord[1:], wexb_rad_mean_amp[1:])
        #ax[1].set_title(r'$R/L_T =$ ' + rlt + ', time interval [' + str(start) + ' ' + str(end) + ']', pad=20)
        ax[1].set_xlabel(r'$x[\rho]$')
        ax[1].set_ylabel(r'Amplitude')
        
        ax_ticks_subplot(ax[1])

        #savefig_subplot(fig, ax[1],'../pictures/'+data+'/'+path+'/'+data+'_'+resolution+'_Amp_Rad_'+str(start)+'_'+str(end)+'.pdf', 0.02)
    elif figuresize == (12,8):
        fig, ax = plt.subplots(figsize=figuresize)
        
        # Plot shearing rate
        ax.plot(rad_coord, wexb_rad_mean)
        ax.plot(rad_coord, wexb_rad_middle, 'black', linestyle='--', linewidth=1)
        ax.plot(rad_coord, np.repeat(wexb_rad_mean_amp_max, len(rad_coord)), 'r', linestyle='--', linewidth=1)
        ax.plot(rad_coord, -np.repeat(wexb_rad_mean_amp_max, len(rad_coord)), 'r', linestyle='--', linewidth=1)
        #ax.set_title(r'$R/L_T =$ ' + rlt + ', time interval [' + str(start) + ' ' + str(end) + ']', pad=20)
        ax.set_xlabel(r'$x[\rho]$')
        ax.set_ylabel(r'$\omega_{\mathrm{E \times B}}$')
        
        ax_ticks_subplot(ax)

# Multisubplot for evolution of shearing rate
def time_interval(interval, stepsize, distance, length_run, repetition):
    while interval[1][-1] < length_run*repetition:
    
        if interval[1][-1] % length_run == 0:
            interval[0].append(interval[1][-1])
            interval[1].append(interval[0][-1]+stepsize)
        else:
            interval[0].append(interval[1][-1]+distance)
            interval[1].append(interval[0][-1]+stepsize)
    
    return interval

def mean_shearingrate_radialcoordinate_subplots_grid(wexb, stepsize, distance, share_x, share_y):
    
    length = len(wexb[0])
    length_run = 10000
    
    x, y = 3, int(length/length_run)
    grid_x, grid_y = np.arange(x), np.arange(y)
    
    interval = [[ 500, 4000,  7000],
                [2000, 6000, 10000]]
        
    interval = time_interval(interval, stepsize, distance, length_run, y)
        
    fig, ax = plt.subplots(y, x,figsize=(12*x, 8*y), sharey=share_y, sharex=share_x, squeeze=True)
    
    return fig, ax, grid_y, grid_x, interval
      
def mean_shearingrate_radialcoordinate_subplot(rad_coord, rad_boxsize, wexb_rad_mean, wexb_rad_middle, wexb_rad_mean_amp_max, 
                                               ax, x, y, start, end):
    if y > 1:
        axis = ax[y,x]
    else:
        axis = ax[x]

    # Plot shearing rate
    axis.plot(rad_coord, wexb_rad_mean)
    #axis.plot(rad_coord, wexb_rad_middle, 'black', linestyle='--', linewidth=1)
    axis.plot(rad_coord, np.repeat(wexb_rad_mean_amp_max, len(rad_coord)), 'r', linestyle='--', linewidth=1)
    axis.plot(rad_coord, -np.repeat(wexb_rad_mean_amp_max, len(rad_coord)), 'r', linestyle='--', linewidth=1)
    #ax.set_title(r'$R/L_T =$ ' + rlt + ', time interval [' + str(start) + ' ' + str(end) + ']', pad=20)
    #ax[y,x].set_xlabel(r'$x[\rho]$')
    #ax[y,x].set_ylabel(r'$\omega_{\mathrm{E \times B}}$')
    
    y_axis_height = 0.45
    
    axis.set_ylim([-y_axis_height, y_axis_height])
    axis.set_xlim([0, rad_boxsize])
    
    ax_ticks_subplot(axis)
    
    #amp_label_height_neg = (y_axis_height - wexb_rad_mean_amp_max)/(2*y_axis_height) - 0.1
    amp_label_height_neg = 0.05
    amp_label_neg=r'$-\,A_{\mathrm{max}}$'
    
    amp_label_height_pos = 1 - amp_label_height_neg
    amp_label_pos=r'$A_{\mathrm{max}}$ = ' + format(round(wexb_rad_mean_amp_max, 2), '.2f')
    
    title = 'time interval [' + str(start) + ', ' + str(end) + ']'
    
    axis.text(0.805, amp_label_height_neg, amp_label_neg, color='r', ha='center', va='center', transform=axis.transAxes)
    axis.text(0.87, amp_label_height_pos, amp_label_pos, color='r', ha='center', va='center', transform=axis.transAxes)
    axis.text(0.5, 0.95, title, ha='center', va='center', transform=axis.transAxes)