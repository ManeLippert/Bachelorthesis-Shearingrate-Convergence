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

def eflux_time(time, eflux, figuresize, data, resolution):
    fig, ax = plt.subplots(figsize=figuresize)
    
    ax.plot(time, eflux)
    ax.set_xlabel(r'$t~[R/ \nu_{\mathrm{th}}]$')
    ax.set_ylabel(r'$\chi~[\rho^2 \nu_{\mathrm{th}} / R]$')
    
    plt.savefig('../pictures/'+data+'/'+resolution+'/'+data+'_'+resolution+'_eflux.pdf', bbox_inches='tight')
    
def max_shearingrate_time(time, wexb_max, figuresize, data, resolution):
    fig, ax = plt.subplots(figsize=figuresize)

    ax.plot(time,wexb_max)    
    ax.set_xlabel(r'$t~[R/ \nu_{\mathrm{th}}]$')
    ax.set_ylabel(r'$|k_x^2 \phi|$')

    plt.savefig('../pictures/'+data+'/'+resolution+'/'+data+'_'+resolution+'_wexb_max.pdf', bbox_inches='tight')
    
def all_shearingrate_radialcoordinate(rad_coord, wexb, figuresize, rlt, data, resolution):
    fig, ax = plt.subplots(figsize=figuresize)

    start, end = 0, 999 
    
    while end <= wexb.shape[1]:
        
        wexb_mean = np.mean(wexb[:,start:end],1)
    
        ax.plot(rad_coord, wexb_mean)
    
        start += 1000
        end += 1000
    
    
    ax.set_title(r'$R/L_T =$ ' + rlt + ', time interval [0 '+str(wexb.shape[1])+']', pad=20)
    ax.set_xlabel(r'$x[\rho]$')
    ax.set_ylabel(r'$\omega_{\mathrm{E \times B}}$')
    
    plt.savefig('../pictures/'+data+'/'+resolution+'/'+data+'_'+resolution+'_wexb_all.pdf', bbox_inches='tight')
    
def mean_shearingrate_radialcoordinate_amplitude(rad_coord, wexb_rad_mean, wexb_rad_middle, wexb_rad_mean_amp, wexb_rad_mean_amp_max, figuresize, start, end, rlt, data, resolution):
    fig, ax = plt.subplots(1, 2, figsize=figuresize)

    # Plot shearing rate
    ax[0].plot(rad_coord, wexb_rad_mean)
    ax[0].plot(rad_coord, wexb_rad_middle, 'black', linestyle='--', linewidth=1)
    ax[0].plot(rad_coord, np.repeat(wexb_rad_mean_amp_max, len(rad_coord)), 'r', linestyle='--', linewidth=1)
    ax[0].plot(rad_coord, -np.repeat(wexb_rad_mean_amp_max, len(rad_coord)), 'r', linestyle='--', linewidth=1)
    ax[0].set_title(r'$R/L_T =$ ' + rlt + ', time interval [' + str(start) + ' ' + str(end) + ']', pad=20)
    ax[0].set_xlabel(r'$x[\rho]$')
    ax[0].set_ylabel(r'$\omega_{\mathrm{E \times B}}$')

    savefig_subplot(fig, ax[0], '../pictures/'+data+'/'+resolution+'/'+data+'_'+resolution+'_wexb_'+str(start)+'_'+str(end)+'.pdf', 0.02)

    # FT{shearing rate}
    ax[1].plot(rad_coord[1:], wexb_rad_mean_amp[1:])
    ax[1].set_title(r'$R/L_T =$ ' + rlt + ', time interval [' + str(start) + ' ' + str(end) + ']', pad=20)
    ax[1].set_xlabel(r'$x[\rho]$')
    ax[1].set_ylabel(r'Amplitude')

    savefig_subplot(fig, ax[1],'../pictures/'+data+'/'+resolution+'/'+data+'_'+resolution+'_Amp_Rad_'+str(start)+'_'+str(end)+'.pdf', 0.02)