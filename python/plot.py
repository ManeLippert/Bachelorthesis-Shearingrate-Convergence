import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.transforms import Bbox
import seaborn as sns
from cycler import cycler


# Plot parameters
def parameters(fontsize, figsize, dpi, usetex = True,
               linewidth = 3, colorpalette = sns.color_palette("colorblind", as_cmap=True),
               tickwidth = 1.5, ticklength = 10, tickspad = 10, ticksdirc = 'in', ticksadditional = True,
               legendpad = -2, legendontop = True):
    
    # TEXT ==========================================================
    
    plt.rcParams['text.usetex'] = usetex
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.size'] = fontsize
    
    # FIGURE ========================================================
    
    plt.rcParams['figure.figsize'] = figsize
    plt.rcParams['figure.dpi'] = dpi
    
    # STYLE =========================================================
    
    plt.rcParams['axes.labelpad'] = 15
    plt.rcParams['lines.linewidth'] = linewidth
    plt.rcParams['axes.prop_cycle'] = cycler('color', colorpalette)
    
    # TICKS AND FRAME ===============================================
    
    plt.rcParams['axes.linewidth'] = tickwidth
    
    plt.rcParams['xtick.top'] = ticksadditional
    plt.rcParams['xtick.bottom'] = ticksadditional
    
    plt.rcParams['ytick.left'] = ticksadditional
    plt.rcParams['ytick.right'] = ticksadditional
    
    plt.rcParams['xtick.major.size'] = ticklength
    plt.rcParams['xtick.major.width'] = tickwidth
    
    plt.rcParams['ytick.major.size'] = ticklength
    plt.rcParams['ytick.major.width'] = tickwidth
    
    plt.rcParams['xtick.minor.size'] = ticklength/2
    plt.rcParams['xtick.minor.width'] = tickwidth
    
    plt.rcParams['ytick.minor.size'] = ticklength/2    
    plt.rcParams['ytick.minor.width'] = tickwidth
    
    plt.rcParams['xtick.direction'] = ticksdirc
    plt.rcParams['ytick.direction'] = ticksdirc
    
    plt.rcParams['xtick.major.pad'] = tickspad
    plt.rcParams['ytick.major.pad'] = tickspad
    
    # LEGEND ========================================================
    
    if legendontop == True:
        plt.rcParams['legend.loc'] = 'upper center'
        plt.rcParams['legend.frameon'] = False
        #plt.rcParams['legend.handleheight'] = 1
        plt.rcParams['legend.borderpad'] = legendpad
    
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

def savefig_subplot(fig, ax, path, pad, bbox_input = None):
    #extent = full_extent(ax, pad).transformed(fig.dpi_scale_trans.inverted())
    if bbox_input == None:
        bbox = ax.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.dpi_scale_trans.inverted())
    else:
        bbox = bbox_input
        
    try:
        bbox = bbox.expanded(1.0 + pad[1], 1.0 + pad[0])
    except TypeError:
        bbox = bbox.expanded(1.0 + pad, 1.0 + pad)
        
    fig.savefig(path, bbox_inches=bbox)

def eflux_time(time, eflux, figuresize = (24,8), xlim = (0, None), ylim = (0, None), label = None, create_plot = True, axis = None):
    
    if create_plot:
        fig, axis = plt.subplots(figsize=figuresize)
    
    axis.plot(time, eflux, label = label)
    axis.set_xlabel(r'$t~[R/ \nu_{\mathrm{th}}]$')
    axis.set_ylabel(r'$\chi~[\rho^2 \nu_{\mathrm{th}} / R]$')
    
    if xlim[1] == None:
        axis.set_xlim((0, max(time)))
    else:
        axis.set_xlim(xlim)
    
    axis.set_ylim(ylim)
    
def max_shearingrate_time(time, wexb_max, fourier_index, figuresize):
    fig, ax = plt.subplots(figsize=figuresize)
    
    if type(fourier_index) == int:
        ax.plot(time,wexb_max[fourier_index], label = r'$k_' + str(fourier_index) + '$')
    else:
        for i in fourier_index:
            ax.plot(time,wexb_max[i], label = r'$k_' + str(i) + '$')

    ax.set_xlabel(r'$t~[R/ \nu_{\mathrm{th}}]$')
    ax.set_ylabel(r'$|k_\psi^2 \phi|$')
    
    ax.set_xlim(xmin=0, xmax=time[-1])
    ax.set_ylim(ymin=0)
    
    ax_ticks_subplot(ax)
    
    plt.legend(loc = 'center right')
    
def all_shearingrate_radialcoordinate(rad_coord, wexb, figuresize, stepsize):
    fig, ax = plt.subplots(figsize=figuresize)

    start, end = 0, stepsize - 1 
    
    while end <= wexb.shape[1]:
        
        wexb_mean = np.mean(wexb[:,start:end],1)
    
        ax.plot(rad_coord, wexb_mean)
    
        start += stepsize
        end += stepsize
    
    #ax.set_title(r'$R/L_T =$ ' + rlt + ', time interval [0 '+str(wexb.shape[1])+']', pad=20)
    ax.set_xlabel(r'$\psi[\rho]$')
    ax.set_ylabel(r'$\omega_{\mathrm{E \times B}}$')
    
    ax.set_xlim(xmin=0, xmax=rad_coord[-1])
    ax.set_ylim(ymin=-0.45, ymax=0.45)
    
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
        ax[0].set_xlabel(r'$\psi[\rho]$')
        ax[0].set_ylabel(r'$\omega_{\mathrm{E \times B}}$')
        
        ax[0].set_xlim(xmin=0)

        ax_ticks_subplot(ax[0])

        #savefig_subplot(fig, ax[0], '../pictures/'+data+'/'+path+'/'+data+'_'+resolution+'_wexb_'+str(start)+'_'+str(end)+'.pdf', 0.02)

        # FT{shearing rate}
        ax[1].plot(rad_coord[1:], wexb_rad_mean_amp[1:])
        #ax[1].set_title(r'$R/L_T =$ ' + rlt + ', time interval [' + str(start) + ' ' + str(end) + ']', pad=20)
        ax[1].set_xlabel(r'$\psi[\rho]$')
        ax[1].set_ylabel(r'Amplitude')
        
        ax[1].set_xlim(xmin=0)
        
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
        ax.set_xlabel(r'$\psi[\rho]$')
        ax.set_ylabel(r'$\omega_{\mathrm{E \times B}}$')
        
        ax.set_xlim(xmin=0)
        
        ax_ticks_subplot(ax)
      
def mean_shearingrate_radialcoordinate_subplot(rad_coord, rad_boxsize, wexb_rad_mean, wexb_rad_middle, wexb_rad_mean_amp_max, 
                                               ax, x, y, x_max, y_max, start_time, end_time):
    
    if y_max > 1:
        axis = ax[y,x]
    else:
        if x_max > 1:
            axis = ax[x]
        else:
            axis = ax

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
    
    title = 'time interval [' + str(start_time) + ', ' + str(end_time) + ']'
    
    axis.text(0.805, amp_label_height_neg, amp_label_neg, color='r', ha='center', va='center', transform=axis.transAxes)
    axis.text(0.87, amp_label_height_pos, amp_label_pos, color='r', ha='center', va='center', transform=axis.transAxes)
    axis.text(0.5, 0.95, title, ha='center', va='center', transform=axis.transAxes)
    
def dim_subplot(interval):
    
    for i in [3,2,1]:
        if interval.shape[1] % i == 0:
            xdim, ydim = i, int(interval[1]/i)
            break