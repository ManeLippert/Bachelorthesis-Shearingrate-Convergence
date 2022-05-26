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