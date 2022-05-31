#Import modules
import numpy as np
import gkw
import derivative
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.transforms import Bbox

# From Florian Rath
def get_eflux_time(hdf5_file):
    
    # find out number of time steps
    time = hdf5_file['diagnostic/diagnos_growth_freq/time'][()]
    nt = time.shape[1]
    
    # name of the dataset
    node_name = 'diagnostic/diagnos_fluxes/eflux_species01'
    
    # load data into array
    data = np.transpose(hdf5_file[node_name][()])
    
    # reshape GKW flux ordering
    flux = np.reshape(data,(nt,2))[:,0]
    
    return flux, time[0]

def get_shearingrate_radialcoordinate_radialboxsize_ddphi_dx_zonalpot(hdf5_file):
    
    # Elektrostatic potencial
    phi = hdf5_file[gkw.find_key(hdf5_file, 'phi')][()]
    nx = phi.shape[0]

    # Mean over y to get a approximation for the zonal potenzial
    zonal_pot = np.mean(phi,1)

    # Finite Difference for shearing rate omega_ExB

    # Stepsize
    rad_boxsize = hdf5_file[gkw.find_key(hdf5_file, 'lxn')][()][0]
    rad_coord = hdf5_file[gkw.find_key(hdf5_file,'xphi')][0,:]
    dx = rad_boxsize/nx

    ddphi= derivative.finite_second_order(zonal_pot[:,:], dx, 'period')
    wexb = 0.5 * ddphi
    
    return wexb, rad_coord, rad_boxsize, ddphi, dx, zonal_pot

def get_max_shearingrate(wexb, time):
    
    wexb_max = []

    for time_point in range(len(time)):
        wexb_time_point = wexb[:,time_point]
        wexb_fft = np.fft.fft(wexb_time_point)
        wexb_fft_amp = 2/len(wexb_fft) * np.abs(wexb_fft)
        wexb_max.append(wexb_fft_amp[1])

    wexb_max = np.array(wexb_max)
    
    return wexb_max

def get_mean_middle_shearingrate(start, end, wexb):
    middle = int((end - start)/2 + start)

    # Shearing rate with mean over time
    wexb_rad_mean = np.mean(wexb[:,start:end],1)
    wexb_rad_middle = wexb[:,middle]
    
    return wexb_rad_mean, wexb_rad_middle

def get_fft_mean_max_shearingrate_amplitude(wexb_mean):
    wexb_rad_mean_fft = np.fft.fft(wexb_mean)
    wexb_rad_mean_amp = 2/ wexb_rad_mean_fft.shape[0] * np.abs(wexb_rad_mean_fft)
    wexb_rad_mean_amp_max = max(wexb_rad_mean_amp)
    
    return wexb_rad_mean_amp, wexb_rad_mean_amp_max