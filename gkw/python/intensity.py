#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import h5py

try:
    # derive the plot object from the class which is part of the
    # gkwgui package, if it is available.
    import gkwgui.plottable
    Plottable = gkwgui.plottable.Plottable
except ImportError:
    Plottable = object

def find_convergence_time(h5filename, data):
    f = h5py.File(h5filename, "r")
    non_linear = f['/'].attrs['control.non_linear'][0] == b'T'
    time = f['/grid/time'].value.flatten()
    import scipy.interpolate
    # obtain a smoothed cubic spline representation of the curve
    smoothing_cond = 1.0e0
    tck = scipy.interpolate.splrep(time,data, k=3, s=smoothing_cond)
    # evaluate the profile
    smoothed_data = scipy.interpolate.splev(time, tck, der=0)
    smoothed_data_1stderiv = scipy.interpolate.splev(time, tck, der=1)
    # approximately determine the time where the step occurs. this is where the first derivative its large bump
    sm_dat_1stderiv_center_of_mass = time.dot(smoothed_data_1stderiv)/smoothed_data_1stderiv.sum()
    data = smoothed_data

    # now fit a Gaussian bell curve to the derivative of the smoothed
    # timetrace, in order to determine its width, and to avoid the
    # linear timeintervall properly.
    
    def gauss_bell(x, mu, sigma, norm): return norm/(sigma*np.sqrt(2*np.pi)) * np.exp(-0.5*((x-mu)/sigma)**2)
    
    import scipy.optimize

    # now fit a gauss bell
    time_crit_0 = sm_dat_1stderiv_center_of_mass
    if(time_crit_0 >= time[-1]):
        time_crit_0 = time[0]
    amp_0 = data[-1]
    sigma_0 = 10.0
    
    popt, pcov = scipy.optimize.curve_fit(gauss_bell, time, smoothed_data_1stderiv,
                                          p0=(time_crit_0, sigma_0, amp_0),
                                          # sigma=data*0.0,
                                          # absolute_sigma=True,
                                          #bounds=([0.0,0], [time[-1],np.max(data)]),
    )
    time_crit_opt, sigma_opt, amp_opt = popt

    time_crit_opt += sigma_opt*3.0

    itime_crit_opt = time.searchsorted(time_crit_opt)
    
    #return critical time, i.e. the time after the step has occurred.
    return time_crit_opt, itime_crit_opt


class IntensityTimetracePlot(Plottable):
    dsetname = '/diagnostic/diagnos_fields/kyspec'
    
    @classmethod
    def isPlottable(cls,h5filename):
        return True

    def plot(self, ax):
        f = h5py.File(self.h5filename, "r")
        time = f['/grid/time'].value.flatten()

        data = self.get_data(self.h5filename)
        p = ax.plot(time, data)
        p = p[0]
        unit = f[self.dsetname].attrs['physical unit'][0].decode('utf-8')
        ax.set(xlabel=r'time $t$ [$R_{ref}/v_{th,ref}$]',
               ylabel=r'intensity $|\phi|^2$ [$'+unit+'$]'
               )
        if(self.arguments['mark_average'].get()):
            # draw a thick blue vline at x=0 that spans the upper quadrant of
            # the yrange
            avg_start_time = self.arguments['avg_start'].get()
            avg_end_time = self.arguments['avg_end'].get()
            # avg_start_time = 5.0
            # avg_end_time = 15.0

            if(avg_start_time > time[-1] or avg_end_time < time[0]):
                avg_end_time = time[-1]
                conv_time, conv_itime = find_convergence_time(self.h5filename, self.get_data(self.h5filename))
                avg_start_time = conv_time
                self.arguments['avg_start'].set(avg_start_time)
                self.arguments['avg_end'].set(avg_end_time)
            
            avg_start_itime = time.searchsorted(self.arguments['avg_start'].get())
            avg_end_itime = time.searchsorted(self.arguments['avg_end'].get())

            # print("avg_start_time:", avg_start_time)
            # print("avg_end_time:", avg_end_time)
            # print("avg_start_itime:", avg_start_itime)
            # print("avg_end_itime:", avg_end_itime)


            
            
            data_mean = data[avg_start_itime:avg_end_itime].mean()
            data_std = data[avg_start_itime:avg_end_itime].std()
            xlim = ax.get_xlim()
            ax.axhline(y=data_mean,
                       xmin=avg_start_time/xlim[1],
                       xmax=avg_end_time/xlim[1],
                       color=p.get_color(),
                       linestyle='--',
            )
            
            ax.axhspan(ymin=data_mean-data_std,
                       ymax=data_mean+data_std,
                       xmin=avg_start_time/xlim[1],
                       xmax=avg_end_time/xlim[1],
                       facecolor=p.get_color(),
                       alpha=0.4)
            
                       
        return data, time, None

    def get_data(self, h5filename):
        f = h5py.File(h5filename, "r")
        data = f[self.dsetname].value
        # to correctly integrate over ky, we have to do what is called
        # the parseval_correction in the gkw srccode.
        total_timetrace = data[0,:]
        total_timetrace += 2*data[1:,:].sum(axis=0)
        data = total_timetrace

        return data

    def export(self):
        fig, ax = plt.subplots()
        desc = "total intensity $|\phi|^2$ timetrace"
        data, x, y = self.plot(ax)
        return data, x, y, desc

    def specifyArguments(self, parser):
        f = h5py.File(self.h5filename, "r")
        time = f['/grid/time'].value.flatten()
        self.parser.add_argument("--mark_average", help="mark the average intervall in the plot",
                                 default=True,
                                 type=bool)
        conv_time, conv_itime = find_convergence_time(self.h5filename, self.get_data(self.h5filename))
        self.parser.add_argument("--avg_start", help="start time for avg.",
                                 default=conv_time,
                                 type=float)
        self.parser.add_argument("--avg_end", help="end time for avg.",
                                 default=time[-1],
                                 type=float)
    
class IntensitySpectrumPlot(Plottable):

    dsetname = '/diagnostic/diagnos_fields/kyspec'
    
    @classmethod
    def isPlottable(cls,h5filename):
        f = h5py.File(h5filename, "r")
        ret = True
        ret = ret and f['/'].attrs['grid.nmod'][0] > 1
        return ret

    def plot(self, ax):
        f = h5py.File(self.h5filename, "r")
        try:
            krho = f["/grid/krho"][:,0]
            xlabel = r'$k_\theta\rho$'
        except:
            krho = f["/grid/kzeta"][:,0]
            xlabel = r'$k_\zeta$'

        if(self.arguments['average'].get()):
            time = f['/grid/time'].value.flatten()
            avg_start_time = self.arguments['avg_start'].get()
            avg_end_time = self.arguments['avg_end'].get()
            
            if(avg_start_time > time[-1] or avg_end_time < time[0]):
                avg_end_time = time[-1]
                conv_time, conv_itime = find_convergence_time(self.h5filename, self.get_data(self.h5filename))
                avg_start_time = conv_time
                self.arguments['avg_start'].set(avg_start_time)
                self.arguments['avg_end'].set(avg_end_time)
                
            avg_start_itime = time.searchsorted(self.arguments['avg_start'].get())
            avg_end_itime = time.searchsorted(self.arguments['avg_end'].get())

                
            # take the time average
            data = f[self.dsetname].value[:,avg_start_itime:avg_end_itime]
            print(data.shape)
            time_axis = data.ndim-1
            data_mean = data.mean(axis=time_axis)
            data_std = data.std(axis=time_axis)
        else:
            # take the last sample
            data_mean = f[self.dsetname].value[:,-1]
            data_std=None
        if(self.arguments['stddev_bars'].get() and data_std is not None):
            print(krho.shape)
            print(data_mean.shape)
            print(data_std.shape)
            ax.errorbar(krho, data_mean, yerr=data_std)
        else:
            ax.plot(krho, data_mean)
        unit = f[self.dsetname].attrs['physical unit'][0].decode('utf-8')
        ax.set(xlabel=xlabel,
               ylabel=r'intensity $|\phi|^2$ [$'+unit+'$]'
        )

    def get_data(self, h5filename):
        f = h5py.File(h5filename, "r")
        data = f[self.dsetname].value
        # to correctly integrate over ky, we have to do what is called
        # the parseval_correction in the gkw srccode.
        total_timetrace = data[0,:]
        total_timetrace += 2*data[1:,:].sum(axis=0)
        data = total_timetrace

        return data

    def specifyArguments(self, parser):
        f = h5py.File(self.h5filename, "r")
        time = f['/grid/time'].value.flatten()
        self.parser.add_argument("--average", help="True: average over interval;False: last sample",
                                 default=True,
                                 type=bool)
        conv_time, conv_itime = find_convergence_time(self.h5filename, self.get_data(self.h5filename))
        self.parser.add_argument("--avg_start", help="start time for avg.",
                                 default=conv_time,
                                 type=float)
        self.parser.add_argument("--avg_end", help="end time for avg.",
                                 default=time[-1],
                                 type=float)
        self.parser.add_argument("--stddev_bars", help="plot stddev bars",
                                 default=True,
                                 type=bool)

    
        
        
 
if(__name__ == '__main__'):
    raise NotImplementedError("Standalone functionality has not yet been added for this, sorry.")
    
