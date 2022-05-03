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


class DominantGrowthRatePlot(Plottable):
    
    @classmethod
    def isPlottable(cls,h5filename):
        return True

    def plot(self, ax):
        f = h5py.File(self.h5filename, "r")
        time = f['/grid/time'].value.flatten()
        data = f['/diagnostic/diagnos_growth_freq/dominant_growth_rate'].value.flatten()
        ax.plot(time, data)
        ax.set(xlabel=r'time $t$ [$R_{ref}/v_{th,ref}$]', ylabel=r'growth rate $\gamma$ [$v_{th,ref}/R_{ref}$]')
        return data, time, None

    def export(self):
        fig, ax = plt.subplots()
        desc = "dominant growth rate $\gamma$ timetrace"
        data, x, y = self.plot(ax)
        return data, x, y, desc
    
        

class DominantFrequencyPlot(Plottable):
    
    @classmethod
    def isPlottable(cls,h5filename):
        return True

    def plot(self, ax):
        f = h5py.File(self.h5filename, "r")
        time = f['/grid/time'].value.flatten()
        data = f['/diagnostic/diagnos_growth_freq/dominant_real_freq'].value.flatten()
        ax.plot(time, data)
        ax.set(xlabel=r'time $t$ [$R_{ref}/v_{th,ref}$]', ylabel=r'frequency $\omega$ [$v_{th,ref}/R_{ref}$]')
        return data, time, None

    def export(self):
        fig, ax = plt.subplots()
        desc = "dominant frequency $\omega$ timetrace"
        data, x, y = self.plot(ax)
        return data, x, y, desc
        

class GrowthRateSpectrumTimeevoPlot(Plottable):
    @classmethod
    def isPlottable(cls,h5filename):
        return True

    def plot(self, ax):
        f = h5py.File(self.h5filename, "r")
        time = f['/grid/time'].value.flatten()
        data = f['/diagnostic/diagnos_growth_freq/growth_rates'].value
        # print("data.shape:", str(data.shape))
        # print("time.shape:", str(time.shape))
        for j in range(self.arguments['imod_min'].get(), self.arguments['imod_max'].get()):
            try:
                line = data[j,:]
                if(np.any(line > 0)):
                    ax.loglog(time[line >0], line[line > 0], label="mode %d" % j, linestyle='solid')
                ## give this the same color:
                #ax.loglog(time[line<0], -line[line < 0], label="mode %d" % j, linestyle='dotted')
            except:
                continue
        ax.set(xlabel=r'time $t$ [$R_{ref}/v_{th,ref}$]', ylabel=r'growth rate $\gamma$ [$v_{th,ref}/R_{ref}$]')
        ax.legend(loc='best')


    def specifyArguments(self, parser):
        self.parser.add_argument("--imod_min", help="binormal mode index (0-based)",
                                 default=0,
                                 type=int)

        f = h5py.File(self.h5filename, "r")
        nmod = f['/'].attrs['grid.nmod'][0]
        self.parser.add_argument("--imod_max", help="binormal mode index (0-based)",
                                 default=nmod,
                                 type=int)


class GrowthRateSpectrumPlot(Plottable):
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
        # take the last sample
        data = f['/diagnostic/diagnos_growth_freq/growth_rates'].value[:,-1]
        ax.plot(krho, data)
        ax.set(xlabel=xlabel, ylabel=r'growth rate $\gamma$ [$v_{th,ref}/R_{ref}$]')
        

class FrequencySpectrumPlot(Plottable):
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
        
        # take the last sample
        data = f['/diagnostic/diagnos_growth_freq/frequencies'].value[:,-1]
        ax.plot(krho, data)
        ax.set(xlabel=xlabel, ylabel=r'frequency $\omega$ [$v_{th,ref}/R_{ref}$]')

 
if(__name__ == '__main__'):
    p = DominantGrowthRatePlot("gkwdata.h5", parse_cmd_line_args=True)
    fig,ax = plt.subplots()
    p.plot(ax)
    fig.savefig("dominant_growth_rate.png")
    plt.close()

    p = DominantFrequencyPlot("gkwdata.h5", parse_cmd_line_args=True)
    fig,ax = plt.subplots()
    p.plot(ax)
    fig.savefig("dominant_real_freq.png")
    plt.close()

    p = GrowthRateSpectrumPlot("gkwdata.h5", parse_cmd_line_args=True)
    fig,ax = plt.subplots()
    p.plot(ax)
    fig.savefig("growthrate_spectr.png")
    plt.close()

    p = FrequencySpectrumPlot("gkwdata.h5", parse_cmd_line_args=True)
    fig,ax = plt.subplots()
    p.plot(ax)
    fig.savefig("frequencies_spectr.png")
    plt.close()
