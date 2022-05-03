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


class TotalFluxesPlot(Plottable):
    
    @classmethod
    def isPlottable(cls,h5filename):
        return True

    def plot(self,ax):
        f = h5py.File(self.h5filename, "r")
        time = f['/grid/time'].value.flatten()
        number_of_species = f['/'].attrs['grid.number_of_species'][0]

        pflux_es = f['/diagnostic/diagnos_fluxes/pflux_es'].value
        eflux_es = f['/diagnostic/diagnos_fluxes/eflux_es'].value
        vflux_es = f['/diagnostic/diagnos_fluxes/vflux_es'].value
        isp = np.clip(self.arguments['species'].get(),0,number_of_species-1)
        #ax = fig.add_subplot(100*number_of_species + 10 + 1+isp)
        charge = f['/'].attrs[('species%02d' % (isp+1)) + '.z']
        if(charge < 0):
            sp_string = 'electron'
        else:
            sp_string = 'ion (species %d)' % (isp+1)

        r_tiny = 1e-10
        
        if(np.any(np.abs(pflux_es[isp,:]) >= r_tiny)):
            ax.plot(time, pflux_es[isp,:], label="particle flux",
                    marker='x')

        if(np.any(np.abs(eflux_es[isp,:]) >= r_tiny)):
            ax.plot(time, eflux_es[isp,:], label="heat flux",
                    marker='o')

        if(np.any(np.abs(vflux_es[isp,:]) >= r_tiny)):
            ax.plot(time, vflux_es[isp,:], label="flux of toroidal momentum",
                    marker='.')

        ax.set(xlabel=r'time $t$ [$R_{ref}/v_{th,ref}$]',
               ylabel='normalised flux ')
        ax.set(title=sp_string + r' radial component of normalised flux ')
        ax.legend(loc='best')

    def specifyArguments(self, parser):
        f = h5py.File(self.h5filename, "r")
        self.parser.add_argument("--species", help="species index (0-based)",
                                 default=0,
                                 type=int)

class AllSpeciesFluxesTotalPlot(Plottable):

    ylabels = {'pflux_es':r'particle flux',
               'eflux_es':r'heat flux',
               'vflux_es':r'flux of tor. momentum',
               }
    
    @classmethod
    def isPlottable(cls,h5filename):
        return True

    def plot(self,ax):
        f = h5py.File(self.h5filename, "r")
        time = f['/grid/time'].value.flatten()
        number_of_species = f['/'].attrs['grid.number_of_species'][0]

        data = f['/diagnostic/diagnos_fluxes/' + self.arguments['dsetname'].get()].value
        for isp in range(number_of_species):
            charge = f['/'].attrs[('species%02d' % (isp+1)) + '.z']
            if(charge < 0):
                sp_string = 'electron'
            else:
                sp_string = 'ion (species %d)' % (isp+1)

            r_tiny = 0

            if(np.any(np.abs(data[isp,:]) >= r_tiny)):
                ax.plot(time, data[isp,:], label=sp_string)

        ax.set(xlabel=r'time $t$ [$R_{ref}/v_{th,ref}$]',
               ylabel='normalised flux')
        ax.set(title=r'radial component of '+self.ylabels[self.arguments['dsetname'].get()])
        ax.legend(loc='best')

    def getPossibleChoices(self):
        f = h5py.File(self.h5filename, "r")
        choices = list(self.ylabels.keys())
        choices.sort()
        dset_exists = lambda dsetname: dsetname in f
        return [c for c in choices if '/diagnostic/diagnos_fluxes/'+c in f]

    def specifyArguments(self, parser):
        f = h5py.File(self.h5filename, "r")
        self.parser.add_argument("--dsetname", help="dataset name",
                                 choices=self.getPossibleChoices(),
                                 type=str)


if(__name__ == '__main__'):
    p = TotalFluxesPlot("gkwdata.h5", parse_cmd_line_args=True)
    fig,ax = plt.subplots()
    p.plot(ax)
    fig.savefig("fluxes.png")
    plt.close()
