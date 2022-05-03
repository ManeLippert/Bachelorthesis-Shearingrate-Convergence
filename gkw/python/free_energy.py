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


class FreeEnergyBalancePlot(Plottable):

    @classmethod
    def isPlottable(cls, h5filename):
        try:
            f = h5py.File(h5filename, "r")
            d = f['/diagnostic/diagnos_energetics/dt_entr_ky']
            return True
        except:
            return False

    def plot(self,h5filename):
        f = h5py.File(h5filename, "r")
        fig = plt.figure()
        time = f['/grid/time'].value.flatten()
        dt_energy = -(f['/diagnostic/diagnos_energetics/dt_entr_ky'].value +
          f['/diagnostic/diagnos_energetics/dt_entr_field_ky'].value +
          f['/diagnostic/diagnos_energetics/dt_entr_field_adiacorr_ky'].value)
        dt_energy = np.sum(dt_energy,axis=0)

        dissipation = -(f['/diagnostic/diagnos_energetics/entr_num_dis_ky'].value +
                        f['/diagnostic/diagnos_energetics/entr_num_vp_ky'].value +
                        f['/diagnostic/diagnos_energetics/entr_num_perp_ky'].value +
                        f['/diagnostic/diagnos_energetics/entr_coll_ky'].value +
                        f['/diagnostic/diagnos_energetics/entr_outflow_ky'].value +
                        f['/diagnostic/diagnos_energetics/entr_temp_src_ky'].value
                        )
        dissipation = np.sum(dissipation,axis=0)

        free_energy_source = -f['/diagnostic/diagnos_energetics/entr_source01_ky'].value
        for i in range(2,7):
          free_energy_source += -(f['/diagnostic/diagnos_energetics/entr_source%02d_ky' % i].value)
        free_energy_source = np.sum(free_energy_source,axis=0)

        total = dt_energy + dissipation + free_energy_source 
        
        plt.plot(time, dt_energy, label='$dF/dt$')
        plt.plot(time, dissipation, label='dissipation $D$')
        plt.plot(time, free_energy_source, label='free energy production $S$')
        plt.plot(time, total, label='total $\Sigma$')
        ax = plt.gca()
        print(f['/diagnostic/diagnos_energetics/dt_entr_field_ky'].attrs['description'][0])
        ax.set(xlabel=r'time $t$ [$R_{ref}/v_{th,ref}$]',
               ylabel=r'%s [$%s$]' % (f['/diagnostic/diagnos_energetics/dt_entr_ky'].attrs['description'][0].decode('utf-8'),
                                    f['/diagnostic/diagnos_energetics/dt_entr_ky'].attrs['physical unit'][0].decode('utf-8')))
        return fig

if(__name__ == '__main__'):
    p = FreeEnergyBalancePlot("gkwdata.h5", parse_cmd_line_args=True)
    fig,ax = plt.subplots()
    p.plot(ax)
    fig.savefig("free_energy_balance.png")

    plt.close()
