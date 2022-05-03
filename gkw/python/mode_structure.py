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

class ModeStructurePlot(Plottable):

    ylabels = {'phi':r'es. potential $\phi$',
               'Apar':r'vector potential component $A_\parallel$',
               'Bpar':r'B-field flutter $B_\parallel$',
               'dens':r'density',
               'T':r"temperature",
               'Tpar':r"parallel temp. $T_\parallel$",
               'Tperp':r"perpendicular temp. $T_\perp$",
               'vpar':r'current',
               }
    
    @classmethod
    def isPlottable(cls,h5filename):
        return True

    def plot(self,ax):
        f = h5py.File(self.h5filename, "r")
        n_x_grid = f['/'].attrs['grid.nx'][0]
        n_s_grid = f['/'].attrs['grid.n_s_grid'][0]
        nmod = f['/'].attrs['grid.nmod'][0]

        flux_tube = f['/'].attrs['control.flux_tube'][0] == b'T'
        spectral_radius = f['/'].attrs['control.spectral_radius'][0] == b'T'
        mode_box = f['/'].attrs['mode.mode_box'][0] == b'T'
        if(mode_box and spectral_radius):
            s = f['/geom/s_grid'].value
            s = np.expand_dims(s,axis=0)
            s = np.repeat(s, repeats=nmod, axis=0)
            s = np.expand_dims(s,axis=2)
            s = np.repeat(s, repeats=n_x_grid, axis=2)
        else:
            s = f['/geom/s_grid'].value

        print("s.shape:", str(s.shape))
        
        def rotate_for_maxval_eq_1(rotate_every_plot, field):
            if(rotate_every_plot):
                i = np.argmax(np.abs(field))
                rotate = np.abs(field.flat[i]) / field.flat[i]
            else:
                rotate = 1
            return field * rotate


        imod = self.arguments['imod'].get()
        if(spectral_radius):
            # ikxspace may be larger than nx/2...
            ikxspace = f['/'].attrs['mode.ikxspace'][0]
            ikxspace_effective = min(ikxspace,np.ceil(n_x_grid/2))
            mode_label = f['/diagnostic/diagnos_grid/mode_label'].value
            print("mode_label.shape:", str(mode_label.shape))
            kxrh = f['/grid/kxrh'];
            ixzero = np.argmin(np.abs(kxrh))
        else:
            ixzero = 0

        ieiv = self.arguments['eigenmode'].get()
        
        ix = self.arguments['delta_ix'].get()
        
        if(mode_box and spectral_radius):
            radially_connected_mask = mode_label[:,imod] == mode_label[ixzero+ix,imod]
            print("radially_connected_mask.shape:", str(radially_connected_mask.shape))
            n_radially_connected = np.count_nonzero(radially_connected_mask)
            if(n_radially_connected % 2 != 1):
                raise Exception("shouldnt there be an odd number of radial modes connected, always?")
            # from -something ...-1 0 +1 .. +something
            #with_ixzero_connected_ix = find(radially_connected_mask)
            staircase = np.arange(1,n_radially_connected+1) - np.ceil(n_radially_connected/2.0)
            
            staircase = np.expand_dims(staircase,axis=1)
            staircase = np.repeat(staircase, repeats=n_s_grid, axis=1)
            s = s[0,:,radially_connected_mask]+staircase
        else:
            # there is no radial connection between modes for nonspectral or if the modes are not equidistant.
            # In this case just pick one line of data (aka. "a profile").
            radially_connected_mask = self.arguments['delta_ix'].get()

        isp = self.arguments['species'].get()
        charge = f['/'].attrs[('species%02d' % (isp+1)) + '.z']
        if(charge < 0):
            sp_string = 'electron'
        else:
            sp_string = 'ion'
            if(f['/'].attrs['grid.number_of_species'][0] > 2 or
               (f['/'].attrs['grid.number_of_species'][0] > 1 and f['/'].attrs['spcgeneral.adiabatic_electrons']==b'T')):
                
                sp_string += ' (species %d)' % (isp+1)

        dsetname = self.arguments['dsetname'].get()
        data = f['/diagnostic/diagnos_mode_struct/'+dsetname+'_real'].value + 1j*f['/diagnostic/diagnos_mode_struct/'+dsetname+'_imag'].value
        rotate_every_plot = True
        if(data.ndim == 5):
            ###### PLOT A MOMENT
            #print("NOTE that radial and parallel dimension are interchanged, comparing fields and moments data!")

            data = np.sum(data, axis=2, keepdims=True);

            # print("s.shape:", str(s.shape))
            # print("data.shape:", str(data.shape))
            # print("data_.shape:", str(data[radially_connected_mask, :,0,isp,ieiv].shape))

            data[radially_connected_mask, :,0,isp,ieiv] = rotate_for_maxval_eq_1(rotate_every_plot,data[radially_connected_mask, :,0,isp,ieiv])

            ax.plot(s.flatten(), np.real(data[radially_connected_mask, :,0,isp,ieiv]).flatten(),marker='x', label="real")
            ax.plot(s.flatten(), np.imag(data[radially_connected_mask, :,0,isp,ieiv]).flatten(),marker='.', label="imag")
            ax.set(xlabel="s",ylabel= sp_string+" " + self.ylabels[dsetname])
            ax.legend(loc='best')

        else:
            ###### PLOT A FIELD
            print("data.shape:", str(data.shape))
            data = np.sum(data, axis=2, keepdims=True);
            data[:,radially_connected_mask,0,ieiv] = rotate_for_maxval_eq_1(rotate_every_plot,data[:,radially_connected_mask,0,ieiv]);

            ax.plot(s.flatten(),
                 np.real(data[:,radially_connected_mask,0,ieiv]).transpose().flatten(),marker='x', label='real')
            ax.plot(s.flatten(),
                 np.imag(data[:,radially_connected_mask,0,ieiv]).transpose().flatten(),marker='.', label='imag')
            ax.set(xlabel="s",ylabel=self.ylabels[dsetname])
            ax.legend(loc='best');
            #print_all_formats(["mode_struct",imodstr,eigmodestr,"line_phi","_ixzero+",num2str(ix,"%1d")])


        
        # ax = plt.subplot(111)
        # staircase = np.repeat(np.arange(-(n_x_grid-1)/2, +(n_x_grid-1)/2+1),repeats=n_s_grid,axis=0)
        # staircase = np.repeat(np.arange(0,n_x_grid),repeats=n_s_grid,axis=0)
        # imod = 1
        # print("kram:")
        # print(staircase.shape)
        # print(staircase)
        # print(data3d[:,:,imod].flatten('F').shape)
        # print(s.flatten('F').shape)
        # print(staircase + s.flatten('F'))
        # plt.plot(staircase + s.flatten('F'),data3d[:,:,imod].flatten('F'), label='raw')
        # plt.plot(r_i[0,:],q_i[0,:], label='interp')
        # ax.set(xlabel='s', ylabel='data')
        # ax.legend(loc='best')

    def getPossibleChoices(self):
        f = h5py.File(self.h5filename, "r")
        choices = list(self.ylabels.keys())
        choices.sort()
        dset_exists = lambda dsetname: dsetname in f
        return [c for c in choices if '/diagnostic/diagnos_mode_struct/'+c+'_real' in f]
        
    def specifyArguments(self, parser):
        f = h5py.File(self.h5filename, "r")
        self.parser.add_argument("--dsetname", help="dataset name",
                                 choices=self.getPossibleChoices(),
                                 type=str)
        self.parser.add_argument("--imod", help="binormal mode index (0-based)",
                                 default=f['/'].attrs['grid.nmod'][0]-1,
                                 type=int)
        self.parser.add_argument("--eigenmode", help="eigenmode index (0-based)",
                                 default=0,
                                 type=int)
        self.parser.add_argument("--species", help="species index (0-based)",
                                 default=0,
                                 type=int)
        self.parser.add_argument("--delta_ix", help="radial index delta (0-based)",
                                 default=0,
                                 type=int)
    
    
if(__name__ == '__main__'):
    p = ModeStructurePlot("gkwdata.h5", parse_cmd_line_args=True)
    fig,ax = plt.subplots()
    p.plot(ax)
    fig.savefig("mode_struct.png")
    plt.close()
