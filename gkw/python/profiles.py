#!/usr/bin/env python3

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np
from numpy import pi
import h5py

try:
    # derive the plot object from the class which is part of the
    # gkwgui package, if it is available.
    import gkwgui.plottable
    Plottable = gkwgui.plottable.Plottable
except ImportError:
    Plottable = object

import gkwlib

class FlutubeInBackgroundMarkerPlot(Plottable):
    quantities = ['dens',
              'temp',
              'rln',
              'rlt',
          ]
    
    @classmethod
    def isPlottable(cls,h5filename):
        return True

    def __init__(self, h5filename):
        super().__init__(h5filename)

    def plot(self,ax):
        import os
        try:
            import namelist_python
        except:
            # if the GKW_HOME/python folder is not contained in the PYTHONPATH
            sys.path.append(os.path.join(os.environ['GKW_HOME'],'python'))
            import namelist_python

        input_dat_filename = os.path.join(os.path.dirname(self.h5filename), 'input.dat')
        nlist = namelist_python.read_namelist_file(input_dat_filename)
        n = nlist.groups

        import re
        with open(input_dat_filename, 'r') as f:
            m = re.search(r'! *x = ([0-9.+-e]+)', f.read());
            if(m is None):
                x = float(self.arguments["x"].get())
            else:
                x = float(m.group(1))
                self.arguments["x"].set("%f (parsed from input.dat)" % x)
                
        for isp in range(n['gridsize']['number_of_species']):
                ax.plot(x, n['species'][isp][self.arguments["dsetname"].get()],
                        marker='.',
                        color='black')

        
    def specifyArguments(self, parser):
        self.parser.add_argument("--dsetname", help="dataset name",
                                 choices=self.quantities,
                                 type=str)
        self.parser.add_argument("--x", help="flux surface label",
                                 default=0.5,
                                 type=float)
        

class BackgroundInputProfPlot(Plottable):
    BACKGROUND_CHARGE_DENS = 'background charge density'
    BACKGROUND_CHARGE_DENS_GRADIENT = 'background charge density gradient'
    QUASI_NEUTRALITY = 'total charge density'
    QUASI_NEUTRALITY_GRADIENT = 'total charge density gradient'
    further_profiles = [
        BACKGROUND_CHARGE_DENS,
        BACKGROUND_CHARGE_DENS_GRADIENT,
        QUASI_NEUTRALITY,
        QUASI_NEUTRALITY_GRADIENT,
    ]
    
    ylabels = {
        'dens':r'background density',
        'temp':r'background temperature',
        'rln':r'background density gradient',
        'rlt':r'background temp. gradient',
        BACKGROUND_CHARGE_DENS:r'background charge density',
        BACKGROUND_CHARGE_DENS_GRADIENT:r'background charge density gradient',
        QUASI_NEUTRALITY:r'total charge density',
        QUASI_NEUTRALITY_GRADIENT:r'total charge density gradient',
    }

    input_prof_columns = ['xgr',
                          'dens',
                          'temp',
                          'rln',
                          'rlt',
    ]


    @classmethod
    def isPlottable(cls,h5filename):
        import os
        input_prof_filename = os.path.join(os.path.dirname(h5filename), 'input.prof')
        return os.path.exists(input_prof_filename)

    def __init__(self, h5filename):
        super().__init__(h5filename)

    def plot(self,ax):
        import os
        import scipy.interpolate

        complete_input_prof_data = self.read_input_prof_data(os.path.join(os.path.dirname(self.h5filename),
                                                      self.arguments["input_prof_filename"].get()))

        ions_species_counter = 0
        total_profile = None
        total_profile_dsets = (self.QUASI_NEUTRALITY, self.QUASI_NEUTRALITY_GRADIENT)
        gradient_dsets = (self.QUASI_NEUTRALITY_GRADIENT, self.BACKGROUND_CHARGE_DENS_GRADIENT)
        
        for iblock in range(len(complete_input_prof_data)):
            if(complete_input_prof_data[iblock][0].startswith("Species")):

                # check if this is indeed the input.prof entry for this species:
                mass = complete_input_prof_data[iblock][1][0,0]
                charge = complete_input_prof_data[iblock][1][0,1]
                if(charge < 0):
                    sp_string = 'electrons'
                else:
                    sp_string = 'ions'
                    ions_species_counter += 1
                    if(ions_species_counter >= 2):
                      sp_string += ' (species %d)' % ions_species_counter
                # due to the not-so-great file format of
                # input.prof there is no header line specified
                # to address the block containing the
                # background profiles. Instead, it is the next block!
                dsetname = self.arguments['dsetname'].get()
                if(dsetname in self.further_profiles):
                    if(dsetname in gradient_dsets):
                        dsetname = "rln"
                    else:
                        dsetname = "dens"
                profile = complete_input_prof_data[iblock+1][1][:,self.input_prof_columns.index(dsetname)]

                xgrid = complete_input_prof_data[iblock+1][1][:,0]
                if(self.arguments['dsetname'].get() in self.further_profiles):
                    profile *= charge
                if(self.arguments['dsetname'].get() in gradient_dsets):
                    profile *= complete_input_prof_data[iblock+1][1][:,self.input_prof_columns.index("dens")]
                    
                if(self.arguments['dsetname'].get() in total_profile_dsets):
                    if(total_profile is None):
                        total_profile = profile
                    else:
                        total_profile += profile
                else:
                    ax.plot(xgrid, profile, label=sp_string)

        if(self.arguments['dsetname'].get() in total_profile_dsets):
            ax.plot(xgrid, total_profile, label=self.ylabels[self.arguments['dsetname'].get()])
        ax.set(xlabel=r'$x$ [$\rho_{ref}$]')
        ax.set(ylabel=self.ylabels[self.arguments['dsetname'].get()])
        ax.legend(loc='best')

    def read_input_prof_data(self, input_prof_filename, blockname=None):
        import re
        import numpy as np
        ret = []

        print("Reading content of %s..." % input_prof_filename)
        with open(input_prof_filename, "r") as f:
            content = f.read()
            for block in re.split(r'^#', content, flags=re.MULTILINE):
                if(len(block) == 0):
                    continue
                header, data = block.split(sep='\n', maxsplit=1)
                # we have to use strip, because notation with and without
                # blanks may occur mixed.
                data_parsed = np.array([[float(v) for v in line.split()] for line in data.split(sep='\n') if len(line)>0])
                if(blockname is not None and header.strip().startswith(blockname.strip())):
                    #parse the whole thing
                    return data_parsed
                elif(blockname is None):
                    ret.append((header, data_parsed))

        if(blockname is not None):
            raise KeyError("A block '%s' could not be found in %s" % (blockname, input_prof_filename))
        elif(blockname is None):
            return ret

    def specifyArguments(self, parser):
        self.parser.add_argument("--dsetname", help="dataset name",
                                 choices=self.input_prof_columns[1:] + self.further_profiles,
                                 type=str)
        self.parser.add_argument("--input_prof_filename", help="alternative input.prof\nfilename",
                                 default="input.prof",
                                 type=str)

class BackgroundProfilePlot(Plottable):

    ylabels = {'background_dens':r'background density',
              'background_temp':r'background temperature',
              'background_uprim':r"background rotation gradient $u'$",
              'background_rln':r'background density gradient',
              'background_rlt':r'background temp. gradient',
              'shatx':r'magn. shear $\hat s$',
              'qx':r'safety factor $q$'}

    @classmethod
    def isPlottable(cls,h5filename):
        f = h5py.File(h5filename, "r")
        return f['/'].attrs['control.flux_tube'] == b'F'

    def plot(self,ax):
        f = h5py.File(self.h5filename, "r")
        dsetname = self.arguments['dsetname'].get()
        data = f['/diagnostic/diagnos_rad/background_profiles/'+dsetname].value
        n_x_grid = f['/'].attrs['grid.nx'][0]
        number_of_species = f['/'].attrs['grid.number_of_species'][0]

        xgrid = f['/grid/xgr'].value
        print("xgrid.shape:", xgrid.shape)
        print("data.shape:", data.shape)
        if(data.ndim == 1):
            ax.plot(xgrid,data)
        else:
            for isp in range(number_of_species):
                charge = f['/'].attrs[('species%02d' % (isp+1)) + '.z']
                if(charge < 0):
                    sp_string = 'electrons'
                else:
                    sp_string = 'ions'
                    if(number_of_species > 2 or
                       (number_of_species > 1 and f['/'].attrs['spcgeneral.adiabatic_electrons']==b'T')):

                        sp_string += ' (species %d)' % (isp+1)
                ax.plot(xgrid,data[isp, :], label=sp_string)

        ax.set(xlabel=r'$x$ [$\rho_{ref}$]',
               ylabel=self.ylabels[dsetname])
        ax.legend(loc='best')

    def specifyArguments(self, parser):
        choices = list(self.ylabels.keys())
        choices.sort()
        self.parser.add_argument("--dsetname", help="profile dataset name",
                                 choices=choices,
                                 #default='background_dens',
                                 type=str)


if(__name__ == "__main__"):

  p = BackgroundProfilePlot("gkwdata.h5", parse_cmd_line_args=True)
  fig,ax = plt.subplots()
  p.plot(ax)
  fig.savefig("background_profile.png")
  plt.close()

    


