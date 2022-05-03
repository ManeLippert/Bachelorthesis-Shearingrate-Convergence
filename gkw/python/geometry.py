#!/usr/bin/env python3

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np
from numpy import pi
import h5py
import argparse


try:
    import namelist_python
except:
    # if the GKW_HOME/python folder is not contained in the PYTHONPATH
    import os
    sys.path.append(os.path.join(os.environ['GKW_HOME'],'python'))
    import namelist_python

try:
    # derive the plot object from the class which is part of the
    # gkwgui package, if it is available.
    import gkwgui.plottable
    Plottable = gkwgui.plottable.Plottable
except ImportError:
    Plottable = object

class SAlphaSurfPlot(Plottable):

    @classmethod
    def isPlottable(cls,h5filename):
        return True

    def plot(self,ax):
      f = h5py.File(self.h5filename, "r")

      ###################################################
      # s - alpha Geometry
      ###################################################
      
      dpsi = 0.3
      psi = np.arange(0.0,1,dpsi)
      
      dvarphi = 2*pi/10
      varphi = np.arange(0,2*pi+dvarphi,dvarphi)

      dtheta = 2*pi*0.05
      theta = np.arange(0,2*pi+dtheta,dtheta)
      
      R = 3.0
    
      ############ 3D flux surfaces

      dvarphi = 2*pi*0.05
      varphi = np.arange(0,2*pi+dvarphi,dvarphi)
      
      ax = plt.gca(projection='3d')
      
      print(psi[2:3])
      for psi_ in psi[2:3]:
        varphi = np.arange((0+psi_),
                           (2*pi-psi_),
                           2*pi*0.01) + 2*pi*0.15
        theta_,varphi_ = np.meshgrid(theta, varphi);
        x = (R + psi_*np.cos(theta_))*np.cos(varphi_);
        y = -(R + psi_*np.cos(theta_))*np.sin(varphi_);
        z = psi_*np.sin(theta_);
        c = theta_ * 0 + psi_;
        surf = ax.plot_surface(x, y, z, rstride=1, cstride=1, cmap=cm.coolwarm,
                               linewidth=0, antialiased=False)


class CircAnalyticPlot(Plottable):

    @classmethod
    def isPlottable(cls,h5filename):
        return True

    def plot(self,ax):
        f = h5py.File(self.h5filename, "r")
        ax = plt.gca(projection='polar')
        ###################################################
        # circ Geometry
        ###################################################
        
        ########### poloidal cut
        
        # transform the psi-s grid
        def circ_stheta(theta, epsilon, s):
            ret = 1.0/(2.0*np.pi) * (theta + epsilon * np.sin(theta)) - s;

        outer_epsilon = 0.3
        delta_eps = outer_epsilon*0.1
        eps = np.arange(0,outer_epsilon,delta_eps);
        delta_s = 0.01
        s = np.arange(-0.5,+0.5,delta_s)
        import scipy.optimize
        theta_start = (s + 0.5) * 2*pi
        
        #np.array(map(f, x))
        import functools
        for eps_ in eps:
            theta = np.fromiter((scipy.optimize.fsolve(func=functools.partial(circ_stheta, epsilon=eps_, s=s[i]),
                                                       x0=theta_start[i])
                                 for i in range(len(theta_start))),                                      
                                theta_start.dtype,
                                count=len(theta_start))
            theta = np.append(theta, theta[0])
            
            ax.plot(theta, np.ones(theta.shape)*eps_,'-r')

        for i, s_ in enumerate(s):
            theta = np.fromiter((scipy.optimize.fsolve(func=functools.partial(circ_stheta, epsilon=eps[ix], s=s_),
                                                       x0=theta_start[i])
                                 for ix in range(len(eps))),                                      
                                theta_start.dtype,
                                count=len(eps))
            
            ax.plot(theta, eps,'-r')

        # for i in range(0,n_s_grid):
        #     ax.plot(theta[i,:], xgrid[i,:],'-r')

        # r = eps[-1]*1.1;
        # if(np.abs(s) <= 0.01):
        #   t = "s = " + str(_s);
        # else:
        #   t = str(_s);
        #text(r*cos(theta_circle(end))-0.01, r*sin(theta_circle(end)),t) 

class FluxSurfaceThetaXPlot(Plottable):

    @classmethod
    def isPlottable(cls,h5filename):
        return True

    def plot(self,ax):
        f = h5py.File(self.h5filename, "r")
        ax = plt.gca(projection='polar')
        
        n_x_grid = f['/'].attrs['grid.nx'][0]
        n_s_grid = f['/'].attrs['grid.n_s_grid'][0]
        flux_tube = f['/'].attrs['control.flux_tube'][0] == b'T'
        if(flux_tube):
            theta = f['/geom/poloidal_angle'].value
            theta = np.pad(theta,((0,1)), mode='wrap')
            theta = theta.reshape([-1,1])
            theta = np.repeat(theta, repeats=n_x_grid, axis=1)
        else:
            theta = f['/geom/poloidal_angle'].value.reshape([n_x_grid,n_s_grid]).transpose()
            theta = np.pad(theta,((0,1),(0,0)), mode='wrap')
        theta[-1,:] += 2*np.pi
            
        #print('theta.shape:',theta.shape)
        #theta = np.append(theta, theta[0])
        #print('theta.shape:',theta.shape)
        
        if(flux_tube):
            xgrid = f['/grid/xphi'].value
            xgrid = np.pad(xgrid,((0,1),(0,0)), mode='wrap')
        else:
            xgrid = f['/grid/xgr'].value
            xgrid = xgrid.reshape([1,-1])
            xgrid = np.repeat(xgrid, repeats=n_s_grid+1, axis=0)
            
        # outer_epsilon = 0.3
        # delta_eps = outer_epsilon*0.1
        # eps = np.arange(0,outer_epsilon,delta_eps)
        # for eps_ in eps:
        #   ax.plot(theta, theta*0 + eps_,'-r')
        for ix in range(0,n_x_grid,max(1,self.arguments['step'].get())):
            ax.plot(theta[:,ix], xgrid[:,ix],'-r')
        if(self.arguments['draw_last'].get()):
            ax.plot(theta[:,-1], xgrid[:,-1],'-r')

        for i in range(0,n_s_grid):
            ax.plot(theta[i,:], xgrid[i,:],'-r')


    def specifyArguments(self, parser):
        f = h5py.File(self.h5filename, "r")
        self.parser.add_argument("--step", help="every n'th surface",
                                 default=1,
                                 type=int)
        self.parser.add_argument("--draw_last", help="draw last surface",
                                 default=True,
                                 type=bool)

class FluxSurfaceRZPlot(Plottable):
    @classmethod
    def isPlottable(cls,h5filename):
        return True

    def plot(self,ax):
        f = h5py.File(self.h5filename, "r")
        n_x_grid = f['/'].attrs['grid.nx'][0]
        n_s_grid = f['/'].attrs['grid.n_s_grid'][0]
        R = f['/geom/R'].value.reshape((n_x_grid, n_s_grid))
        Z = f['/geom/Z'].value.reshape((n_x_grid, n_s_grid))
        R = np.pad(R,((0,0),(0,1)), mode='wrap')
        Z = np.pad(Z,((0,0),(0,1)), mode='wrap')
        for ix in range(0,n_x_grid,max(1,self.arguments['step'].get())):
            ax.plot(R[ix,:],Z[ix,:],'-r')
        if(self.arguments['draw_last'].get()):
            ax.plot(R[-1,:], Z[-1,:],'-r')

        for i in range(0,n_s_grid):
            ax.plot(R[:,i],Z[:,i],'-r')
        
        ax.set(xlabel=r'$R$', ylabel=r'$Z$')
        ax.set_aspect('equal')
    
    def specifyArguments(self, parser):
        f = h5py.File(self.h5filename, "r")
        self.parser.add_argument("--step", help="every n'th surface",
                                 default=1,
                                 type=int)
        self.parser.add_argument("--draw_last", help="draw last surface",
                                 default=True,
                                 type=bool)
    

class GeomProfilePlot(Plottable):
    ylabels = {'q':r'safety factor $q$',
               'shat':r'magn. shear $\hat s$',
               }

    @classmethod
    def isPlottable(cls,h5filename):
        return True

    def plot(self,ax):
        f = h5py.File(self.h5filename, "r")
        q = f['/geom/'+self.arguments['dsetname'].get()].value
        n_x_grid = f['/'].attrs['grid.nx'][0]
        n_s_grid = f['/'].attrs['grid.n_s_grid'][0]
        q = np.pad(q, (0,n_x_grid - q.size),mode='edge')

        flux_tube = f['/'].attrs['control.flux_tube'][0] == b'T'
        if(flux_tube):
            xgrid = f['/grid/xphi'].value[0,:]
        else:
            xgrid = f['/grid/xgr'].value
        plt.plot(xgrid,q)
        plt.gca().set(xlabel=r'$x$ [$\rho_{ref}$]', ylabel=self.ylabels[self.arguments['dsetname'].get()])

    def getPossibleChoices(self):
        f = h5py.File(self.h5filename, "r")
        choices = list(self.ylabels.keys())
        choices.sort()
        dset_exists = lambda dsetname: dsetname in f
        return [c for c in choices if '/geom/'+c in f]

    def specifyArguments(self, parser):
        f = h5py.File(self.h5filename, "r")
        self.parser.add_argument("--dsetname", help="dataset name",
                                 choices=self.getPossibleChoices(),
                                 type=str)

class MillerInputProfPlot(Plottable):
    ylabels = {
        'xgr':r'flux surface label',
        'q':r'safety factor $q$',
        'shat':r'magn. shear $\hat s$',
        'kappax':r'elongation $\kappa$',
        'deltax':r'triangularity $\delta$',
        'squarex':r'squareness $\zeta$',
        'skappax':r'rad. derivative $s_\kappa$',
        'sdeltax':r'rad. derivative $s_\delta$',
        'ssquarex':r'rad. derivative $s_\zeta$',
        'Zmilx':r'$Z_{mil}$',
        'dRmilx':r'rad. derivative $\frac{\mathrm{d}R_{mil}}{\mathrm{d}r}$',
        'dZmilx':r'rad. derivative $\frac{\mathrm{d}Z_{mil}}{\mathrm{d}r}$',
        'gradpx':r'pressure gradient dp/dPsi'
    }

    input_prof_columns = ['xgr',
                          'q',
                          'shat',
                          'kappax',
                          'deltax',
                          'squarex',
                          'skappax',
                          'sdeltax',
                          'ssquarex',
                          'Zmilx',
                          'dRmilx',
                          'dZmilx',
                          'gradpx',
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
        #f = h5py.File(self.h5filename, "r")
        # parse the file into a Namelist object, and take the groups member of that.
        input_dat_filename = os.path.join(os.path.dirname(self.h5filename), 'input.dat')
        nlist = namelist_python.read_namelist_file(input_dat_filename)
        n = nlist.groups
        # that object is an ordered case-insensitive dictionary


        # n_x_grid = f['/'].attrs['grid.nx'][0]
        # n_s_grid = f['/'].attrs['grid.n_s_grid'][0]
        n_x_grid = n['gridsize']['nx']
        n_s_grid = n['gridsize']['n_s_grid']

        geometry_data = self.read_input_prof_data(os.path.join(os.path.dirname(self.h5filename),
                                                               self.arguments["input_prof_filename"].get()),
                                                  "Geometry")
        # the column containing the chosen dataset is
        icol = self.input_prof_columns.index(self.arguments['dsetname'].get())

        #flux_tube = f['/'].attrs['control.flux_tube'][0] == b'T'
        flux_tube = n['control']['flux_tube']
        if(flux_tube):
            f = h5py.File(self.h5filename, "r")
            xgrid = f['/grid/xphi'].value[0,:]
        else:
            #xgrid = f['/grid/xgr'].value
            xgrid = geometry_data[:,0]


        
        if(self.arguments["compute_rad_deriv"].get() and self.arguments['dsetname'].get() in ("sdeltax", "skappax", "ssquarex")):
            # Compute the expression from the manual for each of the
            # s_* quantities, using a computed derivative.  The
            # resulting line must be identical to the s* data from the
            # input.prof file.

            # strip the leading 's'
            corresponding_dsetname = self.arguments['dsetname'].get()[1:]
            # obtain a cubic spline representation of the curve
            icol = self.input_prof_columns.index(corresponding_dsetname)
            smoothing_cond = 0.00
            tck = scipy.interpolate.splrep(xgrid,geometry_data[:, icol], k=3, s=smoothing_cond)
            # compute the 1st derivative
            derivative_profile = scipy.interpolate.splev(xgrid, tck, der=1)
            #r = f['/geom/eps'].value * 1.0
            r = xgrid
            if(corresponding_dsetname == "kappax"):
                kappa = geometry_data[:, icol]
                data = (r/kappa) * derivative_profile
            elif(corresponding_dsetname == "deltax"):
                delta = geometry_data[:, icol]
                data = r * derivative_profile / np.sqrt(1-delta**2)
            elif(corresponding_dsetname == "squarex"):
                data = r * derivative_profile
            plt.plot(xgrid,data, marker='.')
            
        else:
            plt.plot(xgrid,geometry_data[:, icol])
        plt.gca().set(xlabel=r'$x$ [$\rho_{ref}$]', ylabel=self.ylabels[self.arguments['dsetname'].get()])

    def read_input_prof_data(self, input_prof_filename, blockname):
        import re
        with open(input_prof_filename, "r") as f:
            content = f.read()
            for block in re.split(r'^#', content, flags=re.MULTILINE):
                if(len(block) == 0):
                    continue
                header, data = block.split(sep='\n', maxsplit=1)
                # we have to use strip, because notation with and without
                # blanks may occur mixed.
                if(header.strip().startswith(blockname.strip())):
                    #parse the whole thing
                    return np.array([[float(v) for v in line.split()] for line in data.split(sep='\n') if len(line)>0])
        raise KeyError("A block '%s' could not be found in %s" % (blockname, input_prof_filename)) 

    def getPossibleChoices(self):
        import os
        input_prof_filename = os.path.join(os.path.dirname(self.h5filename), 'input.prof')
        geometry_data = self.read_input_prof_data(input_prof_filename, "Geometry")
        columns = geometry_data.shape[1]
        choices = self.input_prof_columns[:columns]
        choices.sort()
        return choices

    def specifyArguments(self, parser):
        self.parser.add_argument("--dsetname", help="dataset name",
                                 choices=self.getPossibleChoices(),
                                 type=str)
        self.parser.add_argument("--input_prof_filename", help="alternative input.prof\nfilename",
                                 default="input.prof",
                                 type=str)
        self.parser.add_argument("--compute_rad_deriv", help="draw computed derivative\nnot data",
                                 default=False,
                                 type=bool)

class PoloidalAnglePlot(Plottable):

    @classmethod
    def isPlottable(cls,h5filename):
        return True

    def plot(self,ax):
        f = h5py.File(self.h5filename, "r")
        
        n_s_grid = f['/'].attrs['grid.n_s_grid'][0]
        n_x_grid = f['/'].attrs['grid.n_x_grid'][0]
        sgrid = f['/diagnostic/diagnos_grid/sgrid'].value
        
        flux_tube = f['/'].attrs['control.flux_tube'][0] == b'T'
        if(flux_tube):
            theta = f['/geom/poloidal_angle'].value.reshape([n_s_grid,1])
        else:
            theta = f['/geom/poloidal_angle'].value.reshape([n_x_grid,n_s_grid])
            theta = theta[self.arguments['ix'].get(),:]

        ax.plot(sgrid.flatten(),theta.flatten())
        ax.set(xlabel=r'$s$', ylabel=r'poloidal angle')

    def specifyArguments(self, parser):
        f = h5py.File(self.h5filename, "r")
        self.parser.add_argument("--ix", help="radial index (0-based)",
                                 default=0,
                                 type=int)
        
class MillerAnalyticPlot(Plottable):

    FLUX_SURFACE_PLOT = "flux surface plot"
    dRdr_PLOT = "dR/dr plot"
    dZdr_PLOT = "dZ/dr plot"
    
    @classmethod
    def isPlottable(cls,h5filename):
        return True

    def plot(self,ax):
        f = h5py.File(self.h5filename, "r")

        num = 100
        theta = np.linspace(0, 2*np.pi, num=num, endpoint=True, dtype=np.float64)
        r = self.arguments['eps'].get() * self.arguments['Rmil'].get()
        if(self.arguments['quantity'].get() == self.FLUX_SURFACE_PLOT):
            try:
                R = self.arguments['Rmil'].get() + r * np.cos(theta + np.arcsin(self.arguments['delta'].get())*np.sin(theta))
            except Exception as err:
                # maybe the arcsin argument is invalid
                traceback.print_exc()
                return
            Z = self.arguments['Zmil'].get() + r * self.arguments['kappa'].get() * np.sin(theta + self.arguments['square'].get() *np.sin(2*theta))
        
            ax.plot(R, Z)
            ax.set(xlabel=r'$R$', ylabel=r'$Z$')
            ax.set_aspect('equal')
        elif(self.arguments['quantity'].get() == self.dRdr_PLOT):
            try:
                dRdr = self.arguments['dRmil'].get() + np.cos(theta + np.arcsin(self.arguments['delta'].get())*np.sin(theta)) - self.arguments['sdelta'].get() * np.cos(theta) * np.sin(theta + np.arcsin(self.arguments['delta'].get())*np.sin(theta))
            except Exception as err:
                #maybe the arcsin argument is invalid
                traceback.print_exc()
                return
            ax.plot(theta, dRdr)
            ax.set_xticks(np.linspace(0,4,num=5, endpoint=True) * np.pi/2)
            ax.set_xticklabels([r'0', r'$\frac{\pi}{2}$',r'$\pi$',r'$\frac{3\pi}{2}$',r'2$\pi$'])
            ax.set(xlabel=r'poloidal angle $\theta$', ylabel=r'rad. derivative $\frac{\mathrm{d}R}{\mathrm{d}r}$')
            
        elif(self.arguments['quantity'].get() == self.dZdr_PLOT):
            dZdr = (self.arguments['dZmil'].get()
                    + self.arguments['kappa'].get() * np.sin(theta + self.arguments['square'].get() *np.sin(2*theta))
                    + self.arguments['skappa'].get() * np.sin(theta + self.arguments['square'].get() *np.sin(2*theta))
                    + 2 * self.arguments['kappa'].get() * self.arguments['skappa'].get()
                        * np.cos(theta + self.arguments['square'].get() *np.sin(2*theta)))
            ax.plot(theta, dZdr)
            ax.set_xticks(np.linspace(0,4,num=5, endpoint=True) * np.pi/2)
            ax.set_xticklabels([r'0', r'$\frac{\pi}{2}$',r'$\pi$',r'$\frac{3\pi}{2}$r',r'2$\pi$'])
            ax.set(xlabel=r'poloidal angle $\theta$', ylabel=r'$Z$')
        
    
    def specifyArguments(self, parser):
        f = h5py.File(self.h5filename, "r")
        self.parser.add_argument("--quantity", help="Which quantity to plot",
                                 choices=[self.FLUX_SURFACE_PLOT, self.dRdr_PLOT, self.dZdr_PLOT,],
                                 type=str)
        self.parser.add_argument("--kappa", help="elongation",
                                 default=f['/'].attrs['geom.kappa'][0],
                                 type=float)
        self.parser.add_argument("--delta", help="triangularity",
                                 default=f['/'].attrs['geom.delta'][0],
                                 type=float)
        self.parser.add_argument("--square", help="squareness zeta",
                                 default=f['/'].attrs['geom.square'][0],
                                 type=float)
        self.parser.add_argument("--Rmil", help="geometric axis R",
                                 default=1.0,
                                 type=float)
        self.parser.add_argument("--Zmil", help="geometric axis Z",
                                 default=f['/'].attrs['geom.Zmil'][0],
                                 type=float)
        self.parser.add_argument("--skappa", help="rad. derivative s_kappa",
                                 default=f['/'].attrs['geom.skappa'][0],
                                 type=float)
        self.parser.add_argument("--sdelta", help="rad. derivative s_delta",
                                 default=f['/'].attrs['geom.sdelta'][0],
                                 type=float)
        self.parser.add_argument("--ssquare", help="rad. derivative s_zeta",
                                 default=f['/'].attrs['geom.ssquare'][0],
                                 type=float)
        self.parser.add_argument("--dRmil", help="rad. derivative R",
                                 default=f['/'].attrs['geom.dRmil'][0],
                                 type=float)
        self.parser.add_argument("--dZmil", help="rad. derivative Z",
                                 default=f['/'].attrs['geom.dZmil'][0],
                                 type=float)
        # self.parser.add_argument("--gradp", help="pressure gradient",
        #                          default=f['/'].attrs['geom.gradp'][0],
        #                          type=float)
        self.parser.add_argument("--q", help="safety factor",
                                 default=f['/'].attrs['geom.q'][0],
                                 type=float)
        self.parser.add_argument("--shat", help="magn. shear",
                                 default=f['/'].attrs['geom.shat'][0],
                                 type=float)
        self.parser.add_argument("--eps", help="inv. aspect ratio",
                                 default=f['/'].attrs['geom.eps'][0],
                                 type=float)

def parse_for_plottable_name():
  parser = argparse.ArgumentParser()
  import sys
  import inspect
  # get a list of all subclasses of Plottable visible here
  choices = [name for name, obj in inspect.getmembers(sys.modules[__name__], inspect.isclass) if obj is not gkwgui.plottable.Plottable and  issubclass(obj, gkwgui.plottable.Plottable)]
  parser.add_argument('--plottable',
                      help="The name of the Plottable class you want to plot.", 
                      choices=choices,
                      type=str,
                      required=True)
  args = parser.parse_known_args()[0]
  return args.plottable


if(__name__ == "__main__"):

  print(parse_for_plottable_name())

  p = SAlphaSurfPlot("gkwdata.h5", parse_cmd_line_args=True)
  fig,ax = plt.subplots()
  p.plot(ax)
  fig.savefig("s-alpha-surf.png")
  plt.close()

  p = FluxSurfaceThetaXPlot("gkwdata.h5", parse_cmd_line_args=True)
  fig,ax = plt.subplots()
  p.plot(ax)
  fig.savefig("psiconst.png")
  plt.close()

  p = PsiConstAnalyticCircPlot("gkwdata.h5", parse_cmd_line_args=True)
  fig,ax = plt.subplots()
  p.plot(ax)
  fig.savefig("circ-psiconst.png")
  plt.close()

  p = QProfilePlot("gkwdata.h5", parse_cmd_line_args=True)
  fig,ax = plt.subplots()
  p.plot(ax)
  fig.savefig("qprofile.png")
  plt.close()

    
