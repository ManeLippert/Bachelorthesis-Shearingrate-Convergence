#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
from numpy import pi

import h5py
import scipy
import scipy.optimize
import scipy.signal
import scipy.fftpack
import scipy.interpolate

try:
    # derive the plot object from the class which is part of the
    # gkwgui package, if it is available.
    import gkwgui.plottable
    Plottable = gkwgui.plottable.Plottable
except ImportError:
    Plottable = object

class XSContourPlot(Plottable):

    ylabels = {'Poten':r'es. potential',
               'Apara':r'vector potential component $A_\parallel$',
               'Bpara':r'B-field flutter $B_\parallel$',
           }

    @classmethod
    def isPlottable(cls,h5filename):
        return len(cls.getPossibleChoices(h5filename)) > 0

    def plot(self,ax):
        f = h5py.File(self.h5filename, "r")
        flux_tube = f['/'].attrs['control.flux_tube'][0] == b'T'
        if(flux_tube):
            xgrid = f['/grid/xphi'].value[0,:]
        else:
            xgrid = f['/grid/xgr'].value
        sgrid = f['/diagnostic/diagnos_grid/sgrid'].value

        itime = self.arguments['itime'].get()
        dsetname = "/diagnostic/diagnos_fields/"+self.arguments['dsetname'].get()+"%08d" % itime
        data3d = f[dsetname].value
        data2d = np.sum(data3d, axis=2)
        print('xgrid.shape:',xgrid.shape)
        print('sgrid.shape:',sgrid.shape)
        print('data2d.shape:',data2d.shape)
        c = ax.contourf(xgrid, sgrid, data2d)
        plt.colorbar(c)
        ax.set(ylabel='s [dimensionless]',
               xlabel=r'x [$\rho_{ref}$]',
               title=self.ylabels[self.arguments['dsetname'].get()]+', timestep %d' %itime
        )

    @classmethod
    def getPossibleChoices(cls, h5filename):
        f = h5py.File(h5filename, "r")
        choices = list(cls.ylabels.keys())
        choices.sort()
        dset_exists = lambda dsetname: dsetname in f
        # FIXME better check if there is any, with any itime!
        return [c for c in choices if ('/diagnostic/diagnos_fields/'+c+'%08d' % 1) in f]
    
    def specifyArguments(self, parser):
        f = h5py.File(self.h5filename, "r")
        ntime_end = f['/grid/time'].value.shape[1]
        self.parser.add_argument("--itime", help="time step (0-based)",
                                 default=ntime_end,
                                 type=int)

        c = self.getPossibleChoices(self.h5filename)
        self.parser.add_argument("--dsetname", help="dataset name",
                                 choices=c,
                                 default=c[0],
                                 type=str)

    def getDefault_Title(self):
        return self.ylabels[self.arguments['dsetname'].get()]+', timestep %d' % itime



class ModeStructurePlot(Plottable):

    DSET_PHI = 'phi'
    ylabels = {DSET_PHI:r'es. potential $\phi$',
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
        return len(cls.getPossibleChoices(h5filename)) > 0

    def plot(self,ax):
        f = h5py.File(self.h5filename, "r")

        dsetname = "/diagnostic/diagnos_mode_struct/" + self.arguments['dsetname'].get()
        flux_tube = f['/'].attrs['control.flux_tube'][0] == b'T'
        if(flux_tube):
            rhostar = self.arguments['rhostar'].get()
            r_boxsize = self.arguments['r_boxsize'].get()
        else:
            rhostar = f['/'].attrs['spcgeneral.rhostar'][0]
            r_boxsize = None
        spectral_radius = f['/'].attrs['control.spectral_radius'][0] == b'T'
        if(spectral_radius):
            delta_ix = self.arguments['delta_ix'].get()
        else:
            delta_ix = None
        ballooning_angle, data, x, y = plot_poloidal_slice(ax, f, dsetname,
                                       normalise=self.arguments['normalise'].get(),
                                       imod=self.arguments['imod'].get(),
                                       isp=self.arguments['species'].get(),
                                       ieiv=self.arguments['eigenmode'].get(),
                                       interpnum=self.arguments['interpnum'].get(),
                                       rhostar=rhostar,
                                       r_boxsize=r_boxsize,
                                       delta_ix=delta_ix,
                                       part=self.arguments['part'].get(),
                                       mark_ballooning_angle=self.arguments['mark_ballooning_angle'].get(),
                                       ballooning_angle=self.arguments['ballooning_angle'].get(),
                                       detect_ballooning_angle=self.arguments['detect_ballooning_angle'].get(),
        )
        ax.set(title=self.ylabels[self.arguments['dsetname'].get()])
        ax.set_aspect('equal')
        self.arguments['ballooning_angle'].set(ballooning_angle)

        return data, x, y

    def export(self):
        fig, ax = plt.subplots()
        desc = "data, as well as associated x and y coordinates for a contour plot of a perpendicular slice (i.e. a slice through a phi=const. plane)."
        data, x, y = self.plot(ax)
        return data, x, y, desc


    @classmethod
    def getPossibleChoices(cls, h5filename):
        f = h5py.File(h5filename, "r")
        choices = list(cls.ylabels.keys())
        choices.sort()
        return [c for c in choices if '/diagnostic/diagnos_mode_struct/'+c+'_real' in f]

    def specifyArguments(self, parser):
        f = h5py.File(self.h5filename, "r")
        flux_tube = f['/'].attrs['control.flux_tube'][0] == b'T'
        non_linear = f['/'].attrs['control.non_linear'][0] == b'T'
        spectral_radius = f['/'].attrs['control.spectral_radius'][0] == b'T'
        c = self.getPossibleChoices(self.h5filename)
        self.parser.add_argument("--dsetname", help="dataset name",
                                 choices=c,
                                 default=self.DSET_PHI,
                                 type=str)
        c = ["real", "imag","abs"]
        self.parser.add_argument("--part", help="part",
                                 choices=c,
                                 default=c[0],
                                 type=str)
        self.parser.add_argument("--imod", help="binormal mode index (0-based)",
                                 default=f['/'].attrs['grid.nmod'][0]-1,
                                 type=int)
        if(f['/'].attrs['control.method'] == b'EIV'):
            # select the first eigenmode, should be the most unstable
            eigenvalue_default = 0
        else:
            # select the last eigenmode, i.e. the mode structure at
            # the end of the last restarted run
            eigenvalue_default = f['/diagnostic/diagnos_mode_struct/phi_real'].value.shape[-1] - 1

        self.parser.add_argument("--eigenmode", help="eigenmode index (0-based)",
                                 default=eigenvalue_default,
                                 type=int)
        self.parser.add_argument("--species", help="species index (0-based)",
                                 default=0,
                                 type=int)
        if(spectral_radius):
            self.parser.add_argument("--delta_ix", help="radial index delta (0-based)",
                                     default=0,
                                     type=int)
        self.parser.add_argument("--interpnum", help="interpolation factor\n(recommended 8-16)",
                                 default=10,
                                 type=int)
        if(flux_tube):
            self.parser.add_argument("--rhostar", help="rhostar",
                                     default=0.001,
                                     type=float)
            self.parser.add_argument("--r_boxsize", help="radial boxsize [rhostar]",
                                     default=80,
                                     type=float)
        self.parser.add_argument("--normalise", help="normalise",
                                 default=non_linear,
                                 type=bool)
        
        self.parser.add_argument("--mark_ballooning_angle", help="mark the ballooning angle with a line",
                                 default=True,
                                 type=bool)
        self.parser.add_argument("--ballooning_angle", help="the ballooning angle",
                                 default=0.0,
                                 type=float)
        self.parser.add_argument("--detect_ballooning_angle", help="detect the ballooning angle automatically",
                                 default=True,
                                 type=bool)
    

        
class RealSpace3dDataPlot(Plottable):

    ylabels = {'Poten':r'es. potential',
               'Apara':r'vector potential component $A_\parallel$',
               'Bpara':r'B-field flutter $B_\parallel$',
           }
    
    @classmethod
    def isPlottable(cls,h5filename):
        return len(cls.getPossibleChoices(h5filename)) > 0

    def plot(self,ax):
        f = h5py.File(self.h5filename, "r")
        spectral_radius = f['/'].attrs['control.spectral_radius'][0] == b'T'
        itime = self.arguments['itime'].get()
        dsetname = "/diagnostic/diagnos_fields/"+self.arguments['dsetname'].get()+"%08d" % itime

        flux_tube = f['/'].attrs['control.flux_tube'][0] == b'T'
        if(flux_tube):
            rhostar = self.arguments['rhostar'].get()
            r_boxsize = self.arguments['r_boxsize'].get()
        else:
            rhostar = f['/'].attrs['spcgeneral.rhostar'][0]
            r_boxsize = None
        if(spectral_radius):
            delta_ix = self.arguments['delta_ix'].get()
        else:
            delta_ix = None
        
        ballooning_angle, data, x, y = plot_poloidal_slice(ax, f, dsetname,
                                       normalise=self.arguments['normalise'].get(),
                                       imod=self.arguments['imod'].get(),
                                       isp=self.arguments['species'].get(),
                                       ieiv=0,
                                       interpnum=self.arguments['interpnum'].get(),
                                       rhostar=rhostar,
                                       r_boxsize=r_boxsize,
                                       delta_ix=delta_ix,
                                       part=self.arguments['part'].get(),
                                       mark_ballooning_angle=self.arguments['mark_ballooning_angle'].get(),
                                       ballooning_angle=self.arguments['ballooning_angle'].get(),
                                       detect_ballooning_angle=self.arguments['detect_ballooning_angle'].get(),
        )
        ax.set(
               title=self.ylabels[self.arguments['dsetname'].get()]+', timestep %d' %itime
        )
        ax.set_aspect('equal')

    @classmethod
    def getPossibleChoices(cls, h5filename):
        f = h5py.File(h5filename, "r")
        choices = list(cls.ylabels.keys())
        choices.sort()
        dset_exists = lambda dsetname: dsetname in f
        # FIXME better check if there is any, with any itime!
        ntime_end = f['/grid/time'].value.shape[1]
        return [c for c in choices if ('/diagnostic/diagnos_fields/'+c+'%08d' % ntime_end) in f]
    
    def specifyArguments(self, parser):
        f = h5py.File(self.h5filename, "r")
        ntime_end = f['/grid/time'].value.shape[1]
        spectral_radius = f['/'].attrs['control.spectral_radius'][0] == b'T'

        self.parser.add_argument("--itime", help="time step (0-based)",
                                 default=ntime_end,
                                 type=int)

        c = self.getPossibleChoices(self.h5filename)
        self.parser.add_argument("--dsetname", help="dataset name",
                                 choices=c,
                                 default=c[0],
                                 type=str)
        c = ["real", "imag","abs"]
        self.parser.add_argument("--part", help="part",
                                 choices=c,
                                 default=c[0],
                                 type=str)
        self.parser.add_argument("--imod", help="binormal mode index (0-based)",
                                 default=f['/'].attrs['grid.nmod'][0]-1,
                                 type=int)
        self.parser.add_argument("--species", help="species index (0-based)",
                                 default=0,
                                 type=int)
        if(spectral_radius):
            self.parser.add_argument("--delta_ix", help="radial index delta (0-based)",
                                     default=0,
                                     type=int)
        self.parser.add_argument("--interpnum", help="interpolation factor",
                                 default=10,
                                 type=int)
        flux_tube = f['/'].attrs['control.flux_tube'][0] == b'T'
        if(flux_tube):
            self.parser.add_argument("--rhostar", help="rhostar",
                                     default=0.001,
                                     type=float)
            self.parser.add_argument("--r_boxsize", help="radial boxsize [rhostar]",
                                     default=80,
                                     type=float)
        self.parser.add_argument("--normalise", help="normalise",
                                 default=True,
                                 type=bool)
        self.parser.add_argument("--mark_ballooning_angle", help="mark the ballooning angle with a line",
                                 default=True,
                                 type=bool)
        self.parser.add_argument("--ballooning_angle", help="the ballooning angle",
                                 default=0.0,
                                 type=float)
        self.parser.add_argument("--detect_ballooning_angle", help="detect the ballooning angle automatically",
                                 default=True,
                                 type=bool)
        



    # % Adds first s-row at the end of each block (compared to version 5) and prof_back
    # % is no function parameter (as it is not used anyway).
    
    # % Input parameters:
    # %   parallel:  Array created from loading the parallel.dat file.
    # %   geom:      Structure with the geometry output (as from read_geom).
    # %              Is used for the gridsize and the safety factor profile.
    # %   rhostar:   The value of rhostar used for the simulation.
    # %   n:         The value  of the modenumber used for the simulation.
    # %   field:     The field that should be extracted from the file.
    # %              1 = perturbed potential (default).
    # %              2 = || component of vector potential.
    # %              3 = perturbed density.
    # %              4 = perturbed parallel temperature.
    # %              5 = perturbed perpendicular temperature.
    # %              6 = perturbed parallel flow velocity.
    # %              7 = perturbed parallel (compressional) magnetic field.
    # %              \note As real and imaginary part are extracted the field
    # %                    number and the column number differ.
    # %   write:   A string, if this is not empty (the default) then is it used
    # %            as a name to save the fields and grids to.
    # %
    # %
    # % Output parameters:
    # %   a structure with the following fields:
    # %   \note : All arrays (fields/grids) that are returned, will have a size of ns times nx.
    # %   theta:     Grid with the theta coordinate. 
    # %   r:         Grid with the r/epsilon/psi coordinate
    # %   q:         Array with the value of the safety factor for each point.
    # %   s:         Grid with the s-coordinate.
    # %   x:         Grid with the x-coordinate.
    # %   y:         Grid with the y-coordinate.
    # %   ind:       Grid with the indices that correspond to the radial direction.
    # %   zeta:      Grid  with the zeta-coordinate/difference.
    # %   z2:        Array with the rescaled absolute values of the field.
    # %   rzc:       Array with the real part of the perturbed field.
    # %   izc:       Array with the imaginary part of the perturbed field.
    # %   rzc2:      Array with the real part of the perturbed rescaled field.
    # %   izc2:      Array with the imgainary part of the perturbed rescaled field.
    # %
    # %   \note For working with (parallel.dat) files that contain multiple
    # %         species we do not use the complete data, but only the first
    # %         ns*nx entries, that correspond to the first species.
    # %   \todo Determine as many parameters as possible itself.
    # %   \todo Use a method for determining thetazero that works.
    # % ----------------------------------------------------------------------

    # function [phi] = GetDataForPhi5b(parallel, geom, rhostar, n, field = 1, write = '')
    #   ns = geom.ns;
    #   nx = geom.nx;

    #   phi.kzeta = 2*pi*rhostar*n;

    #   phi.r = geom.eps;
    #   phi.r = repmat(phi.r.', ns, 1);% r-coordinate does not depend on s, so simply repeat the vector ns times.

    #   phi.q = repmat((geom.q).', ns, 1);

    #   % Get the s-coordinate from the data.
    #   phi.s = reshape(parallel(1:ns*nx, 1), ns, nx);

    #   % Get the real/imaginary part of the perturbed potential.
    #   phi.rzc = reshape(parallel(1:ns*nx, 2*field + 0), ns, nx);
    #   phi.izc = reshape(parallel(1:ns*nx, 2*field + 1), ns, nx);
    #   % Build the complex potential, and calculate the absolute value.
    #   zc  = complex(phi.rzc, phi.izc);
    #   z   = abs(zc);

    #   %kzeta = (2*pi/rhostar)*krho.*(r(1,:))./(q(1,:));
    #   %kzeta = repmatk(zeta,ns,1);

    #   % For x(radial)-direction we use at first a simple index.
    #   phi.ind = reshape(floor((0:(ns*nx-1))/ns), ns, nx);

    #   phi.theta = reshape(geom.poloidal_angle, ns, nx);

    #   % Compute x- and y-(cartesian)coordinate.
    #   phi.x = phi.r.*cos(phi.theta);
    #   phi.y = phi.r.*sin(phi.theta);

    #   phi.zeta = (phi.q/pi).* atan(sqrt((1-phi.r)./(1+phi.r)).*tan(phi.theta/2));
    #   phi.z2 = z.*abs(sin(phi.kzeta.*phi.zeta));
    #   zc2 = zc.*exp(complex(0,phi.kzeta*phi.zeta/rhostar));

    #   phi.izc2 = imag(zc2);
    #   phi.rzc2 = real(zc2);

    # endfunctio
    # contour(phi.x, phi.y, phi.rzc2);

def plot_poloidal_slice(ax, f, dsetname, rhostar = 0.001, interpnum = 16, imod = -1, normalise=True,
                        isp=0,
                        ieiv=-1,
                        r_boxsize=80,
                        delta_ix=0,
                        part="abs",
                        mark_ballooning_angle=False,
                        ballooning_angle=0.0,
                        detect_ballooning_angle=True,
                    ):

    x,y,data = get_poloidal_slice_xy(f,dsetname,rhostar,interpnum,imod,normalise,isp,ieiv,r_boxsize,delta_ix)

    #ax.set_aspect('equal')
    cm = 'RdBu'
    cm = 'PuOr'
    #cm = 'coolwarm'
    cm = 'BrBG'
    #cm = 'rainbow'

    print("part:", part)
    if(part == "real"):
        dat = data.real
    elif(part == "imag"):
        dat = data.imag
    elif(part == "abs"):
        dat = np.abs(data)

    name = 'poloidal slice'
    if(dsetname.startswith("/diagnostic/diagnos_fields")):
        token = dsetname.replace('/diagnostic/diagnos_fields/','')
        filename = name + '_%s_%s.png' % (token, part)
    else:
        token = dsetname.replace('/diagnostic/diagnos_mode_struct/','')
        filename = name + '_mode_struct_%s_%s.png' % (token, part)
    filename = filename.replace(' ','_').replace('/','_')

    if(ax is not None):
        print('x.shape:', x.shape, x.dtype)
        print('y.shape:', y.shape, y.dtype)
        print('dat.shape:', dat.shape, dat.dtype)
        c = ax.pcolor(x,y,dat, cmap=cm)
        plt.colorbar(c)
        ax.set_aspect('equal')

    if(mark_ballooning_angle):
        #determine the center of the innermost flux surface
        x_center = np.mean(x[0,:])
        y_center = np.mean(y[0,:])
        radius = np.min([(np.max(x) - np.min(x))/2.0,
                          (np.max(y) - np.min(y))/2.0])
        if(detect_ballooning_angle):
            # find the center of mass of the matrix
            xx, yy = np.meshgrid(np.arange(data.shape[0]), np.arange(data.shape[1]), indexing="ij", sparse=False)
            com_ix = int(np.round(np.average(xx, weights=np.abs(data))))
            com_iy = int(np.round(np.average(yy, weights=np.abs(data))))
            com_x = x[com_ix,com_iy]
            com_y = y[com_ix,com_iy]
            ballooning_angle = np.arctan2(com_y-y_center, com_x-x_center)

        if(ax is not None):
            ax.plot([x_center, x_center+radius*np.cos(0.0)],
                    [y_center, y_center+radius*np.sin(0.0)],
                    linestyle="solid",
                    marker=None,
                    color="black",
                    linewidth=1,
                    solid_capstyle="round",
                )
            ax.plot([x_center, x_center+radius*np.cos(ballooning_angle)],
                    [y_center, y_center+radius*np.sin(ballooning_angle)],
                    linestyle="solid",
                    marker=None,
                    color="black",
                    linewidth=3,
                    solid_capstyle="round",
                )

    return ballooning_angle, data, x, y


def get_poloidal_slice_xy(f, dsetname, rhostar = 0.001, interpnum = 16,
                          imod = None, normalise=True, isp=0, ieiv=-1, r_boxsize=80,
                          delta_ix=0):
    # This function will extract from the data of a run the grids as well as
    # the perturbed potential.
    # The x,y coordinates are calculated from these, also the potential for
    # constant phi and plot poloidal slices of the real and imag parts

    if(f['/'].attrs['control.io_legacy'][0] == b'T'):
        raise NotImplementedError("For now, this function uses data in a form as produced with the the parameter io_legacy=.false.")

    spectral_radius = f['/'].attrs['control.spectral_radius'][0] == b'T'
    flux_tube = f['/'].attrs['control.flux_tube'][0] == b'T'
    n_s_grid = f['/'].attrs['grid.n_s_grid'][0]
    n_x_grid = f['/'].attrs['grid.n_x_grid'][0]
    nmod = f['/'].attrs['grid.nmod'][0]
    n_y_grid = nmod*2-1

    if(imod is None):
        imod = nmod-1
    imod = np.clip(imod, 0, nmod-1)

    print('n_s_grid:', str(n_s_grid))
    print('n_x_grid:', str(n_x_grid))
    print('nmod:', str(nmod))
    print('n_y_grid:', str(n_y_grid))
    print()

    if(flux_tube):
        # local
        kzeta = f['/grid/kzeta'].value
        kthrho = f['/grid/krho'].value
        kthnorm = f['/grid/kthnorm'].value
        kzeta = kthrho[imod,ix]/kthnorm

        # g_zeta_zeta = f['geom/g_zeta_zeta'].value
        # print("g_zeta_zeta.shape raw:", str(g_zeta_zeta.shape))
        # #promote this array to a 3d matrix.
        # # FIXME there must be a better way to do this in numpy...
        # g_zeta_zeta = np.expand_dims(g_zeta_zeta, axis=1)
        # g_zeta_zeta = np.repeat(g_zeta_zeta, repeats=n_x_grid, axis=1)
        # g_zeta_zeta = np.expand_dims(g_zeta_zeta, axis=2)
        # g_zeta_zeta = np.repeat(g_zeta_zeta, repeats=n_y_grid, axis=2)
        # print("g_zeta_zeta.shape:", str(g_zeta_zeta.shape))

        inv_asp_ratio = f["/geom/eps"].value
        inv_asp_ratio = np.pad(inv_asp_ratio, (0,n_x_grid - inv_asp_ratio.size),mode='edge')
        inv_asp_ratio += np.linspace(-r_boxsize/2*rhostar,+r_boxsize/2*rhostar,inv_asp_ratio.size)
        r = inv_asp_ratio.reshape([1,n_x_grid])
        print("r.shape raw:", str(r.shape))
        # +1 to close the gap. otherwise one pizza slice is missing from the poloidal slice
        r = np.repeat(r, repeats=n_s_grid+1, axis=0)

        #n = kzeta/(2*np.pi*rhostar)
    else:
        # global

        kzeta = f['/grid/kzeta'].value
        print('kzeta.shape:', kzeta.shape)
        if(nmod == 1):
            kzeta_min = np.min(kzeta[0,:])
            n_y_grid = 2*(nmod+1) - 1
        else:
            kzeta_min = np.min(kzeta[1,:])
            n_y_grid = 2*nmod - 1
        boxsize_zeta = 2*np.pi/kzeta_min
        zeta = np.linspace(-boxsize_zeta/2.0, +boxsize_zeta/2.0, num=n_y_grid, endpoint=True)
        
        if(rhostar is None):
            rhostar = f['/'].attrs['spcgeneral.rhostar'][0]
        #kzeta = np.tile(kzeta, [1, n_x_grid])

        r = f['/grid/xgr'].value
        r = r.reshape([1,n_x_grid])
        # r-coordinate does not depend on s, so simply repeat the vector n_s_grid times.
        r = np.repeat(r, repeats=n_s_grid+1, axis=0)
        # +1 to close the gap. otherwise one pizza slice is missing from the poloidal slice


    R = f['/geom/R'].value.reshape((n_x_grid, n_s_grid))
    Z = f['/geom/Z'].value.reshape((n_x_grid, n_s_grid))
    R = np.pad(R,((0,0),(0,1)), mode='wrap')
    Z = np.pad(Z,((0,0),(0,1)), mode='wrap')
              

    #print("kzeta.shape:", str(kzeta.shape))
    print("r.shape:", str(r.shape))

    q = f['/geom/q'].value
    print("q.shape raw:", str(q.shape))
    q = np.pad(q, (0,n_x_grid - q.size),mode='edge')
    print("q.shape raw:", str(q.shape))
    q = q.reshape([1,n_x_grid])
    # +1 to close the gap. otherwise one pizza slice is missing from the poloidal slice
    q = np.repeat(q, repeats=n_s_grid+1, axis=0)
    print("q.shape:", str(q.shape))
    
    #Also read the data, in this
    # case it is the square of the value of the perturbed potential.
    #data, real and imag part, has shape n_s_grid,n_x_grid

    if(dsetname.startswith("/diagnostic/diagnos_mode_struct")):
        if(spectral_radius):
            mode_label = f['/diagnostic/diagnos_grid/mode_label'].value
            kxrh = f['/grid/kxrh'];
            ixzero = np.argmin(np.abs(kxrh))
            radially_connected_mask = mode_label[:,imod] == mode_label[ixzero+delta_ix,imod]
            print(radially_connected_mask)

        data = f[dsetname+"_real"].value + 1j*f[dsetname+"_imag"].value
        binormal_axis = 2
        if(ieiv > data.shape[-1]-1):
            ieiv = data.shape[-1]-1
        if(data.ndim == 5):
            ###### PLOT A MOMENT
            data = data[:,:,:,isp,ieiv]
            data = np.pad(data,((0,0),(0,1),(0,0)), mode='wrap')
            #NOTE that radial and parallel dimension are interchanged, comparing fields and moments data!
            parallel_axis = 1
            radial_axis = 0
            if(spectral_radius):
                #now erase the radial modes which are not connected to
                #the one we have chosen through delta_ix
                data[radially_connected_mask==False,:,:] = 0.0
        else:
            ###### PLOT A FIELD
            data = data[:,:,:,ieiv]
            data = np.pad(data,((0,1),(0,0),(0,0)), mode='wrap')
            parallel_axis = 0
            radial_axis = 1
            if(spectral_radius):
                #now erase the radial modes which are not connected to
                #the one we have chosen through delta_ix
                data[:,radially_connected_mask==False,:] = 0.0

        if(normalise):
            argmax = np.argmax(np.abs(data))
            argmax_tupel = np.unravel_index(argmax, data.shape)
            data /= data[argmax_tupel]


        print("data.shape raw:", str(data.shape))
        # inverse fft the radial direction
        if(spectral_radius):
            data = scipy.fftpack.ifft(scipy.fftpack.ifftshift(data, axes=radial_axis), axis=radial_axis)

    elif(dsetname.startswith("/diagnostic/diagnos_fields")):
        #spc is useless because it does not have s.
        # the realspace datasets like
        #dsetname = "/diagnostic/diagnos_fields/Poten%08d" % ntime
        # are 3d.

        parallel_axis = 0
        radial_axis = 1
        binormal_axis = 2
            
        data = f[dsetname].value

        if(normalise):
            argmax = np.argmax(np.abs(data))
            argmax_tupel = np.unravel_index(argmax, data.shape)
            data /= data[argmax_tupel]

        print("data.shape raw:", str(data.shape))
        data = scipy.fftpack.fft(data, axis=binormal_axis)[:,:,0:nmod]
        
        #[:,imod,-1]
        # sum over y
        #data = np.sum(data, axis = 1)

    print("data.shape:", str(data.shape), data.dtype)

    r_i = r
    q_i = q
    R_i = R
    Z_i = Z
    data_i = data

    interpnum = max(1,interpnum)
    # interpolate with respect to s
    if(interpnum > 1):
        r_i = scipy.signal.resample(r_i,(n_s_grid+1)*interpnum,axis=0);
        q_i = scipy.signal.resample(q_i,(n_s_grid+1)*interpnum,axis=0);

        R_i = scipy.signal.resample(R_i,(n_s_grid+1)*interpnum,axis=1);
        Z_i = scipy.signal.resample(Z_i,(n_s_grid+1)*interpnum,axis=1);
        
        #s_i_fft = scipy.signal.resample(s,(n_s_grid+1)*interpnum,axis=0);
        #downsample_factor = 2
        #s_i_poly = scipy.signal.resample_poly(s, interpnum*downsample_factor, downsample_factor)

        #kzeta_i = scipy.signal.resample(kzeta,(n_s_grid+1)*interpnum,axis=0);
        data_i = scipy.signal.resample(data_i,(n_s_grid+1)*interpnum,axis=parallel_axis);
        
    #print("kzeta_i.shape:", str(kzeta_i.shape))
    print("r_i.shape:", str(r_i.shape))
    print("q_i.shape:", str(q_i.shape))
    print("data_i.shape:", str(data_i.shape),data_i.dtype)

    # if(flux_tube and spectral_radius):
    #     xgrid = f['/grid/xphi'].value[0,:]
    # else:
    #     xgrid = f['/grid/xgr'].value
    sgrid = f['/diagnostic/diagnos_grid/sgrid'].value
    #print(xgrid.shape)
    print(sgrid.shape)
    #print(xgrid)
    print(sgrid)
    # c = ax.contourf(xgrid, sgrid, data2d)
    # plt.colorbar(c)
    # ax.set(ylabel='s [dimensionless]', xlabel=r'x [$\rho_{ref}$]')
    # ax.legend(loc='best')

    # plt.tight_layout()
    # plt.savefig(filename)
    # plt.close()


    # #for debugging, some plots
    # plt.figure()
    # plt.plot(s[:,0], marker='.')
    # plt.plot(s_i_fft[:,0], marker = 'o')
    # plt.plot(s_i_poly[:,0], marker = 'x')
    # plt.plot(s_i_interp[:,0], marker = '*')
    # plt.show()
    # plt.close()

    # plt.figure()
    # plt.plot(kzeta[:,0], marker='.')
    # plt.plot(kzeta_i[:,0], marker = 'o')
    # plt.show()a
    # plt.close()

    # plt.figure()
    # plt.plot(data[:,0], marker='x')
    # plt.plot(data_i[:,0], marker = 'o')
    # plt.show()
    # plt.close()

    try:
        if(flux_tube):
            theta = f['/geom/poloidal_angle'].value
            theta = np.pad(theta,((0,1)), mode='wrap')
            theta = theta.reshape([-1,1])
            theta = np.repeat(theta, repeats=n_x_grid, axis=1)
        else:
            theta = f['/geom/poloidal_angle'].value.reshape([n_x_grid,n_s_grid]).transpose()
            theta = np.pad(theta,((0,1),(0,0)), mode='wrap')
        print('-theta.shape:', theta.shape)
        theta[-1,:] += 2*np.pi
        print('=theta.shape:', theta.shape)

        print("theta.shape before:", str(theta.shape))
        #print(theta)
        if(interpnum > 1):
            #theta_i = scipy.signal.resample(theta_i,(n_s_grid+1)*interpnum,axis=0);
            tmp_x = np.linspace(0,1,num=(n_s_grid+1)*interpnum)
            if(False):
                interp_func = scipy.interpolate.interp1d(np.linspace(0,1,num=(n_s_grid+1)), theta,fill_value='extrapolate' )
                theta_i = interp_func(tmp_x)
            else:
                theta_i = np.zeros(((n_s_grid+1)*interpnum, n_x_grid))
                for ix in range(theta_i.shape[1]):
                    interp_func = scipy.interpolate.interp1d(np.linspace(0,1,num=(n_s_grid+1)), theta[:,ix],fill_value='extrapolate' )
                    theta_i[:,ix] = interp_func(tmp_x)
        else:
            theta_i = theta
            print("theta is not interpolated")
        print("theta was read from /geom/poloidal_angle")
        print('*theta_i.shape:', theta_i.shape)
    except Exception as err:
        # if we cannot read the poloidal angle from the data, we
        # construct it from the analytical expression for circ
        # geometry.
        print(err)
        import traceback
        traceback.print_exc()

        s = f['/geom/s_grid'].value.reshape([n_s_grid,1])
        # +1 to close the gap. otherwise one pizza slice is missing from the poloidal slice
        s = np.pad(s,((0,1),(0,0)), mode='wrap')
        s = np.repeat(s, repeats=n_x_grid, axis=1)
        print("s.shape:", str(s.shape))
        s_i = s
        if(interpnum > 1):
            #s_i = scipy.signal.resample(s_i,(n_s_grid+1)*interpnum,axis=0);
            interp_func = scipy.interpolate.interp1d(np.linspace(0,1,num=(n_s_grid+1)), s_i[:,0],fill_value='extrapolate' )
            tmp_x = np.linspace(0,1,num=(n_s_grid+1)*interpnum)
            s_i = interp_func(tmp_x).reshape([-1,1])
            s_i = np.repeat(s_i, repeats=n_x_grid, axis=1)
        
        theta_i = np.zeros(r_i.shape)
        for ii in range(theta_i.shape[0]):
          for jj in range(theta_i.shape[1]):
            # Compute the angle from the s-coordinate.
            # Note: This assumes circular geometry.
            #theta_i(ii, jj) = fzero(@(x) (x + r_i(ii, jj)*sin(x) - 2*pi*s_i(ii, jj)), 2*pi*s_i(ii, jj));
            func = lambda x: x + r_i[ii, jj]*np.sin(x) - 2*np.pi*s_i[ii, jj]
            a = 2*np.pi*s_i[ii, jj]-r_i[ii, jj]
            b = 2*np.pi*s_i[ii, jj]+r_i[ii, jj]
            if(not (np.sign(func(a)) != np.sign(func(b)))):
                # cannot happen if Im not wrong
                print("func bracketing:" + str((func(a), func(b))))
            theta_i[ii, jj] = scipy.optimize.brentq(func,a,b)
        print("theta was computed from circ geometry formulas!")

    # theta2 = f['/geom/poloidal_angle'].value.reshape([-1,1])
    # theta2 = np.append(theta2,[theta2[0,:]],axis=0)
    # theta2 = np.tile(theta2, [1,n_x_grid])
    # #print(theta2.shape)
    # theta2 = scipy.signal.resample(theta2,n_s_grid*interpnum,axis=0);
    # #print("+"*20)
    # #print(theta-theta2)
    # theta = theta2
    print("theta_i.shape:", str(theta_i.shape))

    # Compute x- and y- cartesian coordinates for each point.
    x = r_i*np.cos(theta_i);
    y = r_i*np.sin(theta_i);

    print("x.shape:", str(x.shape))
    print("y.shape:", str(y.shape))

    # compute the position space value from the complex fourier coefficient
    #print("kzeta:" + str(kzeta_i))

    #kzeta_i = np.repeat(kzeta_i,repeats=n_y_grid, axis=2)

    #\int_{von phi=phi(s,psi,zeta=0); entlang s=const, psi=const}^{phi=0} dzeta = \int sqrt(g_zeta_zeta(s, psi)dzeta *dzeta)
    
    print("data_i.shape:", str(data_i.shape), data_i.dtype)
    if(imod < 0 or imod >= nmod):
        #sum over y
        data_i = np.sum(data_i, axis=binormal_axis)
    else:
        data_i = data_i[:,:,imod]
        # the zeta grid, (for circ geometry only!)
        zeta = (q_i/np.pi) * np.arctan(np.sqrt((1-r_i)/(1+r_i))*np.tan(theta_i/2));
        print("kzeta.shape:", str(kzeta.shape))
        print("zeta.shape:", str(zeta.shape))
        exp_factor = np.exp(1j*kzeta[imod]*zeta/rhostar)
        #exp_factor = 1
        print("exp_factor.shape:",exp_factor.shape)
        #if(dsetname.startswith("/diagnostic/diagnos_mode_struct")):
            #data_i = data_i*np.repeat(exp_factor,repeats=n_y_grid, axis=2)
            #data_i = data_i*exp_factor
            
    print("data_i.shape:", str(data_i.shape))
    # if(radial_axis == 0):
    #     data_i = data_i.transpose()
    #return (x,y,data_i)

    if(radial_axis == 1):
        data_i = data_i.transpose()
    data_i *= exp_factor.transpose()

    print("data_i.shape:", str(data_i.shape), data_i.dtype)
    
    return (R_i,Z_i,data_i)



# def plot_poloidal_slice(f, dsetname):

#     #data = f[dsetname].value
#     data = f[dsetname+"_real"].value + 1j*f[dsetname+"_imag"].value
#     print("Data shape:" + str(data.shape))
#     rhostar = 0.001

#     n_s_grid = f['/'].attrs['grid.n_s_grid'][0]
#     n_x_grid = f['/'].attrs['grid.nx'][0]
#     nmod = f['/'].attrs['grid.nmod'][0]

#     print(n_s_grid)
#     print(n_x_grid)
#     print(nmod)

#     inv_asp_ratio = f["/geom/eps"].value
#     inv_asp_ratio = np.pad(inv_asp_ratio, (0,n_x_grid - inv_asp_ratio.size),mode='edge')
#     print("inv_asp_ratio grid: " + str(inv_asp_ratio.shape))
#     inv_asp_ratio += np.linspace(-80/2*rhostar,+80/2*rhostar,inv_asp_ratio.size)

#     q = f["/geom/q"].value
#     q = np.pad(q, (0,n_s_grid - q.size),mode='edge')

#     sgrid = f['/geom/s_grid'].value

#     data = abs(data[:,:,int((nmod-1)/2 )])
#     #   # Get the real/imaginary part of the perturbed potential.
#     #   phi.rzc = reshape(parallel(1:n_s_grid*nx, 2*field + 0), n_s_grid, nx);
#     #   phi.izc = reshape(parallel(1:n_s_grid*nx, 2*field + 1), n_s_grid, nx);
#     #   # Build the complex potential, and calculate the absolute value.
#     #   zc  = complex(phi.rzc, phi.izc);
#     #   z   = abs(zc);

#     #   #kzeta = (2*pi/rhostar)*krho*(r(0,:))/(q(1,:));
#     #   #kzeta = repmat(kzeta,n_s_grid,1);

#     #   # For x(radial)-direction we use at first a simple index.
#     #   phi.ind = reshape(floor((0:(n_s_grid*nx-1))/n_s_grid), n_s_grid, nx);

#     #   phi.theta = reshape(geom.poloidal_angle, n_s_grid, nx);
#     poloidal_angle = f['/geom/poloidal_angle'].value
#     if(poloidal_angle.size == n_x_grid*n_s_grid):
#         poloidal_angle = poloidal_angle.reshape((n_s_grid,n_x_grid))
#     elif(poloidal_angle.size == n_s_grid):
#         poloidal_angle = poloidal_angle.reshape((-1,1))
#         poloidal_angle = np.tile(poloidal_angle, (1,n_x_grid))
#     else:
#         raise SystemExit('Unexpected shape of poloidal_angle array: ' + str(poloidal_angle.shape))

#     xg = inv_asp_ratio
#     xg = np.tile(xg,(n_s_grid,1))
#     print("xg:" + str(xg.shape))

#     #   Compute the cartesian x and y coordinates, i.e. the coordinates in the poloidal cut.
#     x = xg*np.cos(poloidal_angle)
#     y = xg*np.sin(poloidal_angle)

#     #zeta = (q/np.pi) * np.arctan(np.sqrt((1-inv_asp_ratio)/(1+inv_asp_ratio)) * np.tan(thetag/2));
#     # z2 = data * np.abs(np.sin(kzeta * zeta));
#     # data_rotated = np.real(data * exp(1j * kzeta*zeta/rhostar))
#     import matplotlib.tri as tri
#     triang = tri.Triangulation(x.flatten(), y.flatten())

#     plt.figure()
#     plt.gca().set_aspect('equal')
#     # plot the triangulation lines
#     #plt.triplot(triang, lw=0.5, color='black')
#     #plt.scatter(x.flatten(), y.flatten())
#     plt.contourf(x,y,data[:,:,0])


#     # levels = np.arange(0., 1., 0.025)
#     # cmap = cm.get_cmap(name='terrain', lut=None)
#     # plt.tricontourf(tri_refi, z_test_refi, levels=levels, cmap=cmap)
#     # plt.tricontour(tri_refi, z_test_refi, levels=levels,
#     #                colors=['0.25', '0.5', '0.5', '0.5', '0.5'],
#     #                linewidths=[1.0, 0.5, 0.5, 0.5, 0.5])

#     plt.title(dsetname)
#     #plt.show()
#     plt.savefig(dsetname.replace('/','_')+'.png')
#     plt.close()


if(__name__ == '__main__'):

    f = h5py.File("gkwdata.h5", "r")

    # nmod = f['/'].attrs['grid.nmod'][0]
    p = ModeStructurePlot("gkwdata.h5", parse_cmd_line_args=True)
    fig, ax = plt.subplots()
    p.plot(ax)
    fig.savefig("poloidal_slice_mode_struct_%s_%s.png" % (p.arguments["dsetname"].get(), p.arguments["part"].get()))
    plt.close()

    exit(0)
    
    # for imod in range(1,nmod+1):
    #     print("Plotting %s (_real and _imag)..." % dsetname)
    #     #plot_poloidal_slice(f, dsetname)
    #     import math
    #     plot_poloidal_slice(f, dsetname, interpnum=(math.floor(imod/16.0)+1)*16,imod=imod)
    
    ntime_end = f['/grid/time'].value.shape[1]
    # #for ntime in range(1,f['/grid/time'].shape[-1]):
    # for ntime in (150,):
    ntime = 0
    while(True):
        ntime += 1
        if(ntime > 1 and ntime < ntime_end):
            continue
        print("-"*40)

    
        try:
            dsetname = "/diagnostic/diagnos_fields/Poten%08d" % ntime
            print("Plotting %s..." % dsetname)
            fig, ax = plt.subplots()
            ballooning_angle, data, x, y = plot_poloidal_slice(ax,f, dsetname, normalise=False)
            
            fig.savefig(filename)
            plt.close()
        except KeyError:
            print()
            print("There are no further datasets.")
            exit(0)
