#!/usr/bin/env python3

'''
  Simulation of Landau damping
  using the most simple 1D electrostatic PIC model (one spatial, one velocity dimension)
  after pic.py by Philip Mocs
  https://medium.com/swlh/create-your-own-plasma-pic-simulation-with-python-39145c66578b
'''

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve

class Electrostatic_1d():

    def __init__(self, Nx, boxsize):
        '''
          electrostatic field only
          :Nx:       number of grid points, preferrably a power of 2
          :boxsize:  length of computational domain
        '''
        self.Nx = Nx
        self.boxsize = boxsize
        self.grid = np.linspace(0,boxsize,Nx)

	# Construct matrix G to compute Gradient  (1st derivative)
        dx = boxsize/Nx
        e = np.ones(Nx)
        diags = np.array([-1,1])
        vals  = np.vstack((-e,e))
        Gmtx = sp.spdiags(vals, diags, Nx, Nx);
        Gmtx = sp.lil_matrix(Gmtx)
        Gmtx[0,Nx-1] = -1
        Gmtx[Nx-1,0] = 1
        Gmtx /= (2*dx)
        self.Gmtx = sp.csr_matrix(Gmtx)

        # Construct matrix L to compute Laplacian (2nd derivative)
        diags = np.array([-1,0,1])
        vals  = np.vstack((e,-2*e,e))
        Lmtx = sp.spdiags(vals, diags, Nx, Nx);
        Lmtx = sp.lil_matrix(Lmtx)
        Lmtx[0,Nx-1] = 1
        Lmtx[Nx-1,0] = 1
        Lmtx /= dx**2
        self.Lmtx = sp.csr_matrix(Lmtx)

    def weights(self, pos):
        '''
         weights for particles at positions 'pos'
        '''
        dx         = self.boxsize / self.Nx
        j          = np.floor(pos/dx).astype(int)
        jp1        = j+1
        weight_j   = ( jp1*dx - pos  )/dx
        weight_jp1 = ( pos    - j*dx )/dx
        j          = np.mod(j, self.Nx)
        jp1        = np.mod(jp1, self.Nx)   # periodic BC
        return j, jp1, weight_j, weight_jp1
        
    def clear_charge(self):
        '''
          clear charge density and total accumulated charge
        '''
        self.rho = np.zeros(Nx)
        self.rho_tot = 0.0

    def add_charge(self, pos, rho):
        '''
          add charges at positions 'pos' with magnitude 'rho'
        '''
        j, jp1, weight_j, weight_jp1 = self.weights(pos)
        n  = np.bincount(j,   weights=weight_j,   minlength=self.Nx)
        n += np.bincount(jp1, weights=weight_jp1, minlength=self.Nx)
        N = pos.shape[0]
        dx = self.boxsize / self.Nx
        self.rho += n * rho * self.boxsize / N / dx
        self.rho_tot += rho

    def update(self):
        '''
          solve for potential phi and electrical field E
        '''
        # Solve Poisson's Equation: laplacian(phi) = n-n0
        self.phi = spsolve(self.Lmtx, self.rho - self.rho_tot, \
                           permc_spec="MMD_AT_PLUS_A")
        # Apply Derivative to get the Electric field
        self.E = - self.Gmtx @ self.phi

    def E_field(self, pos):
        '''
          E field at positions 'pos'
        '''
        j, jp1, weight_j, weight_jp1 = self.weights(pos)
        return weight_j * self.E[j] + weight_jp1 * self.E[jp1]


# ------------------------------------------------------------
        
class Particles():

    def __init__(self, field, pos, vel, n0):
        '''
          particle species
        '''
        self.field = field
        self.pos = pos
        self.vel = vel
        self.n0 = n0

    def move(self, dt):
        '''
          move particles for time step 'dt' and apply bc
        '''
        self.pos += self.vel * dt
        # periodic boundary conditions
        self.pos = np.mod(self.pos, self.field.boxsize)

    def add_charge(self):
        self.field.add_charge(self.pos, self.n0)
        
    def accelerate(self, dt):
        '''
          accelerate particles 
        '''
        acc = - self.field.E_field(self.pos)
        self.vel += acc * dt


# ---------------------- "main" program ----------------------------------

if True:
    import matplotlib.pyplot as plt

    N         = 80000    # number of particles
    n0        = 1        # electron number density
    Nx        = 32       # number of grid points
    boxsize   = 6.283185307  # periodic domain [0,boxsize]

# -------------------- Initial Conditions --------------------------------
    np.random.seed(42)            # set the random number generator seed    

    pos  = np.random.rand(N) * boxsize   # random particle positions
    # pos  = np.linspace(0.0, boxsize * (N-1.0)/N, N)  # uniform particle spacing

    if True:  # position perturbation 
        mode = 1
        x1 = 0.1
        pos += x1 * np.cos( 2.0*np.pi* mode* pos / boxsize)  # n=1 preset for initial phi

    v0 = -1.3    # offset velocity
    vth  = 0.4    # thermal velocity 
    vel = v0 + vth * np.random.randn(N)  # random velocity

    field = Electrostatic_1d(Nx, boxsize)
    spec1 = Particles(field, pos, vel, n0)

# ------------------ diagnostic set-up ------------------------
    vmin = np.min(spec1.vel)
    vmax = np.max(spec1.vel)

    nvbins = 50   #  #bins for distribution function
    vbinedges = np.linspace(vmin, vmax, nvbins+1)
    vbincentr = 0.5*(vbinedges[1:]+vbinedges[:-1])

    nmarkers = 1000   # number of marker particles
    markerint = 5   # number of most recent iterations to plot
    markervel = -0.5 + np.random.rand(nmarkers)    # random initial velocities 
    markerpos = boxsize * np.random.rand(nmarkers) # random initial positions
    markers = Particles(field, markerpos, markervel, nmarkers/N)

    f_intrlv = 20    # phase space interleave factor (plot every n-th data point)

    # analytical growth rate
    k2 = (2*np.pi/boxsize)**2   # =1, here
    gamma = - np.sqrt(np.pi / 8) / k2**1.5 / vth**3 * np.exp(-1/(k2*vth**2) - 1.5)
    
# ............ simulation ....................

    plotRealTime = True
    Nt = 256  # number of timesteps
    t  = 0    # current time of the simulation
    dt = 0.1  # timestep

    # prepare figure
    fig = plt.figure(figsize=(14,8), dpi=80)

    # initial charge density and potential 
    field.clear_charge()
    spec1.add_charge()
    # markers are passive, they do not add to the charge density
    field.update()
    
    # simulation main loop
    for i in range(Nt):
        spec1.accelerate (dt/2.0)  # (1/2) kick
        markers.accelerate (dt/2.0)

        # drift (and apply periodic boundary conditions)
        spec1.move(dt)
        markers.move (dt)

        # update field
        field.clear_charge()
        spec1.add_charge()
        field.update()

        spec1.accelerate (dt/2.0)  # (1/2) kick
        markers.accelerate (dt/2.0)
        t += dt

        # if i==Nt//2:  # PIC echo!
        #    spec1.vel = -spec1.vel
        
        # diagnostics         
        vhist = np.histogram(spec1.vel, vbinedges)  # velocity distribution
        phispec = np.abs(np.fft.fft(field.phi))  # mode (k) spectral amplitude
        phi1 = phispec[1]     # fundamental: n=1, k=2*pi/boxsize

        # accumulate time evolution
        if i==0:
            timebase = np.array([t])
            phi_t = np.array([phi1])
            vhist0 = vhist[0]
        else:
            timebase = np.hstack( (timebase, np.array([t])) )
            phi_t = np.hstack( (phi_t, np.array([phi1])) )
            markerpos_t = np.vstack( (markerpos_t, np.reshape(markers.pos, (1,-1))) )
            markervel_t = np.vstack( (markervel_t, np.reshape(markers.vel, (1,-1))) )

        if i==0:
            markerpos_t =  np.reshape(markers.pos, (1,-1))
            markervel_t =  np.reshape(markers.vel, (1,-1))
        else:
            markerpos_t = np.vstack( (markerpos_t, np.reshape(markers.pos, (1,-1))) )
            markervel_t = np.vstack( (markervel_t, np.reshape(markers.vel, (1,-1))) )
            
        # plot in real time
        if plotRealTime or (i == Nt-1) or (i < 33):
            ax1 = fig.add_subplot(2,2,1, xlim=(0,boxsize), ylim=(-6,6))
            ax1.cla()
            plt.title('t = %f' % t)
            ax1.scatter(spec1.pos[::f_intrlv], spec1.vel[::f_intrlv], \
                        marker='.', s=0.4, color='blue')
            mi = min(i, markerint)
            ax1.scatter(markerpos_t[-mi:,:], markervel_t[-mi:,:], \
                        s=0.2, color='red')
            plt.ylabel('v')
            
            ax2 = fig.add_subplot(2,2,3, xlim=(0,boxsize))
            ax2.cla()
            ax2.plot(field.grid, field.phi)
            plt.xlabel('x')
            plt.ylabel('phi')

            ax3 = fig.add_subplot(2,2,2, xlim=(0,boxsize))
            ax3.cla()
            ax3.semilogy(timebase, phi_t)
            plt.xlabel('time')
            plt.ylabel('phi')

            ax4 = fig.add_subplot(2,2,4, xlim=(0,boxsize))
            ax4.cla()
            ax4.semilogy(vbincentr, vhist0, color='grey')
            ax4.semilogy(vbincentr, vhist[0], color='red')
            plt.xlabel('v')
            plt.ylabel('f(v)')

            plt.pause(0.001)
        else:
            print("step: %d   t = %8.3f" % (i+1, t))

    # Save figure
    # plt.savefig('pic.png',dpi=240)
    plt.show()


	

