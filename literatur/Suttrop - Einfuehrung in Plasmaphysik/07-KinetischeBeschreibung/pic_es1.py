#!/usr/bin/env python3

'''
  most simple 1D electrostatic PIC model (one spatial, one velocity dimension)
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


# ------------------------------------------------------------

if True:
    import matplotlib.pyplot as plt

    Nx        = 256     # number of grid points
    N         = 4000    # number of particles
    boxsize   = 1000     # periodic domain [0,boxsize]
    n0        = 1       # electron number density
    A         = 1       # perturbation amplitude
    vth       = 0.0     # thermal velocity 

# -------------------- Initial Conditions --------------------------------
    np.random.seed(42)            # set the random number generator seed
    
    if False:  # two-particle test
        pos = boxsize * np.array([0.4,0.6])
        vel = np.array([0.0, 0.0])
        plotRealTime = True 
        Nt = 64  # number of timesteps

    if False:  # mode excitation test
        pos  = np.linspace(0.0, boxsize * (N-1.0)/N, N)
        pos = np.hstack((pos[3:], pos[0:3]))
        l = 15   # mode number
        vel  = 0.05 * np.cos(2*np.pi*l*pos/boxsize)
        Nt = 256  # number of timesteps
        plotRealTime = False
    
    if True:   # random particle pos, beam with Gaussian broadening
        pos  = np.random.rand(N) * boxsize
        vth  = 1.0    # thermal velocity 
        vel  = vth * np.random.randn(N)
        plotRealTime = False
        Nt = 512  # number of timesteps

    field = Electrostatic_1d(Nx, boxsize)
    spec1 = Particles(field, pos, vel, n0)

# ............ simulation ....................

    t  = 0   # current time of the simulation
    dt = 0.2   # timestep

    Nxh = int(Nx/2)
    Nth = int(Nt/2)
    
    # prepare figure
    fig = plt.figure(figsize=(13,8), dpi=80)

    field.clear_charge()
    spec1.add_charge()
    field.update()
    
    # Simulation Main Loop
    for i in range(Nt):
        spec1.accelerate (dt/2.0)  # (1/2) kick

        # drift (and apply periodic boundary conditions)
        spec1.move(dt)

        # update field
        field.clear_charge()
        spec1.add_charge()
        field.update()

        spec1.accelerate (dt/2.0)  # (1/2) kick
        t += dt

        psispec = np.fft.fft(field.rho)
        if i==0:
            partpos = spec1.pos
            timebase = np.array([t])
            firstspec = psispec
            history = field.rho
        else:
            partpos = np.vstack( (partpos, spec1.pos) )
            timebase = np.hstack( (timebase, np.array([t])) )
            history = np.vstack( (history, field.rho) )
        
        # plot in real time - color 1/2 particles blue, other half red
        if plotRealTime or (i == Nt-1) or (i < 30):
            ax1 = fig.add_subplot(2,2,1, xlim=(0,boxsize), ylim=(-6,6))
            ax1.cla()
            plt.title('t = %f' % t)
            ax1.scatter(spec1.pos, spec1.vel, s=0.8, color='blue')
            plt.ylabel('v')
            
            ax2 = fig.add_subplot(2,2,3, xlim=(0,boxsize))
            ax2.cla()
            ax2.plot(field.grid, field.phi)
            plt.xlabel('x')
            plt.ylabel('phi')

            ax3 = fig.add_subplot(2,2,2)
            ax3.cla()
            if i!=0:
                ax3.plot(timebase, partpos[:,0], color='blue')
                ax3.plot(timebase, partpos[:,1], color='red')
                
            ax4 = fig.add_subplot(2,2,4, ylim=(0,20))
            ax4.cla()

            if i!=0:
                ax4.plot(np.abs(firstspec[1:Nxh]), color='grey')

            ax4.plot(np.abs(psispec[1:Nxh]), color='red')
            plt.pause(0.001)
        else:
            print("step: %d   t = %8.3f" % (i+1, t))

    # Save figure
    # plt.savefig('pic.png',dpi=240)
    plt.show()

    # 2D FFT (frequency vs. wavenumber)
    histspec = np.fft.fft2(history)

    # compute x axis (wavenumber axis)
    ikmax = int(Nxh/2)
    dk = 2*np.pi/boxsize # wave number slot width
    kmax = dk*ikmax
    kax = np.linspace(0, kmax-dk, ikmax)
    kbase = np.linspace(0, kmax-dk, 100)
    
    # compute y-axis (frequency axis)
    df_fp = 2*np.pi/(dt*Nt)  # frequency slot width
    fmax_fp = 2.3   # maximum normalised frequency we want to plot
    fmin_fp = 0.5   # minimum normalised frequency we want to plot
    ifmin = int(np.floor(fmin_fp / df_fp))
    ifmax = int(np.ceil(fmax_fp / df_fp))
    freqax = np.linspace(ifmin*df_fp, ifmax*df_fp, ifmax-ifmin)

    # contour plot of spectrum
    plt.contourf(kax, freqax, np.abs(histspec[ifmin:ifmax,0:ikmax]))
    # omega_p = const line
    plt.plot([0, kmax-dk], [1.0, 1.0], color='red', linestyle='dotted')
    # analytical dispersion relation
    plt.plot(kbase, np.sqrt(0.5 + 0.5*np.sqrt(1 + 12*(kbase**2)*(vth**2))), \
             color='red')
    # Bohm-Gross dispersion relation
    plt.plot(kbase, np.sqrt(1 + 3*(kbase**2)*(vth**2)), \
             color='red', linestyle = 'dotted')
    
    plt.xlabel('k')
    plt.ylabel('omega / omega_p')
    plt.show()

	

