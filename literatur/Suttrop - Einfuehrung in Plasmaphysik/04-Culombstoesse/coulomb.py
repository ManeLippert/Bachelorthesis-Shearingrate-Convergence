#!/usr/bin/env python3

"""
  Electron undergoing a Coulomb collision
"""

import numpy as np
import scipy.integrate as integrate

# natural constants
e = 1.602e-19   # [As] elementary charge
epsilon0 = 8.854e-12   # [As/Vm]  dielectric constant
m_e = 9.11e-31   # [kg] electron mass

# global parameters
E0 = 2*e   # [eV]   incident kinetic energy
vz0 = np.sqrt(2*E0/m_e)   # [m/s] incident velocity

q1 = e          # [As] of scattering centre
q2 = -e         # [As] of scattered particle
r0perp = np.abs(q1*q2)/(8*np.pi*epsilon0*E0)  # [m]  impact parameter for 90 degrees deflection


# ---------------------------------------------------------------------
class Particle:
    '''
      follow one particle in electrostatic potential
      :dt:  time step [s]   (no default value)
      :m:   particle mass [kg]    (default: electron mass)
      :lambda_D:  Debye shielding length [m] (None if no shielding)
    '''
    def __init__(self, dt, m=m_e, lambda_D=None):
        self.m = m
        self.lambda_D = lambda_D
        self.dt = dt
        self.CoulombConst = q1*q2 / (4*np.pi*epsilon0)
        self.tot_weight = 0.0   # accumulated weight
        self.tot_dvz = 0.0      # accumulated
        self.r = np.array([])   # empty arrays to start with
        self.z = np.array([])

    def init(self, init_state):
        '''
        initialise particle position and velocity
        :init_state:  [r0, z0, vr0, vz0]
        '''
        self.init_state = init_state
        self.state = self.init_state
        self.r = np.array([init_state[0]])
        self.z = np.array([init_state[1]])

    def dvz_acc(self, r0perp):
        '''
         evaluate slowing down: accumulated dvz-dvz0
        '''
        weight = (self.init_state[0]/r0perp)**2 # cyl-sym around z -> weight = (r0/r0perp)^2 
        dvz = self.state[3] - self.init_state[3]  # vz - vz0
        self.tot_weight += weight
        self.tot_dvz += (dvz*weight)
        return self.tot_dvz  # accumulated dvz

    def energy(self):
        v2 = self.state[2]**2 + self.state[3]**2
        return self.m*v2/2
        
    def position(self):
        return (self.state[0], self.state[1])

    def force(self, r, z):
        '''
           return force (Fr, Fz) [N] on particle
        '''
        return np.array([0.0, 0.0])   # no force

    def dstate_dt(self, t, state):
        ''' 
           return derivative of a given state
        '''
        a = self.force(state[0], state[1]) / self.m
        return np.array([state[2], state[3], a[0], a[1]])

    def step(self):
        '''
           follow particle for one time step (self.dt)
        '''
        sol = integrate.solve_ivp(self.dstate_dt, [0, self.dt], self.state, \
                                  t_eval = [0, self.dt], \
                                  method='Radau', atol=1e-9, rtol=1e-8,
                                  max_step=self.dt/10)        
        nt = sol.y.shape[1]
        self.state = sol.y[:,nt-1]
        self.r = np.hstack((self.r, np.array([self.state[0]])))
        self.z = np.hstack((self.z, np.array([self.state[1]])))

    def nstep(self):
        '''
           follow particle and record states at a timebase
        '''
        sol = integrate.solve_ivp(self.dstate_dt, [0, self.dt], self.state, \
                                  method='Radau', atol=1e-9, rtol=1e-8,
                                  max_step=self.dt/10)
        # print('solve_ivp status: ', sol.status, ' nfev: ', sol.nfev)
        self.state = sol.y[:,-1]
        self.r = np.hstack((self.r, sol.y[0, 1:]))
        self.z = np.hstack((self.z, sol.y[1, 1:]))
        return sol.t, sol.y


# ---------------------------------------------------------------------
class Particle_Coulomb (Particle):

    def force(self, r, z):
        '''
           return Coulomb force (Fr, Fz) [N] on particle
        '''
        d2 = r**2 + z**2
        d = np.sqrt(d2)
        if self.lambda_D is not None and d > self.lambda_D:
          F0 = 0.0   # cut
        else:  # Coulomb force
          F0 = self.CoulombConst / d2    # Coulomb force
        return F0*np.array([r/d, z/d])

# ---------------------------------------------------------------------
class Particle_shielded (Particle):

    def force(self, r, z):
        '''
           return force (Fr, Fz) [N] on particle in Debye-shielded potential
        '''
        d2 = r**2 + z**2
        d = np.sqrt(d2)
        F0 = self.CoulombConst * np.exp(-d/self.lambda_D) \
            * ( 1.0/d2 + 1.0/self.lambda_D*d)    # shielded Coulomb force
        return F0*np.array([r/d, z/d])

# -------------------- main program -------------------------
if __name__ == '__main__':
    import matplotlib
    import matplotlib.pyplot as plt
    import matplotlib.animation as animation

    zd = 40*r0perp  # domain size
    rd = 30*r0perp
    tmax = 2*zd/vz0
    lambda_D = 5.0*r0perp  # Debye length

    fig = plt.figure()
    ax = fig.add_subplot(111, aspect='equal', autoscale_on=False,
                         xlim=(-zd, zd), ylim=(-rd, rd))
    ax.plot([-zd,0],[r0perp,r0perp],linestyle='dotted',color='grey')
    ax.plot([-zd,zd],[0,0],linestyle='solid',color='black',linewidth=0.2)
    ax.plot([0,0],[-rd,rd],linestyle='solid',color='black',linewidth=0.2)
    DebyeSphere = plt.Circle((0,0),radius=lambda_D,color='#D0D0D0',linestyle='dotted')
    ax.add_artist(DebyeSphere)
    
    rf = 1.0  # impact parameter in units of r0perp
    init_state = [rf*r0perp, -zd, 0.0, vz0]

    if True:
        p1 = Particle_Coulomb  (tmax, m=m_e, lambda_D=None) # Coulomb 
        p1.init(init_state); p1t, p1sol = p1.nstep();
        ax.plot(p1.z, p1.r, color='black')

    if True:
        p2 = Particle_Coulomb  (tmax, m=m_e, lambda_D=lambda_D) # Coulomb 
        p2.init(init_state); p2t, p2sol = p2.nstep();
        ax.plot(p2.z, p2.r, color='red')

    if True:
        p3 = Particle_shielded  (tmax, m=m_e, lambda_D=lambda_D) # Coulomb 
        p3.init(init_state); p3t, p3sol = p3.nstep();
        ax.plot(p3.z, p3.r, color='blue')

    plt.show()

# -------------------- end of coulomb.py --------------------------------
