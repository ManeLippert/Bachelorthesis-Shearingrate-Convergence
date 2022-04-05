#!/usr/bin/env python
# animate a (possibly gyrating, possibly damped) particle in an oscillating electrical field

import numpy as np
import scipy.integrate
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation


epsilon0 = 8.854e-12   # [As/Vm]  free-space permissivity
mu0 = 4*np.pi*1e-7     # [Vs/Am] free-space permeability
c0 = 1.0/np.sqrt(epsilon0 * mu0)   # [m/s] free-space speed of light

# ------------------ plane wave RF field --------------
class RF_PlaneWave ():

  def __init__(self, Ex1, Bz0, omega):
    '''
      plane wave - oscillating E_x field and static homogeneous B_z field 
      E1 = [Ex1, 0, 0] * cos(omega*t)
      B0 = [0, 0, Bz0] 
    '''
    self.E1 = np.array([Ex1, 0.0, 0.0])
    self.B0 = np.array([0.0, 0.0, Bz0])
    self.B1 = np.array([0.0, 0.0, -Ex1/c0])  # from Faraday's law, integrated
    self.omega = omega
    self.ky = omega / c0   # wave vector in y direction
    print('RF wavelength: %f mm' % (1e3/self.ky))
    print('Wave magnetic field amplitude: %f T' % (np.abs(Ex1)/c0))
    if Bz0!=0:
      print('  B1 / B0 = %f' % (np.abs(Ex1/c0/Bz0)))

  def Evec(self, t, r):
    '''
      return electrical field 
      at time t and position vector r.
      r can be a vector, r=[x,y,z], 
      or an array r=[xyz, caseindex]
    '''
    E = self.E1 * np.cos(self.omega * t)
    if r.ndim==1:
      return E
    else:   # same at all positions
      return np.outer(E, np.ones(r.shape[1]))

  def Bvec(self, t, r):
    '''
     return magnetic induction vector
     at time t and position(s) r.
     r can be a vector, r=[x,y,z], 
     or an array r=[xyz, caseindex]
    '''
    B = self.B0 + self.B1*np.sin(self.omega*t)
    if r.ndim==1:
      return self.B0
    else:
      return np.outer(self.B0, np.ones(r.shape[1]))


# ============ damped particle, full orbit ===================

class FullOrbit_Damped():

  def __init__(self, fields, dt, m, q, nu_c):
    '''
     initialise particle
     :fields: class to compute electromagnetic fields
     :dt: time step
     :m: particle mass
     :q: particle charge
     :nu_c: damping frequency (effective 90 deg scattering, by which momentum is lost)
    '''
    self.fields = fields
    self.dt = dt
    self.m = m
    self.q = q
    self.nu_c = nu_c    

  def init(self, init_state):
    '''
       initialise particle state
    '''
    self.state = init_state
    self.solution = np.array([init_state]).reshape(-1,1)
    self.time = 0.0
    self.timebase = np.array([self.time])

  def dstate_dt(self, t, state):
    '''
       time derivative of state (full orbit)
    '''
    r = state[0:3]
    v = state[3:6]
    E = self.fields.Evec(t, r)
    B = self.fields.Bvec(t, r)
    dvdt = (self.q/self.m) * (E + np.cross(v, B)) - v*self.nu_c
    return np.concatenate([v, dvdt])

  def step(self):
    '''
       follow particle and record states at iteration time steps
    '''
    sol = scipy.integrate.solve_ivp(self.dstate_dt, [self.time, self.time+self.dt], \
              self.state, method='Radau', atol=1e-9, rtol=1e-8,
              max_step=self.dt/10)
    self.state = sol.y[:,-1]
    self.solution = np.hstack((self.solution, sol.y[:, 1:]))
    self.timebase = np.hstack((self.timebase, sol.t[1:]))
    self.time += self.dt
    return sol.t, sol.y

   

# ================= main program =============================
   
# constants
e = 1.602E-19   # [As]   elementary charge
m_e = 9.11E-31   # [kg]   electron mass
m_p = 1.67E-27   # [kg]   proton mass

# ------ simulation cases -----------

if True:     # RF heated plasma
   Ex1 = 500   # [V/m]  wave amplitude
   omega = 2.0*np.pi* 13.56e6  # [rad/s] technical RF frequency
   omega_c = 0   # unmagnetised plasma
   nu_c = 0.0*omega   # collision frequency
   Ekin0 = 0*e   # start with resting particle

else:        # magnetised, electron cyclotron-heated plasma
   Ex1 = 14000
   omega = 2.0*np.pi* 2.45e9  # [rad/s] technical microwave frequency
   omega_c = omega   # cyclotron resonance
   nu_c = 0*omega
   Ekin0 = 20*e   # initial kinetic energy (in gyromotion) [J]

# -------- particle parameters ---------
m = m_e
q = -e
col = 'green'

# ------- field parameters -------------
Bz0 = np.abs(omega_c*m/q)

# ------- simulation parameters --------
ncycles = 10
nframes = 100 * ncycles
tmax = ncycles * (2*np.pi)/np.abs(omega)
tstep = tmax/nframes    # time step for each frame

# initialisation 
xstart = e*Ex1/(m*omega**2)
vstart = q/e*np.sqrt(2*Ekin0/m)
if omega_c!=0:
   r_L = np.sqrt(2*Ekin0/m)/np.abs(omega_c)
else:
   r_L = 0   # infinity, rather
   
# the simulation
f = RF_PlaneWave(Ex1, Bz0, omega) # the EM field
p = FullOrbit_Damped(f, tstep, m, q, nu_c)
x0 = [xstart, -r_L * np.sign(vstart), 0.0]
v0 = [vstart, 0.0, 0.0]

# plot
matplotlib.rcParams['figure.figsize'] = [19.0,10.0]
matplotlib.rcParams['figure.subplot.left'] = 0.09
matplotlib.rcParams['figure.subplot.right'] = 0.975
matplotlib.rcParams['figure.subplot.bottom'] = 0.05
matplotlib.rcParams['figure.subplot.top'] = 0.98

fig = plt.figure()
ab = 1.1e3*(2*r_L+2*np.abs(xstart))  # plot boundary
ax = fig.add_subplot(1, 2, 1, aspect='equal', autoscale_on=False,
                     xlim=(-ab, ab), ylim=(-ab, ab))
plt.xlabel('x [mm]')
plt.ylabel('y [mm]')
ax.plot([-ab,ab], [0.0,0.0], color='#D0D0D0')
ax.plot([0.0,0.0], [-ab,ab], color='#D0D0D0')
if omega_c!=0:
   GyroOrbit = plt.Circle((0,0),radius=1e3*r_L, color='grey', linestyle='dotted', fill=False)
   ax.add_artist(GyroOrbit)
point1, = ax.plot([], [], linestyle='none', marker='o', color=col, markersize=20)
line1, = ax.plot([], [], color=col)

# electrical field vs. time
ax2 = fig.add_subplot(3, 2, 2, autoscale_on=False, \
                     xlim=(0, 1e6*tmax), ylim=(-1.1e-3*Ex1, 1.1e-3*Ex1))
plt.ylabel('E_x [kV/m]')
ax2.plot([0,1e6*tmax], [0,0], color='grey', linestyle='dotted')
t_est = np.linspace(0, tmax, nframes)
ax2.plot(1e6*t_est, 1e-3*Ex1*np.cos(t_est*omega), color='grey', linestyle='dotted')
t_hist = np.array([])
Exhist = np.array([])
line2, = ax2.plot([], [], color='black')

# x-velocity vs. time
vscale = np.max((np.abs(vstart), np.abs(xstart*omega)))
ax3 = fig.add_subplot(3, 2, 4, autoscale_on=False,
                      xlim=(0, 1e6*tmax), ylim=(-2e-6*np.abs(vscale), 2e-6*np.abs(vscale)))
plt.ylabel('v_x [10^6 m/s]')
ax3.plot([0,1e6*tmax], [-1e-6*vstart,-1e-6*vstart], color='grey', linestyle='dotted')
ax3.plot([0,1e6*tmax], [1e-6*vstart,1e-6*vstart], color='grey', linestyle='dotted')
ax3.plot([0,1e6*tmax], [-1e-6*xstart*omega,-1e-6*xstart*omega], color='grey', linestyle='dotted')
ax3.plot([0,1e6*tmax], [1e-6*xstart*omega,1e-6*xstart*omega], color='grey', linestyle='dotted')
line3, = ax3.plot([], [], color='black')

# x-position vs. time
ax4 = fig.add_subplot(3, 2, 6, autoscale_on=False,
                      xlim=(0, 1e6*tmax), ylim=(-ab, ab))
plt.xlabel('t [us]')
plt.ylabel('x [mm]')
ax4.plot([0,1e6*tmax], [1e3*xstart,1e3*xstart], color='grey', linestyle='dotted')
ax4.plot([0,1e6*tmax], [-1e3*xstart,-1e3*xstart], color='grey', linestyle='dotted')
line4, = ax4.plot([], [], color='black')


#------------------------------------------------------------
# set up animation

def animate(i):
    global p, t_hist, Exhist
    p.step(); pos=1e3*p.solution[:,-1]
    t = p.timebase[-1]
    t_hist=np.hstack((t_hist, np.array([t])))
    Evec = f.Evec(t, pos[:3])
    Exhist=np.hstack((Exhist, np.array([Evec[0]])))
    point1.set_data([pos[0]], [pos[1]]);
    line1.set_data(1e3*p.solution[0,:], 1e3*p.solution[1,:])
    line2.set_data(1e6*t_hist, 1e-3*Exhist)
    line3.set_data(1e6*p.timebase, 1e-6*p.solution[3,:])
    line4.set_data(1e6*p.timebase, 1e3*p.solution[0,:])
    return point1,line1,line2,line3,line4,

p.init( np.concatenate([x0, v0]) )
ani = animation.FuncAnimation(fig, animate, frames=nframes,
                              interval=40, blit=True, repeat=False)

plt.show()
