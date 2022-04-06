#!/usr/bin/env python3
# animate a sinusoidal magnetic trap

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation

import particle

# --------------------- mirror field -----------------
class SinusoidalTrap(particle.SimpleFields):

  '''
   magnetic double mirror (trap) along z-direction
   Bz = Bmax - (Bmax-Bmin)*cos(k*z), k*L=2*pi
   Br according to div(B)=0
   no E field
  '''
  def __init__(self, Bmin, Bmax, L):
    '''
      :Bmin: [T] minimum field
      :Bmax: [T] maximum field
      :L: [m] mirror length
    '''
    self.Bmid = (Bmax+Bmin)/2.0
    self.Ba = (Bmax-Bmin)/2.0
    self.k = 2*np.pi/L

  def Evec(self, t, r):
    '''
      here: Zero E
    '''
    return [0.0, 0.0, 0.0]

  def Bvec(self, t, r):
    '''
      returns B(r) (time-independent field)
    '''
    kr = self.k*r[2]
    Bz = self.Bmid - self.Ba*np.cos(kr)
    Br = -np.sqrt(r[0]**2 + r[1]**2) * self.k * self.Ba/2.0 * np.sin(kr)
    theta = np.arctan2(r[1], r[0])
    return np.array([Br*np.cos(theta), Br*np.sin(theta), Bz])

  def grad_B(self, t, r):
    '''
      returns grad |B|  (time-independent field gradient)
    '''
    kr = self.k*r[2]
    Basinkr = self.Ba*np.sin(kr)
    Bacoskr = self.Ba*np.cos(kr)
    Bz = self.Bmid - Bacoskr
    rad = np.sqrt(r[0]**2 + r[1]**2)
    Br = -rad * self.k * Basinkr/2.0
    B_ = np.sqrt(Br**2 + Bz**2)    # |B|
    print('x: ', r[0], ' y:', r[1], ' z:', r[2], 'Br: ', Br, ' Bz: ', Bz)
    dBdr = Br**2/rad/2/B_
    dBdz = self.k/2/B_ *Basinkr *(2*self.Bmid + Bacoskr*(kr**2/2 + 2.0))
    theta = np.arctan2(r[1], r[0])
    return np.array([dBdr*np.cos(theta), dBdr*np.sin(theta), dBdz])


# ---------------- main program ----------------------

e = 1.602E-19   # As   elementary charge
me = 9.11E-31   # kg   electron mass
mp = 1.67E-27   # kg   proton mass

Bmin = 0.3  # T  minimum B
Bmax = 0.6  # T  maximum B
L = 0.5 # m  length between maxima
f = SinusoidalTrap (Bmin, Bmax, L)

# particle
q = e   
m = mp
Ekin0 = 200*e   # kinetic energy [J]
# alpha = 0.25*np.pi 
alpha = 0.26*np.pi   # pitch angle tan(alpha) = v_perp / v_par

# initial conditions
v_0 = np.sqrt(2*Ekin0/m) # initial velocity
v_perp0 = v_0 * np.sin(alpha)
v_par0  = v_0 * np.cos(alpha)
omega_c = q*Bmin/m
r_L = v_perp0/np.abs(omega_c)
v0 = [np.sign(q)*v_perp0, 0.0, v_par0]

Bm = ((v_par0/v_perp0)**2 + 1)*Bmin  # mirror field
alpha_tpb = np.arctan(np.sqrt(Bmin/(Bmax-Bmin)))
print('trapped-passing boundary: alpha = ', alpha_tpb)

# simulation parameters
tmax = 10 * L / v_par0 
nframes = 2000
p = particle.FullOrbit(f, tmax/nframes, m, q)

# plot
matplotlib.rcParams['figure.figsize'] = [14.0,7.0]
matplotlib.rcParams['figure.subplot.left'] = 0.05
matplotlib.rcParams['figure.subplot.right'] = 0.98
matplotlib.rcParams['figure.subplot.bottom'] = 0.07
matplotlib.rcParams['figure.subplot.top'] = 0.98

fig = plt.figure()
zl = 0.75*L     # z scale
rl = 1.1e3*r_L  # r scale

# z-y trajectory
ax = fig.add_subplot(2, 1, 1, autoscale_on=False,
                     xlim=(-zl, zl), ylim=(-rl, rl))
plt.xlabel('z [m]')
plt.ylabel('y [mm]')
point1, = ax.plot([], [], linestyle='none', marker='o', color='blue', markersize=5)
line1, = ax.plot([], [], color='grey')

# magnetic field
ax2 = fig.add_subplot(4, 1, 3, autoscale_on=False, \
                     xlim=(-zl, zl), ylim=(0.95*Bmin, 1.05*Bmax))
plt.ylabel('B [T]')
zbase = np.linspace(-zl, zl, 201)
B = f.Bvec(0.0, [0.0, 0.0, zbase])
ax2.plot(zbase, B[2,:], color='black')
ax2.plot([-zl,zl], [Bm,Bm], color='grey', linestyle='dotted')
line2, = ax2.plot([], [], color='black')

# parallel velocity
ax3 = fig.add_subplot(4, 1, 4, autoscale_on=False,
                      xlim=(-zl, zl), ylim=(-1.1e-6*v_par0, 1.1e-6*v_par0))
plt.ylabel('v_par [10^6 m/s]')
plt.xlabel('z [m]')
ax3.plot([-zl,zl], [0.0,0.0], color='grey', linestyle='dotted')
point3, = ax3.plot([], [], linestyle='none', marker='o', color='blue', markersize=5)
line3, = ax3.plot([], [], color='black')

#------------------------------------------------------------
# set up animation
def animate(i):
    global p
    p.step(); pos=p.solution[:,-1]
    point1.set_data([pos[2]], [1e3*pos[1]])
    line1.set_data(p.solution[2,:], 1e3*p.solution[1,:])
    point3.set_data([pos[2]], [1e-6*pos[5]])
    line3.set_data(p.solution[2,:], 1e-6*p.solution[5,:])  # vz to approximate vpar
    return point1,line1,point3,line3,

r0 = [0.0, r_L, 0.0]
p.init( np.concatenate([r0, v0]) )

ani = animation.FuncAnimation(fig, animate, frames=nframes, \
                              interval=20, blit=True, repeat=False)
plt.show()
# ------------------------ end of mirror_anim.py -------------------
