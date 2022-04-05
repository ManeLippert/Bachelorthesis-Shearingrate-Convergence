#!/usr/bin/env python3
# animate gyromotion

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation

import particle

e = 1.602E-19   # As   elementary charge
me = 9.11E-31   # kg   electron mass
mp = 1.67E-27   # kg   proton mass

# the E and B fields
Bz0 = 1.0
B0    = [0.0, 0.0, Bz0]  # [T]
# gradB = 30.0  # [T/m]
gradB = 0.0  # [T/m]
# E0 = np.array([1e4, 0.0, 0.0])
E0 = np.array([1e4, 0.0, 0.0])
f = particle.SimpleFields(E0, B0, gradB)

# particle
q = e   
m = mp
Ekin0 = 200*e   # kinetic energy [J]
v_perp0 = np.sqrt(2*Ekin0/m) # initial (perpendicular) velocity
    
# analytical solution
omega_c = q*np.sqrt(np.dot(B0,B0))/m
r_L = v_perp0/np.abs(omega_c)
    
# simulation parameters
ngyrations = 10
nframes = 200 * ngyrations
tmax = ngyrations * (2*np.pi)/np.abs(omega_c)
tstep = tmax/nframes    # time step for each frame

# the simulation
p = particle.FullOrbit(f, tstep, m, q, track_vperp=True)
x0 = [0.0, r_L, 0.0]
v0 = [np.sign(q)*v_perp0, 0.0, 0.0]

# plot
matplotlib.rcParams['figure.figsize'] = [14.0,7.0]
matplotlib.rcParams['figure.subplot.left'] = 0.05
matplotlib.rcParams['figure.subplot.right'] = 0.98
matplotlib.rcParams['figure.subplot.bottom'] = 0.07
matplotlib.rcParams['figure.subplot.top'] = 0.98

fig = plt.figure()
ab = 4*1.1e3*r_L  # plot boundary (radius)

# x-y trajectory
ax = fig.add_subplot(1, 2, 1, aspect='equal', autoscale_on=False,
                     xlim=(-ab, ab), ylim=(-ab, ab))
plt.xlabel('x [mm]')
plt.ylabel('y [mm]')
ax.plot([-ab,ab], [0.0,0.0], color='#D0D0D0')
ax.plot([0.0,0.0], [-ab,ab], color='#D0D0D0')
GyroOrbit = plt.Circle((0,0),radius=1e3*r_L, color='grey', linestyle='dotted', fill=False)
ax.add_artist(GyroOrbit)
point1, = ax.plot([], [], linestyle='none', marker='o', color='blue', markersize=5)
line1, = ax.plot([], [], color='black')

# magnetic field vs. time
ax2 = fig.add_subplot(3, 2, 2, autoscale_on=False, \
                     xlim=(0, 1e6*tmax), ylim=(0.9*Bz0, 1.1*Bz0))
plt.ylabel('B [T]')
ax2.plot([0,1e6*tmax], [Bz0,Bz0], color='grey', linestyle='dotted')
line2, = ax2.plot([], [], color='black')

# orbit velocity vs. time
ax3 = fig.add_subplot(3, 2, 4, autoscale_on=False,
                      xlim=(0, 1e6*tmax), ylim=(0.9e-6*v_perp0, 1.1e-6*v_perp0))
plt.ylabel('v_perp [10^6 m/s]')
ax3.plot([0,1e6*tmax], [1e-6*v_perp0,1e-6*v_perp0], color='grey', linestyle='dotted')
line3, = ax3.plot([], [], color='black')

# gyroradius vs. time
ax4 = fig.add_subplot(3, 2, 6, autoscale_on=False,
                     xlim=(0, 1e6*tmax), ylim=(0.9e3*r_L, 1.1e3*r_L))
plt.xlabel('t [us]')
plt.ylabel('r_L [mm]')
ax4.plot([0,1e6*tmax], [1e3*r_L,1e3*r_L], color='grey', linestyle='dotted')
line4, = ax4.plot([], [], color='black')


#------------------------------------------------------------
# set up animation
def animate(i):
    global p
    p.step(); pos=1e3*p.solution[:,-1]
    point1.set_data([pos[0]], [pos[1]]);
    line1.set_data(1e3*p.solution[0,:], 1e3*p.solution[1,:])
    line2.set_data(1e6*p.timebase, p.Bnorm)
    line3.set_data(1e6*p.timebase, 1e-6*p.vperp)
    rLt = p.m/np.abs(p.q)*p.vperp/p.Bnorm
    line4.set_data(1e6*p.timebase, 1e3*rLt)
    return point1,line1,line2,line3,line4,

p.init( np.concatenate([x0, v0]) )
ani = animation.FuncAnimation(fig, animate, frames=nframes, \
                              interval=20, blit=True, repeat=False)

plt.show()
