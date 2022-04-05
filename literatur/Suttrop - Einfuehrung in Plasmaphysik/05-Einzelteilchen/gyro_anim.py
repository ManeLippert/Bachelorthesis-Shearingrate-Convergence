#!/usr/bin/env python3
# animate gyromotion

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation

import particle

# constants
e = 1.602E-19   # As   elementary charge
m_e = 9.11E-31   # kg   electron mass
m_p = 1.67E-27   # kg   proton mass

# homogeneous magnetic field, no electrical field
B0    = [0.0, 0.0, 1.0]  # [T]
gradB = [0.0, 0.0, 0.0]  # [T/m]
E0 = np.array([0.0, 0.0, 0.0])
f = particle.SimpleFields(E0, B0, 0.0)  # the EM field

# particle
Ekin0 = 200*e   # kinetic energy [J]
if False:   # electron
    m = m_e
    q = -e
    col = 'green'
else:      # proton
    m = m_p
    q = e
    col = 'red'

v_perp = np.sqrt(2*Ekin0/m) # initial (perpendicular) velocity
    
# analytical solution
omega_c = q*np.sqrt(np.dot(B0,B0))/m
r_L = v_perp/np.abs(omega_c)
    
# simulation parameters
ngyrations = 10
nframes = 100 * ngyrations
tmax = ngyrations * (2*np.pi)/np.abs(omega_c)
tstep = tmax/nframes    # time step for each frame

# the simulation
p = particle.FullOrbit(f, tstep, m, q)
x0 = [0.0, r_L, 0.0]
v0 = [q/e*v_perp, 0.0, 0.0]

# plot
matplotlib.rcParams['figure.figsize'] = [8.0,8.0]
matplotlib.rcParams['figure.subplot.left'] = 0.09
matplotlib.rcParams['figure.subplot.right'] = 0.975
matplotlib.rcParams['figure.subplot.bottom'] = 0.05
matplotlib.rcParams['figure.subplot.top'] = 0.98

fig = plt.figure()
ab = 1.1e3*r_L  # plot boundary
ax = fig.add_subplot(111, aspect='equal', autoscale_on=False,
                     xlim=(-ab, ab), ylim=(-ab, ab))
plt.xlabel('x [mm]')
plt.ylabel('y [mm]')
ax.plot([-ab,ab], [0.0,0.0], color='#D0D0D0')
ax.plot([0.0,0.0], [-ab,ab], color='#D0D0D0')
GyroOrbit = plt.Circle((0,0),radius=1e3*r_L, color='grey', linestyle='dotted', fill=False)
ax.add_artist(GyroOrbit)
point1, = ax.plot([], [], linestyle='none', marker='o', color=col, markersize=10)
line1, = ax.plot([], [], color=col)

#------------------------------------------------------------
# set up animation

def animate(i):
    global p
    p.step(); pos=1e3*p.solution[:,-1]
    point1.set_data([pos[0]], [pos[1]]);
    line1.set_data(1e3*p.solution[0,:], 1e3*p.solution[1,:])
    return point1,line1,

p.init( np.concatenate([x0, v0]) )
ani = animation.FuncAnimation(fig, animate, frames=nframes,
                              interval=40, blit=True)

plt.show()
