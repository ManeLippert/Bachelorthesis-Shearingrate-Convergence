#!/usr/bin/env python3

"""
Animated electron undergoing a Coulomb collision

animation adapted from code at 
http://matplotlib.sourceforge.net/examples/animation/double_pendulum_animated.py
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation

from coulomb import e, epsilon0, m_e, Particle_Coulomb, Particle_shielded

# global parameters
E0 = 2*e   # [eV]   incident kinetic energy
vz0 = np.sqrt(2*E0/m_e)   # [m/s] incident velocity

q1 = e          # [As] of scattering centre
q2 = -e         # [As] of scattered particle
# [m]  impact parameter for 90 degrees deflection
r0perp = np.abs(q1*q2)/(8*np.pi*epsilon0*E0)

# impact parameter in units of r0perp
r0fac = 1.0

Chi = 2*np.arctan(1/r0fac)*np.sign(q1*q2)
print('Chi: ', Chi)

#------------------------------------------------------------
# simulation parameters
nframes = 250  # number  of animation frames
zd = 40*r0perp  # domain size
rd = 30*r0perp
dt = 2*zd/vz0/nframes  # time step (simulated time) for each frame
lambda_D = 10.0*r0perp  # Debye length

# one particle for each potential type
p = Particle_Coulomb  (dt, m=m_e, lambda_D=None) # Coulomb 
#p = Particle_Coulomb  (dt, m=m_e, lambda_D=lambda_D)  # truncated Coulomb
#p = Particle_shielded (dt, m=m_e, lambda_D=lambda_D)  # shielded Coulomb

#------------------------------------------------------------
# set up figure

matplotlib.rcParams['figure.figsize'] = [10.0,8.0]
matplotlib.rcParams['figure.subplot.left'] = 0.09
matplotlib.rcParams['figure.subplot.right'] = 0.975
matplotlib.rcParams['figure.subplot.bottom'] = 0.05
matplotlib.rcParams['figure.subplot.top'] = 0.98

fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal', autoscale_on=False,
                     xlim=(-zd, zd), ylim=(-rd, rd))

plt.xlabel('z [m]')
plt.ylabel('r [m]')
r0 = r0perp*r0fac
ax.plot([-zd,0], [r0perp,r0perp],linestyle='dotted',color='grey')
ax.plot([-zd,zd], [r0-zd*np.tan(Chi), r0+zd*np.tan(Chi)],\
        linestyle='dotted',color='grey')
ax.plot([-zd,zd],[0,0],linestyle='solid',color='black',linewidth=0.2)
ax.plot([0,0],[-rd,rd],linestyle='solid',color='black',linewidth=0.2)
# DebyeSphere = plt.Circle((0,0),radius=lambda_D,color='#D0D0D0',linestyle='dotted')
# ax.add_artist(DebyeSphere)
point1, = ax.plot([], [], linestyle='none', marker='o', color='red', markersize=10)
line1, = ax.plot([], [], color='black')

#------------------------------------------------------------
# set up animation
def init():
    global p
    pos=p.position()
    point1.set_data([pos[1]], [pos[0]]); line1.set_data([], [])
    return point1,line1

def animate(i):
    global p
    p.nstep(); pos=p.position()
    point1.set_data([pos[1]], [pos[0]]); line1.set_data(p.z, p.r)
    return point1,line1,

init_state = [r0fac*r0perp, -zd, 0.0, vz0]
p.init(init_state);
ani = animation.FuncAnimation(fig, animate, frames=nframes,
                              interval=20, blit=True)
ax.text(0.2*zd, 0.85*rd, 'r0/r0perp = '+str(r0fac))
# save the animation as an mp4. See also:
# http://matplotlib.sourceforge.net/api/animation_api.html
#ani.save('coulomb.mp4', fps=30, extra_args=['-vcodec', 'libx264'])
plt.show()

# -------------------- end of coulomb_anim1.py --------------------------------
