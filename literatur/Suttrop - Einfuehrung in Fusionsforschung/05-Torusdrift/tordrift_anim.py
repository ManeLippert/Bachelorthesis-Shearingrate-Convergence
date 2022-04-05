#!/usr/bin/env python

# tordrift_anim.py -
# Gyrocenter motion in toroidal (axisymmetrical) geometry;
# not quite a tokamak equilibrium, though, and div(B)!=0 :-(  but
# nevertheless, this illustrates a few nice properties of particle orbits ...

# blue curves: numerical integration of drift orbit ODE
# red curves: analytical large aspect ratio estimates

# Wolfgang Suttrop  from tordrift.sce

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.integrate import odeint
import numpy as np

# particle parameters
q = 1.602E-19		 # particle charge [As]
m = 1.67E-27 		 # particle mass [kg]
E0 = 60000 * q           # kinetic energy [J]

# geometry
epsilon=0.18
R0 = 1.65	         # [m] major radius
r0 = R0*epsilon  	 # [m] minor radius particle orbit
a=0.5                    # [m] minor radius plasma

# initial v_par / v
pitch = 0.78            # a passing particle
# pitch = 0.4             # a trapped particle
# pitch = 0.5702          # marginally trapped

# field parameters
B0 = 2.0	         # [T] central toroidal field
Bp0 = 0.4	         # [T] poloidal field at r=a
Er0 = 0                  # [kV/m] radial (perpendicular) E field, negative=inward
Ep0 = 0                  # [kV/m] parallel E field

# time base to obtain result on
tg = np.arange(0.0, 3e-4, 1e-7)     # tstart, tend, dt [s]

# Vectors are in cylindrical coordinates (R, phi, z)
x0 = np.array([R0+r0, 0.0, 0.0])    # gc start position

# --------------- subroutines -----------------------------------

# Electric field vector; Bvec and B1 parameters for convenience
def E(x, Bvec, B1):
  Evec = Er0 * np.array([Bvec[2], 0.0, -Bvec[0]]) /B1 + Ep0 * Bvec /B1
  return Evec

# Magnetic induction (B vector), gradient (abs(B)) and abs(B)^2
# at position x = (R, phi, z)
# caution: r>a case not treated
def B(x):
  R = x[0,]; z = x[2,]
  # r = np.sqrt((R-R0)**2 + z*z)
  Bvec = np.array([ -Bp0*z/a, B0*R0/R, Bp0*(R-R0)/a ])
  B2 = np.sum(Bvec * Bvec, axis=0)
  B1 = np.sqrt(B2)
  gradB = [ (-Bvec[1,] * Bvec[1,] /R - (Bp0**2)*(R-R0)/a**2) /B1, \
            0.0, -(Bp0**2) *z /a**2 /B1 ]
#  np.zeros(R.shape[1:])
  return Bvec, gradB, B2

# d[x,vpar]/dt in guiding centre approximation 
def dydt_gc (y, t):
  x = y[0:3]
  vp = y[3]    # parallel velocity
  B_, gradB_, B2 = B(x)
  B1 = np.sqrt(B2)
  E_ = E(x,B_,B1)

# motion along field
  vD = vp * B_ / B1
# add ExB drift
  vD = vD + np.cross(E_, B_) / B2
# add (B x gradB) drift velocity
  vD = vD + (mu/q) * np.cross(B_ , gradB_ ) / B2
# add curvature drift (zero current density approximation)
  vD = vD + m/q * vp * vp * (np.cross(B_, gradB_)) / (B2 * B1)
# parallel acceleration: E_par, momentum conservation
  dvdt = (q/m)*np.sum(E_ * B_) / np.sqrt(B2)  - (mu/m)*np.sum(gradB_ * B_) /B1
#  dvdt = - (mu/m)*sum(gradB_ * B_) /B1
  return np.append(vD, dvdt)


# ----------- main program ----------

# initial values
B_, gradB0, B2 = B(x0)
B1 = np.sqrt(B2)
v0 = np.sqrt(0.5 * E0 / m)
v_par0 = pitch * v0
v_perp2 = v0**2 - v_par0**2
mu = 0.5 * m * v_perp2 / B1  # magnetic moment

print("Trapped/passing boundary: sqrt(2*epsilon/(1-epsilon)) = ", \
      np.sqrt(2*epsilon/(1-epsilon)))
print("v_par / v_perp = ", v_par0 / np.sqrt(v_perp2))

# integrate orbit
y0 = np.append (x0, v_par0)
yg = odeint (dydt_gc, y0, tg)

# .............  plot results ...................
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['text.usetex'] = False
matplotlib.rcParams['figure.figsize'] = [11.0,10.0]
matplotlib.rcParams['figure.subplot.left'] = 0.08
matplotlib.rcParams['figure.subplot.right'] = 0.975
matplotlib.rcParams['figure.subplot.bottom'] = 0.06
matplotlib.rcParams['figure.subplot.top'] = 0.96
matplotlib.rcParams['figure.subplot.wspace'] = 0.26
matplotlib.rcParams['figure.subplot.hspace'] = 0.22

fig1 = plt.figure()
color='blue'

# R,z projection of result gc orbit
ax1 = plt.subplot(221)
# plt.plot (yg[:,0],yg[:,2], color=color)
l1, = plt.plot ([],[], color=color)
l1a, = plt.plot([],[], color='red', marker='o', markersize=5)
plt.xlim(1.1*np.min(yg[:,0]), 1.1*np.max(yg[:,0]))
plt.ylim(1.1*np.min(yg[:,2]), 1.1*np.max(yg[:,2]))
plt.title("R-z projection")
plt.xlabel("R [m]")
plt.ylabel("z [m]")
plt.axis("equal")
# magnetic surface
circ1=plt.Circle((R0,0.0), radius=r0, color='black', fill=False)
ax1.add_patch(circ1)
# passing particles: shifted circle orbit
vD0 = m/q * (v_par0*v_par0 + 0.5*v_perp2) /R0 /B0
B_, gradB0, B2 = B(x0)
Bpt0 = np.sqrt(B_[0]*B_[0] + B_[2]*B_[2])
rshift= vD0/v_par0 * B0/Bpt0*r0
circ2=plt.Circle((R0-rshift,0.0), radius=r0+rshift, color='red', fill=False)
ax1.add_patch(circ2)
# trapped particles: banana width
xbanana = 4 * vD0 / v0 * B0/Bpt0 * np.sqrt(r0*R0)
plt.plot([R0+r0, R0+r0+xbanana], [0.0,0.0], color='red')

# Toroidal/vertical (phi, z) motion 
plt.subplot(222)
# plt.plot (180.0/np.pi*yg[:,1],yg[:,2], color=color)
l2, = plt.plot ([],[], color=color)
l2a, = plt.plot([],[], color='red', marker='o', markersize=5)
plt.xlim(np.min(180.0/np.pi*yg[:,1]), np.max(180.0/np.pi*yg[:,1]))
plt.ylim(np.min(yg[:,2]),np.max(yg[:,2]))
plt.title("Toroidal/vertical")
plt.xlabel("$\phi$ [deg]")
plt.ylabel("z [m]")

# Top view (R, phi)
ax3 = plt.subplot(223)
# plt.plot (yg[:,0]*np.cos(yg[:,1]),yg[:,0]*np.sin(yg[:,1]), color=color)
l3x = yg[:,0]*np.cos(yg[:,1])
l3y = yg[:,0]*np.sin(yg[:,1])
circ3=plt.Circle((0.0,0.0), radius=R0+r0, color='black', fill=False)
ax3.add_patch(circ3)
circ4=plt.Circle((0.0,0.0), radius=R0-r0, color='black', fill=False)
ax3.add_patch(circ4)
l3, = plt.plot([],[], color=color)
l3a, = plt.plot([],[], color='red', marker='o', markersize=5)
plt.xlim(np.min(l3x),np.max(l3x))
plt.ylim(np.min(l3y),np.max(l3y))
plt.title("Top view")
plt.xlabel("x [m]")
plt.ylabel("y [m]")
plt.axis("equal")


# Parallel velocity vs. time
plt.subplot(224)
# plt.plot (1e6*tg, 1e-3*yg[:,3], color=color)
l4, = plt.plot ([], [], color=color)
plt.xlim(np.min(1e6*tg),np.max(1e6*tg))
plt.ylim(np.min(1e-3*yg[:,3]),np.max(1e-3*yg[:,3]))
plt.title("Parallel velocity")
plt.xlabel("time [$\mu$s]")
plt.ylabel("v$_\|$ [km/s]")

# =================== animation ======================================

def update_line(num):
    l1.set_data([yg[:num,0], yg[:num,2]])
    l1a.set_data([yg[num,0], yg[num,2]])
    l2.set_data([180.0/np.pi*yg[:num,1],yg[:num,2]])
    l2a.set_data([180.0/np.pi*yg[num,1],yg[num,2]])
    l3.set_data([l3x[:num],l3y[:num]])
    l3a.set_data([l3x[num],l3y[num]])
    l4.set_data([1e6*tg[:num], 1e-3*yg[:num,3]])
    return l1,l1a,l2,l2a,l3,l3a,l4

ntimes = yg.shape[0]
line_ani = animation.FuncAnimation(fig1, update_line, ntimes,
                        interval=20, blit=True, repeat=False)

# write animation to MP4 file - CAUTION: This takes time!
# line_ani.save("tordrift_anim.mp4")

plt.show()

# -------------- end of tordrift_anim.py -----------------------
  
