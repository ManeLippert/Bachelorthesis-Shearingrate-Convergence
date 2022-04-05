#!/usr/bin/env python3

"""
  statistics for coulomb scattering
"""

import numpy as np
from coulomb import e, epsilon0, m_e, Particle_Coulomb, Particle_shielded

q1 = e
q2 = -e
E0 = 2*e   # [eV]   incident kinetic energy
vz0 = np.sqrt(2*E0/m_e)   # [m/s] incident velocity
r0perp = np.abs(q1*q2)/(8*np.pi*epsilon0*E0)  # [m]  impact parameter for 90 degrees deflection

#------------------------------------------------------------
# simulation parameters
nsteps = 220
zd = 100*r0perp  # domain size
rd = 100*r0perp
dt = 2*zd/vz0
lambda_D = 5.0*r0perp  # Debye length

# one particle for each potential type
p1 = Particle_Coulomb (dt, m=m_e, lambda_D=None) # Coulomb 
p2 = Particle_Coulomb (dt, m=m_e, lambda_D=lambda_D) # truncated Coulomb
p3 = Particle_shielded (dt, m=m_e, lambda_D=lambda_D)  # shielded Coulomb

print("r0/r0perp    weight   Coulomb     Truncated   Shielded")
r0fac = np.linspace(0.1,15,300)
# r0fac = [3.0]
for rf in r0fac:
    init_state = [rf*r0perp, -zd, 0.0, vz0]
    p1.init(init_state); p2.init(init_state); p3.init(init_state);
    p1.nstep(); p2.nstep(); p3.nstep();
    dvz1 = -p1.dvz_acc(r0perp) / vz0;
    dvz2 = -p2.dvz_acc(r0perp) / vz0;
    dvz3 = -p3.dvz_acc(r0perp) / vz0;
    print("%8.3f  %8.2f  %10.3e  %10.3e  %10.3e" % (rf, p1.tot_weight, dvz1, dvz2, dvz3))

# -------------------- end of coulomb_stat.py --------------------------------
