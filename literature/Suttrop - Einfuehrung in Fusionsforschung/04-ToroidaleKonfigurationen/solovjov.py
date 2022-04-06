#!/usr/bin/env python
# solovjov.py - Analytische axisymmetrische Gleichgewichte 
# Wolfgang Suttrop, 01-Dez-2010, 09-Nov-2016

import matplotlib.pyplot as plt
import numpy as np

C1=-1; C2=-8; C3=-20; C4=20; C5=0.2  # shaped plasma
# C1=-10; C2=-0.5; C3=80; C4=0; C5=0.0  # no vertical field, no vacuum field
# C1=-10; C2=-0.5; C3=30; C4=1; C5=0.0  # small vertical
# C1=-10; C2=-0.5; C3=30; C4=3; C5=0.0  # large vertical field

R1 = np.linspace(0.1, 6, 60)
z1 = np.linspace(-4.0, 4.0, 81)
R  = np.outer(R1,np.ones(z1.size))
z  = np.outer(np.ones(R1.size),z1)

Psi = C1/2*(z**2) + C2/8*(R**4) + C3 + C4*(R**2) + C5*((R**4) - (4* (R**2) * (z**2)))
levels = np.linspace(0, 1, 21) * np.max(Psi)

ax = plt.subplot(111)
plt.axis('equal')
plt.xlabel('R')
plt.ylabel('z')
ax.contour(R1, z1, Psi.T, levels=levels)
plt.show()

