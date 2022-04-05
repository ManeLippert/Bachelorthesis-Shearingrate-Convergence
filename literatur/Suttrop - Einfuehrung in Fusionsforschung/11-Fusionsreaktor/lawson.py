#!/usr/bin/python
# Lawson criterion

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

# t1 = 2*np.logspace(0,1.65,33,10.0)
t1 = np.logspace(0,2,num=101,base=10.0)
ntau1 = 1e20*np.logspace(0,2,40,10.0)
t2 = np.outer(t1,np.ones(ntau1.size))
ntau2 = np.outer(np.ones(t1.size),ntau1)
sigmav1 = 1.1e-24 * t1 *t1

def oplot_lcr (rho, fZ, Zeff, ax):
  sigmav_ntau2 = 1.1e-24 * t2 *t2 *ntau2
  Z0 = 0.25 * rho * sigmav_ntau2
  N0 = 1+ (0.5 * rho * sigmav_ntau2) 
  fHe = Z0 / N0
  Z1 = 6 *16020 *(2-fHe+fZ-Zeff*fZ) *t2
  N1= 6200 * np.square(1-2*fHe-2*fZ) * t2 *t2
  N2= 4 * 4850 * np.sqrt(t2) *Zeff
  N3= 4 *1e9 *fZ
  rhs= Z1 / (N1-N2-N3)
  y=ntau2*1e-20 / rhs
  ax.contour(t1, ntau1, y.T, levels=[1.0])


matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['figure.figsize'] = [6,5]
matplotlib.rcParams['figure.subplot.left'] = 0.15
matplotlib.rcParams['figure.subplot.right'] = 0.98
matplotlib.rcParams['figure.subplot.bottom'] = 0.12
matplotlib.rcParams['figure.subplot.top'] = 0.98
matplotlib.rcParams['figure.subplot.wspace'] = 0.25

font = {'family': 'sans-serif',
        'color':  'black',
        'weight': 'normal',
        'size': 12,
        }

matplotlib.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
matplotlib.rc('text', usetex=True)
  
ax2 = plt.subplot(111)
ax2.set_xscale("log")
ax2.set_yscale("log") 
plt.xlabel(r'temperature [keV]')
plt.ylabel(r'n ${\tau}_E$ [$m^{-3}$ s]')

# no He, no impurities
oplot_lcr(0.0, 0.0, 1.0, ax2);

# with He
#oplot_lcr(5.0, 0.0, 1.0, ax2);
#oplot_lcr(10.0, 0.0, 1.0, ax2);

# heavy impurities
#oplot_lcr(0.0, 1e-5, 1.0, ax2);
#oplot_lcr(0.0, 1e-4, 1.0, ax2);

# with He and heavy impurities
oplot_lcr(3.0, 1e-5, 1.6, ax2);
oplot_lcr(6.0, 1e-5, 1.6, ax2);

plt.show()
