# tokscaling.py
# scaling of tokamak parameters with major radius R

import matplotlib
import matplotlib.pyplot as plt
import numpy as np


# "free" design parameters
A=3           # aspect ratio (=1/epsilon)
kappa = 1.7   # elongation
Bmax = 8      # [T] limit for superconducting magnets
qa = 3        # edge safety factor
Tk = 10       # [keV] burn temperature in plasma centre

R = np.linspace(1.0, 10.0, 101)

# radius-dependent quantities
Vol = kappa/(A**2)*2 *np.pi**2 *R**3;      # [m^3] torus volume
Bt  = Bmax*(1-1/A)                         # [T]   (independent of R!)
Ip  = 1.328*5*Bt*R/(A**2) *0.5*(1+kappa**2)/qa;     # [MA] ITER shape
nel = 0.85*Ip/np.pi / (R/A)**2;            # [1e20 m-3], electron density = 0.85 * nGreenwald
nDT = 0.5*nel                              # [1e20 m-3] deuterium or tritium density

# Confinement time, required from Lawson criterion  n.tau = 10^21 m^-3 s
tauEreq = 10.0 / nel             # [s] 
# Alpha particle power assuming T0=10 keV and 1/3 active volume
Palpha = 0.3 *Vol *6200e-6 *Tk**2 *nDT**2  # [MW]
# Kinetic stored energy  (electrons and ions, profile factor 0.5)
Wkin = (3.0/2.0) *2 *0.5*Vol *16020e-6 *nel *Tk   # [MJ]
# Loss power needed for this energy
PL = 365.24 *A**0.6216 /kappa**1.81 /Bt**0.2162 * Wkin**2.703 / Ip**2.622 / nel**1.108 / R**5.22  # [MW] conducted
# Radiation losses
Pbr = 4850e-6*nel *nel *np.sqrt(10.0) *Vol #  [MW] Bremsstrahlung
Prad = 10000e-6 *nel**2 *Vol  # [MW] High Z radiation (fZ=1e-4)
# Heating power needed
Pheat = PL+Prad+Pbr
# energy confinement time obtained
tauEobt = 0.11269 / A**0.23 *kappa**0.67 *Bt**0.08 * (Ip**0.97)/PL**0.63 *nel**0.41 *R**1.93;  # [s]
# beta
bet=1e6*Wkin/Vol/((Bt**2)/(np.pi*8e-7));


# --------------- plot the results ------------------------------------

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['figure.figsize'] = [18,10]
matplotlib.rcParams['figure.subplot.left'] = 0.03
matplotlib.rcParams['figure.subplot.right'] = 0.98
matplotlib.rcParams['figure.subplot.bottom'] = 0.04
matplotlib.rcParams['figure.subplot.top'] = 0.97
matplotlib.rcParams['figure.subplot.wspace'] = 0.15
matplotlib.rcParams['figure.subplot.hspace'] = 0.18

font = {'family': 'sans-serif',
        'color':  'black',
        'weight': 'normal',
        'size': 12,
        }

matplotlib.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
matplotlib.rc('text', usetex=True)

# -------------- plasma current ---------------------------------------
ax1 = plt.subplot(231)
plt.title('Plasma current')
plt.xlabel(r'major radius $R$ [m]')
plt.ylabel(r'plasma current $I_p$ [MA]')
ax1.plot(R,Ip)

# --------------- density limit ---------------------------------------
ax2 = plt.subplot(232)
plt.title('Plasma density')
plt.xlabel(r'major radius $R$ [m]')
plt.ylabel(r'85\% density limit [10$^{20}$ m$^{-3}$]')
ax2.plot(R,nel)

# --------------- confinement time ------------------------------------
ax3 = plt.subplot(233)
plt.title('Energy confinement time')
plt.xlabel(r'major radius $R$ [m]')
plt.ylabel(r'confinement time $\tau$ [s]')
ax3.plot(R,tauEreq,color='blue',label='required')
ax3.plot(R,tauEobt,color='red',label='obtained')
ax3.legend(loc='best')

# --------------- confinement time ------------------------------------
ax4 = plt.subplot(234)
plt.title('Power balance')
plt.xlabel(r'major radius $R$ [m]')
plt.ylabel(r'power $P$ [MW]')
ax4.plot(R,Pheat,color='blue',label='$P_{heat} = P_{L} + P_{brems} + P_{rad}$')
ax4.plot(R,Pbr,color='green',label='$P_{brems}$')
ax4.plot(R,Prad,color='magenta',label='$P_{rad}$')
ax4.plot(R,Palpha,color='red',label='$P_{alpha}$')
ax4.legend(loc='best')

# --------------- plasma beta ------------------------------------
ax5 = plt.subplot(235)
plt.title('Plasma beta')
plt.xlabel(r'major radius $R$ [m]')
plt.ylabel(r'beta [\%]')
ax5.plot(R,100*bet,color='blue')


plt.show()
