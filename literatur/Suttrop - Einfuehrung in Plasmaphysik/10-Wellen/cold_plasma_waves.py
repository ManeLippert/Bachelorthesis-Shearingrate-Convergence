'''
  cold_plasma_waves.py  - dispersion relation for cold waves (P=0, T=0) in a plasma
'''

import numpy as np

# ---------------------------------------------------------------------------------------

class ColdElectronWave(object):
    '''
      Wave in a cold electron fluid
    '''

    def __init__(self, ompe2, omce):
        '''
           Initialize wave in a cold electron fluid
           
           ompe2: (omega_p,e / omega)^2   square of (plasma frequency normalised to frequency)
           omcs:  (omega_c,e / omega)     cyclotron frequency normalised to frequency
        '''
        self.ompe2 = ompe2
        self.omce  = omce
        self.P = 1 - ompe2
        self.R = 1 - ompe2/(1.0 - omce)
        self.L = 1 - ompe2/(1.0 + omce)
        self.S = 0.5*(self.R + self.L)
        self.D = 0.5*(self.R - self.L)

    def calc_n2(self, theta):
        '''
           set progagation direction of wave
           theta: propagation angle with respect to static 
                  background magnetic field \vec{B} (in radians)
               theta = 0: wave vector parallel \vec{B}
               theta = pi/s: wave vector perpendicular \vec{B}
        '''
        s2t = np.square(np.sin(theta))
        c2t = np.square(np.cos(theta))
        a = s2t*self.S + c2t*self.P
        b = -(s2t*self.R*self.L + (1+c2t)*self.P*self.S)
        c = self.R *self.L *self.P
        n2a = (-b + np.sqrt(np.square(b) - 4*a*c)) / 2/a
        n2b = (-b - np.sqrt(np.square(b) - 4*a*c)) / 2/a
        return np.array([n2a, n2b])


# --------------------------- main --------------------------------------------------


import matplotlib.pyplot as plt


# k || B, frequency dependence
if True:
    omfac = 0.9   # omega_p^2 / omega_c^2

    logmin = -2
    logmax = 2    
    omega_c = np.logspace(logmin, logmax, num=1000)
    omega_p2 = omfac * np.square(omega_c)

    omega_R =  0.5 + np.sqrt( (0.25 + omfac) )
    omega_L = -0.5 + np.sqrt( (0.25 + omfac) )
    
    theta = 0  # k\| B
    #  theta = np.pi/2   # k \perp B
    
    cew = ColdElectronWave(omega_p2, omega_c)
    n2 = cew.calc_n2(theta)

    fig = plt.figure(figsize=(15,8), dpi=80)
    
    ax1 = fig.add_subplot(1,2,1, xlim=(10**(-logmax), 10**(-logmin)), ylim=(-2,2))    
    ax1.plot([10**(-logmin), 10**(-logmax)], [0.0, 0.0], color='grey')
    ax1.plot([10**(-logmin), 10**(-logmax)], [1.0, 1.0], linestyle='dotted', color='grey')
    ax1.plot([omega_L, omega_L], [-2, 2], linestyle='dashed', color='blue', label='omega_L')
    ax1.plot([1.0, 1.0], [-2, 2], linestyle='dashed', color='green', label='omega_ce')
    ax1.plot([omega_R, omega_R], [-2, 2], linestyle='dashed', color='red', label='omega_R')    
    ax1.semilogx(1/omega_c, n2[0,:], linestyle='none', marker='.', color='blue', label='L wave')
    ax1.semilogx(1/omega_c, n2[1,:], linestyle='none', marker='.', color='red', label='R wave')
    plt.xlabel('omega / omega_c,e')
    plt.ylabel('n^2')
    plt.legend(loc='best')

    ax2 = fig.add_subplot(1,2,2, xlim=(10**(-logmax), 10**(-logmin)), ylim=(-2,2))    
    ax2.plot([10**(-logmin), 10**(-logmax)], [0.0, 0.0], color='grey')
    ax2.plot([10**(-logmin), 10**(-logmax)], [1.0, 1.0], linestyle='dotted', color='grey')
    ax2.plot([omega_L, omega_L], [-2, 2], linestyle='dashed', color='blue', label='omega_L')
    ax2.plot([1.0, 1.0], [-2, 2], linestyle='dashed', color='green', label='omega_ce')
    ax2.plot([omega_R, omega_R], [-2, 2], linestyle='dashed', color='red', label='omega_R')    
    ax2.semilogx(1/omega_c, 1.0/n2[0,:], linestyle='none', marker='.', color='blue', label='L wave')
    ax2.semilogx(1/omega_c, 1.0/n2[1,:], linestyle='none', marker='.', color='red', label='R wave')
    plt.xlabel('omega / omega_c,e')
    plt.ylabel('(v_ph/c)^2')
    plt.legend(loc='best')

    plt.show()

# ...........................................................................................
    
# k perp B, frequency dependence
if False:
    omfac = 0.7   # omega_p^2 / omega_c^2

    logmin = -1
    logmax = 1   
    omega_c = np.logspace(logmin, logmax, num=1000)
    omega_p2 = omfac * np.square(omega_c)

    omega_p = np.sqrt(omfac)
    omega_R =  0.5 + np.sqrt( (0.25 + omfac) )
    omega_L = -0.5 + np.sqrt( (0.25 + omfac) )
    omega_uh = np.sqrt(1.0 + omfac)
    
    theta = np.pi/2   # k \perp B
    
    cew = ColdElectronWave(omega_p2, omega_c)
    n2 = cew.calc_n2(theta)

    fig = plt.figure(figsize=(15,8), dpi=80)
    
    ax1 = fig.add_subplot(1,2,1, xlim=(10**(-logmax), 10**(-logmin)), ylim=(-5,5))    
    ax1.plot([10**(-logmin), 10**(-logmax)], [0.0, 0.0], color='grey')
    ax1.plot([10**(-logmin), 10**(-logmax)], [1.0, 1.0], linestyle='dotted', color='grey')
    ax1.plot([omega_p, omega_p], [-5, 5], linestyle='dashed', color='black', label='omega_p')
    ax1.plot([omega_L, omega_L], [-5, 5], linestyle='dashed', color='blue', label='omega_L')
    ax1.plot([1.0, 1.0], [-5, 5], linestyle='dashed', color='green', label='omega_ce')
    ax1.plot([omega_R, omega_R], [-5, 5], linestyle='dashed', color='red', label='omega_R')    
    ax1.plot([omega_uh, omega_uh], [-5, 5], linestyle='dashed', color='magenta', label='omega_uh')    
    ax1.semilogx(1/omega_c, n2[0,:], linestyle='none', marker='.', color='black')
    ax1.semilogx(1/omega_c, n2[1,:], linestyle='none', marker='.', color='black')
    plt.xlabel('omega / omega_c,e')
    plt.ylabel('n^2')
    plt.legend(loc='best')

    ax2 = fig.add_subplot(1,2,2, xlim=(10**(-logmax), 10**(-logmin)), ylim=(-5,5))    
    ax2.plot([10**(-logmin), 10**(-logmax)], [0.0, 0.0], color='grey')
    ax2.plot([10**(-logmin), 10**(-logmax)], [1.0, 1.0], linestyle='dotted', color='grey')
    ax2.plot([omega_p, omega_p], [-5, 5], linestyle='dashed', color='black', label='omega_p')
    ax2.plot([omega_L, omega_L], [-5, 5], linestyle='dashed', color='blue', label='omega_L')
    ax2.plot([1.0, 1.0], [-5, 5], linestyle='dashed', color='green', label='omega_ce')
    ax2.plot([omega_R, omega_R], [-5, 5], linestyle='dashed', color='red', label='omega_R')    
    ax2.plot([omega_uh, omega_uh], [-5, 5], linestyle='dashed', color='magenta', label='omega_uh')    
    ax2.semilogx(1/omega_c, 1.0/n2[0,:], linestyle='none', marker='.', color='black')
    ax2.semilogx(1/omega_c, 1.0/n2[1,:], linestyle='none', marker='.', color='black')
    plt.xlabel('omega / omega_c,e')
    plt.ylabel('(v_ph/c)^2')
    plt.legend(loc='best')

    plt.show()


# polar diagram
if False:
    theta = np.linspace(0, 2*np.pi, 500)

    if False:    # high frequency case - no resonance
        omega_c = 0.49
        omega_p2 = 0.3

    if False:    # resonance when k || B
        omega_c = 0.8
        omega_p2 = 0.8
        
    if False:    # whistler mode  R <-> X
        omega_c = 1.1
        omega_p2 = 1.1

    if True:    # whistler mode  R <-> X
        omega_c = 0.9
        omega_p2 = 1.1

        
    omega_R = 0.5*omega_c + np.sqrt( np.square(0.5*omega_c) + omega_p2 )
    omega_L = -0.5*omega_c + np.sqrt( np.square(0.5*omega_c) + omega_p2 )
    omega_p = np.sqrt(omega_p2)
    omega_uh = np.sqrt( omega_p2 + np.square(omega_c) )
    print("Freq  c: %5.3f  R: %5.3f  L: %5.3f  p: %5.3f  uh: %5.3f" \
          % (omega_c, omega_R, omega_L, omega_p, omega_uh) )
    
    cew = ColdElectronWave(omega_p2, omega_c)
    n2 = cew.calc_n2(theta)

    plt.polar(theta, np.zeros(len(theta)), color='black')
    plt.polar(theta, np.ones(len(theta)), color='black', linestyle='dotted')
    plt.polar(theta, n2[0,:], linestyle='none', marker='.', color='red')
    plt.polar(theta, n2[1,:], linestyle='none', marker='.', color='blue')
    ax=plt.gca()
    ax.set_rlim(-0.5,2)
    
    plt.show()
