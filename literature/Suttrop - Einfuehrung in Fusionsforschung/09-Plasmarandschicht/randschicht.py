#!/usr/bin/python

from scipy.integrate import odeint
from pylab import * # for plotting commands
from math import exp, sqrt
import numpy as np

# normalised temperature, 2kTe / (m_i v_o^2)
t = 2.0     # Bohm criterium
# t = 5.03      # violating Bohm criterium

# x = normalised position (x/lambda_D)
# y[0] = eta, normalised potential,  eta' = E
# y[1] = normalised electric field,  E' = rho/epsilon

def deriv(y,x):
    # return derivatives of the array y
    a = 1.0 / sqrt(1 + t*y[0])  # ion density
    b = exp(-y[0])  # electron density
    return array([y[1], a-b])

x = linspace(0.0,37.0,300)
yinit = array([0.0,1e-3]) # initial values
y = odeint(deriv,yinit,x)

# plot parameters
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['figure.figsize'] = [12,9]
matplotlib.rcParams['figure.subplot.left'] = 0.1
matplotlib.rcParams['figure.subplot.right'] = 0.98
matplotlib.rcParams['figure.subplot.bottom'] = 0.07
matplotlib.rcParams['figure.subplot.top'] = 0.98
matplotlib.rcParams['figure.subplot.wspace'] = 0.2
matplotlib.rcParams['figure.subplot.hspace'] = 0.1

figure()
subplot(2,2,1)
plot(x,y[:,0], color='red')
ylabel('norm. potential $\eta$')

subplot(2,2,2)
c = np.abs(1.0 / np.sqrt(1 + t*y[:,0]) - np.exp(-y[:,0]))
semilogy(x, c, color='red')
ylabel('normalised charge density')

subplot(2,2,3)
semilogy(x,np.abs(y[:,0]), color='red')
xlabel('normalised distance')
ylabel('ln |$\eta$|')

subplot(2,2,4)
semilogy(x,np.abs(y[:,1]), color='red')
xlabel('normalised distance')
ylabel('normalised electrical field E')

show()
