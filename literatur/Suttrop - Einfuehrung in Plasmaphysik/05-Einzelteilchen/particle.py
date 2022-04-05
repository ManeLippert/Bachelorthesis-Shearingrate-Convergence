#!/usr/bin/env python3
import numpy as np
import scipy.integrate


# --------------------- field definition ---------------------------
class SimpleFields():
  '''
   Base class for electromagnetic field definitions
   Implements a homogeneous E field and 
   B field || z with gradient in x-direction
  '''
  def __init__(self, E, B, gradB):
    '''
      initialise simple static E and B fields with
      B = [Bx, By, Bz + gradB*x]  (field in z-direction, gradient B in x-direction)

      :E: = [Ex, Ey, Ez]
      :B: = [Bx, By, Bz] homogeneous part of B field
      :gradB: (T/m) gradient of Bz in x-direction
    '''
    self.E = np.array(E)
    self.B = np.array(B)
    self.gradB = gradB

  def Evec(self, t, r):
    '''
      return electrical field 
      at time t and position vector r.
      r can be a vector, r=[x,y,z], 
      or an array r=[xyz, caseindex]
    '''
    if r.ndim==1:
      return self.E
    else:   # same at all positions
      return np.outer(self.E, np.ones(r.shape[1]))

  def Bvec(self, t, r):
    '''
     return magnetic induction vector
     at time t and position(s) r.
     r can be a vector, r=[x,y,z], 
     or an array r=[xyz, caseindex]
    '''
    if r.ndim==1:
      return np.array([self.B[0], self.B[1], self.B[2] + self.gradB*r[0]])
    else:
      B = np.outer(self.B, np.ones(r.shape[1]))
      dBz = self.gradB*r[0]
      B[2,:] += dBz
      return B

  def grad_B(self, t, r):
    '''
     return grad B at time t and position vector r 
    '''
    if r.ndim==1:
      return [self.gradB, 0.0, 0.0]
    else:   # same at all positions
      return np.outer([self.gradB, 0.0, 0.0], np.ones(r.shape[1]))

    
# --------------------- particle integrators -----------------------
class FullOrbit():

  def __init__(self, fields, dt, m, q, track_vperp=False):
    '''
     initialise particle
     :fields: class to compute electromagnetic fields
     :dt: time step
     :m: particle mass
     :q: particle charge
    '''
    self.fields = fields
    self.dt = dt
    self.m = m
    self.q = q
    self.track_vperp = track_vperp

  def init(self, init_state):
    '''
     initialise particle state
    '''
    self.state = init_state
    self.solution = np.array([init_state]).reshape(-1,1)
    self.time = 0.0
    self.timebase = np.array([self.time])
    if self.track_vperp:
      B0 = self.fields.Bvec(0, init_state[0:3])
      B1 = np.linalg.norm(B0)
      if B1>0:   # have magnetic field
        self.Bnorm = np.array([B1])
        v0 = init_state[3:6]
        vperp0 = np.linalg.norm(np.cross(v0, B0)) / B1
        self.vperp = np.array([vperp0])
      else:
        self.track_vperp = False  # no field, no vperp

  def dstate_dt(self, t, state):
    '''
    time derivative of state (full orbit)
    '''
    r = state[0:3]
    v = state[3:6]
    E = self.fields.Evec(t, r)
    B = self.fields.Bvec(t, r)
    dvdt = (self.q/self.m) * (E + np.cross(v, B))
    return np.concatenate([v, dvdt])

  def step(self):
    '''
       follow particle and record states at a timebase
    '''
    sol = scipy.integrate.solve_ivp(self.dstate_dt, [self.time, self.time+self.dt], \
              self.state, method='Radau', atol=1e-9, rtol=1e-8,
              max_step=self.dt/10)
    # print('solve_ivp status: ', sol.status, ' nfev: ', sol.nfev)
    self.state = sol.y[:,-1]
    self.solution = np.hstack((self.solution, sol.y[:, 1:]))
    self.timebase = np.hstack((self.timebase, sol.t[1:]))
    self.time += self.dt
    if self.track_vperp:
      B0 = self.fields.Bvec(sol.t[1:], sol.y[0:3, 1:])
      B1 = np.linalg.norm(B0, axis=0)
      self.Bnorm = np.hstack((self.Bnorm, B1))
      v  = sol.y[3:6, 1:]
      vperp = np.linalg.norm(np.cross(v, B0, axis=0), axis=0) / B1
      self.vperp = np.hstack((self.vperp, vperp))
    return sol.t, sol.y
  
# ..............................................................  
class GuidingCentreOrbit(FullOrbit):
  '''
   Integrate particle orbit in guiding centre approximation
  '''

  def init(self, init_state, mu):
    '''
      initialise guiding centre orbit

      :init_state: [x0, y0, z0, vpar0]

      :mu: magnetic moment (constant of motion)
    '''
    self.mu = mu
    self.state = init_state
    self.solution = np.array([init_state]).reshape(-1,1)
    self.time = 0.0
    self.timebase = np.array([self.time])
  
  def dstate_dt(self, t, state):
     r = state[0:3]
     vpar = state[3]  # parallel velocity
     E = self.fields.Evec(t, r)
     B = self.fields.Bvec(t, r)
     B2 = np.dot(B, B)
     B1 = np.sqrt(B2)
     gradB = self.fields.grad_B(t, r)
     # Calculate guiding centre velocity:
     # - parallel motion
     vD = vpar * B / B1
     # - ExB drift velocity
     vD += np.cross(E, B) / B2
     # - (B x gradB) drift velocity
     vD += (self.mu/self.q) * np.cross(B, gradB) / B2
     # - curvature drift. Assume zero current density
     vD += self.m/self.q * vpar**2  * (np.cross(B, gradB)) / (B2 * B1)
     # Parallel acceleration by E_par and conservation of magnetic moment
     Epar = np.dot(E, B) / B1
     dB_ds = np.dot(gradB, B) / B1
     dvdt = np.array( [ (self.q/self.m)*Epar - (self.mu/self.m)*dB_ds ] )
     return np.concatenate([vD, dvdt])


# ----------------- main program ---------------------
if __name__ == '__main__':

  import matplotlib
  import matplotlib.pyplot as plt

  # constants
  e = 1.602E-19   # As   elementary charge
  me = 9.11E-31   # kg   electron mass
  mp = 1.67E-27   # kg   proton mass

  # fields
  E0 = [0, 0, 0] # V/m
  B0 = [0, 0, 1.0]   # T
  gradB = 10.0 # T/m 
  f = SimpleFields(E0, B0, gradB)

  # particle
  m = mp
  q = e
  Ekin0 = 200*e
  v0 = np.array([0, np.sqrt(2*Ekin0/mp), 0])  # initial velocity

  # simulation parameters
  tmax = 1e-6  # [s]   simulated duration
  x0 = np.array([0.0, 0.0, 0.0])  # guiding centre start position
  B0 = f.Bvec(0, x0)   # magnetic induction vector at start position
  B2 = np.dot(B0,B0)   # |B|^2
  B1 = np.sqrt(B2)     # |B|
  
  print('Solve for full orbit')
  p_fo = FullOrbit(f, tmax, m, q)
  p_fo.init( np.concatenate([x0 - mp/(e*B2)*np.cross(v0, B0), v0]) )
  p_fo.step()

  print('Solve for gyrocentre motion')
  p_gc = GuidingCentreOrbit(f, tmax, m, q)
  vp0 = np.cross(v0, B0) / B1   # velocity perpendicular B
  vp02 = np.dot(vp0, vp0)       # |vp|^2
  mu = 0.5 * mp * vp02 / B1     # magnetic moment
  v_par0 = np.array( [ np.dot(v0, B0) / B1 ] )  # parallel velocity
  p_gc.init(np.concatenate([x0, v_par0]), mu)
  p_gc.step()

  # plot result
  matplotlib.rcParams['figure.figsize'] = [10.0,8.0]
  matplotlib.rcParams['figure.subplot.left'] = 0.09
  matplotlib.rcParams['figure.subplot.right'] = 0.975
  matplotlib.rcParams['figure.subplot.bottom'] = 0.05
  matplotlib.rcParams['figure.subplot.top'] = 0.98
  
  plt.figure()

  plt.subplot(2,2,1)
  plt.plot(1e6*p_fo.timebase, 1e3*p_fo.solution[1,:], color='black')
  plt.plot(1e6*p_gc.timebase, 1e3*p_gc.solution[1,:], color='red')
  plt.xlabel('t [us]')
  plt.ylabel('y [mm]')

  plt.subplot(2,2,2, aspect='equal')
  plt.plot(1e3*p_fo.solution[0,:], 1e3*p_fo.solution[1,:], color='black')
  plt.plot(1e3*p_gc.solution[0,:], 1e3*p_gc.solution[1,:], color='red')
  plt.xlabel('x [mm]')
  plt.ylabel('y [mm]')

  plt.subplot(2,2,3)
  plt.plot(1e6*p_fo.timebase, 1e3*p_fo.solution[0, :], color='black')
  plt.plot(1e6*p_gc.timebase, 1e3*p_gc.solution[0, :], color='red')
  plt.xlabel('t [us]')
  plt.ylabel('x [mm]')

  plt.subplot(2,2,4)
  plt.plot(1e6*p_fo.timebase, 1e3*p_fo.solution[2,:], color='black')
  plt.plot(1e6*p_gc.timebase, 1e3*p_gc.solution[2,:], color='red')
  plt.xlabel('t [us]')
  plt.ylabel('z [mm]')

  plt.show()
