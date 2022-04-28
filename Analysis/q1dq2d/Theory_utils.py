#!/usr/local/opt/python@3.8/bin/python3.8
__author__ = "@Tssp"
__date__ = "02/12/20"
import numpy as np
from math import pi
import os
from numba import njit, prange

#--------------------------------------- Functions ---------------------------------------
def Ecm(wx, wy, wz, nx, ny, nz):
  '''
  Parameters:
  -----------
  wj: oscillation frequencies in a.u on the 
  -Axis.
  nj: excitation levels on the 3-Axis.

  Returns:
  --------
  HO Energy.
  '''
  return wx * (nx + 0.5) + wy * (ny + 0.5) + wz * (nz + 0.5)

def Ecm_n(wx, wy, wz, Vx, Vy, Vz, nx, ny, nz):
  '''
  Parameters:
  -----------
  wj: oscillation frequencies in a.u on the 3-Axis.
  nj: excitation levels on the 3-Axis.
  Vj: potential depths in a.u on the 3-Axis.

  Returns:
  --------
  1st order perturbation energy for the sextic anharmonic potential.
  '''
  Ecm_n_harm = wx * (nx + 1/2) + wy * (ny + 1/2) + wz * (nz + 1/2)
  Ecm_n_anharm = -1/(1152*Vx**2) * (36*(2*nx**2 + 2*nx + 1)*Vx*wx**2 - (4*nx**3 + 6*nx**2 + 8*nx + 3)*wx**3 ) + \
  -1/(1152*Vy**2) * (36*(2*ny**2 + 2*ny + 1)*Vy*wy**2 - (4*ny**3 + 6*ny**2 + 8*ny + 3)*wy**3 ) + \
  -1/(1152*Vz**2) * (36*(2*nz**2 + 2*nz + 1)*Vz*wz**2 - (4*nz**3 + 6*nz**2 + 8*nz + 3)*wz**3 )
  return Ecm_n_harm + Ecm_n_anharm

# Integral
@njit(parallel=True)
def trapezoidal_int(a, b, h, eps, eta_x, eta_z):
    '''
    Parameters:
    -----------
    a:     left limit of the X-axis
    b:     right limit of the X-axis
    h:     step size
    eps:   epsilon
    eta_j: wj/wy
    
    Outputs:
    --------
    Value of the integral applying the trapezoidal integration method: step * f(t)
    '''
    out = 0
    N = int((2*b-2*a)/h)
    t = np.arange(2*a, 2*b, h)
    for j in prange(N):
        out += np.sqrt(eta_x*eta_z) * np.exp(eps*t[j]/2) / \
        np.sqrt((1 - np.exp(-t[j]))*(1 - np.exp(-eta_x*t[j]))*(1 - np.exp(-eta_z*t[j]))) - t[j]**(-3/2)
    out *= h
    return out
#---------------------------------------------------------------------------------------------------------------------------