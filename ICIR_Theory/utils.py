#!/usr/local/opt/python@3.8/bin/python3.8
__author__ = "@Tssp"
__date__ = "12/01/21"
import numpy as np
from numpy import sqrt
from math import pi
import os
from numba import njit, prange
from scipy.special import gamma, hyp2f1

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

@njit(parallel=True)
def trapezoidal_int_v2(a, b, h, eps, eta_x, eta_z):
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
    print(eps + 0.5*(1 + eta_x + eta_z))
    for j in prange(N):
        out += np.sqrt(eta_x*eta_z) * np.exp((eps + 0.5*(1 + eta_x + eta_z))*t[j]) / \
        np.sqrt(np.sinh(t[j])*np.sinh(eta_x*t[j])*np.sinh(eta_z*t[j])) - (t[j])**(-3/2)
    out *= h
    return out

#--------------------------------------- Chen, Zhang 2020 ---------------------------------------

@njit(parallel=True)
def A3D_int(a, b, h, eta_x, eta_y, eta_z, E):
    '''
    Chen et al. (2020)
    Parameters:
    -----------
    a:     left limit of the X-axis
    b:     right limit of the X-axis
    h:     step size
    nj:    levels
    eta_j: wj/wy
    E:     energy
    Lambda:cutoff
    
    Outputs:
    --------
    Value of the integral applying the trapezoidal integration method: step * f(t)
    '''
    out = 0
    N = int((2*b-2*a)/h)
    beta = np.arange(2*a, 2*b, h)
    for j in prange(N):
        out += -np.exp(beta[j]*E) * np.sqrt(eta_x*eta_y*eta_z/((4*pi)**3*np.sinh(eta_x*beta[j])*np.sinh(eta_y*beta[j])*np.sinh(eta_z*beta[j])))\
        + 1/(4*pi*beta[j])**(3/2)
    out*=h
    return out

def en(n, eta):
    return eta*(n + 1/2)

def W3D(nx, ny, nz, etax, etay, etaz, E):
    suma = 0
    for i in range(0, nx+1):
        for j in range(0, ny+1):
            if i%2==0 and j%2==0 and en(i, etax) + en(j, etay) + 1/2 <= E:
                suma += 2**(i + j - 1)*gamma(1/4 - (E - en(i, etax) - en(j, etay))/2)/ \
                        (gamma((1-i)/2)**2 * gamma((1-j)/2)**2 * gamma(1+i) * gamma(1+j) * \
                        gamma(3/4 - (E - en(i, etax) - en(j, etay))/2))
    return -pi/2 * sqrt(etax*etay*etaz/2) * suma

def B1_3D(nx, ny, nz, etax, etay, etaz, E, Lambda):
    #np.seterr('raise')
    suma = 0
    for i in range(0, nx+1):
        for j in range(0, ny+1):
            if i%2==0 and j%2==0 and en(i, etax) + en(j, etay) + 1/2 <= E:
                suma += 2**(i + j - 5/2) * sqrt(pi*etax*etay) * gamma(1/4 - (E - en(i, etax) - en(j, etay))/2) *\
                       np.exp((E - en(i, etax) - en(j, etay) - 3/2)*Lambda) / (gamma((1-i)/2)**2 * gamma((1-j)/2)**2 \
                       * gamma(1+i) * gamma(1+j) * gamma(5/4 - (E - en(i, etax) - en(j, etay))/2)) * sqrt(np.exp(2*Lambda) - 1)*\
                        hyp2f1(1, 3/4 - (E - en(i, etax) - en(j, etay))/2, 5/4 - (E - en(i, etax) - en(j, etay))/2, np.exp(-2*Lambda))
    return (-1) * suma

def B2_3D(nx, ny, nz, etax, etay, etaz, E, Lambda):
    suma = 0
    for i in range(0, nx+1):
        for j in range(0, ny+1):
            if en(i, etax) + en(j, etay) + 1/2 >= E:
                suma += 2**(i + j - 1) * (np.exp(2*Lambda) - 1) * hyp2f1(1, 3/4 - (E - en(i, etax) - en(j, etay))/2,\
                        5/4 - (E - en(i, etax) - en(j, etay))/2, np.exp(-2*Lambda)) / (gamma((1-i)/2)**2 * gamma((1-j)/2)**2 * gamma(1+i) * gamma(1+j) * (2*(E - en(i, etax) - en(j, etay)) - 1))
    return sqrt(pi*etax*etay/np.sinh(Lambda)) * suma

def separate_levels(A, E):
    dic = {}
    count = 0
    level = 0
    for i in range(len(E)-1):
        if np.sign(A[i]) != np.sign(A[i+1]):
            level+=1
            dic[f'a3D_n{level}'] = A[i-count:i+1]
            dic[f'E_n{level}']   = E[i-count:i+1]
            count = 0
        elif np.sign(A[i]) == np.sign(A[i+1]):
            count += 1
    dic[f'a3D_n{level+1}'] = A[i-count+1:i+2]
    dic[f'E_n{level+1}']   = E[i-count+1:i+2]
    return dic, level
#---------------------------------------------------------------------------------------------------------------------------