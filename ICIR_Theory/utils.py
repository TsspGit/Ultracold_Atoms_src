#!/usr/local/opt/python@3.8/bin/python3.8
__author__ = "@Tssp"
__date__ = "12/01/21"
import numpy as np
from numpy import sqrt
from math import pi
import os
from numba import njit, prange
from scipy.special import gamma, hyp2f1
import sys
eps = sys.float_info.epsilon

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

def A3D(eta_x, eta_y, E, beta):
    '''
    Chen et al. (2020)
    Parameters:
    -----------
    eta_j: wj/wy
    E:     energy
    beta:x-axis
    
    Outputs:
    --------
    Array of A3D.
    '''
    
    out = -np.exp(beta*E) * np.sqrt(eta_x*eta_y/((4*pi)**3*np.sinh(eta_x*beta)*np.sinh(eta_y*beta)*np.sinh(beta)))\
        + 1/(4*pi*beta)**(3/2)
    return out

def W3D(nx, ny, etax, etay, E):
    suma = 0
    for i in range(0, nx+1):
        for j in range(0, ny+1):
            if i%2 == 0 and j%2 == 0 and en(i, etax) + en(j, etay) + 1/2 <= E:
                suma += 2**(i + j - 1)*gamma(1/4 - (E - en(i, etax) - en(j, etay))/2)/ \
                        (gamma((1-i)/2)**2 * gamma((1-j)/2)**2 * gamma(1+i) * gamma(1+j) * \
                        gamma(3/4 - (E - en(i, etax) - en(j, etay))/2))
    return -pi/2 * sqrt(etax*etay/2) * suma

def I3D(etax, etay, nx, ny, E, beta):
    suma = 0
    for i in range(0, nx+1, 2):
        for j in range(0, ny+1, 2):
            if en(i, etax) + en(j, etay) + 1/2 <= E:
              suma += 2**(i+j-1/2)*np.exp(beta*(E - etax - etay))/ \
              (gamma((1-i)/2)**2 * gamma((1-j)/2)**2 * gamma(1+i) * gamma(1+j))
    out = sqrt(pi*etax*etay/(8*np.sinh(beta))) * suma
    return A3D(etax, etay, E, beta) + out

def I3D_int(a, b, h, nx, ny, etax, etay, etaz, E):
    N = int((2*b-2*a)/h)
    beta = np.arange(2*a, 2*b, h)
    out = 0
    for k in prange(N):
      suma = 0
      for i in prange(0, nx+1):
          for j in prange(0, ny+1):
              if i%2==0 and j%2==0 and en(i, etax) + en(j, etay) + 1/2 <= E:
                suma += 2**(i+j-1/2)*np.exp(beta[k]*(E - etax - etay))/ \
                (gamma((1-i)/2)**2 * gamma((1-j)/2)**2 * gamma(1+i) * gamma(1+j))
      out += sqrt(pi*etax*etay/(8*np.sinh(beta[k]))) * suma
    out *= h
    return A3D_int(a, b, h, etax, etay, etaz, E) + out

def B1_3D(nx, ny, etax, etay, E, Lambda):
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

def B2_3D(nx, ny, etax, etay, E, Lambda):
    suma = 0
    for i in range(0, nx+1):
        for j in range(0, ny+1):
            if en(i, etax) + en(j, etay) + 1/2 >= E:
                suma += 2**(i + j - 1) * (np.exp(2*Lambda) - 1) * hyp2f1(1, 3/4 - (E - en(i, etax) - en(j, etay))/2,\
                        5/4 - (E - en(i, etax) - en(j, etay))/2, np.exp(-2*Lambda)) / (gamma((1-i)/2)**2 * gamma((1-j)/2)**2 * gamma(1+i) * gamma(1+j) * (2*(E - en(i, etax) - en(j, etay)) - 1))
    return sqrt(pi*etax*etay/np.sinh(Lambda)) * suma

def separate_levels(A, E):
    '''
    Parameters
    ----------
    - A: x-axis (usually a3D)
    - E: y-axis (usually the energy)

    Returns
    -------
    - dic: dictionary that contains each level separately denoted by the key "a3D_n{level}" and "E_n{level}"
    - level: last level separated
    '''
    dic = {}
    count = 0
    level = 0
    for i in range(len(E)-1):
        if np.sign(A[i]) != np.sign(A[i+1]) and A[i] > abs(2):
            level+=1
            dic[f'a3D_n{level}'] = A[i-count:i+1]
            dic[f'E_n{level}']   = E[i-count:i+1]
            count = 0
        elif np.sign(A[i]) == np.sign(A[i+1]):
            count += 1
    dic[f'a3D_n{level+1}'] = A[i-count+1:i+2]
    dic[f'E_n{level+1}']   = E[i-count+1:i+2]
    return dic, level

def gauss(f,n,a,b):
    from scipy.special.orthogonal import p_roots
    [x,w] = p_roots(n+1)
    G=0.5*(b-a)*sum(w*f(0.5*(b-a)*x+0.5*(b+a)))
    return G
#---------------------------------------------------------------------------------------------------------------------------