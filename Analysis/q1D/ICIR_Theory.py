#!/usr/local/opt/python@3.8/bin/python3.8
__author__ = "@Tssp"
__date__ = "18/11/20"
import numpy as np
from math import pi
import os
from numba import njit, prange
from utils.atomic_units import ao, vo, e, hbar, me, Eh, to
from input_ICIR_Theory import *

# Parameters and Constants:
mass  = m * 1.66053873e-27 / me # a.u
wL    = wL * 1e-9 / ao # a.u
ky    = 2*np.pi/wL
# Intensities:
Ix    = Ix * (1e4 / Eh * to * ao**2)
Iy    = Iy * (1e4 / Eh * to * ao**2)
Iz    = Iz * (1e4 / Eh * to * ao**2)
# Potential depths:
Vy    = alpha * Iy # a.u
Vx    = alpha * Ix # a.u
Vz    = alpha * Iz # a.u
# Frequencies:
wx    = np.sqrt( 2 * Vx * ky**2 / mass)       
wy    = np.sqrt( 2 * Vy * ky**2 / mass)  
wz    = np.sqrt( 2 * Vz * ky**2 / mass)
eta_x = wx/wy
eta_z = wz/wy
# Characteristic distances:
dhox  = np.sqrt(2 / (mass * wx))
dhoy  = np.sqrt(2 / (mass * wy))
dhoz  = np.sqrt(2 / (mass * wz))
print(f"""            
                       Parameters
           -------------------------------------------------
           mass (a.u)     =    {mass}
           Ix(mW/cm2)     =    {Ix /(1e4 / Eh * to * ao**2)}
           Iy(mW/cm2)     =    {Iy /(1e4 / Eh * to * ao**2)}
           Iz(mW/cm2)     =    {Iz /(1e4 / Eh * to * ao**2)}
           wL(nm)         =    {wL * 0.0529177249}
           alpha(a.u)     =    {alpha}
           Vx(a.u)        =    {Vx}
           Vy(a.u)        =    {Vy}
           Vz(a.u)        =    {Vz}
           wx(a.u)        =    {wx}
           wy(a.u)        =    {wy}
           wz(a.u)        =    {wz}
           eta_x          =    {eta_x}
           eta_z          =    {eta_z}
           dho_x(a.u)     =    {dhox}
           dho_y(a.u)     =    {dhoy}
           dho_z(a.u)     =    {dhoz}
           h              =    {h}
      """)

# Energies needed for compute the integral:
E_ICIR = wy*np.loadtxt(f'Results/ICIR_positions_{int(Ix/(1e4 / Eh * to * ao**2))}_{int(Iy/(1e4 / Eh * to * ao**2))}_{int(Iz/(1e4 / Eh * to * ao**2))}.txt')[2]
Eo = (wx + wy + wz)/2
En_CM = [1.5460962687803033e-10, 1.7540673119226125e-10] 
print(f"""       
                  Energies
         -------------------------
         E_ICIR:
                -(0,2,0):   {E_ICIR[0]}
                -(2,0,0):   {E_ICIR[1]}
         Eo:                {Eo}
         E_CM:
                -(0,2,0):   {En_CM[0]}
                -(2,0,0):   {En_CM[1]}
      """)

# Load CM and rm energies:
folder_path = f'Simulations/ix{int(Ix/(1e4 / Eh * to * ao**2))}_iy{int(Iy/(1e4 / Eh * to * ao**2))}_iz{int(Iz/(1e4 / Eh * to * ao**2))}/orbitals/eva/'
for file in os.listdir(folder_path + 'rm'):
    if 'noint' in file:
        print(f'rm/ folder, file readed:\n{file}')
        Erm = np.loadtxt(folder_path + 'rm/' + file)[1,2]
for file in os.listdir(folder_path + 'CM'):
    print(f'\nCM/ folder, file readed:\n{file}')
    ECM = np.loadtxt(folder_path + 'CM/' + file)[0,2]
print(f"\nRelative energy of the first trap state: {Erm}\nCM fundamental energy: {ECM} ")

# Calculations
Eref = Erm + ECM
if mode == 'C':
    C = np.abs((E_ICIR - Eref)/wz)
elif mode == 'W':
    C = np.array([C, C])
epsilon = (C*wz + 2 * ECM - En_CM - Eo)/wy
print(f"""
          C:       
            -(0,2,0):   {C[0]}
            -(2,0,0):   {C[1]}
          epsilon:       
            -(0,2,0):   {epsilon[0]}
            -(2,0,0):   {epsilon[1]}
""")

# Integral
@njit(parallel=True)
def trapezoidal_int(a, b, h, eps):
    '''
    Parameters:
    -----------
    a:   left limit of the X-axis
    b:   right limit of the X-axis
    h:   step size
    eps: epsilon
    
    Outputs:
    --------
    Value of the integral applying the trapezoidal integration method: step * f(t)
    '''
    out = 0
    N = int((2*b-2*a)/h)
    t = np.arange(2*a, 2*b, h)
    for j in prange(N):
        out += np.sqrt(eta_x*eta_z) * np.exp(eps*t[j]/2) / np.sqrt((1 - np.exp(-t[j]))*(1 - np.exp(-eta_x*t[j]))*(1 - np.exp(-eta_z*t[j]))) - t[j]**(-3/2)
    out *= h
    return out
# Create the partitions of the integral axis.
a = h*1e-3
b = h*1e3
integral020 = 0
integral200 = 0
for i in range(0, 7):
    print(f'''
          {i+1}th integral
          ----------------
          a:    {2*a}
          b:    {2*b}   
          step: {h*10**(-3+i)}        
          N:    {int((2*b-2*a)/(h*10**(-3+i)))}
    ''')
    integral020 += trapezoidal_int(a, b, h*10**(-3+i), epsilon[0])
    integral200 += trapezoidal_int(a, b, h*10**(-3+i), epsilon[1])
    print(f'integ{i+1}: ', trapezoidal_int(a, b, h*10**(-3+i), epsilon[1]))
    a = b
    b *= 10
print(10*'-' + '\nBingo !' )


a_ICIR020 = -1/(1/(np.sqrt(pi))*integral020)
a_ICIR200 = -1/(1/(np.sqrt(pi))*integral200)
print(f"""     
               Results
        ---------------------
        asc/dy:
               -(0,2,0): {a_ICIR020}
               -(2,0,0): {a_ICIR200}
""")

if mode == 'C':
    np.savetxt(f'Results/ICIR_q1d_ix{int(Ix/(1e4 / Eh * to * ao**2))}_iy{int(Iy/(1e4 / Eh * to * ao**2))}_iz{int(Iz/(1e4 / Eh * to * ao**2))}_exact.dat', 
      [wx/wy, a_ICIR020, a_ICIR200])
elif mode == 'W':
    np.savetxt(f'Results/ICIR_q1d_ix{int(Ix/(1e4 / Eh * to * ao**2))}_iy{int(Iy/(1e4 / Eh * to * ao**2))}_iz{int(Iz/(1e4 / Eh * to * ao**2))}_C{C[0]}.dat', 
      [wx/wy, a_ICIR020, a_ICIR200])




