#!/usr/local/opt/python@3.8/bin/python3.8
__author__ = "@Tssp"
__date__ = "12/01/21"
import numpy as np
from math import pi
import os
from numba import njit, prange
from atomic_units import ao, vo, e, hbar, me, Eh, to
from input_ICIR_Theory_main import *
from utils import *

#--------------------------------------- Parameters and Constants ---------------------------------------
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
           model          =    {model}
           Config         =    {Config}
      """)

#--------------------------------------- Energies ---------------------------------------

if model == 'Numerical':
  if Config == False:
    # Erm and ECM actually readed from input file
    Eref = Erm + ECM
  elif Config == True:
    Eref = Econfig

elif model == 'Perturbation':
  ECM = Ecm_n(wx, wy, wz, Vx, Vy, Vz, 0, 0, 0)
  Erm = ECM
  En_CM = Ecm_n(wx, wy, wz, Vx, Vy, Vz, 0, 4, 0)
  Eref = Erm + ECM

E_ICIR = wy*np.array(E_ICIR)
Eo = (wx + wy + wz)/2

print(f"""
          ECM:
              - (0,0,0): {ECM}
              - n-level: {En_CM}
          Erm (1):       {Erm} 
          E_ICIR:        {E_ICIR}
          Eo:            {Eo}
  """)

#--------------------------------------- C and Epsilon ---------------------------------------

if mode == 'C':
  C = (E_ICIR - Eref)/wz
elif mode == 'W':
  # C readed from the input file
  pass

epsilon = (C*wz + 2 * ECM - En_CM - Eo)/wy
print(f"""
          C:       {C}
          eps:     {epsilon}
""")

#--------------------------------------- Integral ---------------------------------------
a = h*1e-3
b = h*1e5
integral = 0
for i in range(0, 7):
    print(f'''
          {i+1}th integral
          ----------------
          a:    {2*a}
          b:    {2*b}   
          step: {h*10**(-3+i)}        
          N:    {int((2*b-2*a)/(h*10**(-3+i)))}
    ''')
    integral += trapezoidal_int_v2(a, b, h*10**(-3+i), epsilon, eta_x, eta_z)/sqrt(2)# if v2 /sqrt(2)
    print(integral)
    print('\nDone!\n')
    a = b
    b *= 10
print(10*'-' + '\nBingo !' )


a_ICIR = -1/(1/(np.sqrt(pi))*integral)
print(f"""     
               Results
        ---------------------
        asc/dy: {a_ICIR}
""")

#-------------------------------------------------------------------------------------------------