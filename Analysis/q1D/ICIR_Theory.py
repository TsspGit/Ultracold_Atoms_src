#!/usr/local/opt/python@3.8/bin/python3.8
__author__ = "@Tssp"
__date__ = "18/11/20"
import numpy as np
from math import pi
import os
from numba import njit, prange
from utils.atomic_units import ao, vo, e, hbar, me, Eh, to
from input_ICIR_Theory import *
from utils.Theory_utils import *

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

if model == 'Numerical':
  #Load CM and rm energies:
  folder_path = f'Simulations/ix{int(Ix/(1e4 / Eh * to * ao**2))}_iy{int(Iy/(1e4 / Eh * to * ao**2))}_iz{int(Iz/(1e4 / Eh * to * ao**2))}/orbitals/eva/'
  for file in os.listdir(folder_path + 'rm'):
      if 'noint' in file:
          print(f'rm/ folder, file readed:\n{file}')
          Erm = np.loadtxt(folder_path + 'rm/' + file)[0,2]
  for file in os.listdir(folder_path + 'CM'):
      print(f'\nCM/ folder, file readed:\n{file}')
      ECM = np.loadtxt(folder_path + 'CM/' + file)[0,2]

elif model == 'Perturbation':
      ECM = Ecm_n(wx, wy, wz, Vx, Vy, Vz, 0, 0, 0)
      En_CM020 = Ecm_n(wx, wy, wz, Vx, Vy, Vz, 0, 2, 0)
      En_CM200 = Ecm_n(wx, wy, wz, Vx, Vy, Vz, 2, 0, 0)
      En_CM = [En_CM020, En_CM200]
      Erm = Ecm_n(wx, wy, wz, Vx, Vy, Vz, 0, 0, 0)
      #Erm = ECM + wz # Simon C=1
E_ICIR = wy*np.loadtxt(f'Results/ICIR_positions_{int(Ix/(1e4 / Eh * to * ao**2))}_{int(Iy/(1e4 / Eh * to * ao**2))}_{int(Iz/(1e4 / Eh * to * ao**2))}.txt')[2]
Eo = (wx + wy + wz)/2

print(f"""
          ECM:
              - (0,0,0): {ECM}
              - (0,2,0): {En_CM[0]}
              - (2,0,0): {En_CM[1]}
          Erm (1):       {Erm} 
          E_ICIR:
              -(0,2,0):  {E_ICIR[0]}
              -(2,0,0):  {E_ICIR[1]}
          Eo:            {Eo}
  """)

# Calculations
Eref = Erm + ECM
if mode == 'C':
    C = (E_ICIR - Eref)/wz
    if Config == True:
      C = (E_ICIR - Econfig)/wz
    #C = (E_ICIR - 1.134854998141340E-010)/wz # Usando la configuracion
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

# Create the partitions of the integral axis.
a = h*1e-3
b = h*1e5
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
    integral020 += trapezoidal_int(a, b, h*10**(-3+i), epsilon[0], eta_x, eta_z)
    integral200 += trapezoidal_int(a, b, h*10**(-3+i), epsilon[1], eta_x, eta_z)
    print('\nDone!\n')
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

#if mode == 'C':
#    if Config == True:
#      np.savetxt(f'Results/ICIR_q1d_ix{int(Ix/(1e4 / Eh * to * ao**2))}_iy{int(Iy/(1e4 / Eh * to * ao**2))}_iz{int(Iz/(1e4 / Eh * to * ao**2))}_config.dat', 
#      [wx/wy, a_ICIR020, a_ICIR200], header='wxwy, asc020, asc200')
#      print(f'ICIR_q1d_ix{int(Ix/(1e4 / Eh * to * ao**2))}_iy{int(Iy/(1e4 / Eh * to * ao**2))}_iz{int(Iz/(1e4 / Eh * to * ao**2))}_config.dat')
#    else:
#      np.savetxt(f'Results/ICIR_q1d_ix{int(Ix/(1e4 / Eh * to * ao**2))}_iy{int(Iy/(1e4 / Eh * to * ao**2))}_iz{int(Iz/(1e4 / Eh * to * ao**2))}_exact.dat', 
#        [wx/wy, a_ICIR020, a_ICIR200], header='wxwy, asc020, asc200')
#      print(f'ICIR_q1d_ix{int(Ix/(1e4 / Eh * to * ao**2))}_iy{int(Iy/(1e4 / Eh * to * ao**2))}_iz{int(Iz/(1e4 / Eh * to * ao**2))}_exact.dat')
#elif mode == 'W':
#    if Config == True:
#      np.savetxt(f'Results/ICIR_q1d_ix{int(Ix/(1e4 / Eh * to * ao**2))}_iy{int(Iy/(1e4 / Eh * to * ao**2))}_iz{int(Iz/(1e4 / Eh * to * ao**2))}_config_C{C[0]}.dat', 
#        [wx/wy, a_ICIR020, a_ICIR200], header='wxwy, asc020, asc200')
#      print(f'ICIR_q1d_ix{int(Ix/(1e4 / Eh * to * ao**2))}_iy{int(Iy/(1e4 / Eh * to * ao**2))}_iz{int(Iz/(1e4 / Eh * to * ao**2))}_config_C{C[0]}.dat')
#    else:
#      np.savetxt(f'Results/ICIR_q1d_ix{int(Ix/(1e4 / Eh * to * ao**2))}_iy{int(Iy/(1e4 / Eh * to * ao**2))}_iz{int(Iz/(1e4 / Eh * to * ao**2))}_C{C[0]}.dat', 
#        [wx/wy, a_ICIR020, a_ICIR200], header='wxwy, asc020, asc200')
#      print(f'ICIR_q1d_ix{int(Ix/(1e4 / Eh * to * ao**2))}_iy{int(Iy/(1e4 / Eh * to * ao**2))}_iz{int(Iz/(1e4 / Eh * to * ao**2))}_C{C[0]}.dat')




