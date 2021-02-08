#!/usr/local/opt/python@3.8/bin/python3.8

__author__ = "@Tssp"
__date__ = "03/02/21"
import numpy as np
import matplotlib.pyplot as plt
from math import pi
from input_DVR_3D import hbar, m, Ix, Iy, Iz, n, pot, wL, alpha, delta, xmax, xmin, mode
from method import DVR_method
plt.rc('text',usetex=True)
plt.rc('font',family='serif')
plt.rcParams['xtick.labelsize']=16
plt.rcParams['ytick.labelsize']=16
plt.rcParams['axes.labelsize']=18
plt.rcParams['axes.titlesize']=18

# Parameters:
Ix = Ix / 6.436409310e15
Iy = Iy / 6.436409310e15
Iz = Iz / 6.436409310e15
wL = wL / 0.0529177249
kx = 2*pi/wL
Vx = alpha * Ix
Vy = alpha * Iy
Vz = alpha * Iz
if mode == 'CM':
    m = 2*m
    wx = np.sqrt(4 * Vx/m) * kx
elif mode == 'all':
    wx = np.sqrt(2 * Vx/m) * kx


print(f""" 
           Potential      =    {pot}
           order (Taylor) =    {n}
           mass (a.u)     =    {m}
           Ix(mW/cm2)     =    {Ix * 6.436409310e15}
           Iy(mW/cm2)     =    {Iy * 6.436409310e15}
           Iz(mW/cm2)     =    {Iz * 6.436409310e15}
           wL(nm)         =    {wL * 0.0529177249}
           wL(a.u)        =    {wL}
           hbar           =    {hbar}
           alpha(a.u)     =    {alpha}
           delta          =    {delta}
           xmin           =    {xmin}
           xmax           =    {xmax}
           kx(a.u)        =    {kx}
           Vx(a.u)        =    {Vx}
           Vy(a.u)        =    {Vy}
           Vz(a.u)        =    {Vz}
           wx(a.u)        =    {wx}
           
      """)
# xaxis:
########
x = np.arange(xmin, xmax, delta)
N = len(x)
print("\nX Axis: \n")
Ex, _ = DVR_method(N, delta, m, kx, x, Vx, wx, mode)

print("\nY Axis: \n")
Ey, _ = DVR_method(N, delta, m, kx, x, Vy, wx, mode)

print("\nZ Axis: \n")
Ez, _ = DVR_method(N, delta, m, kx, x, Vz, wx, mode)

Et = Ex + Ey + Ez
print(f'The corresponding level is:')
while True:
  E_input = float(input('Introduce the energy of the level: '))
  for i in range(5):
      if i%2 == 0:
          for j in range(5):
              if j%2 == 0:
                  for k in range(30):
                      if k%2 == 0:
                          if mode == 'all':
                            if E_input == 2*(Ex[i] + Ey[j] + Ez[k]):
                                print(f'(nx,ny,nz) = ({i}, {j}, {k})')
                          elif mode == 'CM':
                              if float("{:.4e}".format(E_input)) == float("{:.4e}".format(Ex[i] + Ey[j] + Ez[k])):
                                print(f'(nx,ny,nz) = ({i},{j},{k})')
print('\nIf the program did not show the corresponding level of the inputed energy, please try again.')