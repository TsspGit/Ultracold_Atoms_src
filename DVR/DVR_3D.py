#!/usr/local/opt/python@3.8/bin/python3.8

__author__ = "@Tssp"
__date__ = "01/10/20"
import numpy as np
import matplotlib.pyplot as plt
from math import pi
from input_DVR_3D import hbar, m, Ix, Iy, Iz, n, pot, wL, alpha, delta, xmax, xmin
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
wx = np.sqrt(2 * Vx/m) * kx
print(f""" 
           Potential      =    {pot}
           order (Taylor) =    {n}
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
Ex, _ = DVR_method(N, delta, kx, x, Vx, wx)

print("\nY Axis: \n")
Ey, _ = DVR_method(N, delta, kx, x, Vy, wx)

print("\nZ Axis: \n")
Ez, _ =DVR_method(N, delta, kx, x, Vz, wx)

print("\nEx + Ey + Ez: \n")
Et = Ex + Ey + Ez
print("           n        E                       E[hbar wx]")
for i in range(11):
    print('          ', i, Et[i],'   ', Et[i]/wx)
print("\nBingo !")









