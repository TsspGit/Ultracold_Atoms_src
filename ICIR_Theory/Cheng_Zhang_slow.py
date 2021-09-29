# coding: utf-8
__author__ = 'T. Sánchez-Pastor'
__date__   = '23 de Septiembre de 2021'
# Modules
import numpy as np
import matplotlib.pyplot as plt
plt.rc('text',usetex=True)
plt.rc('font',family='serif')
ref_ticksize = 16
plt.rcParams['xtick.labelsize']=ref_ticksize
plt.rcParams['legend.fontsize']=ref_ticksize
plt.rcParams['ytick.labelsize']=ref_ticksize
plt.rcParams['axes.labelsize']=ref_ticksize * 3/2
plt.rcParams['axes.titlesize']=ref_ticksize * 3/2
from math import pi, sqrt, isnan
from scipy.special import gamma, hyp2f1
import os
from atomic_units import ao, vo, e, hbar, me, Eh, to
from utils import A3D_int, I3D_int, W3D, separate_levels
from input_Chen_Zhang_2020 import eta_x, eta_y, eta_z, nx, ny

# Energy:
Eo = 1/2*(eta_x + eta_y + eta_z)
E  = np.linspace(-7.5, 5, num=100)

print(f'''
         Parameters
      ----------------
      etax:    {eta_x}
      etay:    {eta_y}
      etaz:    {eta_z}
      nx:      {nx}    
      ny:      {ny} 
      Eo:      {Eo}
      E:       {E[0]}-{E[-1]}
''')
# Integration
a3D = []
integral_list = []
for e in E:
    integral = 0
    h = 1e-6
    a = h*1e-3
    b = h*1e-1
    for i in range(0, 10):
        print(f'''
              {i+1}th integral
              ----------------
              a:    {2*a}
              b:    {2*b}   
              step: {h*10**(-3+i)}        
              N:    {int((2*b-2*a)/(h*10**(-3+i)))}
              E:    {e}
        ''')
        integral += I3D_int(a, b, h*10**(-3+i), nx, ny, eta_x, eta_y, eta_z, e)
        print(integral)
        print('\nDone!\n')
        a = b
        b *= 10
        integral_list.append(integral)
        print(10*'-' + '\nBingo !' )

    J3D = sqrt(2)*4*pi*(W3D(nx, ny, eta_x, eta_y, eta_z, e) + integral)
    a3D.append(1/J3D)

print(a3D)
Spectrum, level = separate_levels(a3D, E)
fig, ax = plt.subplots(figsize=(8,6))
for i in range(1, level+1):
    plt.plot(Spectrum[f'a3D_n{i}'], Spectrum[f'E_n{i}'], 'C0')
ax.set_xlabel(r'$a_{3D}/d_y$')
ax.set_ylabel(r'$E/(\hbar \omega_z)$')
ax.set_xlim(-10, 10)
ax.set_ylim(-7.5, 12)
plt.grid()
plt.show()