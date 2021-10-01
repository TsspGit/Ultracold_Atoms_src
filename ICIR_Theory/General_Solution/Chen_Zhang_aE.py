#!/usr/local/opt/python@3.8/bin/python3.8
# coding: utf-8
__author__ = 'T. Sánchez-Pastor'
__date__   = '26 de Julio de 2021'
# Modules
import numpy as np
import matplotlib.pyplot as plt
from math import pi, sqrt, isnan
from scipy.special import gamma, hyp2f1
import os
from atomic_units import ao, vo, e, hbar, me, Eh, to
from utils import A3D_int, A3D, I3D, W3D, B1_3D, B2_3D, separate_levels, gauss
from input_Chen_Zhang_2020 import eta_x, eta_y, nx, ny
from scipy.integrate import quad
plt.rc('text',usetex=True)
plt.rc('font',family='serif')
ref_ticksize = 16
plt.rcParams['xtick.labelsize']=ref_ticksize
plt.rcParams['legend.fontsize']=ref_ticksize
plt.rcParams['ytick.labelsize']=ref_ticksize
plt.rcParams['axes.labelsize']=ref_ticksize * 3/2
plt.rcParams['axes.titlesize']=ref_ticksize * 3/2

# Energy:
Eo = 1/2*(eta_x + eta_y + 1)
E  = np.linspace(-7.5, 13, num=1000)

print(f'''
         Parameters
      ----------------
      etax:    {eta_x}
      etay:    {eta_y}
      nx:      {nx}    
      ny:      {ny} 
      Eo:      {Eo}
      E:       {E[0]}-{E[-1]}
''')
# Integration
a3D = []
for e in E:
    integral = 0
    if e>Eo: #and e<=6:
        Lambda = 0.15/(e - Eo)
        #Lambda = 0.86
        integral = gauss(lambda x: I3D(eta_x, eta_y, nx, ny, e, x), 100, 1e-6, Lambda)
    #elif e>6:
    #    Lambda = 1.1
    #    integral = gauss(lambda x: I3D(eta_x, eta_y, nx, ny, e, x), 100, 1e-6, Lambda)
    else:
        Lambda = 100
        integral = gauss(lambda x: I3D(eta_x, eta_y, nx, ny, e, x), 100, 1e-6, Lambda)
    J3D = sqrt(2)*4*pi*(W3D(nx, ny, eta_x, eta_y, e) + integral)
    a3D.append(1/J3D)

#fig, ax = plt.subplots(figsize=(8,6))
#plt.plot(a3D, E, 'C0')
#ax.set_xlabel(r'$a_{3D}/d_y$')
#ax.set_ylabel(r'$E/(\hbar \omega_z)$')
#ax.set_xlim(-10, 10)
#ax.set_ylim(-4, 8)
#plt.grid()
#plt.savefig('E_as_anisotropic.png', dpi=200)
#plt.show()

Spectrum, level = separate_levels(a3D, E)
fig, ax = plt.subplots(figsize=(8,6))
for i in range(1, level+1):
    plt.plot(Spectrum[f'a3D_n{i}'], Spectrum[f'E_n{i}'], 'C0')
ax.set_xlabel(r'$a_{3D}/d_y$')
ax.set_ylabel(r'$E/(\hbar \omega_z)$')
ax.set_xlim(-10, 10)
ax.set_ylim(-7.5, 12)
plt.grid()
plt.savefig('E_as_isotropic.png', dpi=200)
plt.show()