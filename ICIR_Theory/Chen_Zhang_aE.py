#!/usr/local/opt/python@3.8/bin/python3.8
# coding: utf-8
__author__ = 'T. Sánchez-Pastor'
__date__   = '26 de Julio de 2021'
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
from utils import A3D_int, B1_3D, B2_3D, separate_levels
from input_Chen_Zhang_2020 import eta_x, eta_y, eta_z, nx, ny, nz

# Energy:
Eo = 1/2*(eta_x + eta_y + eta_z)
E  = np.linspace(-7.5, 12, num=400)

print(f'''
         Parameters
      ----------------
      etax:    {eta_x}
      etay:    {eta_y}
      etaz:    {eta_z}
      nx:      {nx}    
      ny:      {ny} 
      nz:      {nz}
      Eo:      {Eo}
      E:       {E[0]}-{E[-1]}
''')
# Integration
a3D = []
for e in E:
    integral = 0
    if e>1/2*(eta_x + eta_y + eta_z):
        Lambda = np.float128(10/(e - Eo))
        beta   = np.linspace(1e-6, Lambda, num=50000, dtype=np.float128)
        A3D    = - np.float128(np.exp(beta*e)) * np.sqrt(eta_x*eta_y*eta_z/((4*pi)**3*np.float128(np.sinh(eta_x*beta))*np.float128(np.sinh(eta_y*beta))*np.float128(np.sinh(eta_z*beta)))) + 1/(4*pi*beta)**(3/2)
        integral = np.trapz(A3D)
    else:
        Lambda = 3
        beta   = np.linspace(1e-6, Lambda, num=50000)
        A3D    = - np.exp(beta*e) * np.sqrt(eta_x*eta_y*eta_z/((4*pi)**3*np.sinh(eta_x*beta)*np.sinh(eta_y*beta)*np.sinh(eta_z*beta))) + 1/(4*pi*beta)**(3/2)
        integral = np.trapz(A3D)

    J3D = sqrt(2)*4*pi*(integral + B1_3D(nx, ny, nz, eta_x, eta_y, eta_z, e, Lambda) + B2_3D(nx, ny, nz, eta_x, eta_y, eta_z, e, Lambda) + (1/(2*pi))**3/2 * 1/sqrt(2*Lambda))
    a3D.append(1/J3D)

Spectrum, level = separate_levels(a3D, E)
fig, ax = plt.subplots(figsize=(8,6))
for i in range(1, level+1):
    plt.plot(Spectrum[f'a3D_n{i}'], Spectrum[f'E_n{i}'], 'C0')
ax.set_xlabel(r'$a_{3D}/d_y$')
ax.set_ylabel(r'$E/(\hbar \omega_z)$')
#ax.set_xlim(-10, 10)
ax.set_ylim(-7.5, 12)
plt.grid()
plt.show()