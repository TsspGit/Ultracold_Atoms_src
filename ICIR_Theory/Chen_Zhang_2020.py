#!/usr/local/opt/python@3.8/bin/python3.8
# coding: utf-8
__author__ = 'T. Sánchez-Pastor'
__date__   = '19 de Julio de 2021'
# Modules
import numpy as np
import matplotlib.pyplot as plt
from math import pi, sqrt, isnan
from scipy.special import gamma, hyp2f1
import os
from atomic_units import ao, vo, e, hbar, me, Eh, to
from utils import A3D_int, B1_3D, B2_3D
from input_Chen_Zhang_2020 import eta_x, eta_y, eta_z, nx, ny, nz

# Energy:
E      = 2.5

print(f'''
         Parameters
      ----------------
      etax:    {eta_x}
      etay:    {eta_y}
      etaz:    {eta_z}
      nx:      {nx}    
      ny:      {ny} 
      nz:      {nz}
      Eo:      {1/2*(eta_x + eta_y + eta_z)}
      E:       {E}
''')
# Integration
h        = 1e-6
a        = h*1e-3 
b        = h*1e4
integral_list = []
for i in range(0, 7):
    print(f'''
          {i+1}th integral
          ----------------
          a:    {2*a}
          b:    {2*b}   
          step: {h*10**(-3+i)}        
          N:    {int((2*b-2*a)/(h*10**(-3+i)))}
    ''')
    integral_list.append(A3D_int(a, b, h*10**(-3+i), eta_x, eta_y, eta_z, E))
    print(integral_list[-1])
    print('\nDone!\n')
    a = b
    b *= 10
    if abs(integral_list[-1]) > 1e3 or isnan(integral_list[-1]):
      print('hola')
      integral_list[-1] = 0
      break
A3D = np.sum(integral_list)
print(10*'-' + '\nBingo !' )


print(f'''B1_3D: {B1_3D(nx, ny, nz, eta_x, eta_y, eta_z, E, 2*b)},
B2_3D: {B2_3D(nx, ny, nz, eta_x, eta_y, eta_z, E, 2*b)}''')

J3D = sqrt(2)*4*pi*(A3D + B1_3D(nx, ny, nz, eta_x, eta_y, eta_z, E, 2*b) + B2_3D(nx, ny, nz, eta_x, eta_y, eta_z, E, 2*b) +             (1/(2*pi))**3/2 * 1/sqrt(2*2*b))
print(f"""     
               Results
        ---------------------
        asc/dy: {1/J3D}
""")