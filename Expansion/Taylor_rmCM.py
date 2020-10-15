#!/usr/local/opt/python@3.8/bin/python3.8

__author__ = "@Tssp"
__date__ = "14/10/20"
import sympy as sp
from sympy import init_session
from math import pi
import numpy as np
from input_Taylor_rmCM import *

# Intensities for both particles in every direction:
####################################################
Ix1_value = Ix1/6.436409310e15 # a.u
Ix2_value = Ix2/6.436409310e15 # a.u
Iy1_value = Iy1/6.436409310e15 # a.u
Iy2_value = Iy2/6.436409310e15 # a.u
Iz1_value = Iz1/6.436409310e15 # a.u
Iz2_value = Iz2/6.436409310e15 # a.u

# Potential depths for both particles in every direction:
####################################################
Vx1_value = alpha1 * Ix1_value
Vx2_value = alpha2 * Ix2_value
Vy1_value = alpha1 * Iy1_value
Vy2_value = alpha2 * Iy2_value
Vz1_value = alpha1 * Iz1_value
Vz2_value = alpha2 * Iz2_value

# Wavelengths for both particles in every direction:
####################################################
wLx1_value  = wLx1 / 0.0529177249 # a.u
wLx2_value  = wLx2 / 0.0529177249 # a.u
wLy1_value  = wLy1 / 0.0529177249 # a.u
wLy2_value  = wLy2 / 0.0529177249 # a.u
wLz1_value  = wLz1 / 0.0529177249 # a.u
wLz2_value  = wLz2 / 0.0529177249 # a.u

# Wavenumber for both particles in every direction:
####################################################
kx1_value   = 2*pi/wLx1_value
kx2_value   = 2*pi/wLx2_value
ky1_value   = 2*pi/wLy1_value
ky2_value   = 2*pi/wLy2_value
kz1_value   = 2*pi/wLz1_value
kz2_value   = 2*pi/wLz2_value

nmax = np.max([nx1, ny1, nz1, nx2, ny2, nz2])
print(f""" 
           Atom 1                  =    {atom1}
           X-axis order (Taylor) 1 =    {nx1}
           Y-axis order (Taylor) 1 =    {ny1}
           Z-axis order (Taylor) 1 =    {nz1}
           Ix(mW/cm2) 1            =    {Ix1}
           Iy(mW/cm2) 1            =    {Iy1}
           Iz(mW/cm2) 1            =    {Iz1}
           wLx(nm) 1               =    {wLx1}
           wLy(nm) 1               =    {wLy1_value}
           wLz(nm) 1               =    {wLz1_value}
           alpha(a.u) 1            =    {alpha1}
           kx(a.u) 1               =    {kx1_value}
           ky(a.u) 1               =    {ky1_value}
           kz(a.u) 1               =    {kz1_value}
           Vx(a.u) 1               =    {Vx1_value}  
           Vy(a.u) 1               =    {Vy1_value}   
           Vz(a.u) 1               =    {Vz1_value}  
           
           Atom 2                  =    {atom2}
           X-axis order (Taylor) 2 =    {nx2}
           Y-axis order (Taylor) 2 =    {ny2}
           Z-axis order (Taylor) 2 =    {nz2}
           Ix(mW/cm2) 2            =    {Ix2}
           Iy(mW/cm2) 2            =    {Iy2}
           Iz(mW/cm2) 2            =    {Iz2}
           wLx(nm) 2               =    {wLx2_value}
           wLy(nm) 2               =    {wLy2_value}
           wLz(nm) 2               =    {wLz2_value}
           alpha(a.u) 2            =    {alpha2}
           kx(a.u) 2               =    {kx2_value}
           ky(a.u) 2               =    {ky2_value}
           kz(a.u) 2               =    {kz2_value}
           Vx(a.u) 2               =    {Vy2_value}
           Vy(a.u) 2               =    {Vx2_value}
           Vz(a.u) 2               =    {Vz2_value}
      """)

# Reduced masses:
#################
mu1_value = m1_value/(m1_value * m2_value)
mu2_value = m2_value/(m1_value * m2_value)

def opt_pot_expansion(atom1, atom2, mu1_value, mu2_value, kx1_value, kx2_value, V1_value, V2_value, nx1, nx2):
    ''' Inputs input_Taylor_rmCM.py parameters and expand the optical lattice on a given axis for a given order:
    
        Parameters
        ----------
        atom1, atom2, n:         directly given by input_Taylor_rmCM.py.
        mu1_value, mu2_value:    reduced mass coefficients of the atoms.
        k_value:                 wavenumber for the given axis.
        V1_value, V2_value:      potential depths
    '''
    x1, x2, mu1, mu2, kx1, kx2, Vx1, Vx2, x, X = sp.symbols('x_1 x_2 \mu_1 \mu_2 k^1_x k^2_x V^1_o_x V^2_o_x x X')
    dic_subs = {mu1: mu1_value,
               mu2: mu2_value,
               kx1: kx1_value,
               kx2: kx2_value,
               Vx1: Vx1_value,
               Vx2: Vx2_value}
    
    VoL = Vx1*sp.sin(kx1*x1)**2 + Vx2*sp.sin(kx2*x2)**2           # Define the sinusoidal lattice
    VoL_expanded = sp.series(VoL, x1, n=nx1+1).removeO()          # Expand in x1 up to nth order
    VoL_expanded = sp.series(VoL_expanded, x2, n=nx2+1).removeO() # Expand in x2 up to nth order
    VoL_rmCM = VoL_expanded.subs(x2, X - x*mu1)                   # Substitute x2 by the CM-rm coordinates
    VoL_rmCM = VoL_rmCM.subs(x1, X + x*mu2)                       # Substitute x1 by the CM-rm coordinates
    VoL_rmCM_expanded = sp.expand(VoL_rmCM)                       # Expand the expression
    coeff = list()
    for i in range(0, nmax+1):
        for j in range(0, nmax+1):
            coeff.append(VoL_rmCM_expanded.coeff(x, i).coeff(X, j).subs(dic_subs))
            print(VoL_rmCM_expanded.coeff(x, i).coeff(X, j).subs(dic_subs))
    return coeff

print(f'Expansion in the X axis up to {nx1}th order for the first atom and {nx2}th order for the second:\n')
coeff_x = opt_pot_expansion(atom1, atom2, mu1_value, mu2_value, kx1_value, kx2_value, Vx1_value, Vx2_value, nx1, nx2)

print(f'Expansion in the Y axis up to {ny1}th order for the first atom and {ny2}th order for the second:\n')
coeff_y = opt_pot_expansion(atom1, atom2, mu1_value, mu2_value, ky1_value, ky2_value, Vy1_value, Vy2_value, ny1, ny2)

print(f'Expansion in the Z axis up to {nz1}th order for the first atom and {nz2}th order for the second:\n')
coeff_z = opt_pot_expansion(atom1, atom2, mu1_value, mu2_value, kz1_value, kz2_value, Vz1_value, Vz2_value, nz1, nz2)

f = open(f"out/{atom1}{atom2}_nx{nx1}ny{ny1}nz{nz1}_ix{Ix1}iy{Iy1}iz{Iz1}_coeff.txt", "a")
f.write(str(atom1) + '\n')
f.write(str(atom2) + '\n')
f.write(str(alpha1) + '\n')
f.write(str(alpha1) + '\n')
f.write(str(len(coeff_x)) + '\n')
f.write(str(nx1) + '\n')
f.write(str(ny1) + '\n')
f.write(str(nz1) + '\n')
for line in coeff_x:
    f.write(str(line) + '\n')
for line in coeff_y:
    f.write(str(line) + '\n')
for line in coeff_z:
    f.write(str(line) + '\n')
f.close()

print(f"\nout/{atom1}{atom2}_nx{nx1}ny{ny1}nz{nz1}_ix{Ix1}iy{Iy1}iz{Iz1}.coeff")



# ### Archivo de salida:
# Nombre del átomo 1
# 
# Nombre del átomo 2
# 
# Polarización 1
# 
# Polarización 2
# 
# Número de términos de la expansión total
# 
# Orden de la expansión en x
# 
# "" en y
# 
# "" en z
# 
# (0, 0) coef (no se sabe si entre 2 o por 2)
# 
# (x0, X1) coef
# 
# (x0, X2) coef
# 
# .
# 
# .
# 
# .
# 
# (x0, Xn) coef
# 
# (x1, X0) coef
# 
# (x1, X1) coef
# 
# .
# 
# .
# 
# .
# 
# (xn, Xn)
# 
# "" en y
# 
# "" en z
# 
