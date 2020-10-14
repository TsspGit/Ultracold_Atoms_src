#!/usr/local/opt/python@3.8/bin/python3.8

__author__ = "@Tssp"
__date__ = "14/10/20"
import sympy as sp
from sympy import init_session
from math import pi
from input_Taylor_rmCM import atom1, atom2, m1_value, m2_value, nx, ny, nz, Ix1, Ix2, Iy1, Iy2, Iz1, Iz2, wL, alpha

Ix1_value = Ix1/6.436409310e15 # a.u
Ix2_value = Ix2/6.436409310e15 # a.u
Iy1_value = Iy1/6.436409310e15 # a.u
Iy2_value = Iy2/6.436409310e15 # a.u
Iz1_value = Iz1/6.436409310e15 # a.u
Iz2_value = Iz2/6.436409310e15 # a.u
Vx1_value = alpha * Ix1_value
Vx2_value = alpha * Ix2_value
Vy1_value = alpha * Iy1_value
Vy2_value = alpha * Iy2_value
Vz1_value = alpha * Iz1_value
Vz2_value = alpha * Iz2_value
wL_value  = wL / 0.0529177249 # a.u
kx_value  = 2*pi/wL
print(f""" 
           atom 1                =    {atom1}
           atom 2                =    {atom2}
           X-axis order (Taylor) =    {nx}
           Y-axis order (Taylor) =    {ny}
           Z-axis order (Taylor) =    {nz}
           Ix(mW/cm2) 1          =    {Ix1}
           Ix(mW/cm2) 2          =    {Ix2}
           Iy(mW/cm2) 1          =    {Iy1}
           Iy(mW/cm2) 2          =    {Iy2}
           Iz(mW/cm2) 1          =    {Iz1}
           Iz(mW/cm2) 2          =    {Iz2}
           wL(nm)                =    {wL}
           alpha(a.u)            =    {alpha}
           kx(a.u)               =    {kx_value}
           Vx(a.u) 1             =    {Vx1_value}   
           Vx(a.u) 2             =    {Vy2_value}
           Vy(a.u) 1             =    {Vy1_value}   
           Vy(a.u) 2             =    {Vx2_value}
           Vz(a.u) 1             =    {Vz1_value}   
           Vz(a.u) 2             =    {Vz2_value}
      """)
mu1_value = m1_value/(m1_value * m2_value)
mu2_value = m2_value/(m1_value * m2_value)

def opt_pot_expansion(atom1, atom2, mu1_value, mu2_value, k_value, V1_value, V2_value, n):
    ''' Inputs input_Taylor_rmCM.py parameters and expand the optical lattice on a given axis for a given order:
    
        Parameters
        ----------
        atom1, atom2, n:         directly given by input_Taylor_rmCM.py.
        mu1_value, mu2_value:    reduced mass coefficients of the atoms.
        k_value:                 wavenumber for the given axis.
        V1_value, V2_value:      potential depths
    '''
    x1, x2, mu1, mu2, kx, Vx1, Vx2, x, X = sp.symbols('x_1 x_2 \mu_1 \mu_2 k_x V^1_o_x V^2_o_x x X')
    x1, x2, mu1, mu2, kx, Vx1, Vx2, x, X = sp.symbols('x_1 x_2 \mu_1 \mu_2 k_x V^1_o_x V^2_o_x x X')
    dic_subs = {mu1: mu1_value,
               mu2: mu2_value,
               kx: kx_value,
               Vx1: Vx1_value,
               Vx2: Vx2_value}
    VoL = Vx1*sp.sin(kx*x1)**2 + Vx2*sp.sin(kx*x2)**2           # Define the sinusoidal lattice
    VoL_expanded = sp.series(VoL, x1, n=n+1).removeO()          # Expand in x1 up to nth order
    VoL_expanded = sp.series(VoL_expanded, x2, n=n+1).removeO() # Expand in x2 up to nth order
    VoL_rmCM = VoL_expanded.subs(x2, X - x*mu1)                 # Substitute x2 by the CM-rm coordinates
    VoL_rmCM = VoL_rmCM.subs(x1, X + x*mu2)                     # Substitute x1 by the CM-rm coordinates
    VoL_rmCM_expanded = sp.expand(VoL_rmCM)                     # Expand the expression
    for i in range(0, 7):
        for j in range(0, 7):
            print(VoL_rmCM_expanded.coeff(x, i).coeff(X, j).subs(dic_subs))

print(f'Expansion in the X axis up to {nx}th order:\n')
opt_pot_expansion(atom1, atom2, mu1_value, mu2_value, kx_value, Vx1_value, Vx2_value, nx)

print(f'Expansion in the Y axis up to {ny}th order:\n')
opt_pot_expansion(atom1, atom2, mu1_value, mu2_value, kx_value, Vy1_value, Vy2_value, ny)

print(f'Expansion in the Z axis up to {nz}th order:\n')
opt_pot_expansion(atom1, atom2, mu1_value, mu2_value, kx_value, Vz1_value, Vz2_value, nz)


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
# (0, 0) coef (no se sabe si entre 2 o por 2) Nota: todos los x^n·X^m n + m != 6, son 0!
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
