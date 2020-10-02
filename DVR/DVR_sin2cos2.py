#!/usr/local/opt/python@3.8/bin/python3.8

__author__ = "@Tssp"
__date__ = "01/10/20"
import numpy as np
import matplotlib.pyplot as plt
from math import pi
from input_DVR_sin2cos2 import hbar, m, Ix, n, pot, wL, alpha, delta, xmax, xmin, JOBZ
plt.rc('text',usetex=True)
plt.rc('font',family='serif')
plt.rcParams['xtick.labelsize']=16
plt.rcParams['ytick.labelsize']=16
plt.rcParams['axes.labelsize']=18
plt.rcParams['axes.titlesize']=18

# Parameters:
Ix = Ix / 6.436409310e15
wL = wL / 0.0529177249
kx = 2*pi/wL
Vx = alpha * Ix
wx = np.sqrt(2 * Vx/m) * kx
print(f""" 
           Potential      =    {pot}
           order (Taylor) =    {n}
           Ix(mW/cm2)     =    {Ix * 6.436409310e15}
           Ix(a.u)        =    {Ix}
           wL(nm)         =    {wL * 0.0529177249}
           wL(a.u)        =    {wL}
           hbar           =    {hbar}
           alpha(a.u)     =    {alpha}
           delta          =    {delta}
           xmin           =    {xmin}
           xmax           =    {xmax}
           kx(a.u)        =    {kx}
           Vx(a.u)        =    {Vx}
           wx(a.u)        =    {wx}
           
      """)
# xaxis:
########
x = np.arange(xmin, xmax, delta)
N = len(x)

# Kinetic Energy:
#################
T = np.zeros((N, N))
for i in range(N):
    i += 1 
    for j in range(N):
        j += 1
        if i == j:
            T[i-1, i-1] = -hbar**2/(2*m) * (-1/3*(pi/delta)**2 + 2/delta**2 * (-1)**(i+j)/(i+j)**2)
        elif i != j:
            T[i-1, j-1] = hbar**2/(delta**2*m) * ((-1)**(i-j)/(i-j)**2 - (-1)**(i+j)/(i+j)**2)
            
# Potential:
############
V = 0
if n >= 6:
    V += Vx * 2/45 * (kx*x)**6
if n >= 4:
    V += - Vx * 1/3 * (kx*x)**4
V += Vx * (kx*x)**2
#V = Vx * (2/45 * (kx*x)**6 - 1/3 * (kx*x)**4 + (kx*x)**2)
if pot == 'cos2':
    V = 1 - V
# Hamiltonian:
##############
H = T + np.diagflat(V)
# Diagonalization:
##################
E, cn = np.linalg.eig(H)
del T, H
inds = E.argsort()
E = E[inds[::1]]
cn = cn[:,inds[::1]]
print("           n        E                       E[hbar wx]")
for i in range(11):
    print('          ', i, E[i],'   ', E[i]/wx)
print("\nBingo !")

if JOBZ == 'V':
    i = int(input("Level to plot: "))
    fig, ax = plt.subplots(figsize=(8,6))
    eigen_bool = int(input("Type [0] to plot the eigenfunctions or [1] probability density: "))
    if eigen_bool == 0:
        ax.plot(x, -cn[:,i], label='$\psi$')
    elif eigen_bool == 1:
        ax.plot(x, cn[:,i] * np.conj(cn[:,i]), label='$|\psi|^2$')
        ax.set_ylim(0, 2*np.max(cn[:,i] * np.conj(cn[:,i])))
    ax.plot(x, V*1e7, 'k', lw=2, label='V(x)')
    ax.set_xlabel("$x$")
    ax.set_ylabel("$\psi(x)$")
    plt.legend()
    plt.grid()
    plt.show()
    print(10*"-"+ "Routine finished"+10*"-")
elif JOBZ == 'N':
    print(10*"-"+ "Routine finished"+10*"-")