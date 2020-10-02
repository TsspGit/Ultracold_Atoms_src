#!/usr/local/opt/python@3.8/bin/python3.8

__author__ = "@Tssp"
__date__ = "01/10/20"
import numpy as np
import matplotlib.pyplot as plt
from math import pi
from input_DVR import hbar, m, w, delta, xmax, xmin, JOBZ
plt.rc('text',usetex=True)
plt.rc('font',family='serif')
plt.rcParams['xtick.labelsize']=16
plt.rcParams['ytick.labelsize']=16
plt.rcParams['axes.labelsize']=18
plt.rcParams['axes.titlesize']=18

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
V = 1/2 * m * w**2 * x**2
# Hamiltonian:
##############
H = T + np.diagflat(V)
# Diagonalization:
##################
E, cn = np.linalg.eig(H)
inds = E.argsort()
E = E[inds[::1]]
cn = cn[:,inds[::1]]
print("n        E")
for i in range(11):
    print(i, E[i])
print("\nBingo !")

if JOBZ == 'V':
    i = int(input("Level to plot: "))
    fig, ax = plt.subplots(figsize=(8,6))
    ax.plot(x, -cn[:,i])
    ax.set_xlabel("$x$")
    ax.set_ylabel("$\psi(x)$")
    plt.grid()
    plt.show()
    print(10*"-"+ "Routine finished"+10*"-")
elif JOBZ == 'N':
    print(10*"-"+ "Routine finished"+10*"-")