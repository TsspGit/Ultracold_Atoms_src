__author__ = "@Tssp"
__date__ = "08/04/21"
import numpy as np
from input_DVR_pot import *
from math import pi
import matplotlib.pyplot as plt

def DVR_method(N, delta, m1, m2, k1, k2, x, coeff, w, mode):
    '''
    Parameters
    ----------
    N: length of the grid
    delta: step
    x: grid points
    k: wavenumbers
    V: Potential depth in the j direction
    coeff: coefficients of the Taylor expansion
    w: frequency in X direction
    mode: all for rm+cm contributions or CM for center of mass computation.
    Returns
    -------
    Energies.
    Wavefunctions.
    '''

    # Kinetic Energy:
    #################
    T = np.zeros((N, N))
    for i in range(N):
        i += 1 
        for j in range(N):
            j += 1
            if mode == 'all':
                if i == j:
                    T[i-1, i-1] = -hbar**2/(2*m1) * (-1/3*(pi/delta)**2 + 2/delta**2 * (-1)**(i+j)/(i+j)**2) - hbar**2/(2*m2) * (-1/3*(pi/delta)**2 + 2/delta**2 * (-1)**(i+j)/(i+j)**2)
                elif i != j:
                    T[i-1, j-1] = hbar**2/(delta**2*m2) * ((-1)**(i-j)/(i-j)**2 - (-1)**(i+j)/(i+j)**2) + hbar**2/(delta**2*m2) * ((-1)**(i-j)/(i-j)**2 - (-1)**(i+j)/(i+j)**2)
            elif mode == 'CM':
                m = m1 + m2
                if i == j:
                    T[i-1, i-1] = -hbar**2/(2*m) * (-1/3*(pi/delta)**2 + 2/delta**2 * (-1)**(i+j)/(i+j)**2)
                elif i != j:
                    T[i-1, j-1] = hbar**2/(delta**2*m) * ((-1)**(i-j)/(i-j)**2 - (-1)**(i+j)/(i+j)**2)
            
    # Potential:
    ############
    V = 0
    if mode == 'all':
        V1 = 0
        V2 = 0
        for i in range(len(coeff[0]), 0, -1):
            V1 += coeff[0][i-1] * (x)**(i-1)
            V2 += coeff[1][i-1] * (x)**(i-1)
            print(coeff[0][i-1], 'x^', i-1)
        H = T + np.diagflat(V1) + np.diagflat(V2)
    elif mode == 'CM':
        print('Computing only CM energies')
        V = 0
        for i in range(len(coeff), 0, -1):
            V += coeff[i-1] * x**(i-1)
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
        print('          ', i, E[i],'   ', E[i]/w)
    print("\nBingo !")
    return E, cn