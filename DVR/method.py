__author__ = "@Tssp"
__date__ = "19/10/20"
import numpy as np
from input_DVR_3D import hbar, m, n, pot
from math import pi

def DVR_method(N, delta, k, x, Vj, wx):
    '''
    Parameters
    ----------
    N: length of the grid
    delta: step
    k: wavenumber in the j direction
    x: grid points
    V: Potential depth in the j direction
    wx: frequency in X direction

    Returns
    -------
    Energies
    '''

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
        V += Vj * 2/45 * (k*x)**6
    if n >= 4:
        V += - Vj * 1/3 * (k*x)**4
    V += Vj * (k*x)**2
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
    return E, cn