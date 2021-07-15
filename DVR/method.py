__author__ = "@Tssp"
__date__ = "19/10/20"
import numpy as np
from input_DVR_3D import hbar, n, pot, mode
from math import pi
import matplotlib.pyplot as plt

def DVR_method(N, delta, m, k, x, Vj, w, mode='all'):
    '''
    Parameters
    ----------
    N: length of the grid
    delta: step
    k: wavenumber in the j direction
    x: grid points
    V: Potential depth in the j direction
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
            if i == j:
                T[i-1, i-1] = -hbar**2/(2*m) * (-1/3*(pi/delta)**2 + 2/delta**2 * (-1)**(i+j)/(i+j)**2)
            elif i != j:
                T[i-1, j-1] = hbar**2/(delta**2*m) * ((-1)**(i-j)/(i-j)**2 - (-1)**(i+j)/(i+j)**2)
                
    #print('Kinetic operator\n', T[-3:, -3:], '\n')
    # Potential:
    ############
    V = 0
    if mode == 'all':
        if n >= 6:
            V += Vj * 2/45 * (k*x)**6
            print(Vj * 2/45 * k**6, 'x^6')
        if n >= 4:
            V += - Vj * 1/3 * (k*x)**4
            print(Vj * 1/3 * k**4, 'x^4')
        V += Vj * (k*x)**2
        print(Vj * k**2, 'x^2')
        #V = Vx * (2/45 * (kx*x)**6 - 1/3 * (kx*x)**4 + (kx*x)**2)
    elif mode == 'CM':
        print('Computing only CM energies')
        if n >= 6:
            V += Vj * 4/45 * (k*x)**6
            print(Vj * 4/45 * k**6, 'x^6')
        if n >= 4:
            V += - Vj * 2/3 * (k*x)**4
            print(- Vj * 2/3 * k**4, 'x^4')
        V += Vj * 2*(k*x)**2
        print(Vj * 2*k**2, 'x^2')
    #V = Vx * (4/45 * (kx*x)**6 - 2/3 * (kx*x)**4 + 2*(kx*x)**2)
    if pot == 'cos2':
        V = 1 - V
    # Hamiltonian:
    ##############
    H = T + np.diagflat(V)
    print(f'\n\nx:\n{x[-3:]}\n\n')
    print(f'\n\nKinetic Energy:\n{T[-3:,-3:]}\n\n')
    print(f'\n\nPotential:\n{np.diagflat(V)[-3:,-3:]}\n\n')
    print(f'\n\nHamiltonian:\n{H[-3:,-3:]}\n\n')
    # Diagonalization:
    ##################
    E, cn = np.linalg.eigh(H)
    del T, H
    inds = E.argsort()
    E = E[inds[::1]]
    cn = cn[:,inds[::1]]
    print("           n        E                       E[hbar wx]")
    for i in range(11):
        print('          ', i, E[i],'   ', E[i]/w)
    print("\nBingo !")
    return E, cn