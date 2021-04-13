__author__ = "@Tssp"
__date__   = "22/10/20"
import sympy as sp
import numpy as np

def opt_pot_expansion(atom1, atom2, mu1_value, mu2_value, k1_value, k2_value, V1_value, V2_value, n1, n2):
    ''' Inputs input_Taylor_rmCM.py parameters and expand the optical lattice on a given axis for a given order:
    
        Parameters
        ----------
        atom1, atom2, n:         directly given by input_Taylor_rmCM.py.
        mu1_value, mu2_value:    reduced mass coefficients of the atoms.
        k_value:                 wavenumber for the given axis.
        V1_value, V2_value:      potential depths
    
        Output
        ------
        coeff:                   rm, CM coefficients of the expansion.
        coeff_CM:                CM coefficients of the expansion.
        coeff_config:            dictionary of configuration coefficients.
    '''
    
    nmax = np.max([n1, n2])
    x1, x2, mu1, mu2, k1, k2, V1, V2, x, X = sp.symbols('x_1 x_2 \mu_1 \mu_2 k^1 k^2 V^1_o V^2_o x X')
    dic_subs = {mu1: mu1_value,
               mu2: mu2_value,
               k1: k1_value,
               k2: k2_value,
               V1: V1_value,
               V2: V2_value}
    
    VoL = V1*sp.sin(k1*x1)**2 + V2*sp.sin(k2*x2)**2               # Define the sinusoidal lattice
    VoL_expanded = sp.series(VoL, x1, n=n1+1).removeO()           # Expand in x1 up to nth order
    VoL_expanded = sp.series(VoL_expanded, x2, n=n2+1).removeO()  # Expand in x2 up to nth order
    VoL_rmCM = VoL_expanded.subs(x2, X - x*mu1)                   # Substitute x2 by the CM-rm coordinates
    VoL_rmCM = VoL_rmCM.subs(x1, X + x*mu2)                       # Substitute x1 by the CM-rm coordinates
    VoL_rmCM_expanded = sp.expand(VoL_rmCM)                       # Expand the expression
    
    print('rm-CM coefficients')
    coeff = list()
    for i in range(0, nmax+1):
        for j in range(0, nmax+1):
            if i + j <= nmax:
                coeff.append(np.float64(VoL_rmCM_expanded.coeff(x, i).coeff(X, j).subs(dic_subs)))
                print(f'rm{i}CM{j}', coeff[-1])
    print('CM coefficients')
    coeff_CM = list()
    for i in range(0, nmax+1):
        coeff_CM.append(np.float64(VoL_rmCM_expanded.coeff(x, 0).coeff(X, i).subs(dic_subs)))
        print(coeff_CM[-1])
    
    print('config coefficients')
    coeff_config = {atom1: list(),
                   atom2: list()}
    #atom 1:
    for i in range(0, n1+1):
        coeff_config[atom1].append(np.float64(VoL_expanded.coeff(x2, 0).coeff(x1, i).subs(dic_subs)))
        print(coeff_config[atom1][-1])
    #atom 2:
    for i in range(0, n2+1):
        coeff_config[atom2].append(np.float64(VoL_expanded.coeff(x1, 0).coeff(x2, i).subs(dic_subs)))
        print(coeff_config[atom2][-1])
    
    return coeff, coeff_CM, coeff_config