__author__ = "@Tssp"
__date__   = "15/04/21"
import sympy as sp
from sympy.functions.special.polynomials import chebyshevt
from sympy.abc import n, x
import numpy as np

def Chebyshev_Expansion(f, order):
    '''
    Parameters
    ----------
    f: function to expand in terms of the Chebyshev polynomials f(x)
    m: order of expansion
    
    Returns
    -------
    C: Expanded function
    '''
    C = 0
    for m in range(order):
        if m == 0: am = 2
        elif m>0: am=1
        cm = 2/(np.pi*am) * sp.integrate(f*chebyshevt(m, x)/sp.sqrt(1 - x**2), (x,-1,1)).evalf()
        C += cm * chebyshevt(m, x)
    print(C)
    return C

def np_Chebyshev(x, m):
    '''
    Parameters
    ----------
    x: [array] X axis.
    order: [int] order of the chebyshev polynomials.
    
    Returns
    -------
    Tn: [array] Chebyshev polynomial of order m.
    '''
    if m == 0:
        T0 = np.ones(x.size)
        return T0
    elif m == 1:
        T1 = x
        return T1
    else:
        T0 = np.ones(x.size)
        T1 = x
        T_list = [T0, T1]
        for m in range(1, m):
            T_list.append(2*x*T_list[-1] - T_list[-2])
        return T_list[-1]