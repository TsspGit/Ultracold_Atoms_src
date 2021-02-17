__author__ = "@Tssp"
import numpy as np
import os
import sympy as sp

def transpose_energies(Data):
    '''Inputs the Data list and outputs a dictionary with the evolution of the energy levels
       with respect the scattering length.
       
       
       Parameters
       ----------
       Data: list containing the N arrays from the N data folders.
       
       Returns
       -------
       dic: dictionary with the N size arrays of all the levels of the system, saved as nivel_0, nivel_1....
       '''
    
    L = Data[0].shape[0]
    l = len(Data)
    dic = {}
    for j in range(L):
        out = list()
        for i in range(l):
            out.append(Data[i][:,2][j])
        dic['nivel_{}'.format(j)] = out
    return dic




def dic_from_least_bound_forward(dic, wx):
    '''Inputs the dictionary of the transposed data, detects the least bound state (LBS) and the first trap state and filters
       the dictionary from the LBS forwards.
       
       
       Parameters
       ----------
       dic: The dictionary of the last function.
       
       
       Returns
       -------
       dic: Filtered dictionary
       i-1: Least bound state position.
       
    '''
    L = len(dic.keys())
    for i in range(L):
        if dic['nivel_{}'.format(i)][0] > 0:
            print("Least bound state: ", i-1, dic['nivel_{}'.format(i-1)][0]/wx)
            print("First trap state: ", i, dic['nivel_{}'.format(i)][0]/wx)
            break
    for j in range(1, i):
        del dic['nivel_{}'.format(j-1)]
    return dic, i-1




def cross_points(f, g):
    ''' This function takes the numpy 1-degree linear polyfit variables that contains the slope in the
        [0] position and the intercept in [1] of two different functions f and g and returns the cross
        point between them.
        
        
        Parameters
        ----------
        f: 1-degree linear polyfit variable.
        g: 1-degree linear polyfit variable.
        
        Returns
        -------
        Cross point between f and g.
    '''
    x = sp.symbols('x')
    return float(sp.solve(f[0]*x + f[1] - (g[0]*x + g[1]))[0])




def trap_fit(dic, x, x_lims, y_lims, levels, wy, tol=10):
    ''' This function takes the dictionary of the energy levels and find the trap line among the levels between
    level_inf and level_sup.
    
    Parameters
    ----------
    dic   : dictionary that contains the energy levels.
    x_lims     : tuple| (x_inf, x_sup)
    level : tuple| (level_inf, level_sup)
    tol   : float| trap state slope tolerance.
    
    Returns
    -------
    output  : fitted trap state
    '''
    m_ref = (np.diff(np.array(dic[f'nivel_{levels[0]}'])[(dic[f'nivel_{levels[0]}']/wy < y_lims[1]) & (dic[f'nivel_{levels[0]}']/wy > y_lims[0])])).min()
    print(m_ref)
    trap_slopes = []
    trap_values = []
    x_values    = []
    for l in range(levels[0], levels[1]+1):
        x_values    = np.concatenate((np.array(x)[(dic[f'nivel_{l}']/wy < y_lims[1]) & (dic[f'nivel_{l}']/wy > y_lims[0])], x_values))
        trap_values = np.concatenate((np.array(dic[f'nivel_{l}'])[(dic[f'nivel_{l}']/wy < y_lims[1]) & (dic[f'nivel_{l}']/wy > y_lims[0])], trap_values))
    trap_slopes = np.diff(trap_values)
    trap_values = np.array(trap_values[:-1])[(trap_slopes > -tol*m_ref) & (trap_slopes < tol*m_ref)]
    x_values    = np.array(x_values[:-1])[(trap_slopes > -tol*m_ref) & (trap_slopes < tol*m_ref)]
    output = np.polyfit(x_values, trap_values, deg=1)
    return output