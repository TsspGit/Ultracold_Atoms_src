__author__ = "@Tssp"
import numpy as np
import os

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

def dic_from_least_bound_forward(dic):
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
    
    for i in range(L):
        if dic['nivel_{}'.format(i)][0] > 0:
            print("Least bound state: ", i-1, dic['nivel_{}'.format(i-1)][0]/wx)
            print("First trap state: ", i, dic['nivel_{}'.format(i)][0]/wx)
            break
    for j in range(1, i):
        del dic['nivel_{}'.format(j-1)]
    return dic, i-1