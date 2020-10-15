#!/usr/local/opt/python@3.8/bin/python3.8
__author__ = '@Tssp'
import numpy as np
from atomic_units import ao, vo, e, hbar, me, Eh, to
import sys

class Units_Transformer:
    def __init__(self, lambd, coefER, alpha, m):
        '''lambd: wavelength in nm
           coefER: V = coef*ER
           alpha: polarizability in u.a
           I1: intensity in one direction in W/cm2
           m: mass'''
        self.lambd = lambd * 1e-9 / ao
        self.k = 2*np.pi/self.lambd
        self.coefER = coefER
        self.alpha = alpha
        self.m = m * 1.66053873e-27 / me
        print(f"\nParameters in Atomic units:\n-mass: {self.m}\n-polarizability: {self.alpha}\n-Wavelength: {self.lambd}\n-k: {self.k}")
        
    def Vj_to_Ij(self):
        eta = float(input("\nRelation between omega2 and omega1: "))
        V1 = self.coefER * self.k**2/(2 * self.m)
        I1 = V1/self.alpha
        I2 = eta**2*I1
        omega2 = self.k * np.sqrt(2*self.alpha*I2/self.m)
        dho = np.sqrt(2/(self.m*omega2))
        print(f"\nIntensities in W/cm2:\n-I1: {I1/(1e4 / Eh * to * ao**2)}\n-I2: {I2/(1e4 / Eh * to * ao**2)}")
        print(f"Trap length in direction 2: ", dho, "u.a")
        
        
if __name__ == '__main__':
    lambd = float(input("Wavelength in nm: "))
    coefER = float(input("Potential in terms of Er: "))
    alpha = float(input("Polarizability in u.a: "))
    m = float(input("mass of the atoms in Da: "))
    UT = Units_Transformer(lambd, coefER, alpha, m)
    UT.Vj_to_Ij()