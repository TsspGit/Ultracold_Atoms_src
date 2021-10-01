__author__ = "Tomás Sánchez Sánchez-Pastor"
__date__   = "01/10/2021"
import sys
import numpy as np
from math import pi
from scipy.special import gamma, hyp2f1
from input_Chen_Zhang_2020 import eta_x, eta_y, nx, ny
eps = sys.float_info.epsilon


class J3D_funs:
	''' This class contains all the methods that finally
	computes the dependence between the rel-Energy and the
	s-wave scattering length a0 for a system of two cold
	atoms confined in anisotropic optical traps (wx!=wy!=wz)
	J3D = 4π√2(W3D(E, ß) + ∫dß(I3D(E, ß), 0, +oo))
	a3D = 1/J3D

	Ref. [Physical Review A, 101(5), 053624] 

	Parameters
	----------
	eta_j: wj/wy
    E:     Energy (self)
    beta:  X-axis (self)
	n_j:   Max levels of the summation

	Methods
	-------
	· en(n, eta): HO energy
	· A3D(etax, etay, E, beta): integrand function of J3D
	· W3D(nx, ny, etax, etay, E): part of J3D
	· I3D(etax, etay, nx, ny, E, beta): int(A3D, 0, +oo) + sumation
	'''
	def __init__(self):
		pass

	def A3D(self, eta_x, eta_y):
		out = -np.exp(self.beta*self.E) * np.sqrt(eta_x*eta_y/((4*pi)**3*np.sinh(eta_x*self.beta)*np.sinh(eta_y*self.beta)*np.sinh(self.beta)))\
        + 1/(4*pi*self.beta)**(3/2)
		return out

	def W3D(self, nx, ny, etax, etay, E):
		pass

	def I3D(self, nx, ny, etax, etay, E):
		pass


