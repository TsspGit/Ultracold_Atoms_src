__author__ = "Tomás Sánchez Sánchez-Pastor"
__date__   = "01/10/2021"
import sys
import numpy as np
from numpy import sqrt
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
	def __init__(self, nx, ny):
		self.nx = nx
		self.ny = ny

	def en(self, n, eta):
		return eta*(n + 1/2)

	def A3D(self, eta_x, eta_y, E, beta):
		out = -np.exp(beta*E) * np.sqrt(eta_x*eta_y/((4*pi)**3*np.sinh(eta_x*beta)*np.sinh(eta_y*beta)*np.sinh(beta)))\
        + 1/(4*pi*beta)**(3/2)
		return out

	def W3D(self, etax, etay, E):
		suma = 0
		for i in range(0, self.nx+1):
			for j in range(0, self.ny+1):
				if i%2 == 0 and j%2 == 0 and self.en(i, etax) + self.en(j, etay) + 1/2 <= E:
					suma += 2**(i + j - 1)*gamma(1/4 - (E - self.en(i, etax) - self.en(j, etay))/2)/ \
					(gamma((1-i)/2)**2 * gamma((1-j)/2)**2 * gamma(1+i) * gamma(1+j) * \
					gamma(3/4 - (E - self.en(i, etax) - self.en(j, etay))/2))
		return -pi/2 * sqrt(etax*etay/2) * suma

	def I3D(self, etax, etay, E, beta):
		suma = 0
		for i in range(0, self.nx+1, 2):
			for j in range(0, self.ny+1, 2):
				if self.en(i, etax) + self.en(j, etay) + 1/2 <= E:
					suma += 2**(i+j-1/2)*np.exp(beta*(E - etax - etay))/ \
					(gamma((1-i)/2)**2 * gamma((1-j)/2)**2 * gamma(1+i) * gamma(1+j))
		out = sqrt(pi*etax*etay/(8*np.sinh(beta))) * suma
		return self.A3D(etax, etay, E, beta) + out