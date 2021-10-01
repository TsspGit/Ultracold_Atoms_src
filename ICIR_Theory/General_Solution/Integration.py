__author__ = "Tomás Sánchez Sánchez-Pastor"
__date__   = "01/10/2021"
import numpy as np
from math import pi, sqrt, isnan

class Integral:
	'''Computes the integral of J3D through a Gauss-Legendre quadrature.
	Each case requires a unique configuration'''

	def __init__(self):
		pass

	def gauss(self, f, n, a, b):
		''' Gauss-Legendre quadrature of f in [a, b] with order n.
		'''
		from scipy.special.orthogonal import p_roots
		[x,w] = p_roots(n+1)
		G=0.5*(b-a)*sum(w*f(0.5*(b-a)*x+0.5*(b+a)))
		return G

	def int(self, f, n, E, Eo):
		if E>Eo:
			Lambda = 0.15/(E - Eo)
			integral = self.gauss(f, n, 1e-6, Lambda)
		else:
			Lambda = 100
			integral = self.gauss(f, n, 1e-6, Lambda)
		return integral