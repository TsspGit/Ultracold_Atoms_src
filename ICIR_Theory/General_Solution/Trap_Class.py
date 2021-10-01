__author__ = "Tomás Sánchez Sánchez-Pastor"
__date__   = "01/10/2021"

class Trap:
	''' Modelating and giving the parameters of the optical trap
	'''
	def __init__(self):
		''' eta_j = wj/wz
		'''
		self.eta_x = 1
		self.eta_y = 1
		self.Eo    = 1/2*(self.eta_x + self.eta_y + 1)
		if self.eta_x == self.eta_y:
			print('Isotropic trap')
		else:
			print(f'Anisotropic trap\neta_x: {self.eta_x}\neta_y: {self.eta_y}')