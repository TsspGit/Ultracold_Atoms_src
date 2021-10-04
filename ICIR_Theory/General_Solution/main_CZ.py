__author__ = "Tomás Sánchez Sánchez-Pastor"
__date__   = "01/10/2021"
import numpy as np
import matplotlib.pyplot as plt
from math import pi, sqrt
from Integration import Integral
from Solution_Class import J3D_funs
from Trap_Class import Trap
plt.rc('text',usetex=True)
plt.rc('font',family='serif')
ref_ticksize = 16
plt.rcParams['xtick.labelsize']=ref_ticksize
plt.rcParams['legend.fontsize']=ref_ticksize
plt.rcParams['ytick.labelsize']=ref_ticksize
plt.rcParams['axes.labelsize']=ref_ticksize * 3/2
plt.rcParams['axes.titlesize']=ref_ticksize * 3/2

def separate_levels(A, E):
    '''
    Parameters
    ----------
    - A: x-axis (usually a3D)
    - E: y-axis (usually the energy)

    Returns
    -------
    - dic: dictionary that contains each level separately denoted by the key "a3D_n{level}" and "E_n{level}"
    - level: last level separated
    '''
    dic = {}
    count = 0
    level = 0
    for i in range(len(E)-1):
        if np.sign(A[i]) != np.sign(A[i+1]) and A[i] > abs(2):
            level+=1
            dic[f'a3D_n{level}'] = A[i-count:i+1]
            dic[f'E_n{level}']   = E[i-count:i+1]
            count = 0
        elif np.sign(A[i]) == np.sign(A[i+1]):
            count += 1
    dic[f'a3D_n{level+1}'] = A[i-count+1:i+2]
    dic[f'E_n{level+1}']   = E[i-count+1:i+2]
    return dic, level

if __name__ == '__main__':
	Trap = Trap()
	E  = np.linspace(-7.5, 13, num=1000)
	J3D_init = J3D_funs(nx=100, ny=100)

	print(f'''
	         Parameters
	      ----------------
	      etax:    {Trap.eta_x}
	      etay:    {Trap.eta_y}
	      nx:      {J3D_init.nx}    
	      ny:      {J3D_init.ny} 
	      Eo:      {Trap.Eo}
	      E:       {E[0]}-{E[-1]}
	''')
	
	# Integration
	a3D = []
	for e in E:
		# Initialize J3D
		I3D = Integral().int(lambda x: J3D_init.I3D(Trap.eta_x, Trap.eta_y, e, x), 100, e, Trap.Eo)
		J3D = sqrt(2)*4*pi*(J3D_init.W3D(Trap.eta_x, Trap.eta_y, e) + I3D)
		a3D.append(1/J3D)
	print(f'a3D: [{min(a3D)} - {max(a3D)}]')
	# Figure
	Spectrum, level = separate_levels(a3D, E)
	fig, ax = plt.subplots(figsize=(8,6))
	for i in range(1, level+1):
	    plt.plot(Spectrum[f'a3D_n{i}'], Spectrum[f'E_n{i}'], 'C0')
	ax.set_xlabel(r'$a_{3D}/d_y$')
	ax.set_ylabel(r'$E/(\hbar \omega_z)$')
	ax.set_xlim(-10, 10)
	ax.set_ylim(-7.5, 12)
	plt.grid()
	#plt.savefig('E_as_isotropic.png', dpi=200)
	plt.show()