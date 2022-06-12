__author__ = 'Tomás Sánchez Sánchez-Pastor'
__date__ = '30 de Mayo de 2022'

import numpy as np
import matplotlib.pyplot as plt
from input_plot_potential import file
plt.rc('text',usetex=True)
plt.rc('font',family='serif')
ref_ticksize = 16
plt.rcParams['xtick.labelsize']=ref_ticksize
plt.rcParams['legend.fontsize']=ref_ticksize
plt.rcParams['ytick.labelsize']=ref_ticksize
plt.rcParams['axes.labelsize']=ref_ticksize * 3/2
plt.rcParams['axes.titlesize']=ref_ticksize * 3/2

def main():
    V = np.genfromtxt(file, skip_header=14, skip_footer=22)
    r = V[:,0]
    V = V[:,1]
    minpos = np.min(np.where(V == min(V))[0])
    print(min(V))
    plt.plot(r[minpos-15:minpos+50], V[minpos-15:minpos+50])
    plt.grid()
    plt.xlabel(r'$r(u.a)$', x=0.9)
    plt.ylabel(r'$V(r)(u.a)$', y=.86)
    plt.title(file[:file.find('_')])
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    main()