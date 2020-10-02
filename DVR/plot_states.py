#!/usr/local/opt/python@3.8/bin/python3.8
__author__ = '@Tssp'
__date__   = '30 Sept 2020'
from numpy import loadtxt
import matplotlib.pyplot as plt
plt.rc('text',usetex=True)
plt.rc('font',family='serif')
plt.rcParams['xtick.labelsize']=16
plt.rcParams['ytick.labelsize']=16
plt.rcParams['axes.labelsize']=18
plt.rcParams['axes.titlesize']=18

n = input("Input the energy level: ")
psi = loadtxt('eig_m1_' + '0'*(4 - len(n)) + f'{n}.plt')
print('eig_m1_' + '0'*(4 - len(n)) + f'{n}.plt')
fig, ax = plt.subplots(figsize=(8, 6))
ax.plot(psi[:,0], -psi[:,1])
ax.set_xlabel('x')
ax.set_ylabel('$\psi$')
plt.show()