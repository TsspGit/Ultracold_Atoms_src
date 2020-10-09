#!/usr/local/opt/python@3.8/bin/python3.8
__author__ = "@Tssp"
__date__   = "9/10/20"
import numpy as np
import matplotlib.pyplot as plt
import os
plt.rc('text',usetex=True)
plt.rc('font',family='serif')
plt.rcParams['xtick.labelsize']=16
plt.rcParams['ytick.labelsize']=16
plt.rcParams['axes.labelsize']=18
plt.rcParams['axes.titlesize']=18

for f in os.listdir():
    if 'LiLi' in f:
        print(f)

file = input("\n\nWhich of the previous orbitals do you want to plot? ")
D1007 = np.loadtxt(file)
x = D1007[:,0]
y = D1007[:,1]
wf = D1007[:,2]

fig, ax = plt.subplots(figsize=(8,6))
cm = ax.tricontourf(x, y, wf, cmap='jet')
fig.colorbar(cm)
ax.set_xlabel('$x_1 (a.u)$')
ax.set_ylabel('$x_2 (a.u)$')
plt.tight_layout()
plt.show()