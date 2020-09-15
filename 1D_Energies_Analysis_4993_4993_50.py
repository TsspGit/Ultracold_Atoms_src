#!/usr/local/opt/python@3.8/bin/python3.8
__author__ = "@Tssp"
import numpy as np
import matplotlib.pyplot as plt
import os
from utils.atomic_units import ao, vo, e, hbar, me, Eh, to
from decimal import Decimal
from utils.Energies_Analysis_utils import transpose_energies, dic_from_least_bound_forward, cross_points
plt.rc('text',usetex=True)
plt.rc('font',family='serif')
plt.rcParams['xtick.labelsize']=16
plt.rcParams['ytick.labelsize']=16
plt.rcParams['axes.labelsize']=18
plt.rcParams['axes.titlesize']=18


# Parameters and Constants:
###########################
delta, asc = np.loadtxt('Data/citaold2h/python/delta_asc.txt')
delta = delta.tolist()
delta[40] = Decimal('0.7100')
asc = asc.tolist()
mass = 7.0160040 * 1.66053873e-27 / me # a.u
alpha = 200 # a.u
lambd = 1000 * 1e-9 / ao # a.u
kx = 2*np.pi/lambd
Ix = 4993 * (1e4 / Eh * to * ao**2)
Vx = alpha * Ix # a.u
Ix = Vx/alpha
wx = np.sqrt( 2 * Vx * kx**2 / mass)       
dho = np.sqrt(2 / (mass * wx))
x = dho / asc
print(f'''Parameters:
          wl  [nm]   : {lambd*ao*1e9}
          mass[a.u]  : {mass}
          Ix  [W/cm2]: {Ix / ((1e4 / Eh * to * ao**2))}
          wx  [a.u]  : {wx/to}
          dho [a.u]  : {dho}''')

# Data mining:
##############
print(f'Current directory:\n {os.getcwd()}')
folder_path = "Data/citaold2h/eva/ix4993_iy4993_iz50/"
print(f'folder path:\n {folder_path}')

Data = list()
print(f'Reading files, please wait...')
for d in delta:
    try:
        Data.append(np.loadtxt(folder_path + f'Li7Li7_b_x20000_y20000_z20000_152rm8g2l50m10_102CM8g1L50M10_Li7a200_Li7a200_kx1000_ky1000_kz1000_ix4993_iy4993_iz50_LiLi_a3Sup_{d}_sinTnx6_sinTny6_sinTnz6/Ag_vsLiLi_500-800_75b.eva'))
    except:
        Data.append(np.loadtxt(folder_path + f'Li7Li7_b_x20000_y20000_z20000_152rm8g2l50m10_102CM8g1L50M10_Li7a200_Li7a200_kx1000_ky1000_kz1000_ix4993_iy4993_iz50_LiLi_a3Sup_{d}0_sinTnx6_sinTny6_sinTnz6/Ag_vsLiLi_500-800_75b.eva'))

print(f'{len(Data)} files has been readed, each one with {Data[0].shape[0]} energy levels')
dic = transpose_energies(Data)
dic, least_bound_pos = dic_from_least_bound_forward(dic, wx)

# Plots/Figures:
################
fig, ax = plt.subplots(figsize=(8,6))
for i in range(least_bound_pos, least_bound_pos + 30):
    ax.plot(x, dic['nivel_{}'.format(i)]/wx, 'b')
ax.set_ylim(0, 3)
ax.set_xlabel('$d_{\perp}/a$')
ax.set_ylabel('$E/\hbar \omega_{\perp}$')
plt.tight_layout()
fig.savefig('General_figures/Ix4993_Iy4993_Iz50/Ix4993_Iy4993_Iz50_Easc.png', dpi=200)

fig3, ax3 = plt.subplots(figsize=(8,5))
# nivel 1056:
f1056 = np.polyfit(x[x < 1.2], np.array(dic['nivel_1056'])[x < 1.2]/wx, deg=1)
ax3.plot(x, np.polyval(f1056, x), 'ro', markersize=4, label='Diabetic')
ax3.plot(x, dic['nivel_1056']/wx, 'b', label='Adiabatic')

# nivel 1057:
f1057 = np.polyfit(x[x < 1.2], np.array(dic['nivel_1057'])[x < 1.2]/wx, deg=1)
x1057 = np.linspace(1, 1.3, num=500)
ax3.plot(x1057, np.polyval(f1057, x1057), 'ro', markersize=4)
ax3.plot(x, dic['nivel_1057']/wx, 'b')

# nivel 1058:
f1058 = np.polyfit(x[(x > 1.2) & (x < 1.3)], np.array(dic['nivel_1058'])[(x > 1.2) & (x < 1.3)]/wx, deg=1)
x1058 = np.linspace(1.2, 1.3, num=500)
ax3.plot(x1058, np.polyval(f1058, x1058), 'ro', markersize=4)
ax3.plot(x, dic['nivel_1058']/wx, 'b')

# nivel 1059:
# same as 1058
ax3.plot(x, dic['nivel_1059']/wx, 'b')

# nivel 1060:
f1060 = np.polyfit(x[(x > 1.3)], np.array(dic['nivel_1060'])[(x > 1.3)]/wx, deg=1)
x1060 = np.linspace(1.3, 1.4, num=500)
ax3.plot(x1060, np.polyval(f1060, x1060), 'ro', markersize=4)
ax3.plot(x, dic['nivel_1060']/wx, 'b')

# nivel 1061:
# same as 1060:
ax3.plot(x, dic['nivel_1061']/wx, 'b')

# Crosses:
##########
ax3.plot(cross_points(f1056, f1057), np.polyval(f1056, cross_points(f1056, f1057)), 'k*', markersize=10, label='Cross')
ax3.plot(cross_points(f1056, f1058), np.polyval(f1056, cross_points(f1056, f1058)), 'k*', markersize=10)
ax3.plot(cross_points(f1056, f1060), np.polyval(f1056, cross_points(f1056, f1060)), 'k*', markersize=10)


ax3.set_xlabel('$d_{\perp}/a$')
ax3.set_ylabel('$E/\hbar \omega_{\perp}$')
ax3.legend(fontsize=14, bbox_to_anchor=(1.03, 1))
ax3.set_ylim(2.07, 2.076)
ax3.set_xlim(1.15, 1.45)
plt.tight_layout()
fig3.savefig('General_figures/Ix4993_Iy4993_Iz50/Ix4993_Iy4993_Iz50_Easc_Interpolation.png', dpi=200)

print(f'''Crosses between (x, E):
* 1056-1057: ({cross_points(f1056, f1057)}, {np.polyval(f1056, cross_points(f1056, f1057))})
* 1056-1058: ({cross_points(f1056, f1058)}, {np.polyval(f1056, cross_points(f1056, f1058))})
* 1056-1060: ({cross_points(f1056, f1060)}, {np.polyval(f1056, cross_points(f1056, f1060))})''')
plt.show()