#!/usr/local/opt/python@3.8/bin/python3.8

__author__ = "@Tssp"
__date__ = "09/04/21"
import numpy as np
import matplotlib.pyplot as plt
from math import pi
from input_DVR_pot import hbar, m1, m2, pot, wLx1, wLy1, wLz1, wLx2, wLy2, wLz2, delta, xmax, xmin, mode, coeff
from method_pot import DVR_method
plt.rc('text',usetex=True)
plt.rc('font',family='serif')
plt.rcParams['xtick.labelsize']=16
plt.rcParams['ytick.labelsize']=16
plt.rcParams['axes.labelsize']=18
plt.rcParams['axes.titlesize']=18

# Parameters:
atoms = np.loadtxt('../Expansion/out/' + coeff, delimiter='\n', max_rows=2, dtype=str)
atom1 = atoms[0]
atom2 = atoms[1]
coeff_data = np.loadtxt('../Expansion/out/' + coeff, delimiter='\n', skiprows=2)
alpha1= coeff_data[0]
alpha2= coeff_data[1]
nx1   = int(coeff_data[3])
ny1   = int(coeff_data[4])
nz1   = int(coeff_data[5])
nx2   = int(coeff_data[6])
ny2   = int(coeff_data[7])
nz2   = int(coeff_data[8])
nxmax = np.max([nx1, nx2])
nymax = np.max([ny1, ny2])
nzmax = np.max([nz1, nz2])
if mode == 'CM':
    coeff_x = list()
    for i in range(nxmax+1):
        coeff_x.append(coeff_data[i + 9])
    coeff_y = list()
    for i in range(nymax+1):
        coeff_y.append(coeff_data[i + 10 + nxmax])
    coeff_z = list()
    for i in range(nzmax+1):
        coeff_z.append(coeff_data[i + 11 + nxmax + nymax])
elif mode == 'all':
    # Atom 1:
    coeff_x1 = list()
    for i in range(nx1+1):
        coeff_x1.append(coeff_data[i + 9])
    coeff_x2 = list()
    for i in range(nx2+1):
        coeff_x2.append(coeff_data[i + 10 + nx1])
    coeff_y1 = list()
    for i in range(ny1+1):
        coeff_y1.append(coeff_data[i + 11 + nx1 + nx2])
    # Atom2
    coeff_y2 = list()
    for i in range(ny2+1):
        coeff_y2.append(coeff_data[i + 12 + nx1 + nx2 + ny1])
    coeff_z1 = list()
    for i in range(nz1+1):
        coeff_z1.append(coeff_data[i + 13 + nx1 + nx2 + ny1 + ny2])
    coeff_z2 = list()
    for i in range(nz2+1):
        coeff_z2.append(coeff_data[i + 14 + nx1 + ny1 + nz1 + nx2 + ny2])
    coeff_x = [coeff_x1, coeff_x2]
    coeff_y = [coeff_y1, coeff_y2]
    coeff_z = [coeff_z1, coeff_z2]
    print([coeff_x1, coeff_x2])
    print([coeff_y1, coeff_y2])
    print([coeff_z1, coeff_z2])

# Intensities for both particles in every direction readed from txt filename:
#############################################################################
Ix1 = np.float(coeff[coeff.find('ix')+2:coeff.find('iy')])/6.436409310e15 # a.u
Iy1 = np.float(coeff[coeff.find('iy')+2:coeff.find('iz')])/6.436409310e15 # a.u
Iz1 = np.float(coeff[coeff.find('iz')+2:coeff.rfind('_ix')])/6.436409310e15 # a.u
Ix2 = np.float(coeff[coeff.rfind('ix')+2:coeff.rfind('iy')])/6.436409310e15 # a.u
Iy2 = np.float(coeff[coeff.rfind('iy')+2:coeff.rfind('iz')])/6.436409310e15 # a.u
Iz2 = np.float(coeff[coeff.rfind('iz')+2:coeff.rfind('_')])/6.436409310e15 # a.u

# Potential depths for both particles in every direction:
####################################################
Vx1 = alpha1 * Ix1
Vx2 = alpha2 * Ix2
Vy1 = alpha1 * Iy1
Vy2 = alpha2 * Iy2
Vz1 = alpha1 * Iz1
Vz2 = alpha2 * Iz2

# Wavelengths for both particles in every direction:
####################################################
wLx1 = wLx1 / 0.0529177249
wLy1 = wLy1 / 0.0529177249
wLz1 = wLz1 / 0.0529177249
wLx2 = wLx2 / 0.0529177249
wLy2 = wLy2 / 0.0529177249
wLz2 = wLz2 / 0.0529177249

# Wavenumber for both particles in every direction:
####################################################
kx1  = 2*pi/wLx1
ky1  = 2*pi/wLy1
kz1  = 2*pi/wLz1
kx2  = 2*pi/wLx2
ky2  = 2*pi/wLy2
kz2  = 2*pi/wLz2

# Frequencies for both particles in every direction:
####################################################
if mode == 'CM':
    m = m1 + m2
    wx = np.sqrt(Vx1/m*kx1**2 + Vx2/m*kx2**2)*np.sqrt(2) # Alomejor *sqrt(2)
elif mode == 'all':
    wx = np.sqrt(Vx1/m1*kx1**2 + Vx2/m2*kx2**2)

print(f""" 
           Atom 1                  =    {atom1}
           X-axis order (Taylor) 1 =    {nx1}
           Y-axis order (Taylor) 1 =    {ny1}
           Z-axis order (Taylor) 1 =    {nz1}
           Ix(mW/cm2) 1            =    {Ix1 * 6.436409310e15}
           Iy(mW/cm2) 1            =    {Iy1 * 6.436409310e15}
           Iz(mW/cm2) 1            =    {Iz1 * 6.436409310e15}
           wLx(nm) 1               =    {wLx1 * 0.0529177249}
           wLy(nm) 1               =    {wLy1 * 0.0529177249}
           wLz(nm) 1               =    {wLz1 * 0.0529177249}
           alpha(a.u) 1            =    {alpha1}
           kx(a.u) 1               =    {kx1}
           ky(a.u) 1               =    {ky1}
           kz(a.u) 1               =    {kz1}
           Vx(a.u) 1               =    {Vx1}  
           Vy(a.u) 1               =    {Vy1}   
           Vz(a.u) 1               =    {Vz1}
           
           Atom 2                  =    {atom2}
           X-axis order (Taylor) 2 =    {nx2}
           Y-axis order (Taylor) 2 =    {ny2}
           Z-axis order (Taylor) 2 =    {nz2}
           Ix(mW/cm2) 2            =    {Ix2 * 6.436409310e15}
           Iy(mW/cm2) 2            =    {Iy2 * 6.436409310e15}
           Iz(mW/cm2) 2            =    {Iz2 * 6.436409310e15}
           wLx(nm) 2               =    {wLx2 * 0.0529177249}
           wLy(nm) 2               =    {wLy2 * 0.0529177249}
           wLz(nm) 2               =    {wLz2 * 0.0529177249}
           alpha(a.u) 2            =    {alpha2}
           kx(a.u) 2               =    {kx2}
           ky(a.u) 2               =    {ky2}
           kz(a.u) 2               =    {kz2}
           Vx(a.u) 2               =    {Vx2}
           Vy(a.u) 2               =    {Vy2}
           Vz(a.u) 2               =    {Vz2}
           wx(a.u)                 =    {wx}

           General:
           delta                   =    {delta}
           xmin                    =    {xmin}
           xmax                    =    {xmax}
           coeff file              =    {coeff}
      """)
# xaxis:
########
x = np.linspace(xmin, xmax, int((xmax-xmin)/delta))
N = len(x)
print("\nX Axis: \n")
Ex, _ = DVR_method(N, delta, m1, m2, kx1, kx2, x, coeff_x, wx, mode)
print("\nY Axis: \n")
Ey, _ = DVR_method(N, delta, m1, m2, ky1, ky2, x, coeff_y, wx, mode)
print("\nZ Axis: \n")
Ez, _ = DVR_method(N, delta, m1, m2, kz1, kz2, x, coeff_z, wx, mode)
print("\nEx + Ey + Ez for two atoms: \n")
Et = Ex + Ey + Ez
print("        (nx,ny,nz)           E                 E[hbar wx]")
for i in range(5):
    if i%2 == 0:
        for j in range(5):
            if j%2 == 0:
                for k in range(12):
                    if k%2 == 0:
                        if mode == 'all':
                            print(f'          ({i},{j},{k}) {(Ex[i] + Ey[j] + Ez[k])}   {(Ex[i] + Ey[j] + Ez[k])/wx}')
                        elif mode == 'CM':
                            print(f'          ({i},{j},{k}) {(Ex[i] + Ey[j] + Ez[k])}   {(Ex[i] + Ey[j] + Ez[k])/wx}')
print("\nBingo !")