#!/user/bin/env python
# -----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------
# PLOTTING PROGRAM OF THE LI7-LI7 SPECTRUM FOR THE STUDY OF INNSBRUCK EXPERIMENT
# -----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------
# Fabio Revuelta
# Grupo de Sistemas Complejos
# Universidad Politecnica de Madrid
# April 2019
# -----------------------------------------------------------------------------------------
# Preliminaries
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os, fnmatch
from matplotlib.ticker import  MultipleLocator, FormatStrFormatter
from decimal import Decimal
# -----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------

# Intensity (mW)
Ix = 4993.

# Mass (Da)
mass_Da = 7.0160040
mass = mass_Da * 1822.888457197206133193229
# Polarizability
alpha = 200
# Trap depth (a.u.)
Vx = alpha * Ix / ( 643640931E7 )
# Wave length (nm)
lambda_nm = 1000
# Wave number (a.u.)
a0 = 0.0529177249
kx = 2. * np.pi * a0 / lambda_nm

# Trapping frequency ( a.u. )
hbar = 1.0

wx = np.sqrt( 2.0 * Vx * kx * kx / mass)

w = wx
# -----------------------------------------------------------------------------------------
# Trap length in a.u. 
dx = np.sqrt( 2.0 * hbar / ( mass * w ) )
#d  = dx
#dx = np.sqrt( hbar / ( mass * w ) ) # It is actually a_h
d  = dx 
print('Ix (mW)       : ',  Ix)
print('Ix (a.u.)     : ',  Ix / ( 643640931E7 ) )
print('mass (Da)     : ',  mass_Da)
print('mass (a.u.)   : ',  mass)
print('alpha         : ',  alpha)
print('Vx (a.u.)     : ',  Vx)
print('lambda (nm)   : ',  lambda_nm)
print('kx (a.u.)     : ',  kx)
print('wx (a.u.)     : ',  wx)
print('ah (a.u.)     : ',  dx)
# -----------------------------------------------------------------------------------------
# Scattering length
# Delta    asc (a.u.)




delta = []
asc   = []

delta.append( Decimal('0.7080')  )
asc.append( 1191.58168082959 )

delta.append( 0.70805 )
asc.append(  1199.73395812338 )

delta.append( 0.7081   )
asc.append(1208.00150021930 )

delta.append( 0.70815  )
asc.append( 1216.38675672344 )

delta.append( 0.7082  )
asc.append( 1224.89227055587 )

delta.append( 0.70825 )
asc.append(  1233.52065585274 )

delta.append( 0.7083   )
asc.append(1242.27458068079 )

delta.append( 0.70835 )
asc.append(  1251.15681153772 )

delta.append( 0.7084  )
asc.append( 1260.17019870601 )

delta.append( 0.70845  )
asc.append( 1269.31764930816 )

delta.append( 0.7085   )
asc.append(1278.60219549329 )

delta.append( 0.70855  )
asc.append( 1288.02693854581 )

delta.append( 0.7086  )
asc.append( 1297.59506832038 )

delta.append( 0.70865  )
asc.append( 1307.30989759343 )

delta.append( 0.7087   )
asc.append(1317.17482778704 )

delta.append( 0.70875  )
asc.append( 1327.19334316745 )

delta.append( 0.7088  )
asc.append( 1337.36908859614 )

delta.append( 0.70885 )
asc.append(  1347.70577729459 )

delta.append( 0.7089  )
asc.append( 1358.20726171416 )

delta.append( 0.70895 )
asc.append(  1368.87751842734 )

delta.append( Decimal('0.7090')  )
asc.append( 1379.72064594415 )

delta.append( 0.70905 )
asc.append(  1390.74087793971 )

delta.append( 0.7091   )
asc.append(1401.94260401253 )

delta.append( 0.70915 )
asc.append(  1413.33032919673 )

delta.append( 0.7092  )
asc.append( 1424.90873625383 )

delta.append( 0.70925 )
asc.append(  1436.68264328870 )

delta.append( 0.7093 )
asc.append(  1448.65706791664 )

delta.append( 0.70935  )
asc.append( 1460.83715304107 )

delta.append( 0.7094 )
asc.append(  1473.22824971113 )

delta.append( 0.70945 )
asc.append(  1485.83590393415 )

delta.append( 0.7095 )
asc.append(  1498.66583449300 )

delta.append( 0.70955  )
asc.append( 1511.72396478758 )

delta.append( 0.7096 )
asc.append(  1525.01645260571 )

delta.append( 0.70965 )
asc.append(  1538.54965309584 )

delta.append( 0.7097  )
asc.append( 1552.33017349492 )

delta.append( 0.70975 )
asc.append(  1566.36484680843 )

delta.append( 0.7098  )
asc.append( 1580.66078690389 )

delta.append( 0.70985 )
asc.append(  1595.22536089131 )

delta.append( 0.7099 )
asc.append(  1610.06618747999 )

delta.append( 0.70995 )
asc.append(  1625.19121772037 )

delta.append( Decimal('0.7100') )
asc.append(  1640.60870394915 )

#delta.append( Decimal('0.7105'))  
#asc.append( 1812.93415170932)

delta.append( 0.7101  )
asc.append( 1672.35561079401 )

delta.append( 0.7102 )
asc.append(  1705.37959037706 )

#delta.append( Decimal('0.7115'))  
#asc.append( 2298.38803849248 )
#
#
L = len(delta)
print
print( 'L : ', L )
#-----------------------------------------------------------------------------------------
# Plotting parameters
c = 1.75
cm_inc = 32.0/80.0
w_fig = 7 * cm_inc * c
h_fig = w_fig / 1.6
fig, ax = plt.subplots(figsize = (w_fig, h_fig ) )
#fig, ax = plt.subplots(figsize = (1, 1 ) )
#ax.set_aspect(0.621)

# Plotting range
ascmin = np.min( np.abs(asc) ) 
ascmax = np.max( np.abs(asc) ) 

xmin = 1.05
xmax = 1.46
emin =-0.3 #1.5 #-0.3
emax = 3.0 #3.6

print( 'min(asc) : ', ascmin)
print( 'max(asc) : ', ascmax)
print
print( 'xmin : ', xmin)
print( 'xmax : ', xmax)
print
print( 'emin : ', emin)
print( 'emax : ', emax)
print
#xmin = xmin - 1
#xmax = 5 #40

# Plotting limits
ax.set_xlim([ xmin, xmax] )
ax.set_ylim([ emin, emax] )

# Point size
scattersize = 2

# Labels
labelsize = 15

#ax.set_xlabel('$a_{sc}$\,(a.u.)')
ax.set_xlabel('$a_h/a_{sc}$', fontsize=labelsize)
ax.set_ylabel('$E/(\hbar\,\omega_x)$', fontsize=labelsize)
#-----------------------------------------------------------------------------------------
# Eigenenergies for the Ag symmetry of the D2h symmetry group
#
file1 = '../eva/ix4993_iy4993_iz50/Li7Li7_b_x20000_y20000_z20000_152rm8g2l50m10_102CM8g1L50M10_Li7a200_Li7a200_kx1000_ky1000_kz1000_ix4993_iy4993_iz50_LiLi_a3Sup_'
file2 = '_sinTnx6_sinTny6_sinTnz6/'
#
# Loop
x_resonance   = []
y_resonance   = []
asc_resonance = []
k=-1
for i in range(0, L):
  print(i, delta[i], asc[i], d/asc[i])  
  
# Input file name
  fileref = file1 + str(delta[i]) + file2
  print(fileref)
  config = os.listdir(fileref)
  
  for entry in config:
    if fnmatch.fnmatch(entry, '*b.eva'):
#      configname = entry
      configname = 'Ag_vsLiLi_500-800_75b.eva'
    
  fileref = fileref + configname
# print 'fileref : ', fileref
  
# Read input file
  data = np.loadtxt(fileref)
# x = data[:, 0]
# y = data[:, 2]
  y = data[:, 2]/w
  x = y * 0. + d/asc[i]
  
  Ly = len(y)
  j = -1
  while j <= Ly:
    j = j + 1
    if y[j] < 0 :
      jref = j
    else:
      break
      
  j = -1
  k = -1
  while j < Ly-1:
    j = j + 1
    if x[j] >=xmin and x[j] <= xmax and y[j] >= emin and y[j] <= emax :
      k=k+1
      x_resonance.append(x[j])
      y_resonance.append(y[j])
      asc_resonance.append(asc[i])

  print( 'Least bound state :', j, ' E : ', y[jref]*w, ' E/w : ', y[jref])
  print( 'First trap  state :', j+1, ' E : ', y[jref+1]*w, ' E/w : ', y[jref+1])
  print
  
  plt.scatter(x, y, s = scattersize, c ='k')
  


print( )
print( )
print( )
print( )
for i in range(0, L):
  print(i, delta[i], asc[i], d/asc[i])  



# Ticks
ticksize = (2/3)*labelsize
ax.xaxis.set_tick_params( labelsize = ticksize)
ax.yaxis.set_tick_params( labelsize = ticksize)
#tick_params(axis='both',which='major',length=8,width=1.2,direction='in')
#tick_params(axis='both',which='minor',length=5,width=1,direction='in')
tick_params(direction='in')

plt.tight_layout()

#plt.plot(xexp, yexp, ls = '--', c = 'k')
#plt.plot(xref, yref, ls = '-', c = 'g')
fig.savefig('spectrum_ci.pdf')
plt.show()
