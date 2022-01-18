__author__ = "@Tssp"
__date__ = "08/04/21"

coeff  = 'Li7Li7_nx6ny6nz6_nx6ny6nz6_ix7802iy4993iz4993_ix7802iy50iz50_CM.coeff'
pot    = 'sin2'               # sin2 or cos2 expansion
mode   = 'CM'                 # Compute mode: all or CM
delta  = 50                   # Grid spacing
xmax   = 18000                # xmax
xmin   = -18000               # xmin

# Atom 1:
#########
wLx1   = 1000                 # Wavelength in nm
wLy1   = 1000                 # Wavelength in nm
wLz1   = 1000                 # Wavelength in nm
hbar   = 1                    # Planck reduced Constant
m1     = 12789.3927072494     # Li mass in a.u

# Atom 2:
#########
pot    = 'sin2'               # sin2 or cos2 expansion
wLx2   = 1000                 # Wavelength in nm
wLy2   = 1000                 # Wavelength in nm
wLz2   = 1000                 # Wavelength in nm
m2     = 12789.3927072494     # Li mass in a.u
