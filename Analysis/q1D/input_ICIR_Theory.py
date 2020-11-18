__author__ = "@Tssp"
__date__ = "01/10/20"
atom   = 'Li7'
m      = 7.0160040            # mass in Da
n      = 6                    # Expansion order
Ix     = 7190                 # Laser intensity in mW/cm2 X
Iy     = 4993                 # Laser intensity in mW/cm2 Y
Iz     = 50                   # Laser intensity in mW/cm2 Z
wL     = 1000                 # Wavelength in nm
hbar   = 1                    # Planck reduced Constant
alpha  = 200                  # Polarization a.u
h      = 1e-6                 # Reference integration step
mode   = 'C'                  # Compute mode: - C computes the constant of C*hbar*wz
							  #               - W write your own C below
C      = 0                    # Write value if mode=W