__author__ = "@Tssp"
__date__ = "01/10/20"
atom   = 'Li7'
m      = 7.0160040                                        # mass in Da
n      = 6                                                # Expansion order
Ix     = 4993                                             # Laser intensity in mW/cm2 X
Iy     = 4993                                             # Laser intensity in mW/cm2 Y
Iz     = 50                                               # Laser intensity in mW/cm2 Z
wL     = 1000                                             # Wavelength in nm
hbar   = 1                                                # Planck reduced Constant
alpha  = 200                                              # Polarization a.u
h      = 1e-6                                             # Reference integration step
model  = 'Numerical'                                      # Energies computation: - Perturbation computes the energies using the 1st order perturbation theory
														  #                       - Numerical reads the energies
En_CM  = [1.4942847525476374e-10, 1.4942847525476374e-10] # (0,2,0) and (2,0.0) free energies if model=Numerical DVR
Config = True                                             # Use config energies with model=Numerical
Econfig = 1.160797099979224e-10                           # DVR input if model=numerical
mode   = 'C'                                              # Compute mode: - C computes the constant of C*hbar*wz
							                              #               - W write your own C below
C      = 2                                                # Write value if mode=W