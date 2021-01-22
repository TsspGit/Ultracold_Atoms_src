__author__ = "@Tssp"
__date__ = "12/01/21"
atom   = 'Li7'
m      = 7.0160040                                        # mass in Da
n      = 6                                                # Expansion order
Ix     = 7802                                             # Laser intensity in mW/cm2 X
Iy     = 4993                                             # Laser intensity in mW/cm2 Y
Iz     = 4993                                             # Laser intensity in mW/cm2 Z
wL     = 1000                                             # Wavelength in nm
hbar   = 1                                                # Planck reduced Constant
alpha  = 200                                              # Polarization a.u
h      = 1e-6                                             # Reference integration step
model  = 'Numerical'                                      # Energies computation: - Perturbation computes the energies using the 1st order perturbation theory
														  #                       - Numerical reads the energies
Erm    = 8.250619644654524e-11                           # (0,0,0) relative motion energy (Config False)
ECM    = 8.250619644654524e-11                           # (0,0,0) center of mass motion energy (Config False)
En_CM  = 3.178263170132861e-10                           # (0,n,0) and (n,0.0) free energies if model=Numerical DVR
Config = True                                             # Use config energies with model=Numerical
Econfig= 1.615698448144664e-10                           # DVR input if model=numerical
E_ICIR = 3.603054845215324                               # Energies Resonances
mode   = 'W'                                              # Compute mode: - C computes the constant of C*hbar*wz
							                              #               - W write your own C below
C      = 0                                                # Write value if mode=W