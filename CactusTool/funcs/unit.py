"""
This module provides a class Units representing unit systems or unit 
conversions. All expressed in SI units.
"""

import math


# The following constants are all given in SI units
c          = 299792458.0           # Vacuum speed of light m/s
G          = 6.673e-11             # Gravitational constant m^3/(kg s^2)
M_sun      = 1.98892e30            # Solar mass kg
eV         = 1.602176
565e-19       # Electron volt
MeV        = 1e6 * eV
KB         = 1.3806488e-23         # Boltzmann constant [J/K]
M_e        = 9.10938291e-31        # Electron mass
M_p        = 1.672621777e-27       # Proton mass
M_n        = 1.674927351e-27       # Neutron mass
Lightyear  = 9460730472580800.0    # Lightyear
pc         = 30.856776e15          # Parsec
NA         = 6.02214129e23         # Avogadro constant 1/mol
h          = 6.62606957e-34        # Planck constant [J s]
h_bar      = 1.054571726e-34       # H_SI / (2 pi)

# ------------------------------------------------------------
# All expeession in geometric units refer to unit mass 
# where G=c=M_sun=1
# ------------------------------------------------------------
# convertion factors
M_to_km   = M_sun*G/(1000*c*c)                  # km
M_to_ms   = (1000*M_sun*G/(c*c*c))              # ms
M_to_dens = c*c*c*c*c*c / (G*G*G * M_sun*M_sun) # kg/m^3
M_to_dens_CGS = M_to_dens *1000/100**3         # g/cm^3
