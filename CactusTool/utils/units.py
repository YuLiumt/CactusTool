"""
This module provides a class Units representing unit systems or unit 
conversions. Units can be used to convert from geometrized units to SI.
"""
import numpy as np
import re


# The following constants are all given in SI units from astropy - codata2018.py
h_SI          = 6.62607015e-34        # Planck constant [J s]
hbar_SI       = h_SI / (2 * np.pi)  # H_SI / (2 pi)
k_B_SI        = 1.380649e-23          # Boltzmann constant [J/K]
c_SI          = 299792458.0           # Vacuum speed of light m/s
G_SI          = 6.67430e-11           # Gravitational constant m^3/(kg s^2)
m_e_SI        = 9.1093837015e-31      # Electron mass [kg]
m_p_SI        = 1.67262192369e-27     # Proton mass [kg]
m_n_SI        = 1.67492749804e-27     # Neutron mass [kg]
u_SI          = 1.66053906660e-27     # Atomic mass [kg]
e_SI          = 1.602176634e-19       # Electron charge [C]
N_A_SI        = 6.02214076e23         # Avogadro constant 1/mol
au_SI         = 1.49597870700e11      # Astronomical Unit [m]
pc_SI         = au_SI / np.radians(1. / 3600.) # Parsec [m] 
kpc_SI        = 1000. * pc_SI         # Kiloparsec [m]
L_bol0_SI     = 3.0128e28             # Luminosity for absolute bolometric magnitude 0 [W]
L_sun_SI      = 3.828e26              # Nominal solar luminosity [W]
GM_sun_SI     = 1.3271244e20          # Nominal solar mass parameter m3 / (s2)
M_sun_SI      = GM_sun_SI / G_SI      # Solar mass [kg]
R_sun_SI      = 6.957e8               # Solar radius [m]
lyr_SI        = 9.46073e+15           # Lightyear [m]
eV_SI         = 1.602176634e-19       # electron volt [J]

# The following constants are all given in CU units c = G = M_sun = 1
c_CU          = 1                     # Vacuum speed of light
G_CU          = 1                     # Gravitational constant
M_CU      = 1                     # Solar mass 
MeV_CU        = 1                     # MeV

# The following convertion factors are all from CU to SI
m      = M_sun_SI * G_SI / c_SI **2 # meter
cm     = 100 * m                    # centimeter
km     = m / 1000                   # kilometer
au     = m / au_SI                  # Astronomical Unit
pc     = m / pc_SI                  # Parsec
kpc    = m / kpc_SI                 # Kiloparsec
lyr    = m / lyr_SI                 # Lightyear
R_sun  = m / R_sun_SI               # Solar radius
s      = m / c_SI                   # second
ms     = 1000 * s                   # millisecond
kg     = M_sun_SI                   # kilogram
g      = 1000 * kg                  # gram
m_e    = kg / m_e_SI                # Electron mass
m_p    = kg / m_p_SI                # Proton mass
m_n    = kg / m_n_SI                # Neutron mass
u      = kg / u_SI                  # Atomic mass


_subs_re = [
    (r"([\w\.\-\+\*\\\^])\s+", r"\1 "),  # merge multiple spaces
    # (
    #     r"\b([0-9]+\.?[0-9]*)(?=[e|E][a-zA-Z]|[a-df-zA-DF-Z])",
    #     r"\1*",
    # ),  # Handle numberLetter for multiplication
    (r"([\w\.\-])\s+(?=\w)", r"\1*"),  # Handle space for multiplication
]

#: Compiles the regex and replace {} by a regex that matches an identifier.
_subs_re = [(re.compile(a.format(r"[_a-zA-Z][_a-zA-Z0-9]*")), b) for a, b in _subs_re]


def UnitConversion(input_string):
    """
    Parse a units expression. 
    
    The expression can only contain products, ratios and powers of units.
    
    :param str input_string: units expression
    """
    # Sanitize input_string with whitespaces.
    input_string = input_string.strip()
    
    for a, b in _subs_re:
        input_string = a.sub(b, input_string)
        
    # Handle caret exponentiation
    input_string = input_string.replace("^", "**")
    
    return float(eval(input_string))