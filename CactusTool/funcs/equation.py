def chirp_mass(m1, m2):
    """
    The chirp mass of the binary
    
    \mathcal{M}_{\mathrm{chirp}}=\left(M_{1} M_{2}\right)^{3 / 5}\left(M_{1}+M_{2}\right)^{-1 / 5}
    
    :param float m1: gravitational mass
    :param float m2: gravitational mass
    """
    return (m1*m2)**(3./5) * (m1+m2)**(-1./5)

def mass_ratio(m1, m2):
    """
    The binary mass ratio
    
    q = M_{2} / M_{1}
    
    :param float m1: gravitational mass
    :param float m2: gravitational mass
    """
    return m2/m1

def gravitational_mass(M, q):
    """
    m_{1}=\mathscr{M}(1+q)^{1 / 5} q^{-3 / 5}
    m_{2}=m_{1}*q
    """ 
    m1 = M * (1 + q)**(1./5) * q**(-3./5)
    return m1, m1*q


def baryonic_mass(m):
    """
    Relations between NS gravitational mass and baryonic mass

    M_{b}=M_{g}+A_{1} \times M_{g}^{2}+A_{2} \times M_{g}^{3}
    
    Putting all EOSs together, the overall best fit A1 and A2 values are 0.0729 and 0.0032.
    ref:arXiv:1905.03784
    
    :param float m: gravitational mass
    """
    return m + 0.0729*m**2 + 0.0032*m**3