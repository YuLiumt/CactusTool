#----------------------------------------------------------
#
# Classes and routines for representing and manipulating 
# discrete wave functions.
#
# Released under the MIT License.
# (C) Christian Reisswig 2009-2011
#
#----------------------------------------------------------

import sys

from math import *
from numpy import *
from .DiscreteFunction import *

debug = 0

class WaveFunction(DiscreteFunction):
    """Class describing discrete 1d wave functions"""

    def __init__(self, x_, f_, g_=None):
        """Initialise a 1D wave function"""
        self.x = array(x_)
        if g_ != None:
            self.f = array(f_) + 1j*array(g_)
        else:
            self.f = array(f_)
        if self.HasUniformSpacing() == True:
            self.dx = self.x[1]-self.x[0]
        else:
            self.dx = 0.0

    def Phase(self):
        """Returns the phase as a DiscreteFunction object """
        
        f = unwrap(angle(self.f))
        res = DiscreteFunction(self.x, f)
        
        #res = DiscreteFunction(zeros(self.Length(), float), zeros(self.Length(), float))
        
        #phasecounter = 0
        #for ii in range(0, self.Length()):
        #    res.x[ii] = self.x[ii]
        #    
        #    phase = 0.0
        #    if (self.f[ii].real != 0):
        #       phase = atan(self.f[ii].imag/self.f[ii].real)
        #       sign = 1
        #       if ii != 0: 
        #           if res.f[ii-1] != 0:
        #               sign = res.f[ii-1]/abs(res.f[ii-1])
        #    else:
        #        phase = pi/2.0
        #    phase += phasecounter*pi
        #    
        #    if (ii > 0):
        #       if (abs(phase-res.f[ii-1]) > pi/2):
        #           phasecounter += 1 * sign
        #           phase += pi * sign
        #    
        #    res.f[ii] = phase
            
        # find jumps and remove them
        #counter = 0
        #f_ = []
        #for ii in range(0, self.Length()):
        #    if ii > 0:
        #       c = int((res.f[ii]-res.f[ii-1])/pi) 
        #       if c != 0:
        #           counter += c
        #    f_.append(res.f[ii]+counter*pi)
        #res.f = unwrap(res.f])
        
        return res
            

    def Amplitude(self):
        """Returns the amplitude as a DiscreteFunction object """
        
        res = DiscreteFunction(zeros(self.Length(), float), zeros(self.Length(), float))
        
        for ii in range(0, self.Length()):
            res.x[ii] = self.x[ii]
            res.f[ii] = sqrt(self.f[ii].real**2 + self.f[ii].imag**2)
        
        return res

    def Real(self):
        """Returns the real part of the wave"""
        res = DiscreteFunction(zeros(self.Length(), float), zeros(self.Length(), float))
        
        for ii in range(0, self.Length()):
            res.x[ii] = self.x[ii]
            res.f[ii] = self.f[ii].real
        
        return res
    
    def Imag(self):
        """Returns the imaginary part of the wave"""
        res = DiscreteFunction(zeros(self.Length(), float), zeros(self.Length(), float))
        
        for ii in range(0, self.Length()):
            res.x[ii] = self.x[ii]
            res.f[ii] = self.f[ii].imag
        
        return res
    
    
    def Frequency(self):
        """Returns the instantenous frequency as a DiscreteFunction object """
        
        res = self.Phase()
        res = res.FirstDerivative()
        
        return res
    





