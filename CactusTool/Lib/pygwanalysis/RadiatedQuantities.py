#!/usr/bin/python

#----------------------------------------------------------
#
# This module computes radiated gravitational-wave energy,
# angular/linear momentum.
#
# This module depends on DataAnalysis and DiscreteFunction
#
# Released under the MIT License.
# (C) Christian Reisswig 2010-2011
#
#----------------------------------------------------------


import sys
import os
import glob
from .DiscreteFunction import *
from .WaveFunction import *
from .harmonics import *
from .DataAnalysis import *


def alm(l, m):
    return sqrt((l-m)*(l+m+1.0))*1.0/(l*(l+1.0))

def blm(l, m):
    return sqrt((l-2.0)*(l+2.0)*(l+m)*(l+m-1.0)*1.0/((2.0*l-1.0)*(2.0*l+1.0)))/(2.0*l)

def clm(l, m):
    return m*2.0/(l*(l+1.0))

def dlm(l, m):
    return 1.0/l * sqrt((l-2.0)*(l+2.0)*(l-m)*(l+m)*1.0/((2.0*l-1.0)*(2.0*l+1.0)))

def flm(l, m):
    return sqrt(l*(l+1.0)-m*(m+1.0))
    
    

def RadiatedEnergyFluxRWZM(Qeven, Qodd, lmax):
    """ Given the time integral of rPsi4 (l,m) modes (a DiscreteFunction mode array),
        this function computes the radiated energy flux
        that is emitted as a function of time in units of [M]."""
    
    npoints = Qeven[2][0].Length()
    Edot = DiscreteFunction(Qeven[2][0].x, zeros(npoints, float))

    dtQeven = InitModeArray(2)
    for l in range(2, lmax+1):
        for m in range(-l, l+1):
            dtQeven[l][m] = Qeven[l][m].DifferentiateFunction()
    
    for ii in range(Edot.Length()):
        e = 0.0
        for l in range(2, lmax+1):
            for m in range(-l, l+1):
                e = e + abs(dtQeven[l][m].f[ii])**2 + abs(Qodd[l][m].f[ii])**2
        Edot.f[ii] = e/(32.*pi)
    
    return Edot



def RadiatedEnergyFlux(int_rPsi4_ModeArray, lmax):
    """ Given the time integral of rPsi4 (l,m) modes (a DiscreteFunction mode array),
        this function computes the radiated energy flux
        that is emitted as a function of time in units of [M]."""
    
    t = []
    Edot = []
    
    for ii in range(int_rPsi4_ModeArray[2][2].Length()):
        e = 0.0
        for l in range(2, lmax+1):
            for m in range(-l, l+1):
                e = e + abs(int_rPsi4_ModeArray[l][m].f[ii])**2
        e = e/(16.*pi)
        Edot.append(e)
        t.append(int_rPsi4_ModeArray[2][2].x[ii])
    
    return DiscreteFunction(t, Edot)



def RadiatedAngMomFlux(int_rPsi4_ModeArray, int_int_rPsi4_ModeArray, lmax):
    """ Given the first and second time integral of rPsi4 (l,m) modes (DiscreteFunction mode arrays),
        this function computes the radiated angular momentum flux (Lx, Ly, Lz)
        that is emitted as a function of time in units of [M]."""
    t = []
    Lxdot = []
    Lydot = []
    Lzdot = []
    
    for ii in range(int_rPsi4_ModeArray[2][2].Length()):
        lx_dt = complex(0,0)
        ly_dt = complex(0,0)
        lz_dt = complex(0,0)
        for l in range(2, lmax+1):
            for m in range(-l, l+1):
                factorx = 0
                factory = 0
                if m+1 <= l:
                    factorx = factorx + flm(l, m) * int_rPsi4_ModeArray[l][m+1].f[ii].conjugate()
                    factory = factory + flm(l, m) * int_rPsi4_ModeArray[l][m+1].f[ii].conjugate()
                if m-1 >= -l:
                    factorx = factorx + flm(l,-m) * int_rPsi4_ModeArray[l][m-1].f[ii].conjugate()
                    factory = factory - flm(l,-m) * int_rPsi4_ModeArray[l][m-1].f[ii].conjugate()
                        
                lx_dt = lx_dt + int_int_rPsi4_ModeArray[l][m].f[ii] * factorx
                ly_dt = ly_dt + int_int_rPsi4_ModeArray[l][m].f[ii] * factory
                lz_dt = lz_dt + m * int_int_rPsi4_ModeArray[l][m].f[ii] * int_rPsi4_ModeArray[l][m].f[ii].conjugate()
            
        Lxdot.append( 1./(32*pi) * (lx_dt).imag)
        Lydot.append(-1./(32*pi) * (ly_dt).real)
        Lzdot.append( 1./(16*pi) * (lz_dt).imag)
        t.append(int_rPsi4_ModeArray[2][2].x[ii])
    
    return DiscreteFunction(t, Lxdot), DiscreteFunction(t, Lydot), DiscreteFunction(t, Lzdot)
    

def RadiatedLinMomFlux(int_rPsi4_ModeArray, lmax):
    """ Given the first time integral of rPsi4 (l,m) modes (DiscreteFunction mode array),
        this function computes the radiated linear momentum flux (Px, Py, Pz)
        that is emitted as a function of time in units of [M]."""
    t = []
    Pxdot = []
    Pydot = []
    Pzdot = []
    
    for ii in range(int_rPsi4_ModeArray[2][2].Length()):
        px_py_dt = complex(0, 0)
        pz_dt = complex(0, 0)
        for l in range(2, lmax+1):
            for m in range(-l, l+1):
                factor = complex(0.0, 0.0)
                factor2 = complex(0.0, 0.0)
                if m+1 <= l: 
                    factor = factor + alm(l,m)*int_rPsi4_ModeArray[l][m+1].f[ii].conjugate()
                #else: continue
                if l-1 >= 2:
                    if m+1 <= l-1: 
                        factor = factor + blm(l,-m)*int_rPsi4_ModeArray[l-1][m+1].f[ii].conjugate()
                    #else: continue
                    factor2 = factor2 + dlm(l,m)*int_rPsi4_ModeArray[l-1][m].f[ii].conjugate()
                #else: continue
                if l+1 <= lmax: 
                    factor = factor - blm(l+1,m+1)*int_rPsi4_ModeArray[l+1][m+1].f[ii].conjugate()
                    factor2 = factor2 + dlm(l+1,m)*int_rPsi4_ModeArray[l+1][m].f[ii].conjugate()
                #else: continue
        
                px_py_dt = px_py_dt + int_rPsi4_ModeArray[l][m].f[ii] * factor
                pz_dt = pz_dt + int_rPsi4_ModeArray[l][m].f[ii] * (clm(l,m)*int_rPsi4_ModeArray[l][m].f[ii].conjugate() + factor2)
    
        px_py_dt = 1./(8*pi) * px_py_dt
        pz_dt = 1./(16*pi) * pz_dt
        
        Pxdot.append(px_py_dt.real)
        Pydot.append(px_py_dt.imag)
        Pzdot.append(pz_dt.real)
        t.append(int_rPsi4_ModeArray[2][2].x[ii])
    
    return DiscreteFunction(t, Pxdot), DiscreteFunction(t, Pydot), DiscreteFunction(t, Pzdot)



