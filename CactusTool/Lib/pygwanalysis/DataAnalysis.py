#!/usr/bin/python

#----------------------------------------------------
#
# Routines and algorithms for GW data analysis.
# 
# Released under the MIT Licence.
# (C) Christian Reisswig 2009-2011
# 
#----------------------------------------------------


import sys
import os
import glob
from .DiscreteFunction import *
from .WaveFunction import *
from .harmonics import *

constant_c = 299792458       # speed of light in [m/sec]
constant_G = 6.67428 * 1e-11 # gravitational constant in [m^3 kg^-1 s^-2]
constant_H0 = 70.0         # Hubble constant

Msol = 1.9891*10**30                                # solar mass[kg]
ConversionFac = constant_G/constant_c**2            # [m/kg] (mass faktor G/c^2 to transform from geometrized units to SI units)
Msec = Msol*ConversionFac/constant_c                #  one solar mass in [sec]
Mm = Msol*ConversionFac                             #  one solar mass in[m]
MmMsec = Msec * Mm
Mpc = 3.08568025*10**16 * 10**6    # 1 Mpc in[m] (~3.26 *10^6 ly)

path = os.path.dirname(__file__)+"/"
#path = "/home/reisswig/scripts/DataAnalysis/"

AdLIGOcurve = DiscreteFunction([],[])
AdLIGOcurve.Load(path+"AdLIGO_ZERO_DET_high_P.dat")

eLIGOcurve = DiscreteFunction([],[])
eLIGOcurve.Load(path+"eligo.dat")

AdVirgocurve = DiscreteFunction([],[])
AdVirgocurve.Load(path+"AdVrefsens.dat")


def FStrainToSI(h, M):
    """Converts the Fourier domain strain to SI units"""
    res = WaveFunction(h.x * 1/(M*Msec), h.f * M**2 * MmMsec)
    return res

def StrainToSI(h, M):
    """ Converts time domain strain to SI units"""
    res = WaveFunction(h.x * M*Msec, h.f * Mm*M)
    return res

def StrainToCGS(h, M):
    """ Converts time domain strain to CGS units"""
    res = WaveFunction(h.x * M*Msec, h.f * 1.4772e5*M)
    return res

def Psi4ToSI(h, M):
    """ Converts time domain psi4 to SI units"""
    res = WaveFunction(h.x * M*Msec, h.f * Mm*M / (M*Msec)**2)
    return res

def NewsToSI(h, M):
    """ Converts time domain news to SI units"""
    res = WaveFunction(h.x * M*Msec, h.f * Mm*M / (M*Msec))
    return res


def TimesToSI(t, M):
    """Converts a list of time values to SI units"""
    res = array(map(lambda x: x * M * Msec, t))
    return res

def FrequenciesToSI(f, M):
    """Converts a list of frequencies to SI units"""
    res = array(map(lambda x: x * 1 / (M*Msec), f))
    return res

def SelecthPlus(h):
    """Given a complex function h=h_+-ih_x, this returns h_+"""
    hplus = WaveFunction(h.x, h.f.real)
    return hplus

def SelecthCross(h):
    """Given a complex function h=h_+-ih_x, this returns h_x"""
    hcross = WaveFunction(h.x, -h.f.imag)
    return hcross

def MinMassFromOmega(omega, fcut):
    """Given a frequency omega=2*pi*f in units 1/M and a lower detector sensitivity cut-off fcut [Hz], 
       this computes the minimum allowable mass such that fSI>=fcut, where fSI is f in SI units."""
    return omega / (2*pi*fcut*Msec)

def GaussSmooth(WF, sigma):
    """Smoothes the input by convolving with a normal distribution"""
    res = WaveFunction(WF.x, WF.f*0)
    def Gauss(t, t0, sigma0):
        return 1.0/(sqrt(2*pi)*sigma)*exp(-(t-t0)**2/(2*sigma0**2))
    #WF.dx = WF.x[1]-WF.[0]
    nkernel_points = int(3*sigma / WF.dx)
    dx = WF.dx
    print(nkernel_points)
    x = zeros(2*nkernel_points+1, float)
    for ii in range(0, len(x)):
        x[ii] = -dx*nkernel_points + ii*dx
    G = map(lambda x: Gauss(x, 0, sigma), x)
    for ii in range(nkernel_points, WF.Length()-nkernel_points):
        for jj in range(-nkernel_points, nkernel_points+1):
            res.f[ii] += WF.f[ii+jj]*G[jj+nkernel_points]  #*Gauss(WF.x[ii-jj], WF.x[ii], sigma)
    return res

def Reconstruction(WF, dx, kernel, tau):
    """Resamples the input to resolution dx by convolving with a reconstruction filter kernel of width tau"""
    npoints = int((WF.x[-1]-WF.x[0]) / dx) + 1
    res = WaveFunction(zeros(npoints, float), zeros(npoints, type(WF.f[0])))
    nkernel_points = int(tau / WF.dx) # number of points for the current resampled point necessary for the reconstruction from the orginal samples
    print("nkernel_points = %d" % nkernel_points)
    
    for ii in range(0, res.Length()):
        res.x[ii] = WF.x[0] + ii*dx
        # find nearest neighbor to current point x in old data
        kk = WF.FindNearestNeighbor(res.x[ii])
        # convolve
        norm = 0
        for jj in range(-nkernel_points, nkernel_points+1):
            if (kk+jj >= 0 and kk+jj < WF.Length()):
                norm += kernel(res.x[ii]-WF.x[kk+jj])
                res.f[ii] += WF.f[kk+jj] * kernel(res.x[ii]-WF.x[kk+jj])
        res.f[ii] /= norm
    return res



def Spectrogram(WF, kernel, tau, sigma, t0, t1, n):
    """ Constructs a spectogram of function WF, i.e. computes the frequencies components 
        for n time bins between t0 and t1. The result is a 2D function along time and frequency
        axis. 'kernel' and 'sigma' denote the window function and half-width that is used for 
        computing the Fourier transform for a given time bin of width 'tau'. 'tau' is usually chosen such
        that the bins overlap."""
    times = []
    spectra = []
    
    # Time resolution of WF
    dx = WF.dx
    # number of kernel points for the current time
    nkernel_points = int(tau / dx / 2) 
    print("Kernel points : %d" % nkernel_points)
    # bin (cell) width as it appears in t-axis of spectogram
    dt = float(t1-t0) / float(n)
    # start at cell midpoints
    t0 = t0 + dt/2
    for ii in range(0,n):
        times.append(t0+ii*dt)
        x = []
        f = []
        # find nearest neighbor to current point x in old data
        kk = WF.FindNearestNeighbor(times[ii])
        # apply window to current time bin of width tau
        print("time = %g" % times[ii])
        for jj in range(-nkernel_points, nkernel_points+1):
            if (kk+jj >= 0 and kk+jj < WF.Length()):
                x.append(WF.x[kk+jj])
                f.append(WF.f[kk+jj] * kernel(WF.x[kk+jj]-times[ii], sigma))
            else:
                print("Error: cannot construct proper window: t = %g" % times[ii])
                x.append(WF.x[0] + WF.dx*(kk+jj))
                f.append(0)
        windowedF = DiscreteFunction(x, f)
        spectra.append(windowedF.FourierTransform(True))
        spectra[ii].f = abs(spectra[ii].f)**2
    res = DiscreteFunction(times, spectra)
    return res
    


def TanhWindow(WF, t0, sigma0, t1, sigma1):
    """Applies a tanh filter to the waveform where t0, t1 denote the blending positions and sigma0 and sigma1 the widths of blending"""
    res = WaveFunction(WF.x, WF.f)
    res.f = WF.f * 0.25 * (1.0 + tanh(4.0*(WF.x-t0)/sigma0)) * (1.0 - tanh(4.0*(WF.x-t1)/sigma1))
    return res

def GaussFilter(WF, t0, sigma0, sign):
    """Applies a (half-)Gaussian filter to the waveform where t0, denotes the blending position and sigma0 the width of blending. 
       sign determines which of the two halfes is chosen."""
    res = WaveFunction(WF.x, WF.f)
    def HalfGauss(t, t0, sigma0, sign):
        if sign > 0 and t >= t0: return 1.0
        if sign < 0 and t <= t0: return 1.0
        return exp(-(t-t0)**2/(2*sigma0**2))
    blend = map(lambda x: HalfGauss(x, t0, sigma0, sign), WF.x)
    res.f = WF.f * blend
    return res

def ButterworthFilter(WF, G0, omega0, n, sign):
    """Applies a "Butterworth" filter to the waveform where omega0 denotes the blending position and n the order of blending. 
       sign determines which of the two halfes is chosen."""
    res = WaveFunction(WF.x, WF.f)
    
    blend = map(lambda x: sqrt(G0**2/(1+(WF.x/omega0)**(2*n))), WF.x)
    res.f = WF.f * blend
    return res


def ReisswigFilter(WF1, f0, doItInFreqDomain=False, k=1, n=0):
    """ This filter takes the k-th derivative (optionally in the Fourier domain) and then applies Fixed-frequency integration
        to hi-pass filter the signal """
    WF = WaveFunction(WF1.x, WF1.f) #.conjugate())
    if doItInFreqDomain==False:
        if k == 1: WF = WF.FirstDerivative()
        if k == 2: WF = WF.SecondDerivative()
        if k == 4: WF = WF.SecondDerivative().SecondDerivative()
        return FixedFrequencyIntegration(WF, f0, k)
        
    fr = WF.Extend(n*WF.Length(), WF.dx, lambda x: 0.0).FourierTransform()

    frange = [f0, f0]

    for ii in range(0, fr.Length()):
        if (fr.x[ii] != 0):
            if fr.x[ii] < 0: fr.f[ii] = 0
            div = 1.0
            if fr.x[ii] < frange[0]: div = frange[0]/fr.x[ii]
            fr.f[ii] = -1.0 * fr.f[ii]/(div)

        else:
            fr.f[ii] = 0
    res = fr.InverseFourierTransform()
    return res



def FixedFrequencyIntegration(WF1, f0, k, n=0):
    """Given some function and frequency f0 (ie omega0/2pi)), we apply the integration scheme proposed by Reisswig&Pollney
       to get the k-th integral"""
    # Put lots of zeros after end of signal in time domain to increase resolution of signal in in frequency domain
    # and afterwards Fourier-transform
    WF = WaveFunction(WF1.x, WF1.f)
    fr = WF.Extend(n*WF.Length(), WF.dx, lambda x: 0.0).FourierTransform()

    for ii in range(0, fr.Length()):
        if (fr.x[ii] != 0):
            #if fr.x[ii] < 0: fr.f[ii] = 0
            div = 2*pi*fr.x[ii]
            if abs(fr.x[ii]) < abs(f0): div = 2*pi*abs(f0)*sign(fr.x[ii])
            fr.f[ii] = complex(0.,1.)**k * fr.f[ii]/(div**k)
            
        else:
            fr.f[ii] = 0

    res = fr.InverseFourierTransform()    
        
    return res



def GethFromPsi4(WF1, f0, n=0):
    """Given psi4 and frequency f0 (ie omega0/2pi)), we apply the integration scheme proposed by Reisswig&Pollney
       to get h. n denotes the number of zeros as n times length of signal"""
    # Put lots of zeros after end of signal in time domain to increase resolution of signal in in frequency domain
    # and afterwards Fourier-transform
    WF = WaveFunction(WF1.x, WF1.f.conjugate())
    fr = WF.Extend(n*WF.Length(), WF.dx, lambda x: 0.0).FourierTransform()

    #frange = [f0, f0+5e-4]
    frange = [f0, f0]

    for ii in range(0, fr.Length()):
        if (fr.x[ii] != 0):
            #if fr.x[ii] < 0: fr.f[ii] = 0
            div = 4*pi**2*fr.x[ii]**2
            if abs(fr.x[ii]) < abs(frange[0]): div = 4*pi**2*frange[0]**2
            #if abs(fr.x[ii]) < abs(frange[0]): div = 4*pi**2*frange[0]**2
            #if abs(fr.x[ii]) > abs(frange[0]) and abs(fr.x[ii]) < abs(frange[1]): #div = (1.0-GaussBlend(DF1fr.x[ii], frange[0], (frange[1]-frange[0])/100)) * div + GaussBlend(DF1fr.x[ii], frange[0], (frange[1]-frange[0])/100) * 4*pi**2*frange[0]**2
            #    div = (1 - (abs(fr.x[ii])-abs(frange[0]))/(abs(frange[1]-frange[0]))) * 4*pi**2*frange[0]**2 + (abs(fr.x[ii])-abs(frange[0]))/(abs(frange[1]-frange[0])) * div
            fr.f[ii] = -fr.f[ii]/(div)
            
        else:
            fr.f[ii] = 0
    
    res = fr.InverseFourierTransform()    
    
    # Shift peak to zero
    res.x -= res.x[res.FindMax()]
    
    return res


def GethFromNews(WF1, f0, n=0):
    """Given the news and frequency f0 (ie omega0/2pi)), we apply the integration scheme proposed by Reisswig&Pollney
       to get h"""
    # Put lots of zeros after end of signal in time domain to increase resolution of signal in in frequency domain
    # and afterwards Fourier-transform
    WF = WaveFunction(WF1.x, WF1.f.conjugate())
    fr = WF.Extend(n*WF.Length(), WF.dx, lambda x: 0.0).FourierTransform()

    #frange = [omega0-1e-4, omega0]
    frange = [f0, f0]

    for ii in range(0, fr.Length()):
        if (fr.x[ii] != 0):
            #if fr.x[ii] < 0: fr.f[ii] = 0
            div = 2*pi*fr.x[ii]
            if abs(fr.x[ii]) < abs(frange[0]): div = 2*pi*abs(frange[0])*sign(fr.x[ii])
            #if fr.x[ii] > frange[1]: div = 4*pi**2*frange[1]**2
            #if DF1fr.x[ii] > frange[0] and DF1fr.x[ii] < frange[1]: #div = (1.0-GaussBlend(DF1fr.x[ii], frange[0], (frange[1]-frange[0])/100)) * div + GaussBlend(DF1fr.x[ii], frange[0], (frange[1]-frange[0])/100) * 4*pi**2*frange[0]**2
            #    div = (1 - (DF1fr.x[ii]-frange[0])/(frange[1]-frange[0])) * 4*pi**2*frange[0]**2 + (DF1fr.x[ii]-frange[0])/(frange[1]-frange[0]) * div
            fr.f[ii] = complex(0.,1.)*fr.f[ii]/(div)
            
        else:
            fr.f[ii] = 0
    #plt.loglog(fr.x, abs(fr.f))

    res = fr.InverseFourierTransform()    
    
    # Shift peak to zero
    #res.x -= res.x[res.FindMax()]
    
    return res


def PolynomialFilterDriver(WFfr, x0, x1, blend=False):
    """Given a waveform in Fourier domain, we fit a polynomial between frequncy range frange as
       described in Reisswig&Pollney 2010"""
    fr = WaveFunction(WFfr.x, WFfr.f)
    
    b = 1.0*log(x1[1]/x0[1])/log(x1[0]/x0[0])
    a = x1[1]/(x1[0]**b)
    print("a = %s" % a)
    print("b = %s" % b)
    
    def GaussBlend(x, x0, s0):
        return exp(-(x-x0)**2/s0**2)
    x_ = []
    f_ = []
    x_.append(0.0)
    f_.append(0.0)
    for ii in range(0, fr.Length()):
        if fr.x[ii]-x1[0] < 0 and fr.x[ii] > 0:
            fac = abs(WFfr.f[ii])/(a*WFfr.x[ii]**b)
            if (blend == False or fr.x[ii] < x0[0]): fr.f[ii] /= fac
            if (blend == True and fr.x[ii] >= x0[0]): fr.f[ii] = (1.0-GaussBlend(fr.x[ii], x1[0], (x1[0]-x0[0])/10)) * fr.f[ii]/fac + GaussBlend(fr.x[ii], x1[0], (x1[0]-x0[0])/10)*fr.f[ii]
            #if (DF2.x[ii] >= 0 and DF2.x[ii] <= x0[0]):
            #    x_.append(DF2.x[ii])
            #    f_.append(fac)
        if (fr.x[ii] <= 0): fr.f[ii] = 0
    #H = WaveFunction(x_, f_)
    #return DF2, H 
    return fr


def GethFromPsi4Filtered(WF1, frange, n=0):
    """Given psi4 and frequency f0 (ie omega0/2pi)), we apply the integration scheme proposed by Reisswig&Pollney by applying a filter
       to get h. n denotes the number of zeros as n times length of signal"""
    # Put lots of zeros after end of signal in time domain to increase resolution of signal in in frequency domain
    # and afterwards Fourier-transform
    WF = WaveFunction(WF1.x, WF1.f.conjugate())
    fr = WF.Extend(n*WF.Length(), WF.dx, lambda x: 0.0).FourierTransform()
    frAbs = WaveFunction(fr.x, abs(fr.f))

    # Find slope
    x0 = [frange[0], frAbs.Interpolate(frange[0])]
    x1 = [frange[1], frAbs.Interpolate(frange[1])]
    
    fr = PolynomialFilterDriver(fr, x0, x1, True)

    for ii in range(0, fr.Length()):
        if (fr.x[ii] != 0):
            if fr.x[ii] < 0: fr.f[ii] = 0
            div = 4*pi**2*fr.x[ii]**2
            fr.f[ii] = -fr.f[ii]/(div)
        else:
            fr.f[ii] = 0
    
    res = fr.InverseFourierTransform()    
    
    # Shift peak to zero
    #res.x -= res.x[res.FindMax()]
    
    return res


def LuminosityDistance(z, Omega_m=0.3175, Omega_lambda=0.6825):
    """Computes the luminosity distance in Mpc for a given redshift z
       Code from Ed Wright http://www.astro.ucla.edu/~wright/CosmoCalc.html"""
    
    H0 = constant_H0
    WM = Omega_m # Omega(matter)
    WV = Omega_lambda # Omega(vacuum) or lambda
    
    WR = 0.        # Omega(radiation)
    WK = 0.        # Omega curvaturve = 1-Omega(total)
    c = 299792.458 # velocity of light in km/sec
    Tyr = 977.8    # coefficent for converting 1/H into Gyr
    DTT = 0.5      # time from z to now in units of 1/H0
    DTT_Gyr = 0.0  # value of DTT in Gyr
    age = 0.5      # age of Universe in units of 1/H0
    age_Gyr = 0.0  # value of age in Gyr
    zage = 0.1     # age of Universe at redshift z in units of 1/H0
    zage_Gyr = 0.0 # value of zage in Gyr
    DCMR = 0.0     # comoving radial distance in units of c/H0
    DCMR_Mpc = 0.0 
    DCMR_Gyr = 0.0
    DA = 0.0       # angular size distance
    DA_Mpc = 0.0
    DA_Gyr = 0.0
    kpc_DA = 0.0
    DL = 0.0       # luminosity distance
    DL_Mpc = 0.0
    DL_Gyr = 0.0   # DL in units of billions of light years
    V_Gpc = 0.0
    a = 1.0        # 1/(1+z), the scale factor of the Universe
    az = 0.5       # 1/(1+z(object))
    
    h = H0/100.
    WR = 4.165E-5/(h*h)   # includes 3 massless neutrino species, T0 = 2.72528
    WK = 1-WM-WR-WV
    az = 1.0/(1+1.0*z)
    age = 0.
    n=1000         # number of points in integrals
    for i in range(n):
        a = az*(i+0.5)/n
        adot = sqrt(WK+(WM/a)+(WR/(a*a))+(WV*a*a))
        age = age + 1./adot

    zage = az*age/n
    zage_Gyr = (Tyr/H0)*zage
    DTT = 0.0
    DCMR = 0.0
    
    def DCMT():
        ratio = 1.00
        x = sqrt(abs(WK))*DCMR
        if (x > 0.1):
            if (WK > 0): ratio = 0.5*(exp(x)-exp(-x))/x 
            else: ratio = sin(x)/x;
            y = ratio*DCMR
            return y
        y = x*x
        # statement below fixed 13-Aug-03 to correct sign error in expansion
        if (WK < 0): y = -y
        ratio = 1 + y/6 + y*y/120
        y= ratio*DCMR
        return y


    def VCM():
        ratio = 1.00
        x = sqrt(abs(WK))*DCMR;
        if (x > 0.1):
            if (WK > 0): ratio = (0.125*(exp(2*x)-exp(-2*x))-x/2)/(x*x*x/3)
            else: ratio = (x/2 - sin(2*x)/4)/(x*x*x/3)
            y = ratio*DCMR*DCMR*DCMR/3
            return y
        y = x*x
        # statement below fixed 13-Aug-03 to correct sign error in expansion
        if (WK < 0): y = -y
        ratio = 1 + y/5 + (2/105)*y*y
        y= ratio*DCMR*DCMR*DCMR/3
        return y
    
    
    h = H0/100
    WR = 4.165E-5/(h*h) # includes 3 massless neutrino species, T0 = 2.72528
    WK = 1-WM-WR-WV
    az = 1.0/(1+1.0*z)
    age = 0;
    for i in range (0, n+1):
        a = az*(i+0.5)/n
        adot = sqrt(WK+(WM/a)+(WR/(a*a))+(WV*a*a))
        age = age + 1/adot
    zage = az*age/n
    # correction for annihilations of particles not present now like e+/e-
    # added 13-Aug-03 based on T_vs_t.f
    lpz = log((1+1.0*z))/log(10.0)
    dzage = 0
    if (lpz >  7.500): dzage = 0.002 * (lpz -  7.500)
    if (lpz >  8.000): dzage = 0.014 * (lpz -  8.000) +  0.001
    if (lpz >  8.500): dzage = 0.040 * (lpz -  8.500) +  0.008
    if (lpz >  9.000): dzage = 0.020 * (lpz -  9.000) +  0.028
    if (lpz >  9.500): dzage = 0.019 * (lpz -  9.500) +  0.039
    if (lpz > 10.000): dzage = 0.048
    if (lpz > 10.775): dzage = 0.035 * (lpz - 10.775) +  0.048
    if (lpz > 11.851): dzage = 0.069 * (lpz - 11.851) +  0.086
    if (lpz > 12.258): dzage = 0.461 * (lpz - 12.258) +  0.114
    if (lpz > 12.382): dzage = 0.024 * (lpz - 12.382) +  0.171
    if (lpz > 13.055): dzage = 0.013 * (lpz - 13.055) +  0.188
    if (lpz > 14.081): dzage = 0.013 * (lpz - 14.081) +  0.201
    if (lpz > 15.107): dzage = 0.214
    zage = zage*10.0**dzage

    zage_Gyr = (Tyr/H0)*zage
    DTT = 0.0
    DCMR = 0.0
    # do integral over a=1/(1+z) from az to 1 in n steps, midpoint rule
    for i in range (0,n+1):
        a = az+(1-az)*(i+0.5)/n
        adot = sqrt(WK+(WM/a)+(WR/(a*a))+(WV*a*a))
        DTT = DTT + 1/adot
        DCMR = DCMR + 1/(a*adot)

    DTT = (1-az)*DTT/n
    DCMR = (1-az)*DCMR/n
    age = DTT+zage
    age_Gyr = age*(Tyr/H0)
    DTT_Gyr = (Tyr/H0)*DTT
    DCMR_Gyr = (Tyr/H0)*DCMR
    DCMR_Mpc = (c/H0)*DCMR
    DA = az*DCMT()
    DA_Mpc = (c/H0)*DA
    kpc_DA = DA_Mpc/206.264806
    DA_Gyr = (Tyr/H0)*DA
    DL = DA/(az*az)
    DL_Mpc = (c/H0)*DL
    DL_Gyr = (Tyr/H0)*DL
    V_Gpc = 4*pi*pow(0.001*c/H0,3)*VCM()
    
    return DL_Mpc


def ShiftToZero(DF, trange):
    """Given a range [.,.] shift the wave according to that intervall to zero"""
    s = 0
    n = 0
    if (trange[1] > DF.x[-1]): print("ShiftToZero: trange[1] out of bounds")
    if (trange[0] > DF.x[-1]): print("ShiftToZero: trange[0] out of bounds")
    for ii in range(0, DF.Length()):
        if DF.x[ii] >= trange[0] and DF.x[ii] <= trange[1]:
            s += DF.f[ii]
            n += 1
    return WaveFunction(DF.x, DF.f - s/n)


def GethFD(WF, wf_type='news'):
    """Given psi4 or the news in Fourier domain this returns h in Fourier domain"""
    res = WaveFunction(WF.x, WF.f)
    if wf_type.lower() == 'psi4':
        for ii in range(0, WF.Length()):
            if (WF.x[ii] != 0):
                res.f[ii] = -WF.f[ii]/(4*pi**2*WF.x[ii]**2)
    else:
        for ii in range(0, WF.Length()):
            if (WF.x[ii] != 0):
                res.f[ii] = complex(0.,1.)*WF.f[ii]/(2*pi*WF.x[ii])
    return res


def GethTD(WF, wf_type='news', shiftToZero=False):
    """Given psi4 or the news in the time domain this returns h in the time domain"""
    res = WaveFunction(WF.x, WF.f).IntegrateFunction()
    if (shiftToZero == True):
        print("Have to implement the shift-to-zero")
    if (wf_type.lower() == 'psi4'):
        res = res.IntegrateFunction()
        if (shiftToZero == True):
            print("Have to implement the shift-to-zero")
    return res
    
def InitModeArray(lmax):
    """Creates an array with empty modes up to lmax"""
    res = []
    for ll in range(0, lmax+1):
        res.append([])
        for mm in range(-ll, ll+1):
            res[-1].append([])
    return res

def LoadModeArray(basename, lmax, col_time=1, col_re=2, col_im=3, filenameformat="cce_conversion_script"):
    """Loads an array of modes up to lmax from disk.
       We assume that files begin with 'basename' and end with _l*m*.dat"""
    res = InitModeArray(lmax)
    
    for ll in range(2, lmax+1):
        for mm in range(-ll, ll+1):
            if (filenameformat == "cce_conversion_script"):
                res[ll][mm] = WaveFunction([], [])
                res[ll][mm].Load(basename+"_l"+str(ll)+"m"+str(mm)+".dat", col_time, col_re, col_im)
                continue
            if (filenameformat == "pittnull"):
                res[ll][-mm] = WaveFunction([], [])
                # Get spaces from Fortran outpout filenames right
                tmp1 = "%02d" % ll
                if mm < 0:
                    tmp2 = "m"
                else:
                    tmp2 = "p"
                tmp2 += "%02d" % abs(mm)
                # Take into account that nullcode has a weird convention for psi4!
                res[ll][-mm].Load(basename+".L"+tmp1+"M"+tmp2+".asc", col_time, col_re, col_im)
                res[ll][-mm].f *= 2.0*(-1.0)**(mm+1)
                res[ll][-mm].f = map(lambda x: x.conjugate(), res[ll][-mm].f)
                continue
            print("Error loading mode array: filenameformat not recognized!")
                
    return res


def ReconFromsYlm(arrayWF, theta, phi, lmax):
    """Given an array of WF modes h_lm reconstruct h at angle theta, phi up to l=lmax"""
    res = WaveFunction(arrayWF[2][0].x, arrayWF[2][0].f*0.0)
    
    for ll in range(2,lmax+1):
        for mm in range(-ll, ll+1):
            res.f += sYlm(-2, ll, mm, theta, phi) * arrayWF[ll][mm].f
            
    return res 
    

def SensitivityLIGO(f):
    f0 = 150.0
    s0 = 9.0e-46
    S = s0 * ((4.49*f/f0)**(-56.0) + 0.16*(f/f0)**(-4.52) + 0.52 + 0.32*(f/f0)**2)
    return S

def SensitivityeLIGO(f):
    return eLIGOcurve.Interpolate(f)

def SensitivityAdLIGO_Old(f):
    """This curve is an old one!"""
    f0 = 215.0
    S0 = 1.0e-49
    S = S0*( (f/f0)**(-4.14) - 5.0*(f0/f)**2 + 111.0*( (1.0 - (f/f0)**2 + 0.5*(f/f0)**4) / (1.0 + 0.5*(f/f0)**2 ) ) )
    return S

def SensitivityAdLIGO(f):
    """This curve is using most likely data points from 25. Jan. 2010, LIGO Document T0900288-v3"""
    # web-adress: https://dcc.ligo.org/cgi-bin/DocDB/ShowDocument?docid=2974 
    return AdLIGOcurve.Interpolate(f)**2


def SensitivityVirgo(f):
    s0 = 10.2e-46
    f0 = 500.0
    S = s0 *( (7.87*f/f0)**(-4.8) + 6.0/17.0 * f0/f + 1 + (f/f0)**2 )    
    return S


def SensitivityAdVirgo(f):
    return AdVirgocurve.Interpolate(f)**2


def UnitStep(x):
    """Step function"""
    return 1.0 if (x > 0) else 0.0

def SensitivityLISA(F):
    # code due to Miquel Trias, see LISA PE Wiki, \
    # http://www.tapir.caltech.edu/dokuwiki/lisape:home
    # We used this curve for the Detection paper.
    
    Fst = 9.542690318473885*10**(-3);
    Larm = 5*10**9;
    SPS = 4*10**(-22);
    SACC = 9*10**(-30);
    red = (2*pi*10**(-4))**2;
    
    n0 = 1.0/(4*Larm**2);
    n1 = 4.0*SPS;
    n2 = (16.0*SACC)/(2.0*pi*F)**4;
    n3 = (16.0*SACC*red)/(2.0*pi*F)**6;
    n4 = 2.0*(F/Fst)**2*SPS;
    
    conf = (10**-44.62*F**-2.3)*(1 - UnitStep(F - 10**-3)) + (10**-50.92*
          F**-4.4)*(UnitStep(F - 10**-3) - 
            UnitStep(F - 10**-2.7)) + (10**-62.8*
          F**-8.8)*(UnitStep(F - 10**-2.7) - 
            UnitStep(F - 10**-2.4)) + (10**-89.68*
          F**-20)*(UnitStep(F - 10**-2.4) - UnitStep(F - 10**-2));
    
    S = n0*(n1 + n2 + n3 + n4) + conf
    return S


def SensitivityLISA_Berti(f):
    # Noise curve from Eq. (2.31) of http://arxiv.org/abs/gr-qc/0411129v2
    # Suggested by Ajith. However, the other curve seems to be somewhat newer...
    # There is no real "official" LISA curve to date.
    
    S_NSA = (9.18*1e-52 * f**4 + 1.59*1e-41 + 9.18*1e-38 * f**2)
    S_gal = 2.1*1e-45 * f**(-7.0/3.0)
    S_ex_gal = 4.2*1e-47 * f**(-7.0/3.0)
    
    #T_mission = 
    #kappa = 4.5
    
    #S = min([S_NSA ]) + 
    S = S_NSA + S_gal + S_ex_gal
    
    return S


def SensitivityeLISA(f):
    # Noise curve for eLISA-NGO.
    # Function from eLISA yellow book, Sec. 4.3
    
    f0 = 150.0e-3
    a = 0.41
    T = sqrt(1.0+(f/(a*f0))**2)
    
    Da0 = 3.0
    fL = 0.1e-3
    fH = 8.0e-3
    da = Da0*1.0e-15 * sqrt(1.0+(f/fH)**4) * sqrt(1.0+fL/f)
    
    dx_DRS = 2.0*da / (2.0*pi*f)**2
    
    Dx0 = 12.0
    f0 = 2.8e-3
    dx_IMS = Dx0*1e-12 * sqrt(1.0+(f0/f)**4)
    
    # Is this correct??? No...
    S_IMS = dx_IMS
    S_DRS = dx_DRS
    
    # arm length
    L = 1e9
    
    #h = sqrt(5.0)*2.0/sqrt(3.0)*T*sqrt(S_IMS + S_DRS) / L
    # Use this in terms of dx_IMS, dx_DRS:
    h = sqrt(5.0)*2.0/sqrt(3.0)*T*(dx_IMS + dx_DRS) / L
    
    S = h**2
    
    return S


def SensitivityDECIGO(f):
    """ From arXiv:1110.2865
        For this sensitivity curve, you don't need to select the + or x polarization component since the detector is a triangle.
        The 60 degree abgle between detector arms is incoproated into the noise curve.
        The antenna pattern is incorporated in the noise curve.
        For SNRs, use an additional factor of sqrt(4) since BBO consists of 4 detectors!"""
    
    # White-dwarf confusion noise
    SWD = 4.2e-47 * f**(-7.0/3.0) * exp(-2.0 * (f/5e-2)**2)
    
    fc = 7.69
    S = 3.30 * 1e-50*1.0/f**4 + 3.09*1e-47 * (1.0 + (f/fc)**2) + SWD
    return S


def SensitivityBBO(f):
    """ From arXiv:1110.2865
        For this sensitivity curve, you don't need to select the + or x polarization component since the detector is a triangle.
        The 60 degree abgle between detector arms is incoproated into the noise curve.
        The antenna pattern is incorporated in the noise curve.
        For SNRs, use an additional factor of sqrt(4) since BBO consists of 4 detectors!"""
    
    # White-dwarf confusion noise
    SWD = 4.2e-47 * f**(-7.0/3.0) * exp(-2.0 * (f/5e-2)**2)
    
    S = 6.15*1e-51*(1.0/f)**4 + 1.95*1e-48 + 1.20*1e-48*f**2 + SWD
    return S


    
def SNR(fhSI, d, Sn, fLow, fHi):
    """Given a reconstructed h in SI units in Fourier domain and distance d (in Mpc) this computes the SNR for sensitivity curve Sn"""
    """Remember to select either h_+ or h_\cross! Otherise the SNR will be twice as large!"""
    S = map(Sn, fhSI.x) # create table of values for noise curve
    
    rho = sqrt(4.0 * WaveFunction(fhSI.x, abs(fhSI.f/(d*Mpc))**2 / S).Integrate(fLow, fHi, False).real)
    return rho


def SNRavg(fhSIModes, lmax, d, Sn, fLow, fHi):
    """Given a modes array in SI units in Fourier domain and distance d (in Mpc) this computes the angle averaged SNR for sensitivity curve Sn"""
    S = map(Sn, fhSIModes[2][2].x) # create table of values for noise curve
    
    rho = 0.0
    for ll in range(2, lmax+1):
        for mm in range(-ll, ll+1):
            rho += 1.0/pi * WaveFunction(fhSIModes[ll][mm].x, abs(fhSIModes[ll][mm].f/(d*Mpc))**2 / S).Integrate(fLow, fHi, False).real
    return sqrt(rho)


def IntCharFreq(fhSI, d, Sn, fLow, fHi):
    """Given a reconstructed h in SI units in Fourier domain and distance d (in Mpc) this computes the integrated characteristic signal frequency for sensitivity curve Sn"""
    """Remember to select either h_+ or h_\cross! Otherise the SNR will be twice as large!"""
    S = map(Sn, fhSI.x) # create table of values for noise curve
    
    fc = WaveFunction(fhSI.x, abs(fhSI.f/(d*Mpc))**2 * fhSI.x/ S).Integrate(fLow, fHi, False).real / WaveFunction(fhSI.x, abs(fhSI.f/(d*Mpc))**2 / S).Integrate(fLow, fHi, False).real
    return fc



def WienerScalarProduct(hSI1, hSI2, d, Sn, fLow, fHi):
    """Given a reconstructed h1 and h2 in SI units in Fourier domain and distance d (in Mpc) this computes the Wiener scalar product for sensitivity curve Sn"""
    S = map(Sn, hSI1.x) # create table of values for noise curve
    
    if (any(hSI1.x != hSI2.x) and len(hSI1.x) != len(hSI2.x)): print("Warning: Waveforms are not stored at same frequencies!")
    
    prod = 4.0 * WaveFunction(hSI1.x, hSI1.f/(d*Mpc)*hSI2.f.conjugate()/(d*Mpc) / S).Integrate(fLow, fHi, False).real
    
    return prod

def WienerDiscreteScalarProduct(hSI1, hSI2, d, Sn, fLow, fHi):
    """Given a reconstructed h1 and h2 in SI units in Fourier domain and distance d (in Mpc) this computes the discrete Wiener scalar product for sensitivity curve Sn"""
    S = map(Sn, hSI1.x[1:]) # create table of values for noise curve
    
    if (any(hSI1.x != hSI2.x) and len(hSI1.x) != len(hSI2.x)): print("Warning: Waveforms are not stored at same frequencies!")
    
    s = sum(WaveFunction(hSI1.x[1:], hSI1.f[1:]/(d*Mpc)*hSI2.f[1:].conjugate()/(d*Mpc) / S).f).real  #.Integrate(fLow, fHi, False).real
    prod = 4*(max(abs(hSI1.x)) - min(abs(hSI1.x))) * s/hSI1.Length()
    return prod


def OrthonormalBasis(hSI_plus, hSI_cross, Sn, fLow, fHi):
    """Given the two linearly independent waves h_plus and h_cross in SI units and in the Fourier domain, we find an orthonormalised basis wrt to the detector inner product"""
    # set scalar product
    ScalProd = WienerDiscreteScalarProduct   # Use discrete product for now. The continuous product yields the same answer to high precision
    
    norm_plus = sqrt(ScalProd(hSI_plus, hSI_plus, 1, Sn, fLow, fHi))
    norm_cross = sqrt(ScalProd(hSI_cross, hSI_cross, 1, Sn, fLow, fHi))
    
    # prevent NANs
    if norm_plus == 0: norm_plus = 1.0
    if norm_cross == 0: norm_cross = 1.0
    
    # First vector is simply the nomralised h_+ 
    e_plus = WaveFunction(hSI_plus.x, hSI_plus.f/norm_plus)
    # Normalise cross mode as well
    e_cross = WaveFunction(hSI_cross.x, hSI_cross.f/norm_cross)
    # Get orthonormal e_ortho that is orthonormal to e_plus
    e_ortho = WaveFunction(e_cross.x, e_cross.f - e_plus.f*ScalProd(e_plus, e_cross, 1, Sn, fLow, fHi)
              / sqrt(1.0 - ScalProd(e_plus, e_cross, 1, Sn, fLow, fHi)**2))
    
    #print("<e_plus, e_ortho> = {0}".fornat(ScalProd(e_plus, e_ortho, 1, Sn, fLow, fHi)))
    #print("<e_plus, e_plus> = {0}".format(ScalProd(e_plus, e_plus, 1, Sn, fLow, fHi)))
    #print("<e_cross, e_cross> = {0}".format(ScalProd(e_cross, e_cross, 1, Sn, fLow, fHi)))
    
    return e_plus, e_ortho

def CorrelationFunction(fhSI1, fhSI2, Sn, fLow, fHi):
    """Computes the correlation function over a given detector by using Fourier-modes in SI units.
       The definition of the correlation function is given in Phys. Rev. D 44, 3819-3834, Eq. 2.6, and especially Eq. 2.8 in terms of Fourier transforms."""
    
    S = map(Sn, abs(fhSI1.x)) # create table of values for noise curve
    S[0] = 1
    
    norm1 = max(WaveFunction(fhSI1.x, fhSI1.f * fhSI1.f.conjugate() / S).InverseFourierTransform(False, True).f.real)
    norm2 = max(WaveFunction(fhSI2.x, fhSI2.f * fhSI2.f.conjugate() / S).InverseFourierTransform(False, True).f.real)
    
    # Prevent NANs
    if norm1 == 0: norm1 = 1.0
    if norm2 == 0: norm2 = 1.0
    
    C = WaveFunction(fhSI1.x, fhSI1.f * fhSI2.f.conjugate() / S).InverseFourierTransform(False, True)
    C.f = C.f.real/sqrt(norm1*norm2)
    #C.f[0] = 0 # get rid of zero-freq component
    return C

def Match(hSI1, hSI2, Sn, fLow, fHi, WFtype='h'):
    """Computes the best match given two waveforms hSI1 and hSI2 in SI units in the time domain.
    The waveforms are assumed to be in the standard from h=h_+-ih_x.
    If WFtype is not set to 'h', but 'psi4' or 'news' we will perform the integration in the Fourier domain to get h"""
    
    if (any(hSI1.x != hSI2.x) and len(hSI1.x) != len(hSI2.x)): print("Warning: Inconsistent waveforms. Interpolate to same gridpoints!")
    
    hSI1plus = SelecthPlus(hSI1)
    hSI1cross = SelecthCross(hSI1)
    hSI2plus = SelecthPlus(hSI2)
    hSI2cross = SelecthCross(hSI2)
    
    fhSI1plus = hSI1plus.FourierTransform(True)
    fhSI1cross = hSI1cross.FourierTransform(True)
    fhSI2plus = hSI2plus.FourierTransform(True)
    fhSI2cross = hSI2cross.FourierTransform(True)
    
    # restrict to requested frequency range
    for ii in range(0, fhSI1plus.Length()):
        if fhSI1plus.x[ii] < fLow or fhSI1plus.x[ii] > fHi: 
            fhSI1plus.f[ii] = 0.0
            fhSI1cross.f[ii] = 0.0
            fhSI2plus.f[ii] = 0.0
            fhSI2cross.f[ii] = 0.0
    
    # user has not provided h but psi4 in SI units.
    if (WFtype == 'psi4'):
        fhSI1plus = GethFD(fhSI1plus, 'psi4')
        fhSI1cross = GethFD(fhSI1cross, 'psi4')
        fhSI2plus = GethFD(fhSI2plus, 'psi4')
        fhSI2cross = GethFD(fhSI2cross, 'psi4')
    # user has not provided h but news in SI units.
    if (WFtype == 'news'):
        fhSI1plus = GethFD(fhSI1plus, 'news')
        fhSI1cross = GethFD(fhSI1cross, 'news')
        fhSI2plus = GethFD(fhSI2plus, 'news')
        fhSI2cross = GethFD(fhSI2cross, 'news')
    
    e1_plus, e1_ortho = OrthonormalBasis(fhSI1plus, fhSI1cross, Sn, fLow, fHi)
    e2_plus, e2_ortho = OrthonormalBasis(fhSI2plus, fhSI2cross, Sn, fLow, fHi)
    
    # We compute the correlation function in order to maximize over time-lag between the two given waveforms.
    p1p2p = CorrelationFunction(e1_plus, e2_plus, Sn, fLow, fHi)
    p1p2n = CorrelationFunction(e1_plus, e2_ortho, Sn, fLow, fHi)
    
    p1n2p = CorrelationFunction(e1_ortho, e2_plus, Sn, fLow, fHi)
    p1n2n = CorrelationFunction(e1_ortho, e2_ortho, Sn, fLow, fHi)
    
    A = p1p2p.f**2 + p1p2n.f**2
    B = p1n2p.f**2 + p1n2n.f**2
    C = p1p2p.f*p1n2p.f + p1p2n.f*p1n2n.f
    
    #A = WienerScalarProduct(e1_plus, e2_plus, 1, Sn, fLow, fHi)**2 + WienerScalarProduct(e1_plus, e2_ortho, 1, Sn, fLow, fHi)**2
    #B = WienerScalarProduct(e1_ortho, e2_plus, 1, Sn, fLow, fHi)**2 + WienerScalarProduct(e1_ortho, e2_ortho, 1, Sn, fLow, fHi)**2
    #C = WienerScalarProduct(e1_plus, e2_plus, 1, Sn, fLow, fHi) * WienerScalarProduct(e1_ortho, e2_plus, 1, Sn, fLow, fHi) 
    #    + WienerScalarProduct(e1_plus, e2_ortho, 1, Sn, fLow, fHi) * WienerScalarProduct(e1_ortho, e2_ortho, 1, Sn, fLow, fHi)
    
    
    best = sqrt( (A+B)/2 + ( sqrt(((A-B)/2)**2 + C**2) ))
    minmax = sqrt( (A+B)/2 - ( sqrt(((A-B)/2)**2 + C**2) ))
    return max(best), max(minmax)  # pick max value according to correlation.


def MisMatch(hSI1, hSI2, Sn, fLow, fHi, WFtype='h'):
    """Same as above but 1-Match"""
    best, minmax = Match(hSI1, hSI2, Sn, fLow, fHi, WFtype='h')
    return 1-best, 1-minmax



def Richardson_error(order, dxLow, dxHi, errorDF):
    """Richardson extrapolates the errorDF assumed as given by the difference on low dxLow and hi resolution dxHi by the order specified."""
    err = DiscreteFunction(errorDF.x, 1.0/((dxLow/dxHi)**order-1.0) * errorDF.f)
    return err


def PhaseAlignOverRange(DF1, DF2, x0, x1, start_xoff, delta, dx, accuracy=1e-2, monotonic=False):
    """Phase-aligns two waveforms by minimizing their difference in the L2 norm and 
       returns the time-shift (and phase shift) that needs to be applied. The procedure is that
       from Boyle et al, arxiv:0804.4184v2.
       DF1: phase of waveform 1
       DF2: phase of waveform 2
       [x0,x1]: The minimization interval.
       start_xoff: An initial guess to speed up computation.
       delta: the shift will be applied in the range [start_xoff-delta, start_xoff+delta]. 
              The smaller delta, the faster (but keep in mind that you need to be sure that the minimum is located in that range!
       dx: The delta spacing of interpolated function values on which the pointwise difference is taken.
       accuracy: The accuracy in finding the shift that minimizes the L2 norm of the error ('the delta spacing of the shift').
       monotonic: Set to 'True' if function is monotonic (this will greatly speed up computation)"""

    # Note: this function is essentially identical to AlignOverRange in DiscreteFunction.py.
    #       But this also return a phase shift!

    # Keep DF1 fixed. Find offset x_off that needs to be applied to DF2.x
    # in order to make the L2-norm of the difference over a range [x0,x1]
    # as small as possible.

    x_off = start_xoff-delta
    minErr = 10e10
    minX = x_off
    trueMinX = x_off
    trueDeltaPhi = 0
    firstTime = True

    # If function is monotonic, we can use a faster algorithm
    if (monotonic==True):
        factor = 2.0*delta/accuracy*0.1
        maxX = start_xoff+delta
        while (factor >= 1.0 or firstTime == True):
            x_off = minX
            minErr = 10e10
            while (x_off < maxX):
                DF2shift = DiscreteFunction(DF2.x+x_off, DF2.f)
                subDF1, subDF2 = MakeConsistent(DF1, DF2shift, x0, x1, dx)
                # first, compute DeltaPhi analytically according to Boyle et al. via midpoint rule
                DeltaPhi = 0
                DeltaPhi += 0.5*(subDF1.f[0]-subDF2.f[0])
                DeltaPhi += 0.5*(subDF1.f[-1]-subDF2.f[-1])
                for ii in range (1, subDF1.Length()-1):
                    DeltaPhi += (subDF1.f[ii]-subDF2.f[ii])
                DeltaPhi /= (subDF1.Length()-1)
                # Now compute error
                err = 0
                for ii in range (0, subDF1.Length()):
                    err += (subDF1.f[ii]-subDF2.f[ii]-DeltaPhi)**2
                #err = sqrt(err)  # sqrt not necessary
                #print(minErr, err)
                if (err < minErr):
                    minErr = err
                    minX = x_off-accuracy*factor
                    trueMinX = x_off
                    trueDeltaPhi = DeltaPhi
                else:
                    maxX = x_off
                    break
                x_off += accuracy*factor
            factor *= 0.1
            firstTime = False
            #print(factor)
        return trueMinX, trueDeltaPhi
        
        print("ERROR! ...using brute force method!")
    
    # Proceed with brute force algorithm...

    x_off = start_xoff-delta

    minErr = 10e10
    minX = x_off
    minDeltaPhi = 0


    # We shift DF2 and compute the error over range [x0,x1]
    while (x_off < start_xoff+delta):
        DF2shift = DiscreteFunction(DF2.x+x_off, DF2.f)
        subDF1, subDF2 = MakeConsistent(DF1, DF2shift, x0, x1, dx)
        # first, compute DeltaPhi analytically according to Boyle et al.
        DeltaPhi = 0
        DeltaPhi += 0.5*(subDF1.f[0]-subDF2.f[0])
        DeltaPhi += 0.5*(subDF1.f[-1]-subDF2.f[-1])
        for ii in range (1, subDF1.Length()-1):
            DeltaPhi += (subDF1.f[ii]-subDF2.f[ii])
        DeltaPhi /= (subDF1.Length()-1)
        err = 0
        for ii in range (0, subDF1.Length()):
            err += (subDF1.f[ii]-subDF2.f[ii]-DeltaPhi)**2
        #err = sqrt(err) # sqrt not necessary
        #print(err)
        if (err < minErr):
            minErr = err
            minX = x_off
            minDeltaPhi = DeltaPhi
        x_off += accuracy

    return minX, minDeltaPhi
