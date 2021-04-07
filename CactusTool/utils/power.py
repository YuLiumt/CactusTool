# Copyright (c) 2017 The Board of Trustees of the University of Illinois
# All rights reserved.
#
# Developed by: Daniel Johnson, E. A. Huerta, Roland Haas
#               NCSA Gravity Group
#               National Center for Supercomputing Applications
#               University of Illinois at Urbana-Champaign
#               http://gravity.ncsa.illinois.edu/

import math
import numpy as np
import scipy

def RadialToTortoise(r, M):
    """
    Convert the radial coordinate to the tortoise coordinate

    :param float r: radial coordinate
    :param float M: ADMMass used to convert coordinate
    :return: tortoise coordinate value
    """
    return r + 2. * M * math.log( r / (2. * M) - 1.)

def psi4ToStrain(mp_psi4, f0):
    """
    Convert the input mp_psi4 data to the strain of the gravitational wave

    :param array mp_psi4: Weyl scalar result from simulation
    :param float f0: cutoff frequency
    :return: = strain (h) of the gravitational wave
    """
    t0 = mp_psi4[:, 0]
    list_len = len(t0)
    complexPsi = np.zeros(list_len, dtype=np.complex_)
    complexPsi = mp_psi4[:, 1]+1.j*mp_psi4[:, 2]

    freq, psif = myFourierTransform(t0, complexPsi)
    dhf = ffi(freq, psif, f0)
    hf = ffi(freq, dhf, f0)

    return myFourierTransformInverse(freq, hf, t0[0])

def ffi(freq, data, f0):
    """
    Fixed frequency integration. Integrates the data according to the input frequency and cutoff frequency. See https://arxiv.org/abs/1508.07250 for method

    freq = fourier transform frequency
    data = input on which ffi is performed
    f0 = cutoff frequency
    """
    f1 = f0/(2*math.pi)
    fs = freq
    gs = data
    mask1 = (np.sign((fs/f1) - 1) + 1)/2.
    mask2 = (np.sign((-fs/f1) - 1) + 1)/2.
    mask = 1 - (1 - mask1) * (1 - mask2)
    fs2 = mask * fs + (1-mask) * f1 * np.sign(fs - np.finfo(float).eps)
    new_gs = gs/(2*math.pi*1.j*fs2)
    return new_gs

def myFourierTransform(t0, complexPsi):
    """
    Transforms the complexPsi data to frequency space

    t0 = time data points
    complexPsi = data points of Psi to be transformed
    """
    psif = np.fft.fft(complexPsi, norm="ortho")
    l = len(complexPsi)
    n = int(math.floor(l/2.))
    newpsif = psif[l-n:]
    newpsif = np.append(newpsif, psif[:l-n])
    T = np.amin(np.diff(t0))*l
    freq = range(-n, l-n)/T
    return freq, newpsif

#Inverse Fourier Transform
def myFourierTransformInverse(freq, hf, t0):
    l = len(hf)
    n = int(math.floor(l/2.))
    newhf = hf[n:]
    newhf = np.append(newhf, hf[:n])
    amp = np.fft.ifft(newhf, norm="ortho")
    df = np.amin(np.diff(freq))
    time = t0 + range(0, l)/(df*l)
    return time, amp

def InitModeArray(lmax):
    """Creates an array with empty modes up to lmax"""
    res = []
    for ll in range(0, lmax+1):
        res.append([])
        for mm in range(-ll, ll+1):
            res[-1].append([])
    return res

def Extrapolate(psi4, distances, f0, ADMMass):
    """
    Extrapolate to infinite.
    """
    phase = []
    amp = []
    t_start = - np.inf
    t_end = np.inf
    for dist in sorted(distances):
        dset = np.loadtxt(psi4[dist], comments="#")
        dset[:, 0] -= RadialToTortoise(dist, ADMMass)
        dset[:, 1] *= dist
        dset[:, 2] *= dist      
        dt = dset[1, 0] - dset[0, 0]
        # Fixed-frequency integration twice to get strain
        time, h = psi4ToStrain(dset, f0)
        # Get phase and amplitude of strain
        phase.append(scipy.interpolate.interp1d(time, np.unwrap(np.angle(h)), kind=9))
        amp.append(scipy.interpolate.interp1d(time, np.absolute(h), kind=9))
        if time[0] > t_start:
            t_start = time[0]
        if time[-1] < t_end:
            t_end = time[-1]
    t = np.arange(t_start, t_end, dt)
    for i in range(len(distances)):
        phase_resampled = phase[i](t)
        amp_resampled = amp[i](t)
    return phase_resampled