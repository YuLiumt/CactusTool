# Copyright (c) 2017 The Board of Trustees of the University of Illinois
# All rights reserved.
#
# Developed by: Daniel Johnson, E. A. Huerta, Roland Haas
#               NCSA Gravity Group
#               National Center for Supercomputing Applications
#               University of Illinois at Urbana-Champaign
#               http://gravity.ncsa.illinois.edu/

import numpy as np
import glob
import math
import scipy.interpolate

#-----Function Definitions-----#

#Convert radial to tortoise coordinates
def RadialToTortoise(r, M):
        """
        Convert the radial coordinate to the tortoise coordinate

        r = radial coordinate
        M = ADMMass used to convert coordinate
        return = tortoise coordinate value
        """
        return r + 2. * M * math.log( r / (2. * M) - 1.)

#Convert modified psi4 to strain
def psi4ToStrain(t0, rpsi4, f0):
        """
        Convert the input mp_psi4 data to the strain of the gravitational wave

        f0 = cutoff frequency
        return = strain (h) of the gravitational wave
        """
        freq, psif = myFourierTransform(t0, rpsi4)
        dhf = ffi(freq, psif, f0)
        hf = ffi(freq, dhf, f0)

        return myFourierTransformInverse(freq, hf, t0[0])

#Fixed frequency integration
# See https://arxiv.org/abs/1508.07250 for method
def ffi(freq, data, f0):
        """
        Integrates the data according to the input frequency and cutoff frequency

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

#Fourier Transform
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

def angular_momentum(x, q, m, chi1, chi2, LInitNR):
        eta = q/(1.+q)**2
        m1 = (1.+math.sqrt(1.-4.*eta))/2.
        m2 = m - m1
        S1 = m1**2. * chi1
        S2 = m2**2. * chi2
        Sl = S1+S2
        Sigmal = S2/m2 - S1/m1
        DeltaM = m1 - m2
        mu = eta
        nu = eta
        GammaE = 0.5772156649;
        e4 = -(123671./5760.)+(9037.* math.pi**2.)/1536.+(896.*GammaE)/15.+(-(498449./3456.)+(3157.*math.pi**2.)/576.)*nu+(301. * nu**2.)/1728.+(77.*nu**3.)/31104.+(1792. *math.log(2.))/15.
        e5 = -55.13
        j4 = -(5./7.)*e4+64./35.
        j5 = -(2./3.)*e5-4988./945.-656./135. * eta;
        a1 = -2.18522;
        a2 = 1.05185;
        a3 = -2.43395;
        a4 = 0.400665;
        a5 = -5.9991;
        CapitalDelta = (1.-4.*eta)**0.5

        l = (eta/x**(1./2.)*(
                1. +
                x*(3./2. + 1./6.*eta) +
        x**2. *(27./8. - 19./8.*eta + 1./24.*eta**2.) + 
        x**3. *(135./16. + (-6889./144. + 41./24. * math.pi**2.)*eta + 31./24.*eta**2. + 7./1296.*eta**3.) + 
        x**4. *((2835./128.) + eta*j4 - (64.*eta*math.log(x)/3.))+ 
        x**5. *((15309./256.) + eta*j5 + ((9976./105.) + (1312.*eta/15.))*eta*math.log(x))+
        x**(3./2.)*(-(35./6.)*Sl - 5./2.*DeltaM* Sigmal) + 
        x**(5./2.)*((-(77./8.) + 427./72.*eta)*Sl + DeltaM* (-(21./8.) + 35./12.*eta)*Sigmal) + 
        x**(7./2.)*((-(405./16.) + 1101./16.*eta - 29./16.*eta**2.)*Sl + DeltaM*(-(81./16.) + 117./4.*eta - 15./16.*eta**2.)*Sigmal) + 
        (1./2. + (m1 - m2)/2. - eta)* chi1**2. * x**2. +
        (1./2. + (m2 - m1)/2. - eta)* chi2**2. * x**2. + 
        2.*eta*chi1*chi2*x**2. +
        ((13.*chi1**2.)/9. +
        (13.*CapitalDelta*chi1**2.)/9. -
        (55.*nu*chi1**2.)/9. - 
        29./9.*CapitalDelta*nu*chi1**2. + 
        (14.*nu**2. *chi1**2.)/9. +
        (7.*nu*chi1*chi2)/3. +
        17./18.* nu**2. * chi1 * chi2 + 
        (13.* chi2**2.)/9. -
        (13.*CapitalDelta*chi2**2.)/9. -
        (55.*nu*chi2**2.)/9. +
        29./9.*CapitalDelta*nu*chi2**2. +
        (14.*nu**2. * chi2**2.)/9.)
        * x**3.))
        return l - LInitNR

# Convert Cacus truth values to python's
def to_bool(s):
        s = s.lower()
        if s in ["true", "yes", "1"]:
            return True
        elif s in ["false", "no", "0"]:
            return False
        else:
            raise(ValueError("Not a boolean values: %s" % s))

#Get Energy
def get_energy(sim):
        """
        Save the energy radiated energy
        sim = string of simulation
        """
        python_strain = np.loadtxt("./Extrapolated_Strain/"+sim+"/"+sim+"_radially_extrapolated_strain_l2_m2.dat")
        val = np.zeros(len(python_strain))
        val = val.astype(np.complex_)
        cur_max_time = python_strain[0][0]
        cur_max_amp = abs(pow(python_strain[0][1], 2))
        # TODO: rewrite as array operations (use numpy.argmax)
        for i in python_strain[:]:
                cur_time = i[0]
                cur_amp = abs(pow(i[1], 2))
                if(cur_amp>cur_max_amp):
                        cur_max_amp = cur_amp
                        cur_max_time = cur_time

        max_idx = 0
        for i in range(0, len(python_strain[:])):
                if(python_strain[i][1] > python_strain[max_idx][1]):
                        max_idx = i

        paths = glob.glob("./Extrapolated_Strain/"+sim+"/"+sim+"_radially_extrapolated_strain_l[2-4]_m*.dat")
        for path in paths:
                python_strain = np.loadtxt(path)

                t = python_strain[:, 0]
                t = t.astype(np.complex_)
                h = python_strain[:, 1] + 1j * python_strain[:, 2]
                dh = np.zeros(len(t), dtype=np.complex_)
                for i in range(0, len(t)-1):
                        dh[i] = ((h[i+1] - h[i])/(t[i+1] - t[i]))
                dh[len(t)-1] = dh[len(t)-2]

                dh_conj = np.conj(dh)
                prod = np.multiply(dh, dh_conj)
                local_val = np.zeros(len(t))
                local_val = local_val.astype(np.complex_)
                # TODO: rewrite as array notation using numpy.cumtrapz
                for i in range(0, len(t)):
                        local_val[i] = np.trapz(prod[:i], x=(t[:i]))
                val += local_val

        val *= 1/(16 * math.pi)
        np.savetxt("./Extrapolated_Strain/"+sim+"/"+sim+"_radially_extrapolated_energy.dat", val)

#Get angular momentum
def get_angular_momentum(python_strain):
        """
        Save the energy radiated angular momentum
        sim = string of simulation
        """
        python_strain = np.loadtxt("./Extrapolated_Strain/"+sim+"/"+sim+"_radially_extrapolated_strain_l2_m2.dat")
        val = np.zeros(len(python_strain))
        val = val.astype(np.complex_)
        cur_max_time = python_strain[0][0]
        cur_max_amp = abs(pow(python_strain[0][1], 2))
        # TODO: rewrite as array operations (use numpy.argmax)
        for i in python_strain[:]:
                cur_time = i[0]
                cur_amp = abs(pow(i[1], 2))
                if(cur_amp>cur_max_amp):
                        cur_max_amp = cur_amp
                        cur_max_time = cur_time

        max_idx = 0
        for i in range(0, len(python_strain[:])):
                if(python_strain[i][1] > python_strain[max_idx][1]):
                        max_idx = i

        paths = glob.glob("./Extrapolated_Strain/"+sim+"/"+sim+"_radially_extrapolated_strain_l[2-4]_m*.dat")
        for path in paths:
                python_strain = np.loadtxt(path)

                t = python_strain[:, 0]
                t = t.astype(np.complex_)
                h = python_strain[:, 1] + 1j * python_strain[:, 2]
                dh = np.zeros(len(t), dtype=np.complex_)
                # TODO: rewrite using array notation
                for i in range(0, len(t)-1):
                        dh[i] = ((h[i+1] - h[i])/(t[i+1] - t[i]))
                dh[len(t)-1] = dh[len(t)-2]

                dh_conj = np.conj(dh)
                prod = np.multiply(h, dh_conj)
                local_val = np.zeros(len(t))
                local_val = local_val.astype(np.complex_)
                # TODO: rewrite as array notation using numpy.cumtrapz. Move atoi call out of inner loop.
                for i in range(0, len(t)):
                        local_val[i] = np.trapz(prod[:i], x=(t[:i])) * int(((path.split("_")[-1]).split("m")[-1]).split(".")[0])
                val += local_val

        val *= 1/(16 * math.pi)
        np.savetxt("./Extrapolated_Strain/"+sim+"/"+sim+"_radially_extrapolated_angular_momentum.dat", val)

# Some adjustments to the source code
def Extrapolate(strain, radii, phase_extrapolation_order, amp_extrapolation_order, interpolation_order):
    tstart = strain[0][0][0]
    dt = strain[0][0][1] - strain[0][0][0]
    tend = strain[-1][0][-1]
    t = np.arange(tstart, tend, dt)
    phase_interp = []
    amp_interp = []
    for i in range(len(radii)):
        time = strain[i][0]
        h = strain[i][1]
        h_phase = np.unwrap(np.angle(h))
        interp_function = scipy.interpolate.interp1d(time, h_phase, kind=interpolation_order)
        resampled_phase_vals = interp_function(t)
        if i > 0:
            phase_shift = round((resampled_phase_vals[0] - phase_interp[0][0])/(2.*math.pi))*2.*math.pi
            resampled_phase_vals -= phase_shift
        phase_interp.append(resampled_phase_vals)
        h_amp = np.absolute(h)
        interp_function = scipy.interpolate.interp1d(time, h_amp, kind=interpolation_order)
        amp_interp.append(interp_function(t))    
    
    A_phase = np.power(radii, 0)
    for i in range(1, phase_extrapolation_order+1):
        A_phase = np.column_stack((A_phase, np.power(radii, -1*i)))
    A_amp = np.power(radii, 0)
    for i in range(1, amp_extrapolation_order+1):
        A_amp = np.column_stack((A_amp, np.power(radii, -1*i)))
    
    extrapolated_phase = np.empty(0)
    extrapolated_amp = np.empty(0)
    for i in range(0, len(t)):
        b_phase = np.empty(0)
        for j in range(len(radii)):
            b_phase = np.append(b_phase, phase_interp[j][i])   
        extrapolated_phase = np.append(extrapolated_phase, np.linalg.lstsq(A_phase, b_phase, rcond=None)[0][0]) 
        b_amp = np.empty(0)
        for j in range(len(radii)):        
            b_amp = np.append(b_amp, amp_interp[j][i])
        extrapolated_amp = np.append(extrapolated_amp, np.linalg.lstsq(A_amp, b_amp, rcond=None)[0][0])
    
    return t, extrapolated_amp*np.exp(1j*extrapolated_phase)