from ..utils.log import logger
from ..utils.file import columns_asc, dataset_asc, header_h5, select_header_h5, dataset_h5, read
from ..Analysis import TimeSeries, GravitationalWave, VectorSeries
from ..Analysis import Trajectory, BlackHoleBinary, BHDiags
from ..Lib.power import RadialToTortoise, psi4ToStrain, Extrapolate, angular_momentum
from ..Lib.pygwanalysis import InitModeArray
from ..Lib.surrkick import coeffs
import numpy as np
import scipy.optimize
import configparser
import re
import os

class multipole:
    """
    A class representing the output of the Multipole thorn in a simulation。We store a dictionary {key: files} to be able to track the location of the data across the multiple files.
    """

    def __init__(self, files, ftype='h5'):
        """
        Output a dict for each mode at each radius.
        """
        self.files = files
        self.psi4 = dict()
        if ftype == 'h5':
            for file in files:
                with read(file) as f:
                    for dset in f:
                        mp = re.match(r'l(\d*)_m(-?\d*)_r(\d*\.\d)', dset)
                        if mp is not None:
                            l = int(mp.group(1))
                            m = int(mp.group(2))
                            radius = float(mp.group(3))
                            if (l, m) in self.psi4:
                                if radius in self.psi4[(l, m)]:
                                    self.psi4[(l, m)][radius] = np.append(self.psi4[(l, m)][radius], np.array(f[dset]).T, axis=1)
                                else:
                                    self.psi4[(l, m)][radius] = np.array(f[dset]).T
                            else:
                                self.psi4[(l, m)] = {radius: np.array(f[dset]).T}
        elif ftype == 'asc':
            pat = re.compile('^mp_([a-zA-Z0-9\[\]_]+)_l(\d+)_m([-]?\d+)_r([0-9.]+).asc$')
            for file in files:
                mp = pat.match(os.path.basename(file))
                if mp is not None:
                    var = mp.group(1).lower()
                    l = int(mp.group(2))
                    m = int(mp.group(3))
                    radius = float(mp.group(4))
                    if var == 'psi4':
                        if (l, m) in self.psi4:
                            if radius in self.psi4[(l, m)]:
                                self.psi4[(l, m)][radius] = np.append(self.psi4[(l, m)][radius], np.loadtxt(file, comments="#", unpack=True), axis=1)
                            else:
                                self.psi4[(l, m)][radius] = np.loadtxt(file, comments="#", unpack=True)
                        else:
                            self.psi4[(l, m)] = {radius: np.loadtxt(file, comments="#", unpack=True)}

        for mode in self.psi4:
            for radius in self.psi4[mode]:
                self.psi4[mode][radius] = np.unique(self.psi4[mode][radius], axis=1)

        self._M_ADM = None
        self._f0 = None
        self._radiu = None 
        self._phase_extrapolation_order = None
        self._amp_extrapolation_order = None
        self._t = None
        self._h = None
        self._hdot = None
        self._dEdt = None
        self._dPdt = None
        self._dJdt = None

    @property
    def lmax(self):
        return max([l for l, m in self.psi4.keys()])

    def rPsi4(self, M_ADM, mode, radiu):
        assert mode in self.psi4.keys()
        assert radiu in self.psi4[mode].keys(), "Radius: {}".format(self.psi4[mode].keys())
        psi4 = self.psi4[mode][radiu] 
        t = psi4[0] - RadialToTortoise(radiu, M_ADM)
        rpsi4 = (psi4[1]+1.j*psi4[2]) * radiu
        if np.absolute(rpsi4).min() <= np.finfo(float).eps and l >= 2:
            logger.warning("The psi4 amplitude of mode ({}, {}) at radiu {} is near zero. The phase is ill-defined.", l, m, radiu) 
        return TimeSeries(t, rpsi4)

    def Strain(self, M_ADM, f0, mode, radiu=-1, phase_extrapolation_order=1, amp_extrapolation_order=2):
        """
        This class is used to obtain GW signal multipole components from :math:`\Psi_4` Weyl scalar multipole components extracted at various distances.

        :param float f0: FFI cutoff frequency (ie :math:`\omega`). This must be choosen smaller than any physically expected frequency.
        :param int lmax: Maximum l-mode to process
        :param int distance: distance to process
        """
        if radiu != -1:
            assert radiu in self.psi4[mode].keys(), "Radius: {}".format(self.psi4[mode].keys())
            psi4 = self.psi4[mode][radiu] 
            t = psi4[0] - RadialToTortoise(radiu, M_ADM)
            rpsi4 = (psi4[1]+1.j*psi4[2]) * radiu
            if np.absolute(rpsi4).min() <= np.finfo(float).eps and mode[0]>= 2:
                logger.warning("The psi4 amplitude of mode {} at radiu {} is near zero. The phase is ill-defined.", mode, radiu) 
            #Fixed-frequency integration twice to get strain
            return GravitationalWave(*psi4ToStrain(t, rpsi4, f0))
        else:
            # Extrapolate
            radii = sorted(self.psi4[mode])
            strains = []
            for radiu in radii: 
                psi4 = self.psi4[mode][radiu] 
                t = psi4[0] - RadialToTortoise(radiu, M_ADM)
                rpsi4 = (psi4[1]+1.j*psi4[2]) * radiu
                if np.absolute(rpsi4).min() <= np.finfo(float).eps and mode[0]>= 2:
                    logger.warning("The psi4 amplitude of mode {} at radiu {} is near zero. The phase is ill-defined.", mode, radiu) 
                #Fixed-frequency integration twice to get strain
                strains.append(psi4ToStrain(t, rpsi4, f0))
            return GravitationalWave(*Extrapolate(strains, radii, phase_extrapolation_order, amp_extrapolation_order, 9))

    def _init_h(self, M_ADM, f0, radiu=-1, phase_extrapolation_order=1, amp_extrapolation_order=2):
        lmax = self.lmax
        self._h = InitModeArray(lmax)
        self._hdot = InitModeArray(lmax)
        self._M_ADM = M_ADM
        self._f0 = f0
        self._radiu = radiu 
        self._phase_extrapolation_order = phase_extrapolation_order
        self.amp_extrapolation_order = amp_extrapolation_order

        for l in np.arange(2, lmax+1):
            for m in np.arange(-l, l+1):
                h = self.Strain(M_ADM, f0, (l,m), radiu, phase_extrapolation_order, amp_extrapolation_order)
                self._t = h.t
                self._h[l][m] = h.y
                self._hdot[l][m] = h.gradient().y

    def h(self, l, m):
        """Correct the strain values to return zero if either l or m are not allowed (this is how the expressions of arXiv:0707.4654 are supposed to be used). Returns a single mode."""

        if l < 2 or l > self.lmax:
            return np.zeros(len(self._t), dtype=complex)
        elif m < -l or m > l:
            return np.zeros(len(self._t),dtype=complex)
        else:
            return self._h[l][m]

    def hdot(self, l, m):
        """Correct the hdot values to return zero if either l or m are not allowed (this is how the expressions of arXiv:0707.4654 are supposed to be used). Returns a single mode."""

        if l < 2 or l > self.lmax:
            return np.zeros(len(self._t),dtype=complex)
        elif m < -l or m > l:
            return np.zeros(len(self._t),dtype=complex)
        else:
            return self._hdot[l][m]

    def RadiatedEnergyFlux(self, lmax=None, **kwargs):
        """Radiated energy as a function of time. Time integral of Eq. (3.8) of arXiv:0707.4654, evaluated at the time nodes.

        :param int lmax: Maximum l-mode to process
        """
        if lmax is None:
            lmax = self.lmax
        assert lmax <= self.lmax

        initialize = False
        if kwargs:
            for i in kwargs:       
                if eval('self._'+i) != kwargs[i]:
                    initialize = True
        if ((self._h is None) or (self._hdot is None)):
            initialize = True

        if initialize:
            self._init_h(**kwargs)

        dEdt = 0
        for l in np.arange(2, lmax+1):
            for m in np.arange(-l, l+1):
                dEdt += (1/(16*np.pi)) * np.abs(self.hdot(l,m))**2
        self._dEdt = TimeSeries(self._t, dEdt)

        return self._dEdt

    def RadiatedEnergy(self, E0, **kwargs):
        self.RadiatedEnergyFlux(**kwargs)
        Eoft = self._dEdt.integrate()
        Eoft.y += E0
        return Eoft

    def RadiatedLinMomFlux(self, lmax=None, **kwargs):
        """Implement Eq. (3.14-3.15) of arXiv:0707.4654 for the three component of the linear momentum momentum flux. Note that the modes provided by the surrogate models are actually h*(r/M) extracted as r=infinity, so the r^2 factor is already in there. Returned array has size len(times)x3.
        """
        if lmax is None:
            lmax = self.lmax
        assert lmax <= self.lmax

        initialize = False
        if kwargs:
            for i in kwargs:       
                if eval('self._'+i) != kwargs[i]:
                    initialize = True
        if ((self._h is None) or (self._hdot is None)):
            initialize = True

        if initialize:
            self._init_h(**kwargs)

        dPpdt = 0
        dPzdt = 0
        for l in np.arange(2, lmax+1):
            for m in np.arange(-l, l+1):
                # Eq. 3.14. dPpdt= dPxdt + i dPydt
                dPpdt += (1/(8*np.pi)) * self.hdot(l,m) * ( coeffs.a(l,m) * np.conj(self.hdot(l,m+1)) + coeffs.b(l,-m) * np.conj(self.hdot(l-1,m+1)) - coeffs.b(l+1,m+1) * np.conj(self.hdot(l+1,m+1)) )
                # Eq. 3.15
                dPzdt += (1/(16*np.pi)) * self.hdot(l,m) * ( coeffs.c(l,m) * np.conj(self.hdot(l,m)) + coeffs.d(l,m) * np.conj(self.hdot(l-1,m)) + coeffs.d(l+1,m) * np.conj(self.hdot(l+1,m)) )
        dPxdt=dPpdt.real
        dPydt=dPpdt.imag
        assert max(dPzdt.imag)<1e-6 # Check...
        dPzdt=dPzdt.real

        self._dPdt = VectorSeries(self._t, np.stack([dPxdt, dPydt, dPzdt], axis=1))
        return self._dPdt

    def RadiatedLinMom(self, **kwargs):
        self.RadiatedLinMomFlux(**kwargs)
        Poft = self._dPdt.integrate()
        # Eliminate unphysical drift due to the starting point of the integration. Integrate for tbuffer and substract the mean.
        tbuffer=1000
        tstart = self._t[0]
        tend= tstart + tbuffer
        P0 = Poft.integrate().clip(tstart, tend).y[-1]/tbuffer
        Poft.y -= P0
        return Poft

    def RadiatedAngMomFlux(self, lmax=None, **kwargs):
        """
        Implement Eq. (3.22-3.24) of arXiv:0707.4654 for the three component of the angular momentum momentum flux. Note that the modes provided by the surrogate models are actually h*(r/M) extracted as r=infinity, so the r^2 factor is already in there.
        """
        if lmax is None:
            lmax = self.lmax
        assert lmax <= self.lmax

        initialize = False
        if kwargs:
            for i in kwargs:       
                if eval('self._'+i) != kwargs[i]:
                    initialize = True
        if ((self._h is None) or (self._hdot is None)):
            initialize = True

        if initialize:
            self._init_h(**kwargs)

        dJxdt = 0
        dJydt = 0
        dJzdt = 0
        for l in np.arange(2, lmax+1):
            for m in np.arange(-l, l+1):
                # Eq. 3.22
                dJxdt += (1/(32*np.pi)) * self.h(l,m) * ( coeffs.f(l,m) * np.conj(self.hdot(l,m+1)) + coeffs.f(l,-m) * np.conj(self.hdot(l,m-1)) )
                # Eq. 3.23
                dJydt += (-1/(32*np.pi)) * self.h(l,m) * ( coeffs.f(l,m) * np.conj(self.hdot(l,m+1)) - coeffs.f(l,-m) * np.conj(self.hdot(l,m-1)) )
                # Eq. 3.24
                dJzdt += (1/(16*np.pi)) * m * self.h(l,m) * np.conj(self.hdot(l,m))

        dJxdt=dJxdt.imag
        dJydt=dJydt.real
        dJzdt=dJzdt.imag

        self._dJdt = VectorSeries(self._t, np.stack([dJxdt, dJydt, dJzdt], axis=1))
        return self._dJdt

    def RadiatedAngMom(self, **kwargs):
        self.RadiatedAngMomFlux(**kwargs)
        Joft = self._dJdt.integrate()
        return Joft

class puncturetracker(Trajectory):
    
    def __init__(self, files):
        columns = columns_asc(files[0])
        index = []
        vars = [
            'pt_loc_t[0]', 
            'pt_loc_x[0]', 
            'pt_loc_y[0]', 
            'pt_loc_z[0]', 
            'pt_loc_x[1]', 
            'pt_loc_y[1]', 
            'pt_loc_z[1]',
        ]
        for var in vars:
            index.append(columns.index(var))
        super().__init__(*dataset_asc(files, index))


class quasilocalmeasures(BlackHoleBinary):

    def __init__(self, files):
        columns = columns_asc(files[0])
        index = []
        vars = [
            'qlm_time[0]',
            'qlm_area[0]',
            'qlm_area[1]',
            'qlm_area[2]',
            'qlm_mass[0]',
            'qlm_mass[1]',
            'qlm_mass[2]',
            'qlm_spin[0]',
            'qlm_spin[1]',
            'qlm_spin[2]',
            'qlm_coordspinx[0]', 
            'qlm_coordspiny[0]', 
            'qlm_coordspinz[0]', 
            'qlm_coordspinx[1]', 
            'qlm_coordspiny[1]',
            'qlm_coordspinz[1]',
            'qlm_coordspinx[2]', 
            'qlm_coordspiny[2]',
            'qlm_coordspinz[2]',
        ]
        for var in vars:
            index.append(columns.index(var))
        super().__init__(*dataset_asc(files, index))

class twopunctures:
    
    def __init__(self, files):
        self.bbh = configparser.ConfigParser()
        self.bbh.read(files[0])

    @property
    def m1(self):
        return round(float(self.bbh['metadata']['initial-bh-puncture-adm-mass1']), 10)

    @property
    def m2(self):
        return round(float(self.bbh['metadata']['initial-bh-puncture-adm-mass2']), 10)

    @property
    def ADMMass(self):
        return float(self.bbh['metadata']['initial-ADM-energy'])

    @property
    def mass_ratio(self):
        return self.m1 / self.m2

    @property
    def eta(self):
        return self.m1*self.m2/(self.m1+self.m2)**2

    @property
    def separation(self):
        return float(self.bbh['metadata']['initial-separation'])

    @property
    def spin1(self):
        spin1x = float(self.bbh['metadata']['initial-bh-spin1x'])
        spin1y = float(self.bbh['metadata']['initial-bh-spin1y'])
        spin1z = float(self.bbh['metadata']['initial-bh-spin1z'])
        return np.asarray([spin1x, spin1y, spin1z])

    @property
    def spin2(self):
        spin2x = float(self.bbh['metadata']['initial-bh-spin2x'])
        spin2y = float(self.bbh['metadata']['initial-bh-spin2y'])
        spin2z = float(self.bbh['metadata']['initial-bh-spin2z'])
        return np.asarray([spin2x, spin2y, spin2z])

    @property
    def chi1(self):    
        return self.spin1 / self.m1**2

    @property
    def chi2(self):    
        return self.spin2 / self.m2**2

    @property
    def chi_eff(self):    
        """
        The mass-weighted sum of spins projected along the direction perpendicular to the orbital plane.

        \chi_{\mathrm{eff}}=\frac{m_{1} \chi_{1}+m_{2} \chi_{2}}{m_{1}+m_{2}}
        """
        return (self.m1*self.chi1[2] + self.m2*self.chi2[2]) / (self.m1 + self.m2)

    @property
    def chi_p(self):    
        """
        The dimensionless precession spin parameter Eq.(3.4) arXiv:1408.1810
        """
        A1 = 2 + 3/(2*self.mass_ratio)
        A2 = 2 + 3*self.mass_ratio/2
        S1 = [self.spin1[0], self.spin1[1], 0]
        S2 = [self.spin2[0], self.spin2[1], 0]
        Sp = max(A1*np.linalg.norm(S1), A2*np.linalg.norm(S2))
        return Sp / (A1*self.m1**2)      

    @property
    def Omega_orb(self):
        x1 = float(self.bbh['metadata']['initial-bh-position1x'])
        x2 = float(self.bbh['metadata']['initial-bh-position2x'])
        py1 = float(self.bbh['metadata']['initial-bh-momentum1y'])
        py2 = float(self.bbh['metadata']['initial-bh-momentum2y'])
        LInitNR = x1*py1 + x2*py2
        # .014 is the initial guess for cutoff frequency
        omOrbPN = scipy.optimize.fsolve(angular_momentum, .014, (self.mass_ratio, self.ADMMass, self.chi1[2], self.chi2[2], LInitNR))[0]
        return omOrbPN**(3./2.)        

    @property
    def CutoffFrequency(self):
        omGWPN = 2. * self.Omega_orb
        return 0.75 * omGWPN

    def __str__(self):
        output  = "m1 = {}\n".format(self.m1)
        output += "m2 = {}\n".format(self.m2)
        output += "ADMMass = {}\n".format(self.ADMMass)
        output += "Mass Ratio = {}\n".format(self.mass_ratio)
        output += "eta = {}\n".format(self.eta)
        output += "separation = {}\n".format(self.separation)
        output += "spin1 = {}\n".format(self.spin1)
        output += "spin2 = {}\n".format(self.spin2)
        output += "chi1 = {}\n".format(self.chi1)
        output += "chi2 = {}\n".format(self.chi2)
        output += "chi_eff = {}\n".format(self.chi_eff)
        output += "chi_p = {}\n".format(self.chi_p)
        output += "Omega_orb = {}\n".format(self.Omega_orb)
        output += "CutoffFrequency = {}\n".format(self.CutoffFrequency)
        return output


class nanchecker:

    def __init__(self, files):
        self.vars = set()
        self.it = set()
        self.hierarchy = {}
        self.header = {}
        for file in files:
            self.header[file] = header_h5(file)
            for h, v in self.header[file].items():
                self.vars.add(v['varname'])
                self.it.add(v['iteration'])
                ml = v['ml']
                if ml not in self.hierarchy:
                    self.hierarchy[ml] = dict()
                rl = v['rl']
                if rl not in self.hierarchy[ml]:
                    self.hierarchy[ml][rl] = set()
                if 'c' in v:
                    self.hierarchy[ml][rl].add(v['c'])
    
    def dsets(self, it=-1, ml=-1, rl=-1, c=-1):
        """
        Select specified datasets.

        :param int it: -1: all iteration number.
        :param int ml: -1: all map.
        :param int rl: -1: all refinement level.
        :param int c: -1: all component.
        :return: a dict of datasets
        """
        if it != -1:
            assert it in self.it, "iteration {} is not exist".format(it)
        if ml != -1:
            assert ml in self.hierarchy, "map {} is not exist".format(ml)
        if rl != -1:
            assert rl in self.hierarchy[0], "refinement level {} is not exist at map 0".format(rl)
        if c != -1:
            assert rl != -1
            assert c in self.hierarchy[0][rl], "component {} is not exist in refinement level {}".format(c, rl)

        header = select_header_h5(self.header, -1, it, ml, rl, c)
        return dataset_h5(header)

class ahfinderdirect(BHDiags):

    def __init__(self, files):
        dsets = {}
        pat = re.compile(r'BH_diagnostics.(ah\d+).gp')
        for file in files:
            mp = pat.search(os.path.basename(file))
            if (mp is not None):
                ahidx = mp.group(1)
                dsets.setdefault(ahidx,[]).append(file)
        idx = [1, 2, 3, 4, 7, 20, 21, 22, 26]
        for k, v in dsets.items():
            dsets[k] = dataset_asc(v, idx)
            
        super().__init__(dsets)

class carpet:

    def __init__(self, files):
        self.files = files

    def coordinates(self, it=0):
        # '^[\t ]*([^\s=:"\'#!\]\[]+)::([^\s=:"\'#!\]\[]+)[\t ]*=[\t ]*([^\s=:"\'#!]+)[\t ]*(?:!|#|\n|\r\n)'
        p = {}
        pat = re.compile('(\d) (\d) (\d+) (\d)[\t ]*\(([\-\d\:\/\,\.\]\[]+)\)')
        match = False
        for file in self.files:
            with read(file) as f:
                for line in f.readlines():
                    if "iteration" in line:
                        if int(line[10:]) == it:
                            match = True
                        else:
                            match = False
                    if match:
                        m = pat.match(line)
                        if m is not None:
                            map = int(m.group(1))
                            if map not in p:
                                p[map] = {}
                            rl = int(m.group(3))
                            if rl not in p[map]:
                                p[map][rl] = {}
                            c = int(m.group(4))
                            v = m.group(5).split('/')[0]
                            p[map][rl][c] = v

        return p

    @property
    def dx(self):
        map0 = self.coordinates()[0]
        p = {}
        for rl in map0:
            p[rl] = eval(map0[rl][0].split(':')[-1]) 
        return p