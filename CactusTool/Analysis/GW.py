from ..funcs import add_three_key, logger
from ..Lib.pygwanalysis import InitModeArray, WaveFunction, GethFromPsi4, SelecthPlus, SelecthCross, StrainToSI, Spectrogram
# , ReconFromsYlm, FStrainToSI
import numpy as np
import re
import os

class multipole:
    """
    A class representing the output of the Multipole thorn in a simulation。We store a dictionary {key: files} to be able to track the location of the data across the multiple files.
    """

    def __init__(self, files):
        self.vars = {}
        pat = re.compile('^mp_([a-zA-Z0-9\[\]_]+)_l(\d+)_m([-]?\d+)_r([0-9.]+).asc$')
        for file in files:
            mp = pat.match(os.path.basename(file))
            if mp is not None:
                var = mp.group(1).lower()
                l = int(mp.group(2))
                m = int(mp.group(3))
                distance = float(mp.group(4))
                add_three_key(self.vars, var, (l, m), distance, file)

    def Psi4ToStrain(self, f0, lmax=None, distance=None):
        """
        This class is used to obtain GW signal multipole components from :math:`\Psi_4` Weyl scalar multipole components extracted at various distances.

        :param float f0: FFI cutoff frequency (ie :math:`\omega/2\pi`). This must be choosen smaller than any physically expected frequency.
        :param int lmax: Maximum l-mode to process
        :param int distance: distance to process
        """
        assert 'psi4' in self.vars, "No psi4 data!"
        psi4 = self.vars['psi4']
        # Maximum l-mode to process
        if lmax is None:
            lmax = max([mode[0] for mode in psi4.keys()])
        # Maximum distance to process
        if distance is None:
            distance = max([dist for dist in psi4[(lmax, 0)].keys()])
        else:
            assert distance in psi4[(lmax, 0)], "No distance {}".format(distance)
        # Initialize a mode array
        Psi4 = InitModeArray(lmax)
        Strain = InitModeArray(lmax)
        # Load each psi4-mode into a WaveFunction object and store it in mode array
        logger.info("Load (l, m) psi4-mode into a WaveFunction.")
        for l in range(2,lmax+1):
            for m in range(-l, l+1):
                Psi4[l][m] = WaveFunction([], []) 
                Psi4[l][m].Load(psi4[(l, m)][distance])
        # Integrate (l,m) mode using FFI
        logger.info("Integrate (l, m) mode using FFI.")
        for l in range(2,lmax+1):
            for m in range(-l, l+1):
                WF = GethFromPsi4(Psi4[l][m], f0)
                Strain[l][m] = GravitationalWave(WF)
        logger.info("Get time of merger from maximum of amplitude and shift WaveFunction according to tmerger")
        return Strain


class GravitationalWave:

    def __init__(self, WF):
        self.WF = WF
    
    # def hplus(self):
    #     return SelecthPlus(self.WF)

    # def hcross(self):
    #     return SelecthCross(self.WF)

    # def Amplitude(self):
    #     return self.WF.Amplitude()

    def Preview(self, M):
        import matplotlib.pyplot as plt

        Mpc = 3.08568025*10**16 * 10**6 # 1 Mpc in[m] (~3.26 *10^6 ly)
        fig, ax = plt.subplots()
        logger.info("Converts time domain strain to SI units with {} M_sun", M)
        WF = StrainToSI(self.WF, M)
        hplus = SelecthPlus(WF)
        hcross = SelecthCross(WF)
        Amplitude = WF.Amplitude()
        ax.plot(hplus.x, hplus.f/(100*Mpc), label=r'$h^+$')
        ax.plot(hcross.x, hcross.f/(100*Mpc), label=r'$h^{\times}$')
        ax.plot(Amplitude.x, Amplitude.f/(100*Mpc), '--')
        logger.info("the amplitude at a hypothetical distance of D_obs = 100 Mpc")
        ax.legend()
        plt.xlabel(r'$t$')
        plt.ylabel(r'$h$')
        plt.show()

    # def Spectrogram(self, tau t0, t1, n):
    #     def hamming(t, sigma):
    #         return 0.54 - 0.46*cos(2.0*pi*n/(M-1))

    #     return Spectrogram(self.WF, kernel, tau, sigma, t0, t1, n)
