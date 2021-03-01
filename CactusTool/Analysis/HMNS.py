from ..funcs.array import list_bounds
import matplotlib.pyplot as plt
import numpy as np
import os

class hmns:

    def __init__(self, files):
        self.files = files
        self.dsets = {}
    
    @property
    def it(self):
        if 'Omega' not in self.dsets:
            self.dsets['Omega'] = self.load('HMNS_Omega.asc')
        it = self.dsets['Omega'].T[0].astype(np.int)
        return sorted(list(set(it)))

    def load(self, fname):
        file = [file for file in self.files if os.path.basename(file) == fname]
        assert len(file) == 1, "Too much file"
        return np.loadtxt(file[0], comments="#")

    def Omega(self, it=0, ax=None):
        if 'Omega' not in self.dsets:
            self.dsets['Omega'] = self.load('HMNS_Omega.asc')

        if ax is None:
            fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))

        dset = self.dsets['Omega'][self.dsets['Omega'][:, 0] == it].T
        time = dset[1][0]
        r = dset[3]
        theta = dset[4]
        data = dset[5]
        nr = len(set(r))
        ntheta = len(set(theta))
        Theta = theta.reshape(ntheta, nr)[:, 0]
        R = r.reshape(ntheta, nr)[0]
        omega = data.reshape(ntheta, nr)
        ax.pcolormesh(list_bounds(Theta), list_bounds(R), omega.T)
        ax.set_title('Time: {}'.format(time), fontsize=12)

    def RotationProfile(self, it=0, r_c=True, ax=None):
        if 'Omega_r' not in self.dsets:
            self.dsets['Omega_r'] = self.load('HMNS_Omega_r.asc')
        if r_c:
            if 'r_c' not in self.dsets:
                self.dsets['r_c'] = self.load('HMNS_r_c.asc')

        dset = self.dsets['Omega_r'][self.dsets['Omega_r'][:, 0] == it].T
        time = dset[1][0]
        r = dset[3]
        data = dset[4]

        if ax is None:
            fig, ax = plt.subplots()

        if r_c:
            dset_r = self.dsets['r_c'][self.dsets['r_c'][:, 0] == it].T
            R_c = dset_r[4]
            ax.plot(R_c, data)
        else:
            ax.plot(r, data)

        ax.set_title('Time: {}'.format(time), fontsize=12) 

    def RotationProfile_2D(self, ax=None):
        if ax is None:
            fig, ax = plt.subplots()

        if 'Omega_r' not in self.dsets:
            self.dsets['Omega_r'] = self.load('HMNS_Omega_r.asc')
        
        dset = self.dsets['Omega_r'].T
        time = dset[1]
        r = dset[3]
        data = dset[4]
        nt = len(set(time))
        nr = len(set(r))
        Time = time.reshape(nt, nr)[:, 0]
        R = r.reshape(nt, nr)[0]
        omega = data.reshape(nt, nr)
        im = ax.pcolormesh(list_bounds(R), list_bounds(Time), omega)
        if fig:
            fig.colorbar(im, ax=ax)
        else:
            return im        
        
