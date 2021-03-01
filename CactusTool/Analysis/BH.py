from ..Visualize.plot3d import plot_trajectory
from ..Carpet.Scalar import Scalar
from ..Lib.power import angular_momentum
import numpy as np
import scipy.optimize
# class QuasiLocalMeasures:

class puncturetracker:
    
    def __init__(self, files):
        p = Scalar(files)
        self.t  = p.dsets('pt_loc_t[0]')[0]
        self.x1 = p.dsets('pt_loc_x[0]')[1]
        self.x2 = p.dsets('pt_loc_x[1]')[1]
        self.y1 = p.dsets('pt_loc_y[0]')[1]
        self.y2 = p.dsets('pt_loc_y[1]')[1]
        self.z1 = p.dsets('pt_loc_z[0]')[1]
        self.z2 = p.dsets('pt_loc_z[1]')[1]

    def Preview(self, tstart=0, tend=-1, view_angle=None, axlim=None, grid=True):
        """
        Plot BH tracks in initial simulation frame
        """
        import matplotlib.pyplot as plt

        fig = plt.figure(figsize=(10,10))
        ax = fig.gca(projection='3d')

        plot_trajectory(ax, self.x1[tstart:tend]-self.x2[tstart:tend], self.y1[tstart:tend]-self.y2[tstart:tend], self.z1[tstart:tend]-self.z2[tstart:tend], color='k')

        ax.set_title("From {:.2f} to {:.2f} [M]".format(self.t[tstart], self.t[tend]))

        if view_angle:
            ax.view_init(view_angle[0],view_angle[1])

        if grid:
            ax.set_xlabel('(x1-x2)/M')
            ax.set_ylabel('(y1-y2)/M')
            ax.set_zlabel('(z1-z2)/M')
        else:
            ax.axis('off')
        if axlim is not None:
            ax.set_xlim(axlim)
            ax.set_ylim(axlim)
        plt.show()

    @property
    def distance(self):
        """
        Evolution of the coordinate distance between the star centers.
        """
        distance = np.sqrt((self.x1-self.x2)**2 + (self.y1-self.y2)**2 + (self.z1-self.z2)**2)
        return np.vstack((self.t, distance)) 

    def eccentricity(self, p0=[0, 0, 0.01, 0.002, 0], show=True):
        """
        $$
        \dot{D}(t)=A_{0}+A_{1} t-e D_{0} \omega_{e} \sin \left(\omega_{e} t+\phi_{e}\right)
        $$
        where $e$ is the eccentricity and $D_{0}$ the initial coordinate interbinary distance.

        arXiv:1605.03424
        """
        from scipy.optimize import curve_fit

        dist = self.distance
        # dist[0] = dist[0]
        # dist[1] = dist[1]

        def orbital_evolution(t, A0, A1, e, we, phie):
            return A0 + A1*t - e*dist[1][0]*we*np.sin(we*t + phie)

        try:
            tmerger = dist[0][np.where(dist[1] == 0.)][0]
        except:
            print("Make sure BNS have merged")

        # The fit is performed in the time interval between $t_{\mathrm{ret}}=50 \mathrm{M}$ and $t_{\mathrm{ret}}=\frac{2}{3} t_{\mathrm{merger}}$ to avoid the initial spurious radiation and the plunge phase but having at least one eccentricity cycle included.
        index = np.where((dist[0] >= 50) & (dist[0] <= 2/3 * tmerger))[0]
        t = dist[0][index]
        y = np.gradient(dist[1], dist[0])[index]

        params, params_covariance = curve_fit(orbital_evolution, t, y, p0=p0)
        print('Orbital Eccentricity =', params[2])

        if show:
            import matplotlib.pyplot as plt

            fig, ax = plt.subplots()
            ax.plot(t, y)
            ax.set_xlabel('time [M]')
            ax.set_ylabel(r'$\dot{D}$')
            ax.plot(t, orbital_evolution(t, params[0], params[1], params[2], params[3], params[4]), label='Fitted function')
            plt.legend(loc='best')
            ax.set_xlim(t.min(), t.max())
            plt.show()

class twopunctures:
    
    def __init__(self, files):
        from ..Lib.power import angular_momentum

        self.file = files[0]
        self.dsets = dict()
        with open(self.file) as f:
            for line in f.readlines():
                if '= ' in line:
                    key, value = line.rstrip("\n").partition("=")[::2]
                    self.dsets[key.strip()] = value.strip()
    
    @property
    def ADMMass(self):
        if "initial-ADM-energy" in self.dsets: 
            return float(self.dsets["initial-ADM-energy"])
    
    @property
    def ADMAngularMomentum(self):
        if "initial-ADM-angular-momentumz" in self.dsets: 
            return [float(self.dsets["initial-ADM-angular-momentumx"]), float(self.dsets["initial-ADM-angular-momentumy"]), float(self.dsets["initial-ADM-angular-momentumz"])]

    @property
    def Mp(self):
        if "initial-bh-puncture-adm-mass1" in self.dsets:  
            return float(self.dsets["initial-bh-puncture-adm-mass1"])

    @property
    def Mm(self):
        if "initial-bh-puncture-adm-mass2" in self.dsets:  
            return float(self.dsets["initial-bh-puncture-adm-mass2"])

    @property
    def xp(self):
        if "initial-bh-position1x" in self.dsets:  
            return float(self.dsets["initial-bh-position1x"])      

    @property
    def xm(self):
        if "initial-bh-position2x" in self.dsets:  
            return float(self.dsets["initial-bh-position2x"])   
    
    @property
    def Spinp(self):    
        if "initial-bh-spin1z" in self.dsets:   
            return [float(self.dsets["initial-bh-spin1x"]), float(self.dsets["initial-bh-spin1y"]), float(self.dsets["initial-bh-spin1z"])]

    @property
    def Spinm(self):    
        if "initial-bh-spin2z" in self.dsets:   
            return [float(self.dsets["initial-bh-spin2x"]), float(self.dsets["initial-bh-spin2y"]), float(self.dsets["initial-bh-spin2z"])]

    @property
    def Pp(self):    
        if "initial-bh-momentum1x" in self.dsets:   
            return [float(self.dsets["initial-bh-momentum1x"]), float(self.dsets["initial-bh-momentum1y"]), float(self.dsets["initial-bh-momentum1z"])]

    @property
    def Pm(self):    
        if "initial-bh-momentum2x" in self.dsets:   
            return [float(self.dsets["initial-bh-momentum2x"]), float(self.dsets["initial-bh-momentum2y"]), float(self.dsets["initial-bh-momentum2z"])]

    @property
    def CutoffFrequency(self):
        m1 = self.Mp
        m2 = self.Mm
        q = m1 / m2
        M = m1 + m2
        chi1 = self.Spinp[2] / m1**2
        chi2 = self.Spinm[2] / m2**2
        LInitNR = self.xp*self.Pp[1] + self.xm*self.Pm[1] 
        # .014 is the initial guess for cutoff frequency
        omOrbPN = scipy.optimize.fsolve(angular_momentum, .014, (q, M, chi1, chi2, LInitNR))[0]
        omOrbPN = omOrbPN**(3./2.)
        omGWPN = 2. * omOrbPN
        return 0.75 * omGWPN