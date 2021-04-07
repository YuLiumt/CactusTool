from ..Visualize.plot3d import plot_trajectory
from ..utils.units import UnitConversion
import numpy as np

class volumeintegrals_grmhd:
    """
    A class representing the output of the VolumeIntegrals_GRMHD thorn in a simulation.
    """
    def __init__(self, files):
        dset = np.loadtxt(files[0], unpack=True, comments="#")
        # Computational units CU to ms or km
        self.t = dset[0]
        self.x1 = dset[2]/dset[5]
        self.y1 = dset[3]/dset[5]
        self.z1 = dset[4]/dset[5]
        self.M_1 = dset[5]
        self.x2 = dset[7]/dset[10]
        self.y2 = dset[8]/dset[10] 
        self.z2 = dset[9]/dset[10]
        self.M_2 = dset[10]
        self.ADM = dset[12]

    def Preview(self, tstart=0, tend=-1, view_angle=None, axlim=None, axis=None, unit=('ms', 'km')):
        """
        Plot NS/BH tracks in initial simulation frame
        """
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D

        fig = plt.figure(figsize=(10,10))
        ax = fig.gca(projection='3d')

        plot_trajectory(ax, self.x1[tstart:tend] * UnitConversion(unit[1]), self.y1[tstart:tend] * UnitConversion(unit[1]), self.z1[tstart:tend] * UnitConversion(unit[1]), color='#0392ff', label=r'$Star_1$')
        plot_trajectory(ax, self.x2[tstart:tend] * UnitConversion(unit[1]), self.y2[tstart:tend] * UnitConversion(unit[1]), self.z2[tstart:tend] * UnitConversion(unit[1]), color='#ff1c03', label=r'$Star_2$')

        ax.set_title("From {:.2f} to {:.2f} [ms]".format(self.t[tstart] * UnitConversion(unit[0]), self.t[tend] * UnitConversion(unit[0])))

        if view_angle:
            ax.view_init(view_angle[0],view_angle[1])
        if not axis:
            ax.axis('off')
        else:
            ax.set_xlabel('X [km]')
            ax.set_ylabel('Y [km]')
            ax.set_zlabel('Z [km]')
        ax.legend()
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

    def eccentricity(self, p0=[0, 0, 0.01, 1, 0], show=True):
        """
        $$
        \dot{D}(t)=A_{0}+A_{1} t-e D_{0} \omega_{e} \sin \left(\omega_{e} t+\phi_{e}\right)
        $$
        where $e$ is the eccentricity and $D_{0}$ the initial coordinate interbinary distance.

        arXiv:1605.03424
        """
        from scipy.optimize import curve_fit

        dist = self.distance
        dist[0] = dist[0] * UnitConversion('ms')
        dist[1] = dist[1] * UnitConversion('km')

        def orbital_evolution(t, A0, A1, e, we, phie):
            return A0 + A1*t - e*dist[1][0]*we*np.sin(we*t + phie)

        try:
            tmerger = dist[0][np.where(dist[1] == 0.)][0]
        except:
            print("Make sure BNS have merged")

        # The fit is performed in the time interval between $t_{\mathrm{ret}}=3 \mathrm{ms}$ and $t_{\mathrm{ret}}=\frac{2}{3} t_{\mathrm{merger}}$ to avoid the initial spurious radiation and the plunge phase but having at least one eccentricity cycle included.
        index = np.where((dist[0] >= 3) & (dist[0] <= 2/3 * tmerger))[0]
        t = dist[0][index]
        y = np.gradient(dist[1], dist[0])[index]

        params, params_covariance = curve_fit(orbital_evolution, t, y, p0=p0)
        print('Orbital Eccentricity =', params[2])

        if show:
            import matplotlib.pyplot as plt

            fig, ax = plt.subplots()
            ax.plot(t, y)
            ax.set_xlabel('time [ms]')
            ax.set_ylabel(r'$\dot{D}$')
            ax.plot(t, orbital_evolution(t, params[0], params[1], params[2], params[3], params[4]), label='Fitted function')
            plt.legend(loc='best')
            ax.set_xlim(t.min(), t.max())
            plt.show()