from ..Visualize.plot3d import plot_trajectory
import numpy as np

class volumeintegrals_grmhd:
    """
    A class representing the output of the VolumeIntegrals_GRMHD thorn in a simulation.
    """
    def __init__(self, files):
        dset = np.loadtxt(files[0], unpack=True, comments="#")
        self.t = dset[0]
        self.CoM_1 = [dset[2], dset[3], dset[4]]
        self.M_1 = dset[5]
        self.CoM_2 = [dset[7], dset[8], dset[9]]
        self.M_2 = dset[10]
        self.ADM = dset[12]

    def Preview(self, tstart=0, tend=-1, view_angle=None, axlim=None, axis=None):
        """
        Plot NS/BH tracks in initial simulation frame
        """
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D

        x1, y1, z1 = self.CoM_1
        x2, y2, z2 = self.CoM_2

        fig = plt.figure(figsize=(10,10))
        ax = fig.gca(projection='3d')

        plot_trajectory(ax, x1[tstart:tend], y1[tstart:tend], z1[tstart:tend], color='#0392ff', label=r'$NR_1$')
        plot_trajectory(ax, x2[tstart:tend], y2[tstart:tend], z1[tstart:tend], color='#ff1c03', label=r'$NR_2$')

        if view_angle:
            ax.view_init(view_angle[0],view_angle[1])
        if not axis:
            ax.axis('off')
        else:
            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            ax.set_zlabel('Z')
        ax.legend()
        if axlim is not None:
            ax.set_xlim(axlim)
            ax.set_ylim(axlim)
        plt.show()