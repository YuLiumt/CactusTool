"""
Makes a 2D color plot
"""
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
from matplotlib.patches import Rectangle

def pcolormesh(dsets, dim='xy', axlim=None, Normalize=None, title=None):
    vmax = -np.inf
    vmin = np.inf
    for rl in sorted(dsets):
        for c in sorted(dsets[rl]):
            if vmax < np.amax(dsets[rl][c]['data']):
                vmax = np.amax(dsets[rl][c]['data'])
            if vmin > np.amin(dsets[rl][c]['data']):
                vmin = np.amin(dsets[rl][c]['data'])

    for rl in sorted(dsets):
        for c in sorted(dsets[rl]):
            if Normalize is None:
                plt.pcolormesh(dsets[rl][c][dim[0]], dsets[rl][c][dim[1]], dsets[rl][c]['data'], vmin=vmin, vmax=vmax)
            elif Normalize == 'LogNorm':
                assert vmin > 0, "values must all be positive"
                plt.pcolormesh(dsets[rl][c][dim[0]], dsets[rl][c][dim[1]], dsets[rl][c]['data'], norm=colors.LogNorm(vmin=vmin, vmax=vmax))

    if axlim is not None:
        plt.xlim(axlim)
        plt.ylim(axlim)
    plt.xlabel(dim[0])
    plt.ylabel(dim[1])
    if title is not None:
        plt.title(title)
    plt.colorbar()
    plt.show()

    

    # def grid(self, facecolor=None, edgecolor=None):
    #     """
    #     Plots grid structure
    #     """
    #     if edgecolor is None:
    #         edgecolor=['k']
    #     if facecolor is None:
    #         facecolor=['w']
        # Rectangle(x0, width=dx[0], height=dx[1], edgecolor=eclr, facecolor=fclr, **kwargs)
        # self.ax.add_patch(patch)
    # def extent():
    #     if extent[0] is not None:
    #         self.ax.set_xlim(xmin=extent[0])
    #     if extent[1] is not None:
    #         self.ax.set_xlim(xmax=extent[1])
    #     if extent[2] is not None:
    #         self.ax.set_ylim(ymin=extent[2])
    #     if extent[3] is not None:
    #         self.ax.set_ylim(ymax=extent[3])