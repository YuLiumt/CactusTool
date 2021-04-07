"""
Makes a 2D color plot
"""
from ..utils.log import logger
from ..utils.units import UnitConversion
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
from matplotlib.patches import Rectangle

def pcolormesh(ax, dsets, scale=None, unit=('1', '1'), label=None, **kwargs):
    """
    Create a pseudocolor 2-D plot. 

    :param ax: Which axes to use for the plot.
    :param dict dsets: AMRgrid data
    :param str scale: log or not
    :param kwargs: Unknown keyword arguments are passed to :py:func:`ax.pcolormesh()`.
    :return: image object
    """
    if label:
        ax.set_xlabel(label[0])
        ax.set_ylabel(label[1])

    vmax = -np.inf
    vmin = np.inf
    for rl in sorted(dsets):
        for c in sorted(dsets[rl]):
            if vmax < np.amax(dsets[rl][c]['data']*UnitConversion(unit[1])):
                vmax = np.amax(dsets[rl][c]['data']*UnitConversion(unit[1]))
            if vmin > np.amin(dsets[rl][c]['data']*UnitConversion(unit[1])):
                vmin = np.amin(dsets[rl][c]['data']*UnitConversion(unit[1]))

    if vmin < 0 and scale == 'log':
        scale = None
        logger.info("values must all be positive. we reset scale as None")

    for rl in sorted(dsets):
        for c in sorted(dsets[rl]):
            origin = dsets[rl][c]['origin']
            delta = dsets[rl][c]['delta']
            mesh = dsets[rl][c]['data']
            size = mesh.shape
            x = np.arange(0, size[1]+1)*delta[0] + origin[0] - delta[0]/2
            y = np.arange(0, size[0]+1)*delta[1] + origin[1] - delta[1]/2
            if scale is None:
                norm = colors.Normalize(vmin=vmin, vmax=vmax)
            elif scale == 'log':
                mesh = np.abs(mesh)
                norm = colors.LogNorm(vmin=vmin, vmax=vmax)
            im = ax.pcolormesh(x*UnitConversion(unit[0]), y*UnitConversion(unit[0]), mesh*UnitConversion(unit[1]), norm=norm, **kwargs)

    return im

def contourf(ax, dsets, scale=None, levels=3, **kwargs):
    """
    Draw contour lines.

    :param ax: Which axes to use for the plot.
    :param dict dsets: AMRgrid data
    :param str scale: log or not
    :param kwargs: Unknown keyword arguments are passed to :py:func:`ax.pcolormesh()`.
    :return: image object
    """
    vmax = -np.inf
    vmin = np.inf
    for rl in sorted(dsets):
        for c in sorted(dsets[rl]):
            if vmax < np.amax(dsets[rl][c]['data']):
                vmax = np.amax(dsets[rl][c]['data'])
            if vmin > np.amin(dsets[rl][c]['data']):
                vmin = np.amin(dsets[rl][c]['data'])

    rl = max(rho.keys())
    for c in sorted(dsets[rl]):
        origin = dsets[rl][c]['origin']
        delta = dsets[rl][c]['delta']
        mesh = dsets[rl][c]['data']
        size = mesh.shape
        x = np.arange(0, size[1])*delta[0] + origin[0]
        y = np.arange(0, size[0])*delta[1] + origin[1]
        if scale is None:
            norm = colors.Normalize(vmin=vmin, vmax=vmax)
        elif scale == 'log':
            if vmin < 0:
                mesh = np.abs(mesh)
                logger.info("Use log_abs replace log")
            norm = colors.LogNorm(vmin=vmin, vmax=vmax)
        im = ax.contourf(x, y, mesh, norm=norm, levels=levels, **kwargs)

    return im

def grid(self, facecolor=None, edgecolor=None):
    """
    Plots grid structure
    """
    pass
    #     if edgecolor is None:
    #         edgecolor=['k']
    #     if facecolor is None:
    #         facecolor=['w']
    #     Rectangle(x0, width=dx[0], height=dx[1], edgecolor=eclr, facecolor=fclr, **kwargs)
    #     self.ax.add_patch(patch)
    # def extent():
    #     if extent[0] is not None:
    #         self.ax.set_xlim(xmin=extent[0])
    #     if extent[1] is not None:
    #         self.ax.set_xlim(xmax=extent[1])
    #     if extent[2] is not None:
    #         self.ax.set_ylim(ymin=extent[2])
    #     if extent[3] is not None:
    #         self.ax.set_ylim(ymax=extent[3])